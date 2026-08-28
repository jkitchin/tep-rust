//! Known-answer tests for the FFT.
//!
//! The trap this file is built around: **a round-trip test and Parseval's
//! identity both pass for a transform whose output is permuted.** Bit-reversal
//! is the easiest part of an FFT to get wrong and the hardest to notice, and
//! those two tests, which are the ones everyone writes, are exactly blind to
//! it. So the transforms below are pinned against hand-computed values and
//! against the definition, and the round trip and Parseval are here as
//! supporting evidence rather than as the case.

#![allow(
    clippy::float_cmp,
    reason = "known-answer tests assert exact values where the answer is exact"
)]
#![allow(
    clippy::suboptimal_flops,
    reason = "fixtures mirror the formulas they stand in for"
)]

use tepsim_stats::{Complex, Fft, dft_naive};

/// How close a transform has to be to its hand-computed value.
///
/// The butterflies accumulate `log2(n)` roundings, so at n = 4096 that is
/// twelve, about 2.6e-15 relative. This asks for 1e-14, which is four times
/// that and thirty times the worst actually measured. It is deliberately not
/// the 1e-13 first written here: once the reference was fixed (see
/// `dft_naive`) the agreement turned out to be sub-ULP, and a tolerance two
/// orders above the truth would let a real regression through unnoticed.
const TIGHT: f64 = 1e-14;

fn assert_close(actual: &[Complex], expected: &[Complex], tolerance: f64, what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: length");
    let scale = expected
        .iter()
        .map(|z| z.abs())
        .fold(0.0_f64, f64::max)
        .max(1.0);
    for (index, (a, e)) in actual.iter().zip(expected).enumerate() {
        let error = (*a - *e).abs() / scale;
        assert!(
            error <= tolerance,
            "{what}: bin {index} is {a:?}, expected {e:?}, error {error:.3e}"
        );
    }
}

fn real(values: &[f64]) -> Vec<Complex> {
    values.iter().copied().map(Complex::real).collect()
}

// ---------------------------------------------------------------------------
// Hand-computed transforms
// ---------------------------------------------------------------------------

/// `n = 1` is the identity, and `n = 2` is `[x0 + x1, x0 - x1]`.
#[test]
fn the_smallest_transforms_are_exact() {
    let fft = Fft::new(1);
    let mut buffer = real(&[7.0]);
    fft.forward(&mut buffer);
    assert_eq!(buffer, vec![Complex::new(7.0, 0.0)]);

    let fft = Fft::new(2);
    for (x0, x1) in [(1.0, 0.0), (0.0, 1.0), (3.0, -5.0), (2.5, 2.5)] {
        let mut buffer = real(&[x0, x1]);
        fft.forward(&mut buffer);
        assert_eq!(
            buffer,
            vec![Complex::new(x0 + x1, 0.0), Complex::new(x0 - x1, 0.0)],
            "n=2 on [{x0}, {x1}]"
        );
    }
}

/// `n = 4`, four cases worked out on paper.
///
/// The twiddles are `1, -i, -1, i`, so every entry is exact in `f64` and the
/// answers below are exact too.
///
/// | input          | transform        | why |
/// |----------------|------------------|-----|
/// | `[1,0,0,0]`    | `[1,1,1,1]`      | an impulse is flat |
/// | `[1,1,1,1]`    | `[4,0,0,0]`      | a constant is a single DC bin |
/// | `[1,0,-1,0]`   | `[0,2,0,2]`      | `1 - (-1)^k` |
/// | `[0,1,0,-1]`   | `[0,-2i,0,2i]`   | the same, shifted by one sample |
#[test]
fn the_four_point_transform_is_exact() {
    let fft = Fft::new(4);

    let cases: &[(&[f64], &[Complex])] = &[
        (
            &[1.0, 0.0, 0.0, 0.0],
            &[
                Complex::new(1.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(1.0, 0.0),
                Complex::new(1.0, 0.0),
            ],
        ),
        (
            &[1.0, 1.0, 1.0, 1.0],
            &[
                Complex::new(4.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 0.0),
            ],
        ),
        (
            &[1.0, 0.0, -1.0, 0.0],
            &[
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(2.0, 0.0),
            ],
        ),
        (
            &[0.0, 1.0, 0.0, -1.0],
            &[
                Complex::new(0.0, 0.0),
                Complex::new(0.0, -2.0),
                Complex::new(0.0, 0.0),
                Complex::new(0.0, 2.0),
            ],
        ),
    ];

    for (input, expected) in cases {
        let mut buffer = real(input);
        fft.forward(&mut buffer);
        assert_close(&buffer, expected, TIGHT, &format!("n=4 on {input:?}"));
    }
}

/// An impulse at sample `d` transforms to `exp(-2 pi i d k / n)`, whose
/// magnitude is one in every bin and whose *phase* winds `d` times around.
///
/// This is the test a permuted output cannot survive: the phase in bin `k`
/// names `k`, so any reordering is visible immediately. A flat magnitude alone
/// would not do it.
#[test]
fn a_shifted_impulse_winds_its_phase_the_right_number_of_times() {
    for size in [8_usize, 16, 64] {
        let fft = Fft::new(size);
        for delay in [0_usize, 1, 3, size / 2, size - 1] {
            let mut signal = vec![0.0; size];
            signal[delay] = 1.0;
            let mut buffer = real(&signal);
            fft.forward(&mut buffer);

            // `(delay * k) % size` before the trig, for the same reason
            // `dft_naive` does it: at size = 64 the unreduced argument reaches
            // 389 radians, which carries 8.6e-14 of error in the argument
            // alone and would make this fixture, not the transform, the least
            // accurate thing in the comparison.
            let expected: Vec<Complex> = (0..size)
                .map(|k| {
                    let reduced = (delay * k) % size;
                    let angle = -2.0 * std::f64::consts::PI * reduced as f64 / size as f64;
                    Complex::new(angle.cos(), angle.sin())
                })
                .collect();
            assert_close(
                &buffer,
                &expected,
                TIGHT,
                &format!("n={size}, impulse at {delay}"),
            );
        }
    }
}

/// A pure cosine at bin `k0` puts `n/2` in bins `k0` and `n - k0` and nothing
/// anywhere else.
#[test]
fn a_pure_cosine_lands_in_exactly_two_bins() {
    let size = 64;
    let fft = Fft::new(size);
    for k0 in [1_usize, 5, 17, 31] {
        let signal: Vec<f64> = (0..size)
            .map(|j| (2.0 * std::f64::consts::PI * k0 as f64 * j as f64 / size as f64).cos())
            .collect();
        let spectrum = fft.forward_real(&signal);

        for (bin, value) in spectrum.iter().enumerate() {
            let expected = if bin == k0 || bin == size - k0 {
                size as f64 / 2.0
            } else {
                0.0
            };
            assert!(
                (value.abs() - expected).abs() < 1e-12 * size as f64,
                "k0={k0}: bin {bin} has magnitude {}, expected {expected}",
                value.abs()
            );
        }
        // And it is real: a cosine has no imaginary part in its spectrum.
        assert!(
            spectrum[k0].im.abs() < 1e-13 * size as f64,
            "k0={k0}: bin {k0} has imaginary part {}",
            spectrum[k0].im
        );
    }
}

/// A sine puts the same magnitude in the same two bins, with the phase a
/// quarter turn away. Together with the cosine test this pins the sign
/// convention, which a magnitude-only check leaves free.
#[test]
fn a_pure_sine_is_imaginary_and_antisymmetric() {
    let size = 32;
    let fft = Fft::new(size);
    let k0 = 5;
    let signal: Vec<f64> = (0..size)
        .map(|j| (2.0 * std::f64::consts::PI * k0 as f64 * j as f64 / size as f64).sin())
        .collect();
    let spectrum = fft.forward_real(&signal);

    // X[k0] = -i n/2 with the exp(-i...) convention; X[n-k0] = +i n/2.
    assert!(
        (spectrum[k0] - Complex::new(0.0, -(size as f64) / 2.0)).abs() < 1e-12 * size as f64,
        "bin {k0} is {:?}",
        spectrum[k0]
    );
    assert!(
        (spectrum[size - k0] - Complex::new(0.0, size as f64 / 2.0)).abs() < 1e-12 * size as f64,
        "bin {} is {:?}",
        size - k0,
        spectrum[size - k0]
    );
}

// ---------------------------------------------------------------------------
// Against the definition
// ---------------------------------------------------------------------------

/// A deterministic scatter, so this needs no random-number crate.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed | 1)
    }
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((self.0 >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
    }
}

#[test]
fn the_fast_transform_equals_the_definition() {
    let mut rng = Lcg::new(0xF17);
    let mut worst = 0.0_f64;
    for size in [2_usize, 4, 8, 16, 32, 64, 128, 256, 512, 1024] {
        let fft = Fft::new(size);
        let input: Vec<Complex> = (0..size)
            .map(|_| Complex::new(rng.next(), rng.next()))
            .collect();
        let mut fast = input.clone();
        fft.forward(&mut fast);
        let slow = dft_naive(&input);

        let scale = slow.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
        for (bin, (a, b)) in fast.iter().zip(&slow).enumerate() {
            let error = (*a - *b).abs() / scale;
            if error > worst {
                worst = error;
            }
            assert!(
                error < 1e-14,
                "n={size} bin {bin}: fast {a:?}, definition {b:?}, error {error:.3e}"
            );
        }
    }
    println!("fast vs definition, n up to 1024: worst relative {worst:.3e}");
}

/// Complex input, not just real: the imaginary path is exercised too.
#[test]
fn the_transform_is_linear() {
    let mut rng = Lcg::new(0x11EA5);
    let size = 128;
    let fft = Fft::new(size);
    let a: Vec<Complex> = (0..size)
        .map(|_| Complex::new(rng.next(), rng.next()))
        .collect();
    let b: Vec<Complex> = (0..size)
        .map(|_| Complex::new(rng.next(), rng.next()))
        .collect();
    let alpha = Complex::new(1.7, -0.4);

    let combined: Vec<Complex> = a.iter().zip(&b).map(|(x, y)| alpha * *x + *y).collect();

    let mut ta = a.clone();
    let mut tb = b.clone();
    let mut tc = combined;
    fft.forward(&mut ta);
    fft.forward(&mut tb);
    fft.forward(&mut tc);

    let expected: Vec<Complex> = ta.iter().zip(&tb).map(|(x, y)| alpha * *x + *y).collect();
    assert_close(&tc, &expected, 1e-14, "linearity");
}

// ---------------------------------------------------------------------------
// Supporting evidence: Parseval and the round trip
// ---------------------------------------------------------------------------

/// Parseval's identity, `sum |x|^2 = (1/n) sum |X|^2`.
///
/// Necessary and nowhere near sufficient: it holds for any unitary map, so a
/// permuted output passes it. It is here because Welch's method divides by the
/// window's power and an error in the transform's overall scale would land
/// straight in the spectrum.
#[test]
fn parseval_holds() {
    let mut rng = Lcg::new(0x9A23);
    for size in [4_usize, 32, 256, 2048] {
        let fft = Fft::new(size);
        let input: Vec<Complex> = (0..size)
            .map(|_| Complex::new(rng.next(), rng.next()))
            .collect();
        let time: f64 = input.iter().map(|z| z.norm_squared()).sum();

        let mut spectrum = input;
        fft.forward(&mut spectrum);
        let frequency: f64 = spectrum.iter().map(|z| z.norm_squared()).sum::<f64>() / size as f64;

        let error = (time - frequency).abs() / time;
        assert!(
            error < 1e-14,
            "n={size}: time domain {time:.17e}, frequency domain {frequency:.17e}"
        );
    }
}

/// The round trip. Also necessary and not sufficient, for the same reason: a
/// permutation composed with its inverse is the identity.
#[test]
fn the_inverse_undoes_the_forward_transform() {
    let mut rng = Lcg::new(0x5017);
    let mut worst = 0.0_f64;
    for size in [1_usize, 2, 8, 64, 512, 4096] {
        let fft = Fft::new(size);
        let input: Vec<Complex> = (0..size)
            .map(|_| Complex::new(rng.next(), rng.next()))
            .collect();
        let mut buffer = input.clone();
        fft.forward(&mut buffer);
        fft.inverse(&mut buffer);

        let scale = input.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
        for (index, (a, b)) in buffer.iter().zip(&input).enumerate() {
            let error = (*a - *b).abs() / scale;
            if error > worst {
                worst = error;
            }
            assert!(error < 1e-14, "n={size}, index {index}: {a:?} vs {b:?}");
        }
    }
    println!("round trip, n up to 4096: worst relative {worst:.3e}");
}

/// The size Welch will use, end to end, so B-0045 starts from a measured
/// accuracy rather than an assumed one.
#[test]
fn the_welch_segment_length_is_accurate_enough() {
    let size = 4096;
    let fft = Fft::new(size);
    // Two tones plus a constant, which is roughly the shape of a plant
    // measurement: a mean, a slow oscillation, and a faster one.
    let signal: Vec<f64> = (0..size)
        .map(|j| {
            let t = j as f64;
            2705.0
                + 0.4 * (2.0 * std::f64::consts::PI * 7.0 * t / size as f64).cos()
                + 0.1 * (2.0 * std::f64::consts::PI * 300.0 * t / size as f64).cos()
        })
        .collect();

    let fast = fft.forward_real(&signal);
    let slow = dft_naive(&real(&signal));
    let scale = slow.iter().map(|z| z.abs()).fold(0.0_f64, f64::max);
    let worst = fast
        .iter()
        .zip(&slow)
        .map(|(a, b)| (*a - *b).abs() / scale)
        .fold(0.0_f64, f64::max);
    println!("n=4096 plant-shaped: worst relative against the definition {worst:.3e}");
    assert!(worst < 1e-15, "worst {worst:.3e}");

    // The three components land where they should.
    assert!(
        (fast[0].abs() / size as f64 - 2705.0).abs() < 1e-9,
        "DC bin is {}",
        fast[0].abs() / size as f64
    );
    for (bin, amplitude) in [(7_usize, 0.4), (300, 0.1)] {
        let found = 2.0 * fast[bin].abs() / size as f64;
        assert!(
            (found - amplitude).abs() < 1e-11,
            "bin {bin}: amplitude {found}, expected {amplitude}"
        );
    }
}

/// The transform is bit-for-bit reproducible, which Tier 9 will rest on.
#[test]
fn the_transform_is_deterministic() {
    let size = 1024;
    let fft = Fft::new(size);
    let mut rng = Lcg::new(0xDEDE);
    let input: Vec<Complex> = (0..size)
        .map(|_| Complex::new(rng.next(), rng.next()))
        .collect();

    let mut first = input.clone();
    fft.forward(&mut first);
    // A second `Fft` of the same size, built independently.
    let mut second = input;
    Fft::new(size).forward(&mut second);

    for (index, (a, b)) in first.iter().zip(&second).enumerate() {
        assert_eq!(a.re.to_bits(), b.re.to_bits(), "bin {index} real part");
        assert_eq!(a.im.to_bits(), b.im.to_bits(), "bin {index} imaginary part");
    }
}

#[test]
#[should_panic(expected = "not a power of two")]
fn a_non_power_of_two_length_is_rejected_loudly() {
    let _ = Fft::new(100);
}
