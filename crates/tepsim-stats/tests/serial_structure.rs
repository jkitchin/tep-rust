//! Known-answer tests for the autocorrelation function and Welch spectra.
//!
//! The one that matters most is the *absolute level* of the spectrum. Tier 5
//! compares two spectra to each other, and a normalisation that is wrong by a
//! constant factor cancels in every such comparison and is therefore invisible
//! to the tests one would naturally write. So a unit-amplitude sinusoid is
//! required here to come back with the power it actually has.

#![allow(
    clippy::float_cmp,
    reason = "known-answer tests assert exact values where the answer is exact"
)]
#![allow(
    clippy::suboptimal_flops,
    reason = "fixtures mirror the formulas they stand in for"
)]

use tepsim_stats::serial::{autocorrelation_standard_error, decorrelation_lag};
use tepsim_stats::{Window, autocorrelation, band_comparison, log_band_edges, welch};
use tepsim_stats::{autocorrelation_direct, bartlett_standard_error};

/// A deterministic generator, so these tests need no random-number crate and
/// give the same answer on every platform.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed | 1)
    }
    /// Uniform on [-0.5, 0.5), variance 1/12.
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        ((self.0 >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
    }
}

// ---------------------------------------------------------------------------
// Autocorrelation
// ---------------------------------------------------------------------------

/// `x = [1, 2, 3, 4]`, worked out by hand.
///
/// Mean 2.5, deviations `[-1.5, -0.5, 0.5, 1.5]`, denominator `5`.
///
/// ```text
/// r(1) = (0.75 - 0.25 + 0.75) / 5 = 0.25
/// r(2) = (-0.75 - 0.75)      / 5 = -0.30
/// r(3) = (-2.25)             / 5 = -0.45
/// ```
#[test]
fn the_autocorrelation_matches_a_hand_computed_case() {
    let acf = autocorrelation(&[1.0, 2.0, 3.0, 4.0], 3);
    assert_eq!(acf.len(), 4);
    assert_eq!(acf[0], 1.0);
    assert!((acf[1] - 0.25).abs() < 1e-15, "r(1) = {}", acf[1]);
    assert!((acf[2] + 0.30).abs() < 1e-15, "r(2) = {}", acf[2]);
    assert!((acf[3] + 0.45).abs() < 1e-15, "r(3) = {}", acf[3]);
}

#[test]
fn lag_zero_is_exactly_one() {
    let mut rng = Lcg::new(0xACF);
    for n in [2_usize, 7, 100, 1000] {
        let x: Vec<f64> = (0..n).map(|_| 3.0 + 10.0 * rng.next()).collect();
        assert_eq!(autocorrelation(&x, 5)[0], 1.0, "n={n}");
    }
}

/// An AR(1) process has theoretical autocorrelation `phi^k`.
///
/// The tolerance is Bartlett's standard error, `1/sqrt(n)`, times a factor
/// stated here rather than tuned: three standard errors is the usual band and
/// is what this asks for.
#[test]
fn the_autocorrelation_of_an_ar1_process_is_phi_to_the_k() {
    let n = 400_000;
    println!(
        "n={n}, white-noise standard error {:.5}",
        autocorrelation_standard_error(n)
    );

    for phi in [0.5_f64, 0.8, 0.95] {
        let mut rng = Lcg::new(0xA51 + (phi * 100.0) as u64);
        let mut x = Vec::with_capacity(n);
        let mut previous = 0.0;
        // Burn in, so the series starts from the stationary distribution
        // rather than from zero.
        for _ in 0..10_000 {
            previous = phi * previous + rng.next();
        }
        for _ in 0..n {
            previous = phi * previous + rng.next();
            x.push(previous);
        }

        // Out to 400 lags, not 20: Bartlett's formula sums over the whole
        // ACF, and at phi = 0.95 the correlation is still 0.36 at lag 20, so a
        // short ACF understates the standard error badly.
        let acf = autocorrelation(&x, 400);
        let mut worst = (0.0_f64, 0);
        for (k, r) in acf.iter().enumerate().take(21) {
            let expected = phi.powi(k as i32);
            let deviations = if k == 0 {
                // r(0) is exactly 1 by construction, so it has no error and
                // Bartlett's formula gives zero there. Dividing would be 0/0.
                0.0
            } else {
                (r - expected).abs() / bartlett_standard_error(&acf, k, n)
            };
            if deviations > worst.0 {
                worst = (deviations, k);
            }
        }
        println!(
            "  phi={phi}: worst deviation {:.2} Bartlett standard errors, at              lag {} (where one standard error is {:.5}, against {:.5} for              white noise)",
            worst.0,
            worst.1,
            bartlett_standard_error(&acf, worst.1.max(1), n),
            autocorrelation_standard_error(n)
        );
        assert!(
            worst.0 < 3.0,
            "phi={phi}: r({}) is {:.2} Bartlett standard errors from phi^k",
            worst.1,
            worst.0
        );
    }
}

/// White noise decorrelates at lag 1, and the sample ACF stays inside the
/// Bartlett band.
#[test]
fn white_noise_has_no_autocorrelation() {
    let n = 200_000;
    let mut rng = Lcg::new(0x0741);
    let x: Vec<f64> = (0..n).map(|_| rng.next()).collect();
    let acf = autocorrelation(&x, 50);
    let band = 3.0 * autocorrelation_standard_error(n);

    let worst = acf[1..].iter().map(|r| r.abs()).fold(0.0_f64, f64::max);
    println!("white noise, n={n}: worst |r(k)| over lags 1..50 = {worst:.5}, band {band:.5}");
    assert!(worst < band, "worst {worst:.5} outside the band {band:.5}");
    assert_eq!(decorrelation_lag(&acf, 0.5), Some(1));
}

/// A pure sinusoid never decorrelates: its ACF oscillates forever.
///
/// This is the property that distinguishes a limit cycle from noise, and it is
/// the whole reason Tier 5 looks at serial structure at all.
#[test]
fn a_sinusoid_never_decorrelates() {
    let n = 100_000;
    let period = 50.0;
    let x: Vec<f64> = (0..n)
        .map(|t| (2.0 * std::f64::consts::PI * t as f64 / period).sin())
        .collect();
    let acf = autocorrelation(&x, 200);

    // At every whole period the ACF returns to nearly one.
    for multiple in 1..=4 {
        let lag = multiple * period as usize;
        assert!(
            acf[lag] > 0.99,
            "r({lag}) = {:.4}, but the signal has period {period}",
            acf[lag]
        );
    }
    // And at every half period it is nearly minus one.
    for multiple in 0..4 {
        let lag = multiple * period as usize + period as usize / 2;
        assert!(
            acf[lag] < -0.99,
            "r({lag}) = {:.4}, half a period out of phase",
            acf[lag]
        );
    }
    assert_eq!(
        decorrelation_lag(&acf, 0.5),
        Some(9),
        "the ACF of a period-50 sine first drops below 0.5 at lag 9, where \
         cos(2 pi 9 / 50) = 0.4258"
    );
}

/// The transform route against the definition it claims to compute.
///
/// The equality is exact in exact arithmetic; this is about the code, and
/// specifically about the zero padding. Padding to less than `2n` makes the
/// correlation circular, and the resulting error is largest for a series with
/// a trend, which is why one of the fixtures below has one.
#[test]
fn the_transform_autocorrelation_equals_the_definition() {
    let mut rng = Lcg::new(0x4C0F);
    let mut worst = 0.0_f64;

    let fixtures: Vec<(&str, Vec<f64>)> = vec![
        ("white noise", (0..5_000).map(|_| rng.next()).collect()),
        (
            "a strong trend, where a circular correlation would show",
            (0..5_000).map(|t| t as f64 * 0.01 + rng.next()).collect(),
        ),
        (
            "a slow sinusoid",
            (0..5_000)
                .map(|t| (2.0 * std::f64::consts::PI * t as f64 / 617.0).sin())
                .collect(),
        ),
        (
            "plant-shaped: a large offset and a small ripple",
            (0..5_000)
                .map(|t| 2705.0 + 0.1 * (t as f64 * 0.37).sin() + 0.01 * rng.next())
                .collect(),
        ),
    ];

    for (what, series) in &fixtures {
        let fast = autocorrelation(series, 300);
        let slow = autocorrelation_direct(series, 300);
        assert_eq!(fast.len(), slow.len(), "{what}");
        let mut here = 0.0_f64;
        for (lag, (a, b)) in fast.iter().zip(&slow).enumerate() {
            let error = (a - b).abs();
            if error > here {
                here = error;
            }
            assert!(error < 1e-11, "{what}, lag {lag}: {a} against {b}");
        }
        println!("  {what}: worst absolute difference {here:.3e}");
        if here > worst {
            worst = here;
        }
    }
    println!("transform vs definition: worst absolute difference {worst:.3e}");
}

/// The cost, measured, because B-0047 has to budget for it.
#[test]
fn the_transform_autocorrelation_is_faster_than_the_definition() {
    // One Tier 5 series: 48 h at one second.
    let mut rng = Lcg::new(0xC057);
    let x: Vec<f64> = (0..172_800).map(|_| rng.next()).collect();

    let start = std::time::Instant::now();
    let fast = autocorrelation(&x, 200);
    let fast_time = start.elapsed();

    let start = std::time::Instant::now();
    let slow = autocorrelation_direct(&x, 200);
    let slow_time = start.elapsed();

    println!(
        "172,800 samples to lag 200: transform {fast_time:?}, definition {slow_time:?}, ratio {:.2}",
        slow_time.as_secs_f64() / fast_time.as_secs_f64()
    );
    for (lag, (a, b)) in fast.iter().zip(&slow).enumerate() {
        assert!((a - b).abs() < 1e-11, "lag {lag}: {a} against {b}");
    }

    // The claim is about *scaling*, not about the wall clock, and the scaling
    // is what is asserted: the transform's cost does not depend on how many
    // lags are asked for, and the direct method's is linear in them.
    //
    // An earlier version asserted `fast_time < slow_time` directly and failed
    // intermittently. It was not measuring the code. At this size the two are
    // within a factor of two of each other, so on a machine busy with other
    // work the ordering flips, and the test becomes a measurement of the load
    // average.
    let start = std::time::Instant::now();
    let _ = autocorrelation(&x, 2_000);
    let fast_deep = start.elapsed();
    let start = std::time::Instant::now();
    let _ = autocorrelation_direct(&x, 2_000);
    let slow_deep = start.elapsed();

    let fast_growth = fast_deep.as_secs_f64() / fast_time.as_secs_f64();
    let slow_growth = slow_deep.as_secs_f64() / slow_time.as_secs_f64();
    println!(
        "  going from 200 lags to 2,000: transform x{fast_growth:.2}, definition x{slow_growth:.2}"
    );
    assert!(
        slow_growth > 3.0 * fast_growth,
        "asking for ten times the lags cost the transform x{fast_growth:.2} and the definition x{slow_growth:.2}, so the transform is not the lag-independent one"
    );
}

#[test]
fn a_degenerate_series_returns_nothing_rather_than_inventing_a_number() {
    assert!(autocorrelation(&[], 5).is_empty());
    assert!(autocorrelation(&[1.0], 5).is_empty());
    // Zero variance: the ratio is 0/0.
    assert!(autocorrelation(&[7.0, 7.0, 7.0, 7.0], 2).is_empty());
}

// ---------------------------------------------------------------------------
// Welch: the absolute level
// ---------------------------------------------------------------------------

/// **The test the rest of Tier 5 depends on.**
///
/// A sinusoid of amplitude `A` carries power `A^2 / 2`. Integrating the
/// one-sided density over all frequency must return exactly that, whatever
/// window is used, whatever the segment length, and whatever the sample rate.
///
/// Nothing that compares two spectra to each other can see an error here: a
/// normalisation wrong by a constant cancels in every ratio.
#[test]
fn a_sinusoid_carries_the_power_it_actually_has() {
    let sample_rate = 100.0;
    let segment = 1024;
    let n = 100 * segment;

    for amplitude in [1.0_f64, 0.25, 7.5] {
        for window in [Window::Rectangular, Window::Hann] {
            // Exactly on a bin centre, so there is no scalloping loss to
            // confuse the level with.
            let bin = 64.0;
            let frequency = bin * sample_rate / segment as f64;
            let signal: Vec<f64> = (0..n)
                .map(|t| {
                    amplitude
                        * (2.0 * std::f64::consts::PI * frequency * t as f64 / sample_rate).cos()
                })
                .collect();

            let spectrum = welch(&signal, sample_rate, segment, segment / 2, window);
            let power = spectrum.total_power();
            let expected = amplitude * amplitude / 2.0;
            let error = (power - expected).abs() / expected;
            println!(
                "A={amplitude}, {window:?}: total power {power:.9}, expected \
                 {expected:.9}, error {error:.2e}"
            );
            assert!(
                error < 1e-9,
                "A={amplitude}, {window:?}: total power {power}, expected {expected}"
            );
        }
    }
}

/// The same claim for noise: the integral of the density is the variance.
#[test]
fn the_integrated_density_is_the_variance() {
    let sample_rate = 1.0;
    let segment = 512;
    let mut rng = Lcg::new(0x5DE1);
    let x: Vec<f64> = (0..200 * segment).map(|_| rng.next()).collect();

    // Uniform on [-0.5, 0.5) has variance 1/12.
    let expected = 1.0 / 12.0;
    for window in [Window::Rectangular, Window::Hann] {
        let spectrum = welch(&x, sample_rate, segment, segment / 2, window);
        let power = spectrum.total_power();
        let error = (power - expected).abs() / expected;
        println!(
            "{window:?}: integrated density {power:.6}, variance {expected:.6}, \
             error {:.2}%",
            error * 100.0
        );
        // Sampling error on a Welch estimate of total power over `segments`
        // averages is about 1/sqrt(segments); with 50% overlap and 200
        // segments' worth of data that is a few percent.
        assert!(error < 0.02, "{window:?}: {power} against {expected}");
    }
}

/// A tone lands in the bin it belongs to, and the Hann window keeps its
/// neighbours far quieter than a rectangular one does.
#[test]
fn a_tone_lands_in_the_right_bin_and_hann_contains_its_leakage() {
    let sample_rate = 256.0;
    let segment = 256;
    let n = 64 * segment;
    // Deliberately *between* bins, which is where leakage is worst and where
    // the choice of window earns its keep.
    let bin = 40.5;
    let frequency = bin * sample_rate / segment as f64;
    let signal: Vec<f64> = (0..n)
        .map(|t| (2.0 * std::f64::consts::PI * frequency * t as f64 / sample_rate).cos())
        .collect();

    let mut leakage = Vec::new();
    for window in [Window::Rectangular, Window::Hann] {
        let spectrum = welch(&signal, sample_rate, segment, segment / 2, window);
        let peak = spectrum
            .density
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .expect("a peak");
        assert!(
            peak.0 == 40 || peak.0 == 41,
            "{window:?}: peak at bin {} for a tone at bin {bin}",
            peak.0
        );
        // Power more than five bins from the tone: pure leakage.
        let far: f64 = spectrum
            .density
            .iter()
            .enumerate()
            .filter(|(k, _)| k.abs_diff(40) > 5)
            .map(|(_, p)| p * spectrum.resolution)
            .sum();
        let fraction = far / spectrum.total_power();
        println!(
            "{window:?}: leakage beyond 5 bins = {:.4}%",
            fraction * 100.0
        );
        leakage.push(fraction);
    }
    assert!(
        leakage[1] < leakage[0] / 100.0,
        "Hann leaked {:.3e} against rectangular {:.3e}, less than the hundredfold \
         improvement that is the reason for using it",
        leakage[1],
        leakage[0]
    );
}

/// Removing each segment's mean is what makes the spectrum usable on plant
/// data, and it is not optional.
#[test]
fn the_segment_mean_is_removed() {
    let sample_rate = 1.0;
    let segment = 256;
    // Reactor pressure: a huge offset and a tiny oscillation.
    let signal: Vec<f64> = (0..40 * segment)
        .map(|t| {
            2705.0 + 0.4 * (2.0 * std::f64::consts::PI * 16.0 * t as f64 / segment as f64).sin()
        })
        .collect();

    let spectrum = welch(&signal, sample_rate, segment, segment / 2, Window::Hann);
    // The oscillation carries 0.4^2/2 = 0.08. The offset carries 7.3e6 and
    // must not appear at all.
    let power = spectrum.total_power();
    println!("offset 2705, oscillation 0.4: total power {power:.9}, expected 0.08");
    assert!(
        (power - 0.08).abs() / 0.08 < 1e-9,
        "total power {power}, expected 0.08. If this is near 7.3e6 the segment \
         mean is not being removed."
    );
    assert!(
        spectrum.density[0] < 1e-6,
        "the DC bin holds {:.3e}",
        spectrum.density[0]
    );
}

#[test]
fn the_spectrum_reports_its_own_shape() {
    let sample_rate = 3600.0;
    let segment = 4096;
    let signal = vec![0.0; 172_800];
    let spectrum = welch(&signal, sample_rate, segment, segment / 2, Window::Hann);

    assert_eq!(spectrum.frequencies.len(), segment / 2 + 1);
    assert_eq!(spectrum.resolution, sample_rate / segment as f64);
    assert_eq!(spectrum.frequencies[0], 0.0);
    assert_eq!(
        *spectrum.frequencies.last().expect("bins"),
        sample_rate / 2.0,
        "the last bin is Nyquist"
    );
    // 172800 samples, 4096-sample segments, 2048 step: floor((172800-4096)/2048)+1
    assert_eq!(spectrum.segments, 83);
    println!(
        "48 h at 1 s, 4096-sample Hann segments: {} averages, {:.4} mHz resolution",
        spectrum.segments,
        spectrum.resolution * 1000.0
    );
}

// ---------------------------------------------------------------------------
// Band comparison
// ---------------------------------------------------------------------------

#[test]
fn identical_spectra_compare_as_one_in_every_band() {
    let sample_rate = 100.0;
    let segment = 512;
    let mut rng = Lcg::new(0xBA4D);
    let x: Vec<f64> = (0..50 * segment).map(|_| rng.next()).collect();
    let spectrum = welch(&x, sample_rate, segment, segment / 2, Window::Hann);

    let edges = log_band_edges(spectrum.resolution, sample_rate / 2.0, 6);
    let bands = band_comparison(&spectrum, &spectrum, &edges);
    assert_eq!(bands.len(), 6);
    for band in &bands {
        println!("{band}");
        assert_eq!(band.ratio, 1.0, "{band}");
    }
    // The bands cover everything except the DC and Nyquist bins: the lowest
    // edge is bin 1, and `band_power` is half-open so a top edge exactly at
    // Nyquist excludes the Nyquist bin. Both exclusions are deliberate and
    // both are asserted, because a set of bands that quietly dropped a bin
    // would make every ratio in Tier 5 slightly wrong.
    let covered: f64 = bands.iter().map(|b| b.power_b).sum();
    let ends =
        (spectrum.density[0] + spectrum.density[spectrum.density.len() - 1]) * spectrum.resolution;
    let total = spectrum.total_power() - ends;
    assert!(
        (covered - total).abs() / total < 1e-12,
        "bands cover {covered:.6e} of {total:.6e}"
    );

    // And extending the top edge past Nyquist does pick it up.
    let mut wider = edges.clone();
    let last = wider.len() - 1;
    wider[last] = sample_rate;
    let all = band_comparison(&spectrum, &spectrum, &wider);
    let covered_all: f64 = all.iter().map(|b| b.power_b).sum();
    let expected = spectrum.total_power() - spectrum.density[0] * spectrum.resolution;
    assert!(
        (covered_all - expected).abs() / expected < 1e-12,
        "with the top edge above Nyquist, bands cover {covered_all:.6e} of          {expected:.6e}"
    );
}

/// A difference confined to one band shows up in that band and nowhere else.
/// This is the whole argument for comparing band by band rather than with one
/// norm.
#[test]
fn a_difference_in_one_band_is_localised_to_it() {
    let sample_rate = 256.0;
    let segment = 512;
    let n = 60 * segment;
    let mut rng = Lcg::new(0x10CA1);
    let noise: Vec<f64> = (0..n).map(|_| rng.next()).collect();

    // The same noise, plus a tone at 32 Hz in one of them.
    let tone = 32.0;
    let with_tone: Vec<f64> = noise
        .iter()
        .enumerate()
        .map(|(t, v)| v + 0.3 * (2.0 * std::f64::consts::PI * tone * t as f64 / sample_rate).sin())
        .collect();

    let a = welch(&with_tone, sample_rate, segment, segment / 2, Window::Hann);
    let b = welch(&noise, sample_rate, segment, segment / 2, Window::Hann);

    let edges = log_band_edges(a.resolution, sample_rate / 2.0, 6);
    let bands = band_comparison(&a, &b, &edges);
    let mut guilty = 0;
    for (index, band) in bands.iter().enumerate() {
        println!("{band}");
        if band.low <= tone && tone < band.high {
            guilty = index;
        }
    }
    assert!(
        bands[guilty].ratio > 1.5,
        "the band containing the tone has ratio {:.4}",
        bands[guilty].ratio
    );
    for (index, band) in bands.iter().enumerate() {
        if index != guilty {
            assert!(
                (band.ratio - 1.0).abs() < 0.15,
                "band {index} moved to {:.4} although the tone is not in it",
                band.ratio
            );
        }
    }
}

#[test]
fn log_band_edges_span_the_range_and_ascend() {
    let edges = log_band_edges(0.01, 100.0, 8);
    assert_eq!(edges.len(), 9);
    assert!((edges[0] - 0.01).abs() < 1e-15);
    assert!((edges[8] - 100.0).abs() < 1e-12);
    for pair in edges.windows(2) {
        assert!(pair[1] > pair[0]);
        // Constant ratio, which is what "logarithmic" means.
        assert!((pair[1] / pair[0] - edges[1] / edges[0]).abs() < 1e-12);
    }
}

#[test]
#[should_panic(expected = "not a power of two")]
fn a_non_power_of_two_segment_is_rejected() {
    let _ = welch(&[0.0; 1000], 1.0, 100, 50, Window::Hann);
}

#[test]
#[should_panic(expected = "shorter than one")]
fn a_signal_shorter_than_a_segment_is_rejected() {
    let _ = welch(&[0.0; 100], 1.0, 256, 128, Window::Hann);
}
