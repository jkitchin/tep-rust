//! A radix-2 Cooley-Tukey FFT, and the naive transform it is checked against.
//!
//! Needed for Welch power spectra (B-0045), which is the Tier 5 test that
//! detects "the Fortran limit-cycles and the port damps". No moment or
//! distribution comparison sees that failure: the two runs can have the same
//! mean, the same variance and the same marginal distribution while one
//! oscillates and the other does not.
//!
//! # Determinism
//!
//! The same rule as everywhere else in this project. The butterfly order is
//! fixed by the algorithm, there is no threading and no runtime planning, and
//! the twiddle factors come from the vendored `libm` rather than the platform
//! one. Two runs on two architectures produce identical bits.
//!
//! Each twiddle is computed as `exp(-2 pi i m / n)` directly rather than by the
//! recurrence `w_{m+1} = w_m * w_1`. The recurrence is faster and accumulates
//! error along the table, which at n = 4096 is visible in the last few digits.
//!
//! # Sizes
//!
//! Powers of two only. That is not a limitation in practice here: Welch chooses
//! its own segment length, and 4096 samples at one second gives 83 averages
//! over a 48-hour run with a resolution of 0.244 mHz. Nothing in this project
//! needs Bluestein.

use alloc::vec;
use alloc::vec::Vec;

/// A complex number.
///
/// Defined here rather than pulled in, to keep this crate's dependency list at
/// exactly one (`libm`, for determinism) and because only three operations are
/// needed.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex {
    /// From real and imaginary parts.
    #[must_use]
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// A real number.
    #[must_use]
    pub const fn real(re: f64) -> Self {
        Self { re, im: 0.0 }
    }

    /// The squared magnitude, `|z|^2`.
    ///
    /// The square rather than the magnitude, because a power spectrum wants
    /// exactly this and taking a square root only to square it again loses
    /// half an ulp for nothing.
    #[must_use]
    pub fn norm_squared(self) -> f64 {
        self.re * self.re + self.im * self.im
    }

    /// The magnitude.
    #[must_use]
    pub fn abs(self) -> f64 {
        libm::hypot(self.re, self.im)
    }

    /// The complex conjugate.
    #[must_use]
    pub const fn conjugate(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }
}

impl core::ops::Add for Complex {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self::new(self.re + other.re, self.im + other.im)
    }
}

impl core::ops::Sub for Complex {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self::new(self.re - other.re, self.im - other.im)
    }
}

impl core::ops::Mul for Complex {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        // Fixed evaluation order; see the module docs on determinism.
        Self::new(
            self.re * other.re - self.im * other.im,
            self.re * other.im + self.im * other.re,
        )
    }
}

/// A transform of one fixed size, with its twiddle table.
///
/// Built once and reused. Welch runs the same length dozens of times per
/// series, and rebuilding the table each time would be `n/2` `sin` and `cos`
/// calls per segment for no reason.
#[derive(Clone, Debug)]
pub struct Fft {
    size: usize,
    /// `exp(-2 pi i m / n)` for `m` in `0..n/2`.
    twiddles: Vec<Complex>,
}

impl Fft {
    /// A transform of length `size`, which must be a power of two.
    ///
    /// # Panics
    ///
    /// If `size` is zero or not a power of two. This is a programming error
    /// rather than a data-dependent condition, so it is caught loudly at
    /// construction instead of returning a `NaN` spectrum later.
    #[must_use]
    pub fn new(size: usize) -> Self {
        assert!(
            size.is_power_of_two(),
            "FFT length {size} is not a power of two"
        );
        let half = size / 2;
        let mut twiddles = Vec::with_capacity(half);
        for m in 0..half {
            // Computed directly rather than by recurrence; see the module docs.
            let angle = -2.0 * core::f64::consts::PI * m as f64 / size as f64;
            twiddles.push(Complex::new(libm::cos(angle), libm::sin(angle)));
        }
        Self { size, twiddles }
    }

    /// The transform length.
    #[must_use]
    pub const fn size(&self) -> usize {
        self.size
    }

    /// The forward transform, in place.
    ///
    /// ```text
    /// X[k] = sum_j x[j] exp(-2 pi i j k / n)
    /// ```
    ///
    /// Unnormalised, which is the usual convention: the `1/n` lives in
    /// [`Fft::inverse`].
    ///
    /// # Panics
    ///
    /// If `data.len()` is not this transform's size.
    pub fn forward(&self, data: &mut [Complex]) {
        assert_eq!(
            data.len(),
            self.size,
            "buffer length does not match the transform"
        );
        bit_reverse_permute(data);

        let n = self.size;
        let mut span = 2;
        while span <= n {
            let half = span / 2;
            // The twiddle table is indexed for the full length, so a stage of
            // width `span` steps through it by `n / span`.
            let stride = n / span;
            let mut start = 0;
            while start < n {
                for k in 0..half {
                    let w = self.twiddles[k * stride];
                    let upper = data[start + k + half] * w;
                    let lower = data[start + k];
                    data[start + k] = lower + upper;
                    data[start + k + half] = lower - upper;
                }
                start += span;
            }
            span *= 2;
        }
    }

    /// The inverse transform, in place, normalised by `1/n`.
    ///
    /// Implemented by conjugation: `ifft(x) = conj(fft(conj(x))) / n`. That
    /// reuses the forward pass exactly, so any error in it shows up in the
    /// round trip rather than cancelling against a second, differently-wrong
    /// implementation.
    ///
    /// # Panics
    ///
    /// If `data.len()` is not this transform's size.
    pub fn inverse(&self, data: &mut [Complex]) {
        for value in data.iter_mut() {
            *value = value.conjugate();
        }
        self.forward(data);
        let scale = 1.0 / self.size as f64;
        for value in data.iter_mut() {
            *value = Complex::new(value.re * scale, -value.im * scale);
        }
    }

    /// The forward transform of a real signal, as a fresh buffer.
    #[must_use]
    pub fn forward_real(&self, signal: &[f64]) -> Vec<Complex> {
        let mut buffer: Vec<Complex> = signal.iter().copied().map(Complex::real).collect();
        self.forward(&mut buffer);
        buffer
    }
}

/// The transform by its definition, `O(n^2)`.
///
/// Present so that [`Fft::forward`] is checked against the thing it claims to
/// compute. A round-trip test and Parseval's identity both pass for a transform
/// whose output is permuted, which is exactly the bug a hand-written
/// bit-reversal invites, so neither is a substitute for this.
///
/// Any length, not just powers of two.
///
/// # Accuracy of the reference
///
/// A reference has to be at least as accurate as the thing it judges, and the
/// obvious transcription of the definition is not. Two things were needed, and
/// they were found in this order by measuring rather than by reasoning:
///
/// **The angle must be reduced on the integers.** Writing
/// `angle = -2 pi j k / n` and calling `cos` on it means evaluating `cos` at up
/// to `2 pi (n-1)^2 / n` radians, which at `n = 4096` is 2.6e4. The argument
/// itself then carries `eps * 2.6e4 = 5.7e-12` radians of error, and that lands
/// directly in the twiddle. Reducing `j * k` modulo `n` first keeps every
/// argument inside one turn and costs nothing, because the reduction is exact
/// integer arithmetic. This was the whole discrepancy: 1.7e-13 became 3e-15.
///
/// **The sum is compensated.** A plain sequential sum of `n` terms accumulates
/// about `sqrt(n) * eps`, while the FFT accumulates only `log2(n) * eps`.
#[must_use]
pub fn dft_naive(data: &[Complex]) -> Vec<Complex> {
    let n = data.len();
    let mut out = vec![Complex::default(); n];
    for (k, slot) in out.iter_mut().enumerate() {
        let mut sum = Complex::default();
        let mut compensation = Complex::default();
        for (j, value) in data.iter().enumerate() {
            // Exact on the integers, so the trig argument never leaves one turn.
            let reduced = (j * k) % n;
            let angle = -2.0 * core::f64::consts::PI * reduced as f64 / n as f64;
            let term = *value * Complex::new(libm::cos(angle), libm::sin(angle));
            // Kahan, on each component.
            let adjusted = term - compensation;
            let next = sum + adjusted;
            compensation = (next - sum) - adjusted;
            sum = next;
        }
        *slot = sum;
    }
    out
}

/// Reorder in place so that each element sits at its bit-reversed index.
///
/// Decimation in time needs the input in this order for the butterflies to
/// read contiguous pairs. The swap is done only when `i < reversed`, so each
/// pair is exchanged once rather than twice back to where it started.
fn bit_reverse_permute(data: &mut [Complex]) {
    let n = data.len();
    if n <= 2 {
        return;
    }
    let bits = n.trailing_zeros();
    for i in 0..n {
        let reversed = (i as u64).reverse_bits() >> (64 - bits);
        let reversed = reversed as usize;
        if i < reversed {
            data.swap(i, reversed);
        }
    }
}
