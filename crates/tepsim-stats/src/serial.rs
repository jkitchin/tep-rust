//! Serial structure: the autocorrelation function, and Welch power spectra.
//!
//! These are the Tier 5 tests that see *time*. Two runs can agree on every
//! moment and every marginal distribution while one oscillates and the other
//! does not: shuffling a series leaves its mean, variance, histogram and
//! Kolmogorov-Smirnov statistic untouched and destroys its spectrum entirely.
//! `PLAN.org` names the failure this catches directly, "the Fortran
//! limit-cycles, the port damps".

use alloc::vec;
use alloc::vec::Vec;

use crate::fft::Fft;
use crate::special::sqrt;

/// The sample autocorrelation function, `r(0) .. r(max_lag)`.
///
/// ```text
///          sum_{t=0}^{n-1-k} (x_t - xbar)(x_{t+k} - xbar)
/// r(k) =  ------------------------------------------------
///              sum_{t=0}^{n-1} (x_t - xbar)^2
/// ```
///
/// The denominator is the full sum at every lag, not the `n - k` terms the
/// numerator has. That is the *biased* estimator, and it is the right one:
/// dividing by `n - k` gives an unbiased estimate at each lag separately but a
/// sequence that need not be positive semidefinite, so its implied spectrum can
/// be negative. The bias is a factor `(n - k) / n`, which at the lags Tier 5
/// uses (200 out of 172,800) is 0.999.
///
/// `r(0)` is exactly 1 by construction.
///
/// Returns an empty vector if the series is shorter than two points or has zero
/// variance, since the ratio is then `0/0` and any answer would be invented.
///
/// # Cost
///
/// Computed through the FFT, by the Wiener-Khinchin theorem: the
/// autocovariance is the inverse transform of the power spectrum. That makes
/// it `O(n log n)` *independent of `max_lag`*, where the direct double loop is
/// `O(n * max_lag)`.
///
/// At Tier 5's own size the two are close, and the scaling is the point rather
/// than the constant. Measured in release on one 172,800-sample series to lag
/// 200: 41 ms here against 65 ms direct, a factor of 1.6. Going to lag 1000
/// leaves this unchanged and makes the direct form five times worse, so the
/// choice buys headroom rather than a big win today.
///
/// The padding to at least `2n` is what makes it exact rather than circular.
/// Without it, lag `k` would pick up the wrap-around term `x[n-k..] * x[..k]`,
/// which for a series with any trend is not small.
///
/// [`autocorrelation_direct`] is the definition, kept as the reference.
#[must_use]
pub fn autocorrelation(x: &[f64], max_lag: usize) -> Vec<f64> {
    if x.len() < 2 {
        return Vec::new();
    }
    let n = x.len();
    let mean = kahan_sum(x) / n as f64;
    let centred: Vec<f64> = x.iter().map(|v| v - mean).collect();

    let lags = max_lag.min(n - 1);

    // At least 2n, so the circular correlation the transform computes has no
    // wrap-around inside the lags being read.
    let padded = (2 * n).next_power_of_two();
    let fft = Fft::new(padded);
    let mut buffer: Vec<crate::fft::Complex> = centred
        .iter()
        .map(|v| crate::fft::Complex::real(*v))
        .chain(core::iter::repeat_n(
            crate::fft::Complex::default(),
            padded - n,
        ))
        .collect();
    fft.forward(&mut buffer);
    for value in &mut buffer {
        // |X|^2, as a real number: the power spectrum.
        *value = crate::fft::Complex::real(value.norm_squared());
    }
    fft.inverse(&mut buffer);

    let denominator = buffer[0].re;
    if denominator <= 0.0 {
        return Vec::new();
    }
    (0..=lags).map(|k| buffer[k].re / denominator).collect()
}

/// The autocorrelation by its definition, `O(n * max_lag)`.
///
/// Present so that [`autocorrelation`] is checked against what it claims to
/// compute rather than against itself. The transform route is exact in exact
/// arithmetic, but "exact in exact arithmetic" is a claim about the algebra
/// and not about the code.
#[must_use]
pub fn autocorrelation_direct(x: &[f64], max_lag: usize) -> Vec<f64> {
    if x.len() < 2 {
        return Vec::new();
    }
    let n = x.len();
    let mean = kahan_sum(x) / n as f64;
    let centred: Vec<f64> = x.iter().map(|v| v - mean).collect();

    let denominator: f64 = centred.iter().map(|v| v * v).sum();
    if denominator == 0.0 {
        return Vec::new();
    }

    let lags = max_lag.min(n - 1);
    let mut out = Vec::with_capacity(lags + 1);
    for k in 0..=lags {
        let numerator: f64 = (0..n - k).map(|t| centred[t] * centred[t + k]).sum();
        out.push(numerator / denominator);
    }
    out
}

/// The lag at which the autocorrelation first falls below `threshold`.
///
/// A single number summarising how long the series remembers, for the report
/// table. `None` if it has not fallen that far by `max_lag`, which is itself
/// worth reporting: it means the series is more correlated than the window can
/// see.
#[must_use]
pub fn decorrelation_lag(acf: &[f64], threshold: f64) -> Option<usize> {
    acf.iter().position(|r| *r < threshold)
}

/// Which window a spectrum was taken with.
///
/// Named rather than passed as a closure so that a report can say which one was
/// used, and so that the noise bandwidth travels with it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Window {
    /// No taper. Maximum frequency resolution, terrible sidelobes: a strong
    /// tone leaks across the whole band. Present mainly so the Hann case has
    /// something to be compared against.
    Rectangular,
    /// The periodic (DFT-even) Hann window, `0.5 (1 - cos(2 pi i / L))`.
    ///
    /// Periodic and not symmetric: the symmetric form `2 pi i / (L - 1)`
    /// belongs to filter design, and using it here makes the window
    /// non-periodic over the segment, which puts a small error in every bin.
    /// Noise bandwidth 1.5 bins.
    Hann,
}

impl Window {
    /// The window's samples for a segment of length `size`.
    #[must_use]
    pub fn samples(self, size: usize) -> Vec<f64> {
        match self {
            Self::Rectangular => vec![1.0; size],
            Self::Hann => (0..size)
                .map(|i| {
                    0.5 * (1.0 - libm::cos(2.0 * core::f64::consts::PI * i as f64 / size as f64))
                })
                .collect(),
        }
    }

    /// Equivalent noise bandwidth, in bins.
    ///
    /// The factor by which the window widens a spectral line. 1 for
    /// rectangular, 1.5 for Hann. Needed to turn a peak height back into a
    /// power.
    #[must_use]
    pub const fn noise_bandwidth(self) -> f64 {
        match self {
            Self::Rectangular => 1.0,
            Self::Hann => 1.5,
        }
    }
}

/// A one-sided power spectral density.
#[derive(Clone, Debug, PartialEq)]
pub struct Spectrum {
    /// Bin centre frequencies, in the same units as the sample rate.
    pub frequencies: Vec<f64>,
    /// Power spectral density, in signal units squared per unit frequency.
    pub density: Vec<f64>,
    /// The spacing between bins.
    pub resolution: f64,
    /// How many segments were averaged.
    pub segments: usize,
    /// Which window was used.
    pub window: Window,
}

impl Spectrum {
    /// The power in a frequency band, by integrating the density over it.
    ///
    /// Half-open, `[low, high)`, so adjacent bands sharing an edge do not
    /// double-count the bin on it.
    ///
    /// The consequence worth knowing: a band whose upper edge is *exactly*
    /// Nyquist excludes the Nyquist bin. To cover the whole spectrum, the top
    /// edge has to sit above it. That is why
    /// [`crate::serial::band_comparison`]'s callers should pass a high edge
    /// past `sample_rate / 2`, and why a set of bands from the resolution to
    /// Nyquist accounts for the total power minus the DC and Nyquist bins
    /// rather than all of it.
    #[must_use]
    pub fn band_power(&self, low: f64, high: f64) -> f64 {
        self.frequencies
            .iter()
            .zip(&self.density)
            .filter(|(f, _)| **f >= low && **f < high)
            .map(|(_, p)| p * self.resolution)
            .sum()
    }

    /// The total power, which for a mean-removed signal is its variance.
    #[must_use]
    pub fn total_power(&self) -> f64 {
        self.density.iter().map(|p| p * self.resolution).sum()
    }
}

/// Welch's method: split, detrend, window, transform, average.
///
/// Welch, P. D. (1967), "The use of fast Fourier transform for the estimation
/// of power spectra", *IEEE Transactions on Audio and Electroacoustics*
/// 15(2), 70-73.
///
/// `segment` must be a power of two. `overlap` is the number of samples two
/// consecutive segments share; 50% is the usual choice with Hann, because
/// that is where the windows sum to a constant and no data is under-weighted.
///
/// # Normalisation
///
/// ```text
/// P_k = |X_k|^2 / (fs * sum_i w_i^2)
/// ```
///
/// with bins strictly between DC and Nyquist doubled to make the spectrum
/// one-sided. Dividing by `sum w^2` rather than by `(sum w)^2` or by `L` is
/// what makes the result a *density* whose integral is the signal's variance.
/// Get it wrong and every spectrum is off by a constant factor, which no test
/// that compares two spectra to each other would ever notice; hence the
/// absolute-level tests beside this.
///
/// **Each segment has its mean removed.** Plant measurements sit at 2705 kPa
/// and vary by tenths, so the DC bin would otherwise be seven orders above
/// everything else and its leakage would swamp the low-frequency bins that
/// matter. The consequence is that `density[0]` is not meaningful and
/// [`Spectrum::total_power`] is the variance rather than the mean square.
///
/// # Panics
///
/// If `segment` is not a power of two, if `overlap >= segment`, or if the
/// signal is shorter than one segment. All three are programming errors rather
/// than data-dependent conditions.
#[must_use]
pub fn welch(
    signal: &[f64],
    sample_rate: f64,
    segment: usize,
    overlap: usize,
    window: Window,
) -> Spectrum {
    assert!(
        segment.is_power_of_two(),
        "segment length {segment} is not a power of two"
    );
    assert!(
        overlap < segment,
        "overlap {overlap} is not less than {segment}"
    );
    assert!(
        signal.len() >= segment,
        "signal of {} samples is shorter than one {segment}-sample segment",
        signal.len()
    );

    let taper = window.samples(segment);
    let window_power: f64 = taper.iter().map(|w| w * w).sum();
    let fft = Fft::new(segment);
    let bins = segment / 2 + 1;

    let step = segment - overlap;
    let mut accumulated = vec![0.0_f64; bins];
    let mut segments = 0_usize;

    let mut start = 0;
    while start + segment <= signal.len() {
        let slice = &signal[start..start + segment];
        let mean = kahan_sum(slice) / segment as f64;
        let tapered: Vec<f64> = slice
            .iter()
            .zip(&taper)
            .map(|(x, w)| (x - mean) * w)
            .collect();

        let spectrum = fft.forward_real(&tapered);
        for (bin, slot) in accumulated.iter_mut().enumerate() {
            // Bins strictly inside carry the power of their negative-frequency
            // twin as well.
            let fold = if bin == 0 || bin == segment / 2 {
                1.0
            } else {
                2.0
            };
            *slot += fold * spectrum[bin].norm_squared() / (sample_rate * window_power);
        }
        segments += 1;
        start += step;
    }

    let scale = 1.0 / segments as f64;
    let resolution = sample_rate / segment as f64;
    Spectrum {
        frequencies: (0..bins).map(|k| k as f64 * resolution).collect(),
        density: accumulated.iter().map(|p| p * scale).collect(),
        resolution,
        segments,
        window,
    }
}

/// Band-by-band comparison of two spectra.
///
/// `PLAN.org` asks for this rather than a single scalar: a port whose spectrum
/// is right everywhere except in one band would pass any norm over the whole
/// band and fail here, and the band it fails in says what is wrong.
///
/// Each entry is `(low, high, power_a, power_b, ratio)`. The ratio is
/// `power_a / power_b`, or `NaN` where `power_b` is zero.
#[must_use]
pub fn band_comparison(a: &Spectrum, b: &Spectrum, edges: &[f64]) -> Vec<BandComparison> {
    edges
        .windows(2)
        .map(|pair| {
            let (low, high) = (pair[0], pair[1]);
            let power_a = a.band_power(low, high);
            let power_b = b.band_power(low, high);
            BandComparison {
                low,
                high,
                power_a,
                power_b,
                ratio: if power_b == 0.0 {
                    f64::NAN
                } else {
                    power_a / power_b
                },
            }
        })
        .collect()
}

/// One band's worth of [`band_comparison`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BandComparison {
    /// Lower edge, inclusive.
    pub low: f64,
    /// Upper edge, exclusive.
    pub high: f64,
    /// Power in the band for the first spectrum.
    pub power_a: f64,
    /// And for the second.
    pub power_b: f64,
    /// `power_a / power_b`.
    pub ratio: f64,
}

impl core::fmt::Display for BandComparison {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "[{:.5e}, {:.5e})  a={:.6e}  b={:.6e}  ratio={:.6}",
            self.low, self.high, self.power_a, self.power_b, self.ratio
        )
    }
}

/// Logarithmically spaced band edges from `low` to `high`.
///
/// Log rather than linear because process dynamics live on a log frequency
/// axis: the difference between a 1-hour and a 2-hour oscillation matters as
/// much as the one between a 1-minute and a 2-minute one, and linear bands
/// would put both in the same bucket.
#[must_use]
pub fn log_band_edges(low: f64, high: f64, count: usize) -> Vec<f64> {
    let ratio = crate::special::ln(high / low) / count as f64;
    (0..=count)
        .map(|i| low * crate::special::exp(ratio * i as f64))
        .collect()
}

fn kahan_sum(v: &[f64]) -> f64 {
    let mut sum = 0.0_f64;
    let mut compensation = 0.0_f64;
    for value in v {
        let adjusted = value - compensation;
        let next = sum + adjusted;
        compensation = (next - sum) - adjusted;
        sum = next;
    }
    sum
}

/// The standard error of a sample autocorrelation at lag `k` under the null
/// that the series is white noise: `1 / sqrt(n)`.
///
/// The band everyone draws on a correlogram. It is only valid for white noise,
/// and a plant measurement is not white noise, so use
/// [`bartlett_standard_error`] to judge one.
#[must_use]
pub fn autocorrelation_standard_error(n: usize) -> f64 {
    1.0 / sqrt(n as f64)
}

/// The standard error of `r(k)` for a series with the given autocorrelation.
///
/// Bartlett's formula:
///
/// ```text
/// Var[r(k)] ~ (1/n) sum_j [ rho(j+k) + rho(j-k) - 2 rho(k) rho(j) ]^2
/// ```
///
/// Bartlett, M. S. (1946), "On the theoretical specification and sampling
/// properties of autocorrelated time-series", *Supplement to the Journal of
/// the Royal Statistical Society* 8(1), 27-41.
///
/// For white noise this collapses to `1 / sqrt(n)`. For anything correlated it
/// is larger, often much larger: an AR(1) with `phi = 0.95` has a standard
/// error at large lag of about `4.4 / sqrt(n)`, so judging it against the
/// white-noise band would reject a perfectly good estimate four times out of
/// five.
///
/// `acf` supplies `rho`, and is taken as zero beyond its end. It therefore has
/// to reach far enough out that the tail it drops is genuinely negligible: at
/// `phi = 0.95`, twenty lags leaves `rho = 0.36` and truncating there
/// understates the error badly. Passing an ACF that extends several times
/// further than the lag of interest is the fix, and this asserts nothing about
/// that because it cannot know the series.
#[must_use]
pub fn bartlett_standard_error(acf: &[f64], lag: usize, n: usize) -> f64 {
    if acf.is_empty() || n == 0 {
        return f64::NAN;
    }
    let rho = |j: isize| -> f64 {
        let index = j.unsigned_abs();
        if index < acf.len() { acf[index] } else { 0.0 }
    };
    let reach = (acf.len() - 1) as isize;
    let k = lag as isize;
    let rho_k = rho(k);
    let mut sum = 0.0;
    for j in -reach..=reach {
        let term = rho(j + k) + rho(j - k) - 2.0 * rho_k * rho(j);
        sum += term * term;
    }
    sqrt(sum / n as f64)
}
