//! The two-sample Kolmogorov-Smirnov statistic and its asymptotic p-value.
//!
//! KS compares whole marginal distributions rather than their first two
//! moments, which is why `PLAN.org` puts it in Tier 5 alongside TOST: two
//! samples can agree in mean and variance and still have visibly different
//! shapes, and a simulator that produced the right moments with the wrong
//! distribution would pass everything else.

use alloc::vec::Vec;

use crate::special::{exp, sqrt};

/// The two-sample Kolmogorov-Smirnov statistic,
/// `D = sup_x |F_n(x) - G_m(x)|`.
///
/// Both empirical CDFs are right-continuous step functions, so the supremum is
/// attained at one of the observed values and the walk below finds it exactly.
/// Ties are handled by advancing *both* samples past the whole tied value
/// before measuring, which is what makes `D` well defined on discrete data.
///
/// Returns `NaN` if either sample is empty or contains a `NaN`. A silent zero
/// there would read as perfect agreement.
#[must_use]
pub fn ks_statistic(a: &[f64], b: &[f64]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return f64::NAN;
    }
    if a.iter().chain(b).any(|v| v.is_nan()) {
        return f64::NAN;
    }

    let x = sorted(a);
    let y = sorted(b);
    let n = x.len();
    let m = y.len();

    let mut i = 0;
    let mut j = 0;
    let mut d = 0.0_f64;
    while i < n && j < m {
        // The next value either sample steps at.
        let value = if x[i] < y[j] { x[i] } else { y[j] };
        // Past every copy of it, in both. Measuring mid-tie would report a gap
        // the CDFs never actually have.
        while i < n && x[i] <= value {
            i += 1;
        }
        while j < m && y[j] <= value {
            j += 1;
        }
        let gap = libm::fabs(i as f64 / n as f64 - j as f64 / m as f64);
        if gap > d {
            d = gap;
        }
    }
    // No tail to check: once one sample is exhausted its CDF is 1 and the
    // other's only rises toward 1, so the gap can only shrink from here.
    d
}

fn sorted(v: &[f64]) -> Vec<f64> {
    let mut out = v.to_vec();
    // `total_cmp` rather than `partial_cmp().unwrap()`: NaN is already
    // rejected above, and this cannot panic even if that changes.
    out.sort_unstable_by(f64::total_cmp);
    out
}

/// The asymptotic p-value of a two-sample KS statistic.
///
/// The limiting distribution of `sqrt(n_e) * D` is Kolmogorov's, with the
/// effective sample size
///
/// ```text
/// n_e = n * m / (n + m)
/// ```
///
/// The argument carries Stephens' finite-sample correction,
///
/// ```text
/// lambda = (sqrt(n_e) + 0.12 + 0.11 / sqrt(n_e)) * D
/// ```
///
/// Stephens, M. A. (1970), "Use of the Kolmogorov-Smirnov, Cramer-von Mises and
/// related statistics without extensive tables", *Journal of the Royal
/// Statistical Society B* 32(1), 115-122. The correction matters at the sample
/// sizes Tier 5 uses: at `n = m = 100` it moves `lambda` by about 2%.
///
/// Asymptotic, so it is a good approximation and not an exact test. Tier 5
/// uses it to rank and report, never as the sole basis of a verdict.
#[must_use]
pub fn ks_two_sample_p(d: f64, n: usize, m: usize) -> f64 {
    if n == 0 || m == 0 || d.is_nan() {
        return f64::NAN;
    }
    let effective = (n as f64 * m as f64) / (n + m) as f64;
    let root = sqrt(effective);
    kolmogorov_q((root + 0.12 + 0.11 / root) * d)
}

/// `Q(x) = P(K > x)` for Kolmogorov's distribution: the complement of its CDF.
///
/// Two series, because neither works everywhere:
///
/// ```text
/// Q(x) = 2 sum_{k>=1} (-1)^(k-1) exp(-2 k^2 x^2)                     (1)
/// Q(x) = 1 - (sqrt(2 pi) / x) sum_{k>=1} exp(-(2k-1)^2 pi^2 / (8x^2)) (2)
/// ```
///
/// (1) converges in three or four terms for `x > 1` and is useless below it:
/// at `x = 0.2` the answer is 5e-13 and the terms are of order 1, so every
/// significant digit is lost to cancellation before the series even converges.
/// (2) is the theta-function transform of the same distribution and converges
/// faster the *smaller* `x` is. Switching at `x = 1` puts each on the side
/// where it is well conditioned; the two agree to 1e-15 across the overlap,
/// which is asserted rather than assumed.
#[must_use]
pub fn kolmogorov_q(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 1.0;
    }
    if x < 1.0 {
        // The theta-transformed series.
        const PI_SQUARED_OVER_8: f64 = core::f64::consts::PI * core::f64::consts::PI / 8.0;
        let mut sum = 0.0;
        for k in 1..=20_u32 {
            let odd = f64::from(2 * k - 1);
            let term = exp(-odd * odd * PI_SQUARED_OVER_8 / (x * x));
            sum += term;
            if term < 1e-20 * sum {
                break;
            }
        }
        1.0 - sqrt(2.0 * core::f64::consts::PI) / x * sum
    } else {
        // The alternating series.
        let mut sum = 0.0;
        for k in 1..=100_u32 {
            let k = f64::from(k);
            let term = exp(-2.0 * k * k * x * x);
            sum += if (k as u32) % 2 == 1 { term } else { -term };
            if term < 1e-20 {
                break;
            }
        }
        let q = 2.0 * sum;
        // The series can overshoot by less than an ulp at large x.
        q.clamp(0.0, 1.0)
    }
}
