//! Student's t distribution: the CDF, and its quantile by bisection.

use crate::special::{regularized_incomplete_beta, sqrt};

/// The CDF of Student's t distribution with `df` degrees of freedom.
///
/// Expressed through the regularised incomplete beta:
///
/// ```text
/// P(T <= t) = 1 - I_{df / (df + t^2)}(df/2, 1/2) / 2      for t >= 0
/// P(T <= t) =     I_{df / (df + t^2)}(df/2, 1/2) / 2      for t <  0
/// ```
///
/// `df` need not be an integer, which matters because Welch's test produces a
/// fractional one.
///
/// Returns `NaN` for `df <= 0`.
#[must_use]
pub fn student_t_cdf(t: f64, df: f64) -> f64 {
    if df.is_nan() || df <= 0.0 || t.is_nan() {
        return f64::NAN;
    }
    if t.is_infinite() {
        return if t > 0.0 { 1.0 } else { 0.0 };
    }
    let x = df / (df + t * t);
    let tail = 0.5 * regularized_incomplete_beta(0.5 * df, 0.5, x);
    if t >= 0.0 { 1.0 - tail } else { tail }
}

/// The two-sided p-value for a t statistic: `P(|T| >= |t|)`.
#[must_use]
pub fn student_t_two_sided_p(t: f64, df: f64) -> f64 {
    if df.is_nan() || df <= 0.0 || t.is_nan() {
        return f64::NAN;
    }
    // Computed from the upper tail directly rather than as `2 * (1 - cdf)`,
    // which loses every significant digit when `cdf` is close to one.
    let x = df / (df + t * t);
    regularized_incomplete_beta(0.5 * df, 0.5, x)
}

/// The `p`-quantile of Student's t: the `t` with `student_t_cdf(t, df) == p`.
///
/// Found by bisection on a bracket widened until it straddles. Bisection
/// rather than Newton because the CDF is monotone, bisection cannot diverge,
/// and it terminates in a fixed number of steps, which keeps the answer
/// reproducible bit for bit across platforms. Sixty halvings exhaust an `f64`
/// bracket.
///
/// Returns `NaN` for `df <= 0` or `p` outside `(0, 1)`; infinities at the ends
/// would be correct but a caller almost certainly has a bug.
#[must_use]
pub fn student_t_quantile(p: f64, df: f64) -> f64 {
    if df.is_nan() || df <= 0.0 || p.is_nan() || p <= 0.0 || p >= 1.0 {
        return f64::NAN;
    }
    // Exact: the t distribution is symmetric, so the median is exactly zero
    // and there is nothing for bisection to find.
    if p.to_bits() == 0.5_f64.to_bits() {
        return 0.0;
    }

    // Widen until the target is bracketed. Doubling from 1 reaches the tails of
    // even a df = 1 Cauchy in well under sixty steps.
    let mut hi = 1.0_f64;
    while student_t_cdf(hi, df) < p && hi < 1e300 {
        hi *= 2.0;
    }
    let mut lo = -1.0_f64;
    while student_t_cdf(lo, df) > p && lo > -1e300 {
        lo *= 2.0;
    }

    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        // Bit equality, not a tolerance: the bracket has collapsed to two
        // adjacent representable numbers and there is no midpoint left.
        if mid.to_bits() == lo.to_bits() || mid.to_bits() == hi.to_bits() {
            break;
        }
        if student_t_cdf(mid, df) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// The standard normal CDF, for the `df -> infinity` limit and for sanity
/// checks against it.
#[must_use]
pub fn normal_cdf(z: f64) -> f64 {
    0.5 * libm::erfc(-z / sqrt(2.0))
}
