//! Student's t and Snedecor's F: their CDFs, and their quantiles by bisection.
//!
//! Every quantile here is found the same way: bracket, then bisect a fixed
//! number of times until the bracket collapses onto two adjacent representable
//! numbers. Bisection rather than Newton because the CDFs are monotone,
//! bisection cannot diverge on a flat tail, and it takes the same number of
//! steps on every platform, which is what makes the answer reproducible bit for
//! bit. That is the same argument `CLAUDE.md` makes about reordered reductions.

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
/// Found by bisection on a bracket widened until it straddles, as described in
/// the module documentation.
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

    bisect(lo, hi, |t| student_t_cdf(t, df) < p)
}

/// The standard normal CDF, for the `df -> infinity` limit and for sanity
/// checks against it.
#[must_use]
pub fn normal_cdf(z: f64) -> f64 {
    0.5 * libm::erfc(-z / sqrt(2.0))
}

/// The `p`-quantile of the standard normal.
///
/// Bisection on [`normal_cdf`], for the same reason [`student_t_quantile`]
/// bisects: a fixed step count is a reproducible answer.
///
/// The Jackson-Mudholkar SPE limit (see [`crate::pca::spe_limit`]) needs this
/// and nothing else does, which is why there is no separate rational
/// approximation here: one call per control limit does not justify a second
/// implementation to keep correct.
///
/// Returns `NaN` for `p` outside `(0, 1)`.
///
/// # Resolution in the far tail
///
/// `libm::erfc` underflows to zero near `z = -38`, so below about
/// `p = 1e-316` the CDF is flat and the bisection returns the edge of the
/// region where it still carries information rather than the true quantile.
/// No control limit in this crate goes anywhere near there.
#[must_use]
pub fn normal_quantile(p: f64) -> f64 {
    if p.is_nan() || p <= 0.0 || p >= 1.0 {
        return f64::NAN;
    }
    // Exact: the normal is symmetric about zero.
    if p.to_bits() == 0.5_f64.to_bits() {
        return 0.0;
    }

    let mut hi = 1.0_f64;
    while normal_cdf(hi) < p && hi < 1e300 {
        hi *= 2.0;
    }
    let mut lo = -1.0_f64;
    while normal_cdf(lo) > p && lo > -1e300 {
        lo *= 2.0;
    }

    bisect(lo, hi, |z| normal_cdf(z) < p)
}

/// The CDF of Snedecor's F distribution with `df1` and `df2` degrees of
/// freedom.
///
/// Expressed through the regularised incomplete beta, which is the standard
/// identity:
///
/// ```text
/// P(F <= x) = I_y(df1/2, df2/2),    y = df1 x / (df1 x + df2)
/// ```
///
/// `df1` and `df2` need not be integers.
///
/// Returns `NaN` for a non-positive degrees-of-freedom argument, `0` for
/// `x <= 0`, and `1` for `x` infinite.
#[must_use]
pub fn f_cdf(x: f64, df1: f64, df2: f64) -> f64 {
    if x.is_nan() || df1.is_nan() || df2.is_nan() || df1 <= 0.0 || df2 <= 0.0 {
        return f64::NAN;
    }
    if x <= 0.0 {
        return 0.0;
    }
    let u = df1 * x;
    // `u` overflowing means `y` would be computed as `inf / inf = NaN`. The
    // limit is `y = 1`, and `I_1(a, b)` is exactly 1, so say so directly rather
    // than letting the arithmetic decide.
    if !u.is_finite() {
        return 1.0;
    }
    regularized_incomplete_beta(0.5 * df1, 0.5 * df2, u / (u + df2))
}

/// The `p`-quantile of Snedecor's F: the `x` with `f_cdf(x, df1, df2) == p`.
///
/// This is the value a published F table prints. `f_quantile(0.95, 5, 20)` is
/// the 2.71 in the `alpha = 0.05` table at 5 and 20 degrees of freedom.
///
/// Bisection, bracketed from `[0, 1]` and doubled upward until it straddles.
/// The support is `[0, inf)` so the lower end of the bracket is exact and needs
/// no widening.
///
/// Returns `NaN` for `p` outside `(0, 1)` or a non-positive degrees-of-freedom
/// argument.
#[must_use]
pub fn f_quantile(p: f64, df1: f64, df2: f64) -> f64 {
    if p.is_nan()
        || p <= 0.0
        || p >= 1.0
        || df1.is_nan()
        || df2.is_nan()
        || df1 <= 0.0
        || df2 <= 0.0
    {
        return f64::NAN;
    }

    let mut hi = 1.0_f64;
    while f_cdf(hi, df1, df2) < p && hi < 1e300 {
        hi *= 2.0;
    }
    // Every doubling proved the previous `hi` was below the quantile, so the
    // last one is a valid lower bound. Starting the bisection from zero would
    // throw that away and waste the bracketing work.
    let lo = if hi > 1.0 { 0.5 * hi } else { 0.0 };

    bisect(lo, hi, |x| f_cdf(x, df1, df2) < p)
}

/// Bisect `[lo, hi]` for the point where `below` stops holding.
///
/// `below(lo)` is assumed true and `below(hi)` false; the caller establishes
/// that by bracketing. Two hundred halvings is more than an `f64` bracket can
/// absorb, and the loop stops early once the midpoint lands on an endpoint,
/// which is the bit-level statement that there is nothing between them left to
/// try. A tolerance would be a tuned number; this is not.
fn bisect(mut lo: f64, mut hi: f64, below: impl Fn(f64) -> bool) -> f64 {
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid.to_bits() == lo.to_bits() || mid.to_bits() == hi.to_bits() {
            break;
        }
        if below(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}
