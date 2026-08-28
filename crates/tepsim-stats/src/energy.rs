//! The energy distance between two samples.
//!
//! ```text
//! E(X, Y) = 2 E|X - Y| - E|X - X'| - E|Y - Y'|
//! ```
//!
//! Zero exactly when the two distributions are equal, and positive otherwise.
//! Unlike KS it uses the whole sample rather than the single worst point of the
//! CDF gap, so it is sensitive to differences spread thinly across the support,
//! which is the shape a systematically slightly-wrong simulator produces.
//!
//! Szekely, G. J. and Rizzo, M. L. (2013), "Energy statistics: a class of
//! statistics based on distances", *Journal of Statistical Planning and
//! Inference* 143(8), 1249-1272.
//!
//! # Cost
//!
//! The definition is a double sum: `O(n*m + n^2 + m^2)`. Tier 5 compares runs
//! of 172,800 samples, where that is about 9e10 distance evaluations per
//! comparison per measurement, which is not affordable.
//!
//! In one dimension it collapses. For sorted `x`,
//!
//! ```text
//! sum_{i,j} |x_i - x_j| = 2 sum_j (2j - n - 1) x_j        (j one-based)
//! ```
//!
//! because each pair contributes once with each sign, and the cross term falls
//! out of a single merge walk. So the whole thing is `O(n log n)` for the sort
//! and `O(n + m)` after it. [`energy_distance`] does that;
//! [`energy_distance_naive`] is the definition, kept so the fast path has
//! something to be checked against.

use alloc::vec::Vec;

/// The energy distance between two samples, in `O((n + m) log(n + m))`.
///
/// Returns `NaN` if either sample is empty or contains a `NaN`.
#[must_use]
pub fn energy_distance(a: &[f64], b: &[f64]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return f64::NAN;
    }
    if a.iter().chain(b).any(|v| v.is_nan()) {
        return f64::NAN;
    }
    // Centred first. Every term is a distance, so the whole statistic is
    // invariant under a shift applied to both samples, but the *sums* that
    // compute it are not conditioned that way.
    //
    // `self_absolute_sum` accumulates `(2j - n - 1) * x_j`, whose weights sum
    // to zero. On plant data at 2705 kPa with a spread of tenths, those terms
    // reach 1e6 and the answer is 1.6e4: four digits gone before the outer
    // subtraction removes three more. Measured against the definition, the
    // uncentred version was wrong in the eighth significant figure. Centred,
    // the terms are the size of the answer and nothing cancels.
    let offset = midpoint(a, b);
    let x = sorted_centred(a, offset);
    let y = sorted_centred(b, offset);
    let n = x.len() as f64;
    let m = y.len() as f64;

    let cross = 2.0 * cross_absolute_sum(&x, &y) / (n * m);
    let within_x = self_absolute_sum(&x) / (n * n);
    let within_y = self_absolute_sum(&y) / (m * m);
    cross - within_x - within_y
}

/// A shift that brings both samples near zero.
///
/// The mean of the two sample means rather than the pooled mean: it needs no
/// weighting by size, and any value inside the data's range serves. What
/// matters is only that it is the *same* value for both samples, since a
/// different shift for each would change the cross term.
fn midpoint(a: &[f64], b: &[f64]) -> f64 {
    0.5 * (kahan_sum(a) / a.len() as f64 + kahan_sum(b) / b.len() as f64)
}

fn sorted_centred(v: &[f64], offset: f64) -> Vec<f64> {
    let mut out: Vec<f64> = v.iter().map(|x| x - offset).collect();
    out.sort_unstable_by(f64::total_cmp);
    out
}

/// The definition, `O(n*m + n^2 + m^2)`.
///
/// Present so that [`energy_distance`] is checked against the thing it claims
/// to compute rather than against itself. Not for production use: see the
/// module docs for what it costs.
#[must_use]
pub fn energy_distance_naive(a: &[f64], b: &[f64]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return f64::NAN;
    }
    let mean_pairwise = |p: &[f64], q: &[f64]| {
        let mut total = 0.0;
        for u in p {
            for v in q {
                total += libm::fabs(u - v);
            }
        }
        total / (p.len() as f64 * q.len() as f64)
    };
    2.0 * mean_pairwise(a, b) - mean_pairwise(a, a) - mean_pairwise(b, b)
}

/// `sum_{i,j} |x_i - x_j|` for sorted `x`, in one pass.
///
/// Every ordered pair `(i, j)` with `i < j` contributes `x_j - x_i` twice, once
/// for each order. Collecting by index, `x_j` appears `+(j-1)` times as the
/// larger and `-(n-j)` times as the smaller, one-based.
fn self_absolute_sum(x: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mut total = 0.0;
    for (index, value) in x.iter().enumerate() {
        let j = (index + 1) as f64;
        total += (2.0 * j - n - 1.0) * value;
    }
    2.0 * total
}

/// `sum_{i,j} |x_i - y_j|` for sorted `x` and `y`, by a merge walk.
///
/// For a fixed `y_j`, splitting `x` at `y_j` gives
///
/// ```text
/// sum_i |x_i - y_j| = (below * y_j - sum_below) + (sum_above - above * y_j)
/// ```
///
/// and walking `y` in order means the split point only moves forward.
fn cross_absolute_sum(x: &[f64], y: &[f64]) -> f64 {
    let total_x: f64 = kahan_sum(x);
    let mut below_count = 0.0_f64;
    let mut below_sum = 0.0_f64;
    let mut index = 0;
    let mut total = 0.0;
    for value in y {
        while index < x.len() && x[index] <= *value {
            below_sum += x[index];
            below_count += 1.0;
            index += 1;
        }
        let above_count = x.len() as f64 - below_count;
        let above_sum = total_x - below_sum;
        total += (below_count * value - below_sum) + (above_sum - above_count * value);
    }
    total
}

/// Compensated summation.
///
/// `total_x` is subtracted from a running partial sum in [`cross_absolute_sum`],
/// so an error in it does not cancel out: it biases every term. Plant data runs
/// to 172,800 samples around 2705 kPa, where a plain sum loses about five
/// digits.
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
