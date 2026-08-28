//! The cross-correlation matrix, and the distance between two of them.
//!
//! `PLAN.org` singles this out because it is what PCA-based fault detectors
//! actually consume: a detector trained on one simulator and tested on another
//! sees the two only through their correlation structure. Two runs can agree on
//! every marginal distribution and disagree here, if the port has the right
//! variables with the wrong relationships between them.

use alloc::vec;
use alloc::vec::Vec;

use crate::special::sqrt;

/// A symmetric correlation matrix, stored whole.
///
/// `variables` by `variables`, row-major. Small enough to store densely: Tier 5
/// uses 53 variables, so this is 22 KB.
#[derive(Clone, Debug, PartialEq)]
pub struct CorrelationMatrix {
    variables: usize,
    entries: Vec<f64>,
}

impl CorrelationMatrix {
    /// The Pearson correlation matrix of a set of series.
    ///
    /// `series[v][t]` is variable `v` at time `t`. Every series must have the
    /// same length.
    ///
    /// A variable with zero variance correlates with nothing: its row and
    /// column are `NaN` off the diagonal and 1 on it. `NaN` rather than zero
    /// because "no linear relationship" and "this variable never moved" are
    /// different findings, and a constant column in a Tier 5 run is itself
    /// worth seeing.
    ///
    /// # Panics
    ///
    /// If the series have different lengths, or if there are fewer than two
    /// samples.
    #[must_use]
    pub fn of(series: &[Vec<f64>]) -> Self {
        let variables = series.len();
        if variables == 0 {
            return Self {
                variables: 0,
                entries: Vec::new(),
            };
        }
        let samples = series[0].len();
        assert!(samples >= 2, "a correlation needs at least two samples");
        for (index, s) in series.iter().enumerate() {
            assert_eq!(
                s.len(),
                samples,
                "series {index} has {} samples, not {samples}",
                s.len()
            );
        }

        // Centre once. Doing it inside the double loop would repeat the work
        // and, worse, would compute each mean twice with no guarantee the two
        // agree in the last bits.
        let centred: Vec<Vec<f64>> = series
            .iter()
            .map(|s| {
                let mean = kahan_sum(s) / samples as f64;
                s.iter().map(|v| v - mean).collect()
            })
            .collect();
        // The *sums of squares*, not the norms. The correlation is then
        // `dot / sqrt(ss_i * ss_j)`, one square root of a product rather than
        // the product of two square roots.
        //
        // That is not a micro-optimisation, it is the difference between right
        // and nearly right. For `x = [1,2,3,4]` and `y = 2x` the sums of
        // squares are 5 and 20, whose product is exactly 100 and whose square
        // root is exactly 10, so the correlation comes out exactly 1. Taking
        // `sqrt(5) * sqrt(20)` instead gives 10.000000000000002 and a
        // correlation one ulp below 1, which then has to be clamped and still
        // fails an exact test.
        //
        // Overflow is not reachable: a sum of squares would have to exceed
        // 1e154 for the product to leave `f64`.
        let sum_squares: Vec<f64> = centred
            .iter()
            .map(|c| c.iter().map(|v| v * v).sum::<f64>())
            .collect();

        let mut entries = vec![0.0; variables * variables];
        for i in 0..variables {
            for j in i..variables {
                let value = if i == j {
                    // Exactly one, not the computed ratio. A variable is
                    // perfectly correlated with itself, including a constant
                    // one, and rounding must not be allowed to say otherwise.
                    1.0
                } else if sum_squares[i] == 0.0 || sum_squares[j] == 0.0 {
                    f64::NAN
                } else {
                    let dot: f64 = centred[i].iter().zip(&centred[j]).map(|(a, b)| a * b).sum();
                    // Clamped: the ratio can still exceed one by an ulp on
                    // near-collinear data, and a correlation above one is a
                    // number no downstream consumer expects.
                    (dot / sqrt(sum_squares[i] * sum_squares[j])).clamp(-1.0, 1.0)
                };
                entries[i * variables + j] = value;
                entries[j * variables + i] = value;
            }
        }
        Self { variables, entries }
    }

    /// How many variables.
    #[must_use]
    pub const fn variables(&self) -> usize {
        self.variables
    }

    /// One entry.
    ///
    /// # Panics
    ///
    /// If either index is out of range.
    #[must_use]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        assert!(
            i < self.variables && j < self.variables,
            "({i}, {j}) is outside a {}-variable matrix",
            self.variables
        );
        self.entries[i * self.variables + j]
    }

    /// The entries, row-major.
    #[must_use]
    pub fn entries(&self) -> &[f64] {
        &self.entries
    }
}

/// The Frobenius norm of the difference between two correlation matrices.
///
/// ```text
/// ||A - B||_F = sqrt( sum_ij (a_ij - b_ij)^2 )
/// ```
///
/// The single number `PLAN.org` asks for. Reported alongside
/// [`worst_correlation_difference`], because a norm says how much the structure
/// moved and not *where*, and the pair of variables that moved is what a reader
/// needs in order to act.
///
/// `NaN` entries are skipped on both sides and counted; a matrix with a
/// constant variable would otherwise poison the whole norm.
///
/// Returns `(norm, skipped)`.
///
/// # Panics
///
/// If the two matrices have different sizes.
#[must_use]
pub fn frobenius_distance(a: &CorrelationMatrix, b: &CorrelationMatrix) -> (f64, usize) {
    assert_eq!(
        a.variables(),
        b.variables(),
        "cannot compare a {}-variable matrix with a {}-variable one",
        a.variables(),
        b.variables()
    );
    let mut total = 0.0;
    let mut skipped = 0;
    for (x, y) in a.entries().iter().zip(b.entries()) {
        if x.is_nan() || y.is_nan() {
            skipped += 1;
            continue;
        }
        let difference = x - y;
        total += difference * difference;
    }
    (sqrt(total), skipped)
}

/// The pair of variables whose correlation moved most between two matrices.
///
/// Returns `(i, j, a_ij, b_ij)`, or `None` if every entry is `NaN` or the
/// matrices are empty.
///
/// # Panics
///
/// If the two matrices have different sizes.
#[must_use]
pub fn worst_correlation_difference(
    a: &CorrelationMatrix,
    b: &CorrelationMatrix,
) -> Option<(usize, usize, f64, f64)> {
    assert_eq!(a.variables(), b.variables(), "size mismatch");
    let n = a.variables();
    let mut worst: Option<(f64, usize, usize)> = None;
    for i in 0..n {
        // Upper triangle only: the matrix is symmetric, so reporting the lower
        // half would just duplicate every answer.
        for j in (i + 1)..n {
            let (x, y) = (a.get(i, j), b.get(i, j));
            if x.is_nan() || y.is_nan() {
                continue;
            }
            let difference = libm::fabs(x - y);
            if worst.is_none_or(|(seen, _, _)| difference > seen) {
                worst = Some((difference, i, j));
            }
        }
    }
    worst.map(|(_, i, j)| (i, j, a.get(i, j), b.get(i, j)))
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
