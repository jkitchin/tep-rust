//! First and second moments, computed so that they survive the data.
//!
//! # Why not the textbook formula
//!
//! The computational formula for the variance,
//!
//! ```text
//! s^2 = (sum(x^2) - n * mean^2) / (n - 1)
//! ```
//!
//! subtracts two nearly equal large numbers whenever the data's mean is large
//! relative to its spread. Plant measurements are exactly that shape: reactor
//! pressure sits near 2705 kPa and moves by tenths. In `f64` the formula can
//! return a *negative* variance.
//!
//! This uses Welford's recurrence instead, which never forms `sum(x^2)`:
//!
//! ```text
//! n   <- n + 1
//! d   <- x - mean
//! mean <- mean + d / n
//! M2  <- M2 + d * (x - mean)     // the new mean, deliberately
//! ```
//!
//! Welford, B. P. (1962), "Note on a method for calculating corrected sums of
//! squares and products", *Technometrics* 4(3), 419-420.
//!
//! [`Summary::of`] folds the recurrence over a slice in index order, so the
//! result depends on the order of the data but not on anything else. That is
//! the determinism rule: a reordered reduction is a different number.

use crate::special::sqrt;

/// A running count, mean and sum of squared deviations.
///
/// Combine two with [`Summary::merge`] (Chan, Golub and LeVeque 1979), which is
/// how a statistic over several runs is accumulated without keeping the data.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Summary {
    n: usize,
    mean: f64,
    /// The sum of squared deviations from the current mean, `M2`.
    m2: f64,
}

impl Summary {
    /// An empty summary.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            n: 0,
            mean: 0.0,
            m2: 0.0,
        }
    }

    /// Summarise a slice, folding in index order.
    #[must_use]
    pub fn of(xs: &[f64]) -> Self {
        let mut summary = Self::new();
        for x in xs {
            summary.push(*x);
        }
        summary
    }

    /// Add one observation.
    pub fn push(&mut self, x: f64) {
        self.n += 1;
        let delta = x - self.mean;
        // Dividing by the *new* n, which is why the count is incremented first.
        self.mean += delta / self.n as f64;
        // The second deviation uses the *updated* mean. Using `delta` twice is
        // the classic transcription error and it biases M2 upward.
        self.m2 += delta * (x - self.mean);
    }

    /// How many observations.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// The arithmetic mean. Zero for an empty summary.
    #[must_use]
    pub const fn mean(&self) -> f64 {
        self.mean
    }

    /// The **sample** variance, dividing by `n - 1`.
    ///
    /// `NaN` for fewer than two observations, which is the honest answer: a
    /// single point carries no information about spread, and returning zero
    /// would let a degenerate sample pass an equivalence test silently.
    #[must_use]
    pub fn variance(&self) -> f64 {
        if self.n < 2 {
            return f64::NAN;
        }
        self.m2 / (self.n - 1) as f64
    }

    /// The **population** variance, dividing by `n`.
    ///
    /// `NaN` for an empty summary.
    #[must_use]
    pub fn population_variance(&self) -> f64 {
        if self.n == 0 {
            return f64::NAN;
        }
        self.m2 / self.n as f64
    }

    /// The sample standard deviation.
    #[must_use]
    pub fn sd(&self) -> f64 {
        sqrt(self.variance())
    }

    /// The standard error of the mean, `s / sqrt(n)`.
    #[must_use]
    pub fn standard_error(&self) -> f64 {
        sqrt(self.variance() / self.n as f64)
    }

    /// Combine two independent summaries.
    ///
    /// Chan, T. F., Golub, G. H. and LeVeque, R. J. (1979), "Updating formulae
    /// and a pairwise algorithm for computing sample variances", Stanford
    /// technical report STAN-CS-79-773.
    ///
    /// The `delta^2 * n_a * n_b / n` term is what the naive "add the M2s"
    /// version omits, and omitting it understates the variance of the
    /// combination by exactly the between-group contribution.
    #[must_use]
    pub fn merge(self, other: Self) -> Self {
        if self.n == 0 {
            return other;
        }
        if other.n == 0 {
            return self;
        }
        let n_a = self.n as f64;
        let n_b = other.n as f64;
        let n = n_a + n_b;
        let delta = other.mean - self.mean;
        Self {
            n: self.n + other.n,
            mean: self.mean + delta * n_b / n,
            m2: self.m2 + other.m2 + delta * delta * n_a * n_b / n,
        }
    }
}

impl core::fmt::Display for Summary {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "n={} mean={:.6e} sd={:.6e}",
            self.n,
            self.mean(),
            self.sd()
        )
    }
}
