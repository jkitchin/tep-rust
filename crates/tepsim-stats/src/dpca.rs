//! Dynamic PCA: the same detector, on a lag-augmented data matrix.
//!
//! Static PCA treats each sample as independent of the last. Plant data is not:
//! reactor level at 3 p.m. is nearly reactor level at 2.57 p.m., and a control
//! loop's whole job is to create relationships between a measurement now and a
//! manipulated variable a few minutes ago. A static model sees that serial
//! structure as unexplained variance, which inflates SPE on ordinary data and
//! costs the detector sensitivity to the faults that matter.
//!
//! The fix is embarrassingly simple and is the reason this module is thin.
//! Stack each sample with its predecessors into one wide row, and fit exactly
//! the same PCA to that. The correlations the model then learns include
//! `variable i now` against `variable j three samples ago`, so a fault has to
//! break a *dynamic* relationship to raise SPE, not merely be a moving plant.
//!
//! Ku, W., Storer, R. H. and Georgakis, C. (1995), "Disturbance detection and
//! isolation by dynamic principal component analysis", *Chemometrics and
//! Intelligent Laboratory Systems* 30(1), 179-196.
//!
//! # Choosing the lag count
//!
//! Ku et al. give a procedure for it, based on counting the linear relations the
//! augmented matrix acquires at each order. This module does not implement that
//! procedure and does not pick a lag count for the caller: the TEP literature
//! almost universally uses two lags on this process, and a number chosen by the
//! caller and written in the log is worth more than one chosen here and
//! forgotten.

use alloc::vec;
use alloc::vec::Vec;

use crate::pca::{ControlLimits, Pca, Retention, Statistics};

/// Stack each sample with the `lags` samples before it.
///
/// `data` is `samples` by `variables`, row-major. The result is
/// `samples - lags` by `variables * (lags + 1)`, row-major, with row `r`
/// holding
///
/// ```text
/// [ x(r + lags) | x(r + lags - 1) | ... | x(r) ]
/// ```
///
/// so block zero is the present sample and block `l` is the one `l` steps
/// before it. Present first, because that is the convention Ku et al. use and
/// because it makes `lags = 0` the identity: the augmented matrix is then the
/// original, entry for entry.
///
/// The first `lags` samples have no complete history and are dropped, which is
/// why the row count falls. That is a real cost on short records and it is
/// visible rather than hidden: with the TEP training file's 500 samples and two
/// lags, 498 rows survive.
///
/// # Panics
///
/// If `data.len()` is not `samples * variables`, or if `lags >= samples`.
#[must_use]
pub fn augment_with_lags(
    data: &[f64],
    samples: usize,
    variables: usize,
    lags: usize,
) -> (Vec<f64>, usize, usize) {
    assert_eq!(
        data.len(),
        samples * variables,
        "a {samples}-by-{variables} matrix needs {} entries, not {}",
        samples * variables,
        data.len()
    );
    assert!(
        lags < samples,
        "{lags} lags need more than {samples} samples to leave a single row"
    );

    let rows = samples - lags;
    let columns = variables * (lags + 1);
    let mut out = vec![0.0; rows * columns];
    for r in 0..rows {
        for l in 0..=lags {
            // Block `l` is the sample `l` steps before the row's present, which
            // is at index `r + lags`. Off by one here is the classic dynamic-PCA
            // bug and it is invisible in the fitted spectrum, so the tests check
            // the layout directly rather than checking a statistic computed from
            // it.
            let source = r + lags - l;
            for v in 0..variables {
                out[r * columns + l * variables + v] = data[source * variables + v];
            }
        }
    }
    (out, rows, columns)
}

/// A fitted dynamic PCA model.
///
/// A [`Pca`] over the augmented matrix, plus the lag count needed to build a
/// row to score against it.
#[derive(Clone, Debug, PartialEq)]
pub struct Dpca {
    lags: usize,
    variables: usize,
    pca: Pca,
}

impl Dpca {
    /// Fit a dynamic model.
    ///
    /// `data` is `samples` by `variables`, row-major, in time order. With
    /// `lags = 0` this is exactly [`Pca::fit`] on the same matrix, bit for bit,
    /// which is the identity the tests use to pin the wrapper to the thing it
    /// wraps.
    ///
    /// # Panics
    ///
    /// The conditions [`augment_with_lags`] and [`Pca::fit`] panic on: a
    /// mismatched length, a lag count that consumes the record, fewer than two
    /// augmented rows.
    #[must_use]
    pub fn fit(
        data: &[f64],
        samples: usize,
        variables: usize,
        lags: usize,
        rule: Retention,
    ) -> Self {
        let (augmented, rows, columns) = augment_with_lags(data, samples, variables, lags);
        Self {
            lags,
            variables,
            pca: Pca::fit(&augmented, rows, columns, rule),
        }
    }

    /// The lag count.
    #[must_use]
    pub const fn lags(&self) -> usize {
        self.lags
    }

    /// How many variables per sample, before augmentation.
    #[must_use]
    pub const fn variables(&self) -> usize {
        self.variables
    }

    /// The underlying model over the augmented variables.
    ///
    /// Its `variables()` is `variables * (lags + 1)` and its `samples()` is the
    /// augmented row count, which is `lags` short of the record it was fitted
    /// to. Both matter when reading a control limit, which is why they are
    /// reachable rather than hidden behind the wrapper.
    #[must_use]
    pub const fn pca(&self) -> &Pca {
        &self.pca
    }

    /// Turn a window of consecutive samples into one augmented row.
    ///
    /// `window` is `lags + 1` samples of `variables` values each, row-major, in
    /// **time order, oldest first**, the same order they appear in the source
    /// matrix. The result reverses the blocks, so it is present-first like the
    /// rows of [`augment_with_lags`].
    ///
    /// # Panics
    ///
    /// If `window.len()` is not `variables * (lags + 1)`.
    #[must_use]
    pub fn augment_window(&self, window: &[f64]) -> Vec<f64> {
        let expected = self.variables * (self.lags + 1);
        assert_eq!(
            window.len(),
            expected,
            "a {}-lag window over {} variables is {expected} values, not {}",
            self.lags,
            self.variables,
            window.len()
        );
        let mut row = vec![0.0; expected];
        for l in 0..=self.lags {
            let source = self.lags - l;
            for v in 0..self.variables {
                row[l * self.variables + v] = window[source * self.variables + v];
            }
        }
        row
    }

    /// Hotelling's T-squared for one window.
    ///
    /// # Panics
    ///
    /// If the window is the wrong length.
    #[must_use]
    pub fn t_squared(&self, window: &[f64]) -> f64 {
        self.pca.t_squared(&self.augment_window(window))
    }

    /// The squared prediction error for one window.
    ///
    /// # Panics
    ///
    /// If the window is the wrong length.
    #[must_use]
    pub fn spe(&self, window: &[f64]) -> f64 {
        self.pca.spe(&self.augment_window(window))
    }

    /// Both statistics for one window.
    ///
    /// # Panics
    ///
    /// If the window is the wrong length.
    #[must_use]
    pub fn statistics(&self, window: &[f64]) -> Statistics {
        self.pca.statistics(&self.augment_window(window))
    }

    /// Both control limits at one confidence level.
    #[must_use]
    pub fn limits(&self, confidence: f64) -> ControlLimits {
        self.pca.limits(confidence)
    }
}
