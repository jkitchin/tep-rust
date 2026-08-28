//! The statistics behind TEP-Rust's validation ladder, implemented in Rust.
//!
//! **Development only.** Nothing shipped depends on this crate, and
//! `cargo xtask ci` asserts it. Unlike [`tepsim-oracle`] it needs no gfortran,
//! so it builds and runs everywhere, wasm included.
//!
//! [`tepsim-oracle`]: https://github.com/jkitchin/tep-rust
//!
//! # Why not scipy
//!
//! Tiers 5 through 7 are the project's actual equivalence claim. Delegating
//! them to numpy and scipy would rest that claim on a dependency's behaviour at
//! a version nobody recorded, in a language the rest of the validation does not
//! use. Every statistic here is small, and every one of them is checked against
//! a published value or an exact identity before it is pointed at plant data.
//! See the decision entry of 2026-08-28 in `BACKLOG.org`.
//!
//! # Determinism
//!
//! The same rules as the model: `f64` only, no `f32`, no reordered reductions,
//! no parallelism, and the vendored [`libm`] rather than the platform one. A
//! statistic computed here rounds identically on x86-64, aarch64 and wasm32.
//!
//! # What is here
//!
//! | Module | Tier 5 role |
//! |---|---|
//! | [`summary`] | first and second moments, stably |
//! | [`special`] | `ln_gamma` and the regularised incomplete beta |
//! | [`distribution`] | Student's t and Snedecor's F, CDFs and quantiles |
//! | [`equivalence`] | Welch's t-test, and TOST on top of it |
//! | [`ks`] | the two-sample Kolmogorov-Smirnov test |
//! | [`energy`] | the energy distance between two samples |
//! | [`fft`] | a deterministic radix-2 transform, for Welch spectra |
//! | [`serial`] | autocorrelation and Welch power spectra |
//! | [`correlation`] | the cross-correlation matrix, which PCA consumes |
//!
//! and, for Tier 6, the fault detectors that turn those into a downstream task:
//!
//! | Module | Tier 6 role |
//! |---|---|
//! | [`eigen`] | cyclic Jacobi, the deterministic symmetric eigensolver |
//! | [`pca`] | PCA, Hotelling's T-squared, SPE, and both control limits |
//! | [`dpca`] | the same on a lag-augmented matrix |
//! | [`cva`] | canonical variate analysis, over a generalised eigenproblem |
//! | [`detection`] | detection rate, false alarm rate, detection delay |
//!
//! # Reporting
//!
//! Every result type carries the numbers that produced it, not only the
//! verdict: [`equivalence::Tost`] reports the confidence interval and both
//! one-sided p-values, and [`equivalence::WelchT`] reports the statistic, the
//! degrees of freedom and the standard error. That is the same rule
//! `CLAUDE.md` applies to log entries, for the same reason: a verdict cannot
//! be compared against the previous run and a number can.

#![no_std]
#![forbid(unsafe_code)]

extern crate alloc;

pub mod correlation;
pub mod cva;
pub mod detection;
pub mod distribution;
pub mod dpca;
pub mod eigen;
pub mod energy;
pub mod equivalence;
pub mod fft;
pub mod ks;
pub mod pca;
pub mod serial;
pub mod special;
pub mod summary;

pub use correlation::{CorrelationMatrix, frobenius_distance, worst_correlation_difference};
pub use cva::{Cva, GeneralizedSymmetricEigen, cholesky, generalized_symmetric_eigen, past_future};
pub use detection::{
    DetectionReport, alarms_above, detection_delay, detection_report, false_alarm_rate,
    fault_detection_rate,
};
pub use distribution::{f_cdf, f_quantile, normal_quantile, student_t_cdf, student_t_quantile};
pub use dpca::{Dpca, augment_with_lags};
pub use eigen::{SymmetricEigen, symmetric_eigen};
pub use energy::{energy_distance, energy_distance_naive};
pub use equivalence::{OneSampleT, Tost, WelchT, one_sample_t, tost, tost_paired, welch_t};
pub use fft::{Complex, Fft, dft_naive};
pub use ks::{kolmogorov_q, ks_statistic, ks_two_sample_p};
pub use pca::{ControlLimits, Pca, Retention, Statistics, spe_limit, t_squared_limit};
pub use serial::{
    BandComparison, Spectrum, Window, autocorrelation, autocorrelation_direct, band_comparison,
    bartlett_standard_error, log_band_edges, welch,
};
pub use special::{ln_gamma, regularized_incomplete_beta};
pub use summary::Summary;
