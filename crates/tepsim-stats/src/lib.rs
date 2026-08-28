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
//! | [`distribution`] | Student's t CDF and quantile |
//! | [`equivalence`] | Welch's t-test, and TOST on top of it |
//! | [`ks`] | the two-sample Kolmogorov-Smirnov test |
//! | [`energy`] | the energy distance between two samples |
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

pub mod distribution;
pub mod energy;
pub mod equivalence;
pub mod ks;
pub mod special;
pub mod summary;

pub use distribution::{student_t_cdf, student_t_quantile};
pub use energy::{energy_distance, energy_distance_naive};
pub use equivalence::{Tost, WelchT, tost, welch_t};
pub use ks::{kolmogorov_q, ks_statistic, ks_two_sample_p};
pub use special::{ln_gamma, regularized_incomplete_beta};
pub use summary::Summary;
