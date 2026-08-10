//! Regulatory control for the Tennessee Eastman Process.
//!
//! The `Controller` trait, the multi-rate scheduler, and the classic
//! decentralized PI suite. The original Fortran spells this out as nineteen
//! near-identical copy-pasted subroutines on three different sample periods
//! (3 s, 360 s and 900 s); here it is one velocity-form implementation plus a
//! table of parameters.
//!
//! # Status
//!
//! Skeleton. Lands in phase 4; see `BACKLOG.org`.

#![forbid(unsafe_code)]
