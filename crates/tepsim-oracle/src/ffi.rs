//! Raw bindings to the gfortran-compiled original.
//!
//! gfortran lowers a subroutine `TEFUNC` to the symbol `tefunc_` and a named
//! `COMMON /PV/` block to `pv_`, with every argument passed by reference. That
//! is the entire calling convention here, which is why this crate needs no C
//! shim: Rust can name those symbols directly.
//!
//! Nothing in this module is public. [`super::Oracle`] is the only way in, and
//! it exists to concentrate every soundness argument in one place.

#![allow(non_upper_case_globals)]

/// `COMMON/PV/ XMEAS(41), XMV(12)`, from `teprob.f:209-210`.
///
/// Measurements and manipulated variables. Fortran arrays are one-based, so
/// `xmeas[0]` here is `XMEAS(1)` there.
#[repr(C)]
pub(crate) struct Pv {
    pub xmeas: [f64; 41],
    pub xmv: [f64; 12],
}

/// `COMMON/DVEC/ IDV(20)`, from `teprob.f:211-212`.
///
/// The twenty disturbance flags. `INTEGER` in gfortran is 32-bit by default,
/// and no flag here is ever compiled with `-fdefault-integer-8`.
#[repr(C)]
pub(crate) struct Dvec {
    pub idv: [i32; 20],
}

/// `COMMON/RANDSD/ G`, from `teprob.f:836`.
///
/// The entire state of the simulator's random number generator: one `f64`.
/// Being able to set this is what makes the noise sequence reproducible across
/// implementations, and reading it after each step is what Tier 3 diffs.
#[repr(C)]
pub(crate) struct Randsd {
    pub g: f64,
}

unsafe extern "C" {
    /// `SUBROUTINE TEFUNC(NN, TIME, YY, YP)`, from `teprob.f:194`.
    ///
    /// Evaluates derivatives. Not a pure function: it also advances the
    /// disturbance walks, draws measurement noise, ticks the sampled analysers
    /// and latches valve sticking. That impurity is why the Rust port splits
    /// the two roles apart; see `PLAN.org`.
    pub(crate) fn tefunc_(nn: *const i32, time: *const f64, yy: *mut f64, yp: *mut f64);

    /// `SUBROUTINE TEINIT(NN, TIME, YY, YP)`, from `teprob.f:817`.
    ///
    /// Loads the nominal steady state and sets `TIME` to zero.
    pub(crate) fn teinit_(nn: *const i32, time: *mut f64, yy: *mut f64, yp: *mut f64);

    pub(crate) static mut pv_: Pv;
    pub(crate) static mut dvec_: Dvec;
    pub(crate) static mut randsd_: Randsd;
}
