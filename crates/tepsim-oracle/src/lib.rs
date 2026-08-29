//! Development-only harness that drives the *original* Tennessee Eastman
//! Fortran so the Rust port can be differentially tested against it.
//!
//! This is the instrument the whole validation ladder rests on. It compiles the
//! unmodified `reference/fortran/teprob.f` with gfortran and lets a test force
//! the Fortran into a chosen state, evaluate once, and compare against the Rust
//! port under identical conditions.
//!
//! # This crate is never shipped
//!
//! `publish = false`, and nothing in `tepsim`, `tepsim-py` or `tepsim-wasm` may
//! depend on it. `cargo xtask ci` asserts that. The `oracle` feature is off by
//! default; without it this crate is an empty shell that needs no Fortran
//! compiler, so a contributor with no gfortran can still build the workspace.
//!
//! # The Fortran has exactly one global state
//!
//! `COMMON` blocks are process-wide mutable globals, so the original supports
//! precisely one simulation per process and is not reentrant. That is a large
//! part of *why* the Rust port exists, and it is also a hazard here: Rust runs
//! tests on multiple threads by default. `Oracle` is therefore a process-wide
//! singleton guarded by a mutex, obtained through `Oracle::lock`. There is no
//! way to touch the Fortran without holding it.
//!
//! (Those two are code spans rather than links because the type only exists
//! under the `oracle` feature, and this page is built without it.)
//!
//! # Reference numbers depend on the compiler
//!
//! Tier 1 and Tier 2 tolerances are measured against a specific gfortran build
//! with pinned flags, asserted by [`build_info`]. Changing either invalidates
//! every number recorded in `LOG.org`, so it is a deliberate, logged
//! re-baseline. Never add `-ffast-math` or anything else permitting
//! reassociation.

// FFI to a Fortran library cannot be written in safe Rust. The unsafety is
// confined to `ffi` and to the accessors below, each of which states the
// argument for why it is sound. See CLAUDE.md for the project-wide rule.
#![deny(unsafe_op_in_unsafe_fn)]

pub mod build_info;
// Deliberately not gated on `oracle`: `cargo xtask fidelity` parses the golden
// trace on machines with no Fortran toolchain at all.
pub mod golden;
// Likewise: the Tier 1 pools and the ULP report are pure Rust, so they can be
// proved correct without a Fortran compiler present.
pub mod tier1;
// Gated: forcing the Fortran into a chosen state needs the Fortran.
#[cfg(feature = "oracle")]
pub mod tier2;

/// Tier 3: the generator draw trace and its differ.
pub mod tier3;

// Gated for the same reason `tier2` is: driving the Fortran's controllers
// needs the Fortran.
#[cfg(feature = "oracle")]
pub mod tier5;

// Tier 6 sits on Tier 5's runs and on the Fortran side of them, so it carries
// the same gate.
#[cfg(feature = "oracle")]
pub mod tier6;

// Tier 7 needs no Fortran of its own: it compares the port against files
// vendored under `reference/data/`. It is gated all the same, because the
// battery it judges with lives in `tier5` and that does need the Fortran.
#[cfg(feature = "oracle")]
pub mod tier7;

// Tier 8 forces the Fortran into a randomly generated state and reads its
// derivative back, so it needs the Fortran for the same reason `tier2` does.
#[cfg(feature = "oracle")]
pub mod tier8;

#[cfg(feature = "oracle")]
mod ffi;

#[cfg(feature = "oracle")]
pub use ffi::{Const, Ctrlall, Flag6, TRACE_CAPACITY, Teproc, Wlk};
#[cfg(feature = "oracle")]
pub use ffi::{
    Ctrl1, Ctrl2, Ctrl3, Ctrl4, Ctrl5, Ctrl6, Ctrl7, Ctrl8, Ctrl9, Ctrl10, Ctrl11, Ctrl13, Ctrl14,
    Ctrl15, Ctrl16, Ctrl17, Ctrl18, Ctrl19, Ctrl20, Ctrl22,
};
#[cfg(feature = "oracle")]
pub use oracle::{N_STATES, Oracle, WalkSegment, WalkSegmentStart, WalkSpans};

#[cfg(feature = "oracle")]
mod oracle {
    use std::sync::{Mutex, MutexGuard, OnceLock};

    use crate::ffi::{self, Const, Ctrlall, Teproc, Wlk};
    use crate::ffi::{
        Ctrl1, Ctrl2, Ctrl3, Ctrl4, Ctrl5, Ctrl6, Ctrl7, Ctrl8, Ctrl9, Ctrl10, Ctrl11, Ctrl13,
        Ctrl14, Ctrl15, Ctrl16, Ctrl17, Ctrl18, Ctrl19, Ctrl20, Ctrl22,
    };

    /// The original integrates 50 states. `teprob.f:24-26` invites callers to
    /// append their own beyond that; we never do, so this is exact.
    pub const N_STATES: usize = 50;

    /// Where a walk segment starts: `S`, `SP` and `TLAST` at `teprob.f:1506`.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct WalkSegmentStart {
        /// `S`, the value the previous segment ended at.
        pub value: f64,
        /// `SP`, the slope it ended with.
        pub slope: f64,
        /// `TLAST`, the time it ended.
        pub tlast: f64,
    }

    /// One channel's five span parameters and its disturbance flag.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct WalkSpans {
        /// `HSPAN`, the half-range of the segment duration.
        pub hspan: f64,
        /// `HZERO`, its centre.
        pub hzero: f64,
        /// `SSPAN`, the half-range of the endpoint value.
        pub sspan: f64,
        /// `SZERO`, its centre.
        pub szero: f64,
        /// `SPSPAN`, the half-range of the endpoint slope.
        pub spspan: f64,
        /// `IDVWLK(I)`: zero unless this channel's disturbance is active. It
        /// multiplies *both* endpoint draws (`teprob.f:1529-1530`), so an
        /// inactive channel lands exactly on `SZERO` with zero slope and only
        /// its duration stays random. The draws still happen either way.
        pub idvflag: i32,
    }

    /// The cubic `TESUB5` produces, and when it ends.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct WalkSegment {
        /// `ADIST`, the constant term.
        pub adist: f64,
        /// `BDIST`, the linear term.
        pub bdist: f64,
        /// `CDIST`, the quadratic term.
        pub cdist: f64,
        /// `DDIST`, the cubic term.
        pub ddist: f64,
        /// `TNEXT`, when this segment ends and the next is built.
        pub tnext: f64,
    }

    /// Exclusive access to the Fortran's process-wide `COMMON` state.
    ///
    /// Held for as long as you need a consistent view. Every accessor takes
    /// `&mut self`, so the borrow checker prevents interleaving two logical
    /// operations even within one thread.
    #[derive(Debug)]
    pub struct Oracle {
        _private: (),
    }

    fn mutex() -> &'static Mutex<Oracle> {
        static LOCK: OnceLock<Mutex<Oracle>> = OnceLock::new();
        LOCK.get_or_init(|| Mutex::new(Oracle { _private: () }))
    }

    impl Oracle {
        /// Take exclusive control of the Fortran.
        ///
        /// Blocks until any other holder is done. Recovers from a poisoned
        /// mutex rather than panicking: a test that failed while holding the
        /// lock left the Fortran in an arbitrary state, not an unsound one, and
        /// every entry point here overwrites the state it depends on.
        pub fn lock() -> MutexGuard<'static, Oracle> {
            mutex()
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner())
        }

        /// `TEINIT`: load the nominal steady state. Returns `(time, states)`.
        pub fn init(&mut self) -> (f64, [f64; N_STATES]) {
            let nn = N_STATES as i32;
            let mut time = 0.0_f64;
            let mut yy = [0.0_f64; N_STATES];
            let mut yp = [0.0_f64; N_STATES];
            // SAFETY: `nn` matches the length of both arrays, so the Fortran
            // cannot index past them. We hold the lock, so nothing else is
            // touching the COMMON blocks concurrently.
            unsafe { ffi::teinit_(&nn, &mut time, yy.as_mut_ptr(), yp.as_mut_ptr()) };
            (time, yy)
        }

        /// `TEINIT` from a cold `COMMON`, so the result does not depend on
        /// what ran before it.
        ///
        /// `TCR`, `TCS`, `TCC` and `TCV` are the Newton warm starts for the
        /// four vessel temperatures, and they are **never initialised**.
        /// `TESUB2` takes each as both guess and answer (`teprob.f:460`), and
        /// nothing else in `teprob.f` ever assigns them. On a freshly loaded
        /// process they are the loader's zeros; after any run they are
        /// wherever that run left them, and `TEINIT` does not put them back.
        ///
        /// So two identical runs in one process give different answers in the
        /// last bits, and which answer you get depends on the order the tests
        /// happened to run in. Zeroing the four before calling `TEINIT`
        /// reproduces the freshly-loaded-process result exactly, and does so
        /// whatever the history: verified against runs of 1, 137, 5,000 and
        /// 20,000 steps.
        ///
        /// Use this for anything that starts from the nominal state and must
        /// be comparable against another run or against a recorded number.
        /// Tier 2 does not need it, because
        /// [`crate::tier2::Scenario::force`] restores the whole of
        /// `COMMON/TEPROC/` and so specifies the warm start as part of the
        /// scenario.
        pub fn init_cold(&mut self) -> (f64, [f64; N_STATES]) {
            let mut common = self.teproc();
            common.tcr = 0.0;
            common.tcs = 0.0;
            common.tcc = 0.0;
            common.tcv = 0.0;
            self.set_teproc(&common);
            self.init()
        }

        /// `TEFUNC`: evaluate derivatives at `time` for state `yy`.
        ///
        /// Also advances the disturbance walks, draws measurement noise and
        /// ticks the sampled analysers, so calling it twice with the same
        /// arguments does not generally give the same answer. Reset
        /// [`Oracle::set_rng`] first if you need repeatability.
        pub fn derivatives(&mut self, time: f64, yy: &[f64; N_STATES]) -> [f64; N_STATES] {
            let nn = N_STATES as i32;
            let mut state = *yy;
            let mut yp = [0.0_f64; N_STATES];
            // SAFETY: as above. `state` is a private copy, so the Fortran
            // writing through `yy` cannot surprise the caller.
            unsafe { ffi::tefunc_(&nn, &time, state.as_mut_ptr(), yp.as_mut_ptr()) };
            yp
        }

        /// `XMEAS(1..41)`, the measurement vector, zero-based here.
        pub fn measurements(&mut self) -> [f64; 41] {
            // SAFETY: we hold the lock, and `Pv` mirrors the Fortran layout of
            // `COMMON/PV/`, so this reads initialised f64s written by Fortran.
            unsafe { (&raw const ffi::pv_.xmeas).read() }
        }

        /// Set `XMEAS(1..41)`.
        ///
        /// Needed because the analysers *hold* `XMEAS(23..41)` between
        /// samples (`teprob.f:744` writes them only inside the schedule
        /// check), so those nineteen are part of what `TEFUNC` reads and a
        /// scenario is not fully specified without them.
        pub fn set_measurements(&mut self, xmeas: &[f64; 41]) {
            // SAFETY: as above, writing a correctly typed and sized value.
            unsafe { (&raw mut ffi::pv_.xmeas).write(*xmeas) };
        }

        /// `XMV(1..12)`, the manipulated variables, zero-based here.
        pub fn manipulated(&mut self) -> [f64; 12] {
            // SAFETY: as above.
            unsafe { (&raw const ffi::pv_.xmv).read() }
        }

        /// Overwrite `XMV(1..12)`. This is how a controller is simulated.
        pub fn set_manipulated(&mut self, xmv: &[f64; 12]) {
            // SAFETY: as above, writing a correctly typed and sized value.
            unsafe { (&raw mut ffi::pv_.xmv).write(*xmv) };
        }

        /// The twenty `IDV` disturbance flags.
        pub fn disturbances(&mut self) -> [i32; 20] {
            // SAFETY: as above.
            unsafe { (&raw const ffi::dvec_.idv).read() }
        }

        /// Set the twenty `IDV` flags. `TEFUNC` clamps anything positive to 1.
        pub fn set_disturbances(&mut self, idv: &[i32; 20]) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::dvec_.idv).write(*idv) };
        }

        /// The RNG state: `COMMON/RANDSD/ G`, a single `f64`.
        pub fn rng(&mut self) -> f64 {
            // SAFETY: as above.
            unsafe { (&raw const ffi::randsd_.g).read() }
        }

        /// Set the RNG state. Necessary before any comparison that must be
        /// reproducible, and the basis of the Tier 3 call-order diff.
        pub fn set_rng(&mut self, g: f64) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::randsd_.g).write(g) };
        }

        /// A snapshot of `COMMON/TEPROC/`, the plant's whole working set.
        pub fn teproc(&mut self) -> Teproc {
            // SAFETY: as above.
            unsafe { (&raw const ffi::teproc_).read() }
        }

        /// Overwrite `COMMON/TEPROC/`. This is how a differential test forces
        /// the Fortran into a chosen state before evaluating.
        pub fn set_teproc(&mut self, v: &Teproc) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::teproc_).write(*v) };
        }

        /// A snapshot of `COMMON/WLK/`, the disturbance random-walk state.
        pub fn wlk(&mut self) -> Wlk {
            // SAFETY: as above.
            unsafe { (&raw const ffi::wlk_).read() }
        }

        /// Overwrite `COMMON/WLK/`. Needed to reproduce a disturbance
        /// trajectory, since the walks carry state between calls.
        pub fn set_wlk(&mut self, v: &Wlk) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::wlk_).write(*v) };
        }

        /// `COMMON/CONST/`, the thermodynamic coefficients. Read-only in
        /// practice: `TEINIT` sets them and nothing else writes them.
        pub fn constants(&mut self) -> Const {
            // SAFETY: as above.
            unsafe { (&raw const ffi::const_).read() }
        }

        /// The shutdown flag, non-zero when the plant has tripped.
        ///
        /// Only reachable because `build.rs` hoists `ISD` into
        /// `COMMON/SHUTDN/`; the original keeps it as a local. A test proves
        /// that rewrite changes no numbers.
        pub fn shutdown_flag(&mut self) -> i32 {
            // SAFETY: as above.
            unsafe { (&raw const ffi::shutdn_.isd).read() }
        }

        /// `TESUB1`: mixture enthalpy of composition `z` at `t` degrees C.
        ///
        /// `ity` selects the correlation: 0 for liquid, 1 for vapour, 2 for
        /// vapour with the ideal-gas correction subtracted.
        pub fn tesub1(&mut self, z: &[f64; 8], t: f64, ity: i32) -> f64 {
            let mut h = 0.0;
            // SAFETY: `z` is exactly the 8 elements the Fortran indexes.
            unsafe { ffi::tesub1_(z.as_ptr(), &t, &mut h, &ity) };
            h
        }

        /// `TESUB2`: temperature from enthalpy, by Newton iteration.
        ///
        /// `t` is both the initial guess and, on failure to converge in 100
        /// iterations, the value returned unchanged. The original reports no
        /// error in that case; quantifying how often it happens is B-0011.
        pub fn tesub2(&mut self, z: &[f64; 8], t_guess: f64, h: f64, ity: i32) -> f64 {
            let mut t = t_guess;
            // SAFETY: as above.
            unsafe { ffi::tesub2_(z.as_ptr(), &mut t, &h, &ity) };
            t
        }

        /// `TESUB3`: heat capacity, the temperature derivative of `TESUB1`.
        pub fn tesub3(&mut self, z: &[f64; 8], t: f64, ity: i32) -> f64 {
            let mut dh = 0.0;
            // SAFETY: as above.
            unsafe { ffi::tesub3_(z.as_ptr(), &t, &mut dh, &ity) };
            dh
        }

        /// `TESUB4`: liquid density of composition `x` at `t` degrees C.
        pub fn tesub4(&mut self, x: &[f64; 8], t: f64) -> f64 {
            let mut r = 0.0;
            // SAFETY: as above.
            unsafe { ffi::tesub4_(x.as_ptr(), &t, &mut r) };
            r
        }

        /// `TESUB5`: build the next cubic walk segment for one channel.
        ///
        /// Consumes three draws (`teprob.f:1528-1530`).
        pub fn tesub5(&mut self, start: WalkSegmentStart, spans: WalkSpans) -> WalkSegment {
            let mut out = WalkSegment {
                adist: 0.0,
                bdist: 0.0,
                cdist: 0.0,
                ddist: 0.0,
                tnext: 0.0,
            };
            // SAFETY: every pointer is to a live local of the right type, and
            // `TESUB5` writes only the six it is declared to write. No arrays
            // are involved, so there is no length to get wrong.
            unsafe {
                ffi::tesub5_(
                    &start.value,
                    &start.slope,
                    &mut out.adist,
                    &mut out.bdist,
                    &mut out.cdist,
                    &mut out.ddist,
                    &start.tlast,
                    &mut out.tnext,
                    &spans.hspan,
                    &spans.hzero,
                    &spans.sspan,
                    &spans.szero,
                    &spans.spspan,
                    &spans.idvflag,
                );
            }
            out
        }

        /// `TESUB6`: one measurement-noise sample of standard deviation `std`.
        ///
        /// Consumes twelve draws from the generator, summing them and
        /// subtracting 6 to approximate a Gaussian.
        pub fn tesub6(&mut self, std: f64) -> f64 {
            let mut x = 0.0;
            // SAFETY: no arrays involved.
            unsafe { ffi::tesub6_(&std, &mut x) };
            x
        }

        /// `TESUB7`: one draw. Negative `i` gives [-1,1), otherwise [0,1).
        pub fn tesub7(&mut self, i: i32) -> f64 {
            // SAFETY: no arrays involved.
            unsafe { ffi::tesub7_(&i) }
        }

        /// `TESUB8`: evaluate disturbance walk channel `i` (1-based) at time
        /// `t`, as a cubic in the time since the channel's last knot.
        pub fn tesub8(&mut self, i: i32, t: f64) -> f64 {
            // SAFETY: no arrays involved.
            unsafe { ffi::tesub8_(&i, &t) }
        }

        /// `COMMON/CTRLALL/`: setpoints and the controller sample time.
        pub fn ctrlall(&mut self) -> Ctrlall {
            // SAFETY: we hold the lock; `Ctrlall` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrlall_).read() }
        }

        /// Set `COMMON/CTRLALL/`.
        pub fn set_ctrlall(&mut self, value: &Ctrlall) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::ctrlall_).write(*value) };
        }

        /// `COMMON/FLAG6/`: the purge override's latch.
        pub fn flag6(&mut self) -> i32 {
            // SAFETY: as above.
            unsafe { (&raw const ffi::flag6_.flag).read() }
        }

        /// Set the purge override's latch.
        pub fn set_flag6(&mut self, flag: i32) {
            // SAFETY: as above.
            unsafe { (&raw mut ffi::flag6_.flag).write(flag) };
        }

        /// Run `CONTRL1` once. `temain_mod.f`.
        pub fn contrl1(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl1_() };
        }

        /// `COMMON/CTRL1/`.
        pub fn ctrl1(&mut self) -> Ctrl1 {
            // SAFETY: as above; `Ctrl1` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl1_).read() }
        }

        /// Set `COMMON/CTRL1/`.
        pub fn set_ctrl1(&mut self, value: &Ctrl1) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl1_).write(*value) };
        }

        /// Run `CONTRL2` once. `temain_mod.f`.
        pub fn contrl2(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl2_() };
        }

        /// `COMMON/CTRL2/`.
        pub fn ctrl2(&mut self) -> Ctrl2 {
            // SAFETY: as above; `Ctrl2` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl2_).read() }
        }

        /// Set `COMMON/CTRL2/`.
        pub fn set_ctrl2(&mut self, value: &Ctrl2) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl2_).write(*value) };
        }

        /// Run `CONTRL3` once. `temain_mod.f`.
        pub fn contrl3(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl3_() };
        }

        /// `COMMON/CTRL3/`.
        pub fn ctrl3(&mut self) -> Ctrl3 {
            // SAFETY: as above; `Ctrl3` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl3_).read() }
        }

        /// Set `COMMON/CTRL3/`.
        pub fn set_ctrl3(&mut self, value: &Ctrl3) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl3_).write(*value) };
        }

        /// Run `CONTRL4` once. `temain_mod.f`.
        pub fn contrl4(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl4_() };
        }

        /// `COMMON/CTRL4/`.
        pub fn ctrl4(&mut self) -> Ctrl4 {
            // SAFETY: as above; `Ctrl4` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl4_).read() }
        }

        /// Set `COMMON/CTRL4/`.
        pub fn set_ctrl4(&mut self, value: &Ctrl4) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl4_).write(*value) };
        }

        /// Run `CONTRL5` once. `temain_mod.f`.
        pub fn contrl5(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl5_() };
        }

        /// `COMMON/CTRL5/`.
        pub fn ctrl5(&mut self) -> Ctrl5 {
            // SAFETY: as above; `Ctrl5` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl5_).read() }
        }

        /// Set `COMMON/CTRL5/`.
        pub fn set_ctrl5(&mut self, value: &Ctrl5) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl5_).write(*value) };
        }

        /// Run `CONTRL6` once. `temain_mod.f`.
        pub fn contrl6(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl6_() };
        }

        /// `COMMON/CTRL6/`.
        pub fn ctrl6(&mut self) -> Ctrl6 {
            // SAFETY: as above; `Ctrl6` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl6_).read() }
        }

        /// Set `COMMON/CTRL6/`.
        pub fn set_ctrl6(&mut self, value: &Ctrl6) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl6_).write(*value) };
        }

        /// Run `CONTRL7` once. `temain_mod.f`.
        pub fn contrl7(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl7_() };
        }

        /// `COMMON/CTRL7/`.
        pub fn ctrl7(&mut self) -> Ctrl7 {
            // SAFETY: as above; `Ctrl7` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl7_).read() }
        }

        /// Set `COMMON/CTRL7/`.
        pub fn set_ctrl7(&mut self, value: &Ctrl7) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl7_).write(*value) };
        }

        /// Run `CONTRL8` once. `temain_mod.f`.
        pub fn contrl8(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl8_() };
        }

        /// `COMMON/CTRL8/`.
        pub fn ctrl8(&mut self) -> Ctrl8 {
            // SAFETY: as above; `Ctrl8` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl8_).read() }
        }

        /// Set `COMMON/CTRL8/`.
        pub fn set_ctrl8(&mut self, value: &Ctrl8) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl8_).write(*value) };
        }

        /// Run `CONTRL9` once. `temain_mod.f`.
        pub fn contrl9(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl9_() };
        }

        /// `COMMON/CTRL9/`.
        pub fn ctrl9(&mut self) -> Ctrl9 {
            // SAFETY: as above; `Ctrl9` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl9_).read() }
        }

        /// Set `COMMON/CTRL9/`.
        pub fn set_ctrl9(&mut self, value: &Ctrl9) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl9_).write(*value) };
        }

        /// Run `CONTRL10` once. `temain_mod.f`.
        pub fn contrl10(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl10_() };
        }

        /// `COMMON/CTRL10/`.
        pub fn ctrl10(&mut self) -> Ctrl10 {
            // SAFETY: as above; `Ctrl10` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl10_).read() }
        }

        /// Set `COMMON/CTRL10/`.
        pub fn set_ctrl10(&mut self, value: &Ctrl10) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl10_).write(*value) };
        }

        /// Run `CONTRL11` once. `temain_mod.f`.
        pub fn contrl11(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl11_() };
        }

        /// `COMMON/CTRL11/`.
        pub fn ctrl11(&mut self) -> Ctrl11 {
            // SAFETY: as above; `Ctrl11` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl11_).read() }
        }

        /// Set `COMMON/CTRL11/`.
        pub fn set_ctrl11(&mut self, value: &Ctrl11) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl11_).write(*value) };
        }

        /// Run `CONTRL13` once. `temain_mod.f`.
        pub fn contrl13(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl13_() };
        }

        /// `COMMON/CTRL13/`.
        pub fn ctrl13(&mut self) -> Ctrl13 {
            // SAFETY: as above; `Ctrl13` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl13_).read() }
        }

        /// Set `COMMON/CTRL13/`.
        pub fn set_ctrl13(&mut self, value: &Ctrl13) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl13_).write(*value) };
        }

        /// Run `CONTRL14` once. `temain_mod.f`.
        pub fn contrl14(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl14_() };
        }

        /// `COMMON/CTRL14/`.
        pub fn ctrl14(&mut self) -> Ctrl14 {
            // SAFETY: as above; `Ctrl14` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl14_).read() }
        }

        /// Set `COMMON/CTRL14/`.
        pub fn set_ctrl14(&mut self, value: &Ctrl14) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl14_).write(*value) };
        }

        /// Run `CONTRL15` once. `temain_mod.f`.
        pub fn contrl15(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl15_() };
        }

        /// `COMMON/CTRL15/`.
        pub fn ctrl15(&mut self) -> Ctrl15 {
            // SAFETY: as above; `Ctrl15` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl15_).read() }
        }

        /// Set `COMMON/CTRL15/`.
        pub fn set_ctrl15(&mut self, value: &Ctrl15) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl15_).write(*value) };
        }

        /// Run `CONTRL16` once. `temain_mod.f`.
        pub fn contrl16(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl16_() };
        }

        /// `COMMON/CTRL16/`.
        pub fn ctrl16(&mut self) -> Ctrl16 {
            // SAFETY: as above; `Ctrl16` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl16_).read() }
        }

        /// Set `COMMON/CTRL16/`.
        pub fn set_ctrl16(&mut self, value: &Ctrl16) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl16_).write(*value) };
        }

        /// Run `CONTRL17` once. `temain_mod.f`.
        pub fn contrl17(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl17_() };
        }

        /// `COMMON/CTRL17/`.
        pub fn ctrl17(&mut self) -> Ctrl17 {
            // SAFETY: as above; `Ctrl17` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl17_).read() }
        }

        /// Set `COMMON/CTRL17/`.
        pub fn set_ctrl17(&mut self, value: &Ctrl17) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl17_).write(*value) };
        }

        /// Run `CONTRL18` once. `temain_mod.f`.
        pub fn contrl18(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl18_() };
        }

        /// `COMMON/CTRL18/`.
        pub fn ctrl18(&mut self) -> Ctrl18 {
            // SAFETY: as above; `Ctrl18` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl18_).read() }
        }

        /// Set `COMMON/CTRL18/`.
        pub fn set_ctrl18(&mut self, value: &Ctrl18) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl18_).write(*value) };
        }

        /// Run `CONTRL19` once. `temain_mod.f`.
        pub fn contrl19(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl19_() };
        }

        /// `COMMON/CTRL19/`.
        pub fn ctrl19(&mut self) -> Ctrl19 {
            // SAFETY: as above; `Ctrl19` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl19_).read() }
        }

        /// Set `COMMON/CTRL19/`.
        pub fn set_ctrl19(&mut self, value: &Ctrl19) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl19_).write(*value) };
        }

        /// Run `CONTRL20` once. `temain_mod.f`.
        pub fn contrl20(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl20_() };
        }

        /// `COMMON/CTRL20/`.
        pub fn ctrl20(&mut self) -> Ctrl20 {
            // SAFETY: as above; `Ctrl20` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl20_).read() }
        }

        /// Set `COMMON/CTRL20/`.
        pub fn set_ctrl20(&mut self, value: &Ctrl20) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl20_).write(*value) };
        }

        /// Run `CONTRL22` once. `temain_mod.f`.
        pub fn contrl22(&mut self) {
            // SAFETY: we hold the lock, and the routine takes no arguments;
            // everything it touches is a `COMMON` block this crate mirrors.
            unsafe { ffi::contrl22_() };
        }

        /// `COMMON/CTRL22/`.
        pub fn ctrl22(&mut self) -> Ctrl22 {
            // SAFETY: as above; `Ctrl22` mirrors the Fortran layout.
            unsafe { (&raw const ffi::ctrl22_).read() }
        }

        /// Set `COMMON/CTRL22/`.
        pub fn set_ctrl22(&mut self, value: &Ctrl22) {
            // SAFETY: as above, writing a correctly typed value.
            unsafe { (&raw mut ffi::ctrl22_).write(*value) };
        }
    }
}
