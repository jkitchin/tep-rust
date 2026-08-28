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

// ---------------------------------------------------------------------------
// The closed-loop driver, `reference/fortran/temain_mod.f`
//
// The nineteen active controllers and the twentieth that is never called. Each
// is a parameterless `SUBROUTINE` communicating entirely through `COMMON`,
// which makes the binding trivial and the *state* the whole difficulty.
// ---------------------------------------------------------------------------

/// `COMMON/CTRLALL/`: the twenty setpoints and the controller sample time.
///
/// `SETPT` is indexed one-based in the Fortran; slot 12 belongs to the dead
/// `CONTRL22` and slot 21 is unused.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ctrlall {
    /// `SETPT(1..20)`, one-based in the Fortran.
    pub setpt: [f64; 20],
    /// `DELTAT`: the plant step, in hours.
    pub deltat: f64,
}

/// `COMMON/FLAG6/`: the purge override's latch.
///
/// Zero when the PI loop is running, 1 while the valve is latched open, 2
/// while it is latched shut. See `temain_mod.f:710-731`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Flag6 {
    /// `FLAG`: 0 running, 1 latched open, 2 latched shut.
    pub flag: i32,
}

/// `COMMON/CTRL1/`: controller 1's tuning and error history.
///
/// proportional-only; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl1 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL2/`: controller 2's tuning and error history.
///
/// proportional-only; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl2 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL3/`: controller 3's tuning and error history.
///
/// proportional-only; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl3 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL4/`: controller 4's tuning and error history.
///
/// proportional-only; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl4 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL5/`: controller 5's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl5 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL6/`: controller 6's tuning and error history.
///
/// proportional-only; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl6 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL7/`: controller 7's tuning and error history.
///
/// proportional-only; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl7 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL8/`: controller 8's tuning and error history.
///
/// proportional-only; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl8 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL9/`: controller 9's tuning and error history.
///
/// proportional-only; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl9 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL10/`: controller 10's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl10 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL11/`: controller 11's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl11 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL13/`: controller 13's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl13 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL14/`: controller 14's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl14 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL15/`: controller 15's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl15 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL16/`: controller 16's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl16 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL17/`: controller 17's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl17 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL18/`: controller 18's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl18 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL19/`: controller 19's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl19 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL20/`: controller 20's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl20 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
}

/// `COMMON/CTRL22/`: controller 22's tuning and error history.
///
/// PI; the presence of `taui` is what distinguishes the two shapes.
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ctrl22 {
    /// `GAINn`: controller gain, in output units per percent of span.
    pub gain: f64,
    /// `TAUIn`: reset time, in hours.
    pub taui: f64,
    /// `ERROLDn`: the error at the previous call.
    pub errold: f64,
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

/// `COMMON/TEPROC/`, the plant's entire working set, from `teprob.f:243-271`.
///
/// 580 doubles followed by 12 integers, laid out exactly in declaration order.
/// The trailing `IVST` integers after a long run of doubles are the part most
/// likely to be got wrong, so [`Teproc::LEN_BYTES`] pins the total and the
/// layout tests read fields from both ends of the block.
///
/// # Two-dimensional arrays are column-major
///
/// Fortran stores `FCM(8,13)` with the first index varying fastest, so
/// `FCM(i, j)` lives at `fcm[j - 1][i - 1]` here: the *outer* Rust index is the
/// stream and the *inner* one is the component. Getting this backwards would
/// still typecheck and still have the right size.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
// Every field is the Fortran name from the COMMON declaration cited above,
// lowercased. Per-field doc comments would restate the name and nothing else;
// the declaration and the unit-operation docs are the real documentation.
#[allow(missing_docs)]
pub struct Teproc {
    pub uclr: [f64; 8],
    pub ucvr: [f64; 8],
    pub utlr: f64,
    pub utvr: f64,
    pub xlr: [f64; 8],
    pub xvr: [f64; 8],
    pub etr: f64,
    pub esr: f64,
    pub tcr: f64,
    pub tkr: f64,
    pub dlr: f64,
    pub vlr: f64,
    pub vvr: f64,
    pub vtr: f64,
    pub ptr: f64,
    pub ppr: [f64; 8],
    pub crxr: [f64; 8],
    pub rr: [f64; 4],
    pub rh: f64,
    pub fwr: f64,
    pub twr: f64,
    pub qur: f64,
    pub hwr: f64,
    pub uar: f64,
    pub ucls: [f64; 8],
    pub ucvs: [f64; 8],
    pub utls: f64,
    pub utvs: f64,
    pub xls: [f64; 8],
    pub xvs: [f64; 8],
    pub ets: f64,
    pub ess: f64,
    pub tcs: f64,
    pub tks: f64,
    pub dls: f64,
    pub vls: f64,
    pub vvs: f64,
    pub vts: f64,
    pub pts: f64,
    pub pps: [f64; 8],
    pub fws: f64,
    pub tws: f64,
    pub qus: f64,
    pub hws: f64,
    pub uclc: [f64; 8],
    pub utlc: f64,
    pub xlc: [f64; 8],
    pub etc: f64,
    pub esc: f64,
    pub tcc: f64,
    pub dlc: f64,
    pub vlc: f64,
    pub vtc: f64,
    pub quc: f64,
    pub ucvv: [f64; 8],
    pub utvv: f64,
    pub xvv: [f64; 8],
    pub etv: f64,
    pub esv: f64,
    pub tcv: f64,
    pub tkv: f64,
    pub vtv: f64,
    pub ptv: f64,
    pub vcv: [f64; 12],
    pub vrng: [f64; 12],
    pub vtau: [f64; 12],
    pub ftm: [f64; 13],
    /// `FCM(8,13)`: `fcm[stream][component]`.
    pub fcm: [[f64; 8]; 13],
    /// `XST(8,13)`: `xst[stream][component]`.
    pub xst: [[f64; 8]; 13],
    pub xmws: [f64; 13],
    pub hst: [f64; 13],
    pub tst: [f64; 13],
    pub sfr: [f64; 8],
    pub cpflmx: f64,
    pub cpprmx: f64,
    pub cpdh: f64,
    pub tcwr: f64,
    pub tcws: f64,
    pub htr: [f64; 3],
    pub agsp: f64,
    pub xdel: [f64; 41],
    pub xns: [f64; 41],
    pub tgas: f64,
    pub tprod: f64,
    pub vst: [f64; 12],
    pub ivst: [i32; 12],
}

impl Teproc {
    /// 580 doubles then 12 integers, counted from the Fortran declaration.
    pub const LEN_BYTES: usize = 580 * 8 + 12 * 4;
}

/// `COMMON/WLK/`, the disturbance random-walk state, from `teprob.f:285-297`.
///
/// Eleven arrays of 12 doubles, then 12 integers. Channels 1-9 are driven by
/// `TESUB5`; channels 10-12 use the different, spikier rule in the `DO 910`
/// loop at `teprob.f:372`.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
// Every field is the Fortran name from the COMMON declaration cited above,
// lowercased. Per-field doc comments would restate the name and nothing else;
// the declaration and the unit-operation docs are the real documentation.
#[allow(missing_docs)]
pub struct Wlk {
    pub adist: [f64; 12],
    pub bdist: [f64; 12],
    pub cdist: [f64; 12],
    pub ddist: [f64; 12],
    pub tlast: [f64; 12],
    pub tnext: [f64; 12],
    pub hspan: [f64; 12],
    pub hzero: [f64; 12],
    pub sspan: [f64; 12],
    pub szero: [f64; 12],
    pub spspan: [f64; 12],
    pub idvwlk: [i32; 12],
}

impl Wlk {
    /// Eleven arrays of 12 doubles, then 12 integers.
    pub const LEN_BYTES: usize = 11 * 12 * 8 + 12 * 4;
}

/// `COMMON/CONST/`, the thermodynamic coefficients, from `teprob.f:305-311`.
///
/// Fourteen arrays of 8, one entry per component A through H. Antoine vapour
/// pressure (`AVP`, `BVP`, `CVP`), liquid enthalpy (`AH`, `BH`, `CH`), vapour
/// enthalpy (`AG`, `BG`, `CG`), heat of vaporisation (`AV`), liquid density
/// (`AD`, `BD`, `CD`), and molecular weight (`XMW`).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq)]
// Every field is the Fortran name from the COMMON declaration cited above,
// lowercased. Per-field doc comments would restate the name and nothing else;
// the declaration and the unit-operation docs are the real documentation.
#[allow(missing_docs)]
pub struct Const {
    pub avp: [f64; 8],
    pub bvp: [f64; 8],
    pub cvp: [f64; 8],
    pub ah: [f64; 8],
    pub bh: [f64; 8],
    pub ch: [f64; 8],
    pub ag: [f64; 8],
    pub bg: [f64; 8],
    pub cg: [f64; 8],
    pub av: [f64; 8],
    pub ad: [f64; 8],
    pub bd: [f64; 8],
    pub cd: [f64; 8],
    pub xmw: [f64; 8],
}

impl Const {
    /// Fourteen arrays of 8 doubles.
    pub const LEN_BYTES: usize = 14 * 8 * 8;
}

/// `COMMON/SHUTDN/ ISD`, present only after build-time instrumentation.
///
/// The original keeps `ISD` as a local in `TEFUNC`, so there is nothing to link
/// against until `build.rs` hoists it. See `instrument.rs`.
#[repr(C)]
pub(crate) struct Shutdn {
    pub isd: i32,
}

/// How many draws `COMMON/RNGTRC/` can hold. Must match `TRCCAP` in
/// `instrument.rs`.
///
/// B-0027 measured at most 522 draws in one evaluation, so 4096 is eight times
/// the worst case seen. The counter keeps counting past it, so an overflow is
/// reported rather than silently truncating.
pub const TRACE_CAPACITY: usize = 4096;

/// `COMMON/RNGTRC/`: the Tier 3 draw trace.
///
/// Does not exist in the vendored Fortran; `instrument.rs` adds it. See that
/// file for why recording cannot change the numbers.
#[repr(C)]
pub(crate) struct Rngtrc {
    pub value: [f64; TRACE_CAPACITY],
    pub sign: [i32; TRACE_CAPACITY],
    pub count: i32,
}

unsafe extern "C" {
    pub(crate) static mut teproc_: Teproc;
    pub(crate) static mut wlk_: Wlk;
    pub(crate) static mut const_: Const;
    pub(crate) static mut shutdn_: Shutdn;
    pub(crate) static mut rngtrc_: Rngtrc;
    pub(crate) static mut ctrlall_: Ctrlall;
    pub(crate) static mut flag6_: Flag6;
    pub(crate) static mut ctrl1_: Ctrl1;
    /// `SUBROUTINE CONTRL1`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl1_();
    pub(crate) static mut ctrl2_: Ctrl2;
    /// `SUBROUTINE CONTRL2`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl2_();
    pub(crate) static mut ctrl3_: Ctrl3;
    /// `SUBROUTINE CONTRL3`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl3_();
    pub(crate) static mut ctrl4_: Ctrl4;
    /// `SUBROUTINE CONTRL4`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl4_();
    pub(crate) static mut ctrl5_: Ctrl5;
    /// `SUBROUTINE CONTRL5`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl5_();
    pub(crate) static mut ctrl6_: Ctrl6;
    /// `SUBROUTINE CONTRL6`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl6_();
    pub(crate) static mut ctrl7_: Ctrl7;
    /// `SUBROUTINE CONTRL7`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl7_();
    pub(crate) static mut ctrl8_: Ctrl8;
    /// `SUBROUTINE CONTRL8`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl8_();
    pub(crate) static mut ctrl9_: Ctrl9;
    /// `SUBROUTINE CONTRL9`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl9_();
    pub(crate) static mut ctrl10_: Ctrl10;
    /// `SUBROUTINE CONTRL10`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl10_();
    pub(crate) static mut ctrl11_: Ctrl11;
    /// `SUBROUTINE CONTRL11`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl11_();
    pub(crate) static mut ctrl13_: Ctrl13;
    /// `SUBROUTINE CONTRL13`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl13_();
    pub(crate) static mut ctrl14_: Ctrl14;
    /// `SUBROUTINE CONTRL14`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl14_();
    pub(crate) static mut ctrl15_: Ctrl15;
    /// `SUBROUTINE CONTRL15`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl15_();
    pub(crate) static mut ctrl16_: Ctrl16;
    /// `SUBROUTINE CONTRL16`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl16_();
    pub(crate) static mut ctrl17_: Ctrl17;
    /// `SUBROUTINE CONTRL17`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl17_();
    pub(crate) static mut ctrl18_: Ctrl18;
    /// `SUBROUTINE CONTRL18`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl18_();
    pub(crate) static mut ctrl19_: Ctrl19;
    /// `SUBROUTINE CONTRL19`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl19_();
    pub(crate) static mut ctrl20_: Ctrl20;
    /// `SUBROUTINE CONTRL20`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl20_();
    pub(crate) static mut ctrl22_: Ctrl22;
    /// `SUBROUTINE CONTRL22`: no arguments; reads and writes `COMMON`.
    pub(crate) fn contrl22_();

    /// `SUBROUTINE TESUB1(Z, T, H, ITY)`: mixture enthalpy. `teprob.f:1376`.
    pub(crate) fn tesub1_(z: *const f64, t: *const f64, h: *mut f64, ity: *const i32);
    /// `SUBROUTINE TESUB2(Z, T, H, ITY)`: temperature from enthalpy, by Newton.
    /// `teprob.f:1416`. Silently returns `T` unchanged if it fails to converge.
    pub(crate) fn tesub2_(z: *const f64, t: *mut f64, h: *const f64, ity: *const i32);
    /// `SUBROUTINE TESUB3(Z, T, DH, ITY)`: heat capacity. `teprob.f:1444`.
    pub(crate) fn tesub3_(z: *const f64, t: *const f64, dh: *mut f64, ity: *const i32);
    /// `SUBROUTINE TESUB4(X, T, R)`: liquid density. `teprob.f:1482`.
    pub(crate) fn tesub4_(x: *const f64, t: *const f64, r: *mut f64);
    /// `SUBROUTINE TESUB5(S, SP, ADIST, BDIST, CDIST, DDIST, TLAST, TNEXT,
    /// HSPAN, HZERO, SSPAN, SZERO, SPSPAN, IDVFLAG)`: build the next cubic
    /// walk segment. `teprob.f:1506-1537`.
    ///
    /// Fourteen arguments, and Fortran passes every one by reference, so the
    /// six that the routine writes (`ADIST` through `TNEXT`) are `*mut` here
    /// while the rest are `*const`. `S` and `SP` are the segment's starting
    /// value and slope; the five span parameters and the flag are read only.
    ///
    /// Consumes *three* draws, at `teprob.f:1528-1530`.
    #[allow(
        clippy::too_many_arguments,
        reason = "the Fortran signature has fourteen; wrapping it is `tesub5`"
    )]
    pub(crate) fn tesub5_(
        s: *const f64,
        sp: *const f64,
        adist: *mut f64,
        bdist: *mut f64,
        cdist: *mut f64,
        ddist: *mut f64,
        tlast: *const f64,
        tnext: *mut f64,
        hspan: *const f64,
        hzero: *const f64,
        sspan: *const f64,
        szero: *const f64,
        spspan: *const f64,
        idvflag: *const i32,
    );
    /// `SUBROUTINE TESUB6(STD, X)`: Gaussian-ish noise, twelve uniforms summed.
    /// `teprob.f:1539`.
    pub(crate) fn tesub6_(std: *const f64, x: *mut f64);
    /// `DOUBLE PRECISION FUNCTION TESUB7(I)`: the generator. `teprob.f:1547-1555`.
    /// Negative `i` gives [-1,1), non-negative gives [0,1).
    pub(crate) fn tesub7_(i: *const i32) -> f64;
    /// `DOUBLE PRECISION FUNCTION TESUB8(I, T)`: evaluate walk channel `i` at
    /// time `t` as a cubic. `teprob.f:1557`.
    pub(crate) fn tesub8_(i: *const i32, t: *const f64) -> f64;
}
