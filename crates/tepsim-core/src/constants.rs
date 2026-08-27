//! The thermodynamic coefficient tables, `COMMON/CONST/`, and the vessel
//! geometry.
//!
//! Fourteen arrays of eight, one entry per species, set once in `TEINIT` and
//! never written again, plus the four fixed vessel volumes at the end of the
//! file. Everything here is a number the original writes down once and then
//! only reads.
//!
//! # Precision is per-constant and is not a detail
//!
//! Fixed-form Fortran treats a literal with no exponent letter as **single**
//! precision. Assigning it to a `DOUBLE PRECISION` variable widens a value that
//! has already been rounded to 24 bits of mantissa, so `XMW(2) = 25.4` at
//! `teprob.f:942` does not store the nearest double to 25.4. It stores
//! 25.399999618530273.
//!
//! Of the 112 assignments in this block, **62 are single precision and 50 carry
//! a `D` suffix**. The error from getting one wrong is about 1.5e-8 relative,
//! against a Tier 1 tolerance of 1e-13: five orders of magnitude, and it would
//! present as a deep numerical problem rather than the clerical one it is.
//!
//! So every single-precision literal is wrapped in [`single`], including the
//! ones that happen to be exactly representable. Uniformity is the point: the
//! wrapper marks what the Fortran says, not what rounding happens to occur, and
//! a reader can check any line against the source without arithmetic.
//!
//! All 24 Antoine coefficients are single precision, which matters most: they
//! feed `exp` at `teprob.f:484` and `teprob.f:487`, so the error is amplified
//! rather than merely carried.
//!
//! # How this table is checked
//!
//! Twice, and neither check trusts the transcription. `tests/constants.rs`
//! reparses `reference/fortran/teprob.f` and independently derives what each
//! value should be, and the oracle test compares every entry bit for bit
//! against the values gfortran actually stored in `COMMON/CONST/`. The second
//! is the decisive one: it is ground truth from the compiler itself.

use crate::component::ByComponent;

/// Widen a single-precision Fortran literal the way the original does.
///
/// Marks every constant whose Fortran literal carries no `D` suffix. See the
/// module documentation for why this is load-bearing rather than cosmetic.
#[must_use]
pub const fn single(value: f32) -> f64 {
    value as f64
}

/// Fold a quotient of two single-precision Fortran literals, in single
/// precision, the way the original does.
///
/// # Why this is not `single(a) / single(b)`
///
/// Fortran types an expression from its operands, not from its destination.
/// `40000.0/1.987` at `teprob.f:503` is `REAL(4) / REAL(4)`, so the *division*
/// happens at single precision and only the result is widened. Writing it as a
/// double division of two widened literals gives a different number:
///
/// | form | value |
/// |---|---|
/// | what gfortran stores | `0x40D3A8B680000000` |
/// | `single(40000.0) / single(1.987)` | `0x40D3A8B670F51E20` |
/// | `40000.0_f64 / 1.987_f64` | `0x40D3A8B66F0ED0F9` |
///
/// That is 4e-9 relative, four orders past the Tier 2 gate, and it lands
/// inside a `DEXP` argument where it is amplified rather than merely carried.
/// The hex above was read out of a gfortran program compiled with this
/// project's pinned flags, not derived from the standard.
///
/// The tell is the mantissa: 29 trailing zero bits is what a widened `f32`
/// looks like, and neither double-precision form has them.
#[must_use]
pub const fn single_quotient(numerator: f32, denominator: f32) -> f64 {
    (numerator / denominator) as f64
}

/// Antoine vapour pressure, constant term.
///
/// `ln(Pvap) = AVP + BVP / (T + CVP)` with `T` in degrees Celsius, evaluated
/// at `teprob.f:484` and `teprob.f:487`. Zero for A, B and C: those three are
/// non-condensible here and are handled as ideal gases instead.
///
/// 0 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:949-956
pub const AVP: ByComponent<f64> = ByComponent::new([
    single(0.0),   // AVP(1) = 0.0
    single(0.0),   // AVP(2) = 0.0
    single(0.0),   // AVP(3) = 0.0
    single(15.92), // AVP(4) = 15.92
    single(16.35), // AVP(5) = 16.35
    single(16.35), // AVP(6) = 16.35
    single(16.43), // AVP(7) = 16.43
    single(17.21), // AVP(8) = 17.21
]);

/// Antoine vapour pressure, reciprocal-temperature term.
///
/// See [`AVP`].
///
/// 0 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:957-964
pub const BVP: ByComponent<f64> = ByComponent::new([
    single(0.0),     // BVP(1) = 0.0
    single(0.0),     // BVP(2) = 0.0
    single(0.0),     // BVP(3) = 0.0
    single(-1444.0), // BVP(4) = -1444.0
    single(-2114.0), // BVP(5) = -2114.0
    single(-2114.0), // BVP(6) = -2114.0
    single(-2748.0), // BVP(7) = -2748.0
    single(-3318.0), // BVP(8) = -3318.0
]);

/// Antoine vapour pressure, temperature offset.
///
/// See [`AVP`].
///
/// 0 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:965-972
pub const CVP: ByComponent<f64> = ByComponent::new([
    single(0.0),   // CVP(1) = 0.0
    single(0.0),   // CVP(2) = 0.0
    single(0.0),   // CVP(3) = 0.0
    single(259.0), // CVP(4) = 259.0
    single(265.5), // CVP(5) = 265.5
    single(265.5), // CVP(6) = 265.5
    single(232.9), // CVP(7) = 232.9
    single(249.6), // CVP(8) = 249.6
]);

/// Liquid heat capacity, constant term.
///
/// Liquid enthalpy is `1.8 * T * (AH + BH*T/2 + CH*T^2/3)` per unit mass,
/// weighted by mole fraction and molecular weight in `TESUB1` with `ITY = 0`.
///
/// 8 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:997-1004
pub const AH: ByComponent<f64> = ByComponent::new([
    1.0e-6,   // AH(1) = 1.0D-6
    1.0e-6,   // AH(2) = 1.0D-6
    1.0e-6,   // AH(3) = 1.0D-6
    0.960e-6, // AH(4) = 0.960D-6
    0.573e-6, // AH(5) = 0.573D-6
    0.652e-6, // AH(6) = 0.652D-6
    0.515e-6, // AH(7) = 0.515D-6
    0.471e-6, // AH(8) = 0.471D-6
]);

/// Liquid heat capacity, linear term.
///
/// See [`AH`].
///
/// 5 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:1005-1012
pub const BH: ByComponent<f64> = ByComponent::new([
    single(0.0), // BH(1) = 0.0
    single(0.0), // BH(2) = 0.0
    single(0.0), // BH(3) = 0.0
    8.70e-9,     // BH(4) = 8.70D-9
    2.41e-9,     // BH(5) = 2.41D-9
    2.18e-9,     // BH(6) = 2.18D-9
    5.65e-10,    // BH(7) = 5.65D-10
    8.70e-10,    // BH(8) = 8.70D-10
]);

/// Liquid heat capacity, quadratic term.
///
/// See [`AH`].
///
/// 5 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:1013-1020
pub const CH: ByComponent<f64> = ByComponent::new([
    single(0.0), // CH(1) = 0.0
    single(0.0), // CH(2) = 0.0
    single(0.0), // CH(3) = 0.0
    4.81e-11,    // CH(4) = 4.81D-11
    1.82e-11,    // CH(5) = 1.82D-11
    1.94e-11,    // CH(6) = 1.94D-11
    3.82e-12,    // CH(7) = 3.82D-12
    2.62e-12,    // CH(8) = 2.62D-12
]);

/// Vapour heat capacity, constant term.
///
/// Vapour enthalpy is `1.8 * T * (AG + BG*T/2 + CG*T^2/3) + AV`, the same
/// shape as the liquid correlation plus the heat of vaporisation, in `TESUB1`
/// with `ITY` non-zero.
///
/// 8 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:1029-1036
pub const AG: ByComponent<f64> = ByComponent::new([
    3.411e-6,  // AG(1) = 3.411D-6
    0.3799e-6, // AG(2) = 0.3799D-6
    0.2491e-6, // AG(3) = 0.2491D-6
    0.3567e-6, // AG(4) = 0.3567D-6
    0.3463e-6, // AG(5) = 0.3463D-6
    0.3930e-6, // AG(6) = 0.3930D-6
    0.170e-6,  // AG(7) = 0.170D-6
    0.150e-6,  // AG(8) = 0.150D-6
]);

/// Vapour heat capacity, linear term.
///
/// See [`AG`].
///
/// 8 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:1037-1044
pub const BG: ByComponent<f64> = ByComponent::new([
    7.18e-10, // BG(1) = 7.18D-10
    1.08e-9,  // BG(2) = 1.08D-9
    1.36e-11, // BG(3) = 1.36D-11
    8.51e-10, // BG(4) = 8.51D-10
    8.96e-10, // BG(5) = 8.96D-10
    1.02e-9,  // BG(6) = 1.02D-9
    0.0e0,    // BG(7) = 0.D0
    0.0e0,    // BG(8) = 0.D0
]);

/// Vapour heat capacity, quadratic term.
///
/// See [`AG`].
///
/// 8 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:1045-1052
pub const CG: ByComponent<f64> = ByComponent::new([
    6.0e-13,   // CG(1) = 6.0D-13
    -3.98e-13, // CG(2) = -3.98D-13
    -3.93e-14, // CG(3) = -3.93D-14
    -3.12e-13, // CG(4) = -3.12D-13
    -3.27e-13, // CG(5) = -3.27D-13
    -3.12e-13, // CG(6) = -3.12D-13
    0.0e0,     // CG(7) = 0.D0
    0.0e0,     // CG(8) = 0.D0
]);

/// Heat of vaporisation.
///
/// Added to the vapour enthalpy in `TESUB1`. Zero for A, B and C.
///
/// 8 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:1021-1028
pub const AV: ByComponent<f64> = ByComponent::new([
    1.0e-6,   // AV(1) = 1.0D-6
    1.0e-6,   // AV(2) = 1.0D-6
    1.0e-6,   // AV(3) = 1.0D-6
    86.7e-6,  // AV(4) = 86.7D-6
    160.0e-6, // AV(5) = 160.D-6
    160.0e-6, // AV(6) = 160.D-6
    225.0e-6, // AV(7) = 225.D-6
    209.0e-6, // AV(8) = 209.D-6
]);

/// Liquid density, constant term.
///
/// Specific volume is `sum over i of X(i) * XMW(i) / (AD + (BD + CD*T) * T)`
/// and density is its reciprocal, in `TESUB4`.
///
/// 0 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:973-980
pub const AD: ByComponent<f64> = ByComponent::new([
    single(1.0),  // AD(1) = 1.0
    single(1.0),  // AD(2) = 1.0
    single(1.0),  // AD(3) = 1.0
    single(23.3), // AD(4) = 23.3
    single(33.9), // AD(5) = 33.9
    single(32.8), // AD(6) = 32.8
    single(49.9), // AD(7) = 49.9
    single(50.5), // AD(8) = 50.5
]);

/// Liquid density, linear term.
///
/// See [`AD`].
///
/// 0 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:981-988
pub const BD: ByComponent<f64> = ByComponent::new([
    single(0.0),     // BD(1) = 0.0
    single(0.0),     // BD(2) = 0.0
    single(0.0),     // BD(3) = 0.0
    single(-0.0700), // BD(4) = -0.0700
    single(-0.0957), // BD(5) = -0.0957
    single(-0.0995), // BD(6) = -0.0995
    single(-0.0191), // BD(7) = -0.0191
    single(-0.0541), // BD(8) = -0.0541
]);

/// Liquid density, quadratic term.
///
/// See [`AD`].
///
/// 0 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:989-996
pub const CD: ByComponent<f64> = ByComponent::new([
    single(0.0),       // CD(1) = 0.0
    single(0.0),       // CD(2) = 0.0
    single(0.0),       // CD(3) = 0.0
    single(-0.0002),   // CD(4) = -0.0002
    single(-0.000152), // CD(5) = -0.000152
    single(-0.000233), // CD(6) = -0.000233
    single(-0.000425), // CD(7) = -0.000425
    single(-0.000150), // CD(8) = -0.000150
]);

/// Molecular weight.
///
/// Used everywhere a molar quantity becomes a mass quantity.
///
/// 0 of 8 carry a `D` suffix; the rest are single precision.
//
// @port teprob.f:941-948
pub const XMW: ByComponent<f64> = ByComponent::new([
    single(2.0),  // XMW(1) = 2.0
    single(25.4), // XMW(2) = 25.4
    single(28.0), // XMW(3) = 28.0
    single(32.0), // XMW(4) = 32.0
    single(46.0), // XMW(5) = 46.0
    single(48.0), // XMW(6) = 48.0
    single(62.0), // XMW(7) = 62.0
    single(76.0), // XMW(8) = 76.0
]);

// ---------------------------------------------------------------------------
// Vessel geometry
//
// Not part of `COMMON/CONST/`: these live in `COMMON/TEPROC/` and are set once
// in `TEINIT`, alongside quantities the model recomputes on every call. They
// are constants in every sense that matters here, and a reader looking for a
// vessel volume looks in this file.
//
// All four literals are written without a `D` suffix and so are single
// precision, like most of the table above. All four happen to be exactly
// representable in 24 bits, so `single` changes nothing numerically; it is
// applied anyway, for the reason given in the module documentation.
// ---------------------------------------------------------------------------

/// Reactor total volume, cubic feet.
///
/// The vapour space is this minus the liquid volume (`teprob.f:473`).
//
// @port teprob.f:1118
pub const VTR: f64 = single(1300.0);

/// Separator total volume, cubic feet.
///
/// The vapour space is this minus the liquid volume (`teprob.f:474`).
//
// @port teprob.f:1119
pub const VTS: f64 = single(3500.0);

/// Stripper total volume, cubic feet.
///
/// The stripper is liquid-filled, so this appears only in the level
/// measurement (`teprob.f:693`), never in a vapour balance.
//
// @port teprob.f:1120
pub const VTC: f64 = single(156.5);

/// Mixing zone total volume, cubic feet.
///
/// The mixing zone is all vapour, so its volume is fixed and its pressure
/// follows the ideal gas law directly (`teprob.f:492`).
//
// @port teprob.f:1121
pub const VTV: f64 = single(5000.0);

// ---------------------------------------------------------------------------
// The nominal initial state
// ---------------------------------------------------------------------------

/// `YY(1..50)` as `TEINIT` sets them: the plant's nominal operating point.
///
/// This is the condition every published `d00`-`d21` run starts from, and the
/// only state in the file that is a *plant* rather than an arbitrary vector.
/// Without it `tepsim-core` cannot start a simulation at all.
///
/// # Precision
///
/// 3 of the 50 carry a `D` suffix: `YY(20)`, `YY(22)` and `YY(24)`, the three
/// written in exponential form because they are near 1e-2. The other 47 are
/// ordinary decimals with no suffix, so they are single precision and are
/// widened.
///
/// That split is not arbitrary: the three exponential literals needed a suffix
/// to be double, and whoever wrote them supplied one, while the 47 plain
/// decimals did not. `YY(43) = 22.21000000` is written to eight decimal places
/// and stored to seven significant figures.
///
/// Getting one wrong shifts the starting point rather than the model, so it
/// would show up as a trajectory that diverges from the published data while
/// every Tier 1 and Tier 2 number stayed perfect. `tests/nominal_state.rs`
/// reparses the Fortran and checks all fifty.
//
// Transcribed digit for digit; `single` does the rounding. See the note on the
// pre-exponentials in `crate::kinetics` for why clippy's advice to shorten
// these is exactly backwards here.
#[allow(
    clippy::excessive_precision,
    reason = "transcribed verbatim from teprob.f; `single` does the rounding"
)]
// @port teprob.f:1053-1102
pub const NOMINAL_STATE: [f64; 50] = [
    single(10.40491389),  // YY(1), teprob.f:1053
    single(4.363996017),  // YY(2), teprob.f:1054
    single(7.570059737),  // YY(3), teprob.f:1055
    single(0.4230042431), // YY(4), teprob.f:1056
    single(24.15513437),  // YY(5), teprob.f:1057
    single(2.942597645),  // YY(6), teprob.f:1058
    single(154.3770655),  // YY(7), teprob.f:1059
    single(159.1865960),  // YY(8), teprob.f:1060
    single(2.808522723),  // YY(9), teprob.f:1061
    single(63.75581199),  // YY(10), teprob.f:1062
    single(26.74026066),  // YY(11), teprob.f:1063
    single(46.38532432),  // YY(12), teprob.f:1064
    single(0.2464521543), // YY(13), teprob.f:1065
    single(15.20484404),  // YY(14), teprob.f:1066
    single(1.852266172),  // YY(15), teprob.f:1067
    single(52.44639459),  // YY(16), teprob.f:1068
    single(41.20394008),  // YY(17), teprob.f:1069
    single(0.5699317760), // YY(18), teprob.f:1070
    single(0.4306056376), // YY(19), teprob.f:1071
    7.9906200783e-03,     // YY(20), teprob.f:1072
    single(0.9056036089), // YY(21), teprob.f:1073
    1.6054258216e-02,     // YY(22), teprob.f:1074
    single(0.7509759687), // YY(23), teprob.f:1075
    8.8582855955e-02,     // YY(24), teprob.f:1076
    single(48.27726193),  // YY(25), teprob.f:1077
    single(39.38459028),  // YY(26), teprob.f:1078
    single(0.3755297257), // YY(27), teprob.f:1079
    single(107.7562698),  // YY(28), teprob.f:1080
    single(29.77250546),  // YY(29), teprob.f:1081
    single(88.32481135),  // YY(30), teprob.f:1082
    single(23.03929507),  // YY(31), teprob.f:1083
    single(62.85848794),  // YY(32), teprob.f:1084
    single(5.546318688),  // YY(33), teprob.f:1085
    single(11.92244772),  // YY(34), teprob.f:1086
    single(5.555448243),  // YY(35), teprob.f:1087
    single(0.9218489762), // YY(36), teprob.f:1088
    single(94.59927549),  // YY(37), teprob.f:1089
    single(77.29698353),  // YY(38), teprob.f:1090
    single(63.05263039),  // YY(39), teprob.f:1091
    single(53.97970677),  // YY(40), teprob.f:1092
    single(24.64355755),  // YY(41), teprob.f:1093
    single(61.30192144),  // YY(42), teprob.f:1094
    single(22.21000000),  // YY(43), teprob.f:1095
    single(40.06374673),  // YY(44), teprob.f:1096
    single(38.10034370),  // YY(45), teprob.f:1097
    single(46.53415582),  // YY(46), teprob.f:1098
    single(47.44573456),  // YY(47), teprob.f:1099
    single(41.10581288),  // YY(48), teprob.f:1100
    single(18.11349055),  // YY(49), teprob.f:1101
    single(50.00000000),  // YY(50), teprob.f:1102
];
