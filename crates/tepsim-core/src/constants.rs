//! The thermodynamic coefficient tables, `COMMON/CONST/`.
//!
//! Fourteen arrays of eight, one entry per species, set once in `TEINIT` and
//! never written again.
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
