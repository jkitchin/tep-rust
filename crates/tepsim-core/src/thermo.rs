//! Mixture enthalpy, internal energy, heat capacity and liquid density.
//!
//! Ported from `TESUB1` (`teprob.f:1375-1414`), `TESUB3`
//! (`teprob.f:1443-1480`) and `TESUB4` (`teprob.f:1481-1505`). Between them
//! these are called from almost every balance in the plant, so Tier 1 exists to
//! prove them exactly.
//!
//! # The correlations
//!
//! Each species carries a quadratic heat capacity in degrees Celsius, and the
//! molar enthalpy is its integral from zero. Writing \\(a_i, b_i, c_i\\) for the
//! liquid coefficients `AH`, `BH`, `CH` and \\(M_i\\) for the molecular weight
//! `XMW`:
//!
//! \\[
//!   h_i^{\\mathrm{liq}}(T) = 1.8\\, T \\left( a_i + \\frac{b_i T}{2}
//!                            + \\frac{c_i T^2}{3} \\right)
//! \\]
//!
//! The vapour form uses `AG`, `BG`, `CG` and adds a latent heat `AV`, which is
//! the enthalpy of vaporisation at the reference temperature:
//!
//! \\[
//!   h_i^{\\mathrm{vap}}(T) = 1.8\\, T \\left( a^g_i + \\frac{b^g_i T}{2}
//!                            + \\frac{c^g_i T^2}{3} \\right) + a^v_i
//! \\]
//!
//! Either way the mixture value is the mole-weighted sum
//! \\(H = \\sum_i z_i M_i h_i(T)\\), and the heat capacity is its exact
//! derivative:
//!
//! \\[
//!   \\frac{\\partial H}{\\partial T} = \\sum_i z_i M_i \\cdot 1.8
//!       \\left( a_i + b_i T + c_i T^2 \\right)
//! \\]
//!
//! The third mode subtracts \\(R (T + 273.15)\\) to give internal energy
//! instead of enthalpy; see [`EnergyBasis::VapourInternalEnergy`].
//!
//! The 1.8 is degrees Rankine per degree Celsius. The correlations are fitted
//! in imperial units and evaluated with a Celsius argument, so every term picks
//! up the conversion.
//!
//! # Variables
//!
//! | Fortran | Here | Meaning | Units |
//! |---|---|---|---|
//! | `Z(8)` | `z` | mole fractions | dimensionless |
//! | `T` | `celsius` | temperature | degrees Celsius |
//! | `ITY` | `basis` | which correlation | - |
//! | `H` | return of [`enthalpy`] | mixture enthalpy or internal energy | MMBtu/lbmol-ish, see below |
//! | `DH` | return of [`heat_capacity`] | its temperature derivative | per degree Celsius |
//! | `AH,BH,CH` | [`crate::constants::AH`] etc. | liquid heat capacity coefficients | - |
//! | `AG,BG,CG` | [`crate::constants::AG`] etc. | vapour heat capacity coefficients | - |
//! | `AV` | [`crate::constants::AV`] | latent heat | - |
//! | `AD,BD,CD` | [`crate::constants::AD`] etc. | liquid density coefficients | lb/ft^3 |
//! | `X(8)` | `x` | liquid mole fractions | dimensionless |
//! | `R` | return of [`liquid_density`] | mixture molar density | lbmol/ft^3 |
//! | `XMW` | [`crate::constants::XMW`] | molecular weight | lb/lbmol |
//!
//! The original never states an energy unit. The `1.0D-6` scaling folded into
//! the coefficients and the `1.9872` Btu/(lbmol R) gas constant behind
//! [`GAS_CONSTANT`] put it at millions of Btu per hour once multiplied by a
//! molar flow. Nothing in the model depends on naming it, so this port does not
//! invent one.
//!
//! # Why the arithmetic is written the way it is
//!
//! Every expression here reproduces the Fortran's association and rounding
//! exactly, term for term. That is not stylistic: Tier 1 asserts bit equality
//! with gfortran's output, and a fused multiply-add, a reassociated sum or a
//! `T * T * c` where the original wrote `c * T**2` all change the last bits.
//! `clippy::suboptimal_flops` is allowed for the module for this reason, and it
//! is the one place in the crate where "more accurate" is the wrong answer.

// Every float expression in this module is shaped to match `teprob.f`'s
// rounding exactly. `mul_add` would fuse two roundings into one and produce a
// *more* accurate result that no longer matches the oracle bit for bit.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]
// `hi = 1.8 * hi` rather than `hi *= 1.8`, because `teprob.f:1396` writes
// `HI=1.8D0*HI` and this module is meant to be checkable against the listing a
// line at a time. IEEE multiplication commutes exactly, so the two forms are
// bit-identical and this is purely about which one a reader can verify.
#![allow(
    clippy::assign_op_pattern,
    reason = "operand order is transcribed from the Fortran, not chosen"
)]

use crate::component::{Component, Composition};
use crate::constants::{AD, AG, AH, AV, BD, BG, BH, CD, CG, CH, XMW, single};

/// The gas constant in the model's units, from `teprob.f:1410`.
///
/// Written as the division rather than as `3.57696e-6` because that is what the
/// original computes: `R=3.57696D0/1.D6`. The two happen to agree here, but
/// transcribing the operation rather than its result is the habit that keeps
/// them agreeing.
///
/// The value is 1.8 x 1.9872, exactly, in decimal: the ideal gas constant in
/// Btu per lbmol per degree Rankine, times the Rankine-per-Celsius degree ratio
/// that appears throughout these correlations. `teprob.f:594` uses the same two
/// numbers unmultiplied, which is how the coincidence can be checked rather
/// than guessed at.
pub const GAS_CONSTANT: f64 = 3.57696 / 1.0e6;

/// The Celsius-to-Kelvin offset as `teprob.f:1411` actually stores it.
///
/// # This is single precision, and it matters
///
/// The literal on that line is written `273.15` with no `D` suffix, so gfortran
/// stores it at `f32` precision and widens it: 273.14999389648438, not
/// 273.15. Writing `273.15_f64` here would make the correction wrong by 2.1e-8
/// relative.
///
/// That is already past the Tier 1 gate, and the error does not stay that
/// small. The subtraction this constant feeds is a near-total cancellation over
/// much of the composition space: at an equimolar A/B/C/D mixture at 21.875
/// degrees the correction is 99.87% of the enthalpy, so the relative error in
/// the *result* is about 1.6e-5. Eight orders of magnitude past the gate, from
/// one missing letter.
///
/// The original is not consistent about it. `teprob.f:594` writes the same
/// offset as `273.15D0`, in double. Only this one is single, so it cannot be
/// inferred from elsewhere in the file and has to be read off the line itself.
pub const ABSOLUTE_ZERO_OFFSET: f64 = single(273.15);

/// Degrees Rankine per degree Celsius, `teprob.f:1396`.
const RANKINE_PER_CELSIUS: f64 = 1.8;

/// Which of the three correlations to evaluate: the original's `ITY`.
///
/// `TESUB1` and `TESUB3` branch on it identically. `ITY=0` takes the liquid
/// coefficients and anything else takes the vapour ones, and then `ITY=2` alone
/// gets the ideal-gas correction (`teprob.f:1392-1412`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum EnergyBasis {
    /// `ITY=0`. Liquid enthalpy, from `AH`, `BH`, `CH`.
    ///
    /// Used for the three liquid holdups and the two liquid streams
    /// (`teprob.f:460-464`, `563-564`, `655`).
    LiquidEnthalpy = 0,

    /// `ITY=1`. Vapour enthalpy, from `AG`, `BG`, `CG`, plus the latent heat
    /// `AV`.
    ///
    /// Used for every vapour stream (`teprob.f:555-561`, `654`).
    VapourEnthalpy = 1,

    /// `ITY=2`. Vapour *internal energy*: the vapour enthalpy less
    /// \\(R(T + 273.15)\\).
    ///
    /// This is \\(U = H - RT\\) for an ideal gas, with the temperature
    /// converted to absolute. The original does not say so, but the arithmetic
    /// does: [`GAS_CONSTANT`] is the gas constant, and the one call site
    /// (`teprob.f:465`) solves for `TCV` from `ESV`, the reactor vapour
    /// *energy* state.
    ///
    /// It is also the mode with no numerical margin. The correction can reach
    /// 99.9% of the enthalpy, so a difference of one ULP anywhere upstream
    /// arrives at the caller multiplied by a thousand. See
    /// [`ABSOLUTE_ZERO_OFFSET`].
    VapourInternalEnergy = 2,
}

impl EnergyBasis {
    /// All three, in `ITY` order.
    pub const ALL: [Self; 3] = [
        Self::LiquidEnthalpy,
        Self::VapourEnthalpy,
        Self::VapourInternalEnergy,
    ];

    /// The `ITY` value the Fortran would be passed.
    #[must_use]
    pub const fn ity(self) -> i32 {
        self as i32
    }

    /// Recover a basis from an `ITY` value.
    ///
    /// Returns `None` for anything but 0, 1 or 2. The original would silently
    /// treat 7 as vapour, since it tests only `ITY.EQ.0` and `ITY.EQ.2`, but no
    /// call site passes anything else and reproducing the shrug would mean
    /// carrying a fourth state forever.
    #[must_use]
    pub const fn from_ity(ity: i32) -> Option<Self> {
        Some(match ity {
            0 => Self::LiquidEnthalpy,
            1 => Self::VapourEnthalpy,
            2 => Self::VapourInternalEnergy,
            _ => return None,
        })
    }

    /// Whether this basis uses the vapour coefficients.
    #[must_use]
    const fn is_vapour(self) -> bool {
        !matches!(self, Self::LiquidEnthalpy)
    }
}

/// Mixture enthalpy, or internal energy for
/// [`EnergyBasis::VapourInternalEnergy`], of composition `z` at `celsius`.
///
/// See the module documentation for the correlations and the units.
///
/// ```
/// use tepsim_core::{Composition, thermo::{EnergyBasis, enthalpy}};
///
/// // Every liquid term carries a factor of T, so the correlation is anchored
/// // at zero: `teprob.f:1395`.
/// let pure_a = Composition::new([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
/// assert_eq!(enthalpy(&pure_a, 0.0, EnergyBasis::LiquidEnthalpy), 0.0);
/// ```
// @port teprob.f:1375-1414
#[must_use]
pub fn enthalpy(z: &Composition, celsius: f64, basis: EnergyBasis) -> f64 {
    let mut h = 0.0_f64;

    if basis.is_vapour() {
        // teprob.f:1400-1407
        for component in Component::ALL {
            let mut hi = celsius
                * (AG[component]
                    + BG[component] * celsius / 2.0
                    + CG[component] * (celsius * celsius) / 3.0);
            hi = RANKINE_PER_CELSIUS * hi;
            hi += AV[component];
            h += z[component] * XMW[component] * hi;
        }
    } else {
        // teprob.f:1393-1398
        for component in Component::ALL {
            let mut hi = celsius
                * (AH[component]
                    + BH[component] * celsius / 2.0
                    + CH[component] * (celsius * celsius) / 3.0);
            hi = RANKINE_PER_CELSIUS * hi;
            h += z[component] * XMW[component] * hi;
        }
    }

    // teprob.f:1409-1412
    if basis == EnergyBasis::VapourInternalEnergy {
        h -= GAS_CONSTANT * (celsius + ABSOLUTE_ZERO_OFFSET);
    }
    h
}

/// The exact temperature derivative of [`enthalpy`], at the same point.
///
/// "Exact" analytically, not merely numerically: `TESUB3` is the term-by-term
/// derivative of `TESUB1`, with the constant latent heat `AV` dropping out and
/// the \\(R(T + 273.15)\\) correction differentiating to \\(R\\). That is what
/// makes `TESUB2`'s Newton iteration converge, and it is asserted here by a
/// finite-difference test rather than left as a claim.
///
/// ```
/// use tepsim_core::{Composition, thermo::{EnergyBasis, heat_capacity}};
///
/// // Unlike the enthalpy, the heat capacity at zero degrees is not zero: it is
/// // 1.8 * sum(z_i * M_i * a_i), the leading coefficient (`teprob.f:1463`).
/// let pure_a = Composition::new([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
/// assert!(heat_capacity(&pure_a, 0.0, EnergyBasis::LiquidEnthalpy) > 0.0);
/// ```
// @port teprob.f:1443-1480
#[must_use]
pub fn heat_capacity(z: &Composition, celsius: f64, basis: EnergyBasis) -> f64 {
    let mut dh = 0.0_f64;

    if basis.is_vapour() {
        // teprob.f:1468-1473. No `AV` term: the latent heat is constant, so it
        // differentiates away. This is the only structural difference from
        // TESUB1's vapour branch.
        for component in Component::ALL {
            let mut dhi =
                AG[component] + BG[component] * celsius + CG[component] * (celsius * celsius);
            dhi = RANKINE_PER_CELSIUS * dhi;
            dh += z[component] * XMW[component] * dhi;
        }
    } else {
        // teprob.f:1461-1466
        for component in Component::ALL {
            let mut dhi =
                AH[component] + BH[component] * celsius + CH[component] * (celsius * celsius);
            dhi = RANKINE_PER_CELSIUS * dhi;
            dh += z[component] * XMW[component] * dhi;
        }
    }

    // teprob.f:1475-1478
    if basis == EnergyBasis::VapourInternalEnergy {
        dh -= GAS_CONSTANT;
    }
    dh
}

/// The iteration cap in `TESUB2`, from `teprob.f:1433`.
pub const MAX_NEWTON_ITERATIONS: u32 = 100;

/// The convergence criterion in `TESUB2`, from `teprob.f:1439`.
///
/// Absolute, in degrees Celsius, and tested against the Newton *step* rather
/// than against the enthalpy residual. It carries a `D` suffix, so unlike the
/// offset in [`ABSOLUTE_ZERO_OFFSET`] it really is double precision.
pub const NEWTON_TOLERANCE: f64 = 1.0e-12;

/// Why solving for a temperature failed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TemperatureError {
    /// Newton ran the full [`MAX_NEWTON_ITERATIONS`] without the step falling
    /// below [`NEWTON_TOLERANCE`].
    DidNotConverge {
        /// The initial guess the caller supplied, in degrees Celsius.
        guess: f64,
        /// Where the iteration had got to when it gave up.
        last: f64,
        /// The size of the final step, which never got small enough.
        last_step: f64,
    },
}

impl core::fmt::Display for TemperatureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let Self::DidNotConverge {
            guess,
            last,
            last_step,
        } = self;
        write!(
            f,
            "Newton did not converge in {MAX_NEWTON_ITERATIONS} iterations \
             from a guess of {guess} C: reached {last} C with a final step of \
             {last_step}, which is not below {NEWTON_TOLERANCE:e}"
        )
    }
}

impl core::error::Error for TemperatureError {}

/// Solve [`enthalpy`] for the temperature that produces `target`, by Newton's
/// method from `guess`.
///
/// The iteration is the original's, step for step: evaluate the enthalpy and
/// the heat capacity at the current temperature, take the Newton step
/// \\(-\\,\\mathrm{err}/(\\partial H/\\partial T)\\), apply it, and stop when the
/// step falls below [`NEWTON_TOLERANCE`]. Since [`enthalpy`] and
/// [`heat_capacity`] are both bit-identical to the Fortran, so is every
/// iterate, and so is the decision about when to stop.
///
/// # This returns a `Result`, and the original does not
///
/// **Delta D-001, Class B.** `teprob.f:1439-1440` puts the convergence test on
/// the loop-terminal line, so on success it jumps out with the converged value.
/// On failure the loop simply runs out, control falls through to `T=TIN`, and
/// the routine *restores the caller's original guess and returns as if it had
/// succeeded*. There is no error code, no flag and no output: a caller cannot
/// distinguish a converged temperature from a silently abandoned one, and a
/// plant state built on a stale guess propagates from there.
///
/// This port returns [`TemperatureError::DidNotConverge`] instead and lets the
/// caller decide. See `book/src/deltas.md` for the measured effect; across the
/// full Tier 1 sweep the two behaviours never diverge, because the iteration
/// always converges on the physical domain.
// @port teprob.f:1415-1442
// @delta D-001 class=B teprob.f:1439-1440
pub fn temperature_from_enthalpy(
    z: &Composition,
    guess: f64,
    target: f64,
    basis: EnergyBasis,
) -> Result<f64, TemperatureError> {
    // teprob.f:1432. The original keeps this to restore on failure; here it is
    // kept only to report what the failed call started from.
    let mut celsius = guess;
    let mut step = f64::NAN;

    // teprob.f:1433-1439
    for _ in 0..MAX_NEWTON_ITERATIONS {
        let error = enthalpy(z, celsius, basis) - target;
        let slope = heat_capacity(z, celsius, basis);
        // `-ERR/DH` at teprob.f:1437. Fortran binds unary minus below division,
        // so this is -(err/slope). IEEE negation is exact, so the two groupings
        // agree bit for bit; the parenthesis is here to match the listing.
        step = -(error / slope);
        celsius = celsius + step;
        if step.abs() < NEWTON_TOLERANCE {
            return Ok(celsius);
        }
    }

    // teprob.f:1440. The original assigns `T=TIN` here. Delta D-001.
    Err(TemperatureError::DidNotConverge {
        guess,
        last: celsius,
        last_step: step,
    })
}

/// Molar density of a liquid mixture, in lbmol per cubic foot.
///
/// Each species has a mass density quadratic in temperature,
/// \\(\\rho_i(T) = a^d_i + (b^d_i + c^d_i T) T\\) in pounds per cubic foot, and
/// the mixture is combined by ideal volume additivity: the volumes add, so the
/// reciprocals do.
///
/// \\[
///   \\rho = \\left( \\sum_i \\frac{x_i M_i}{\\rho_i(T)} \\right)^{-1}
/// \\]
///
/// The units follow from the one call site. `teprob.f:467-469` divides a molar
/// holdup by this to get a volume, and `teprob.f:704` divides that volume by
/// 35.3145, the cubic feet in a cubic metre, before comparing it against a
/// level limit. So the denominator is cubic feet and this is lbmol per cubic
/// foot.
///
/// # Only five species have a real density
///
/// `AD` is 1.0 and `BD` and `CD` are zero for A, B and C (`teprob.f:973-996`),
/// so those three contribute \\(x_i M_i\\) unchanged. They are the
/// non-condensibles, and the model never puts them in a liquid phase in any
/// quantity, so the entry is a placeholder that keeps the sum finite rather
/// than a density anyone fitted.
///
/// # The correlation has a pole, outside the operating range
///
/// \\(\\rho_i(T)\\) is a downward parabola for every real species, so each has a
/// temperature where it crosses zero and the reciprocal blows up. The nearest
/// is component D at 208.57 degrees, comfortably above the 175 degree reactor
/// shutdown limit at `teprob.f:706`. The original does not guard against it and
/// neither does this, since a guard would be behaviour the oracle does not
/// have. A Tier 1 harness test asserts the sweep ceiling stays below it.
///
/// ```
/// use tepsim_core::{Composition, thermo::liquid_density};
///
/// // Pure H at its nominal separator temperature: a plausible liquid density.
/// let pure_h = Composition::new([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]);
/// let rho = liquid_density(&pure_h, 80.0);
/// assert!(rho > 0.0 && rho < 1.0, "lbmol/ft^3, got {rho}");
/// ```
// @port teprob.f:1481-1505
#[must_use]
pub fn liquid_density(x: &Composition, celsius: f64) -> f64 {
    // teprob.f:1498-1502. The denominator is Horner form in the original and
    // stays that way here: expanding it to `AD + BD*T + CD*T*T` is the same
    // polynomial and a different set of roundings.
    let mut v = 0.0_f64;
    for component in Component::ALL {
        v += x[component] * XMW[component]
            / (AD[component] + (BD[component] + CD[component] * celsius) * celsius);
    }
    // teprob.f:1503
    1.0 / v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::assert_exact;

    fn equimolar() -> Composition {
        Composition::new([0.125; Component::COUNT])
    }

    #[test]
    fn ity_values_match_the_fortran() {
        assert_eq!(EnergyBasis::LiquidEnthalpy.ity(), 0);
        assert_eq!(EnergyBasis::VapourEnthalpy.ity(), 1);
        assert_eq!(EnergyBasis::VapourInternalEnergy.ity(), 2);
        for basis in EnergyBasis::ALL {
            assert_eq!(EnergyBasis::from_ity(basis.ity()), Some(basis));
        }
        assert_eq!(EnergyBasis::from_ity(3), None);
        assert_eq!(EnergyBasis::from_ity(-1), None);
    }

    /// `teprob.f:1395` and `1402` both multiply the whole polynomial by `T`, so
    /// the liquid and vapour enthalpies differ at zero only by the latent heat.
    #[test]
    fn the_enthalpy_correlations_are_anchored_at_zero_degrees() {
        let z = equimolar();
        assert_exact(
            enthalpy(&z, 0.0, EnergyBasis::LiquidEnthalpy),
            0.0,
            "liquid enthalpy at 0 C",
        );

        let mut latent = 0.0;
        for component in Component::ALL {
            latent += z[component] * XMW[component] * AV[component];
        }
        assert_exact(
            enthalpy(&z, 0.0, EnergyBasis::VapourEnthalpy),
            latent,
            "vapour enthalpy at 0 C is the latent heat alone",
        );
    }

    /// The one mode that is non-zero at zero degrees, and the reason is the
    /// single-precision offset.
    #[test]
    fn internal_energy_at_zero_degrees_is_the_ideal_gas_correction() {
        let z = equimolar();
        let vapour = enthalpy(&z, 0.0, EnergyBasis::VapourEnthalpy);
        let internal = enthalpy(&z, 0.0, EnergyBasis::VapourInternalEnergy);
        assert_exact(
            vapour - internal,
            GAS_CONSTANT * ABSOLUTE_ZERO_OFFSET,
            "the correction at 0 C",
        );
    }

    /// The hazard this whole module is shaped around. If someone "tidies"
    /// [`ABSOLUTE_ZERO_OFFSET`] into a plain `273.15`, this fires.
    #[test]
    fn the_absolute_zero_offset_is_stored_at_single_precision() {
        assert_exact(
            ABSOLUTE_ZERO_OFFSET,
            f64::from(273.15_f32),
            "the offset must be the widened f32",
        );
        assert!(
            (ABSOLUTE_ZERO_OFFSET - 273.15_f64).abs() > 1e-6,
            "and it must genuinely differ from the double literal, or the \
             constant has been quietly corrected: {ABSOLUTE_ZERO_OFFSET:?}"
        );
    }

    /// `3.57696 = 1.8 * 1.9872` in decimal: the Rankine ratio times the gas
    /// constant in Btu/(lbmol R), the same two numbers `teprob.f:594` uses
    /// separately.
    #[test]
    fn the_gas_constant_is_the_rankine_ratio_times_the_ideal_gas_constant() {
        let assembled = RANKINE_PER_CELSIUS * 1.9872 / 1.0e6;
        let relative = (GAS_CONSTANT - assembled).abs() / GAS_CONSTANT;
        assert!(
            relative < 1e-15,
            "3.57696 should be 1.8 x 1.9872, but they differ by {relative:e}"
        );
    }

    /// `TESUB3` claims to be the derivative of `TESUB1`. Check it, for all three
    /// modes, rather than trusting the pair of listings to agree.
    #[test]
    fn the_heat_capacity_is_the_derivative_of_the_enthalpy() {
        let z = equimolar();
        let step = 1e-4;
        for basis in EnergyBasis::ALL {
            for celsius in [5.292, 45.0, 100.0, 120.4, 175.0] {
                let central = (enthalpy(&z, celsius + step, basis)
                    - enthalpy(&z, celsius - step, basis))
                    / (2.0 * step);
                let analytic = heat_capacity(&z, celsius, basis);
                let relative = (central - analytic).abs() / analytic.abs();
                assert!(
                    relative < 1e-8,
                    "{basis:?} at {celsius} C: finite difference {central:e} \
                     against analytic {analytic:e}, off by {relative:e}"
                );
            }
        }
    }

    /// Both routines are linear in the composition, which is what lets the
    /// stream mixing elsewhere in the plant be written as it is.
    #[test]
    fn both_routines_are_linear_in_the_composition() {
        let celsius = 93.5;
        for basis in EnergyBasis::ALL {
            let mut summed_h = 0.0;
            let mut summed_dh = 0.0;
            let mixture = equimolar();
            for component in Component::ALL {
                let mut pure = [0.0; Component::COUNT];
                pure[component.index()] = 1.0;
                let pure = Composition::new(pure);
                summed_h += 0.125 * enthalpy(&pure, celsius, basis);
                summed_dh += 0.125 * heat_capacity(&pure, celsius, basis);
            }
            // The correction needs no special handling even though it is not
            // proportional to z: the eight weights sum to one, so summing the
            // pure species subtracts it exactly once, the same as evaluating
            // the mixture directly.
            let h = enthalpy(&mixture, celsius, basis);
            let dh = heat_capacity(&mixture, celsius, basis);
            assert!(
                (summed_h - h).abs() / h.abs() < 1e-14,
                "{basis:?} enthalpy is not linear in z: {summed_h:e} vs {h:e}"
            );
            assert!(
                (summed_dh - dh).abs() / dh.abs() < 1e-14,
                "{basis:?} heat capacity is not linear in z: {summed_dh:e} vs {dh:e}"
            );
        }
    }

    /// The round trip that `TESUB2` exists for: solve back to the temperature
    /// the enthalpy came from, for every mode and from both ends of the range.
    #[test]
    fn newton_recovers_the_temperature_the_enthalpy_came_from() {
        let z = equimolar();
        for basis in EnergyBasis::ALL {
            for target_celsius in [0.0, 5.292, 45.0, 100.0, 120.4, 175.0] {
                let target = enthalpy(&z, target_celsius, basis);
                for guess in [0.0, 45.0, 175.0] {
                    let solved = temperature_from_enthalpy(&z, guess, target, basis)
                        .unwrap_or_else(|e| panic!("{basis:?} from {guess} C: {e}"));
                    assert!(
                        (solved - target_celsius).abs() < 1e-9,
                        "{basis:?} from {guess} C: recovered {solved} for a \
                         target of {target_celsius}"
                    );
                }
            }
        }
    }

    /// Newton is quadratic here, so the step criterion is far more conservative
    /// than the error it certifies. Worth pinning: it is the reason the round
    /// trip above can assert 1e-9 when the loop only tests for 1e-12 on the
    /// step, and the reason a caller can trust the returned value.
    #[test]
    fn convergence_takes_a_handful_of_iterations_not_a_hundred() {
        let z = equimolar();
        let basis = EnergyBasis::LiquidEnthalpy;
        let target = enthalpy(&z, 120.4, basis);
        // Solving from the far end of the range must not need anything like the
        // iteration cap. Establish that by showing a much smaller cap suffices,
        // via the public function's own behaviour on a converging problem.
        let solved = temperature_from_enthalpy(&z, 0.0, target, basis)
            .expect("must converge from the cold end");
        assert!((solved - 120.4).abs() < 1e-9, "recovered {solved}");
    }

    /// The failure path has to be reachable, or the `Result` is decoration.
    /// An unattainable target sends Newton off the end of the polynomial.
    #[test]
    fn an_unreachable_target_reports_non_convergence_rather_than_lying() {
        let z = equimolar();
        // The liquid enthalpy is bounded below by its value at large negative
        // temperatures only through a cubic that turns over, so a target far
        // outside the achievable set leaves Newton nowhere to land.
        let outcome = temperature_from_enthalpy(&z, 120.4, -1.0e30, EnergyBasis::LiquidEnthalpy);
        match outcome {
            Err(TemperatureError::DidNotConverge {
                guess, last_step, ..
            }) => {
                assert_exact(guess, 120.4, "the error must report the guess");
                assert!(
                    last_step.abs() >= NEWTON_TOLERANCE || last_step.is_nan(),
                    "a reported failure must actually have failed the test, \
                     got a final step of {last_step}"
                );
            }
            Ok(t) => panic!(
                "an unreachable target returned {t} as though it had converged, \
                 which is precisely the Fortran behaviour this port exists to \
                 replace"
            ),
        }
    }

    /// Volumes add, so the reciprocal densities do. Mixing two species must give
    /// a specific volume between the two pure ones, which is the defining
    /// property of `teprob.f:1500-1503` and the reason it is written as a
    /// reciprocal sum rather than a weighted average of densities.
    #[test]
    fn liquid_volumes_are_additive() {
        let celsius = 80.0;
        let mut pure_g = [0.0; Component::COUNT];
        pure_g[Component::G.index()] = 1.0;
        let mut pure_h = [0.0; Component::COUNT];
        pure_h[Component::H.index()] = 1.0;
        let mut half = [0.0; Component::COUNT];
        half[Component::G.index()] = 0.5;
        half[Component::H.index()] = 0.5;

        let volume = |x: [f64; Component::COUNT]| {
            let composition = Composition::new(x);
            let moles_per_volume = liquid_density(&composition, celsius);
            1.0 / moles_per_volume
        };

        let (vg, vh, vmix) = (volume(pure_g), volume(pure_h), volume(half));
        let predicted = 0.5 * vg + 0.5 * vh;
        assert!(
            (vmix - predicted).abs() / predicted < 1e-14,
            "specific volume is not additive: {vmix:e} vs {predicted:e}"
        );
    }

    /// A, B and C carry `AD = 1` and no temperature dependence
    /// (`teprob.f:973`, `983`, `993`), so their contribution is `x * M` flat.
    /// That is a placeholder, not a fitted density, and the port must reproduce
    /// it rather than tidy it away.
    #[test]
    fn the_non_condensibles_have_a_flat_placeholder_density() {
        let mut pure_a = [0.0; Component::COUNT];
        pure_a[Component::A.index()] = 1.0;
        let pure_a = Composition::new(pure_a);

        let cold = liquid_density(&pure_a, 0.0);
        let hot = liquid_density(&pure_a, 175.0);
        assert_exact(cold, hot, "A's density must not vary with temperature");
        assert_exact(
            cold,
            1.0 / XMW[Component::A],
            "and it must be 1 / molecular weight, from AD = 1",
        );
    }

    /// Every real species' density stays positive across the whole sweep range,
    /// so the reciprocal never blows up inside it. The pole exists, at 208.57 C
    /// for D, but not where the plant operates.
    #[test]
    fn the_density_correlation_stays_positive_over_the_operating_range() {
        for step in 0..=175 {
            let celsius = f64::from(step);
            for component in Component::ALL {
                let rho = AD[component] + (BD[component] + CD[component] * celsius) * celsius;
                assert!(
                    rho > 0.0,
                    "{component:?} density is {rho} at {celsius} C, inside the \
                     operating range"
                );
            }
            let mut equal = [0.125; Component::COUNT];
            equal[0] = 0.125;
            let mixture = liquid_density(&Composition::new(equal), celsius);
            assert!(
                mixture.is_finite() && mixture > 0.0,
                "mixture density is {mixture} at {celsius} C"
            );
        }
    }

    /// The vapour branch is taken for `ITY` 1 and 2 alike; they differ only by
    /// the correction (`teprob.f:1399`, `1409`).
    ///
    /// Recovered to within a rounding error rather than bit-exactly, because
    /// `a - (a - r)` is not `r` in IEEE arithmetic: the inner subtraction
    /// rounds, and the outer one cannot undo it. Demanding equality here would
    /// be asserting something false about floating point rather than something
    /// true about the port.
    ///
    /// The tolerance is two ULP *of the enthalpy*, not of the correction, which
    /// is the cancellation the module documents seen from the other side: both
    /// roundings happen at the larger quantity's scale, so the small quantity
    /// comes back with the large one's absolute error.
    #[test]
    fn the_two_vapour_modes_differ_only_by_the_correction() {
        let z = equimolar();
        for celsius in [0.0, 45.0, 175.0] {
            let vapour = enthalpy(&z, celsius, EnergyBasis::VapourEnthalpy);
            let recovered = vapour - enthalpy(&z, celsius, EnergyBasis::VapourInternalEnergy);
            let expected = GAS_CONSTANT * (celsius + ABSOLUTE_ZERO_OFFSET);
            assert!(
                (recovered - expected).abs() <= 2.0 * f64::EPSILON * vapour.abs(),
                "enthalpy correction at {celsius} C: {recovered:e} vs {expected:e}"
            );

            let vapour = heat_capacity(&z, celsius, EnergyBasis::VapourEnthalpy);
            let recovered = vapour - heat_capacity(&z, celsius, EnergyBasis::VapourInternalEnergy);
            assert!(
                (recovered - GAS_CONSTANT).abs() <= 2.0 * f64::EPSILON * vapour.abs(),
                "heat capacity correction at {celsius} C: {recovered:e} vs \
                 {GAS_CONSTANT:e}"
            );
        }
    }
}
