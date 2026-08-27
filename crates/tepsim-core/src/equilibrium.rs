//! Vapour-liquid equilibrium and the three vessel pressures.
//!
//! Ported from `teprob.f:473-502`. Given the unpacked state, this produces the
//! partial pressures in the reactor and the separator, the totals, the vapour
//! mole fractions, and the vapour holdups that the balances and the reaction
//! rates consume.
//!
//! # The equations
//!
//! The vapour space of a two-phase vessel is what the liquid does not occupy:
//!
//! \\[ V_v = V_\\text{total} - V_l \\]
//!
//! A, B and C are non-condensible and are treated as ideal gases, so their
//! partial pressures come from the holdup directly:
//!
//! \\[ p_i = \\frac{n_i R T}{V_v}, \\qquad i \\in \\{A, B, C\\} \\]
//!
//! D through H are condensible, so their partial pressures come from Raoult's
//! law with an Antoine vapour pressure, in degrees Celsius:
//!
//! \\[
//!   p_i = x_i \\exp\\!\\left(A_i + \\frac{B_i}{T + C_i}\\right),
//!   \\qquad i \\in \\{D \\ldots H\\}
//! \\]
//!
//! The total is their sum, the vapour composition is
//! \\(y_i = p_i / P\\), and the vapour holdup follows from the ideal gas law
//! applied to the mixture:
//!
//! \\[ N_v = \\frac{P V_v}{R T}, \\qquad n_i = N_v y_i \\]
//!
//! The mixing zone is all vapour at a fixed volume, so it needs only the
//! middle step: \\(P = N R T / V\\).
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `VVR`, `VVS` | [`VapourSpace::volume`] | vapour volume |
//! | `PPR`, `PPS` | [`VapourSpace::partial`] | partial pressures |
//! | `PTR`, `PTS` | [`VapourSpace::pressure`] | total pressure |
//! | `PTV` | [`Equilibrium::mixing_pressure`] | mixing zone pressure |
//! | `XVR`, `XVS` | [`VapourSpace::fractions`] | vapour mole fractions |
//! | `UTVR`, `UTVS` | [`VapourSpace::total`] | total vapour moles |
//! | `UCVR`, `UCVS` | [`VapourSpace::moles`] | vapour component holdups |
//! | `RG` | [`GAS_CONSTANT`] | gas constant in mmHg ft^3 / lbmol K |
//! | `VTR`, `VTS`, `VTV` | [`crate::constants::VTR`] etc. | vessel volumes |
//!
//! # `UCVR` is both an input and an output
//!
//! `teprob.f:418` fills `UCVR(1..3)` from the state and `teprob.f:500` fills
//! `UCVR(4..8)` from the equilibrium computed here. It is one Fortran array
//! written by two different mechanisms, and the halves are not
//! interchangeable: the non-condensibles are integrated, the condensibles are
//! derived. [`crate::vessels::Unpacked`] leaves the D-H slots zero for exactly
//! this reason, and [`VapourSpace::moles`] is where they are filled in.
//!
//! Read in that order it also explains why the A/B/C partial pressures come
//! first at `teprob.f:478-483`: they need `UCVR` before it is overwritten.
//!
//! # This is where bit equality with gfortran ends
//!
//! `DEXP` at `teprob.f:485` and `teprob.f:488` is the model's first
//! transcendental call, and the port answers it from the vendored pure-Rust
//! `libm` rather than the platform's, for the determinism reason set out in
//! [`crate::math`]. Measured over the whole Antoine range this model reaches,
//! the two disagree on 9.945% of arguments, by exactly one ULP.
//!
//! Everything downstream of a D-H partial pressure therefore carries about
//! 1.1e-16 of relative difference from the Fortran that no amount of care in
//! the algebra removes. Tier 2's gate is 1e-12. The bit-exactness claim does
//! not disappear, it moves: under the `libm-system` feature the transcendental
//! agrees with gfortran exactly, and then the port is bit-identical again, so
//! the algebra is still held to zero ULP rather than to a tolerance.
//!
//! # Precision hazards in this range
//!
//! [`GAS_CONSTANT`] is written `RG=998.9` at `teprob.f:475` with no `D`
//! suffix, so it is single precision and multiplies every ideal-gas partial
//! pressure and divides every vapour holdup. All 24 Antoine coefficients are
//! single as well ([`crate::constants::AVP`]), and they sit inside the
//! exponent argument, where an absolute error of 5e-7 becomes a relative error
//! of the same size in the result.

// Every float expression here reproduces `teprob.f`'s association and rounding
// exactly. `mul_add` would fuse two roundings into one and give a more
// accurate answer that no longer matches the oracle.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]

use crate::component::{ByComponent, Component, Composition};
use crate::constants::{AVP, BVP, CVP, VTR, VTS, VTV, single};
use crate::math::exp;
use crate::vessels::Unpacked;

/// The gas constant in the units this block works in, from `teprob.f:475`.
///
/// Pressures are mmHg, volumes cubic feet, holdups lbmol and temperatures
/// kelvin, which puts it at 998.9 rather than at any familiar value. Nothing
/// in the original says so; it is recoverable from the units of the equation.
///
/// # This is single precision
///
/// The literal carries no `D` suffix, so gfortran stores 998.90002441406250,
/// not 998.9. The difference is 2.4e-8 relative, which is 4 orders of
/// magnitude past the Tier 2 gate, and it multiplies six of the eight partial
/// pressures in every vessel. See [`crate::constants`].
///
/// Unrelated to [`crate::thermo::GAS_CONSTANT`], which is the same physical
/// quantity in the energy units the enthalpy correlations use.
pub const GAS_CONSTANT: f64 = single(998.9);

/// The non-condensibles, which are ideal gases here (`teprob.f:478`).
const NON_CONDENSIBLE: [Component; 3] = [Component::A, Component::B, Component::C];

/// The condensibles, which get an Antoine vapour pressure (`teprob.f:484`).
const CONDENSIBLE: [Component; 5] = [
    Component::D,
    Component::E,
    Component::F,
    Component::G,
    Component::H,
];

/// The vapour phase of a two-phase vessel, in equilibrium with its liquid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct VapourSpace {
    /// `VVR`/`VVS`: vessel volume less liquid volume, cubic feet.
    pub volume: f64,
    /// `PPR`/`PPS`: partial pressure of each species, mmHg.
    pub partial: ByComponent<f64>,
    /// `PTR`/`PTS`: total pressure, mmHg, summed in Fortran order.
    pub pressure: f64,
    /// `XVR`/`XVS`: vapour mole fractions.
    pub fractions: Composition,
    /// `UTVR`/`UTVS`: total vapour moles, from the ideal gas law.
    pub total: f64,
    /// `UCVR`/`UCVS`: vapour component holdups.
    ///
    /// A, B and C are carried through from the state; D through H are computed
    /// here. See the module documentation.
    pub moles: ByComponent<f64>,
}

/// Everything `teprob.f:473-502` produces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Equilibrium {
    /// The reactor's vapour space.
    pub reactor: VapourSpace,
    /// The separator's vapour space.
    pub separator: VapourSpace,
    /// `PTV`: the mixing zone pressure, mmHg.
    ///
    /// The mixing zone is single-phase and fixed-volume, so it has no vapour
    /// space of its own to describe: the pressure is the whole answer.
    pub mixing_pressure: f64,
}

/// Solve the vapour-liquid equilibrium in the reactor and the separator, and
/// compute the three pressures.
///
/// Total pressures accumulate in the original's order, A through H, because
/// reassociating a sum of eight terms spanning four orders of magnitude
/// changes the last bits.
// @port teprob.f:473-502
#[must_use]
pub fn equilibrium(unpacked: &Unpacked) -> Equilibrium {
    // teprob.f:473-474. Nothing stops this going negative on an adversarial
    // state, and the original does not check either; a negative vapour volume
    // gives negative pressures, which is a true statement about a state that
    // cannot physically occur.
    let vvr = VTR - unpacked.reactor.volume;
    let vvs = VTS - unpacked.separator.volume;

    let tcr = unpacked.reactor.celsius;
    let tcs = unpacked.separator.celsius;
    let tkr = unpacked.reactor.kelvin();
    let tks = unpacked.separator.kelvin();

    let mut ppr = ByComponent::new([0.0; Component::COUNT]);
    let mut pps = ByComponent::new([0.0; Component::COUNT]);

    // teprob.f:476-477
    let mut ptr = 0.0;
    let mut pts = 0.0;

    // teprob.f:478-483. Both vessels share the loop in the original; they
    // accumulate into separate totals, so the interleaving is presentation
    // rather than arithmetic.
    for component in NON_CONDENSIBLE {
        ppr[component] = unpacked.reactor_vapour[component] * GAS_CONSTANT * tkr / vvr;
        ptr += ppr[component];
        pps[component] = unpacked.separator_vapour[component] * GAS_CONSTANT * tks / vvs;
        pts += pps[component];
    }

    // teprob.f:484-491. Antoine in degrees Celsius, not kelvin: the offset is
    // folded into `CVP`.
    for component in CONDENSIBLE {
        let vpr = exp(AVP[component] + BVP[component] / (tcr + CVP[component]));
        ppr[component] = vpr * unpacked.reactor.fractions[component];
        ptr += ppr[component];
        let vps = exp(AVP[component] + BVP[component] / (tcs + CVP[component]));
        pps[component] = vps * unpacked.separator.fractions[component];
        pts += pps[component];
    }

    // teprob.f:492
    let mixing_pressure = unpacked.mixing.total * GAS_CONSTANT * unpacked.mixing.kelvin() / VTV;

    // teprob.f:493-496
    let mut xvr = [0.0; Component::COUNT];
    let mut xvs = [0.0; Component::COUNT];
    for component in Component::ALL {
        xvr[component.index()] = ppr[component] / ptr;
        xvs[component.index()] = pps[component] / pts;
    }
    // Unchecked, for the same reason as in `vessels`: an adversarial state can
    // drive a total pressure to zero, and the resulting NaN is the honest
    // answer rather than something to assert away.
    let xvr = Composition::new_unchecked(xvr);
    let xvs = Composition::new_unchecked(xvs);

    // teprob.f:497-498
    let utvr = ptr * vvr / GAS_CONSTANT / tkr;
    let utvs = pts * vvs / GAS_CONSTANT / tks;

    // teprob.f:499-502. Only D through H: A, B and C keep the values the state
    // gave them. See the module documentation.
    let mut ucvr = unpacked.reactor_vapour;
    let mut ucvs = unpacked.separator_vapour;
    for component in CONDENSIBLE {
        ucvr[component] = utvr * xvr[component];
        ucvs[component] = utvs * xvs[component];
    }

    Equilibrium {
        reactor: VapourSpace {
            volume: vvr,
            partial: ppr,
            pressure: ptr,
            fractions: xvr,
            total: utvr,
            moles: ucvr,
        },
        separator: VapourSpace {
            volume: vvs,
            partial: pps,
            pressure: pts,
            fractions: xvs,
            total: utvs,
            moles: ucvs,
        },
        mixing_pressure,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::State;
    use crate::testing::assert_exact;
    use crate::vessels::{TemperatureSeeds, unpack};

    /// A state near the nominal operating point, enough to exercise the whole
    /// block without an oracle.
    fn plausible() -> Unpacked {
        let mut y = State::default();
        for component in Component::ALL {
            y.reactor.moles[component] = 10.0 + component.index() as f64;
            y.separator.moles[component] = 20.0 + component.index() as f64;
            y.stripper.moles[component] = 5.0;
            y.mixing.moles[component] = 40.0;
        }
        y.reactor.energy = 200.0;
        y.separator.energy = 180.0;
        y.stripper.energy = 30.0;
        y.mixing.energy = 320.0;
        unpack(&y, TemperatureSeeds::default()).expect("converges")
    }

    /// The single-precision gas constant, spelled out. Getting this wrong is
    /// worth 2.4e-8 relative on most of the plant.
    #[test]
    fn the_gas_constant_is_the_widened_single_precision_literal() {
        assert_exact(GAS_CONSTANT, 998.9_f32 as f64, "RG is single precision");
        // 998.9 rounds *up* in binary32, so a double literal would be smaller.
        // Constant-folded, hence the const block.
        const { assert!(GAS_CONSTANT > 998.9_f64) };
    }

    /// The total is the sum of the parts, in the original's order. A quiet
    /// reassociation would still pass this; the oracle differential is what
    /// catches that. This catches a dropped term.
    #[test]
    fn each_total_pressure_is_the_sum_of_its_partials() {
        let e = equilibrium(&plausible());
        for space in [e.reactor, e.separator] {
            let mut sum = 0.0;
            for component in NON_CONDENSIBLE {
                sum += space.partial[component];
            }
            for component in CONDENSIBLE {
                sum += space.partial[component];
            }
            assert_exact(space.pressure, sum, "the total is the sum of the partials");
            assert!(space.fractions.sums_to_one());
        }
    }

    /// The non-condensibles keep the holdups the state gave them and the
    /// condensibles are replaced. Getting this backwards would silently
    /// discard three integrated states.
    #[test]
    fn only_the_condensible_holdups_are_replaced() {
        let unpacked = plausible();
        let e = equilibrium(&unpacked);
        for component in NON_CONDENSIBLE {
            assert_exact(
                e.reactor.moles[component],
                unpacked.reactor_vapour[component],
                "A, B and C come from the state",
            );
        }
        for component in CONDENSIBLE {
            assert_exact(
                unpacked.reactor_vapour[component],
                0.0,
                "unpack leaves D-H empty",
            );
            assert!(
                e.reactor.moles[component] > 0.0,
                "D-H must be filled in from the equilibrium"
            );
        }
    }

    /// The mixing zone is fixed-volume and single-phase, so its pressure is
    /// the ideal gas law with nothing subtracted.
    #[test]
    fn the_mixing_zone_pressure_is_the_ideal_gas_law_at_fixed_volume() {
        let unpacked = plausible();
        let e = equilibrium(&unpacked);
        assert_exact(
            e.mixing_pressure,
            unpacked.mixing.total * GAS_CONSTANT * unpacked.mixing.kelvin() / VTV,
            "PTV",
        );
    }

    /// Antoine is evaluated in Celsius. Passing kelvin would be a plausible
    /// mistake that no unit check catches, and it lands about 5% low on the
    /// reactor pressure rather than obviously wrong.
    #[test]
    fn the_antoine_argument_is_in_celsius() {
        let unpacked = plausible();
        let e = equilibrium(&unpacked);
        let d = Component::D;
        let expected = exp(AVP[d] + BVP[d] / (unpacked.reactor.celsius + CVP[d]))
            * unpacked.reactor.fractions[d];
        assert_exact(e.reactor.partial[d], expected, "PPR(4)");
        let kelvin = exp(AVP[d] + BVP[d] / (unpacked.reactor.kelvin() + CVP[d]))
            * unpacked.reactor.fractions[d];
        assert!(
            kelvin.to_bits() != e.reactor.partial[d].to_bits(),
            "Celsius and kelvin gave the same answer, so this test proves nothing"
        );
    }

    /// The vapour volume is the vessel less the liquid, and the two vessels do
    /// not share a volume. Crossing them is a copy-paste error the oracle
    /// would catch, but only after the whole block was written.
    #[test]
    fn the_vapour_volume_is_the_vessel_less_the_liquid() {
        let unpacked = plausible();
        let e = equilibrium(&unpacked);
        assert_exact(e.reactor.volume, VTR - unpacked.reactor.volume, "VVR");
        assert_exact(e.separator.volume, VTS - unpacked.separator.volume, "VVS");
    }
}
