//! The four gas-phase reactions, their rates, and the reaction enthalpy.
//!
//! Ported from `teprob.f:503-528`. Everything here happens in the reactor's
//! vapour space, on the partial pressures [`mod@crate::equilibrium`] produced.
//!
//! # The reactions
//!
//! ```text
//! 1:  A + C + D -> G          2:  A + C + E -> H
//! 3:  A + E     -> F          4:  3 D       -> 2 F
//! ```
//!
//! B appears in none of them. It enters with the feed and leaves only through
//! the purge, which is the entire reason the plant has one.
//!
//! # The equations
//!
//! Rates 1 and 2 are Arrhenius in reactor temperature with fractional pressure
//! orders on A and C, multiplied by a disturbance drift factor:
//!
//! \\[
//!   r_1 = f_1 \\, e^{\\,a_1 - E_1/T} \; p_A^{1.1544} \\, p_C^{0.3735} \\, p_D \\, V_v
//! \\]
//!
//! \\[
//!   r_2 = f_2 \\, e^{\\,a_2 - E_2/T} \; p_A^{1.1544} \\, p_C^{0.3735} \\, p_E \\, V_v
//! \\]
//!
//! Rates 3 and 4 are first order in each reactant, and rate 4 shares rate 3's
//! exponential rather than having one of its own:
//!
//! \\[
//!   r_3 = e^{\\,a_3 - E_3/T} \\, p_A \\, p_E \\, V_v,
//!   \\qquad
//!   r_4 = 0.767488334 \; e^{\\,a_3 - E_3/T} \\, p_A \\, p_D \\, V_v
//! \\]
//!
//! Net production per species follows the stoichiometry, and the heat release
//! comes from reactions 1 and 2 only:
//!
//! \\[ Q = r_1 h_1 + r_2 h_2 \\]
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `RR(1..4)` | [`Kinetics::rates`] | extent of each reaction |
//! | `CRXR(1..8)` | [`Kinetics::production`] | net production per species |
//! | `RH` | [`Kinetics::heat`] | heat of reaction |
//! | `R1F`, `R2F` (at `403`) | [`ReactionDrift`] | disturbance drift factors |
//! | `HTR(1)`, `HTR(2)` | [`HEAT_OF_REACTION`] | heats of reaction 1 and 2 |
//!
//! # `R1F` and `R2F` are two different quantities under one name
//!
//! They arrive from `TESUB8(7)` and `TESUB8(8)` at `teprob.f:415-416` as the
//! IDV(13) kinetics-drift multipliers, are consumed at `503-504`, and are then
//! **reassigned in place** at `508-509` to hold the fractional pressure
//! powers. The two meanings share nothing but the storage.
//!
//! Reading `510` as though `R1F` were still the drift factor gives a plausible
//! and completely wrong rate law. This port gives the two roles separate names
//! ([`ReactionDrift`] and a local `pressure_order`), which is delta D-002: no
//! numerical effect, and the only defence against a misreading that no test
//! would catch, because a wrong-but-consistent reading still reproduces
//! itself.
//!
//! The visible consequence is that after `TEFUNC` returns, `COMMON`'s `R1F`
//! holds a pressure power and not the drift factor a caller might expect. Tier
//! 2 has to fetch the drift from `TESUB8` directly rather than from the block.
//!
//! # `CRXR(2)` is never assigned
//!
//! Seven of the eight slots are written at `teprob.f:521-527`. `CRXR(2)`, the
//! inert, is not, and it is read anyway at `teprob.f:763`. It works because
//! `COMMON` is zero-initialised and nothing ever writes it, so B's net
//! production is zero by static initialisation rather than by statement. That
//! is delta D-003, class A: the value is right, the mechanism is an accident.
//! Here the slot is explicitly zero and a test asserts the oracle agrees.
//!
//! # Precision hazards in this range
//!
//! Two, and both are silent.
//!
//! Every literal except `0.767488334D0` and `1.5D0` is **single precision**:
//! the three pre-exponentials, the three activation energies, the gas constant
//! `1.987`, and both fractional exponents. See [`crate::constants::single`].
//!
//! Worse, `40000.0/1.987` is a quotient of *two* single-precision literals, so
//! Fortran evaluates the division itself in single precision. Widening the
//! operands first and dividing in double is wrong by 4e-9 relative, inside a
//! `DEXP` argument. That is what [`crate::constants::single_quotient`] is for,
//! and its documentation carries the measured bit patterns.

// Every float expression here reproduces `teprob.f`'s association and rounding
// exactly; see `crate::thermo` for why that forbids the "better" forms.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]
// `RR(I)=RR(I)*VVR` at teprob.f:519, written the way the listing writes it so
// the module can be checked against the source a line at a time. IEEE
// multiplication makes the two forms bit-identical; this is about readability
// against the original, as in `crate::thermo`.
#![allow(
    clippy::assign_op_pattern,
    reason = "operand order is transcribed from the Fortran, not chosen"
)]

use crate::component::{ByComponent, Component};
use crate::constants::{single, single_quotient};
use crate::equilibrium::VapourSpace;
use crate::math::{exp, pow};

/// The heats of reaction for reactions 1 and 2, `HTR` (`teprob.f:1122-1123`).
///
/// Reactions 3 and 4 contribute no heat: the original simply does not include
/// them in `RH`, and `HTR(3)` is declared but never assigned or read.
///
/// Both literals carry a `D` suffix, so unlike almost everything else in this
/// module they are full double precision.
//
// @port teprob.f:1122-1123
pub const HEAT_OF_REACTION: [f64; 2] = [0.06899381054, 0.05];

// Transcribed digit for digit from the listing and rounded by `single`, not
// pre-rounded by hand. `clippy::excessive_precision` is exactly backwards
// here: the digits the f32 cannot hold are the evidence that the constant was
// copied rather than retyped, and dropping them is how a transcription error
// gets in.
#[allow(
    clippy::excessive_precision,
    reason = "transcribed verbatim from teprob.f; `single` does the rounding"
)]
/// Reaction 1: pre-exponential, `teprob.f:503`. Single precision.
const LN_PRE_EXPONENTIAL_1: f64 = single(31.5859536);
#[allow(
    clippy::excessive_precision,
    reason = "transcribed verbatim from teprob.f; `single` does the rounding"
)]
/// Reaction 2: pre-exponential, `teprob.f:504`. Single precision.
const LN_PRE_EXPONENTIAL_2: f64 = single(3.00094014);
#[allow(
    clippy::excessive_precision,
    reason = "transcribed verbatim from teprob.f; `single` does the rounding"
)]
/// Reaction 3: pre-exponential, `teprob.f:505`. Single precision.
const LN_PRE_EXPONENTIAL_3: f64 = single(53.4060443);

/// Reaction 1: `40000.0/1.987`, folded in single precision (`teprob.f:503`).
const ACTIVATION_1: f64 = single_quotient(40000.0, 1.987);
/// Reaction 2: `20000.0/1.987`, folded in single precision (`teprob.f:504`).
const ACTIVATION_2: f64 = single_quotient(20000.0, 1.987);
/// Reaction 3: `60000.0/1.987`, folded in single precision (`teprob.f:505`).
const ACTIVATION_3: f64 = single_quotient(60000.0, 1.987);

/// Reaction 4's rate relative to reaction 3's, `teprob.f:506`.
///
/// One of only two double-precision literals in this range.
const REACTION_4_RATIO: f64 = 0.767488334;

/// Pressure order on A in reactions 1 and 2, `teprob.f:508`. Single precision.
const ORDER_A: f64 = single(1.1544);
/// Pressure order on C in reactions 1 and 2, `teprob.f:509`. Single precision.
const ORDER_C: f64 = single(0.3735);

/// Moles of D consumed per unit of reaction 4, `teprob.f:523`.
///
/// The other double-precision literal here. Reaction 4 is `3D -> 2F`, so the
/// stoichiometric coefficient is 1.5 relative to the rate as defined.
const D_PER_REACTION_4: f64 = 1.5;

/// The IDV(13) kinetics drift multipliers, `R1F` and `R2F` as they arrive.
///
/// Produced by `TESUB8(7, t)` and `TESUB8(8, t)` at `teprob.f:415-416`, which
/// is disturbance-walk machinery and therefore Phase 3. Until then they are an
/// input, so that this module can be ported and validated without it.
///
/// See the module documentation: the same two Fortran variables are reused for
/// something else eight lines later.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactionDrift {
    /// `R1F` at `teprob.f:503`: multiplies reaction 1.
    pub first: f64,
    /// `R2F` at `teprob.f:504`: multiplies reaction 2.
    pub second: f64,
}

impl Default for ReactionDrift {
    /// No drift. `TESUB8` returns `SZERO`, which is 1 for both channels, when
    /// the walk has not moved.
    fn default() -> Self {
        Self {
            first: 1.0,
            second: 1.0,
        }
    }
}

/// Everything `teprob.f:503-528` produces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Kinetics {
    /// `RR(1..4)`: the extent of each reaction, already multiplied by the
    /// vapour volume (`teprob.f:518-520`).
    pub rates: [f64; 4],
    /// `CRXR(1..8)`: net production per species, negative for consumption.
    ///
    /// B is always exactly zero; see the module documentation.
    pub production: ByComponent<f64>,
    /// `RH`: heat released by reactions 1 and 2.
    pub heat: f64,
}

/// Evaluate the four reaction rates, the net production, and the heat release.
///
/// `reactor` is the reactor's vapour space and `reactor_kelvin` is `TKR`.
// @port  teprob.f:503-528
// @delta D-002 class=A teprob.f:503-511
// @delta D-003 class=A teprob.f:521-527
#[must_use]
pub fn kinetics(reactor: &VapourSpace, reactor_kelvin: f64, drift: ReactionDrift) -> Kinetics {
    let partial = &reactor.partial;

    // teprob.f:503-506. Reaction 4 shares reaction 3's exponential.
    let mut rates = [
        exp(LN_PRE_EXPONENTIAL_1 - ACTIVATION_1 / reactor_kelvin) * drift.first,
        exp(LN_PRE_EXPONENTIAL_2 - ACTIVATION_2 / reactor_kelvin) * drift.second,
        exp(LN_PRE_EXPONENTIAL_3 - ACTIVATION_3 / reactor_kelvin),
        0.0,
    ];
    rates[3] = rates[2] * REACTION_4_RATIO;

    // teprob.f:507-515. `pow` of a non-positive base is not defined the way
    // this rate law needs, so the original guards it and zeroes both rates
    // instead. The guard is on A and C only: `PPR(4)` and `PPR(5)` enter
    // linearly and are allowed to be negative.
    if partial[Component::A] > 0.0 && partial[Component::C] > 0.0 {
        // These two are `R1F` and `R2F` in the original, reusing the names the
        // drift factors arrived under. See the module documentation.
        let order_a = pow(partial[Component::A], ORDER_A);
        let order_c = pow(partial[Component::C], ORDER_C);
        rates[0] = rates[0] * order_a * order_c * partial[Component::D];
        rates[1] = rates[1] * order_a * order_c * partial[Component::E];
    } else {
        rates[0] = 0.0;
        rates[1] = 0.0;
    }

    // teprob.f:516-517. Outside the branch, so reactions 3 and 4 still run
    // when 1 and 2 have been zeroed.
    rates[2] = rates[2] * partial[Component::A] * partial[Component::E];
    rates[3] = rates[3] * partial[Component::A] * partial[Component::D];

    // teprob.f:518-520. Rates become extents by scaling with the vapour
    // volume, which is why everything downstream uses the scaled values.
    for rate in &mut rates {
        *rate = *rate * reactor.volume;
    }
    let [r1, r2, r3, r4] = rates;

    // teprob.f:521-527. Seven assignments for eight species: B is inert and is
    // never written. See the module documentation.
    let mut production = ByComponent::new([0.0; Component::COUNT]);
    production[Component::A] = -r1 - r2 - r3;
    production[Component::C] = -r1 - r2;
    production[Component::D] = -r1 - D_PER_REACTION_4 * r4;
    production[Component::E] = -r2 - r3;
    production[Component::F] = r3 + r4;
    production[Component::G] = r1;
    production[Component::H] = r2;

    Kinetics {
        rates,
        production,
        // teprob.f:528
        heat: r1 * HEAT_OF_REACTION[0] + r2 * HEAT_OF_REACTION[1],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::equilibrium::equilibrium;
    use crate::state::State;
    use crate::testing::assert_exact;
    use crate::vessels::{TemperatureSeeds, Unpacked, unpack};

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

    fn solved() -> Kinetics {
        let unpacked = plausible();
        let eq = equilibrium(&unpacked);
        kinetics(
            &eq.reactor,
            unpacked.reactor.kelvin(),
            ReactionDrift::default(),
        )
    }

    /// The activation energies are folded in single precision. This is the
    /// hazard of the whole module, worth 4e-9 relative if missed.
    #[test]
    fn the_activation_energies_are_folded_in_single_precision() {
        assert_exact(
            ACTIVATION_1,
            f64::from(40000.0_f32 / 1.987_f32),
            "40000.0/1.987",
        );
        assert!(
            ACTIVATION_1.to_bits() != (40000.0_f64 / 1.987_f64).to_bits(),
            "the single-precision fold agreed with the double one, so this \
             constant proves nothing"
        );
        // 29 trailing zero bits is the signature of a widened f32.
        for (name, value) in [
            ("ACTIVATION_1", ACTIVATION_1),
            ("ACTIVATION_2", ACTIVATION_2),
            ("ACTIVATION_3", ACTIVATION_3),
        ] {
            assert_eq!(
                value.to_bits() & 0x1FFF_FFFF,
                0,
                "{name} has mantissa bits an f32 cannot hold, so it was not \
                 folded in single precision"
            );
        }
    }

    /// Stoichiometry, read off the four reactions rather than off the Fortran.
    /// An independent statement of the same thing, so a transcription slip in
    /// `teprob.f:521-527` has somewhere to fail.
    #[test]
    fn net_production_follows_the_stoichiometry() {
        let k = solved();
        let [r1, r2, r3, r4] = k.rates;
        // A + C + D -> G, A + C + E -> H, A + E -> F, 3D -> 2F.
        assert_exact(k.production[Component::G], r1, "G comes only from 1");
        assert_exact(k.production[Component::H], r2, "H comes only from 2");
        assert_exact(k.production[Component::F], r3 + r4, "F from 3 and 4");
        assert_exact(
            k.production[Component::B],
            0.0,
            "B is inert and takes part in nothing",
        );
        // A is consumed by 1, 2 and 3, and not by 4.
        assert_exact(k.production[Component::A], -r1 - r2 - r3, "A");
        // C is consumed by 1 and 2 only.
        assert_exact(k.production[Component::C], -r1 - r2, "C");
        // E is consumed by 2 and 3.
        assert_exact(k.production[Component::E], -r2 - r3, "E");
    }

    /// Reaction 4 consumes three D per two F, so D's coefficient is 1.5 times
    /// the rate. Getting it wrong is a plausible slip that leaves every other
    /// species right.
    #[test]
    fn reaction_four_consumes_three_d_for_every_two_f() {
        let k = solved();
        let [r1, _, _, r4] = k.rates;
        assert_exact(k.production[Component::D], -r1 - 1.5 * r4, "D");
        assert!(
            k.production[Component::D] < -r1 - r4,
            "D is consumed faster than one per unit of reaction 4"
        );
    }

    /// Only reactions 1 and 2 release heat. Including 3 or 4 would be a
    /// reasonable-looking physical assumption and is not what the model says.
    #[test]
    fn only_the_first_two_reactions_release_heat() {
        let k = solved();
        let [r1, r2, _, _] = k.rates;
        assert_exact(
            k.heat,
            r1 * HEAT_OF_REACTION[0] + r2 * HEAT_OF_REACTION[1],
            "RH",
        );
    }

    /// The guard is on A and C, and it zeroes reactions 1 and 2 while leaving
    /// 3 and 4 running. Both halves of that matter.
    #[test]
    fn a_non_positive_pressure_on_a_or_c_zeroes_only_the_first_two_reactions() {
        let unpacked = plausible();
        let mut eq = equilibrium(&unpacked);
        eq.reactor.partial[Component::A] = 0.0;

        let k = kinetics(
            &eq.reactor,
            unpacked.reactor.kelvin(),
            ReactionDrift::default(),
        );
        assert_exact(k.rates[0], 0.0, "reaction 1 is off");
        assert_exact(k.rates[1], 0.0, "reaction 2 is off");
        // 3 and 4 are outside the branch, but both are linear in A, so with A
        // at zero they vanish too. Use C instead to separate the two effects.
        let mut eq = equilibrium(&unpacked);
        eq.reactor.partial[Component::C] = -1.0;
        let k = kinetics(
            &eq.reactor,
            unpacked.reactor.kelvin(),
            ReactionDrift::default(),
        );
        assert_exact(k.rates[0], 0.0, "reaction 1 is off");
        assert_exact(k.rates[1], 0.0, "reaction 2 is off");
        assert!(
            k.rates[2] != 0.0 && k.rates[3] != 0.0,
            "reactions 3 and 4 are outside the guard and must still run"
        );
        assert_exact(k.heat, 0.0, "no heat without reactions 1 and 2");
    }

    /// The drift factors multiply reactions 1 and 2 and nothing else. This is
    /// the half of the `R1F` reuse that a misreading would get wrong.
    #[test]
    fn the_drift_factors_scale_only_the_first_two_reactions() {
        let unpacked = plausible();
        let eq = equilibrium(&unpacked);
        let kelvin = unpacked.reactor.kelvin();

        let base = kinetics(&eq.reactor, kelvin, ReactionDrift::default());
        let drifted = kinetics(
            &eq.reactor,
            kelvin,
            ReactionDrift {
                first: 2.0,
                second: 3.0,
            },
        );
        assert_exact(drifted.rates[0], base.rates[0] * 2.0, "reaction 1 drifts");
        assert_exact(drifted.rates[1], base.rates[1] * 3.0, "reaction 2 drifts");
        assert_exact(drifted.rates[2], base.rates[2], "reaction 3 does not");
        assert_exact(drifted.rates[3], base.rates[3], "reaction 4 does not");
    }

    /// Reaction 4 shares reaction 3's exponential rather than having its own,
    /// so their ratio is a constant that does not depend on temperature.
    #[test]
    fn reaction_four_is_a_fixed_multiple_of_reaction_three() {
        let k = solved();
        // Both were scaled by the same volume and by A, and differ only in the
        // ratio and in D versus E.
        let unpacked = plausible();
        let eq = equilibrium(&unpacked);
        let expected = k.rates[2] / eq.reactor.partial[Component::E]
            * REACTION_4_RATIO
            * eq.reactor.partial[Component::D];
        let relative = (k.rates[3] - expected).abs() / expected.abs();
        assert!(
            relative < 1e-14,
            "reaction 4 is not 0.767488334 times reaction 3's exponential: \
             {relative:e}"
        );
    }
}
