//! The fifty derivatives: component and energy balances, and the valve lags.
//!
//! Ported from `teprob.f:762-811`. Everything the rest of the model computes
//! exists to feed these fifty lines.
//!
//! # The balances
//!
//! Four vessels, each with eight component balances and one energy balance,
//! plus two cooling-water wall temperatures and twelve valve positions.
//!
//! \\[
//!   \\frac{dn_i}{dt} = \\sum_{\\text{in}} \\dot n_i - \\sum_{\\text{out}} \\dot n_i + r_i
//! \\]
//!
//! The reactor is the only vessel with a reaction term. The mixing zone has
//! five inlets and one outlet; the separator has one inlet and three outlets;
//! the stripper is a straight pass-through of the two streams the column
//! produced.
//!
//! Energy balances have the same shape with enthalpy in place of moles, plus
//! the heat duty:
//!
//! \\[
//!   \\frac{dE}{dt} = \\sum_{\\text{in}} h F - \\sum_{\\text{out}} h F + Q
//! \\]
//!
//! and the reactor's carries the heat of reaction as well.
//!
//! # The coolant walls
//!
//! \\[
//!   \\frac{dT_w}{dt} = \\frac{F_w \\times 500.53 \\times (T_{in} - T_w) - Q \\times 10^6 / 1.8}{H_w}
//! \\]
//!
//! `500.53` converts a cooling water flow to a heat capacity rate, and the
//! `10^6 / 1.8` undoes the scaling the enthalpy correlations carry.
//!
//! ## `1.8` here is single precision, and it is the only one in the file
//!
//! `teprob.f:790` and `792` write `1.8` with no `D` suffix. Every other
//! occurrence in `teprob.f` is `1.8D0`: lines `1396`, `1404`, `1464` and
//! `1471`. So this constant alone is stored as 1.7999999523162842, and it
//! divides the entire heat term of both wall balances.
//!
//! The error from getting it wrong is 2.6e-8 relative on that term. This is
//! the same shape as the `273.15` against `273.15D0` split that
//! [`crate::thermo::ABSOLUTE_ZERO_OFFSET`] documents: the original is not
//! consistent, so the precision cannot be inferred from elsewhere in the file
//! and has to be read off the line itself.
//!
//! # The valve lags
//!
//! \\[ \\frac{dv_i}{dt} = \\frac{c_i - v_i}{\\tau_i} \\]
//!
//! with `c` the latched command. The latch itself is *hoisted into the
//! pre-phase*; see [`mod@crate::plant`] and
//! `crates/tepsim-oracle/tests/hoist_valve_latch.rs`, which proves mechanically
//! that moving `teprob.f:793-804` changes no number. Only `teprob.f:805` stays
//! here.
//!
//! # The shutdown freeze is a Class C quirk and is **on** by default
//!
//! `teprob.f:807-811` zeroes all fifty derivatives whenever any shutdown
//! condition holds. That freezes the plant rather than stopping it: the state
//! stops moving, the clock keeps running, and nothing says so.
//!
//! `PLAN.org` classes this C, "behaviour-defining and benchmark-relevant", and
//! says each such entry "gets a full Tier 5 and Tier 6 delta report and an
//! explicit sign-off *before it becomes the default*". The thing needing
//! sign-off is the **fix**, so the faithful freeze is what ships until then and
//! [`QuirkFixes::trip_ends_the_run`] is off by default.
//!
//! That is also the only reading Tier 2 can live with: the adversarial pool
//! contains states that trip, and a port that did not freeze would disagree
//! with the oracle on all fifty components for each of them.
//!
//! The trip is still *reported*, through [`Balances::shutdown`], so a caller
//! never has to infer a freeze from a vector of zeros.
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `YP(1..8)`, `YP(9)` | reactor | component and energy balances |
//! | `YP(10..17)`, `YP(18)` | separator | likewise |
//! | `YP(19..26)`, `YP(27)` | stripper | likewise |
//! | `YP(28..35)`, `YP(36)` | mixing zone | likewise |
//! | `YP(37)`, `YP(38)` | walls | cooling water outlet temperatures |
//! | `YP(39..50)` | valves | first-order lags |

// Every float expression here reproduces `teprob.f`'s association and
// rounding exactly; see `crate::thermo`. The balances are long sums whose
// order is load-bearing.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]

use crate::component::Component;
use crate::constants::single;
use crate::flows::Flows;
use crate::heat::HeatTransfer;
use crate::kinetics::Kinetics;
use crate::measurements::Shutdown;
use crate::state::{Derivative, State};
use crate::stream::Stream;
use crate::streams::Streams;

/// Cooling water heat capacity rate factor (`teprob.f:789`). Single precision.
const COOLANT_HEAT_CAPACITY: f64 = single(500.53);

/// The energy-unit scale the enthalpy correlations carry (`teprob.f:790`).
/// Double precision.
const ENERGY_SCALE: f64 = 1.0e6;

/// Rankine per Celsius degree, as `teprob.f:790` actually stores it.
///
/// **Single precision**, unlike every other `1.8` in the file. See the module
/// documentation: this is 1.7999999523162842, not 1.8.
const RANKINE_PER_CELSIUS: f64 = single(1.8);

/// Reactor coolant wall heat capacity, `HWR` (`teprob.f:1124`). Single.
//
// @port teprob.f:1124-1125
pub const REACTOR_WALL_CAPACITY: f64 = single(7060.);
/// Condenser coolant wall heat capacity, `HWS` (`teprob.f:1125`). Single.
pub const CONDENSER_WALL_CAPACITY: f64 = single(11138.);

/// Valve time constants, `VTAU`, in hours (`teprob.f:1172-1186`).
///
/// The original writes twelve values in *seconds* and then divides the whole
/// array by 3600 in a loop, so the constant here is the quotient rather than
/// the seconds figure. Every literal, including the 3600, is single precision;
/// all thirteen are exactly representable, so the division is exact.
//
// @port teprob.f:1172-1186
pub const VALVE_TIME_CONSTANT: [f64; 12] = [
    single(8.) / single(3600.),
    single(8.) / single(3600.),
    single(6.) / single(3600.),
    single(9.) / single(3600.),
    single(7.) / single(3600.),
    single(5.) / single(3600.),
    single(5.) / single(3600.),
    single(5.) / single(3600.),
    single(120.) / single(3600.),
    single(5.) / single(3600.),
    single(5.) / single(3600.),
    single(5.) / single(3600.),
];

/// The sticking threshold, `VST` (`teprob.f:1106`). Double precision.
///
/// Multiplied by `IVST`, which is zero unless a valve-sticking disturbance is
/// active, so the effective threshold is normally zero and the valves track
/// their command exactly.
pub const VALVE_STICTION: f64 = 2.0;

/// Which Class C quirks are fixed rather than reproduced.
///
/// Every field is off by default, which is the faithful configuration.
/// `PLAN.org` requires a Tier 5 and Tier 6 delta report and an explicit
/// sign-off before any of these becomes the default.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct QuirkFixes {
    /// When `false` (the default), a shutdown zeroes all fifty derivatives,
    /// exactly as `teprob.f:807-811` does.
    ///
    /// When `true`, the derivatives are returned un-zeroed and the caller is
    /// expected to end the run on [`Balances::shutdown`]. That is `PLAN.org`'s
    /// `SimulationOutcome::Trip`, and it is a genuine behaviour change: it is
    /// **blocked on sign-off**, see B-0025b.
    pub trip_ends_the_run: bool,
}

impl QuirkFixes {
    /// Every quirk reproduced rather than fixed, which is the default.
    ///
    /// A `const` constructor as well as `Default`, so a scenario can be built
    /// in a `const fn`.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            trip_ends_the_run: false,
        }
    }
}

/// The fifty derivatives, and whether the plant is down.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Balances {
    /// `YP(1..50)`.
    pub derivative: Derivative,
    /// The magnitude of the largest term entering each balance.
    ///
    /// # Why a derivative carries its own scale
    ///
    /// A balance is inflow minus outflow. Near steady state those nearly
    /// agree, so the result is a small difference of large numbers and its own
    /// magnitude says nothing about how accurately it can be computed.
    /// `YP(2)`, the inert's reactor balance, is a difference of two flows
    /// around 660 whose value is a few parts in ten thousand of either: an
    /// error of 1e-16 *of the flows* is 1e-12 *of the answer*.
    ///
    /// So the error budget of a balance is set by its terms, not by its value,
    /// and this is that budget. Tier 2's gate is the difference from the
    /// Fortran divided by this, which is the decision of 2026-08-27 in
    /// `BACKLOG.org`.
    ///
    /// It is computed here rather than in the harness because only this
    /// function knows what the terms are. Recovering them from the outside
    /// would mean a second implementation of every balance.
    pub scale: Derivative,
    /// The trip, carried alongside rather than inferred from a zero vector.
    pub shutdown: Shutdown,
    /// Whether the derivative was frozen by the trip.
    ///
    /// True only when the plant is down *and* the quirk is being reproduced.
    /// Without this, a caller cannot tell a frozen plant from a steady one.
    pub frozen: bool,
}

/// The two walk-driven coolant inlet temperatures.
///
/// `TCWR` and `TCWS` at `teprob.f:413-414`, from `TESUB8(5)` and `TESUB8(6)`
/// plus `IDV(4)` and `IDV(5)`. Phase 3 machinery, so an input until then.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoolantInlet {
    /// `TCWR`, the reactor coolant inlet temperature.
    pub reactor: f64,
    /// `TCWS`, the condenser coolant inlet temperature.
    pub condenser: f64,
}

impl Default for CoolantInlet {
    /// The nominal operating point, `SZERO(5)` and `SZERO(6)`.
    fn default() -> Self {
        Self {
            reactor: 35.0,
            condenser: 40.0,
        }
    }
}

/// Assemble all fifty derivatives.
///
/// `valve_command` is `VCV`, latched in the pre-phase; see the module
/// documentation.
// @port  teprob.f:762-811
// @delta D-007 class=C teprob.f:807-811
#[allow(
    clippy::too_many_arguments,
    reason = "the balances read the whole model; grouping them into a struct \
              would add a type whose only job is to be unpacked again here"
)]
#[must_use]
pub fn balances(
    y: &State,
    stream_table: &Streams,
    flow: &Flows,
    kinetics: &Kinetics,
    heat: &HeatTransfer,
    shutdown: Shutdown,
    coolant: CoolantInlet,
    valve_command: &[f64; 12],
    fixes: QuirkFixes,
) -> Balances {
    let n = |stream: Stream, c: Component| flow.component[stream][c];
    // Stream 9's enthalpy is the *post-compressor* one. `teprob.f:601` adds
    // the compressor work to `HST(9)` after the stream table is built, and
    // both energy balances that read stream 9 come after that line. Stream 10
    // keeps the pre-bump value, because `teprob.f:562` copied it first.
    //
    // Reading the stream table for stream 9 here gives a mixing-zone energy
    // balance wrong by a factor of 6e8. That happened; see the log for
    // B-0025. `crate::streams` documents the snapshot, and documenting it was
    // not enough.
    let h = |stream: Stream| {
        let enthalpy = if stream == Stream::Recycle {
            flow.recycle_enthalpy
        } else {
            stream_table.enthalpy[stream]
        };
        enthalpy * flow.molar[stream]
    };

    let mut yp = State::default();
    // The largest term entering each balance; see `Balances::scale`.
    let mut scale = State::default();
    /// The magnitude of the largest of a set of terms.
    fn largest(terms: &[f64]) -> f64 {
        let mut worst = 0.0_f64;
        for term in terms {
            let magnitude = term.abs();
            if magnitude > worst {
                worst = magnitude;
            }
        }
        worst
    }

    // teprob.f:762-770. Four component balances per species, in the original's
    // order: reassociating a five-term sum would change the last bits.
    for c in Component::ALL {
        // Reactor: inlet less outlet, plus reaction.
        yp.reactor.moles[c] =
            n(Stream::ReactorInlet, c) - n(Stream::ReactorOutlet, c) + kinetics.production[c];
        // Separator: reactor effluent in, recycle, purge and underflow out.
        yp.separator.moles[c] = n(Stream::ReactorOutlet, c)
            - n(Stream::Recycle, c)
            - n(Stream::Purge, c)
            - n(Stream::SeparatorUnderflow, c);
        // Stripper: the column's downflow in, product out.
        yp.stripper.moles[c] = n(Stream::StripperDownflow, c) - n(Stream::Product, c);
        // Mixing zone: three feeds, the stripper overhead and the recycle in,
        // and its own outlet out.
        yp.mixing.moles[c] = n(Stream::DFeed, c)
            + n(Stream::EFeed, c)
            + n(Stream::AFeed, c)
            + n(Stream::StripperOverhead, c)
            + n(Stream::Recycle, c)
            - n(Stream::MixingZoneOutlet, c);

        scale.reactor.moles[c] = largest(&[
            n(Stream::ReactorInlet, c),
            n(Stream::ReactorOutlet, c),
            kinetics.production[c],
        ]);
        scale.separator.moles[c] = largest(&[
            n(Stream::ReactorOutlet, c),
            n(Stream::Recycle, c),
            n(Stream::Purge, c),
            n(Stream::SeparatorUnderflow, c),
        ]);
        scale.stripper.moles[c] = largest(&[n(Stream::StripperDownflow, c), n(Stream::Product, c)]);
        scale.mixing.moles[c] = largest(&[
            n(Stream::DFeed, c),
            n(Stream::EFeed, c),
            n(Stream::AFeed, c),
            n(Stream::StripperOverhead, c),
            n(Stream::Recycle, c),
            n(Stream::MixingZoneOutlet, c),
        ]);
    }

    // teprob.f:771-772. The reactor is the only vessel with a reaction
    // enthalpy, and `QUR` is heat *removed*, so it enters positive here
    // because `UAR*(TWR-TCR)` is already negative when the coil is cooling.
    yp.reactor.energy =
        h(Stream::ReactorInlet) - h(Stream::ReactorOutlet) + kinetics.heat + heat.reactor_duty;

    // teprob.f:773-777
    yp.separator.energy = h(Stream::ReactorOutlet)
        - h(Stream::Recycle)
        - h(Stream::Purge)
        - h(Stream::SeparatorUnderflow)
        + heat.condenser_duty;

    // teprob.f:778-782. The mixed feed enters the stripper directly, not
    // through the mixing zone.
    yp.stripper.energy = h(Stream::AcFeed) + h(Stream::SeparatorUnderflow)
        - h(Stream::StripperOverhead)
        - h(Stream::Product)
        + heat.stripper_duty;

    // teprob.f:783-788. No heat duty: the mixing zone is adiabatic.
    yp.mixing.energy = h(Stream::DFeed)
        + h(Stream::EFeed)
        + h(Stream::AFeed)
        + h(Stream::StripperOverhead)
        + h(Stream::Recycle)
        - h(Stream::MixingZoneOutlet);

    // teprob.f:789-792. `RANKINE_PER_CELSIUS` is single precision here and
    // nowhere else in the file; see the module documentation.
    yp.reactor_cw_out_c =
        (flow.reactor_coolant * COOLANT_HEAT_CAPACITY * (coolant.reactor - y.reactor_cw_out_c)
            - heat.reactor_duty * ENERGY_SCALE / RANKINE_PER_CELSIUS)
            / REACTOR_WALL_CAPACITY;
    yp.condenser_cw_out_c = (flow.condenser_coolant
        * COOLANT_HEAT_CAPACITY
        * (coolant.condenser - y.condenser_cw_out_c)
        - heat.condenser_duty * ENERGY_SCALE / RANKINE_PER_CELSIUS)
        / CONDENSER_WALL_CAPACITY;

    scale.reactor.energy = largest(&[
        h(Stream::ReactorInlet),
        h(Stream::ReactorOutlet),
        kinetics.heat,
        heat.reactor_duty,
    ]);
    scale.separator.energy = largest(&[
        h(Stream::ReactorOutlet),
        h(Stream::Recycle),
        h(Stream::Purge),
        h(Stream::SeparatorUnderflow),
        heat.condenser_duty,
    ]);
    scale.stripper.energy = largest(&[
        h(Stream::AcFeed),
        h(Stream::SeparatorUnderflow),
        h(Stream::StripperOverhead),
        h(Stream::Product),
        heat.stripper_duty,
    ]);
    scale.mixing.energy = largest(&[
        h(Stream::DFeed),
        h(Stream::EFeed),
        h(Stream::AFeed),
        h(Stream::StripperOverhead),
        h(Stream::Recycle),
        h(Stream::MixingZoneOutlet),
    ]);

    // The wall balances divide by their heat capacity, so the terms are taken
    // after that division: the scale has to be in the same units as the answer.
    scale.reactor_cw_out_c = largest(&[
        flow.reactor_coolant * COOLANT_HEAT_CAPACITY * coolant.reactor / REACTOR_WALL_CAPACITY,
        flow.reactor_coolant * COOLANT_HEAT_CAPACITY * y.reactor_cw_out_c / REACTOR_WALL_CAPACITY,
        heat.reactor_duty * ENERGY_SCALE / RANKINE_PER_CELSIUS / REACTOR_WALL_CAPACITY,
    ]);
    scale.condenser_cw_out_c = largest(&[
        flow.condenser_coolant * COOLANT_HEAT_CAPACITY * coolant.condenser
            / CONDENSER_WALL_CAPACITY,
        flow.condenser_coolant * COOLANT_HEAT_CAPACITY * y.condenser_cw_out_c
            / CONDENSER_WALL_CAPACITY,
        heat.condenser_duty * ENERGY_SCALE / RANKINE_PER_CELSIUS / CONDENSER_WALL_CAPACITY,
    ]);

    // teprob.f:805. The latch at 799-804 is hoisted into the pre-phase; only
    // the derivative stays here.
    for i in 0..12 {
        yp.valve_pos[i] = (valve_command[i] - y.valve_pos[i]) / VALVE_TIME_CONSTANT[i];
        scale.valve_pos[i] = largest(&[
            valve_command[i] / VALVE_TIME_CONSTANT[i],
            y.valve_pos[i] / VALVE_TIME_CONSTANT[i],
        ]);
    }

    // teprob.f:807-811. Class C; see the module documentation for why this is
    // the default.
    let frozen = shutdown.is_tripped() && !fixes.trip_ends_the_run;
    if frozen {
        yp = State::default();
        // A frozen derivative is exactly zero on both sides, so its budget is
        // zero too: there is nothing to be accurate about.
        scale = State::default();
    }

    Balances {
        derivative: Derivative::new(yp),
        scale: Derivative::new(scale),
        shutdown,
        frozen,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::equilibrium::equilibrium;
    use crate::flows::{FlowDrift, flows};
    use crate::heat::{HeatDrift, heat_transfer};
    use crate::kinetics::{ReactionDrift, kinetics};
    use crate::measurements::{ShutdownCause, measurements};
    use crate::state::N_STATES;
    use crate::streams::{FeedConditions, streams};
    use crate::stripper::stripper;
    use crate::testing::assert_exact;
    use crate::vessels::{TemperatureSeeds, unpack};

    struct Fixture {
        y: State,
        table: Streams,
        flow: Flows,
        kinetics: Kinetics,
        heat: HeatTransfer,
    }

    fn plant() -> Fixture {
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
        y.valve_pos = [50.0; 12];
        y.reactor_cw_out_c = 94.6;
        y.condenser_cw_out_c = 77.3;
        let unpacked = unpack(&y, TemperatureSeeds::default()).expect("converges");
        let eq = equilibrium(&unpacked);
        let mut table = streams(&unpacked, &eq, &FeedConditions::default());
        let mut flow = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        let _ = stripper(&mut table, &mut flow, unpacked.stripper.celsius);
        let heat = heat_transfer(&y, &unpacked, &table, &flow, HeatDrift::default());
        let kin = kinetics(
            &eq.reactor,
            unpacked.reactor.kelvin(),
            ReactionDrift::default(),
        );
        Fixture {
            y,
            table,
            flow,
            kinetics: kin,
            heat,
        }
    }

    fn solve(f: &Fixture, shutdown: Shutdown, fixes: QuirkFixes) -> Balances {
        balances(
            &f.y,
            &f.table,
            &f.flow,
            &f.kinetics,
            &f.heat,
            shutdown,
            CoolantInlet::default(),
            &[50.0; 12],
            fixes,
        )
    }

    /// `1.8` at `teprob.f:790` is single precision and every other `1.8` in
    /// the file is not. This is the whole hazard of the module in one
    /// assertion.
    #[test]
    fn the_rankine_ratio_is_single_precision_here_and_double_elsewhere() {
        assert_exact(RANKINE_PER_CELSIUS, 1.8_f32 as f64, "teprob.f:790");
        const { assert!(RANKINE_PER_CELSIUS < 1.8_f64) };
        // `crate::thermo` carries the double one, from `teprob.f:1396`. If
        // these ever became equal, one of the two was typed wrong.
        assert!(
            RANKINE_PER_CELSIUS.to_bits() != 1.8_f64.to_bits(),
            "the single and double 1.8 agree, which cannot be right"
        );
    }

    /// A frozen plant returns exactly fifty zeros, and says so rather than
    /// leaving the caller to infer it.
    #[test]
    fn a_trip_freezes_every_derivative_by_default() {
        let f = plant();
        // Any cause will do; the freeze does not look at which.
        let measured = measurements(
            &f.y,
            &unpack(&f.y, TemperatureSeeds::default()).expect("converges"),
            &f.table,
            &f.flow,
            &f.heat,
            (1.0e6, 1.0e6, 1.0e6),
        );
        let tripped = measured.shutdown;
        assert!(tripped.is_tripped());

        let b = solve(&f, tripped, QuirkFixes::default());
        assert!(b.frozen);
        assert_eq!(b.shutdown, tripped);
        for (slot, value) in b.derivative.to_flat().iter().enumerate() {
            assert_eq!(
                value.to_bits(),
                0.0_f64.to_bits(),
                "YP({}) is not zero on a frozen plant",
                slot + 1
            );
        }
    }

    /// The fix returns the un-frozen derivative and is not the default.
    /// `PLAN.org` requires sign-off before it becomes one.
    #[test]
    fn the_fix_is_off_by_default_and_changes_the_answer_when_on() {
        let f = plant();
        let tripped = measurements(
            &f.y,
            &unpack(&f.y, TemperatureSeeds::default()).expect("converges"),
            &f.table,
            &f.flow,
            &f.heat,
            (1.0e6, 1.0e6, 1.0e6),
        )
        .shutdown;

        assert!(
            !QuirkFixes::default().trip_ends_the_run,
            "the default is faithful"
        );

        let fixed = solve(
            &f,
            tripped,
            QuirkFixes {
                trip_ends_the_run: true,
            },
        );
        assert!(!fixed.frozen);
        assert!(fixed.shutdown.is_tripped(), "the trip is still reported");
        let moving = fixed
            .derivative
            .to_flat()
            .iter()
            .filter(|v| v.to_bits() != 0.0_f64.to_bits())
            .count();
        assert!(
            moving > 20,
            "only {moving} of 50 derivatives are non-zero with the fix on, so \
             the flag is not actually changing anything"
        );
    }

    /// A healthy plant is never frozen, whatever the flag says.
    #[test]
    fn a_healthy_plant_is_never_frozen() {
        let f = plant();
        let healthy = Shutdown::default();
        assert!(!healthy.is_tripped());
        for fixes in [
            QuirkFixes::default(),
            QuirkFixes {
                trip_ends_the_run: true,
            },
        ] {
            let b = solve(&f, healthy, fixes);
            assert!(!b.frozen);
        }
    }

    /// The inert accumulates in the mixing zone and leaves only through the
    /// purge, which is the entire reason the plant has one. Its reactor
    /// balance therefore has no reaction term.
    #[test]
    fn the_inert_has_no_reaction_term_in_the_reactor_balance() {
        let f = plant();
        let b = solve(&f, Shutdown::default(), QuirkFixes::default());
        let flat = b.derivative.to_flat();
        assert_exact(
            flat[Component::B.index()],
            f.flow.component[Stream::ReactorInlet][Component::B]
                - f.flow.component[Stream::ReactorOutlet][Component::B],
            "YP(2) is inlet less outlet, with nothing added",
        );
        assert_exact(f.kinetics.production[Component::B], 0.0, "CRXR(2)");
    }

    /// Each balance uses its own vessel's streams. Crossing two would leave
    /// the plant running and wrong, and is the mistake this layer is most
    /// exposed to.
    #[test]
    fn each_vessel_balances_its_own_streams() {
        let f = plant();
        let b = solve(&f, Shutdown::default(), QuirkFixes::default());
        let flat = b.derivative.to_flat();
        let c = Component::D;
        // Reactor: YP(4).
        assert_exact(
            flat[3],
            f.flow.component[Stream::ReactorInlet][c] - f.flow.component[Stream::ReactorOutlet][c]
                + f.kinetics.production[c],
            "YP(4)",
        );
        // Stripper: YP(22) is component 4 of the stripper block, YP(19..26).
        assert_exact(
            flat[21],
            f.flow.component[Stream::StripperDownflow][c] - f.flow.component[Stream::Product][c],
            "YP(22)",
        );
    }

    /// The valve lags drive each position toward its command with the valve's
    /// own time constant. Valve 9 is twenty-four times slower than the rest,
    /// which is easy to lose in a table of twelve.
    #[test]
    fn the_valve_lags_use_per_valve_time_constants() {
        let f = plant();
        let b = solve(&f, Shutdown::default(), QuirkFixes::default());
        let flat = b.derivative.to_flat();
        for i in 0..12 {
            assert_exact(
                flat[38 + i],
                (50.0 - f.y.valve_pos[i]) / VALVE_TIME_CONSTANT[i],
                "valve lag",
            );
        }
        // Valve 9 is the slow one, at 120 seconds against 5 to 9.
        let slowest = VALVE_TIME_CONSTANT
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert_exact(VALVE_TIME_CONSTANT[8], slowest, "valve 9 is the slowest");
        assert_exact(VALVE_TIME_CONSTANT[8], 120.0 / 3600.0, "120 seconds");
    }

    /// A valve already at its command does not move.
    #[test]
    fn a_valve_at_its_command_has_no_rate() {
        let f = plant();
        let b = balances(
            &f.y,
            &f.table,
            &f.flow,
            &f.kinetics,
            &f.heat,
            Shutdown::default(),
            CoolantInlet::default(),
            &f.y.valve_pos,
            QuirkFixes::default(),
        );
        for value in &b.derivative.to_flat()[38..N_STATES] {
            assert_exact(*value, 0.0, "a valve at its command should be still");
        }
    }

    /// A shutdown cause is still reported on a frozen plant. Without that, a
    /// caller sees fifty zeros and cannot tell a trip from a steady state.
    #[test]
    fn a_frozen_plant_still_reports_why() {
        let f = plant();
        let tripped = measurements(
            &f.y,
            &unpack(&f.y, TemperatureSeeds::default()).expect("converges"),
            &f.table,
            &f.flow,
            &f.heat,
            (1.0e6, 1.0e6, 1.0e6),
        )
        .shutdown;
        let b = solve(&f, tripped, QuirkFixes::default());
        assert!(b.frozen);
        assert!(b.shutdown.first().is_some());
        assert!(
            b.shutdown.holds(ShutdownCause::ReactorPressureHigh),
            "the cause survives the freeze"
        );
    }
}
