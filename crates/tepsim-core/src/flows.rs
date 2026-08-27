//! Valve-lagged flows, pressure-driven flows, and the compressor.
//!
//! Ported from `teprob.f:565-613`, with the valve ranges and compressor
//! constants from `teprob.f:1109-1117` and `1170-1171`. This is where the
//! plant's hydraulics live: everything before it describes vessel contents,
//! everything after it moves them around.
//!
//! # Three kinds of flow
//!
//! **Valve-lagged.** Five flows are proportional to a lagged valve position:
//!
//! \\[ F = \\frac{v \\, R}{100} \\]
//!
//! with `v` the valve position (a state, `VPOS`) and `R` the valve's range
//! ([`VALVE_RANGE`]). Streams 1, 2, 3, 4, 11 and 13, plus the two cooling
//! water flows and the stripper steam coefficient.
//!
//! **Pressure-driven.** Three flows come from a square-root resistance across
//! a pressure difference, converted from mass to moles by the stream's mean
//! molecular weight:
//!
//! \\[ F = \\frac{k \\sqrt{\\Delta P}}{\\overline{M}} \\]
//!
//! Every one of them clamps \\(\\Delta P\\) at zero first, so a reversed
//! pressure gradient gives no flow rather than a `NaN` from the square root.
//!
//! **The compressor.** A fixed-speed centrifugal machine with a cubic
//! pressure-ratio curve, a recycle valve bleeding flow back, and a floor:
//!
//! \\[
//!   \\dot m = \\max\\!\\left(
//!     F_{\\max}\\left(1 + \\frac{1 - r^3}{1.197}\\right)
//!     - v_5 \\, 53.349 \\sqrt{\\Delta P},
//!     \; 10^{-3}
//!   \\right)
//! \\]
//!
//! # Five clamps, and one of them cannot be reached
//!
//! | Clamp | Line | Reachable from a physical state |
//! |---|---|---|
//! | `PTV-PTR` at zero | 577 | yes, adversarial |
//! | `PTR-PTS` at zero | 581 | yes, adversarial |
//! | `PTS-760` at zero | 586 | **no** |
//! | `PTV-PTS` at zero | 597 | yes, adversarial |
//! | `PR` low at 1, high at [`MAX_PRESSURE_RATIO`] | 590-591 | yes, adversarial |
//! | compressor flow floor at 1e-3 | 599 | yes, adversarial |
//!
//! The purge clamp is the interesting one. `PTS` is a sum of eight partial
//! pressures in a vessel that always holds material, and on the whole sampled
//! domain it floors around 811 mmHg against a threshold of 760. No trajectory
//! state, perturbation or adversarial boundary will ever take that branch, so
//! it is covered by a unit test at a composition chosen to reach it rather than
//! by the differential. A branch that no test enters is indistinguishable from
//! a branch that is wrong.
//!
//! # `PR**3` is an integer power
//!
//! `teprob.f:593` writes `PR**3` with an integer exponent, and gfortran emits
//! multiplications for it rather than a call to `pow`. Routing it through
//! [`crate::math::pow`] would change the last bits, so it is written out as
//! `pr * pr * pr`. The same applies to `AGSP**2` and `(FTM(8)/3528.73)**4` in
//! B-0023.
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `FTM(13)` | [`Flows::molar`] | molar flow per stream |
//! | `FCM(8,13)` | [`Flows::component`] | component molar flow |
//! | `UAC` | [`Flows::steam_coefficient`] | stripper steam heat transfer |
//! | `FWR`, `FWS` | [`Flows::reactor_coolant`] etc. | cooling water flows |
//! | `AGSP` | [`Flows::agitator`] | agitator speed, fraction of nominal |
//! | `CPDH` | [`Flows::compressor_work`] | compressor enthalpy bump |
//! | `VRNG` | [`VALVE_RANGE`] | valve full-travel capacity |
//! | `CPFLMX`, `CPPRMX` | [`MAX_COMPRESSOR_FLOW`] etc. | compressor limits |

// Every float expression here reproduces `teprob.f`'s association and
// rounding exactly; see `crate::thermo`.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]
// `FLMS=FLMS-...` at teprob.f:598, written the way the listing writes it so
// this module can be checked against the source a line at a time. See
// `crate::thermo`.
#![allow(
    clippy::assign_op_pattern,
    reason = "operand order is transcribed from the Fortran, not chosen"
)]

use crate::component::{ByComponent, Component};
use crate::constants::single;
use crate::equilibrium::Equilibrium;
use crate::math::sqrt;
use crate::state::State;
use crate::stream::{ByStream, Stream};
use crate::streams::Streams;
use crate::vessels::Unpacked;

/// Valve full-travel capacities, `VRNG` (`teprob.f:1109-1117`).
///
/// Indexed by valve number minus one. Valves 5, 6 and 12 have no range: the
/// recycle and purge valves carry their coefficient inline at
/// `teprob.f:587` and `598`, and valve 12 drives the agitator, which is a
/// speed rather than a flow. Those three slots are never written by the
/// original and never read, so they are zero here for the same reason.
///
/// All nine literals are single precision; none has a `D` suffix.
//
// @port teprob.f:1109-1117
pub const VALVE_RANGE: [f64; 12] = [
    single(400.00),  // VRNG(1),  D feed
    single(400.00),  // VRNG(2),  E feed
    single(100.00),  // VRNG(3),  A feed
    single(1500.00), // VRNG(4),  A and C feed
    0.0,             // VRNG(5),  never written
    0.0,             // VRNG(6),  never written
    single(1500.00), // VRNG(7),  separator underflow
    single(1000.00), // VRNG(8),  product
    single(0.03),    // VRNG(9),  stripper steam
    single(1000.),   // VRNG(10), reactor coolant
    single(1200.0),  // VRNG(11), condenser coolant
    0.0,             // VRNG(12), never written
];

/// `CPFLMX`, the compressor's maximum mass flow (`teprob.f:1170`).
///
/// Single precision, though 280275 is exactly representable in 24 bits.
//
// @port teprob.f:1170
pub const MAX_COMPRESSOR_FLOW: f64 = single(280275.);

/// `CPPRMX`, the compressor's maximum pressure ratio (`teprob.f:1171`).
///
/// Single precision, and 1.3 is *not* exactly representable, so this is
/// 1.2999999523162842. It appears in a comparison and as the clamped value, so
/// the difference is observable in both.
//
// @port teprob.f:1171
pub const MAX_PRESSURE_RATIO: f64 = single(1.3);

/// The compressor curve's slope divisor, `teprob.f:592`. Double precision.
const COMPRESSOR_CURVE_DIVISOR: f64 = 1.197;

/// Minimum compressor mass flow, `teprob.f:599`. Double precision.
///
/// Prevents a fully-open recycle valve from driving the flow to zero and
/// producing a division by zero at `teprob.f:600-601`.
const MIN_COMPRESSOR_FLOW: f64 = 1.0e-3;

/// Resistance coefficient, mixing zone to reactor (`teprob.f:578`).
const RESISTANCE_MIXING_TO_REACTOR: f64 = 1937.6;
/// Resistance coefficient, reactor to separator (`teprob.f:582`).
const RESISTANCE_REACTOR_TO_SEPARATOR: f64 = 4574.21;
/// Resistance coefficient per unit purge valve opening (`teprob.f:587`).
const RESISTANCE_PURGE: f64 = 0.151169;
/// Resistance coefficient per unit recycle valve opening (`teprob.f:598`).
const RESISTANCE_RECYCLE: f64 = 53.349;

/// Atmospheric pressure the purge discharges against, mmHg
/// (`teprob.f:585`). Single precision, and exactly representable.
const ATMOSPHERIC: f64 = single(760.0);

/// The agitator speed offset and scale, `teprob.f:575`. Both single precision.
const AGITATOR_OFFSET: f64 = single(150.0);
/// Percent-of-travel divisor, used throughout this range. Single precision.
const PERCENT: f64 = single(100.0);

/// The two disturbance walk channels this range reads.
///
/// `TESUB8(9, t)` and `TESUB8(12, t)` at `teprob.f:572` and `583`. Walk
/// machinery is Phase 3, so until then they are inputs, exactly as
/// [`crate::kinetics::ReactionDrift`] and [`crate::streams::FeedConditions`]
/// are.
///
/// [`Default`] is no drift: both channels sit at zero, which is `SZERO` for
/// these two.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FlowDrift {
    /// `TESUB8(9, t)`: fractional drift on the stripper steam coefficient.
    ///
    /// Enters as `1 + drift`, so zero means nominal.
    pub steam_capacity: f64,
    /// `TESUB8(12, t)`: fractional restriction on the reactor outlet.
    ///
    /// Enters as `1 - 0.25 * drift`, so zero means unrestricted and the
    /// quarter-scaling caps the worst case at a 25% loss.
    pub reactor_outlet: f64,
}

/// Everything `teprob.f:565-613` produces.
///
/// Streams 5, 7 and 12 are left at zero: they are the stripper's (B-0022).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Flows {
    /// `FTM`: molar flow leaving each stream.
    pub molar: ByStream<f64>,
    /// `FCM`: molar flow of each component in each stream.
    pub component: ByStream<ByComponent<f64>>,
    /// `UAC`: the stripper's steam-side heat transfer coefficient.
    pub steam_coefficient: f64,
    /// `FWR`: reactor cooling water flow.
    pub reactor_coolant: f64,
    /// `FWS`: condenser cooling water flow.
    pub condenser_coolant: f64,
    /// `AGSP`: agitator speed as a fraction of nominal.
    pub agitator: f64,
    /// `CPDH`: the compressor's enthalpy contribution, before division by the
    /// recycle molar flow.
    pub compressor_work: f64,
    /// `HST(9)` after `teprob.f:601` adds the compressor work.
    ///
    /// Kept here rather than written back into [`Streams`] because the
    /// pre-bump value is still needed: `HST(10)` is a copy of it, taken at
    /// `teprob.f:562`. See [`mod@crate::streams`].
    pub recycle_enthalpy: f64,
    /// Whether the compressor's pressure ratio was clamped, and which way.
    pub pressure_ratio_clamp: RatioClamp,
}

/// Which way the compressor pressure ratio hit its limits at
/// `teprob.f:590-591`.
///
/// Reported rather than merely applied, because Tier 2 has no other way to
/// tell whether a sampled state exercised the branch, and an unexercised
/// branch is not evidence of anything.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RatioClamp {
    /// `PTV/PTS` was inside `[1, CPPRMX]` and passed through.
    None,
    /// Below 1, clamped up. The separator is at higher pressure than the
    /// mixing zone, so the compressor would be running backwards.
    Low,
    /// Above `CPPRMX`, clamped down to the machine's limit.
    High,
}

/// Clamp a pressure difference at zero, as the four `IF(DLP.LT.0.0)DLP=0.0`
/// lines do.
///
/// Written as a function rather than four copies because the differential
/// cannot tell them apart, and a reader checking one against the listing
/// should not have to check four.
fn no_reverse_flow(difference: f64) -> f64 {
    if difference < 0.0 { 0.0 } else { difference }
}

/// Compute every flow this range produces.
///
/// `streams` supplies `XMWS` and `XST`; the vapour pressures come from `eq`,
/// the valve positions from `y`, and `disturbances` is the twenty `IDV` flags
/// already clamped to 0 or 1.
// @port  teprob.f:565-613
// @delta D-004 class=B teprob.f:568-569
#[must_use]
pub fn flows(
    y: &State,
    unpacked: &Unpacked,
    eq: &Equilibrium,
    stream_table: &Streams,
    disturbances: &[f64; 20],
    drift: FlowDrift,
) -> Flows {
    let valve = |n: usize| y.valve_pos[n - 1];
    let range = |n: usize| VALVE_RANGE[n - 1];
    let idv = |n: usize| disturbances[n - 1];
    let mut molar = ByStream::new([0.0; Stream::COUNT]);

    // teprob.f:565-566, 570-571. Straight valve-lagged flows.
    molar[Stream::DFeed] = valve(1) * range(1) / PERCENT;
    molar[Stream::EFeed] = valve(2) * range(2) / PERCENT;
    molar[Stream::SeparatorUnderflow] = valve(7) * range(7) / PERCENT;
    molar[Stream::Product] = valve(8) * range(8) / PERCENT;

    // teprob.f:567. IDV(6) is A feed loss: the flag shuts the stream off
    // entirely rather than reducing it.
    molar[Stream::AFeed] = valve(3) * (1.0 - idv(6)) * range(3) / PERCENT;

    // teprob.f:568-569. IDV(7) is C header pressure loss, a 20% reduction.
    //
    // The trailing `+1.D-10` is delta D-004: it keeps `FTM(4)` non-zero so
    // that nothing downstream divides by it, and it is a Class B quirk because
    // it silently adds a flow the plant does not have.
    molar[Stream::AcFeed] =
        valve(4) * (1.0 - idv(7) * 0.2) * range(4) / PERCENT + FEED_FLOW_EPSILON;

    // teprob.f:572-574.
    let steam_coefficient = valve(9) * range(9) * (1.0 + drift.steam_capacity) / PERCENT;
    let reactor_coolant = valve(10) * range(10) / PERCENT;
    let condenser_coolant = valve(11) * range(11) / PERCENT;

    // teprob.f:575. Not a flow: the agitator runs between 1.5 and 2.5 times
    // nominal over the valve's travel, and never reaches zero.
    let agitator = (valve(12) + AGITATOR_OFFSET) / PERCENT;

    let ptr = eq.reactor.pressure;
    let pts = eq.separator.pressure;
    let ptv = eq.mixing_pressure;
    let weight = |stream: Stream| stream_table.molar_mass[stream];

    // teprob.f:576-579. Mixing zone to reactor.
    let mass = RESISTANCE_MIXING_TO_REACTOR * sqrt(no_reverse_flow(ptv - ptr));
    molar[Stream::MixingZoneOutlet] = mass / weight(Stream::MixingZoneOutlet);

    // teprob.f:580-584. Reactor to separator, throttled by IDV(20).
    let mass = RESISTANCE_REACTOR_TO_SEPARATOR
        * sqrt(no_reverse_flow(ptr - pts))
        * (1.0 - 0.25 * drift.reactor_outlet);
    molar[Stream::ReactorOutlet] = mass / weight(Stream::ReactorOutlet);

    // teprob.f:585-588. Purge to atmosphere. See the module documentation for
    // why this clamp is unreachable from a sampled state.
    let mass = valve(6) * RESISTANCE_PURGE * sqrt(no_reverse_flow(pts - ATMOSPHERIC));
    molar[Stream::Purge] = mass / weight(Stream::Purge);

    // teprob.f:589-591. The compressor's operating point on its curve.
    let raw_ratio = ptv / pts;
    let (ratio, pressure_ratio_clamp) = if raw_ratio < 1.0 {
        (1.0, RatioClamp::Low)
    } else if raw_ratio > MAX_PRESSURE_RATIO {
        (MAX_PRESSURE_RATIO, RatioClamp::High)
    } else {
        (raw_ratio, RatioClamp::None)
    };

    // teprob.f:592-593. `PR**3` is an integer power and expands to
    // multiplication; see the module documentation.
    let slope = MAX_COMPRESSOR_FLOW / COMPRESSOR_CURVE_DIVISOR;
    let mut mass = MAX_COMPRESSOR_FLOW + slope * (1.0 - ratio * ratio * ratio);

    // teprob.f:594-595. Note `273.15D0` here is *double* precision, unlike the
    // one in `TESUB2`; the original is not consistent about it and this line
    // has to be read rather than inferred.
    let compressor_work =
        mass * (unpacked.separator.celsius + 273.15) * 1.8e-6 * 1.9872 * (ptv - pts)
            / (weight(Stream::Recycle) * pts);

    // teprob.f:596-599. The recycle valve bleeds flow back, and the result has
    // a floor so that the division at 600-601 cannot blow up.
    mass = mass - valve(5) * RESISTANCE_RECYCLE * sqrt(no_reverse_flow(ptv - pts));
    if mass < MIN_COMPRESSOR_FLOW {
        mass = MIN_COMPRESSOR_FLOW;
    }
    // teprob.f:600-601.
    molar[Stream::Recycle] = mass / weight(Stream::Recycle);
    let recycle_enthalpy =
        stream_table.enthalpy[Stream::Recycle] + compressor_work / molar[Stream::Recycle];

    // teprob.f:602-613. The same ten streams the stream table covers.
    let mut component = ByStream::new([ByComponent::new([0.0; Component::COUNT]); Stream::COUNT]);
    for stream in crate::streams::ASSEMBLED_STREAMS {
        for c in Component::ALL {
            component[stream][c] = stream_table.composition[stream][c] * molar[stream];
        }
    }

    Flows {
        molar,
        component,
        steam_coefficient,
        reactor_coolant,
        condenser_coolant,
        agitator,
        compressor_work,
        recycle_enthalpy,
        pressure_ratio_clamp,
    }
}

/// The `1.D-10` added to the mixed feed at `teprob.f:569`.
///
/// Delta D-004. It exists so that `FTM(4)` is never exactly zero, which
/// matters because a closed valve 4 would otherwise give a zero flow that
/// downstream code divides by. It is Class B rather than Class A because it is
/// numerically observable: the plant receives 1e-10 lbmol/h of A/B/C feed that
/// nothing accounts for, forever.
pub const FEED_FLOW_EPSILON: f64 = 1.0e-10;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::component::Composition;
    use crate::equilibrium::equilibrium;
    use crate::streams::{FeedConditions, streams};
    use crate::testing::assert_exact;
    use crate::vessels::{TemperatureSeeds, unpack};

    fn plant() -> (State, Unpacked, Equilibrium, Streams) {
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
        let unpacked = unpack(&y, TemperatureSeeds::default()).expect("converges");
        let eq = equilibrium(&unpacked);
        let table = streams(&unpacked, &eq, &FeedConditions::default());
        (y, unpacked, eq, table)
    }

    fn solved() -> Flows {
        let (y, unpacked, eq, table) = plant();
        flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default())
    }

    /// The purge clamp cannot be reached from any state the differential
    /// samples, so it is exercised here instead. Without this the branch is
    /// untested, and an untested branch is indistinguishable from a wrong one.
    ///
    /// A separator holding almost nothing has a total pressure below
    /// atmospheric, so the purge cannot discharge.
    #[test]
    fn the_purge_stops_when_the_separator_falls_below_atmospheric() {
        let (y, unpacked, eq, table) = plant();
        let mut eq = eq;
        // Below the 760 mmHg threshold, so `PTS - 760` is negative.
        eq.separator.pressure = 700.0;
        assert!(eq.separator.pressure < ATMOSPHERIC);

        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        assert_exact(
            f.molar[Stream::Purge],
            0.0,
            "a reversed gradient gives no purge flow, not a NaN",
        );
        assert!(
            f.molar[Stream::Purge].is_finite(),
            "the clamp exists to stop sqrt of a negative"
        );

        // And just above it, the purge flows again.
        eq.separator.pressure = 900.0;
        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        assert!(
            f.molar[Stream::Purge] > 0.0,
            "above atmospheric the purge must flow"
        );
    }

    /// All four `DLP` clamps, checked as a group: a reversed gradient gives
    /// zero flow and never a `NaN`.
    #[test]
    fn every_reversed_pressure_gradient_gives_zero_flow_rather_than_a_nan() {
        let (y, unpacked, eq, table) = plant();
        let mut eq = eq;
        // Invert every gradient at once.
        eq.mixing_pressure = 100.0;
        eq.reactor.pressure = 200.0;
        eq.separator.pressure = 300.0;

        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        for stream in [
            Stream::MixingZoneOutlet,
            Stream::ReactorOutlet,
            Stream::Purge,
        ] {
            assert_exact(f.molar[stream], 0.0, stream.name());
        }
        assert!(
            f.molar[Stream::Recycle].is_finite(),
            "the compressor must survive a fully reversed plant"
        );
    }

    /// The pressure ratio clamps both ways, and the port reports which.
    #[test]
    fn the_compressor_pressure_ratio_clamps_at_both_ends() {
        let (y, unpacked, eq, table) = plant();
        let mut eq = eq;

        eq.separator.pressure = 1000.0;
        eq.mixing_pressure = 500.0; // ratio 0.5
        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        assert_eq!(f.pressure_ratio_clamp, RatioClamp::Low);

        eq.mixing_pressure = 5000.0; // ratio 5.0
        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        assert_eq!(f.pressure_ratio_clamp, RatioClamp::High);

        eq.mixing_pressure = 1100.0; // ratio 1.1, inside the band
        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        assert_eq!(f.pressure_ratio_clamp, RatioClamp::None);
    }

    /// `CPPRMX` is single precision, so it is just below 1.3. A state at
    /// exactly 1.3 must therefore clamp, which the double literal would not
    /// do.
    #[test]
    fn the_maximum_pressure_ratio_is_the_single_precision_value() {
        assert_exact(MAX_PRESSURE_RATIO, 1.3_f32 as f64, "CPPRMX");
        // 1.3 rounds *down* in binary32, so a ratio of exactly 1.3 is above
        // `CPPRMX` and clamps. Constant-folded, hence the const block.
        const { assert!(MAX_PRESSURE_RATIO < 1.3_f64) };
    }

    /// `IDV(6)` shuts the A feed off completely; `IDV(7)` takes 20% off the
    /// mixed feed. Swapping the two is a plausible slip.
    #[test]
    fn the_two_feed_disturbances_do_different_things() {
        let (y, unpacked, eq, table) = plant();
        let base = solved();

        let mut idv = [0.0; 20];
        idv[5] = 1.0; // IDV(6)
        let f = flows(&y, &unpacked, &eq, &table, &idv, FlowDrift::default());
        assert_exact(f.molar[Stream::AFeed], 0.0, "IDV(6) stops the A feed dead");
        assert_exact(
            f.molar[Stream::AcFeed],
            base.molar[Stream::AcFeed],
            "IDV(6) does not touch the mixed feed",
        );

        let mut idv = [0.0; 20];
        idv[6] = 1.0; // IDV(7)
        let f = flows(&y, &unpacked, &eq, &table, &idv, FlowDrift::default());
        assert_exact(
            f.molar[Stream::AFeed],
            base.molar[Stream::AFeed],
            "IDV(7) does not touch the A feed",
        );
        let expected = (base.molar[Stream::AcFeed] - FEED_FLOW_EPSILON) * 0.8 + FEED_FLOW_EPSILON;
        let relative = (f.molar[Stream::AcFeed] - expected).abs() / expected;
        assert!(relative < 1e-15, "IDV(7) takes 20% off: {relative:e}");
    }

    /// The mixed feed never reaches exactly zero, even with the valve shut.
    /// That is the whole point of delta D-004.
    #[test]
    fn the_mixed_feed_never_reaches_exactly_zero() {
        let (mut y, unpacked, eq, table) = plant();
        y.valve_pos[3] = 0.0;
        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        assert_exact(f.molar[Stream::AcFeed], FEED_FLOW_EPSILON, "FTM(4)");
        assert!(f.molar[Stream::AcFeed] > 0.0);
    }

    /// The compressor flow floor stops a fully-open recycle valve from
    /// producing a division by zero at `teprob.f:600-601`.
    #[test]
    fn the_compressor_flow_has_a_floor() {
        let (mut y, unpacked, eq, table) = plant();
        let mut eq = eq;
        // Both pressures set explicitly: the fixture's separator pressure is
        // far above the mixing zone's, so raising `PTV` alone clamps the ratio
        // *low* and leaves the bleed at zero. That is how this test failed the
        // first time it was written.
        eq.separator.pressure = 1000.0;
        eq.mixing_pressure = 1.0e6;
        y.valve_pos[4] = 100.0;
        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        assert!(
            f.molar[Stream::Recycle] > 0.0 && f.molar[Stream::Recycle].is_finite(),
            "the floor must keep the recycle flow positive and finite"
        );
        assert_exact(
            f.molar[Stream::Recycle],
            MIN_COMPRESSOR_FLOW / table.molar_mass[Stream::Recycle],
            "the floor is applied to the mass flow, before the division",
        );
    }

    /// The agitator is a speed, not a flow: it runs from 1.5 to 2.5 over the
    /// valve's travel and never reaches zero.
    #[test]
    fn the_agitator_never_stops() {
        let (mut y, unpacked, eq, table) = plant();
        y.valve_pos[11] = 0.0;
        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        assert_exact(f.agitator, 1.5, "AGSP at a shut valve");
        y.valve_pos[11] = 100.0;
        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        assert_exact(f.agitator, 2.5, "AGSP at a fully open valve");
    }

    /// Component flows are the stream's composition times its total flow, and
    /// the three stripper streams are left alone.
    #[test]
    fn component_flows_split_the_total_by_composition() {
        let (_, _, _, table) = plant();
        let f = solved();
        for stream in crate::streams::ASSEMBLED_STREAMS {
            let mut sum = 0.0;
            for c in Component::ALL {
                sum += f.component[stream][c];
                assert_exact(
                    f.component[stream][c],
                    table.composition[stream][c] * f.molar[stream],
                    "FCM",
                );
            }
            let expected = f.molar[stream] * table.composition[stream].sum();
            assert!((sum - expected).abs() <= 1e-9 * expected.abs().max(1.0));
        }
        for stream in [
            Stream::StripperOverhead,
            Stream::ReactorInlet,
            Stream::StripperDownflow,
        ] {
            for c in Component::ALL {
                assert_exact(f.component[stream][c], 0.0, "left for B-0022");
            }
        }
    }

    /// `PR**3` must expand to multiplication. If it ever went through `pow`,
    /// the cube of an exactly-representable ratio would still agree, so the
    /// check is on a value where the two differ.
    #[test]
    fn the_cube_in_the_compressor_curve_is_not_a_pow_call() {
        let ratio = 1.2345678901234567_f64;
        let by_multiplication = ratio * ratio * ratio;
        let by_pow = crate::math::pow(ratio, 3.0);
        assert!(
            by_multiplication.to_bits() != by_pow.to_bits(),
            "multiplication and pow agree at this ratio, so this test cannot \
             distinguish them and needs a different value"
        );
    }

    /// A stream whose molecular weight is zero divides by zero, and the
    /// original does not check either.
    ///
    /// The IEEE answer is what it should be: a non-finite flow that propagates
    /// visibly, rather than a panic or a silent zero that would look like a
    /// shut valve. An adversarial state can drive a vessel empty, so this is
    /// reachable, not hypothetical.
    #[test]
    fn a_zero_molecular_weight_divides_rather_than_panicking() {
        let (y, unpacked, eq, table) = plant();
        let mut table = table;
        table.molar_mass[Stream::Purge] = 0.0;
        table.composition[Stream::Purge] = Composition::new_unchecked([0.0; 8]);

        let f = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        let purge = f.molar[Stream::Purge];
        assert!(
            !purge.is_finite(),
            "a zero molecular weight must give infinity or NaN, not a finite              number that reads as a real flow: got {purge}"
        );
        assert!(
            purge != 0.0,
            "a silent zero would be indistinguishable from a shut valve"
        );
    }
}
