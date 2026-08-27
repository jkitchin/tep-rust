//! The twenty-two continuous measurements, and the shutdown detector.
//!
//! Ported from `teprob.f:679-710`. This is the last of the *pure* right-hand
//! side: everything after line 710 draws noise or ticks an analyser and
//! belongs to [`crate::plant::Plant::sample_measurements`].
//!
//! # What this layer is
//!
//! Nothing here is physics. Every one of the twenty-two lines takes a quantity
//! the model already computed and converts it into the unit an operator's
//! instrument would read: molar flow to standard cubic metres per hour, molar
//! flow to kilograms per hour, absolute pressure to kPa gauge, volume to
//! percent of a span. Getting a conversion wrong changes no state and no
//! derivative; it changes only what the controller sees, which is worse,
//! because the plant then runs correctly and is controlled wrongly.
//!
//! # The conversion factors, and what they mean
//!
//! | Factor | Where | Meaning |
//! |---|---|---|
//! | `0.359` | `679`, `682`-`684`, `688` | standard cubic feet per lbmol |
//! | `35.3145` | same, and `692`, `695`, `704`-`710` | cubic feet per cubic metre |
//! | `0.454` | `680`, `681`, `697` | kilograms per pound |
//! | `760` | `685`, `691`, `694` | mmHg per atmosphere |
//! | `101.325` | same | kPa per atmosphere |
//!
//! So `FTM * 0.359 / 35.3145` is lbmol/h to standard m^3/h, and
//! `(P - 760)/760 * 101.325` is mmHg absolute to kPa gauge.
//!
//! # `XMEAS(20)` is assigned twice
//!
//! ```fortran
//!       XMEAS(20)=CPDH*0.0003927D6
//!       XMEAS(20)=CPDH*0.29307D3
//! ```
//!
//! The first is dead, and the two factors are not equal: 392.7 against 293.07,
//! a third apart. So this is not a harmless duplicate but a superseded
//! conversion, and a port that took the first line would report compressor
//! work 34% high. Delta D-006.
//!
//! # The shutdown detector
//!
//! Eight limits, checked at `teprob.f:702-710`. The original records only that
//! *something* tripped, in a single integer `ISD`; this port reports which,
//! because "the plant tripped" without a reason is nearly useless to a caller
//! and the information is free.
//!
//! All eight comparisons are strict, so a state exactly on a limit does *not*
//! trip. That is the same shape B-0022 found in the stripper and B-0021 in the
//! compressor, and it means the adversarial states placed *on* the limits
//! exercise the not-tripped side; the tripping side needs states past them.
//! B-0016 built both, which is why `tier2_adversarial` reports four of its
//! boundaries as tripping.
//!
//! Note that two of the eight are phrased in terms of `XMEAS` rather than the
//! underlying quantity: `teprob.f:703` tests the *converted* reactor pressure
//! against 3000 kPa gauge, and `706` tests `XMEAS(9)`, which is `TCR`
//! unconverted. Testing `PTR` against an equivalent mmHg threshold instead
//! would be arithmetically different in the last bits.
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `XMEAS(1..22)` | [`Measured::continuous`] | the instrument readings |
//! | `ISD` | [`Measured::shutdown`] | the trip, with its cause |

// Every float expression here reproduces `teprob.f`'s association and
// rounding exactly; see `crate::thermo`.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]

use crate::constants::{VTC, single};
use crate::flows::Flows;
use crate::heat::HeatTransfer;
use crate::plant::N_CONTINUOUS;
use crate::state::State;
use crate::stream::Stream;
use crate::streams::Streams;
use crate::vessels::Unpacked;

/// Standard cubic feet per lbmol (`teprob.f:679`). Single precision.
const SCF_PER_LBMOL: f64 = single(0.359);
/// Cubic feet per cubic metre (`teprob.f:679`). Single precision.
const CUBIC_FEET_PER_CUBIC_METRE: f64 = single(35.3145);
/// Kilograms per pound (`teprob.f:680`). Single precision.
const KG_PER_LB: f64 = single(0.454);
/// mmHg per atmosphere (`teprob.f:685`). Single precision.
const MMHG_PER_ATM: f64 = single(760.0);
/// kPa per atmosphere (`teprob.f:685`). Single precision.
const KPA_PER_ATM: f64 = single(101.325);

/// Reactor level span offset and range (`teprob.f:686`). Single precision.
const REACTOR_LEVEL: (f64, f64) = (single(84.6), single(666.7));
/// Separator level span offset and range (`teprob.f:690`). Single precision.
const SEPARATOR_LEVEL: (f64, f64) = (single(27.5), single(290.0));
/// Stripper level span offset (`teprob.f:693`). The range is [`VTC`].
const STRIPPER_LEVEL_OFFSET: f64 = single(78.25);
/// Percent, the span every level is reported against. Single precision.
const PERCENT: f64 = single(100.0);

/// Stripper steam duty to kilograms per hour (`teprob.f:697`). Double.
const STEAM_DUTY_TO_KG: f64 = 1.04e3;
/// Compressor duty to kilowatts (`teprob.f:699`). Double.
const COMPRESSOR_DUTY_TO_KW: f64 = 0.29307e3;

/// Reactor pressure trip, kPa gauge (`teprob.f:703`). Single precision.
const REACTOR_PRESSURE_LIMIT: f64 = single(3000.0);
/// Reactor temperature trip, degrees Celsius (`teprob.f:706`).
const REACTOR_TEMPERATURE_LIMIT: f64 = single(175.0);
/// Reactor liquid volume trips, cubic metres (`teprob.f:704-705`).
const REACTOR_VOLUME_LIMITS: (f64, f64) = (single(2.0), single(24.0));
/// Separator liquid volume trips, cubic metres (`teprob.f:707-708`).
const SEPARATOR_VOLUME_LIMITS: (f64, f64) = (single(1.0), single(12.0));
/// Stripper liquid volume trips, cubic metres (`teprob.f:709-710`).
const STRIPPER_VOLUME_LIMITS: (f64, f64) = (single(1.0), single(8.0));

/// One of the eight conditions that shuts the plant down.
///
/// The original has no such type: `teprob.f:702-710` sets a single flag and
/// throws the reason away. Naming them costs nothing and is the difference
/// between a scenario report that says "tripped" and one that says which
/// vessel overfilled.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum ShutdownCause {
    /// `teprob.f:703`: reactor pressure above 3000 kPa gauge.
    ReactorPressureHigh = 0,
    /// `teprob.f:704`: reactor liquid above 24 cubic metres.
    ReactorLevelHigh = 1,
    /// `teprob.f:705`: reactor liquid below 2 cubic metres.
    ReactorLevelLow = 2,
    /// `teprob.f:706`: reactor temperature above 175 C.
    ReactorTemperatureHigh = 3,
    /// `teprob.f:707`: separator liquid above 12 cubic metres.
    SeparatorLevelHigh = 4,
    /// `teprob.f:708`: separator liquid below 1 cubic metre.
    SeparatorLevelLow = 5,
    /// `teprob.f:709`: stripper liquid above 8 cubic metres.
    StripperLevelHigh = 6,
    /// `teprob.f:710`: stripper liquid below 1 cubic metre.
    StripperLevelLow = 7,
}

impl ShutdownCause {
    /// All eight, in the order `teprob.f:703-710` tests them.
    pub const ALL: [ShutdownCause; 8] = [
        ShutdownCause::ReactorPressureHigh,
        ShutdownCause::ReactorLevelHigh,
        ShutdownCause::ReactorLevelLow,
        ShutdownCause::ReactorTemperatureHigh,
        ShutdownCause::SeparatorLevelHigh,
        ShutdownCause::SeparatorLevelLow,
        ShutdownCause::StripperLevelHigh,
        ShutdownCause::StripperLevelLow,
    ];

    /// A short human-readable description.
    #[must_use]
    pub const fn describe(self) -> &'static str {
        match self {
            Self::ReactorPressureHigh => "reactor pressure high",
            Self::ReactorLevelHigh => "reactor level high",
            Self::ReactorLevelLow => "reactor level low",
            Self::ReactorTemperatureHigh => "reactor temperature high",
            Self::SeparatorLevelHigh => "separator level high",
            Self::SeparatorLevelLow => "separator level low",
            Self::StripperLevelHigh => "stripper level high",
            Self::StripperLevelLow => "stripper level low",
        }
    }
}

/// Which of the eight shutdown conditions hold.
///
/// More than one can hold at once, and the original would report that as the
/// same single flag it reports one with, so the set is carried rather than
/// just the first.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Shutdown(u8);

impl Shutdown {
    /// Whether any condition holds. This is the original's `ISD != 0`.
    #[must_use]
    pub const fn is_tripped(self) -> bool {
        self.0 != 0
    }

    /// Whether a particular condition holds.
    #[must_use]
    pub const fn holds(self, cause: ShutdownCause) -> bool {
        self.0 & (1 << (cause as u8)) != 0
    }

    /// The first condition that holds, in the order the original tests them.
    #[must_use]
    pub fn first(self) -> Option<ShutdownCause> {
        ShutdownCause::ALL.into_iter().find(|c| self.holds(*c))
    }

    /// Every condition that holds.
    pub fn causes(self) -> impl Iterator<Item = ShutdownCause> {
        ShutdownCause::ALL
            .into_iter()
            .filter(move |c| self.holds(*c))
    }

    /// Record a condition.
    fn set(&mut self, cause: ShutdownCause, tripped: bool) {
        if tripped {
            self.0 |= 1 << (cause as u8);
        }
    }
}

/// The three vessel pressures this layer converts, in mmHg absolute:
/// `(PTR, PTS, PTV)`.
///
/// A tuple rather than a struct because it is passed straight through from
/// [`crate::equilibrium::Equilibrium`] and naming three fields twice would add
/// nothing.
pub type Pressures = (f64, f64, f64);

/// Everything `teprob.f:679-710` produces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Measured {
    /// `XMEAS(1..22)`, before any noise is added.
    pub continuous: [f64; N_CONTINUOUS],
    /// The trip, with its causes.
    pub shutdown: Shutdown,
}

/// Assemble the noise-free measurements and test the eight shutdown limits.
// @port  teprob.f:679-710
// @delta D-006 class=A teprob.f:698
#[must_use]
pub fn measurements(
    y: &State,
    unpacked: &Unpacked,
    stream_table: &Streams,
    flow: &Flows,
    heat: &HeatTransfer,
    pressures: Pressures,
) -> Measured {
    let (reactor_pressure, separator_pressure, mixing_pressure) = pressures;
    // lbmol/h to standard cubic metres per hour.
    let gas_flow = |stream: Stream| flow.molar[stream] * SCF_PER_LBMOL / CUBIC_FEET_PER_CUBIC_METRE;
    // lbmol/h to kilograms per hour, through the stream's molecular weight.
    let mass_flow =
        |stream: Stream| flow.molar[stream] * stream_table.molar_mass[stream] * KG_PER_LB;
    // mmHg absolute to kPa gauge.
    let gauge = |p: f64| (p - MMHG_PER_ATM) / MMHG_PER_ATM * KPA_PER_ATM;
    // Volume to percent of a span.
    let level = |v: f64, span: (f64, f64)| (v - span.0) / span.1 * PERCENT;
    // lbmol/h to cubic metres per hour, through a liquid density.
    let liquid_flow =
        |stream: Stream, density: f64| flow.molar[stream] / density / CUBIC_FEET_PER_CUBIC_METRE;

    let mut continuous = [0.0; N_CONTINUOUS];
    continuous[0] = gas_flow(Stream::AFeed); // teprob.f:679
    continuous[1] = mass_flow(Stream::DFeed); // 680
    continuous[2] = mass_flow(Stream::EFeed); // 681
    continuous[3] = gas_flow(Stream::AcFeed); // 682
    continuous[4] = gas_flow(Stream::Recycle); // 683
    continuous[5] = gas_flow(Stream::MixingZoneOutlet); // 684
    continuous[6] = gauge(reactor_pressure); // 685
    continuous[7] = level(unpacked.reactor.volume, REACTOR_LEVEL); // 686
    continuous[8] = unpacked.reactor.celsius; // 687
    continuous[9] = gas_flow(Stream::Purge); // 688
    continuous[10] = unpacked.separator.celsius; // 689
    continuous[11] = level(unpacked.separator.volume, SEPARATOR_LEVEL); // 690
    continuous[12] = gauge(separator_pressure); // 691
    continuous[13] = liquid_flow(Stream::SeparatorUnderflow, unpacked.separator.density); // 692
    // teprob.f:693. The stripper's span is the vessel volume itself, unlike
    // the other two, which carry their own hard-coded ranges.
    continuous[14] = (unpacked.stripper.volume - STRIPPER_LEVEL_OFFSET) / VTC * PERCENT;
    continuous[15] = gauge(mixing_pressure); // 694
    continuous[16] = liquid_flow(Stream::Product, unpacked.stripper.density); // 695
    continuous[17] = unpacked.stripper.celsius; // 696
    continuous[18] = heat.stripper_duty * STEAM_DUTY_TO_KG * KG_PER_LB; // 697
    // teprob.f:698-699. The first assignment is dead; see delta D-006.
    continuous[19] = flow.compressor_work * COMPRESSOR_DUTY_TO_KW;
    continuous[20] = y.reactor_cw_out_c; // 700
    continuous[21] = y.condenser_cw_out_c; // 701

    // teprob.f:702-710. Every comparison is strict, so a state exactly on a
    // limit does not trip.
    let cubic_metres = |v: f64| v / CUBIC_FEET_PER_CUBIC_METRE;
    let reactor = cubic_metres(unpacked.reactor.volume);
    let separator = cubic_metres(unpacked.separator.volume);
    let stripper = cubic_metres(unpacked.stripper.volume);

    let mut shutdown = Shutdown::default();
    // Tested against the *converted* measurement, not against `PTR`; an
    // equivalent threshold in mmHg would round differently.
    shutdown.set(
        ShutdownCause::ReactorPressureHigh,
        continuous[6] > REACTOR_PRESSURE_LIMIT,
    );
    shutdown.set(
        ShutdownCause::ReactorLevelHigh,
        reactor > REACTOR_VOLUME_LIMITS.1,
    );
    shutdown.set(
        ShutdownCause::ReactorLevelLow,
        reactor < REACTOR_VOLUME_LIMITS.0,
    );
    shutdown.set(
        ShutdownCause::ReactorTemperatureHigh,
        continuous[8] > REACTOR_TEMPERATURE_LIMIT,
    );
    shutdown.set(
        ShutdownCause::SeparatorLevelHigh,
        separator > SEPARATOR_VOLUME_LIMITS.1,
    );
    shutdown.set(
        ShutdownCause::SeparatorLevelLow,
        separator < SEPARATOR_VOLUME_LIMITS.0,
    );
    shutdown.set(
        ShutdownCause::StripperLevelHigh,
        stripper > STRIPPER_VOLUME_LIMITS.1,
    );
    shutdown.set(
        ShutdownCause::StripperLevelLow,
        stripper < STRIPPER_VOLUME_LIMITS.0,
    );

    Measured {
        continuous,
        shutdown,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::component::Component;
    use crate::equilibrium::equilibrium;
    use crate::flows::{FlowDrift, flows};
    use crate::heat::{HeatDrift, heat_transfer};
    use crate::streams::{FeedConditions, streams};
    use crate::stripper::stripper;
    use crate::testing::assert_exact;
    use crate::variables::MeasIndex;
    use crate::vessels::{TemperatureSeeds, unpack};

    fn plant() -> (
        State,
        Unpacked,
        Streams,
        Flows,
        HeatTransfer,
        (f64, f64, f64),
    ) {
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
        let pressures = (
            eq.reactor.pressure,
            eq.separator.pressure,
            eq.mixing_pressure,
        );
        (y, unpacked, table, flow, heat, pressures)
    }

    fn solved() -> Measured {
        let (y, unpacked, table, flow, heat, pressures) = plant();
        measurements(&y, &unpacked, &table, &flow, &heat, pressures)
    }

    /// The three temperatures are reported raw, with no conversion at all.
    /// Anything else would be an invented unit change.
    #[test]
    fn the_temperatures_are_reported_unconverted() {
        let (_, unpacked, _, _, _, _) = plant();
        let m = solved();
        assert_exact(m.continuous[8], unpacked.reactor.celsius, "XMEAS(9)");
        assert_exact(m.continuous[10], unpacked.separator.celsius, "XMEAS(11)");
        assert_exact(m.continuous[17], unpacked.stripper.celsius, "XMEAS(18)");
    }

    /// The three gauge pressures read zero at one atmosphere, which is what
    /// makes them gauge rather than absolute.
    #[test]
    fn gauge_pressures_read_zero_at_atmospheric() {
        let (y, unpacked, table, flow, heat, _) = plant();
        let m = measurements(
            &y,
            &unpacked,
            &table,
            &flow,
            &heat,
            (MMHG_PER_ATM, MMHG_PER_ATM, MMHG_PER_ATM),
        );
        for slot in [6, 12, 15] {
            assert_exact(m.continuous[slot], 0.0, "gauge at one atmosphere");
        }
    }

    /// `XMEAS(20)` uses the *second* factor. The first line is dead and the
    /// two differ by a third, so taking it would be a large silent error.
    #[test]
    fn the_compressor_work_uses_the_second_conversion_factor() {
        let (_, _, _, flow, _, _) = plant();
        let m = solved();
        assert_exact(
            m.continuous[19],
            flow.compressor_work * COMPRESSOR_DUTY_TO_KW,
            "XMEAS(20)",
        );
        // The dead first assignment, for contrast.
        let superseded = flow.compressor_work * 0.0003927e6;
        assert!(
            (superseded - m.continuous[19]).abs() > 0.3 * m.continuous[19].abs(),
            "the two factors are close enough that the dead line would not \
             show up, so delta D-006 is not worth its entry"
        );
    }

    /// Every shutdown comparison is strict, so a state exactly on a limit does
    /// not trip. Same shape as B-0021's `CPPRMX` and B-0022's `TCC`.
    #[test]
    fn a_state_exactly_on_a_limit_does_not_trip() {
        let (y, unpacked, table, flow, heat, _) = plant();
        let mut at_limit = unpacked;
        // The fixture runs well over 3000 kPa gauge, so the pressure trip
        // fires on its own and would mask the one under test.
        let pressures = (20_000.0, 20_000.0, 20_000.0);
        at_limit.reactor.celsius = 120.0;
        // Exactly 24 cubic metres of reactor liquid.
        at_limit.reactor.volume = REACTOR_VOLUME_LIMITS.1 * CUBIC_FEET_PER_CUBIC_METRE;
        let m = measurements(&y, &at_limit, &table, &flow, &heat, pressures);
        assert!(
            !m.shutdown.holds(ShutdownCause::ReactorLevelHigh),
            "a level of exactly 24 m^3 tripped, so teprob.f:704 was read as \
             `>=` rather than the `.GT.` it is"
        );

        // And a hair above it does.
        at_limit.reactor.volume *= 1.000_001;
        let m = measurements(&y, &at_limit, &table, &flow, &heat, pressures);
        assert!(m.shutdown.holds(ShutdownCause::ReactorLevelHigh));
        assert!(m.shutdown.is_tripped());
        assert_eq!(m.shutdown.first(), Some(ShutdownCause::ReactorLevelHigh));
    }

    /// Each of the eight limits can be reached independently, and each names
    /// itself. A transposed pair would be invisible in a bare flag.
    #[test]
    fn each_of_the_eight_limits_reports_its_own_cause() {
        let (y, unpacked, table, flow, heat, _) = plant();
        let ft3 = CUBIC_FEET_PER_CUBIC_METRE;

        /// One condition, and the perturbation that provokes it.
        type Case<'a> = (ShutdownCause, &'a dyn Fn(&mut Unpacked, &mut Pressures));

        let cases: [Case<'_>; 8] = [
            (ShutdownCause::ReactorPressureHigh, &|_, p| {
                p.0 = 3.0e5;
            }),
            (ShutdownCause::ReactorLevelHigh, &|u, _| {
                u.reactor.volume = 30.0 * ft3;
            }),
            (ShutdownCause::ReactorLevelLow, &|u, _| {
                u.reactor.volume = 1.0 * ft3;
            }),
            (ShutdownCause::ReactorTemperatureHigh, &|u, _| {
                u.reactor.celsius = 200.0;
            }),
            (ShutdownCause::SeparatorLevelHigh, &|u, _| {
                u.separator.volume = 20.0 * ft3;
            }),
            (ShutdownCause::SeparatorLevelLow, &|u, _| {
                u.separator.volume = 0.5 * ft3;
            }),
            (ShutdownCause::StripperLevelHigh, &|u, _| {
                u.stripper.volume = 10.0 * ft3;
            }),
            (ShutdownCause::StripperLevelLow, &|u, _| {
                u.stripper.volume = 0.5 * ft3;
            }),
        ];

        for (expected, apply) in cases {
            let mut u = unpacked;
            // Put every level and pressure comfortably inside its band first,
            // so only the condition under test fires.
            let mut p: Pressures = (20_000.0, 20_000.0, 20_000.0);
            u.reactor.volume = 10.0 * ft3;
            u.separator.volume = 5.0 * ft3;
            u.stripper.volume = 4.0 * ft3;
            u.reactor.celsius = 120.0;
            apply(&mut u, &mut p);

            let m = measurements(&y, &u, &table, &flow, &heat, p);
            let fired: alloc::vec::Vec<_> = m.shutdown.causes().collect();
            assert_eq!(
                fired,
                alloc::vec![expected],
                "expected only {expected:?} to fire, got {fired:?}"
            );
        }
    }

    /// A plant well inside every limit does not trip.
    #[test]
    fn a_healthy_plant_does_not_trip() {
        let (y, unpacked, table, flow, heat, _) = plant();
        let ft3 = CUBIC_FEET_PER_CUBIC_METRE;
        let mut u = unpacked;
        u.reactor.volume = 10.0 * ft3;
        u.separator.volume = 5.0 * ft3;
        u.stripper.volume = 4.0 * ft3;
        u.reactor.celsius = 120.0;
        let m = measurements(&y, &u, &table, &flow, &heat, (20_000.0, 20_000.0, 20_000.0));
        assert!(!m.shutdown.is_tripped());
        assert_eq!(m.shutdown.first(), None);
        assert_eq!(m.shutdown.causes().count(), 0);
    }

    /// Several conditions can hold at once, and all of them are reported. The
    /// original would report exactly the same single flag for one as for four.
    #[test]
    fn more_than_one_cause_can_hold_at_once() {
        let (y, unpacked, table, flow, heat, _) = plant();
        let mut u = unpacked;
        u.reactor.volume = 0.0;
        u.separator.volume = 0.0;
        u.stripper.volume = 0.0;
        u.reactor.celsius = 300.0;
        let m = measurements(&y, &u, &table, &flow, &heat, (1.0e6, 1.0e6, 1.0e6));
        let fired: alloc::vec::Vec<_> = m.shutdown.causes().collect();
        assert!(
            fired.len() >= 4,
            "only {fired:?} fired on a plant that is empty, overheated and \
             overpressured"
        );
        assert_eq!(m.shutdown.first(), Some(ShutdownCause::ReactorPressureHigh));
    }

    /// The measurement slots line up with `MeasIndex`, so a caller indexing by
    /// name reaches the value this module put there.
    #[test]
    fn the_slots_line_up_with_the_measurement_indices() {
        let m = solved();
        assert_eq!(m.continuous.len(), N_CONTINUOUS);
        // XMEAS(9) is the reactor temperature; index 9 one-based is slot 8.
        let index = MeasIndex::new(9).expect("in range");
        assert_eq!(index.zero_based(), 8);
        assert!(m.continuous[index.zero_based()] > 0.0);
    }
}
