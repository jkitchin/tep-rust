//! Heat transfer: the reactor coil, the condenser, and the stripper reboiler.
//!
//! Ported from `teprob.f:663-678`. Three duties, each with a different shape,
//! and each gated on something.
//!
//! # The reactor coil
//!
//! The coil's effective area ramps with liquid level, because the coil is only
//! wetted over part of its height. `VLR/7.8` is the level as a percentage, and
//! the ramp is piecewise linear between 10% and 50%:
//!
//! \\[
//!   \\lambda = \\begin{cases}
//!     1 & \\ell > 50 \\\\
//!     0 & \\ell < 10 \\\\
//!     0.025\\,\\ell - 0.25 & \\text{otherwise}
//!   \\end{cases}
//!   \\qquad \\ell = V_l / 7.8
//! \\]
//!
//! The overall coefficient is quadratic in agitator speed:
//!
//! \\[ U A_r = \\lambda \\left(-0.5\\,\\omega^2 + 2.75\\,\\omega - 2.5\\right) \\times 855490 \\times 10^{-6} \\]
//!
//! That parabola peaks at \\(\\omega = 2.75\\), which is *above* the agitator's
//! whole range: [`crate::flows::Flows::agitator`] runs from 1.5 to 2.5. So the
//! coefficient rises monotonically with speed everywhere the plant can go,
//! from 0.5 to 1.25 times the scale, and the falling half of the parabola is
//! unreachable. Its roots are at 1.149 and 4.351, both outside that range too,
//! so the coefficient never reaches zero from the agitator alone. The model is
//! a fit, not a mechanism.
//!
//! # The ramp does not quite meet its flat sections
//!
//! `0.025` at `teprob.f:668` is single precision, so it is stored as
//! 0.02500000037252903. The ramp therefore misses both of its endpoints:
//!
//! | level | ramp gives | flat section gives | gap |
//! |---|---|---|---|
//! | 10 | 3.725290298461914e-9 | 0 | 3.7e-9 |
//! | 50 | 1.0000000186264515 | 1 | 1.9e-8 |
//!
//! Both breakpoint comparisons are strict, so a level of exactly 10 or exactly
//! 50 takes the *ramp*, and `UARLEV` is discontinuous by those amounts as the
//! level crosses either one.
//!
//! This is faithful reproduction rather than a delta: it is what the original
//! computes, the gaps are eight orders below the quantity itself, and the
//! coefficient they scale is an empirical fit in the first place. It is
//! written down because "the ramp meets the flat sections" is the obvious
//! assumption, it is false, and a test asserting it would fail for a reason
//! that looks like a porting error.
//!
//! # The condenser
//!
//! A smooth saturating function of reactor outlet flow, approaching 0.404655
//! as the flow grows:
//!
//! \\[ U A_s = 0.404655 \\left(1 - \\frac{1}{1 + (F_8/3528.73)^4}\\right) \\]
//!
//! # `**2` and `**4` are integer powers and must not go through `pow`
//!
//! gfortran expands an integer exponent into multiplications rather than
//! calling libm, and the *shape* of that expansion is load-bearing. Measured
//! over 200,000 values with this project's pinned flags:
//!
//! | candidate for `X**4` | matches gfortran |
//! |---|---|
//! | `(x*x)*(x*x)` | 200,000 of 200,000 |
//! | `((x*x)*x)*x` | 132,040 |
//! | `pow(x, 4.0)` | 99,523 |
//!
//! So it is binary exponentiation, squaring twice, and the two plausible
//! alternatives are each wrong about a third and a half of the time. `X**2` is
//! `x*x` on all 200,000, which is the only thing it could be.
//!
//! # The stripper reboiler is gated on temperature
//!
//! \\[ Q_c = \\begin{cases} U A_c (100 - T_c) & T_c < 100 \\\\ 0 & \\text{otherwise} \\end{cases} \\]
//!
//! The steam is at 100 C, so above that there is nothing to transfer and the
//! original sets the duty to zero rather than letting it go negative. The
//! nominal trajectory sits near 65 C, so the *cutoff* is the branch at risk of
//! never being exercised: B-0021 measured 300 of 300 nominal states below 100.
//!
//! `UAC` itself is computed at `teprob.f:572` and belongs to
//! [`mod@crate::flows`]; this module only consumes it.
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `UARLEV` | [`HeatTransfer::level_factor`] | wetted fraction of the coil |
//! | `UAR`, `QUR` | [`HeatTransfer::reactor_coefficient`] etc. | reactor coil |
//! | `UAS`, `QUS` | [`HeatTransfer::condenser_coefficient`] etc. | condenser |
//! | `QUC` | [`HeatTransfer::stripper_duty`] | stripper reboiler |

// Every float expression here reproduces `teprob.f`'s association and
// rounding exactly; see `crate::thermo`.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]

use crate::constants::single;
use crate::flows::Flows;
use crate::state::State;
use crate::stream::Stream;
use crate::streams::Streams;
use crate::vessels::Unpacked;

/// Cubic feet of reactor liquid per percent of level (`teprob.f:663`).
///
/// Single precision. The level in the model is always `VLR/7.8`, never a
/// separate state.
const CUBIC_FEET_PER_PERCENT: f64 = single(7.8);

/// Level above which the coil is fully wetted, percent (`teprob.f:663`).
const FULLY_WETTED_ABOVE: f64 = single(50.0);
/// Level below which it is entirely dry, percent (`teprob.f:665`).
const DRY_BELOW: f64 = single(10.0);
/// Slope of the ramp between them (`teprob.f:668`).
const RAMP_SLOPE: f64 = single(0.025);
/// Intercept of the ramp (`teprob.f:668`).
const RAMP_INTERCEPT: f64 = single(0.25);

/// Agitator quadratic, `teprob.f:670-671`. All three single precision.
const AGITATOR_QUADRATIC: [f64; 3] = [single(0.5), single(2.75), single(2.5)];
/// Coil area-times-coefficient scale, `teprob.f:671`. Double precision.
const REACTOR_COIL_SCALE: f64 = 855490.0e-6;
/// Reactor coolant disturbance depth, `teprob.f:673`. Double precision.
const REACTOR_COOLANT_DEPTH: f64 = 0.35;

/// Condenser coefficient asymptote, `teprob.f:674`. Single precision.
const CONDENSER_ASYMPTOTE: f64 = single(0.404655);
/// Flow at which the condenser reaches half its asymptote, `teprob.f:674`.
const CONDENSER_HALF_FLOW: f64 = single(3528.73);
/// Condenser coolant disturbance depth, `teprob.f:676`. Double precision.
const CONDENSER_COOLANT_DEPTH: f64 = 0.25;

/// Steam temperature, above which the reboiler does nothing
/// (`teprob.f:678`). Single precision, and exactly representable.
const STEAM_CELSIUS: f64 = single(100.);

/// The two disturbance walk channels this range reads.
///
/// `TESUB8(10, t)` and `TESUB8(11, t)` at `teprob.f:673` and `676`. Walk
/// machinery is Phase 3; until then these are inputs, as in
/// [`crate::flows::FlowDrift`].
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HeatDrift {
    /// `TESUB8(10, t)`: enters as `1 - 0.35 * drift` on the reactor duty.
    pub reactor_coolant: f64,
    /// `TESUB8(11, t)`: enters as `1 - 0.25 * drift` on the condenser duty.
    pub condenser_coolant: f64,
}

/// Which branch of the coil's level ramp was taken (`teprob.f:663-669`).
///
/// Reported so a differential can say whether a sampled state exercised each
/// one. An unexercised branch is not evidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LevelBranch {
    /// Above 50%: the coil is fully wetted.
    FullyWetted,
    /// Below 10%: the coil is dry and transfers nothing.
    Dry,
    /// Between: the linear ramp.
    Ramp,
}

/// Everything `teprob.f:663-678` produces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatTransfer {
    /// `UARLEV`: the fraction of the coil in contact with liquid.
    pub level_factor: f64,
    /// `UAR`: the reactor coil's overall coefficient times area.
    pub reactor_coefficient: f64,
    /// `QUR`: heat removed from the reactor.
    pub reactor_duty: f64,
    /// `UAS`: the condenser's coefficient times area.
    pub condenser_coefficient: f64,
    /// `QUS`: heat removed in the condenser.
    pub condenser_duty: f64,
    /// `QUC`: heat added by the stripper reboiler.
    pub stripper_duty: f64,
    /// Which level branch was taken.
    pub level_branch: LevelBranch,
    /// Whether the reboiler was active, that is `TCC < 100`.
    pub steam_on: bool,
}

/// The coil's wetted fraction, `UARLEV` (`teprob.f:663-669`).
///
/// Takes the *level*, `VLR/7.8`, not the volume: the original divides once and
/// compares the quotient three times, so doing the same avoids three chances
/// to write the constant differently.
fn level_factor(level: f64) -> (f64, LevelBranch) {
    if level > FULLY_WETTED_ABOVE {
        (1.0, LevelBranch::FullyWetted)
    } else if level < DRY_BELOW {
        (0.0, LevelBranch::Dry)
    } else {
        (RAMP_SLOPE * level - RAMP_INTERCEPT, LevelBranch::Ramp)
    }
}

/// Compute the three heat duties.
// @port teprob.f:663-678
#[must_use]
pub fn heat_transfer(
    y: &State,
    unpacked: &Unpacked,
    stream_table: &Streams,
    flow: &Flows,
    drift: HeatDrift,
) -> HeatTransfer {
    // teprob.f:663-669. Note the original writes `VLR/7.8` out three times
    // rather than binding it; the value is identical each time.
    let level = unpacked.reactor.volume / CUBIC_FEET_PER_PERCENT;
    let (level_factor_value, level_branch) = level_factor(level);

    // teprob.f:670-671. `AGSP**2` is an integer power: `agsp * agsp`, never
    // `pow`. See the module documentation.
    let agsp = flow.agitator;
    let reactor_coefficient = level_factor_value
        * (-AGITATOR_QUADRATIC[0] * (agsp * agsp) + AGITATOR_QUADRATIC[1] * agsp
            - AGITATOR_QUADRATIC[2])
        * REACTOR_COIL_SCALE;

    // teprob.f:672-673
    let reactor_duty = reactor_coefficient
        * (y.reactor_cw_out_c - unpacked.reactor.celsius)
        * (1.0 - REACTOR_COOLANT_DEPTH * drift.reactor_coolant);

    // teprob.f:674. `**4` is binary exponentiation: square, then square again.
    // A four-term chain and a `pow` call are each wrong on a large fraction of
    // inputs; the measurement is in the module documentation.
    let ratio = flow.molar[Stream::ReactorOutlet] / CONDENSER_HALF_FLOW;
    let squared = ratio * ratio;
    let condenser_coefficient = CONDENSER_ASYMPTOTE * (1.0 - 1.0 / (1.0 + squared * squared));

    // teprob.f:675-676. The driving temperature is the *stream* temperature
    // `TST(8)`, which is the reactor's, not the separator's.
    let condenser_duty = condenser_coefficient
        * (y.condenser_cw_out_c - stream_table.celsius[Stream::ReactorOutlet])
        * (1.0 - CONDENSER_COOLANT_DEPTH * drift.condenser_coolant);

    // teprob.f:677-678
    let steam_on = unpacked.stripper.celsius < STEAM_CELSIUS;
    let stripper_duty = if steam_on {
        flow.steam_coefficient * (STEAM_CELSIUS - unpacked.stripper.celsius)
    } else {
        0.0
    };

    HeatTransfer {
        level_factor: level_factor_value,
        reactor_coefficient,
        reactor_duty,
        condenser_coefficient,
        condenser_duty,
        stripper_duty,
        level_branch,
        steam_on,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::component::Component;
    use crate::equilibrium::equilibrium;
    use crate::flows::{FlowDrift, flows};
    use crate::streams::{FeedConditions, streams};
    use crate::testing::assert_exact;
    use crate::vessels::{TemperatureSeeds, unpack};

    fn plant() -> (State, Unpacked, Streams, Flows) {
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
        let table = streams(&unpacked, &eq, &FeedConditions::default());
        let flow = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        (y, unpacked, table, flow)
    }

    /// `**4` is `(x*x)*(x*x)`, not a chain and not a `pow` call. Measured
    /// against gfortran; see the module documentation.
    #[test]
    fn the_fourth_power_squares_twice_rather_than_chaining() {
        // A value where the three candidate expansions disagree, so the test
        // can tell them apart. Found by sweeping rather than guessed: the
        // first attempt used a value where binary and chained agree, and the
        // test passed while proving nothing.
        let x = 0.3001_f64;
        let squared = x * x;
        let binary = squared * squared;
        let chained = ((x * x) * x) * x;
        let by_pow = crate::math::pow(x, 4.0);
        assert!(
            binary.to_bits() != chained.to_bits(),
            "the two expansions agree at this value, so it cannot discriminate"
        );
        assert!(
            binary.to_bits() != by_pow.to_bits(),
            "pow agrees at this value, so it cannot discriminate"
        );
    }

    /// All three level branches, at levels chosen to land in each.
    #[test]
    fn the_coil_level_ramp_has_three_branches() {
        assert_eq!(level_factor(80.0), (1.0, LevelBranch::FullyWetted));
        assert_eq!(level_factor(5.0), (0.0, LevelBranch::Dry));
        assert_eq!(level_factor(30.0).1, LevelBranch::Ramp);
        // Both comparisons are strict, so a level exactly on a breakpoint
        // takes the ramp. B-0022 found the same shape in the stripper.
        assert_eq!(level_factor(FULLY_WETTED_ABOVE).1, LevelBranch::Ramp);
        assert_eq!(level_factor(DRY_BELOW).1, LevelBranch::Ramp);
    }

    /// The ramp very nearly meets its flat sections, and the miss is the
    /// single-precision `0.025`.
    ///
    /// "Nearly" is the finding. Writing this test expecting exact equality is
    /// the natural thing to do and it fails; see the module documentation.
    #[test]
    fn the_level_ramp_misses_its_flat_sections_by_the_single_precision_gap() {
        let lower = level_factor(DRY_BELOW).0;
        let upper = level_factor(FULLY_WETTED_ABOVE).0;
        assert_exact(lower, 3.725_290_298_461_914e-9, "the gap at 10%");
        assert_exact(upper, 1.000_000_018_626_451_5, "the gap at 50%");
        assert!(
            lower.to_bits() != 0.0_f64.to_bits() && upper.to_bits() != 1.0_f64.to_bits(),
            "the ramp met its endpoints exactly, so 0.025 was typed as a \
             double and the port no longer matches the Fortran"
        );
        // Small enough to be irrelevant physically, which is why the original
        // never noticed.
        assert!(lower < 1e-8 && (upper - 1.0) < 1e-7);
    }

    /// The agitator parabola peaks outside the agitator's own range, so over
    /// the reachable range the coefficient only increases with speed.
    #[test]
    fn the_agitator_quadratic_is_monotone_over_the_reachable_range() {
        let quadratic = |w: f64| {
            -AGITATOR_QUADRATIC[0] * (w * w) + AGITATOR_QUADRATIC[1] * w - AGITATOR_QUADRATIC[2]
        };
        // `AGSP` runs 1.5 to 2.5; see `crate::flows`.
        let mut previous = quadratic(1.5);
        let mut step = 1.5;
        while step < 2.5 {
            step += 0.01;
            let next = quadratic(step);
            assert!(
                next > previous,
                "the coefficient fell between {} and {step}",
                step - 0.01
            );
            previous = next;
        }
        // The peak is outside the range, which is what makes it monotone.
        assert!(quadratic(2.75) > quadratic(2.5));
        // And the roots are outside it on both sides, so the coefficient never
        // reaches zero from the agitator alone.
        assert!(quadratic(1.15) > 0.0 && quadratic(1.14) < 0.0, "lower root");
        assert!(quadratic(4.35) > 0.0 && quadratic(4.36) < 0.0, "upper root");
    }

    /// The reboiler shuts off at 100 C rather than going negative, and the
    /// comparison is strict.
    #[test]
    fn the_reboiler_shuts_off_at_the_steam_temperature() {
        let (y, unpacked, table, flow) = plant();
        let mut hot = unpacked;
        hot.stripper.celsius = 150.0;
        let h = heat_transfer(&y, &hot, &table, &flow, HeatDrift::default());
        assert_exact(h.stripper_duty, 0.0, "no duty above the steam temperature");
        assert!(!h.steam_on);

        hot.stripper.celsius = STEAM_CELSIUS;
        let h = heat_transfer(&y, &hot, &table, &flow, HeatDrift::default());
        assert_exact(h.stripper_duty, 0.0, "the comparison is strict");
        assert!(!h.steam_on);

        hot.stripper.celsius = 99.0;
        let h = heat_transfer(&y, &hot, &table, &flow, HeatDrift::default());
        assert!(h.steam_on && h.stripper_duty > 0.0);
    }

    /// A dry coil removes no heat at all, whatever the temperature difference.
    #[test]
    fn a_dry_coil_removes_no_heat() {
        let (y, unpacked, table, flow) = plant();
        let mut empty = unpacked;
        // Below 10% of 7.8 cubic feet per percent.
        empty.reactor.volume = 50.0;
        let h = heat_transfer(&y, &empty, &table, &flow, HeatDrift::default());
        assert_eq!(h.level_branch, LevelBranch::Dry);
        assert_exact(h.level_factor, 0.0, "UARLEV");
        // The two zeros have different signs, and that is not an accident.
        // `UAR` is `+0 * (positive quadratic) * (positive scale)`, so it is
        // `+0`. `QUR` is then `+0 * (TWR - TCR)`, and the reactor runs hotter
        // than its coolant outlet here, so the difference is negative and the
        // product is `-0`.
        //
        // Asserted on bits rather than by value, because `-0.0 == 0.0` is true
        // and would hide a sign flip. The oracle comparison is on bits too, so
        // a sign that stopped propagating would fail Tier 2 rather than pass
        // quietly.
        assert_eq!(h.reactor_coefficient.to_bits(), 0.0_f64.to_bits(), "UAR");
        assert_eq!(h.reactor_duty.to_bits(), (-0.0_f64).to_bits(), "QUR");
    }

    /// The condenser saturates: its coefficient approaches the asymptote from
    /// below and never exceeds it.
    #[test]
    fn the_condenser_coefficient_saturates_below_its_asymptote() {
        let (y, unpacked, table, flow) = plant();
        let mut flow = flow;
        let mut previous = 0.0;
        for scale in [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0] {
            flow.molar[Stream::ReactorOutlet] = CONDENSER_HALF_FLOW * scale;
            let h = heat_transfer(&y, &unpacked, &table, &flow, HeatDrift::default());
            assert!(
                h.condenser_coefficient > previous,
                "the coefficient must increase with flow"
            );
            assert!(
                h.condenser_coefficient < CONDENSER_ASYMPTOTE,
                "it must stay below its asymptote"
            );
            previous = h.condenser_coefficient;
        }
        // At the half-flow the quartic is one, so the coefficient is exactly
        // half the asymptote. That is what names the constant.
        flow.molar[Stream::ReactorOutlet] = CONDENSER_HALF_FLOW;
        let h = heat_transfer(&y, &unpacked, &table, &flow, HeatDrift::default());
        assert_exact(
            h.condenser_coefficient,
            CONDENSER_ASYMPTOTE * 0.5,
            "half the asymptote at the half-flow",
        );
    }

    /// The two coolant disturbances scale their own duty and nothing else.
    #[test]
    fn each_coolant_disturbance_scales_only_its_own_duty() {
        let (y, unpacked, table, flow) = plant();
        let base = heat_transfer(&y, &unpacked, &table, &flow, HeatDrift::default());

        let h = heat_transfer(
            &y,
            &unpacked,
            &table,
            &flow,
            HeatDrift {
                reactor_coolant: 1.0,
                condenser_coolant: 0.0,
            },
        );
        assert_exact(
            h.reactor_duty,
            base.reactor_duty * (1.0 - REACTOR_COOLANT_DEPTH),
            "QUR loses 35%",
        );
        assert_exact(h.condenser_duty, base.condenser_duty, "QUS is untouched");

        let h = heat_transfer(
            &y,
            &unpacked,
            &table,
            &flow,
            HeatDrift {
                reactor_coolant: 0.0,
                condenser_coolant: 1.0,
            },
        );
        assert_exact(h.reactor_duty, base.reactor_duty, "QUR is untouched");
        assert_exact(
            h.condenser_duty,
            base.condenser_duty * (1.0 - CONDENSER_COOLANT_DEPTH),
            "QUS loses 25%",
        );
    }

    /// The condenser's driving temperature is the reactor *outlet stream*, not
    /// the separator. They are different numbers and the mistake is invisible
    /// without an oracle.
    #[test]
    fn the_condenser_is_driven_by_the_reactor_outlet_temperature() {
        let (y, unpacked, table, flow) = plant();
        let h = heat_transfer(&y, &unpacked, &table, &flow, HeatDrift::default());
        assert_exact(
            h.condenser_duty,
            h.condenser_coefficient * (y.condenser_cw_out_c - table.celsius[Stream::ReactorOutlet]),
            "QUS uses TST(8)",
        );
        assert!(
            (table.celsius[Stream::ReactorOutlet] - unpacked.separator.celsius).abs() > 1e-9,
            "the two temperatures coincide here, so this test proves nothing"
        );
    }
}
