//! The stripper, and the reactor-inlet alias.
//!
//! Ported from `teprob.f:614-662`. Steam strips light components out of the
//! separator liquid; what leaves overhead rejoins the mixing zone as stream 5,
//! and what does not becomes the stripper's own liquid, stream 12.
//!
//! # The equations
//!
//! The feed is the mixed A/C feed plus the separator underflow:
//!
//! \\[ f_i = \\dot n_{i,4} + \\dot n_{i,11} \\]
//!
//! A vapour-to-liquid ratio sets how hard the column strips, scaled by a
//! temperature factor that is piecewise in the stripper temperature:
//!
//! \\[ \\Lambda = \\frac{F_4}{F_{11}} \\, \\tau(T_c) \\]
//!
//! Each condensible then strips according to a Langmuir-shaped saturating
//! function of that ratio:
//!
//! \\[ s_i = \\frac{k_i \\Lambda}{1 + k_i \\Lambda}, \\qquad i \\in \\{D \\ldots H\\} \\]
//!
//! and the split is simply
//!
//! \\[ \\dot n_{i,5} = s_i f_i, \\qquad \\dot n_{i,12} = f_i - \\dot n_{i,5} \\]
//!
//! # `SFR(1..3)` are never recomputed
//!
//! `teprob.f:623-627` and `629-633` both write slots 4 through 8 only. Slots 1,
//! 2 and 3 are set once in `TEINIT` (`teprob.f:1126-1128`) and are read at
//! `teprob.f:643` on every evaluation, so A, B and C strip at a *fixed* 99.5%,
//! 99.1% and 99.0% no matter what the column is doing.
//!
//! That is the intended physics rather than an oversight: the non-condensibles
//! are gases, they leave overhead essentially completely, and no temperature
//! or flow ratio in the plant's range would change that. It is worth stating
//! because the loop at `643` runs `I=1,8` and looks as though all eight
//! factors come from the branch above it.
//!
//! # The temperature factor has a pole at 177 C
//!
//! \\[
//!   \\tau(T) = \\begin{cases}
//!     T - 120.262 & T > 170 \\\\
//!     0.1 & T < 5.292 \\\\
//!     \\dfrac{363.744}{177 - T} - 2.22579488 & \\text{otherwise}
//!   \\end{cases}
//! \\]
//!
//! The middle branch diverges at 177 C, which is *inside* the range the two
//! outer branches leave for it only if `TCC` exceeds 170, and it does not: the
//! `T > 170` branch takes over first. So the pole is unreachable by seven
//! degrees, and the two branches are continuous to within 0.1% at 170.
//!
//! B-0016 built an adversarial state at 176 C anyway, to sit near the pole and
//! confirm it stays on the linear branch. That state is coverage of the
//! *guard*, not of the pole.
//!
//! # `FTM(11) > 0.1` switches the whole block
//!
//! Below that threshold the column is not really running, and
//! `teprob.f:629-633` substitutes five fixed stripping factors rather than
//! evaluating the correlation. The reason is visible in the arithmetic:
//! \\(\\Lambda = F_4/F_{11}\\) diverges as \\(F_{11} \\to 0\\).
//!
//! Both sides are covered by the adversarial catalogue, which places a state
//! exactly on `FTM(11) = 0.1`. Since the test is `.GT.`, that state takes the
//! *fixed-factor* branch.
//!
//! # The reactor inlet is an alias, and this is where it is made
//!
//! `teprob.f:656-661` copies flow, enthalpy, temperature, composition and
//! component flows from stream 6 to stream 7 wholesale. There is no mixing,
//! no pressure drop and no heat loss between them: stream 7 exists so that the
//! reactor's balance at `teprob.f:763-772` can name its own inlet. See
//! [`crate::stream::Stream::ReactorInlet`].
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `SFR(8)` | [`Stripper::factors`] | fraction of each species stripped |
//! | `FIN(8)` | [`Stripper::feed`] | combined feed to the column |
//! | `VOVRL` | [`Stripper::vapour_to_liquid`] | vapour-to-liquid ratio |
//! | `TMPFAC` | [`Stripper::temperature_factor`] | temperature scaling |

// Every float expression here reproduces `teprob.f`'s association and
// rounding exactly; see `crate::thermo`.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]
// Every literal in this module is transcribed digit for digit from the listing
// and rounded by `single`, never pre-rounded by hand. The digits an f32 cannot
// hold are the evidence that the constant was copied rather than retyped, and
// dropping them is exactly how a transcription error gets in: B-0019 found one
// that way, where a hand-shortened 31.5859536 gave a different f32.
#![allow(
    clippy::excessive_precision,
    reason = "transcribed verbatim from teprob.f; `single` does the rounding"
)]

use crate::component::{ByComponent, Component, Composition};
use crate::constants::single;
use crate::flows::Flows;
use crate::stream::Stream;
use crate::streams::Streams;
use crate::thermo::{EnergyBasis, enthalpy};

/// Fixed stripping factors for the three non-condensibles
/// (`teprob.f:1126-1128`).
///
/// Set once in `TEINIT` and never recomputed; see the module documentation.
/// All three literals are single precision.
///
/// `TEINIT` also sets slots 4 through 8 (`teprob.f:1129-1133`), and those are
/// dead: `teprob.f:614-634` overwrites all five on every evaluation, through
/// whichever branch it takes, before anything reads them. That is delta D-005.
//
// @port  teprob.f:1126-1133
// @delta D-005 class=A teprob.f:1129-1133
pub const NON_CONDENSIBLE_STRIPPING: [f64; 3] = [
    single(0.99500), // SFR(1), A
    single(0.99100), // SFR(2), B
    single(0.99000), // SFR(3), C
];

/// Langmuir coefficients for the five condensibles (`teprob.f:623-627`).
///
/// All single precision. D and E strip readily, F slightly more so, and G and
/// H barely at all, which is what makes the column a product separator.
const STRIPPING_COEFFICIENT: [f64; 5] = [
    single(8.5010), // D
    single(11.402), // E
    single(11.795), // F
    single(0.0480), // G
    single(0.0242), // H
];

/// Fixed factors used when the column is barely flowing (`teprob.f:629-633`).
///
/// All single precision.
const IDLE_STRIPPING: [f64; 5] = [
    single(0.9999), // D
    single(0.999),  // E
    single(0.999),  // F
    single(0.99),   // G
    single(0.98),   // H
];

/// Below this separator underflow the correlation is not evaluated at all
/// (`teprob.f:614`). Single precision.
///
/// The test is `.GT.`, so a flow of exactly 0.1 takes the fixed-factor branch.
const MINIMUM_UNDERFLOW: f64 = single(0.1);

/// Above this stripper temperature the factor is linear (`teprob.f:615`).
const LINEAR_ABOVE: f64 = single(170.);
/// Offset of the linear branch (`teprob.f:616`). Single precision.
const LINEAR_OFFSET: f64 = single(120.262);
/// Below this the factor pins to its floor (`teprob.f:617`).
const PINNED_BELOW: f64 = single(5.292);
/// The floor itself (`teprob.f:618`).
const PINNED_VALUE: f64 = single(0.1);
/// Numerator of the hyperbolic branch (`teprob.f:620`).
const HYPERBOLIC_SCALE: f64 = single(363.744);
/// The pole of the hyperbolic branch (`teprob.f:620`). Unreachable; see the
/// module documentation.
const HYPERBOLIC_POLE: f64 = single(177.);
/// Offset of the hyperbolic branch (`teprob.f:620`).
const HYPERBOLIC_OFFSET: f64 = single(2.22579488);

/// Which way `teprob.f:614-634` went.
///
/// Reported rather than merely applied, so that a differential can say whether
/// a sampled state exercised each branch. An unexercised branch is not
/// evidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StripperBranch {
    /// `FTM(11) <= 0.1`: the column is idle and the fixed factors are used.
    Idle,
    /// `TCC > 170`: the linear temperature branch.
    Linear,
    /// `TCC < 5.292`: the temperature factor is pinned to 0.1.
    Pinned,
    /// Between the two: the hyperbolic branch.
    Hyperbolic,
}

/// Everything `teprob.f:614-662` produces that is not already a stream.
///
/// The stream quantities are written back into [`Streams`] and [`Flows`],
/// filling in the three slots those two left empty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stripper {
    /// `SFR`: the fraction of each species that leaves overhead.
    pub factors: ByComponent<f64>,
    /// `FIN`: the combined feed to the column.
    pub feed: ByComponent<f64>,
    /// `VOVRL`: the vapour-to-liquid ratio. Meaningless in the idle branch,
    /// where the original leaves it at whatever the last evaluation left.
    pub vapour_to_liquid: f64,
    /// `TMPFAC`: the temperature scaling. Likewise.
    pub temperature_factor: f64,
    /// Which branch was taken.
    pub branch: StripperBranch,
}

/// The temperature factor, `TMPFAC` (`teprob.f:615-621`).
fn temperature_factor(celsius: f64) -> (f64, StripperBranch) {
    if celsius > LINEAR_ABOVE {
        (celsius - LINEAR_OFFSET, StripperBranch::Linear)
    } else if celsius < PINNED_BELOW {
        (PINNED_VALUE, StripperBranch::Pinned)
    } else {
        (
            HYPERBOLIC_SCALE / (HYPERBOLIC_POLE - celsius) - HYPERBOLIC_OFFSET,
            StripperBranch::Hyperbolic,
        )
    }
}

/// Run the stripper and fill in streams 5, 12 and 7.
///
/// `stripper_celsius` is `TCC`. `stream_table` and `flow` are updated in
/// place: they carry the ten streams B-0020 and B-0021 filled in, and this
/// completes all thirteen.
// @port teprob.f:614-662
pub fn stripper(stream_table: &mut Streams, flow: &mut Flows, stripper_celsius: f64) -> Stripper {
    // teprob.f:614-634. The branch decides all five condensible factors; the
    // three non-condensibles never move.
    let underflow = flow.molar[Stream::SeparatorUnderflow];
    let mut factors = ByComponent::new([0.0; Component::COUNT]);
    for (slot, value) in NON_CONDENSIBLE_STRIPPING.iter().enumerate() {
        factors.as_mut_array()[slot] = *value;
    }

    let (vapour_to_liquid, temperature_scaling, branch) = if underflow > MINIMUM_UNDERFLOW {
        let (scaling, branch) = temperature_factor(stripper_celsius);
        // teprob.f:622
        let ratio = flow.molar[Stream::AcFeed] / underflow * scaling;
        // teprob.f:623-627. Written as the listing writes it: the product
        // appears twice rather than being hoisted, which is bit-identical and
        // checkable line by line.
        for (offset, k) in STRIPPING_COEFFICIENT.iter().enumerate() {
            factors.as_mut_array()[offset + 3] = k * ratio / (1.0 + k * ratio);
        }
        (ratio, scaling, branch)
    } else {
        // teprob.f:629-633
        for (offset, fixed) in IDLE_STRIPPING.iter().enumerate() {
            factors.as_mut_array()[offset + 3] = *fixed;
        }
        // `VOVRL` and `TMPFAC` are locals the original does not write on this
        // path. Nothing reads them either, so reporting NaN says "not
        // computed" rather than inventing a number that looks meaningful.
        (f64::NAN, f64::NAN, StripperBranch::Idle)
    };

    // teprob.f:635-639. The accumulation is written out in three statements in
    // the original; the order is `(0 + FCM4) + FCM11`.
    let mut feed = ByComponent::new([0.0; Component::COUNT]);
    for component in Component::ALL {
        let mut total = 0.0;
        total += flow.component[Stream::AcFeed][component];
        total += flow.component[Stream::SeparatorUnderflow][component];
        feed[component] = total;
    }

    // teprob.f:640-647. Both totals accumulate in component order.
    let mut overhead = ByComponent::new([0.0; Component::COUNT]);
    let mut downflow = ByComponent::new([0.0; Component::COUNT]);
    let mut overhead_total = 0.0;
    let mut downflow_total = 0.0;
    for component in Component::ALL {
        overhead[component] = factors[component] * feed[component];
        downflow[component] = feed[component] - overhead[component];
        overhead_total += overhead[component];
        downflow_total += downflow[component];
    }

    // teprob.f:648-651. Unchecked, as everywhere: an idle column can drive a
    // total to zero, and the original does not check either.
    let mut overhead_fractions = [0.0; Component::COUNT];
    let mut downflow_fractions = [0.0; Component::COUNT];
    for component in Component::ALL {
        overhead_fractions[component.index()] = overhead[component] / overhead_total;
        downflow_fractions[component.index()] = downflow[component] / downflow_total;
    }

    flow.component[Stream::StripperOverhead] = overhead;
    flow.component[Stream::StripperDownflow] = downflow;
    flow.molar[Stream::StripperOverhead] = overhead_total;
    flow.molar[Stream::StripperDownflow] = downflow_total;
    stream_table.composition[Stream::StripperOverhead] =
        Composition::new_unchecked(overhead_fractions);
    stream_table.composition[Stream::StripperDownflow] =
        Composition::new_unchecked(downflow_fractions);

    // teprob.f:652-655. The overhead leaves as vapour, the downflow as liquid.
    stream_table.celsius[Stream::StripperOverhead] = stripper_celsius;
    stream_table.celsius[Stream::StripperDownflow] = stripper_celsius;
    stream_table.enthalpy[Stream::StripperOverhead] = enthalpy(
        &stream_table.composition[Stream::StripperOverhead],
        stripper_celsius,
        EnergyBasis::VapourEnthalpy,
    );
    stream_table.enthalpy[Stream::StripperDownflow] = enthalpy(
        &stream_table.composition[Stream::StripperDownflow],
        stripper_celsius,
        EnergyBasis::LiquidEnthalpy,
    );

    // teprob.f:656-661. The reactor inlet is stream 6, copied.
    flow.molar[Stream::ReactorInlet] = flow.molar[Stream::MixingZoneOutlet];
    flow.component[Stream::ReactorInlet] = flow.component[Stream::MixingZoneOutlet];
    stream_table.enthalpy[Stream::ReactorInlet] = stream_table.enthalpy[Stream::MixingZoneOutlet];
    stream_table.celsius[Stream::ReactorInlet] = stream_table.celsius[Stream::MixingZoneOutlet];
    stream_table.composition[Stream::ReactorInlet] =
        stream_table.composition[Stream::MixingZoneOutlet];

    Stripper {
        factors,
        feed,
        vapour_to_liquid,
        temperature_factor: temperature_scaling,
        branch,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::equilibrium::equilibrium;
    use crate::flows::{FlowDrift, flows};
    use crate::state::State;
    use crate::streams::{FeedConditions, streams};
    use crate::testing::assert_exact;
    use crate::vessels::{TemperatureSeeds, unpack};

    fn plant() -> (Streams, Flows, f64) {
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
        let flow = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        (table, flow, unpacked.stripper.celsius)
    }

    /// A, B and C strip at the `TEINIT` factors, whatever the column is doing.
    /// The loop at `teprob.f:643` runs `I=1,8` and looks as though all eight
    /// come from the branch above it.
    #[test]
    fn the_non_condensibles_strip_at_fixed_factors_in_both_branches() {
        let (mut table, mut flow, tcc) = plant();
        let running = stripper(&mut table, &mut flow, tcc);

        let (mut table, mut flow, tcc) = plant();
        flow.molar[Stream::SeparatorUnderflow] = 0.05; // below the threshold
        let idle = stripper(&mut table, &mut flow, tcc);
        assert_eq!(idle.branch, StripperBranch::Idle);

        for (offset, expected) in NON_CONDENSIBLE_STRIPPING.iter().enumerate() {
            assert_exact(running.factors.as_array()[offset], *expected, "running");
            assert_exact(idle.factors.as_array()[offset], *expected, "idle");
        }
        // And the condensibles do differ between the two, or this proves
        // nothing.
        assert!(
            running.factors[Component::D].to_bits() != idle.factors[Component::D].to_bits(),
            "the two branches gave the same D factor"
        );
    }

    /// All three temperature branches, at values chosen to land in each.
    #[test]
    fn the_temperature_factor_has_three_branches() {
        assert_eq!(temperature_factor(200.0).1, StripperBranch::Linear);
        assert_exact(
            temperature_factor(200.0).0,
            200.0 - LINEAR_OFFSET,
            "linear branch",
        );

        assert_eq!(temperature_factor(0.0).1, StripperBranch::Pinned);
        assert_exact(temperature_factor(0.0).0, PINNED_VALUE, "pinned branch");

        assert_eq!(temperature_factor(65.0).1, StripperBranch::Hyperbolic);

        // The boundaries themselves. Both tests are strict, so a temperature
        // exactly on one goes to the hyperbolic branch.
        assert_eq!(
            temperature_factor(LINEAR_ABOVE).1,
            StripperBranch::Hyperbolic
        );
        assert_eq!(
            temperature_factor(PINNED_BELOW).1,
            StripperBranch::Hyperbolic
        );
    }

    /// The pole at 177 C is unreachable: the linear branch takes over at 170,
    /// seven degrees before it. Getting the guard backwards would produce an
    /// infinity in a perfectly ordinary plant.
    #[test]
    fn the_pole_at_177_is_guarded_by_the_linear_branch() {
        for celsius in [170.1, 175.0, 176.999, 177.0, 177.001, 200.0] {
            let (value, branch) = temperature_factor(celsius);
            assert_eq!(
                branch,
                StripperBranch::Linear,
                "{celsius} C must not reach the hyperbolic branch"
            );
            assert!(value.is_finite(), "{celsius} C gave {value}");
        }
        // The two branches nearly agree at the handover, which is the evidence
        // that 170 and 120.262 belong together.
        let hyperbolic = HYPERBOLIC_SCALE / (HYPERBOLIC_POLE - 170.0) - HYPERBOLIC_OFFSET;
        let linear = 170.0 - LINEAR_OFFSET;
        let gap = (hyperbolic - linear).abs() / linear;
        assert!(gap < 0.01, "the branches disagree by {gap:e} at 170 C");
    }

    /// The threshold is `.GT.`, so a flow of exactly 0.1 is *idle*.
    #[test]
    fn a_flow_of_exactly_the_threshold_takes_the_idle_branch() {
        let (mut table, mut flow, tcc) = plant();
        flow.molar[Stream::SeparatorUnderflow] = MINIMUM_UNDERFLOW;
        assert_eq!(
            stripper(&mut table, &mut flow, tcc).branch,
            StripperBranch::Idle
        );

        let (mut table, mut flow, tcc) = plant();
        flow.molar[Stream::SeparatorUnderflow] = MINIMUM_UNDERFLOW * 1.000_001;
        assert_ne!(
            stripper(&mut table, &mut flow, tcc).branch,
            StripperBranch::Idle
        );
    }

    /// What goes in comes out: the two product streams partition the feed
    /// exactly, component by component.
    #[test]
    fn the_two_outlets_partition_the_feed() {
        let (mut table, mut flow, tcc) = plant();
        let s = stripper(&mut table, &mut flow, tcc);
        for component in Component::ALL {
            assert_exact(
                flow.component[Stream::StripperOverhead][component]
                    + flow.component[Stream::StripperDownflow][component],
                s.feed[component],
                "the split is exact by construction",
            );
        }
        // And the feed is the two inlets.
        for component in Component::ALL {
            assert_exact(
                s.feed[component],
                flow.component[Stream::AcFeed][component]
                    + flow.component[Stream::SeparatorUnderflow][component],
                "FIN",
            );
        }
    }

    /// Every stripping factor lies in [0, 1]: it is a fraction. The Langmuir
    /// form guarantees it for any non-negative ratio, and a sign error would
    /// break it.
    #[test]
    fn every_stripping_factor_is_a_fraction() {
        let (mut table, mut flow, tcc) = plant();
        let s = stripper(&mut table, &mut flow, tcc);
        for component in Component::ALL {
            let f = s.factors[component];
            assert!(
                (0.0..=1.0).contains(&f),
                "{} strips at {f}, which is not a fraction",
                component.name()
            );
        }
        // The heavies must strip less than the lights, or the column is not a
        // product separator.
        assert!(s.factors[Component::D] > s.factors[Component::H]);
        assert!(s.factors[Component::E] > s.factors[Component::G]);
    }

    /// Stream 7 is stream 6, in every field. Anything else is a copy that
    /// missed one.
    #[test]
    fn the_reactor_inlet_copies_every_field_of_the_mixing_zone_outlet() {
        let (mut table, mut flow, tcc) = plant();
        let _ = stripper(&mut table, &mut flow, tcc);
        assert_exact(
            flow.molar[Stream::ReactorInlet],
            flow.molar[Stream::MixingZoneOutlet],
            "FTM(7)",
        );
        assert_exact(
            table.enthalpy[Stream::ReactorInlet],
            table.enthalpy[Stream::MixingZoneOutlet],
            "HST(7)",
        );
        assert_exact(
            table.celsius[Stream::ReactorInlet],
            table.celsius[Stream::MixingZoneOutlet],
            "TST(7)",
        );
        assert_eq!(
            table.composition[Stream::ReactorInlet],
            table.composition[Stream::MixingZoneOutlet],
            "XST(*,7)"
        );
        assert_eq!(
            flow.component[Stream::ReactorInlet],
            flow.component[Stream::MixingZoneOutlet],
            "FCM(*,7)"
        );
    }

    /// After this module runs, every one of the thirteen streams has been
    /// filled in. That is the property B-0025 depends on.
    #[test]
    fn all_thirteen_streams_are_populated_afterwards() {
        let (mut table, mut flow, tcc) = plant();
        let _ = stripper(&mut table, &mut flow, tcc);
        for stream in Stream::ALL {
            assert!(
                table.enthalpy[stream].is_finite() && table.enthalpy[stream] != 0.0,
                "{} still has no enthalpy",
                stream.name()
            );
        }
    }
}
