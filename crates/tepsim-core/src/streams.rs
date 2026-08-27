//! Stream compositions, molecular weights, temperatures and enthalpies.
//!
//! Ported from `teprob.f:529-564`, with the feed data from
//! `teprob.f:1134-1169`. This is the bookkeeping layer between the vessels and
//! the flow network: it takes what [`mod@crate::vessels`] and
//! [`mod@crate::equilibrium`] know about the four vessels and writes it into
//! the per-stream arrays the balances read.
//!
//! # Which streams this covers, and which it does not
//!
//! Ten of the thirteen. The four feeds are constants or walk-driven; six more
//! are the vessels' own contents leaving them:
//!
//! | Stream | Composition | Temperature |
//! |---|---|---|
//! | 1, 2, 3 | fixed feed ([`FEED_COMPOSITION`]) | fixed, 45 C |
//! | 4 | fixed, except A/B/C from the walks | from the walks |
//! | 6 | `XVV`, the mixing zone vapour | `TCV` |
//! | 8 | `XVR`, the reactor vapour | `TCR` |
//! | 9, 10 | `XVS`, the separator vapour | `TCS` |
//! | 11 | `XLS`, the separator liquid | `TCS` |
//! | 13 | `XLC`, the stripper liquid | `TCC` |
//!
//! Streams 5, 7 and 12 belong to the stripper and are set at
//! `teprob.f:654-661`, which is B-0022.
//!
//! # `XMWS` exists for six streams only
//!
//! `teprob.f:529-534` zeroes exactly six slots and `542-547` accumulates into
//! the same six. The other seven are never written anywhere in the file. They
//! are also never *read* anywhere in the file, so this is dead storage rather
//! than a latent zero-divide, and it is left at zero here for the same reason
//! it is zero there.
//!
//! The six that exist are the ones a mass flow has to be converted through:
//! streams 1 and 2 for the feed-rate measurements, and 6, 8, 9 and 10 for the
//! four pressure-driven flows.
//!
//! # `HST(10) = HST(9)` is a snapshot, not an alias
//!
//! `teprob.f:562` copies the separator vapour enthalpy into the purge. Both
//! streams leave the separator vapour space, so at that moment they are the
//! same fluid at the same temperature and the copy is exact.
//!
//! Then `teprob.f:601` does this:
//!
//! ```fortran
//!       HST(9)=HST(9)+CPDH/FTM(9)
//! ```
//!
//! The recycle gains the compressor work; the purge does not, because it was
//! copied first. Reading line 562 as an alias rather than as a copy would give
//! the purge a share of compressor work it never receives, and the two lines
//! are seventy apart, so the ordering is easy to miss.
//!
//! The consequence for validation is immediate and useful: after `TEFUNC`
//! returns, `COMMON`'s `HST(9)` carries the compressor bump and `HST(10)` does
//! not. So `HST(10)` is exactly the value this module computes for stream 9,
//! and Tier 2 checks it there. `HST(9)` itself becomes checkable in B-0021,
//! once the compressor exists.
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `XST(8,13)` | [`Streams::composition`] | mole fractions per stream |
//! | `XMWS(13)` | [`Streams::molar_mass`] | mean molecular weight |
//! | `TST(13)` | [`Streams::celsius`] | temperature |
//! | `HST(13)` | [`Streams::enthalpy`] | specific enthalpy |

// Every float expression here reproduces `teprob.f`'s association and
// rounding exactly; see `crate::thermo`.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]

use crate::component::{ByComponent, Component, Composition};
use crate::constants::{XMW, single};
use crate::equilibrium::Equilibrium;
use crate::stream::{ByStream, Stream};
use crate::thermo::{EnergyBasis, enthalpy};
use crate::vessels::Unpacked;

/// The four feed compositions as `TEINIT` sets them (`teprob.f:1134-1168`).
///
/// Indexed as `FEED_COMPOSITION[stream - 1][component]`, covering streams 1
/// through 4 only.
///
/// # These do not sum to one
///
/// Every literal is single precision, so `0.9999` and `0.0001` widen to values
/// whose sum is 1.0000000116860974, not 1. That residual of 1.2e-8 is the
/// original's own data and is the reason [`Composition`] checks its sum with a
/// tolerance rather than exactly; see that type's documentation.
///
/// Stream 4's A, B and C entries are overwritten on every evaluation by the
/// disturbance walks at `teprob.f:407-410`, so the values here are only the
/// starting point. D through H stay zero: the mixed feed carries no heavies.
//
// The claim covers the whole feed block, `1134-1169`, compositions and the
// four `TST` lines interleaved with them. `FEED_CELSIUS` claims one of those
// lines again, more finely; overlapping claims merge.
//
// @port teprob.f:1134-1169
pub const FEED_COMPOSITION: [[f64; Component::COUNT]; 4] = [
    // Stream 1, the D feed. teprob.f:1134-1141.
    [
        single(0.0),
        single(0.0001),
        single(0.0),
        single(0.9999),
        single(0.0),
        single(0.0),
        single(0.0),
        single(0.0),
    ],
    // Stream 2, the E feed. teprob.f:1143-1150.
    [
        single(0.0),
        single(0.0),
        single(0.0),
        single(0.0),
        single(0.9999),
        single(0.0001),
        single(0.0),
        single(0.0),
    ],
    // Stream 3, the A feed. teprob.f:1152-1159.
    [
        single(0.9999),
        single(0.0001),
        single(0.0),
        single(0.0),
        single(0.0),
        single(0.0),
        single(0.0),
        single(0.0),
    ],
    // Stream 4, the A and C feed. teprob.f:1161-1168. A, B and C are replaced
    // by the walks on every evaluation.
    [
        single(0.4850),
        single(0.0050),
        single(0.5100),
        single(0.0),
        single(0.0),
        single(0.0),
        single(0.0),
        single(0.0),
    ],
];

/// The feed temperature `TEINIT` sets for all four feeds, degrees Celsius.
///
/// `teprob.f:1142`, `1151`, `1160` and `1169` all write `45.`, without a `D`
/// suffix. Streams 1 and 4 are then overwritten by the walks at
/// `teprob.f:411-412`, so only streams 2 and 3 keep this value.
///
/// The other three lines are inside [`FEED_COMPOSITION`]'s claim.
//
// @port teprob.f:1142
pub const FEED_CELSIUS: f64 = single(45.);

/// The feed conditions the disturbance walks supply.
///
/// `teprob.f:407-412` overwrites four numbers on every evaluation from
/// `TESUB8` and the `IDV` flags. That is disturbance machinery and therefore
/// Phase 3, so until then it is an input, exactly as
/// [`crate::kinetics::ReactionDrift`] is.
///
/// [`Default`] is the nominal operating point: the walks sit at `SZERO` and no
/// `IDV` is active, so the values are what `TEINIT` left behind.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FeedConditions {
    /// `XST(1..3, 4)`: A, B and C in the mixed feed, after `TESUB8(1..2)` and
    /// the two composition step disturbances. C is the remainder
    /// (`teprob.f:410`), never an independent draw.
    pub ac_feed_light: [f64; 3],
    /// `TST(1)`: the D feed temperature, after `TESUB8(3)` and `IDV(3)`.
    pub d_feed_celsius: f64,
    /// `TST(4)`: the mixed feed temperature, after `TESUB8(4)`.
    pub ac_feed_celsius: f64,
}

impl Default for FeedConditions {
    fn default() -> Self {
        Self {
            ac_feed_light: [
                FEED_COMPOSITION[3][0],
                FEED_COMPOSITION[3][1],
                FEED_COMPOSITION[3][2],
            ],
            d_feed_celsius: FEED_CELSIUS,
            ac_feed_celsius: FEED_CELSIUS,
        }
    }
}

/// Everything `teprob.f:529-564` produces.
///
/// Streams 5, 7 and 12 are left at their defaults; the stripper fills them in
/// (B-0022).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Streams {
    /// `XST`: mole fractions leaving each stream.
    pub composition: ByStream<Composition>,
    /// `XMWS`: mean molecular weight. Only six slots are meaningful; see the
    /// module documentation.
    pub molar_mass: ByStream<f64>,
    /// `TST`: temperature, degrees Celsius.
    pub celsius: ByStream<f64>,
    /// `HST`: specific enthalpy.
    ///
    /// Stream 9 is the value *before* the compressor work is added at
    /// `teprob.f:601`; see the module documentation.
    pub enthalpy: ByStream<f64>,
}

/// Mean molecular weight of a stream, `teprob.f:542-547`.
///
/// Accumulated in component order from zero, because that is the order the
/// `DO 2010` loop runs in and a reassociated sum would change the last bits.
fn molar_mass(composition: &Composition) -> f64 {
    let mut total = 0.0;
    for component in Component::ALL {
        total += composition[component] * XMW[component];
    }
    total
}

/// Assemble the stream table for everything known at this point in the
/// evaluation.
// @port teprob.f:529-564
#[must_use]
pub fn streams(unpacked: &Unpacked, eq: &Equilibrium, feed: &FeedConditions) -> Streams {
    let mut composition = ByStream::new([Composition::default(); Stream::COUNT]);
    let mut celsius = ByStream::new([0.0; Stream::COUNT]);

    // teprob.f:1134-1169, as modified by teprob.f:407-412. The feeds are not
    // recomputed here; they are what `TEINIT` and the walks left in `COMMON`.
    for (offset, stream) in [Stream::DFeed, Stream::EFeed, Stream::AFeed, Stream::AcFeed]
        .into_iter()
        .enumerate()
    {
        composition[stream] = Composition::new_unchecked(FEED_COMPOSITION[offset]);
    }
    // teprob.f:407-410. C is the remainder of one, not an independent value.
    let mixed = composition[Stream::AcFeed].fractions().as_array();
    let mut mixed = *mixed;
    mixed[0] = feed.ac_feed_light[0];
    mixed[1] = feed.ac_feed_light[1];
    mixed[2] = feed.ac_feed_light[2];
    composition[Stream::AcFeed] = Composition::new_unchecked(mixed);

    celsius[Stream::DFeed] = feed.d_feed_celsius;
    celsius[Stream::EFeed] = FEED_CELSIUS;
    celsius[Stream::AFeed] = FEED_CELSIUS;
    celsius[Stream::AcFeed] = feed.ac_feed_celsius;

    // teprob.f:536-541. Six streams take a vessel's composition wholesale.
    // Unchecked for the same reason as elsewhere: an adversarial state can
    // leave a vessel's fractions un-normalised, and that is a true answer.
    let vessel_compositions = [
        (Stream::MixingZoneOutlet, unpacked.mixing.fractions),
        (Stream::ReactorOutlet, eq.reactor.fractions),
        (Stream::Recycle, eq.separator.fractions),
        (Stream::Purge, eq.separator.fractions),
        (Stream::SeparatorUnderflow, unpacked.separator.fractions),
        (Stream::Product, unpacked.stripper.fractions),
    ];
    for (stream, fractions) in vessel_compositions {
        composition[stream] = fractions;
    }

    // teprob.f:549-554.
    celsius[Stream::MixingZoneOutlet] = unpacked.mixing.celsius;
    celsius[Stream::ReactorOutlet] = unpacked.reactor.celsius;
    celsius[Stream::Recycle] = unpacked.separator.celsius;
    celsius[Stream::Purge] = unpacked.separator.celsius;
    celsius[Stream::SeparatorUnderflow] = unpacked.separator.celsius;
    celsius[Stream::Product] = unpacked.stripper.celsius;

    // teprob.f:529-534 and 542-547. Six slots, and the rest stay zero.
    let mut molar = ByStream::new([0.0; Stream::COUNT]);
    for stream in [
        Stream::DFeed,
        Stream::EFeed,
        Stream::MixingZoneOutlet,
        Stream::ReactorOutlet,
        Stream::Recycle,
        Stream::Purge,
    ] {
        molar[stream] = molar_mass(&composition[stream]);
    }

    // teprob.f:555-564. Vapour basis for the feeds and the vapour streams,
    // liquid basis for the two liquid ones.
    let mut hst = ByStream::new([0.0; Stream::COUNT]);
    for stream in [
        Stream::DFeed,
        Stream::EFeed,
        Stream::AFeed,
        Stream::AcFeed,
        Stream::MixingZoneOutlet,
        Stream::ReactorOutlet,
        Stream::Recycle,
    ] {
        hst[stream] = enthalpy(
            &composition[stream],
            celsius[stream],
            EnergyBasis::VapourEnthalpy,
        );
    }
    // teprob.f:562. A copy taken before the compressor bump, not an alias; see
    // the module documentation.
    hst[Stream::Purge] = hst[Stream::Recycle];
    for stream in [Stream::SeparatorUnderflow, Stream::Product] {
        hst[stream] = enthalpy(
            &composition[stream],
            celsius[stream],
            EnergyBasis::LiquidEnthalpy,
        );
    }

    Streams {
        composition,
        molar_mass: molar,
        celsius,
        enthalpy: hst,
    }
}

/// The six streams whose molecular weight the original computes.
///
/// Public so that a validation harness can compare exactly these and leave the
/// dead slots alone rather than asserting something about storage nobody
/// reads.
pub const WEIGHED_STREAMS: [Stream; 6] = [
    Stream::DFeed,
    Stream::EFeed,
    Stream::MixingZoneOutlet,
    Stream::ReactorOutlet,
    Stream::Recycle,
    Stream::Purge,
];

/// The ten streams this module fills in. Streams 5, 7 and 12 are the
/// stripper's.
pub const ASSEMBLED_STREAMS: [Stream; 10] = [
    Stream::DFeed,
    Stream::EFeed,
    Stream::AFeed,
    Stream::AcFeed,
    Stream::MixingZoneOutlet,
    Stream::ReactorOutlet,
    Stream::Recycle,
    Stream::Purge,
    Stream::SeparatorUnderflow,
    Stream::Product,
];

/// Unused, but kept honest: the mole fractions of a feed as a plain array.
#[must_use]
pub fn feed_composition(stream: Stream) -> Option<ByComponent<f64>> {
    let index = stream.fortran_index();
    (1..=4)
        .contains(&index)
        .then(|| ByComponent::new(FEED_COMPOSITION[index - 1]))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::equilibrium::equilibrium;
    use crate::state::State;
    use crate::testing::assert_exact;
    use crate::vessels::{TemperatureSeeds, unpack};

    fn solved() -> (Unpacked, Equilibrium, Streams) {
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
        let unpacked = unpack(&y, TemperatureSeeds::default()).expect("converges");
        let eq = equilibrium(&unpacked);
        let s = streams(&unpacked, &eq, &FeedConditions::default());
        (unpacked, eq, s)
    }

    /// The purge takes the recycle's enthalpy *before* the compressor work is
    /// added, which happens seventy lines later. Both are the separator vapour
    /// at the separator temperature, so at this point they must be identical.
    #[test]
    fn the_purge_enthalpy_is_a_snapshot_of_the_recycle_enthalpy() {
        let (_, _, s) = solved();
        assert_exact(
            s.enthalpy[Stream::Purge],
            s.enthalpy[Stream::Recycle],
            "HST(10) = HST(9)",
        );
        assert_exact(
            s.celsius[Stream::Purge],
            s.celsius[Stream::Recycle],
            "same fluid, same temperature",
        );
    }

    /// Only six streams get a molecular weight, and the other seven stay
    /// exactly zero. Filling them in would look like an improvement and would
    /// diverge from the original.
    #[test]
    fn only_six_streams_carry_a_molecular_weight() {
        let (_, _, s) = solved();
        for stream in Stream::ALL {
            if WEIGHED_STREAMS.contains(&stream) {
                assert!(
                    s.molar_mass[stream] > 0.0,
                    "{} should have a molecular weight",
                    stream.name()
                );
            } else {
                assert_exact(
                    s.molar_mass[stream],
                    0.0,
                    "the original never writes this slot",
                );
            }
        }
    }

    /// The feeds are what they are named after. Getting the internal stream
    /// numbering backwards is the mistake `crate::stream` exists to prevent,
    /// and this is where it would first show up as a wrong number.
    #[test]
    fn each_feed_carries_the_component_it_is_named_for() {
        let (_, _, s) = solved();
        assert!(s.composition[Stream::DFeed][Component::D] > 0.99, "D feed");
        assert!(s.composition[Stream::EFeed][Component::E] > 0.99, "E feed");
        assert!(s.composition[Stream::AFeed][Component::A] > 0.99, "A feed");
        let mixed = &s.composition[Stream::AcFeed];
        assert!(
            mixed[Component::A] > 0.4 && mixed[Component::C] > 0.5,
            "A/C feed"
        );
        assert_exact(mixed[Component::D], 0.0, "the mixed feed carries no D");
    }

    /// The feed literals are single precision, so they do not sum to one.
    /// Enforcing an exact sum would reject the original's own data.
    #[test]
    fn the_feed_compositions_miss_one_by_the_single_precision_residual() {
        for (index, feed) in FEED_COMPOSITION.iter().enumerate() {
            let composition = Composition::new_unchecked(*feed);
            let residual = (composition.sum() - 1.0).abs();
            assert!(
                residual < Composition::SUM_TOLERANCE,
                "feed {} sums to {}",
                index + 1,
                composition.sum()
            );
            if index != 3 {
                assert!(
                    residual > 0.0,
                    "feed {} summed to exactly one, so the literals were not \
                     transcribed as single precision",
                    index + 1
                );
            }
        }
    }

    /// The six vessel-fed streams take their composition and temperature from
    /// the vessel they leave. Crossing two of them is a copy-paste error that
    /// leaves the plant running and wrong.
    #[test]
    fn the_vessel_streams_carry_their_vessels_conditions() {
        let (unpacked, eq, s) = solved();
        assert_exact(
            s.celsius[Stream::ReactorOutlet],
            unpacked.reactor.celsius,
            "stream 8 leaves the reactor",
        );
        assert_exact(
            s.celsius[Stream::Product],
            unpacked.stripper.celsius,
            "stream 13 leaves the stripper",
        );
        assert_eq!(
            s.composition[Stream::ReactorOutlet],
            eq.reactor.fractions,
            "stream 8 carries the reactor *vapour*, not its liquid",
        );
        assert_eq!(
            s.composition[Stream::SeparatorUnderflow],
            unpacked.separator.fractions,
            "stream 11 carries the separator *liquid*, not its vapour",
        );
    }

    /// The two liquid streams use the liquid enthalpy correlation and the rest
    /// use the vapour one. Using the wrong basis is a large error that still
    /// produces plausible numbers.
    #[test]
    fn the_liquid_streams_use_the_liquid_enthalpy_basis() {
        let (_, _, s) = solved();
        for stream in [Stream::SeparatorUnderflow, Stream::Product] {
            assert_exact(
                s.enthalpy[stream],
                enthalpy(
                    &s.composition[stream],
                    s.celsius[stream],
                    EnergyBasis::LiquidEnthalpy,
                ),
                "liquid basis",
            );
            assert!(
                s.enthalpy[stream].to_bits()
                    != enthalpy(
                        &s.composition[stream],
                        s.celsius[stream],
                        EnergyBasis::VapourEnthalpy,
                    )
                    .to_bits(),
                "the two bases agreed, so this test proves nothing"
            );
        }
    }

    /// Streams 5, 7 and 12 belong to the stripper and must be left alone here.
    #[test]
    fn the_stripper_streams_are_not_touched() {
        let (_, _, s) = solved();
        for stream in [
            Stream::StripperOverhead,
            Stream::ReactorInlet,
            Stream::StripperDownflow,
        ] {
            assert_exact(s.enthalpy[stream], 0.0, "left for B-0022");
            assert_exact(s.celsius[stream], 0.0, "left for B-0022");
        }
    }
}
