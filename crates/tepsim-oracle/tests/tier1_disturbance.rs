//! Tier 1 for the two generator-consuming utilities: `TESUB5` and `TESUB6`,
//! `teprob.f:1506-1546`.
//!
//! # Why this is not an input sweep
//!
//! `TESUB1` through `TESUB4` are pure functions of their arguments, so Tier 1
//! for them means sweeping ten million inputs. These two *draw*, so their
//! answer depends on the generator word as much as on their arguments, and a
//! sweep would be comparing two different streams.
//!
//! So each case pins `G` on both sides first. What is then swept is the seed
//! together with the arguments, which covers the same ground: a different seed
//! is a different set of draws.
//!
//! # The draw count is part of the contract
//!
//! Tier 3 will check that the port and the Fortran make the *same draws in the
//! same order* across a whole run. That rests on each routine consuming a
//! fixed, known number, so the count is asserted here at the source rather
//! than inferred later from a trace that has already gone wrong.
//!
//! Both counts are recovered without instrumenting anything: `TESUB7` is a
//! pure function of its own word, so stepping a port-side generator from the
//! word before a call to the word after gives the count exactly. B-0027
//! established the method.

#![cfg(feature = "oracle")]

use tepsim_core::disturbance::{
    ChannelSpans, NOISE_DRAWS, SEGMENT_DRAWS, SegmentStart, noise, walk_segment,
};
use tepsim_core::{TepRng, constants};
use tepsim_oracle::tier1::{Comparison, Sampler};
use tepsim_oracle::{Oracle, WalkSegmentStart, WalkSpans};

/// `PLAN.org`, "Tier 1".
const TIER1_TOLERANCE: f64 = 1e-13;

/// Identifies one compared number.
#[derive(Clone, Copy, Debug)]
struct Case {
    seed: u64,
    what: &'static str,
}

impl core::fmt::Display for Case {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "seed#{}[{}]", self.seed, self.what)
    }
}

/// How many draws took `before` to `after`.
fn draws_between(before: f64, after: f64) -> usize {
    let mut probe = TepRng::new(before);
    for step in 0..100 {
        if probe.state() == after {
            return step;
        }
        let _ = probe.unit();
    }
    panic!("more than a hundred draws between {before} and {after}");
}

/// A uniform draw on `[lo, hi)`. `Sampler` offers only `[0, 1)`.
fn between(sampler: &mut Sampler, lo: f64, hi: f64) -> f64 {
    lo + (hi - lo) * sampler.unit()
}

/// A generator word to pin both sides to.
///
/// Spread across the modulus rather than clustered near the default seed, so
/// the sweep sees the whole state space.
fn seed_for(index: u64) -> f64 {
    // The modulus is 2^32; step by a large odd stride so the seeds do not
    // fall into a short cycle of their own.
    ((index.wrapping_mul(2_654_435_761) % 4_294_967_291) + 1) as f64
}

#[test]
fn tesub6_matches_the_fortran_over_the_state_space() {
    let mut oracle = Oracle::lock();
    // `XNS` and the walk spans are `TEINIT`'s, and a bare lock has them at
    // zero. Reading them before this call gives a sweep whose "real" half is
    // all zeros, which agrees with the Fortran and proves nothing.
    let _ = oracle.init();
    let mut comparison: Comparison<Case> = Comparison::new("TESUB6 noise sample");
    let mut sampler = Sampler::new(0x7E2_0028);

    // Every one of the 41 published noise magnitudes, plus drawn ones.
    for index in 0..20_000_u64 {
        let seed = seed_for(index);
        let std = if index % 2 == 0 {
            // XNS(1..41) are the magnitudes the model actually uses. They are
            // B-0030's to transcribe, so they are read from the oracle here
            // rather than duplicated.
            oracle.teproc().xns[(index as usize / 2) % 41]
        } else {
            between(&mut sampler, 0.0, 25.0)
        };

        oracle.set_rng(seed);
        let theirs = oracle.tesub6(std);
        let after_theirs = oracle.rng();

        let mut rng = TepRng::new(seed);
        let ours = noise(&mut rng, std);

        comparison.observe(
            Case {
                seed: index,
                what: "X",
            },
            ours,
            theirs,
        );
        assert_eq!(
            rng.state().to_bits(),
            after_theirs.to_bits(),
            "the generator ended in a different place at seed {seed}, so the \
             draw sequence differs even where the value happens to agree"
        );
        assert_eq!(draws_between(seed, rng.state()), NOISE_DRAWS);
    }

    println!("{comparison}");
    comparison.assert_within(TIER1_TOLERANCE);
    assert_eq!(
        comparison.max_ulp(),
        0,
        "TESUB6 is twelve additions and a multiply, with no transcendental \
         anywhere, so anything but bit equality is a porting error"
    );
}

#[test]
fn tesub5_matches_the_fortran_over_the_state_space() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();
    let mut fields = [
        Comparison::<Case>::new("TESUB5 ADIST"),
        Comparison::<Case>::new("TESUB5 BDIST"),
        Comparison::<Case>::new("TESUB5 CDIST"),
        Comparison::<Case>::new("TESUB5 DDIST"),
        Comparison::<Case>::new("TESUB5 TNEXT"),
    ];
    let mut sampler = Sampler::new(0x7E2_0028 + 1);

    for index in 0..20_000_u64 {
        let seed = seed_for(index);
        // Half the cases use a real channel's spans, half are drawn, so the
        // sweep covers both the operating point and the wider space.
        let (spans, flag) = if index % 2 == 0 {
            let wlk = oracle.wlk();
            let channel = (index as usize / 2) % 12;
            (
                ChannelSpans {
                    duration_span: wlk.hspan[channel],
                    duration_centre: wlk.hzero[channel],
                    value_span: wlk.sspan[channel],
                    value_centre: wlk.szero[channel],
                    slope_span: wlk.spspan[channel],
                },
                i32::from(index % 4 == 0),
            )
        } else {
            (
                ChannelSpans {
                    // `HZERO` must stay clear of zero: `H` divides the cubic
                    // coefficients at `teprob.f:1533-1534`, and the original
                    // has no guard. The published spans keep it well away.
                    duration_span: between(&mut sampler, 0.0, 2.0),
                    duration_centre: between(&mut sampler, 2.5, 5.0),
                    value_span: between(&mut sampler, 0.0, 1.0),
                    value_centre: between(&mut sampler, -2.0, 2.0),
                    slope_span: between(&mut sampler, 0.0, 0.5),
                },
                i32::from(index % 3 == 0),
            )
        };
        let start = SegmentStart {
            value: between(&mut sampler, -2.0, 2.0),
            slope: between(&mut sampler, -1.0, 1.0),
            since: between(&mut sampler, 0.0, 50.0),
        };

        oracle.set_rng(seed);
        let theirs = oracle.tesub5(
            WalkSegmentStart {
                value: start.value,
                slope: start.slope,
                tlast: start.since,
            },
            WalkSpans {
                hspan: spans.duration_span,
                hzero: spans.duration_centre,
                sspan: spans.value_span,
                szero: spans.value_centre,
                spspan: spans.slope_span,
                idvflag: flag,
            },
        );
        let after_theirs = oracle.rng();

        let mut rng = TepRng::new(seed);
        let ours = walk_segment(&mut rng, start, &spans, flag);

        let pairs = [
            (ours.constant, theirs.adist, "ADIST"),
            (ours.linear, theirs.bdist, "BDIST"),
            (ours.quadratic, theirs.cdist, "CDIST"),
            (ours.cubic, theirs.ddist, "DDIST"),
            (ours.until, theirs.tnext, "TNEXT"),
        ];
        assert!(
            spans.duration_centre > 0.0,
            "a case has a zero-centred duration, so `H` can reach zero and \
             the cubic coefficients divide by it"
        );
        for (comparison, (a, b, what)) in fields.iter_mut().zip(pairs) {
            comparison.observe(Case { seed: index, what }, a, b);
        }

        assert_eq!(
            rng.state().to_bits(),
            after_theirs.to_bits(),
            "the generator ended in a different place at seed {seed}"
        );
        assert_eq!(draws_between(seed, rng.state()), SEGMENT_DRAWS);
    }

    for comparison in &fields {
        println!("{comparison}");
    }
    for comparison in &fields {
        comparison.assert_within(TIER1_TOLERANCE);
        assert_eq!(
            comparison.max_ulp(),
            0,
            "TESUB5 is three draws and some arithmetic, with no transcendental \
             anywhere, so anything but bit equality is a porting error"
        );
    }
}

/// The disturbance flag must scale the endpoint draws and not the duration
/// draw, and both flag values must consume three draws.
///
/// Reading `teprob.f:1528-1530` the other way round, so that an inactive
/// channel skipped its draws, would give identical *values* and a completely
/// different stream for everything downstream.
#[test]
fn the_flag_changes_the_segment_but_never_the_draw_count() {
    let mut oracle = Oracle::lock();
    let spans = WalkSpans {
        hspan: 0.2,
        hzero: 0.5,
        sspan: 0.03,
        szero: 0.485,
        spspan: 0.1,
        idvflag: 0,
    };
    let start = WalkSegmentStart {
        value: 0.485,
        slope: 0.0,
        tlast: 0.0,
    };

    for seed in [1.0, 12345.0, TepRng::DEFAULT_SEED] {
        let mut ends = Vec::new();
        let mut durations = Vec::new();
        for flag in [0, 1] {
            oracle.set_rng(seed);
            let segment = oracle.tesub5(
                start,
                WalkSpans {
                    idvflag: flag,
                    ..spans
                },
            );
            assert_eq!(
                draws_between(seed, oracle.rng()),
                SEGMENT_DRAWS,
                "flag {flag} at seed {seed} did not consume three draws"
            );
            ends.push(segment.cdist);
            durations.push(segment.tnext);
        }
        // The duration draw comes first and is not scaled by the flag, so the
        // segment ends at the same time either way.
        assert_eq!(
            durations[0].to_bits(),
            durations[1].to_bits(),
            "the flag changed the segment duration at seed {seed}, so it is \
             being applied to the first draw rather than the last two"
        );
        // The curvature does change, because the endpoint did.
        assert_ne!(
            ends[0].to_bits(),
            ends[1].to_bits(),
            "the flag changed nothing at seed {seed}"
        );
    }
}

/// An inactive channel lands exactly on `SZERO`, in the Fortran too.
///
/// That is what makes a disturbance-free run a flat line rather than a small
/// wander, and it is the property `crate::disturbance` documents.
#[test]
fn an_inactive_channel_lands_exactly_on_its_centre_in_the_fortran() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();
    let wlk = oracle.wlk();
    assert!(
        wlk.hzero.iter().all(|h| *h > 0.0),
        "the walk spans are still zero, so TEINIT did not run and every \
         segment here divides by zero"
    );

    for channel in 0..12 {
        let centre = wlk.szero[channel];
        let spans = WalkSpans {
            hspan: wlk.hspan[channel],
            hzero: wlk.hzero[channel],
            sspan: wlk.sspan[channel],
            szero: centre,
            spspan: wlk.spspan[channel],
            idvflag: 0,
        };
        oracle.set_rng(seed_for(channel as u64 + 7));
        let segment = oracle.tesub5(
            WalkSegmentStart {
                value: centre,
                slope: 0.0,
                tlast: 0.0,
            },
            spans,
        );
        // Starting at the centre with zero slope and ending there too, the
        // cubic is identically constant.
        assert_eq!(
            segment.adist.to_bits(),
            centre.to_bits(),
            "channel {channel}"
        );
        assert_eq!(
            segment.bdist.to_bits(),
            0.0_f64.to_bits(),
            "channel {channel}"
        );
        assert!(
            segment.cdist.abs() < 1e-15 && segment.ddist.abs() < 1e-15,
            "channel {channel} is not flat with its disturbance off: \
             c={} d={}",
            segment.cdist,
            segment.ddist
        );
    }
    // The nominal state must itself be disturbance-free, or the whole idea of
    // a nominal run is wrong.
    assert_eq!(constants::NOMINAL_STATE.len(), 50);
}
