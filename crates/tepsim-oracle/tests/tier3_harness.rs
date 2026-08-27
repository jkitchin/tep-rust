//! The Tier 3 harness itself: does the trace record what actually happened?
//!
//! B-0029. This file tests the *measuring instrument*, not the port. The port
//! has nothing to trace yet: `Plant` does not consume the generator until
//! B-0031 ports the walk advance and B-0024b the noise. Building the
//! instrument first is deliberate, and B-0027 gives the reason: the walk state
//! is readable from `COMMON/WLK/` so values stay checkable per item, but draw
//! *order* is not a value in any `COMMON` block and cannot be recovered after
//! the fact.
//!
//! So what is checked here is that the trace is faithful: that it records
//! every draw, in order, with the right scaling, and that it agrees with the
//! generator word it is supposed to explain.

#![cfg(feature = "oracle")]

use tepsim_core::disturbance::{Draw, NOISE_DRAWS, SEGMENT_DRAWS};
use tepsim_core::{TepRng, TracingRng};
use tepsim_oracle::tier2::{Pools, Scenario};
use tepsim_oracle::{Oracle, TRACE_CAPACITY, tier3};

const DT: f64 = 1.0 / 3600.0;

/// Force a scenario with the trace cleared first, and return what it drew.
fn traced(oracle: &mut Oracle, scenario: &Scenario) -> Vec<Draw> {
    oracle.set_teproc(&scenario.common);
    oracle.set_wlk(&scenario.walk);
    oracle.set_rng(scenario.rng);
    oracle.set_manipulated(&scenario.manipulated);
    oracle.set_disturbances(&scenario.disturbances);
    tier3::clear(oracle);
    let _ = oracle.derivatives(scenario.time, &scenario.state);
    tier3::trace(oracle)
}

/// The trace must explain the generator word: replaying its draws from the
/// starting word has to land exactly where the evaluation left it.
///
/// This is the strongest single check available on the instrument. A trace
/// that dropped a draw, recorded one twice, or reordered any of them would
/// land somewhere else.
#[test]
fn replaying_the_trace_reproduces_the_generator_word() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);

    let mut total = 0_usize;
    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let before = scenario.rng;
        let trace = traced(&mut oracle, &scenario);
        let after = oracle.rng();
        total += trace.len();

        let mut replay = TepRng::new(before);
        for (position, draw) in trace.iter().enumerate() {
            let value = if draw.signed {
                replay.signed()
            } else {
                replay.unit()
            };
            assert_eq!(
                value.to_bits(),
                draw.value.to_bits(),
                "trace entry {position} at nominal#{index} does not match a \
                 replay of the same sequence"
            );
        }
        assert_eq!(
            replay.state().to_bits(),
            after.to_bits(),
            "replaying {} traced draws from {before} did not reach {after} at \
             nominal#{index}: the trace does not account for the whole stream",
            trace.len()
        );
    }
    println!("{total} draws traced and replayed across 400 evaluations");
    assert!(total > 100_000, "only {total} draws; the pool is too quiet");
}

/// The counts the trace reports must match B-0027's census, which was measured
/// a completely different way: by stepping a generator from the word before to
/// the word after, with no instrumentation at all.
///
/// Two independent methods agreeing is what makes either trustworthy.
#[test]
fn the_trace_length_agrees_with_the_uninstrumented_census() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);
    let base = pools.nominal_case(30);

    let cases = [
        (0.0, 0_usize, "TIME=0: noise skipped, walks reset"),
        (1.0e-6, 264, "noise only"),
        (0.15, 462, "noise, walk advance and the gas analysers"),
        (0.30, 522, "and the product analyser"),
    ];
    for (time, expected, what) in cases {
        let scenario = Scenario {
            time,
            ..base.clone()
        };
        let trace = traced(&mut oracle, &scenario);
        println!("t={time}: {} draws, {what}", trace.len());
        assert_eq!(
            trace.len(),
            expected,
            "the trace disagrees with B-0027's census at t={time}"
        );
    }
}

/// The sign flag is recorded, and both forms occur.
///
/// If every draw came back as the same form the field would be dead weight and
/// the differ's most useful message would never fire.
#[test]
fn both_scalings_appear_in_a_real_evaluation() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 60, DT);

    // A step where the walks advance: `TESUB5` draws signed, `TESUB6` unit.
    let trace = traced(
        &mut oracle,
        &Scenario {
            time: 0.15,
            ..pools.nominal_case(30)
        },
    );
    let signed = trace.iter().filter(|d| d.signed).count();
    let unit = trace.len() - signed;
    println!("{signed} signed draws, {unit} unit draws");

    assert!(signed > 0 && unit > 0, "only one scaling appeared");
    // The walk advance is the only signed consumer: 9 channels through TESUB5
    // at 3 draws each, plus 1 each from the 3 spike channels at teprob.f:388.
    assert_eq!(
        signed,
        9 * SEGMENT_DRAWS + 3,
        "the signed draws are the walk's"
    );
    assert_eq!(
        unit % NOISE_DRAWS,
        0,
        "the unit draws should all come from TESUB6, twelve at a time"
    );
    // And the signed ones are on [-1, 1) while the unit ones are on [0, 1).
    for draw in &trace {
        if draw.signed {
            assert!((-1.0..1.0).contains(&draw.value), "{:?}", draw.value);
        } else {
            assert!((0.0..1.0).contains(&draw.value), "{:?}", draw.value);
        }
    }
}

/// The differ finds the first divergence, and says something useful about it.
#[test]
fn the_differ_reports_the_first_divergence_and_its_kind() {
    let a = |v: f64| Draw {
        value: v,
        signed: false,
    };
    let s = |v: f64| Draw {
        value: v,
        signed: true,
    };

    assert_eq!(tier3::diff(&[a(0.1), a(0.2)], &[a(0.1), a(0.2)]), None);

    // A wrong value, at index 1 rather than 0.
    let d = tier3::diff(&[a(0.1), a(0.3)], &[a(0.1), a(0.2)]).expect("differs");
    assert_eq!(d.index, 1);
    assert!(format!("{d}").contains("Same scaling, different value"));

    // The wrong scaling: same position, same count, different form.
    let d = tier3::diff(&[a(0.1), s(0.2)], &[a(0.1), a(0.2)]).expect("differs");
    assert_eq!(d.index, 1);
    let message = format!("{d}");
    assert!(message.contains("wrong sign flag"), "{message}");

    // Too few draws, which is B-0028's failure mode.
    let d = tier3::diff(&[a(0.1)], &[a(0.1), a(0.2)]).expect("differs");
    assert_eq!(d.index, 1);
    assert!(format!("{d}").contains("fewer"));

    // Too many.
    let d = tier3::diff(&[a(0.1), a(0.2)], &[a(0.1)]).expect("differs");
    assert_eq!(d.index, 1);
    assert!(format!("{d}").contains("more"));
}

/// `TracingRng` records what `TepRng` produces, and nothing else.
///
/// The port side of the comparison has to be as faithful as the Fortran side,
/// and it is easy to write a wrapper that perturbs the thing it observes.
#[test]
fn the_tracing_generator_records_without_disturbing() {
    let mut plain = TepRng::with_default_seed();
    let mut traced = TracingRng::new(TepRng::with_default_seed());

    for i in 0..1_000 {
        // Alternate the two forms, so both paths are exercised.
        let (a, b) = if i % 3 == 0 {
            (plain.signed(), traced.signed())
        } else {
            (plain.unit(), traced.unit())
        };
        assert_eq!(a.to_bits(), b.to_bits(), "draw {i}");
    }
    assert_eq!(plain.state().to_bits(), traced.state().to_bits());
    assert_eq!(traced.draws().len(), 1_000);
    for (i, draw) in traced.draws().iter().enumerate() {
        assert_eq!(draw.signed, i % 3 == 0, "draw {i} recorded the wrong form");
    }

    // Clearing keeps the generator where it is, which is what makes per-step
    // comparison possible.
    let state = traced.state();
    traced.clear();
    assert!(traced.draws().is_empty());
    assert_eq!(traced.state().to_bits(), state.to_bits());
}

/// The trace buffer is comfortably larger than any evaluation needs, and an
/// overflow would be reported rather than truncating.
#[test]
fn the_trace_capacity_has_real_headroom() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let mut worst = 0;
    for index in 0..pools.trajectory.len() {
        worst = worst.max(traced(&mut oracle, &pools.nominal_case(index)).len());
    }
    println!("worst evaluation: {worst} draws against a capacity of {TRACE_CAPACITY}");
    assert!(
        worst * 4 < TRACE_CAPACITY,
        "the worst evaluation used {worst} of {TRACE_CAPACITY}, which is less \
         than four times of headroom. Longer runs and active disturbances draw \
         more; raise TRCCAP."
    );
}
