//! Tier 9: this platform against the committed digests.
//!
//! Every one of these runs on whatever machine `cargo test` was invoked on, so
//! the native half of Tier 9 is complete for any architecture that runs the
//! suite. No comparison here is between two values computed in this process;
//! the constants in [`tepsim::tier9::CASES`] were measured once and written
//! down, and that is what makes an x86-64 run of this file extend the claim to
//! x86-64 rather than merely observe that x86-64 agrees with itself.
//!
//! The wasm32 half needs a WebAssembly runtime and lives in `cargo xtask
//! tier9`, which builds `tepsim-wasm` for `wasm32-unknown-unknown` and has the
//! runtime evaluate the same table against the same constants.

use tepsim::run::Sample;
use tepsim::tier9::{self, CASES, Case, Fnv1a64, SUITE_DIGEST};
use tepsim::{Recorder, Scenario, Simulation};

/// The headline. If this fails, the platform, the compiler or the model
/// changed the numbers the simulator produces, and no constant should be
/// touched until it is known which.
#[test]
fn every_case_reproduces_its_committed_digest() {
    if let Err(mismatch) = tier9::check() {
        panic!(
            "{mismatch}\n\n\
             This is the most important failure this project can report. Two \
             platforms, or two builds, disagree about what the simulator \
             computes. Find out which changed before editing the constant: \
             moving it is a logged re-baseline, not a fix."
        );
    }
}

/// The one-number summary, which is what goes in a log entry.
#[test]
fn the_suite_digest_is_the_committed_one() {
    assert_eq!(
        tier9::suite_digest(),
        SUITE_DIGEST,
        "the suite digest moved. Either a case digest moved, which \
         `every_case_reproduces_its_committed_digest` will say, or the table \
         was renamed or reordered."
    );
}

/// A digest that did not depend on the numbers would pass every test above
/// while proving nothing. Perturbing one input has to move it.
#[test]
fn the_digest_responds_to_the_run() {
    let base = tier9::digest(Scenario::baseline().with_hours(1.0));
    assert_eq!(base, CASES[0].digest, "the same scenario, the same digest");

    // A different seed is a different noise stream from the first sample on.
    let reseeded = tier9::digest(Scenario::baseline().with_hours(1.0).with_seed(1.0));
    assert_ne!(base, reseeded, "the seed must reach the digest");

    // A different length is a different number of rows.
    let longer = tier9::digest(Scenario::baseline().with_hours(2.0));
    assert_ne!(base, longer, "the duration must reach the digest");

    // And a fault must, or the fault cases in the table are decoration.
    let faulted = tier9::digest(Scenario::fault(1).with_hours(1.0));
    assert_ne!(base, faulted, "a disturbance must reach the digest");
}

/// The table must actually cover different code, or six cases are one case
/// written six times.
#[test]
fn the_cases_are_distinct() {
    for (i, case) in CASES.iter().enumerate() {
        assert!(!case.name.is_empty(), "case {i} has no name");
        assert!(!case.covers.is_empty(), "case `{}` says nothing", case.name);
        for other in &CASES[i + 1..] {
            assert_ne!(case.name, other.name, "duplicate case name");
            assert_ne!(
                case.digest, other.digest,
                "cases `{}` and `{}` have the same digest, so one of them is \
                 not exercising what it claims to",
                case.name, other.name
            );
        }
    }
}

/// Every case has to be a valid, whole number of sampling intervals, or the
/// chunked consumer in `tepsim-wasm` and the batch one here record different
/// numbers of rows and the shared constant means two different things.
#[test]
fn every_case_is_a_whole_number_of_sampling_intervals() {
    for case in CASES {
        let scenario = case.scenario();
        scenario
            .validate()
            .unwrap_or_else(|e| panic!("case `{}` is invalid: {e:?}", case.name));
        assert_eq!(
            scenario.steps() % scenario.sample_every,
            0,
            "case `{}` ends part way through a sampling interval",
            case.name
        );
        assert!(
            scenario.samples() > 0,
            "case `{}` records nothing",
            case.name
        );
    }
}

/// [`Case::compute`] must be a pure function of the case. Running it twice in
/// one process proves nothing about another platform, but it does catch a
/// digest that accumulated global state, which would make every cross-platform
/// comparison meaningless.
///
/// One case rather than all six: `every_case_reproduces_its_committed_digest`
/// already walks the table, and the cases share every line of the code this is
/// checking. Repeating the nine-hour run to learn the same thing would double
/// what this file costs on every `cargo test`.
#[test]
fn computing_a_case_twice_gives_the_same_answer() {
    let case = &CASES[0];
    assert_eq!(case.compute(), case.compute(), "case `{}`", case.name);
    assert!(case.agrees(), "case `{}`", case.name);
}

/// The streaming digest and a digest taken over a collected run must agree, or
/// [`Recorder`] is doing something to the samples on the way past.
#[test]
fn streaming_and_collected_digests_agree() {
    let scenario = Scenario::baseline().with_hours(1.0);
    let run = Simulation::new(scenario).run();

    let mut collected = Fnv1a64::new();
    for sample in &run.samples {
        collected.record(sample);
    }

    assert_eq!(collected.finish(), tier9::digest(scenario));
    assert_eq!(collected.finish(), CASES[0].digest);
}

/// The row layout is load-bearing: `tepsim-wasm` hashes the packed row it hands
/// a browser, and this hashes the sample. They must be the same bytes in the
/// same order or the shared constant is a coincidence.
#[test]
fn the_recorder_hashes_time_then_the_row() {
    let sample = Sample {
        step: 7,
        hours: 0.25,
        measurements: core::array::from_fn(|i| i as f64),
        manipulated: core::array::from_fn(|i| 100.0 + i as f64),
        labels: tepsim::Labels::none(),
    };

    let mut by_recorder = Fnv1a64::new();
    by_recorder.record(&sample);

    let mut by_hand = Fnv1a64::new();
    by_hand.write_f64(sample.hours);
    by_hand.write_slice(&sample.measurements);
    by_hand.write_slice(&sample.manipulated);

    assert_eq!(by_recorder.finish(), by_hand.finish());

    // And the step number is not in there. It is redundant with the time and
    // is a `usize`, whose width differs between wasm32 and a 64-bit host.
    let mut renumbered = sample;
    renumbered.step = 9;
    let mut other = Fnv1a64::new();
    other.record(&renumbered);
    assert_eq!(by_recorder.finish(), other.finish());
}

/// FNV-1a is a published algorithm with published test vectors. Checking one
/// pins the constants, so a mistyped prime cannot make every digest in the
/// project self-consistently wrong.
#[test]
fn the_hash_matches_the_published_fnv1a_vector() {
    let mut hash = Fnv1a64::new();
    for byte in b"foobar" {
        hash.write_u8(*byte);
    }
    assert_eq!(
        hash.finish(),
        0x8594_4171_f739_67e8,
        "FNV-1a 64 of \"foobar\""
    );
    assert_eq!(
        Fnv1a64::new().finish(),
        0xcbf2_9ce4_8422_2325,
        "the offset basis"
    );
}

/// Bit patterns, not values. A signed zero appearing on one architecture and
/// not another is exactly the finding Tier 9 exists to surface, so the digest
/// must not normalise it away the way `tepsim_scenario::Digest` does.
#[test]
fn the_digest_distinguishes_signed_zero_and_is_order_sensitive() {
    let of = |values: &[f64]| {
        let mut hash = Fnv1a64::new();
        hash.write_slice(values);
        hash.finish()
    };
    assert_ne!(of(&[0.0]), of(&[-0.0]));
    assert_ne!(of(&[1.0, 2.0]), of(&[2.0, 1.0]));
    assert_ne!(of(&[]), of(&[0.0]));
    assert_eq!(of(&[]), Fnv1a64::new().finish());

    // The scenario digest does normalise, on purpose. Asserting the contrast
    // here means a future edit that unified the two has to face the question.
    let mut scenario_digest = tepsim::Digest::new();
    scenario_digest.push_f64(-0.0);
    let mut positive = tepsim::Digest::new();
    positive.push_f64(0.0);
    assert_eq!(scenario_digest.value(), positive.value());
}

/// `write_str` has to be length-terminated, or two case names that concatenate
/// the same way collide in the suite digest.
#[test]
fn the_suite_digest_cannot_be_confused_by_concatenation() {
    let of = |parts: &[&str]| {
        let mut hash = Fnv1a64::new();
        for part in parts {
            hash.write_str(part);
        }
        hash.finish()
    };
    assert_ne!(of(&["ab", "c"]), of(&["a", "bc"]));
}

/// A `Case` is data, so it is easy to add one and forget the constant. Zero is
/// the value a copy-paste leaves behind, and no real run produces it.
#[test]
fn no_case_has_a_placeholder_digest() {
    for Case { name, digest, .. } in CASES {
        assert_ne!(*digest, 0, "case `{name}` still has a placeholder digest");
        assert_ne!(
            *digest,
            Fnv1a64::new().finish(),
            "case `{name}` digests as the empty stream, so it recorded nothing"
        );
    }
}
