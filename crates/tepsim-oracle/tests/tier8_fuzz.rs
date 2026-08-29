//! Tier 8: differential fuzzing of the fifty derivatives against the oracle.
//!
//! `PLAN.org` asks for "random (state, inputs, IDV mask, RNG word) tuples
//! compared against the oracle, with automatic shrinking of any counterexample
//! to a minimal reproducer that lands in the regression suite. This is how we
//! find the branch we did not think to test." That is the whole of it: Tier 2
//! compares states somebody chose, and this compares states nobody did.
//!
//! `tepsim_oracle::tier8` holds the machinery and the reasoning. This file is
//! the gate, the teeth, and the corpus.
//!
//! # What is asserted, and what is only counted
//!
//! Asserted: over every tuple on which both implementations produce a running,
//! converged, finite answer, every one of the fifty components agrees to 1e-12
//! of the scale of the terms. That is the Tier 2 gate, asked of unguided
//! states.
//!
//! Counted, not asserted: how many tuples never got that far. A state drawn
//! over twelve orders of magnitude is usually not physical, `TESUB2` does not
//! converge on it, and the original has no way to say so (delta D-001). The
//! census is printed on every run because the fraction of tuples that carry the
//! claim is the number the next session needs, and a generator that quietly
//! stopped producing physical states would otherwise pass forever.
//!
//! # What the first full run found
//!
//! Five million tuples at seed `0x7E2_0062`: 2,551,618 compared, 2,263,855
//! frozen, 184,526 lost to `TESUB2` falling through, and **one
//! counterexample**. Shrunk from 22 differing knobs to 3, it is the nominal
//! state with `IDV(13)` on, one particular generator word, and a time of 973
//! hours, which leaves the kinetic drift at 3e7 and the reaction rates at
//! 1.1e10. `YP(9)` misses by 4.61e-12 of the scale of its terms, against a gate
//! of 1e-12.
//!
//! It is a libm finding, not an algebra finding, and
//! `the_open_finding_is_the_vendored_libm_and_not_the_algebra` proves that
//! rather than asserting it: under `--features oracle,libm-system` the same
//! tuple is bit-identical in all fifty components. The gate has not been
//! touched, so `TEP_TIER8=full` reports the counterexample. That is the
//! intended behaviour of a tier whose job is to find things.
//!
//! # Run twice
//!
//! Default and `--features oracle,libm-system`, like every Phase 2
//! differential. See `tier2_equilibrium.rs`. The two configurations make
//! genuinely different claims here: 1e-12 of the scale of the terms against the
//! build that ships, and 0 ULP against the build that shares gfortran's `exp`.
//!
//! # Size
//!
//! `TEP_TIER8` unset is 400 tuples, and the whole file runs in under a second.
//! `TEP_TIER8=full` is five million, about six minutes. A decimal count
//! selects a specific size.
//!
//! Two tests run at fixed sizes regardless, because what they measure is not a
//! statistical claim but a property of the harness: the branch demonstration
//! needs twenty thousand tuples to reach `teprob.f:615` reliably, and the
//! corpus is however many entries it has.

// Differential tests: exact comparison against the Fortran is the property
// under test, and the tuples in the corpus are transcribed literals.
#![allow(
    clippy::float_cmp,
    reason = "bit equality against the Fortran is the property under test"
)]
#![allow(
    clippy::unreadable_literal,
    reason = "corpus tuples are machine-emitted and pasted verbatim"
)]
#![cfg(feature = "oracle")]

use tepsim_core::math;
use tepsim_oracle::Oracle;
use tepsim_oracle::tier2::Scenario;
use tepsim_oracle::tier8::{
    Branch, Budget, Disagreement, Generator, Mutation, Outcome, TOLERANCE, Tuple, check, run,
    shrink,
};

/// The seed the gate runs at. Fixed, so that a failure is reproducible from the
/// log entry alone.
const SEED: u64 = 0x7E2_0062;

/// The nominal starting condition, which every tuple is an overlay on.
///
/// `tier8::nominal_scenario` rather than `Pools::collect`, and that is not a
/// stylistic preference: `TEINIT` leaves the four Newton warm starts wherever
/// the previous run in this process put them. See that function for the two
/// ways it bit this file.
fn base(oracle: &mut Oracle) -> Scenario {
    tepsim_oracle::tier8::nominal_scenario(oracle)
}

// ---------------------------------------------------------------------------
// The generator
// ---------------------------------------------------------------------------

/// A tuple must be a function of its seed and its index alone.
///
/// Without this the shrinker's output is not a reproducer and a logged
/// counterexample cannot be replayed, which would make every other test in this
/// file decorative.
#[test]
fn a_tuple_is_reproducible_from_its_seed_and_index() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    let generator = Generator::new(SEED, &base);

    for index in [0, 1, 7, 4_242, 999_983] {
        let first = generator.tuple(index);
        let second = Generator::new(SEED, &base).tuple(index);
        assert_eq!(
            first, second,
            "tuple {index} is not a function of the seed and the index"
        );
    }

    // And a different seed must give different tuples, or the seed is
    // decoration and a nightly at a fresh seed would re-run the same search.
    let other = Generator::new(SEED ^ 1, &base);
    let differing = (0..64)
        .filter(|i| generator.tuple(*i) != other.tuple(*i))
        .count();
    assert_eq!(differing, 64, "two seeds produced overlapping tuples");
}

/// The generator must actually move the state, over the range it claims to.
///
/// A generator that drew everything within a percent of nominal would pass the
/// gate forever and search nothing. So this measures the spread it produces and
/// requires it, rather than trusting the shapes to be wired up.
#[test]
fn the_generator_covers_orders_of_magnitude_and_every_slot() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    let generator = Generator::new(SEED, &base);
    let nominal = *generator.nominal();

    let mut moved = [false; 50];
    let mut widest = 0.0_f64;
    let mut zero_time = 0;
    let mut faulted = 0;
    for index in 0..2_000 {
        let tuple = generator.tuple(index);
        for (slot, seen) in moved.iter_mut().enumerate() {
            if tuple.state[slot].to_bits() != nominal.state[slot].to_bits() {
                *seen = true;
            }
            // Only the extensive slots are scaled multiplicatively; the
            // temperatures and the valve positions are drawn on their own
            // scales and a ratio of them means nothing.
            if slot < 36 && nominal.state[slot] != 0.0 {
                let ratio = (tuple.state[slot] / nominal.state[slot]).abs();
                if ratio > 0.0 && ratio.is_finite() {
                    widest = widest.max(ratio.log10().abs());
                }
            }
        }
        if tuple.time == 0.0 {
            zero_time += 1;
        }
        if tuple.disturbances.iter().any(|d| *d != 0) {
            faulted += 1;
        }
    }

    let still: Vec<usize> = (0..50).filter(|s| !moved[*s]).map(|s| s + 1).collect();
    println!("widest excursion: 10^{widest:.1} of nominal");
    println!("{zero_time} of 2000 tuples at TIME = 0, {faulted} carrying a fault");
    assert!(
        still.is_empty(),
        "YY{still:?} are never moved, so nothing about them is being searched"
    );
    assert!(
        widest > 5.0,
        "the widest excursion is only 10^{widest:.1}; the generator is not \
         searching the range it documents"
    );
    // `teprob.f:397-406` resets the whole walk state at TIME = 0, and the
    // ordinary path is everywhere else. Both have to be in the sample.
    assert!(
        zero_time > 100 && zero_time < 1_000,
        "{zero_time} of 2000 tuples at TIME = 0 leaves one of the two walk \
         paths effectively unsampled"
    );
    assert!(faulted > 100, "only {faulted} of 2000 tuples carry a fault");
}

/// A tuple must be a complete specification of one evaluation.
///
/// `TEFUNC` is not reproducible unless it is made so: it advances the walks,
/// draws noise and ticks the analysers as side effects, and it warm-starts four
/// Newton solves from the previous call's answers. `Scenario::force` restores
/// all of that, and `tier2::reproducible` asserts it. Tier 8 rests on the same
/// property twice over: the shrinker's predicate has to be a function of the
/// tuple, and a corpus entry has to give the same answer in ten years as it
/// does now. So it is checked directly on generated tuples, not inherited from
/// Tier 2's pools.
#[test]
fn forcing_the_same_tuple_twice_gives_the_same_answer() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    let generator = Generator::new(SEED, &base);

    for index in 0..64 {
        let tuple = generator.tuple(index);
        tepsim_oracle::tier2::reproducible(&mut oracle, &tuple.scenario(&base));
    }

    // And the verdict, not merely the derivative, has to be stable: the
    // shrinker asks "does this still fail" hundreds of times per finding, in an
    // order nothing controls.
    for index in 0..64 {
        let tuple = generator.tuple(index);
        let first = check(
            &mut oracle,
            &base,
            &tuple,
            TOLERANCE,
            Mutation::None,
            0,
            None,
        );
        // An unrelated evaluation in between, so a stale-state bug cannot pass
        // by the Fortran simply not having moved.
        let _ = check(
            &mut oracle,
            &base,
            &generator.tuple(index + 1),
            TOLERANCE,
            Mutation::None,
            0,
            None,
        );
        let second = check(
            &mut oracle,
            &base,
            &tuple,
            TOLERANCE,
            Mutation::None,
            0,
            None,
        );
        assert_eq!(
            core::mem::discriminant(&first),
            core::mem::discriminant(&second),
            "tuple {index} produced {first:?} and then {second:?}"
        );
        if let (Outcome::Agreed { worst: a, .. }, Outcome::Agreed { worst: b, .. }) =
            (first, second)
        {
            assert_eq!(a, b, "tuple {index} changed its worst error between runs");
        }
    }
}

// ---------------------------------------------------------------------------
// The gate
// ---------------------------------------------------------------------------

/// The tier itself: unguided tuples, compared against the oracle.
///
/// Passes at the default size. At `TEP_TIER8=full` it reports the one
/// counterexample described in this file's header, which is attributed to the
/// vendored `exp` by
/// `the_open_finding_is_the_vendored_libm_and_not_the_algebra` and is left
/// failing on purpose: the gate is 1e-12 of the scale of the terms and nothing
/// here moves it.
#[test]
fn random_tuples_agree_with_the_fortran() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    let budget = Budget::selected();
    let generator = Generator::new(SEED, &base);

    let report = run(
        &mut oracle,
        &base,
        &generator,
        budget,
        TOLERANCE,
        Mutation::None,
        8,
    );

    println!(
        "transcendentals come from the {} libm",
        if math::USES_SYSTEM_LIBM {
            "platform"
        } else {
            "vendored"
        }
    );
    println!("{report}");

    // A run that compared nothing would pass every tolerance ever set.
    assert!(
        report.compared > budget.tuples / 20,
        "only {} of {} tuples reached a comparison; the generator has drifted \
         off the physical domain and the gate is measuring almost nothing",
        report.compared,
        budget.tuples
    );
    // The freeze has to be exercised too, or `teprob.f:807-811` is never
    // compared against the port's copy of it.
    assert!(
        report.frozen > 0,
        "no tuple tripped the plant, so the shutdown freeze went uncompared"
    );

    if report.counterexamples.is_empty() {
        assert_eq!(report.disagreed, 0);
        return;
    }

    // A counterexample is a finding, not something to tune away. Shrink it,
    // print the literal, and fail with the reproducer in the output.
    let (index, tuple, finding) = report.counterexamples[0];
    let shrunk = shrink(
        &mut oracle,
        &base,
        generator.nominal(),
        &tuple,
        TOLERANCE,
        Mutation::None,
    );
    println!("counterexample at fuzz#{index}: {finding}");
    println!(
        "shrunk from {} knobs to {}",
        shrunk.knobs_before, shrunk.knobs
    );
    for line in shrunk.tuple.differences(generator.nominal()) {
        println!("  {line}");
    }
    println!(
        "{}",
        shrunk
            .tuple
            .as_rust_literal("FOUND", "found by tier8; see LOG.org")
    );
    panic!(
        "tier8 found {} counterexample(s); the smallest is {} and fails with {}",
        report.disagreed, shrunk.knobs, shrunk.finding
    );
}

/// With the transcendentals taken out of the comparison, every compared
/// component must be bit-identical.
///
/// The same claim `tier2_balances.rs` makes, over states nobody chose. If this
/// holds and the default run does not, the difference is `exp` and `pow` and
/// not the algebra.
#[test]
#[cfg(feature = "libm-system")]
fn every_compared_component_is_bit_identical_under_the_platform_libm() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    let generator = Generator::new(SEED, &base);
    let report = run(
        &mut oracle,
        &base,
        &generator,
        Budget::selected(),
        TOLERANCE,
        Mutation::None,
        8,
    );
    println!("{report}");
    assert_eq!(
        report.comparison.max_ulp(),
        0,
        "the derivative is not bit-identical under the platform libm, so the \
         difference is in the algebra somewhere between `teprob.f:407` and \
         `811` and not in `exp` or `pow`"
    );
}

/// The non-finite policy, stated as tests rather than only as prose.
///
/// Two hundred and forty million comparisons in the first full run produced
/// zero non-finite values on either side, so every branch of `slot_verdict`
/// below the finite one is unexercised by the search itself. An unexercised
/// policy is a policy nobody has checked, and this one decides whether a
/// wild-state overflow is a counterexample or not, so it is checked directly.
#[test]
fn the_non_finite_policy_is_what_the_documentation_says() {
    use tepsim_oracle::tier8::slot_verdict;

    let tol = TOLERANCE;
    // Two NaNs agree: IEEE-754 does not say which payload arithmetic produces,
    // and gfortran and LLVM need not pick the same one.
    assert!(slot_verdict(f64::NAN, f64::NAN, 1.0, tol).is_none());
    assert!(slot_verdict(f64::NAN, -f64::NAN, 1.0, tol).is_none());

    // A NaN against a number is a disagreement about whether an answer exists.
    for (ours, theirs) in [(f64::NAN, 1.0), (1.0, f64::NAN)] {
        let (error, kind) =
            slot_verdict(ours, theirs, 1.0, tol).expect("a NaN against a number must be a finding");
        assert_eq!(kind, Disagreement::Existence);
        assert!(error.is_infinite());
    }

    // Matching infinities agree; opposite ones do not.
    assert!(slot_verdict(f64::INFINITY, f64::INFINITY, 1.0, tol).is_none());
    assert_eq!(
        slot_verdict(f64::INFINITY, f64::NEG_INFINITY, 1.0, tol).map(|v| v.1),
        Some(Disagreement::Existence)
    );
    assert_eq!(
        slot_verdict(f64::INFINITY, 1e308, 1.0, tol).map(|v| v.1),
        Some(Disagreement::Existence)
    );

    // The finite path is the gate: the error is against the supplied scale and
    // not against the value, which is the decision of 2026-08-27. The same two
    // numbers pass or fail depending only on the scale they are handed, which
    // is the whole content of that decision.
    assert!(slot_verdict(1.0, 1.0 + 1e-11, 1.0, tol).is_some());
    assert!(slot_verdict(1.0, 1.0 + 1e-11, 1e3, tol).is_none());
    // And a difference below the gate is not a finding at either scale.
    assert!(slot_verdict(1.0, 1.0 + 1e-13, 1.0, tol).is_none());
    // A zero scale means the balance has no terms, so only bit equality passes.
    assert!(slot_verdict(0.0, 0.0, 0.0, tol).is_none());
    assert!(slot_verdict(1e-300, 0.0, 0.0, tol).is_some());
}

// ---------------------------------------------------------------------------
// Teeth
// ---------------------------------------------------------------------------

/// The component the mutations corrupt: `YP(9)`, the reactor energy balance.
///
/// Chosen because every tuple that runs at all reaches it, so the demonstration
/// does not depend on the search finding a rare state.
const MUTATED: usize = 9;

/// Ten times the gate. A term this size must be caught.
const ABOVE_THE_GATE: f64 = 1e-11;

/// A tenth of the gate. A term this size must not be.
const BELOW_THE_GATE: f64 = 1e-13;

/// A wrong constant, everywhere, is found.
///
/// This is the weakest thing a differential fuzzer has to be able to do, and it
/// is asserted rather than assumed because "no counterexamples" is not a result
/// unless the search could have produced one.
#[test]
fn the_search_finds_a_term_of_the_wrong_size() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    let generator = Generator::new(SEED, &base);
    let mutation = Mutation::WrongConstant {
        component: MUTATED,
        relative: ABOVE_THE_GATE,
    };

    let report = run(
        &mut oracle,
        &base,
        &generator,
        Budget::SMOKE,
        TOLERANCE,
        mutation,
        4,
    );
    println!("{report}");
    assert!(
        report.disagreed > 0,
        "a term {ABOVE_THE_GATE:e} of the scale of YP({MUTATED}), on every \
         state, went unnoticed over {} tuples. The harness has no teeth.",
        Budget::SMOKE.tuples
    );
    assert_eq!(
        report.counterexamples[0].2.component, MUTATED,
        "the search blamed the wrong component"
    );
}

/// A term below the gate is not reported.
///
/// The other half of the teeth demonstration, and the one that keeps the first
/// half honest: a harness that failed on everything would also "find" the
/// mutation above. The gate is 1e-12 of the scale of the terms, so 1e-13 of it
/// is inside the budget and must pass.
#[test]
fn a_term_below_the_gate_is_not_reported() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    let generator = Generator::new(SEED, &base);
    let report = run(
        &mut oracle,
        &base,
        &generator,
        Budget::SMOKE,
        TOLERANCE,
        Mutation::WrongConstant {
            component: MUTATED,
            relative: BELOW_THE_GATE,
        },
        4,
    );
    println!("{report}");
    assert_eq!(
        report.disagreed, 0,
        "a term {BELOW_THE_GATE:e} of the scale of YP({MUTATED}) was reported \
         as a counterexample, so the harness fails on noise and its findings \
         mean nothing"
    );
}

/// How many tuples the branch demonstration needs.
///
/// Measured, not guessed. `TCC > 170` is reached by 25 of every 10,222 compared
/// tuples at this seed, so 400 would find it about once and would be a coin
/// flip on any change to the generator. Twenty thousand takes 1.4 seconds and
/// finds it 25 times.
const BRANCH_TUPLES: Budget = Budget { tuples: 20_000 };

/// A fault behind a branch is found, and shrinking isolates the branch.
///
/// This is the case `PLAN.org` names: a mistake that only fires when the plant
/// is somewhere the hand-built pools do not go. Most tuples pass; the search
/// has to reach the branch on its own, and the shrinker has to hand back a
/// tuple that still crosses it.
#[test]
fn the_search_finds_a_fault_behind_a_branch_and_shrinks_to_it() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    let generator = Generator::new(SEED, &base);
    let mutation = Mutation::DroppedTermBehindBranch {
        component: MUTATED,
        relative: ABOVE_THE_GATE,
        branch: Branch::StripperAbove170C,
    };

    let report = run(
        &mut oracle,
        &base,
        &generator,
        BRANCH_TUPLES,
        TOLERANCE,
        mutation,
        4,
    );
    println!("{report}");
    assert!(
        report.disagreed > 0,
        "the branch at teprob.f:615 was never reached in {} tuples",
        BRANCH_TUPLES.tuples
    );
    // Most tuples must still pass, or the "branch" is not a branch and the
    // demonstration is the same as the unconditional one.
    assert!(
        report.disagreed < report.compared,
        "every compared tuple took the branch, so this is not a test of \
         reaching a rare state"
    );

    let (index, tuple, finding) = report.counterexamples[0];
    let shrunk = shrink(
        &mut oracle,
        &base,
        generator.nominal(),
        &tuple,
        TOLERANCE,
        mutation,
    );
    println!("found at fuzz#{index}: {finding}");
    println!(
        "shrunk {} knobs -> {} in {} evaluations, still failing with {}",
        shrunk.knobs_before, shrunk.knobs, shrunk.evaluations, shrunk.finding
    );
    for line in shrunk.tuple.differences(generator.nominal()) {
        println!("  {line}");
    }
    println!(
        "{}",
        shrunk.tuple.as_rust_literal(
            "BRANCH_REPRODUCER",
            "teprob.f:615, TCC above 170 C, reached by scaling the stripper"
        )
    );

    assert!(
        shrunk.knobs < shrunk.knobs_before,
        "shrinking removed nothing: {} knobs in, {} out",
        shrunk.knobs_before,
        shrunk.knobs
    );
    // The whole point of shrinking: what is left has to be about the branch.
    // `TCC` comes from `ESC = ETC/UTLC` (`teprob.f:458`), so the eight stripper
    // holdups and the stripper energy, slots 19 to 27, are the only knobs that
    // can put a state above 170 C. A minimal reproducer must still differ in
    // one of them.
    let touches_stripper =
        (18..27).any(|slot| !shrunk.tuple.knob_matches(generator.nominal(), slot));
    assert!(
        touches_stripper,
        "the shrunk tuple no longer differs anywhere in the stripper, so it \
         cannot be crossing teprob.f:615 and the shrink lost the fault"
    );
    assert_eq!(
        shrunk.finding.kind,
        Disagreement::Magnitude,
        "the shrink wandered off the fault it started from"
    );
}

/// Shrinking must terminate and must be idempotent.
///
/// A shrinker that could be run again and get smaller is not reporting a
/// minimum, and one that does not terminate is worse than none at all.
#[test]
fn shrinking_reaches_a_fixpoint() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    let generator = Generator::new(SEED, &base);
    let mutation = Mutation::WrongConstant {
        component: MUTATED,
        relative: ABOVE_THE_GATE,
    };
    let report = run(
        &mut oracle,
        &base,
        &generator,
        Budget::SMOKE,
        TOLERANCE,
        mutation,
        1,
    );
    let (_, tuple, _) = report.counterexamples[0];

    let once = shrink(
        &mut oracle,
        &base,
        generator.nominal(),
        &tuple,
        TOLERANCE,
        mutation,
    );
    let twice = shrink(
        &mut oracle,
        &base,
        generator.nominal(),
        &once.tuple,
        TOLERANCE,
        mutation,
    );
    println!(
        "{} knobs -> {} -> {}",
        once.knobs_before, once.knobs, twice.knobs
    );
    assert_eq!(
        twice.knobs, once.knobs,
        "a second shrink got further, so the first did not reach a fixpoint"
    );
}

// ---------------------------------------------------------------------------
// The regression corpus
// ---------------------------------------------------------------------------

/// One tuple that is checked forever after, and why.
///
/// A counterexample found by the search stops being a fuzz result the moment it
/// lands here: it becomes a fixed input with a fixed expectation, replayable
/// with no generator and no seed. [`Tuple::as_rust_literal`] emits these, so
/// what is pasted in is exactly the bits the search ran.
struct Case {
    /// What it is, for the failure message.
    name: &'static str,
    /// Why it is in the corpus.
    why: &'static str,
    /// The input.
    tuple: Tuple,
    /// A corruption this tuple is known to detect.
    ///
    /// The corpus's own test. An entry that passed under every mutation would
    /// be checking nothing, and would go on checking nothing forever without
    /// anybody noticing.
    detects: Mutation,
}

/// Every tuple in the corpus.
///
/// One entry, and it is deliberately seeded rather than a real counterexample:
/// it is what `the_search_finds_a_fault_behind_a_branch_and_shrinks_to_it`
/// shrank a mutant's counterexample down to. So it is a real product of the
/// search and not a hand-written state, and it comes with a pair of facts that
/// is what makes the mechanism worth having before there is anything to put in
/// it: with the port as it ships it agrees with the Fortran, and with the
/// mutation that produced it, it does not.
///
/// The one *real* counterexample B-0062 turned up is not here, because a corpus
/// entry is a tuple that agrees and this one does not. It has its own section
/// at the end of the file, with the evidence that it is the vendored `exp`
/// rather than the algebra.
fn corpus() -> Vec<Case> {
    vec![Case {
        name: "stripper above the 170 C stripping-factor branch",
        why: "teprob.f:615, TMPFAC goes linear above 170 C. Shrunk from the \
              tuple the branch mutant found at seed 0x7E2_0062, index 688: \
              eighteen knobs in, one out.",
        tuple: BRANCH_REPRODUCER,
        detects: Mutation::DroppedTermBehindBranch {
            component: MUTATED,
            relative: ABOVE_THE_GATE,
            branch: Branch::StripperAbove170C,
        },
    }]
}

/// The one corpus tuple, emitted verbatim by `Tuple::as_rust_literal`.
///
/// Only `YY(27)`, the stripper's internal energy `ETC`, differs from the
/// nominal condition: 1.1875 against 0.3755, which is what carries `TCC` past
/// the 170 C branch at `teprob.f:615`. Everything else is the nominal starting
/// condition, and the single-precision-looking values in it are not a
/// transcription error: `teprob.f:1053` writes `YY(1)=10.40491389` with no `D0`
/// exponent, so gfortran reads it as a single and promotes.
const BRANCH_REPRODUCER: Tuple = Tuple {
    time: 0.0,
    state: [
        10.404913902282715,
        4.3639960289001465,
        7.570059776306152,
        0.42300423979759216,
        24.155134201049805,
        2.9425976276397705,
        154.37705993652344,
        159.1865997314453,
        2.8085227012634277,
        63.75581359863281,
        26.74026107788086,
        46.38532257080078,
        0.24645215272903442,
        15.20484447479248,
        1.8522661924362183,
        52.44639587402344,
        41.203941345214844,
        0.5699318051338196,
        0.4306056499481201,
        0.0079906200783,
        0.9056035876274109,
        0.016054258216,
        0.7509759664535522,
        0.088582855955,
        48.27726364135742,
        39.38459014892578,
        1.1875036682231526,
        107.75627136230469,
        29.772504806518555,
        88.32481384277344,
        23.039295196533203,
        62.85848617553711,
        5.546318531036377,
        11.92244815826416,
        5.555448055267334,
        0.9218489527702332,
        94.59927368164063,
        77.29698181152344,
        63.05263137817383,
        53.979705810546875,
        24.643558502197266,
        61.30192184448242,
        22.209999084472656,
        40.06374740600586,
        38.100341796875,
        46.534156799316406,
        47.445735931396484,
        41.105812072753906,
        18.11349105834961,
        50.0,
    ],
    manipulated: [
        63.05263137817383,
        53.979705810546875,
        24.643558502197266,
        61.30192184448242,
        22.209999084472656,
        40.06374740600586,
        38.100341796875,
        46.534156799316406,
        47.445735931396484,
        41.105812072753906,
        18.11349105834961,
        50.0,
    ],
    disturbances: [0; 20],
    rng: 4651207995.0,
};

/// Every corpus tuple still agrees with the Fortran.
#[test]
fn the_corpus_agrees_with_the_fortran() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    for case in corpus() {
        let outcome = check(
            &mut oracle,
            &base,
            &case.tuple,
            TOLERANCE,
            Mutation::None,
            0,
            None,
        );
        match outcome {
            Outcome::Agreed { worst, component } => {
                println!("{}: worst {worst:.3e} at YP({component})", case.name);
            }
            Outcome::Disagreed(finding) => {
                panic!("{} regressed: {finding}. {}", case.name, case.why)
            }
            other => panic!(
                "{} no longer reaches a comparison ({other:?}), so it is not \
                 checking anything. {}",
                case.name, case.why
            ),
        }
    }
}

/// Every corpus tuple detects the fault it was recorded for.
///
/// Without this the corpus is a list of inputs that pass, which is what an
/// empty list also is. An entry has to be able to fail.
#[test]
fn every_corpus_tuple_can_still_fail() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);
    for case in corpus() {
        let outcome = check(
            &mut oracle,
            &base,
            &case.tuple,
            TOLERANCE,
            case.detects,
            0,
            None,
        );
        assert!(
            matches!(outcome, Outcome::Disagreed(_)),
            "{} does not detect the fault it was recorded for ({outcome:?}), \
             so it is an inert entry in the corpus. {}",
            case.name,
            case.why
        );
        println!("{}: detects its fault", case.name);
    }
}

// ---------------------------------------------------------------------------
// The open finding
// ---------------------------------------------------------------------------

/// What the five-million-tuple run measured on the one counterexample.
///
/// Recorded so a regression is visible rather than inferred, and so the next
/// session compares against a number instead of a verdict.
const OPEN_FINDING_ERROR: f64 = 4.606_900_216_695_249e-12;

/// The counterexample the five-million-tuple run found, shrunk to three knobs.
///
/// Seed `0x7E2_0062`, tuple 863,105, shrunk from 22 differing knobs to 3. The
/// fifty states are *exactly* nominal; what is left is `IDV(13)`, one generator
/// word, and a time of 973 hours.
///
/// `IDV(13)` is the slow drift in the reaction kinetics. Nine hundred and
/// seventy-three hours of that walk, from this generator word, leaves
/// `TESUB8(7)` at -3.2e7 and `TESUB8(8)` at +5.4e7, so the two reaction rates
/// come out at -8.07e9 and +1.11e10 against a nominal of order one. The plant
/// is nowhere near physical; it is also nowhere near tripping, because none of
/// the eight shutdown conditions looks at a reaction rate.
///
/// See `the_open_finding_is_the_vendored_libm_and_not_the_algebra` for the
/// attribution and the numbers.
const OPEN_FINDING: Tuple = Tuple {
    time: 973.308134169278,
    state: [
        10.404913902282715,
        4.3639960289001465,
        7.570059776306152,
        0.42300423979759216,
        24.155134201049805,
        2.9425976276397705,
        154.37705993652344,
        159.1865997314453,
        2.8085227012634277,
        63.75581359863281,
        26.74026107788086,
        46.38532257080078,
        0.24645215272903442,
        15.20484447479248,
        1.8522661924362183,
        52.44639587402344,
        41.203941345214844,
        0.5699318051338196,
        0.4306056499481201,
        0.0079906200783,
        0.9056035876274109,
        0.016054258216,
        0.7509759664535522,
        0.088582855955,
        48.27726364135742,
        39.38459014892578,
        0.3755297362804413,
        107.75627136230469,
        29.772504806518555,
        88.32481384277344,
        23.039295196533203,
        62.85848617553711,
        5.546318531036377,
        11.92244815826416,
        5.555448055267334,
        0.9218489527702332,
        94.59927368164063,
        77.29698181152344,
        63.05263137817383,
        53.979705810546875,
        24.643558502197266,
        61.30192184448242,
        22.209999084472656,
        40.06374740600586,
        38.100341796875,
        46.534156799316406,
        47.445735931396484,
        41.105812072753906,
        18.11349105834961,
        50.0,
    ],
    manipulated: [
        63.05263137817383,
        53.979705810546875,
        24.643558502197266,
        61.30192184448242,
        22.209999084472656,
        40.06374740600586,
        38.100341796875,
        46.534156799316406,
        47.445735931396484,
        41.105812072753906,
        18.11349105834961,
        50.0,
    ],
    disturbances: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    rng: 985491900.0,
};

/// The one counterexample five million tuples produced, and where it comes
/// from.
///
/// **This tuple does not meet the Tier 8 gate under the vendored libm**, and
/// nothing here is written to make it. `YP(9)`, the reactor energy balance,
/// comes out 4.61e-12 of the scale of its terms away from the Fortran, against
/// a gate of 1e-12. That is a finding and it is recorded as one.
///
/// What the finding *is* is settled here rather than left open, because the
/// evidence is one feature flag away. Under `--features oracle,libm-system`,
/// where `exp` and `pow` are the ones gfortran itself calls, all fifty
/// components are **bit-identical**. The algebra between `teprob.f:407` and
/// `811` is therefore exact, and the entire difference is the vendored `exp`.
///
/// Why it shows up here and nowhere in Tier 2: the reaction rates on this state
/// are 1.1e10, because `IDV(13)`'s kinetic drift has run for 973 hours. `exp`
/// is accurate to a relative ULP, so a 1.1e-16 relative difference on a rate of
/// 1.1e10 is 1.2e-6 absolute, and the reactor energy balance's own terms are
/// only 5.2e4. Amplification by fourteen orders of magnitude turns the one ULP
/// the project has always accepted into 4.6e-12 of scale. No state in any of
/// the three Tier 2 pools has a drift factor above about one, which is exactly
/// why an unguided search was worth building.
///
/// This test asserts the two halves of that, so it fails if either changes: bit
/// equality under the platform libm, and a bounded, attributed error under the
/// vendored one.
#[test]
fn the_open_finding_is_the_vendored_libm_and_not_the_algebra() {
    let mut oracle = Oracle::lock();
    let base = base(&mut oracle);

    // The drift is the whole explanation, so it is measured rather than
    // asserted from the prose.
    let scenario = OPEN_FINDING.scenario(&base);
    let snapshot = scenario.force(&mut oracle);
    let drift = (
        oracle.tesub8(7, scenario.time),
        oracle.tesub8(8, scenario.time),
    );
    println!(
        "IDV(13) drift after {:.1} h: TESUB8(7) = {:.4e}, TESUB8(8) = {:.4e}",
        scenario.time, drift.0, drift.1
    );
    println!("reaction rates: {:?}", snapshot.common.rr);
    assert!(
        drift.0.abs() > 1e6 && drift.1.abs() > 1e6,
        "the drift is no longer extreme, so this tuple is no longer the case \
         the finding is about and the attribution below does not apply"
    );

    let outcome = check(
        &mut oracle,
        &base,
        &OPEN_FINDING,
        TOLERANCE,
        Mutation::None,
        0,
        None,
    );
    println!("{outcome:?}");

    if math::USES_SYSTEM_LIBM {
        // The strong claim, and the one that settles the attribution.
        assert!(
            matches!(outcome, Outcome::Agreed { worst: 0.0, .. }),
            "the port is no longer bit-identical to the Fortran on this state \
             under the platform libm, so the difference has moved into the \
             algebra and this is no longer a libm finding"
        );
        return;
    }

    let Outcome::Disagreed(finding) = outcome else {
        panic!(
            "this tuple now agrees under the vendored libm ({outcome:?}). That \
             is better than the recorded state of the world, not worse, but it \
             means the recorded number is stale and the finding needs \
             re-measuring."
        );
    };
    assert_eq!(finding.component, 9, "the finding moved component");
    println!(
        "YP(9): {:.4e} of scale, recorded {OPEN_FINDING_ERROR:.4e}",
        finding.error
    );
    // A bound, not an acceptance. The gate is still 1e-12; this pins how far
    // past it the one known state sits, so a real regression cannot hide
    // behind a known divergence.
    assert!(
        finding.error <= 2.0 * OPEN_FINDING_ERROR,
        "the known libm divergence on this state grew from \
         {OPEN_FINDING_ERROR:e} to {:e} of scale",
        finding.error
    );
}
