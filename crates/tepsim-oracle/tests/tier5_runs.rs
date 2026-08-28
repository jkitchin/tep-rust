//! The Tier 5 run harness: does it produce comparable data from both sources?
//!
//! B-0047a. Nothing here judges the plant. These check that a *run* is a
//! well-defined object: reproducible from its scenario and seed, the same
//! shape from either source, and actually different when the seed is.
//!
//! The last one matters more than it looks. A harness whose seed did nothing
//! would produce a hundred identical runs, every statistic would have zero
//! variance, and every equivalence test in B-0047b would pass trivially.

// Differential tests: exact comparisons are the property under test, and
// arithmetic is transcribed from `teprob.f` so that a reader can check it
// against the listing line by line. Rearranging either would defeat the point.
#![allow(
    clippy::float_cmp,
    reason = "bit equality against the Fortran is the property under test"
)]
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions are transcribed to be checkable against teprob.f"
)]
#![cfg(feature = "oracle")]

use std::collections::HashSet;

use tepsim_core::math;
use tepsim_oracle::Oracle;
use tepsim_oracle::tier5::{
    Battery, FAULTS, SAMPLE_EVERY, SCENARIOS, Scenario, VARIABLES, run_fortran, run_port, seed,
    start,
};

/// Long enough to cross a few analyser samples and a composition-loop fire,
/// short enough to run in a test.
const HOURS: usize = 2;

// ---------------------------------------------------------------------------
// Seeds
// ---------------------------------------------------------------------------

/// `TESUB7` is `G <- DMOD(G * 9228907, 2^32)`. Zero is a fixed point, and a
/// generator modulo a power of two keeps whatever factors of two its seed has,
/// so an even seed has a shorter period and low bits that never move.
///
/// The seeds this harness generates are therefore odd. They are *not* required
/// to be below 2^32: the first `DMOD` reduces them, and the original's own
/// published seeds routinely exceed it. Index 0 is the seed compiled into
/// `teprob.f:1187`, which is 4,651,207,995 and so above the modulus itself.
#[test]
fn every_seed_is_a_valid_generator_word() {
    let mut seen = HashSet::new();
    for index in 0..100 {
        let g = seed(index);
        assert!(g > 0.0, "seed {index} is {g}");
        assert_eq!(g.fract(), 0.0, "seed {index} is not an integer: {g}");
        assert!(seen.insert(g as u64), "seed {index} = {g} is a duplicate");
        if index > 0 {
            assert_eq!(
                (g as u64) % 2,
                1,
                "seed {index} is {g}, which is even; the generator's low bits                  would never move"
            );
            assert!(g < 4_294_967_296.0, "seed {index} is {g}");
        }
    }
    // Index 0 is the golden seed every other test in this repository uses.
    assert_eq!(seed(0), tepsim_oracle::golden::SEED);
    assert!(
        seed(0) > 4_294_967_296.0,
        "the compiled-in seed no longer exceeds the modulus, so the note about          published seeds doing so needs revisiting"
    );
}

/// The seed table `teprob.f` carries in comments, which Tier 7 will need.
///
/// Transcribed rather than derived, so this checks the transcription: the
/// counts, that none is zero, and the two properties the module documentation
/// claims about them.
#[test]
fn the_published_seed_table_is_transcribed_correctly() {
    use tepsim_oracle::tier5::published_seeds as p;

    assert_eq!(p::COMPILED_IN, tepsim_oracle::golden::SEED);
    assert_eq!(p::TRAINING.len(), 27, "d00_tr through d26_tr");
    assert_eq!(p::TESTING.len(), 27, "d00_te through d26_te");

    for (label, table) in [("training", &p::TRAINING), ("testing", &p::TESTING)] {
        for (index, g) in table.iter().enumerate() {
            assert!(*g > 0.0, "{label} seed {index} is {g}");
            assert_eq!(g.fract(), 0.0, "{label} seed {index} is not an integer");
        }
    }

    // The two claims the docs make about the original's seeds.
    let above_modulus = p::TRAINING
        .iter()
        .chain(&p::TESTING)
        .filter(|g| **g > 4_294_967_296.0)
        .count();
    let even = p::TRAINING
        .iter()
        .chain(&p::TESTING)
        .filter(|g| (**g as u64) % 2 == 0)
        .count();
    println!("published seeds: {above_modulus} of 54 exceed 2^32, {even} of 54 are          even");
    assert!(
        above_modulus > 0,
        "no published seed exceeds the modulus, so the note saying they do is          wrong"
    );
    assert!(
        even > 0,
        "no published seed is even, so the note saying some are is wrong"
    );

    // `d19_te` is written `9090909232.DO` in the source, with the letter O.
    assert_eq!(p::TESTING[p::MALFORMED_TESTING_INDEX], 9_090_909_232.0);
}

// ---------------------------------------------------------------------------
// Scenarios
// ---------------------------------------------------------------------------

/// The scenario count comes from the vendored Fortran, not from the
/// literature.
///
/// `teprob.f:340` is `DO 500 I=1,20`, and `tepsim_core::FAULTS` has twenty
/// entries. `PLAN.org` says twenty-two scenarios, which counts an `IDV(21)`
/// that later versions of the model have and this source does not.
#[test]
fn there_are_twenty_faults_and_twenty_one_scenarios() {
    assert_eq!(FAULTS, tepsim_core::FAULTS.len());
    assert_eq!(SCENARIOS, FAULTS + 1);
    assert_eq!(Scenario::all().count(), SCENARIOS);

    let nominal = Scenario::NOMINAL;
    assert_eq!(nominal.disturbances(), [0.0; 20]);
    assert_eq!(nominal.label(), "nominal");

    for n in 1..=FAULTS {
        let s = Scenario::fault(n);
        let idv = s.disturbances();
        assert_eq!(idv[n - 1], 1.0, "IDV({n}) is not set");
        assert_eq!(idv.iter().sum::<f64>(), 1.0, "IDV({n}) set more than one");
        assert_eq!(s.label(), format!("IDV({n})"));
    }
}

// ---------------------------------------------------------------------------
// Shape
// ---------------------------------------------------------------------------

#[test]
fn both_sources_produce_the_same_shape() {
    let mut oracle = Oracle::lock();
    let start = start(&mut oracle);

    let ours = run_port(&start, Scenario::NOMINAL, seed(0), HOURS);
    let theirs = run_fortran(&mut oracle, &start, Scenario::NOMINAL, seed(0), HOURS);

    let expected = HOURS * 3_600 / SAMPLE_EVERY;
    assert_eq!(ours.samples.len(), expected, "port sample count");
    assert_eq!(theirs.samples.len(), expected, "Fortran sample count");
    assert_eq!(ours.series(0).len(), expected);
    assert_eq!(ours.all_series().len(), VARIABLES);
    println!("{HOURS} h gives {expected} samples of {VARIABLES} variables");
}

// ---------------------------------------------------------------------------
// Reproducibility
// ---------------------------------------------------------------------------

/// The same scenario and seed give bit-identical runs, every time.
///
/// This is what lets B-0047b compare a hundred runs against a hundred others
/// and attribute every difference to the source rather than to the harness.
#[test]
fn a_run_is_reproducible_from_its_scenario_and_seed() {
    let mut oracle = Oracle::lock();
    let start = start(&mut oracle);

    let first = run_port(&start, Scenario::fault(4), seed(7), HOURS);
    let second = run_port(&start, Scenario::fault(4), seed(7), HOURS);
    assert_eq!(first, second, "two port runs of the same scenario differ");

    // And on the Fortran side, where `COMMON` survives between runs and the
    // warm start is not reset by `TEINIT` (see `Oracle::init_cold`).
    let theirs = run_fortran(&mut oracle, &start, Scenario::fault(4), seed(7), HOURS);
    // Deliberately run something else in between, to disturb `COMMON`.
    let _ = run_fortran(&mut oracle, &start, Scenario::fault(13), seed(3), HOURS);
    let again = run_fortran(&mut oracle, &start, Scenario::fault(4), seed(7), HOURS);
    assert_eq!(
        theirs, again,
        "two Fortran runs of the same scenario differ, so COMMON is leaking \
         between runs and the whole battery is order-dependent"
    );
}

/// Different seeds give different runs.
///
/// A harness whose seed did nothing would produce a hundred identical runs,
/// every statistic would have zero variance, and every equivalence test in
/// B-0047b would pass trivially and mean nothing.
#[test]
fn different_seeds_give_different_runs() {
    let mut oracle = Oracle::lock();
    let start = start(&mut oracle);

    let runs: Vec<_> = (0..5)
        .map(|i| run_port(&start, Scenario::NOMINAL, seed(i), HOURS))
        .collect();

    let mut spread = 0.0_f64;
    for i in 0..runs.len() {
        for j in (i + 1)..runs.len() {
            assert_ne!(runs[i].samples, runs[j].samples, "seeds {i} and {j}");
            // And the difference is a real one, not a last-bit one: the
            // disturbance walks are different realisations.
            let worst = runs[i]
                .samples
                .iter()
                .zip(&runs[j].samples)
                .flat_map(|(a, b)| a.iter().zip(b))
                .filter(|(_, y)| **y != 0.0)
                .map(|(x, y)| (x - y).abs() / y.abs())
                .fold(0.0_f64, f64::max);
            spread = spread.max(worst);
        }
    }
    println!("five nominal seeds: worst relative spread {spread:.3e}");
    assert!(
        spread > 1e-4,
        "the largest difference between two seeds is {spread:.3e}, which is \
         rounding rather than a different disturbance realisation"
    );
}

// ---------------------------------------------------------------------------
// The two sources agree
// ---------------------------------------------------------------------------

/// Under `libm-system` the port and the Fortran are the same run, bit for bit.
///
/// This is the harness's real correctness check: it says the two sides are
/// driven identically, seeded identically and sampled identically, which is
/// the whole job of B-0047a. Under the vendored libm the transcendentals
/// differ and only the magnitude of the divergence is recorded.
#[test]
fn the_two_sources_agree_over_a_whole_run() {
    let mut oracle = Oracle::lock();
    let start = start(&mut oracle);

    for scenario in [Scenario::NOMINAL, Scenario::fault(1), Scenario::fault(13)] {
        let ours = run_port(&start, scenario, seed(0), HOURS);
        let theirs = run_fortran(&mut oracle, &start, scenario, seed(0), HOURS);

        let mut worst = (0.0_f64, 0usize, 0usize);
        for (sample, (a, b)) in ours.samples.iter().zip(&theirs.samples).enumerate() {
            for (variable, (x, y)) in a.iter().zip(b).enumerate() {
                if *y == 0.0 {
                    continue;
                }
                let relative = (x - y).abs() / y.abs();
                if relative > worst.0 {
                    worst = (relative, sample, variable);
                }
            }
        }
        println!(
            "{}: worst {:.3e} at sample {} variable {} ({} libm)",
            scenario.label(),
            worst.0,
            worst.1,
            worst.2 + 1,
            if math::USES_SYSTEM_LIBM {
                "platform"
            } else {
                "vendored"
            }
        );
        assert_eq!(ours.tripped, theirs.tripped, "{}", scenario.label());

        if math::USES_SYSTEM_LIBM {
            assert_eq!(
                worst.0,
                0.0,
                "{}: the harness does not drive the two sources identically",
                scenario.label()
            );
        } else {
            assert!(
                worst.0 < 1e-6,
                "{}: {:.3e} is past amplification of a one-ULP exp difference",
                scenario.label(),
                worst.0
            );
        }
    }
}

/// A disturbance actually changes the plant.
///
/// Twenty scenarios that all produced the nominal trajectory would make the
/// battery a hundred-fold repetition of one comparison.
#[test]
fn every_disturbance_moves_the_plant() {
    let mut oracle = Oracle::lock();
    let start = start(&mut oracle);
    let nominal = run_port(&start, Scenario::NOMINAL, seed(0), HOURS);

    let mut quiet = Vec::new();
    for n in 1..=FAULTS {
        let scenario = Scenario::fault(n);
        let faulted = run_port(&start, scenario, seed(0), HOURS);
        let worst = nominal
            .samples
            .iter()
            .zip(&faulted.samples)
            .flat_map(|(a, b)| a.iter().zip(b))
            .filter(|(x, _)| **x != 0.0)
            .map(|(x, y)| (x - y).abs() / x.abs())
            .fold(0.0_f64, f64::max);
        println!(
            "  {:<8} worst departure from nominal {:.3e}{}",
            scenario.label(),
            worst,
            match faulted.tripped {
                Some(step) => format!("  TRIPPED at {step}"),
                None => String::new(),
            }
        );
        if worst < 1e-6 {
            quiet.push(n);
        }
    }
    assert!(
        quiet.is_empty(),
        "IDV {quiet:?} left the plant within 1e-6 of nominal over {HOURS} h. \
         Either the disturbance is not wired up, or it needs longer than \
         {HOURS} h to show and the battery's horizon has to say so."
    );
}

// ---------------------------------------------------------------------------
// The battery's own shape
// ---------------------------------------------------------------------------

#[test]
fn the_battery_sizes_are_what_they_claim() {
    let smoke = Battery::SMOKE;
    let full = Battery::FULL;

    assert_eq!(full.scenarios, SCENARIOS);
    assert_eq!(full.samples(), 960, "48 h at 180-step sampling");
    assert_eq!(full.runs(), 2_100);
    assert!(
        smoke.runs() < full.runs() / 100,
        "the smoke battery is not small"
    );
    assert!(smoke.scenarios <= SCENARIOS);

    // 766 ms per 48-hour port run, measured in release.
    println!(
        "full battery: {} runs per source, {} samples each, about {:.0} min \
         per source at 766 ms per 48 h run",
        full.runs(),
        full.samples(),
        full.runs() as f64 * 0.766 * (full.hours as f64 / 48.0) / 60.0
    );

    // Selected from the environment, defaulting to smoke.
    assert_eq!(Battery::selected(), smoke, "TEP_TIER5 should be unset here");
}
