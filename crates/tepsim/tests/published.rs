//! Every constant in `tepsim::published`, read back out of the vendored source.
//!
//! B-0064. The project's rule is that a constant is asserted and never
//! retyped, because a digit read off a listing is a silent failure with no
//! test that catches it. Fifty-two of the numbers in that module are seeds
//! from a *comment block*, which no `COMMON` block can be interrogated for, so
//! the assertion is a parse of `reference/fortran/teprob.f` itself.

// Exact comparisons throughout: a transcribed constant either is the number in
// the vendored source or it is not, and a published value either survives
// five-digit rounding unchanged or it does not. A tolerance here would defeat
// the point of the file.
#![allow(
    clippy::float_cmp,
    reason = "exact equality with the vendored source is the property under test"
)]

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use tepsim::published::{
    self, COLUMNS, COMPILED_IN_SEED, D00_TRAINING_ALTERNATIVE, D00_TRAINING_ROWS,
    D18_TRAINING_ALTERNATIVE, DRIVER_STEPS, FILES, ORIGINAL_SEED, PUBLISHED_DIGITS, SAMPLE_EVERY,
    SETTLING_STEPS, Split, TESTING_ROWS, TRAINING_ROWS, Unavailable,
};

fn reference(rest: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../reference")
        .join(rest)
}

/// Every `C   dNN_tr: G=1234.D0` line, keyed by the label.
///
/// `d19_te` is written `.DO` rather than `.D0` in the source, an `O` for a
/// zero. It is a comment, so the compiler never saw it and nothing ever
/// complained. The parser accepts both endings on purpose; tightening it to
/// `D0` would silently drop that one seed.
fn seed_comments() -> BTreeMap<String, f64> {
    let text = fs::read_to_string(reference("fortran/teprob.f")).expect("teprob.f");
    let mut seeds = BTreeMap::new();
    for line in text.lines() {
        let Some(rest) = line.strip_prefix('C') else {
            continue;
        };
        let Some((label, value)) = rest.split_once(": G=") else {
            continue;
        };
        let value = value.trim();
        let Some(digits) = value
            .strip_suffix(".D0")
            .or_else(|| value.strip_suffix(".DO"))
        else {
            continue;
        };
        let seed: f64 = digits.trim().parse().expect("a seed is a number");
        seeds.insert(label.trim().to_string(), seed);
    }
    seeds
}

#[test]
fn the_parser_finds_the_whole_seed_block() {
    let seeds = seed_comments();
    // 27 training plus 27 testing plus `d00_tr_new`, `original` and `dd18_tr`.
    assert_eq!(
        seeds.len(),
        57,
        "found {} labelled seeds, so the parser is reading the wrong thing and \
         every assertion below is vacuous",
        seeds.len()
    );
    assert!(seeds.contains_key("d19_te"), "the `.DO` line was dropped");
}

#[test]
fn every_files_seed_is_the_one_the_source_records() {
    let seeds = seed_comments();
    for file in &FILES {
        let label = format!(
            "d{:02}_{}",
            file.fault,
            match file.split {
                Split::Training => "tr",
                Split::Testing => "te",
            }
        );
        let recorded = seeds
            .get(&label)
            .unwrap_or_else(|| panic!("teprob.f has no seed labelled {label}"));
        assert_eq!(
            file.seed, *recorded,
            "{label}: the module says {} and teprob.f says {recorded}",
            file.seed
        );
    }
}

#[test]
fn the_ambiguous_and_compiled_in_seeds_match_too() {
    let seeds = seed_comments();
    assert_eq!(seeds["d00_tr_new"], D00_TRAINING_ALTERNATIVE);
    assert_eq!(seeds["dd18_tr"], D18_TRAINING_ALTERNATIVE);
    assert_eq!(seeds["original"], ORIGINAL_SEED);

    // The compiled-in one is code, not a comment, so it is matched separately.
    let text = fs::read_to_string(reference("fortran/teprob.f")).expect("teprob.f");
    assert!(
        text.contains(&format!("G={COMPILED_IN_SEED:.0}.D0")),
        "teprob.f does not assign G={COMPILED_IN_SEED:.0}"
    );
}

#[test]
fn the_run_geometry_is_the_drivers() {
    let text = fs::read_to_string(reference("fortran/temain_mod.f")).expect("temain_mod.f");
    let stripped: String = text
        .lines()
        .filter(|l| !l.starts_with('C') && !l.starts_with('c'))
        .collect::<Vec<_>>()
        .join("\n");

    assert!(
        stripped.contains(&format!("NPTS = {DRIVER_STEPS}")),
        "temain_mod.f does not set NPTS = {DRIVER_STEPS}"
    );
    assert!(
        stripped.contains(&format!("SSPTS = 3600 * {}", SETTLING_STEPS / 3600)),
        "temain_mod.f does not set SSPTS = 3600 * {}",
        SETTLING_STEPS / 3600
    );
    assert!(
        stripped.contains(&format!("MOD(I,{SAMPLE_EVERY})")),
        "temain_mod.f does not write every {SAMPLE_EVERY} steps"
    );
    assert!(
        stripped.contains(&format!("FORMAT(1X,E13.{PUBLISHED_DIGITS})")),
        "temain_mod.f does not write E13.{PUBLISHED_DIGITS}"
    );
}

/// The shipped files, read for their actual shape.
fn shape(stem: &str) -> (usize, usize) {
    let text = fs::read_to_string(reference(&format!("data/{stem}.dat")))
        .unwrap_or_else(|_| panic!("{stem}.dat"));
    let rows: Vec<usize> = text
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().count())
        .collect();
    let width = rows.first().copied().unwrap_or(0);
    assert!(
        rows.iter().all(|w| *w == width),
        "{stem}.dat is ragged, so its shape is not a number"
    );
    (rows.len(), width)
}

#[test]
fn every_files_row_count_is_the_shipped_files_row_count() {
    for file in &FILES {
        let stem = file.stem();
        // `d00.dat` alone is stored transposed: 52 rows of 500 columns.
        let (rows, width) = shape(&stem);
        let (rows, width) = if stem == "d00" {
            (width, rows)
        } else {
            (rows, width)
        };

        assert_eq!(
            rows, file.rows,
            "{stem}.dat has {rows} rows and the module says {}",
            file.rows
        );
        assert_eq!(
            width, COLUMNS,
            "{stem}.dat has {width} columns and the module says {COLUMNS}"
        );
    }
}

#[test]
fn the_row_counts_are_the_ones_named() {
    assert_eq!(TESTING_ROWS, 960);
    assert_eq!(TRAINING_ROWS, 480);
    assert_eq!(D00_TRAINING_ROWS, 500);
    assert_eq!(FILES.len(), 44);
}

/// Every value in a published file is already five significant digits, so
/// rounding it again changes nothing.
///
/// This is what makes [`published::round_as_published`] a reproduction of the
/// original's output stage rather than a guess about it. If the rounding were
/// wrong in the last digit, or applied to the wrong number of digits, this
/// would fail on the first row.
#[test]
fn published_values_are_already_rounded_to_five_digits() {
    let mut checked = 0_usize;
    let mut idempotent = 0_usize;
    for stem in ["d00_te", "d01_te", "d06", "d14_te", "d21"] {
        let text = fs::read_to_string(reference(&format!("data/{stem}.dat"))).expect(stem);
        for field in text.split_whitespace() {
            let value: f64 = field.parse().expect("a published field is a number");
            checked += 1;
            if published::round_as_published(value) == value {
                idempotent += 1;
            }
        }
    }
    assert!(
        checked > 100_000,
        "only {checked} fields, too few to mean much"
    );
    assert_eq!(
        idempotent,
        checked,
        "{} of {checked} published values carry more than {PUBLISHED_DIGITS} \
         significant digits",
        checked - idempotent
    );
}

/// And the same check has teeth: a sixth digit does not survive.
#[test]
fn rounding_actually_discards_a_sixth_digit() {
    assert_eq!(published::round_as_published(1.234_567_8), 1.2346);
    assert_eq!(published::round_as_published(-9.876_543_21e-3), -9.8765e-3);
    assert_eq!(published::round_as_published(0.0), 0.0);
    assert!(published::round_as_published(f64::NAN).is_nan());
}

#[test]
fn d21_is_not_available_and_everything_else_is() {
    for file in &FILES {
        let stem = file.stem();
        match file.scenario() {
            Ok(scenario) => {
                assert_ne!(file.fault, 21, "{stem} built a scenario for IDV(21)");
                assert_eq!(scenario.seed, file.seed);
                assert_eq!(scenario.sample_every, SAMPLE_EVERY);
                assert!(
                    !scenario.driver_forces_idv12,
                    "{stem}: the published bytes say IDV(12) was not forced"
                );
                let kept = scenario.samples() - file.discarded_rows();
                assert_eq!(
                    kept,
                    file.rows,
                    "{stem}: the scenario yields {kept} rows after discarding {}, \
                     and the file has {}",
                    file.discarded_rows(),
                    file.rows
                );
            }
            Err(Unavailable::FaultNotInThisRevision { fault }) => {
                assert_eq!(fault, 21, "{stem} was refused for IDV({fault})");
            }
        }
    }
}

#[test]
fn the_published_column_order_is_measurements_then_valves() {
    let names = published::column_names();
    assert_eq!(names.len(), COLUMNS);
    let all = tepsim::channel_names();
    assert_eq!(names[0], all[0], "the first column is XMEAS(1)");
    assert_eq!(names[40], all[40], "the forty-first is XMEAS(41)");
    assert_eq!(names[41], all[41], "then the valves start");
    assert_eq!(names[51], all[51], "and stop at XMV(11)");
    assert!(
        names.iter().filter(|n| n.starts_with("XMEAS")).count() == 41,
        "41 measurements, got {:?}",
        names.iter().filter(|n| n.starts_with("XMEAS")).count()
    );
    // The agitator is the last channel a `Run` carries and is in no published
    // file, so it must be exactly the one that got dropped.
    assert!(
        !names.contains(all.last().expect("a last channel")),
        "the agitator {} is in no published file",
        all.last().expect("a last channel")
    );
}

// ---------------------------------------------------------------------------
// The Rieth ensemble
// ---------------------------------------------------------------------------

use tepsim::published::rieth;

/// The geometry is the one the dataset documents, and the arithmetic closes.
#[test]
fn the_rieth_geometry_is_self_consistent() {
    // 25 h and 48 h at three minutes.
    assert_eq!(
        rieth::TRAINING_ROWS,
        (rieth::TRAINING_HOURS * 3600.0 / (SAMPLE_EVERY as f64)) as usize
    );
    assert_eq!(
        rieth::TESTING_ROWS,
        (rieth::TESTING_HOURS * 3600.0 / (SAMPLE_EVERY as f64)) as usize
    );
    // The onset row is where the documented normal/faulty split falls.
    assert_eq!(
        rieth::TRAINING_NORMAL_ROWS,
        (rieth::TRAINING_ONSET_HOURS * 3600.0 / (SAMPLE_EVERY as f64)) as usize
    );
    assert_eq!(
        rieth::TESTING_NORMAL_ROWS,
        (rieth::TESTING_ONSET_HOURS * 3600.0 / (SAMPLE_EVERY as f64)) as usize
    );
    // 20 + 480 and 160 + 800, as the dataset's description states.
    assert_eq!(rieth::TRAINING_ROWS - rieth::TRAINING_NORMAL_ROWS, 480);
    assert_eq!(rieth::TESTING_ROWS - rieth::TESTING_NORMAL_ROWS, 800);
}

/// The training protocol is *not* the original's, and the test says so out
/// loud, because both are 25 hours and both split 20 against 480.
#[test]
fn the_rieth_training_run_is_not_the_original_training_run() {
    let original = FILES
        .iter()
        .find(|f| f.fault == 1 && f.split == Split::Training)
        .expect("d01")
        .scenario()
        .expect("a scenario");
    let theirs = rieth::scenario(1, Split::Training, 0).expect("a scenario");

    assert_eq!(original.hours, theirs.hours, "both are 25 hours");
    // Ours has the fault from step zero and discards the first hour; theirs
    // schedules it at hour one and keeps every row.
    assert!(original.schedule.is_empty(), "the original has no schedule");
    assert!(
        !theirs.schedule.is_empty(),
        "Rieth's fault arrives on a schedule"
    );
    assert_eq!(original.samples(), rieth::TRAINING_ROWS);
    assert_eq!(theirs.samples(), rieth::TRAINING_ROWS);
    // The rows kept differ: 480 against 500.
    assert_eq!(
        original.samples()
            - FILES
                .iter()
                .find(|f| f.fault == 1 && f.split == Split::Training)
                .expect("d01")
                .discarded_rows(),
        480
    );
}

/// Every run in the ensemble has its own seed, and the same coordinates always
/// give the same one.
#[test]
fn every_run_has_a_distinct_and_reproducible_seed() {
    use std::collections::BTreeSet;
    let mut seen = BTreeSet::new();
    for fault in 0..=20 {
        for split in [Split::Training, Split::Testing] {
            for run in 0..rieth::RUNS {
                let seed = rieth::seed(fault, split, run);
                assert!(seed > 0.0 && seed.is_finite(), "seed {seed}");
                assert!(
                    seen.insert(seed.to_bits()),
                    "seed {seed} collides at fault {fault} run {run}"
                );
                assert_eq!(seed, rieth::seed(fault, split, run), "not reproducible");
            }
        }
    }
    assert_eq!(seen.len(), 21 * 2 * rieth::RUNS);
}

/// Two runs of the ensemble really are different runs, and the fault really
/// does arrive where the geometry says.
#[test]
fn the_ensemble_produces_distinct_runs_with_the_fault_where_it_belongs() {
    let a = tepsim::Simulation::new(
        rieth::scenario(1, Split::Training, 0)
            .expect("a scenario")
            .with_hours(3.0),
    )
    .run();
    let b = tepsim::Simulation::new(
        rieth::scenario(1, Split::Training, 1)
            .expect("a scenario")
            .with_hours(3.0),
    )
    .run();
    assert_ne!(
        a.samples.first().map(tepsim::Sample::row),
        b.samples.first().map(tepsim::Sample::row),
        "two runs of the ensemble are identical, so the seeds did not take"
    );

    // Fault-free before the onset, faulted after it.
    let onset = rieth::TRAINING_NORMAL_ROWS;
    assert!(
        !a.samples[onset - 1].labels.faulted(),
        "the fault arrived before its onset row"
    );
    assert!(
        a.samples[onset + 1].labels.faulted(),
        "the fault had not arrived after its onset row"
    );
}

#[test]
fn the_ensemble_refuses_idv21_like_the_rest_of_the_module() {
    assert!(matches!(
        rieth::scenario(21, Split::Testing, 0),
        Err(Unavailable::FaultNotInThisRevision { fault: 21 })
    ));
}
