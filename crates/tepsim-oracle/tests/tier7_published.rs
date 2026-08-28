//! Tier 7: can the port regenerate the published `d00`-`d21` datasets?
//!
//! B-0051. The answer this file arrives at is *no, and here is exactly why
//! not*, which is what `PLAN.org` asks for:
//!
//! > Where reproduction is imperfect, document exactly which protocol detail
//! > is unknown rather than tuning until the numbers match.
//!
//! Accordingly, **nothing in this file or in [`tepsim_oracle::tier7`] is a
//! knob.** Every protocol quantity is read off `temain_mod.f` or derived from
//! a row count, and which of the two it is is recorded in
//! [`tier7::Provenance`]. The findings are printed rather than gated, because
//! a gate would have to be a number nobody measured.
//!
//! # How agreement is judged without inventing a threshold
//!
//! A raw Kolmogorov-Smirnov statistic of 0.06 between a published file and a
//! port run is not interpretable on its own. Tier 5 solves this by calibrating
//! against the Fortran's own run-to-run spread; Tier 7 cannot, because there
//! is exactly one published file per scenario and no way to make another.
//!
//! So the null is built the other way round. For each file the port is run at
//! the seed `teprob.f` records for it, and also at several *decoy* seeds:
//! wrong words, right protocol. Comparing the published file against each of
//! those gives a distribution of "same simulator, same protocol, unrelated
//! noise realisation". The recorded seed's value is then one more draw from
//! that distribution, and the question becomes a rank:
//!
//! > Does knowing the published seed make the port's output closer to the
//! > published file than an arbitrary seed would?
//!
//! If reproduction worked, the recorded seed would rank first on every
//! statistic by a wide margin. If the protocol or the revision is wrong, it
//! ranks where chance puts it. No tolerance appears anywhere in that sentence,
//! which is why it is the test.
//!
//! # What was found
//!
//! Trajectory reproduction fails completely and distributional agreement is
//! good, and the two are separate results.
//!
//! - **The recorded seeds carry no information.** They rank first on 10 of 42
//!   files against 8.4 expected by chance. Whatever generated these files was
//!   not this `teprob.f` stepped from those words, and no protocol choice
//!   recovers it: five significant figures of output and two decades of
//!   unrecorded `exp` and `pow` see to that on their own.
//! - **The operating point matches to a twentieth of a standard deviation.**
//!   `d00`'s 52 channels sit a median 0.037 published standard deviations
//!   apart, worst 0.16, with a median spread ratio of 1.08. The port is the
//!   same plant; it is not the same trajectory.
//! - **`temain_mod.f:367` was replaced, not kept.** The published files carry
//!   `IDV(12)` only in `d12` and `d12_te`. Reading the line the other way
//!   inflates the port's spread against `d00` by up to 9.9x and triples the
//!   median KS across every file in the table.
//! - **Four files shut down and the fifth nearly does.** `d06`, `d06_te`,
//!   `d18` and `d18_te` reach the 3000 kPa reactor-pressure limit and freeze;
//!   the port reaches it too, within ten minutes on three of the four. After a
//!   freeze the comparison is a point mass against a point mass and no
//!   distributional statistic can be small, which is why those are the worst
//!   four rows. `d13` is the mirror image: the port crosses where the
//!   published file did not.
//!
//! # Cost
//!
//! Forty-two runs at the recorded seed plus `TEP_TIER7_DECOYS` (default 4) per
//! file, half of them 25 simulated hours and half 48. About five minutes in
//! release, which is why this is not in `xtask ci --fast`.

#![cfg(feature = "oracle")]

use std::collections::BTreeMap;

use tepsim_oracle::tier5::{Run, VARIABLES, battery, published_seeds};
use tepsim_oracle::tier7::{
    self, D00_TRAINING_ROWS, DRIVER_SETTLING_STEPS, DRIVER_STEPS, PUBLISHED_COLUMNS, Protocol,
    Provenance, Published, SIGNIFICANT_DIGITS, Split, TESTING_ROWS, TRAINING_ROWS,
    UNRECORDED_AGITATOR, Unknown,
};
use tepsim_stats::Summary;

/// How many wrong-seed runs build the null. Four is enough to see a rank of
/// one out of five; more only sharpens a conclusion that is already obvious.
fn decoys() -> usize {
    std::env::var("TEP_TIER7_DECOYS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or_else(|| if full() { 4 } else { 1 })
}

/// Whether to run every published file, or a representative few.
///
/// The full sweep is 42 files by five seeds by 172,800 steps. That is 321
/// seconds in release and *thirty-five minutes* in the debug build
/// `cargo xtask ci` uses, which is not a per-commit cost. `TEP_TIER7=full`
/// asks for all of it; anything else, including the variable being absent,
/// runs [`SMOKE_FILES`].
///
/// The same idiom as `TEP_TIER1_SWEEP`, `TEP_TIER4_HOURS` and `TEP_TIER5`.
pub fn full() -> bool {
    std::env::var("TEP_TIER7").as_deref() == Ok("full")
}

/// The files the smoke sweep covers, chosen to span the outcomes rather than
/// to be the first few.
///
/// `d00` is the fault-free operating point and the only transposed file.
/// `d01` is an ordinary step fault. `d06` freezes on the reactor-pressure
/// shutdown, which is the case no distributional statistic can be small on.
/// `d12` is the one pair that genuinely carries `IDV(12)`, which is what the
/// forced-disturbance finding turns on.
pub const SMOKE_FILES: &[&str] = &["d00", "d01", "d06", "d12", "d12_te"];

/// Whether this file is in the current sweep.
pub fn included(name: impl AsRef<str>) -> bool {
    full() || SMOKE_FILES.contains(&name.as_ref())
}

/// Seeds no published file was made with.
///
/// Deliberately *not* other files' published seeds: those were used to
/// generate real datasets, and reusing one would make the null share a noise
/// realisation with some other published file. These are odd words inside
/// 2^32, which is what `TESUB7` wants, and nothing else.
fn decoy_seed(k: usize) -> f64 {
    let mut z = (k as u64 + 1).wrapping_mul(0xD1B5_4A32_D192_ED03);
    z = (z ^ (z >> 33)).wrapping_mul(0xFF51_AFD7_ED55_8CCD);
    z ^= z >> 29;
    ((z & 0xFFFF_FFFF) | 1) as f64
}

// ---------------------------------------------------------------------------
// The inventory: what is actually in reference/data.
// ---------------------------------------------------------------------------

#[test]
fn the_inventory_is_forty_four_files_and_only_d00_is_transposed() {
    let files = tier7::files();
    assert_eq!(files.len(), tier7::FILES);

    println!("\n=== Tier 7 file inventory ===");
    println!(
        "{:<9} {:>6} {:>5} {:>7}  shape as stored",
        "file", "rows", "cols", "raw"
    );
    for file in &files {
        let text = std::fs::read_to_string(file.path()).expect("a published file");
        let raw_rows = text.lines().filter(|l| !l.trim().is_empty()).count();
        let raw_cols = text
            .lines()
            .find(|l| !l.trim().is_empty())
            .map_or(0, |l| l.split_whitespace().count());
        let rows = file.load();

        assert!(
            rows.iter().all(|r| r.len() == PUBLISHED_COLUMNS),
            "{}: not 52 wide",
            file.name()
        );
        assert_eq!(
            rows.len(),
            file.protocol().rows,
            "{}: row count disagrees with the protocol",
            file.name()
        );

        let transposed = raw_cols != PUBLISHED_COLUMNS;
        assert_eq!(
            transposed,
            file.name() == "d00",
            "{}: only d00.dat is stored transposed",
            file.name()
        );
        println!(
            "{:<9} {:>6} {:>5} {:>3}x{:<3}  {}",
            file.name(),
            rows.len(),
            PUBLISHED_COLUMNS,
            raw_rows,
            raw_cols,
            if transposed {
                "TRANSPOSED (52 x 500)"
            } else {
                "samples x 52"
            }
        );
    }

    // The three row counts, and the one that upstream's own README gets wrong.
    let training: Vec<usize> = files
        .iter()
        .filter(|f| f.split == Split::Training)
        .map(|f| f.load().len())
        .collect();
    assert_eq!(training[0], D00_TRAINING_ROWS, "d00.dat holds 500, not 480");
    assert!(
        training[1..].iter().all(|n| *n == TRAINING_ROWS),
        "every other training file holds 480"
    );
    assert!(
        files
            .iter()
            .filter(|f| f.split == Split::Testing)
            .all(|f| f.load().len() == TESTING_ROWS)
    );
}

#[test]
fn the_published_values_carry_five_significant_digits() {
    // `temain_mod.f:1358` is `FORMAT(1X,E13.5)`. The claim is checkable from
    // the bytes: no published value should carry a sixth significant figure.
    let mut histogram = BTreeMap::new();
    let mut worst = (0usize, String::new());
    for file in tier7::files() {
        let text = std::fs::read_to_string(file.path()).expect("a published file");
        for token in text.split_whitespace() {
            let mantissa = token
                .split(['e', 'E'])
                .next()
                .expect("a mantissa")
                .replace(['-', '+', '.'], "");
            let digits = mantissa.trim_end_matches('0').trim_start_matches('0').len();
            *histogram.entry(digits).or_insert(0usize) += 1;
            if digits > worst.0 {
                worst = (digits, token.to_string());
            }
        }
    }
    println!("\n=== significant digits per published value ===");
    for (digits, count) in &histogram {
        println!("  {digits} digits: {count}");
    }
    println!("  widest value seen: {} ({} digits)", worst.1, worst.0);
    assert!(
        worst.0 <= SIGNIFICANT_DIGITS as usize,
        "a published value carries more than {SIGNIFICANT_DIGITS} figures: {}",
        worst.1
    );
}

// ---------------------------------------------------------------------------
// The protocol, as far as the source states it.
// ---------------------------------------------------------------------------

#[test]
fn the_testing_protocol_is_the_drivers_and_the_data_agrees() {
    // What the driver states.
    assert_eq!(DRIVER_STEPS, 172_800, "temain_mod.f:220, NPTS");
    assert_eq!(DRIVER_SETTLING_STEPS, 3600 * 8, "temain_mod.f:226, SSPTS");
    assert_eq!(DRIVER_STEPS / 180, TESTING_ROWS, "temain_mod.f:401");
    assert_eq!(
        Protocol::testing(true).source,
        Provenance::Stated,
        "every testing number is read off the driver"
    );

    // What the data says, independently. The row where each `_te` file
    // changes is found by a plain change-point search over the file's
    // standardised distance from the `d00` operating point: the split that
    // maximises the two-sample separation. There is no threshold in that, so
    // the detected row is not steered towards the answer the driver implies.
    println!("\n=== change point in the testing files, found without assuming one ===");
    println!(
        "{:<9} {:>11} {:>12} {:>8} {:>9}  note",
        "file", "change row", "change hours", "channel", "effect"
    );
    let mut detected = Vec::new();
    for file in tier7::files()
        .into_iter()
        .filter(|f| f.split == Split::Testing)
    {
        let rows = file.load();
        // The channel that changes most, rather than an average over 52: a
        // fault that moves two variables is invisible in a mean over all of
        // them, and which channel it is is itself worth printing.
        let mut best = (0usize, 0usize, f64::NEG_INFINITY, 0.0);
        for column in 0..PUBLISHED_COLUMNS {
            let series: Vec<f64> = rows.iter().map(|r| r[column]).collect();
            let (row, statistic, effect) = change_point(&series);
            if statistic > best.2 {
                best = (column, row, statistic, effect);
            }
        }
        let (column, row, _, effect) = best;
        println!(
            "{:<9} {:>11} {:>12.2} {:>8} {:>9.2}  {}",
            file.name(),
            row,
            row as f64 * 180.0 / 3600.0,
            channel_name(column),
            effect,
            if file.fault == 0 {
                "no fault requested"
            } else {
                ""
            }
        );
        if file.fault != 0 {
            detected.push((file.name(), row));
        }
    }

    // The driver switches the fault on at step 28,800, which is the boundary
    // between rows 159 and 160 zero-based. The claim is asserted on the modal
    // *hour* across the twenty-one fault files rather than on any single file:
    // a fault whose signature is weak has a change point set by noise, and
    // requiring every file to agree would be requiring the data to be
    // something it is not.
    let mut counts: BTreeMap<usize, usize> = BTreeMap::new();
    for (_, row) in &detected {
        *counts.entry(row / 20).or_default() += 1;
    }
    println!("\n  detected change point by hour: ");
    for (hour, count) in &counts {
        println!("    hour {hour:>2}: {}", "#".repeat(*count));
    }
    let (modal_hour, modal_count) = counts
        .iter()
        .max_by_key(|(hour, count)| (**count, std::cmp::Reverse(**hour)))
        .map(|(hour, count)| (*hour, *count))
        .expect("at least one fault file");
    println!(
        "  modal hour: {modal_hour}, in {modal_count} of {} fault files",
        detected.len()
    );
    assert_eq!(
        modal_hour, 8,
        "the testing files' modal change point is not the driver's eight hours"
    );
}

#[test]
fn the_training_protocol_is_not_documented_anywhere() {
    // The driver cannot produce a training file. It is 960 rows with the
    // fault at row 160; the training files are 480 rows with the fault at row
    // zero. Two constants were edited and the edit is not recorded.
    assert_ne!(DRIVER_STEPS / 180, TRAINING_ROWS);
    assert_eq!(
        Protocol::training(true).source,
        Provenance::Inferred,
        "the training protocol is a hypothesis, and the type says so"
    );

    // The evidence for the hypothesis, printed rather than asserted, because
    // it is evidence and not proof.
    let d00 = Published {
        fault: 0,
        split: Split::Training,
    };
    println!("\n=== the training protocol, reconstructed ===");
    println!("  d00.dat rows              : {}", d00.load().len());
    println!("  every other training file : {TRAINING_ROWS}");
    println!(
        "  difference                : {} rows = {} h",
        D00_TRAINING_ROWS - TRAINING_ROWS,
        (D00_TRAINING_ROWS - TRAINING_ROWS) as f64 * 180.0 / 3600.0
    );
    println!(
        "  500 rows at 180 s         : {} h",
        D00_TRAINING_ROWS as f64 * 180.0 / 3600.0
    );
    println!(
        "  hypothesis                : {} h run, fault at {} h, first {} rows dropped",
        Protocol::training(true).hours,
        Protocol::training(true).onset_hours.expect("a fault"),
        Protocol::training(true).discarded
    );
    println!("  competing hypothesis      : 24 h run, fault from step 1, nothing dropped");
    println!("  both put the first published row 3 minutes after onset, so the");
    println!("  fault response cannot separate them; only d00's 500 rows can, and");
    println!("  d00 has no fault to place.");

    // The faults are live in the first published training row, which is what
    // rules out a settling period inside the published section.
    let baseline = d00.load();
    let (centre, spread) = operating_point(&baseline);
    println!(
        "\n{:<9} {:>10} {:>10} {:>10}",
        "file", "row 0 z", "row 20 z", "row 479 z"
    );
    for file in tier7::files()
        .into_iter()
        .filter(|f| f.split == Split::Training && f.fault != 0)
    {
        let rows = file.load();
        println!(
            "{:<9} {:>10.2} {:>10.2} {:>10.2}",
            file.name(),
            distance(&rows[0], &centre, &spread),
            distance(&rows[20], &centre, &spread),
            distance(&rows[rows.len() - 1], &centre, &spread)
        );
    }
}

#[test]
// The seeds are integers held in `f64`, so exact equality is the only
// comparison that means anything about them and a margin would be nonsense.
#[allow(clippy::float_cmp)]
fn every_published_file_has_a_recorded_seed_and_two_have_more_than_one() {
    println!("\n=== seeds, teprob.f:1187-1256 ===");
    println!("{:<9} {:>14}  note", "file", "seed");
    for file in tier7::files() {
        let note = match (file.fault, file.split) {
            (0, Split::Training) => "also d00_tr_new = 5687912315",
            (18, Split::Training) => "also dd18_tr = 1234567890",
            (19, Split::Testing) => "source reads 9090909232.DO, letter O",
            (21, _) => "IDV(21) is not in this teprob.f",
            _ => "",
        };
        println!("{:<9} {:>14.0}  {note}", file.name(), file.seed());
        assert!(file.seed() > 0.0, "{}: no seed", file.name());
    }
    // The two alternatives are recorded and distinct from the primaries.
    assert_ne!(
        published_seeds::D00_TRAINING_NEW,
        published_seeds::TRAINING[0]
    );
    assert_ne!(
        published_seeds::D18_TRAINING_ALTERNATIVE,
        published_seeds::TRAINING[18]
    );
    assert_eq!(published_seeds::MALFORMED_TESTING_INDEX, 19);
    // d21 has a seed but no fault to apply it to.
    let d21 = Published {
        fault: 21,
        split: Split::Testing,
    };
    assert!(d21.seed() > 0.0);
    assert!(!d21.is_representable());
}

#[test]
fn the_agitator_is_constant_in_the_port() {
    // `UNRECORDED_AGITATOR` supplies the 53rd channel on the published side.
    // It is only legitimate because the port holds it there for the whole
    // run, which is checked here rather than assumed.
    let run = tier7::generate(
        &Published {
            fault: 0,
            split: Split::Training,
        },
        published_seeds::TRAINING[0],
        &Protocol {
            hours: 2.0,
            onset_hours: None,
            discarded: 0,
            rows: 40,
            source: Provenance::Stated,
        },
    );
    let agitator = run.series(PUBLISHED_COLUMNS);
    assert!(
        agitator.iter().all(|v| *v == UNRECORDED_AGITATOR),
        "XMV(12) moved, so supplying it as a constant would be fabrication"
    );
}

// ---------------------------------------------------------------------------
// The reproduction attempt.
// ---------------------------------------------------------------------------

/// The headline numbers for one published-versus-generated comparison.
#[derive(Clone, Copy, Debug)]
struct Agreement {
    ks_median: f64,
    ks_max: f64,
    energy_max: f64,
    acf_max: f64,
    spectrum_max: f64,
    frobenius: f64,
    z_mean_max: f64,
    z_mean_median: f64,
    logvar_max: f64,
    /// The channel behind `ks_max`, which is what says *where* a file fails.
    ks_worst: usize,
}

fn agreement(file: &Published, published: &Run, generated: &Run) -> Agreement {
    let report = battery::compare(
        if file.is_representable() {
            file.scenario()
        } else {
            tepsim_oracle::tier5::Scenario::NOMINAL
        },
        std::slice::from_ref(published),
        std::slice::from_ref(generated),
    );

    let mut ks = Vec::new();
    let mut energy = 0.0_f64;
    let mut acf = 0.0_f64;
    let mut spectrum = 0.0_f64;
    let mut z = Vec::new();
    let mut logvar = 0.0_f64;

    // Only the 52 channels a published file actually holds. The 53rd is the
    // agitator, which is a constant on both sides by construction.
    for variable in 0..PUBLISHED_COLUMNS {
        let entry = &report.variables[variable];
        ks.push(entry.ks.cross);
        energy = energy.max(entry.energy.cross);
        if entry.autocorrelation.cross.is_finite() {
            acf = acf.max(entry.autocorrelation.cross);
        }
        if entry.spectrum.cross.is_finite() {
            spectrum = spectrum.max(entry.spectrum.cross);
        }

        // TOST needs many runs per side and there is one published file, so
        // the moments are reported directly instead: the mean gap in units of
        // the published spread, and the log variance ratio.
        let a = Summary::of(&published.series(variable));
        let b = Summary::of(&generated.series(variable));
        if a.sd() > 0.0 {
            z.push((b.mean() - a.mean()).abs() / a.sd());
        }
        if a.variance() > 0.0 && b.variance() > 0.0 {
            logvar = logvar.max((b.variance() / a.variance()).ln().abs());
        }
    }

    let ks_worst = ks
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.partial_cmp(b.1).expect("no NaN"))
        .map_or(0, |(index, _)| index);

    Agreement {
        ks_worst,
        ks_median: median(&ks),
        ks_max: ks.iter().copied().fold(0.0_f64, f64::max),
        energy_max: energy,
        acf_max: acf,
        spectrum_max: spectrum,
        frobenius: report.structure.frobenius.cross,
        z_mean_max: z.iter().copied().fold(0.0_f64, f64::max),
        z_mean_median: median(&z),
        logvar_max: logvar,
    }
}

#[test]
fn tier7_per_file_agreement() {
    let decoys = decoys();
    let files = tier7::files();

    println!("\n=== Tier 7: published file against the port, per file ===");
    println!(
        "protocol: testing = 48 h / fault at 8 h (stated); training = 25 h / fault at 1 h / drop 20 (inferred)"
    );
    println!(
        "line 367 read as REPLACED, which is what the published bytes say; see\n\
         `the_published_files_were_not_generated_with_the_forced_idv12`. The `kept`\n\
         column is the same comparison with the driver's shipped IDV(12), for contrast."
    );
    println!("null: {decoys} decoy seeds per file, same protocol, wrong generator word\n");
    println!(
        "{:<9} {:>8} {:>8} {:>9} {:>8} {:>8} {:>8} {:>7} {:>6} {:>8} {:>9}",
        "file",
        "KS med",
        "KS max",
        "worst ch",
        "z mean",
        "logvar",
        "acf",
        "frob",
        "rank",
        "KS kept",
        "port trip"
    );

    let mut ranks = Vec::new();
    let mut rows = Vec::new();
    for file in &files {
        if !included(file.name()) {
            continue;
        }
        if !file.is_representable() {
            println!(
                "{:<9} {:>8} {:>8} {:>9} {:>8} {:>8} {:>8} {:>7} {:>6} {:>8} {:>9}   IDV(21) absent from this teprob.f",
                file.name(),
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-",
                "-"
            );
            continue;
        }
        let protocol = file.protocol();
        let published = file.run();
        let generated = tier7::generate_without_forced_idv12(file, file.seed(), &protocol);
        let trip = generated.tripped;
        let truth = agreement(file, &published, &generated);
        let kept = agreement(
            file,
            &published,
            &tier7::generate(file, file.seed(), &protocol),
        );
        let null: Vec<Agreement> = (0..decoys)
            .map(|k| {
                agreement(
                    file,
                    &published,
                    &tier7::generate_without_forced_idv12(file, decoy_seed(k), &protocol),
                )
            })
            .collect();

        // The rank of the recorded seed among itself and the decoys, by
        // median KS. One means the recorded seed is the best of the set,
        // which is what reproduction would look like.
        let rank = 1 + null
            .iter()
            .filter(|d| d.ks_median < truth.ks_median)
            .count();
        ranks.push(rank);

        for value in [
            truth.ks_median,
            truth.ks_max,
            truth.z_mean_max,
            truth.logvar_max,
            truth.frobenius,
        ] {
            assert!(
                value.is_finite(),
                "{}: a statistic came back non-finite, so the comparison did nothing",
                file.name()
            );
        }

        println!(
            "{:<9} {:>8.4} {:>8.4} {:>9} {:>8.3} {:>8.3} {:>8.4} {:>7.3} {:>4}/{} {:>8.4} {:>9}",
            file.name(),
            truth.ks_median,
            truth.ks_max,
            channel_name(truth.ks_worst),
            truth.z_mean_max,
            truth.logvar_max,
            truth.acf_max,
            truth.frobenius,
            rank,
            decoys + 1,
            kept.ks_median,
            trip.map_or_else(
                || "-".to_string(),
                |step| format!("{:.1} h", step as f64 / 3600.0)
            )
        );
        rows.push((file.name(), truth, null));
    }

    // What the decoys say the same statistics look like when the seed is
    // known to be wrong. If the recorded-seed column above sits inside this
    // band, the seed carries no information and reproduction did not happen.
    println!("\n=== the recorded seed against the {decoys}-seed null ===");
    println!(
        "{:<9} {:>9} {:>9} {:>9} | {:>9} {:>9} | {:>9} {:>9} {:>9}",
        "file", "KS seed", "KS lo", "KS hi", "z med", "z lo", "energy", "spectrum", "e null lo"
    );
    for (name, truth, null) in &rows {
        let lo = null
            .iter()
            .map(|d| d.ks_median)
            .fold(f64::INFINITY, f64::min);
        let hi = null
            .iter()
            .map(|d| d.ks_median)
            .fold(f64::NEG_INFINITY, f64::max);
        let zlo = null
            .iter()
            .map(|d| d.z_mean_median)
            .fold(f64::INFINITY, f64::min);
        let elo = null
            .iter()
            .map(|d| d.energy_max)
            .fold(f64::INFINITY, f64::min);
        println!(
            "{name:<9} {:>9.4} {lo:>9.4} {hi:>9.4} | {:>9.3} {zlo:>9.3} | {:>9.3} {:>9.3} {elo:>9.3}",
            truth.ks_median, truth.z_mean_median, truth.energy_max, truth.spectrum_max
        );
    }

    let firsts = ranks.iter().filter(|r| **r == 1).count();
    let expected = ranks.len() as f64 / (decoys + 1) as f64;
    println!(
        "\nrecorded seed ranked first on {firsts} of {} files; chance alone gives {expected:.1}",
        ranks.len()
    );
    println!(
        "A reproduction would put every file at rank 1 with the statistic near zero.\n\
         It does not, and the reasons are enumerated by `the_unknowns_are_enumerated`."
    );
}

#[test]
fn the_published_files_were_not_generated_with_the_forced_idv12() {
    // `temain_mod.f:367` is literally `IDV(12)=1`, inside the settling
    // conditional, and `temain_mod.f:101-102` tells the user to "go to line
    // 367" and type the disturbance wanted. Two readings:
    //
    //   replaced : the line is an example. Only d12 and d12_te carry IDV(12),
    //             and every other file is free of it.
    //   added to : the line stays. Every file carries IDV(12) from eight
    //             hours, on top of whatever was requested.
    //
    // Both are readings of the source, so neither can be picked by argument.
    // They make different predictions about the published bytes, and the
    // predictions are what decide it. Nothing here is fitted: the question is
    // binary and the source states both halves of it.
    let signature = [6usize, 12, 15, 19, 10, 21];
    let name = |c: usize| channel_name(c);

    // Prediction 1. If the line stayed, `d00_te` would carry IDV(12) from row
    // 160 and its spread would jump there the way `d12_te`'s does.
    let d00_te = Published {
        fault: 0,
        split: Split::Testing,
    }
    .load();
    let d12_te = Published {
        fault: 12,
        split: Split::Testing,
    }
    .load();
    println!("\n=== was IDV(12) left forced on? ===");
    println!("\nprediction 1: spread after row 160 over spread before it\n");
    println!("{:<10} {:>12} {:>12}", "channel", "d00_te", "d12_te");
    // Paired channel by channel, because the six differ by an order of
    // magnitude in how much IDV(12) moves them and a worst-against-best
    // comparison across channels would be comparing two different things.
    let mut weakest_pair = f64::INFINITY;
    for column in signature {
        let split = |rows: &[[f64; PUBLISHED_COLUMNS]]| {
            let pre: Vec<f64> = rows[..160].iter().map(|r| r[column]).collect();
            let post: Vec<f64> = rows[160..].iter().map(|r| r[column]).collect();
            Summary::of(&post).sd() / Summary::of(&pre).sd()
        };
        let (a, b) = (split(&d00_te), split(&d12_te));
        weakest_pair = weakest_pair.min(b / a);
        println!("{:<10} {a:>12.2} {b:>12.2}", name(column));
    }

    // Prediction 2. The port under each reading, against `d00`, which has no
    // fault of its own to confound it.
    let d00 = Published {
        fault: 0,
        split: Split::Training,
    };
    let published = d00.run();
    let protocol = d00.protocol();
    let forced = tier7::generate(&d00, d00.seed(), &protocol);
    let unforced = tier7::generate_without_forced_idv12(&d00, d00.seed(), &protocol);
    println!("\nprediction 2: port spread over published d00 spread\n");
    println!(
        "{:<10} {:>16} {:>16}",
        "channel", "line kept", "line replaced"
    );
    let mut worst_forced = 0.0_f64;
    let mut worst_unforced = 0.0_f64;
    for column in signature {
        let reference = Summary::of(&published.series(column)).sd();
        let a = Summary::of(&forced.series(column)).sd() / reference;
        let b = Summary::of(&unforced.series(column)).sd() / reference;
        worst_forced = worst_forced.max(a);
        worst_unforced = worst_unforced.max(b);
        println!("{:<10} {a:>16.2} {b:>16.2}", name(column));
    }

    println!(
        "\n  On every one of those channels d12_te's spread jumps at row 160 by at least\n  \
         {weakest_pair:.1}x what d00_te's does. Keeping line 367 inflates the port's spread\n  \
         against d00 by up to {worst_forced:.1}x; replacing it, by {worst_unforced:.2}x. Both\n  \
         predictions point the same way: the published files were generated with\n  \
         line 367 REPLACED, so only d12 and d12_te carry IDV(12).\n"
    );
    println!(
        "  This does not change the port's default, which reproduces the driver as\n  \
         shipped (delta D-011, a Class C quirk needing sign-off). It says the\n  \
         published datasets were not made with the driver as shipped."
    );

    assert!(
        weakest_pair > 4.0,
        "d12_te's spread does not jump distinguishably more than d00_te's at row 160, \
         so this diagnostic has no discriminating power and its conclusion is unsupported"
    );
    assert!(
        worst_unforced < worst_forced,
        "replacing line 367 does not bring the port closer to d00's spread, which \
         contradicts the conclusion printed above"
    );
}

#[test]
fn the_nominal_operating_point_side_by_side() {
    // Due diligence before reporting a non-reproduction. `d00` is the one file
    // with no fault to place and no onset to guess, so if the port and the
    // published data disagree *here* the disagreement is either the operating
    // point itself or something wrong with this harness, and either way it is
    // not a fault-protocol question. Printed per channel so a reader can see
    // which it is rather than take a summary statistic on trust.
    let d00 = Published {
        fault: 0,
        split: Split::Training,
    };
    let published = d00.run();
    let protocol = d00.protocol();

    // The two readings of `temain_mod.f:367`, side by side. Neither is tuned:
    // one is the line as shipped, the other is the line replaced as the
    // instructions at `temain_mod.f:101-102` describe.
    let forced = tier7::generate(&d00, d00.seed(), &protocol);
    let unforced = tier7::generate_without_forced_idv12(&d00, d00.seed(), &protocol);

    println!("\n=== d00: the nominal operating point, published against the port ===");
    println!(
        "25 h, no fault, line 367 replaced. The `kept` column is the same run with\n\
         the driver's shipped IDV(12), which is the size of that one decision.\n"
    );
    println!(
        "{:<9} {:>12} {:>12} {:>8} | {:>10} {:>10} {:>9} | {:>9}",
        "channel", "pub mean", "port mean", "z", "pub sd", "port sd", "sd ratio", "kept"
    );
    let mut worst = (0usize, 0.0_f64);
    let mut worst_ratio = (0usize, 1.0_f64);
    for column in 0..PUBLISHED_COLUMNS {
        let a = Summary::of(&published.series(column));
        let b = Summary::of(&unforced.series(column));
        let c = Summary::of(&forced.series(column));
        let z = if a.sd() > 0.0 {
            (b.mean() - a.mean()).abs() / a.sd()
        } else {
            0.0
        };
        let ratio = if a.sd() > 0.0 { b.sd() / a.sd() } else { 1.0 };
        let kept_ratio = if a.sd() > 0.0 { c.sd() / a.sd() } else { 1.0 };
        if z > worst.1 {
            worst = (column, z);
        }
        if (ratio.ln()).abs() > worst_ratio.1.ln().abs() {
            worst_ratio = (column, ratio);
        }
        println!(
            "{:<9} {:>12.4} {:>12.4} {:>8.2} | {:>10.4} {:>10.4} {:>9.3} | {:>9.3}",
            channel_name(column),
            a.mean(),
            b.mean(),
            z,
            a.sd(),
            b.sd(),
            ratio,
            kept_ratio
        );
    }
    println!(
        "\n  worst mean gap : {} at {:.2} published standard deviations",
        channel_name(worst.0),
        worst.1
    );
    println!(
        "  worst sd ratio : {} at {:.3}",
        channel_name(worst_ratio.0),
        worst_ratio.1
    );

    // The operating point itself has to agree, or nothing downstream means
    // anything. A mean sitting one published standard deviation away is a
    // different plant, not a different noise realisation, and would say this
    // harness is wrong rather than that the protocol is unknown.
    let gaps: Vec<f64> = (0..PUBLISHED_COLUMNS)
        .filter_map(|c| {
            let a = Summary::of(&published.series(c));
            let b = Summary::of(&unforced.series(c));
            (a.sd() > 0.0).then(|| (b.mean() - a.mean()).abs() / a.sd())
        })
        .collect();
    let ratios: Vec<f64> = (0..PUBLISHED_COLUMNS)
        .filter_map(|c| {
            let a = Summary::of(&published.series(c));
            let b = Summary::of(&unforced.series(c));
            (a.sd() > 0.0).then(|| b.sd() / a.sd())
        })
        .collect();
    println!("  median mean gap: {:.3} published sd", median(&gaps));
    println!("  median sd ratio: {:.3}", median(&ratios));

    // Is the residual spread difference bigger than the port's own
    // seed-to-seed wander? Four more runs at decoy seeds say what that wander
    // is, so `1.55` can be read against something instead of eyeballed.
    let others: Vec<Run> = (0..4)
        .map(|k| tier7::generate_without_forced_idv12(&d00, decoy_seed(k), &protocol))
        .collect();
    println!("\n  the same ratio at four other seeds, on the channels that moved most:");
    println!(
        "{:<10} {:>10} {:>26}",
        "channel", "at d00_tr", "at four other seeds"
    );
    for column in [6usize, 12, 15, 10, 19] {
        let reference = Summary::of(&published.series(column)).sd();
        let here = Summary::of(&unforced.series(column)).sd() / reference;
        let elsewhere: Vec<String> = others
            .iter()
            .map(|r| format!("{:.2}", Summary::of(&r.series(column)).sd() / reference))
            .collect();
        println!(
            "{:<10} {here:>10.2} {:>26}",
            channel_name(column),
            elsewhere.join(" ")
        );
    }
    assert!(
        median(&gaps) < 1.0,
        "the port's nominal operating point is more than a published standard \
         deviation from d00's on half its channels, which is a harness fault \
         rather than a protocol unknown"
    );
}

#[test]
// The freeze test compares consecutive samples for *bit* equality on purpose:
// `teprob.f:807-811` zeroes every derivative, so a frozen plant repeats its
// state exactly, and a running one never does. A margin would blur the very
// distinction being detected.
#[allow(clippy::float_cmp)]
fn the_files_that_shut_down_disagree_on_when() {
    // Four published files reach the reactor-pressure shutdown at
    // `teprob.f:703`, `IF(XMEAS(7).GT.3000.0) ISD=1`. On that condition the
    // Fortran zeroes all fifty derivatives and stops adding measurement noise
    // (`teprob.f:711`), so the file's tail is frozen and the freeze row is
    // readable off the bytes to the sample.
    //
    // That makes the trip a *threshold crossing on a trajectory*, and Tier 4
    // already measured trajectories separating within hours on `exp` and `pow`
    // rounding alone. Two simulators that agree everywhere except in the last
    // bits will cross 3000 kPa at different times, and after that they are
    // comparing a frozen plant against a running one. This is the specific
    // reason d06 and d18 are the worst rows in the agreement table, and it is
    // not a protocol unknown: it is arithmetic.
    println!("\n=== the four files that hit the 3000 kPa shutdown ===");
    println!(
        "{:<9} {:>14} {:>14} {:>16} {:>14}",
        "file", "published row", "published h", "port trip h", "port peak kPa"
    );
    let mut compared = 0;
    for file in tier7::files() {
        if !included(file.name()) || !file.is_representable() {
            continue;
        }
        let rows = file.load();
        let pressure: Vec<f64> = rows.iter().map(|r| r[6]).collect();
        // The freeze row: the first sample from which reactor pressure never
        // changes again. A running plant's pressure is never bit-identical
        // twice in a row, so this is unambiguous.
        let Some(freeze) =
            (0..pressure.len() - 1).find(|r| pressure[*r..].windows(2).all(|w| w[0] == w[1]))
        else {
            continue;
        };
        compared += 1;
        let protocol = file.protocol();
        let generated = tier7::generate_without_forced_idv12(&file, file.seed(), &protocol);
        let port_peak = generated
            .series(6)
            .into_iter()
            .fold(f64::NEG_INFINITY, f64::max);
        println!(
            "{:<9} {freeze:>14} {:>14.2} {:>16} {port_peak:>14.1}",
            file.name(),
            (freeze + 1 + protocol.discarded) as f64 * 180.0 / 3600.0,
            generated.tripped.map_or_else(
                || "did not trip".to_string(),
                |step| format!("{:.2}", step as f64 / 3600.0)
            )
        );
    }
    println!(
        "\n  A trip is a threshold crossing, so the two simulators need only differ in\n  \
         the last bits for it to land in a different hour. After it, one plant is\n  \
         frozen and the other is not, and no distributional statistic can be small."
    );
    assert!(
        compared > 0,
        "no published file freezes, so this diagnostic found nothing to compare and \
         the explanation it prints is unsupported"
    );

    // The other half of the same story, and the one that is not explained by a
    // threshold crossing at all: a file the port shuts down on and the
    // published data does not. `d13` is the case, and it is the reason that
    // row of the agreement table has a mean 10 published standard deviations
    // out while `d13_te` agrees as well as any file does.
    //
    // Whether the trip is luck or the model's normal behaviour under a 25-hour
    // `IDV(13)` is answerable by running other seeds, which is what this does.
    println!("\n=== where the port shuts down and the published file does not ===");
    println!(
        "{:<9} {:>12} {:>34}",
        "file", "published", "port trip at 5 seeds"
    );
    for file in tier7::files().into_iter().filter(|f| f.fault == 13) {
        let rows = file.load();
        let pressure: Vec<f64> = rows.iter().map(|r| r[6]).collect();
        let freezes =
            (0..pressure.len() - 1).any(|r| pressure[r..].windows(2).all(|w| w[0] == w[1]));
        let protocol = file.protocol();
        let trips: Vec<String> = std::iter::once(file.seed())
            .chain((0..4).map(decoy_seed))
            .map(|seed| {
                tier7::generate_without_forced_idv12(&file, seed, &protocol)
                    .tripped
                    .map_or_else(
                        || "-".to_string(),
                        |step| format!("{:.1}h", step as f64 / 3600.0),
                    )
            })
            .collect();
        println!(
            "{:<9} {:>12} {:>34}",
            file.name(),
            if freezes { "freezes" } else { "runs on" },
            trips.join(" ")
        );
    }
    println!(
        "\n  The port's 25-hour IDV(13) run crosses the reactor-pressure limit at the\n  \
         seed recorded for d13 and at none of four other seeds, and its 48-hour run\n  \
         crosses at none of the five. So the crossing is rare rather than routine,\n  \
         and it is what puts d13's mean 10 published standard deviations out while\n  \
         d13_te agrees as well as any file does. Same mechanism as d06 and d18: a\n  \
         threshold crossing on a diverging trajectory, not a separate unknown.\n  \
         What it is NOT is evidence that the port trips where the Fortran does not.\n  \
         That comparison needs the oracle driven under this protocol, which Tier 7\n  \
         does not do and Tier 5 does not cover, because Tier 5 has no delayed onset."
    );
}

#[test]
fn neither_training_duration_hypothesis_beats_the_other() {

    // The conclusion is a count across the training files: neither hypothesis
    // wins on enough of them to separate. The smoke sweep has four files, and
    // "9 against 11" is not a statement two files can make.
    if !full() {
        println!(
            "skipped: comparing two training-duration hypotheses needs the \
             whole set of training files. Run with TEP_TIER7=full."
        );
        return;
    }
    // `Unknown::TrainingProtocolIsUndocumented`, measured rather than argued.
    // Two readings produce 480 rows with the fault three minutes before the
    // first of them:
    //
    //   25 h : run 25 hours, fault at 1 h, drop the first 20 rows. Explains
    //          d00.dat's 500 rows as the same run with nothing dropped.
    //   24 h : run 24 hours with the fault live from the first step. Explains
    //          nothing about d00.dat.
    //
    // Both are printed. Neither is adopted: they are two candidate protocols,
    // and picking one because it scored better on the data would be exactly
    // the fitting `PLAN.org` forbids. What this test can do is say whether
    // the choice matters at all.
    let twenty_four = Protocol {
        hours: 24.0,
        onset_hours: Some(0.0),
        discarded: 0,
        rows: 480,
        source: Provenance::Inferred,
    };

    println!("\n=== the two training-duration hypotheses ===");
    println!(
        "{:<9} {:>14} {:>14} {:>14} {:>14}",
        "file", "KS med 25 h", "KS med 24 h", "frob 25 h", "frob 24 h"
    );
    let mut twenty_five_wins = 0;
    let mut twenty_four_wins = 0;
    let mut worst_gap = 0.0_f64;
    for file in tier7::files().into_iter().filter(|f| {
        f.split == Split::Training
            && f.fault != 0
            && f.is_representable()
            && included(f.name())
    }) {
        let published = file.run();
        let a = agreement(
            &file,
            &published,
            &tier7::generate_without_forced_idv12(&file, file.seed(), &file.protocol()),
        );
        let b = agreement(
            &file,
            &published,
            &tier7::generate_without_forced_idv12(&file, file.seed(), &twenty_four),
        );
        if a.ks_median < b.ks_median {
            twenty_five_wins += 1;
        } else {
            twenty_four_wins += 1;
        }
        worst_gap = worst_gap.max((a.ks_median - b.ks_median).abs());
        println!(
            "{:<9} {:>14.4} {:>14.4} {:>14.3} {:>14.3}",
            file.name(),
            a.ks_median,
            b.ks_median,
            a.frobenius,
            b.frobenius
        );
    }
    println!(
        "\n  25 h closer on {twenty_five_wins} files, 24 h closer on {twenty_four_wins}; \
         largest KS gap {worst_gap:.4}."
    );
    println!(
        "  The two hypotheses are not separated by this data. That is the finding:\n  \
         the training protocol stays an unknown, and d00.dat's 500 rows remain the\n  \
         only argument for the 25-hour reading."
    );

    // The split has to be close to even, or one hypothesis really is better
    // and the unknown would be overstated. Stated as a two-sided check on the
    // count, with no tuning: 21 files, so a clean win would be near 21-0.
    assert!(
        twenty_five_wins > 0 && twenty_four_wins > 0,
        "one training-duration hypothesis fits every file and the other fits none, \
         which would settle Unknown::TrainingProtocolIsUndocumented and this test \
         claims it does not"
    );
}

#[test]
fn the_unknowns_are_enumerated() {
    println!("\n=== Tier 7 unknowns: what stops reproduction ===");
    let mut open = 0;
    for (n, unknown) in Unknown::all().iter().enumerate() {
        println!("{:>2}. {unknown:?}\n    scope: {}", n + 1, unknown.scope());
        match unknown.settled_by_evidence() {
            Some(answer) => println!("    SETTLED by the published data: {answer}"),
            None => {
                open += 1;
                println!("    OPEN: neither the source nor the data decides it");
            }
        }
    }
    println!("\n  {open} of {} remain open.", Unknown::all().len());
    // The enumeration is the deliverable of B-0051, so its size is asserted:
    // a variant removed without a decision would otherwise vanish silently.
    assert_eq!(Unknown::all().len(), 8);
}

// ---------------------------------------------------------------------------
// Helpers.
// ---------------------------------------------------------------------------

/// Column means and standard deviations of a published file.
fn operating_point(
    rows: &[[f64; PUBLISHED_COLUMNS]],
) -> ([f64; PUBLISHED_COLUMNS], [f64; PUBLISHED_COLUMNS]) {
    let mut centre = [0.0; PUBLISHED_COLUMNS];
    let mut spread = [0.0; PUBLISHED_COLUMNS];
    for column in 0..PUBLISHED_COLUMNS {
        let series: Vec<f64> = rows.iter().map(|r| r[column]).collect();
        let summary = Summary::of(&series);
        centre[column] = summary.mean();
        spread[column] = summary.sd();
    }
    (centre, spread)
}

/// Mean standardised distance of one row from an operating point.
fn distance(
    row: &[f64; PUBLISHED_COLUMNS],
    centre: &[f64; PUBLISHED_COLUMNS],
    spread: &[f64; PUBLISHED_COLUMNS],
) -> f64 {
    let mut total = 0.0;
    let mut counted = 0;
    for column in 0..PUBLISHED_COLUMNS {
        if spread[column] > 0.0 {
            total += ((row[column] - centre[column]) / spread[column]).abs();
            counted += 1;
        }
    }
    if counted == 0 {
        0.0
    } else {
        total / f64::from(counted)
    }
}

/// The split of a series that best separates it into two levels.
///
/// Maximises `r (n - r) / n * (mean_before - mean_after)^2`, the usual
/// single-change-point statistic, standardised by the series' variance so that
/// the statistic is comparable across channels of wildly different scale.
/// Returns the row the change begins at (zero-based), that statistic, and the
/// level difference in units of the series' own standard deviation.
///
/// No threshold: a series with no change still returns its best split, and the
/// effect size is what says whether to believe it.
fn change_point(series: &[f64]) -> (usize, f64, f64) {
    let n = series.len();
    let total: f64 = series.iter().sum();
    let variance = Summary::of(series).variance();
    if variance.is_nan() || variance <= 0.0 {
        return (0, f64::NEG_INFINITY, 0.0);
    }
    let mut prefix = 0.0;
    let mut best = (0usize, f64::NEG_INFINITY, 0.0);
    // Both halves need enough rows for a mean to mean anything.
    for r in 1..n {
        prefix += series[r - 1];
        if r < 20 || n - r < 20 {
            continue;
        }
        let before = prefix / r as f64;
        let after = (total - prefix) / (n - r) as f64;
        let gap = after - before;
        let statistic = (r * (n - r)) as f64 / n as f64 * gap * gap / variance;
        if statistic > best.1 {
            best = (r, statistic, gap / variance.sqrt());
        }
    }
    best
}

/// `XMEAS(n)` or `XMV(n)`, for a published column index.
fn channel_name(column: usize) -> String {
    if column < 41 {
        format!("XMEAS{}", column + 1)
    } else {
        format!("XMV{}", column - 40)
    }
}

fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    sorted[sorted.len() / 2]
}

/// Keeps `VARIABLES` referenced: the published files are 52 wide and a `Run`
/// is 53, and the difference is the whole reason `UNRECORDED_AGITATOR` exists.
#[test]
fn a_run_is_one_channel_wider_than_a_published_file() {
    assert_eq!(VARIABLES, PUBLISHED_COLUMNS + 1);
}
