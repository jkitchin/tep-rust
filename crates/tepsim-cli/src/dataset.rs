//! `tep dataset`: generate a `d00`-`d21` shaped dataset.
//!
//! The geometry, the seeds and the column order all come from
//! [`tepsim::published`], which asserts them against the vendored Fortran. This
//! file is the part that writes bytes.
//!
//! # What this produces is not the published data
//!
//! It has the same shape and it is not the same numbers, and the subcommand
//! says so on stderr every time it runs rather than leaving that in a manual.
//! `crates/tepsim-oracle/src/tier7.rs` measures the gap. The short version is
//! that the toolchain that made the shipped files is unrecorded, its `exp` was
//! not this one's, the training protocol was never written down, and `d21`
//! names a disturbance this revision of the model does not contain.

use std::fs;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use tepsim::published::{self, COLUMNS, File, Split, Unavailable, rieth};
use tepsim::{Outcome, Run, Simulation};

/// How the rows are written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Format {
    /// The published layout: 52 space-separated fields at five significant
    /// digits, no header.
    Dat,
    /// Comma-separated with a header row, at full precision.
    Csv,
}

impl Format {
    fn parse(raw: &str) -> Option<Self> {
        match raw {
            "dat" | "published" => Some(Self::Dat),
            "csv" => Some(Self::Csv),
            _ => None,
        }
    }

    const fn extension(self) -> &'static str {
        match self {
            Self::Dat => "dat",
            Self::Csv => "csv",
        }
    }
}

/// Which distribution's shape to produce.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Layout {
    /// The original `d00`-`d21` files: one run each, at the seeds
    /// `teprob.f:1187-1256` records.
    Published,
    /// `rieth-2017-addit`: many runs per fault, in one long-format file per
    /// split, with `faultNumber`, `simulationRun` and `sample` columns.
    ///
    /// This is the shape most current machine-learning work on the process is
    /// trained against, and it is *not* the original distribution. Its testing
    /// protocol matches; its training protocol does not. See
    /// `tepsim::published::rieth`.
    Rieth,
}

impl Layout {
    fn parse(raw: &str) -> Option<Self> {
        match raw {
            "published" | "original" => Some(Self::Published),
            "rieth" => Some(Self::Rieth),
            _ => None,
        }
    }
}

/// Which half of the dataset to generate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Which {
    /// `dNN_te.dat` only. The default, because it is the half the driver
    /// documents.
    Testing,
    /// `dNN.dat` only.
    Training,
    /// Both.
    All,
}

impl Which {
    fn parse(raw: &str) -> Option<Self> {
        match raw {
            "te" | "test" | "testing" => Some(Self::Testing),
            "tr" | "train" | "training" => Some(Self::Training),
            "all" | "both" => Some(Self::All),
            _ => None,
        }
    }

    const fn wants(self, split: Split) -> bool {
        match self {
            Self::All => true,
            Self::Testing => matches!(split, Split::Testing),
            Self::Training => matches!(split, Split::Training),
        }
    }
}

/// What `tep dataset` was asked to do.
#[derive(Debug)]
pub(crate) struct Options {
    out: PathBuf,
    which: Which,
    faults: Option<Vec<usize>>,
    format: Format,
    seed: Option<f64>,
    force_idv12: bool,
    list_only: bool,
    layout: Layout,
    runs: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            out: PathBuf::from("tep-data"),
            which: Which::Testing,
            faults: None,
            format: Format::Dat,
            seed: None,
            force_idv12: false,
            list_only: false,
            layout: Layout::Published,
            runs: 1,
        }
    }
}

/// Parse the flags after `tep dataset`.
///
/// # Errors
///
/// A message for the user, on an unknown flag or an unusable value.
pub(crate) fn parse(args: &[String]) -> Result<Options, String> {
    let mut options = Options::default();
    let mut rest = args.iter();

    while let Some(flag) = rest.next() {
        let mut value = || {
            rest.next()
                .ok_or_else(|| format!("`{flag}` needs a value"))
                .cloned()
        };
        match flag.as_str() {
            "--out" => options.out = PathBuf::from(value()?),
            "--set" => {
                let raw = value()?;
                options.which = Which::parse(&raw)
                    .ok_or_else(|| format!("`--set {raw}` is not one of te, training, all"))?;
            }
            "--format" => {
                let raw = value()?;
                options.format = Format::parse(&raw)
                    .ok_or_else(|| format!("`--format {raw}` is not one of dat, csv"))?;
            }
            "--faults" => {
                let raw = value()?;
                if raw == "all" {
                    options.faults = None;
                } else {
                    let mut wanted = Vec::new();
                    for part in raw.split(',') {
                        let n: usize = part
                            .trim()
                            .parse()
                            .map_err(|_| format!("`{part}` in `--faults` is not a number"))?;
                        if n > 21 {
                            return Err(format!(
                                "`--faults {n}`: the published dataset stops at 21"
                            ));
                        }
                        wanted.push(n);
                    }
                    if wanted.is_empty() {
                        return Err("`--faults` with nothing in it".to_string());
                    }
                    options.faults = Some(wanted);
                }
            }
            "--seed" => {
                let raw = value()?;
                let seed: f64 = raw
                    .parse()
                    .map_err(|_| format!("`--seed {raw}` is not a number"))?;
                if seed.is_nan() || seed <= 0.0 {
                    return Err("`--seed` must be positive".to_string());
                }
                options.seed = Some(seed);
            }
            "--layout" => {
                let raw = value()?;
                options.layout = Layout::parse(&raw)
                    .ok_or_else(|| format!("`--layout {raw}` is not one of published, rieth"))?;
            }
            "--runs" => {
                let raw = value()?;
                let runs: usize = raw
                    .parse()
                    .map_err(|_| format!("`--runs {raw}` is not a number"))?;
                if runs == 0 {
                    return Err("`--runs 0` would generate nothing".to_string());
                }
                options.runs = runs;
            }
            "--force-idv12" => options.force_idv12 = true,
            "--list" => options.list_only = true,
            other => return Err(format!("unknown option `{other}`")),
        }
    }
    Ok(options)
}

/// Whether a file is one the caller asked for.
fn wanted(options: &Options, file: &File) -> bool {
    options.which.wants(file.split)
        && options
            .faults
            .as_ref()
            .is_none_or(|list| list.contains(&file.fault))
}

/// Width of one published field. 52 of these is [`PUBLISHED_LINE_BYTES`].
const FIELD_WIDTH: usize = 16;

/// Bytes in one row of a published file, newline excluded.
///
/// Every line of every shipped `.dat` file is exactly this long, which is how
/// the layout was determined: a value is right-justified in sixteen columns
/// rather than separated by three spaces. The two look identical until a
/// negative appears, and negatives do appear, 25 of them in `d00_te.dat`
/// alone. `the_generated_layout_is_the_published_layout` pins it.
const PUBLISHED_LINE_BYTES: usize = COLUMNS * FIELD_WIDTH;

/// One published field, in the layout the shipped files use.
///
/// The shipped files are not raw `E13.5`: they were reformatted at some point
/// into `%16.7e`, which is why a value reads `2.5025000e-01` with three
/// trailing zeros that carry no information. The five significant digits are
/// real and the rest is padding, so the value is rounded first and then
/// printed wide. Rust writes a one-character exponent where C writes two, so
/// the exponent is rebuilt here.
fn field(value: f64) -> String {
    let raw = format!("{:.7e}", published::round_as_published(value));
    let (mantissa, exponent) = raw
        .split_once('e')
        .expect("`{:e}` always writes an exponent");
    let exponent: i32 = exponent.parse().expect("`{:e}` writes an integer exponent");
    let sign = if exponent < 0 { '-' } else { '+' };
    format!("{mantissa}e{sign}{:02}", exponent.abs())
}

/// One row, right-justified into fixed-width columns.
fn published_line(row: &[f64; COLUMNS]) -> String {
    let mut line = String::with_capacity(PUBLISHED_LINE_BYTES);
    for value in row {
        let text = field(*value);
        for _ in text.len()..FIELD_WIDTH {
            line.push(' ');
        }
        line.push_str(&text);
    }
    line
}

/// Write one run out, already trimmed to the rows the file keeps.
fn write_file(path: &PathBuf, run: &Run, skip: usize, format: Format) -> io::Result<usize> {
    let file = fs::File::create(path)?;
    let mut out = BufWriter::new(file);
    let mut rows = 0;

    if format == Format::Csv {
        writeln!(out, "{}", published::column_names().join(","))?;
    }
    for sample in run.samples.iter().skip(skip) {
        let row = published::row(sample);
        match format {
            Format::Dat => writeln!(out, "{}", published_line(&row))?,
            Format::Csv => {
                let cells: Vec<String> = row.iter().map(|v| format!("{v:.17e}")).collect();
                writeln!(out, "{}", cells.join(","))?;
            }
        }
        rows += 1;
    }
    out.flush()?;
    Ok(rows)
}

/// Write the Rieth-shaped ensemble: one long-format file per split.
///
/// # Long format, not one file per run
///
/// Their distribution ships four files, each stacking every run of every fault
/// with `faultNumber`, `simulationRun` and `sample` identifying the row. That
/// is what a dataframe wants and it is what the papers built on this expect, so
/// it is what this writes, rather than twenty-one thousand small files.
fn write_rieth(options: &Options) -> Result<(), String> {
    let faults: Vec<usize> = options.faults.clone().unwrap_or_else(|| (0..=20).collect());
    let splits: Vec<Split> = match options.which {
        Which::Testing => vec![Split::Testing],
        Which::Training => vec![Split::Training],
        Which::All => vec![Split::Training, Split::Testing],
    };

    // The cost, before starting rather than after. A 48 hour run is about
    // three quarters of a second in release, and this is the one subcommand
    // that can be asked for tens of thousands of them.
    // Seconds per run, measured in release on the development machine: a
    // 25 hour training run and a 48 hour testing one. Only ever used to print
    // an estimate, so being wrong by a factor of two costs a sentence.
    const TRAINING_SECONDS: f64 = 0.42;
    const TESTING_SECONDS: f64 = 0.80;
    let runs = splits.len() * faults.len() * options.runs;
    let estimate: f64 = splits
        .iter()
        .map(|split| {
            let per = match split {
                Split::Training => TRAINING_SECONDS,
                Split::Testing => TESTING_SECONDS,
            };
            per * (faults.len() * options.runs) as f64
        })
        .sum();
    eprintln!(
        "tep: {runs} run(s): {} fault(s) x {} run(s) x {} split(s), about {}",
        faults.len(),
        options.runs,
        splits.len(),
        if estimate < 90.0 {
            format!("{estimate:.0} s")
        } else {
            format!("{:.0} min", estimate / 60.0)
        }
    );

    for split in splits {
        let name = match split {
            Split::Training => "training",
            Split::Testing => "testing",
        };
        let path = options.out.join(format!("{name}.csv"));
        let file =
            fs::File::create(&path).map_err(|e| format!("creating {}: {e}", path.display()))?;
        let mut out = BufWriter::new(file);

        // The identifying columns their files carry, then the 52 variables.
        write!(out, "faultNumber,simulationRun,sample").map_err(io_err(&path))?;
        for name in published::column_names() {
            write!(out, ",{name}").map_err(io_err(&path))?;
        }
        writeln!(out).map_err(io_err(&path))?;

        let mut written = 0_usize;
        let mut skipped = Vec::new();
        for &fault in &faults {
            for run in 0..options.runs {
                let scenario = match rieth::scenario(fault, split, run) {
                    Ok(scenario) => scenario,
                    Err(Unavailable::FaultNotInThisRevision { fault }) => {
                        if run == 0 {
                            skipped.push(format!("IDV({fault}) is not in this revision"));
                        }
                        continue;
                    }
                };
                let finished = Simulation::new(scenario).run();
                for (index, sample) in finished.samples.iter().enumerate() {
                    // One-based, as their `sample` column is.
                    write!(out, "{fault},{},{}", run + 1, index + 1).map_err(io_err(&path))?;
                    for value in published::row(sample) {
                        write!(out, ",{value:.17e}").map_err(io_err(&path))?;
                    }
                    writeln!(out).map_err(io_err(&path))?;
                }
                written += 1;
            }
            eprintln!("tep: {name}: IDV({fault}) done, {written} run(s) so far");
        }
        out.flush().map_err(io_err(&path))?;
        eprintln!(
            "tep: {} {written} run(s) x {} rows",
            path.display(),
            rieth::rows(split)
        );
        for reason in &skipped {
            eprintln!("tep: skipped {reason}");
        }
    }
    Ok(())
}

/// A closure that turns an IO error into a message naming the file.
fn io_err(path: &std::path::Path) -> impl Fn(io::Error) -> String + '_ {
    move |e| format!("writing {}: {e}", path.display())
}

/// Run the subcommand.
///
/// # Errors
///
/// A message for the user, on any filesystem failure.
pub(crate) fn run(options: &Options) -> Result<(), String> {
    let files: Vec<&File> = published::FILES
        .iter()
        .filter(|f| wanted(options, f))
        .collect();
    if files.is_empty() {
        return Err("no files match that selection".to_string());
    }

    if options.list_only {
        println!("{:<10} {:>7} {:>7} {:>14}", "file", "rows", "cols", "seed");
        for file in &files {
            println!(
                "{:<10} {:>7} {:>7} {:>14.0}",
                format!("{}.{}", file.stem(), options.format.extension()),
                file.rows,
                COLUMNS,
                options.seed.unwrap_or(file.seed)
            );
        }
        return Ok(());
    }

    fs::create_dir_all(&options.out)
        .map_err(|e| format!("creating {}: {e}", options.out.display()))?;

    eprintln!(
        "tep: this is published-SHAPED data, not the published data. The \
         toolchain that made the shipped files is unrecorded and its `exp` was \
         not this one's; see `tepsim::published` and Tier 7."
    );

    if options.layout == Layout::Rieth {
        eprintln!(
            "tep: rieth-2017 layout. Their runs use seeds that are not recorded \
             anywhere this repository can read, so these are distinct and \
             reproducible seeds of our own, not theirs."
        );
        eprintln!(
            "tep: note their training protocol is NOT the original's: 25 h with \
             the fault at hour one and all 500 rows kept, where d01-d21 discard \
             the first hour. Both are 25 h and both split 20 against 480."
        );
        return write_rieth(options);
    }
    if options.which != Which::Testing {
        eprintln!(
            "tep: the training protocol is a hypothesis, not a record: 25 h with \
             the fault from step zero, first hour discarded except for d00. It is \
             the only hypothesis that explains d00.dat's 500 rows and the other \
             files' 480 at the same time, and it is not written down anywhere."
        );
    }
    if options.seed.is_some() {
        eprintln!("tep: --seed overrides the per-file seeds teprob.f records");
    }
    if options.force_idv12 {
        eprintln!(
            "tep: --force-idv12 adds IDV(12) at hour eight to every run, which \
             the published bytes say was NOT done"
        );
    }

    let mut written = 0_usize;
    let mut skipped = Vec::new();
    for file in &files {
        let stem = file.stem();
        let mut scenario = match file.scenario() {
            Ok(scenario) => scenario,
            Err(Unavailable::FaultNotInThisRevision { fault }) => {
                skipped.push(format!(
                    "{stem}: IDV({fault}) is not in this revision of teprob.f"
                ));
                continue;
            }
        };
        if let Some(seed) = options.seed {
            scenario = scenario.with_seed(seed);
        }
        scenario.driver_forces_idv12 = options.force_idv12;

        let run = Simulation::new(scenario).run();
        let note = match run.outcome {
            Outcome::Completed => String::new(),
            Outcome::Tripped { hours, cause, .. } => {
                format!("  (tripped at {hours:.2} h on {cause:?}, then frozen)")
            }
            Outcome::SolveFailed { step } => {
                skipped.push(format!("{stem}: temperature solve failed at step {step}"));
                continue;
            }
        };

        let path = options
            .out
            .join(format!("{stem}.{}", options.format.extension()));
        let rows = write_file(&path, &run, file.discarded_rows(), options.format)
            .map_err(|e| format!("writing {}: {e}", path.display()))?;
        if rows != file.rows {
            return Err(format!(
                "{stem}: wrote {rows} rows and the published file has {}",
                file.rows
            ));
        }
        eprintln!("tep: {} {rows} x {COLUMNS}{note}", path.display());
        written += 1;
    }

    eprintln!("tep: {written} file(s) in {}", options.out.display());
    for reason in &skipped {
        eprintln!("tep: skipped {reason}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        FIELD_WIDTH, Format, Layout, Options, PUBLISHED_LINE_BYTES, Which, field, published_line,
    };
    use tepsim::published::COLUMNS;

    fn published(stem: &str) -> String {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../reference/data")
            .join(format!("{stem}.dat"));
        std::fs::read_to_string(path).expect("a published file")
    }

    /// Take a shipped line apart into numbers and put it back together with
    /// this module's formatter. The bytes have to come back identical.
    ///
    /// This is the whole layout claim in one assertion. It fails if the field
    /// width is wrong, if the separator is wrong, if the exponent is not
    /// zero-padded to two digits, if the mantissa carries the wrong number of
    /// digits, or if a negative is handled differently from a positive. Each
    /// of those was a plausible reading of the files before this test existed,
    /// and the three-spaces reading was the one actually written first: it
    /// agrees on every positive value and is one byte too wide on every line
    /// that happens to contain a negative.
    #[test]
    fn the_generated_layout_is_the_published_layout() {
        let mut lines = 0_usize;
        let mut negatives = 0_usize;
        for stem in ["d00_te", "d01_te", "d06_te", "d14_te"] {
            let text = published(stem);
            for line in text.lines().filter(|l| !l.trim().is_empty()) {
                let values: Vec<f64> = line
                    .split_whitespace()
                    .map(|f| f.parse().expect("a published field is a number"))
                    .collect();
                assert_eq!(values.len(), COLUMNS, "{stem}: {} fields", values.len());
                negatives += values.iter().filter(|v| **v < 0.0).count();

                let mut row = [0.0; COLUMNS];
                row.copy_from_slice(&values);
                assert_eq!(
                    published_line(&row),
                    line,
                    "{stem} line {lines}: re-emitting the file's own values does \
                     not reproduce its bytes"
                );
                lines += 1;
            }
        }
        assert_eq!(lines, 4 * 960);
        // Without a negative in the sample the test cannot see the difference
        // between fixed-width columns and a three-space separator, which is
        // the error it exists to catch.
        assert!(
            negatives > 0,
            "no negative value in the sample, so this test proves nothing"
        );
        println!("{lines} lines reproduced exactly, {negatives} of them negative");
    }

    #[test]
    fn a_field_is_always_the_column_width_or_narrower() {
        for value in [0.0, 1.0, -1.0, 0.25025, -9.8765e-3, 2704.5, -0.000_012_345] {
            let text = field(value);
            assert!(
                text.len() <= FIELD_WIDTH,
                "{value} formats to {} characters",
                text.len()
            );
        }
        assert_eq!(field(0.250_25), "2.5025000e-01");
        assert_eq!(field(-0.010_865), "-1.0865000e-02");
        assert_eq!(field(0.0), "0.0000000e+00");
        // Five significant digits, and the sixth is gone.
        assert_eq!(field(1.234_567_8), "1.2346000e+00");
    }

    #[test]
    fn a_row_is_exactly_the_published_line_width() {
        let row = [1.0; COLUMNS];
        assert_eq!(published_line(&row).len(), PUBLISHED_LINE_BYTES);
        assert_eq!(PUBLISHED_LINE_BYTES, 832);
    }

    #[test]
    fn the_flags_parse() {
        let args: Vec<String> = ["--set", "all", "--faults", "0,4", "--format", "csv"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        let options = super::parse(&args).expect("valid flags");
        assert_eq!(options.which, Which::All);
        assert_eq!(options.faults, Some(vec![0, 4]));
        assert_eq!(options.format, Format::Csv);

        assert!(super::parse(&["--set".to_string(), "sideways".to_string()]).is_err());
        assert!(super::parse(&["--faults".to_string(), "22".to_string()]).is_err());
        assert!(super::parse(&["--out".to_string()]).is_err());
        assert!(super::parse(&["--nope".to_string()]).is_err());
    }

    #[test]
    fn the_rieth_flags_parse_and_are_refused_when_wrong() {
        let args: Vec<String> = ["--layout", "rieth", "--runs", "500"]
            .iter()
            .map(|s| (*s).to_string())
            .collect();
        let options = super::parse(&args).expect("valid flags");
        assert_eq!(options.layout, Layout::Rieth);
        assert_eq!(options.runs, 500);

        assert!(super::parse(&["--layout".to_string(), "nope".to_string()]).is_err());
        assert!(super::parse(&["--runs".to_string(), "0".to_string()]).is_err());
        assert!(super::parse(&["--runs".to_string(), "many".to_string()]).is_err());
    }

    /// The default layout stays the original distribution, and one run each.
    ///
    /// `--runs` without `--layout rieth` would otherwise look like it did
    /// something.
    #[test]
    fn the_default_layout_is_the_original_distribution() {
        let options = Options::default();
        assert_eq!(options.layout, Layout::Published);
        assert_eq!(options.runs, 1);
    }

    #[test]
    fn the_default_is_the_documented_half_only() {
        let options = Options::default();
        assert_eq!(options.which, Which::Testing);
        assert_eq!(options.format, Format::Dat);
        assert!(
            !options.force_idv12,
            "the published bytes say IDV(12) was not forced"
        );
    }
}
