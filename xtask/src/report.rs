//! Generating `book/src/validation/` from what the validation suite measures.
//!
//! The validation chapter used to be written by hand from `LOG.org`, which made
//! it a transcription of a transcription. Both copies rot, and the rot is
//! invisible: a page that says "worst 6.093e-14" reads exactly the same whether
//! the suite still measures that or has not measured it since B-0026.
//!
//! So the pages here are written by the command that runs the suite, from the
//! suite's own stdout, and from nothing else. Three rules keep that honest:
//!
//! 1. **Verbatim first.** Each chapter carries the tests' output as they
//!    printed it. A parser that summarises can be wrong; a transcript cannot.
//! 2. **Only what ran.** A tier that was not selected gets no chapter written,
//!    and the chapters that do exist each say which commit and which command
//!    produced them. No page carries a number from prose or from memory.
//! 3. **Loud, not clever.** The summary table is built from a block shape the
//!    suite already prints. Where a target prints no such block the page says
//!    so, rather than leaving a row silently blank.
//!
//! The block shape, which `tepsim_oracle::tier1::compare` emits and most of
//! Tiers 1 and 2 use:
//!
//! ```text
//! tier1 kinetics RR
//!   cases          : 9700
//!   max rel err    : 8.492e-16 at perturbed#618[2]
//!   max ulp        : 7 at perturbed#1279[1]
//!   ulp percentiles: p50=0 p90=2 p99=3 p100=7
//! ```
//!
//! A flush-left label, then indented `key : value` lines. The summary table's
//! columns are the union of the keys actually seen, in the order the run first
//! printed them, so a new field in the reporter becomes a new column rather
//! than data this module quietly drops.

use std::fmt::Write as _;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Command, Stdio};

/// Where the generated chapters live, relative to the workspace root.
pub(crate) const DIR: &str = "book/src/validation";

// ---------------------------------------------------------------------------
// running a target and keeping what it said
// ---------------------------------------------------------------------------

/// Which `libm` the port was built against for a run.
///
/// The distinction is load-bearing rather than cosmetic. The vendored `libm`
/// differs from gfortran's by an ULP on about a tenth of `exp` and `pow` calls,
/// so a differential against the default build can only assert 1e-12, and only
/// the platform build can assert bit equality. A page that did not say which
/// one produced a zero would be reporting the wrong claim.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Libm {
    Vendored,
    Platform,
}

impl Libm {
    fn features(self) -> &'static str {
        match self {
            Libm::Vendored => "oracle",
            Libm::Platform => "oracle,libm-system",
        }
    }

    pub(crate) fn label(self) -> &'static str {
        match self {
            Libm::Vendored => "vendored",
            Libm::Platform => "platform",
        }
    }
}

/// libtest's summary line, which every test binary prints in the same shape.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Tally {
    pub(crate) passed: usize,
    pub(crate) failed: usize,
    pub(crate) ignored: usize,
}

/// One test binary's run: the command, the transcript, and libtest's tally.
pub(crate) struct TargetRun {
    /// The test target, e.g. `tier2_balances`.
    pub(crate) target: String,
    /// The exact command, so the page says how to reproduce the numbers.
    pub(crate) command: String,
    pub(crate) libm: Libm,
    /// Every line the run printed, in order, verbatim.
    pub(crate) transcript: Vec<String>,
    pub(crate) tally: Option<Tally>,
}

/// Run one test target, echoing its output and keeping it.
///
/// Streamed rather than collected with `output()`: a Tier 2 target can run for
/// a minute at full volume, and a session watching a silent terminal cannot
/// tell a slow run from a hung one.
pub(crate) fn run_target(
    root: &Path,
    target: &str,
    libm: Libm,
    extra: &[&str],
    env: &[(&str, &str)],
) -> Result<TargetRun, String> {
    let mut args: Vec<String> = ["test", "-p", "tepsim-oracle", "--features"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    args.push(libm.features().to_string());
    // Release always: the full sweeps are minutes in release and hours in
    // debug, and the numbers a page records should come from the configuration
    // anyone would actually reproduce them in.
    args.push("--release".to_string());
    args.push("--test".to_string());
    args.push(target.to_string());
    args.push("--".to_string());
    args.push("--nocapture".to_string());
    // Single-threaded so the transcript is in the order the tests ran. With
    // several threads the printed lines interleave and the page becomes a
    // record of nothing in particular.
    args.push("--test-threads".to_string());
    args.push("1".to_string());
    args.extend(extra.iter().map(|s| s.to_string()));

    let prefix: String = env.iter().map(|(k, v)| format!("{k}={v} ")).collect();
    let command = format!("{prefix}cargo {}", args.join(" "));
    println!("\n$ {command}");

    let mut child = Command::new("cargo")
        .args(&args)
        .current_dir(root)
        .envs(env.iter().copied())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .spawn()
        .map_err(|e| format!("failed to run `{command}`: {e}"))?;

    let stdout = child.stdout.take().expect("piped");
    let mut transcript = Vec::new();
    for line in BufReader::new(stdout).lines().map_while(Result::ok) {
        println!("{line}");
        transcript.push(line);
    }

    let status = child
        .wait()
        .map_err(|e| format!("waiting on `{command}`: {e}"))?;
    if !status.success() {
        // Nothing is written from a red run, deliberately. A page generated
        // from a failing suite would record numbers nobody should quote, and it
        // would look exactly like a page generated from a green one.
        return Err(format!(
            "`{command}` failed.\nThe validation report is generated only from a \
             green run, so no chapter was written."
        ));
    }

    let tally = transcript.iter().find_map(|line| parse_tally(line));
    Ok(TargetRun {
        target: target.to_string(),
        command,
        libm,
        transcript,
        tally,
    })
}

/// `test result: ok. 12 passed; 0 failed; 1 ignored; 0 measured; 0 filtered out`
pub(crate) fn parse_tally(line: &str) -> Option<Tally> {
    let rest = line.trim().strip_prefix("test result:")?;
    let mut tally = Tally {
        passed: 0,
        failed: 0,
        ignored: 0,
    };
    // Every `<count> <word>` pair anywhere in the line, rather than the first
    // word of each clause: the first clause is `ok. 12 passed`, so the count is
    // not where a naive reading puts it, and the whole tally silently came back
    // as zero passed.
    let words: Vec<&str> = rest.split_whitespace().collect();
    for pair in words.windows(2) {
        let Ok(count) = pair[0].parse::<usize>() else {
            continue;
        };
        // The separator rides along on the word: `passed;`, not `passed`.
        match pair[1].trim_end_matches([';', ',', '.']) {
            "passed" => tally.passed = count,
            "failed" => tally.failed = count,
            "ignored" => tally.ignored = count,
            _ => {}
        }
    }
    Some(tally)
}

// ---------------------------------------------------------------------------
// what the tests printed
// ---------------------------------------------------------------------------

/// One `label` + indented `key : value` block lifted out of a transcript.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub(crate) struct Block {
    /// The flush-left label line, e.g. `tier1 kinetics RR`.
    pub(crate) label: String,
    /// The test function it was printed from, when libtest's line said.
    pub(crate) test: String,
    /// The fields, in the order printed. Values are kept as text: reformatting
    /// a number is how a report starts disagreeing with the run it came from.
    pub(crate) fields: Vec<(String, String)>,
}

/// How many words a `key : value` key may have before it is read as prose.
///
/// The widest key the suite actually prints is four words (`worst at the end`),
/// so five leaves room without admitting sentences.
const MAX_KEY_WORDS: usize = 5;

/// Split a transcript into the measurement blocks it contains.
///
/// With `--nocapture` libtest writes `test <name> ... ` without its newline and
/// the test's first printed line lands on the same line. So the prefix is
/// stripped and the remainder is treated as output, which is also what makes
/// the block attributable to a test function.
pub(crate) fn blocks(transcript: &[String]) -> Vec<Block> {
    let mut out: Vec<Block> = Vec::new();
    let mut test = String::new();
    let mut pending_label: Option<String> = None;

    for raw in transcript {
        let mut line = raw.as_str();
        if let Some((name, rest)) = split_libtest_prefix(line) {
            test = name.to_string();
            line = rest;
        }
        if line.trim().is_empty() {
            continue;
        }

        let indented = line.starts_with(' ') || line.starts_with('\t');
        if !indented {
            // A flush-left line either closes a block or offers a label for the
            // next one. libtest's own bookkeeping is neither.
            pending_label = (!is_bookkeeping(line)).then(|| line.trim().to_string());
            continue;
        }

        let Some((key, value)) = line.split_once(':') else {
            continue;
        };
        let (key, value) = (key.trim(), value.trim());
        if key.is_empty() || value.is_empty() {
            continue;
        }
        // A field's value carries a digit. An indented sentence with a colon in
        // it usually does not, and the transcript keeps it either way.
        if !value.chars().any(|c| c.is_ascii_digit()) {
            continue;
        }
        // A key is a label, not a sentence. Without this, an indented line
        // reading `accumulation 0.0613 lb/h against a residual of about
        // 2.10e-11 lb/h: a factor of 2.9e9` splits at its late colon and turns
        // the whole clause into a column heading, which is what the Tier 2 page
        // did before the guard came back. The widest real key is four words.
        if key.split_whitespace().count() > MAX_KEY_WORDS {
            continue;
        }

        match pending_label.take() {
            Some(label) => out.push(Block {
                label,
                test: test.clone(),
                fields: vec![(key.to_string(), value.to_string())],
            }),
            None => match out.last_mut() {
                Some(block) => block.fields.push((key.to_string(), value.to_string())),
                // An indented field with no label above it anywhere. Recorded
                // under the test's name rather than dropped.
                None => out.push(Block {
                    label: format!("(unlabelled, in `{test}`)"),
                    test: test.clone(),
                    fields: vec![(key.to_string(), value.to_string())],
                }),
            },
        }
    }
    out
}

/// `test some_name ... <rest>` split into the name and whatever followed.
fn split_libtest_prefix(line: &str) -> Option<(&str, &str)> {
    let rest = line.strip_prefix("test ")?;
    let (name, tail) = rest.split_once(" ... ")?;
    (!name.contains(char::is_whitespace)).then_some((name, tail))
}

fn is_bookkeeping(line: &str) -> bool {
    let trimmed = line.trim();
    trimmed.starts_with("running ")
        || trimmed.starts_with("test result:")
        || matches!(trimmed, "ok" | "FAILED" | "ignored")
}

// ---------------------------------------------------------------------------
// writing pages
// ---------------------------------------------------------------------------

/// The banner every generated page opens with.
///
/// It names the command and the commit, because a generated page whose
/// provenance is not on the page is indistinguishable from a hand-written one
/// the moment someone reads it in isolation.
pub(crate) fn header(root: &Path, title: &str, command: &str, claim: &str) -> String {
    let commit = describe_commit(root);
    format!(
        "<!-- {PROVENANCE} `{command}` from commit `{commit}`. Do not edit by \
         hand: the next run overwrites it. -->\n\n\
         # {title}\n\n\
         > **This page is generated.** `{command}` wrote it from commit \
         `{commit}`.\n\
         > {claim}\n"
    )
}

/// What a generated validation chapter claims about its own numbers.
pub(crate) const MEASURED: &str = "Every number on it was captured from that \
    run's own output. To change what it says,\n> change what the suite measures \
    and run the command again.";

/// The marker the provenance comment opens with, so a page can be read back.
const PROVENANCE: &str = "GENERATED by";

/// The command and commit a generated page records, read back from the page.
///
/// The index reports what each chapter on disk says about itself rather than
/// what this run did, because chapters are written a tier at a time and one
/// index claiming they all came from the same run would be false.
pub(crate) fn read_provenance(root: &Path, relative: &str) -> Option<(String, String)> {
    let text = fs::read_to_string(root.join(relative)).ok()?;
    let line = text.lines().next()?;
    let rest = line.split_once(PROVENANCE)?.1;
    let mut quoted = rest.split('`').skip(1).step_by(2);
    Some((quoted.next()?.to_string(), quoted.next()?.to_string()))
}

/// `HEAD`'s short hash, with `-dirty` when the tree has uncommitted changes.
///
/// The dirty marker matters: a page generated from a modified tree records
/// numbers that no commit reproduces, and saying so is the difference between
/// a provenance line and a decoration.
pub(crate) fn describe_commit(root: &Path) -> String {
    let hash = Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .current_dir(root)
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty());
    let Some(hash) = hash else {
        return "unknown (not a git checkout)".to_string();
    };
    let dirty = Command::new("git")
        .args(["status", "--porcelain"])
        .current_dir(root)
        .output()
        .ok()
        .is_some_and(|o| !o.stdout.is_empty());
    if dirty { format!("{hash}-dirty") } else { hash }
}

pub(crate) fn write_generated(root: &Path, relative: &str, body: &str) -> Result<(), String> {
    let path = root.join(relative);
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| format!("creating {}: {e}", parent.display()))?;
    }
    fs::write(&path, body).map_err(|e| format!("writing {}: {e}", path.display()))?;
    println!("[write] {relative}");
    Ok(())
}

/// A fenced transcript, with the fence widened if the text contains one.
pub(crate) fn fenced(body: &[String]) -> String {
    let mut fence = "```".to_string();
    while body
        .iter()
        .any(|line| line.trim_start().starts_with(&fence))
    {
        fence.push('`');
    }
    let mut out = format!("{fence}text\n");
    for line in body {
        out.push_str(line);
        out.push('\n');
    }
    out.push_str(&fence);
    out.push('\n');
    out
}

/// The toolchain the numbers were produced with.
///
/// On every page, because the Tier 1 and Tier 2 numbers are only meaningful
/// against a fixed compiler and fixed flags. `CLAUDE.md` treats a `gfortran`
/// change as a logged re-baseline rather than a regression hunt, and that is
/// only possible if the pages say which one they used.
pub(crate) fn toolchain() -> String {
    let rustc = tool_version("rustc", &["--version"]).unwrap_or_else(|| "rustc unknown".into());
    let gfortran = tool_version("gfortran", &["-dumpfullversion"])
        .map(|v| format!("gfortran {v}"))
        .unwrap_or_else(|| "no gfortran".into());
    format!("{rustc}, {gfortran}")
}

fn tool_version(program: &str, args: &[&str]) -> Option<String> {
    let output = Command::new(program).args(args).output().ok()?;
    let text = String::from_utf8(output.stdout).ok()?;
    let first = text.lines().next()?.trim().to_string();
    (!first.is_empty()).then_some(first)
}

/// Trim the blank lines a transcript opens and closes with.
///
/// Nothing else is removed. libtest's own lines stay, because they are part of
/// the record and they are what the tally table above is derived from: a reader
/// who distrusts the table can check it against the transcript on the same page.
pub(crate) fn trimmed(transcript: &[String]) -> Vec<String> {
    let start = transcript
        .iter()
        .position(|l| !l.trim().is_empty())
        .unwrap_or(transcript.len());
    let end = transcript
        .iter()
        .rposition(|l| !l.trim().is_empty())
        .map_or(start, |i| i + 1);
    transcript[start..end].to_vec()
}

/// Render one tier's chapter from the targets that ran.
pub(crate) fn render_tier(
    root: &Path,
    tier: u8,
    title: &str,
    lead: &str,
    command: &str,
    runs: &[TargetRun],
) -> String {
    let mut page = header(root, &format!("Tier {tier}: {title}"), command, MEASURED);
    let _ = write!(page, "\n{lead}\n\n");
    let _ = writeln!(
        page,
        "Produced with {}. The oracle's compiler flags are fixed in\n\
         `crates/tepsim-oracle/build.rs` and asserted by a test; changing them \
         invalidates\nevery number on this page, which is why it is a logged \
         re-baseline and not an edit.\n",
        toolchain()
    );

    let (passed, failed, ignored) = runs.iter().fold((0, 0, 0), |acc, run| {
        run.tally.map_or(acc, |t| {
            (acc.0 + t.passed, acc.1 + t.failed, acc.2 + t.ignored)
        })
    });
    let _ = writeln!(
        page,
        "## What ran\n\n\
         {} test binar{}: **{passed} test(s) passed**, {failed} failed, \
         {ignored} ignored.\n",
        runs.len(),
        if runs.len() == 1 { "y" } else { "ies" }
    );
    page.push_str("| target | `libm` | passed | failed | ignored |\n|---|---|---|---|---|\n");
    for run in runs {
        let (p, f, i) = run.tally.map_or(("?".into(), "?".into(), "?".into()), |t| {
            (
                t.passed.to_string(),
                t.failed.to_string(),
                t.ignored.to_string(),
            )
        });
        let _ = writeln!(
            page,
            "| `{}` | {} | {p} | {f} | {i} |",
            run.target,
            run.libm.label()
        );
    }

    page.push_str(&render_summary(runs));

    page.push_str(
        "\n## Transcripts\n\n\
         Each block below is a test binary's own output, verbatim, with the \
         command that\nproduced it. The summary above is derived from these; \
         they are not derived from it.\n",
    );
    for run in runs {
        let _ = writeln!(
            page,
            "\n### `{}`, {} `libm`\n\n```text\n{}\n```\n",
            run.target,
            run.libm.label(),
            run.command
        );
        page.push_str(&fenced(&trimmed(&run.transcript)));
    }
    page
}

/// The summary table: one row per measurement block, columns from the keys.
///
/// The columns are the union of the keys the run actually printed, in
/// first-seen order, so a field added to the reporter shows up as a new column
/// instead of being silently dropped by a hard-coded list.
fn render_summary(runs: &[TargetRun]) -> String {
    let mut rows: Vec<(&str, Block)> = Vec::new();
    for run in runs {
        for block in blocks(&run.transcript) {
            rows.push((run.target.as_str(), block));
        }
    }
    if rows.is_empty() {
        return "\n## Measurements\n\nNo target in this tier printed the \
                `label` then indented `key : value` shape the\nsummary table is \
                built from, so there is no table: the transcripts below are \
                the\nmeasurement. That is stated rather than left as an empty \
                table, because an empty\ntable and a table nobody filled in look \
                the same.\n"
            .to_string();
    }

    let mut keys: Vec<&str> = Vec::new();
    for (_, block) in &rows {
        for (key, _) in &block.fields {
            if !keys.contains(&key.as_str()) {
                keys.push(key);
            }
        }
    }

    let mut out = String::from("\n## Measurements\n\n");
    let _ = write!(
        out,
        "\
         {} block(s), lifted from the transcripts below. The columns are \
         whatever fields the\nrun printed, so a new field in the reporter \
         becomes a new column here rather than\ndata this page drops.\n\n\
         The test column matters as much as the numbers, because not every row \
         is the port\nbeing measured. Some tests deliberately mis-type a \
         constant, or solve from the wrong\nguess, to show what that would cost; \
         a row from one of those is supposed to be\nenormous, and its test name \
         says so. The `what` column carries a `tier1` prefix in\nevery tier, \
         because it is the shared comparison reporter's own label rather than \
         a\nclaim about which tier printed it.\n\n",
        rows.len()
    );
    out.push_str("| target | from test | what | ");
    out.push_str(&keys.join(" | "));
    out.push_str(" |\n|---|---|---|");
    for _ in &keys {
        out.push_str("---|");
    }
    out.push('\n');
    for (target, block) in &rows {
        let test = if block.test.is_empty() {
            "-".to_string()
        } else {
            format!("`{}`", escape(&block.test))
        };
        let _ = write!(out, "| `{target}` | {test} | {} |", escape(&block.label));
        for key in &keys {
            let value = block
                .fields
                .iter()
                .find(|(k, _)| k == key)
                .map_or(String::new(), |(_, v)| format!("`{}`", escape(v)));
            let _ = write!(out, " {value} |");
        }
        out.push('\n');
    }
    out
}

/// A pipe inside a cell ends the cell. Nothing else needs escaping here.
fn escape(text: &str) -> String {
    text.replace('|', "\\|")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lines(text: &str) -> Vec<String> {
        text.lines().map(str::to_string).collect()
    }

    #[test]
    fn reads_libtests_tally() {
        let line = "test result: ok. 12 passed; 0 failed; 1 ignored; 0 measured; \
                    0 filtered out; finished in 3.21s";
        assert_eq!(
            parse_tally(line),
            Some(Tally {
                passed: 12,
                failed: 0,
                ignored: 1
            })
        );
    }

    /// The count is not the first word of its clause: libtest writes
    /// `ok. 12 passed`, and reading the first word gave zero passed on every
    /// run, which the summary table reported as a straight face.
    #[test]
    fn the_count_is_found_past_the_verdict_word() {
        assert_eq!(
            parse_tally("test result: FAILED. 3 passed; 1 failed; 0 ignored"),
            Some(Tally {
                passed: 3,
                failed: 1,
                ignored: 0
            })
        );
    }

    #[test]
    fn a_line_that_is_not_a_tally_is_not_one() {
        assert_eq!(parse_tally("  worst result: 3e-14"), None);
        assert_eq!(parse_tally("running 4 tests"), None);
    }

    /// The index reads each chapter's provenance back off the page, so the
    /// comment it writes has to be the comment it can parse.
    #[test]
    fn the_provenance_comment_round_trips() {
        let root = std::env::temp_dir().join(format!("xtask-report-{}", std::process::id()));
        let _ = fs::remove_dir_all(&root);
        let page = header(
            &root,
            "Tier 1: the utility routines",
            "cargo xtask validate --tiers 1",
            MEASURED,
        );
        write_generated(&root, "book/src/validation/tier1.md", &page).expect("write");

        let (by, commit) =
            read_provenance(&root, "book/src/validation/tier1.md").expect("provenance");
        assert_eq!(by, "cargo xtask validate --tiers 1");
        assert!(!commit.is_empty(), "the commit field must not be blank");
        assert!(page.contains(&commit), "{page}");

        assert_eq!(read_provenance(&root, "book/src/validation/tier9.md"), None);
        let _ = fs::remove_dir_all(&root);
    }

    /// The real shape, copied from a `tier2_kinetics` run: libtest's prefix
    /// with no newline, so the first printed line is glued to it.
    #[test]
    fn lifts_blocks_out_of_a_real_transcript() {
        let transcript = lines(
            "running 3 tests\n\
             test the_kinetics_match ... exp and pow come from the vendored libm\n\
             tier1 kinetics RR\n\
             \x20 cases          : 9700\n\
             \x20 max rel err    : 8.492e-16 at perturbed#618[2]\n\
             \x20 max ulp        : 7 at perturbed#1279[1]\n\
             tier1 kinetics CRXR\n\
             \x20 cases          : 19400\n\
             \x20 max rel err    : 8.492e-16 at perturbed#618[8]\n\
             ok\n\
             \n\
             test result: ok. 3 passed; 0 failed; 0 ignored\n",
        );
        let found = blocks(&transcript);
        assert_eq!(found.len(), 2, "{found:#?}");
        assert_eq!(found[0].label, "tier1 kinetics RR");
        assert_eq!(found[0].test, "the_kinetics_match");
        assert_eq!(found[0].fields.len(), 3);
        assert_eq!(found[0].fields[0], ("cases".into(), "9700".into()));
        assert_eq!(
            found[0].fields[1],
            ("max rel err".into(), "8.492e-16 at perturbed#618[2]".into())
        );
        assert_eq!(found[1].label, "tier1 kinetics CRXR");
        assert_eq!(found[1].fields.len(), 2);
    }

    /// libtest's own lines must never become a block label: `running 3 tests`
    /// sits flush left immediately above the first field.
    #[test]
    fn libtest_bookkeeping_is_never_a_label() {
        let transcript = lines(
            "running 1 test\n\
             \x20 cases : 5\n\
             ok\n\
             test result: ok. 1 passed; 0 failed; 0 ignored\n",
        );
        let found = blocks(&transcript);
        assert_eq!(found.len(), 1);
        assert!(found[0].label.starts_with("(unlabelled"), "{found:#?}");
    }

    /// A tier whose tests print prose rather than blocks yields none, and the
    /// page then says so instead of showing an empty table.
    #[test]
    fn prose_output_yields_no_blocks() {
        let transcript = lines(
            "running 2 tests\n\
             test both_scalings_appear ... 30 signed draws, 432 unit draws\n\
             ok\n\
             test the_capacity_has_headroom ... worst evaluation: 462 draws\n\
             ok\n\
             test result: ok. 2 passed; 0 failed; 0 ignored\n",
        );
        assert!(blocks(&transcript).is_empty());
    }

    /// The exact line from `tier5_invariants` that turned a whole clause into
    /// a column heading on the Tier 2 page. Its colon falls near the end, so
    /// the "key" was thirteen words of prose.
    #[test]
    fn an_indented_sentence_with_a_late_colon_is_not_a_field() {
        let transcript = lines(
            "a label\n\
             \x20 accumulation 0.0613 lb/h against a residual of about \
             2.10e-11 lb/h: a factor of 2.9e9\n\
             \x20 worst at the end : 7.882e-11\n",
        );
        let found = blocks(&transcript);
        assert_eq!(found.len(), 1);
        // The four-word key survives; the thirteen-word sentence does not.
        assert_eq!(
            found[0].fields,
            vec![("worst at the end".into(), "7.882e-11".into())]
        );
    }

    #[test]
    fn an_indented_sentence_without_a_number_is_not_a_field() {
        let transcript = lines(
            "a label\n\
             \x20 every line of this file is claimed by some function\n\
             \x20 cases : 5\n",
        );
        let found = blocks(&transcript);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].fields, vec![("cases".into(), "5".into())]);
    }

    #[test]
    fn the_summary_columns_are_the_union_of_the_keys_printed() {
        let run = TargetRun {
            target: "tier2_kinetics".into(),
            command: "cargo test".into(),
            libm: Libm::Vendored,
            transcript: lines(
                "one\n\x20 cases : 5\n\x20 max ulp : 1\ntwo\n\x20 cases : 7\n\x20 skipped : 2\n",
            ),
            tally: None,
        };
        let table = render_summary(&[run]);
        assert!(
            table.contains("| target | from test | what | cases | max ulp | skipped |"),
            "{table}"
        );
        // The block that never printed `max ulp` leaves that cell empty rather
        // than borrowing the other block's value.
        assert!(table.contains("| two | `7` |  | `2` |"), "{table}");
    }

    #[test]
    fn a_pipe_in_a_value_does_not_break_the_table() {
        assert_eq!(escape("a|b"), "a\\|b");
    }

    /// A test that prints a fenced block of its own must not break the page.
    #[test]
    fn a_fence_in_the_output_widens_the_fence() {
        let body = vec!["```".to_string(), "inner".to_string()];
        let out = fenced(&body);
        assert!(out.starts_with("````text\n"), "{out}");
        assert!(out.trim_end().ends_with("````"), "{out}");
    }

    #[test]
    fn a_plain_transcript_uses_a_plain_fence() {
        assert!(fenced(&["a".to_string()]).starts_with("```text\n"));
    }

    #[test]
    fn trimming_removes_only_the_surrounding_blank_lines() {
        let transcript = lines("\n\na\n\nb\n\n");
        assert_eq!(trimmed(&transcript), vec!["a", "", "b"]);
        assert!(trimmed(&lines("\n\n")).is_empty());
    }
}
