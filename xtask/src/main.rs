//! Build, gate and validation automation for the tep-rust workspace.
//!
//! Deliberately dependency-free: this runs at the top of every development
//! session, so its own compile time is a tax on every iteration.
//!
//! # Subcommands
//!
//! - `ci [--fast]` is *the* gate. Everything else is advisory; this is not.
//! - `provenance` reports which line ranges of the original `teprob.f` no Rust
//!   function claims. It is how a silently dropped term gets caught, since no
//!   differential test can find a term that is never evaluated.
//! - `fidelity` diffs a short run against a committed golden oracle trace. It is
//!   a stub until B-0004 produces that trace, and it fails loudly rather than
//!   passing vacuously.
//!
//! See `CLAUDE.md` for how these fit into the session protocol.

#![forbid(unsafe_code)]

use std::fmt;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitCode};

/// Crates that are published or shipped, and therefore must never depend on the
/// development-only Fortran oracle.
const SHIPPED_CRATES: &[&str] = &[
    "tepsim",
    "tepsim-core",
    "tepsim-control",
    "tepsim-scenario",
    "tepsim-cli",
    "tepsim-py",
    "tepsim-wasm",
];

/// The marker that opens a provenance claim. Must be anchored at the start
/// of a comment line so prose and test fixtures cannot masquerade as coverage.
const PROVENANCE_MARKER: &str = "@port";

/// The file a claim refers to, following the marker.
const PROVENANCE_TAG: &str = "teprob.f:";

/// Path to the vendored original, relative to the workspace root.
const REFERENCE_FORTRAN: &str = "reference/fortran/teprob.f";

/// Path to the committed golden trace, relative to the workspace root.
const GOLDEN_TRACE: &str = "reference/golden/nominal-100-steps.json";

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (cmd, flags) = match args.split_first() {
        Some((c, rest)) => (c.as_str(), rest),
        None => {
            usage();
            return ExitCode::from(2);
        }
    };

    let root = workspace_root();
    let result = match cmd {
        "ci" => cmd_ci(&root, flags.iter().any(|f| f == "--fast")),
        "provenance" => cmd_provenance(&root),
        "fidelity" => cmd_fidelity(&root),
        "help" | "--help" | "-h" => {
            usage();
            return ExitCode::SUCCESS;
        }
        other => Err(format!("unknown subcommand `{other}`")),
    };

    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("\nxtask: {msg}");
            ExitCode::from(2)
        }
    }
}

fn usage() {
    eprintln!(
        "\
usage: cargo xtask <command>

  ci [--fast]   the gate: fmt, clippy, test, doc, deny, oracle isolation
                --fast skips the Fortran oracle differential job
  provenance    teprob.f line ranges not claimed by any Rust function
  fidelity      diff a short run against the committed golden oracle trace
  help          this message"
    );
}

/// The workspace root, derived from this crate's manifest directory.
fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .to_path_buf()
}

// ---------------------------------------------------------------------------
// ci
// ---------------------------------------------------------------------------

fn cmd_ci(root: &Path, fast: bool) -> Result<(), String> {
    check_oracle_isolation(root)?;

    step(root, "cargo", &["fmt", "--all", "--check"])?;
    step(
        root,
        "cargo",
        &[
            "clippy",
            "--workspace",
            "--all-targets",
            "--",
            "-D",
            "warnings",
        ],
    )?;
    step(root, "cargo", &["test", "--workspace"])?;
    step_with_env(
        root,
        "cargo",
        &["doc", "--workspace", "--no-deps"],
        &[("RUSTDOCFLAGS", "-D warnings")],
    )?;
    step(root, "cargo", &["deny", "check"])?;

    if fast {
        println!("\n[skip] oracle differential job (--fast)");
    } else if which("gfortran").is_none() {
        println!(
            "\n[skip] oracle differential job: gfortran not found.\n\
             Per CLAUDE.md, a session without gfortran must not do model work."
        );
    } else {
        step(
            root,
            "cargo",
            &["test", "-p", "tepsim-oracle", "--features", "oracle"],
        )?;
    }

    println!("\nci: green");
    Ok(())
}

/// The oracle links the original Fortran and must never reach a shipped
/// artifact. Cheap to check, expensive to discover at publish time.
fn check_oracle_isolation(root: &Path) -> Result<(), String> {
    for krate in SHIPPED_CRATES {
        let manifest = root.join("crates").join(krate).join("Cargo.toml");
        let text = fs::read_to_string(&manifest)
            .map_err(|e| format!("reading {}: {e}", manifest.display()))?;
        if text.contains("tepsim-oracle") {
            return Err(format!(
                "{krate} references tepsim-oracle. The oracle is development-only \
                 and must never be reachable from a shipped crate."
            ));
        }
    }
    println!("[ok] oracle isolation: no shipped crate depends on tepsim-oracle");
    Ok(())
}

fn step(root: &Path, program: &str, args: &[&str]) -> Result<(), String> {
    step_with_env(root, program, args, &[])
}

fn step_with_env(
    root: &Path,
    program: &str,
    args: &[&str],
    env: &[(&str, &str)],
) -> Result<(), String> {
    println!("\n$ {program} {}", args.join(" "));
    let mut cmd = Command::new(program);
    cmd.args(args).current_dir(root);
    for (k, v) in env {
        cmd.env(k, v);
    }
    let status = cmd
        .status()
        .map_err(|e| format!("failed to run `{program}`: {e}"))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("`{program} {}` failed", args.join(" ")))
    }
}

fn which(program: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    std::env::split_paths(&path)
        .map(|dir| dir.join(program))
        .find(|candidate| candidate.is_file())
}

// ---------------------------------------------------------------------------
// provenance
// ---------------------------------------------------------------------------

/// An inclusive, one-based line range in the original Fortran.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
struct LineRange {
    start: usize,
    end: usize,
}

impl LineRange {
    fn len(self) -> usize {
        self.end.saturating_sub(self.start) + 1
    }
}

impl fmt::Display for LineRange {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.start == self.end {
            write!(f, "{}", self.start)
        } else {
            write!(f, "{}-{}", self.start, self.end)
        }
    }
}

fn cmd_provenance(root: &Path) -> Result<(), String> {
    let fortran = root.join(REFERENCE_FORTRAN);
    let Ok(source) = fs::read_to_string(&fortran) else {
        return Err(format!(
            "{REFERENCE_FORTRAN} not present, so coverage cannot be computed.\n\
             Vendor the reference material first (B-0003). Failing rather than \
             reporting a vacuous all-clear."
        ));
    };
    let total = source.lines().count();

    let mut claimed = Vec::new();
    let mut files_with_claims = 0usize;
    for file in rust_sources(root)? {
        let text =
            fs::read_to_string(&file).map_err(|e| format!("reading {}: {e}", file.display()))?;
        let found = parse_annotations(&text);
        if !found.is_empty() {
            files_with_claims += 1;
        }
        claimed.extend(found);
    }

    let merged = merge(claimed);
    let claimed_lines: usize = merged.iter().map(|r| r.len()).sum();
    let unclaimed = gaps(&merged, total);

    println!("provenance against {REFERENCE_FORTRAN}");
    println!("  total lines:     {total}");
    println!(
        "  claimed:         {claimed_lines} ({:.1}%) across {files_with_claims} file(s)",
        percent(claimed_lines, total)
    );
    println!("  unclaimed spans: {}", unclaimed.len());
    for range in &unclaimed {
        println!("    {range}  ({} lines)", range.len());
    }
    if unclaimed.is_empty() {
        println!("  every line of the original is claimed by some Rust function");
    }
    Ok(())
}

fn percent(part: usize, whole: usize) -> f64 {
    if whole == 0 {
        0.0
    } else {
        100.0 * part as f64 / whole as f64
    }
}

/// Every `.rs` file in the workspace, skipping build output.
fn rust_sources(root: &Path) -> Result<Vec<PathBuf>, String> {
    let mut out = Vec::new();
    // Only the porting crates. xtask is tooling; scanning it is exactly how
    // the tool came to count its own tests as coverage.
    let mut stack = vec![root.join("crates")];
    while let Some(dir) = stack.pop() {
        if !dir.exists() {
            continue;
        }
        let entries = fs::read_dir(&dir).map_err(|e| format!("reading {}: {e}", dir.display()))?;
        for entry in entries {
            let entry = entry.map_err(|e| format!("reading {}: {e}", dir.display()))?;
            let path = entry.path();
            if path.is_dir() {
                if path.file_name().is_some_and(|n| n == "target") {
                    continue;
                }
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "rs") {
                out.push(path);
            }
        }
    }
    out.sort();
    Ok(out)
}

/// Extract every provenance claim from a source file.
///
/// A claim is one whole comment line. Anchoring is the whole point: an earlier
/// version scanned for the bare `teprob.f:` tag anywhere in the text and
/// happily counted its own doc comments and test fixtures, reporting 2.6%
/// coverage of a file that nothing had ported. A tool that overstates coverage
/// is worse than no tool, so a claim must look exactly like a claim.
fn parse_annotations(text: &str) -> Vec<LineRange> {
    text.lines().filter_map(parse_claim_line).collect()
}

/// Parse one line as `// @port teprob.f:505-522`, in any comment style.
///
/// Returns `None` for prose that merely mentions the convention, for string
/// fixtures, and for a bare tag with no marker.
fn parse_claim_line(line: &str) -> Option<LineRange> {
    let body = strip_comment_prefix(line.trim_start())?;
    let after_marker = body.trim_start().strip_prefix(PROVENANCE_MARKER)?;
    // Require whitespace after the marker so `@ported` is not a claim.
    if !after_marker.starts_with(char::is_whitespace) {
        return None;
    }
    let after_tag = after_marker.trim_start().strip_prefix(PROVENANCE_TAG)?;
    let (start, tail) = take_usize(after_tag)?;
    let end = match tail.strip_prefix('-').and_then(take_usize) {
        Some((e, _)) => e,
        None => start,
    };
    Some(LineRange {
        start,
        end: end.max(start),
    })
}

/// Strip `//`, `///` or `//!` from the front of an already-trimmed line.
fn strip_comment_prefix(s: &str) -> Option<&str> {
    Some(s.strip_prefix("//")?.trim_start_matches(['/', '!']))
}

fn take_usize(s: &str) -> Option<(usize, &str)> {
    let end = s.find(|c: char| !c.is_ascii_digit()).unwrap_or(s.len());
    if end == 0 {
        return None;
    }
    s[..end].parse().ok().map(|n| (n, &s[end..]))
}

/// Sort and coalesce ranges, joining ones that touch or abut.
fn merge(mut ranges: Vec<LineRange>) -> Vec<LineRange> {
    ranges.sort_unstable();
    let mut out: Vec<LineRange> = Vec::with_capacity(ranges.len());
    for range in ranges {
        match out.last_mut() {
            Some(last) if range.start <= last.end.saturating_add(1) => {
                last.end = last.end.max(range.end);
            }
            _ => out.push(range),
        }
    }
    out
}

/// The complement of `merged` within `1..=total`.
fn gaps(merged: &[LineRange], total: usize) -> Vec<LineRange> {
    let mut out = Vec::new();
    let mut cursor = 1usize;
    for range in merged {
        if range.start > cursor {
            out.push(LineRange {
                start: cursor,
                end: range.start - 1,
            });
        }
        cursor = cursor.max(range.end.saturating_add(1));
    }
    if cursor <= total {
        out.push(LineRange {
            start: cursor,
            end: total,
        });
    }
    out
}

// ---------------------------------------------------------------------------
// fidelity
// ---------------------------------------------------------------------------

fn cmd_fidelity(root: &Path) -> Result<(), String> {
    if !root.join(GOLDEN_TRACE).exists() {
        return Err(format!(
            "no golden trace at {GOLDEN_TRACE}.\n\
             The fidelity preflight is not available until B-0004 generates one \
             from the Fortran oracle. Failing deliberately: a preflight that \
             passed with nothing to check would be worse than no preflight."
        ));
    }
    Err("golden trace exists but the comparison is not implemented (B-0004)".into())
}

// ---------------------------------------------------------------------------
// tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn r(start: usize, end: usize) -> LineRange {
        LineRange { start, end }
    }

    #[test]
    fn parses_a_line_comment_claim() {
        assert_eq!(
            parse_annotations("// @port teprob.f:505-522"),
            vec![r(505, 522)]
        );
    }

    #[test]
    fn parses_doc_and_module_doc_claims() {
        assert_eq!(parse_annotations("/// @port teprob.f:1"), vec![r(1, 1)]);
        assert_eq!(
            parse_annotations("//! @port teprob.f:1552-1560"),
            vec![r(1552, 1560)]
        );
    }

    #[test]
    fn indentation_and_trailing_text_are_allowed() {
        assert_eq!(
            parse_annotations("        // @port teprob.f:100-200 (TEFUNC prologue)"),
            vec![r(100, 200)]
        );
    }

    #[test]
    fn parses_several_across_lines() {
        let text = "// @port teprob.f:1-10\ncode();\n/// @port teprob.f:20-30\n";
        assert_eq!(parse_annotations(text), vec![r(1, 10), r(20, 30)]);
    }

    // The next four are the regression tests for the false-positive bug: the
    // first version of this scanner reported 2.6% coverage of a file nothing
    // had touched, by counting the tool's own documentation and fixtures.

    #[test]
    fn prose_mentioning_the_convention_is_not_a_claim() {
        let text = "//! Annotations look like `@port teprob.f:505-522`.";
        assert_eq!(parse_annotations(text), vec![]);
    }

    #[test]
    fn a_string_fixture_is_not_a_claim() {
        let text = r#"assert_eq!(parse("@port teprob.f:1-10"), vec![]);"#;
        assert_eq!(parse_annotations(text), vec![]);
    }

    #[test]
    fn a_bare_tag_without_the_marker_is_not_a_claim() {
        assert_eq!(parse_annotations("// see teprob.f:505-522"), vec![]);
    }

    #[test]
    fn a_marker_glued_to_a_word_is_not_a_claim() {
        assert_eq!(parse_annotations("// @ported teprob.f:1-10"), vec![]);
    }

    #[test]
    fn ignores_a_claim_with_no_number() {
        assert_eq!(parse_annotations("// @port teprob.f:foo"), vec![]);
    }

    #[test]
    fn a_reversed_range_is_clamped_not_panicked() {
        assert_eq!(
            parse_annotations("// @port teprob.f:90-10"),
            vec![r(90, 90)]
        );
    }

    #[test]
    fn merge_coalesces_overlapping_and_abutting() {
        let merged = merge(vec![r(10, 20), r(15, 25), r(26, 30), r(40, 50)]);
        assert_eq!(merged, vec![r(10, 30), r(40, 50)]);
    }

    #[test]
    fn merge_handles_full_containment() {
        assert_eq!(merge(vec![r(1, 100), r(20, 30)]), vec![r(1, 100)]);
    }

    #[test]
    fn merge_of_nothing_is_nothing() {
        assert_eq!(merge(vec![]), vec![]);
    }

    #[test]
    fn gaps_finds_head_middle_and_tail() {
        let merged = merge(vec![r(10, 20), r(30, 40)]);
        assert_eq!(gaps(&merged, 50), vec![r(1, 9), r(21, 29), r(41, 50)]);
    }

    #[test]
    fn nothing_claimed_means_the_whole_file_is_a_gap() {
        assert_eq!(gaps(&[], 1594), vec![r(1, 1594)]);
    }

    #[test]
    fn full_coverage_leaves_no_gaps() {
        assert_eq!(gaps(&[r(1, 1594)], 1594), vec![]);
    }

    #[test]
    fn claims_past_the_end_do_not_wrap_or_panic() {
        assert_eq!(gaps(&[r(1, 5000)], 1594), vec![]);
    }

    #[test]
    fn range_length_is_inclusive() {
        assert_eq!(r(1, 1).len(), 1);
        assert_eq!(r(10, 20).len(), 11);
    }

    #[test]
    fn percent_of_an_empty_file_is_zero_not_nan() {
        // 0/0 would be NaN, which formats as "NaN%" and silently poisons any
        // arithmetic downstream. Compared by bits so the assertion is exact
        // without an approximate-equality fudge.
        let p = percent(0, 0);
        assert!(p.is_finite(), "percent of an empty file must not be NaN");
        assert_eq!(p.to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn percent_is_a_plain_ratio() {
        let p = percent(1, 4);
        assert_eq!(p.to_bits(), 25.0_f64.to_bits());
    }
}
