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

use tepsim_oracle::golden::{self, Trace};

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
/// The vendored files a `@port` claim can name.
///
/// A claim is `@port <file>:<range>`, and the file must be one of these. An
/// unrecognised name is *not* silently ignored: `cmd_provenance` reports it,
/// because a typo in a claim is a claim that stops counting and nothing else
/// would notice.
const PROVENANCE_FILES: &[(&str, &str)] = &[
    ("teprob.f:", "reference/fortran/teprob.f"),
    ("temain_mod.f:", "reference/fortran/temain_mod.f"),
];

/// The pinned toolchain manifest, relative to the workspace root.
const TOOLCHAIN_FILE: &str = "rust-toolchain.toml";

/// Path to the committed golden trace, relative to the workspace root.
const GOLDEN_TRACE: &str = tepsim_oracle::golden::PATH;

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
        "validate" => cmd_validate(&root, flags),
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
  validate [--tiers 1,2,3] [--compare-to-log]
                run the validation ladder at full volume, in release
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

/// Test files re-run with the port on the platform libm, asserting 0 ULP.
///
/// Every Phase 2 differential belongs here. See `tepsim_core::math`: the
/// vendored libm differs from gfortran's on about a tenth of `exp` and `pow`
/// calls, so the default run can only assert 1e-12, and this run is what still
/// holds the algebra to bit equality. Append the new file when an item lands.
/// The Tier 2 differentials, in the order the model computes them.
///
/// A failure in an early one explains failures in the later ones, so running
/// them in this order means the first red line is the informative one.
const TIER2_TESTS: &[&str] = &[
    "tier2_unpack",
    "tier2_equilibrium",
    "tier2_kinetics",
    "tier2_streams",
    "tier2_flows",
    "tier2_stripper",
    "tier2_heat",
    "tier2_measurements",
    "tier2_balances",
];

/// The Tier 3 differentials, in the order the stream is consumed.
const TIER3_TESTS: &[&str] = &[
    "rng_call_order",
    "tier3_harness",
    "tier1_disturbance",
    "tier3_walk",
    "tier3_walk_inputs",
    "tier3_analysers",
    "fault_table",
];

const LIBM_SYSTEM_TESTS: &[&str] = &[
    "tier4_closed_loop",
    "tier5_invariants",
    "tier5_runs",
    "tier2_equilibrium",
    "tier2_kinetics",
    "tier2_streams",
    "tier2_flows",
    "tier2_stripper",
    "tier2_heat",
    "tier2_measurements",
    "tier2_balances",
];

fn cmd_ci(root: &Path, fast: bool) -> Result<(), String> {
    check_toolchain(root)?;
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
        // The same Tier 2 differential again, with the port using the platform
        // libm instead of the vendored one. That is not a shipping
        // configuration; it exists so the comparison can be made with the
        // transcendental removed, leaving bit equality as the assertion about
        // the algebra. Without it, every item from B-0018 onward would be held
        // only to 1e-12, which is four orders of magnitude of room to hide a
        // reassociation in. Scoped to the one test that needs it, since it
        // rebuilds the workspace under a different feature set.
        let mut bit_exact = vec![
            "test",
            "-p",
            "tepsim-oracle",
            "--features",
            "oracle,libm-system",
        ];
        for name in LIBM_SYSTEM_TESTS {
            bit_exact.push("--test");
            bit_exact.push(name);
        }
        step(root, "cargo", &bit_exact)?;
    }

    println!("\nci: green");
    Ok(())
}

/// The running compiler must be the one `rust-toolchain.toml` pins.
///
/// `rust-toolchain.toml` only takes effect when `cargo` is rustup's. On a
/// machine whose PATH reaches a distribution or Homebrew toolchain first, the
/// file is *silently ignored* and every validation number is produced by an
/// unknown compiler. That is precisely the kind of quiet mismatch this project
/// refuses to tolerate, so it is checked rather than assumed.
fn check_toolchain(root: &Path) -> Result<(), String> {
    let manifest = root.join(TOOLCHAIN_FILE);
    let text =
        fs::read_to_string(&manifest).map_err(|e| format!("reading {TOOLCHAIN_FILE}: {e}"))?;

    let pinned = parse_pinned_channel(&text)
        .ok_or_else(|| format!("{TOOLCHAIN_FILE} has no `channel` entry"))?;

    let output = Command::new("rustc")
        .arg("--version")
        .output()
        .map_err(|e| format!("running rustc: {e}"))?;
    let version = String::from_utf8_lossy(&output.stdout);
    let running = parse_rustc_version(&version)
        .ok_or_else(|| format!("could not parse `{}`", version.trim()))?;

    if running != pinned {
        return Err(format!(
            "toolchain mismatch: {TOOLCHAIN_FILE} pins {pinned}, but `rustc` here \
             is {running}.\n\
             The pin only applies when cargo comes from rustup. Check `command -v \
             cargo`;\n if it is not under ~/.cargo/bin, another toolchain is \
             shadowing it and\n every number this gate produces would come from an \
             unrecorded compiler."
        ));
    }
    println!("[ok] toolchain: rustc {running} matches {TOOLCHAIN_FILE}");
    Ok(())
}

/// The `channel` value from a `rust-toolchain.toml`, ignoring comments.
fn parse_pinned_channel(text: &str) -> Option<String> {
    text.lines()
        .map(str::trim)
        .filter(|line| !line.starts_with('#'))
        .find_map(|line| line.strip_prefix("channel"))
        .and_then(|rest| rest.trim().strip_prefix('='))
        .map(|rest| {
            // Strip a trailing inline comment before unquoting.
            let value = rest.split('#').next().unwrap_or(rest).trim();
            value.trim_matches('"').to_string()
        })
        .filter(|value| !value.is_empty())
}

/// The version from `rustc --version` output, e.g. `1.97.1`.
fn parse_rustc_version(output: &str) -> Option<&str> {
    let mut parts = output.split_whitespace();
    (parts.next()? == "rustc").then(|| parts.next())?
}

/// The development-only crates, which must never reach a shipped artifact.
///
/// `tepsim-oracle` links the original Fortran. `tepsim-stats` needs no gfortran
/// but exists to *judge* the model, and a shipped crate that depended on its
/// own judge would be a strange object. Both are cheap to check here and
/// expensive to discover at publish time.
const DEV_ONLY_CRATES: &[&str] = &["tepsim-oracle", "tepsim-stats"];

/// No shipped crate may name a development-only one.
fn check_oracle_isolation(root: &Path) -> Result<(), String> {
    for krate in SHIPPED_CRATES {
        let manifest = root.join("crates").join(krate).join("Cargo.toml");
        let text = fs::read_to_string(&manifest)
            .map_err(|e| format!("reading {}: {e}", manifest.display()))?;
        for dev in DEV_ONLY_CRATES {
            if text.contains(dev) {
                return Err(format!(
                    "{krate} references {dev}, which is development-only and \
                     must never be reachable from a shipped crate."
                ));
            }
        }
    }
    // The check has teeth only if the crates it names exist. A typo in
    // SHIPPED_CRATES would otherwise turn this into a loop over nothing.
    for dev in DEV_ONLY_CRATES {
        let manifest = root.join("crates").join(dev).join("Cargo.toml");
        if !manifest.exists() {
            return Err(format!(
                "{dev} is listed as development-only but {} does not exist",
                manifest.display()
            ));
        }
    }
    println!(
        "[ok] isolation: none of {} shipped crates depends on {:?}",
        SHIPPED_CRATES.len(),
        DEV_ONLY_CRATES
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// validate
// ---------------------------------------------------------------------------

/// The validation ladder at full volume.
///
/// Split out of `ci` deliberately. The Tier 1 sweep is ten million cases per
/// mode, about three minutes per routine in a debug build and growing with
/// every routine ported, which is the wrong thing to put on every commit. `ci`
/// runs the same assertions over the smoke sweep; this runs the gate, in
/// release, and the session protocol invokes it at every preflight.
fn cmd_validate(root: &Path, flags: &[String]) -> Result<(), String> {
    let tiers = parse_tiers(flags)?;
    if flags.iter().any(|f| f == "--compare-to-log") {
        println!(
            "[note] --compare-to-log is accepted but does nothing yet: it needs\n\
             the recorded numbers parsed out of LOG.org, which lands with the\n\
             validation report. Compare against the previous log entry by hand."
        );
    }

    for tier in &tiers {
        match tier {
            1 => {
                if which("gfortran").is_none() {
                    return Err(
                        "tier 1 needs gfortran, which is not on PATH. Per CLAUDE.md, \
                         a session without it must not do model work."
                            .to_string(),
                    );
                }
                println!("\n=== tier 1: utility routines vs the Fortran ===");
                // One target per invocation, so a failure names which routine
                // family broke rather than which binary did.
                for target in TIER1_TESTS {
                    step_with_env(
                        root,
                        "cargo",
                        &[
                            "test",
                            "-p",
                            "tepsim-oracle",
                            "--features",
                            "oracle",
                            "--release",
                            "--test",
                            target,
                            "--",
                            "--nocapture",
                            "--test-threads",
                            "1",
                        ],
                        &[("TEP_TIER1_SWEEP", "full")],
                    )?;
                }
            }
            2 => {
                if which("gfortran").is_none() {
                    return Err(
                        "tier 2 needs gfortran, which is not on PATH. Per CLAUDE.md, \
                         a session without it must not do model work."
                            .to_string(),
                    );
                }
                println!("\n=== tier 2: the plant model vs the Fortran ===");
                // One file per invocation, so a failure names which part of
                // the model broke rather than which binary did.
                for target in TIER2_TESTS {
                    step(
                        root,
                        "cargo",
                        &[
                            "test",
                            "-p",
                            "tepsim-oracle",
                            "--features",
                            "oracle",
                            "--release",
                            "--test",
                            target,
                            "--",
                            "--nocapture",
                            "--test-threads",
                            "1",
                        ],
                    )?;
                }
                // And again with the transcendentals taken out, where the
                // claim is bit equality rather than a tolerance. See
                // `tepsim_core::math`.
                println!("\n--- tier 2 again, on the platform libm ---");
                let mut exact = vec![
                    "test",
                    "-p",
                    "tepsim-oracle",
                    "--features",
                    "oracle,libm-system",
                    "--release",
                ];
                for name in LIBM_SYSTEM_TESTS {
                    exact.push("--test");
                    exact.push(name);
                }
                step(root, "cargo", &exact)?;
            }
            3 => {
                if which("gfortran").is_none() {
                    return Err(
                        "tier 3 needs gfortran, which is not on PATH. Per CLAUDE.md, \
                         a session without it must not do model work."
                            .to_string(),
                    );
                }
                println!("\n=== tier 3: the generator stream vs the Fortran ===");
                for target in TIER3_TESTS {
                    step(
                        root,
                        "cargo",
                        &[
                            "test",
                            "-p",
                            "tepsim-oracle",
                            "--features",
                            "oracle",
                            "--release",
                            "--test",
                            target,
                            "--",
                            "--nocapture",
                            "--test-threads",
                            "1",
                        ],
                    )?;
                }
            }
            4 => {
                if which("gfortran").is_none() {
                    return Err(
                        "tier 4 needs gfortran, which is not on PATH. Per CLAUDE.md, \
                         a session without it must not do model work."
                            .to_string(),
                    );
                }
                // Diagnostic, not a gate: `PLAN.org` is explicit that
                // long-horizon divergence is expected. Both configurations are
                // run, because the *contrast* is the result.
                println!("\n=== tier 4: trajectories (diagnostic, not a gate) ===");
                for features in ["oracle", "oracle,libm-system"] {
                    println!(
                        "\n--- with the {} libm ---",
                        if features.contains("system") {
                            "platform"
                        } else {
                            "vendored"
                        }
                    );
                    for target in TIER4_TESTS {
                        step_with_env(
                            root,
                            "cargo",
                            &[
                                "test",
                                "-p",
                                "tepsim-oracle",
                                "--features",
                                features,
                                "--release",
                                "--test",
                                target,
                                "--",
                                "--nocapture",
                                "--include-ignored",
                                "--test-threads",
                                "1",
                            ],
                            // `NPTS = 172800` at a one-second step: the run
                            // `temain_mod.f` was written to do. `ci` runs ten
                            // hours, which is already past the driver's forced
                            // `IDV(12)`; this runs the whole thing.
                            &[(TIER4_HOURS_ENV, "48")],
                        )?;
                    }
                }
            }
            5 => {
                if which("gfortran").is_none() {
                    return Err(
                        "tier 5 needs gfortran, which is not on PATH. Per CLAUDE.md, \
                         a session without it must not do model work."
                            .to_string(),
                    );
                }
                // 21 scenarios by 100 seeds by 48 h, on both sources: about
                // an hour of simulation. `ci` runs a smoke battery of twelve
                // runs over the same code.
                println!("\n=== tier 5: statistical equivalence ===");
                step_with_env(
                    root,
                    "cargo",
                    &[
                        "test",
                        "-p",
                        "tepsim-oracle",
                        "--features",
                        "oracle",
                        "--release",
                        "--test",
                        "tier5_battery",
                        "--",
                        "--nocapture",
                        "--test-threads",
                        "1",
                    ],
                    &[(TIER5_ENV, "full")],
                )?;
                // The invariants are Tier 5 too, and they are cheap.
                for features in ["oracle", "oracle,libm-system"] {
                    step(
                        root,
                        "cargo",
                        &[
                            "test",
                            "-p",
                            "tepsim-oracle",
                            "--features",
                            features,
                            "--release",
                            "--test",
                            "tier5_invariants",
                            "--",
                            "--nocapture",
                            "--test-threads",
                            "1",
                        ],
                    )?;
                }
            }
            other => println!(
                "\n[skip] tier {other}: not implemented yet. Tiers 6-10 land \
                 with their phases; see BACKLOG.org."
            ),
        }
    }

    println!("\nvalidate: green for tier(s) {tiers:?}");
    Ok(())
}

/// The integration tests that make up Tier 1, run at full sweep volume.
const TIER1_TESTS: [&str; 2] = ["tier1_enthalpy", "tier1_temperature"];

/// Selects Tier 4's horizon, as `TEP_TIER1_SWEEP` selects Tier 1's volume.
const TIER4_HOURS_ENV: &str = "TEP_TIER4_HOURS";

/// Selects Tier 5's battery size, likewise.
const TIER5_ENV: &str = "TEP_TIER5";

/// Tier 4. The open-loop trajectory is diagnostic; the closed-loop one is not
/// quite, because the controllers hold the plant at a setpoint and so remove
/// the amplification the open-loop run exists to characterise.
const TIER4_TESTS: [&str; 2] = ["tier4_trajectory", "tier4_closed_loop"];

/// `--tiers 1,2,3`, defaulting to every tier.
fn parse_tiers(flags: &[String]) -> Result<Vec<u8>, String> {
    let Some(position) = flags.iter().position(|f| f == "--tiers") else {
        return Ok((1..=10).collect());
    };
    let list = flags
        .get(position + 1)
        .ok_or_else(|| "--tiers needs a value, for example `--tiers 1,2`".to_string())?;
    list.split(',')
        .map(|t| {
            t.trim()
                .parse::<u8>()
                .map_err(|_| format!("`{t}` is not a tier number"))
                .and_then(|n| {
                    if (1..=10).contains(&n) {
                        Ok(n)
                    } else {
                        Err(format!("tier {n} does not exist; they run 1 to 10"))
                    }
                })
        })
        .collect()
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
    // Every vendored file, with its claims collected across the whole tree.
    let mut claims: Vec<Vec<LineRange>> = vec![Vec::new(); PROVENANCE_FILES.len()];
    let mut files_with_claims = 0usize;
    let mut unknown = Vec::new();

    for file in rust_sources(root)? {
        let text =
            fs::read_to_string(&file).map_err(|e| format!("reading {}: {e}", file.display()))?;
        let (found, bad) = parse_annotations(&text);
        if !found.is_empty() {
            files_with_claims += 1;
        }
        for claim in found {
            claims[claim.file].push(claim.range);
        }
        for entry in bad {
            unknown.push((file.clone(), entry));
        }
    }

    // A mistyped file name is a claim that stopped counting. Coverage would
    // fall with no explanation, so it is an error rather than a warning.
    if !unknown.is_empty() {
        let known: Vec<&str> = PROVENANCE_FILES.iter().map(|(tag, _)| *tag).collect();
        let lines: Vec<String> = unknown
            .iter()
            .map(|(path, entry)| format!("  {}: @port {}", path.display(), entry.tag))
            .collect();
        return Err(format!(
            "{} claim(s) name a file provenance does not know about:\n{}\n\
             Known prefixes: {known:?}. A mistyped name is a claim that stops \
             counting, which shows up as coverage falling for no reason.",
            unknown.len(),
            lines.join("\n")
        ));
    }

    println!(
        "provenance across {} vendored file(s)",
        PROVENANCE_FILES.len()
    );
    let mut any_missing = false;
    let (mut total_all, mut claimed_all) = (0usize, 0usize);

    for (index, (_, path)) in PROVENANCE_FILES.iter().enumerate() {
        let full = root.join(path);
        let Ok(source) = fs::read_to_string(&full) else {
            return Err(format!(
                "{path} not present, so coverage cannot be computed.\n\
                 Vendor the reference material first (B-0003). Failing rather \
                 than reporting a vacuous all-clear."
            ));
        };
        let total = source.lines().count();
        let merged = merge(core::mem::take(&mut claims[index]));
        let claimed_lines: usize = merged.iter().map(|r| r.len()).sum();
        let unclaimed = gaps(&merged, total);
        total_all += total;
        claimed_all += claimed_lines;

        println!("\n{path}");
        println!("  total lines:     {total}");
        println!(
            "  claimed:         {claimed_lines} ({:.1}%)",
            percent(claimed_lines, total)
        );
        println!("  unclaimed spans: {}", unclaimed.len());
        for range in &unclaimed {
            println!("    {range}  ({} lines)", range.len());
        }
        if unclaimed.is_empty() {
            println!("  every line of this file is claimed by some Rust function");
        } else {
            any_missing = true;
        }
    }

    println!(
        "\ntotal: {claimed_all} of {total_all} ({:.1}%) across {files_with_claims} Rust file(s)",
        percent(claimed_all, total_all)
    );
    if !any_missing {
        println!("every line of every vendored file is claimed");
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
fn parse_annotations(text: &str) -> (Vec<Claim>, Vec<UnknownClaim>) {
    let mut claims = Vec::new();
    let mut unknown = Vec::new();
    for line in text.lines() {
        match parse_claim_line(line) {
            Some(Ok(claim)) => claims.push(claim),
            Some(Err(bad)) => unknown.push(bad),
            None => {}
        }
    }
    (claims, unknown)
}

/// Parse one line as `// @port teprob.f:505-522`, in any comment style.
///
/// Returns `None` for prose that merely mentions the convention, for string
/// fixtures, and for a bare tag with no marker.
/// One `@port` claim: which file, and which lines of it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Claim {
    /// Index into [`PROVENANCE_FILES`].
    file: usize,
    range: LineRange,
}

/// A claim whose file name is not one of [`PROVENANCE_FILES`].
///
/// Reported rather than dropped. A mistyped file name is a claim that silently
/// stops counting, and the coverage number would go *down* with no explanation.
#[derive(Clone, Debug, PartialEq, Eq)]
struct UnknownClaim {
    tag: String,
}

fn parse_claim_line(line: &str) -> Option<Result<Claim, UnknownClaim>> {
    let body = strip_comment_prefix(line.trim_start())?;
    let after_marker = body.trim_start().strip_prefix(PROVENANCE_MARKER)?;
    // Require whitespace after the marker so `@ported` is not a claim.
    if !after_marker.starts_with(char::is_whitespace) {
        return None;
    }
    let rest = after_marker.trim_start();

    let Some((file, after_tag)) = PROVENANCE_FILES
        .iter()
        .enumerate()
        .find_map(|(index, (tag, _))| rest.strip_prefix(tag).map(|r| (index, r)))
    else {
        // It is a `@port` claim; it just names something unrecognised. Take
        // the first word so the report can quote it.
        let tag: String = rest.split_whitespace().next().unwrap_or(rest).to_string();
        return Some(Err(UnknownClaim { tag }));
    };

    let (start, tail) = take_usize(after_tag)?;
    let end = match tail.strip_prefix('-').and_then(take_usize) {
        Some((e, _)) => e,
        None => start,
    };
    Some(Ok(Claim {
        file,
        range: LineRange {
            start,
            end: end.max(start),
        },
    }))
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
    let path = root.join(GOLDEN_TRACE);
    let text = fs::read_to_string(&path).map_err(|e| {
        format!(
            "cannot read the golden trace at {GOLDEN_TRACE}: {e}\n\
             Regenerate it with:\n  \
             cargo run -p tepsim-oracle --features oracle --bin gen-golden-trace"
        )
    })?;

    let trace = Trace::parse(&text).map_err(|e| format!("{GOLDEN_TRACE}: {e}"))?;
    trace
        .require_full_length()
        .map_err(|e| format!("{GOLDEN_TRACE}: {e}"))?;

    println!("golden trace {GOLDEN_TRACE}");
    println!("  steps          : {}", trace.steps.len());
    println!("  values/step    : {}", golden::VALUES_PER_STEP);
    println!("  recorded with  : gfortran {}", trace.gfortran);
    println!("  fflags         : {}", trace.fflags);
    println!("  seed           : {}", trace.seed);
    println!("  dt             : {} h", trace.dt_hours);

    // The trace is only meaningful against the compiler that produced it. A
    // mismatch is a re-baseline, not a regression, and confusing the two costs
    // a whole session.
    match local_gfortran_version() {
        Some(local) if local != trace.gfortran => {
            println!(
                "\n  WARNING: local gfortran is {local}, the trace was recorded with \
                 {}.\n  Numbers measured against this trace may not reproduce here. \
                 Per CLAUDE.md\n  that is a logged re-baseline, not a regression hunt.",
                trace.gfortran
            );
        }
        Some(local) => println!("  local gfortran : {local} (matches)"),
        None => println!("  local gfortran : absent (fine; this check needs none)"),
    }

    // The other half of the preflight lives in tepsim-oracle, where the Fortran
    // can actually be re-run. This half runs the *port* against the trace, so
    // it needs no gfortran and takes about a second.
    diff_port_against(&trace)
}

/// Run the port forward and compare it against the recorded trace.
///
/// The port reproduces states, derivatives, measurements and the generator
/// word, so all four are compared. The tolerance is Tier 2's, applied the way
/// Tier 2 applies it: the derivatives are a cancelling quantity and are
/// measured against the scale of their own terms, everything else against its
/// own value. See the decision of 2026-08-27 in `BACKLOG.org`.
fn diff_port_against(trace: &Trace) -> Result<(), String> {
    use tepsim_core::{Inputs, Plant, SimTime, State, constants};

    /// `PLAN.org`, "Tier 2".
    const TOLERANCE: f64 = 1e-12;

    let mut plant = Plant::new();
    plant.set_rng(trace.seed);
    let mut state = State::from_flat(&constants::NOMINAL_STATE);
    // `TEINIT` sets `XMV(I) = YY(I+38)` at `teprob.f:1104`, so the manipulated
    // variables are the nominal valve positions. Taken from the state rather
    // than retyped: the first attempt here wrote them out by hand and
    // transposed two digits of `XMV(1)`, which this preflight caught at step
    // 79 as a 1.6e-5 divergence in `XMEAS(2)`.
    let inputs = Inputs {
        manipulated: core::array::from_fn(|i| constants::NOMINAL_STATE[38 + i]),
        disturbances: [0.0; 20],
    };

    // The first recorded state must be the one the port starts from, or the
    // whole comparison is against a different plant.
    for (slot, (ours, theirs)) in state
        .to_flat()
        .iter()
        .zip(trace.steps[0].states)
        .enumerate()
    {
        if ours.to_bits() != theirs.to_bits() {
            return Err(format!(
                "the port's nominal state disagrees with the trace at YY({}): \
                 {ours} against {theirs}.\n  \
                 constants::NOMINAL_STATE and the trace must come from the \
                 same TEINIT.",
                slot + 1
            ));
        }
    }

    let mut worst = (0.0_f64, String::new());
    let mut time = 0.0;
    for (index, step) in trace.steps.iter().enumerate() {
        let t = SimTime(time);
        plant.advance_discrete(t, &inputs);
        let (derivative, scale, signals) = plant
            .derivatives_with_scale(t, &state, &inputs)
            .map_err(|e| format!("step {index}: {e}"))?;
        let measurements = plant.sample_measurements(t, &signals);

        let mut check = |what: &str, slot: usize, ours: f64, theirs: f64, against: f64| {
            let error = if against == 0.0 {
                if ours.to_bits() == theirs.to_bits() {
                    0.0
                } else {
                    f64::INFINITY
                }
            } else {
                (ours - theirs).abs() / against.abs()
            };
            if error > worst.0 {
                worst = (error, format!("{what}({}) at step {index}", slot + 1));
            }
        };

        for (slot, (ours, theirs)) in state.to_flat().iter().zip(step.states).enumerate() {
            check("YY", slot, *ours, theirs, theirs);
        }
        for (slot, ((ours, theirs), budget)) in derivative
            .to_flat()
            .iter()
            .zip(step.derivatives)
            .zip(scale.to_flat())
            .enumerate()
        {
            check("YP", slot, *ours, theirs, budget);
        }
        for (slot, (ours, theirs)) in measurements
            .as_array()
            .iter()
            .zip(step.measurements)
            .enumerate()
        {
            check("XMEAS", slot, *ours, theirs, theirs);
        }
        check("G", 0, plant.rng(), step.rng, step.rng);

        plant
            .step_seeds(&state)
            .map_err(|e| format!("step {index}: {e}"))?;
        state = state.step(trace.dt_hours, &derivative);
        time += trace.dt_hours;
    }

    println!(
        "\n  port vs trace  : {} of {} steps diffed",
        trace.steps.len(),
        trace.steps.len()
    );
    println!("  worst          : {:e} at {}", worst.0, worst.1);
    println!("  gate           : {TOLERANCE:e}");
    if worst.0 > TOLERANCE {
        return Err(format!(
            "the port diverges from the golden trace: {:e} at {}, past the \
             {TOLERANCE:e} gate.\n  \
             This is the whole point of the preflight. Do not regenerate the \
             trace to make it pass.",
            worst.0, worst.1
        ));
    }
    println!("  fidelity: green");
    Ok(())
}

/// The local gfortran version, or `None` if there is no Fortran compiler here.
fn local_gfortran_version() -> Option<String> {
    let output = Command::new("gfortran")
        .arg("-dumpfullversion")
        .output()
        .ok()?;
    let version = String::from_utf8(output.stdout).ok()?.trim().to_string();
    (!version.is_empty()).then_some(version)
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

    /// The ranges a text claims against `teprob.f`, which is what the
    /// pre-B-0036 tests below were written against.
    fn teprob_ranges(text: &str) -> Vec<LineRange> {
        let (claims, unknown) = parse_annotations(text);
        assert!(unknown.is_empty(), "unexpected unknown claim: {unknown:?}");
        claims
            .into_iter()
            .filter(|c| c.file == 0)
            .map(|c| c.range)
            .collect()
    }

    #[test]
    fn parses_a_line_comment_claim() {
        assert_eq!(
            teprob_ranges("// @port teprob.f:505-522"),
            vec![r(505, 522)]
        );
    }

    #[test]
    fn parses_doc_and_module_doc_claims() {
        assert_eq!(teprob_ranges("/// @port teprob.f:1"), vec![r(1, 1)]);
        assert_eq!(
            teprob_ranges("//! @port teprob.f:1552-1560"),
            vec![r(1552, 1560)]
        );
    }

    #[test]
    fn indentation_and_trailing_text_are_allowed() {
        assert_eq!(
            teprob_ranges("        // @port teprob.f:100-200 (TEFUNC prologue)"),
            vec![r(100, 200)]
        );
    }

    #[test]
    fn parses_several_across_lines() {
        let text = "// @port teprob.f:1-10\ncode();\n/// @port teprob.f:20-30\n";
        assert_eq!(teprob_ranges(text), vec![r(1, 10), r(20, 30)]);
    }

    // The next four are the regression tests for the false-positive bug: the
    // first version of this scanner reported 2.6% coverage of a file nothing
    // had touched, by counting the tool's own documentation and fixtures.

    #[test]
    fn prose_mentioning_the_convention_is_not_a_claim() {
        let text = "//! Annotations look like `@port teprob.f:505-522`.";
        assert_eq!(teprob_ranges(text), vec![]);
    }

    #[test]
    fn a_string_fixture_is_not_a_claim() {
        let text = r#"assert_eq!(parse("@port teprob.f:1-10"), vec![]);"#;
        assert_eq!(teprob_ranges(text), vec![]);
    }

    #[test]
    fn a_bare_tag_without_the_marker_is_not_a_claim() {
        assert_eq!(teprob_ranges("// see teprob.f:505-522"), vec![]);
    }

    #[test]
    fn a_marker_glued_to_a_word_is_not_a_claim() {
        assert_eq!(teprob_ranges("// @ported teprob.f:1-10"), vec![]);
    }

    #[test]
    fn ignores_a_claim_with_no_number() {
        assert_eq!(teprob_ranges("// @port teprob.f:foo"), vec![]);
    }

    #[test]
    fn a_claim_names_its_file() {
        let teprob = parse_claim_line("// @port teprob.f:100-200").expect("a claim");
        let driver = parse_claim_line("// @port temain_mod.f:477-514").expect("a claim");
        assert_eq!(
            teprob,
            Ok(Claim {
                file: 0,
                range: r(100, 200)
            })
        );
        assert_eq!(
            driver,
            Ok(Claim {
                file: 1,
                range: r(477, 514)
            })
        );
        // A single line is a range of one, as before.
        assert_eq!(
            parse_claim_line("// @port teprob.f:42"),
            Some(Ok(Claim {
                file: 0,
                range: r(42, 42)
            }))
        );
    }

    /// A claim naming an unknown file is reported, not dropped.
    ///
    /// Silently ignoring it would show up as coverage falling with no
    /// explanation, which is the worst way for a typo to present.
    #[test]
    fn an_unknown_file_is_an_error_and_not_a_shrug() {
        let bad = parse_claim_line("// @port teprob.for:1-9").expect("recognised as a claim");
        assert!(bad.is_err(), "an unknown file should not parse as a claim");
        // And it quotes what it saw, so the typo is visible in the message.
        let Err(entry) = bad else { unreachable!() };
        assert!(entry.tag.starts_with("teprob.for"), "{}", entry.tag);
    }

    /// Prose about the convention must still not count as coverage.
    #[test]
    fn prose_mentioning_the_marker_is_still_not_a_claim() {
        assert!(parse_claim_line("// a line reading @ported teprob.f:1-9").is_none());
        assert!(parse_claim_line("//! `@port teprob.f:1-9` above a function").is_none());
        assert!(parse_claim_line("let x = 1; // @port").is_none());
    }

    #[test]
    fn a_reversed_range_is_clamped_not_panicked() {
        assert_eq!(teprob_ranges("// @port teprob.f:90-10"), vec![r(90, 90)]);
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

    /// The oracle links the original Fortran and is `publish = false`. If a
    /// shipped crate ever depended on it, wheels and the wasm bundle would try
    /// to drag a Fortran toolchain along. `xtask ci` checks this as a gate
    /// step; this makes plain `cargo test` catch it too.
    #[test]
    fn reads_the_pinned_channel_past_the_comment_block() {
        let text = "# a comment mentioning channel = \"9.9.9\"\n\
                    [toolchain]\n\
                    channel = \"1.97.1\"\n\
                    profile = \"minimal\"\n";
        assert_eq!(parse_pinned_channel(text).as_deref(), Some("1.97.1"));
    }

    #[test]
    fn tolerates_spacing_and_inline_comments() {
        assert_eq!(
            parse_pinned_channel("channel   =   \"1.90.0\"   # pinned").as_deref(),
            Some("1.90.0")
        );
    }

    #[test]
    fn a_file_with_no_channel_yields_nothing() {
        assert_eq!(
            parse_pinned_channel("[toolchain]\nprofile = \"minimal\"\n"),
            None
        );
    }

    #[test]
    fn reads_the_running_rustc_version() {
        assert_eq!(
            parse_rustc_version("rustc 1.97.1 (8bab26f4f 2026-07-14)\n"),
            Some("1.97.1")
        );
    }

    #[test]
    fn rejects_output_that_is_not_rustc() {
        assert_eq!(parse_rustc_version("cargo 1.97.1 (c980f4866)"), None);
        assert_eq!(parse_rustc_version(""), None);
        assert_eq!(parse_rustc_version("rustc"), None);
    }

    /// The failure this guard exists for: a pinned file plus a shadowing
    /// toolchain, which otherwise passes silently.
    #[test]
    fn a_mismatch_is_detectable_from_the_two_parsers_alone() {
        let pinned = parse_pinned_channel("channel = \"1.97.1\"").expect("pinned");
        let running =
            parse_rustc_version("rustc 1.89.0 (29483883e 2025-08-04) (Homebrew)").expect("running");
        assert_ne!(running, pinned);
    }

    #[test]
    fn no_shipped_crate_depends_on_the_oracle() {
        check_oracle_isolation(&workspace_root()).expect("oracle isolation");
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
