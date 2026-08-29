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
//! - `validate` runs the ladder at full volume and writes `book/src/validation/`
//!   from what the tests printed. The chapter used to be transcribed by hand
//!   out of `LOG.org`, which is a copy of a copy and rots invisibly.
//! - `deltas` cross-checks the `@delta` markers in the source against the
//!   entries in `book/src/deltas.md` and fails if either half is missing.
//! - `python` builds the wheel, installs it into a throwaway virtualenv and runs
//!   the pytest suite against it. Part of `ci`, skipped when maturin is absent.
//! - `licences` checks that the licence texts the wheel ships are the texts and
//!   not the symlink targets a checkout without symlink support leaves behind.
//!
//! See `CLAUDE.md` for how these fit into the session protocol.

#![forbid(unsafe_code)]

mod deltas;
mod plot;
mod report;
mod tier9;

use std::fmt;
use std::fmt::Write as _;
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
    "tepsim-operations",
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
        "fidelity" => cmd_fidelity(&root).map(|_| ()),
        "validate" => cmd_validate(&root, flags),
        "tier9" => tier9::cmd_tier9(&root, flags),
        "deltas" => deltas::cmd_deltas(&root),
        "python" => cmd_python(&root, true),
        "licences" => check_wheel_licences(&root),
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

  ci [--fast]   the gate: fmt, clippy, test, doc, deny, oracle isolation,
                the wheel's licence texts, and its pytest suite
                --fast skips the Fortran oracle differential job
  provenance    teprob.f line ranges not claimed by any Rust function
  fidelity      diff a short run against the committed golden oracle trace
  validate [--tiers 1,2,3] [--smoke] [--compare-to-log]
                run the validation ladder at full volume, in release, and
                write book/src/validation/ from what the tests printed
                --smoke runs the reduced sweeps, and says so on the page
  tier9 [--profiles release-wasm,release]
                cross-platform determinism: run the committed digest table
                natively, then again from a wasm32 build in a WebAssembly
                runtime, and compare both against the committed constants.
                Skipped, with the reason, without the target or a runtime
  deltas        cross-check the @delta markers against book/src/deltas.md,
                and write book/src/validation/delta-index.md
  python        build the wheel, install it into a throwaway virtualenv, and
                run the pytest suite against it. Skipped without maturin
  licences      the licence texts the wheel ships are texts, not symlink
                targets left by a checkout that could not follow them
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
    // Up here with the other structural checks rather than down in the python
    // job, because it costs nothing and a checkout that mangled the licence
    // symlinks is worth hearing about before five minutes of clippy, not after.
    check_wheel_licences(root)?;

    // The delta register's two halves have to agree: every `## D-0NN` heading
    // in `book/src/deltas.md` needs a `@delta` marker in the source, and every
    // marker needs a heading. Prose and code drift apart otherwise, and this
    // register is the project's record of every place the port knowingly
    // differs from the Fortran, which is exactly the document that must not
    // rot.
    //
    // When this check was written it failed on a real, pre-existing finding:
    // D-008 was documented and no marker named it, although the doc comment
    // above `const FAST` in `tepsim-control` said in prose that `CONTRL22` is
    // absent and that this is D-008. The marker was added and the check is a
    // gate step from that same commit.
    deltas::cmd_deltas(root)?;

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

    // Again with the oracle feature on. Without this pass, everything behind
    // `#[cfg(feature = "oracle")]` is never linted at all: the tier 2, 5, 6
    // and 7 machinery, the differentials, and the harness. That gap went
    // unnoticed until B-0049a reported ten errors nobody had ever seen, and
    // the fix was ten real lints in the library plus the observation that the
    // differential tests need `float_cmp` and `suboptimal_flops` allowed,
    // because exact comparison against the Fortran is what they are for.
    //
    // Skipped where the oracle is unavailable, like the rest of the oracle work.
    if oracle_supported() {
        step(
            root,
            "cargo",
            &[
                "clippy",
                "--workspace",
                "--all-targets",
                "--features",
                "tepsim-oracle/oracle",
                "--",
                "-D",
                "warnings",
            ],
        )?;
    } else {
        println!(
            "[skip] clippy with the oracle feature: {}",
            oracle_unavailable()
        );
    }
    // And a third pass under `libm-system`. Not redundant: that configuration
    // is `std`, so `f64::mul_add` exists and `clippy::suboptimal_flops` starts
    // firing on arithmetic that is silent under `no_std`. Nothing linted it
    // until B-0067, and a lint nothing runs is a lint that is not on.
    step(
        root,
        "cargo",
        &[
            "clippy",
            "--workspace",
            "--all-targets",
            "--features",
            "tepsim-core/libm-system",
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
    } else if !oracle_supported() {
        println!(
            "\n[skip] oracle differential job: {}.\n\
             Per CLAUDE.md, a session without the oracle must not do model work.",
            oracle_unavailable()
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

    // Last, because it builds a release wheel and the fast signals should come
    // first. The licences were already checked at the top of this function, so
    // it does not repeat that here.
    cmd_python(root, false)?;

    // TEP Studio's Node suite. Nothing ran it until B-0071, and the cost of
    // that was exactly what this project keeps rediscovering: its runner had
    // broken silently on a Node upgrade, and behind the broken runner sat a
    // test still asserting the pre-sign-off D-007 behaviour. Both were found
    // the moment something ran it.
    //
    // Node is optional, like gfortran and maturin: a checkout without it can
    // still gate the Rust. The suite itself skips its deployed-artifact check
    // when `apps/studio/dist` has not been built.
    let studio = root.join("apps/studio");
    if which("node").is_none() {
        println!("\n[skip] TEP Studio's Node suite: node is not on PATH");
    } else if !studio.join("dist").is_dir() {
        // The suite imports from `dist/` on purpose, not from `js/`:
        // `page.test.mjs` says so outright, because reading the built output is
        // how it proves `build.sh` copied what it should. That makes a build a
        // precondition rather than a convenience, and building it here would
        // mean a wasm toolchain on every `xtask ci`.
        //
        // `.github/workflows/pages.yml` builds the studio and then runs this
        // suite, so CI has the coverage unconditionally; locally it costs one
        // `apps/studio/build.sh`.
        println!(
            "\n[skip] TEP Studio's Node suite: apps/studio/dist does not exist.\n\
             Run `apps/studio/build.sh` once; the suite reads the built output \
             by design."
        );
    } else {
        step(&studio, "npm", &["test", "--silent"])?;
    }

    println!("\nci: green");
    Ok(())
}

// ---------------------------------------------------------------------------
// python
// ---------------------------------------------------------------------------

/// The maturin project: the wheel, the pytest suite, and the licence texts.
const PY_CRATE_DIR: &str = "crates/tepsim-py";

/// Where the throwaway wheel and virtualenv are built.
///
/// Under `target/`, which is already ignored, so a gate run leaves the working
/// tree clean.
/// Where the Python venv and the built wheel live.
///
/// Deliberately *not* under `target/`, though it was until B-0075b and that is
/// where you would expect scratch work to go.
///
/// `Swatinem/rust-cache` walks `target/` looking for build directories to
/// prune, and a Python venv underneath it is both a waste of that walk and a
/// hazard: numpy ships a fixture directory literally named
/// `numpy/testing/tests/target`, which the cache action tried to open and
/// failed on, annotating every CI run with an `ENOENT` that had nothing to do
/// with the build. Keeping the venv out of `target/` removes the whole class.
const PY_WORK_DIR: &str = ".xtask-python";

/// Build the wheel, install it into a throwaway virtualenv, run pytest.
///
/// The pytest suite tests the *binding*, and a binding only exists once it is
/// compiled and installed, so there is no way to run these tests against the
/// source tree. That is why they went unrun by the gate until B-0058a: a
/// `cargo test` cannot reach them.
///
/// Release, not debug. The suite asserts that `run()` releases the GIL by
/// timing four twelve-hour runs against each other, and a debug build turns
/// that from a second into a minute. Release is also the configuration a wheel
/// actually ships in, so this is the artifact under test rather than a proxy
/// for it.
///
/// Absent maturin this prints why and returns success. A machine that cannot
/// build a wheel cannot ship one either, and failing the Rust gate over it
/// would only teach people to stop running the gate. The same reasoning the
/// oracle job already uses for gfortran.
fn cmd_python(root: &Path, check_licences: bool) -> Result<(), String> {
    // Before the wheel is built, never after: a wheel whose licence file holds
    // `../../LICENSE` is a licence violation that has already been packaged.
    if check_licences {
        check_wheel_licences(root)?;
    }

    let Some(maturin) = which("maturin") else {
        println!(
            "\n[skip] python job: maturin is not on PATH.\n\
             The pytest suite in {PY_CRATE_DIR}/tests runs against an installed \
             wheel,\nso without maturin there is nothing to run it against. \
             `pip install maturin`."
        );
        return Ok(());
    };
    let Some(python) = which_python() else {
        println!(
            "\n[skip] python job: neither python3 nor python is on PATH.\n\
             The suite needs an interpreter to make a virtualenv with."
        );
        return Ok(());
    };

    // Throwaway means throwaway: a virtualenv left over from a previous run
    // could hold a stale `tepsim`, and pip would then have nothing to do and
    // report success. Remove first, ask questions never.
    let work = root.join(PY_WORK_DIR);
    if work.exists() {
        fs::remove_dir_all(&work).map_err(|e| format!("clearing {}: {e}", work.display()))?;
    }
    let dist = work.join("dist");
    let venv = work.join("venv");
    fs::create_dir_all(&dist).map_err(|e| format!("creating {}: {e}", dist.display()))?;

    let maturin = maturin.to_string_lossy().into_owned();
    let python = python.to_string_lossy().into_owned();
    let manifest = format!("{PY_CRATE_DIR}/Cargo.toml");
    let dist_arg = dist.to_string_lossy().into_owned();
    let venv_arg = venv.to_string_lossy().into_owned();

    step(
        root,
        &maturin,
        &[
            "build",
            "--release",
            "--out",
            &dist_arg,
            "--manifest-path",
            &manifest,
            // Explicit, so the wheel is built for the interpreter the tests are
            // about to run on. With `abi3-py39` the tag is `cp39-abi3` whatever
            // this is, but pinning it keeps the two halves from drifting.
            "--interpreter",
            &python,
        ],
    )?;

    let wheel = one_wheel(&dist)?;
    step(root, &python, &["-m", "venv", &venv_arg])?;

    let venv_python = venv_python(&venv).to_string_lossy().into_owned();
    // `[test]` pulls pytest from `[project.optional-dependencies]`, and numpy
    // comes from the wheel's own `dependencies`. Naming them here instead would
    // be a second place for the requirements to live.
    let requirement = format!("{}[test]", wheel.to_string_lossy());
    step(
        root,
        &venv_python,
        &[
            "-m",
            "pip",
            "install",
            "--disable-pip-version-check",
            &requirement,
        ],
    )?;

    // Prove the tests are about to import the wheel just built. Without this
    // the job passes just as happily against a `tepsim` that leaked in from
    // somewhere else, which is the one way a green Python gate could lie.
    // Written to a file rather than passed with `-c`, so the echoed command
    // stays one readable line.
    let probe = work.join("check_import_site.py");
    fs::write(&probe, IMPORT_SITE_CHECK)
        .map_err(|e| format!("writing {}: {e}", probe.display()))?;
    let probe_arg = probe.to_string_lossy().into_owned();
    step(root, &venv_python, &[&probe_arg, &venv_arg])?;

    step(
        root,
        &venv_python,
        &[
            "-m",
            "pytest",
            &format!("{PY_CRATE_DIR}/tests"),
            "-q",
            // No `.pytest_cache` in the working tree: the gate must leave the
            // tree clean enough to commit from.
            "-p",
            "no:cacheprovider",
        ],
    )?;

    println!("\n[ok] python: {} passed its pytest suite", wheel.display());
    Ok(())
}

/// Asserts that `import tepsim` resolves inside the throwaway virtualenv.
const IMPORT_SITE_CHECK: &str = r#"# Fail unless `tepsim` comes from the virtualenv named in argv[1].

import pathlib
import sys

import tepsim

here = pathlib.Path(tepsim.__file__).resolve()
venv = pathlib.Path(sys.argv[1]).resolve()
if venv not in here.parents:
    sys.exit(
        f"tepsim imported from {here}, which is outside {venv}. The tests would "
        "be checking some other installation, not the wheel just built."
    )
print(f"tepsim {tepsim.__version__} from {here}")
"#;

/// The single wheel maturin produced, or an explanation of why there is not one.
fn one_wheel(dist: &Path) -> Result<PathBuf, String> {
    let entries = fs::read_dir(dist).map_err(|e| format!("reading {}: {e}", dist.display()))?;
    let mut wheels = Vec::new();
    for entry in entries {
        let path = entry
            .map_err(|e| format!("reading {}: {e}", dist.display()))?
            .path();
        if path.extension().is_some_and(|e| e == "whl") {
            wheels.push(path);
        }
    }
    wheels.sort();
    match wheels.len() {
        1 => Ok(wheels.remove(0)),
        // Not "take the first": two wheels in a directory this command emptied
        // a moment ago means maturin built something unexpected, and installing
        // an arbitrary one of them would test an arbitrary artifact.
        n => Err(format!(
            "expected exactly one wheel in {}, found {n}: {wheels:?}",
            dist.display()
        )),
    }
}

/// The interpreter inside a virtualenv, which Windows puts somewhere else.
fn venv_python(venv: &Path) -> PathBuf {
    if cfg!(windows) {
        venv.join("Scripts").join("python.exe")
    } else {
        venv.join("bin").join("python")
    }
}

/// The Python to build the throwaway environment with.
///
/// `python3` first. A bare `python` is still Python 2 on a few systems, and on
/// Windows it may be the App Execution Alias stub that opens the Store.
fn which_python() -> Option<PathBuf> {
    which("python3").or_else(|| which("python"))
}

// ---------------------------------------------------------------------------
// licences
// ---------------------------------------------------------------------------

/// The licence texts the wheel ships must be texts, not paths.
///
/// `crates/tepsim-py/{LICENSE,LICENSE-NCSA,NOTICE.md}` are stored in git as
/// symlinks (mode 120000) to the ones at the repository root, so there is one
/// copy of each text and nothing can drift. That is fine everywhere git can
/// create a symlink.
///
/// Where it cannot -- a Windows checkout without `core.symlinks=true`, which is
/// the default when the installer could not enable them -- git writes the
/// *target path* as the file's content. `crates/tepsim-py/LICENSE` becomes a
/// 13-byte file reading `../../LICENSE`, `maturin build` copies it into the
/// wheel because `pyproject.toml` names it in `license-files`, and the wheel
/// ships with `../../LICENSE` where the NCSA licence should be. Nothing else in
/// the build would notice: the wheel is well-formed, the metadata is valid, and
/// the attribution the NCSA licence requires is simply gone.
///
/// So the check is: every text `license-files` declares must be the text, and
/// where the repository root holds the same name, the two must be byte for
/// byte the same file.
fn check_wheel_licences(root: &Path) -> Result<(), String> {
    let crate_dir = root.join(PY_CRATE_DIR);
    let manifest = crate_dir.join("pyproject.toml");
    let text = fs::read_to_string(&manifest)
        .map_err(|e| format!("reading {}: {e}", manifest.display()))?;
    // Read from pyproject rather than hard-coded here, so that a fourth licence
    // text added to the wheel cannot silently escape the check.
    let names = parse_license_files(&text).ok_or_else(|| {
        format!(
            "{} declares no `license-files`. The wheel is derived from the \
             original Tennessee Eastman Fortran, whose NCSA licence requires \
             attribution in binary distributions, so the texts have to ship \
             inside it.",
            manifest.display()
        )
    })?;

    verify_license_texts(&crate_dir, root, &names)?;
    println!(
        "[ok] licences: {} text(s) in {PY_CRATE_DIR} are texts and match the root",
        names.len()
    );
    Ok(())
}

/// The check itself, over explicit directories so a test can point it at one.
fn verify_license_texts(crate_dir: &Path, root: &Path, names: &[String]) -> Result<(), String> {
    for name in names {
        let path = crate_dir.join(name);
        let content = fs::read(&path).map_err(|e| {
            format!(
                "{} is declared in license-files but cannot be read: {e}",
                path.display()
            )
        })?;
        if content.is_empty() {
            return Err(format!("{} is empty", path.display()));
        }

        // The signature of a symlink written as text: one line, no newline, and
        // it names a file that is actually there. A licence is never that.
        let body = String::from_utf8_lossy(&content);
        let trimmed = body.trim();
        if !trimmed.contains('\n') && crate_dir.join(trimmed).exists() {
            return Err(format!(
                "{} is not a licence text. It holds the single line `{trimmed}`, \
                 which is the\n  *target* of a symlink rather than its contents.\n\n  \
                 git stores this file as a symlink (mode 120000) so there is one \
                 copy of the\n  text. A checkout that cannot create symlinks -- \
                 Windows without\n  `core.symlinks=true` -- writes the target path \
                 as the file body instead, and\n  a wheel built here would ship \
                 `{trimmed}` as its licence. The NCSA terms the\n  original Fortran \
                 comes under require that attribution in binary\n  distributions, so \
                 that is not cosmetic.\n\n  Fix the checkout, not the file. Replacing \
                 it with a copy of the text would\n  pass this check and reintroduce \
                 the drift the symlinks exist to prevent:\n    \
                 git config --global core.symlinks true\n    \
                 git checkout-index -f -a     # rewrite the working tree from the \
                 index\n  On Windows that needs Developer Mode or an elevated shell. \
                 A fresh clone\n  made after setting the config works too.",
                path.display(),
            ));
        }

        // One copy of each text is the whole point of the symlinks, so where
        // the root has the same name the bytes must agree. This is what catches
        // someone "fixing" a broken checkout by pasting in a stale copy.
        let counterpart = root.join(name);
        if counterpart.exists() {
            let expected = fs::read(&counterpart)
                .map_err(|e| format!("reading {}: {e}", counterpart.display()))?;
            if expected != content {
                return Err(format!(
                    "{} and {} differ ({} bytes against {}).\n  They are meant to \
                     be the same file: the one in the crate is a symlink to the \
                     one at\n  the root, so that there is a single copy of the \
                     text. Two copies drift.",
                    path.display(),
                    counterpart.display(),
                    content.len(),
                    expected.len(),
                ));
            }
        }
    }
    Ok(())
}

/// The `license-files` array from a `pyproject.toml`.
///
/// A deliberately small reader rather than a TOML dependency: xtask compiles at
/// the top of every session, and this is one flat array of strings.
fn parse_license_files(text: &str) -> Option<Vec<String>> {
    let mut lines = text
        .lines()
        .map(str::trim)
        .filter(|line| !line.starts_with('#'));

    let mut body = lines.find_map(|line| {
        let rest = line.strip_prefix("license-files")?;
        let rest = rest.trim_start().strip_prefix('=')?;
        Some(rest.trim_start().strip_prefix('[')?.to_string())
    })?;
    // Multi-line arrays are legal TOML. Truncating one silently would leave the
    // later entries unchecked, which is the failure this whole function exists
    // to prevent, so keep reading until the array closes.
    while !body.contains(']') {
        let line = lines.next()?;
        body.push(' ');
        body.push_str(line);
    }
    let body = body.split(']').next()?;

    let mut out = Vec::new();
    for piece in body.split(',') {
        let piece = piece.trim();
        if piece.is_empty() {
            continue;
        }
        let quote = piece.chars().next()?;
        if quote != '"' && quote != '\'' {
            return None;
        }
        let value = piece.strip_prefix(quote)?.strip_suffix(quote)?;
        if value.is_empty() {
            return None;
        }
        out.push(value.to_string());
    }
    (!out.is_empty()).then_some(out)
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
    // Reduced sweeps, and the pages say so. Without this the only way to run
    // the ladder at all is at full volume, which is minutes per Tier 1 target
    // and is why the report never got generated from a real run before now.
    let smoke = flags.iter().any(|f| f == "--smoke");
    if flags.iter().any(|f| f == "--compare-to-log") {
        println!(
            "[note] --compare-to-log is accepted but does nothing yet: it needs\n\
             the recorded numbers parsed out of LOG.org. The generated chapters\n\
             under book/src/validation/ are the diffable record in the meantime."
        );
    }

    // Which tiers wrote a chapter this run, for the index page. A tier that was
    // not selected, or that has no generator yet, must not be listed as if it
    // had been measured.
    let mut written: Vec<u8> = Vec::new();
    // The command actually being run, recorded on every page it writes. Not a
    // per-chapter command that would reproduce that one page: the header claims
    // "this wrote it", and that has to be literally what happened.
    let invocation = format!(
        "cargo xtask validate --tiers {}{}",
        tiers
            .iter()
            .map(u8::to_string)
            .collect::<Vec<_>>()
            .join(","),
        if smoke { " --smoke" } else { "" }
    );

    for tier in &tiers {
        match tier {
            1 => {
                require_gfortran(1)?;
                println!("\n=== tier 1: utility routines vs the Fortran ===");
                // One target per invocation, so a failure names which routine
                // family broke rather than which binary did.
                let env: &[(&str, &str)] = if smoke {
                    &[]
                } else {
                    &[("TEP_TIER1_SWEEP", "full")]
                };
                let mut runs = Vec::new();
                for target in TIER1_TESTS {
                    runs.push(report::run_target(
                        root,
                        target,
                        report::Libm::Vendored,
                        &[],
                        env,
                    )?);
                }
                write_chapter(
                    root,
                    1,
                    "the utility routines",
                    TIER1_LEAD,
                    smoke,
                    &invocation,
                    &runs,
                )?;
                written.push(1);
            }
            2 => {
                require_gfortran(2)?;
                println!("\n=== tier 2: the plant model vs the Fortran ===");
                // One file per invocation, so a failure names which part of
                // the model broke rather than which binary did.
                let mut runs = Vec::new();
                for target in TIER2_TESTS {
                    runs.push(report::run_target(
                        root,
                        target,
                        report::Libm::Vendored,
                        &[],
                        &[],
                    )?);
                }
                // And again with the transcendentals taken out, where the claim
                // is bit equality rather than a tolerance. See
                // `tepsim_core::math`. Per target rather than one invocation,
                // as above, and because a chapter needs one libtest tally per
                // row to report.
                println!("\n--- tier 2 again, on the platform libm ---");
                for target in LIBM_SYSTEM_TESTS {
                    runs.push(report::run_target(
                        root,
                        target,
                        report::Libm::Platform,
                        &[],
                        &[],
                    )?);
                }
                write_chapter(
                    root,
                    2,
                    "the plant model",
                    TIER2_LEAD,
                    smoke,
                    &invocation,
                    &runs,
                )?;
                written.push(2);
            }
            3 => {
                require_gfortran(3)?;
                println!("\n=== tier 3: the generator stream vs the Fortran ===");
                let mut runs = Vec::new();
                for target in TIER3_TESTS {
                    runs.push(report::run_target(
                        root,
                        target,
                        report::Libm::Vendored,
                        &[],
                        &[],
                    )?);
                }
                write_chapter(
                    root,
                    3,
                    "the generator stream",
                    TIER3_LEAD,
                    smoke,
                    &invocation,
                    &runs,
                )?;
                written.push(3);
            }
            4 => {
                require_gfortran(4)?;
                // Diagnostic, not a gate: `PLAN.org` is explicit that
                // long-horizon divergence is expected. Both configurations are
                // run, because the *contrast* is the result.
                println!("\n=== tier 4: trajectories (diagnostic, not a gate) ===");
                // Captured rather than merely echoed, because the contrast
                // between the two `libm` builds is the whole result and a
                // picture of it is worth more than two tables a page apart.
                let mut runs = Vec::new();
                for libm in [report::Libm::Vendored, report::Libm::Platform] {
                    println!("\n--- with the {} libm ---", libm.label());
                    for target in TIER4_TESTS {
                        runs.push(report::run_target(
                            root,
                            target,
                            libm,
                            &["--include-ignored"],
                            // `NPTS = 172800` at a one-second step: the run
                            // `temain_mod.f` was written to do. `ci` runs ten
                            // hours, which is already past the driver's forced
                            // `IDV(12)`; this runs the whole thing.
                            &[(TIER4_HOURS_ENV, "48")],
                        )?);
                    }
                }
                write_tier4_figure(root, &invocation, &runs)?;
            }
            5 => {
                require_gfortran(5)?;
                // 21 scenarios by 100 seeds by 48 h, on both sources: about
                // an hour of simulation. `ci` runs a smoke battery of twelve
                // runs over the same code.
                println!("\n=== tier 5: statistical equivalence ===");
                let battery = report::run_target(
                    root,
                    "tier5_battery",
                    report::Libm::Vendored,
                    &[],
                    if smoke { &[] } else { &[(TIER5_ENV, "full")] },
                )?;
                write_tier5_figure(root, &invocation, &[battery])?;
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
            6 => {
                if !oracle_supported() {
                    return Err(format!("tier 6 cannot run: {}.", oracle_unavailable()));
                }
                // The cross-source detector experiment, plus the detection
                // rates against the published files. Both take a size
                // selector; `validate` asks for the full one.
                println!("\n=== tier 6: downstream-task equivalence ===");
                for (target, size) in [
                    ("tier6_cross_source", TIER6_ENV),
                    ("tier6_published_rates", TIER6_RATES_ENV),
                ] {
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
                        &[(size, "full")],
                    )?;
                }
            }
            7 => {
                if !oracle_supported() {
                    return Err(format!("tier 7 cannot run: {}.", oracle_unavailable()));
                }
                // Every published file rather than the four-file smoke set.
                // About five minutes in release; thirty-five in debug, which
                // is why it is not in `ci`.
                println!("\n=== tier 7: published dataset reproduction ===");
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
                        "tier7_published",
                        "--",
                        "--nocapture",
                        "--test-threads",
                        "1",
                    ],
                    &[(TIER7_ENV, "full")],
                )?;
            }
            9 => {
                // The one tier that needs no Fortran at all. What it compares
                // is this build against constants committed to the repository,
                // so it is meaningful on a machine that could not build the
                // oracle, and it is the only tier a Windows runner could run.
                println!("\n=== tier 9: cross-platform determinism ===");
                tier9::cmd_tier9(root, &[])?;
            }
            10 => {
                if !oracle_supported() {
                    return Err(format!("tier 10 cannot run: {}.", oracle_unavailable()));
                }
                println!("\n=== tier 10: measured deltas for every Class C quirk ===");
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
                        "tier10_quirk_deltas",
                        "--",
                        "--nocapture",
                        "--test-threads",
                        "1",
                    ],
                )?;
            }
            other => println!("\n[skip] tier {other}: no harness in this arm; see BACKLOG.org."),
        }
    }

    // The book's tutorial transcripts. Not a tier, but the same shape of check
    // and the same reason it is here rather than in `ci`: re-running the three
    // examples costs about a minute in a debug build, and the cheap half, which
    // pins each listing to its example file, already runs on every commit.
    println!("\n--- book tutorials ---");
    step_with_env(
        root,
        "cargo",
        &["test", "-p", "tepsim", "--test", "book_examples"],
        &[("TEP_BOOK", "1")],
    )?;

    // The index last, so it describes a run that finished. It re-reads the
    // chapters on disk rather than trusting `written`, which is what keeps it
    // honest about a chapter some earlier run left behind.
    write_validation_index(root, &tiers, &written)?;

    println!("\nvalidate: green for tier(s) {tiers:?}");
    if written.is_empty() {
        println!(
            "no chapter was generated: only tiers {GENERATED_TIERS:?} have a \
             generator today."
        );
    }
    Ok(())
}

/// Whether this machine can build and run the Fortran oracle.
///
/// Two conditions, and the second is easy to forget. gfortran has to be on
/// PATH, and the platform has to be one the oracle is supported on, which per
/// `CLAUDE.md` is Linux and macOS and never Windows.
///
/// # Why the platform check is not redundant
///
/// It was, until CI ran on Windows for the first time in 180 commits. The
/// `windows-latest` runner image ships gfortran, so a bare `which("gfortran")`
/// says yes and the oracle build script runs. It then fails, because
/// `build.rs` hands gfortran a path from `fs::canonicalize`, which on Windows
/// is an extended-length `\\?\` path that MinGW's `f951.exe` cannot parse:
/// `Fatal Error: Cannot open file '\\teprob.f'`.
///
/// The fix is the policy, not the path. Every Tier 1 and Tier 2 number in
/// `LOG.org` was baselined against a specific gfortran on Linux or macOS, and
/// silently admitting a third toolchain would invalidate them without anyone
/// deciding to. See the note in `crates/tepsim-oracle/build.rs` for what a
/// deliberate Windows port would have to fix first.
fn oracle_supported() -> bool {
    !cfg!(windows) && which("gfortran").is_some()
}

/// Why the oracle is unavailable, for a skip message.
fn oracle_unavailable() -> &'static str {
    if cfg!(windows) {
        "the oracle is not supported on Windows; see `oracle_supported`"
    } else {
        "gfortran is not on PATH"
    }
}

/// No oracle, no model work. `CLAUDE.md` is explicit about this.
fn require_gfortran(tier: u8) -> Result<(), String> {
    if oracle_supported() {
        return Ok(());
    }
    Err(format!(
        "tier {tier} cannot run: {}. Per CLAUDE.md, a session without the \
         oracle must not do model work.",
        oracle_unavailable()
    ))
}

/// The tiers `validate` writes a chapter for.
///
/// Tiers 4 to 7, 9 and 10 still *run*; they just do not write a page yet.
/// Generating a chapter nobody has ever seen rendered, at the end of an
/// hour-long Tier 5 battery, is the wrong place to find out the renderer was
/// wrong, so they wait for a session that can afford to run them. Tier 8 has no
/// harness at all. The index says both of these on the page rather than leaving
/// a reader to infer it from an absence.
const GENERATED_TIERS: &[u8] = &[1, 2, 3];

/// The tiers `validate` runs but writes no chapter for.
const RUNNABLE_TIERS: &[u8] = &[4, 5, 6, 7, 9, 10];

/// Write one tier's chapter, and say which volume produced it.
fn write_chapter(
    root: &Path,
    tier: u8,
    title: &str,
    lead: &str,
    smoke: bool,
    command: &str,
    runs: &[report::TargetRun],
) -> Result<(), String> {
    let volume = if smoke {
        "**Reduced volume.** This run passed `--smoke`, so the sweeps are the \
         short ones the CI\ngate uses rather than the full ones `PLAN.org` \
         specifies. The case counts in the\ntables below are what actually ran. \
         Drop `--smoke` for the gate volume."
    } else {
        "Full volume: the sweeps `PLAN.org` specifies, not the reduced ones the \
         CI gate runs."
    };
    let lead = format!("{lead}\n\n{volume}");
    let figures = tier_figures(root, tier, command, runs)?;
    let page = report::render_tier(root, tier, title, &lead, command, runs, &figures);
    report::write_generated(root, &format!("{}/tier{tier}.md", report::DIR), &page)
}

/// The gate a tier's comparisons are drawn against.
///
/// `PLAN.org`: Tier 1 is a maximum relative error below 1e-13, Tier 2 below
/// 1e-12 of the scale of the terms. Tier 3 is a trace diff with no numeric
/// gate of its own, but the blocks it prints are `TESUB5` and `TESUB6`
/// comparisons, so it is drawn against Tier 1's line and the caption says so.
fn tier_gate(tier: u8) -> f64 {
    match tier {
        2 => 1e-12,
        _ => 1e-13,
    }
}

/// Which `LOG.org` iteration first measured what a tier's figures draw.
fn tier_provenance(tier: u8) -> &'static str {
    match tier {
        1 => "B-0009, B-0010 and B-0011",
        2 => "B-0026",
        3 => "B-0028 and B-0029",
        _ => "LOG.org",
    }
}

/// Draw a tier's figures, and return the markdown that inlines them.
///
/// Every caption states the failure condition, because a picture with no
/// stated way to be wrong is decoration. Nothing is drawn from a number that
/// was typed in: the two functions below read the same transcripts the tables
/// on the page are built from.
fn tier_figures(
    root: &Path,
    tier: u8,
    command: &str,
    runs: &[report::TargetRun],
) -> Result<String, String> {
    let gate = tier_gate(tier);
    let from = tier_provenance(tier);
    let mut out = String::new();

    let points = plot::error_points(runs);
    // Named from the data rather than described in prose. A tier reports some
    // comparisons it deliberately does not gate, and Tier 1 mis-types a
    // constant on purpose, so the orange dots are expected; which ones they
    // are is the thing worth pinning down, and a caption that lists them
    // changes the moment a different one appears.
    let outside = plot::outside_the_gate(&points, gate);
    let named = if outside.is_empty() {
        "Nothing is beyond it in this run.".to_string()
    } else {
        format!(
            "Beyond it in this run: {}, and nothing else. Any other dot \
             crossing the line is a failure.",
            outside
                .iter()
                .map(|label| format!("`{label}`"))
                .collect::<Vec<_>>()
                .join("; ")
        )
    };
    let strip_caption = format!(
        "Every comparison this tier made, at its own maximum relative error, \
         in a lane named for the target that ran it. Hovering a dot names the \
         comparison. A dot is orange because its value is past the {gate:e} \
         gate and for no other reason: no test is recognised by name, so a \
         regression and a deliberate positive control are drawn identically. \
         {named} The bit-equality half of the claim is falsified separately, \
         by any dot leaving the `= 0` lane under the platform `libm`, where \
         the two sides call the same `exp`."
    );
    let strip_id = format!("tier{tier}-errors");
    let ulp_id = format!("tier{tier}-ulp");

    match plot::strip(tier, gate, &points) {
        Some(svg) => {
            plot::write_figure(root, &svg, command)?;
            let caption = plot::Caption {
                id: &strip_id,
                title: "Every comparison, against the gate",
                caption: &strip_caption,
                from,
            };
            if let Some(markdown) = plot::include(root, "figures/", &caption) {
                out.push_str("\n## Figures\n\n");
                out.push_str(&markdown);
            }
        }
        // Said rather than left to be inferred from an absence. Every
        // comparison being exactly zero is the strongest result a tier can
        // have, and it is the one result a logarithmic error axis cannot draw.
        None if !points.is_empty() => out.push_str(
            "\n## Figures\n\nEvery one of this tier's comparisons is \
             bit-identical to the Fortran, so there is\nnothing to place on a \
             logarithmic error axis and no error figure is drawn. The \
             figure\nbelow states the same result in the units that suit it.\n",
        ),
        None => {}
    }
    if let Some(svg) = plot::ulp(tier, gate, runs) {
        plot::write_figure(root, &svg, command)?;
        // Only Tier 2 runs both builds, so only Tier 2's caption may name the
        // bit-equality condition. A caption that stated a failure condition
        // the figure cannot show would be exactly the decoration these are
        // supposed not to be.
        let both = runs.iter().any(|run| run.libm == report::Libm::Platform);
        let caption = plot::Caption {
            id: &ulp_id,
            title: "How many bits actually differ",
            caption: if both { ULP_CAPTION_BOTH } else { ULP_CAPTION },
            from,
        };
        if let Some(markdown) = plot::include(root, "figures/", &caption) {
            if out.is_empty() {
                out.push_str("\n## Figures\n\n");
            }
            out.push('\n');
            out.push_str(&markdown);
        }
    }
    Ok(out)
}

/// Tier 4's figure, from the scenario sweep both `libm` builds printed.
///
/// Tiers 4 and 5 write no chapter, for the reason `GENERATED_TIERS` gives, so
/// their figures land on the generated index instead. The file carries its own
/// provenance, so the index can say which run drew it even when that run was
/// not this one.
fn write_tier4_figure(
    root: &Path,
    command: &str,
    runs: &[report::TargetRun],
) -> Result<(), String> {
    let points = plot::noise_points(runs);
    let hours = plot::sweep_hours(runs).unwrap_or(0.0);
    match plot::noise_band(&points, hours) {
        Some(svg) => plot::write_figure(root, &svg, command),
        None => {
            println!(
                "[note] tier 4 printed no scenario sweep, so no figure was \
                 written. The sweep is `#[ignore]`d, and `validate` passes \
                 `--include-ignored` to run it."
            );
            Ok(())
        }
    }
}

/// Tier 5's figure: the cross-source value against the within-source null.
fn write_tier5_figure(
    root: &Path,
    command: &str,
    runs: &[report::TargetRun],
) -> Result<(), String> {
    let points = plot::calibration_points(runs);
    let size = plot::battery_size(runs).unwrap_or_else(|| "size not reported".to_string());
    match plot::calibration(&points, &size) {
        Some(svg) => plot::write_figure(root, &svg, command),
        None => {
            println!("[note] tier 5 printed no calibrated statistic, so no figure was written.");
            Ok(())
        }
    }
}

/// The two figures the generated index carries, with their failure conditions.
///
/// They are constants rather than built at the point of drawing because the
/// index is written on every run, including runs that did not touch Tier 4 or
/// Tier 5. What varies between runs, the command and the commit, is read back
/// off the figure itself.
const STANDALONE_FIGURES: &[plot::Caption<'static>] = &[
    plot::Caption {
        id: "tier4-noise-band",
        title: "How far apart the trajectories get, in units of instrument noise",
        caption: "\
Every disturbance scenario, run on both implementations, with the worst \
disagreement anywhere in the run divided by the noise standard deviation \
`XNS(i)` of the channel it happened on. The reference line is therefore not \
zero but the point at which the plant's own instruments could resolve the \
difference. The claim is false if any marker reaches the shaded band, and the \
explanation is false if the platform `libm` markers ever leave the `= 0` lane: \
that would mean the divergence is something other than transcendental \
rounding.",
        from: "B-0034",
    },
    plot::Caption {
        id: "tier5-calibration",
        title: "The two sources against the reference's own run-to-run spread",
        caption: "\
There is no absolute scale on which a Frobenius distance or a \
Kolmogorov-Smirnov statistic is small, so the battery measures one: it splits \
the Fortran's own runs in half and computes the same statistic Fortran against \
Fortran. The horizontal axis is the statistic across the two sources, the \
vertical axis is that null, and both axes share a scale so the diagonal is \
equality. The claim is false if a cross-source point crosses to the low side \
of the diagonal; the two that already sit there are the battery's own positive \
control, where one variable of the reference was shifted by ten standard \
deviations before comparing.",
        from: "B-0047b",
    },
];

/// The ULP figure says the same thing in the units that cannot be flattered.
const ULP_CAPTION: &str = "\
The same comparisons counted by how many bits differ rather than by how much. \
A relative error is a ratio, and dividing by a large number makes any \
difference look small; a count of differing units in the last place cannot be \
flattered that way. This tier is straight-line arithmetic over constants \
already proved bit-identical, so the claim is that the only bar is the one at \
zero. A second bar appearing anywhere is a regression to chase rather than a \
tolerance to widen.";

/// The same figure for a tier that ran both `libm` builds, where the contrast
/// between them is the result and the caption is allowed to say so.
const ULP_CAPTION_BOTH: &str = "\
The same comparisons counted by how many bits differ rather than by how much. \
A relative error is a ratio, and dividing by a large number makes any \
difference look small; a count of differing units in the last place cannot be \
flattered that way. The two groups are the whole argument: under the platform \
`libm`, where both sides call the same `exp`, the distribution is one bar at \
zero, and every comparison in the tier is identical to the last bit. Under the \
vendored `libm` a tail appears, and it is the transcendentals rather than the \
algebra. The claim is false the moment the platform group grows a second bar.";

const TIER1_LEAD: &str = "\
`TESUB1` (enthalpy), `TESUB2` (temperature from enthalpy by Newton), `TESUB3`
(heat capacity) and `TESUB4` (liquid density) are swept against the Fortran over
a simplex grid, a Dirichlet sample and a boundary pool, at every temperature in
the physical range, for each of the three `ITY` modes. `PLAN.org` sets the gate
at a maximum relative error below 1e-13, with a ULP histogram reported rather
than a verdict.";

const TIER2_LEAD: &str = "\
Both implementations are forced into an identical state, evaluated once, and
compared on all fifty derivative components. Sampling is from three pools: states
along the nominal closed-loop trajectory, random perturbations of those states,
and adversarial states placed at every discontinuity and clamp in the model.

Every comparison runs twice. The vendored `libm` disagrees with gfortran's by an
ULP on about a tenth of `exp` and `pow` calls, so the default build can only be
held to 1e-12; the `libm-system` build removes the transcendental from the
comparison and is held to bit equality. Both runs are below, and the `libm`
column says which is which.

The tolerance is relative to the scale of the terms rather than to the result. A
balance is inflow minus outflow, and near steady state those nearly cancel, so an
error that is 1e-16 of either term can be 1e-4 of their difference.";

/// The ladder, as `PLAN.org` defines it. No numbers here: the numbers are in
/// the chapters, and each of those came from a run.
const LADDER: &[(u8, &str)] = &[
    (1, "`TESUB1` to `TESUB8` match the oracle"),
    (2, "single-step derivatives match, over all three pools"),
    (3, "the generator call *order* matches, draw for draw"),
    (
        4,
        "trajectories stay inside the measurement noise (diagnostic)",
    ),
    (
        5,
        "statistical equivalence: TOST, KS, ACF, spectra, correlations",
    ),
    (6, "downstream detectors cannot tell the two sources apart"),
    (7, "the published `d00` to `d21` files are reproduced"),
    (8, "differential fuzzing finds no counterexample"),
    (9, "identical digests across platforms, wasm included"),
    (10, "every quirk fix ships with a measured delta"),
];

/// The index page: what has a generated chapter, and the preflight number.
fn write_validation_index(root: &Path, ran: &[u8], written: &[u8]) -> Result<(), String> {
    let command = "cargo xtask validate";
    let mut page = report::header(root, "Validation, measured", command, report::MEASURED);
    page.push_str(
        "\nThe narrative version of this material, with the reasoning behind \
         each tier and the\nhistory of what it caught, is in \
         [Validation](../validation.md). This section is the\nother half: the \
         numbers, written by the command that ran the suite, from the suite's\n\
         own output. Nothing here was transcribed.\n\n\
         A tier with no chapter has not been generated. That is stated rather \
         than left to be\ninferred from an absence, because a missing page and \
         a page nobody updated look the\nsame from the table of contents.\n\n",
    );
    let _ = writeln!(
        page,
        "Toolchain on the machine that ran this: {}.\n",
        report::toolchain()
    );

    // The preflight is cheap and needs no Fortran, so the index always carries
    // it. It is also the one number `CLAUDE.md` asks every session to look at.
    page.push_str("## Fidelity preflight\n\n");
    match cmd_fidelity(root) {
        Ok(f) => {
            let _ = writeln!(
                page,
                "`cargo xtask fidelity` runs the port forward from the nominal \
                 state and diffs states,\nderivatives, measurements and the \
                 generator word against a golden oracle trace\ncommitted to the \
                 repository. It needs no Fortran toolchain, so it runs \
                 everywhere in\nabout a second.\n\n\
                 | steps diffed | worst | where | gate | trace recorded with |\n\
                 |---|---|---|---|---|\n\
                 | {} of {} | {:e} | `{}` | {:e} | gfortran {} |\n",
                f.steps, f.steps, f.worst, f.worst_at, f.tolerance, f.recorded_with
            );
        }
        Err(e) => {
            // Reported on the page, not swallowed. An index that silently omits
            // a failing preflight is worse than one that has no preflight.
            let _ = writeln!(
                page,
                "**The preflight did not complete**, so this page carries no \
                 fidelity number:\n\n```text\n{e}\n```\n"
            );
        }
    }

    // Tiers 4 and 5 write no chapter, so their figures live here. Each is
    // included only if some run has actually drawn it, and the caption says
    // which run that was: a figure nobody has generated is reported as absent
    // rather than left out, because a missing picture and a picture nobody
    // updated look the same on a page.
    page.push_str(
        "## What the long tiers look like\n\nTiers 4 and 5 run but write no \
         chapter yet, for the reason `GENERATED_TIERS`\ngives in `xtask`. \
         Their figures are here, each drawn by the run named under it.\n",
    );
    for figure in STANDALONE_FIGURES {
        match plot::include(root, "figures/", figure) {
            Some(markdown) => {
                let _ = writeln!(page, "\n{markdown}");
            }
            None => {
                let _ = writeln!(
                    page,
                    "\n**{}.** No run has drawn this yet. `cargo xtask validate \
                     --tiers 4,5` writes it.",
                    figure.title
                );
            }
        }
    }
    page.push('\n');

    page.push_str("## The ladder\n\n| tier | what it proves | chapter |\n|---|---|---|\n");
    for (tier, proves) in LADDER {
        let relative = format!("{}/tier{tier}.md", report::DIR);
        let state = match report::read_provenance(root, &relative) {
            Some((by, commit)) => {
                format!("[tier {tier}](tier{tier}.md), from `{commit}` by `{by}`")
            }
            None if GENERATED_TIERS.contains(tier) => {
                format!("none yet: `cargo xtask validate --tiers {tier}`")
            }
            // Two different absences, and conflating them would misread the
            // second as the first. Tiers 4 and 5 have a suite and no generator;
            // 6 to 10 have neither.
            None if RUNNABLE_TIERS.contains(tier) => "runs, but writes no chapter yet".to_string(),
            None => "no harness yet".to_string(),
        };
        let _ = writeln!(page, "| {tier} | {proves} | {state} |");
    }

    let _ = writeln!(
        page,
        "\nThis run selected tier(s) {ran:?} and wrote chapter(s) for {written:?}. \
         Tiers 4 to 7, 9\nand 10 run but do not write a chapter yet; tier 8 has no \
         harness. The [delta\nindex](delta-index.md) is generated separately, by \
         `cargo xtask deltas`.\n"
    );
    report::write_generated(root, &format!("{}/index.md", report::DIR), &page)
}

const TIER3_LEAD: &str = "\
Both sides are instrumented to emit every generator draw, and the traces are
diffed. This is the tier that catches a port whose arithmetic is right and whose
*call order* is not, which no statistical comparison would find until after a
48-hour run.";

/// The integration tests that make up Tier 1, run at full sweep volume.
const TIER1_TESTS: [&str; 2] = ["tier1_enthalpy", "tier1_temperature"];

/// Selects Tier 4's horizon, as `TEP_TIER1_SWEEP` selects Tier 1's volume.
const TIER4_HOURS_ENV: &str = "TEP_TIER4_HOURS";

/// Selects Tier 5's battery size, likewise.
const TIER5_ENV: &str = "TEP_TIER5";

/// Selects Tier 6's cross-source battery size.
const TIER6_ENV: &str = "TEP_TIER6";

/// Selects how much of the published-rate sweep runs.
const TIER6_RATES_ENV: &str = "TEP_TIER6_RATES";

/// Selects whether Tier 7 covers every published file or the smoke set.
const TIER7_ENV: &str = "TEP_TIER7";

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
    // On Windows an executable on PATH is `maturin.exe` or `maturin.bat`, never
    // a bare `maturin`. Probing only the plain name reports every tool as
    // absent, which for a job that *skips* when a tool is missing means it
    // skips forever and says so convincingly.
    let names: Vec<String> = if cfg!(windows) {
        [".exe", ".cmd", ".bat", ""]
            .iter()
            .map(|extension| format!("{program}{extension}"))
            .collect()
    } else {
        vec![program.to_string()]
    };
    std::env::split_paths(&path)
        .flat_map(|dir| names.iter().map(move |name| dir.join(name)))
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

/// What the preflight measured, so the report can quote it.
///
/// Returned rather than only printed: the validation index carries this number,
/// and re-deriving it by parsing this command's own stdout would be a second
/// implementation of the same measurement.
pub(crate) struct Fidelity {
    steps: usize,
    worst: f64,
    worst_at: String,
    tolerance: f64,
    /// The gfortran that recorded the trace, which is not necessarily the one
    /// on this machine. A mismatch is a re-baseline, not a regression.
    recorded_with: String,
}

fn cmd_fidelity(root: &Path) -> Result<Fidelity, String> {
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
fn diff_port_against(trace: &Trace) -> Result<Fidelity, String> {
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
    Ok(Fidelity {
        steps: trace.steps.len(),
        worst: worst.0,
        worst_at: worst.1,
        tolerance: TOLERANCE,
        recorded_with: trace.gfortran.clone(),
    })
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

    // -----------------------------------------------------------------------
    // licences
    // -----------------------------------------------------------------------

    #[test]
    fn reads_the_license_files_array() {
        let text = "[project]\nname = \"tepsim\"\n\
                    license-files = [\"LICENSE\", \"LICENSE-NCSA\", \"NOTICE.md\"]\n";
        assert_eq!(
            parse_license_files(text),
            Some(vec![
                "LICENSE".to_string(),
                "LICENSE-NCSA".to_string(),
                "NOTICE.md".to_string(),
            ])
        );
    }

    /// A multi-line array must not be read as its first line only: the entries
    /// past the truncation would go unchecked, which is the exact failure the
    /// guard exists to prevent.
    #[test]
    fn reads_an_array_spread_over_lines() {
        let text = "license-files = [\n  \"LICENSE\",\n  \"NOTICE.md\",\n]\n";
        assert_eq!(
            parse_license_files(text),
            Some(vec!["LICENSE".to_string(), "NOTICE.md".to_string()])
        );
        // Unterminated is not "everything up to here", it is unreadable.
        assert_eq!(
            parse_license_files("license-files = [\n  \"LICENSE\",\n"),
            None
        );
    }

    #[test]
    fn a_pyproject_with_no_license_files_yields_nothing() {
        assert_eq!(parse_license_files("[project]\nname = \"tepsim\"\n"), None);
        assert_eq!(parse_license_files("license-files = []\n"), None);
        // A bare word is not a string, and guessing at what was meant would be
        // worse than saying so.
        assert_eq!(parse_license_files("license-files = [LICENSE]\n"), None);
    }

    /// A scratch tree shaped like the repository: a root, and a crate under it.
    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("xtask-licence-{}-{name}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(dir.join(PY_CRATE_DIR)).expect("scratch tree");
        fs::write(dir.join("LICENSE"), "MIT License\n\nCopyright ...\n").expect("root licence");
        dir
    }

    fn names() -> Vec<String> {
        vec!["LICENSE".to_string()]
    }

    #[test]
    fn a_real_licence_text_passes() {
        let dir = scratch("good");
        let crate_dir = dir.join(PY_CRATE_DIR);
        fs::copy(dir.join("LICENSE"), crate_dir.join("LICENSE")).expect("copy");

        verify_license_texts(&crate_dir, &dir, &names()).expect("the texts agree");
        let _ = fs::remove_dir_all(&dir);
    }

    /// The hazard itself, reproduced without a Windows machine: git writing the
    /// symlink target as the file body.
    #[test]
    fn a_symlink_written_as_text_is_caught() {
        let dir = scratch("symlink-as-text");
        let crate_dir = dir.join(PY_CRATE_DIR);
        fs::write(crate_dir.join("LICENSE"), "../../LICENSE").expect("write");

        let error = verify_license_texts(&crate_dir, &dir, &names())
            .expect_err("a path is not a licence text");
        assert!(error.contains("../../LICENSE"), "{error}");
        // The message has to name the fix, because the file looks fine in an
        // editor and the wheel that ships it is well-formed.
        assert!(error.contains("core.symlinks"), "{error}");
        let _ = fs::remove_dir_all(&dir);
    }

    /// Only when it really is a path. A licence whose first line happens to be
    /// short must not trip the check.
    #[test]
    fn a_short_first_line_is_not_a_symlink() {
        let dir = scratch("short-line");
        let crate_dir = dir.join(PY_CRATE_DIR);
        let text = "LICENSE\n\nThe text of a licence that opens with its own name.\n";
        fs::write(dir.join("LICENSE"), text).expect("root");
        fs::write(crate_dir.join("LICENSE"), text).expect("crate");

        verify_license_texts(&crate_dir, &dir, &names()).expect("multi-line text is a text");
        let _ = fs::remove_dir_all(&dir);
    }

    /// Two copies drift. Replacing the symlink with a stale paste is the
    /// obvious wrong fix for the failure above, so it is caught too.
    #[test]
    fn a_second_copy_that_drifted_is_caught() {
        let dir = scratch("drift");
        let crate_dir = dir.join(PY_CRATE_DIR);
        fs::write(crate_dir.join("LICENSE"), "MIT License\n\nCopyright 1993\n").expect("write");

        let error =
            verify_license_texts(&crate_dir, &dir, &names()).expect_err("the copies disagree");
        assert!(error.contains("differ"), "{error}");
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_missing_or_empty_licence_is_caught() {
        let dir = scratch("missing");
        let crate_dir = dir.join(PY_CRATE_DIR);

        let error = verify_license_texts(&crate_dir, &dir, &names()).expect_err("nothing to ship");
        assert!(error.contains("cannot be read"), "{error}");

        fs::write(crate_dir.join("LICENSE"), "").expect("write");
        let error = verify_license_texts(&crate_dir, &dir, &names()).expect_err("empty");
        assert!(error.contains("is empty"), "{error}");
        let _ = fs::remove_dir_all(&dir);
    }

    /// And the real checkout this session is running in.
    #[test]
    fn this_checkout_ships_the_licence_texts() {
        check_wheel_licences(&workspace_root()).expect("licence texts");
    }

    #[test]
    fn a_virtualenv_puts_python_where_the_platform_does() {
        let venv = Path::new("/tmp/venv");
        let python = venv_python(venv);
        assert!(python.starts_with(venv));
        let tail = if cfg!(windows) {
            "Scripts/python.exe"
        } else {
            "bin/python"
        };
        assert!(python.ends_with(tail), "{}", python.display());
    }
}
