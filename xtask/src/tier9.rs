//! `cargo xtask tier9`: the same table, on this architecture and on wasm32.
//!
//! # What this proves and what it does not
//!
//! [`tepsim::tier9::CASES`] pairs six fixed scenarios with the digest each is
//! committed to produce. The constants are in the source tree, so a platform
//! checks itself: this command computes them natively, then builds
//! `tepsim-wasm` for `wasm32-unknown-unknown` and has a WebAssembly runtime
//! compute them again inside the module, and compares both against the
//! committed values.
//!
//! What that establishes is "this architecture and wasm32 agree with what was
//! committed". It says nothing about an architecture nobody has run it on.
//! Completing the claim for x86-64 means running this command on an x86-64
//! machine, and it needs no edit here to do that: the constants are committed,
//! not computed twice in one process, so the second machine is comparing
//! against the first machine's measurement rather than against itself.
//!
//! # Why a runtime and not just `cargo test --target wasm32`
//!
//! `cargo test` for `wasm32-unknown-unknown` needs a test runner configured for
//! the target, which means `wasm-bindgen-test` and a browser or a Node harness
//! driven by `wasm-bindgen-cli` at a version pinned to the crate. That is a
//! reasonable thing to ask of someone building the browser app and too much to
//! put between a session and the check that says the numbers are the same
//! everywhere. `tepsim-wasm`'s `tier9` module exports the table through plain
//! `extern "C"` functions taking and returning integers, so any runtime that
//! can instantiate a module reads them back with no code generator involved.
//!
//! # The stubbed imports
//!
//! The rest of `tepsim-wasm` is `#[wasm_bindgen]`, so the compiled module still
//! *declares* imports from `__wbindgen_placeholder__`. Nothing on the Tier 9
//! path calls them. The Node driver supplies a throwing function for each,
//! which is deliberate: if the digest path ever did reach into JavaScript, that
//! is a determinism bug, and it should fail loudly rather than be quietly
//! satisfied by a stub that returns zero.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

use tepsim::tier9::{self, CASES};

use crate::{step, which};

/// Profiles the wasm module is built and checked in, unless `--profiles` says
/// otherwise.
///
/// Two rather than one because optimisation level is exactly the kind of thing
/// that could reassociate an expression, and a determinism claim that only
/// holds at `opt-level = "s"` is not the claim. `release-wasm` is what the
/// browser app ships (`apps/studio/build.sh`); `release` is the plain one.
/// `dev` is available through the flag and is left out of the default because
/// an unoptimised wasm build runs the eighteen simulated hours in this table
/// slowly enough to be noticed.
const DEFAULT_PROFILES: &[&str] = &["release-wasm", "release"];

/// The target the wasm half is built for.
const WASM_TARGET: &str = "wasm32-unknown-unknown";

/// Where the generated Node driver is written.
///
/// Under `target/` next to `xtask-python`, so it is a build artifact and not
/// something a reader has to wonder about in the source tree.
const WORK_DIR: &str = "target/xtask-tier9";

pub(crate) fn cmd_tier9(root: &Path, flags: &[String]) -> Result<(), String> {
    let profiles = parse_profiles(flags)?;

    println!("=== tier 9: cross-platform determinism ===\n");
    println!("host      : {}", host_triple());
    println!("toolchain : {}", crate::report::toolchain());
    println!(
        "table     : {} cases, committed in crates/tepsim/src/tier9.rs",
        CASES.len()
    );

    // Native first. If this machine cannot reproduce its own committed
    // constants, nothing the wasm half says is interesting.
    let native = native_pass()?;

    check_browser_copies(root)?;

    let wasm = wasm_pass(root, &profiles)?;

    summary(&native, &wasm, &profiles);

    if !native.agrees() {
        return Err(mismatch_message("native", &native));
    }
    for run in &wasm {
        if !run.digests.agrees() {
            return Err(mismatch_message(&run.label, &run.digests));
        }
    }
    Ok(())
}

/// One target's answers for the whole table.
#[derive(Debug)]
struct Digests {
    /// Per case, in table order.
    cases: Vec<u64>,
    /// The suite digest as that target computed it.
    suite: u64,
}

impl Digests {
    fn agrees(&self) -> bool {
        self.cases.len() == CASES.len()
            && self
                .cases
                .iter()
                .zip(CASES)
                .all(|(computed, case)| *computed == case.digest)
            && self.suite == tier9::SUITE_DIGEST
    }
}

/// A wasm build that was actually evaluated.
struct WasmRun {
    label: String,
    runtime: String,
    bytes: u64,
    digests: Digests,
}

// ---------------------------------------------------------------------------
// native
// ---------------------------------------------------------------------------

fn native_pass() -> Result<Digests, String> {
    println!("\n--- native: {} ---\n", host_triple());
    let digests = Digests {
        cases: CASES.iter().map(tier9::Case::compute).collect(),
        suite: tier9::suite_digest(),
    };
    print_table(&digests);
    Ok(digests)
}

// ---------------------------------------------------------------------------
// the browser app's copies
// ---------------------------------------------------------------------------

/// JavaScript and HTML that hard-code case 0's digest.
///
/// Each of these compares a browser's `selfCheckDigest()` against a string
/// typed into the file. Nothing in the Rust build sees those strings, so a
/// legitimate re-baseline of the table would leave them behind and the browser
/// would report a determinism failure that is really a stale constant, or a
/// change made here and not there would go unnoticed until somebody opened the
/// page. Checking them is three file reads.
const BROWSER_COPIES: &[&str] = &[
    "apps/studio/js/app.js",
    "apps/studio/node/protocol.test.mjs",
    "crates/tepsim-wasm/www/index.html",
];

/// The identifier the three files assign the digest to.
const BROWSER_CONSTANT: &str = "EXPECTED_SELF_CHECK";

fn check_browser_copies(root: &Path) -> Result<(), String> {
    println!("\n--- the browser app's copy of case 0 ---\n");
    let expected = format!("{:016x}", CASES[0].digest);

    for relative in BROWSER_COPIES {
        let path = root.join(relative);
        let Ok(text) = std::fs::read_to_string(&path) else {
            // Absent is a note, not a failure. A checkout without the browser
            // app is a legitimate thing to run tier 9 in.
            println!("  {relative:<38} not present, skipped");
            continue;
        };
        match find_browser_digest(&text) {
            Some(found) if found == expected => println!("  {relative:<38} {found} match"),
            Some(found) => {
                return Err(format!(
                    "{relative} has {BROWSER_CONSTANT} = \"{found}\", and tier 9 \
                     case `{}` is committed to {expected}.\n\nThe browser \
                     compares its own digest against that string, so as it \
                     stands the page would report a determinism failure that is \
                     really a stale constant. Fix whichever is wrong; do not \
                     assume it is the Rust one.",
                    CASES[0].name
                ));
            }
            None => {
                return Err(format!(
                    "{relative} no longer assigns {BROWSER_CONSTANT} a \
                     sixteen-digit hex string. Either the browser stopped \
                     checking its digest, which is a regression, or it was \
                     renamed and BROWSER_COPIES in xtask/src/tier9.rs needs to \
                     follow it."
                ));
            }
        }
    }
    Ok(())
}

/// Pull the digest out of `const EXPECTED_SELF_CHECK = "c8a2...";`.
///
/// A search for the identifier and then for the next quoted run of sixteen hex
/// digits, rather than a regex or a JavaScript parse. The three files write the
/// assignment three slightly different ways and all that matters is the string
/// next to the name.
fn find_browser_digest(text: &str) -> Option<String> {
    for (offset, _) in text.match_indices(BROWSER_CONSTANT) {
        let after = &text[offset + BROWSER_CONSTANT.len()..];
        // Only an assignment, not a use. `=== EXPECTED_SELF_CHECK` and
        // `${EXPECTED_SELF_CHECK}` are comparisons and interpolations, and
        // taking the next quoted hex run after one of those would find some
        // other string entirely.
        let Some(rest) = after.trim_start().strip_prefix('=') else {
            continue;
        };
        // `continue`, not `?`: `EXPECTED_SELF_CHECK == "..."` also survives the
        // `=` above, and giving up on the whole file at the first name that is
        // not an assignment would miss the assignment further down.
        let Some(quoted) = rest.trim_start().strip_prefix(['"', '\'']) else {
            continue;
        };
        let digest: String = quoted.chars().take(16).collect();
        if digest.len() == 16 && digest.chars().all(|c| c.is_ascii_hexdigit()) {
            return Some(digest);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// wasm32
// ---------------------------------------------------------------------------

fn wasm_pass(root: &Path, profiles: &[String]) -> Result<Vec<WasmRun>, String> {
    println!("\n--- wasm32-unknown-unknown ---");

    if !wasm_target_installed() {
        println!(
            "\n[skip] the {WASM_TARGET} target is not installed, so the wasm \
             half of tier 9\n       did not run. Install it with:\n\n         \
             rustup target add {WASM_TARGET}\n\n       The native half above \
             still stands on its own."
        );
        return Ok(Vec::new());
    }

    let Some(runtime) = pick_runtime() else {
        println!(
            "\n[skip] no WebAssembly runtime was found on PATH, so the wasm \
             half of tier 9\n       did not run. Any one of these is enough:\n\n\
             \x20        node       (what this command drives; V8, the engine \
             the browser app runs on)\n\
             \x20        wasmtime   (see the note below)\n\
             \x20        wasmer     (see the note below)\n\n\
             {RUNTIME_NOTE}\n\n       The native half above still stands on \
             its own."
        );
        return Ok(Vec::new());
    };

    let Runtime::Node(node) = &runtime;
    println!(
        "\nruntime   : node {} at {}",
        node_version(node).unwrap_or_else(|| "unknown".into()),
        node.display()
    );

    let driver = write_driver(root)?;

    let mut runs = Vec::new();
    for profile in profiles {
        println!("\n--- wasm32, profile `{profile}` ---");
        let module = build_wasm(root, profile)?;
        let bytes = std::fs::metadata(&module).map(|m| m.len()).unwrap_or(0);
        let digests = run_driver(node, &driver, &module)?;
        print_table(&digests);
        runs.push(WasmRun {
            label: format!("wasm32 `{profile}`"),
            runtime: format!("node {}", node_version(node).unwrap_or_default()),
            bytes,
            digests,
        });
    }
    Ok(runs)
}

/// The runtimes this command knows how to drive.
///
/// One variant today. It is an enum rather than a `PathBuf` because the shape
/// of the problem is "which of several runtimes did we find", and flattening
/// that to a path would lose the reason a `wasmtime` on PATH is not usable
/// here.
enum Runtime {
    Node(PathBuf),
}

/// Why `wasmtime` and `wasmer` are found but not used.
///
/// Written once and printed wherever it is relevant, because a reader who has
/// `wasmtime` installed and sees it skipped deserves the actual reason rather
/// than a shrug.
const RUNTIME_NOTE: &str = "\
       Note on wasmtime and wasmer: the module carries import declarations from\n\
       `__wbindgen_placeholder__`, because the rest of tepsim-wasm is\n\
       #[wasm_bindgen]. Nothing on the tier 9 path calls them, but an embedder\n\
       must still supply one function per import before it can instantiate, and\n\
       neither CLI can be told to do that from the command line. Driving them\n\
       would mean linking their embedding API into xtask, which is a large\n\
       dependency for a job node already does. If you have one of them and no\n\
       node, that is the change to make.";

fn pick_runtime() -> Option<Runtime> {
    // node first, and not only because it is the one implemented: it is V8,
    // which is the engine the browser app actually runs on, so agreement here
    // is the closest thing to the claim `PLAN.org` makes about a real browser.
    if let Some(node) = which("node") {
        return Some(Runtime::Node(node));
    }
    // Found, reported, and not usable. See RUNTIME_NOTE.
    for other in ["wasmtime", "wasmer"] {
        if let Some(path) = which(other) {
            println!(
                "\n[note] found {other} at {}, but this command cannot drive \
                 it.\n\n{RUNTIME_NOTE}",
                path.display()
            );
        }
    }
    None
}

fn node_version(node: &Path) -> Option<String> {
    let out = Command::new(node).arg("--version").output().ok()?;
    let text = String::from_utf8(out.stdout).ok()?;
    Some(text.trim().to_string()).filter(|s| !s.is_empty())
}

/// Whether the standard library for `wasm32-unknown-unknown` is present.
///
/// Asked of the sysroot rather than of `rustup target list`, so the answer is
/// right on a toolchain that rustup did not install. `cargo` would otherwise
/// fail several seconds into a build with a message about a missing crate,
/// which reads like a code error rather than a missing target.
fn wasm_target_installed() -> bool {
    let Ok(out) = Command::new("rustc").arg("--print").arg("sysroot").output() else {
        return false;
    };
    let Ok(sysroot) = String::from_utf8(out.stdout) else {
        return false;
    };
    Path::new(sysroot.trim())
        .join("lib/rustlib")
        .join(WASM_TARGET)
        .join("lib")
        .is_dir()
}

/// Build the module and return the path to it.
fn build_wasm(root: &Path, profile: &str) -> Result<PathBuf, String> {
    step(
        root,
        "cargo",
        &[
            "build",
            "-p",
            "tepsim-wasm",
            "--target",
            WASM_TARGET,
            "--profile",
            profile,
        ],
    )?;

    // Cargo puts the `dev` profile's output in `debug/` and every other
    // profile's under its own name. Special-casing the one exception is
    // shorter than asking cargo where it put things.
    let directory = if profile == "dev" { "debug" } else { profile };
    let module = root
        .join("target")
        .join(WASM_TARGET)
        .join(directory)
        .join("tepsim_wasm.wasm");
    if !module.is_file() {
        return Err(format!(
            "the build reported success but {} does not exist",
            module.display()
        ));
    }
    Ok(module)
}

/// Write the Node driver and return its path.
///
/// Generated rather than committed so it cannot drift from the export names it
/// calls, which live in `crates/tepsim-wasm/src/tier9.rs`.
fn write_driver(root: &Path) -> Result<PathBuf, String> {
    let dir = root.join(WORK_DIR);
    std::fs::create_dir_all(&dir)
        .map_err(|e| format!("could not create {}: {e}", dir.display()))?;
    let path = dir.join("tier9.mjs");
    std::fs::write(&path, DRIVER)
        .map_err(|e| format!("could not write {}: {e}", path.display()))?;
    Ok(path)
}

/// The Node driver.
///
/// Deliberately small and deliberately dependency-free: `node` and nothing
/// else, no `package.json`, no install step. It prints one line per case in a
/// format `parse_driver_output` reads with `split_whitespace`, because adding a
/// JSON parser to xtask for six lines of output would be worse than the
/// format.
///
/// There is no `Math.random`, no `Date`, and no floating-point arithmetic in
/// here. The values cross the boundary as wasm `i64`, which JavaScript sees as
/// `BigInt`, so nothing passes through a double on the way out.
const DRIVER: &str = r#"// Generated by `cargo xtask tier9`. Edits here are overwritten.
//
// Instantiates a tepsim-wasm module and reads the tier 9 table out of it
// through the glue-free `extern "C"` exports. No wasm-bindgen, no npm.
import { readFileSync } from 'node:fs';

const path = process.argv[2];
if (!path) {
  console.error('usage: node tier9.mjs <module.wasm>');
  process.exit(2);
}

const mod = await WebAssembly.compile(readFileSync(path));

// The module declares wasm-bindgen's placeholder imports because the rest of
// the crate is #[wasm_bindgen]. Nothing on the tier 9 path calls them, and a
// stub that threw is how we find out if that ever stops being true.
const imports = {};
for (const i of WebAssembly.Module.imports(mod)) {
  imports[i.module] ??= {};
  imports[i.module][i.name] = () => {
    throw new Error(`tier 9 reached JavaScript through ${i.module}.${i.name}`);
  };
}

const { exports: e } = await WebAssembly.instantiate(mod, imports);

const hex = (value) => BigInt.asUintN(64, value).toString(16).padStart(16, '0');

const count = e.tepsim_wasm_tier9_case_count();
const lines = [];
for (let i = 0; i < count; i += 1) {
  lines.push(`case ${i} ${hex(e.tepsim_wasm_tier9_case_digest(i))} ` +
             `${hex(e.tepsim_wasm_tier9_expected_digest(i))}`);
}
lines.push(`suite ${hex(e.tepsim_wasm_tier9_suite_digest())} ` +
           `${hex(e.tepsim_wasm_tier9_expected_suite_digest())}`);
console.log(lines.join('\n'));
"#;

fn run_driver(node: &Path, driver: &Path, module: &Path) -> Result<Digests, String> {
    println!("\n$ node {} {}", driver.display(), module.display());
    let out = Command::new(node)
        .arg(driver)
        .arg(module)
        .output()
        .map_err(|e| format!("failed to run node: {e}"))?;
    if !out.status.success() {
        return Err(format!(
            "the node driver failed:\n{}",
            String::from_utf8_lossy(&out.stderr)
        ));
    }
    parse_driver_output(&String::from_utf8_lossy(&out.stdout))
}

/// Read the driver's lines back.
///
/// Each `case` line carries the digest the module computed *and* the constant
/// compiled into it. Comparing the second against this process's own copy of
/// the table is what catches a stale `.wasm` left in `target/` by an older
/// commit, which would otherwise pass as agreement.
fn parse_driver_output(text: &str) -> Result<Digests, String> {
    let mut cases = Vec::new();
    let mut suite = None;

    for line in text.lines() {
        let fields: Vec<&str> = line.split_whitespace().collect();
        match fields.as_slice() {
            ["case", index, computed, expected] => {
                let index: usize = index
                    .parse()
                    .map_err(|_| format!("bad case index in `{line}`"))?;
                if index != cases.len() {
                    return Err(format!(
                        "the driver reported case {index} where case {} was \
                         expected, so the table is out of order",
                        cases.len()
                    ));
                }
                let case = CASES
                    .get(index)
                    .ok_or_else(|| format!("the module has a case {index} this build does not"))?;
                let expected = parse_hex(expected, line)?;
                if expected != case.digest {
                    return Err(format!(
                        "case `{}`: the module was built with {expected:016x} \
                         as the committed digest and this build has \
                         {:016x}. The .wasm is from a different commit; \
                         rebuild it.",
                        case.name, case.digest
                    ));
                }
                cases.push(parse_hex(computed, line)?);
            }
            ["suite", computed, expected] => {
                let expected = parse_hex(expected, line)?;
                if expected != tier9::SUITE_DIGEST {
                    return Err(format!(
                        "the module was built with {expected:016x} as the \
                         suite digest and this build has {:016x}. The .wasm \
                         is from a different commit; rebuild it.",
                        tier9::SUITE_DIGEST
                    ));
                }
                suite = Some(parse_hex(computed, line)?);
            }
            [] => {}
            _ => return Err(format!("unrecognised driver output: `{line}`")),
        }
    }

    if cases.len() != CASES.len() {
        return Err(format!(
            "the module reported {} cases and this build has {}",
            cases.len(),
            CASES.len()
        ));
    }
    let suite = suite.ok_or_else(|| "the driver reported no suite digest".to_string())?;
    Ok(Digests { cases, suite })
}

fn parse_hex(field: &str, line: &str) -> Result<u64, String> {
    u64::from_str_radix(field, 16).map_err(|_| format!("`{field}` is not a digest, in `{line}`"))
}

// ---------------------------------------------------------------------------
// reporting
// ---------------------------------------------------------------------------

fn print_table(digests: &Digests) {
    println!("| case | committed | computed | |");
    println!("|---|---|---|---|");
    for (case, computed) in CASES.iter().zip(&digests.cases) {
        println!(
            "| {:<18} | {:016x} | {:016x} | {} |",
            case.name,
            case.digest,
            computed,
            verdict(*computed == case.digest)
        );
    }
    println!(
        "| {:<18} | {:016x} | {:016x} | {} |",
        "suite",
        tier9::SUITE_DIGEST,
        digests.suite,
        verdict(digests.suite == tier9::SUITE_DIGEST)
    );
}

/// Words, not symbols. A reader scanning for a failure should not have to
/// distinguish two similar glyphs in a column of hex.
fn verdict(ok: bool) -> &'static str {
    if ok { "match" } else { "DIFFERS" }
}

fn mismatch_message(what: &str, digests: &Digests) -> String {
    let mut message = format!(
        "TIER 9 FAILED on {what}.\n\n\
         This is the most important result this project can produce: a \
         platform does not\nreproduce the digests committed in \
         crates/tepsim/src/tier9.rs. Do not update a\nconstant. Find out which \
         of the platform, the compiler or the model changed the\nnumbers.\n\n"
    );
    for (case, computed) in CASES.iter().zip(&digests.cases) {
        if *computed != case.digest {
            let _ = writeln!(
                message,
                "  {:<18} committed {:016x}, computed {:016x}   ({})",
                case.name, case.digest, computed, case.covers
            );
        }
    }
    if digests.suite != tier9::SUITE_DIGEST {
        let _ = writeln!(
            message,
            "  {:<18} committed {:016x}, computed {:016x}",
            "suite",
            tier9::SUITE_DIGEST,
            digests.suite
        );
    }
    message
}

/// What was compared, and what was not.
///
/// The second half matters more than the first. A reader who sees only "green"
/// will believe the claim covers every platform in `PLAN.org`, and it does not
/// until somebody runs this on each of them.
fn summary(native: &Digests, wasm: &[WasmRun], profiles: &[String]) {
    println!("\n--- what this run covers ---\n");
    println!(
        "  {:<24} {}",
        host_triple(),
        if native.agrees() {
            "matches the committed table"
        } else {
            "DOES NOT MATCH the committed table"
        }
    );
    for run in wasm {
        println!(
            "  {:<24} {} (under {}, {} bytes)",
            run.label,
            if run.digests.agrees() {
                "matches the committed table"
            } else {
                "DOES NOT MATCH the committed table"
            },
            run.runtime,
            run.bytes
        );
    }
    if wasm.is_empty() {
        println!(
            "  {:<24} not run: see the skip message above",
            "wasm32-unknown-unknown"
        );
    } else if wasm.len() < profiles.len() {
        println!("  (fewer profiles ran than were asked for)");
    }

    println!("\n--- what this run does not cover ---\n");
    println!(
        "  Any architecture other than {}. The digests are committed \
         constants, so running\n  this same command on an x86-64 machine \
         extends the claim to x86-64 with no code\n  change: it compares that \
         machine against what was measured here, not against\n  itself.",
        std::env::consts::ARCH
    );
    println!(
        "\n  A real browser. This drives node, which is V8, the same engine \
         Chrome runs; a\n  browser additionally exercises the wasm-bindgen \
         glue and the Web Worker\n  boundary. `PLAN.org` asks for a headless \
         browser and that is still owed. The\n  browser app's own page prints \
         the case 0 digest, so the check is available by\n  eye today."
    );
}

// ---------------------------------------------------------------------------
// flags and environment
// ---------------------------------------------------------------------------

/// `aarch64-apple-darwin`, or the best description available.
///
/// Taken from `rustc -vV`, which names the triple exactly as the target
/// directory does, rather than assembled from `std::env::consts`, which would
/// give `aarch64` and `macos` and leave a reader to guess the vendor.
fn host_triple() -> String {
    let fallback = || format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS);
    let Ok(out) = Command::new("rustc").arg("-vV").output() else {
        return fallback();
    };
    let Ok(text) = String::from_utf8(out.stdout) else {
        return fallback();
    };
    text.lines()
        .find_map(|line| line.strip_prefix("host: "))
        .map(str::to_string)
        .unwrap_or_else(fallback)
}

fn parse_profiles(flags: &[String]) -> Result<Vec<String>, String> {
    let Some(position) = flags.iter().position(|f| f == "--profiles") else {
        return Ok(DEFAULT_PROFILES.iter().map(|s| (*s).to_string()).collect());
    };
    let list = flags.get(position + 1).ok_or_else(|| {
        "--profiles needs a value, for example `--profiles release-wasm,dev`".to_string()
    })?;
    let profiles: Vec<String> = list
        .split(',')
        .map(|p| p.trim().to_string())
        .filter(|p| !p.is_empty())
        .collect();
    if profiles.is_empty() {
        return Err("--profiles was given no profile names".to_string());
    }
    Ok(profiles)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(index: usize, computed: u64) -> String {
        format!("case {index} {computed:016x} {:016x}", CASES[index].digest)
    }

    fn good_output() -> String {
        let mut text = String::new();
        for (index, case) in CASES.iter().enumerate() {
            text.push_str(&line(index, case.digest));
            text.push('\n');
        }
        let _ = writeln!(
            text,
            "suite {:016x} {:016x}",
            tier9::SUITE_DIGEST,
            tier9::SUITE_DIGEST
        );
        text
    }

    #[test]
    fn a_matching_transcript_parses_and_agrees() {
        let digests = parse_driver_output(&good_output()).expect("parses");
        assert!(digests.agrees());
        assert_eq!(digests.cases.len(), CASES.len());
    }

    /// The point of the whole command: a differing digest must survive parsing
    /// and be reported, not be rejected as malformed.
    #[test]
    fn a_differing_digest_parses_and_disagrees() {
        let text = good_output().replace(
            &format!("{:016x}", CASES[0].digest),
            // Both occurrences on the line would change if this replaced the
            // expected column too, so the constant chosen differs from every
            // committed digest and the assertion below checks the effect.
            "0000000000000001",
        );
        // The expected column changed as well, which is the stale-module case.
        assert!(parse_driver_output(&text).is_err());

        // Now only the computed column.
        let mut lines: Vec<String> = good_output().lines().map(str::to_string).collect();
        lines[0] = line(0, 0x0123_4567_89ab_cdef);
        let digests = parse_driver_output(&lines.join("\n")).expect("parses");
        assert!(!digests.agrees());
        assert_eq!(digests.cases[0], 0x0123_4567_89ab_cdef);
        assert!(mismatch_message("wasm32", &digests).contains(CASES[0].name));
    }

    /// A `.wasm` built from another commit carries that commit's constants.
    /// Accepting it would report agreement between two different tables.
    #[test]
    fn a_stale_module_is_rejected() {
        let text = format!(
            "case 0 {:016x} 0000000000000002\nsuite {:016x} {:016x}\n",
            CASES[0].digest,
            tier9::SUITE_DIGEST,
            tier9::SUITE_DIGEST
        );
        let error = parse_driver_output(&text).expect_err("stale");
        assert!(error.contains("different commit"), "{error}");
    }

    #[test]
    fn a_short_transcript_is_rejected() {
        let text = format!("case 0 {:016x} {:016x}\n", CASES[0].digest, CASES[0].digest);
        assert!(parse_driver_output(&text).is_err());
    }

    #[test]
    fn out_of_order_cases_are_rejected() {
        let text = format!("case 1 {:016x} {:016x}\n", CASES[1].digest, CASES[1].digest);
        let error = parse_driver_output(&text).expect_err("out of order");
        assert!(error.contains("out of order"), "{error}");
    }

    #[test]
    fn profiles_default_and_parse() {
        assert_eq!(parse_profiles(&[]).expect("the default"), DEFAULT_PROFILES);
        let flags = vec!["--profiles".to_string(), "dev, release".to_string()];
        assert_eq!(
            parse_profiles(&flags).expect("a list"),
            vec!["dev", "release"]
        );
        let flags = vec!["--profiles".to_string()];
        assert!(parse_profiles(&flags).is_err());
    }

    #[test]
    fn the_browser_digest_is_read_from_an_assignment_only() {
        assert_eq!(
            find_browser_digest("const EXPECTED_SELF_CHECK = \"c8a26889992f1719\";"),
            Some("c8a26889992f1719".to_string())
        );
        // Single quotes, and no space around the equals.
        assert_eq!(
            find_browser_digest("      const EXPECTED_SELF_CHECK='0123456789abcdef';"),
            Some("0123456789abcdef".to_string())
        );
        // A use is not an assignment, and must not stop the search.
        assert_eq!(
            find_browser_digest(
                "if (d === EXPECTED_SELF_CHECK) {}\n\
                 const EXPECTED_SELF_CHECK = \"c8a26889992f1719\";"
            ),
            Some("c8a26889992f1719".to_string())
        );
        assert_eq!(find_browser_digest("`${EXPECTED_SELF_CHECK}`"), None);
        // Too short, and not hex.
        assert_eq!(
            find_browser_digest("const EXPECTED_SELF_CHECK = \"c8a2\";"),
            None
        );
        assert_eq!(
            find_browser_digest("const EXPECTED_SELF_CHECK = \"not-a-digest-xx\";"),
            None
        );
        assert_eq!(find_browser_digest("nothing here"), None);
    }

    /// The browser app's hard-coded digest against the table, without needing a
    /// wasm build. This is the check that would otherwise only happen when
    /// somebody opened the page.
    #[test]
    fn this_checkouts_browser_copies_match_the_table() {
        check_browser_copies(&crate::workspace_root()).expect("the browser copies agree");
    }

    /// The generated driver has to call the exports the crate actually has.
    /// Nothing else checks the two against each other, because one is a string
    /// and the other is Rust.
    #[test]
    fn the_driver_calls_the_exports_that_exist() {
        for name in [
            "tepsim_wasm_tier9_case_count",
            "tepsim_wasm_tier9_case_digest",
            "tepsim_wasm_tier9_expected_digest",
            "tepsim_wasm_tier9_suite_digest",
            "tepsim_wasm_tier9_expected_suite_digest",
        ] {
            assert!(DRIVER.contains(name), "the driver never calls {name}");
        }
        // And nothing that would make it non-deterministic.
        for banned in ["Math.random", "Date.now", "new Date"] {
            assert!(!DRIVER.contains(banned), "the driver uses {banned}");
        }
    }
}
