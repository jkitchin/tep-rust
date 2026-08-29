//! Compiles the original Tennessee Eastman Fortran and links it into the
//! development-only oracle crate.
//!
//! Two rules govern this file.
//!
//! The vendored source under `reference/` is ground truth and is never copied
//! or edited. It is compiled in place, straight out of the repository.
//!
//! The flags are pinned and asserted by a test. Reference numbers for Tier 1
//! and Tier 2 are measured against a specific compiler and a specific flag set;
//! changing either invalidates every number recorded in `LOG.org`, so it is a
//! deliberate re-baseline rather than a casual edit.
//!
//! # This does not build on Windows, and that is the policy
//!
//! `xtask::oracle_supported` refuses to build the oracle on Windows, so nothing
//! in CI reaches this file there. The policy is `CLAUDE.md`'s: the oracle runs
//! on Linux and macOS runners only. Every Tier 1 and Tier 2 number was
//! baselined against a specific gfortran on one of those two, and admitting a
//! third toolchain is a re-baseline, not a portability fix.
//!
//! For whoever does decide to port it, the first thing in the way is a path.
//! `fs::canonicalize` on Windows returns an extended-length `\\?\` path, and
//! MinGW's `f951.exe` cannot parse one: it reported
//! `Fatal Error: Cannot open file '\\teprob.f'` the first time CI ran this on
//! `windows-latest`. Stripping the prefix before handing the path to gfortran
//! is the fix, and it is the beginning of the work rather than the end of it.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

#[path = "instrument.rs"]
mod instrument;

/// Flags used to compile the Fortran. Deliberately conservative.
///
/// `-O0` keeps the generated code as close to the source as possible, which
/// matters because the oracle's whole purpose is to be an unambiguous
/// statement of what the original computes. There is no `-ffast-math`, no
/// `-funsafe-math-optimizations`, and no `-fassociative-math`: any of those
/// would let the compiler reassociate floating-point expressions and quietly
/// change the reference values we validate against.
///
/// Kept in sync with `FORTRAN_FLAGS` in `src/build_info.rs`, which a test
/// asserts against this list.
const FORTRAN_FLAGS: &[&str] = &[
    "-c",
    "-O0",
    "-fno-fast-math",
    "-fno-unsafe-math-optimizations",
    "-fPIC",
    "-std=legacy",
];

fn main() {
    // Without the feature, this crate is an empty shell and must not require a
    // Fortran toolchain to build.
    if env::var_os("CARGO_FEATURE_ORACLE").is_none() {
        println!("cargo:rerun-if-changed=build.rs");
        return;
    }

    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR"));
    let reference = manifest
        .join("../../reference/fortran")
        .canonicalize()
        .expect("reference/fortran must exist; see B-0003");
    let source = reference.join("teprob.f");
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR"));

    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-changed={}", source.display());

    // The vendored file is ground truth and is checksummed, so the shutdown
    // flag is exposed by rewriting into OUT_DIR rather than by editing it. Each
    // rewrite asserts its own pre-image; see instrument.rs.
    let pristine = std::fs::read_to_string(&source)
        .unwrap_or_else(|e| panic!("reading {}: {e}", source.display()));
    let instrumented = instrument::instrument(&pristine);
    let instrumented_path = out_dir.join("teprob_instrumented.f");
    std::fs::write(&instrumented_path, &instrumented).expect("writing instrumented source");

    let object = out_dir.join("teprob.o");
    run(
        Command::new(gfortran())
            .args(FORTRAN_FLAGS)
            .arg(&instrumented_path)
            .arg("-o")
            .arg(&object),
        "compiling the instrumented teprob.f",
    );

    // The closed-loop driver, for Phase 4. Its nineteen `CONTRLn` subroutines
    // are what the control layer is ported from, and they are only linkable
    // once the file's unnamed main program stops defining `main`; see
    // instrument.rs. `temain.f` stays out: it is the older single-loop driver
    // and nothing in this project ports from it.
    let driver_source = reference.join("temain_mod.f");
    println!("cargo:rerun-if-changed={}", driver_source.display());
    let driver_pristine = std::fs::read_to_string(&driver_source)
        .unwrap_or_else(|e| panic!("reading {}: {e}", driver_source.display()));
    let driver_path = out_dir.join("temain_mod_instrumented.f");
    std::fs::write(
        &driver_path,
        instrument::instrument_driver(&driver_pristine),
    )
    .expect("writing instrumented driver");
    let driver_object = out_dir.join("temain_mod.o");
    run(
        Command::new(gfortran())
            .args(FORTRAN_FLAGS)
            .arg(&driver_path)
            .arg("-o")
            .arg(&driver_object),
        "compiling the instrumented temain_mod.f",
    );

    let archive = out_dir.join("libteoracle.a");
    let _ = std::fs::remove_file(&archive);
    run(
        Command::new("ar")
            .arg("crs")
            .arg(&archive)
            .arg(&object)
            .arg(&driver_object),
        "archiving teprob.o and temain_mod.o",
    );

    // Build a standalone probe from the *pristine* source so a test can prove
    // the instrumentation changed no numbers. Linking both copies into one
    // process is impossible, since they define the same symbols, so the
    // pristine one runs as a separate executable and reports its answers as
    // raw bit patterns. Comparing bits rather than decimal text keeps the
    // comparison exact.
    let probe_src = out_dir.join("pristine_probe.f");
    std::fs::write(&probe_src, PRISTINE_PROBE).expect("writing probe driver");
    let probe_bin = out_dir.join("pristine_probe");
    run(
        Command::new(gfortran())
            .args(FORTRAN_FLAGS.iter().filter(|f| **f != "-c"))
            .arg(&source)
            .arg(&probe_src)
            .arg("-o")
            .arg(&probe_bin),
        "building the pristine probe",
    );
    println!(
        "cargo:rustc-env=TEP_ORACLE_PRISTINE_PROBE={}",
        probe_bin.display()
    );

    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static=teoracle");

    // The Fortran runtime. gfortran knows where its own library lives; asking
    // it beats guessing at Homebrew and distribution layouts.
    if let Some(dir) = gfortran_lib_dir() {
        println!("cargo:rustc-link-search=native={}", dir.display());
    }
    println!("cargo:rustc-link-lib=dylib=gfortran");

    // Record what actually built this, so the golden trace in B-0004c can carry
    // it and a mismatch can be reported rather than silently tolerated.
    println!(
        "cargo:rustc-env=TEP_ORACLE_FORTRAN_FLAGS={}",
        FORTRAN_FLAGS.join(" ")
    );
    println!(
        "cargo:rustc-env=TEP_ORACLE_GFORTRAN_VERSION={}",
        gfortran_version()
    );
}

fn gfortran() -> String {
    env::var("FC").unwrap_or_else(|_| "gfortran".to_string())
}

fn run(cmd: &mut Command, what: &str) {
    let status = cmd.status().unwrap_or_else(|e| {
        panic!(
            "{what}: could not run the command: {e}.\n\
             The oracle feature requires a Fortran compiler. Set FC to override \
             the default `gfortran`."
        )
    });
    assert!(status.success(), "{what}: command failed with {status}");
}

/// Ask gfortran where its runtime lives, then take the containing directory.
fn gfortran_lib_dir() -> Option<PathBuf> {
    for name in ["libgfortran.dylib", "libgfortran.so"] {
        let output = Command::new(gfortran())
            .arg(format!("-print-file-name={name}"))
            .output()
            .ok()?;
        let path = String::from_utf8(output.stdout).ok()?;
        let path = Path::new(path.trim());
        // gfortran echoes the bare name back when it cannot find the file.
        if path.is_absolute() {
            if let Some(parent) = path.parent() {
                return parent
                    .canonicalize()
                    .ok()
                    .or_else(|| Some(parent.to_path_buf()));
            }
        }
    }
    None
}

fn gfortran_version() -> String {
    Command::new(gfortran())
        .arg("-dumpfullversion")
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| "unknown".to_string())
}

/// A driver compiled against the *pristine* Fortran, used to prove the
/// instrumentation is behaviour-preserving.
///
/// It runs the same sequence the equivalence test runs in-process against the
/// instrumented build: initialise, pin the generator, take one step past zero.
/// Results are emitted as raw IEEE bit patterns via `TRANSFER`, because decimal
/// formatting would not round-trip exactly and this comparison has to be exact.
const PRISTINE_PROBE: &str = r#"
      PROGRAM PROBE
      INTEGER NN,I
      DOUBLE PRECISION TIME, YY(50), YP(50)
      DOUBLE PRECISION G
      COMMON/RANDSD/G
      INTEGER*8 IB
      NN=50
      CALL TEINIT(NN,TIME,YY,YP)
      G=4651207995.D0
      TIME=1.D0/3600.D0
      CALL TEFUNC(NN,TIME,YY,YP)
      DO 10 I=1,50
      IB=TRANSFER(YP(I),IB)
      WRITE(*,*) IB
   10 CONTINUE
      IB=TRANSFER(G,IB)
      WRITE(*,*) IB
      END
"#;
