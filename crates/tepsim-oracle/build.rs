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

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

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

    // Only teprob.f. temain.f and temain_mod.f each declare a PROGRAM, which
    // would collide with the Rust test harness's own entry point.
    let object = out_dir.join("teprob.o");
    run(
        Command::new(gfortran())
            .args(FORTRAN_FLAGS)
            .arg(&source)
            .arg("-o")
            .arg(&object),
        "compiling teprob.f",
    );

    let archive = out_dir.join("libteoracle.a");
    let _ = std::fs::remove_file(&archive);
    run(
        Command::new("ar").arg("crs").arg(&archive).arg(&object),
        "archiving teprob.o",
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
