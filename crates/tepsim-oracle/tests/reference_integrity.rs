//! The vendored reference material is ground truth for every validation tier.
//!
//! An accidental edit to `teprob.f` or to a `.dat` file would silently
//! invalidate every recorded validation number while leaving the build green,
//! which is the worst possible failure mode for this project. So the bytes are
//! pinned and checked on every gate run.
//!
//! If this test fails, restore the file from upstream. Do not regenerate
//! `CHECKSUMS.sha256` to make it pass: that converts a caught error into an
//! accepted one.

use std::fs;
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

fn reference_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../reference")
        .canonicalize()
        .expect("reference/ must exist; see B-0003")
}

/// One `<hex>  <relative path>` line of the manifest.
fn parse_manifest(text: &str) -> Vec<(String, String)> {
    text.lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .filter_map(|line| {
            let (hash, path) = line.split_once("  ")?;
            Some((hash.to_ascii_lowercase(), path.trim().to_string()))
        })
        .collect()
}

fn sha256_of(path: &Path) -> String {
    let bytes = fs::read(path).unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    let digest = Sha256::digest(&bytes);
    digest.iter().map(|b| format!("{b:02x}")).collect()
}

#[test]
fn every_vendored_file_matches_its_checksum() {
    let root = reference_dir();
    let manifest_path = root.join("CHECKSUMS.sha256");
    let manifest = fs::read_to_string(&manifest_path).expect("CHECKSUMS.sha256 must exist");
    let entries = parse_manifest(&manifest);

    assert!(
        !entries.is_empty(),
        "CHECKSUMS.sha256 lists nothing. An empty manifest would let this test \
         pass while checking absolutely nothing."
    );

    let mut mismatched = Vec::new();
    for (expected, rel) in &entries {
        let path = root.join(rel);
        assert!(
            path.is_file(),
            "{rel} is listed in the manifest but missing"
        );
        let actual = sha256_of(&path);
        if &actual != expected {
            mismatched.push(format!(
                "  {rel}\n    expected {expected}\n    actual   {actual}"
            ));
        }
    }
    assert!(
        mismatched.is_empty(),
        "vendored reference material has been modified:\n{}\n\nRestore these files. \
         Do not regenerate the manifest.",
        mismatched.join("\n")
    );
}

/// The manifest must cover everything present, not just the files it happens to
/// list. Otherwise a new or renamed file slips in unchecked.
#[test]
fn the_manifest_covers_every_file_present() {
    let root = reference_dir();
    let manifest = fs::read_to_string(root.join("CHECKSUMS.sha256")).expect("manifest");
    let listed: Vec<String> = parse_manifest(&manifest)
        .into_iter()
        .map(|(_, path)| path)
        .collect();

    let mut unlisted = Vec::new();
    for sub in ["fortran", "data"] {
        let dir = root.join(sub);
        for entry in fs::read_dir(&dir).unwrap_or_else(|e| panic!("reading {sub}: {e}")) {
            let entry = entry.expect("dir entry");
            if !entry.path().is_file() {
                continue;
            }
            let name = entry.file_name().to_string_lossy().to_string();
            let rel = format!("{sub}/{name}");
            if !listed.contains(&rel) {
                unlisted.push(rel);
            }
        }
    }
    unlisted.sort();
    assert!(
        unlisted.is_empty(),
        "present in reference/ but absent from CHECKSUMS.sha256: {unlisted:?}"
    );
}

/// Guards the two traps documented in `reference/README.org`. A loader written
/// against upstream's `data/README.md` gets both of these wrong.
#[test]
fn d00_training_set_is_transposed_with_500_samples() {
    let root = reference_dir();
    let text = fs::read_to_string(root.join("data/d00.dat")).expect("d00.dat");
    let rows = text.lines().filter(|l| !l.trim().is_empty()).count();
    let cols = text
        .lines()
        .find(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().count())
        .unwrap_or(0);

    assert_eq!(
        rows, 52,
        "d00.dat is stored transposed: 52 rows of variables"
    );
    assert_eq!(
        cols, 500,
        "d00.dat holds 500 samples, not the 480 that upstream's data/README.md claims"
    );
}

#[test]
fn every_other_dataset_is_samples_by_52() {
    let root = reference_dir();
    for fault in 1..=21 {
        for (suffix, expected_rows) in [("", 480), ("_te", 960)] {
            let name = format!("data/d{fault:02}{suffix}.dat");
            let text =
                fs::read_to_string(root.join(&name)).unwrap_or_else(|e| panic!("{name}: {e}"));
            let rows = text.lines().filter(|l| !l.trim().is_empty()).count();
            let cols = text
                .lines()
                .find(|l| !l.trim().is_empty())
                .map(|l| l.split_whitespace().count())
                .unwrap_or(0);
            assert_eq!(rows, expected_rows, "{name} row count");
            assert_eq!(cols, 52, "{name} column count");
        }
    }
    // d00_te follows the convention even though d00 does not.
    let text = fs::read_to_string(root.join("data/d00_te.dat")).expect("d00_te.dat");
    assert_eq!(text.lines().filter(|l| !l.trim().is_empty()).count(), 960);
}
