//! Checks [`TepRng`] against output recorded from the original Fortran.
//!
//! This is the half of the RNG check that needs no Fortran toolchain: the
//! vectors in `golden/tesub7-vectors.txt` were recorded from the original and
//! are committed, so a contributor with no gfortran still gets a bit-exact
//! check. The complementary half, `crates/tepsim-oracle/tests/rng_differential.rs`,
//! compares draw by draw against the live Fortran.
//!
//! An integration test rather than a unit test because `tepsim-core` is
//! `no_std` and cannot read files; integration tests are separate crates with
//! std available.

use std::collections::BTreeMap;
use std::path::PathBuf;

use tepsim_core::TepRng;

/// One seed's recorded expectations.
#[derive(Debug)]
struct Vectors {
    seed: f64,
    draws: usize,
    fold: u64,
    final_state: u64,
    unit: Vec<u64>,
    signed: Vec<u64>,
}

fn vectors_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../golden/tesub7-vectors.txt")
}

fn hex(field: &str) -> u64 {
    u64::from_str_radix(field, 16).unwrap_or_else(|e| panic!("bad hex {field:?}: {e}"))
}

fn load() -> Vec<Vectors> {
    let path = vectors_path();
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "reading {}: {e}\nRegenerate with: cargo run -p tepsim-oracle \
             --features oracle --bin gen-rng-vectors",
            path.display()
        )
    });

    let mut all = Vec::new();
    let mut current: BTreeMap<&str, &str> = BTreeMap::new();
    let push = |fields: &BTreeMap<&str, &str>, all: &mut Vec<Vectors>| {
        if fields.is_empty() {
            return;
        }
        all.push(Vectors {
            seed: f64::from_bits(hex(fields["seed"])),
            draws: fields["draws"].parse().expect("draw count"),
            fold: hex(fields["fold"]),
            final_state: hex(fields["final"]),
            unit: fields["unit"].split_whitespace().map(hex).collect(),
            signed: fields["signed"].split_whitespace().map(hex).collect(),
        });
    };

    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (key, value) = line.split_once(' ').expect("key and value");
        if key == "seed" && current.contains_key("seed") {
            push(&current, &mut all);
            current.clear();
        }
        current.insert(key, value.trim());
    }
    push(&current, &mut all);

    assert!(!all.is_empty(), "the vector file records no seeds at all");
    all
}

#[test]
fn the_state_sequence_matches_the_fortran_for_every_recorded_seed() {
    for v in load() {
        let mut rng = TepRng::new(v.seed);
        let mut fold = 0_u64;
        let mut last = 0_u64;
        for _ in 0..v.draws {
            rng.unit();
            last = rng.state().to_bits();
            fold ^= last;
        }
        assert_eq!(
            fold, v.fold,
            "seed {}: the XOR fold over {} draws differs from the Fortran. \
             Some state in the run is wrong, though not necessarily the last.",
            v.seed, v.draws
        );
        assert_eq!(
            last, v.final_state,
            "seed {}: the state after {} draws differs from the Fortran",
            v.seed, v.draws
        );
    }
}

#[test]
fn both_output_modes_match_the_fortran() {
    for v in load() {
        let mut rng = TepRng::new(v.seed);
        for (i, want) in v.unit.iter().enumerate() {
            let got = rng.unit().to_bits();
            assert_eq!(
                got,
                *want,
                "seed {}: unit draw {i} is {:?}, Fortran gives {:?}",
                v.seed,
                f64::from_bits(got),
                f64::from_bits(*want)
            );
        }

        let mut rng = TepRng::new(v.seed);
        for (i, want) in v.signed.iter().enumerate() {
            let got = rng.signed().to_bits();
            assert_eq!(
                got,
                *want,
                "seed {}: signed draw {i} is {:?}, Fortran gives {:?}",
                v.seed,
                f64::from_bits(got),
                f64::from_bits(*want)
            );
        }
    }
}

/// The compiled-in seed is the one every published dataset and the golden trace
/// start from, so it gets the deepest check.
#[test]
fn the_compiled_in_seed_is_checked_to_ten_million_draws() {
    let v = load()
        .into_iter()
        .find(|v| v.seed.to_bits() == TepRng::DEFAULT_SEED.to_bits())
        .expect("the compiled-in seed must be among the recorded vectors");
    assert_eq!(
        v.draws, 10_000_000,
        "the primary seed should be recorded to 10^7 draws"
    );
}
