//! Records `TESUB7` output from the original Fortran, so the Rust port can be
//! checked on machines with no Fortran toolchain.
//!
//! ```text
//! cargo run -p tepsim-oracle --features oracle --bin gen-rng-vectors
//! ```
//!
//! Like the golden trace, regenerating this is a deliberate re-baseline.
//!
//! # What is recorded, and why it is enough
//!
//! The state sequence *is* the algorithm: both output modes are pure functions
//! of the state. So per seed we record the state after every draw folded into
//! one word by XOR, plus the final state. The fold catches a divergence
//! anywhere in the run, and the final state pins the endpoint. Sixteen raw
//! outputs of each mode then pin the two output formulas.
//!
//! That is a few hundred bytes per seed instead of the 160 MB a full 10^7-draw
//! dump would take, and it fails just as loudly.

use std::fmt::Write as _;
use std::path::PathBuf;

use tepsim_oracle::{Oracle, build_info};

/// Seeds to record, from `reference/README.org`.
///
/// The compiled-in seed gets ten million draws; the dataset seeds get a million
/// each, which is ample to catch a divergence that the fold would see at all.
const SEEDS: &[(&str, f64, usize)] = &[
    ("compiled-in, teprob.f:1187", 4_651_207_995.0, 10_000_000),
    ("original", 1_431_655_765.0, 1_000_000),
    ("d00_tr", 4_243_534_565.0, 1_000_000),
    ("d00_tr_new", 5_687_912_315.0, 1_000_000),
    ("d01_tr", 7_854_912_354.0, 1_000_000),
    ("d00_te", 1_254_545_354.0, 1_000_000),
];

/// How many raw outputs of each mode to record per seed.
const SAMPLES: usize = 16;

fn main() {
    let mut oracle = Oracle::lock();
    let mut out = String::new();

    out.push_str("# tep-rust TESUB7 vectors, format 1\n");
    out.push_str("#\n");
    out.push_str("# Recorded from the ORIGINAL Fortran. The Rust TepRng is checked against\n");
    out.push_str("# these, so they must never be regenerated from the Rust implementation:\n");
    out.push_str("# that would make the test compare the port with itself.\n");
    out.push_str("#\n");
    out.push_str("# Regenerate with:\n");
    out.push_str("#   cargo run -p tepsim-oracle --features oracle --bin gen-rng-vectors\n");
    out.push_str("#\n");
    let _ = writeln!(out, "# gfortran: {}", build_info::GFORTRAN_VERSION);
    let _ = writeln!(out, "# fflags: {}", build_info::FORTRAN_FLAGS);
    out.push_str("# fold is the XOR of the bit patterns of every state in the run.\n");
    out.push_str("# All values are hexadecimal IEEE-754 f64 bit patterns.\n");

    for (label, seed, draws) in SEEDS {
        let _ = writeln!(out, "\n# {label}");
        let _ = writeln!(out, "seed {:016x}", seed.to_bits());
        let _ = writeln!(out, "draws {draws}");

        // The state sequence. `unit` advances the state exactly once, and the
        // state after the call is what the fold accumulates.
        oracle.set_rng(*seed);
        let mut fold: u64 = 0;
        let mut final_state = 0.0_f64;
        for _ in 0..*draws {
            let _ = oracle.tesub7(1);
            final_state = oracle.rng();
            fold ^= final_state.to_bits();
        }
        let _ = writeln!(out, "fold {fold:016x}");
        let _ = writeln!(out, "final {:016x}", final_state.to_bits());

        // The two output scalings, from a fresh seed each time.
        oracle.set_rng(*seed);
        let _ = write!(out, "unit");
        for _ in 0..SAMPLES {
            let _ = write!(out, " {:016x}", oracle.tesub7(1).to_bits());
        }
        out.push('\n');

        oracle.set_rng(*seed);
        let _ = write!(out, "signed");
        for _ in 0..SAMPLES {
            let _ = write!(out, " {:016x}", oracle.tesub7(-1).to_bits());
        }
        out.push('\n');

        println!("{label:28} {draws:>9} draws  fold {fold:016x}");
    }

    let path = workspace_root().join("golden/tesub7-vectors.txt");
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("creating golden/");
    }
    std::fs::write(&path, &out).expect("writing vectors");
    println!("\nwrote {} bytes to {}", out.len(), path.display());
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}
