//! The committed golden trace must still be what the Fortran produces.
//!
//! This is the check that tells a genuine regression apart from a toolchain
//! change. It needs gfortran, so it runs in the non-fast gate and in CI on
//! Linux and macOS. The complementary half, `cargo xtask fidelity`, validates
//! the same file with no Fortran toolchain at all.

#![cfg(feature = "oracle")]

use std::path::PathBuf;

use tepsim_oracle::golden::{self, Step, Trace};
use tepsim_oracle::{Oracle, build_info};

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn committed() -> Trace {
    let path = workspace_root().join(golden::PATH);
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        panic!(
            "reading {}: {e}. Generate it with the gen-golden-trace bin.",
            path.display()
        )
    });
    Trace::parse(&text).unwrap_or_else(|e| panic!("parsing {}: {e}", path.display()))
}

/// Re-run exactly what `gen-golden-trace` runs.
fn regenerate() -> Vec<Step> {
    let mut oracle = Oracle::lock();
    let (mut time, mut yy) = oracle.init();
    oracle.set_disturbances(&[0; 20]);
    oracle.set_rng(golden::SEED);

    let mut steps = Vec::with_capacity(golden::STEPS);
    for _ in 0..golden::STEPS {
        let states = yy;
        let derivatives = oracle.derivatives(time, &yy);
        steps.push(Step {
            states,
            derivatives,
            measurements: oracle.measurements(),
            rng: oracle.rng(),
        });
        time += golden::DT_HOURS;
        for (y, dy) in yy.iter_mut().zip(&derivatives) {
            *y += dy * golden::DT_HOURS;
        }
    }
    steps
}

#[test]
fn the_fortran_still_reproduces_the_committed_trace() {
    let trace = committed();
    trace.require_full_length().expect("trace length");

    if trace.gfortran != build_info::GFORTRAN_VERSION {
        panic!(
            "toolchain changed: the trace was recorded with gfortran {} but this \
             machine has {}.\n\nThis is a re-baseline, not a regression. Per \
             CLAUDE.md, regenerate the trace and re-record the affected \
             validation numbers in LOG.org in the same commit, as its own logged \
             event.",
            trace.gfortran,
            build_info::GFORTRAN_VERSION
        );
    }
    assert_eq!(
        trace.fflags,
        build_info::FORTRAN_FLAGS,
        "the Fortran flag set changed since the trace was recorded, which \
         invalidates it the same way a compiler change would"
    );

    let fresh = regenerate();
    assert_eq!(fresh.len(), trace.steps.len());

    let mut mismatches = Vec::new();
    for (i, (got, want)) in fresh.iter().zip(&trace.steps).enumerate() {
        if got == want {
            continue;
        }
        let mut which = Vec::new();
        for (j, (a, b)) in got.states.iter().zip(&want.states).enumerate() {
            if a.to_bits() != b.to_bits() {
                which.push(format!("state[{j}] {a:?} vs {b:?}"));
            }
        }
        for (j, (a, b)) in got.derivatives.iter().zip(&want.derivatives).enumerate() {
            if a.to_bits() != b.to_bits() {
                which.push(format!("deriv[{j}] {a:?} vs {b:?}"));
            }
        }
        for (j, (a, b)) in got.measurements.iter().zip(&want.measurements).enumerate() {
            if a.to_bits() != b.to_bits() {
                which.push(format!("meas[{j}] {a:?} vs {b:?}"));
            }
        }
        if got.rng.to_bits() != want.rng.to_bits() {
            which.push(format!("rng {:?} vs {:?}", got.rng, want.rng));
        }
        which.truncate(4);
        mismatches.push(format!("  step {i}: {}", which.join(", ")));
        if mismatches.len() >= 5 {
            break;
        }
    }

    assert!(
        mismatches.is_empty(),
        "the Fortran no longer reproduces the committed golden trace:\n{}\n\n\
         With the same compiler and the same flags this is deterministic, so a \
         difference means something changed that should not have. Investigate \
         before regenerating.",
        mismatches.join("\n")
    );
}

/// The trace has to be a plausible run, not just a well-formed file.
#[test]
fn the_trace_records_a_plausible_run() {
    let trace = committed();

    assert_eq!(trace.seed.to_bits(), golden::SEED.to_bits());
    assert_eq!(trace.dt_hours.to_bits(), golden::DT_HOURS.to_bits());

    // Reactor temperature, XMEAS(9), should sit near the nominal 120.4 C for a
    // 100-second open-loop run with no disturbances and no control.
    for (i, step) in trace.steps.iter().enumerate() {
        let temp = step.measurements[8];
        assert!(
            (119.0..=122.0).contains(&temp),
            "step {i}: reactor temperature {temp} C left the nominal band"
        );
    }

    // The generator must advance every step: each one draws measurement noise.
    for pair in trace.steps.windows(2) {
        assert_ne!(
            pair[0].rng.to_bits(),
            pair[1].rng.to_bits(),
            "the generator must advance on every step past t=0"
        );
    }
}

/// Round-trips through the text format without losing a bit.
#[test]
fn the_trace_format_round_trips_exactly() {
    let trace = committed();
    let reparsed = Trace::parse(&trace.to_text()).expect("re-parse");
    assert_eq!(
        reparsed, trace,
        "hex bit patterns must round-trip exactly; if they do not, the format \
         is unfit for anchoring a validation ladder"
    );
}
