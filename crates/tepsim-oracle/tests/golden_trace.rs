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

/// Agreement required when the local gfortran is *not* the one that recorded
/// the trace.
///
/// Bit-equality is only meaningful against the recording compiler. Different
/// gfortran versions, and different platform `libm` implementations behind
/// them, evaluate `exp` and `pow` to different last bits, and 100 Euler steps
/// amplify that.
///
/// This bound is a measurement, not a preference. Observed over this 100-step
/// trace, against a trace recorded with gfortran 15.2.0 on macOS:
///
/// | Running | Platform libm | Values differing | Worst relative deviation |
/// |---------|---------------|------------------|--------------------------|
/// | 15.2.0  | Apple         | 0 of 14200       | 0, bit exact             |
/// | 16.1.0  | Apple         | under the bound  | below 1e-9               |
/// | 13.3.0  | glibc         | 16 of 14200      | 3.307e-8 at step 91      |
///
/// The pattern is the point. A different compiler on the *same* platform barely
/// moves anything, because it calls the same `libm`. A different platform moves
/// the last bits of `exp` and `pow`, and 91 Euler steps amplify that to 3.3e-8.
/// That is the mechanism `PLAN.org` predicts will limit cross-platform
/// agreement, measured rather than assumed, and it is why Tier 4 is diagnostic
/// while Tiers 5 and 6 are the gates.
///
/// Set at roughly three times the largest observed figure. Loosening it hides
/// exactly the regressions this test exists to catch, so it is only ever raised
/// alongside a recorded measurement.
const CROSS_COMPILER_TOLERANCE: f64 = 1e-7;

/// Relative difference, falling back to absolute near zero.
fn deviation(a: f64, b: f64) -> f64 {
    let scale = a.abs().max(b.abs());
    if scale < 1e-300 {
        0.0
    } else {
        (a - b).abs() / scale
    }
}

/// Where the largest disagreement was, for reporting.
#[derive(Default)]
struct Worst {
    deviation: f64,
    what: String,
}

impl Worst {
    fn consider(&mut self, a: f64, b: f64, what: impl FnOnce() -> String) {
        let d = deviation(a, b);
        if d > self.deviation {
            self.deviation = d;
            self.what = what();
        }
    }
}

#[test]
fn the_fortran_still_reproduces_the_committed_trace() {
    let trace = committed();
    trace.require_full_length().expect("trace length");

    let same_compiler = trace.gfortran == build_info::GFORTRAN_VERSION;
    assert_eq!(
        trace.fflags,
        build_info::FORTRAN_FLAGS,
        "the Fortran flag set changed since the trace was recorded, which \
         invalidates it the same way a compiler change would"
    );

    let fresh = regenerate();
    assert_eq!(fresh.len(), trace.steps.len());

    let mut worst = Worst::default();
    let mut exact_mismatches = 0_usize;
    for (i, (got, want)) in fresh.iter().zip(&trace.steps).enumerate() {
        for (j, (a, b)) in got.states.iter().zip(&want.states).enumerate() {
            worst.consider(*a, *b, || format!("step {i} state[{j}]"));
            exact_mismatches += usize::from(a.to_bits() != b.to_bits());
        }
        for (j, (a, b)) in got.derivatives.iter().zip(&want.derivatives).enumerate() {
            worst.consider(*a, *b, || format!("step {i} deriv[{j}]"));
            exact_mismatches += usize::from(a.to_bits() != b.to_bits());
        }
        for (j, (a, b)) in got.measurements.iter().zip(&want.measurements).enumerate() {
            worst.consider(*a, *b, || format!("step {i} meas[{j}]"));
            exact_mismatches += usize::from(a.to_bits() != b.to_bits());
        }
        exact_mismatches += usize::from(got.rng.to_bits() != want.rng.to_bits());
    }

    // Always report, pass or fail. A number in the CI log is what lets the next
    // session tell drift from noise.
    println!(
        "golden trace: recorded with gfortran {}, running {}; \
         {} of {} values differ in bits, worst relative deviation {:.3e} at {}",
        trace.gfortran,
        build_info::GFORTRAN_VERSION,
        exact_mismatches,
        trace.steps.len() * 142,
        worst.deviation,
        if worst.what.is_empty() {
            "nowhere"
        } else {
            &worst.what
        }
    );

    if same_compiler {
        assert_eq!(
            exact_mismatches, 0,
            "the same gfortran must reproduce the trace bit for bit; {} values \
             differ, worst {:.3e} at {}. This is deterministic, so a difference \
             means something changed that should not have.",
            exact_mismatches, worst.deviation, worst.what
        );
    } else {
        // The generator is pure f64 arithmetic with no library calls, so it
        // must agree exactly regardless of compiler. Only the transcendental
        // paths may differ.
        for (i, (got, want)) in fresh.iter().zip(&trace.steps).enumerate() {
            assert_eq!(
                got.rng.to_bits(),
                want.rng.to_bits(),
                "step {i}: the generator is exact f64 arithmetic and must match \
                 across compilers even when the physics does not"
            );
        }
        assert!(
            worst.deviation <= CROSS_COMPILER_TOLERANCE,
            "gfortran {} disagrees with the recording compiler {} by {:.3e} at \
             {}, beyond the {:.0e} allowed for a compiler difference. Either a \
             real regression, or transcendental differences larger than \
             expected; find out which before touching the bound.",
            build_info::GFORTRAN_VERSION,
            trace.gfortran,
            worst.deviation,
            worst.what,
            CROSS_COMPILER_TOLERANCE
        );
    }
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
