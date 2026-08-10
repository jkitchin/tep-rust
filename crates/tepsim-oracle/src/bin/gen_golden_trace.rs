//! Regenerates the committed golden trace from the original Fortran.
//!
//! ```text
//! cargo run -p tepsim-oracle --features oracle --bin gen-golden-trace
//! ```
//!
//! Running this is a **deliberate re-baseline**. Every validation number in
//! `LOG.org` was measured against the current trace, so replacing it without
//! also re-recording those numbers destroys the ability to tell a regression
//! from a toolchain change. It exists so that re-baselining is one obvious
//! command rather than an improvised script.

use std::path::PathBuf;

use tepsim_oracle::golden::{self, Step, Trace};
use tepsim_oracle::{Oracle, build_info};

fn main() {
    let mut oracle = Oracle::lock();

    // TEINIT loads the nominal steady state. It also calls TEFUNC once itself,
    // at t=0, which is why the generator word is pinned afterwards rather than
    // before: we want the trace to start from a known seed regardless.
    let (mut time, mut yy) = oracle.init();
    oracle.set_disturbances(&[0; 20]);
    oracle.set_rng(golden::SEED);

    let mut steps = Vec::with_capacity(golden::STEPS);
    for _ in 0..golden::STEPS {
        let states = yy;
        // Same order as the original's INTGTR at temain_mod.f:1372: evaluate at
        // the current point, then advance time and state.
        let derivatives = oracle.derivatives(time, &yy);
        let measurements = oracle.measurements();
        let rng = oracle.rng();

        steps.push(Step {
            states,
            derivatives,
            measurements,
            rng,
        });

        time += golden::DT_HOURS;
        for (y, dy) in yy.iter_mut().zip(&derivatives) {
            *y += dy * golden::DT_HOURS;
        }
    }

    let trace = Trace {
        gfortran: build_info::GFORTRAN_VERSION.to_string(),
        fflags: build_info::FORTRAN_FLAGS.to_string(),
        seed: golden::SEED,
        dt_hours: golden::DT_HOURS,
        steps,
    };

    let path = workspace_root().join(golden::PATH);
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("creating golden/");
    }
    std::fs::write(&path, trace.to_text()).expect("writing the trace");

    println!("wrote {} steps to {}", trace.steps.len(), path.display());
    println!("  gfortran : {}", trace.gfortran);
    println!("  fflags   : {}", trace.fflags);
    println!(
        "  reactor T: {:.6} C at step 0, {:.6} C at step {}",
        trace.steps[0].measurements[8],
        trace.steps[trace.steps.len() - 1].measurements[8],
        trace.steps.len() - 1
    );
}

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}
