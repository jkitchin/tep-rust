//! Tier 4: trajectory equivalence. **Diagnostic, not a gate.**
//!
//! `PLAN.org` is explicit that long-horizon divergence is expected and is not
//! a bug: the plant amplifies one-ULP differences in `exp` and `pow`, so two
//! implementations that agree to the last bit on every single step will still
//! part company given enough steps.
//!
//! What Tier 4 owes is therefore not a pass mark but two numbers and an
//! explanation:
//!
//! - **How long agreement lasts.** Error must stay below the corresponding
//!   measurement noise `XNS(i)` for at least the first several hours. A
//!   difference smaller than the instrument's own noise is not observable.
//! - **Why it ends.** Switching between the vendored `libm` and the platform
//!   one must *move* the divergence point. That converts an unexplained
//!   mismatch into a characterised property of floating-point arithmetic, and
//!   it is the whole reason the `libm-system` feature exists.
//!
//! Under `libm-system` the port and the Fortran call the same `exp` and `pow`,
//! so there is nothing left to diverge and the trajectories stay together
//! indefinitely. That is the demonstration.

#![cfg(feature = "oracle")]

use tepsim_core::constants::{MEASUREMENT_NOISE, NOMINAL_STATE};
use tepsim_core::{Inputs, Plant, SimTime, State, math};
use tepsim_oracle::Oracle;

const DT: f64 = 1.0 / 3600.0;
/// Steps per simulated hour at a one-second step.
const PER_HOUR: usize = 3_600;

/// How far the two trajectories ran together, measured against `XNS`.
struct Divergence {
    /// The last hour at which every measurement was still within its own
    /// noise magnitude.
    hours_within_noise: f64,
    /// The measurement that gave out first, one-based.
    first_out: usize,
    /// Worst error over noise at the end of the run.
    worst_at_end: f64,
}

/// Run both from the nominal state and report where they part.
fn run(oracle: &mut Oracle, hours: usize, fault: usize) -> Divergence {
    let (_, mut fortran) = oracle.init_cold();
    oracle.set_disturbances(&core::array::from_fn(|i| i32::from(i + 1 == fault)));
    oracle.set_rng(tepsim_oracle::golden::SEED);

    let mut plant = Plant::new();
    plant.set_rng(tepsim_oracle::golden::SEED);
    // `TEINIT` calls `TEFUNC` once itself at t=0 (`teprob.f:1369`), so by the
    // time it returns the four Newton warm-start temperatures have already
    // converged away from the nominal literals `TemperatureSeeds::default()`
    // carries. B-0017 measured that seeding them differently moves up to 21 of
    // the 50 derivatives in the last bits, so a trajectory started from the
    // defaults is a *different* trajectory, not a rounding of the same one.
    let after_init = oracle.teproc();
    plant.set_seeds(tepsim_core::TemperatureSeeds {
        reactor: after_init.tcr,
        separator: after_init.tcs,
        stripper: after_init.tcc,
        mixing: after_init.tcv,
    });
    let mut state = State::from_flat(&NOMINAL_STATE);
    let inputs = Inputs {
        manipulated: core::array::from_fn(|i| NOMINAL_STATE[38 + i]),
        disturbances: core::array::from_fn(|i| f64::from(u8::from(i + 1 == fault))),
    };

    let mut hours_within_noise = 0.0;
    let mut first_out = 0;
    let mut worst_at_end = 0.0_f64;
    let mut t = 0.0;

    for step in 0..hours * PER_HOUR {
        let time = SimTime(t);
        plant.advance_discrete(time, &inputs);
        let Ok((derivative, signals)) = plant.derivatives(time, &state, &inputs) else {
            break;
        };
        let ours = plant.sample_measurements(time, &signals);
        let yp = oracle.derivatives(t, &fortran);
        let theirs = oracle.measurements();

        // Error relative to each instrument's own noise. A difference smaller
        // than `XNS(i)` is not observable on that instrument.
        let mut worst = 0.0_f64;
        let mut culprit = 0;
        for (index, (a, b)) in ours.as_array().iter().zip(theirs).enumerate() {
            let noise = MEASUREMENT_NOISE[index];
            if noise == 0.0 {
                continue;
            }
            let ratio = (a - b).abs() / noise;
            if ratio > worst {
                worst = ratio;
                culprit = index + 1;
            }
        }
        if worst < 1.0 {
            hours_within_noise = t;
        } else if first_out == 0 {
            first_out = culprit;
        }
        worst_at_end = worst;

        plant.step_seeds(&state).expect("converges");
        state = state.step(DT, &derivative);
        for (slot, rate) in fortran.iter_mut().zip(yp) {
            *slot += DT * rate;
        }
        t += DT;
        let _ = step;
    }

    Divergence {
        hours_within_noise,
        first_out,
        worst_at_end,
    }
}

/// The nominal trajectory, recorded rather than gated.
#[test]
fn tier4_nominal_trajectory() {
    let mut oracle = Oracle::lock();
    let hours = 8;
    let d = run(&mut oracle, hours, 0);

    println!(
        "Tier 4, nominal, {} libm, {hours} h:",
        if math::USES_SYSTEM_LIBM {
            "platform"
        } else {
            "vendored"
        }
    );
    println!("  within XNS for : {:.3} h", d.hours_within_noise);
    println!(
        "  first out      : {}",
        if d.first_out == 0 {
            "none".to_string()
        } else {
            format!("XMEAS({})", d.first_out)
        }
    );
    println!("  worst at end   : {:.3e} of XNS", d.worst_at_end);

    if math::USES_SYSTEM_LIBM {
        // Same `exp`, same `pow`, nothing left to diverge.
        assert!(
            d.first_out == 0,
            "on the platform libm the two call the same transcendentals, so \
             XMEAS({}) diverging means something other than rounding is \
             wrong",
            d.first_out
        );
        assert!(
            d.worst_at_end < 1e-6,
            "worst error was {:.3e} of XNS after {hours} h with identical \
             transcendentals",
            d.worst_at_end
        );
    } else {
        // `PLAN.org`: "error must stay below the corresponding measurement
        // noise standard deviation for at least the first several hours".
        assert!(
            d.hours_within_noise >= 1.0,
            "the trajectories parted within {:.3} h, which is less than the \
             'several hours' PLAN.org asks for. That is short enough to be a \
             bug rather than amplification.",
            d.hours_within_noise
        );
    }
}

/// The divergence is caused by the transcendentals, demonstrated rather than
/// asserted.
///
/// This is the test `PLAN.org` describes when it says the onset "must be
/// *explained*, by showing that switching between the vendored libm and the
/// system libm moves the divergence point". It runs under the default build
/// and reports; the `libm-system` run of the same file is the other half, and
/// `xtask validate --tiers 4` runs both.
#[test]
fn the_divergence_is_transcendental_rounding_and_nothing_else() {
    let mut oracle = Oracle::lock();
    let d = run(&mut oracle, 4, 0);

    if math::USES_SYSTEM_LIBM {
        println!(
            "platform libm: still together after 4 h, worst {:.3e} of XNS",
            d.worst_at_end
        );
        assert!(
            d.worst_at_end < 1e-6,
            "the trajectories parted even with identical transcendentals, so \
             the explanation is wrong and something else differs: {:.3e}",
            d.worst_at_end
        );
    } else {
        println!(
            "vendored libm: within XNS for {:.3} h, worst {:.3e} of XNS at 4 h",
            d.hours_within_noise, d.worst_at_end
        );
        assert!(
            d.worst_at_end > 1e-6,
            "the vendored libm produced no measurable divergence in 4 h, so \
             this test cannot demonstrate the contrast and the claim is \
             untested"
        );
    }
}

/// A handful of fault scenarios, recorded.
///
/// `PLAN.org` asks for all twenty-one eventually; this runs a representative
/// four short, because the point is the shape of the divergence and not its
/// exhaustive enumeration. `xtask validate --tiers 4` runs the full set.
#[test]
#[ignore = "long; run through `cargo xtask validate --tiers 4`"]
fn tier4_fault_scenarios() {
    let mut oracle = Oracle::lock();
    println!("Tier 4, fault scenarios, 4 h each:");
    for fault in 0..=20 {
        let d = run(&mut oracle, 4, fault);
        let name = if fault == 0 {
            "nominal".to_string()
        } else {
            format!("IDV({fault})")
        };
        println!(
            "  {name:10} within XNS {:.3} h, first out {}, worst {:.2e}",
            d.hours_within_noise,
            if d.first_out == 0 {
                "none".to_string()
            } else {
                format!("XMEAS({})", d.first_out)
            },
            d.worst_at_end
        );
    }
}
