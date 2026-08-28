//! The facade against the loop that is proven against the Fortran.
//!
//! B-0052. `tier5::run_port` is bit-identical to `teprob.f` driven by
//! `temain_mod.f` over 172,800 closed-loop steps, which makes it the most
//! thoroughly validated forty lines in this repository. The facade's only
//! claim is that it is the same loop, so that is what is tested: not "close",
//! not "within a tolerance", *the same bits*.
//!
//! If this passes, everything Tier 4 and Tier 5 proved about `run_port`
//! transfers to `tepsim::Simulation` unchanged, and by extension to the Python
//! and wasm bindings built on it.

#![cfg(feature = "oracle")]

use tepsim::{Scenario, Simulation};
use tepsim_core::TemperatureSeeds;
use tepsim_oracle::Oracle;
use tepsim_oracle::tier5::{Scenario as TierScenario, run_port, seed, start};

/// Long enough to cross analyser samples, composition-loop fires and several
/// disturbance-walk segments.
const HOURS: usize = 3;

/// The constants the facade starts from are the Fortran's, bit for bit.
///
/// `TEINIT` calls `TEFUNC` once before returning, so the four Newton warm
/// starts a run begins from are not the nominal literals. `tepsim` cannot
/// depend on the oracle, so they are constants in `tepsim-core`; this is what
/// keeps them honest.
#[test]
fn the_facade_starts_from_the_fortrans_warm_start() {
    let mut oracle = Oracle::lock();
    let start = start(&mut oracle);
    let ours = TemperatureSeeds::after_initialisation();

    assert_eq!(
        ours.reactor.to_bits(),
        start.seeds.reactor.to_bits(),
        "TCR: {} against {}",
        ours.reactor,
        start.seeds.reactor
    );
    assert_eq!(ours.separator.to_bits(), start.seeds.separator.to_bits());
    assert_eq!(ours.stripper.to_bits(), start.seeds.stripper.to_bits());
    assert_eq!(ours.mixing.to_bits(), start.seeds.mixing.to_bits());

    // And they are *not* the nominal literals, which is the whole point.
    let naive = TemperatureSeeds::default();
    assert_ne!(
        ours.reactor.to_bits(),
        naive.reactor.to_bits(),
        "the converged warm start equals the literal, so TEINIT's own TEFUNC \
         call no longer moves it and this constant is unnecessary"
    );
    println!(
        "warm start: reactor {:.10} against the literal {:.10}",
        ours.reactor, naive.reactor
    );
}

/// The facade and the validated loop produce identical samples.
#[test]
fn the_facade_reproduces_the_validated_loop() {
    let mut oracle = Oracle::lock();
    let start = start(&mut oracle);

    for fault in [0_usize, 1, 6, 13, 20] {
        let tier = if fault == 0 {
            TierScenario::NOMINAL
        } else {
            TierScenario::fault(fault)
        };
        let theirs = run_port(&start, tier, seed(0), HOURS);

        let scenario = if fault == 0 {
            Scenario::baseline()
        } else {
            Scenario::fault(fault)
        }
        .with_hours(HOURS as f64)
        .with_seed(seed(0));
        let ours = Simulation::new(scenario).run();

        assert_eq!(
            ours.samples.len(),
            theirs.samples.len(),
            "{} sample count",
            tier.label()
        );
        for (index, (a, b)) in ours.samples.iter().zip(&theirs.samples).enumerate() {
            let row = a.row();
            for (channel, (x, y)) in row.iter().zip(b).enumerate() {
                assert_eq!(
                    x.to_bits(),
                    y.to_bits(),
                    "{}: sample {index}, channel {}: {x} against {y}",
                    tier.label(),
                    channel + 1
                );
            }
        }
        println!(
            "  {:<8} {} samples, identical to run_port bit for bit",
            tier.label(),
            ours.samples.len()
        );
    }
}

/// The priming measurements are never read, which is why the facade can start
/// from zeros without needing the oracle.
///
/// The fastest control loop is `MOD(I,3)`, so nothing fires until step 3, and
/// by then step 2 has produced real measurements. That reasoning is load
/// bearing, so it is checked rather than trusted: a run primed with `TEINIT`'s
/// measurements must be identical to one primed with zeros.
#[test]
fn the_priming_measurements_are_never_read() {
    let mut oracle = Oracle::lock();
    let start = start(&mut oracle);

    // `run_port` primes from `TEINIT`; the facade primes from zeros.
    let primed = run_port(&start, TierScenario::NOMINAL, seed(3), 1);
    let unprimed = Simulation::new(Scenario::baseline().with_hours(1.0).with_seed(seed(3))).run();

    assert!(
        start.measurements.iter().any(|m| *m != 0.0),
        "TEINIT left all measurements at zero, so priming from them and from \
         zeros is the same thing and this test proves nothing"
    );
    for (index, (a, b)) in unprimed.samples.iter().zip(&primed.samples).enumerate() {
        assert_eq!(a.row(), *b, "sample {index}");
    }
    println!(
        "priming from TEINIT's {} non-zero measurements changes nothing",
        start.measurements.iter().filter(|m| **m != 0.0).count()
    );
}

/// Stepping and running are the same computation.
#[test]
fn stepping_by_hand_matches_running() {
    let scenario = Scenario::baseline().with_hours(1.0).with_seed(seed(11));
    let batch = Simulation::new(scenario).run();

    let mut simulation = Simulation::new(scenario);
    let mut collected = Vec::new();
    while !simulation.is_halted() {
        if let Some(sample) = simulation.step() {
            collected.push(sample);
        }
    }

    assert_eq!(collected.len(), batch.samples.len());
    for (a, b) in collected.iter().zip(&batch.samples) {
        assert_eq!(a, b);
    }
    assert_eq!(simulation.outcome().unwrap_or(batch.outcome), batch.outcome);
    println!(
        "{} samples, stepped and batched identically",
        collected.len()
    );
}
