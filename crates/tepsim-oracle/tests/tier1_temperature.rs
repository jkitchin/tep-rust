//! Tier 1 for `TESUB2`, and the measurement behind delta D-001.
//!
//! Two questions, answered by the same sweep. Does the ported Newton iteration
//! land on the same temperature as the Fortran's, bit for bit? And how often
//! does the Fortran's silent non-convergence path actually fire, which is what
//! decides whether returning a `Result` is a free improvement or a behaviour
//! change?
//!
//! # How the delta is measured without reimplementing the quirk
//!
//! On failure the original restores the caller's guess (`teprob.f:1440`) and
//! returns as though it had converged. The port returns an error. So the two
//! are reconciled case by case: where the port returns `Ok`, the value must
//! match the Fortran's exactly; where it returns `Err`, the Fortran must have
//! returned *precisely the guess it was given*, which is the quirk's signature
//! and is checkable without a second implementation of it.

#![cfg(feature = "oracle")]

use std::time::Instant;

use tepsim_core::thermo::{EnergyBasis, temperature_from_enthalpy};
use tepsim_oracle::{
    Oracle,
    tier1::{Case, Comparison, Sweep},
};

/// `PLAN.org`, "Tier 1".
const TIER1_TOLERANCE: f64 = 1e-13;

/// Where the iteration is started from.
///
/// The plant always has a warm start: every `TESUB2` call site passes the
/// previous step's temperature (`teprob.f:460-465`), which after one second of
/// simulated time is a very good guess. That is the case the model actually
/// exercises, so it is measured; but a warm start alone would never test the
/// iteration, so the cold start is measured too.
#[derive(Clone, Copy)]
struct Start {
    name: &'static str,
    mirror: bool,
}

const STARTS: [Start; 2] = [
    Start {
        name: "warm",
        mirror: false,
    },
    Start {
        name: "cold",
        mirror: true,
    },
];

impl Start {
    fn guess(self, case: &Case, low: f64, high: f64) -> f64 {
        if self.mirror {
            // The far end of the range: as bad a guess as the domain allows.
            low + high - case.celsius
        } else {
            case.celsius
        }
    }
}

/// The Tier 1 gate for `TESUB2`, and the D-001 count in the same pass.
#[test]
fn tesub2_matches_the_fortran_and_the_silent_failure_never_fires() {
    let sweep = Sweep::from_env();
    println!("{}", sweep.provenance_note());
    let (low, high) = (sweep.range.low, sweep.range.high);
    let mut oracle = Oracle::lock();
    oracle.init();

    let mut total_abandoned = 0_u64;

    for basis in EnergyBasis::ALL {
        for start in STARTS {
            let began = Instant::now();
            let mut comparison: Comparison<Case> =
                Comparison::new(format!("TESUB2 ity={} start={}", basis.ity(), start.name));
            let mut abandoned = 0_u64;
            let mut worst_round_trip = 0.0_f64;
            let mut worst_round_trip_case = None;

            for case in sweep.cases() {
                let target = oracle.tesub1(&case.z(), case.celsius, basis.ity());
                let guess = start.guess(&case, low, high);
                let fortran = oracle.tesub2(&case.z(), guess, target, basis.ity());

                match temperature_from_enthalpy(&case.composition, guess, target, basis) {
                    Ok(solved) => {
                        comparison.observe(case, solved, fortran);
                        // Diagnostic, not the gate: how close the solve landed
                        // to the temperature the target enthalpy was built from.
                        let deviation = (solved - case.celsius).abs();
                        if deviation > worst_round_trip {
                            worst_round_trip = deviation;
                            worst_round_trip_case = Some(case);
                        }
                    }
                    Err(error) => {
                        abandoned += 1;
                        assert_eq!(
                            fortran.to_bits(),
                            guess.to_bits(),
                            "the port gave up at {case} but the Fortran did not \
                             restore the guess, so the two disagree about \
                             convergence rather than about what to do about it.\
                             \n  port: {error}\n  Fortran returned {fortran:?} \
                             for a guess of {guess:?}"
                        );
                    }
                }
            }

            println!(
                "{comparison}\n  abandoned      : {abandoned} (delta D-001)\
                 \n  round trip     : max |solved - true| = {worst_round_trip:e} C{}\
                 \n  elapsed        : {:.1} s",
                worst_round_trip_case
                    .map(|c| format!(" at {c}"))
                    .unwrap_or_default(),
                began.elapsed().as_secs_f64(),
            );

            assert_eq!(
                comparison.cases() + abandoned,
                sweep.len() as u64,
                "every case must be accounted for as either compared or abandoned"
            );
            comparison.assert_within(TIER1_TOLERANCE);

            // Newton is quadratic here, so once the step is under 1e-12 the
            // error is far smaller still. Three orders of margin, asserted so
            // that a future change which merely *passes* the step test without
            // actually solving the problem is caught.
            assert!(
                worst_round_trip < 1e-9,
                "the solve satisfied its step criterion without landing on the \
                 right temperature: worst deviation {worst_round_trip:e} C"
            );

            total_abandoned += abandoned;
        }
    }

    // The Class B claim: the fix is a no-op on the physical domain. If this
    // ever fires it is a genuine finding about the original, not a bug here,
    // and it belongs in book/src/deltas.md with the cases that produced it.
    assert_eq!(
        total_abandoned, 0,
        "delta D-001 fired {total_abandoned} times across the sweep, so \
         returning a Result is a behaviour change rather than a free \
         improvement. Record the cases in book/src/deltas.md before deciding."
    );
}

/// The `Result` has to be reachable through the real code path, or D-001 is a
/// story rather than a fix.
///
/// A target enthalpy no temperature produces sends Newton away from the domain.
/// The Fortran answers it by handing back the guess with no indication that
/// anything happened; the port answers it with an error. Both behaviours are
/// asserted here, side by side, because the delta register's whole claim is
/// about the difference between them.
#[test]
fn the_fortran_silently_returns_the_guess_where_the_port_reports_failure() {
    let mut oracle = Oracle::lock();
    oracle.init();

    let z = [0.125_f64; 8];
    let composition = tepsim_core::Composition::new(z);
    let guess = 120.4;
    let unreachable = -1.0e30;

    let fortran = oracle.tesub2(&z, guess, unreachable, 0);
    let port = temperature_from_enthalpy(
        &composition,
        guess,
        unreachable,
        EnergyBasis::LiquidEnthalpy,
    );

    assert_eq!(
        fortran.to_bits(),
        guess.to_bits(),
        "the Fortran should have restored the guess verbatim (teprob.f:1440), \
         got {fortran:?}"
    );
    let error = port.expect_err("the port must not report success here");
    println!("D-001 demonstrated:\n  Fortran: returned {fortran} silently\n  port:    {error}");
}
