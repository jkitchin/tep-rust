//! Drives the Tier 1 harness against the real Fortran, end to end.
//!
//! No routine has been ported yet, so there is nothing here to *validate*. What
//! there is to establish is that the instrument works before it is trusted with
//! a result: that a case from the sweep survives the trip through the FFI
//! unchanged, that the report reads zero when the two sides genuinely agree,
//! and, the part that actually matters, that it does not read zero when they do
//! not.
//!
//! The discrepancy is injected by perturbing the *input* rather than by writing
//! a stand-in implementation. A stand-in would be a port of `TESUB1` in all but
//! name, and the next session would be tempted to copy it instead of reading
//! `teprob.f`, which is exactly the failure mode `CLAUDE.md` exists to prevent.

#![cfg(feature = "oracle")]

use std::time::Instant;

use tepsim_oracle::{
    Oracle,
    tier1::{Case, Comparison, EXACT_BUCKETS, Sweep},
};

/// Every `(routine, ity)` pair Tier 1 covers, evaluated the way the harness
/// will drive it in B-0009 and after.
fn evaluate(oracle: &mut Oracle, routine: &str, case: &Case) -> f64 {
    match routine {
        "TESUB1 ity=0" => oracle.tesub1(&case.z(), case.celsius, 0),
        "TESUB1 ity=1" => oracle.tesub1(&case.z(), case.celsius, 1),
        "TESUB1 ity=2" => oracle.tesub1(&case.z(), case.celsius, 2),
        "TESUB3 ity=0" => oracle.tesub3(&case.z(), case.celsius, 0),
        "TESUB3 ity=1" => oracle.tesub3(&case.z(), case.celsius, 1),
        "TESUB3 ity=2" => oracle.tesub3(&case.z(), case.celsius, 2),
        "TESUB4" => oracle.tesub4(&case.z(), case.celsius),
        other => unreachable!("no such routine: {other}"),
    }
}

const ROUTINES: [&str; 7] = [
    "TESUB1 ity=0",
    "TESUB1 ity=1",
    "TESUB1 ity=2",
    "TESUB3 ity=0",
    "TESUB3 ity=1",
    "TESUB3 ity=2",
    "TESUB4",
];

/// Run the whole sweep through every routine, comparing the Fortran with
/// itself. The report must be all zeros, and every case must produce a finite
/// number: a sweep that quietly generated infinities would make the eventual
/// Tier 1 result meaningless without failing anything.
#[test]
fn the_sweep_drives_every_routine_and_the_fortran_agrees_with_itself() {
    let sweep = Sweep::SMOKE;
    let mut oracle = Oracle::lock();
    oracle.init();

    let started = Instant::now();
    let mut evaluations = 0_u64;

    for routine in ROUTINES {
        let mut comparison: Comparison<Case> = Comparison::new(routine);
        for case in sweep.cases() {
            let value = evaluate(&mut oracle, routine, &case);
            let again = evaluate(&mut oracle, routine, &case);
            comparison.observe(case, again, value);
            evaluations += 2;
        }
        println!("{comparison}");

        assert_eq!(comparison.cases(), sweep.len() as u64);
        assert_eq!(
            comparison.non_finite(),
            (0, 0),
            "{routine} produced a non-finite value somewhere in the sweep"
        );
        assert_eq!(
            comparison.max_ulp(),
            0,
            "{routine} is not a pure function of its arguments"
        );
        comparison.assert_within(0.0);
    }

    // Not an assertion, a datapoint: `Sweep::FULL` is about 10^7 cases, and
    // B-0009 needs to know whether that fits in a test.
    let elapsed = started.elapsed();
    println!(
        "sweep: {} cases x {} routines, {evaluations} evaluations in {:.2} s \
         ({:.0} evaluations/s, debug build)",
        sweep.len(),
        ROUTINES.len(),
        elapsed.as_secs_f64(),
        evaluations as f64 / elapsed.as_secs_f64(),
    );
}

/// The teeth. Nudge the temperature by one ULP and the report must stop reading
/// zero, for every routine.
///
/// One ULP in is the smallest perturbation that exists, so a harness that sees
/// this will see anything a real porting mistake produces.
#[test]
fn a_one_ulp_change_in_the_input_shows_up_in_the_report() {
    let sweep = Sweep::SMOKE;
    let mut oracle = Oracle::lock();
    oracle.init();

    for routine in ROUTINES {
        let mut comparison: Comparison<Case> = Comparison::new(routine);
        for case in sweep.cases() {
            let reference = evaluate(&mut oracle, routine, &case);
            let nudged = Case {
                celsius: case.celsius.next_up(),
                ..case
            };
            comparison.observe(case, evaluate(&mut oracle, routine, &nudged), reference);
        }
        println!("{comparison}");

        assert!(
            comparison.max_ulp() > 0,
            "{routine}: a one-ULP change in the input left the report reading \
             zero, so the harness is blind\n{comparison}"
        );
        assert!(
            comparison.max_relative_error() > 0.0,
            "{routine}: relative error stayed at zero\n{comparison}"
        );
    }
}

/// `TESUB1` with `ITY=2` is ill-conditioned, and the harness says so the first
/// time it is pointed at it.
///
/// `teprob.f:1410-1411` subtracts `R*(T+273.15)` from the vapour enthalpy, with
/// `R = 3.57696e-6`. Over much of the composition space that correction is
/// nearly the whole of `H`: at an equimolar A/B/C/D mixture at 21.875 degrees,
/// `H` is 1.0567e-3 and the correction is 1.0553e-3, so the subtraction throws
/// away about ten bits. A one-ULP change in the temperature comes out the other
/// side as a thousand.
///
/// This does not make a correct port fail. Identical arithmetic in an identical
/// order still gives identical bits, however badly conditioned. It does mean
/// `ITY=2` is where any difference at all shows up magnified by three orders of
/// magnitude, so it is the mode to suspect first when Tier 1 moves, and the one
/// with the least margin under the 1e-13 gate.
#[test]
fn ity_2_amplifies_a_one_ulp_input_change_by_three_orders_of_magnitude() {
    let sweep = Sweep::SMOKE;
    let mut oracle = Oracle::lock();
    oracle.init();

    let mut amplification = [0_u64; 3];
    for ity in 0..3_i32 {
        let mut comparison: Comparison<Case> = Comparison::new(format!("TESUB1 ity={ity}"));
        for case in sweep.cases() {
            let reference = oracle.tesub1(&case.z(), case.celsius, ity);
            let nudged = oracle.tesub1(&case.z(), case.celsius.next_up(), ity);
            comparison.observe(case, nudged, reference);
        }
        println!("{comparison}");
        amplification[ity as usize] = comparison.max_ulp();
    }

    for (ity, worst) in amplification.iter().enumerate().take(2) {
        assert!(
            *worst < EXACT_BUCKETS as u64,
            "ity={ity} is well conditioned and should stay inside the counted \
             buckets, got {worst} ULP"
        );
    }
    assert!(
        amplification[2] > 100 * amplification[1],
        "ity=2 lost its cancellation, or ity=1 acquired one: {} vs {} ULP. \
         Either way the conditioning of the port has changed and the Tier 1 \
         numbers need re-reading, not re-baselining",
        amplification[2],
        amplification[1]
    );
}

/// The `273.15` at `teprob.f:1411` is a single-precision literal.
///
/// It carries no `D` suffix, so gfortran stores 273.14999389648438 and widens
/// that. Writing `273.15_f64` in the port would make the correction wrong by
/// 2.1e-8 relative, which the cancellation above then multiplies into roughly
/// 1.6e-5 in the returned enthalpy: eight orders of magnitude past the Tier 1
/// gate. This is the same hazard B-0006 found in the constants table, in the
/// one place a constants test cannot see it, because this literal is not in
/// `COMMON/CONST/` at all.
///
/// The equality below is exact rather than approximate, because it can be: at a
/// cancelling case the subtraction is exact by Sterbenz's lemma, so `H(ity=1) -
/// H(ity=2)` recovers the correction term bit for bit.
#[test]
fn the_ity_2_correction_uses_a_single_precision_273_15() {
    use tepsim_core::constants::single;

    let mut oracle = Oracle::lock();
    oracle.init();

    // Equimolar A/B/C/D at 21.875 C, where the correction is 99.87% of H, so
    // the difference of the two modes is exact.
    let z = [0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0];
    let celsius = 21.875;
    let correction = oracle.tesub1(&z, celsius, 1) - oracle.tesub1(&z, celsius, 2);

    let r = 3.57696_f64 / 1.0e6;
    let as_single = r * (celsius + single(273.15));
    let as_double = r * (celsius + 273.15_f64);

    assert_eq!(
        correction.to_bits(),
        as_single.to_bits(),
        "expected the single-precision literal: Fortran {correction:?}, \
         single {as_single:?}, double {as_double:?}"
    );
    let discrepancy = (as_double - correction).abs() / correction.abs();
    assert!(
        discrepancy > 1e-8,
        "and the double-precision reading must be detectably wrong, not a \
         rounding away; it differs by {discrepancy:e}"
    );
}

/// A composition reaches the Fortran in the order it was generated in.
///
/// Every Tier 1 number depends on this and it is one array-layout mistake away
/// from being false, but checking it must not mean writing a second
/// implementation to compare against. It does not have to: `TESUB1` is exactly
/// linear in the composition for `ITY` 0, so the Fortran can supply its own
/// reference. Evaluate the eight pure species, then predict a mixture from
/// them. If the array arrived reversed or shifted, the prediction would be
/// built from the wrong species and would miss by a wide margin.
///
/// The tolerance is not a Tier 1 gate. It is the cost of summing eight positive
/// terms in a different association than the Fortran uses, at most a few ULP.
#[test]
fn a_composition_arrives_in_the_order_it_was_generated_in() {
    let mut oracle = Oracle::lock();
    oracle.init();

    let mut comparison: Comparison<Case> = Comparison::new("TESUB1 ity=0 linearity");
    for case in Sweep::SMOKE.cases().take(2_000) {
        let mut predicted = 0.0_f64;
        for (index, fraction) in case.z().into_iter().enumerate() {
            let mut pure = [0.0_f64; 8];
            pure[index] = 1.0;
            predicted = fraction.mul_add(oracle.tesub1(&pure, case.celsius, 0), predicted);
        }
        comparison.observe(case, predicted, oracle.tesub1(&case.z(), case.celsius, 0));
    }
    println!("{comparison}");
    comparison.assert_within(1e-14);
}
