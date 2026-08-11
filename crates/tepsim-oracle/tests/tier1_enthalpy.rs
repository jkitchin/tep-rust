//! Tier 1 for `TESUB1`, `TESUB3` and `TESUB4`: the ported enthalpy, heat
//! capacity and liquid density against the original Fortran.
//!
//! The gate is a maximum relative error below 1e-13, from `PLAN.org`. What is
//! actually expected is zero: both routines are straight-line arithmetic over
//! constants that B-0006 already proved bit-identical, so if the association
//! and the literal precisions are right there is nothing left to differ. The
//! tolerance is there for the case where that reasoning is wrong, and the ULP
//! histogram is there to say by how much.
//!
//! Run with `--nocapture` to see the reports; they are what goes in `LOG.org`.

#![cfg(feature = "oracle")]

use std::time::Instant;

use tepsim_core::thermo::{EnergyBasis, enthalpy, heat_capacity, liquid_density};
use tepsim_oracle::{
    Oracle,
    tier1::{Case, Comparison, Sweep},
};

/// `PLAN.org`, "Tier 1": maximum relative error below 1e-13.
const TIER1_TOLERANCE: f64 = 1e-13;

#[test]
fn tesub1_matches_the_fortran_over_the_full_sweep() {
    let sweep = Sweep::from_env();
    println!("{}", sweep.provenance_note());
    let mut oracle = Oracle::lock();
    oracle.init();

    for basis in EnergyBasis::ALL {
        let started = Instant::now();
        let mut comparison: Comparison<Case> =
            Comparison::new(format!("TESUB1 ity={}", basis.ity()));
        for case in sweep.cases() {
            comparison.observe(
                case,
                enthalpy(&case.composition, case.celsius, basis),
                oracle.tesub1(&case.z(), case.celsius, basis.ity()),
            );
        }
        println!("{comparison}  [{:.1} s]", started.elapsed().as_secs_f64());
        assert_eq!(comparison.cases(), sweep.len() as u64);
        comparison.assert_within(TIER1_TOLERANCE);
    }
}

#[test]
fn tesub3_matches_the_fortran_over_the_full_sweep() {
    let sweep = Sweep::from_env();
    println!("{}", sweep.provenance_note());
    let mut oracle = Oracle::lock();
    oracle.init();

    for basis in EnergyBasis::ALL {
        let started = Instant::now();
        let mut comparison: Comparison<Case> =
            Comparison::new(format!("TESUB3 ity={}", basis.ity()));
        for case in sweep.cases() {
            comparison.observe(
                case,
                heat_capacity(&case.composition, case.celsius, basis),
                oracle.tesub3(&case.z(), case.celsius, basis.ity()),
            );
        }
        println!("{comparison}  [{:.1} s]", started.elapsed().as_secs_f64());
        assert_eq!(comparison.cases(), sweep.len() as u64);
        comparison.assert_within(TIER1_TOLERANCE);
    }
}

#[test]
fn tesub4_matches_the_fortran_over_the_full_sweep() {
    let sweep = Sweep::from_env();
    println!("{}", sweep.provenance_note());
    let mut oracle = Oracle::lock();
    oracle.init();

    let started = Instant::now();
    let mut comparison: Comparison<Case> = Comparison::new("TESUB4");
    for case in sweep.cases() {
        comparison.observe(
            case,
            liquid_density(&case.composition, case.celsius),
            oracle.tesub4(&case.z(), case.celsius),
        );
    }
    println!("{comparison}  [{:.1} s]", started.elapsed().as_secs_f64());
    assert_eq!(comparison.cases(), sweep.len() as u64);
    let (non_finite, _) = comparison.non_finite();
    assert_eq!(
        non_finite, 0,
        "the density correlation went singular inside the sweep range, which \
         it must not: the nearest pole is 208.57 C\n{comparison}"
    );
    comparison.assert_within(TIER1_TOLERANCE);
}

/// The tolerance above is a backstop, not the claim. The claim is bit equality,
/// which is much stronger and would be invisible under a 1e-13 gate: a port
/// that drifted to 1e-15 would still pass every day while something real had
/// changed. Assert what is actually true.
#[test]
fn both_routines_are_bit_identical_to_the_fortran() {
    let sweep = Sweep::from_env();
    println!("{}", sweep.provenance_note());
    let mut oracle = Oracle::lock();
    oracle.init();

    let mut differing = 0_u64;
    let mut first: Option<String> = None;
    for basis in EnergyBasis::ALL {
        for case in sweep.cases() {
            let pairs = [
                (
                    "TESUB1",
                    enthalpy(&case.composition, case.celsius, basis),
                    oracle.tesub1(&case.z(), case.celsius, basis.ity()),
                ),
                (
                    "TESUB3",
                    heat_capacity(&case.composition, case.celsius, basis),
                    oracle.tesub3(&case.z(), case.celsius, basis.ity()),
                ),
                // TESUB4 takes no ITY, so it is re-checked once per mode.
                // Redundant, and cheaper than a second pass over the sweep.
                (
                    "TESUB4",
                    liquid_density(&case.composition, case.celsius),
                    oracle.tesub4(&case.z(), case.celsius),
                ),
            ];
            for (name, ours, theirs) in pairs {
                if ours.to_bits() != theirs.to_bits() {
                    differing += 1;
                    first.get_or_insert_with(|| {
                        format!(
                            "{name} ity={} at {case}: Rust {ours:?}, Fortran \
                             {theirs:?}, z = {:?}",
                            basis.ity(),
                            case.z()
                        )
                    });
                }
            }
        }
    }

    assert_eq!(
        differing,
        0,
        "{differing} of {} evaluations differ; first: {}",
        9 * sweep.len(),
        first.unwrap_or_default()
    );
}

/// Prove the differential test has teeth, by breaking the one thing most likely
/// to be got wrong and checking the sweep notices.
///
/// The `273.15` at `teprob.f:1411` carries no `D` suffix, so it is single
/// precision. This evaluates the `ITY=2` correction both ways and confirms the
/// double-precision reading is not merely different but catastrophically so,
/// which is what justifies [`tepsim_core::thermo::ABSOLUTE_ZERO_OFFSET`] being
/// written the way it is.
#[test]
fn reading_the_offset_as_double_precision_would_fail_the_gate_by_orders() {
    use tepsim_core::thermo::{ABSOLUTE_ZERO_OFFSET, GAS_CONSTANT};

    let sweep = Sweep::SMOKE;
    let mut oracle = Oracle::lock();
    oracle.init();

    let mut sabotaged: Comparison<Case> = Comparison::new("TESUB1 ity=2, offset as f64");
    for case in sweep.cases() {
        let correct = enthalpy(&case.composition, case.celsius, EnergyBasis::VapourEnthalpy);
        // Exactly what the port would compute with `273.15_f64` in place of the
        // widened f32, and nothing else changed.
        let wrong = correct - GAS_CONSTANT * (case.celsius + 273.15_f64);
        sabotaged.observe(
            case,
            wrong,
            oracle.tesub1(
                &case.z(),
                case.celsius,
                EnergyBasis::VapourInternalEnergy.ity(),
            ),
        );
    }
    println!("{sabotaged}");

    assert!(
        sabotaged.max_relative_error() > 1e-6,
        "the single/double distinction has stopped mattering, which means \
         either the sweep no longer reaches the cancelling region or the \
         correction is no longer being applied\n{sabotaged}"
    );
    assert!(
        ABSOLUTE_ZERO_OFFSET < 273.15_f64,
        "the widened f32 rounds down, so this ordering is a cheap check that \
         the constant has not been silently corrected"
    );
}
