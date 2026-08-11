//! Tier 2 for the state unpacking: `tepsim_core::vessels::unpack` against
//! `teprob.f:417-472`, field by field, over all three sampling pools.
//!
//! The first physics comparison in the project. It calls `TESUB2` four times
//! and `TESUB4` three times, all of which Tier 1 already proved bit-exact over
//! ten million cases, so anything that differs here differs in the unpacking
//! and nowhere else.

#![cfg(feature = "oracle")]

use tepsim_core::{State, vessels};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::Comparison;
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, Snapshot, adversarial, compare_field};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2".
const TIER2_TOLERANCE: f64 = 1e-12;

/// Every field this item computes, paired with where to read it from each side.
struct Field {
    name: &'static str,
    ours: fn(&vessels::Unpacked) -> Vec<f64>,
    theirs: fn(&Snapshot) -> Vec<f64>,
}

fn fields() -> Vec<Field> {
    vec![
        Field {
            name: "UCLR",
            ours: |u| u.reactor.moles.as_array().to_vec(),
            theirs: |s| s.common.uclr.to_vec(),
        },
        Field {
            name: "UCLS",
            ours: |u| u.separator.moles.as_array().to_vec(),
            theirs: |s| s.common.ucls.to_vec(),
        },
        Field {
            name: "UCLC",
            ours: |u| u.stripper.moles.as_array().to_vec(),
            theirs: |s| s.common.uclc.to_vec(),
        },
        Field {
            name: "UCVV",
            ours: |u| u.mixing.moles.as_array().to_vec(),
            theirs: |s| s.common.ucvv.to_vec(),
        },
        Field {
            name: "UTLR/UTLS/UTLC/UTVV",
            ours: |u| {
                vec![
                    u.reactor.total,
                    u.separator.total,
                    u.stripper.total,
                    u.mixing.total,
                ]
            },
            theirs: |s| vec![s.common.utlr, s.common.utls, s.common.utlc, s.common.utvv],
        },
        Field {
            name: "XLR",
            ours: |u| u.reactor.fractions.fractions().as_array().to_vec(),
            theirs: |s| s.common.xlr.to_vec(),
        },
        Field {
            name: "XLS",
            ours: |u| u.separator.fractions.fractions().as_array().to_vec(),
            theirs: |s| s.common.xls.to_vec(),
        },
        Field {
            name: "XLC",
            ours: |u| u.stripper.fractions.fractions().as_array().to_vec(),
            theirs: |s| s.common.xlc.to_vec(),
        },
        Field {
            name: "XVV",
            ours: |u| u.mixing.fractions.fractions().as_array().to_vec(),
            theirs: |s| s.common.xvv.to_vec(),
        },
        Field {
            name: "ESR/ESS/ESC/ESV",
            ours: |u| {
                vec![
                    u.reactor.specific_energy,
                    u.separator.specific_energy,
                    u.stripper.specific_energy,
                    u.mixing.specific_energy,
                ]
            },
            theirs: |s| vec![s.common.esr, s.common.ess, s.common.esc, s.common.esv],
        },
        Field {
            name: "TCR/TCS/TCC/TCV",
            ours: |u| {
                vec![
                    u.reactor.celsius,
                    u.separator.celsius,
                    u.stripper.celsius,
                    u.mixing.celsius,
                ]
            },
            theirs: |s| vec![s.common.tcr, s.common.tcs, s.common.tcc, s.common.tcv],
        },
        Field {
            name: "TKR/TKS/TKV",
            ours: |u| vec![u.reactor.kelvin(), u.separator.kelvin(), u.mixing.kelvin()],
            theirs: |s| vec![s.common.tkr, s.common.tks, s.common.tkv],
        },
        Field {
            name: "DLR/DLS/DLC",
            ours: |u| vec![u.reactor.density, u.separator.density, u.stripper.density],
            theirs: |s| vec![s.common.dlr, s.common.dls, s.common.dlc],
        },
        Field {
            name: "VLR/VLS/VLC",
            ours: |u| vec![u.reactor.volume, u.separator.volume, u.stripper.volume],
            theirs: |s| vec![s.common.vlr, s.common.vls, s.common.vlc],
        },
    ]
}

/// Force the scenario, unpack the same state in Rust *with the same warm-start
/// seeds*, and record every field.
///
/// The seeds come out of the scenario's own `COMMON` block, which is exactly
/// what the Fortran will use, so the two Newton solves start from the same
/// place. Without that the comparison would be measuring path dependence
/// rather than the port.
fn observe(
    oracle: &mut Oracle,
    scenario: &Scenario,
    pool: Pool,
    index: usize,
    comparisons: &mut [Comparison<Case>],
) -> bool {
    let snapshot = scenario.force(oracle);
    let seeds = vessels::TemperatureSeeds {
        reactor: scenario.common.tcr,
        separator: scenario.common.tcs,
        stripper: scenario.common.tcc,
        mixing: scenario.common.tcv,
    };
    let state = State::from_flat(&scenario.state);
    let Ok(unpacked) = vessels::unpack(&state, seeds) else {
        return false;
    };
    for (field, comparison) in fields().iter().zip(comparisons.iter_mut()) {
        compare_field(
            comparison,
            pool,
            index,
            &(field.ours)(&unpacked),
            &(field.theirs)(&snapshot),
        );
    }
    true
}

fn fresh() -> Vec<Comparison<Case>> {
    fields()
        .iter()
        .map(|f| Comparison::new(format!("unpack {}", f.name)))
        .collect()
}

#[test]
fn the_unpacking_matches_the_fortran_over_all_three_pools() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let mut comparisons = fresh();
    let mut skipped = 0;

    // Nominal trajectory.
    for index in 0..pools.trajectory.len() {
        if !observe(
            &mut oracle,
            &pools.nominal_case(index),
            Pool::Nominal,
            index,
            &mut comparisons,
        ) {
            skipped += 1;
        }
    }

    // Perturbations.
    let mut sampler = tepsim_oracle::tier1::Sampler::new(0x7E2_0017);
    for index in 0..2_000 {
        if !observe(
            &mut oracle,
            &pools.perturbed_case(index, &mut sampler),
            Pool::Perturbed,
            index,
            &mut comparisons,
        ) {
            skipped += 1;
        }
    }

    // The seventeen constructed boundaries.
    let (boundaries, missed) = adversarial::build(&mut oracle, &pools.nominal_case(0));
    assert!(missed.is_empty(), "the catalogue regressed: {missed:?}");
    for (index, boundary) in boundaries.iter().enumerate() {
        if !observe(
            &mut oracle,
            &boundary.scenario,
            Pool::Adversarial,
            index,
            &mut comparisons,
        ) {
            skipped += 1;
        }
    }

    for comparison in &comparisons {
        println!("{comparison}");
    }
    println!("states skipped because the port would not converge: {skipped}");

    assert_eq!(
        skipped, 0,
        "{skipped} states failed to converge in the port. The Fortran would \
         have returned its guess instead (delta D-001), so these are not \
         disagreements, but they are also not evidence, and they need naming."
    );
    for comparison in &comparisons {
        comparison.assert_within(TIER2_TOLERANCE);
    }
}

/// The tolerance is a backstop. The claim is bit equality, as everywhere else
/// in this port, and it is asserted directly so that a drift to 1e-15 fails
/// instead of passing quietly for the rest of the phase.
#[test]
fn the_unpacking_is_bit_identical_to_the_fortran() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);
    let mut comparisons = fresh();

    for index in 0..pools.trajectory.len() {
        observe(
            &mut oracle,
            &pools.nominal_case(index),
            Pool::Nominal,
            index,
            &mut comparisons,
        );
    }
    let mut sampler = tepsim_oracle::tier1::Sampler::new(1);
    for index in 0..500 {
        observe(
            &mut oracle,
            &pools.perturbed_case(index, &mut sampler),
            Pool::Perturbed,
            index,
            &mut comparisons,
        );
    }

    let mut wrong = Vec::new();
    for comparison in &comparisons {
        if comparison.max_ulp() != 0 {
            wrong.push(format!("{comparison}"));
        }
    }
    assert!(
        wrong.is_empty(),
        "{} field(s) are not bit-identical:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// The warm start is what makes bit equality possible, so check that dropping
/// it actually breaks the comparison. Otherwise the seeds are ceremony.
#[test]
fn solving_from_a_fixed_guess_instead_of_the_seed_breaks_bit_equality() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 100, DT);

    let mut seeded: Comparison<Case> = Comparison::new("TCV with the carried seed");
    let mut fixed: Comparison<Case> = Comparison::new("TCV from a fixed guess");

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let snapshot = scenario.force(&mut oracle);
        let state = State::from_flat(&scenario.state);

        let carried = vessels::TemperatureSeeds {
            reactor: scenario.common.tcr,
            separator: scenario.common.tcs,
            stripper: scenario.common.tcc,
            mixing: scenario.common.tcv,
        };
        if let Ok(u) = vessels::unpack(&state, carried) {
            seeded.observe(
                Case {
                    pool: Pool::Nominal,
                    index,
                    component: 4,
                },
                u.mixing.celsius,
                snapshot.common.tcv,
            );
        }
        // A perfectly reasonable fixed starting guess, which is what a port
        // that had not noticed the warm start would use.
        if let Ok(u) = vessels::unpack(&state, vessels::TemperatureSeeds::default()) {
            fixed.observe(
                Case {
                    pool: Pool::Nominal,
                    index,
                    component: 4,
                },
                u.mixing.celsius,
                snapshot.common.tcv,
            );
        }
    }

    println!("{seeded}");
    println!("{fixed}");
    assert_eq!(seeded.max_ulp(), 0, "the carried seed must be bit-exact");
    assert!(
        fixed.max_ulp() > 0,
        "solving from a fixed guess gave bit-identical answers too, so the \
         warm start does not matter after all and TemperatureSeeds could be \
         dropped. If this ever becomes true, simplify and record why."
    );
}
