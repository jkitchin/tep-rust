//! Tier 2 for the reaction kinetics: `tepsim_core::kinetics` against
//! `teprob.f:503-528`, over all three sampling pools.
//!
//! Run twice, default and `--features oracle,libm-system`, for the reason set
//! out in `tier2_equilibrium.rs`. This range is where that matters most: it
//! adds `pow` to `exp`, and `RR(1)` is a six-term product chain, which is
//! precisely the shape where a reassociation is invisible to 1e-12.
//!
//! # The drift factors have to come from `TESUB8`, not from `COMMON`
//!
//! `R1F` and `R2F` reach `teprob.f:503-504` as the IDV(13) kinetics-drift
//! multipliers and are **reassigned in place** at `508-509` to hold the
//! fractional pressure powers. So after `TEFUNC` returns, `COMMON`'s `R1F` is
//! not the number the rate law used, and a harness that read it there would be
//! comparing against the wrong input while looking perfectly reasonable.
//!
//! They are fetched from `TESUB8(7, t)` and `TESUB8(8, t)` instead, called
//! after the evaluation. That is sound because nothing between `teprob.f:406`
//! and the `RETURN` writes any of `ADIST`, `BDIST`, `CDIST`, `DDIST` or
//! `TLAST`, so the walk state on return is the one line 415 read, and `TESUB8`
//! is a pure Horner evaluation of it (`teprob.f:1585-1586`).

#![cfg(feature = "oracle")]

use tepsim_core::{Component, State, equilibrium, kinetics, math, vessels};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::Comparison;
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, Snapshot, adversarial, compare_field};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2".
const TIER2_TOLERANCE: f64 = 1e-12;

struct Field {
    name: &'static str,
    ours: fn(&kinetics::Kinetics) -> Vec<f64>,
    theirs: fn(&Snapshot) -> Vec<f64>,
}

fn fields() -> Vec<Field> {
    vec![
        Field {
            name: "RR",
            ours: |k| k.rates.to_vec(),
            theirs: |s| s.common.rr.to_vec(),
        },
        Field {
            name: "CRXR",
            ours: |k| k.production.as_array().to_vec(),
            theirs: |s| s.common.crxr.to_vec(),
        },
        Field {
            name: "RH",
            ours: |k| vec![k.heat],
            theirs: |s| vec![s.common.rh],
        },
    ]
}

fn observe(
    oracle: &mut Oracle,
    scenario: &Scenario,
    pool: Pool,
    index: usize,
    comparisons: &mut [Comparison<Case>],
) -> bool {
    let snapshot = scenario.force(oracle);
    // See the module documentation: after the call, not from `COMMON`.
    let drift = kinetics::ReactionDrift {
        first: oracle.tesub8(7, scenario.time),
        second: oracle.tesub8(8, scenario.time),
    };
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
    let eq = equilibrium::equilibrium(&unpacked);
    let solved = kinetics::kinetics(&eq.reactor, unpacked.reactor.kelvin(), drift);
    for (field, comparison) in fields().iter().zip(comparisons.iter_mut()) {
        compare_field(
            comparison,
            pool,
            index,
            &(field.ours)(&solved),
            &(field.theirs)(&snapshot),
        );
    }
    true
}

fn fresh() -> Vec<Comparison<Case>> {
    fields()
        .iter()
        .map(|f| Comparison::new(format!("kinetics {}", f.name)))
        .collect()
}

fn sweep(
    oracle: &mut Oracle,
    steps: usize,
    perturbations: usize,
    seed: u64,
) -> (Vec<Comparison<Case>>, usize) {
    let pools = Pools::collect(oracle, steps, DT);
    let mut comparisons = fresh();
    let mut skipped = 0;

    for index in 0..pools.trajectory.len() {
        if !observe(
            oracle,
            &pools.nominal_case(index),
            Pool::Nominal,
            index,
            &mut comparisons,
        ) {
            skipped += 1;
        }
    }

    let mut sampler = tepsim_oracle::tier1::Sampler::new(seed);
    for index in 0..perturbations {
        if !observe(
            oracle,
            &pools.perturbed_case(index, &mut sampler),
            Pool::Perturbed,
            index,
            &mut comparisons,
        ) {
            skipped += 1;
        }
    }

    let (boundaries, missed) = adversarial::build(oracle, &pools.nominal_case(0));
    assert!(missed.is_empty(), "the catalogue regressed: {missed:?}");
    for (index, boundary) in boundaries.iter().enumerate() {
        if !observe(
            oracle,
            &boundary.scenario,
            Pool::Adversarial,
            index,
            &mut comparisons,
        ) {
            skipped += 1;
        }
    }

    (comparisons, skipped)
}

#[test]
fn the_kinetics_match_the_fortran_over_all_three_pools() {
    let mut oracle = Oracle::lock();
    let (comparisons, skipped) = sweep(&mut oracle, 400, 2_000, 0x7E2_0019);

    println!(
        "exp and pow come from the {} libm",
        if math::USES_SYSTEM_LIBM {
            "platform"
        } else {
            "vendored"
        }
    );
    for comparison in &comparisons {
        println!("{comparison}");
    }
    println!("states skipped because the port would not converge: {skipped}");

    assert_eq!(
        skipped, 0,
        "{skipped} states failed to converge in the port"
    );
    for comparison in &comparisons {
        comparison.assert_within(TIER2_TOLERANCE);
    }
}

/// With both transcendentals taken out of the comparison, the algebra must be
/// bit-identical. See `tier2_equilibrium.rs` for what a failure here means and
/// what it does not.
#[test]
#[cfg(feature = "libm-system")]
fn the_algebra_is_bit_identical_once_exp_and_pow_agree() {
    let mut oracle = Oracle::lock();
    let (comparisons, _) = sweep(&mut oracle, 200, 500, 2);

    let mut wrong = Vec::new();
    for comparison in &comparisons {
        println!("{comparison}");
        if comparison.max_ulp() != 0 {
            wrong.push(format!("{comparison}"));
        }
    }
    assert!(
        wrong.is_empty(),
        "{} field(s) are not bit-identical under the platform libm.\n\
         Either the algebra in `tepsim_core::kinetics` does not associate the \
         way `teprob.f:503-528` does, which is a port bug; or gfortran no \
         longer resolves `DEXP`/`**` to the same code as the platform `exp` \
         and `pow`, in which case this configuration no longer removes the \
         transcendentals and the claim it supports has to be restated.\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// `CRXR(2)` is never assigned by the original, and is read anyway. Delta
/// D-003.
///
/// The oracle must report exactly zero for it, on every state, or the whole
/// reading is wrong and B's balance at `teprob.f:763` is picking up garbage.
#[test]
fn the_inert_has_no_net_production_in_either_implementation() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);
    let mut comparison: Comparison<Case> = Comparison::new("CRXR(2), never assigned");

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let snapshot = scenario.force(&mut oracle);
        assert_eq!(
            snapshot.common.crxr[Component::B.index()].to_bits(),
            0.0_f64.to_bits(),
            "the Fortran reported a non-zero CRXR(2) at nominal#{index}. \
             Nothing in teprob.f writes that slot, so this would mean it is \
             picking up whatever was left in COMMON, and delta D-003 is not \
             the benign reading it is recorded as."
        );
        comparison.observe(
            Case {
                pool: Pool::Nominal,
                index,
                component: Component::B.fortran_index(),
            },
            0.0,
            snapshot.common.crxr[Component::B.index()],
        );
    }
    println!("{comparison}");
    assert_eq!(comparison.max_ulp(), 0);
}

/// Measure the price of the vendored `pow`, over the pressure range the
/// fractional orders at `teprob.f:508-509` actually see.
///
/// The companion to `the_vendored_and_platform_exp_differ_only_by_rounding` in
/// `tier2_equilibrium.rs`, and the reason it is worth having separately: `pow`
/// is the harder function to round, so its disagreement is not predictable
/// from `exp`'s.
#[test]
fn the_vendored_and_platform_pow_differ_only_by_rounding() {
    use tepsim_core::constants::single;

    let orders = [single(1.1544), single(0.3735)];
    let (mut cases, mut differing, mut worst) = (0u64, 0u64, 0i64);

    // Partial pressures of A and C run to a few hundred mmHg on the nominal
    // trajectory; sweep two decades either side of that.
    for step in 1..=20_000u32 {
        let base = f64::from(step) * 0.25;
        for order in orders {
            cases += 1;
            let delta =
                (math::pow(base, order).to_bits() as i64) - (base.powf(order).to_bits() as i64);
            if delta != 0 {
                differing += 1;
            }
            if delta.abs() > worst.abs() {
                worst = delta;
            }
        }
    }

    println!(
        "pow over [0.25, 5000] at orders 1.1544 and 0.3735: {cases} cases, \
         {differing} differ ({:.3}%), worst {worst} ulp",
        100.0 * differing as f64 / cases as f64
    );
    if math::USES_SYSTEM_LIBM {
        assert_eq!(differing, 0, "this build *is* the platform libm");
    } else {
        assert!(
            worst.abs() <= 1,
            "the vendored pow is {worst} ulp from the platform, not at most 1. \
             Every Phase 2 error budget is written against a one-ulp difference."
        );
    }
}
