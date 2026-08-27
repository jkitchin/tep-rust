//! Tier 2 for heat transfer: `tepsim_core::heat` against `teprob.f:663-678`,
//! over all three sampling pools.
//!
//! Run twice, default and `--features oracle,libm-system`; see
//! `tier2_equilibrium.rs`.
//!
//! # `UAS` and `UARLEV` are locals
//!
//! `COMMON/TEPROC/` carries `UAR`, `QUR`, `QUS` and `QUC`, but `UAS` and
//! `UARLEV` are declared among `TEFUNC`'s locals at `teprob.f:325` and `328`,
//! so neither can be read back directly.
//!
//! `UAR` covers `UARLEV`: it is `UARLEV` times a quadratic in `AGSP` times a
//! constant, and B-0021 validated `AGSP` at 0 ULP, so a wrong `UARLEV` cannot
//! hide in it. `UAS` is covered by `QUS` the same way `UAC` was covered by
//! `QUC` in B-0021, except that here no state needs excluding: `QUS` has no
//! conditional.
//!
//! # The branch-coverage test is the point of this file
//!
//! B-0022 found that a boundary state does not exercise the branch it bounds,
//! because every comparison in `TEFUNC` is strict. Both `UARLEV` breakpoints
//! have boundary states from B-0016, and by that rule *both take the ramp*.
//! The two flat branches would be implemented and unvalidated if nothing said
//! so, which is exactly the failure B-0022 caught one item too late.

#![cfg(feature = "oracle")]

use std::collections::BTreeSet;

use tepsim_core::{State, Stream, equilibrium, flows, heat, math, streams, stripper, vessels};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::Comparison;
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, Snapshot, adversarial, compare_field};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2".
const TIER2_TOLERANCE: f64 = 1e-12;

struct Solved {
    heat: heat::HeatTransfer,
    condenser_delta: f64,
}

struct Field {
    name: &'static str,
    ours: fn(&Solved) -> Vec<f64>,
    theirs: fn(&Snapshot) -> Vec<f64>,
}

fn fields() -> Vec<Field> {
    vec![
        Field {
            name: "UAR",
            ours: |s| vec![s.heat.reactor_coefficient],
            theirs: |s| vec![s.common.uar],
        },
        Field {
            name: "QUR",
            ours: |s| vec![s.heat.reactor_duty],
            theirs: |s| vec![s.common.qur],
        },
        Field {
            name: "QUS",
            ours: |s| vec![s.heat.condenser_duty],
            theirs: |s| vec![s.common.qus],
        },
        Field {
            name: "QUC",
            ours: |s| vec![s.heat.stripper_duty],
            theirs: |s| vec![s.common.quc],
        },
    ]
}

/// Run the port through everything up to and including this item's range.
fn solve(oracle: &mut Oracle, scenario: &Scenario) -> Option<Solved> {
    let t = scenario.time;
    let raw = |n: usize| f64::from(scenario.disturbances[n - 1]);
    // teprob.f:407-408 subtracts *two* terms from A, on two source lines.
    // Every one of these harnesses dropped the second until B-0032, and it
    // never showed because no pooled scenario has a disturbance active.
    //
    // These stay independent of `Plant::advance_discrete` on purpose. If they
    // asked the plant for their own inputs, a bug in the plant would feed both
    // sides of the comparison and Tier 2 would pass on wrong-against-wrong.
    // `tier3_walk_inputs.rs` is what checks the plant's version, against the
    // oracle, with all twenty faults switched on.
    let a = oracle.tesub8(1, t) - raw(1) * 0.03 - raw(2) * 2.43719e-3;
    let b = oracle.tesub8(2, t) + raw(2) * 0.005;
    let feed = streams::FeedConditions {
        ac_feed_light: [a, b, 1.0 - a - b],
        d_feed_celsius: oracle.tesub8(3, t) + raw(3) * 5.0,
        ac_feed_celsius: oracle.tesub8(4, t),
    };
    let flow_drift = flows::FlowDrift {
        steam_capacity: oracle.tesub8(9, t),
        reactor_outlet: oracle.tesub8(12, t),
    };
    // teprob.f:673 and 676.
    let heat_drift = heat::HeatDrift {
        reactor_coolant: oracle.tesub8(10, t),
        condenser_coolant: oracle.tesub8(11, t),
    };
    let seeds = vessels::TemperatureSeeds {
        reactor: scenario.common.tcr,
        separator: scenario.common.tcs,
        stripper: scenario.common.tcc,
        mixing: scenario.common.tcv,
    };
    let state = State::from_flat(&scenario.state);
    let unpacked = vessels::unpack(&state, seeds).ok()?;
    let eq = equilibrium::equilibrium(&unpacked);
    let mut table = streams::streams(&unpacked, &eq, &feed);
    let mut idv = [0.0; 20];
    for (slot, r) in idv.iter_mut().zip(scenario.disturbances) {
        *slot = f64::from(r.clamp(0, 1));
    }
    let mut flow = flows::flows(&state, &unpacked, &eq, &table, &idv, flow_drift);
    let _ = stripper::stripper(&mut table, &mut flow, unpacked.stripper.celsius);
    let result = heat::heat_transfer(&state, &unpacked, &table, &flow, heat_drift);
    Some(Solved {
        heat: result,
        condenser_delta: state.condenser_cw_out_c - table.celsius[Stream::ReactorOutlet],
    })
}

fn observe(
    oracle: &mut Oracle,
    scenario: &Scenario,
    pool: Pool,
    index: usize,
    comparisons: &mut [Comparison<Case>],
) -> bool {
    let snapshot = scenario.force(oracle);
    let Some(solved) = solve(oracle, scenario) else {
        return false;
    };
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
        .map(|f| Comparison::new(format!("heat {}", f.name)))
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
fn heat_transfer_matches_the_fortran_over_all_three_pools() {
    let mut oracle = Oracle::lock();
    let (comparisons, skipped) = sweep(&mut oracle, 400, 2_000, 0x7E2_0023);

    println!(
        "transcendentals come from the {} libm",
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

/// With the transcendentals taken out of the comparison, the algebra must be
/// bit-identical. See `tier2_equilibrium.rs`.
#[test]
#[cfg(feature = "libm-system")]
fn the_algebra_is_bit_identical_once_exp_and_pow_agree() {
    let mut oracle = Oracle::lock();
    let (comparisons, _) = sweep(&mut oracle, 200, 500, 6);

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
         Either the algebra in `tepsim_core::heat` does not associate the way \
         `teprob.f:663-678` does, which is a port bug; or gfortran no longer \
         resolves the transcendentals to the same code as the platform \
         libm.\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// Both flat branches of the level ramp must be entered by the pool, not just
/// the ramp between them.
///
/// This is the check B-0022's finding says to write. B-0016 placed states at
/// `VLR = 78` and `VLR = 390`, the two breakpoints, and both comparisons at
/// `teprob.f:663` and `665` are strict, so *both boundary states take the
/// ramp*. Without states past each threshold, `UARLEV = 1` and `UARLEV = 0`
/// are implemented and never evaluated.
#[test]
fn every_level_branch_is_exercised_by_the_pool() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let mut seen: BTreeSet<&'static str> = BTreeSet::new();
    let name = |b: heat::LevelBranch| match b {
        heat::LevelBranch::FullyWetted => "fully wetted",
        heat::LevelBranch::Dry => "dry",
        heat::LevelBranch::Ramp => "ramp",
    };

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let _ = scenario.force(&mut oracle);
        if let Some(s) = solve(&mut oracle, &scenario) {
            seen.insert(name(s.heat.level_branch));
        }
    }
    let mut sampler = tepsim_oracle::tier1::Sampler::new(0x7E2_0023);
    for index in 0..1_000 {
        let scenario = pools.perturbed_case(index, &mut sampler);
        let _ = scenario.force(&mut oracle);
        if let Some(s) = solve(&mut oracle, &scenario) {
            seen.insert(name(s.heat.level_branch));
        }
    }
    let (boundaries, missed) = adversarial::build(&mut oracle, &pools.nominal_case(0));
    assert!(missed.is_empty(), "the catalogue regressed: {missed:?}");
    for boundary in &boundaries {
        let _ = boundary.scenario.force(&mut oracle);
        if let Some(s) = solve(&mut oracle, &boundary.scenario) {
            println!("{:52} {}", boundary.target.name, name(s.heat.level_branch));
            seen.insert(name(s.heat.level_branch));
        }
    }

    println!("level branches reached: {seen:?}");
    for branch in ["fully wetted", "dry", "ramp"] {
        assert!(
            seen.contains(branch),
            "no state in the pool reached the '{branch}' branch of UARLEV, so \
             it is implemented but unvalidated. Both breakpoint states take \
             the ramp, because teprob.f:663 and 665 are strict; a state past \
             each threshold is needed."
        );
    }
}

/// Both sides of the `TCC < 100` steam cutoff must be entered.
///
/// The nominal trajectory sits near 65 C, so the *cutoff* is the branch at
/// risk, not the active one.
#[test]
fn both_sides_of_the_steam_cutoff_are_exercised() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let (mut on, mut off) = (0, 0);

    let (boundaries, _) = adversarial::build(&mut oracle, &pools.nominal_case(0));
    let mut count = |s: &Solved, on: &mut i32, off: &mut i32| {
        if s.heat.steam_on {
            *on += 1;
        } else {
            *off += 1;
        }
    };
    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let _ = scenario.force(&mut oracle);
        if let Some(s) = solve(&mut oracle, &scenario) {
            count(&s, &mut on, &mut off);
        }
    }
    for boundary in &boundaries {
        let _ = boundary.scenario.force(&mut oracle);
        if let Some(s) = solve(&mut oracle, &boundary.scenario) {
            count(&s, &mut on, &mut off);
        }
    }

    println!("steam on: {on} states, off: {off} states");
    assert!(on > 0, "the reboiler was never active");
    assert!(
        off > 0,
        "no state reached TCC >= 100, so the steam cutoff at teprob.f:678 is \
         implemented and never evaluated. B-0016's 'TCC at the upper \
         stripping-factor branch' state sits at 170 C and should reach it."
    );
}

/// `QUS` must be driven by `TST(8)`, the reactor outlet, and the two are far
/// enough apart on real states that using the separator temperature instead
/// would show up immediately.
///
/// Confirms against the oracle that the reading of `teprob.f:675` is right,
/// rather than just internally consistent.
#[test]
fn the_condenser_driving_difference_is_large_enough_to_discriminate() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);
    let mut smallest = f64::INFINITY;

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let snapshot = scenario.force(&mut oracle);
        let Some(s) = solve(&mut oracle, &scenario) else {
            continue;
        };
        // The separator temperature is the plausible wrong choice.
        let wrong = snapshot.common.tws - snapshot.common.tcs;
        smallest = smallest.min((s.condenser_delta - wrong).abs());
    }
    println!("smallest gap between TST(8) and TCS as the driving temperature: {smallest:.4} C");
    assert!(
        smallest > 1.0,
        "the two candidate driving temperatures came within {smallest} C, so \
         the differential could not distinguish them"
    );
}
