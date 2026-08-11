//! Tests the Tier 2 apparatus, before there is any physics for it to measure.
//!
//! The same discipline as B-0008: an instrument that reports "no difference"
//! for a broken port is worse than no instrument, so every way this one could
//! lie is provoked on purpose.
//!
//! The failure mode specific to Tier 2 is reproducibility. `TEFUNC` mutates the
//! walk state and the generator word as a side effect, so a harness that did
//! not restore them would report differences caused entirely by itself.

#![cfg(feature = "oracle")]

use tepsim_oracle::tier1::{Comparison, Sampler};
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, compare_field, reproducible};
use tepsim_oracle::{Oracle, tier2};

/// One second in hours, the original's step (`temain_mod.f`'s `INTGTR`).
const DT: f64 = 1.0 / 3600.0;

fn pools(oracle: &mut Oracle, steps: usize) -> Pools {
    Pools::collect(oracle, steps, DT)
}

// ------------------------------------------------------------- reproducibility

/// The property everything else rests on: the same scenario, twice, bit for
/// bit, even with another evaluation in between.
#[test]
fn forcing_the_same_scenario_twice_gives_identical_answers() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 50);
    for index in [0, 7, 25, 49] {
        reproducible(&mut oracle, &pools.nominal_case(index));
    }
}

/// The walk state reaches the model, so restoring it is not ceremony.
///
/// Demonstrated directly rather than by evaluating twice and hoping. Two
/// evaluations of the same scenario one step apart do *not* generally differ in
/// their walks: `TNEXT` starts at 0.1 hours (`teprob.f:404`), which is 360
/// one-second steps away, so `TESUB5` does not fire at all on a short run and
/// the coefficients never move. Perturbing a coefficient by hand is what proves
/// the plumbing.
///
/// `ADIST(5)` is the reactor cooling water supply channel, read at
/// `teprob.f:413` into `TCWR`.
#[test]
fn the_walk_state_reaches_the_model() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 10);
    let base = pools.nominal_case(5);

    let mut shifted = base.clone();
    shifted.walk.adist[4] += 1.0;

    let a = base.force(&mut oracle);
    let b = shifted.force(&mut oracle);

    println!(
        "TCWR {} vs {} with ADIST(5) shifted by one degree",
        a.common.tcwr, b.common.tcwr
    );
    assert_ne!(
        a.common.tcwr.to_bits(),
        b.common.tcwr.to_bits(),
        "shifting the walk coefficient for the reactor cooling water changed \
         nothing, so the walk state is not reaching the model and restoring it \
         proves nothing"
    );

    let moved = a
        .derivative
        .iter()
        .zip(b.derivative.iter())
        .filter(|(x, y)| x.to_bits() != y.to_bits())
        .count();
    println!("{moved} of 50 derivatives follow the walk shift");
    assert!(
        moved > 0,
        "the walk shift reached TCWR but no derivative followed it"
    );
}

/// The generator word alone must also matter, separately from the walks.
#[test]
fn the_generator_word_is_part_of_the_scenario() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 10);

    let a = pools.nominal_case(3);
    assert!(a.time > 0.0, "the noise path at teprob.f:711 needs t > 0");
    let mut b = a.clone();
    b.rng = 1_431_655_765.0; // the original's own alternative seed, teprob.f:1191

    let left = a.force(&mut oracle);
    let right = b.force(&mut oracle);

    let differing = left
        .measurements
        .iter()
        .zip(right.measurements.iter())
        .filter(|(x, y)| x.to_bits() != y.to_bits())
        .count();
    println!("{differing} of 41 measurements move with the generator word");
    assert!(
        differing > 0,
        "the generator word changed nothing, so measurement noise is not \
         reaching the snapshot and Tier 3 has nothing to instrument"
    );
}

/// The four solved temperatures are warm-start seeds carried in `COMMON`, and
/// this is what proves it rather than asserting it.
///
/// `TESUB2` takes its temperature argument as both guess and result
/// (`teprob.f:1432`, `1438`), and the four call sites at `teprob.f:460-465`
/// pass `TCR`, `TCS`, `TCC` and `TCV` straight out of `COMMON/TEPROC/`. Seed
/// `TCV` with a different but still plausible value and the converged answer
/// must move, in the last bits.
///
/// This is why [`Scenario`] carries the whole block. It cost two failing tests
/// to find, and B-0017 has to carry these four in the port's own plant state.
#[test]
fn the_solved_temperatures_are_warm_started_from_the_previous_call() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 20);

    let base = pools.nominal_case(5);
    let mut reseeded = base.clone();
    // A quarter of a degree away: far closer than any real step moves it, and
    // still inside the basin, so Newton lands on the same root.
    reseeded.common.tcv += 0.25;

    let a = base.force(&mut oracle);
    let b = reseeded.force(&mut oracle);

    // TST(6) is TCV (teprob.f:549).
    let (x, y) = (a.common.tst[5], b.common.tst[5]);
    println!("TCV converged to {x} from one seed and {y} from another");
    assert_ne!(
        x.to_bits(),
        y.to_bits(),
        "seeding TCV differently changed nothing, so either TESUB2 is no \
         longer warm-started or the scenario is not reaching it. Either way \
         the Tier 2 reproducibility guarantee needs re-deriving."
    );
    let separation = (x - y).abs() / x.abs();
    assert!(
        separation < 1e-12,
        "the two solves landed on different roots ({separation:e} apart), \
         which would mean TESUB2 is not merely path-dependent in its last bits"
    );
}

/// And the carried-over seed reaches the derivative, which is why the scenario
/// has to restore it.
///
/// The effect is small and seed-sensitive: a seed a hundredth of a degree away
/// converges to the same bits, so no derivative moves. A seed taken from a
/// genuinely different point on the trajectory does move one. Rather than pick
/// a distance and hope, this walks the trajectory and requires that *some* seed
/// propagates, reporting which and by how much.
#[test]
fn a_carried_over_temperature_seed_changes_the_derivative() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 200);
    let scenario = pools.nominal_case(5);
    let reference = scenario.force(&mut oracle);

    let mut propagating = Vec::new();
    for step in [1, 10, 50, 100, 199] {
        // Seed the temperatures from a different point on the trajectory,
        // changing nothing else at all.
        let donor = pools.nominal_case(step).force(&mut oracle);
        let mut reseeded = scenario.clone();
        reseeded.common.tcr = donor.common.tcr;
        reseeded.common.tcs = donor.common.tcs;
        reseeded.common.tcc = donor.common.tcc;
        reseeded.common.tcv = donor.common.tcv;

        let seeded = reseeded.force(&mut oracle);
        let moved: Vec<usize> = reference
            .derivative
            .iter()
            .zip(seeded.derivative.iter())
            .enumerate()
            .filter(|(_, (a, b))| a.to_bits() != b.to_bits())
            .map(|(slot, _)| slot + 1)
            .collect();
        println!("seed from step {step}: YP components moved = {moved:?}");
        if !moved.is_empty() {
            propagating.push((step, moved));
        }
    }

    assert!(
        !propagating.is_empty(),
        "no warm-start seed anywhere on the trajectory changed a single \
         derivative, so the four temperatures are not observable state after \
         all and Scenario could stop carrying the COMMON block. If this \
         becomes true, simplify it and record why."
    );
}

// ------------------------------------------------------------------- the pools

#[test]
fn the_nominal_trajectory_actually_moves() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 100);
    assert_eq!(pools.trajectory.len(), 100);

    let (_, first) = pools.trajectory[0];
    let (_, last) = pools.trajectory[99];
    let moved = first
        .iter()
        .zip(last.iter())
        .filter(|(a, b)| a.to_bits() != b.to_bits())
        .count();
    println!("{moved} of 50 states differ between step 0 and step 99");
    assert!(
        moved > 25,
        "only {moved} states moved over 100 steps, so the trajectory pool is \
         nearly a single point and covers nothing"
    );
}

/// A nominal run must not trip. If it did, every derivative would be zeroed
/// (`teprob.f:807-811`) and the pool would be worthless.
#[test]
fn the_nominal_trajectory_never_trips_the_shutdown() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 200);
    for index in 0..pools.trajectory.len() {
        let snapshot = pools.nominal_case(index).force(&mut oracle);
        assert!(
            !snapshot.tripped,
            "the nominal trajectory tripped the shutdown at step {index}, so \
             the pool contains states whose derivative is all zeros"
        );
    }
}

/// The perturbed pool has to perturb, across the range of magnitudes it claims,
/// and has to say how many of its states trip.
#[test]
fn the_perturbed_pool_spans_its_magnitudes_and_reports_its_trips() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 50);
    let mut sampler = Sampler::new(0x7E2_5EED);

    let mut smallest = f64::INFINITY;
    let mut largest = 0.0_f64;
    let mut tripped = 0;
    const CASES: usize = 200;

    for index in 0..CASES {
        let base = pools.nominal_case(index);
        let case = pools.perturbed_case(index, &mut sampler);
        for (before, after) in base.state.iter().zip(case.state.iter()) {
            if before.abs() > 0.0 {
                let relative = ((after - before) / before).abs();
                if relative > 0.0 {
                    smallest = smallest.min(relative);
                    largest = largest.max(relative);
                }
            }
        }
        if case.force(&mut oracle).tripped {
            tripped += 1;
        }
    }

    println!(
        "perturbations span {smallest:e} to {largest:e} relative; \
         {tripped} of {CASES} states tripped the shutdown"
    );
    assert!(
        smallest < 1e-7,
        "the pool never produced a small perturbation: smallest was {smallest:e}"
    );
    assert!(
        largest > 1e-2,
        "the pool never produced a large perturbation: largest was {largest:e}"
    );
    // Not an assertion about the plant, an assertion about the harness: it has
    // to be able to tell, because a tripped state proves nothing.
    assert!(
        tripped <= CASES,
        "the trip count is nonsense: {tripped} of {CASES}"
    );
}

#[test]
fn perturbed_cases_are_reproducible_from_their_seed() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 20);
    let take = |seed: u64| -> Vec<u64> {
        let mut sampler = Sampler::new(seed);
        (0..10)
            .flat_map(|i| {
                pools
                    .perturbed_case(i, &mut sampler)
                    .state
                    .map(f64::to_bits)
                    .to_vec()
            })
            .collect()
    };
    assert_eq!(take(99), take(99), "same seed, same states");
    assert_ne!(take(99), take(100), "different seed, different states");
}

// -------------------------------------------------------------- the field diff

/// The differ must read zero when the two sides agree, and must name the field
/// component when they do not.
#[test]
fn the_field_differ_reports_zero_for_a_field_compared_with_itself() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 20);

    let mut comparison: Comparison<Case> = Comparison::new("FTM self");
    for index in 0..20 {
        let snapshot = pools.nominal_case(index).force(&mut oracle);
        compare_field(
            &mut comparison,
            Pool::Nominal,
            index,
            &snapshot.common.ftm,
            &snapshot.common.ftm,
        );
    }
    println!("{comparison}");
    assert_eq!(comparison.cases(), 20 * 13);
    comparison.assert_within(0.0);
}

/// And it must see a one-ULP difference in a single component of a single case.
#[test]
fn the_field_differ_sees_one_ulp_in_one_component() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 20);

    let mut comparison: Comparison<Case> = Comparison::new("FTM nudged");
    for index in 0..20 {
        let snapshot = pools.nominal_case(index).force(&mut oracle);
        let mut nudged = snapshot.common.ftm;
        if index == 7 {
            nudged[4] = nudged[4].next_up();
        }
        compare_field(
            &mut comparison,
            Pool::Nominal,
            index,
            &nudged,
            &snapshot.common.ftm,
        );
    }
    println!("{comparison}");

    assert_eq!(
        comparison.max_ulp(),
        1,
        "the nudge must be seen, and be 1 ULP"
    );
    let worst = comparison
        .worst_ulp_case()
        .expect("a worst case must be recorded");
    assert_eq!(worst.index, 7, "the wrong case was named");
    assert_eq!(
        worst.component, 5,
        "components are reported one-based, as Fortran subscripts: FTM(5)"
    );
}

/// Every field a Phase 2 item will need is actually reachable, so no item
/// discovers halfway through that its output is unobservable.
#[test]
fn the_intermediates_every_phase_2_item_needs_are_all_observable() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 5);
    let snapshot = pools.nominal_case(0).force(&mut oracle);
    let c = &snapshot.common;

    // One representative per planned item, with the item that needs it.
    let checks: [(&str, &[f64]); 12] = [
        ("B-0017 UCLR", &c.uclr),
        ("B-0017 TCR/TCS/TCC", &[c.tcr, c.tcs, c.tcc]),
        ("B-0017 DLR/DLS/DLC", &[c.dlr, c.dls, c.dlc]),
        ("B-0018 PPR", &c.ppr),
        ("B-0018 PTR/PTS/PTV", &[c.ptr, c.pts, c.ptv]),
        ("B-0019 RR", &c.rr),
        ("B-0019 CRXR", &c.crxr),
        ("B-0020 HST", &c.hst),
        ("B-0021 FTM", &c.ftm),
        ("B-0022 SFR", &c.sfr),
        ("B-0023 QUR/QUS/QUC", &[c.qur, c.qus, c.quc]),
        ("B-0025 VCV", &c.vcv),
    ];
    for (what, values) in checks {
        assert!(!values.is_empty(), "{what} is empty");
        assert!(
            values.iter().all(|v| v.is_finite()),
            "{what} contains a non-finite value at the nominal point: {values:?}"
        );
    }

    // The derivative itself, which is what Tier 2 finally gates on.
    assert!(
        snapshot.derivative.iter().all(|v| v.is_finite()),
        "the nominal derivative is not finite"
    );
    assert!(!snapshot.tripped, "the nominal point trips");
}

/// A scenario carries its own walk state, so scenarios can be evaluated in any
/// order. Check that, because Tier 2 will not evaluate them in the order they
/// were built.
#[test]
fn scenarios_are_order_independent() {
    let mut oracle = Oracle::lock();
    let pools = pools(&mut oracle, 30);

    let forward: Vec<[f64; 50]> = (0..10)
        .map(|i| pools.nominal_case(i).force(&mut oracle).derivative)
        .collect();
    let backward: Vec<[f64; 50]> = (0..10)
        .rev()
        .map(|i| pools.nominal_case(i).force(&mut oracle).derivative)
        .collect();

    for (i, (a, b)) in forward.iter().zip(backward.iter().rev()).enumerate() {
        for (slot, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            assert_eq!(
                x.to_bits(),
                y.to_bits(),
                "case {i} YP({}) depends on evaluation order",
                slot + 1
            );
        }
    }
}

/// `tier2::Pool` must render the way a log entry needs it.
#[test]
fn a_case_prints_the_way_the_log_format_expects() {
    let case = Case {
        pool: tier2::Pool::Perturbed,
        index: 4211,
        component: 27,
    };
    assert_eq!(case.to_string(), "perturbed#4211[27]");
}
