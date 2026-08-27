//! Tier 2 for the flow network: `tepsim_core::flows` against
//! `teprob.f:565-613`, over all three sampling pools.
//!
//! Run twice, default and `--features oracle,libm-system`; see
//! `tier2_equilibrium.rs`.
//!
//! # `HST(9)` becomes checkable here
//!
//! B-0020 could only compare its stream 9 enthalpy against `COMMON`'s
//! `HST(10)`, because `teprob.f:601` adds the compressor work to `HST(9)`
//! after the stream table is built and B-0020 had no compressor. It does now,
//! so this file compares the *post-bump* value against slot 9 directly. The
//! two checks together pin both sides of that line.
//!
//! # `UAC` is a local, so it is checked through `QUC`
//!
//! `teprob.f:326` declares `UAC` among `TEFUNC`'s locals rather than in
//! `COMMON/TEPROC/`, so there is no symbol to read it back from and the
//! oracle cannot expose it. Hoisting it the way `ISD` was hoisted would work
//! and is not worth an edit to the instrumented source for one number.
//!
//! `teprob.f:678` computes `QUC = UAC*(100.0-TCC)` whenever `TCC < 100`, and
//! `QUC` *is* in `COMMON`. So the port's `UAC` is multiplied by the same
//! factor and compared against `QUC`, which is exact: it is the identical
//! expression, not a division that would introduce rounding of its own. States
//! at or above 100 C are counted and excluded rather than silently passed,
//! since there `QUC` is zero regardless of `UAC` and proves nothing.

#![cfg(feature = "oracle")]

use tepsim_core::constants::single;
use tepsim_core::{Component, State, Stream, equilibrium, flows, math, streams, vessels};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::Comparison;
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, Snapshot, adversarial, compare_field};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2".
const TIER2_TOLERANCE: f64 = 1e-12;

/// What the port produced for one scenario, alongside what it needs from the
/// rest of the model to be comparable.
struct Solved {
    flows: flows::Flows,
    stripper_celsius: f64,
}

struct Field {
    name: &'static str,
    ours: fn(&Solved) -> Vec<f64>,
    theirs: fn(&Snapshot) -> Vec<f64>,
}

fn fields() -> Vec<Field> {
    vec![
        Field {
            name: "FTM (the 10 assembled streams)",
            ours: |s| {
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .map(|stream| s.flows.molar[*stream])
                    .collect()
            },
            theirs: |s| {
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .map(|stream| s.common.ftm[stream.index()])
                    .collect()
            },
        },
        Field {
            name: "FCM (10 streams x 8)",
            ours: |s| {
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .flat_map(|stream| *s.flows.component[*stream].as_array())
                    .collect()
            },
            theirs: |s| {
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .flat_map(|stream| s.common.fcm[stream.index()])
                    .collect()
            },
        },
        Field {
            name: "FWR/FWS/AGSP",
            ours: |s| {
                vec![
                    s.flows.reactor_coolant,
                    s.flows.condenser_coolant,
                    s.flows.agitator,
                ]
            },
            theirs: |s| vec![s.common.fwr, s.common.fws, s.common.agsp],
        },
        Field {
            name: "CPDH",
            ours: |s| vec![s.flows.compressor_work],
            theirs: |s| vec![s.common.cpdh],
        },
        Field {
            name: "HST(9) after the compressor bump",
            ours: |s| vec![s.flows.recycle_enthalpy],
            theirs: |s| vec![s.common.hst[Stream::Recycle.index()]],
        },
    ]
}

/// Reproduce `teprob.f:407-412`, as `tier2_streams.rs` does.
fn feed_conditions(oracle: &mut Oracle, scenario: &Scenario) -> streams::FeedConditions {
    let t = scenario.time;
    let idv = |n: usize| f64::from(scenario.disturbances[n - 1]);
    let a = oracle.tesub8(1, t) - idv(1) * 0.03;
    let b = oracle.tesub8(2, t) + idv(2) * 0.005;
    streams::FeedConditions {
        ac_feed_light: [a, b, 1.0 - a - b],
        d_feed_celsius: oracle.tesub8(3, t) + idv(3) * 5.0,
        ac_feed_celsius: oracle.tesub8(4, t),
    }
}

/// Run the port to the end of this item's range.
fn solve(oracle: &mut Oracle, scenario: &Scenario) -> Option<Solved> {
    let feed = feed_conditions(oracle, scenario);
    // teprob.f:572 and 583.
    let drift = flows::FlowDrift {
        steam_capacity: oracle.tesub8(9, scenario.time),
        reactor_outlet: oracle.tesub8(12, scenario.time),
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
    let table = streams::streams(&unpacked, &eq, &feed);
    // `teprob.f:341-344` clamps every IDV to exactly 0 or 1 before use.
    let mut idv = [0.0; 20];
    for (slot, raw) in idv.iter_mut().zip(scenario.disturbances) {
        *slot = f64::from(raw.clamp(0, 1));
    }
    Some(Solved {
        flows: flows::flows(&state, &unpacked, &eq, &table, &idv, drift),
        stripper_celsius: unpacked.stripper.celsius,
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
        .map(|f| Comparison::new(format!("flows {}", f.name)))
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
fn the_flow_network_matches_the_fortran_over_all_three_pools() {
    let mut oracle = Oracle::lock();
    let (comparisons, skipped) = sweep(&mut oracle, 400, 2_000, 0x7E2_0021);

    println!(
        "exp, pow and sqrt come from the {} libm",
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
    let (comparisons, _) = sweep(&mut oracle, 200, 500, 4);

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
         Either the algebra in `tepsim_core::flows` does not associate the way \
         `teprob.f:565-613` does, which is a port bug; or gfortran no longer \
         resolves the transcendentals to the same code as the platform \
         libm.\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// `UAC` through `QUC`; see the module documentation for why it is indirect.
#[test]
fn the_steam_coefficient_matches_through_the_condenser_duty() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 300, DT);
    let mut comparison: Comparison<Case> = Comparison::new("UAC, via QUC");
    let (mut checked, mut too_hot) = (0, 0);

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let snapshot = scenario.force(&mut oracle);
        let Some(solved) = solve(&mut oracle, &scenario) else {
            continue;
        };
        // teprob.f:677-678. At or above 100 C the duty is zero whatever UAC
        // is, so those states are not evidence.
        if solved.stripper_celsius >= 100.0 {
            too_hot += 1;
            continue;
        }
        checked += 1;
        comparison.observe(
            Case {
                pool: Pool::Nominal,
                index,
                component: 1,
            },
            solved.flows.steam_coefficient * (single(100.0) - solved.stripper_celsius),
            snapshot.common.quc,
        );
    }

    println!("{comparison}");
    println!("{checked} states below 100 C, {too_hot} at or above and excluded");
    assert!(
        checked > 100,
        "only {checked} states could check UAC, which is not enough to call it \
         validated"
    );
    comparison.assert_within(TIER2_TOLERANCE);
}

/// The two compressor pressure-ratio clamps must actually be exercised by the
/// pool, and the port must agree with the Fortran about which states hit them.
///
/// An unexercised branch is not evidence, and B-0016 could not build these
/// states because `CPPRMX` did not exist yet.
#[test]
fn the_compressor_ratio_clamps_are_exercised_by_the_adversarial_pool() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 100, DT);
    let (boundaries, missed) = adversarial::build(&mut oracle, &pools.nominal_case(0));
    assert!(missed.is_empty(), "the catalogue regressed: {missed:?}");

    let mut seen = Vec::new();
    for boundary in &boundaries {
        let _ = boundary.scenario.force(&mut oracle);
        let Some(solved) = solve(&mut oracle, &boundary.scenario) else {
            continue;
        };
        seen.push((boundary.target.name, solved.flows.pressure_ratio_clamp));
    }

    for (name, clamp) in &seen {
        println!("{name:52} {clamp:?}");
    }
    let any = |want: flows::RatioClamp| seen.iter().any(|(_, c)| *c == want);
    assert!(
        any(flows::RatioClamp::Low),
        "no adversarial state clamped the ratio up to 1"
    );
    assert!(
        any(flows::RatioClamp::High),
        "no adversarial state clamped the ratio down to CPPRMX"
    );
    assert!(
        any(flows::RatioClamp::None),
        "every state clamped, so the unclamped path is untested"
    );

    // `teprob.f:591` is `IF(PR.GT.CPPRMX)`, strictly greater, so the state
    // sitting exactly on `CPPRMX` must *not* clamp. Reading that line as
    // `.GE.` would be invisible everywhere except here.
    let at_boundary = seen
        .iter()
        .find(|(name, _)| *name == "PR = CPPRMX, the compressor maximum-ratio clamp")
        .expect("the boundary state is in the catalogue");
    assert_eq!(
        at_boundary.1,
        flows::RatioClamp::None,
        "a ratio of exactly CPPRMX clamped, so the comparison was read as          `>=` rather than the `.GT.` teprob.f:591 actually writes"
    );
}

/// The purge clamp cannot be reached from any sampled state. That is a
/// measured fact from B-0016, and it is re-measured here because this is the
/// item that implements the branch.
#[test]
fn no_sampled_state_reaches_the_purge_clamp() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let mut lowest = f64::INFINITY;
    let mut sampler = tepsim_oracle::tier1::Sampler::new(0x7E2_0021);

    let mut check = |oracle: &mut Oracle, scenario: &Scenario, lowest: &mut f64| {
        let snapshot = scenario.force(oracle);
        *lowest = lowest.min(snapshot.common.pts);
    };
    for index in 0..pools.trajectory.len() {
        check(&mut oracle, &pools.nominal_case(index), &mut lowest);
    }
    for index in 0..1_000 {
        let scenario = pools.perturbed_case(index, &mut sampler);
        check(&mut oracle, &scenario, &mut lowest);
    }
    let (boundaries, _) = adversarial::build(&mut oracle, &pools.nominal_case(0));
    for boundary in &boundaries {
        check(&mut oracle, &boundary.scenario, &mut lowest);
    }

    println!("lowest PTS over the whole pool: {lowest:.2} mmHg against a 760 threshold");
    assert!(
        lowest > 760.0,
        "PTS reached {lowest}, below the purge clamp threshold. The clamp is \
         reachable after all, so it should be in the adversarial catalogue \
         rather than covered only by a unit test."
    );
}

/// Component flows must sum back to the total, up to the composition's own
/// residual. A cheap independent statement of `teprob.f:602-613`.
#[test]
fn component_flows_sum_to_the_stream_total() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 100, DT);

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let _ = scenario.force(&mut oracle);
        let Some(solved) = solve(&mut oracle, &scenario) else {
            continue;
        };
        for stream in streams::ASSEMBLED_STREAMS {
            let mut sum = 0.0;
            for c in Component::ALL {
                sum += solved.flows.component[stream][c];
            }
            let total = solved.flows.molar[stream];
            if total.abs() < 1e-30 {
                continue;
            }
            let relative = (sum - total).abs() / total.abs();
            assert!(
                relative < 1e-6,
                "{} component flows sum to {sum} against a total of {total}",
                stream.name()
            );
        }
    }
}
