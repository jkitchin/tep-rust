//! Tier 2 for the stripper: `tepsim_core::stripper` against
//! `teprob.f:614-662`, over all three sampling pools.
//!
//! Run twice, default and `--features oracle,libm-system`; see
//! `tier2_equilibrium.rs`.
//!
//! # The coverage existed before the code did
//!
//! This is the first Phase 2 item whose branches were all in the adversarial
//! catalogue before it was written. B-0016 built states at `TCC = 5.292`, at
//! `TCC = 170`, near the 177 C pole, and exactly on `FTM(11) = 0.1`, all from
//! reading the listing. So no new boundaries were needed here, and the
//! branch-coverage test below asserts that the pool really does reach every
//! branch rather than assuming the catalogue still does what it claims.

#![cfg(feature = "oracle")]

use std::collections::BTreeSet;

use tepsim_core::{Component, State, Stream, equilibrium, flows, math, streams, stripper, vessels};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::Comparison;
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, Snapshot, adversarial, compare_field};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2".
const TIER2_TOLERANCE: f64 = 1e-12;

/// The three streams this item fills in.
const STRIPPER_STREAMS: [Stream; 3] = [
    Stream::StripperOverhead,
    Stream::ReactorInlet,
    Stream::StripperDownflow,
];

/// The two the column actually computes, as opposed to the one it copies.
const COLUMN: [Stream; 2] = [Stream::StripperOverhead, Stream::StripperDownflow];

struct Solved {
    stripper: stripper::Stripper,
    streams: streams::Streams,
    flows: flows::Flows,
}

struct Field {
    name: &'static str,
    ours: fn(&Solved) -> Vec<f64>,
    theirs: fn(&Snapshot) -> Vec<f64>,
}

fn fields() -> Vec<Field> {
    vec![
        Field {
            name: "SFR",
            ours: |s| s.stripper.factors.as_array().to_vec(),
            theirs: |s| s.common.sfr.to_vec(),
        },
        // The column's own two outlets and the reactor-inlet alias are
        // reported separately, because they behave completely differently.
        // Everything the stripper computes is built from valve-lagged flows
        // and unpacked compositions, all of which are bit-exact, so streams 5
        // and 12 are bit-exact too. Stream 7 is a copy of stream 6, whose flow
        // came through a square-root resistance, so it inherits that error and
        // nothing here can remove it. Reporting them together would average a
        // clean result with a borrowed one and hide both facts.
        Field {
            name: "FTM (5, 12): the column's own outlets",
            ours: |s| vec![s.flows.molar[COLUMN[0]], s.flows.molar[COLUMN[1]]],
            theirs: |s| {
                vec![
                    s.common.ftm[COLUMN[0].index()],
                    s.common.ftm[COLUMN[1].index()],
                ]
            },
        },
        Field {
            name: "FTM (7): the reactor-inlet alias",
            ours: |s| vec![s.flows.molar[Stream::ReactorInlet]],
            theirs: |s| vec![s.common.ftm[Stream::ReactorInlet.index()]],
        },
        Field {
            name: "FCM (5, 12): the column's own outlets",
            ours: |s| {
                COLUMN
                    .iter()
                    .flat_map(|stream| *s.flows.component[*stream].as_array())
                    .collect()
            },
            theirs: |s| {
                COLUMN
                    .iter()
                    .flat_map(|stream| s.common.fcm[stream.index()])
                    .collect()
            },
        },
        Field {
            name: "FCM (7): the reactor-inlet alias",
            ours: |s| s.flows.component[Stream::ReactorInlet].as_array().to_vec(),
            theirs: |s| s.common.fcm[Stream::ReactorInlet.index()].to_vec(),
        },
        Field {
            name: "XST (5, 12): the column's own outlets",
            ours: |s| {
                COLUMN
                    .iter()
                    .flat_map(|stream| *s.streams.composition[*stream].fractions().as_array())
                    .collect()
            },
            theirs: |s| {
                COLUMN
                    .iter()
                    .flat_map(|stream| s.common.xst[stream.index()])
                    .collect()
            },
        },
        Field {
            name: "XST (7): the reactor-inlet alias",
            ours: |s| {
                s.streams.composition[Stream::ReactorInlet]
                    .fractions()
                    .as_array()
                    .to_vec()
            },
            theirs: |s| s.common.xst[Stream::ReactorInlet.index()].to_vec(),
        },
        // `TST` needs no split: 5 and 12 are `TCC` and 7 is `TCV`, and both
        // came out of `TESUB2` bit-exact in B-0017. It is the one field here
        // the alias does not contaminate.
        Field {
            name: "TST (5, 7, 12)",
            ours: |s| {
                STRIPPER_STREAMS
                    .iter()
                    .map(|stream| s.streams.celsius[*stream])
                    .collect()
            },
            theirs: |s| {
                STRIPPER_STREAMS
                    .iter()
                    .map(|stream| s.common.tst[stream.index()])
                    .collect()
            },
        },
        Field {
            name: "HST (5, 12): the column's own outlets",
            ours: |s| {
                COLUMN
                    .iter()
                    .map(|stream| s.streams.enthalpy[*stream])
                    .collect()
            },
            theirs: |s| {
                COLUMN
                    .iter()
                    .map(|stream| s.common.hst[stream.index()])
                    .collect()
            },
        },
        Field {
            name: "HST (7): the reactor-inlet alias",
            ours: |s| vec![s.streams.enthalpy[Stream::ReactorInlet]],
            theirs: |s| vec![s.common.hst[Stream::ReactorInlet.index()]],
        },
    ]
}

/// Run the port through everything up to and including this item's range.
fn solve(oracle: &mut Oracle, scenario: &Scenario) -> Option<Solved> {
    let t = scenario.time;
    let idv_raw = |n: usize| f64::from(scenario.disturbances[n - 1]);
    let a = oracle.tesub8(1, t) - idv_raw(1) * 0.03;
    let b = oracle.tesub8(2, t) + idv_raw(2) * 0.005;
    let feed = streams::FeedConditions {
        ac_feed_light: [a, b, 1.0 - a - b],
        d_feed_celsius: oracle.tesub8(3, t) + idv_raw(3) * 5.0,
        ac_feed_celsius: oracle.tesub8(4, t),
    };
    let drift = flows::FlowDrift {
        steam_capacity: oracle.tesub8(9, t),
        reactor_outlet: oracle.tesub8(12, t),
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
    for (slot, raw) in idv.iter_mut().zip(scenario.disturbances) {
        *slot = f64::from(raw.clamp(0, 1));
    }
    let mut flow = flows::flows(&state, &unpacked, &eq, &table, &idv, drift);
    let result = stripper::stripper(&mut table, &mut flow, unpacked.stripper.celsius);
    Some(Solved {
        stripper: result,
        streams: table,
        flows: flow,
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
        .map(|f| Comparison::new(format!("stripper {}", f.name)))
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
fn the_stripper_matches_the_fortran_over_all_three_pools() {
    let mut oracle = Oracle::lock();
    let (comparisons, skipped) = sweep(&mut oracle, 400, 2_000, 0x7E2_0022);

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

    // The stripper introduces no error of its own, and exactly two of its
    // eleven reported fields may carry any at all.
    //
    // Everything the column computes descends from valve-lagged flows and
    // unpacked compositions, and no transcendental reaches either, so `SFR`,
    // `FIN` and both outlet streams are bit-exact even under the vendored
    // libm. The alias inherits whatever stream 6 has, and stream 6's *flow*
    // came through a square-root resistance while its *composition* and
    // *temperature* came from the unpack. So `FTM(7)` and `FCM(7)` may differ
    // and `XST(7)`, `TST(7)` and `HST(7)` may not.
    //
    // Asserted rather than observed, and named field by field rather than by a
    // blanket exemption, so that error appearing anywhere it structurally
    // cannot fails here instead of merging into a 1e-12 tolerance.
    if !math::USES_SYSTEM_LIBM {
        const MAY_INHERIT: [&str; 2] = [
            "FTM (7): the reactor-inlet alias",
            "FCM (7): the reactor-inlet alias",
        ];
        for (field, comparison) in fields().iter().zip(comparisons.iter()) {
            if MAY_INHERIT.contains(&field.name) {
                continue;
            }
            assert_eq!(
                comparison.max_ulp(),
                0,
                "`{}` descends only from bit-exact quantities, so a difference \
                 here is the algebra and not libm:\n{comparison}",
                field.name
            );
        }
    }
}

/// With the transcendentals taken out of the comparison, the algebra must be
/// bit-identical. See `tier2_equilibrium.rs`.
#[test]
#[cfg(feature = "libm-system")]
fn the_algebra_is_bit_identical_once_exp_and_pow_agree() {
    let mut oracle = Oracle::lock();
    let (comparisons, _) = sweep(&mut oracle, 200, 500, 5);

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
         Either the algebra in `tepsim_core::stripper` does not associate the \
         way `teprob.f:614-662` does, which is a port bug; or gfortran no \
         longer resolves the transcendentals to the same code as the platform \
         libm.\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// All four branches must be reached by the pool, using the states B-0016
/// built before this code existed. An unexercised branch is not evidence.
#[test]
fn every_stripper_branch_is_exercised_by_the_pool() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let mut seen: BTreeSet<&'static str> = BTreeSet::new();
    let mut record = |branch: stripper::StripperBranch, seen: &mut BTreeSet<&'static str>| {
        seen.insert(match branch {
            stripper::StripperBranch::Idle => "idle",
            stripper::StripperBranch::Linear => "linear",
            stripper::StripperBranch::Pinned => "pinned",
            stripper::StripperBranch::Hyperbolic => "hyperbolic",
        });
    };

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let _ = scenario.force(&mut oracle);
        if let Some(s) = solve(&mut oracle, &scenario) {
            record(s.stripper.branch, &mut seen);
        }
    }
    let (boundaries, missed) = adversarial::build(&mut oracle, &pools.nominal_case(0));
    assert!(missed.is_empty(), "the catalogue regressed: {missed:?}");
    for boundary in &boundaries {
        let _ = boundary.scenario.force(&mut oracle);
        if let Some(s) = solve(&mut oracle, &boundary.scenario) {
            record(s.stripper.branch, &mut seen);
            println!("{:52} {:?}", boundary.target.name, s.stripper.branch);
        }
    }

    println!("branches reached: {seen:?}");
    for branch in ["idle", "linear", "pinned", "hyperbolic"] {
        assert!(
            seen.contains(branch),
            "no state in the pool reached the {branch} branch, so it is \
             implemented but unvalidated. B-0016 built states for all four; if \
             one no longer lands, the catalogue regressed."
        );
    }
}

/// `SFR(1..3)` must be bit-identical to the `TEINIT` constants on every state,
/// in both branches. If the Fortran ever reported something else, the reading
/// that they are never recomputed would be wrong.
#[test]
fn the_non_condensible_factors_never_move_in_the_fortran_either() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 300, DT);
    let expected = stripper::NON_CONDENSIBLE_STRIPPING;

    for index in 0..pools.trajectory.len() {
        let snapshot = pools.nominal_case(index).force(&mut oracle);
        for (slot, want) in expected.iter().enumerate() {
            assert_eq!(
                snapshot.common.sfr[slot].to_bits(),
                want.to_bits(),
                "SFR({}) moved to {} at nominal#{index}. Nothing in \
                 teprob.f:614-634 writes slots 1-3, so this would mean they \
                 are recomputed somewhere and the port is missing it.",
                slot + 1,
                snapshot.common.sfr[slot]
            );
        }
    }
    let (boundaries, _) = adversarial::build(&mut oracle, &pools.nominal_case(0));
    for boundary in &boundaries {
        let snapshot = boundary.scenario.force(&mut oracle);
        for (slot, want) in expected.iter().enumerate() {
            assert_eq!(
                snapshot.common.sfr[slot].to_bits(),
                want.to_bits(),
                "SFR({}) moved at boundary {}",
                slot + 1,
                boundary.target.name
            );
        }
    }
    println!("SFR(1..3) fixed at {expected:?} across 300 nominal and 20 adversarial states");
}

/// Stream 7 must equal stream 6 in the *Fortran*, not just in the port. That
/// is what makes reproducing it as a copy correct rather than convenient.
#[test]
fn the_reactor_inlet_is_an_alias_in_the_fortran_too() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);

    for index in 0..pools.trajectory.len() {
        let s = pools.nominal_case(index).force(&mut oracle);
        let (six, seven) = (
            Stream::MixingZoneOutlet.index(),
            Stream::ReactorInlet.index(),
        );
        assert_eq!(s.common.ftm[seven].to_bits(), s.common.ftm[six].to_bits());
        assert_eq!(s.common.hst[seven].to_bits(), s.common.hst[six].to_bits());
        assert_eq!(s.common.tst[seven].to_bits(), s.common.tst[six].to_bits());
        for c in Component::ALL {
            assert_eq!(
                s.common.xst[seven][c.index()].to_bits(),
                s.common.xst[six][c.index()].to_bits()
            );
            assert_eq!(
                s.common.fcm[seven][c.index()].to_bits(),
                s.common.fcm[six][c.index()].to_bits()
            );
        }
    }
}
