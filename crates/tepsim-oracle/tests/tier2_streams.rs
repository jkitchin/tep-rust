//! Tier 2 for the stream table: `tepsim_core::streams` against
//! `teprob.f:529-564`, over all three sampling pools.
//!
//! Run twice, default and `--features oracle,libm-system`; see
//! `tier2_equilibrium.rs`.
//!
//! # `HST(9)` cannot be compared where it lives
//!
//! `teprob.f:562` copies the separator vapour enthalpy into the purge, and
//! `teprob.f:601` then adds the compressor work to the recycle *only*. So
//! after `TEFUNC` returns, `COMMON`'s `HST(9)` is not what this module
//! computes for stream 9; it is that value plus `CPDH/FTM(9)`.
//!
//! `HST(10)` is untouched by anything after line 562, so it holds exactly the
//! pre-compressor value. Both of this module's stream 9 and stream 10
//! enthalpies are therefore compared against `HST(10)`, and `HST(9)` itself
//! becomes checkable in B-0021 once the compressor exists. Comparing against
//! `HST(9)` here would have failed, and the obvious reading of that failure is
//! that the enthalpy is wrong rather than that it is early.
//!
//! # The walk-driven feed values are derived, not read back
//!
//! `teprob.f:407-412` overwrites `XST(1..3, 4)`, `TST(1)` and `TST(4)` from
//! `TESUB8` and the `IDV` flags. Those are inputs to this module, so they could
//! simply be read out of `COMMON` after the call, and the comparison on those
//! slots would then be trivially true.
//!
//! They are recomputed from `TESUB8` and the scenario's `IDV` instead. It costs
//! four lines and it turns six trivially-equal comparisons into a real check on
//! the reading of `teprob.f:407-412`, which B-0024 will otherwise have to
//! establish from scratch.

#![cfg(feature = "oracle")]

use tepsim_core::{State, Stream, equilibrium, math, streams, vessels};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::Comparison;
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, Snapshot, adversarial, compare_field};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2".
const TIER2_TOLERANCE: f64 = 1e-12;

/// Which oracle slot each stream's enthalpy must be compared against.
///
/// The identity everywhere except stream 9; see the module documentation.
fn enthalpy_slot(stream: Stream) -> usize {
    match stream {
        Stream::Recycle => Stream::Purge.index(),
        other => other.index(),
    }
}

struct Field {
    name: &'static str,
    ours: fn(&streams::Streams) -> Vec<f64>,
    theirs: fn(&Snapshot) -> Vec<f64>,
}

fn fields() -> Vec<Field> {
    vec![
        Field {
            name: "XST (10 streams x 8)",
            ours: |s| {
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .flat_map(|stream| s.composition[*stream].fractions().as_array().to_vec())
                    .collect()
            },
            theirs: |s| {
                // `xst` mirrors `XST(8,13)`, so the outer index is the stream
                // and the inner one the component: the Fortran's column-major
                // layout already puts the eight components together.
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .flat_map(|stream| s.common.xst[stream.index()])
                    .collect()
            },
        },
        Field {
            name: "XMWS (the 6 that exist)",
            ours: |s| {
                streams::WEIGHED_STREAMS
                    .iter()
                    .map(|stream| s.molar_mass[*stream])
                    .collect()
            },
            theirs: |s| {
                streams::WEIGHED_STREAMS
                    .iter()
                    .map(|stream| s.common.xmws[stream.index()])
                    .collect()
            },
        },
        Field {
            name: "TST",
            ours: |s| {
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .map(|stream| s.celsius[*stream])
                    .collect()
            },
            theirs: |s| {
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .map(|stream| s.common.tst[stream.index()])
                    .collect()
            },
        },
        Field {
            name: "HST",
            ours: |s| {
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .map(|stream| s.enthalpy[*stream])
                    .collect()
            },
            theirs: |s| {
                streams::ASSEMBLED_STREAMS
                    .iter()
                    .map(|stream| s.common.hst[enthalpy_slot(*stream)])
                    .collect()
            },
        },
    ]
}

/// Reproduce `teprob.f:407-412`, the walk-driven feed conditions.
fn feed_conditions(oracle: &mut Oracle, scenario: &Scenario) -> streams::FeedConditions {
    let t = scenario.time;
    let idv = |n: usize| f64::from(scenario.disturbances[n - 1]);
    // teprob.f:407-408
    // teprob.f:407-408 subtracts *two* terms from A, on two source lines.
    // Every one of these harnesses dropped the second until B-0032, and it
    // never showed because no pooled scenario has a disturbance active.
    //
    // These stay independent of `Plant::advance_discrete` on purpose. If they
    // asked the plant for their own inputs, a bug in the plant would feed both
    // sides of the comparison and Tier 2 would pass on wrong-against-wrong.
    // `tier3_walk_inputs.rs` is what checks the plant's version, against the
    // oracle, with all twenty faults switched on.
    let a = oracle.tesub8(1, t) - idv(1) * 0.03 - idv(2) * 2.43719e-3;
    // teprob.f:409
    let b = oracle.tesub8(2, t) + idv(2) * 0.005;
    // teprob.f:410
    let c = 1.0 - a - b;
    streams::FeedConditions {
        ac_feed_light: [a, b, c],
        // teprob.f:411
        d_feed_celsius: oracle.tesub8(3, t) + idv(3) * 5.0,
        // teprob.f:412
        ac_feed_celsius: oracle.tesub8(4, t),
    }
}

fn observe(
    oracle: &mut Oracle,
    scenario: &Scenario,
    pool: Pool,
    index: usize,
    comparisons: &mut [Comparison<Case>],
) -> bool {
    let snapshot = scenario.force(oracle);
    let feed = feed_conditions(oracle, scenario);
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
    let solved = streams::streams(&unpacked, &eq, &feed);
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
        .map(|f| Comparison::new(format!("streams {}", f.name)))
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
fn the_stream_table_matches_the_fortran_over_all_three_pools() {
    let mut oracle = Oracle::lock();
    let (comparisons, skipped) = sweep(&mut oracle, 400, 2_000, 0x7E2_0020);

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

/// With the transcendentals taken out of the comparison, the algebra must be
/// bit-identical. See `tier2_equilibrium.rs` for what a failure means.
#[test]
#[cfg(feature = "libm-system")]
fn the_algebra_is_bit_identical_once_exp_and_pow_agree() {
    let mut oracle = Oracle::lock();
    let (comparisons, _) = sweep(&mut oracle, 200, 500, 3);

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
         Either the algebra in `tepsim_core::streams` does not associate the \
         way `teprob.f:529-564` does, which is a port bug; or gfortran no \
         longer resolves the transcendentals to the same code as the platform \
         libm, in which case this configuration no longer removes them.\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// The seven `XMWS` slots the original never writes must be zero on both
/// sides, and none of them may be read anywhere.
///
/// If the Fortran ever reported a non-zero one, this module's decision to
/// leave them at zero would be hiding a real assignment somewhere in the file.
#[test]
fn the_unweighed_streams_are_zero_in_both_implementations() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);
    let unweighed: Vec<Stream> = Stream::ALL
        .into_iter()
        .filter(|s| !streams::WEIGHED_STREAMS.contains(s))
        .collect();
    assert_eq!(unweighed.len(), 7, "six of thirteen are weighed");

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let snapshot = scenario.force(&mut oracle);
        for stream in &unweighed {
            assert_eq!(
                snapshot.common.xmws[stream.index()].to_bits(),
                0.0_f64.to_bits(),
                "the Fortran reported a non-zero XMWS({}) at nominal#{index}. \
                 Nothing in teprob.f writes that slot, so either it does after \
                 all or it is picking up stale COMMON.",
                stream.fortran_index()
            );
        }
    }
    println!(
        "{} unweighed streams, all zero over 200 states",
        unweighed.len()
    );
}

/// `HST(9)` and `HST(10)` are equal on return only if the compressor adds
/// nothing. They are not, which is what makes comparing stream 9 against
/// `HST(10)` a real choice rather than a cosmetic one.
///
/// This pins the reason, so that a future session does not "simplify" the slot
/// mapping away and get a failure it cannot explain.
#[test]
fn the_compressor_makes_the_recycle_enthalpy_differ_from_the_purge() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);
    let mut differing = 0;

    for index in 0..pools.trajectory.len() {
        let snapshot = pools.nominal_case(index).force(&mut oracle);
        if snapshot.common.hst[Stream::Recycle.index()].to_bits()
            != snapshot.common.hst[Stream::Purge.index()].to_bits()
        {
            differing += 1;
        }
    }
    println!("HST(9) differs from HST(10) on {differing} of 200 nominal states");
    assert_eq!(
        differing, 200,
        "HST(9) equalled HST(10) on some states, so teprob.f:601 added nothing \
         there and this test no longer demonstrates why stream 9 is compared \
         against slot 10"
    );
}
