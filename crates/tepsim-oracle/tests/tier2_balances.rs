//! Tier 2 for the balances: all fifty derivatives, `tepsim_core::plant`
//! against `teprob.f:762-811`, over all three sampling pools.
//!
//! This is the comparison `PLAN.org` describes when it says "Tier 2": the
//! whole right-hand side, end to end. Every earlier `tier2_*` file compares an
//! intermediate; this one compares the answer.
//!
//! Run twice, default and `--features oracle,libm-system`; see
//! `tier2_equilibrium.rs`.
//!
//! # The 1e-12 relative gate is not met, and not because of a porting error
//!
//! **Tier 2 acceptance is BLOCKED on this; see B-0026a.** The numbers are
//! here rather than in the log alone because this file is where a future
//! session will look.
//!
//! Under the platform libm the whole right-hand side is *bit-identical* to the
//! Fortran, all fifty components, every state. So the algebra is exactly
//! right. Under the vendored libm, 28 of the 50 exceed 1e-12 **relative**, the
//! worst being `YP(2)` at 1.393e-4.
//!
//! The reason is cancellation, not error. A balance is
//! \(\sum_\text{in} - \sum_\text{out}\), and at anything near steady state
//! those two sums nearly agree: `YP(2)`, the inert's reactor balance, is a
//! difference of two flows around 660 whose result is a few parts in ten
//! thousand of either. A one-ULP difference in each term, which is all the
//! vendored `exp` costs, is 1e-16 of the *terms* and 1e-4 of the *result*.
//!
//! The twenty-two components that do meet 1e-12 are exactly the ones that do
//! not cancel: the stripper balances, whose inlet and outlet are unrelated in
//! magnitude, and the twelve valve lags, which are a single subtraction of two
//! independent numbers.
//!
//! So `PLAN.org`'s "rel err < 1e-12" cannot be read as relative-to-result for
//! a balance equation. Choosing what it *should* be measured against changes
//! what Tier 2 means, so it is a decision rather than a tolerance to nudge,
//! and per `CLAUDE.md` this file does not make it. What it asserts instead:
//!
//! - **0 ULP under `libm-system`, on all fifty.** The correctness claim, and
//!   stronger than any tolerance.
//! - **1e-12 relative on the twenty-two non-cancelling components** in the
//!   default build. Real, and it would catch a regression in them.
//! - The other twenty-eight are measured and reported, not gated.
//!
//! Nothing here is a relaxed threshold: the 1e-12 is still applied wherever it
//! is a meaningful question, and the components it is not applied to are named
//! individually rather than excluded by a widened margin.
//!
//! # Tripping states are counted, not compared
//!
//! `teprob.f:807-811` zeroes all fifty derivatives when the plant is down, and
//! the port reproduces that (delta D-007, and see `crate::balances` for why it
//! is the default). So a tripping state produces fifty zeros on both sides and
//! *any* port matches, correct or not.
//!
//! Those states are therefore reported separately rather than folded into the
//! headline number, which would otherwise be diluted by however many of them
//! the pool happens to contain. What they do check is the freeze itself, and
//! B-0024a already checks the detector that drives it on every state.

#![cfg(feature = "oracle")]

use tepsim_core::state::N_STATES;
use tepsim_core::{Inputs, Plant, SimTime, State, math, plant, vessels};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::Comparison;
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, adversarial};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2".
const TIER2_TOLERANCE: f64 = 1e-12;

/// Put a `Plant` into exactly the condition a scenario describes.
///
/// The walk-driven inputs come from `TESUB8` after the evaluation, as in every
/// earlier file; the seeds and the latched valve commands come out of the
/// scenario's own `COMMON`, which is what the Fortran will use.
fn configure(oracle: &mut Oracle, scenario: &Scenario) -> (Plant, Inputs) {
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

    let walks = plant::WalkInputs {
        feed: tepsim_core::streams::FeedConditions {
            ac_feed_light: [a, b, 1.0 - a - b],
            d_feed_celsius: oracle.tesub8(3, t) + raw(3) * 5.0,
            ac_feed_celsius: oracle.tesub8(4, t),
        },
        flow: tepsim_core::flows::FlowDrift {
            steam_capacity: oracle.tesub8(9, t),
            reactor_outlet: oracle.tesub8(12, t),
        },
        heat: tepsim_core::heat::HeatDrift {
            reactor_coolant: oracle.tesub8(10, t),
            condenser_coolant: oracle.tesub8(11, t),
        },
        reaction: tepsim_core::kinetics::ReactionDrift {
            first: oracle.tesub8(7, t),
            second: oracle.tesub8(8, t),
        },
        // teprob.f:413-414.
        coolant: tepsim_core::balances::CoolantInlet {
            reactor: oracle.tesub8(5, t) + raw(4) * 5.0,
            condenser: oracle.tesub8(6, t) + raw(5) * 5.0,
        },
    };

    let mut p = Plant::new();
    p.set_walk_inputs(walks);
    p.set_seeds(vessels::TemperatureSeeds {
        reactor: scenario.common.tcr,
        separator: scenario.common.tcs,
        stripper: scenario.common.tcc,
        mixing: scenario.common.tcv,
    });
    // `VCV` after the call is what `teprob.f:805` used, because the latch at
    // 799-804 runs immediately before it and nothing else writes `VCV`.
    p.set_valve_command(scenario.common.vcv);

    let inputs = Inputs {
        manipulated: scenario.manipulated,
        disturbances: core::array::from_fn(|i| f64::from(scenario.disturbances[i])),
    };
    (p, inputs)
}

/// The components whose balance does not catastrophically cancel, and which
/// therefore meet 1e-12 relative under the vendored libm.
///
/// Measured, not chosen: `per_component_relative_error` prints the full table,
/// and these are the complement of what it reports as over. They are the nine
/// stripper balances, the condenser wall, and the twelve valve lags.
const NON_CANCELLING: [usize; 22] = [
    19, 20, 21, 22, 23, 24, 25, 26, 27, // stripper components and energy
    37, // reactor coolant wall
    39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, // valve lags
];

/// One comparison of all fifty derivatives, split by whether the plant tripped.
struct Split {
    healthy: Comparison<Case>,
    gated: Comparison<Case>,
    tripped: Comparison<Case>,
    healthy_states: usize,
    tripped_states: usize,
    skipped: usize,
}

impl Split {
    fn new() -> Self {
        Self {
            healthy: Comparison::new("YP(1..50), plant running (reported)"),
            gated: Comparison::new("YP, the 22 non-cancelling components (gated)"),
            tripped: Comparison::new("YP(1..50), plant frozen"),
            healthy_states: 0,
            tripped_states: 0,
            skipped: 0,
        }
    }

    fn observe(&mut self, oracle: &mut Oracle, scenario: &Scenario, pool: Pool, index: usize) {
        let snapshot = scenario.force(oracle);
        let (p, u) = configure(oracle, scenario);
        let state = State::from_flat(&scenario.state);
        let Ok((derivative, signals)) = p.derivatives(SimTime(scenario.time), &state, &u) else {
            self.skipped += 1;
            return;
        };
        assert_eq!(
            signals.shutdown.is_tripped(),
            snapshot.tripped,
            "the port and the Fortran disagree about the trip at {pool}#{index}"
        );

        let tripped = snapshot.tripped;
        if tripped {
            self.tripped_states += 1;
        } else {
            self.healthy_states += 1;
        }
        for (slot, (ours, theirs)) in derivative
            .to_flat()
            .iter()
            .zip(snapshot.derivative.iter())
            .enumerate()
        {
            assert!(
                ours.is_finite(),
                "YP({}) is {ours} at {pool}#{index}",
                slot + 1
            );
            let case = Case {
                pool,
                index,
                component: slot + 1,
            };
            if tripped {
                self.tripped.observe(case, *ours, *theirs);
            } else {
                self.healthy.observe(case, *ours, *theirs);
                if NON_CANCELLING.contains(&(slot + 1)) {
                    self.gated.observe(case, *ours, *theirs);
                }
            }
        }
    }
}

fn sweep(oracle: &mut Oracle, steps: usize, perturbations: usize, seed: u64) -> Split {
    let pools = Pools::collect(oracle, steps, DT);
    let mut split = Split::new();

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        split.observe(oracle, &scenario, Pool::Nominal, index);
    }
    let mut sampler = tepsim_oracle::tier1::Sampler::new(seed);
    for index in 0..perturbations {
        let scenario = pools.perturbed_case(index, &mut sampler);
        split.observe(oracle, &scenario, Pool::Perturbed, index);
    }
    let (boundaries, missed) = adversarial::build(oracle, &pools.nominal_case(0));
    assert!(missed.is_empty(), "the catalogue regressed: {missed:?}");
    for (index, boundary) in boundaries.iter().enumerate() {
        split.observe(oracle, &boundary.scenario, Pool::Adversarial, index);
    }
    split
}

#[test]
fn all_fifty_derivatives_match_the_fortran_over_all_three_pools() {
    let mut oracle = Oracle::lock();
    let split = sweep(&mut oracle, 400, 2_000, 0x7E2_0025);

    println!(
        "transcendentals come from the {} libm",
        if math::USES_SYSTEM_LIBM {
            "platform"
        } else {
            "vendored"
        }
    );
    println!("{}", split.healthy);
    println!("{}", split.gated);
    println!("{}", split.tripped);
    println!(
        "{} states running, {} frozen, {} skipped",
        split.healthy_states, split.tripped_states, split.skipped
    );

    assert_eq!(split.skipped, 0, "states failed to converge in the port");
    assert!(
        split.healthy_states > 2_000,
        "only {} states are evidence; the rest are frozen and prove nothing",
        split.healthy_states
    );
    assert!(
        split.tripped_states > 0,
        "no state froze, so the freeze is implemented and never checked"
    );

    // The twenty-two components where relative error is a meaningful question.
    // See the module documentation, and B-0026a for the other twenty-eight.
    split.gated.assert_within(TIER2_TOLERANCE);

    // A frozen plant is fifty zeros on both sides. That is not evidence about
    // the model, but it *is* evidence that the freeze fires in the same places.
    assert_eq!(split.tripped.max_ulp(), 0, "the freeze disagrees");
}

/// With the transcendentals taken out of the comparison, all fifty must be
/// bit-identical. See `tier2_equilibrium.rs`.
#[test]
#[cfg(feature = "libm-system")]
fn the_whole_right_hand_side_is_bit_identical_once_exp_and_pow_agree() {
    let mut oracle = Oracle::lock();
    let split = sweep(&mut oracle, 200, 500, 8);
    println!("{}", split.healthy);
    assert_eq!(
        split.healthy.max_ulp(),
        0,
        "the assembled derivative is not bit-identical under the platform \
         libm, so the difference is in the algebra somewhere between \
         `teprob.f:407` and `811` and not in `exp` or `pow`"
    );
}

/// The Class C freeze is off when the fix is on, and the difference is real.
///
/// The oracle cannot check this: the fix is a deliberate divergence. What can
/// be checked is that it changes exactly the states it should and no others.
#[test]
fn the_quirk_fix_changes_only_the_frozen_states() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);
    let (boundaries, _) = adversarial::build(&mut oracle, &pools.nominal_case(0));

    let mut changed = 0;
    let mut unchanged = 0;
    for boundary in &boundaries {
        let scenario = &boundary.scenario;
        let _ = scenario.force(&mut oracle);
        let (mut faithful, u) = configure(&mut oracle, scenario);
        let mut fixed = faithful.clone();
        fixed.quirks.trip_ends_the_run = true;

        let state = State::from_flat(&scenario.state);
        faithful.quirks.trip_ends_the_run = false;
        let Ok((a, signals)) = faithful.derivatives(SimTime(scenario.time), &state, &u) else {
            continue;
        };
        let Ok((b, _)) = fixed.derivatives(SimTime(scenario.time), &state, &u) else {
            continue;
        };

        let differs = a.to_flat() != b.to_flat();
        if signals.shutdown.is_tripped() {
            assert!(
                differs,
                "{} trips, so the fix must change its derivative",
                boundary.target.name
            );
            changed += 1;
        } else {
            assert!(
                !differs,
                "{} does not trip, so the fix must change nothing",
                boundary.target.name
            );
            unchanged += 1;
        }
    }
    println!("the fix changes {changed} tripping boundaries and leaves {unchanged} alone");
    assert!(changed > 0 && unchanged > 0, "the test saw only one kind");
}

/// Per-component relative error, to find which slots cancel.
///
/// Diagnostic, not a gate. Printed so the numbers are in the run log.
#[test]
fn per_component_relative_error() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 400, DT);
    let mut worst = [0.0_f64; N_STATES];
    let mut scale = [0.0_f64; N_STATES];

    let mut record = |oracle: &mut Oracle,
                      scenario: &Scenario,
                      worst: &mut [f64; N_STATES],
                      scale: &mut [f64; N_STATES]| {
        let snapshot = scenario.force(oracle);
        if snapshot.tripped {
            return;
        }
        let (p, u) = configure(oracle, scenario);
        let state = State::from_flat(&scenario.state);
        let Ok((derivative, _)) = p.derivatives(SimTime(scenario.time), &state, &u) else {
            return;
        };
        for (slot, (ours, theirs)) in derivative
            .to_flat()
            .iter()
            .zip(snapshot.derivative.iter())
            .enumerate()
        {
            if *theirs != 0.0 {
                let rel = (ours - theirs).abs() / theirs.abs();
                if rel > worst[slot] {
                    worst[slot] = rel;
                }
            }
            scale[slot] = scale[slot].max(theirs.abs());
        }
    };

    for index in 0..pools.trajectory.len() {
        record(
            &mut oracle,
            &pools.nominal_case(index),
            &mut worst,
            &mut scale,
        );
    }
    let mut sampler = tepsim_oracle::tier1::Sampler::new(0x7E2_0025);
    for index in 0..2_000 {
        let scenario = pools.perturbed_case(index, &mut sampler);
        record(&mut oracle, &scenario, &mut worst, &mut scale);
    }

    println!("slot  max-rel-err   max|YP|");
    for slot in 0..N_STATES {
        if worst[slot] > 1e-12 {
            println!(
                "YP({:2})  {:.3e}   {:.4e}   OVER",
                slot + 1,
                worst[slot],
                scale[slot]
            );
        }
    }
    let over: Vec<usize> = (0..N_STATES)
        .filter(|i| worst[*i] > 1e-12)
        .map(|i| i + 1)
        .collect();
    println!("components over 1e-12 relative: {over:?}");
}

/// Every one of the fifty slots must be exercised: a derivative that is always
/// zero on both sides would match perfectly and mean nothing.
#[test]
fn every_derivative_slot_actually_moves_somewhere_in_the_pool() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);
    let mut moved = [false; N_STATES];

    let mut record = |oracle: &mut Oracle, scenario: &Scenario, moved: &mut [bool; N_STATES]| {
        let snapshot = scenario.force(oracle);
        if snapshot.tripped {
            return;
        }
        for (slot, value) in snapshot.derivative.iter().enumerate() {
            if value.to_bits() != 0.0_f64.to_bits() {
                moved[slot] = true;
            }
        }
    };

    for index in 0..pools.trajectory.len() {
        record(&mut oracle, &pools.nominal_case(index), &mut moved);
    }
    // The twelve valve lags are *still* on the nominal trajectory: `VCV`
    // equals `VPOS`, so `(VCV - VPOS)/VTAU` is exactly zero and stays that way
    // as long as nothing moves a command. Only the perturbed pool, which
    // scales the valve-position states away from their commands, exercises
    // them. Without this second loop the twelve comparisons are zero against
    // zero and mean nothing.
    let mut sampler = tepsim_oracle::tier1::Sampler::new(0x7E2_0025);
    for index in 0..500 {
        let scenario = pools.perturbed_case(index, &mut sampler);
        record(&mut oracle, &scenario, &mut moved);
    }

    let still: Vec<usize> = (0..N_STATES)
        .filter(|i| !moved[*i])
        .map(|i| i + 1)
        .collect();
    println!(
        "{} of {N_STATES} slots move somewhere in the pool",
        N_STATES - still.len()
    );
    assert!(
        still.is_empty(),
        "YP{still:?} never move anywhere in the pool, so those comparisons \
         are zero against zero and prove nothing"
    );
}
