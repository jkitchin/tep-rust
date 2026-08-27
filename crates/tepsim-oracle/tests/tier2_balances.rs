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
//! # The gate is 1e-12 of the scale of the terms
//!
//! Not of the derivative. B-0025 measured that relative-to-result cannot be
//! met and that the port is not the reason: under the platform libm the whole
//! right-hand side is *bit-identical* to the Fortran, all fifty components,
//! every state.
//!
//! The reason is cancellation. A balance is
//! \(\sum_\text{in} - \sum_\text{out}\), and near steady state those sums
//! nearly agree: `YP(2)`, the inert's reactor balance, is a difference of two
//! flows around 660 whose result is a few parts in ten thousand of either. A
//! one-ULP difference in each term, which is all the vendored `exp` costs, is
//! 1e-16 of the *terms* and 1e-4 of the *result*.
//!
//! So each balance reports the magnitude of its largest term alongside its
//! value, which is its error budget, and the gate is the error over that.
//! `tepsim_core::balances::Balances::scale` is where it comes from; the
//! decision of 2026-08-27 in `BACKLOG.org` is why.
//!
//! This is not a weakened threshold. It is the same 1e-12, asked of the
//! quantity that can answer it, and it applies to all fifty components rather
//! than to the twenty-two that happened not to cancel.
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

/// One comparison of all fifty derivatives, split by whether the plant tripped.
struct Split {
    healthy: Comparison<Case>,
    gated: Comparison<Case>,
    tripped: Comparison<Case>,
    healthy_states: usize,
    tripped_states: usize,
    skipped: usize,
    /// Worst scaled error seen per component, for the acceptance table.
    worst_by_component: [f64; N_STATES],
    /// Worst error relative to the derivative itself, for the same table. Kept
    /// so the contrast between the two measures is visible rather than
    /// asserted.
    worst_relative_by_component: [f64; N_STATES],
}

impl Split {
    fn new() -> Self {
        Self {
            healthy: Comparison::new("YP(1..50), relative to the derivative (reported)"),
            gated: Comparison::new("YP(1..50), relative to the scale of the terms (the gate)"),
            tripped: Comparison::new("YP(1..50), plant frozen"),
            healthy_states: 0,
            tripped_states: 0,
            skipped: 0,
            worst_by_component: [0.0; N_STATES],
            worst_relative_by_component: [0.0; N_STATES],
        }
    }

    fn observe(&mut self, oracle: &mut Oracle, scenario: &Scenario, pool: Pool, index: usize) {
        let snapshot = scenario.force(oracle);
        let (p, u) = configure(oracle, scenario);
        let state = State::from_flat(&scenario.state);
        let Ok((derivative, scale, signals)) =
            p.derivatives_with_scale(SimTime(scenario.time), &state, &u)
        else {
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
        for (slot, ((ours, theirs), budget)) in derivative
            .to_flat()
            .iter()
            .zip(snapshot.derivative.iter())
            .zip(scale.to_flat())
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
                // The reported figure, relative to the derivative, kept so the
                // log can show what the old reading would have said.
                self.healthy.observe(case, *ours, *theirs);
                // The gate, relative to the scale of the terms.
                self.gated.observe_against(case, *ours, *theirs, budget);
                self.worst_by_component[slot] =
                    self.worst_by_component[slot].max(scaled_error(*ours, *theirs, budget));
                self.worst_relative_by_component[slot] = self.worst_relative_by_component[slot]
                    .max(scaled_error(*ours, *theirs, *theirs));
            }
        }
    }
}

/// The error of one comparison, relative to the balance's own scale.
fn scaled_error(ours: f64, theirs: f64, scale: f64) -> f64 {
    if scale == 0.0 {
        if ours.to_bits() == theirs.to_bits() {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        (ours - theirs).abs() / scale.abs()
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

    // All fifty, against the scale of their own terms.
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

/// The Tier 2 acceptance table: every component, both measures, and the worst
/// case named.
///
/// This is what B-0026 owes the log. It prints rather than asserts; the gate
/// itself is `all_fifty_derivatives_match_the_fortran_over_all_three_pools`.
#[test]
fn tier2_acceptance_table() {
    let mut oracle = Oracle::lock();
    let split = sweep(&mut oracle, 400, 2_000, 0x7E2_0025);

    println!("Tier 2 acceptance, {} running states", split.healthy_states);
    println!("gate: error / scale-of-terms < {TIER2_TOLERANCE:e}");
    println!();
    println!("  YP   err/scale    err/value   ratio");
    let mut worst = (0.0_f64, 0_usize);
    for slot in 0..N_STATES {
        let scaled = split.worst_by_component[slot];
        let plain = split.worst_relative_by_component[slot];
        let ratio = if scaled > 0.0 { plain / scaled } else { 1.0 };
        println!(
            "  {:2}   {scaled:.3e}   {plain:.3e}   {ratio:8.0}x",
            slot + 1
        );
        if scaled > worst.0 {
            worst = (scaled, slot + 1);
        }
    }
    println!();
    println!(
        "worst component: YP({}) at {:.3e} of its own scale",
        worst.1, worst.0
    );
    println!("{}", split.gated);

    // The two measures differ by orders of magnitude on the cancelling
    // components and not at all on the rest. That contrast is the whole reason
    // the gate is written the way it is, so it is asserted rather than left to
    // be read off the table.
    let cancelling = (0..N_STATES)
        .filter(|s| {
            split.worst_by_component[*s] > 0.0
                && split.worst_relative_by_component[*s] / split.worst_by_component[*s] > 100.0
        })
        .count();
    println!("{cancelling} of {N_STATES} components cancel by more than 100x");
    if math::USES_SYSTEM_LIBM {
        // Every component is bit-identical here, so both measures are zero and
        // there is nothing to contrast. That is the point of this
        // configuration, not a weakness of the table.
        assert_eq!(
            split.gated.max_ulp(),
            0,
            "on the platform libm the derivative must be bit-identical"
        );
        return;
    }
    assert!(
        cancelling > 20,
        "only {cancelling} components show cancellation, so relative-to-result \
         would have been nearly as good and the decision of 2026-08-27 needs \
         revisiting"
    );
}

/// Every one of the fifty slots must be exercised: a derivative that is always
/// zero on both sides would match perfectly and mean nothing.
#[test]
fn every_derivative_slot_actually_moves_somewhere_in_the_pool() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);
    let mut moved = [false; N_STATES];

    let record = |oracle: &mut Oracle, scenario: &Scenario, moved: &mut [bool; N_STATES]| {
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
