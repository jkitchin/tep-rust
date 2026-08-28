//! Physics invariants: the only tests in the ladder that can catch an error
//! the port faithfully inherited.
//!
//! Every other tier compares the port against the Fortran, so an error *in the
//! Fortran* is invisible to all of them: the port reproduces it, the
//! differential passes, and the number is wrong. These do not consult the
//! Fortran's answer at all. They consult conservation.
//!
//! Each runs against **both** implementations, and a failure is attributed
//! explicitly. A failure on the port alone is a porting bug. A failure on both
//! is a finding about the original, which is a Class B delta and belongs in
//! `book/src/deltas.md`, not in a fix.
//!
//! # The invariants
//!
//! **I-1, reaction mass conservation.** `teprob.f:521-527` writes the net
//! molar production of each component from the four reaction rates:
//!
//! ```text
//! CRXR(1) = -RR1 - RR2 - RR3        CRXR(5) = -RR2 - RR3
//! CRXR(3) = -RR1 - RR2              CRXR(6) =  RR3 + RR4
//! CRXR(4) = -RR1 - 1.5 RR4          CRXR(7) =  RR1
//!                                   CRXR(8) =  RR2
//! ```
//!
//! With the molecular weights at `teprob.f:1035` (A 2.0, B 25.4, C 28.0,
//! D 32.0, E 46.0, F 48.0, G 62.0, H 76.0) each reaction balances exactly:
//!
//! ```text
//! A + C + D -> G     2 + 28 + 32 = 62
//! A + C + E -> H     2 + 28 + 46 = 76
//! A + E     -> F     2 + 46      = 48
//! 3 D       -> 2 F   3 * 32      = 2 * 48
//! ```
//!
//! so the invariant is
//!
//! ```text
//! sum_i XMW(i) * CRXR(i) = 0     exactly, in exact arithmetic
//! ```
//!
//! **Tolerance.** Relative to the scale of the terms, not to the result: the
//! result is zero by construction and the terms are not. Same convention as
//! Tier 2, decided 2026-08-27. Gate 1e-14 of `max_i |XMW(i) * CRXR(i)|`.
//!
//! **I-2, total mass balance over the whole plant.** The four vessels hold
//! `YY(1..8)`, `YY(10..17)`, `YY(19..26)` and `YY(28..35)` moles of each
//! component (`teprob.f:762-789`). Streams 1 to 4 enter and streams 10 and 11
//! leave (`crate` `Stream`); every other stream is internal and must cancel
//! between two vessels. Since I-1 makes the reaction mass-neutral,
//!
//! ```text
//! sum_i XMW(i) * (YP(i) + YP(i+9) + YP(i+18) + YP(i+27))
//!     = FTM(1) XMWS(1) + FTM(2) XMWS(2) + FTM(3) XMWS(3) + FTM(4) XMWS(4)
//!     - FTM(10) XMWS(10) - FTM(13) XMWS(13)
//! ```
//!
//! This is what catches a dropped term in a balance equation. `xtask
//! provenance` cannot: a term that is never evaluated is also never claimed,
//! and a term evaluated into the wrong equation is claimed correctly.
//!
//! **I-4, per-component molar balance.** Strictly stronger than I-2, and the
//! reason it exists is worth recording: I-2 is the `XMW`-weighted *sum* of
//! I-4, and because I-1 makes the reaction mass-neutral, the reaction term
//! cancels out of that sum entirely. Deleting `CRXR` from the reactor balance
//! therefore passes I-2 unchanged. Mutation testing found that; the invariant
//! that catches it is the unweighted one:
//!
//! ```text
//! YP(i) + YP(i+9) + YP(i+18) + YP(i+27)
//!     = FCM(i,1) + FCM(i,2) + FCM(i,3) + FCM(i,4)
//!     - FCM(i,10) - FCM(i,13) + CRXR(i)
//! ```
//!
//! for each of the eight components.
//!
//! **I-3, the inert has nowhere to go but the purge.** B takes part in no
//! reaction, so `CRXR(2)` is zero (`teprob.f:521-527` never assigns it; see
//! delta D-003) and B's only exit is stream 10. Its balance therefore has no
//! reaction term at all.

// Differential tests: exact comparisons are the property under test, and
// arithmetic is transcribed from `teprob.f` so a reader can check it against
// the listing line by line. Rearranging either would defeat the point.
#![allow(
    clippy::float_cmp,
    reason = "bit equality against the Fortran is the property under test"
)]
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions are transcribed to be checkable against teprob.f"
)]
#![cfg(feature = "oracle")]

use tepsim_core::component::Component;
use tepsim_core::constants::{NOMINAL_STATE, XMW};
use tepsim_core::{
    Inputs, Plant, SimTime, State, Stream, equilibrium, flows, kinetics, streams, stripper, unpack,
};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier2::{Pools, Scenario};

const DT: f64 = 1.0 / 3600.0;

/// How many trajectory states the pool tests sweep.
const POOL_STEPS: usize = 200;

/// I-1's gate, relative to the scale of the terms.
const REACTION_MASS_TOLERANCE: f64 = 1e-14;

/// I-2's gate, relative to the total mass throughput.
///
/// Both implementations measure 2.1e-16 to 2.3e-16, which is one to two ulp of
/// a 32,000 lb/h throughput: the balance closes as exactly as `f64` allows.
/// This asks for 1e-13, three orders of headroom, rather than the 1e-11 first
/// written here. A term dropped from a balance equation would be an *entire
/// stream*, of order 1e-5 relative at the smallest, so there is no risk in
/// asking for what the code actually delivers.
const MASS_BALANCE_TOLERANCE: f64 = 1e-13;

/// The residual of I-1, and the scale it should be judged against.
fn reaction_residual(crxr: &[f64; 8]) -> (f64, f64) {
    let mut total = 0.0;
    let mut scale = 0.0_f64;
    for component in Component::ALL {
        let term = XMW[component] * crxr[component.index()];
        total += term;
        scale = scale.max(term.abs());
    }
    (total.abs(), scale)
}

/// Rebuild the port's reaction rates from a scenario's state.
///
/// The pieces are rebuilt here rather than reached for through `Plant`,
/// following the Tier 2 tests: every input the invariant depends on is then
/// visible at the call site.
fn port_production(scenario: &Scenario) -> Option<[f64; 8]> {
    let state = State::from_flat(&scenario.state);
    let unpacked = unpack(&state, Default::default()).ok()?;
    let eq = equilibrium(&unpacked);
    let kin = kinetics(&eq.reactor, unpacked.reactor.kelvin(), Default::default());
    Some(*kin.production.as_array())
}

/// The mean molecular weight of a stream, from its composition.
///
/// Computed here rather than read from `XMWS`, because `XMWS` is only filled
/// for streams 1, 2, 6, 8, 9 and 10 (`teprob.f:529-548`). Streams 3, 4, 11 and
/// 13 keep whatever was in the slot, and two of those are exactly the streams
/// this invariant needs. Reading `XMWS(13)` gives a mass balance that fails by
/// 95%, identically on both implementations, which is how this was found.
fn mean_molecular_weight(composition: &[f64; 8]) -> f64 {
    let mut total = 0.0;
    for component in Component::ALL {
        total += composition[component.index()] * XMW[component];
    }
    total
}

/// The residual of I-2, and the throughput it should be judged against.
fn mass_balance_residual(
    derivative: &[f64; 50],
    ftm: &[f64; 13],
    composition: &[[f64; 8]; 13],
) -> (f64, f64) {
    // Accumulation, summed over the four vessels.
    let mut accumulation = 0.0;
    for component in Component::ALL {
        let slot = component.index();
        let moles =
            derivative[slot] + derivative[slot + 9] + derivative[slot + 18] + derivative[slot + 27];
        accumulation += XMW[component] * moles;
    }

    // Streams 1..4 in, 10 and 13 out. One-based in the Fortran.
    let mass = |stream: usize| ftm[stream - 1] * mean_molecular_weight(&composition[stream - 1]);
    let inflow = mass(1) + mass(2) + mass(3) + mass(4);
    let outflow = mass(10) + mass(13);

    let residual = accumulation - (inflow - outflow);
    let scale = inflow.abs().max(outflow.abs()).max(accumulation.abs());
    (residual.abs(), scale)
}

/// I-1 on the oracle, over the whole Tier 2 sampling pool.
#[test]
fn the_fortran_conserves_mass_across_the_reaction() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, POOL_STEPS, DT);

    let mut worst = (0.0_f64, 0usize);
    let mut cases = 0;
    for index in 0..POOL_STEPS {
        let scenario = pools.nominal_case(index);
        let snapshot = scenario.force(&mut oracle);
        let (residual, scale) = reaction_residual(&snapshot.common.crxr);
        if scale == 0.0 {
            // No reaction at all: the invariant is 0 = 0 and says nothing.
            continue;
        }
        cases += 1;
        let relative = residual / scale;
        if relative > worst.0 {
            worst = (relative, index);
        }
    }

    println!(
        "I-1 on the Fortran: {cases} states, worst residual {:.3e} of the term \
         scale, at pool index {}",
        worst.0, worst.1
    );
    assert!(cases > 0, "no state had a non-zero reaction rate");
    assert!(
        worst.0 < REACTION_MASS_TOLERANCE,
        "the Fortran's reaction does not conserve mass: residual {:.3e} of the \
         term scale. That is a finding about teprob.f, not about the port, and \
         it belongs in book/src/deltas.md as a Class B delta.",
        worst.0
    );
}

/// I-1 on the port, over the same pool.
#[test]
fn the_port_conserves_mass_across_the_reaction() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, POOL_STEPS, DT);

    let mut worst = (0.0_f64, 0usize);
    let mut cases = 0;
    for index in 0..POOL_STEPS {
        let scenario = pools.nominal_case(index);
        // The port's own kinetics, from the same state.
        let Some(rates) = port_production(&scenario) else {
            continue;
        };
        let (residual, scale) = reaction_residual(&rates);
        if scale == 0.0 {
            continue;
        }
        cases += 1;
        let relative = residual / scale;
        if relative > worst.0 {
            worst = (relative, index);
        }
    }

    println!(
        "I-1 on the port: {cases} states, worst residual {:.3e} of the term \
         scale, at pool index {}",
        worst.0, worst.1
    );
    assert!(cases > 0, "no state had a non-zero reaction rate");
    assert!(
        worst.0 < REACTION_MASS_TOLERANCE,
        "the port's reaction does not conserve mass: residual {:.3e}",
        worst.0
    );
}

/// The four reactions balance exactly on paper, which is why I-1 is an equality
/// rather than a tolerance.
///
/// Asserted from the molecular weights rather than asserted about them, so that
/// a change to `XMW` fails here with the arithmetic on show.
#[test]
fn each_reaction_balances_on_paper() {
    let w = XMW.as_array();
    let (a, b, c, d, e, f, g, h) = (w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7]);

    assert_eq!(a + c + d, g, "A + C + D -> G");
    assert_eq!(a + c + e, h, "A + C + E -> H");
    assert_eq!(a + e, f, "A + E -> F");
    assert_eq!(3.0 * d, 2.0 * f, "3 D -> 2 F");
    // B is inert, and its weight enters no reaction. Recorded so that the
    // absence is deliberate rather than an omission.
    assert!(b > 0.0, "B still has a molecular weight: {b}");
}

/// I-2 on both implementations, from the nominal state.
#[test]
fn the_whole_plant_balances_mass() {
    let mut oracle = Oracle::lock();
    let (_, yy) = oracle.init_cold();
    oracle.set_disturbances(&[0; 20]);
    oracle.set_rng(tepsim_oracle::golden::SEED);

    let yp = oracle.derivatives(0.0, &yy);
    let common = oracle.teproc();
    let (residual, scale) = mass_balance_residual(&yp, &common.ftm, &common.xst);

    println!(
        "I-2 on the Fortran at the nominal state: residual {:.6e} lb/h against \
         a throughput of {scale:.1} lb/h, {:.3e} relative",
        residual,
        residual / scale
    );
    assert!(
        residual / scale < MASS_BALANCE_TOLERANCE,
        "the Fortran's plant does not balance mass: {residual:.6e} lb/h out of \
         {scale:.1} lb/h. A term is missing from a balance equation, or a \
         stream is counted twice."
    );
}

/// I-2 along a trajectory, not just at the nominal point.
///
/// A balance can close at steady state and fail away from it, because at steady
/// state the accumulation term is near zero and a missing term hides inside it.
#[test]
fn the_plant_balances_mass_away_from_steady_state() {
    let mut oracle = Oracle::lock();
    let (_, mut yy) = oracle.init_cold();
    oracle.set_disturbances(&[0; 20]);
    oracle.set_rng(tepsim_oracle::golden::SEED);

    // Open loop from the nominal valves, which drifts hard: B-0041 measured it
    // tripping in three hours. Two hours of that is plenty of "away".
    let mut worst = (0.0_f64, 0usize, 0.0_f64);
    let mut t = 0.0;
    for step in 1..=7_200 {
        let yp = oracle.derivatives(t, &yy);
        let common = oracle.teproc();
        let (residual, scale) = mass_balance_residual(&yp, &common.ftm, &common.xst);
        let accumulation: f64 = Component::ALL
            .into_iter()
            .map(|c| {
                let s = c.index();
                XMW[c] * (yp[s] + yp[s + 9] + yp[s + 18] + yp[s + 27])
            })
            .sum();
        if residual / scale > worst.0 {
            worst = (residual / scale, step, accumulation);
        }
        for (slot, rate) in yy.iter_mut().zip(yp) {
            *slot += DT * rate;
        }
        t += DT;
    }

    println!(
        "I-2 along 2 h open loop: worst {:.3e} relative at step {}, where the \
         accumulation term is {:.1} lb/h",
        worst.0, worst.1, worst.2
    );
    // Teeth: the residual has to be small *relative to a non-zero
    // accumulation*, or the test is only saying that zero equals zero. The
    // plant's mass inventory drifts slowly even while it is heading for a
    // pressure trip, so the bar is not that the accumulation is large but that
    // it is many orders above the residual.
    let residual_absolute = worst.0 * 32_000.0;
    println!(
        "  accumulation {:.4} lb/h against a residual of about {:.2e} lb/h: a \
         factor of {:.1e}",
        worst.2.abs(),
        residual_absolute,
        worst.2.abs() / residual_absolute
    );
    assert!(
        worst.2.abs() > 1e6 * residual_absolute,
        "the accumulation term reached only {:.4} lb/h, not far enough above \
         the {:.2e} lb/h residual for this test to distinguish a closed \
         balance from an empty one",
        worst.2.abs(),
        residual_absolute
    );
    assert!(
        worst.0 < MASS_BALANCE_TOLERANCE,
        "mass balance fails by {:.3e} relative at step {}",
        worst.0,
        worst.1
    );
}

/// I-2 on the port, from the same state.
#[test]
fn the_ported_plant_balances_mass() {
    let mut oracle = Oracle::lock();
    let (_, _) = oracle.init_cold();
    oracle.set_disturbances(&[0; 20]);
    let after_init = oracle.teproc();

    let mut plant = Plant::new();
    plant.set_rng(tepsim_oracle::golden::SEED);
    plant.set_seeds(tepsim_core::TemperatureSeeds {
        reactor: after_init.tcr,
        separator: after_init.tcs,
        stripper: after_init.tcc,
        mixing: after_init.tcv,
    });
    let state = State::from_flat(&NOMINAL_STATE);
    let inputs = Inputs {
        manipulated: core::array::from_fn(|i| NOMINAL_STATE[38 + i]),
        disturbances: [0.0; 20],
    };
    let time = SimTime(0.0);
    plant.advance_discrete(time, &inputs);
    let (derivative, signals) = plant.derivatives(time, &state, &inputs).expect("converges");
    let _ = signals;

    // The port's own stream table and flows, rebuilt from the same state and
    // *the same walk inputs the derivative used*. Rebuilding with the defaults
    // instead leaves the residual at 6.9e-9 rather than 1e-16, because the A/C
    // feed composition then differs from the one the balance saw.
    let walks = plant.walk_inputs();
    let unpacked = unpack(
        &state,
        tepsim_core::TemperatureSeeds {
            reactor: after_init.tcr,
            separator: after_init.tcs,
            stripper: after_init.tcc,
            mixing: after_init.tcv,
        },
    )
    .expect("converges");
    let eq = equilibrium(&unpacked);
    let mut table = streams(&unpacked, &eq, &walks.feed);
    let mut flow = flows(&state, &unpacked, &eq, &table, &[0.0; 20], walks.flow);
    let _ = stripper(&mut table, &mut flow, unpacked.stripper.celsius);

    let ftm: [f64; 13] = core::array::from_fn(|i| flow.molar[Stream::ALL[i]]);
    let composition: [[f64; 8]; 13] = core::array::from_fn(|i| {
        core::array::from_fn(|c| table.composition[Stream::ALL[i]][Component::ALL[c]])
    });
    let flat = derivative.to_flat();
    let (residual, scale) = mass_balance_residual(&flat, &ftm, &composition);

    println!(
        "I-2 on the port at the nominal state: residual {residual:.6e} lb/h \
         against {scale:.1} lb/h, {:.3e} relative",
        residual / scale
    );
    assert!(
        residual / scale < MASS_BALANCE_TOLERANCE,
        "the port's plant does not balance mass: {residual:.6e} lb/h"
    );
}

/// The residual of I-4 for one component, and the scale to judge it against.
fn component_balance_residual(
    component: Component,
    derivative: &[f64; 50],
    fcm: &[[f64; 8]; 13],
    crxr: &[f64; 8],
) -> (f64, f64) {
    let slot = component.index();
    let accumulation =
        derivative[slot] + derivative[slot + 9] + derivative[slot + 18] + derivative[slot + 27];

    let flow = |stream: usize| fcm[stream - 1][slot];
    let inflow = flow(1) + flow(2) + flow(3) + flow(4);
    let outflow = flow(10) + flow(13);
    let reaction = crxr[slot];

    let residual = accumulation - (inflow - outflow + reaction);
    let scale = inflow
        .abs()
        .max(outflow.abs())
        .max(reaction.abs())
        .max(accumulation.abs());
    (residual.abs(), scale)
}

/// I-4 on both implementations, from the nominal state and along a trajectory.
///
/// This is the invariant that sees the reaction. I-2 cannot: it is I-4's
/// mass-weighted sum, and the reaction is mass-neutral by I-1, so the reaction
/// term vanishes from it. Deleting `CRXR` from the reactor's balance leaves
/// I-2 passing and fails this.
#[test]
fn every_component_balances_on_moles() {
    let mut oracle = Oracle::lock();
    let (_, mut yy) = oracle.init_cold();
    oracle.set_disturbances(&[0; 20]);
    oracle.set_rng(tepsim_oracle::golden::SEED);

    let mut worst = (0.0_f64, Component::A, 0usize);
    let mut biggest_reaction = 0.0_f64;
    let mut t = 0.0;
    for step in 1..=3_600 {
        let yp = oracle.derivatives(t, &yy);
        let common = oracle.teproc();
        for component in Component::ALL {
            let (residual, scale) =
                component_balance_residual(component, &yp, &common.fcm, &common.crxr);
            if scale == 0.0 {
                continue;
            }
            if residual / scale > worst.0 {
                worst = (residual / scale, component, step);
            }
            biggest_reaction = biggest_reaction.max(common.crxr[component.index()].abs());
        }
        for (slot, rate) in yy.iter_mut().zip(yp) {
            *slot += DT * rate;
        }
        t += DT;
    }

    println!(
        "I-4 on the Fortran over 1 h: worst {:.3e} relative, on {:?} at step          {}; largest reaction term seen {:.1} lbmol/h",
        worst.0, worst.1, worst.2, biggest_reaction
    );
    assert!(
        biggest_reaction > 1.0,
        "the largest reaction term over the whole run was {biggest_reaction:.3e}          lbmol/h, so this test never actually exercised the reaction and could          not distinguish a balance that includes it from one that does not"
    );
    assert!(
        worst.0 < MASS_BALANCE_TOLERANCE,
        "component {:?} does not balance: {:.3e} relative at step {}",
        worst.1,
        worst.0,
        worst.2
    );
}

/// I-4 on the port, at the nominal state.
#[test]
fn every_component_balances_on_moles_in_the_port() {
    let mut oracle = Oracle::lock();
    let (_, _) = oracle.init_cold();
    oracle.set_disturbances(&[0; 20]);
    let after_init = oracle.teproc();
    let seeds = tepsim_core::TemperatureSeeds {
        reactor: after_init.tcr,
        separator: after_init.tcs,
        stripper: after_init.tcc,
        mixing: after_init.tcv,
    };

    let mut plant = Plant::new();
    plant.set_rng(tepsim_oracle::golden::SEED);
    plant.set_seeds(seeds);
    let state = State::from_flat(&NOMINAL_STATE);
    let inputs = Inputs {
        manipulated: core::array::from_fn(|i| NOMINAL_STATE[38 + i]),
        disturbances: [0.0; 20],
    };
    let time = SimTime(0.0);
    plant.advance_discrete(time, &inputs);
    let (derivative, _) = plant.derivatives(time, &state, &inputs).expect("converges");

    let walks = plant.walk_inputs();
    let unpacked = unpack(&state, seeds).expect("converges");
    let eq = equilibrium(&unpacked);
    let mut table = streams(&unpacked, &eq, &walks.feed);
    let mut flow = flows(&state, &unpacked, &eq, &table, &[0.0; 20], walks.flow);
    let _ = stripper(&mut table, &mut flow, unpacked.stripper.celsius);
    let kin = kinetics(&eq.reactor, unpacked.reactor.kelvin(), walks.reaction);

    let fcm: [[f64; 8]; 13] = core::array::from_fn(|i| {
        core::array::from_fn(|c| flow.component[Stream::ALL[i]][Component::ALL[c]])
    });
    let crxr = *kin.production.as_array();
    let flat = derivative.to_flat();

    let mut worst = (0.0_f64, Component::A);
    for component in Component::ALL {
        let (residual, scale) = component_balance_residual(component, &flat, &fcm, &crxr);
        if scale == 0.0 {
            continue;
        }
        if residual / scale > worst.0 {
            worst = (residual / scale, component);
        }
    }
    println!(
        "I-4 on the port at the nominal state: worst {:.3e} relative, on {:?}",
        worst.0, worst.1
    );
    assert!(
        worst.0 < MASS_BALANCE_TOLERANCE,
        "component {:?} does not balance in the port: {:.3e} relative",
        worst.1,
        worst.0
    );
}

/// I-3: the inert takes part in no reaction.
#[test]
fn the_inert_has_no_reaction_term() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 100, DT);
    let mut seen = 0;
    for index in 0..100 {
        let scenario = pools.nominal_case(index);
        let snapshot = scenario.force(&mut oracle);
        assert_eq!(
            snapshot.common.crxr[1], 0.0,
            "CRXR(2) is {} but B is inert",
            snapshot.common.crxr[1]
        );
        if snapshot.common.crxr[0] != 0.0 {
            seen += 1;
        }
    }
    println!("I-3: CRXR(2) is zero across {seen} reacting states");
    assert!(seen > 0, "no state was reacting, so this proves nothing");
}
