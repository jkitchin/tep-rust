//! Pins the flat state layout against the Fortran's own `COMMON` variables.
//!
//! `tepsim-core`'s own tests prove the mapping is a bijection: fifty slots in,
//! fifty distinct slots out. That is necessary and nowhere near sufficient. A
//! layout that swapped the stripper and the mixing zone would round trip
//! perfectly and mis-attribute every Tier 2 component for the rest of Phase 2,
//! silently, because both are eight moles and an energy.
//!
//! So the semantics are checked, not just the shape. `TEINIT` ends by calling
//! `TEFUNC` (`teprob.f:1369`), which unpacks `YY` into the named `COMMON/TEPROC/`
//! variables at `teprob.f:417-440`. Reading both back and comparing them asks
//! the Fortran which slot is which instead of trusting a transcription.

#![cfg(feature = "oracle")]

use tepsim_core::{Component, State, state::N_STATES};
use tepsim_oracle::Oracle;

/// Every slot that `teprob.f:417-440` copies straight out of `YY`.
///
/// `UCVR(4..8)` and `UCLR(1..3)` are deliberately absent: the first is derived
/// from the vapour-liquid equilibrium later in the same call
/// (`teprob.f:500-501`) and the second is forced to zero (`teprob.f:420-421`),
/// so neither still holds what the state vector put there. The separator pair
/// is the same. Comparing them would fail for a reason that has nothing to do
/// with the layout.
#[test]
fn every_named_field_holds_what_the_fortran_unpacked_into_it() {
    let mut oracle = Oracle::lock();
    let (_, yy) = oracle.init();
    let common = oracle.teproc();
    let state = State::from_flat(&yy);

    let mut checked = 0;
    let mut wrong = Vec::new();
    let mut check = |what: &str, ours: f64, theirs: f64| {
        checked += 1;
        if ours.to_bits() != theirs.to_bits() {
            wrong.push(format!("  {what}: state {ours:?} vs COMMON {theirs:?}"));
        }
    };

    // Reactor: vapour A/B/C, liquid D-H, energy.
    for component in [Component::A, Component::B, Component::C] {
        let i = component.index();
        check(
            &format!("reactor.moles[{component:?}] / UCVR({})", i + 1),
            state.reactor.moles[component],
            common.ucvr[i],
        );
        check(
            &format!("separator.moles[{component:?}] / UCVS({})", i + 1),
            state.separator.moles[component],
            common.ucvs[i],
        );
    }
    for component in [
        Component::D,
        Component::E,
        Component::F,
        Component::G,
        Component::H,
    ] {
        let i = component.index();
        check(
            &format!("reactor.moles[{component:?}] / UCLR({})", i + 1),
            state.reactor.moles[component],
            common.uclr[i],
        );
        check(
            &format!("separator.moles[{component:?}] / UCLS({})", i + 1),
            state.separator.moles[component],
            common.ucls[i],
        );
    }

    // Stripper and mixing zone: all eight slots survive unpacking.
    for component in Component::ALL {
        let i = component.index();
        check(
            &format!("stripper.moles[{component:?}] / UCLC({})", i + 1),
            state.stripper.moles[component],
            common.uclc[i],
        );
        check(
            &format!("mixing.moles[{component:?}] / UCVV({})", i + 1),
            state.mixing.moles[component],
            common.ucvv[i],
        );
    }

    check("reactor.energy / ETR", state.reactor.energy, common.etr);
    check("separator.energy / ETS", state.separator.energy, common.ets);
    check("stripper.energy / ETC", state.stripper.energy, common.etc);
    check("mixing.energy / ETV", state.mixing.energy, common.etv);
    check("reactor_cw_out_c / TWR", state.reactor_cw_out_c, common.twr);
    check(
        "condenser_cw_out_c / TWS",
        state.condenser_cw_out_c,
        common.tws,
    );

    assert!(
        wrong.is_empty(),
        "{} of {checked} slots are mapped to the wrong COMMON variable:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
    assert_eq!(
        checked, 38,
        "expected 38 directly-unpacked slots: 32 moles, 4 energies, 2 temperatures"
    );
}

/// The four vessels are distinguishable at the nominal point, so the test above
/// could not have passed by accident.
///
/// If the reactor and the separator happened to hold identical inventories, a
/// swapped mapping would sail through. They do not, and this says so with
/// numbers rather than assuming it.
#[test]
fn the_four_vessels_are_actually_distinguishable_at_the_nominal_point() {
    let mut oracle = Oracle::lock();
    let (_, yy) = oracle.init();
    let state = State::from_flat(&yy);

    let energies = [
        ("reactor", state.reactor.energy),
        ("separator", state.separator.energy),
        ("stripper", state.stripper.energy),
        ("mixing", state.mixing.energy),
    ];
    println!("nominal vessel energies: {energies:?}");
    for (i, (name, a)) in energies.iter().enumerate() {
        for (other, b) in &energies[i + 1..] {
            let separation = (a - b).abs() / a.abs().max(b.abs());
            assert!(
                separation > 1e-3,
                "{name} and {other} differ by only {separation:e} relative, so \
                 a swapped layout would not be caught"
            );
        }
    }
}

/// The valve positions are the last twelve slots, and the Fortran's own valve
/// dynamics are what say so.
///
/// `VPOS` is not in `COMMON/TEPROC/`, so it cannot be read back directly.
/// `teprob.f:806` gives `YP(I+38) = (VCV(I) - VPOS(I)) / VTAU(I)`, and `VCV`
/// and `VTAU` both are. At the nominal steady state the valve derivative is
/// zero, which pins `VPOS(I)` to `VCV(I)`; perturbing one slot and reading the
/// derivative back pins *which* slot it was.
#[test]
fn the_valve_slots_are_pinned_by_the_valve_dynamics() {
    let mut oracle = Oracle::lock();
    let (time, nominal) = oracle.init();
    let common = oracle.teproc();

    // At the nominal point every valve is at its commanded position.
    let state = State::from_flat(&nominal);
    for (i, (position, commanded)) in state.valve_pos.iter().zip(common.vcv).enumerate() {
        assert_eq!(
            position.to_bits(),
            commanded.to_bits(),
            "valve {} sits at {position} but is commanded to {commanded}",
            i + 1
        );
    }

    // Now move one valve at a time and check the derivative that responds is
    // the one the layout says it should be.
    for valve in 0..12 {
        let mut perturbed = State::from_flat(&nominal);
        perturbed.valve_pos[valve] -= 1.0;
        let yp = oracle.derivatives(time, &perturbed.to_flat());

        for (slot, rate) in yp.iter().enumerate().take(N_STATES).skip(38) {
            let expected_valve = slot - 38;
            if expected_valve == valve {
                let predicted = 1.0 / common.vtau[valve];
                assert!(
                    (rate - predicted).abs() / predicted < 1e-12,
                    "moving valve {} should give YP({}) = 1/VTAU = {predicted}, \
                     got {rate}",
                    valve + 1,
                    slot + 1
                );
            } else {
                assert_eq!(
                    rate.to_bits(),
                    0.0_f64.to_bits(),
                    "moving valve {} disturbed YP({}), so the valve slots are \
                     not where the layout claims",
                    valve + 1,
                    slot + 1
                );
            }
        }
    }
}
