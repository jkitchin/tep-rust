//! The closed-loop driver is linkable, and its controllers can be called.
//!
//! B-0036. Nothing is ported here; this establishes that the *next* six items
//! have an oracle to compare against.
//!
//! # Why `temain_mod.f` needed an edit at all
//!
//! It is an unnamed main program: no `PROGRAM` statement, so gfortran compiles
//! it to `MAIN__` and emits a `main` that calls it. Linked into a Rust test
//! binary that would collide with the harness's own entry point, and if it did
//! run it would open fifteen files under `~/` and simulate 48 hours on load.
//!
//! `instrument.rs` turns the main program into `SUBROUTINE TEMAIN_UNUSED`,
//! which nothing calls. The nineteen `CONTRLn` subroutines below it are
//! untouched and become linkable, which is the entire purpose.

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

use tepsim_oracle::{Ctrl1, Ctrl6, Ctrlall, Oracle};

/// A controller can be called, and it does what its source says.
///
/// `CONTRL1` is the simplest: proportional-only, on the D feed.
///
/// ```fortran
///       ERR1 = (SETPT(1) - XMEAS(2)) * 100. / 5811.
///       DXMV = GAIN1 * ( ( ERR1 - ERROLD1 ) )
///       XMV(1) = XMV(1) + DXMV
///       ERROLD1 = ERR1
/// ```
#[test]
fn a_controller_can_be_called_and_moves_its_valve() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();

    let mut setpoints = oracle.ctrlall();
    setpoints.setpt[0] = 4000.0;
    setpoints.deltat = 1.0 / 3600.0;
    oracle.set_ctrlall(&setpoints);
    oracle.set_ctrl1(&Ctrl1 {
        gain: 1.0,
        errold: 0.0,
    });

    // Put a known reading on XMEAS(2) and a known command on XMV(1).
    let mut measurements = oracle.measurements();
    measurements[1] = 3664.0;
    oracle.set_measurements(&measurements);
    let mut manipulated = oracle.manipulated();
    manipulated[0] = 63.0;
    oracle.set_manipulated(&manipulated);

    oracle.contrl1();

    let expected_error: f64 = (4000.0 - 3664.0) * 100.0 / 5811.0;
    let after = oracle.ctrl1();
    assert_eq!(
        after.errold.to_bits(),
        expected_error.to_bits(),
        "CONTRL1 stored {} as its error, not {expected_error}",
        after.errold
    );
    let moved = oracle.manipulated()[0];
    assert_eq!(
        moved.to_bits(),
        (63.0_f64 + expected_error).to_bits(),
        "XMV(1) went to {moved}"
    );
}

/// The velocity form is *incremental*: with the error unchanged, a
/// proportional-only loop makes no further move.
///
/// That is what "velocity form" means and it is easy to lose in a port, which
/// would give a controller that ramps its valve without limit at a constant
/// offset.
#[test]
fn a_proportional_loop_stops_moving_once_the_error_settles() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();

    let mut setpoints = oracle.ctrlall();
    setpoints.setpt[0] = 4000.0;
    oracle.set_ctrlall(&setpoints);
    oracle.set_ctrl1(&Ctrl1 {
        gain: 1.0,
        errold: 0.0,
    });
    let mut measurements = oracle.measurements();
    measurements[1] = 3664.0;
    oracle.set_measurements(&measurements);
    let mut manipulated = oracle.manipulated();
    manipulated[0] = 63.0;
    oracle.set_manipulated(&manipulated);

    oracle.contrl1();
    let after_first = oracle.manipulated()[0];
    oracle.contrl1();
    let after_second = oracle.manipulated()[0];
    oracle.contrl1();
    let after_third = oracle.manipulated()[0];

    assert_eq!(
        after_second.to_bits(),
        after_first.to_bits(),
        "the valve moved again on an unchanged error, so this is a positional \
         controller and not a velocity-form one"
    );
    assert_eq!(after_third.to_bits(), after_first.to_bits());
}

/// The purge override latches, and while it is latched the PI does not run.
///
/// `temain_mod.f:710-731`. This is the behaviour that makes `CONTRL6`
/// different from every other loop, and the reason B-0037b exists.
#[test]
fn the_purge_override_latches_and_suppresses_the_pi() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();

    let mut setpoints = oracle.ctrlall();
    setpoints.setpt[5] = 0.33712;
    oracle.set_ctrlall(&setpoints);
    oracle.set_ctrl6(&Ctrl6 {
        gain: 1.22,
        errold: 0.0,
    });
    oracle.set_flag6(0);

    // Separator pressure above 2950 kPa: the valve is thrown wide open and the
    // latch is set.
    let mut measurements = oracle.measurements();
    measurements[12] = 3000.0;
    measurements[9] = 0.33712;
    oracle.set_measurements(&measurements);
    oracle.contrl6();

    assert_eq!(oracle.flag6(), 1, "the override did not latch");
    assert_eq!(
        oracle.manipulated()[5].to_bits(),
        100.0_f64.to_bits(),
        "the purge valve was not thrown open"
    );
    assert_eq!(
        oracle.ctrl6().errold.to_bits(),
        0.0_f64.to_bits(),
        "ERROLD6 advanced while the override was latched, so the PI ran after \
         all and the branch structure has been read wrongly"
    );

    // Falling back through 2633.7 releases the latch and resets the loop.
    let mut measurements = oracle.measurements();
    measurements[12] = 2600.0;
    oracle.set_measurements(&measurements);
    oracle.contrl6();

    assert_eq!(oracle.flag6(), 0, "the override did not release");
    // `XMV(6)=40.060` and `SETPT(6)=0.33712` at `temain_mod.f:716-717` carry
    // no `D` suffix, so both are *single precision*: 40.060001373291016 and
    // 0.337119996547699. The thresholds are mixed -- 2950 and 2300 are exactly
    // representable, 2633.7 is not -- so the precision has to be read off each
    // line, exactly as it did throughout `teprob.f`.
    assert_eq!(
        oracle.manipulated()[5].to_bits(),
        f64::from(40.060_f32).to_bits(),
        "the valve was not reset to its nominal position"
    );
    assert_eq!(
        oracle.ctrlall().setpt[5].to_bits(),
        f64::from(0.33712_f32).to_bits(),
        "the setpoint was not reset with it"
    );
}

/// The low-pressure branch latches the other way.
#[test]
fn the_purge_override_also_latches_shut() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();
    oracle.set_ctrl6(&Ctrl6 {
        gain: 1.22,
        errold: 0.0,
    });
    oracle.set_flag6(0);

    let mut measurements = oracle.measurements();
    measurements[12] = 2000.0;
    oracle.set_measurements(&measurements);
    oracle.contrl6();

    assert_eq!(oracle.flag6(), 2, "the low-pressure branch did not latch");
    assert_eq!(oracle.manipulated()[5].to_bits(), 0.0_f64.to_bits());

    // And rising back through 2633.7 releases it.
    let mut measurements = oracle.measurements();
    measurements[12] = 2700.0;
    oracle.set_measurements(&measurements);
    oracle.contrl6();
    assert_eq!(oracle.flag6(), 0);
    assert_eq!(
        oracle.manipulated()[5].to_bits(),
        f64::from(40.060_f32).to_bits()
    );
}

/// A cascade controller writes a *setpoint*, not a valve.
///
/// `CONTRL13` drives `SETPT(3)`, which is `CONTRL3`'s. This is the structure
/// B-0035 found and that `PLAN.org`'s "nineteen near-identical loops" hides.
#[test]
fn a_cascade_controller_moves_a_setpoint_and_no_valve() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();

    let mut setpoints = oracle.ctrlall();
    setpoints.setpt[12] = 32.188;
    setpoints.setpt[2] = 0.25052;
    setpoints.deltat = 1.0 / 3600.0;
    oracle.set_ctrlall(&setpoints);
    oracle.set_ctrl13(&tepsim_oracle::Ctrl13 {
        gain: 0.0018,
        taui: 5000.0 / 3600.0,
        errold: 0.0,
    });
    let mut measurements = oracle.measurements();
    measurements[22] = 30.0;
    oracle.set_measurements(&measurements);
    let valves_before = oracle.manipulated();

    oracle.contrl13();

    assert_eq!(
        oracle.manipulated(),
        valves_before,
        "CONTRL13 moved a valve; it should only move SETPT(3)"
    );
    assert_ne!(
        oracle.ctrlall().setpt[2].to_bits(),
        0.25052_f64.to_bits(),
        "CONTRL13 did not move SETPT(3), which is the loop it cascades onto"
    );
}

/// All twenty controllers are callable, including the one the driver never
/// calls.
///
/// `CONTRL22` is dead code in the driver and perfectly live as a subroutine.
/// Binding it costs nothing and means B-0037's differential can cover it,
/// which is worth doing precisely because nothing else exercises it.
#[test]
fn every_controller_including_the_dead_one_is_callable() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();
    let mut setpoints = oracle.ctrlall();
    setpoints.deltat = 1.0 / 3600.0;
    oracle.set_ctrlall(&setpoints);

    oracle.contrl1();
    oracle.contrl2();
    oracle.contrl3();
    oracle.contrl4();
    oracle.contrl5();
    oracle.contrl6();
    oracle.contrl7();
    oracle.contrl8();
    oracle.contrl9();
    oracle.contrl10();
    oracle.contrl11();
    oracle.contrl13();
    oracle.contrl14();
    oracle.contrl15();
    oracle.contrl16();
    oracle.contrl17();
    oracle.contrl18();
    oracle.contrl19();
    oracle.contrl20();
    oracle.contrl22();

    // Nothing to assert beyond "none of that crashed or hung". The point is
    // that the symbols resolve and the driver's `main` is gone: if
    // `TEMAIN_UNUSED` were still a main program this binary would not link,
    // and if it ran it would open fifteen files and simulate 48 hours.
    assert!(oracle.ctrlall().deltat > 0.0);
}

/// The setpoint array has twenty slots and the dead controller owns one.
#[test]
fn the_setpoint_array_has_room_for_the_dead_controller() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();
    let all: Ctrlall = oracle.ctrlall();
    assert_eq!(all.setpt.len(), 20);
    // `SETPT(12)` is `CONTRL22`'s, which the driver initialises to 2633.7 --
    // exactly the purge override's reset threshold. See B-0035.
    let mut with_value = all;
    with_value.setpt[11] = 2633.7;
    oracle.set_ctrlall(&with_value);
    assert_eq!(oracle.ctrlall().setpt[11].to_bits(), 2633.7_f64.to_bits());
}
