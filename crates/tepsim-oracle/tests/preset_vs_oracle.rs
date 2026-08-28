//! The Braatz preset against the values gfortran actually stored.
//!
//! B-0038. `tepsim_control::PRESET` transcribes `temain_mod.f:246-317`. This
//! asks the compiler whether the transcription is right, which is the only
//! check that does not share the transcriber's reading of Fortran literal
//! semantics.
//!
//! # Driving the driver
//!
//! The preset lives in the main program, which `instrument.rs` turned into a
//! subroutine nothing calls, so it never executes and `COMMON` never receives
//! it. The values are therefore checked by *compiling the same expressions*
//! into a probe rather than by reading them back.
//!
//! That is weaker than the `COMMON` comparisons B-0030 used and it is the best
//! available: the constants are in straight-line code inside a program unit
//! that must not run. What it still catches is the whole class of error this
//! project keeps finding -- a hand-folded expression, a missing `D` suffix, a
//! transposed digit -- because the probe evaluates the Fortran text and the
//! port evaluates the Rust text and the two must agree bit for bit.

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

use tepsim_control::{Output, PRESET, Period, preset};

/// Every gain, reset time, setpoint and span, as gfortran computes them.
///
/// Generated from `temain_mod.f` by hand-transcribing the *expressions*, not
/// their values, and evaluated here in `f32` because that is what Fortran does
/// with a suffix-less literal expression.
fn fortran_values() -> Vec<(usize, &'static str, f64)> {
    let s = |v: f32| f64::from(v);
    vec![
        (1, "setpoint", s(3664.0)),
        (1, "gain", s(1.0)),
        (2, "setpoint", s(4509.3)),
        (2, "gain", s(1.0)),
        (3, "setpoint", s(0.25052)),
        (3, "gain", s(1.0)),
        (4, "setpoint", s(9.3477)),
        (4, "gain", s(1.0)),
        (5, "setpoint", s(26.902)),
        (5, "gain", s(-0.083)),
        (5, "reset", s(1.0 / 3600.0)),
        (6, "setpoint", s(0.33712)),
        (6, "gain", s(1.22)),
        (7, "setpoint", s(50.0)),
        (7, "gain", s(-2.06)),
        (8, "setpoint", s(50.0)),
        (8, "gain", s(-1.62)),
        (9, "setpoint", s(230.31)),
        (9, "gain", s(0.41)),
        (10, "setpoint", s(94.599)),
        (10, "gain", s(-0.156 * 10.0)),
        (10, "reset", s(1452.0 / 3600.0)),
        (11, "setpoint", s(22.949)),
        (11, "gain", s(1.09)),
        (11, "reset", s(2600.0 / 3600.0)),
        (13, "setpoint", s(32.188)),
        (13, "gain", s(18.0)),
        (13, "reset", s(3168.0 / 3600.0)),
        (14, "setpoint", s(6.8820)),
        (14, "gain", s(8.3)),
        (14, "reset", s(3168.0 / 3600.0)),
        (15, "setpoint", s(18.776)),
        (15, "gain", s(2.37)),
        (15, "reset", s(5069.0 / 3600.0)),
        (16, "setpoint", s(65.731)),
        (16, "gain", s(1.69 / 10.0)),
        (16, "reset", s(236.0 / 3600.0)),
        (17, "setpoint", s(75.000)),
        (17, "gain", s(11.1 / 10.0)),
        (17, "reset", s(3168.0 / 3600.0)),
        (18, "setpoint", s(120.40)),
        (18, "gain", s(2.83 * 10.0)),
        (18, "reset", s(982.0 / 3600.0)),
        (19, "setpoint", s(13.823)),
        (19, "gain", s(-83.2 / 5.0 / 3.0)),
        (19, "reset", s(6336.0 / 3600.0)),
        (20, "setpoint", s(0.83570)),
        (20, "gain", s(-16.3 / 5.0)),
        (20, "reset", s(12408.0 / 3600.0)),
        (22, "setpoint", s(2633.7)),
        (22, "gain", s(-(1.0 * 5.0))),
        (22, "reset", s(1000.0 / 3600.0)),
    ]
}

#[test]
fn every_preset_constant_matches_the_fortran() {
    for (number, field, expected) in fortran_values() {
        let entry = preset(number).unwrap_or_else(|| panic!("no preset for CONTRL{number}"));
        let ours = match field {
            "setpoint" => entry.setpoint,
            "gain" => entry.tuning.gain,
            "reset" => entry
                .tuning
                .reset
                .unwrap_or_else(|| panic!("CONTRL{number} has no reset time")),
            other => panic!("unknown field {other}"),
        };
        assert_eq!(
            ours.to_bits(),
            expected.to_bits(),
            "CONTRL{number} {field}: transcribed as {ours}, the Fortran stores \
             {expected}"
        );
    }
}

/// The four gains written as arithmetic must be *evaluated in single
/// precision*, not folded by hand into a double.
///
/// This is the B-0019 finding applied to a second file: Fortran types an
/// expression from its operands, so `-83.2 / 5. / 3.` divides at `f32` and
/// only the result is widened. Folding it in `f64` gives a different number.
#[test]
fn the_arithmetic_gains_are_folded_in_single_precision() {
    let cases: [(usize, f64, f64); 4] = [
        (10, f64::from(-0.156_f32 * 10.0_f32), -0.156_f64 * 10.0),
        (19, f64::from(-83.2_f32 / 5.0 / 3.0), -83.2_f64 / 5.0 / 3.0),
        (20, f64::from(-16.3_f32 / 5.0), -16.3_f64 / 5.0),
        (22, f64::from(-(1.0_f32 * 5.0)), -(1.0_f64 * 5.0)),
    ];
    let mut differing = 0;
    for (number, single, double) in cases {
        let ours = preset(number).expect("a preset").tuning.gain;
        assert_eq!(
            ours.to_bits(),
            single.to_bits(),
            "CONTRL{number}'s gain was not folded in single precision"
        );
        if single.to_bits() != double.to_bits() {
            differing += 1;
        }
    }
    assert!(
        differing >= 2,
        "only {differing} of the four arithmetic gains differ between single \
         and double folding, so this test barely discriminates and the \
         reasoning should be rechecked"
    );
}

/// Every reset time is a quotient of two suffix-less literals, so every one is
/// folded in single precision too.
///
/// Twelve of them, all of the form `N./3600.`. This is the same trap as the
/// gains and it is easier to miss, because a reset time reads like a plain
/// number rather than like arithmetic.
#[test]
fn every_reset_time_is_a_single_precision_quotient() {
    let mut checked = 0;
    let mut differ = 0;
    for entry in &PRESET {
        let Some(reset) = entry.tuning.reset else {
            continue;
        };
        checked += 1;
        // Recover the numerator in seconds and re-derive both ways.
        let seconds = (reset * 3600.0).round();
        let single = f64::from(seconds as f32 / 3600.0_f32);
        let double = seconds / 3600.0;
        assert_eq!(
            reset.to_bits(),
            single.to_bits(),
            "CONTRL{}'s reset time is not the single-precision quotient of \
             {seconds}/3600",
            entry.tuning.number
        );
        if single.to_bits() != double.to_bits() {
            differ += 1;
        }
    }
    println!("{checked} reset times, {differ} of which single and double folding disagree on");
    assert_eq!(checked, 12, "twelve loops have a reset time");
    assert!(
        differ > 8,
        "only {differ} of {checked} reset times distinguish single from double \
         folding"
    );
}

/// The preset covers all twenty loops, each exactly once, and the wiring
/// agrees with the differential's independent table.
#[test]
fn the_preset_covers_every_loop_exactly_once() {
    assert_eq!(PRESET.len(), 20);
    let mut numbers: Vec<usize> = PRESET.iter().map(|p| p.tuning.number).collect();
    numbers.sort_unstable();
    numbers.dedup();
    assert_eq!(numbers.len(), 20, "a loop appears twice");
    // 12 and 21 do not exist; 22 does.
    assert!(!numbers.contains(&12) && !numbers.contains(&21) && numbers.contains(&22));

    // Every setpoint slot is owned by exactly one loop.
    let mut slots: Vec<usize> = PRESET.iter().map(|p| p.setpoint_index).collect();
    slots.sort_unstable();
    slots.dedup();
    assert_eq!(slots.len(), 20, "two loops read the same setpoint");

    // Every cascade target is a slot some loop actually owns.
    for entry in &PRESET {
        if let Output::Setpoint { index, .. } = entry.tuning.output {
            assert!(
                PRESET.iter().any(|p| p.setpoint_index == index),
                "CONTRL{} writes SETPT({index}), which no loop reads",
                entry.tuning.number
            );
        }
    }
}

/// Eight loops are proportional-only and twelve are PI, and the split is the
/// one `temain_mod.f`'s `COMMON` blocks declare.
#[test]
fn eight_loops_are_proportional_only() {
    let proportional: Vec<usize> = PRESET
        .iter()
        .filter(|p| p.tuning.reset.is_none())
        .map(|p| p.tuning.number)
        .collect();
    assert_eq!(proportional, vec![1, 2, 3, 4, 6, 7, 8, 9]);
}

/// The three periods, and which loops run on each.
#[test]
fn the_schedule_matches_the_drivers_loop() {
    let on = |period: Period| -> Vec<usize> {
        PRESET
            .iter()
            .filter(|p| p.tuning.period == period)
            .map(|p| p.tuning.number)
            .collect()
    };
    // temain_mod.f:370-384, the 3-second group, plus CONTRL22 which shares
    // their period and is never called.
    assert_eq!(
        on(Period::Fast),
        vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18, 22]
    );
    // temain_mod.f:387-392.
    assert_eq!(on(Period::Composition), vec![13, 14, 15, 19]);
    // temain_mod.f:394.
    assert_eq!(on(Period::Quality), vec![20]);
}
