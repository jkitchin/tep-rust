//! Every `CONTRLn` against `tepsim_control::Loop`, over a swept range of
//! setpoints, measurements and tunings.
//!
//! B-0037. Twenty subroutines, one implementation. The structural facts about
//! each loop -- which measurement, which output, which span, which period --
//! are written out here from `temain_mod.f` and are what the comparison
//! *tests*; the tuning constants are read out of the oracle's `COMMON` rather
//! than transcribed, so this file does not depend on B-0038 having got them
//! right.
//!
//! # What each case does
//!
//! Set the loop's gain, reset time and error history in `COMMON`. Set the
//! setpoint it reads and the measurement it reads. Snapshot its output. Call
//! the subroutine. Then run the port from the same numbers and compare both
//! the increment that was applied and the error the loop stored.
//!
//! Comparing the stored error as well as the output matters: the two shapes
//! differ only in the increment, so a port that computed the error correctly
//! and the increment wrongly, or the reverse, would be caught by one and not
//! the other.
//!
//! # The spans are single-precision literals
//!
//! Every one is written without a `D` suffix. Sixteen of the eighteen are
//! exactly representable so it makes no difference, and *two are not*:
//! `1.017` is stored as 1.0169999599456787 and `1.6` as 1.600000023841858.
//! Writing them as `f64` gives an error of 5.4e-8 on `CONTRL3`, which is five
//! orders past the Tier 1 gate. Same trap as everywhere else in this project;
//! the first version of this file had it.

// Differential tests: exact comparisons are the property under test, and
// arithmetic is transcribed from `teprob.f` so that a reader can check it
// against the listing line by line. Rearranging either would defeat the point.
#![allow(
    clippy::float_cmp,
    reason = "bit equality against the Fortran is the property under test"
)]
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions are transcribed to be checkable against teprob.f"
)]
#![cfg(feature = "oracle")]

use tepsim_control::{Loop, Output, Period, Tuning};
use tepsim_core::constants::single;
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::{Comparison, Sampler};

/// `PLAN.org`, "Tier 1".
const TOLERANCE: f64 = 1e-13;

/// One loop's structure, read from `temain_mod.f`.
struct Wiring {
    number: usize,
    /// `SETPT(n)`, one-based.
    setpoint: usize,
    /// `XMEAS(n)`, one-based.
    measurement: usize,
    output: Output,
    span: Option<f64>,
    period: Period,
    /// Whether the loop has a `TAUI`.
    integral: bool,
    /// Whether the sweep can drive it. False for `CONTRL6`, which is an
    /// override state machine rather than a PI loop; it is in the table only
    /// so the cascade check can find its span, since `CONTRL19` writes its
    /// setpoint.
    sweepable: bool,
}

/// The twenty, from `temain_mod.f`. Extracted mechanically; see B-0035.
fn wiring() -> Vec<Wiring> {
    use Output::{Setpoint, Valve};
    use Period::{Composition, Fast, Quality};
    let w = |number, setpoint, measurement, output, span, period, integral| Wiring {
        number,
        setpoint,
        measurement,
        output,
        span,
        period,
        integral,
        sweepable: number != 6,
    };
    vec![
        w(1, 1, 2, Valve(1), Some(single(5811.)), Fast, false),
        w(2, 2, 3, Valve(2), Some(single(8354.)), Fast, false),
        w(3, 3, 1, Valve(3), Some(single(1.017)), Fast, false),
        w(4, 4, 4, Valve(4), Some(single(15.25)), Fast, false),
        w(5, 5, 5, Valve(5), Some(single(53.)), Fast, true),
        // CONTRL6 is in the table but not swept: it is an override state
        // machine rather than a PI loop, and B-0037b covers it. It is here
        // because CONTRL19 cascades onto its setpoint, so the cascade check
        // needs its span.
        w(6, 6, 10, Valve(6), Some(single(1.)), Fast, false),
        w(7, 7, 12, Valve(7), Some(single(70.)), Fast, false),
        w(8, 8, 15, Valve(8), Some(single(70.)), Fast, false),
        w(9, 9, 19, Valve(9), Some(single(460.)), Fast, false),
        w(10, 10, 21, Valve(10), Some(single(150.)), Fast, true),
        w(11, 11, 17, Valve(11), Some(single(46.)), Fast, true),
        // The cascade: each writes another loop's setpoint, rescaled by *that*
        // loop's span. Compare against the table above and the pairs line up.
        w(
            13,
            13,
            23,
            Setpoint {
                index: 3,
                span: single(1.017),
            },
            Some(single(100.)),
            Composition,
            true,
        ),
        w(
            14,
            14,
            26,
            Setpoint {
                index: 1,
                span: single(5811.),
            },
            Some(single(100.)),
            Composition,
            true,
        ),
        w(
            15,
            15,
            27,
            Setpoint {
                index: 2,
                span: single(8354.),
            },
            Some(single(100.)),
            Composition,
            true,
        ),
        w(
            16,
            16,
            18,
            Setpoint {
                index: 9,
                span: single(460.),
            },
            Some(single(130.)),
            Fast,
            true,
        ),
        w(
            17,
            17,
            8,
            Setpoint {
                index: 4,
                span: single(15.25),
            },
            Some(single(50.)),
            Fast,
            true,
        ),
        w(
            18,
            18,
            9,
            Setpoint {
                index: 10,
                span: single(150.),
            },
            Some(single(150.)),
            Fast,
            true,
        ),
        w(
            19,
            19,
            30,
            Setpoint {
                index: 6,
                span: single(1.),
            },
            Some(single(26.)),
            Composition,
            true,
        ),
        w(
            20,
            20,
            38,
            Setpoint {
                index: 16,
                span: single(130.),
            },
            Some(single(1.6)),
            Quality,
            true,
        ),
        // The dead one. No span at all, which is consistent with it predating
        // the others; see B-0035.
        w(22, 12, 13, Valve(6), None, Fast, true),
    ]
}

/// Set one loop's tuning in `COMMON`, and read back what it stored.
///
/// A macro rather than a trait: the twenty `COMMON` blocks are twenty
/// different types by construction, since eight of them have no `taui` field,
/// and that difference is the thing worth preserving.
macro_rules! dispatch {
    ($oracle:expr, $n:expr, $gain:expr, $taui:expr, $errold:expr, $call:tt) => {{
        match $n {
            1 => {
                $oracle.set_ctrl1(&tepsim_oracle::Ctrl1 {
                    gain: $gain,
                    errold: $errold,
                });
                $oracle.contrl1();
                $oracle.ctrl1().errold
            }
            2 => {
                $oracle.set_ctrl2(&tepsim_oracle::Ctrl2 {
                    gain: $gain,
                    errold: $errold,
                });
                $oracle.contrl2();
                $oracle.ctrl2().errold
            }
            3 => {
                $oracle.set_ctrl3(&tepsim_oracle::Ctrl3 {
                    gain: $gain,
                    errold: $errold,
                });
                $oracle.contrl3();
                $oracle.ctrl3().errold
            }
            4 => {
                $oracle.set_ctrl4(&tepsim_oracle::Ctrl4 {
                    gain: $gain,
                    errold: $errold,
                });
                $oracle.contrl4();
                $oracle.ctrl4().errold
            }
            5 => {
                $oracle.set_ctrl5(&tepsim_oracle::Ctrl5 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl5();
                $oracle.ctrl5().errold
            }
            7 => {
                $oracle.set_ctrl7(&tepsim_oracle::Ctrl7 {
                    gain: $gain,
                    errold: $errold,
                });
                $oracle.contrl7();
                $oracle.ctrl7().errold
            }
            8 => {
                $oracle.set_ctrl8(&tepsim_oracle::Ctrl8 {
                    gain: $gain,
                    errold: $errold,
                });
                $oracle.contrl8();
                $oracle.ctrl8().errold
            }
            9 => {
                $oracle.set_ctrl9(&tepsim_oracle::Ctrl9 {
                    gain: $gain,
                    errold: $errold,
                });
                $oracle.contrl9();
                $oracle.ctrl9().errold
            }
            10 => {
                $oracle.set_ctrl10(&tepsim_oracle::Ctrl10 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl10();
                $oracle.ctrl10().errold
            }
            11 => {
                $oracle.set_ctrl11(&tepsim_oracle::Ctrl11 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl11();
                $oracle.ctrl11().errold
            }
            13 => {
                $oracle.set_ctrl13(&tepsim_oracle::Ctrl13 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl13();
                $oracle.ctrl13().errold
            }
            14 => {
                $oracle.set_ctrl14(&tepsim_oracle::Ctrl14 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl14();
                $oracle.ctrl14().errold
            }
            15 => {
                $oracle.set_ctrl15(&tepsim_oracle::Ctrl15 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl15();
                $oracle.ctrl15().errold
            }
            16 => {
                $oracle.set_ctrl16(&tepsim_oracle::Ctrl16 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl16();
                $oracle.ctrl16().errold
            }
            17 => {
                $oracle.set_ctrl17(&tepsim_oracle::Ctrl17 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl17();
                $oracle.ctrl17().errold
            }
            18 => {
                $oracle.set_ctrl18(&tepsim_oracle::Ctrl18 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl18();
                $oracle.ctrl18().errold
            }
            19 => {
                $oracle.set_ctrl19(&tepsim_oracle::Ctrl19 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl19();
                $oracle.ctrl19().errold
            }
            20 => {
                $oracle.set_ctrl20(&tepsim_oracle::Ctrl20 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl20();
                $oracle.ctrl20().errold
            }
            22 => {
                $oracle.set_ctrl22(&tepsim_oracle::Ctrl22 {
                    gain: $gain,
                    taui: $taui,
                    errold: $errold,
                });
                $oracle.contrl22();
                $oracle.ctrl22().errold
            }
            other => panic!("no controller {other}"),
        }
    }};
}

/// Identifies one compared number.
#[derive(Clone, Copy, Debug)]
struct Case {
    loop_number: usize,
    trial: usize,
    what: &'static str,
}

impl core::fmt::Display for Case {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "CONTRL{}#{}[{}]",
            self.loop_number, self.trial, self.what
        )
    }
}

fn between(sampler: &mut Sampler, lo: f64, hi: f64) -> f64 {
    lo + (hi - lo) * sampler.unit()
}

#[test]
fn every_controller_matches_the_fortran_over_a_sweep() {
    let mut oracle = Oracle::lock();
    let _ = oracle.init();
    let mut increment: Comparison<Case> = Comparison::new("controller increment");
    let mut stored: Comparison<Case> = Comparison::new("stored error");
    let mut sampler = Sampler::new(0x7E2_0037);
    let dt = 1.0 / 3600.0;

    for wiring in wiring().into_iter().filter(|w| w.sweepable) {
        for trial in 0..500 {
            let gain = between(&mut sampler, -5.0, 5.0);
            let taui = between(&mut sampler, 0.1, 4.0);
            let errold = between(&mut sampler, -50.0, 50.0);
            let setpoint = between(&mut sampler, -100.0, 5000.0);
            let measurement = between(&mut sampler, -100.0, 5000.0);

            // Put both sides on the same footing.
            let mut all = oracle.ctrlall();
            all.deltat = dt;
            all.setpt[wiring.setpoint - 1] = setpoint;
            // Zero whatever this loop writes, so that after the call the
            // output *is* the increment. Starting from a realistic value like
            // 50 and taking a difference loses most of the precision when the
            // increment is small, which showed up as 1.2e-11 of measurement
            // error and nothing to do with the port.
            match wiring.output {
                Output::Valve(v) => {
                    let mut xmv = oracle.manipulated();
                    xmv[v - 1] = 0.0;
                    oracle.set_manipulated(&xmv);
                }
                Output::Setpoint { index, .. } => {
                    all.setpt[index - 1] = 0.0;
                }
            }
            oracle.set_ctrlall(&all);
            let mut xmeas = oracle.measurements();
            xmeas[wiring.measurement - 1] = measurement;
            oracle.set_measurements(&xmeas);

            let their_error = dispatch!(oracle, wiring.number, gain, taui, errold, ());
            let after = match wiring.output {
                Output::Valve(v) => oracle.manipulated()[v - 1],
                Output::Setpoint { index, .. } => oracle.ctrlall().setpt[index - 1],
            };
            let their_move = after;

            // The port, from the same numbers.
            let tuning = Tuning {
                number: wiring.number,
                measurement: wiring.measurement,
                output: wiring.output,
                gain,
                reset: wiring.integral.then_some(taui),
                span: wiring.span,
                period: wiring.period,
            };
            let mut state = Loop {
                previous_error: errold,
            };
            let raw = state.increment(&tuning, setpoint, measurement, dt);
            // A cascade output is rescaled into the inner loop's units.
            let our_move = match wiring.output {
                Output::Valve(_) => raw,
                Output::Setpoint { span, .. } => raw * span / 100.,
            };

            let case = |what| Case {
                loop_number: wiring.number,
                trial,
                what,
            };
            increment.observe(case("move"), our_move, their_move);
            stored.observe(case("errold"), state.previous_error, their_error);
        }
    }

    println!("{increment}");
    println!("{stored}");
    increment.assert_within(TOLERANCE);
    stored.assert_within(TOLERANCE);
    assert_eq!(
        stored.max_ulp(),
        0,
        "the stored error is one subtraction, one multiply and one divide, \
         with no transcendental anywhere, so anything but bit equality is a \
         porting error"
    );
}

/// The cascade rescaling is the *inner* loop's span, and getting it from the
/// outer loop would be a plausible mistake.
///
/// Read off the source: `CONTRL13` writes `SETPT(3)` scaled by 1.017, and
/// 1.017 is `CONTRL3`'s span, not `CONTRL13`'s (which is 100).
#[test]
fn a_cascade_rescales_by_the_inner_loops_span() {
    let table = wiring();
    let span_of = |n: usize| {
        table
            .iter()
            .find(|w| w.number == n)
            .and_then(|w| w.span)
            .expect("a loop with a span")
    };
    for entry in &table {
        if let Output::Setpoint { index, span } = entry.output {
            // The inner loop is the one whose *setpoint index* matches.
            let inner = table
                .iter()
                .find(|w| w.setpoint == index)
                .unwrap_or_else(|| panic!("nothing owns SETPT({index})"));
            assert_eq!(
                span.to_bits(),
                span_of(inner.number).to_bits(),
                "CONTRL{} rescales by {span}, but SETPT({index}) belongs to \
                 CONTRL{}, whose span is {}",
                entry.number,
                inner.number,
                span_of(inner.number)
            );
            assert_ne!(
                entry.number, inner.number,
                "a loop cannot cascade onto itself"
            );
        }
    }
}

/// The purge override against `CONTRL6`, over a pressure walk that enters and
/// leaves both latches many times.
///
/// A single-step comparison would miss the whole point: the machine's answer
/// depends on where it has been, so the sequence is what has to agree. This
/// walks a pressure trajectory that crosses every threshold in both directions
/// and compares the valve, the setpoint, the latch and the error history at
/// every step.
#[test]
fn the_purge_override_matches_the_fortran_over_a_pressure_walk() {
    use tepsim_control::{Override, Purge};

    let mut oracle = Oracle::lock();
    let _ = oracle.init();

    let gain = 1.22;
    let mut all = oracle.ctrlall();
    all.deltat = 1.0 / 3600.0;
    all.setpt[5] = 0.33712;
    oracle.set_ctrlall(&all);
    oracle.set_ctrl6(&tepsim_oracle::Ctrl6 { gain, errold: 0.0 });
    oracle.set_flag6(0);
    let mut xmv = oracle.manipulated();
    xmv[5] = 40.0;
    oracle.set_manipulated(&xmv);

    // The port's mirror of the same condition.
    let mut state = Override::Released;
    let mut valve = 40.0_f64;
    let mut setpoint = 0.33712_f64;
    let mut pi = tepsim_control::Loop::default();
    let tuning = tepsim_control::Tuning {
        number: 6,
        measurement: 10,
        output: tepsim_control::Output::Valve(6),
        gain,
        reset: None,
        span: Some(single(1.)),
        period: tepsim_control::Period::Fast,
    };

    // A walk that crosses 2950, 2633.7 and 2300 repeatedly, in both
    // directions, including landing exactly on each.
    let mut pressures = Vec::new();
    for cycle in 0..6 {
        let base = 2600.0 + f64::from(cycle) * 7.0;
        pressures.extend([
            base, 2960.0, 2900.0, 2700.0, 2633.7, 2600.0, 2400.0, 2250.0, 2280.0, 2500.0, 2633.7,
            2700.0, 2950.0, 3100.0, 2634.0, 2300.0, 2900.0, base,
        ]);
    }

    let purge_reading = 0.25_f64;
    for (step, pressure) in pressures.iter().copied().enumerate() {
        // Both sides see the same measurements.
        let mut xmeas = oracle.measurements();
        xmeas[12] = pressure;
        xmeas[9] = purge_reading;
        oracle.set_measurements(&xmeas);
        oracle.contrl6();

        // The port.
        match state.step(pressure) {
            Purge::Hold(position) => valve = position,
            Purge::Release {
                valve: v,
                setpoint: s,
            } => {
                valve = v;
                setpoint = s;
                pi.previous_error = 0.0;
            }
            Purge::Run => {
                valve += pi.increment(&tuning, setpoint, purge_reading, 1.0 / 3600.0);
            }
        }

        let expected_flag = match state {
            Override::Released => 0,
            Override::Open => 1,
            Override::Shut => 2,
        };
        assert_eq!(
            oracle.flag6(),
            expected_flag,
            "step {step} at {pressure} kPa: the latch disagrees"
        );
        assert_eq!(
            oracle.manipulated()[5].to_bits(),
            valve.to_bits(),
            "step {step} at {pressure} kPa: XMV(6) is {}, the port has {valve}",
            oracle.manipulated()[5]
        );
        assert_eq!(
            oracle.ctrlall().setpt[5].to_bits(),
            setpoint.to_bits(),
            "step {step} at {pressure} kPa: SETPT(6) disagrees"
        );
        assert_eq!(
            oracle.ctrl6().errold.to_bits(),
            pi.previous_error.to_bits(),
            "step {step} at {pressure} kPa: ERROLD6 disagrees. While an \
             override is latched the PI does not run, so this catches a port \
             that advanced the error history anyway."
        );
    }
    println!("{} steps of pressure walk, all exact", pressures.len());
}
