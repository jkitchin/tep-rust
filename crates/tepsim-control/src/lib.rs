//! Regulatory control for the Tennessee Eastman Process.
//!
//! Ported from `reference/fortran/temain_mod.f`. The original spells this out
//! as twenty near-identical copy-pasted subroutines, each with its own
//! `COMMON` block; here it is one velocity-form implementation and a table.
//!
//! # The velocity form, and why it matters
//!
//! Each loop computes an *increment* to its output rather than the output
//! itself:
//!
//! \\[ \\Delta u = K \\left[(e - e_{-1}) + \\frac{e \\, \\Delta t \\, P}{\\tau_I}\\right] \\]
//!
//! and then `u += Δu`. Written positionally instead, a controller with a
//! steady non-zero error would ramp its valve without limit; in velocity form
//! a *proportional* loop with a settled error stops moving entirely, which is
//! measurable and is asserted.
//!
//! # Two shapes, told apart by the type
//!
//! `temain_mod.f` has two: proportional-only, which is the bracket's first
//! term alone, and PI, which adds the second. Which one a loop uses is visible
//! in the original only in whether its `COMMON` block declares a `TAUIn`.
//! Eight loops are proportional-only and twelve are PI.
//!
//! [`Tuning`] carries `reset` as an `Option`, so the shape is a property of
//! the data rather than a flag to get wrong.
//!
//! # Three things the source settles that a summary does not
//!
//! **The error is normalised by a per-loop span** before the gain is applied:
//! `(SETPT - XMEAS) * 100 / span`. So every gain is in percent-of-span, and
//! the spans differ by four orders of magnitude across the twenty loops (1.017
//! for the A feed, 8354 for the E feed).
//!
//! **`Δt · P` is the loop's own sample time**, not the plant's step. `P` is 3,
//! 360 or 900 depending on which schedule the loop runs on, and B-0035
//! confirmed it matches the scheduler loop for loop. Writing 3 everywhere
//! would scale the slow loops' integral action by 120 and 300.
//!
//! **One loop has no normalisation at all.** `CONTRL22` computes
//! `SETPT(12) - XMEAS(13)` raw. It is also the one loop the driver never
//! calls, and the two facts are probably the same fact: it looks like an
//! earlier controller left behind. [`Tuning::span`] is `None` for it.

#![forbid(unsafe_code)]
#![no_std]

use tepsim_core::constants::single;

/// How often a loop runs, as a multiple of the plant's one-second step.
///
/// `temain_mod.f:369-394` schedules on `MOD(I, period)`, and the same number
/// appears again inside each loop as the `Δt · P` of its integral term. That
/// the two agree is checked rather than assumed; see [`Tuning::period`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[repr(usize)]
pub enum Period {
    /// Every three seconds: the flow, level and pressure loops.
    Fast = 3,
    /// Every 360 seconds: the composition loops.
    Composition = 360,
    /// Every 900 seconds: the product quality loop.
    Quality = 900,
}

impl Period {
    /// The multiplier the integral term uses, as a `f64`.
    ///
    /// Written as the Fortran writes it: `3.`, `360.` and `900.` are all
    /// single-precision literals, and all three are exactly representable, so
    /// `single` changes nothing numerically here. It is applied for the reason
    /// [`tepsim_core::constants`] gives: uniformity, so a reader can check the
    /// line without arithmetic.
    #[must_use]
    pub const fn seconds(self) -> f64 {
        match self {
            Period::Fast => single(3.),
            Period::Composition => single(360.),
            Period::Quality => single(900.),
        }
    }

    /// How many plant steps between runs.
    #[must_use]
    pub const fn steps(self) -> usize {
        self as usize
    }
}

/// What a loop writes when it acts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Output {
    /// A manipulated variable, `XMV(n)`, one-based. Eleven loops do this.
    Valve(usize),
    /// Another loop's setpoint, `SETPT(n)`, one-based, rescaled by that loop's
    /// own span. Eight loops do this, and they are the cascade.
    Setpoint {
        /// Which setpoint, one-based.
        index: usize,
        /// The inner loop's span, which converts percent back to engineering
        /// units: `SETPT(n) += DXMV * span / 100`.
        span: f64,
    },
}

/// One loop's constants: everything that does not change while it runs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tuning {
    /// The `CONTRLn` number, for cross-referencing the source.
    pub number: usize,
    /// Which measurement it reads, `XMEAS(n)`, one-based.
    pub measurement: usize,
    /// What it writes.
    pub output: Output,
    /// Controller gain, in output units per percent of span.
    pub gain: f64,
    /// Reset time in hours, or `None` for a proportional-only loop.
    pub reset: Option<f64>,
    /// The span the error is normalised by, or `None` for the one loop that
    /// does not normalise. See the module documentation.
    pub span: Option<f64>,
    /// How often it runs.
    pub period: Period,
}

/// A loop's error history: the only thing it carries between calls.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Loop {
    /// `ERROLDn`: the error at the previous call.
    pub previous_error: f64,
}

impl Loop {
    /// Compute this loop's increment and advance its history.
    ///
    /// `dt` is the plant's step in hours, `DELTAT` in the original.
    ///
    /// Returns the increment; the caller applies it, because where it goes
    /// depends on [`Tuning::output`] and only the caller has the valves and
    /// setpoints.
    // @port temain_mod.f:499-511
    #[must_use]
    pub fn increment(&mut self, tuning: &Tuning, setpoint: f64, measurement: f64, dt: f64) -> f64 {
        // temain_mod.f:501. Normalised to percent of span, except for the one
        // loop with no span; see the module documentation.
        let error = match tuning.span {
            Some(span) => (setpoint - measurement) * single(100.) / span,
            None => setpoint - measurement,
        };

        // temain_mod.f:507 for the proportional form, and the PI form adds the
        // integral term. Written as two expressions rather than one with a
        // zero term, because a zero `TAUI` is a division by zero and not an
        // absent integral: the proportional loops have no `TAUI` at all.
        let increment = match tuning.reset {
            Some(reset) => {
                tuning.gain
                    * ((error - self.previous_error) + error * dt * tuning.period.seconds() / reset)
            }
            None => tuning.gain * (error - self.previous_error),
        };

        // temain_mod.f:511
        self.previous_error = error;
        increment
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn proportional() -> Tuning {
        // CONTRL1, the D feed. temain_mod.f:477-514.
        Tuning {
            number: 1,
            measurement: 2,
            output: Output::Valve(1),
            gain: 1.0,
            reset: None,
            span: Some(5811.),
            period: Period::Fast,
        }
    }

    fn integral() -> Tuning {
        // CONTRL10, the reactor coolant. temain_mod.f:883-921.
        Tuning {
            number: 10,
            measurement: 21,
            output: Output::Valve(10),
            gain: -0.156 * 10.,
            reset: Some(1452. / 3600.),
            span: Some(150.),
            period: Period::Fast,
        }
    }

    /// A proportional loop with a settled error stops moving. That is what the
    /// velocity form means, and a positional implementation would ramp.
    #[test]
    fn a_proportional_loop_stops_moving_once_the_error_settles() {
        let tuning = proportional();
        let mut state = Loop::default();
        let first = state.increment(&tuning, 4000.0, 3664.0, 1.0 / 3600.0);
        let second = state.increment(&tuning, 4000.0, 3664.0, 1.0 / 3600.0);
        let third = state.increment(&tuning, 4000.0, 3664.0, 1.0 / 3600.0);
        assert!(first != 0.0, "the first call should move the valve");
        assert_eq!(second.to_bits(), 0.0_f64.to_bits());
        assert_eq!(third.to_bits(), 0.0_f64.to_bits());
    }

    /// A PI loop keeps moving on a settled error. That is the integral action,
    /// and it is the whole difference between the two shapes.
    #[test]
    fn a_pi_loop_keeps_moving_on_a_settled_error() {
        let tuning = integral();
        let mut state = Loop::default();
        let _ = state.increment(&tuning, 100.0, 94.6, 1.0 / 3600.0);
        let second = state.increment(&tuning, 100.0, 94.6, 1.0 / 3600.0);
        let third = state.increment(&tuning, 100.0, 94.6, 1.0 / 3600.0);
        assert!(second != 0.0, "the integral term should keep it moving");
        assert_eq!(
            second.to_bits(),
            third.to_bits(),
            "with the error constant the integral contribution is constant too"
        );
    }

    /// The period enters the integral term, so the same loop on a slower
    /// schedule integrates proportionally faster per call.
    #[test]
    fn the_period_scales_the_integral_action() {
        let mut fast = integral();
        fast.period = Period::Fast;
        let mut slow = integral();
        slow.period = Period::Composition;

        let mut a = Loop::default();
        let mut b = Loop::default();
        let _ = a.increment(&fast, 100.0, 94.6, 1.0 / 3600.0);
        let _ = b.increment(&slow, 100.0, 94.6, 1.0 / 3600.0);
        let after_fast = a.increment(&fast, 100.0, 94.6, 1.0 / 3600.0);
        let after_slow = b.increment(&slow, 100.0, 94.6, 1.0 / 3600.0);

        let ratio = after_slow / after_fast;
        assert!(
            (ratio - 120.0).abs() < 1e-9,
            "a 360-second loop should integrate 120 times as much per call as \
             a 3-second one, got {ratio}"
        );
    }

    /// The error is normalised by the span, so the same absolute offset on two
    /// loops with different spans gives different increments.
    #[test]
    fn the_span_normalises_the_error() {
        let mut narrow = proportional();
        narrow.span = Some(1.017);
        let mut wide = proportional();
        wide.span = Some(8354.);

        let mut a = Loop::default();
        let mut b = Loop::default();
        let narrow_move = a.increment(&narrow, 1.0, 0.0, 1.0 / 3600.0);
        let wide_move = b.increment(&wide, 1.0, 0.0, 1.0 / 3600.0);
        assert!(
            narrow_move.abs() > wide_move.abs() * 1000.0,
            "the narrow-span loop should react far more strongly to the same \
             absolute error: {narrow_move} against {wide_move}"
        );
    }

    /// The unnormalised loop takes the raw difference.
    #[test]
    fn a_loop_with_no_span_uses_the_raw_error() {
        let tuning = Tuning {
            number: 22,
            measurement: 13,
            output: Output::Valve(6),
            // `GAIN22 = -1.0 * 5.` at temain_mod.f:315, transcribed as the
            // product it is written as rather than as -5. Clippy would rather
            // it were simplified; the whole file's discipline is the other
            // way, and B-0019 measured what happens when a written expression
            // is folded by hand instead of by the compiler.
            #[allow(clippy::neg_multiply, reason = "transcribed from temain_mod.f:315")]
            gain: -1.0 * 5.,
            reset: Some(1000. / 3600.),
            span: None,
            period: Period::Fast,
        };
        let mut state = Loop::default();
        let _ = state.increment(&tuning, 2633.7, 2600.0, 1.0 / 3600.0);
        assert_eq!(
            state.previous_error.to_bits(),
            (2633.7_f64 - 2600.0).to_bits(),
            "CONTRL22 normalises by nothing"
        );
    }

    /// The three periods are the three the scheduler uses.
    #[test]
    fn the_periods_match_the_scheduler() {
        assert_eq!(Period::Fast.steps(), 3);
        assert_eq!(Period::Composition.steps(), 360);
        assert_eq!(Period::Quality.steps(), 900);
        // And the integral multiplier is the same number.
        for period in [Period::Fast, Period::Composition, Period::Quality] {
            assert_eq!(
                period.seconds().to_bits(),
                (period.steps() as f64).to_bits(),
                "the integral multiplier and the schedule have parted company"
            );
        }
    }
}

/// The purge valve's pressure override, `CONTRL6`.
///
/// Ported from `temain_mod.f:710-753`. Alone among the twenty, `CONTRL6` is
/// not a PI loop with a wrapper: it is a latching state machine that *replaces*
/// the loop while it is engaged.
///
/// # What it does
///
/// | Separator pressure | Action |
/// |---|---|
/// | above 2950 | purge wide open, latch [`Override::Open`] |
/// | latched open, still above 2633.7 | hold open |
/// | latched open, back below 2633.7 | reset and release |
/// | below 2300 | purge shut, latch [`Override::Shut`] |
/// | latched shut, still below 2633.7 | hold shut |
/// | latched shut, back above 2633.7 | reset and release |
/// | otherwise | run the PI loop |
///
/// The PI is inside the final `ELSE`, so **while an override is latched the
/// controller does not run and its error history does not advance**. That is
/// checked against the Fortran in `tests/driver_binding.rs`, not inferred from
/// the indentation.
///
/// # Why this exists, and what it replaced
///
/// `CONTRL22` is a separator-pressure controller writing the same valve, with
/// its setpoint initialised to 2633.7 -- exactly this override's release
/// threshold. It is defined, initialised, and never called. The override
/// appears to be what replaced it; see delta D-008.
///
/// # Precision
///
/// `40.060` and `0.33712` at `temain_mod.f:716-717` are **single precision**,
/// and the thresholds are mixed: 2950 and 2300 are exactly representable and
/// 2633.7 is not. Read off each line, not inferred from its neighbour.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Override {
    /// The PI loop is running. `FLAG = 0`.
    #[default]
    Released,
    /// Latched wide open on high pressure. `FLAG = 1`.
    Open,
    /// Latched shut on low pressure. `FLAG = 2`.
    Shut,
}

/// Pressure above which the purge is thrown open (`temain_mod.f:710`).
const OVERRIDE_HIGH: f64 = single(2950.0);
/// Pressure below which it is shut (`temain_mod.f:720`).
const OVERRIDE_LOW: f64 = single(2300.);
/// The pressure both latches release through (`temain_mod.f:713`).
///
/// Not exactly representable, unlike the two limits above.
const OVERRIDE_RELEASE: f64 = single(2633.7);
/// The valve position the loop is reset to (`temain_mod.f:716`).
const OVERRIDE_RESET_VALVE: f64 = single(40.060);
/// The setpoint it is reset to (`temain_mod.f:717`).
const OVERRIDE_RESET_SETPOINT: f64 = single(0.33712);

/// What the override decided this call.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Purge {
    /// Drive the valve to this position, and do not run the loop.
    Hold(f64),
    /// Reset the valve, the setpoint and the error history, and do not run
    /// the loop.
    Release {
        /// Where the valve goes.
        valve: f64,
        /// Where `SETPT(6)` goes.
        setpoint: f64,
    },
    /// Run the PI loop as usual.
    Run,
}

impl Override {
    /// Step the state machine on the current separator pressure.
    ///
    /// `pressure` is `XMEAS(13)`, in kPa gauge.
    // @port temain_mod.f:710-731
    #[must_use]
    pub fn step(&mut self, pressure: f64) -> Purge {
        // The branch order is the original's, and it matters: the `>= 2950`
        // test comes first, so a pressure above 2950 re-latches open even if
        // the machine was latched shut.
        if pressure >= OVERRIDE_HIGH {
            *self = Override::Open;
            Purge::Hold(100.0)
        } else if *self == Override::Open && pressure >= OVERRIDE_RELEASE {
            Purge::Hold(100.0)
        } else if *self == Override::Open {
            // `<= 2633.7`, which the two branches above have already narrowed
            // this to.
            *self = Override::Released;
            Purge::Release {
                valve: OVERRIDE_RESET_VALVE,
                setpoint: OVERRIDE_RESET_SETPOINT,
            }
        } else if pressure <= OVERRIDE_LOW {
            *self = Override::Shut;
            Purge::Hold(0.0)
        } else if *self == Override::Shut && pressure <= OVERRIDE_RELEASE {
            Purge::Hold(0.0)
        } else if *self == Override::Shut {
            *self = Override::Released;
            Purge::Release {
                valve: OVERRIDE_RESET_VALVE,
                setpoint: OVERRIDE_RESET_SETPOINT,
            }
        } else {
            *self = Override::Released;
            Purge::Run
        }
    }
}

#[cfg(test)]
mod override_tests {
    use super::*;

    /// The high branch latches, holds, and releases through 2633.7.
    #[test]
    fn the_high_override_latches_holds_and_releases() {
        let mut state = Override::default();
        assert_eq!(state.step(3000.0), Purge::Hold(100.0));
        assert_eq!(state, Override::Open);
        // Still above the release threshold: hold.
        assert_eq!(state.step(2800.0), Purge::Hold(100.0));
        assert_eq!(state, Override::Open);
        // Through it: reset and release.
        assert_eq!(
            state.step(2600.0),
            Purge::Release {
                valve: OVERRIDE_RESET_VALVE,
                setpoint: OVERRIDE_RESET_SETPOINT
            }
        );
        assert_eq!(state, Override::Released);
        // And now the loop runs again.
        assert_eq!(state.step(2600.0), Purge::Run);
    }

    /// The low branch does the same in the other direction.
    #[test]
    fn the_low_override_latches_holds_and_releases() {
        let mut state = Override::default();
        assert_eq!(state.step(2000.0), Purge::Hold(0.0));
        assert_eq!(state, Override::Shut);
        assert_eq!(state.step(2500.0), Purge::Hold(0.0));
        assert_eq!(
            state.step(2700.0),
            Purge::Release {
                valve: OVERRIDE_RESET_VALVE,
                setpoint: OVERRIDE_RESET_SETPOINT
            }
        );
        assert_eq!(state, Override::Released);
    }

    /// The hysteresis is real: between 2300 and 2950 with nothing latched, the
    /// loop simply runs.
    #[test]
    fn the_band_between_the_limits_runs_the_loop() {
        let mut state = Override::default();
        for pressure in [2301.0, 2500.0, 2633.7, 2800.0, 2949.0] {
            assert_eq!(state.step(pressure), Purge::Run, "at {pressure}");
            assert_eq!(state, Override::Released);
        }
    }

    /// The branch order matters: the high test comes first, so a pressure
    /// above 2950 re-latches open even from the shut state.
    #[test]
    fn a_high_pressure_relatches_open_from_the_shut_state() {
        let mut state = Override::Shut;
        assert_eq!(state.step(3000.0), Purge::Hold(100.0));
        assert_eq!(
            state,
            Override::Open,
            "the high branch must be tested before the latch branches"
        );
    }

    /// The two limits are exactly representable and the release threshold is
    /// not, which is the kind of inconsistency that has to be read off each
    /// line.
    #[test]
    fn the_thresholds_have_mixed_precision() {
        assert_eq!(OVERRIDE_HIGH.to_bits(), 2950.0_f64.to_bits());
        assert_eq!(OVERRIDE_LOW.to_bits(), 2300.0_f64.to_bits());
        assert_ne!(
            OVERRIDE_RELEASE.to_bits(),
            2633.7_f64.to_bits(),
            "2633.7 is not exactly representable in binary32, so the single \
             and double forms must differ"
        );
        assert_ne!(OVERRIDE_RESET_VALVE.to_bits(), 40.060_f64.to_bits());
        assert_ne!(OVERRIDE_RESET_SETPOINT.to_bits(), 0.33712_f64.to_bits());
    }

    /// Exactly on a threshold, the comparisons are `>=` and `<=`, so the
    /// override engages rather than the loop running.
    #[test]
    fn the_limits_are_inclusive() {
        let mut state = Override::default();
        assert_eq!(state.step(OVERRIDE_HIGH), Purge::Hold(100.0));
        let mut state = Override::default();
        assert_eq!(state.step(OVERRIDE_LOW), Purge::Hold(0.0));
    }
}
