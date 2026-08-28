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

/// A loop's tuning together with the setpoint the driver starts it at.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Preset {
    /// The loop's constants.
    pub tuning: Tuning,
    /// Which `SETPT(n)` it reads, one-based.
    pub setpoint_index: usize,
    /// The value the driver initialises that setpoint to.
    pub setpoint: f64,
}

/// The twenty loops as `temain_mod.f` tunes them: the Braatz closed-loop
/// preset.
///
/// # The gains are written as arithmetic, and are transcribed that way
///
/// Four of them are expressions rather than numbers:
///
/// ```fortran
///       GAIN10= -0.156     * 10.
///       GAIN19=-83.2    / 5. /3.
///       GAIN20=-16.3   / 5.
///       GAIN22=-1.0    * 5.
/// ```
///
/// These look like tuning that was adjusted by a factor and left showing its
/// working. They are transcribed as written, not folded by hand, for the
/// reason `crate::single`'s documentation gives and that B-0019 measured: a
/// quotient of single-precision literals is evaluated at *single* precision by
/// Fortran's typing rules, and `-83.2 / 5. / 3.` is exactly that shape.
///
/// Every literal here carries no `D` suffix, so every one is single precision.
///
/// # Slot 12 belongs to the loop that never runs
///
/// `SETPT(12)` is `CONTRL22`'s, and the driver initialises it to 2633.7 --
/// exactly `CONTRL6`'s override release threshold. See delta D-008.
//
// @port temain_mod.f:246-317
pub const PRESET: [Preset; 20] = [
    // CONTRL1. temain_mod.f:246-247.
    Preset {
        tuning: Tuning {
            number: 1,
            measurement: 2,
            output: Output::Valve(1),
            gain: single(1.0),
            reset: None,
            span: Some(single(5811.0)),
            period: Period::Fast,
        },
        setpoint_index: 1,
        setpoint: single(3664.0),
    },
    // CONTRL2. temain_mod.f:249-250.
    Preset {
        tuning: Tuning {
            number: 2,
            measurement: 3,
            output: Output::Valve(2),
            gain: single(1.0),
            reset: None,
            span: Some(single(8354.0)),
            period: Period::Fast,
        },
        setpoint_index: 2,
        setpoint: single(4509.3),
    },
    // CONTRL3. temain_mod.f:252-253.
    Preset {
        tuning: Tuning {
            number: 3,
            measurement: 1,
            output: Output::Valve(3),
            gain: single(1.0),
            reset: None,
            span: Some(single(1.017)),
            period: Period::Fast,
        },
        setpoint_index: 3,
        setpoint: single(0.25052),
    },
    // CONTRL4. temain_mod.f:255-256.
    Preset {
        tuning: Tuning {
            number: 4,
            measurement: 4,
            output: Output::Valve(4),
            gain: single(1.0),
            reset: None,
            span: Some(single(15.25)),
            period: Period::Fast,
        },
        setpoint_index: 4,
        setpoint: single(9.3477),
    },
    // CONTRL5. temain_mod.f:258-260.
    Preset {
        tuning: Tuning {
            number: 5,
            measurement: 5,
            output: Output::Valve(5),
            gain: single(-0.083),
            reset: Some(single(1.0 / 3600.0)),
            span: Some(single(53.0)),
            period: Period::Fast,
        },
        setpoint_index: 5,
        setpoint: single(26.902),
    },
    // CONTRL6. temain_mod.f:262-263.
    Preset {
        tuning: Tuning {
            number: 6,
            measurement: 10,
            output: Output::Valve(6),
            gain: single(1.22),
            reset: None,
            span: Some(single(1.0)),
            period: Period::Fast,
        },
        setpoint_index: 6,
        setpoint: single(0.33712),
    },
    // CONTRL7. temain_mod.f:265-266.
    Preset {
        tuning: Tuning {
            number: 7,
            measurement: 12,
            output: Output::Valve(7),
            gain: single(-2.06),
            reset: None,
            span: Some(single(70.0)),
            period: Period::Fast,
        },
        setpoint_index: 7,
        setpoint: single(50.0),
    },
    // CONTRL8. temain_mod.f:268-269.
    Preset {
        tuning: Tuning {
            number: 8,
            measurement: 15,
            output: Output::Valve(8),
            gain: single(-1.62),
            reset: None,
            span: Some(single(70.0)),
            period: Period::Fast,
        },
        setpoint_index: 8,
        setpoint: single(50.0),
    },
    // CONTRL9. temain_mod.f:271-272.
    Preset {
        tuning: Tuning {
            number: 9,
            measurement: 19,
            output: Output::Valve(9),
            gain: single(0.41),
            reset: None,
            span: Some(single(460.0)),
            period: Period::Fast,
        },
        setpoint_index: 9,
        setpoint: single(230.31),
    },
    // CONTRL10. temain_mod.f:274-276.
    Preset {
        tuning: Tuning {
            number: 10,
            measurement: 21,
            output: Output::Valve(10),
            gain: single(-0.156 * 10.0),
            reset: Some(single(1452.0 / 3600.0)),
            span: Some(single(150.0)),
            period: Period::Fast,
        },
        setpoint_index: 10,
        setpoint: single(94.599),
    },
    // CONTRL11. temain_mod.f:278-280.
    Preset {
        tuning: Tuning {
            number: 11,
            measurement: 17,
            output: Output::Valve(11),
            gain: single(1.09),
            reset: Some(single(2600.0 / 3600.0)),
            span: Some(single(46.0)),
            period: Period::Fast,
        },
        setpoint_index: 11,
        setpoint: single(22.949),
    },
    // CONTRL13. temain_mod.f:282-284. Cascades onto CONTRL3.
    Preset {
        tuning: Tuning {
            number: 13,
            measurement: 23,
            output: Output::Setpoint {
                index: 3,
                span: single(1.017),
            },
            gain: single(18.0),
            reset: Some(single(3168.0 / 3600.0)),
            span: Some(single(100.0)),
            period: Period::Composition,
        },
        setpoint_index: 13,
        setpoint: single(32.188),
    },
    // CONTRL14. temain_mod.f:286-288. Cascades onto CONTRL1.
    Preset {
        tuning: Tuning {
            number: 14,
            measurement: 26,
            output: Output::Setpoint {
                index: 1,
                span: single(5811.0),
            },
            gain: single(8.3),
            reset: Some(single(3168.0 / 3600.0)),
            span: Some(single(100.0)),
            period: Period::Composition,
        },
        setpoint_index: 14,
        setpoint: single(6.8820),
    },
    // CONTRL15. temain_mod.f:290-292. Cascades onto CONTRL2.
    Preset {
        tuning: Tuning {
            number: 15,
            measurement: 27,
            output: Output::Setpoint {
                index: 2,
                span: single(8354.0),
            },
            gain: single(2.37),
            reset: Some(single(5069.0 / 3600.0)),
            span: Some(single(100.0)),
            period: Period::Composition,
        },
        setpoint_index: 15,
        setpoint: single(18.776),
    },
    // CONTRL16. temain_mod.f:294-296. Cascades onto CONTRL9.
    Preset {
        tuning: Tuning {
            number: 16,
            measurement: 18,
            output: Output::Setpoint {
                index: 9,
                span: single(460.0),
            },
            gain: single(1.69 / 10.0),
            reset: Some(single(236.0 / 3600.0)),
            span: Some(single(130.0)),
            period: Period::Fast,
        },
        setpoint_index: 16,
        setpoint: single(65.731),
    },
    // CONTRL17. temain_mod.f:298-300. Cascades onto CONTRL4.
    Preset {
        tuning: Tuning {
            number: 17,
            measurement: 8,
            output: Output::Setpoint {
                index: 4,
                span: single(15.25),
            },
            gain: single(11.1 / 10.0),
            reset: Some(single(3168.0 / 3600.0)),
            span: Some(single(50.0)),
            period: Period::Fast,
        },
        setpoint_index: 17,
        setpoint: single(75.000),
    },
    // CONTRL18. temain_mod.f:302-304. Cascades onto CONTRL10.
    Preset {
        tuning: Tuning {
            number: 18,
            measurement: 9,
            output: Output::Setpoint {
                index: 10,
                span: single(150.0),
            },
            gain: single(2.83 * 10.0),
            reset: Some(single(982.0 / 3600.0)),
            span: Some(single(150.0)),
            period: Period::Fast,
        },
        setpoint_index: 18,
        setpoint: single(120.40),
    },
    // CONTRL19. temain_mod.f:306-308. Cascades onto CONTRL6.
    Preset {
        tuning: Tuning {
            number: 19,
            measurement: 30,
            output: Output::Setpoint {
                index: 6,
                span: single(1.0),
            },
            gain: single(-83.2 / 5.0 / 3.0),
            reset: Some(single(6336.0 / 3600.0)),
            span: Some(single(26.0)),
            period: Period::Composition,
        },
        setpoint_index: 19,
        setpoint: single(13.823),
    },
    // CONTRL20. temain_mod.f:310-312. Cascades onto CONTRL16.
    Preset {
        tuning: Tuning {
            number: 20,
            measurement: 38,
            output: Output::Setpoint {
                index: 16,
                span: single(130.0),
            },
            gain: single(-16.3 / 5.0),
            reset: Some(single(12408.0 / 3600.0)),
            span: Some(single(1.6)),
            period: Period::Quality,
        },
        setpoint_index: 20,
        setpoint: single(0.83570),
    },
    // CONTRL22. temain_mod.f:314-316.
    Preset {
        tuning: Tuning {
            number: 22,
            measurement: 13,
            output: Output::Valve(6),
            // `-1.0 * 5.` at temain_mod.f:315, transcribed as written; see the
            // note on `PRESET`.
            #[allow(clippy::neg_multiply, reason = "transcribed from temain_mod.f:315")]
            gain: single(-1.0 * 5.0),
            reset: Some(single(1000.0 / 3600.0)),
            span: None,
            period: Period::Fast,
        },
        setpoint_index: 12,
        setpoint: single(2633.7),
    },
];

/// Look up a loop's preset by its `CONTRLn` number.
#[must_use]
pub fn preset(number: usize) -> Option<&'static Preset> {
    PRESET.iter().find(|p| p.tuning.number == number)
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

/// The closed-loop scheme: nineteen loops on three schedules, plus the purge
/// override.
///
/// Ported from `temain_mod.f:365-412`, the driver's main loop.
///
/// # Order is load-bearing
///
/// The loops run in *source order* within a period and communicate through
/// shared setpoints, so a cascade's outer loop and its inner loop can both run
/// on the same tick and the outer one's effect is felt immediately. Reordering
/// them changes the plant.
///
/// # The phase of each period is off by one from the obvious
///
/// `temain_mod.f:369` is `MOD(I,3)` with `I` starting at **1**, so the fast
/// loops first fire on step 3, not step 1. Starting a Rust loop at zero and
/// testing `step % 3 == 0` fires on step 0 and shifts every controller action
/// by two seconds for the entire run. [`Scheme::step`] takes the one-based
/// step number for that reason.
///
/// # `CONSHAND` clamps eleven valves, not twelve
///
/// `temain_mod.f:1401` loops `I = 1, 11`. Valve 12 is the agitator, and it is
/// never written by a controller, so clamping it would change nothing today
/// and would be wrong.
#[derive(Clone, Debug, PartialEq)]
pub struct Scheme {
    /// The twenty loops' error histories, in [`PRESET`] order.
    loops: [Loop; 20],
    /// The twenty setpoints, one-based in the Fortran and zero-based here.
    setpoints: [f64; 20],
    /// The purge override's latch.
    purge: Override,
}

impl Default for Scheme {
    /// The Braatz preset's starting condition.
    fn default() -> Self {
        let mut setpoints = [0.0; 20];
        for entry in &PRESET {
            setpoints[entry.setpoint_index - 1] = entry.setpoint;
        }
        Self {
            loops: [Loop::default(); 20],
            setpoints,
            purge: Override::Released,
        }
    }
}

/// Which loops run on a given step. `temain_mod.f:369-394`.
///
/// The fast group in source order, then the composition group, then the
/// quality loop. `CONTRL22` is absent: it is defined and never called, which
/// is delta D-008.
const FAST: [usize; 14] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 16, 17, 18];
const COMPOSITION: [usize; 4] = [13, 14, 15, 19];
const QUALITY: [usize; 1] = [20];

impl Scheme {
    /// A fresh scheme at the preset's setpoints.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// The current setpoints, one-based indexing as `SETPT(n)`.
    #[must_use]
    pub fn setpoint(&self, index: usize) -> f64 {
        self.setpoints[index - 1]
    }

    /// Set one, for a scenario that moves a setpoint.
    pub fn set_setpoint(&mut self, index: usize, value: f64) {
        self.setpoints[index - 1] = value;
    }

    /// The purge override's latch.
    #[must_use]
    pub const fn purge_override(&self) -> Override {
        self.purge
    }

    /// Run whichever loops are due on this step.
    ///
    /// `step` is one-based, as the driver's `I` is. `dt` is the plant step in
    /// hours. `measurements` must be the ones the *previous* step produced;
    /// see [`Driver`] for why, and use it rather than calling this directly.
    ///
    /// Does **not** clamp. `CONSHAND` runs after the integration, not after
    /// the controllers, and [`Driver::step`] is what puts it there.
    // @port temain_mod.f:369-394
    pub fn step(&mut self, step: usize, measurements: &[f64], valves: &mut [f64; 12], dt: f64) {
        if step % Period::Fast.steps() == 0 {
            for number in FAST {
                self.run(number, measurements, valves, dt);
            }
        }
        if step % Period::Composition.steps() == 0 {
            for number in COMPOSITION {
                self.run(number, measurements, valves, dt);
            }
        }
        if step % Period::Quality.steps() == 0 {
            for number in QUALITY {
                self.run(number, measurements, valves, dt);
            }
        }
    }

    /// `CONSHAND`, `temain_mod.f:1401-1404`. Eleven valves, not twelve.
    ///
    /// Written as the listing's two guards rather than as [`f64::clamp`],
    /// which is *not* equivalent here. `CONSHAND` tests `.LE. 0.0`, so it
    /// replaces a negative zero with a positive one; `f64::clamp` tests
    /// `<`, so it leaves negative zero alone. `teprob.f:803-804` does the
    /// same job with strict comparisons and there `clamp` is exact, which is
    /// why [`tepsim_core::Plant`] uses it and this does not.
    // @port temain_mod.f:1391-1407
    #[allow(clippy::manual_clamp, reason = "CONSHAND normalises negative zero")]
    pub fn clamp(&self, valves: &mut [f64; 12]) {
        for valve in valves.iter_mut().take(11) {
            if *valve <= 0.0 {
                *valve = 0.0;
            }
            if *valve >= 100.0 {
                *valve = 100.0;
            }
        }
    }

    /// One loop.
    fn run(&mut self, number: usize, measurements: &[f64], valves: &mut [f64; 12], dt: f64) {
        let slot = PRESET
            .iter()
            .position(|p| p.tuning.number == number)
            .expect("a scheduled loop is in the preset");
        let entry = &PRESET[slot];
        let measurement = measurements[entry.tuning.measurement - 1];

        // CONTRL6 is the override, and the loop runs only when it releases.
        if number == 6 {
            match self.purge.step(measurements[12]) {
                Purge::Hold(position) => {
                    valves[5] = position;
                    return;
                }
                Purge::Release { valve, setpoint } => {
                    valves[5] = valve;
                    self.setpoints[5] = setpoint;
                    self.loops[slot].previous_error = 0.0;
                    return;
                }
                Purge::Run => {}
            }
        }

        let setpoint = self.setpoints[entry.setpoint_index - 1];
        let increment = self.loops[slot].increment(&entry.tuning, setpoint, measurement, dt);
        match entry.tuning.output {
            Output::Valve(v) => valves[v - 1] += increment,
            Output::Setpoint { index, span } => {
                self.setpoints[index - 1] += increment * span / single(100.);
            }
        }
    }
}

/// The valve commands the driver starts a closed-loop run from.
///
/// `temain_mod.f:322-332`. **These are not `TEINIT`'s.** The driver overwrites
/// eleven of the twelve with values rounded to five significant figures:
/// `63.053` where `TEINIT` left `63.05263039`. Valve 12 is not overwritten and
/// keeps `TEINIT`'s 50.
///
/// So a closed-loop run does not start where an open-loop one does, and the
/// difference is in the fourth decimal place of ten valves. That is delta
/// D-009.
///
/// Ten, not eleven: `XMV(5)` is written `22.210` and `TEINIT` leaves
/// `YY(43) = 22.21000000`, which is already five significant figures, so the
/// two agree exactly. That coincidence is the clearest evidence that these
/// numbers *are* `TEINIT`'s rounded rather than an independent set.
///
/// Every literal is single precision, and each is written `value + 0.` in the
/// original, which appears to be a placeholder for a perturbation.
//
// @port  temain_mod.f:322-332
// @delta D-009 class=A temain_mod.f:322-332
pub const DRIVER_INITIAL_VALVES: [f64; 12] = [
    single(63.053 + 0.),
    single(53.980 + 0.),
    single(24.644 + 0.),
    single(61.302 + 0.),
    single(22.210 + 0.),
    single(40.064 + 0.),
    single(38.100 + 0.),
    single(46.534 + 0.),
    single(47.446 + 0.),
    single(41.106 + 0.),
    single(18.114 + 0.),
    // XMV(12) is not among the eleven the driver writes; it keeps TEINIT's
    // value, which is YY(50) = 50 exactly.
    50.0,
];

/// The closed-loop driver: one simulated second of `temain_mod.f`'s main loop.
///
/// [`Scheme`] knows *which* loops fire; the driver knows *when* things happen
/// relative to the plant, and that ordering is as load-bearing as the loop
/// order is.
///
/// # The controllers read stale measurements
///
/// `temain_mod.f:369-394` calls the controllers, and only then, at line 409,
/// calls `INTGTR` (and through it `TEFUNC`, which is what writes `XMEAS`). So
/// on iteration `I` every controller sees the measurements iteration `I - 1`
/// produced. There is one plant step of dead time in every loop, built into
/// the driver rather than into any controller.
///
/// Getting this backwards is not a small error. Feeding a controller the
/// measurements of the step it is about to cause makes the loop tighter than
/// the original by one sample, and B-0039 measured the result: `XMV(7)` parts
/// from the Fortran by 1.5% of range on the very first controller fire, and
/// `XMEAS(14)` is 23% out four hours later.
///
/// # `CONSHAND` runs after the integration, not after the controllers
///
/// Line 411, after line 409. `TEFUNC` clamps its own copy of the valve
/// positions (`teprob.f:803-804`), so with no sticking fault active the
/// placement is unobservable. It stops being unobservable under `IDV(14)`,
/// `IDV(15)` or `IDV(19)`: `teprob.f:801` only moves `VCV` toward `XMV` when
/// they differ by more than the stick threshold, and an unclamped `XMV` of 105
/// crosses that threshold at a different moment than a clamped 100 does.
///
/// # Example
///
/// ```
/// # use tepsim_control::{Driver, DRIVER_INITIAL_VALVES};
/// let mut driver = Driver::new();
/// assert_eq!(driver.valves(), &DRIVER_INITIAL_VALVES);
/// // Nothing fires until step 3.
/// driver.step(&[0.0; 41], 1.0 / 3600.0);
/// assert_eq!(driver.valves(), &DRIVER_INITIAL_VALVES);
/// ```
//
// @port temain_mod.f:365-412
// @delta D-010 class=A temain_mod.f:366-411
#[derive(Clone, Debug, PartialEq)]
pub struct Driver {
    scheme: Scheme,
    valves: [f64; 12],
    /// One-based, as `I` is. Incremented by [`Driver::step`].
    step: usize,
    /// The scenario as asked for, before the driver forces `IDV(12)` on.
    requested: [f64; 20],
    /// What the plant is actually handed. See [`DriverQuirks`].
    disturbances: [f64; 20],
    /// Which Class C driver quirks are fixed rather than reproduced.
    pub quirks: DriverQuirks,
}

impl Default for Driver {
    fn default() -> Self {
        Self {
            scheme: Scheme::new(),
            valves: DRIVER_INITIAL_VALVES,
            step: 0,
            requested: [0.0; 20],
            disturbances: [0.0; 20],
            quirks: DriverQuirks::default(),
        }
    }
}

/// After how many steps the driver forces `IDV(12)` on: `SSPTS = 3600 * 8`,
/// eight simulated hours at a one-second step.
///
/// `temain_mod.f:226`. The comparison is `I .GE. SSPTS` with `I` one-based, so
/// the disturbance is live for the integration of step 28,800 itself, which
/// happens at `TIME = 28799/3600 = 7.99972 h`.
pub const STEADY_STATE_STEPS: usize = 3600 * 8;

/// Which Class C quirks of the *driver* are fixed rather than reproduced.
///
/// All off by default, so the default driver reproduces `temain_mod.f`. This
/// mirrors [`tepsim_core::QuirkFixes`], which does the same job for the plant.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct DriverQuirks {
    /// When `false` (the default), the driver switches `IDV(12)` on at step
    /// [`STEADY_STATE_STEPS`] whatever scenario was asked for, exactly as
    /// `temain_mod.f:366-368` does.
    ///
    /// When `true`, the requested scenario is the scenario. That is a genuine
    /// behaviour change and it is **blocked on sign-off**; see B-0040a and
    /// delta D-011.
    pub only_the_requested_disturbances: bool,
}

impl Driver {
    /// A driver at the preset's setpoints and the driver's initial valves.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// The valve commands to hand the plant on the next step.
    #[must_use]
    pub const fn valves(&self) -> &[f64; 12] {
        &self.valves
    }

    /// The scheme, for reading and moving setpoints.
    #[must_use]
    pub const fn scheme(&self) -> &Scheme {
        &self.scheme
    }

    /// The scheme, mutably.
    pub const fn scheme_mut(&mut self) -> &mut Scheme {
        &mut self.scheme
    }

    /// How many steps have run.
    #[must_use]
    pub const fn steps(&self) -> usize {
        self.step
    }

    /// The scenario as asked for, one-based as `IDV(n)` is.
    ///
    /// This is what the caller set, not what the plant is handed. See
    /// [`Driver::disturbances`].
    pub fn request_disturbance(&mut self, n: usize, on: bool) {
        self.requested[n - 1] = f64::from(u8::from(on));
        self.disturbances[n - 1] = self.requested[n - 1];
    }

    /// Set the whole scenario at once.
    pub fn request(&mut self, idv: &[f64; 20]) {
        self.requested = *idv;
        self.disturbances = *idv;
    }

    /// What the plant is actually handed on this step.
    ///
    /// Equal to the requested scenario except for `IDV(12)`, which the driver
    /// forces on at [`STEADY_STATE_STEPS`] unless
    /// [`DriverQuirks::only_the_requested_disturbances`] is set.
    #[must_use]
    pub const fn disturbances(&self) -> &[f64; 20] {
        &self.disturbances
    }

    /// Whether the driver has forced `IDV(12)` on beyond what was requested.
    #[must_use]
    pub fn scenario_is_overridden(&self) -> bool {
        // Bit comparison, not an epsilon: both are exactly 0.0 or exactly
        // 1.0, and the question is which.
        self.disturbances[11].to_bits() != self.requested[11].to_bits()
    }

    /// Run the controllers due on the next step, given the measurements the
    /// *previous* step produced.
    ///
    /// Returns the valve commands to integrate the plant with. Call
    /// [`Driver::settle`] afterwards, with the plant already advanced.
    ///
    /// `dt` is the plant step in hours.
    // @port temain_mod.f:366-394
    pub fn control(&mut self, previous: &[f64], dt: f64) -> &[f64; 12] {
        self.step += 1;
        // temain_mod.f:366-368. Unconditional, and ahead of the controllers.
        // @delta D-011 class=C temain_mod.f:366-368
        if !self.quirks.only_the_requested_disturbances && self.step >= STEADY_STATE_STEPS {
            self.disturbances[11] = 1.0;
        }
        self.scheme.step(self.step, previous, &mut self.valves, dt);
        &self.valves
    }

    /// `CONSHAND`, run after the plant has been advanced.
    // @port temain_mod.f:411
    pub fn settle(&mut self) {
        self.scheme.clamp(&mut self.valves);
    }

    /// [`Driver::control`] then [`Driver::settle`], for a caller that
    /// integrates in between and does not need the two separated.
    ///
    /// Advancing the plant between them is the *point*; this convenience form
    /// exists for tests that only exercise the control side, and for the
    /// doctest above.
    pub fn step(&mut self, previous: &[f64], dt: f64) -> [f64; 12] {
        let valves = *self.control(previous, dt);
        self.settle();
        valves
    }
}

#[cfg(test)]
mod scheme_tests {
    extern crate std;

    use super::*;
    use std::println;

    /// The fast loops first fire on step 3, not step 1.
    ///
    /// `MOD(I,3)` with `I` from 1. A zero-based loop testing `step % 3 == 0`
    /// fires on step 0 and shifts every controller action by two seconds for
    /// the whole run.
    #[test]
    fn the_fast_loops_first_fire_on_step_three() {
        let mut scheme = Scheme::new();
        let measurements = [1.0; 41];
        let mut valves = [50.0; 12];

        let bits = |v: &[f64; 12]| v.map(f64::to_bits);
        for step in 1..3 {
            let before = bits(&valves);
            scheme.step(step, &measurements, &mut valves, 1.0 / 3600.0);
            assert_eq!(bits(&valves), before, "step {step} should run nothing");
        }
        let before = bits(&valves);
        scheme.step(3, &measurements, &mut valves, 1.0 / 3600.0);
        assert_ne!(bits(&valves), before, "step 3 should run the fast loops");
    }

    /// The composition loops fire on step 360 and the quality loop on 900.
    #[test]
    fn the_slow_loops_fire_on_their_own_multiples() {
        let scheme = Scheme::new();
        let _ = scheme;
        for (period, first) in [
            (Period::Fast, 3),
            (Period::Composition, 360),
            (Period::Quality, 900),
        ] {
            assert_eq!(first % period.steps(), 0);
            assert_ne!((first - 1) % period.steps(), 0);
        }
        // 900 is a multiple of 3 but not of 360, so the quality loop and the
        // composition loops do not always coincide.
        assert_eq!(900 % 3, 0);
        assert_ne!(900 % 360, 0);
        // And 1080 is the first step where the fast and composition groups
        // both fire along with nothing else.
        assert_eq!(1080 % 360, 0);
    }

    /// `CONSHAND` clamps eleven valves and leaves the twelfth alone.
    #[test]
    fn the_clamp_leaves_the_agitator_valve_alone() {
        let scheme = Scheme::new();
        let mut valves = [-5.0; 12];
        valves[11] = -5.0;
        scheme.clamp(&mut valves);
        for (index, valve) in valves.iter().enumerate().take(11) {
            assert_eq!(valve.to_bits(), 0.0_f64.to_bits(), "valve {}", index + 1);
        }
        assert_eq!(
            valves[11].to_bits(),
            (-5.0_f64).to_bits(),
            "valve 12 was clamped; temain_mod.f:1401 loops I = 1, 11"
        );

        let mut valves = [150.0; 12];
        scheme.clamp(&mut valves);
        for valve in valves.iter().take(11) {
            assert_eq!(valve.to_bits(), 100.0_f64.to_bits());
        }
        assert_eq!(valves[11].to_bits(), 150.0_f64.to_bits());

        // The negative-zero normalisation the branch form exists for.
        let mut valves = [-0.0; 12];
        scheme.clamp(&mut valves);
        assert_eq!(
            valves[0].to_bits(),
            0.0_f64.to_bits(),
            "CONSHAND tests .LE. 0.0, so -0.0 becomes +0.0"
        );
        assert_eq!(
            valves[11].to_bits(),
            (-0.0_f64).to_bits(),
            "and valve 12 is not touched at all"
        );
    }

    /// The driver's initial valves are *not* `TEINIT`'s: they are rounded.
    /// The measured gaps, so D-009 quotes numbers rather than an impression.
    #[test]
    fn the_rounding_gaps_are_recorded() {
        let mut largest = (0.0_f64, 0);
        let mut smallest = (f64::INFINITY, 0);
        for (index, driver) in DRIVER_INITIAL_VALVES.iter().enumerate() {
            let teinit = tepsim_core::constants::NOMINAL_STATE[38 + index];
            let gap = (driver - teinit).abs();
            println!("XMV({:2}) gap {gap:.3e}", index + 1);
            if gap > largest.0 {
                largest = (gap, index + 1);
            }
            if gap > 0.0 && gap < smallest.0 {
                smallest = (gap, index + 1);
            }
        }
        println!(
            "largest {:.3e} at XMV({}), smallest non-zero {:.3e} at XMV({})",
            largest.0, largest.1, smallest.0, smallest.1
        );
        assert!(
            largest.0 < 1e-2,
            "XMV({}) is {:.3e} from TEINIT, which is too far to be a rounding \
             to five significant figures",
            largest.1,
            largest.0
        );
    }

    #[test]
    fn the_driver_starts_from_rounded_valve_positions() {
        use tepsim_core::constants::NOMINAL_STATE;
        let mut differing = 0;
        for (index, driver) in DRIVER_INITIAL_VALVES.iter().enumerate().take(11) {
            let teinit = NOMINAL_STATE[38 + index];
            assert!(
                (driver - teinit).abs() < 0.01,
                "valve {} differs by more than rounding: {driver} against {teinit}",
                index + 1
            );
            if driver.to_bits() != teinit.to_bits() {
                differing += 1;
            }
        }
        // Ten of the eleven, not all eleven. `XMV(5)` is the exception:
        // `TEINIT` leaves `YY(43) = 22.21000000`, which is already five
        // significant figures, so the driver's `22.210` rounds to the same
        // `f32`. That is what makes it clear the driver's values *are* the
        // `TEINIT` ones rounded, rather than an independent set of numbers.
        assert_eq!(
            differing, 10,
            "ten of the eleven should differ from TEINIT's values, with XMV(5) \
             the exception because TEINIT already had it at five significant \
             figures. If this changes, delta D-009 needs rereading."
        );
        assert_eq!(
            DRIVER_INITIAL_VALVES[4].to_bits(),
            NOMINAL_STATE[42].to_bits(),
            "XMV(5) is the one that rounds to itself"
        );
        // The twelfth is not overwritten and matches exactly.
        assert_eq!(
            DRIVER_INITIAL_VALVES[11].to_bits(),
            NOMINAL_STATE[49].to_bits(),
            "valve 12 should keep TEINIT's value"
        );
    }

    /// A cascade's outer and inner loop can fire on the same tick, and the
    /// outer one's effect is felt immediately because it runs first.
    #[test]
    fn a_cascade_takes_effect_within_the_same_tick() {
        // CONTRL17 (fast, writes SETPT(4)) runs before CONTRL4 would... except
        // it does not: FAST is [1,2,3,4,...,16,17,18], so CONTRL4 runs *before*
        // CONTRL17. The order is the source's and this pins it.
        let position = |n: usize| FAST.iter().position(|x| *x == n).expect("scheduled");
        assert!(
            position(4) < position(17),
            "CONTRL4 runs before CONTRL17, so a setpoint move by 17 is felt on \
             the *next* tick, not this one. That is the source order and \
             changing it changes the plant."
        );
        assert!(position(9) < position(16));
        assert!(position(10) < position(18));
    }
}
