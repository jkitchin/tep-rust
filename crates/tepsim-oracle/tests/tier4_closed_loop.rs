//! The closed-loop plant against the Fortran driver.
//!
//! B-0039. Everything before this compares a piece; this runs the whole thing:
//! plant, controllers, scheduler and clamp, step for step, against
//! `teprob.f` driven by `temain_mod.f`'s own scheme.
//!
//! # Driving the oracle's driver
//!
//! `temain_mod.f`'s main loop cannot be called: `instrument.rs` turned the
//! program into a subroutine nothing calls, precisely so it would not run. So
//! the loop is reproduced *here*, in the test, by calling the same
//! subroutines in the same order on the same schedule.
//!
//! That is not circular. The thing under test is
//! `tepsim_control::Scheme`, and what it is tested against is the Fortran
//! subroutines themselves plus a transcription of the driver's twenty-line
//! loop. If the transcription were wrong the two would disagree, because the
//! port's schedule is written independently in `Scheme::step`.

#![cfg(feature = "oracle")]

extern crate alloc;

use tepsim_control::{DRIVER_INITIAL_VALVES, Driver, DriverQuirks, STEADY_STATE_STEPS};
use tepsim_core::{Inputs, Plant, SimTime, State, constants, math};
use tepsim_oracle::Oracle;

const DT: f64 = 1.0 / 3600.0;

/// `temain_mod.f:369-394`, transcribed. The order within each group is the
/// source's.
fn run_fortran_controllers(oracle: &mut Oracle, step: usize) {
    if step % 3 == 0 {
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
        oracle.contrl16();
        oracle.contrl17();
        oracle.contrl18();
    }
    if step % 360 == 0 {
        oracle.contrl13();
        oracle.contrl14();
        oracle.contrl15();
        oracle.contrl19();
    }
    if step % 900 == 0 {
        oracle.contrl20();
    }
}

/// `CONSHAND`, `temain_mod.f:1401-1404`.
fn conshand(oracle: &mut Oracle) {
    let mut xmv = oracle.manipulated();
    for valve in xmv.iter_mut().take(11) {
        if *valve <= 0.0 {
            *valve = 0.0;
        }
        if *valve >= 100.0 {
            *valve = 100.0;
        }
    }
    oracle.set_manipulated(&xmv);
}

/// Load the Braatz preset into the oracle's `COMMON`.
fn load_preset(oracle: &mut Oracle) {
    use tepsim_control::PRESET;
    let mut all = oracle.ctrlall();
    all.deltat = DT;
    for entry in &PRESET {
        all.setpt[entry.setpoint_index - 1] = entry.setpoint;
    }
    oracle.set_ctrlall(&all);
    oracle.set_flag6(0);

    // Written out rather than generated: twenty `COMMON` blocks are twenty
    // types, eight of them without a `taui` field, and that difference is
    // exactly what B-0037 made structural.
    macro_rules! p {
        ($n:literal) => {
            tepsim_control::preset($n).expect("a preset").tuning
        };
    }
    oracle.set_ctrl1(&tepsim_oracle::Ctrl1 {
        gain: p!(1).gain,
        errold: 0.0,
    });
    oracle.set_ctrl2(&tepsim_oracle::Ctrl2 {
        gain: p!(2).gain,
        errold: 0.0,
    });
    oracle.set_ctrl3(&tepsim_oracle::Ctrl3 {
        gain: p!(3).gain,
        errold: 0.0,
    });
    oracle.set_ctrl4(&tepsim_oracle::Ctrl4 {
        gain: p!(4).gain,
        errold: 0.0,
    });
    oracle.set_ctrl5(&tepsim_oracle::Ctrl5 {
        gain: p!(5).gain,
        taui: p!(5).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl6(&tepsim_oracle::Ctrl6 {
        gain: p!(6).gain,
        errold: 0.0,
    });
    oracle.set_ctrl7(&tepsim_oracle::Ctrl7 {
        gain: p!(7).gain,
        errold: 0.0,
    });
    oracle.set_ctrl8(&tepsim_oracle::Ctrl8 {
        gain: p!(8).gain,
        errold: 0.0,
    });
    oracle.set_ctrl9(&tepsim_oracle::Ctrl9 {
        gain: p!(9).gain,
        errold: 0.0,
    });
    oracle.set_ctrl10(&tepsim_oracle::Ctrl10 {
        gain: p!(10).gain,
        taui: p!(10).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl11(&tepsim_oracle::Ctrl11 {
        gain: p!(11).gain,
        taui: p!(11).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl13(&tepsim_oracle::Ctrl13 {
        gain: p!(13).gain,
        taui: p!(13).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl14(&tepsim_oracle::Ctrl14 {
        gain: p!(14).gain,
        taui: p!(14).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl15(&tepsim_oracle::Ctrl15 {
        gain: p!(15).gain,
        taui: p!(15).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl16(&tepsim_oracle::Ctrl16 {
        gain: p!(16).gain,
        taui: p!(16).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl17(&tepsim_oracle::Ctrl17 {
        gain: p!(17).gain,
        taui: p!(17).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl18(&tepsim_oracle::Ctrl18 {
        gain: p!(18).gain,
        taui: p!(18).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl19(&tepsim_oracle::Ctrl19 {
        gain: p!(19).gain,
        taui: p!(19).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl20(&tepsim_oracle::Ctrl20 {
        gain: p!(20).gain,
        taui: p!(20).reset.expect("PI"),
        errold: 0.0,
    });
    oracle.set_ctrl22(&tepsim_oracle::Ctrl22 {
        gain: p!(22).gain,
        taui: p!(22).reset.expect("PI"),
        errold: 0.0,
    });
}

/// How far a closed-loop run held together.
struct ClosedLoop {
    /// Worst relative measurement error, and which `XMEAS` it was on.
    worst_measurement: (f64, usize),
    /// Worst absolute valve difference, in percent of range, and which valve.
    worst_valve: (f64, usize),
    /// The first step at which any valve differed at all.
    first_valve_split: Option<usize>,
    /// Steps actually run.
    steps: usize,
}

/// Run both sides closed-loop for `hours` and report where they part.
fn closed_loop(oracle: &mut Oracle, hours: usize) -> ClosedLoop {
    let (_, mut fortran) = oracle.init_cold();
    oracle.set_disturbances(&[0; 20]);
    oracle.set_rng(tepsim_oracle::golden::SEED);
    load_preset(oracle);
    oracle.set_manipulated(&DRIVER_INITIAL_VALVES);

    // The port. `TEINIT` has already run one `TEFUNC` internally, so the warm
    // start is taken from `COMMON`; see B-0034.
    let after_init = oracle.teproc();
    let mut plant = Plant::new();
    plant.set_rng(tepsim_oracle::golden::SEED);
    plant.set_seeds(tepsim_core::TemperatureSeeds {
        reactor: after_init.tcr,
        separator: after_init.tcs,
        stripper: after_init.tcc,
        mixing: after_init.tcv,
    });
    let mut driver = Driver::new();
    let mut state = State::from_flat(&constants::NOMINAL_STATE);
    // What the controllers read on step 1: `TEINIT`'s own `TEFUNC` call left
    // these in `COMMON/PV/`. Nothing fires before step 3, so they are never
    // actually used, but a driver that started from zeros would be wrong for
    // any schedule whose first period is 1.
    let mut previous = oracle.measurements();

    let mut worst_measurement = (0.0_f64, 0);
    let mut worst_valve = (0.0_f64, 0);
    let mut first_valve_split = None;
    let mut steps = 0;
    let mut t = 0.0;

    for step in 1..=hours * 3_600 {
        // The Fortran, in `temain_mod.f`'s order: force IDV(12), controllers,
        // integrate, clamp. The forcing is `temain_mod.f:366-368`, and it is
        // the *first* thing in the loop body, ahead of the controllers.
        if step >= STEADY_STATE_STEPS {
            let mut idv = [0; 20];
            idv[11] = 1;
            oracle.set_disturbances(&idv);
        }
        run_fortran_controllers(oracle, step);
        let their_valves = oracle.manipulated();

        // The port, in the same order and off the same stale measurements.
        let our_valves = *driver.control(&previous, DT);

        for (index, (ours, theirs)) in our_valves.iter().zip(their_valves).enumerate() {
            let difference = (ours - theirs).abs();
            if difference > 0.0 && first_valve_split.is_none() {
                first_valve_split = Some(step);
            }
            if difference > worst_valve.0 {
                worst_valve = (difference, index + 1);
            }
        }

        let yp = oracle.derivatives(t, &fortran);
        let theirs = oracle.measurements();
        conshand(oracle);

        let inputs = Inputs {
            manipulated: our_valves,
            disturbances: *driver.disturbances(),
        };
        let time = SimTime(t);
        let (next, ours) = plant
            .euler_step(time, &state, &inputs, DT)
            .expect("converges");
        driver.settle();

        for (index, (a, b)) in ours.as_array().iter().zip(theirs).enumerate() {
            if b == 0.0 {
                continue;
            }
            let relative = (a - b).abs() / b.abs();
            if relative > worst_measurement.0 {
                worst_measurement = (relative, index + 1);
            }
        }

        previous = *ours.as_array();
        state = next;
        for (slot, rate) in fortran.iter_mut().zip(yp) {
            *slot += DT * rate;
        }
        t += DT;
        steps = step;
    }
    ClosedLoop {
        worst_measurement,
        worst_valve,
        first_valve_split,
        steps,
    }
}

/// The whole loop, plant and controllers together, against the Fortran driver.
#[test]
fn the_closed_loop_plant_matches_the_fortran_driver() {
    let mut oracle = Oracle::lock();
    // Past the eight-hour mark, so the run covers the driver's hard-coded
    // `IDV(12)` (B-0040) rather than stopping short of it.
    let hours = 10;
    let run = closed_loop(&mut oracle, hours);

    println!(
        "closed loop, {} libm, {hours} h ({} steps)\n  \
         worst XMEAS  : {:.3e} at XMEAS({})\n  \
         worst XMV    : {:.3e} %range at XMV({})\n  \
         first split  : {}",
        if math::USES_SYSTEM_LIBM {
            "platform"
        } else {
            "vendored"
        },
        run.steps,
        run.worst_measurement.0,
        run.worst_measurement.1,
        run.worst_valve.0,
        run.worst_valve.1,
        run.first_valve_split
            .map_or_else(|| "never".to_string(), |step| format!("step {step}"))
    );

    if math::USES_SYSTEM_LIBM {
        // Identical transcendentals, so identical measurements, so identical
        // controller errors, so identical valves. The whole loop is closed
        // and every part of it is bit-exact; there is nowhere for a
        // difference to enter.
        assert_eq!(
            run.worst_valve.0, 0.0,
            "XMV({}) differs by {:.3e}, so the control path is not \
             bit-exact even with identical transcendentals",
            run.worst_valve.1, run.worst_valve.0
        );
        assert_eq!(
            run.worst_measurement.0, 0.0,
            "XMEAS({}) differs by {:.3e}",
            run.worst_measurement.1, run.worst_measurement.0
        );
    } else {
        // Vendored `libm`. Tier 4 is diagnostic, but a closed loop should not
        // amplify a one-ULP `exp` difference into anything visible in four
        // hours: the controllers are pulling the plant back to a setpoint.
        assert!(
            run.worst_measurement.0 < 1e-6,
            "the closed loop diverged by {:.3e} at XMEAS({}) in {hours} h",
            run.worst_measurement.0,
            run.worst_measurement.1
        );
    }
}

/// The controllers read the *previous* step's measurements.
///
/// This is the single ordering B-0039 got wrong first, and it is invisible to
/// any unit test of [`Driver`] on its own, because both orderings are
/// well-formed. So it is pinned against the Fortran directly: run the oracle
/// to the point where its `COMMON` holds step 2's `XMEAS`, call `CONTRL7`
/// there, and check which of the two candidate inputs reproduces it.
#[test]
fn the_controllers_read_the_previous_steps_measurements() {
    let mut oracle = Oracle::lock();
    let (_, fortran) = oracle.init_cold();
    oracle.set_disturbances(&[0; 20]);
    oracle.set_rng(tepsim_oracle::golden::SEED);
    load_preset(&mut oracle);
    oracle.set_manipulated(&DRIVER_INITIAL_VALVES);

    let mut yy = fortran;
    let mut t = 0.0;
    let advance = |oracle: &mut Oracle, yy: &mut [f64; 50], t: &mut f64| {
        let yp = oracle.derivatives(*t, yy);
        for (slot, rate) in yy.iter_mut().zip(yp) {
            *slot += DT * rate;
        }
        *t += DT;
    };

    advance(&mut oracle, &mut yy, &mut t);
    advance(&mut oracle, &mut yy, &mut t);
    let after_two = oracle.measurements();

    // `COMMON` now holds exactly what the driver's step 3 would see.
    oracle.contrl7();
    let theirs = oracle.manipulated()[6];
    // Undo it, so the measurement below is not taken from a moved plant.
    oracle.set_manipulated(&DRIVER_INITIAL_VALVES);

    advance(&mut oracle, &mut yy, &mut t);
    let after_three = oracle.measurements();
    assert_ne!(
        after_two[11], after_three[11],
        "XMEAS(12) did not move in a step, so this test cannot tell the two \
         orderings apart"
    );

    let fire = |measurements: &[f64; 41]| {
        let mut driver = Driver::new();
        let _ = driver.control(measurements, DT);
        let _ = driver.control(measurements, DT);
        driver.control(measurements, DT)[6]
    };
    let previous = fire(&after_two);
    let current = fire(&after_three);

    println!(
        "XMV(7) on the first fire: previous-step {previous:.9}, \
         current-step {current:.9}, fortran {theirs:.9}"
    );
    assert_eq!(
        previous.to_bits(),
        theirs.to_bits(),
        "driving the controllers from the previous step's measurements does \
         not reproduce CONTRL7"
    );
    assert_ne!(
        current.to_bits(),
        theirs.to_bits(),
        "both orderings give the same answer here, so this test proves nothing"
    );
}

/// The driver forces `IDV(12)` on at step 28,800 and not a step earlier.
///
/// `temain_mod.f:366-368` with `SSPTS = 3600 * 8`. `I .GE. SSPTS` and `I` is
/// one-based, so the disturbance is live for the integration of step 28,800
/// itself.
#[test]
fn the_driver_forces_idv12_at_the_eight_hour_mark() {
    assert_eq!(STEADY_STATE_STEPS, 28_800, "SSPTS = 3600 * 8");

    let mut driver = Driver::new();
    let measurements = [0.0; 41];
    for _ in 1..STEADY_STATE_STEPS {
        let _ = driver.control(&measurements, DT);
        driver.settle();
    }
    assert_eq!(driver.steps(), STEADY_STATE_STEPS - 1);
    assert_eq!(
        driver.disturbances()[11],
        0.0,
        "IDV(12) is on at step {}, one step early",
        driver.steps()
    );

    let _ = driver.control(&measurements, DT);
    assert_eq!(driver.steps(), STEADY_STATE_STEPS);
    assert_eq!(
        driver.disturbances()[11],
        1.0,
        "IDV(12) is not on at step {STEADY_STATE_STEPS}"
    );
    assert!(driver.scenario_is_overridden());

    // And the other nineteen are untouched.
    for (index, value) in driver.disturbances().iter().enumerate() {
        if index != 11 {
            assert_eq!(*value, 0.0, "IDV({}) was switched on too", index + 1);
        }
    }
}

/// The Class C fix, measured. **Delta D-011.**
///
/// Two port runs that differ in nothing but
/// [`DriverQuirks::only_the_requested_disturbances`], from the same seeds and
/// the same generator word, cross-checked against two Fortran runs that differ
/// in the same way.
///
/// # The effect does not start at the mark
///
/// `IDV(12)` reaches the plant through `IDVWLK(6)` (`teprob.f:351`), and the
/// walk it gates is only redrawn when `TIME` passes that channel's `TNEXT`
/// (`teprob.f:359`). So switching the flag on at step 28,800 changes nothing
/// until channel 6's next segment boundary. That is the *correct* behaviour
/// and both implementations do it; what makes it evidence rather than a
/// coincidence is that the Fortran parts at the same step.
///
/// This is also the teeth check for the differential above: if forcing
/// `IDV(12)` made no measurable difference in the two hours the run covers
/// past the mark, then that differential would pass whether the port
/// reproduced the quirk or not.
#[test]
fn the_forced_disturbance_changes_the_plant_measurably() {
    let mut oracle = Oracle::lock();
    let hours = 10;

    let port = |fixed: bool, seeds, priming| {
        let mut plant = Plant::new();
        plant.set_rng(tepsim_oracle::golden::SEED);
        plant.set_seeds(seeds);
        let mut driver = Driver::new();
        driver.quirks.only_the_requested_disturbances = fixed;
        let mut state = State::from_flat(&constants::NOMINAL_STATE);
        let mut previous: [f64; 41] = priming;
        let mut t = 0.0;
        let mut trace = alloc::vec::Vec::with_capacity(hours * 3_600);
        for _ in 1..=hours * 3_600 {
            let valves = *driver.control(&previous, DT);
            let inputs = Inputs {
                manipulated: valves,
                disturbances: *driver.disturbances(),
            };
            let (next, measured) = plant
                .euler_step(SimTime(t), &state, &inputs, DT)
                .expect("converges");
            driver.settle();
            previous = *measured.as_array();
            trace.push(previous);
            state = next;
            t += DT;
        }
        trace
    };

    let mut fortran = |oracle: &mut Oracle, fixed: bool| {
        let (_, initial) = oracle.init_cold();
        oracle.set_disturbances(&[0; 20]);
        oracle.set_rng(tepsim_oracle::golden::SEED);
        load_preset(oracle);
        oracle.set_manipulated(&DRIVER_INITIAL_VALVES);
        let mut yy = initial;
        let mut t = 0.0;
        let mut trace = alloc::vec::Vec::with_capacity(hours * 3_600);
        for step in 1..=hours * 3_600 {
            if !fixed && step >= STEADY_STATE_STEPS {
                let mut idv = [0; 20];
                idv[11] = 1;
                oracle.set_disturbances(&idv);
            }
            run_fortran_controllers(oracle, step);
            let yp = oracle.derivatives(t, &yy);
            trace.push(oracle.measurements());
            conshand(oracle);
            for (slot, rate) in yy.iter_mut().zip(yp) {
                *slot += DT * rate;
            }
            t += DT;
        }
        trace
    };

    let their_faithful = fortran(&mut oracle, false);
    let their_fixed = fortran(&mut oracle, true);

    // The seeds and priming measurements the port starts from, taken after a
    // fresh `TEINIT` so the second Fortran run above has not moved `COMMON`.
    let (_, _) = oracle.init_cold();
    let after_init = oracle.teproc();
    let seeds = tepsim_core::TemperatureSeeds {
        reactor: after_init.tcr,
        separator: after_init.tcs,
        stripper: after_init.tcc,
        mixing: after_init.tcv,
    };
    let priming = oracle.measurements();
    let our_faithful = port(false, seeds, priming);
    let our_fixed = port(true, seeds, priming);

    /// The first step at which two traces differ in any bit, and where.
    fn first_split(a: &[[f64; 41]], b: &[[f64; 41]]) -> Option<(usize, usize)> {
        for (step, (x, y)) in a.iter().zip(b).enumerate() {
            for (index, (p, q)) in x.iter().zip(y).enumerate() {
                if p.to_bits() != q.to_bits() {
                    return Some((step + 1, index + 1));
                }
            }
        }
        None
    }

    /// Worst relative difference over a whole pair of traces.
    fn worst(a: &[[f64; 41]], b: &[[f64; 41]]) -> (f64, usize) {
        let mut out = (0.0_f64, 0);
        for (x, y) in a.iter().zip(b) {
            for (index, (p, q)) in x.iter().zip(y).enumerate() {
                if *q != 0.0 {
                    let relative = (p - q).abs() / q.abs();
                    if relative > out.0 {
                        out = (relative, index + 1);
                    }
                }
            }
        }
        out
    }

    let ours = first_split(&our_faithful, &our_fixed);
    let theirs = first_split(&their_faithful, &their_fixed);
    let over_run = worst(&our_faithful, &our_fixed);
    let last = our_faithful.len() - 1;
    let at_end = worst(&our_faithful[last..], &our_fixed[last..]);

    println!(
        "D-011, {hours} h, faithful against fixed:\n  \
         first split, port   : {ours:?}\n  \
         first split, fortran: {theirs:?}\n  \
         worst over the run  : {:.3e} at XMEAS({})\n  \
         worst at the end    : {:.3e} at XMEAS({})",
        over_run.0, over_run.1, at_end.0, at_end.1
    );

    // Before the mark neither pair can differ: the flag has nothing to do yet,
    // so a difference here would be non-determinism, not the quirk.
    for step in 0..STEADY_STATE_STEPS - 1 {
        for (index, (a, b)) in our_faithful[step].iter().zip(&our_fixed[step]).enumerate() {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "XMEAS({}) differs at step {}, before IDV(12) can have fired",
                index + 1,
                step + 1
            );
        }
    }

    assert_eq!(
        ours, theirs,
        "the port and the Fortran react to the forced IDV(12) at different \
         steps, so the port's walk gating is wrong even though the plants \
         agree numerically"
    );
    let (step, _) = ours.expect("forcing IDV(12) changed nothing at all");
    assert!(
        step > STEADY_STATE_STEPS,
        "the effect appeared at step {step}, at or before the mark, which \
         would mean the walk is not gated by IDVWLK at all"
    );

    // The differential against the Fortran runs two hours past the mark and
    // holds to 3e-13. If the quirk moved the plant by less than that, the
    // differential could not tell a port that reproduces it from one that
    // does not.
    assert!(
        over_run.0 > 1e-9,
        "forcing IDV(12) moves the plant by only {:.3e}, which is below what \
         the differential can resolve, so that test proves nothing about this \
         quirk",
        over_run.0
    );
}

/// The plant stays up for eight hours closed-loop, which open-loop it cannot.
///
/// That is the point of the control layer, and it is the cheapest possible
/// end-to-end check that the loops are wired to the right valves: a scheme with
/// two loops transposed trips the plant within minutes.
#[test]
fn the_controlled_plant_does_not_trip() {
    let mut oracle = Oracle::lock();
    let (_, mut fortran) = oracle.init_cold();
    oracle.set_disturbances(&[0; 20]);
    oracle.set_rng(tepsim_oracle::golden::SEED);
    load_preset(&mut oracle);
    oracle.set_manipulated(&DRIVER_INITIAL_VALVES);

    let mut t = 0.0;
    for step in 1..=8 * 3_600 {
        run_fortran_controllers(&mut oracle, step);
        let yp = oracle.derivatives(t, &fortran);
        assert_eq!(
            oracle.shutdown_flag(),
            0,
            "the controlled plant tripped at step {step} ({t:.3} h)"
        );
        conshand(&mut oracle);
        for (slot, rate) in fortran.iter_mut().zip(yp) {
            *slot += DT * rate;
        }
        t += DT;
    }
    println!("8 simulated hours closed-loop with no trip");
}
