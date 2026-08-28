//! The run loop.

use alloc::vec::Vec;

use tepsim_control::{DRIVER_INITIAL_VALVES, Driver, STEADY_STATE_STEPS};
use tepsim_core::{Inputs, Plant, SimTime, State, TemperatureSeeds, constants};

use crate::run::{Labels, Outcome, Run, Sample};
use crate::scenario::{DISTURBANCES, Scenario};

/// A simulation in progress.
///
/// Owns everything it needs, is `Send` and `Clone`, and performs no I/O. The
/// original keeps its whole working set in six `COMMON` blocks, which allows
/// exactly one simulation per process and no reentrancy; this allows as many as
/// there are threads.
///
/// # Stepping
///
/// [`Simulation::run`] is the whole scenario in one call. [`Simulation::step`]
/// advances one integrator step and hands back a sample when one is due, which
/// is what an online detector, a live browser display or a reinforcement
/// learning loop needs.
///
/// # The loop's order is the original's
///
/// `temain_mod.f`'s main loop calls the controllers, *then* integrates, *then*
/// clamps. So the controllers on step `I` read the measurements step `I - 1`
/// produced: one plant step of dead time in every loop, present in no
/// controller. Delta D-010, and getting it backwards leaves `XMEAS(14)` 23%
/// out after four hours.
///
/// # Example
///
/// ```
/// use tepsim::{Scenario, Simulation};
///
/// let scenario = Scenario::baseline().with_hours(0.5);
/// let run = Simulation::new(scenario).run();
///
/// assert!(run.outcome.is_completed());
/// assert_eq!(run.samples.len(), 10);
///
/// // Reactor pressure sits near 2705 kPa at the nominal operating point.
/// let pressure = run.measurement(7);
/// assert!(pressure.iter().all(|p| (2600.0..2800.0).contains(p)));
/// ```
#[derive(Clone, Debug)]
pub struct Simulation {
    scenario: Scenario,
    plant: Plant,
    driver: Driver,
    state: State,
    /// What the controllers read on the next step: the measurements the last
    /// one produced. See the type documentation.
    previous: [f64; 41],
    /// One-based, as `temain_mod.f`'s `I` is.
    step: usize,
    /// Simulated time in hours, accumulated by addition rather than computed
    /// from the step number, because that is what the original does and the
    /// two differ in the last bits.
    hours: f64,
    /// Hours at which each disturbance came on.
    onset: [Option<f64>; DISTURBANCES],
    /// The first terminal event, if one has happened. Recording a trip does
    /// *not* stop the run: `teprob.f:807-811` freezes the plant and the frozen
    /// samples are part of the behaviour under test.
    event: Option<Outcome>,
    /// Whether to stop producing samples.
    halted: bool,
}

impl Simulation {
    /// A simulation ready to run the given scenario.
    #[must_use]
    pub fn new(scenario: Scenario) -> Self {
        let mut plant = Plant::new();
        plant.set_rng(scenario.seed);
        // Not `TemperatureSeeds::default()`: `TEINIT` runs one evaluation of
        // its own before returning, so a run starts from the converged values
        // rather than from the nominal literals. See B-0017 and B-0034.
        plant.set_seeds(TemperatureSeeds::after_initialisation());
        plant.quirks = scenario.quirks;

        let mut driver = Driver::new();
        driver.request(&scenario.disturbance_vector());

        let mut onset = [None; DISTURBANCES];
        for (slot, on) in onset.iter_mut().zip(scenario.disturbances) {
            if on {
                *slot = Some(0.0);
            }
        }

        Self {
            scenario,
            plant,
            driver,
            state: State::from_flat(&constants::NOMINAL_STATE),
            // Zeros, not `TEINIT`'s measurements. Nothing reads this: the
            // fastest control loop is `MOD(I,3)`, which first fires on step 3,
            // by which time step 2 has produced real measurements. Asserted in
            // `the_priming_measurements_are_never_read`.
            previous: [0.0; 41],
            step: 0,
            hours: 0.0,
            onset,
            event: None,
            halted: false,
        }
    }

    /// What this simulation is running.
    #[must_use]
    pub const fn scenario(&self) -> &Scenario {
        &self.scenario
    }

    /// How many steps have run.
    #[must_use]
    pub const fn steps_taken(&self) -> usize {
        self.step
    }

    /// Simulated time so far, in hours.
    #[must_use]
    pub const fn hours(&self) -> f64 {
        self.hours
    }

    /// The integrable state, for a caller that wants to inspect or checkpoint
    /// it.
    #[must_use]
    pub const fn state(&self) -> &State {
        &self.state
    }

    /// How the run ended, or `None` while it is still going.
    ///
    /// A trip appears here as soon as it happens, even though the run
    /// continues past it.
    #[must_use]
    pub const fn outcome(&self) -> Option<Outcome> {
        self.event
    }

    /// Whether the run has stopped producing samples.
    #[must_use]
    pub const fn is_halted(&self) -> bool {
        self.halted
    }

    /// Advance one integrator step.
    ///
    /// Returns a [`Sample`] on the steps where one is due and `None`
    /// otherwise, so a caller can drive the loop without deciding when to
    /// record. Returns `None` once the run has finished.
    ///
    /// The order is `temain_mod.f`'s: force `IDV(12)` if due, run the
    /// controllers on the *previous* step's measurements, integrate, then
    /// clamp.
    pub fn step(&mut self) -> Option<Sample> {
        if self.halted || self.step >= self.scenario.steps() {
            self.halted = true;
            return None;
        }
        self.step += 1;
        let step = self.step;
        let dt = self.scenario.step_hours;

        // `temain_mod.f:366-368`. The driver switches IDV(12) on at eight
        // hours whatever the scenario asked for. Delta D-011.
        self.driver.quirks.only_the_requested_disturbances = !self.scenario.driver_forces_idv12;

        let valves = if self.scenario.controlled {
            *self.driver.control(&self.previous, dt)
        } else {
            // Open loop: the driver still advances, so `IDV(12)` still fires
            // on schedule, but nothing moves the valves.
            let _ = self.driver.control(&self.previous, dt);
            DRIVER_INITIAL_VALVES
        };

        // Record the onset of anything now active that was not before,
        // including whatever the driver switched on unasked.
        for (index, active) in self.driver.disturbances().iter().enumerate() {
            if *active != 0.0 && self.onset[index].is_none() {
                self.onset[index] = Some(self.hours);
            }
        }

        let inputs = Inputs {
            manipulated: valves,
            disturbances: *self.driver.disturbances(),
        };
        let time = SimTime(self.hours);

        self.plant.advance_discrete(time, &inputs);
        let Ok((derivative, signals)) = self.plant.derivatives(time, &self.state, &inputs) else {
            // A solve failure is not a trip and not a result. The original
            // cannot report it at all: `TESUB2` returns its guess and claims
            // success (delta D-001). Stopping is the only honest answer.
            self.event = Some(Outcome::SolveFailed { step });
            self.halted = true;
            return None;
        };

        if signals.shutdown.is_tripped() && self.event.is_none() {
            self.event = Some(Outcome::Tripped {
                step,
                hours: self.hours,
                cause: signals.shutdown.first(),
            });
            // Whether it *ends* the run is delta D-007, off by default. With
            // the default the plant freezes and keeps reporting, which is what
            // the original does and what every published dataset contains.
            self.halted = self.scenario.quirks.trip_ends_the_run;
        }

        let measured = self.plant.sample_measurements(time, &signals);
        self.previous = *measured.as_array();

        if self.plant.step_seeds(&self.state).is_err() {
            self.event = Some(Outcome::SolveFailed { step });
            self.halted = true;
            return None;
        }
        self.state = self.state.step(dt, &derivative);
        self.driver.settle();

        let due = step % self.scenario.sample_every == 0;
        let sample = due.then(|| Sample {
            step,
            hours: self.hours,
            measurements: self.previous,
            manipulated: valves,
            labels: self.labels(),
        });

        self.hours += dt;
        sample
    }

    /// Run the whole scenario and collect every sample.
    ///
    /// # Panics
    ///
    /// Never. A plant that trips or fails to converge is reported through
    /// [`Run::outcome`], not by panicking: a run that ended early is data, and
    /// throwing it away would hide the difference between a port that trips
    /// where the original does and one that does not.
    #[must_use]
    pub fn run(mut self) -> Run {
        let mut samples = Vec::with_capacity(self.scenario.samples());
        while !self.halted {
            if let Some(sample) = self.step() {
                samples.push(sample);
            }
        }
        Run {
            scenario: self.scenario,
            samples,
            outcome: self.event.unwrap_or(Outcome::Completed),
        }
    }

    /// Ground truth as of now.
    fn labels(&self) -> Labels {
        let active = core::array::from_fn(|i| self.driver.disturbances()[i] != 0.0);
        let since_onset = core::array::from_fn(|i| self.onset[i].map(|at| self.hours - at));
        Labels {
            active,
            since_onset,
        }
    }
}

/// The step at which the driver forces `IDV(12)` on, for a caller that wants to
/// reason about it. `temain_mod.f:226`.
#[must_use]
pub const fn forced_disturbance_step() -> usize {
    STEADY_STATE_STEPS
}
