//! Driving a [`Simulation`] in chunks, and packing its samples into flat rows.
//!
//! Plain Rust with no JavaScript types in it, so the host-side tests in
//! `tests/` exercise every line of it without a browser. [`crate::sim`] is the
//! `wasm-bindgen` shell over this.
//!
//! # Why a chunked runner rather than [`Simulation::run`]
//!
//! [`Simulation::run`] is the right call for a native program: it consumes the
//! simulation, loops to the end, and hands back a [`tepsim::Run`]. A browser
//! cannot use it. A 48-hour run is 172,800 integrator steps, and whichever
//! thread makes that call is unavailable until every one of them is done: on
//! the main thread the page freezes, and in a worker the run cannot be stopped,
//! paused, or have a fault toggled part-way.
//!
//! [`Simulation::step`] is the seam that fixes it. It advances exactly one
//! integrator step and returns a [`Sample`] on the steps where one is due, so
//! this type can stop after any number of samples, hand them over, and let the
//! event loop run before continuing. [`Runner::step_chunk`] is that loop, and
//! it is the whole reason these bindings are not a thin wrapper over `run`.

use tepsim::run::{CHANNELS, Labels, Outcome, Sample};
use tepsim::scenario::DISTURBANCES;
use tepsim::{Scenario, Simulation};

use crate::digest::Fnv1a64;

/// Values in one packed sample row: the time, then the 53 channels.
///
/// `1 + CHANNELS`, taken from the facade rather than written out, so the stride
/// a browser slices buffers with cannot drift away from the row the simulator
/// produces.
pub const ROW_WIDTH: usize = 1 + CHANNELS;

/// The largest number of integrator steps one scenario may ask for.
///
/// A guard rail, not a modelling limit. 10 million steps is about 2,778
/// simulated hours at the original's one-second step, well past any run the
/// published datasets describe, and far enough short of "wait forever" that a
/// mistyped duration in a browser comes back as an error message instead of a
/// hung tab.
pub const MAX_STEPS: usize = 10_000_000;

/// The largest number of samples one scenario may record.
///
/// [`Runner::run_to_end`] materialises the whole run in one allocation. At this
/// cap that is 432 MB, which a 64-bit host tolerates and a 32-bit wasm heap
/// does not, so a browser hits an allocation failure rather than swapping.
pub const MAX_SAMPLES: usize = 1_000_000;

/// Why a [`Scenario`] was rejected before it could be run.
///
/// A closed set rather than a string, so the browser can render the message and
/// the tests can match on the cause.
///
/// This exists because [`Scenario`] is a plain struct of public fields with no
/// validation of its own, which is right for Rust callers: the type system and
/// the constructors keep them honest. A browser has neither. Every field
/// arrives as a JavaScript number, which may be `NaN`, `Infinity`, negative, or
/// zero, and `Scenario::samples` divides by `sample_every`. In wasm a panic is
/// an abort that takes the module with it, so these are checked up front.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConfigError {
    /// The seed was not finite, or was not positive. Zero is a fixed point of
    /// the generator's recurrence: every draw after it would also be zero.
    Seed,
    /// The duration was not finite, or was not positive.
    Hours,
    /// The integrator step was not finite, or was not positive.
    StepHours,
    /// The sampling cadence was zero. `Scenario::samples` divides by it.
    SampleEvery,
    /// The run works out to fewer than one integrator step.
    NoSteps,
    /// The run would record no samples: the cadence is longer than the run.
    NoSamples,
    /// The run would exceed [`MAX_STEPS`] or [`MAX_SAMPLES`].
    TooLong,
    /// A one-based `IDV` index outside `1..=20`.
    FaultIndex,
}

impl ConfigError {
    /// A message suitable for showing to a user.
    #[must_use]
    pub const fn message(self) -> &'static str {
        match self {
            ConfigError::Seed => "seed must be finite and greater than zero",
            ConfigError::Hours => "duration must be finite and greater than zero",
            ConfigError::StepHours => "integrator step must be finite and greater than zero",
            ConfigError::SampleEvery => "sample every must be at least one step",
            ConfigError::NoSteps => "the run works out to fewer than one integrator step",
            ConfigError::NoSamples => {
                "the run would record no samples: shorten the sampling interval"
            }
            ConfigError::TooLong => "the run exceeds the step or sample limit",
            ConfigError::FaultIndex => "fault index must be in 1..=20",
        }
    }
}

impl core::fmt::Display for ConfigError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(self.message())
    }
}

/// What a validated scenario works out to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Plan {
    /// Integrator steps the run will take.
    pub steps: usize,
    /// Samples it will record, if nothing ends it early.
    pub samples: usize,
}

/// Check a scenario before anything tries to run it.
///
/// # Errors
///
/// The first [`ConfigError`] that applies, in field order.
pub fn validate(scenario: &Scenario) -> Result<Plan, ConfigError> {
    if !scenario.seed.is_finite() || scenario.seed <= 0.0 {
        return Err(ConfigError::Seed);
    }
    if !scenario.hours.is_finite() || scenario.hours <= 0.0 {
        return Err(ConfigError::Hours);
    }
    if !scenario.step_hours.is_finite() || scenario.step_hours <= 0.0 {
        return Err(ConfigError::StepHours);
    }
    if scenario.sample_every == 0 {
        return Err(ConfigError::SampleEvery);
    }

    // Computed here rather than through `Scenario::steps`, which casts an `f64`
    // to `usize`. That cast saturates on infinity and gives zero for `NaN`,
    // both of which would turn a typo into a very long run instead of an error.
    let steps = (scenario.hours / scenario.step_hours).round();
    if !steps.is_finite() {
        return Err(ConfigError::NoSteps);
    }
    if steps < 1.0 {
        return Err(ConfigError::NoSteps);
    }
    if steps > MAX_STEPS as f64 {
        return Err(ConfigError::TooLong);
    }
    let steps = steps as usize;

    let samples = steps / scenario.sample_every;
    if samples == 0 {
        return Err(ConfigError::NoSamples);
    }
    if samples > MAX_SAMPLES {
        return Err(ConfigError::TooLong);
    }

    Ok(Plan { steps, samples })
}

/// Set or clear a one-based `IDV` flag on a scenario.
///
/// # Errors
///
/// [`ConfigError::FaultIndex`] if `id` is outside `1..=20`. `Scenario`'s own
/// `with_fault` panics instead, which is correct for a Rust caller and fatal in
/// wasm.
pub fn set_fault(scenario: &mut Scenario, id: u8, active: bool) -> Result<(), ConfigError> {
    let slot = usize::from(id)
        .checked_sub(1)
        .and_then(|i| scenario.disturbances.get_mut(i))
        .ok_or(ConfigError::FaultIndex)?;
    *slot = active;
    Ok(())
}

/// Whether a one-based `IDV` flag is set. `false` for an out-of-range index.
#[must_use]
pub fn fault(scenario: &Scenario, id: u8) -> bool {
    usize::from(id)
        .checked_sub(1)
        .and_then(|i| scenario.disturbances.get(i))
        .copied()
        .unwrap_or(false)
}

/// A digest of everything a run's output depends on.
///
/// Identifies the scenario, not the samples. Kept separate from
/// [`Runner::checksum`] on purpose: a scenario digest that changed because a
/// run got longer, or a run digest that changed because a label got edited,
/// would waste an afternoon apiece.
#[must_use]
pub fn scenario_digest(scenario: &Scenario) -> u64 {
    let mut hash = Fnv1a64::new();
    hash.write_f64(scenario.seed);
    hash.write_f64(scenario.hours);
    hash.write_f64(scenario.step_hours);
    hash.write_u64(scenario.sample_every as u64);
    for active in &scenario.disturbances {
        hash.write_bool(*active);
    }
    hash.write_bool(scenario.controlled);
    hash.write_bool(scenario.driver_forces_idv12);
    hash.write_bool(scenario.quirks.trip_ends_the_run);
    // The integrator changes every number in the run, so it has to be in here.
    // Hashed by name rather than by discriminant, because a name is stable
    // across a reordering of the enum and a discriminant is not.
    for byte in scenario.integrator.name().as_bytes() {
        hash.write_u8(*byte);
    }
    hash.finish()
}

/// A simulation being driven a chunk at a time.
#[derive(Clone, Debug)]
pub struct Runner {
    sim: Simulation,
    plan: Plan,
    emitted: usize,
    checksum: Fnv1a64,
    /// Ground truth as of the most recent sample. The browser needs it to say
    /// what was actually wrong with the plant while it was drawing the trace.
    labels: Labels,
    /// The scenario as most recently asked for, which can be ahead of the one
    /// running. See [`Runner::request_fault`].
    requested: Scenario,
}

impl Runner {
    /// Validate a scenario and start a run at step zero.
    ///
    /// # Errors
    ///
    /// A [`ConfigError`] if the scenario cannot produce a well-defined run.
    pub fn new(scenario: Scenario) -> Result<Self, ConfigError> {
        let plan = validate(&scenario)?;
        Ok(Self {
            sim: Simulation::new(scenario),
            plan,
            emitted: 0,
            checksum: Fnv1a64::new(),
            labels: Labels::none(),
            requested: scenario,
        })
    }

    /// What this run was asked to do.
    #[must_use]
    pub const fn scenario(&self) -> &Scenario {
        self.sim.scenario()
    }

    /// Steps and samples the scenario works out to.
    #[must_use]
    pub const fn plan(&self) -> Plan {
        self.plan
    }

    /// How many samples have been emitted.
    #[must_use]
    pub const fn emitted_samples(&self) -> usize {
        self.emitted
    }

    /// How many integrator steps have run.
    #[must_use]
    pub const fn steps_taken(&self) -> usize {
        self.sim.steps_taken()
    }

    /// Simulated time so far, in hours.
    #[must_use]
    pub const fn hours(&self) -> f64 {
        self.sim.hours()
    }

    /// Whether the run has stopped producing samples.
    ///
    /// True once every planned sample is out, and also when the plant tripped
    /// under [`tepsim_core::QuirkFixes::trip_ends_the_run`] or a temperature
    /// solve failed to converge. A caller that loops on this cannot spin.
    ///
    /// [`tepsim_core::QuirkFixes::trip_ends_the_run`]: tepsim::tepsim_core::QuirkFixes::trip_ends_the_run
    #[must_use]
    pub const fn is_finished(&self) -> bool {
        self.sim.is_halted() || self.emitted >= self.plan.samples
    }

    /// How the run ended, or `None` while it is still going.
    ///
    /// A trip appears here as soon as it happens even though the run continues
    /// past it, because `teprob.f:807-811` freezes the plant rather than
    /// stopping it and the frozen samples are part of the behaviour under test.
    #[must_use]
    pub const fn outcome(&self) -> Option<Outcome> {
        self.sim.outcome()
    }

    /// How the run ended, as a stable name: `completed`, `tripped`,
    /// `solve_failed`, or `None` while it is still going.
    ///
    /// A completed run reports `completed` once the last sample is out, without
    /// waiting for a further step: [`Simulation::outcome`] holds `None` until
    /// something goes wrong, which is right for it and unhelpful for a progress
    /// display.
    #[must_use]
    pub fn outcome_name(&self) -> Option<&'static str> {
        match self.outcome() {
            Some(Outcome::Tripped { .. }) => Some("tripped"),
            Some(Outcome::SolveFailed { .. }) => Some("solve_failed"),
            Some(Outcome::Completed) => Some("completed"),
            None if self.is_finished() => Some("completed"),
            None => None,
        }
    }

    /// Ground truth as of the most recent sample.
    #[must_use]
    pub const fn labels(&self) -> &Labels {
        &self.labels
    }

    /// Advance until `max_samples` samples are out or the run ends.
    ///
    /// Returns whatever completed, row-major, [`ROW_WIDTH`] values per row:
    /// fewer rows than asked for at the end of a run, an empty vector once it
    /// is over. The caller never has to work out how many are left.
    ///
    /// One sample costs `sample_every` integrator steps, so the chunk size is a
    /// latency control: it is how long the calling thread is unavailable. A
    /// worker that wants to answer a "stop" within a frame should pick a chunk
    /// it can finish in a frame.
    #[must_use]
    pub fn step_chunk(&mut self, max_samples: usize) -> Vec<f64> {
        let capacity = max_samples.min(self.remaining_samples()) * ROW_WIDTH;
        let mut out = Vec::with_capacity(capacity);
        let mut taken = 0;
        while taken < max_samples && !self.is_finished() {
            if let Some(sample) = self.sim.step() {
                self.push_row(&mut out, &sample);
                taken += 1;
            }
        }
        out
    }

    /// Run to the end in one call, returning every remaining row.
    ///
    /// The same numbers [`Simulation::run`] produces, in the same order,
    /// asserted in `tests/determinism.rs`. For tests and short runs. A worker
    /// that calls this cannot answer a message until the run is over, which is
    /// what [`Runner::step_chunk`] exists to avoid.
    #[must_use]
    pub fn run_to_end(&mut self) -> Vec<f64> {
        self.step_chunk(self.remaining_samples())
    }

    /// Set or clear `IDV(id)` while the run is in progress.
    ///
    /// # Errors
    ///
    /// [`ConfigError::FaultIndex`] if `id` is outside `1..=20`.
    ///
    /// # What this does and does not do
    ///
    /// Nothing, yet, to a run already under way. [`Simulation`] takes its
    /// scenario at construction and hands the disturbance vector to the driver
    /// there, so there is no seam in the facade for changing it mid-run and
    /// these bindings must not invent one: a disturbance that switched on
    /// through a path the native API does not have would produce a run no
    /// native caller could reproduce, which is exactly the reproducibility this
    /// project is built on.
    ///
    /// So this records the change and reports it through [`Runner::pending_restart`].
    /// The browser's honest move is to rebuild the run, which is one
    /// constructor call and, at a hundred times real time, quick enough to feel
    /// immediate. `www/worker.js` does that.
    pub fn request_fault(&mut self, id: u8, active: bool) -> Result<(), ConfigError> {
        let mut scenario = *self.sim.scenario();
        set_fault(&mut scenario, id, active)?;
        self.requested = scenario;
        Ok(())
    }

    /// Whether a requested fault change needs the run rebuilt to take effect.
    #[must_use]
    pub fn pending_restart(&self) -> bool {
        self.requested != *self.sim.scenario()
    }

    /// The scenario as requested, including changes not yet in effect.
    #[must_use]
    pub const fn requested_scenario(&self) -> &Scenario {
        &self.requested
    }

    /// Digest of every value emitted so far, in emission order.
    ///
    /// See [`crate::digest`] for what this is and is not.
    #[must_use]
    pub const fn checksum(&self) -> u64 {
        self.checksum.finish()
    }

    /// Samples still to come, if nothing ends the run early.
    #[must_use]
    const fn remaining_samples(&self) -> usize {
        self.plan.samples.saturating_sub(self.emitted)
    }

    /// Pack one sample as `[hours, XMEAS(1..41), XMV(1..12)]` and absorb it.
    ///
    /// The digest is taken over the packed row rather than over the sample it
    /// was packed from, and that is not a stylistic choice. Hashing the sample
    /// would leave the packing itself unchecked: a transposition or an off-by-one
    /// between here and the buffer would send wrong numbers to the browser under
    /// a digest that still matched the native run. Hashing what actually leaves
    /// the module closes that gap, and a deliberately corrupted `push_row` is
    /// how it was found.
    fn push_row(&mut self, out: &mut Vec<f64>, sample: &Sample) {
        let start = out.len();
        out.push(sample.hours);
        out.extend_from_slice(&sample.row());
        debug_assert_eq!(out.len() - start, ROW_WIDTH);
        self.checksum.write_slice(&out[start..]);
        self.labels = sample.labels;
        self.emitted += 1;
    }
}

/// The one-based faults active as of the most recent sample.
#[must_use]
pub fn active_faults(labels: &Labels) -> Vec<u8> {
    labels
        .faults()
        .filter_map(|id| u8::try_from(id).ok())
        .collect()
}

/// Hours since `IDV(id)` came on, or `None` if it never did.
///
/// Not simply the time since the run began: the driver switches `IDV(12)` on at
/// hour eight whatever the scenario asked for, so one onset can be later than
/// the others and is not the caller's doing. Delta D-011.
#[must_use]
pub fn hours_since_onset(labels: &Labels, id: u8) -> Option<f64> {
    let index = usize::from(id).checked_sub(1)?;
    if index >= DISTURBANCES {
        return None;
    }
    labels.since_onset[index]
}

/// The scenario the cross-platform self-check runs.
///
/// One hour of the fault-free plant, closed loop, from the seed `teprob.f`
/// ships with, under the faithful Euler integrator. 3,600 steps and 20
/// samples: long enough that the controllers have fired and the plant has moved
/// off its initial condition, short enough to run on page load without a pause.
///
/// Taken from [`tepsim::tier9::CASES`] rather than written out again, so the
/// number a browser prints is the number the library is committed to and not a
/// second scenario that merely looks the same.
#[must_use]
pub fn self_check_scenario() -> Scenario {
    tepsim::tier9::CASES[0].scenario()
}

/// Run [`self_check_scenario`] and return the digest of its output.
///
/// The point of the whole exercise. This number must be identical on x86-64,
/// aarch64 and wasm32, or the determinism invariant is already broken and every
/// trajectory built on it is suspect.
///
/// Computed through [`Runner`] rather than through [`tepsim::tier9::digest`],
/// on purpose: this is the browser's transport path, chunking and packing
/// included, and asserting that it lands on the library's committed constant is
/// worth more than calling the library and reporting what it said.
/// `tests/determinism.rs` makes that assertion.
#[must_use]
pub fn self_check_digest() -> u64 {
    let Ok(mut runner) = Runner::new(self_check_scenario()) else {
        // Unreachable: the scenario is a constant and `tests/determinism.rs`
        // asserts it validates. Returning zero rather than panicking keeps a
        // self-check from being the thing that aborts a browser tab.
        return 0;
    };
    let _ = runner.run_to_end();
    runner.checksum()
}

/// The self-check digest, for a WebAssembly runtime with no JavaScript glue.
///
/// `wasm-bindgen` generates the glue that makes the rest of this crate usable
/// from a browser, and generating it needs a tool. This export needs none: any
/// runtime that can instantiate the module and call an exported function can
/// read the digest back and compare it against a native run. That is what makes
/// the wasm half of Tier 9 runnable in continuous integration before the
/// browser app exists.
///
/// Named with a prefix because a `no_mangle` symbol is global to the linked
/// artifact.
//
// SAFETY: `no_mangle` is unsafe because the symbol could collide with another
// of the same name and silently redirect calls. The prefix makes a collision
// implausible. The function takes no arguments, touches no shared state, and
// returns a plain `u64`, so there is no aliasing or lifetime obligation for a
// caller to uphold.
#[unsafe(no_mangle)]
pub extern "C" fn tepsim_wasm_self_check_digest() -> u64 {
    self_check_digest()
}
