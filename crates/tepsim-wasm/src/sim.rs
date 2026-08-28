//! The JavaScript-facing classes: [`Scenario`], [`Sim`] and [`Fault`].
//!
//! Thin shells. The arithmetic is in [`tepsim`], the chunking and validation in
//! [`crate::runner`], which the host-side tests exercise without a browser.
//! This module moves values across the boundary and turns errors into
//! exceptions, and does nothing else.
//!
//! # Why `Float64Array` and not JSON
//!
//! A returned `Vec<f64>` becomes a JavaScript `Array` of boxed numbers, built
//! one element at a time. A returned `Float64Array` is one `ArrayBuffer` copied
//! out of the wasm heap in a single memcpy, and, more to the point, its buffer
//! is *transferable*: a worker hands it to the main thread with
//! `postMessage(chunk, [chunk.buffer])`, which moves ownership rather than
//! structured-cloning it. At 54 values a row and a 48-hour run, that is the
//! difference between a copy per chunk and a pointer per chunk.
//!
//! It is also why [`Sim::step_chunk`] returns a fresh array per call rather
//! than a view into wasm memory. A view would avoid even the memcpy, and it
//! would alias the heap, dangle the moment an allocation grew it, and be
//! untransferable. `SharedArrayBuffer` would avoid the copy properly and is
//! ruled out on purpose: it needs COOP and COEP response headers, which neither
//! GitHub Pages nor a Hugging Face Static Space can set, and free static
//! hosting is the reason the browser app exists (`PLAN.org`, "The browser
//! application").

use js_sys::Float64Array;
use wasm_bindgen::prelude::*;

use tepsim::Integrator;
use tepsim::run::Outcome;
use tepsim::tepsim_core::fault::{FAULTS, Shape};

use crate::channels;
use crate::digest::hex64;
use crate::runner::{self, ConfigError, Runner};

/// Render a configuration error as a JavaScript exception.
fn js_error(error: ConfigError) -> JsError {
    JsError::new(error.message())
}

/// One of the twenty `IDV` disturbances, as a value object.
///
/// Reading it does not touch a simulation. Everything on it comes from
/// `tepsim_core::FAULTS`, which `tepsim-oracle/tests/fault_table.rs` compares
/// against the original's header table and against the code that implements
/// each fault.
#[wasm_bindgen]
#[derive(Clone, Debug)]
pub struct Fault {
    id: u8,
    published: String,
    effect: String,
    kind: String,
    line: String,
    channels: Vec<u32>,
    valves: Vec<u32>,
    spiking: bool,
    affects_the_plant: bool,
}

#[wasm_bindgen]
impl Fault {
    /// The one-based `IDV` index, `1..=20`.
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn id(&self) -> u8 {
        self.id
    }

    /// The description from `teprob.f:172-191`, verbatim.
    ///
    /// Five of them say only "Unknown"; [`Fault::effect`] is what the source
    /// actually does in those cases.
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn published(&self) -> String {
        self.published.clone()
    }

    /// What the fault does, where the header is silent or vaguer than the code.
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn effect(&self) -> String {
        self.effect.clone()
    }

    /// How it enters: `step`, `random`, `spike` or `sticking`.
    ///
    /// `spike` is a random-variation fault whose channels are spike trains
    /// rather than walks. The header table calls both "Unknown" or "Random
    /// Variation" and the distinction is why three of them are reported in the
    /// literature as the hardest to detect: they are intermittent.
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }

    /// The `teprob.f` line or lines it acts on.
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn line(&self) -> String {
        self.line.clone()
    }

    /// The one-based walk channels it drives, for a random or spike fault.
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn channels(&self) -> Vec<u32> {
        self.channels.clone()
    }

    /// The one-based valves it sticks, for a sticking fault.
    #[wasm_bindgen(getter, js_name = valves)]
    #[must_use]
    pub fn valves(&self) -> Vec<u32> {
        self.valves.clone()
    }

    /// Whether its channels are spike trains rather than walks.
    #[wasm_bindgen(getter, js_name = spiking)]
    #[must_use]
    pub fn spiking(&self) -> bool {
        self.spiking
    }

    /// Whether the fault reaches the plant model at all.
    ///
    /// False for the three sticking faults, which widen a valve dead band and
    /// touch no equation. In an open-loop run they do nothing, and that is not
    /// a bug in the scenario.
    #[wasm_bindgen(getter, js_name = affectsThePlant)]
    #[must_use]
    pub fn affects_the_plant(&self) -> bool {
        self.affects_the_plant
    }

    /// A non-empty string for a fault panel: the published description, or the
    /// effect where the original says only "Unknown".
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn label(&self) -> String {
        if self.published == "Unknown" {
            self.effect.clone()
        } else {
            self.published.clone()
        }
    }
}

impl Fault {
    fn from_core(fault: &tepsim::tepsim_core::Fault) -> Option<Self> {
        let (kind, channels, valves, spiking) = match fault.shape {
            Shape::Step => ("step", Vec::new(), Vec::new(), false),
            Shape::Random { channels, spiking } => (
                if spiking { "spike" } else { "random" },
                channels
                    .iter()
                    .filter_map(|c| u32::try_from(*c).ok())
                    .collect(),
                Vec::new(),
                spiking,
            ),
            Shape::Sticking { valves } => (
                "sticking",
                Vec::new(),
                valves
                    .iter()
                    .filter_map(|v| u32::try_from(*v).ok())
                    .collect(),
                false,
            ),
        };
        Some(Self {
            id: u8::try_from(fault.index).ok()?,
            published: fault.published.to_string(),
            effect: fault.effect.to_string(),
            kind: kind.to_string(),
            line: fault.line.to_string(),
            channels,
            valves,
            spiking,
            affects_the_plant: fault.affects_the_plant(),
        })
    }
}

/// A complete, reproducible description of a run.
///
/// Mutable, because a settings panel edits it in place. [`Sim`] takes a copy at
/// construction, so editing a scenario never disturbs a simulation already
/// running.
///
/// Every field a run's output depends on is here, which is what makes a
/// recorded dataset reproducible from its scenario alone. [`Scenario::digest`]
/// is a 16-character summary of all of it.
#[wasm_bindgen]
#[derive(Clone, Copy, Debug, Default)]
pub struct Scenario {
    inner: tepsim::Scenario,
}

#[wasm_bindgen]
impl Scenario {
    /// The fault-free plant for 48 hours, closed loop, from the seed
    /// `teprob.f` ships with, integrated by the faithful Euler method.
    ///
    /// 48 hours is `NPTS = 172800` at a one-second step, the run
    /// `temain_mod.f` was written to do.
    #[wasm_bindgen(constructor)]
    #[must_use]
    pub fn new() -> Scenario {
        Scenario::default()
    }

    /// The baseline with one disturbance switched on.
    ///
    /// # Errors
    ///
    /// Throws if `id` is outside `1..=20`.
    #[wasm_bindgen(js_name = withFault)]
    pub fn with_fault(id: u8) -> Result<Scenario, JsError> {
        let mut scenario = Scenario::default();
        scenario.set_fault(id, true)?;
        Ok(scenario)
    }

    /// The generator word. `teprob.f:1187` compiles in 4651207995.
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn seed(&self) -> f64 {
        self.inner.seed
    }

    /// Set the generator word. Must be finite and positive.
    #[wasm_bindgen(setter)]
    pub fn set_seed(&mut self, seed: f64) {
        self.inner.seed = seed;
    }

    /// Run length in simulated hours.
    #[wasm_bindgen(getter, js_name = hours)]
    #[must_use]
    pub fn hours(&self) -> f64 {
        self.inner.hours
    }

    /// Set the run length in simulated hours.
    #[wasm_bindgen(setter, js_name = hours)]
    pub fn set_hours(&mut self, hours: f64) {
        self.inner.hours = hours;
    }

    /// Integrator step in hours. One second by default, as `INTGTR` uses.
    #[wasm_bindgen(getter, js_name = stepHours)]
    #[must_use]
    pub fn step_hours(&self) -> f64 {
        self.inner.step_hours
    }

    /// Set the integrator step in hours.
    ///
    /// Changing it changes every number in the run: the original's one-second
    /// fixed Euler step is part of what is being reproduced, not an
    /// implementation detail. Delta D-005.
    #[wasm_bindgen(setter, js_name = stepHours)]
    pub fn set_step_hours(&mut self, hours: f64) {
        self.inner.step_hours = hours;
    }

    /// Record one sample every this many integrator steps.
    ///
    /// 180 by default, which is the 180-second spacing `temain_mod.f:401`
    /// writes at and the spacing of the published `d00` to `d21` files.
    #[wasm_bindgen(getter, js_name = sampleEvery)]
    #[must_use]
    pub fn sample_every(&self) -> usize {
        self.inner.sample_every
    }

    /// Set the sampling cadence in integrator steps. Must be at least one.
    #[wasm_bindgen(setter, js_name = sampleEvery)]
    pub fn set_sample_every(&mut self, steps: usize) {
        self.inner.sample_every = steps;
    }

    /// Whether the plant runs closed loop under the published control scheme.
    ///
    /// Open loop is a diagnostic mode, not a useful operating one: the plant
    /// trips on reactor pressure after about three hours. It is here because
    /// the difference between the two is the clearest single statement of what
    /// the control layer does.
    #[wasm_bindgen(getter, js_name = controlled)]
    #[must_use]
    pub fn controlled(&self) -> bool {
        self.inner.controlled
    }

    /// Run closed loop, or open loop with the valves held.
    #[wasm_bindgen(setter, js_name = controlled)]
    pub fn set_controlled(&mut self, controlled: bool) {
        self.inner.controlled = controlled;
    }

    /// Whether the driver forces `IDV(12)` on at eight hours regardless of what
    /// was asked for, as `temain_mod.f:366-368` does.
    ///
    /// `true` by default, because every published dataset longer than eight
    /// hours carries it. Delta D-011. Turning it off is a departure from every
    /// number in the literature.
    #[wasm_bindgen(getter, js_name = driverForcesIdv12)]
    #[must_use]
    pub fn driver_forces_idv12(&self) -> bool {
        self.inner.driver_forces_idv12
    }

    /// Set whether the driver forces `IDV(12)` on at eight hours.
    #[wasm_bindgen(setter, js_name = driverForcesIdv12)]
    pub fn set_driver_forces_idv12(&mut self, forced: bool) {
        self.inner.driver_forces_idv12 = forced;
    }

    /// Whether a trip ends the run rather than freezing the plant.
    ///
    /// `false` by default, which is what `teprob.f:807-811` does: the plant
    /// freezes and keeps reporting, and those frozen samples are in every
    /// published dataset. Delta D-007, blocked on sign-off.
    #[wasm_bindgen(getter, js_name = tripEndsTheRun)]
    #[must_use]
    pub fn trip_ends_the_run(&self) -> bool {
        self.inner.quirks.trip_ends_the_run
    }

    /// Set whether a trip ends the run.
    #[wasm_bindgen(setter, js_name = tripEndsTheRun)]
    pub fn set_trip_ends_the_run(&mut self, ends: bool) {
        self.inner.quirks.trip_ends_the_run = ends;
    }

    /// The integrator: `euler`, `rk4` or `dopri5`.
    #[wasm_bindgen(getter, js_name = integrator)]
    #[must_use]
    pub fn integrator(&self) -> String {
        self.inner.integrator.name().to_string()
    }

    /// Choose the integrator by name.
    ///
    /// # Errors
    ///
    /// Throws if the name is not one of `euler`, `rk4`, `dopri5` or
    /// `dormand-prince`.
    #[wasm_bindgen(js_name = setIntegrator)]
    pub fn set_integrator(&mut self, name: &str) -> Result<(), JsError> {
        let integrator = Integrator::parse(name)
            .ok_or_else(|| JsError::new("unknown integrator: expected euler, rk4 or dopri5"))?;
        self.inner.integrator = integrator;
        Ok(())
    }

    /// Whether this scenario's integrator reproduces the original.
    ///
    /// Only Euler does. Everything the validation ladder claims is a claim
    /// about Euler, so a run using anything else is a *better* integration of
    /// the same equations and not a reproduction of the same numbers. A user
    /// interface that offers the choice should say so.
    #[wasm_bindgen(getter, js_name = isFaithful)]
    #[must_use]
    pub fn is_faithful(&self) -> bool {
        self.inner.integrator.is_faithful()
    }

    /// Whether `IDV(id)` is on. `false` for an out-of-range index.
    #[wasm_bindgen(js_name = fault)]
    #[must_use]
    pub fn fault(&self, id: u8) -> bool {
        runner::fault(&self.inner, id)
    }

    /// Set or clear `IDV(id)`.
    ///
    /// # Errors
    ///
    /// Throws if `id` is outside `1..=20`. `Scenario::with_fault` in the
    /// facade panics instead, which is right for a Rust caller and fatal in
    /// wasm, where a panic aborts the module.
    #[wasm_bindgen(js_name = setFault)]
    pub fn set_fault(&mut self, id: u8, active: bool) -> Result<(), JsError> {
        runner::set_fault(&mut self.inner, id, active).map_err(js_error)
    }

    /// The one-based indices of every requested disturbance.
    ///
    /// A getter, to match [`Sim::active_faults`]. That one is ground truth and
    /// this one is configuration, and they differ: the driver switches
    /// `IDV(12)` on at hour eight whatever was asked for.
    #[wasm_bindgen(getter, js_name = activeFaults)]
    #[must_use]
    pub fn active_faults(&self) -> Vec<u8> {
        self.inner
            .active_faults()
            .filter_map(|id| u8::try_from(id).ok())
            .collect()
    }

    /// Clear every disturbance flag.
    #[wasm_bindgen(js_name = clearFaults)]
    pub fn clear_faults(&mut self) {
        self.inner.disturbances = [false; tepsim::scenario::DISTURBANCES];
    }

    /// How many integrator steps this scenario is, or `0` if it is invalid.
    #[wasm_bindgen(getter, js_name = steps)]
    #[must_use]
    pub fn steps(&self) -> usize {
        runner::validate(&self.inner).map_or(0, |plan| plan.steps)
    }

    /// How many samples it will record, or `0` if it is invalid.
    ///
    /// Cheap enough to call on every keystroke in a settings panel, which is
    /// what it is for.
    #[wasm_bindgen(getter, js_name = sampleCount)]
    #[must_use]
    pub fn sample_count(&self) -> usize {
        runner::validate(&self.inner).map_or(0, |plan| plan.samples)
    }

    /// Why this scenario cannot be run, or `undefined` if it can.
    #[wasm_bindgen(js_name = validationError)]
    #[must_use]
    pub fn validation_error(&self) -> Option<String> {
        runner::validate(&self.inner)
            .err()
            .map(|e| e.message().to_string())
    }

    /// A 16-character hex digest of the scenario.
    ///
    /// Identifies the run's description, not its output. [`Sim::checksum`] is
    /// the digest of the samples. Two scenarios with the same digest produce
    /// the same numbers; that is the whole claim, and it is what makes a run
    /// shareable as a link rather than as a file.
    #[wasm_bindgen(js_name = digest)]
    #[must_use]
    pub fn digest(&self) -> String {
        hex64(runner::scenario_digest(&self.inner))
    }
}

/// A run in progress, advanced in chunks.
///
/// The shape this is built for is one `Sim` per Web Worker, driven by a loop
/// that calls [`Sim::step_chunk`], posts the result with its buffer in the
/// transfer list, and yields to the event loop so the worker can answer control
/// messages between chunks. `crates/tepsim-wasm/www/worker.js` is that loop, in
/// about forty lines.
#[wasm_bindgen]
#[derive(Clone, Debug)]
pub struct Sim {
    runner: Runner,
}

#[wasm_bindgen]
impl Sim {
    /// Validate a scenario and start a run at step zero.
    ///
    /// # Errors
    ///
    /// Throws if the scenario cannot produce a well-defined run: a seed that is
    /// not finite and positive, a non-positive duration or step, a zero
    /// sampling cadence, or a run past the step or sample limit.
    #[wasm_bindgen(constructor)]
    pub fn new(scenario: &Scenario) -> Result<Sim, JsError> {
        Runner::new(scenario.inner)
            .map(|runner| Sim { runner })
            .map_err(js_error)
    }

    /// The scenario this run is executing.
    #[wasm_bindgen(getter, js_name = scenario)]
    #[must_use]
    pub fn scenario(&self) -> Scenario {
        Scenario {
            inner: *self.runner.scenario(),
        }
    }

    /// Values per packed row: the time, then the 53 channels.
    ///
    /// The stride into the `Float64Array` [`Sim::step_chunk`] returns.
    #[wasm_bindgen(getter, js_name = rowWidth)]
    #[must_use]
    pub fn row_width(&self) -> usize {
        crate::runner::ROW_WIDTH
    }

    /// Short stable identifiers for the columns, in order.
    #[wasm_bindgen(getter, js_name = columnIds)]
    #[must_use]
    pub fn column_ids(&self) -> Vec<String> {
        channels::column_ids()
    }

    /// Human-readable labels for the columns, in order.
    #[wasm_bindgen(getter, js_name = columnLabels)]
    #[must_use]
    pub fn column_labels(&self) -> Vec<String> {
        channels::column_labels()
    }

    /// How many samples the run will record if nothing ends it early.
    #[wasm_bindgen(getter, js_name = totalSamples)]
    #[must_use]
    pub fn total_samples(&self) -> usize {
        self.runner.plan().samples
    }

    /// How many integrator steps the run will take.
    #[wasm_bindgen(getter, js_name = totalSteps)]
    #[must_use]
    pub fn total_steps(&self) -> usize {
        self.runner.plan().steps
    }

    /// How many samples have been emitted.
    #[wasm_bindgen(getter, js_name = emittedSamples)]
    #[must_use]
    pub fn emitted_samples(&self) -> usize {
        self.runner.emitted_samples()
    }

    /// How many integrator steps have run.
    #[wasm_bindgen(getter, js_name = stepsTaken)]
    #[must_use]
    pub fn steps_taken(&self) -> usize {
        self.runner.steps_taken()
    }

    /// Simulated time so far, in hours.
    #[wasm_bindgen(getter, js_name = hours)]
    #[must_use]
    pub fn hours(&self) -> f64 {
        self.runner.hours()
    }

    /// Fraction of the planned samples emitted, in `0.0..=1.0`.
    #[wasm_bindgen(getter, js_name = progress)]
    #[must_use]
    pub fn progress(&self) -> f64 {
        let total = self.runner.plan().samples;
        if total == 0 {
            1.0
        } else {
            self.runner.emitted_samples() as f64 / total as f64
        }
    }

    /// Whether the run has stopped producing samples.
    #[wasm_bindgen(getter, js_name = isFinished)]
    #[must_use]
    pub fn is_finished(&self) -> bool {
        self.runner.is_finished()
    }

    /// How the run ended: `completed`, `tripped`, `solve_failed`, or
    /// `undefined` while it is still going.
    ///
    /// A trip appears here as soon as it happens even though the run continues
    /// past it: `teprob.f:807-811` freezes the plant rather than stopping it,
    /// and the frozen samples are part of the behaviour under test.
    #[wasm_bindgen(getter, js_name = outcome)]
    #[must_use]
    pub fn outcome(&self) -> Option<String> {
        self.runner.outcome_name().map(str::to_string)
    }

    /// Simulated hours at which the plant tripped, or `undefined`.
    #[wasm_bindgen(getter, js_name = tripHours)]
    #[must_use]
    pub fn trip_hours(&self) -> Option<f64> {
        match self.runner.outcome() {
            Some(Outcome::Tripped { hours, .. }) => Some(hours),
            _ => None,
        }
    }

    /// The first shutdown condition that fired, or `undefined`.
    ///
    /// One of the eight `teprob.f:703-710` tests, as prose.
    #[wasm_bindgen(getter, js_name = tripCause)]
    #[must_use]
    pub fn trip_cause(&self) -> Option<String> {
        match self.runner.outcome() {
            Some(Outcome::Tripped { cause, .. }) => cause.map(|c| c.describe().to_string()),
            _ => None,
        }
    }

    /// Ground truth: the disturbances active as of the most recent sample.
    ///
    /// Not the same as the scenario's list. The driver switches `IDV(12)` on at
    /// eight hours whatever was asked for, so this list can grow mid-run
    /// without the caller doing anything. Delta D-011.
    #[wasm_bindgen(getter, js_name = activeFaults)]
    #[must_use]
    pub fn active_faults(&self) -> Vec<u8> {
        runner::active_faults(self.runner.labels())
    }

    /// Hours since `IDV(id)` came on, or `undefined` if it never did.
    ///
    /// The number a detection-delay figure should be measured against. The
    /// original records nothing of the sort, which is why delays in the
    /// literature are computed against whatever onset the author assumed.
    #[wasm_bindgen(js_name = hoursSinceOnset)]
    #[must_use]
    pub fn hours_since_onset(&self, id: u8) -> Option<f64> {
        runner::hours_since_onset(self.runner.labels(), id)
    }

    /// Advance the run by at most `max_samples` samples.
    ///
    /// Returns a fresh `Float64Array`, row-major, `rowWidth` values per row,
    /// holding **whatever completed**: fewer rows than asked for at the end of
    /// a run, an empty array once it is over. The caller never has to work out
    /// how many are left.
    ///
    /// One sample costs `sampleEvery` integrator steps, so the chunk size is a
    /// latency control: it is how long this thread is unavailable. A worker
    /// that wants to answer a "stop" within a frame should pick a chunk it can
    /// finish in a frame.
    ///
    /// The returned array owns its `ArrayBuffer`, so it goes straight into a
    /// `postMessage` transfer list. After transfer it is detached on this side,
    /// which is correct: this side is done with it.
    #[wasm_bindgen(js_name = stepChunk)]
    #[must_use]
    pub fn step_chunk(&mut self, max_samples: usize) -> Float64Array {
        Float64Array::from(self.runner.step_chunk(max_samples).as_slice())
    }

    /// Run to the end and return every remaining row in one array.
    ///
    /// The same numbers `tepsim::Simulation::run` produces, in the same order.
    /// For tests and short runs: a worker that calls this cannot answer a
    /// message until the run is over, which is what [`Sim::step_chunk`] exists
    /// to avoid.
    #[wasm_bindgen(js_name = runToEnd)]
    #[must_use]
    pub fn run_to_end(&mut self) -> Float64Array {
        Float64Array::from(self.runner.run_to_end().as_slice())
    }

    /// Ask for `IDV(id)` to be set or cleared.
    ///
    /// # Errors
    ///
    /// Throws if `id` is outside `1..=20`.
    ///
    /// # This does not change a run already under way
    ///
    /// `tepsim::Simulation` takes its scenario at construction and hands the
    /// disturbance vector to the driver there. There is no seam in the facade
    /// for changing it mid-run, and these bindings must not invent one: a
    /// disturbance that switched on through a path the native API does not have
    /// would produce a run no native caller could reproduce, and reproducibility
    /// is the entire point of this project.
    ///
    /// So the request is recorded and [`Sim::pending_restart`] goes `true`. The
    /// browser's honest move is to rebuild: `new Sim(sim.requestedScenario)`.
    /// At a hundred times real time that is quick enough to feel immediate, and
    /// the run it produces is one a native caller can reproduce exactly.
    #[wasm_bindgen(js_name = setFault)]
    pub fn set_fault(&mut self, id: u8, active: bool) -> Result<(), JsError> {
        self.runner.request_fault(id, active).map_err(js_error)
    }

    /// Whether a requested change needs the run rebuilt to take effect.
    #[wasm_bindgen(getter, js_name = pendingRestart)]
    #[must_use]
    pub fn pending_restart(&self) -> bool {
        self.runner.pending_restart()
    }

    /// The scenario as requested, including changes not yet in effect.
    ///
    /// Pass it to the [`Sim`] constructor to apply them.
    #[wasm_bindgen(getter, js_name = requestedScenario)]
    #[must_use]
    pub fn requested_scenario(&self) -> Scenario {
        Scenario {
            inner: *self.runner.requested_scenario(),
        }
    }

    /// A 16-character hex digest of every value emitted so far.
    ///
    /// FNV-1a over IEEE 754 bit patterns, integer arithmetic throughout, so the
    /// digest cannot itself become a source of the cross-platform disagreement
    /// it exists to detect. Print it beside a native run of the same scenario:
    /// equal digests mean the wasm build reproduced the native one bit for bit.
    #[wasm_bindgen(js_name = checksum)]
    #[must_use]
    pub fn checksum(&self) -> String {
        hex64(self.runner.checksum())
    }
}

/// The twenty `IDV` disturbances, from `tepsim_core::FAULTS`.
#[wasm_bindgen(js_name = faults)]
#[must_use]
pub fn faults() -> Vec<Fault> {
    FAULTS.iter().filter_map(Fault::from_core).collect()
}

/// One disturbance by its one-based index, or `undefined` if out of range.
#[wasm_bindgen(js_name = faultInfo)]
#[must_use]
pub fn fault_info(id: u8) -> Option<Fault> {
    tepsim::tepsim_core::fault::fault(usize::from(id)).and_then(Fault::from_core)
}

/// How many disturbances this model has.
///
/// Twenty, not the twenty-one of the later literature: `teprob.f:340` loops
/// `DO 500 I=1,20`.
#[wasm_bindgen(js_name = faultCount)]
#[must_use]
pub fn fault_count() -> usize {
    tepsim::scenario::DISTURBANCES
}

/// Short stable identifiers for the columns of a packed row.
#[wasm_bindgen(js_name = columnIds)]
#[must_use]
pub fn column_ids() -> Vec<String> {
    channels::column_ids()
}

/// Human-readable labels for the columns of a packed row.
#[wasm_bindgen(js_name = columnLabels)]
#[must_use]
pub fn column_labels() -> Vec<String> {
    channels::column_labels()
}

/// Units for the columns of a packed row.
#[wasm_bindgen(js_name = columnUnits)]
#[must_use]
pub fn column_units() -> Vec<String> {
    channels::column_units()
}

/// Zero-based offsets of the analyser-sampled measurements, `XMEAS(23..=41)`.
///
/// These hold their value between analyser reports, so a chart should draw them
/// as steps rather than interpolating.
#[wasm_bindgen(js_name = sampledColumns)]
#[must_use]
pub fn sampled_columns() -> Vec<u32> {
    channels::sampled_columns()
}

/// Values per packed row: the time, then the 53 channels.
#[wasm_bindgen(js_name = rowWidth)]
#[must_use]
pub fn row_width() -> usize {
    crate::runner::ROW_WIDTH
}

/// How many measurements the plant reports. `XMEAS(41)`.
#[wasm_bindgen(js_name = measurementCount)]
#[must_use]
pub fn measurement_count() -> usize {
    channels::measurement_count()
}

/// How many manipulated variables the plant accepts. `XMV(12)`.
#[wasm_bindgen(js_name = manipulatedCount)]
#[must_use]
pub fn manipulated_count() -> usize {
    channels::manipulated_count()
}

/// The integrators a scenario can choose: `euler`, `rk4`, `dopri5`.
#[wasm_bindgen(js_name = integrators)]
#[must_use]
pub fn integrators() -> Vec<String> {
    [
        Integrator::Euler,
        Integrator::Rk4,
        Integrator::DormandPrince,
    ]
    .iter()
    .map(|i| i.name().to_string())
    .collect()
}

/// The largest number of integrator steps one scenario may ask for.
#[wasm_bindgen(js_name = maxSteps)]
#[must_use]
pub fn max_steps() -> usize {
    crate::runner::MAX_STEPS
}

/// The largest number of samples one scenario may record.
#[wasm_bindgen(js_name = maxSamples)]
#[must_use]
pub fn max_samples() -> usize {
    crate::runner::MAX_SAMPLES
}

/// The crate version, from `Cargo.toml`.
#[wasm_bindgen(js_name = version)]
#[must_use]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// The algorithm behind [`Sim::checksum`] and [`Scenario::digest`].
#[wasm_bindgen(js_name = checksumAlgorithm)]
#[must_use]
pub fn checksum_algorithm() -> String {
    "fnv1a-64".to_string()
}

/// Run a fixed one-hour baseline scenario and return the digest of its output.
///
/// The cross-platform determinism check, reduced to one string. A browser that
/// prints a different value from a native run of the same commit has broken the
/// invariant the whole validation ladder rests on, and the page should say so
/// loudly rather than plot the numbers.
///
/// 3,600 integrator steps and 20 samples: quick enough to run on page load.
#[wasm_bindgen(js_name = selfCheckDigest)]
#[must_use]
pub fn self_check_digest() -> String {
    hex64(crate::runner::self_check_digest())
}
