//! What to simulate: how long, from which seed, with which disturbances.

use tepsim_core::{Extensions, FAULTS, QuirkFixes};
use tepsim_scenario::{Invalid, Schedule};

use crate::integrator::Integrator;

/// One second, in hours. The step the original's `INTGTR` uses.
pub const DEFAULT_STEP_HOURS: f64 = 1.0 / 3600.0;

/// Output cadence, in steps. `temain_mod.f:401` writes every 180 seconds, and
/// the published `d00`-`d21` files are at that spacing.
pub const DEFAULT_SAMPLE_EVERY: usize = 180;

/// The generator word compiled into `teprob.f:1187`.
pub const DEFAULT_SEED: f64 = 4_651_207_995.0;

/// How many disturbances this model has.
///
/// Twenty, not the twenty-one of the later literature: `teprob.f:340` loops
/// `DO 500 I=1,20`. See [`tepsim_core::FAULTS`].
pub const DISTURBANCES: usize = FAULTS.len();

/// A complete description of a run.
///
/// Cheap to clone and to compare, so a caller can build one, tweak it, and
/// keep both. Everything a run's output depends on is in here, which is what
/// makes a recorded dataset reproducible from its scenario alone.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Scenario {
    /// The generator word. Any positive integer; odd values give the
    /// generator its full period. See `tepsim_core::TepRng`.
    pub seed: f64,
    /// How long to run, in hours.
    pub hours: f64,
    /// Integrator step, in hours.
    pub step_hours: f64,
    /// Record one sample every this many steps.
    pub sample_every: usize,
    /// Which disturbances are on, one-based as `IDV(n)` is: index 0 is
    /// `IDV(1)`.
    pub disturbances: [bool; DISTURBANCES],
    /// Whether to run the plant closed loop under the published control
    /// scheme, or open loop with the valves held.
    ///
    /// Open loop is not a useful operating mode, it is a *diagnostic* one: the
    /// plant trips on reactor pressure after about three hours. It exists here
    /// because the difference between the two is the clearest single statement
    /// of what the control layer does.
    pub controlled: bool,
    /// Which Class C quirks are fixed rather than reproduced.
    ///
    /// All off by default, so a default scenario reproduces the original. See
    /// `book/src/deltas.md`.
    pub quirks: QuirkFixes,
    /// Whether the driver forces `IDV(12)` on at eight hours, as
    /// `temain_mod.f:366-368` does, regardless of what was asked for.
    ///
    /// `true` by default, because every published dataset longer than eight
    /// hours carries it. Delta D-011.
    pub driver_forces_idv12: bool,
    /// Events on a schedule: faults that arrive, clear, or change magnitude
    /// part way through, and setpoint moves.
    ///
    /// Empty by default, which is the only shape the original admits. A
    /// non-empty schedule is what lets a study ask what happens when a fault
    /// arrives at hour twelve and clears at hour twenty, which no published
    /// Tennessee Eastman dataset contains.
    pub schedule: Schedule,
    /// Capabilities beyond the original. All off by default.
    pub extensions: Extensions,
    /// How to advance the state.
    ///
    /// [`Integrator::Euler`] by default, which is what the original does and
    /// the only choice under which this port reproduces it bit for bit.
    /// Anything else is a better integration of the same equations and a
    /// different set of numbers; see [`Integrator::is_faithful`].
    pub integrator: Integrator,
}

impl Default for Scenario {
    fn default() -> Self {
        Self::baseline()
    }
}

impl Scenario {
    /// The fault-free plant for 48 hours, closed loop, from the seed
    /// `teprob.f` ships with.
    ///
    /// 48 hours is `NPTS = 172800` at a one-second step, which is the run
    /// `temain_mod.f` was written to do.
    #[must_use]
    pub const fn baseline() -> Self {
        Self {
            seed: DEFAULT_SEED,
            hours: 48.0,
            step_hours: DEFAULT_STEP_HOURS,
            sample_every: DEFAULT_SAMPLE_EVERY,
            disturbances: [false; DISTURBANCES],
            controlled: true,
            quirks: QuirkFixes::new(),
            driver_forces_idv12: true,
            schedule: Schedule::new(),
            extensions: Extensions::none(),
            integrator: Integrator::Euler,
        }
    }

    /// Add an event to the schedule.
    #[must_use]
    pub fn with_event(mut self, event: tepsim_scenario::Event) -> Self {
        self.schedule.add(event);
        self
    }

    /// Allow a disturbance to be partly on.
    ///
    /// See [`Extensions::continuous_disturbances`]. A run with this on is not
    /// comparable to any published dataset.
    #[must_use]
    pub const fn with_continuous_disturbances(mut self) -> Self {
        self.extensions.continuous_disturbances = true;
        self
    }

    /// Check the scenario is runnable.
    ///
    /// # Errors
    ///
    /// The first problem found. A fractional magnitude without the extension
    /// is rejected rather than rounded: silently turning a request for half a
    /// fault into a whole one would produce a run that does not match its own
    /// description, which is exactly what the content hash exists to prevent.
    pub fn validate(&self) -> Result<(), Invalid> {
        self.schedule
            .validate(self.extensions.continuous_disturbances)
    }

    /// A content hash over everything that affects the run.
    ///
    /// Two scenarios describing the same experiment hash the same, and two
    /// that differ in any respect do not. A dataset carrying this says what
    /// produced it, rather than relying on a filename and a memory.
    #[must_use]
    pub fn digest(&self) -> tepsim_scenario::Digest {
        let mut digest = tepsim_scenario::Digest::new();
        // Versioned, so a later change to what a scenario contains cannot make
        // an old digest silently mean something else.
        digest.push_str("tepsim.scenario.v1");
        digest.push_f64(self.seed);
        digest.push_f64(self.hours);
        digest.push_f64(self.step_hours);
        digest.push_usize(self.sample_every);
        for on in self.disturbances {
            digest.push_bool(on);
        }
        digest.push_bool(self.controlled);
        digest.push_bool(self.quirks.trip_ends_the_run);
        digest.push_bool(self.driver_forces_idv12);
        digest.push_bool(self.extensions.continuous_disturbances);
        digest.push_str(self.integrator.name());
        self.schedule.hash_into(&mut digest);
        digest
    }

    /// The digest as sixteen hex characters, for a filename.
    #[must_use]
    pub fn digest_hex(&self) -> [u8; 16] {
        tepsim_scenario::short_hex(self.digest())
    }

    /// Choose the integrator.
    ///
    /// Leaving this alone keeps the run faithful to the original. See
    /// [`Integrator`].
    #[must_use]
    pub const fn with_integrator(mut self, integrator: Integrator) -> Self {
        self.integrator = integrator;
        self
    }

    /// The baseline with one disturbance on.
    ///
    /// # Panics
    ///
    /// If `n` is not in `1..=20`.
    #[must_use]
    pub const fn fault(n: usize) -> Self {
        assert!(
            n >= 1 && n <= DISTURBANCES,
            "IDV index out of range: this model has twenty disturbances"
        );
        let mut scenario = Self::baseline();
        scenario.disturbances[n - 1] = true;
        scenario
    }

    /// Set the seed.
    #[must_use]
    pub const fn with_seed(mut self, seed: f64) -> Self {
        self.seed = seed;
        self
    }

    /// Set the duration in hours.
    #[must_use]
    pub const fn with_hours(mut self, hours: f64) -> Self {
        self.hours = hours;
        self
    }

    /// Set the output cadence in steps.
    #[must_use]
    pub const fn sampling_every(mut self, steps: usize) -> Self {
        self.sample_every = steps;
        self
    }

    /// Run open loop, with the valves held at their initial positions.
    #[must_use]
    pub const fn open_loop(mut self) -> Self {
        self.controlled = false;
        self
    }

    /// Turn a disturbance on.
    ///
    /// # Panics
    ///
    /// If `n` is not in `1..=20`.
    #[must_use]
    pub const fn with_fault(mut self, n: usize) -> Self {
        assert!(n >= 1 && n <= DISTURBANCES, "IDV index out of range");
        self.disturbances[n - 1] = true;
        self
    }

    /// How many integrator steps this scenario is.
    #[must_use]
    pub fn steps(&self) -> usize {
        // Rounded rather than truncated. `f64::round` is a `std` method and
        // this crate is `no_std`, so the rounding is written out: the quotient
        // is positive by construction, so adding a half and truncating is
        // round-half-up.
        let quotient = self.hours / self.step_hours;
        (quotient + 0.5) as usize
    }

    /// How many samples it will record.
    #[must_use]
    pub fn samples(&self) -> usize {
        self.steps() / self.sample_every
    }

    /// The disturbance vector the plant reads, as the model's `f64` flags.
    #[must_use]
    pub fn disturbance_vector(&self) -> [f64; DISTURBANCES] {
        core::array::from_fn(|i| f64::from(u8::from(self.disturbances[i])))
    }

    /// Which faults are on, one-based, for a label or a report.
    pub fn active_faults(&self) -> impl Iterator<Item = usize> + '_ {
        self.disturbances
            .iter()
            .enumerate()
            .filter(|(_, on)| **on)
            .map(|(index, _)| index + 1)
    }
}
