//! What to simulate: how long, from which seed, with which disturbances.

use alloc::string::String;

use tepsim_core::{Extensions, FAULTS, QuirkFixes};
use tepsim_scenario::{Invalid, Schedule};

use crate::integrator::Integrator;
use crate::text::TextError;

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

/// What a [`Scenario`] contains, as a version string.
///
/// Absorbed by [`Scenario::digest`] and written as the leading tag of
/// [`Scenario::to_text`], deliberately the same string in both places. The
/// digest and the text describe the same set of fields, so a change to that
/// set has to move both at once, and sharing the constant is what makes that
/// automatic rather than remembered.
pub const SCENARIO_VERSION: &str = "tepsim.scenario.v1";

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
    /// [`QuirkFixes::default`] applies every signed-off fix;
    /// [`QuirkFixes::faithful`] reproduces the original instead, and is what
    /// [`Scenario::faithful`] uses. See `book/src/deltas.md`.
    pub quirks: QuirkFixes,
    /// Whether the driver forces `IDV(12)` on at eight hours, as
    /// `temain_mod.f:366-368` does, regardless of what was asked for.
    ///
    /// **`false` by default.** Delta D-011.
    ///
    /// # Why the default is off when the driver's source says on
    ///
    /// `temain_mod.f:367` is literally `IDV(12)=1`, inside
    /// `IF (I.GE.SSPTS)`. Read alone it says every run past eight hours carries
    /// the condenser cooling-water disturbance. The prose at
    /// `temain_mod.f:101-102` says something else: *"Go to line 367, implement
    /// any of the 21 programmed disturbances"*, which reads as an instruction
    /// to replace that line rather than to add beside it.
    ///
    /// The source does not decide between the two readings and the published
    /// bytes do. Every `dNN_te.dat` except `d12_te.dat` sits at the nominal
    /// operating point straight across row 160, which is hour eight, and it
    /// could not if `IDV(12)` had switched on there. Tier 7's
    /// `the_published_files_were_not_generated_with_the_forced_idv12` is that
    /// measurement.
    ///
    /// So the shipped line is an example, the files were made with it replaced,
    /// and the default follows the files. Signed off 2026-08-28.
    /// [`Scenario::faithful`] and `tep run --force-idv12` reproduce the driver
    /// as shipped.
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
            driver_forces_idv12: false,
            schedule: Schedule::new(),
            extensions: Extensions::none(),
            integrator: Integrator::Euler,
        }
    }

    /// The baseline with every Class C quirk reproduced rather than fixed.
    ///
    /// This is what a differential against `teprob.f` runs, and what
    /// regenerating anything published-shaped runs. Two flags separate it from
    /// [`Scenario::baseline`]: the plant freezes on a trip instead of the run
    /// ending ([`QuirkFixes::trip_ends_the_run`], delta D-007), and the driver
    /// forces `IDV(12)` on at eight hours ([`Scenario::driver_forces_idv12`],
    /// delta D-011). Both defaults were decided against on 2026-08-28, on
    /// evidence recorded at each field.
    ///
    /// It is a separate constructor rather than the default because the two
    /// configurations answer different questions. A user of the simulator wants
    /// a trip to stop the run; a test asking whether this port *is*
    /// `teprob.f` cannot afford a single deliberate difference.
    #[must_use]
    pub const fn faithful() -> Self {
        let mut scenario = Self::baseline();
        scenario.quirks = QuirkFixes::faithful();
        scenario.driver_forces_idv12 = true;
        scenario
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
        digest.push_str(SCENARIO_VERSION);
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

    /// This scenario as its canonical text.
    ///
    /// Every field, written out, tagged with [`SCENARIO_VERSION`]. Round-trips
    /// through [`Scenario::from_text`] bit for bit, so the digest survives it,
    /// which is what lets a run travel as a line of text rather than as a file.
    /// The format is described in [`crate::text`].
    ///
    /// ```
    /// use tepsim::Scenario;
    ///
    /// let scenario = Scenario::fault(4).with_hours(8.0);
    /// let text = scenario.to_text();
    /// assert!(text.starts_with("tepsim.scenario.v1;"));
    /// assert_eq!(Scenario::from_text(&text), Ok(scenario));
    /// ```
    #[must_use]
    pub fn to_text(&self) -> String {
        crate::text::to_text(self)
    }

    /// Parse a canonical scenario text.
    ///
    /// Strict: every field must be present, unknown fields are named, and a
    /// value outside its range is rejected rather than clamped. See
    /// [`TextError`] for the reasons and [`crate::text`] for the format.
    ///
    /// # Errors
    ///
    /// The first problem found, as a [`TextError`] whose `Display` says what
    /// was wrong.
    pub fn from_text(text: &str) -> Result<Self, TextError> {
        crate::text::from_text(text)
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
