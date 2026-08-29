//! `tepsim.Scenario`: what to simulate.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use tepsim::DISTURBANCES;

/// The keyword arguments every constructor here shares.
///
/// A struct rather than seven more parameters on the shared builder: the three
/// Python-facing constructors already carry them one by one, because that is
/// how keyword arguments and type stubs work, and repeating that shape a fourth
/// time inside the crate buys nothing.
struct Options {
    seed: f64,
    hours: f64,
    step_hours: f64,
    sample_every: usize,
    controlled: bool,
    driver_forces_idv12: bool,
    trip_ends_the_run: bool,
}

/// A complete description of a run.
///
/// Immutable and cheap to copy, so a caller can build one, derive variants from
/// it, and keep all of them. Everything a run's output depends on is in here,
/// which is what makes a recorded dataset reproducible from its scenario alone
/// rather than from a file.
///
/// `repr()` round-trips: `eval(repr(s))` reconstructs an equal scenario.
///
/// ```python
/// import tepsim as tep
///
/// tep.Scenario.baseline(seed=42, hours=48)
/// tep.Scenario.fault(4, hours=8)
/// tep.Scenario(hours=2, faults=[1, 8], sample_every=60)
/// ```
#[pyclass(frozen, module = "tepsim", eq, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Scenario {
    pub(crate) inner: tepsim::Scenario,
}

impl Scenario {
    /// Validate the options, resolve the fault numbers, and wrap the result.
    fn build(options: Options, faults: &[usize]) -> PyResult<Self> {
        // Rejected rather than clamped: an out-of-range fault number in a
        // scenario is a typo, and quietly running the fault-free plant instead
        // is the worst possible answer, because the resulting dataset looks
        // fine. The finiteness checks are there for the same reason: a NaN
        // `hours` would give a zero-step run that reports success.
        if !options.step_hours.is_finite() || options.step_hours <= 0.0 {
            return Err(PyValueError::new_err(format!(
                "step_hours must be a positive finite number, got {}",
                options.step_hours
            )));
        }
        if !options.hours.is_finite() || options.hours < 0.0 {
            return Err(PyValueError::new_err(format!(
                "hours must be a finite number and not negative, got {}",
                options.hours
            )));
        }
        if options.sample_every == 0 {
            return Err(PyValueError::new_err(
                "sample_every must be at least 1: a run has to record something",
            ));
        }

        let mut active = [false; DISTURBANCES];
        for &n in faults {
            active[check_fault_number(n)? - 1] = true;
        }

        let mut inner = tepsim::Scenario::baseline();
        inner.seed = options.seed;
        inner.hours = options.hours;
        inner.step_hours = options.step_hours;
        inner.sample_every = options.sample_every;
        inner.disturbances = active;
        inner.controlled = options.controlled;
        inner.driver_forces_idv12 = options.driver_forces_idv12;
        inner.quirks.trip_ends_the_run = options.trip_ends_the_run;
        Ok(Self { inner })
    }

    /// Wrap a facade scenario that is already known to be valid.
    pub(crate) const fn wrap(inner: tepsim::Scenario) -> Self {
        Self { inner }
    }
}

/// A one-based `IDV` number, or a `ValueError` naming the range.
fn check_fault_number(n: usize) -> PyResult<usize> {
    if (1..=DISTURBANCES).contains(&n) {
        Ok(n)
    } else {
        Err(PyValueError::new_err(format!(
            "fault number {n} is out of range: this model has {DISTURBANCES} \
             disturbances, numbered 1 to {DISTURBANCES} as IDV(n) is"
        )))
    }
}

#[pymethods]
impl Scenario {
    /// A scenario with every field given explicitly.
    ///
    /// The defaults are the baseline: the fault-free plant for 48 hours, closed
    /// loop, from the generator word compiled into `teprob.f:1187`, sampled
    /// every 180 seconds as the published `d00`-`d21` files are.
    ///
    /// Args:
    ///     seed: The generator word. Any positive integer; odd values give the
    ///         generator its full period.
    ///     hours: How long to run.
    ///     step_hours: Integrator step. The original's is one second.
    ///     sample_every: Record one sample every this many steps.
    ///     faults: Which disturbances are on, one-based as `IDV(n)` is.
    ///     controlled: Closed loop under the published control scheme, or open
    ///         loop with the valves held. Open loop is a diagnostic, not an
    ///         operating mode: the plant trips on reactor pressure after about
    ///         three hours.
    ///     driver_forces_idv12: Whether to reproduce `temain_mod.f:366-368`,
    ///         which switches `IDV(12)` on at eight hours whatever the scenario
    ///         asked for. False by default: the published `dNN_te` files sit at
    ///         the nominal operating point straight across hour eight, so they
    ///         were not made with it. Delta D-011, signed off 2026-08-28.
    ///     trip_ends_the_run: Whether a trip stops the run. True by default.
    ///         `teprob.f:807-811` instead freezes the plant and keeps
    ///         reporting, which is 75.6% of `d06.dat`: 363 of its 480 rows are
    ///         a stopped plant that reads as a very steady one. Pass
    ///         `trip_ends_the_run=False, driver_forces_idv12=True` to get the
    ///         configuration every oracle comparison runs, which Rust names
    ///         `Scenario::faithful`. Delta D-007, signed off 2026-08-28.
    ///
    /// Raises:
    ///     ValueError: If a fault number is outside 1 to 20, or if `hours`,
    ///         `step_hours` or `sample_every` is not a usable size.
    #[new]
    #[pyo3(signature = (
        *,
        seed = tepsim::scenario::DEFAULT_SEED,
        hours = 48.0,
        step_hours = tepsim::scenario::DEFAULT_STEP_HOURS,
        sample_every = tepsim::scenario::DEFAULT_SAMPLE_EVERY,
        faults = Vec::new(),
        controlled = true,
        driver_forces_idv12 = false,
        trip_ends_the_run = true,
    ))]
    #[expect(
        clippy::too_many_arguments,
        reason = "one keyword per field of the scenario. Collapsing them into a \
                  dict would lose the defaults, the type stubs and help()."
    )]
    fn new(
        seed: f64,
        hours: f64,
        step_hours: f64,
        sample_every: usize,
        faults: Vec<usize>,
        controlled: bool,
        driver_forces_idv12: bool,
        trip_ends_the_run: bool,
    ) -> PyResult<Self> {
        Self::build(
            Options {
                seed,
                hours,
                step_hours,
                sample_every,
                controlled,
                driver_forces_idv12,
                trip_ends_the_run,
            },
            &faults,
        )
    }

    /// The fault-free plant, closed loop.
    ///
    /// Identical to the constructor with no faults; it exists because
    /// `Scenario.baseline(seed=42, hours=48)` says what it is at the call site
    /// and `Scenario(seed=42, hours=48)` does not.
    #[staticmethod]
    #[pyo3(signature = (
        *,
        seed = tepsim::scenario::DEFAULT_SEED,
        hours = 48.0,
        step_hours = tepsim::scenario::DEFAULT_STEP_HOURS,
        sample_every = tepsim::scenario::DEFAULT_SAMPLE_EVERY,
        controlled = true,
        driver_forces_idv12 = false,
        trip_ends_the_run = true,
    ))]
    fn baseline(
        seed: f64,
        hours: f64,
        step_hours: f64,
        sample_every: usize,
        controlled: bool,
        driver_forces_idv12: bool,
        trip_ends_the_run: bool,
    ) -> PyResult<Self> {
        Self::build(
            Options {
                seed,
                hours,
                step_hours,
                sample_every,
                controlled,
                driver_forces_idv12,
                trip_ends_the_run,
            },
            &[],
        )
    }

    /// The baseline with one disturbance on, one-based as `IDV(n)` is.
    ///
    /// Raises:
    ///     ValueError: If `n` is outside 1 to 20.
    #[staticmethod]
    #[pyo3(signature = (
        n,
        *,
        seed = tepsim::scenario::DEFAULT_SEED,
        hours = 48.0,
        step_hours = tepsim::scenario::DEFAULT_STEP_HOURS,
        sample_every = tepsim::scenario::DEFAULT_SAMPLE_EVERY,
        controlled = true,
        driver_forces_idv12 = false,
        trip_ends_the_run = true,
    ))]
    #[expect(
        clippy::too_many_arguments,
        reason = "the fault number plus the same seven keywords the other two \
                  constructors take; they have to stay in step."
    )]
    fn fault(
        n: usize,
        seed: f64,
        hours: f64,
        step_hours: f64,
        sample_every: usize,
        controlled: bool,
        driver_forces_idv12: bool,
        trip_ends_the_run: bool,
    ) -> PyResult<Self> {
        Self::build(
            Options {
                seed,
                hours,
                step_hours,
                sample_every,
                controlled,
                driver_forces_idv12,
                trip_ends_the_run,
            },
            &[n],
        )
    }

    /// The generator word.
    #[getter]
    const fn seed(&self) -> f64 {
        self.inner.seed
    }

    /// How long the run is, in hours.
    #[getter]
    const fn hours(&self) -> f64 {
        self.inner.hours
    }

    /// Integrator step, in hours.
    #[getter]
    const fn step_hours(&self) -> f64 {
        self.inner.step_hours
    }

    /// One sample every this many integrator steps.
    #[getter]
    const fn sample_every(&self) -> usize {
        self.inner.sample_every
    }

    /// Whether the run is closed loop.
    #[getter]
    const fn controlled(&self) -> bool {
        self.inner.controlled
    }

    /// Whether the driver forces `IDV(12)` on at eight hours. Delta D-011.
    #[getter]
    const fn driver_forces_idv12(&self) -> bool {
        self.inner.driver_forces_idv12
    }

    /// Whether a trip stops the run. Delta D-007.
    #[getter]
    const fn trip_ends_the_run(&self) -> bool {
        self.inner.quirks.trip_ends_the_run
    }

    /// The active disturbances, one-based and ascending.
    #[getter]
    fn faults<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        // Collected because `PyTuple::new` sizes the tuple up front and
        // `active_faults` is a filter, which cannot say how long it is.
        let active: Vec<usize> = self.inner.active_faults().collect();
        PyTuple::new(py, active)
    }

    /// How many integrator steps this scenario is.
    #[getter]
    fn steps(&self) -> usize {
        self.inner.steps()
    }

    /// How many samples it will record, if it does not stop early.
    #[getter]
    fn samples(&self) -> usize {
        self.inner.samples()
    }

    /// A content hash over everything that affects the run, as 16 hex
    /// characters.
    ///
    /// Two scenarios describing the same experiment share this and two that
    /// differ do not. Written beside a dataset it says what produced the file,
    /// rather than leaving that to a filename and a memory. It survives
    /// `to_text`/`from_text`, which is what makes a serialised scenario worth
    /// trusting.
    #[getter]
    fn digest(&self) -> String {
        let hex = self.inner.digest_hex();
        // ASCII by construction: `digest_hex` writes from a 16-character table.
        String::from_utf8_lossy(&hex).into_owned()
    }

    /// This scenario as one line of canonical text.
    ///
    /// Every field written out and tagged with the format version, so the text
    /// says what it runs rather than what it leaves to a default. It is the
    /// same string the Rust, wasm and browser sides read and write, which is
    /// what lets a run move between them as a line in a file or a fragment in
    /// a URL.
    ///
    /// ```python
    /// import tepsim as tep
    ///
    /// s = tep.Scenario.fault(4, hours=8)
    /// assert tep.Scenario.from_text(s.to_text()) == s
    /// ```
    #[expect(
        clippy::wrong_self_convention,
        reason = "a `#[pymethods]` instance method takes `&self`; the name is \
                  the Python-facing one and pairs with `from_text`"
    )]
    fn to_text(&self) -> String {
        self.inner.to_text()
    }

    /// Parse a scenario back from `to_text`.
    ///
    /// Strict on purpose. A missing field, an unknown field, a malformed
    /// number, a value out of range or a version this build does not know is a
    /// `ValueError` naming what was wrong. Nothing is defaulted quietly,
    /// because a scenario that silently differs from its description is the
    /// one failure the digest exists to prevent.
    ///
    /// Raises:
    ///     ValueError: With a message saying what the text got wrong.
    #[staticmethod]
    fn from_text(text: &str) -> PyResult<Self> {
        tepsim::Scenario::from_text(text)
            .map(Self::wrap)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// This scenario with a different generator word.
    fn with_seed(&self, seed: f64) -> Self {
        Self::wrap(self.inner.with_seed(seed))
    }

    /// This scenario run for a different length of time.
    ///
    /// Raises:
    ///     ValueError: If `hours` is negative or not finite.
    fn with_hours(&self, hours: f64) -> PyResult<Self> {
        if !hours.is_finite() || hours < 0.0 {
            return Err(PyValueError::new_err(format!(
                "hours must be a finite number and not negative, got {hours}"
            )));
        }
        Ok(Self::wrap(self.inner.with_hours(hours)))
    }

    /// This scenario with one more disturbance switched on.
    ///
    /// Raises:
    ///     ValueError: If `n` is outside 1 to 20.
    fn with_fault(&self, n: usize) -> PyResult<Self> {
        Ok(Self::wrap(self.inner.with_fault(check_fault_number(n)?)))
    }

    /// This scenario at a different output cadence, in integrator steps.
    ///
    /// Raises:
    ///     ValueError: If `steps` is zero.
    fn sampling_every(&self, steps: usize) -> PyResult<Self> {
        if steps == 0 {
            return Err(PyValueError::new_err(
                "sample_every must be at least 1: a run has to record something",
            ));
        }
        Ok(Self::wrap(self.inner.sampling_every(steps)))
    }

    /// This scenario run open loop, with the valves held where they started.
    ///
    /// A diagnostic, not an operating mode. The plant trips on reactor pressure
    /// after about three hours, and the difference between this and the closed
    /// loop is the clearest single statement of what the control layer does.
    fn open_loop(&self) -> Self {
        Self::wrap(self.inner.open_loop())
    }

    /// `repr` that round-trips, in both of the two shapes it needs.
    ///
    /// # Why this is not always a constructor call
    ///
    /// The docstring on this class promises `eval(repr(s)) == s`, and the
    /// obvious implementation quietly broke that promise. `Scenario(...)`
    /// cannot express three of the fields a scenario carries: the schedule, the
    /// `continuous_disturbances` extension and the integrator. A `repr` that
    /// printed only the constructor's arguments produced valid Python that
    /// evaluated to a *different* scenario, with a different digest and a
    /// different run, and nothing said so.
    ///
    /// So a scenario the constructor can express is printed as a constructor
    /// call, because that is what is readable at a prompt, and anything else is
    /// printed as `Scenario.from_text(...)`, which round-trips by construction:
    /// the text carries every field and `tests/test_bindings.py` asserts the
    /// pair on scenarios of both shapes.
    fn __repr__(&self) -> String {
        let baseline = tepsim::Scenario::baseline();
        let expressible = self.inner.schedule.is_empty()
            && self.inner.extensions == baseline.extensions
            && self.inner.integrator == baseline.integrator;

        if !expressible {
            return format!("Scenario.from_text({:?})", self.inner.to_text());
        }
        format!(
            "Scenario(seed={:?}, hours={:?}, step_hours={:?}, sample_every={}, \
             faults=({}), controlled={}, driver_forces_idv12={}, trip_ends_the_run={})",
            self.inner.seed,
            self.inner.hours,
            self.inner.step_hours,
            self.inner.sample_every,
            tuple_body(self.inner.active_faults()),
            python_bool(self.inner.controlled),
            python_bool(self.inner.driver_forces_idv12),
            python_bool(self.inner.quirks.trip_ends_the_run),
        )
    }
}

/// `True`/`False` rather than Rust's lowercase, so `repr` stays evaluable.
const fn python_bool(value: bool) -> &'static str {
    if value { "True" } else { "False" }
}

/// The inside of a Python tuple literal, trailing comma included when there is
/// exactly one element.
///
/// Without that comma `faults=(1)` is the integer 1, and `eval(repr(scenario))`
/// would quietly produce a scenario with no faults rather than raising.
fn tuple_body(faults: impl Iterator<Item = usize>) -> String {
    let listed: Vec<String> = faults.map(|n| n.to_string()).collect();
    match listed.as_slice() {
        [] => String::new(),
        [only] => format!("{only},"),
        many => many.join(", "),
    }
}
