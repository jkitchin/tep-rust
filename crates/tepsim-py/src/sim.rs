//! `tepsim.Simulation`: the thing that runs a scenario.

use pyo3::prelude::*;

use crate::run::Run;
use crate::scenario::Scenario;

/// A simulation ready to run.
///
/// Holds a plant, a controller stack and an integrator state, all of it owned.
/// The original keeps its whole working set in six Fortran `COMMON` blocks,
/// which allows exactly one simulation per process and no reentrancy. This
/// allows as many as there are threads, which is the point of the next
/// paragraph.
///
/// ```python
/// run = tep.Simulation(tep.Scenario.baseline(seed=42, hours=48)).run()
/// ```
///
/// # Threads
///
/// `run()` releases the GIL for its whole duration, so an ensemble is a
/// `ThreadPoolExecutor` and not a `ProcessPoolExecutor`. Nothing has to be
/// pickled, and the plant does no I/O and touches no Python object while the
/// GIL is down.
///
/// ```python
/// from concurrent.futures import ThreadPoolExecutor
///
/// sims = [tep.Simulation(tep.Scenario.fault(n, hours=8)) for n in range(1, 21)]
/// with ThreadPoolExecutor() as pool:
///     runs = list(pool.map(tep.Simulation.run, sims))
/// ```
///
/// # `run()` does not consume the simulation
///
/// It runs a copy, so calling it twice gives two equal runs rather than one run
/// and one empty one. A run is a pure function of its scenario: no clock, no
/// global state, no randomness outside the seeded generator. That is what makes
/// a recorded dataset reproducible from its description rather than from a
/// file, and it is asserted in the test suite.
#[pyclass(frozen, module = "tepsim")]
#[derive(Debug)]
pub struct Simulation {
    inner: tepsim::Simulation,
}

#[pymethods]
impl Simulation {
    /// A simulation ready to run the given scenario.
    ///
    /// Args:
    ///     scenario: What to simulate. Defaults to `Scenario.baseline()`.
    #[new]
    #[pyo3(signature = (scenario = None))]
    fn new(scenario: Option<Scenario>) -> Self {
        let scenario = scenario.map_or_else(tepsim::Scenario::baseline, |s| s.inner);
        Self {
            inner: tepsim::Simulation::new(scenario),
        }
    }

    /// What this simulation is running.
    #[getter]
    fn scenario(&self) -> Scenario {
        Scenario::wrap(*self.inner.scenario())
    }

    /// Run the whole scenario and collect every sample.
    ///
    /// Releases the GIL for the integration, which is where all the time goes:
    /// a 48-hour closed-loop run is 172,800 integrator steps.
    ///
    /// Never raises for a plant that misbehaves. A trip or a failed temperature
    /// solve is reported through `Run.outcome`, because a run that ended early
    /// is data: throwing it away would hide the difference between a port that
    /// trips where the original does and one that does not.
    fn run(&self, py: Python<'_>) -> PyResult<Run> {
        // Cloned so the simulation stays usable, and so nothing Python owns is
        // borrowed across the GIL release.
        let simulation = self.inner.clone();
        let finished = py.detach(move || simulation.run());
        Run::build(py, finished)
    }

    fn __repr__(&self) -> String {
        let scenario = self.inner.scenario();
        format!(
            "<tepsim.Simulation {:.1} h, {} steps, {}>",
            scenario.hours,
            scenario.steps(),
            if scenario.controlled {
                "closed loop"
            } else {
                "open loop"
            }
        )
    }
}
