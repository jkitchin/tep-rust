//! `tepsim.faults()`: the twenty disturbances, as what they actually are.

use pyo3::prelude::*;
use pyo3::types::PyTuple;

use tepsim::tepsim_core::{FAULTS, Shape};

/// One of the twenty disturbances.
///
/// `published` is the description at `teprob.f:172-191` verbatim, five of which
/// say only "Unknown". `effect` is what the source actually does, which for
/// those five is perfectly explicit: only the physical interpretation was
/// withheld. `line` names the `teprob.f` line the fault acts on, so a sceptical
/// reader can check.
///
/// ```python
/// for f in tep.faults():
///     print(f.index, f.shape, f.published)
/// ```
#[pyclass(frozen, module = "tepsim", skip_from_py_object)]
#[derive(Clone, Copy, Debug)]
pub struct Fault {
    inner: tepsim::tepsim_core::Fault,
}

#[pymethods]
impl Fault {
    /// The `IDV` index, one-based.
    #[getter]
    const fn index(&self) -> usize {
        self.inner.index
    }

    /// The description from `teprob.f:172-191`, verbatim.
    #[getter]
    const fn published(&self) -> &'static str {
        self.inner.published
    }

    /// What it does, where the header says "Unknown" or where the prose is less
    /// specific than the source.
    #[getter]
    const fn effect(&self) -> &'static str {
        self.inner.effect
    }

    /// How it enters the model: `'step'`, `'random'` or `'sticking'`.
    ///
    /// A step fault changes a feed condition the moment it is switched on and
    /// holds it. A random fault enables a walk channel, which then wanders on
    /// its own schedule. A sticking fault touches no equation at all: it widens
    /// the dead band a valve command must cross before the valve follows, which
    /// is a disturbance to the controller's authority over the plant rather
    /// than to the plant. In an open-loop run, where the command never moves, a
    /// sticking fault does nothing whatever.
    #[getter]
    const fn shape(&self) -> &'static str {
        match self.inner.shape {
            Shape::Step => "step",
            Shape::Random { .. } => "random",
            Shape::Sticking { .. } => "sticking",
        }
    }

    /// The walk channels this fault enables, one-based; empty unless `shape` is
    /// `'random'`. `IDV(8)` and `IDV(13)` drive two each.
    #[getter]
    fn channels<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        match self.inner.shape {
            Shape::Random { channels, .. } => PyTuple::new(py, channels),
            _ => Ok(PyTuple::empty(py)),
        }
    }

    /// Whether this fault's channels are spike trains rather than walks.
    ///
    /// The three spiking faults are intermittent rather than sustained, which
    /// is why the literature reports them as the hardest to detect.
    #[getter]
    const fn spiking(&self) -> bool {
        matches!(self.inner.shape, Shape::Random { spiking: true, .. })
    }

    /// The valves this fault sticks, one-based; empty unless `shape` is
    /// `'sticking'`.
    #[getter]
    fn valves<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        match self.inner.shape {
            Shape::Sticking { valves } => PyTuple::new(py, valves),
            _ => Ok(PyTuple::empty(py)),
        }
    }

    /// The `teprob.f` line this fault acts on.
    #[getter]
    const fn line(&self) -> &'static str {
        self.inner.line
    }

    /// Whether this fault reaches the plant model at all.
    ///
    /// False for the three sticking faults. A scenario engine that treated
    /// those as plant faults would report an injected disturbance with no
    /// effect and look broken.
    #[getter]
    const fn affects_the_plant(&self) -> bool {
        self.inner.affects_the_plant()
    }

    fn __repr__(&self) -> String {
        format!(
            "<tepsim.Fault IDV({}) {} '{}'>",
            self.inner.index,
            self.shape(),
            self.inner.published
        )
    }
}

/// The twenty disturbances, in `IDV` order, with what each one does.
///
/// Twenty, not the twenty-one of the later literature: `teprob.f:340` loops
/// `DO 500 I=1,20`.
#[pyfunction]
pub fn faults(py: Python<'_>) -> PyResult<Bound<'_, PyTuple>> {
    PyTuple::new(py, FAULTS.map(|inner| Fault { inner }))
}
