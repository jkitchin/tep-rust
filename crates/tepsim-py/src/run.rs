//! `tepsim.Run`: what a run produced, as NumPy arrays.
//!
//! # One copy, then views
//!
//! [`tepsim::Run`] holds a `Vec<Sample>`, and a `Sample` interleaves the 41
//! measurements and 12 manipulated variables with a step number, a time, and
//! the ground truth for all twenty disturbances. NumPy needs one contiguous
//! block per array, so the interleaving is undone exactly once. That happens in
//! `Run::build`, which fills five owned `ndarray` buffers and hands each to
//! NumPy through `into_pyarray`, a move rather than a copy: the block this
//! crate filled *is* the block Python sees.
//!
//! Everything after that is a view. `to_numpy` hands back the same object every
//! time. `measurement`, `manipulated` and `column` index a cached transpose, so
//! a column is a strided view of the matrix rather than a fresh allocation, and
//! `columns` is fifty-three of those in a dict. Nothing here copies.
//!
//! # Why the arrays are read-only
//!
//! Because the same objects are handed to every caller. Two callers of
//! `to_numpy` get one array, and a write through either would rewrite the
//! other's data and every column view besides. A run is a record of what the
//! plant did; `.copy()` is the way to get something writable.
//!
//! `make_nonwriteable` is `numpy`'s safe wrapper for clearing the `WRITEABLE`
//! flag, which is why this module contains no `unsafe` block. NumPy lets Python
//! set that flag back only on an array that owns its data, and an array built
//! by `into_pyarray` does not: its data belongs to the Rust `Vec` in its base
//! object. So read-only here means read-only.

use numpy::ndarray::{Array1, Array2, Dimension};
use numpy::{Element, IntoPyArray, PyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use tepsim::{CHANNELS, DISTURBANCES, MANIPULATED, MEASUREMENTS, Outcome, Sample};

use crate::scenario::Scenario;

/// Move an owned array into NumPy and mark it read-only.
///
/// `into_pyarray` transfers the backing `Vec` rather than copying it; the
/// `Vec` becomes the new array's base object and is freed when NumPy drops it.
/// The guard returned by `make_nonwriteable` is discarded on purpose: clearing
/// the flag is a permanent change to the NumPy object, not something the guard
/// holds open.
fn frozen<'py, T, D>(
    py: Python<'py>,
    data: numpy::ndarray::Array<T, D>,
) -> Bound<'py, PyArray<T, D>>
where
    T: Element,
    D: Dimension,
{
    let array = data.into_pyarray(py);
    let _read_only = array.readwrite().make_nonwriteable();
    array
}

/// A finished run: its scenario, its samples as arrays, and how it ended.
///
/// Construct one by calling `Simulation.run`. There is no useful way to build
/// one from Python: its contents are the output of the plant.
///
/// ```python
/// run = tep.Simulation(tep.Scenario.baseline(hours=2)).run()
/// run.to_numpy().shape        # (40, 53)
/// run.measurement(7)          # XMEAS(7), reactor pressure
/// run.outcome                 # 'completed'
/// ```
#[pyclass(frozen, module = "tepsim")]
#[derive(Debug)]
pub struct Run {
    scenario: tepsim::Scenario,
    outcome: Outcome,
    samples: usize,
    /// `(n_samples, 53)`, measurements then manipulated variables.
    matrix: Py<PyArray2<f64>>,
    /// `(53, n_samples)`, a NumPy view of `matrix`, held so that a single
    /// channel is one index away rather than one allocation away.
    transposed: Py<PyArray2<f64>>,
    /// Simulated time at each sample, in hours.
    hours: Py<PyArray1<f64>>,
    /// The integrator step each sample was taken at, one-based as
    /// `temain_mod.f`'s `I` is.
    steps: Py<PyArray1<i64>>,
    /// `(n_samples, 20)`, which disturbances were active.
    active: Py<PyArray2<bool>>,
    /// `(n_samples, 20)`, hours since each disturbance came on, NaN where it
    /// never did.
    since_onset: Py<PyArray2<f64>>,
}

impl Run {
    /// Unpack a finished run into NumPy arrays. The one copy happens here.
    pub(crate) fn build(py: Python<'_>, run: tepsim::Run) -> PyResult<Self> {
        let samples: &[Sample] = &run.samples;
        let n = samples.len();

        // Row-major and filled in row-major order, so the pass over the sample
        // structs is sequential in both source and destination.
        let matrix = Array2::from_shape_fn((n, CHANNELS), |(row, channel)| {
            let sample = &samples[row];
            if channel < MEASUREMENTS {
                sample.measurements[channel]
            } else {
                sample.manipulated[channel - MEASUREMENTS]
            }
        });
        let hours = Array1::from_shape_fn(n, |row| samples[row].hours);
        // i64 rather than u64: a step count never approaches the sign bit, and
        // signed is what pandas and scikit-learn expect to do arithmetic on.
        let steps = Array1::from_shape_fn(n, |row| samples[row].step as i64);
        let active = Array2::from_shape_fn((n, DISTURBANCES), |(row, idv)| {
            samples[row].labels.active[idv]
        });
        // NaN, not zero: a disturbance that never came on has no elapsed time,
        // and zero would read as "it came on this instant".
        let since_onset = Array2::from_shape_fn((n, DISTURBANCES), |(row, idv)| {
            samples[row].labels.since_onset[idv].unwrap_or(f64::NAN)
        });

        let matrix = frozen(py, matrix);
        // A view, so it inherits the read-only flag and shares the buffer.
        let transposed = matrix.transpose()?;

        Ok(Self {
            scenario: run.scenario,
            outcome: run.outcome,
            samples: n,
            matrix: matrix.unbind(),
            transposed: transposed.unbind(),
            hours: frozen(py, hours).unbind(),
            steps: frozen(py, steps).unbind(),
            active: frozen(py, active).unbind(),
            since_onset: frozen(py, since_onset).unbind(),
        })
    }

    /// One row of the cached transpose, which is one channel of the run.
    ///
    /// The index is checked by the callers, which know which of the three
    /// numbering schemes the caller used and can say so in the message. NumPy
    /// would raise its own `IndexError` here, naming a bound of 53 that means
    /// nothing to someone who asked for `XMV(13)`.
    fn channel<'py>(&self, py: Python<'py>, channel: usize) -> PyResult<Bound<'py, PyArray1<f64>>> {
        // Indexing a 2-D float64 array with an integer always gives a 1-D
        // float64 view, so the cast cannot fail; it is checked rather than
        // asserted because there is nothing to gain from panicking.
        self.transposed
            .bind(py)
            .as_any()
            .get_item(channel)?
            .cast_into::<PyArray1<f64>>()
            .map_err(PyErr::from)
    }
}

#[pymethods]
impl Run {
    /// What was asked for.
    #[getter]
    const fn scenario(&self) -> Scenario {
        Scenario::wrap(self.scenario)
    }

    /// How the run ended: `'completed'`, `'tripped'` or `'solve_failed'`.
    ///
    /// A trip does not by itself stop the run: `teprob.f:807-811` freezes the
    /// plant and keeps reporting, so a tripped run still has samples all the
    /// way to the end of its scenario unless `trip_ends_the_run` was set.
    ///
    /// `'solve_failed'` has no counterpart in the original, which cannot report
    /// it: `TESUB2` returns its guess and claims success. Delta D-001.
    #[getter]
    const fn outcome(&self) -> &'static str {
        match self.outcome {
            Outcome::Completed => "completed",
            Outcome::Tripped { .. } => "tripped",
            Outcome::SolveFailed { .. } => "solve_failed",
        }
    }

    /// The integrator step the plant tripped at, or `None` if it did not.
    #[getter]
    const fn tripped_at(&self) -> Option<usize> {
        match self.outcome {
            Outcome::Tripped { step, .. } => Some(step),
            _ => None,
        }
    }

    /// Simulated hours at the trip, or `None` if there was none.
    #[getter]
    const fn tripped_hours(&self) -> Option<f64> {
        match self.outcome {
            Outcome::Tripped { hours, .. } => Some(hours),
            _ => None,
        }
    }

    /// The first shutdown condition that fired, in words, or `None`.
    ///
    /// One of the eight `teprob.f:703-710` tests, named as it is there.
    #[getter]
    const fn trip_cause(&self) -> Option<&'static str> {
        match self.outcome {
            Outcome::Tripped {
                cause: Some(cause), ..
            } => Some(cause.describe()),
            _ => None,
        }
    }

    /// The step a solve failed at, or `None`. See `outcome`.
    #[getter]
    const fn solve_failed_at(&self) -> Option<usize> {
        match self.outcome {
            Outcome::SolveFailed { step } => Some(step),
            _ => None,
        }
    }

    /// The whole run as one `(n_samples, 53)` float64 array, read-only.
    ///
    /// Columns 0 to 40 are `XMEAS(1..41)` and 41 to 52 are `XMV(1..12)`, which
    /// is the layout `channel_names()` names and every downstream consumer
    /// uses. The array is C-contiguous.
    ///
    /// Returns the same object on every call and copies nothing. Use `.copy()`
    /// for something writable.
    fn to_numpy(&self, py: Python<'_>) -> Py<PyArray2<f64>> {
        self.matrix.clone_ref(py)
    }

    /// Simulated time at each sample, in hours: `(n_samples,)` float64,
    /// read-only.
    #[getter]
    fn hours(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.hours.clone_ref(py)
    }

    /// The integrator step each sample was taken at: `(n_samples,)` int64,
    /// read-only. One-based, as `temain_mod.f`'s loop counter is.
    #[getter]
    fn steps(&self, py: Python<'_>) -> Py<PyArray1<i64>> {
        self.steps.clone_ref(py)
    }

    /// One channel across the run, zero-based over the 53.
    ///
    /// A strided view into the matrix `to_numpy()` returns, not a copy. Pass it
    /// to `numpy.ascontiguousarray` if a contiguous buffer matters.
    ///
    /// Raises:
    ///     ValueError: If `channel` is not below 53.
    fn column<'py>(&self, py: Python<'py>, channel: usize) -> PyResult<Bound<'py, PyArray1<f64>>> {
        if channel >= CHANNELS {
            return Err(PyValueError::new_err(format!(
                "channel {channel} is out of range: a row has {CHANNELS} channels, \
                 numbered 0 to {last}",
                last = CHANNELS - 1
            )));
        }
        self.channel(py, channel)
    }

    /// One measurement across the run, one-based as `XMEAS(n)` is.
    ///
    /// A strided view, as `column` is.
    ///
    /// Raises:
    ///     ValueError: If `n` is outside 1 to 41.
    fn measurement<'py>(&self, py: Python<'py>, n: usize) -> PyResult<Bound<'py, PyArray1<f64>>> {
        if !(1..=MEASUREMENTS).contains(&n) {
            return Err(PyValueError::new_err(format!(
                "XMEAS index {n} is out of range: the plant reports \
                 {MEASUREMENTS} measurements, numbered 1 to {MEASUREMENTS}"
            )));
        }
        self.channel(py, n - 1)
    }

    /// One manipulated variable across the run, one-based as `XMV(n)` is.
    ///
    /// A strided view, as `column` is.
    ///
    /// Raises:
    ///     ValueError: If `n` is outside 1 to 12.
    fn manipulated<'py>(&self, py: Python<'py>, n: usize) -> PyResult<Bound<'py, PyArray1<f64>>> {
        if !(1..=MANIPULATED).contains(&n) {
            return Err(PyValueError::new_err(format!(
                "XMV index {n} is out of range: the plant has {MANIPULATED} \
                 manipulated variables, numbered 1 to {MANIPULATED}"
            )));
        }
        self.channel(py, MEASUREMENTS + n - 1)
    }

    /// Every channel by name, in row order.
    ///
    /// A dict of 53 strided views into the same buffer `to_numpy()` returns, so
    /// building it costs one dict and no array data. Insertion order is channel
    /// order, which makes it a `pandas.DataFrame` constructor argument as it
    /// stands.
    fn columns<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let out = PyDict::new(py);
        for (channel, name) in tepsim::channel_names().into_iter().enumerate() {
            out.set_item(name, self.channel(py, channel)?)?;
        }
        Ok(out)
    }

    /// Ground truth: what was actually wrong with the plant at each sample.
    ///
    /// A dict of two read-only `(n_samples, 20)` arrays, both indexed by
    /// `IDV(n) - 1`:
    ///
    /// * `'active'`, bool: whether that disturbance was on at that sample.
    /// * `'since_onset'`, float64: hours since it came on, NaN where it never
    ///   did.
    ///
    /// The original records nothing of the sort; a published dataset is a
    /// matrix and a filename, so detection-delay figures in the literature are
    /// computed against whatever onset the author assumed. `since_onset` is not
    /// simply time since the run began: the driver switches `IDV(12)` on at
    /// hour eight whatever the scenario asked for, so one onset can be later
    /// than the others and is not the caller's doing.
    fn labels<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let out = PyDict::new(py);
        out.set_item("active", self.active.clone_ref(py))?;
        out.set_item("since_onset", self.since_onset.clone_ref(py))?;
        Ok(out)
    }

    /// How many samples the run recorded.
    const fn __len__(&self) -> usize {
        self.samples
    }

    fn __repr__(&self) -> String {
        let ending = match self.outcome {
            Outcome::Completed => String::from("completed"),
            Outcome::Tripped { hours, cause, .. } => format!(
                "tripped at {hours:.3} h ({})",
                cause.map_or("no condition recorded", tepsim_core_describe)
            ),
            Outcome::SolveFailed { step } => format!("solve failed at step {step}"),
        };
        format!(
            "<tepsim.Run {} samples x {CHANNELS} channels, {:.1} h, {ending}>",
            self.samples, self.scenario.hours
        )
    }
}

/// `ShutdownCause::describe` as a function, so `Option::map_or` can name it.
fn tepsim_core_describe(cause: tepsim::tepsim_core::ShutdownCause) -> &'static str {
    cause.describe()
}
