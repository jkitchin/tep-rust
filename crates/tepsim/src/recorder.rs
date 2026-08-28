//! Where samples go.
//!
//! A [`Recorder`] receives every sample as it is produced. That is deliberately
//! not the same shape as "return a `Vec` at the end": a 48-hour run at every
//! step is 172,800 rows of 53 channels, which is 73 MB, and a browser rendering
//! a live trace wants the last few hundred rather than all of them.
//!
//! # Choosing between them
//!
//! | Sink | For |
//! |---|---|
//! | [`Columnar`] | analysis. One `Vec<f64>` per channel, which is what a numpy or Arrow view wants |
//! | [`Csv`] | teaching, and anything that has to be read by eye |
//! | [`Ring`] | a live display: keeps the last `n` samples and nothing else |
//! | [`Decimating`] | wrapping any of the above to keep one sample in `k` |
//! | [`Selecting`] | wrapping any of the above to keep only some channels |
//! | `()` | counting the run without keeping it |
//!
//! The wrappers compose: `Decimating::new(Selecting::new(Columnar::new(), ..), 10)`.
//!
//! # Why the trait is `no_std`
//!
//! `tepsim` is `no_std + alloc`, because the same code has to compile to
//! wasm32 for the browser. [`Csv`] therefore writes through
//! [`core::fmt::Write`] rather than `std::io::Write`; the CLI adapts one to the
//! other in a few lines.

use alloc::string::String;
use alloc::vec::Vec;

use crate::run::{CHANNELS, Sample, channel_names};

/// Something that receives samples.
///
/// Implemented for `()`, which counts them and keeps nothing: useful for
/// timing a run, and as the identity when a caller has no sink.
pub trait Recorder {
    /// Called once per recorded sample, in order.
    fn record(&mut self, sample: &Sample);

    /// Called once when the run ends, so a buffered sink can flush.
    ///
    /// Default does nothing.
    fn finish(&mut self) {}
}

impl Recorder for () {
    fn record(&mut self, _sample: &Sample) {}
}

impl<R: Recorder + ?Sized> Recorder for &mut R {
    fn record(&mut self, sample: &Sample) {
        (**self).record(sample);
    }
    fn finish(&mut self) {
        (**self).finish();
    }
}

/// One `Vec<f64>` per channel.
///
/// The layout every numerical consumer wants: a numpy view, an Arrow column,
/// the correlation matrix, a detector. Storing by channel rather than by sample
/// means a caller that wants one variable does not stride across 52 others.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Columnar {
    columns: Vec<Vec<f64>>,
    steps: Vec<usize>,
    hours: Vec<f64>,
}

impl Columnar {
    /// An empty recorder.
    #[must_use]
    pub fn new() -> Self {
        Self {
            columns: (0..CHANNELS).map(|_| Vec::new()).collect(),
            steps: Vec::new(),
            hours: Vec::new(),
        }
    }

    /// Pre-allocate for a known number of samples.
    ///
    /// Worth doing: a 48-hour run at every step would otherwise grow and copy
    /// 53 vectors seventeen times each.
    #[must_use]
    pub fn with_capacity(samples: usize) -> Self {
        Self {
            columns: (0..CHANNELS).map(|_| Vec::with_capacity(samples)).collect(),
            steps: Vec::with_capacity(samples),
            hours: Vec::with_capacity(samples),
        }
    }

    /// How many samples were recorded.
    #[must_use]
    pub fn len(&self) -> usize {
        self.steps.len()
    }

    /// Whether nothing was recorded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.steps.is_empty()
    }

    /// One channel, zero-based over the 53.
    ///
    /// # Panics
    ///
    /// If `channel` is not below [`CHANNELS`].
    #[must_use]
    pub fn column(&self, channel: usize) -> &[f64] {
        assert!(channel < CHANNELS, "channel {channel} is out of range");
        &self.columns[channel]
    }

    /// Every channel.
    #[must_use]
    pub fn columns(&self) -> &[Vec<f64>] {
        &self.columns
    }

    /// The step number of each sample.
    #[must_use]
    pub fn steps(&self) -> &[usize] {
        &self.steps
    }

    /// The simulated time of each sample, in hours.
    #[must_use]
    pub fn hours(&self) -> &[f64] {
        &self.hours
    }
}

impl Recorder for Columnar {
    fn record(&mut self, sample: &Sample) {
        let row = sample.row();
        for (column, value) in self.columns.iter_mut().zip(row) {
            column.push(value);
        }
        self.steps.push(sample.step);
        self.hours.push(sample.hours);
    }
}

/// Comma-separated values, written through [`core::fmt::Write`].
///
/// # Precision
///
/// Seventeen significant digits, which is what round-trips an `f64` exactly.
/// Fewer would make a recorded dataset only approximately reproducible, and the
/// whole point of a deterministic simulator is that it is exactly so.
#[derive(Debug)]
pub struct Csv<W> {
    out: W,
    labels: bool,
    header_written: bool,
    error: Option<core::fmt::Error>,
}

impl<W: core::fmt::Write> Csv<W> {
    /// A CSV sink writing the 53 channels.
    pub const fn new(out: W) -> Self {
        Self {
            out,
            labels: false,
            header_written: false,
            error: None,
        }
    }

    /// Also write the ground-truth columns.
    #[must_use]
    pub const fn with_labels(mut self) -> Self {
        self.labels = true;
        self
    }

    /// The first write error, if there was one.
    ///
    /// [`Recorder::record`] cannot return a result, so errors are held and
    /// reported here. A sink that silently dropped rows would be worse than
    /// one that stopped.
    #[must_use]
    pub const fn error(&self) -> Option<core::fmt::Error> {
        self.error
    }

    /// The writer back.
    pub fn into_inner(self) -> W {
        self.out
    }

    fn write_header(&mut self) -> core::fmt::Result {
        write!(self.out, "step,hours")?;
        for name in channel_names() {
            write!(self.out, ",{name}")?;
        }
        if self.labels {
            write!(self.out, ",fault,hours_since_onset")?;
        }
        writeln!(self.out)
    }

    fn write_sample(&mut self, sample: &Sample) -> core::fmt::Result {
        write!(self.out, "{},{:.6}", sample.step, sample.hours)?;
        for value in sample.row() {
            write!(self.out, ",{value:.17e}")?;
        }
        if self.labels {
            write!(self.out, ",")?;
            let mut first = true;
            for fault in sample.labels.faults() {
                if !first {
                    write!(self.out, " ")?;
                }
                write!(self.out, "{fault}")?;
                first = false;
            }
            write!(self.out, ",")?;
            if let Some(hours) = sample
                .labels
                .faults()
                .next()
                .and_then(|n| sample.labels.since_onset[n - 1])
            {
                write!(self.out, "{hours:.6}")?;
            }
        }
        writeln!(self.out)
    }
}

impl<W: core::fmt::Write> Recorder for Csv<W> {
    fn record(&mut self, sample: &Sample) {
        if self.error.is_some() {
            return;
        }
        let result = (|| {
            if !self.header_written {
                self.write_header()?;
                self.header_written = true;
            }
            self.write_sample(sample)
        })();
        if let Err(error) = result {
            self.error = Some(error);
        }
    }
}

/// A `String` sink, for tests and for the browser.
pub type CsvString = Csv<String>;

/// The last `capacity` samples and nothing else.
///
/// What a live display wants: bounded memory however long the run, and the
/// recent history is the only part on screen. The browser app runs the
/// simulation in a Web Worker at many times real time, so keeping everything
/// would exhaust memory long before the run ended.
#[derive(Clone, Debug, PartialEq)]
pub struct Ring {
    samples: Vec<Sample>,
    capacity: usize,
    /// Where the next write goes once the buffer is full.
    next: usize,
    /// How many samples have been offered in total, including those overwritten.
    seen: usize,
}

impl Ring {
    /// A ring holding at most `capacity` samples.
    ///
    /// # Panics
    ///
    /// If `capacity` is zero, which would silently discard everything.
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        assert!(capacity > 0, "a ring of zero capacity records nothing");
        Self {
            samples: Vec::with_capacity(capacity),
            capacity,
            next: 0,
            seen: 0,
        }
    }

    /// How many are held.
    #[must_use]
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Whether nothing is held.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// How many samples were offered, including those overwritten.
    #[must_use]
    pub const fn seen(&self) -> usize {
        self.seen
    }

    /// The held samples, oldest first.
    pub fn iter(&self) -> impl Iterator<Item = &Sample> {
        let split = if self.samples.len() < self.capacity {
            0
        } else {
            self.next
        };
        self.samples[split..].iter().chain(&self.samples[..split])
    }
}

impl Recorder for Ring {
    fn record(&mut self, sample: &Sample) {
        self.seen += 1;
        if self.samples.len() < self.capacity {
            self.samples.push(*sample);
        } else {
            self.samples[self.next] = *sample;
            self.next = (self.next + 1) % self.capacity;
        }
    }
}

/// Keep one sample in every `factor`, and pass the rest on to nothing.
///
/// Downsampling belongs at the sink rather than in the run loop: the loop's
/// cadence is part of the scenario and changing it changes the *simulation*,
/// whereas this changes only what is kept. A caller wanting an hourly summary
/// of a 48-hour run should not have to run the plant differently to get one.
#[derive(Clone, Debug, PartialEq)]
pub struct Decimating<R> {
    inner: R,
    factor: usize,
    seen: usize,
}

impl<R: Recorder> Decimating<R> {
    /// Keep one sample in `factor`, starting with the first.
    ///
    /// # Panics
    ///
    /// If `factor` is zero.
    #[must_use]
    pub const fn new(inner: R, factor: usize) -> Self {
        assert!(factor > 0, "a factor of zero would keep nothing");
        Self {
            inner,
            factor,
            seen: 0,
        }
    }

    /// The wrapped sink.
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// The wrapped sink, by reference.
    pub const fn inner(&self) -> &R {
        &self.inner
    }
}

impl<R: Recorder> Recorder for Decimating<R> {
    fn record(&mut self, sample: &Sample) {
        // The first sample is kept, not the `factor`-th. A caller asking for
        // one in ten wants ten percent of the data starting at the beginning,
        // not a series that starts nine samples late.
        if self.seen % self.factor == 0 {
            self.inner.record(sample);
        }
        self.seen += 1;
    }

    fn finish(&mut self) {
        self.inner.finish();
    }
}

/// Keep only some channels, zeroing the rest.
///
/// The sample type is a fixed 53 channels, so this blanks what was not asked
/// for rather than changing the shape. That keeps every downstream index
/// meaning the same thing, which matters more than the memory: a consumer that
/// selected channels and then indexed by position would be silently reading the
/// wrong variable.
#[derive(Clone, Debug, PartialEq)]
pub struct Selecting<R> {
    inner: R,
    keep: [bool; CHANNELS],
}

impl<R: Recorder> Selecting<R> {
    /// Keep the given channels, zero-based over the 53.
    ///
    /// # Panics
    ///
    /// If any index is out of range, or if none is given.
    #[must_use]
    pub fn new(inner: R, channels: &[usize]) -> Self {
        assert!(
            !channels.is_empty(),
            "selecting no channels records nothing"
        );
        let mut keep = [false; CHANNELS];
        for channel in channels {
            assert!(*channel < CHANNELS, "channel {channel} is out of range");
            keep[*channel] = true;
        }
        Self { inner, keep }
    }

    /// The wrapped sink.
    pub fn into_inner(self) -> R {
        self.inner
    }

    /// The wrapped sink, by reference.
    pub const fn inner(&self) -> &R {
        &self.inner
    }
}

impl<R: Recorder> Recorder for Selecting<R> {
    fn record(&mut self, sample: &Sample) {
        let mut masked = *sample;
        for (index, keep) in self.keep.iter().enumerate() {
            if !keep {
                if index < crate::run::MEASUREMENTS {
                    masked.measurements[index] = 0.0;
                } else {
                    masked.manipulated[index - crate::run::MEASUREMENTS] = 0.0;
                }
            }
        }
        self.inner.record(&masked);
    }

    fn finish(&mut self) {
        self.inner.finish();
    }
}
