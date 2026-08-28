//! Disturbances on a schedule, and a content hash that makes a dataset
//! self-describing.
//!
//! The original admits one shape of experiment: set some `IDV` flags before the
//! run and leave them. Every published Tennessee Eastman dataset is that, which
//! is why the literature's fault onsets are all at the same place and why
//! nobody studies a fault that arrives, persists and clears.
//!
//! This crate is the schedule. A `tepsim::Scenario` carries a [`Schedule`] of
//! [`Event`]s, each
//! at a time, each doing one thing. The simulation applies them as it passes
//! their times.
//!
//! # The content hash
//!
//! A [`Digest`] over the scenario's canonical form. Two scenarios that describe
//! the same experiment hash the same and two that do not, do not. A generated
//! dataset carries its scenario's digest, so a file says what produced it
//! rather than relying on a filename and a memory. `PLAN.org` calls this "the
//! piece that makes anomaly-detection results comparable across papers, which
//! the original never offered".
//!
//! # What this is not
//!
//! Not a general event system. The actions are the ones the plant can actually
//! take, and adding one means adding it to the plant too. That is deliberate: a
//! scenario that can express something the simulator cannot do would be a
//! scenario that sometimes silently does nothing.

#![no_std]
#![forbid(unsafe_code)]

pub mod digest;

pub use digest::Digest;

/// How many disturbances this model has. See `tepsim_core::FAULTS`.
pub const DISTURBANCES: usize = 20;

/// One thing that happens at one time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Action {
    /// Turn a disturbance on at full magnitude. One-based, as `IDV(n)` is.
    Start {
        /// Which disturbance, `1..=20`.
        fault: usize,
    },
    /// Turn a disturbance off.
    Stop {
        /// Which disturbance, `1..=20`.
        fault: usize,
    },
    /// Set a disturbance's magnitude.
    ///
    /// Anything other than exactly 0 or 1 needs the plant's
    /// `continuous_disturbances` extension, because
    /// `teprob.f:341-346` forces every flag to one of those two. A scenario
    /// that uses a fractional magnitude without the extension is rejected by
    /// [`Schedule::validate`] rather than silently rounded.
    SetMagnitude {
        /// Which disturbance, `1..=20`.
        fault: usize,
        /// How much of it, `0.0..=1.0`.
        magnitude: f64,
    },
    /// Move a controller setpoint. One-based, as `SETPT(n)` is.
    Setpoint {
        /// Which loop, `1..=20`.
        loop_index: usize,
        /// The new value, in that loop's own units.
        value: f64,
    },
}

impl Action {
    /// The disturbance this action touches, if any.
    #[must_use]
    pub const fn fault(&self) -> Option<usize> {
        match self {
            Self::Start { fault } | Self::Stop { fault } => Some(*fault),
            Self::SetMagnitude { fault, .. } => Some(*fault),
            Self::Setpoint { .. } => None,
        }
    }

    /// Whether this action needs the continuous-disturbance extension.
    #[must_use]
    pub fn needs_continuous_disturbances(&self) -> bool {
        match self {
            Self::SetMagnitude { magnitude, .. } => {
                // Exact comparisons on purpose. The question is whether this
                // magnitude is one of the two the original admits, and
                // 0.9999999999 is not one of them: it needs the extension just
                // as much as 0.5 does. An epsilon here would let a value the
                // plant cannot represent through the guard.
                magnitude.to_bits() != 0.0_f64.to_bits() && magnitude.to_bits() != 1.0_f64.to_bits()
            }
            _ => false,
        }
    }
}

/// An action and when it happens.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Event {
    /// Simulated hours from the start of the run.
    pub at_hours: f64,
    /// What happens.
    pub action: Action,
}

impl Event {
    /// An event at a time.
    #[must_use]
    pub const fn new(at_hours: f64, action: Action) -> Self {
        Self { at_hours, action }
    }

    /// Start a disturbance at a time.
    #[must_use]
    pub const fn start(at_hours: f64, fault: usize) -> Self {
        Self::new(at_hours, Action::Start { fault })
    }

    /// Stop one.
    #[must_use]
    pub const fn stop(at_hours: f64, fault: usize) -> Self {
        Self::new(at_hours, Action::Stop { fault })
    }
}

/// Why a scenario is not runnable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Invalid {
    /// A fault index outside `1..=20`.
    FaultOutOfRange {
        /// The index given.
        fault: usize,
    },
    /// A control loop index outside `1..=20`.
    LoopOutOfRange {
        /// The index given.
        loop_index: usize,
    },
    /// A magnitude outside `0.0..=1.0`, or not a number.
    MagnitudeOutOfRange,
    /// A fractional magnitude without the plant's `continuous_disturbances`
    /// extension.
    ///
    /// Rejected rather than rounded. Silently turning a request for half a
    /// fault into a whole one would produce a run that does not match its own
    /// description, which is exactly what the content hash exists to prevent.
    ContinuousDisturbancesNotEnabled,
    /// An event at a negative time, or at a time that is not a number.
    TimeNotFinite,
}

/// How many events a schedule can hold.
///
/// Fixed rather than growable, and the reason is `Copy`. A `Scenario` is a
/// small configuration value that callers build by chaining `with_*` methods
/// and pass around freely, and both the Python and the WebAssembly bindings
/// rely on that. A `Vec` inside it would take `Copy` away and turn every
/// `base.with_hours(2.0)` into a move or a clone, throughout this crate, both
/// binding layers, and every test.
///
/// Thirty-two is far more than any experiment in the literature uses: the
/// published datasets have exactly one event each, the fault being switched on
/// before the run.
pub const MAX_EVENTS: usize = 32;

/// A schedule of events.
///
/// Kept sorted by time, and stably: two events at the same instant happen in
/// the order they were added. That matters because `Stop` then `Start` on the
/// same fault at the same time is a different scenario from the reverse, and
/// the digest must distinguish them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Schedule {
    events: [Option<Event>; MAX_EVENTS],
    len: usize,
}

impl Default for Schedule {
    fn default() -> Self {
        Self::new()
    }
}

impl Schedule {
    /// An empty schedule.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            events: [None; MAX_EVENTS],
            len: 0,
        }
    }

    /// Add an event, keeping the schedule sorted.
    ///
    /// # Panics
    ///
    /// If the schedule already holds [`MAX_EVENTS`]. Silently dropping an
    /// event would produce a run that does not match its own description,
    /// which is the one thing the content hash exists to make impossible.
    pub fn add(&mut self, event: Event) {
        assert!(
            self.len < MAX_EVENTS,
            "a schedule holds at most {MAX_EVENTS} events"
        );
        // Insert after every event at or before this one's time, so events at
        // the same instant keep their insertion order.
        let mut at = self.len;
        for index in 0..self.len {
            if let Some(existing) = self.events[index]
                && existing.at_hours > event.at_hours
            {
                at = index;
                break;
            }
        }
        let mut index = self.len;
        while index > at {
            self.events[index] = self.events[index - 1];
            index -= 1;
        }
        self.events[at] = Some(event);
        self.len += 1;
    }

    /// Add an event and return the schedule, for building.
    #[must_use]
    pub fn with(mut self, event: Event) -> Self {
        self.add(event);
        self
    }

    /// The events, in time order.
    pub fn events(&self) -> impl Iterator<Item = Event> + '_ {
        self.events[..self.len].iter().filter_map(|e| *e)
    }

    /// How many.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.len
    }

    /// Whether there are none.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The events due in `(previous, now]`.
    ///
    /// Half-open at the start so an event is applied exactly once however the
    /// step size divides its time, and closed at the end so an event at time
    /// zero is applied on the first step rather than never.
    pub fn due(&self, previous_hours: f64, now_hours: f64) -> impl Iterator<Item = Event> + '_ {
        self.events()
            .filter(move |e| e.at_hours > previous_hours && e.at_hours <= now_hours)
    }

    /// Check every event.
    ///
    /// # Errors
    ///
    /// The first problem found, in time order.
    pub fn validate(&self, continuous_disturbances: bool) -> Result<(), Invalid> {
        for event in self.events() {
            if !event.at_hours.is_finite() || event.at_hours < 0.0 {
                return Err(Invalid::TimeNotFinite);
            }
            if let Some(fault) = event.action.fault()
                && !(1..=DISTURBANCES).contains(&fault)
            {
                return Err(Invalid::FaultOutOfRange { fault });
            }
            match event.action {
                Action::SetMagnitude { magnitude, .. } => {
                    if !magnitude.is_finite() || !(0.0..=1.0).contains(&magnitude) {
                        return Err(Invalid::MagnitudeOutOfRange);
                    }
                    if event.action.needs_continuous_disturbances() && !continuous_disturbances {
                        return Err(Invalid::ContinuousDisturbancesNotEnabled);
                    }
                }
                Action::Setpoint { loop_index, value } => {
                    if !(1..=20).contains(&loop_index) {
                        return Err(Invalid::LoopOutOfRange { loop_index });
                    }
                    if !value.is_finite() {
                        return Err(Invalid::MagnitudeOutOfRange);
                    }
                }
                Action::Start { .. } | Action::Stop { .. } => {}
            }
        }
        Ok(())
    }

    /// Feed this schedule into a digest.
    ///
    /// Every field of every event, in order, so nothing about the schedule can
    /// change without the digest changing.
    pub fn hash_into(&self, digest: &mut Digest) {
        digest.push_usize(self.len);
        for event in self.events() {
            digest.push_f64(event.at_hours);
            match event.action {
                Action::Start { fault } => {
                    digest.push_usize(1);
                    digest.push_usize(fault);
                }
                Action::Stop { fault } => {
                    digest.push_usize(2);
                    digest.push_usize(fault);
                }
                Action::SetMagnitude { fault, magnitude } => {
                    digest.push_usize(3);
                    digest.push_usize(fault);
                    digest.push_f64(magnitude);
                }
                Action::Setpoint { loop_index, value } => {
                    digest.push_usize(4);
                    digest.push_usize(loop_index);
                    digest.push_f64(value);
                }
            }
        }
    }
}

/// A scenario's digest, as the sixteen hex characters a filename can carry.
#[must_use]
pub fn short_hex(digest: Digest) -> [u8; 16] {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let value = digest.value();
    let mut out = [b'0'; 16];
    for (index, slot) in out.iter_mut().enumerate() {
        let shift = 60 - 4 * index;
        *slot = HEX[((value >> shift) & 0xf) as usize];
    }
    out
}
