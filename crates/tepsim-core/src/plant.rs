//! The plant, and the three-phase split of the original's right-hand side.
//!
//! # Why three phases
//!
//! `TEFUNC` (`teprob.f:196-816`) is not a pure function of `(t, y)`. It
//! advances the disturbance walks, draws measurement noise, ticks the sampled
//! analysers and latches valve commands, all inside what presents itself as a
//! derivative evaluation. That is harmless for fixed-step Euler, which
//! evaluates the right-hand side exactly once per step, and wrong for anything
//! else: RK4 would advance the walks four times per step and draw four sets of
//! noise.
//!
//! The impure work does not sit in one place, which is the part that is not
//! obvious from reading the plan. Some of it must happen *before* the
//! derivative and some can only happen *after* it:
//!
//! | Phase | `teprob.f` | What |
//! |---|---|---|
//! | [`Plant::advance_discrete`] | 341-406, 793-804 | `IDV` clamp, `IDVWLK` mapping, walk advance and spike draws, the `TIME=0` initialisation, and the valve-command latch |
//! | [`Plant::derivatives`] | 407-710, 762-792, 805 | the entire physical model, the noise-free measurements, the shutdown test, and the balances |
//! | [`Plant::sample_measurements`] | 711-761 | additive measurement noise, and the three sampled analysers with their dead time |
//!
//! The measurement vector is assembled at `teprob.f:679-701` out of flows the
//! pure evaluation computes, so the noise cannot be added before it; the walks
//! are read at `teprob.f:407-416`, so they must be advanced before it. One
//! impure call cannot sit on both sides. The pure phase therefore returns
//! [`Signals`] alongside the derivative, so the post-phase does not have to
//! recompute the model to find out what to add noise to.
//!
//! See `PLAN.org`, "Splitting the right-hand side", and the decision entry of
//! 2026-08-11 in `BACKLOG.org`.
//!
//! # The valve-command latch is hoisted
//!
//! `teprob.f:793-798` sets `IVST` from `IDV` and `799-804` latches `VCV` from
//! `XMV`, at the very end of the routine. It reads only `XMV`, `VST`, `IVST`,
//! `IDV` and `TIME`. Nothing in `345-792` writes any of them, so moving the
//! block into the pre-phase changes no number. That is a claim about four
//! hundred and fifty lines of Fortran, so it is checked mechanically rather
//! than by eye; see `crates/tepsim-oracle/tests/hoist_valve_latch.rs`.
//!
//! The latch shares the `DO 9020` loop with the valve derivative at
//! `teprob.f:805`, so the port splits that loop: the latch goes to the
//! pre-phase and `YP(I+38)` stays in the pure one.
//!
//! # Status
//!
//! Structure only. [`Plant::derivatives`] returns
//! [`PlantError::NotImplemented`] until the physics lands over B-0017 to
//! B-0025. It deliberately does *not* return zeros: an all-zero derivative is a
//! perfectly valid answer describing a frozen plant, and would let a later item
//! pass a test it should have failed.

use crate::state::{Derivative, State};
use crate::thermo::TemperatureError;
use crate::variables::MeasIndex;

/// How many continuous measurements the plant produces, `XMEAS(1..22)`.
pub const N_CONTINUOUS: usize = 22;

/// How many sampled composition measurements, `XMEAS(23..41)`.
pub const N_SAMPLED: usize = 19;

/// Everything the controller drives: `XMV(1..12)` and `IDV(1..20)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Inputs {
    /// The twelve manipulated variables, as percentages of full travel.
    pub manipulated: [f64; 12],
    /// The twenty disturbance magnitudes.
    ///
    /// The original clamps each to exactly 0 or 1 on every call
    /// (`teprob.f:341-344`), which forbids partial-magnitude faults. That is a
    /// Class C quirk: this type carries a magnitude so the eventual fix does
    /// not need a new type, and the faithful path clamps.
    pub disturbances: [f64; 20],
}

impl Default for Inputs {
    fn default() -> Self {
        Self {
            manipulated: [0.0; 12],
            disturbances: [0.0; 20],
        }
    }
}

/// The noise-free signals the pure phase produces for the post-phase.
///
/// Not measurements yet: no noise has been added and the analysers have not
/// been sampled. Kept separate so that [`Plant::sample_measurements`] does not
/// have to re-run the model to find out what to add noise to.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Signals {
    /// `XMEAS(1..22)` before noise (`teprob.f:679-701`).
    pub continuous: [f64; N_CONTINUOUS],
    /// `XCMP(23..41)`, the analyser inputs (`teprob.f:717-735`).
    pub compositions: [f64; N_SAMPLED],
}

impl Default for Signals {
    fn default() -> Self {
        Self {
            continuous: [0.0; N_CONTINUOUS],
            compositions: [0.0; N_SAMPLED],
        }
    }
}

/// `XMEAS(1..41)`, with noise added and the analysers sampled.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Measurements([f64; N_CONTINUOUS + N_SAMPLED]);

impl Measurements {
    /// All forty-one, in `XMEAS` order.
    #[must_use]
    pub const fn as_array(&self) -> &[f64; N_CONTINUOUS + N_SAMPLED] {
        &self.0
    }
}

impl Default for Measurements {
    fn default() -> Self {
        Self([0.0; N_CONTINUOUS + N_SAMPLED])
    }
}

impl core::ops::Index<MeasIndex> for Measurements {
    type Output = f64;

    fn index(&self, index: MeasIndex) -> &f64 {
        &self.0[index.zero_based()]
    }
}

/// Why a derivative evaluation could not produce an answer.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlantError {
    /// The physics has not been ported yet. Removed when B-0025 lands.
    NotImplemented,
    /// A temperature solve failed. See [`TemperatureError`] and delta D-001.
    Temperature(TemperatureError),
}

impl core::fmt::Display for PlantError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NotImplemented => f.write_str("the plant model is not ported yet"),
            Self::Temperature(e) => write!(f, "temperature solve failed: {e}"),
        }
    }
}

impl core::error::Error for PlantError {}

impl From<TemperatureError> for PlantError {
    fn from(error: TemperatureError) -> Self {
        Self::Temperature(error)
    }
}

/// Simulation time, in hours, matching the original's `TIME`.
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Default)]
pub struct SimTime(pub f64);

impl SimTime {
    /// The elapsed hours.
    #[must_use]
    pub const fn hours(self) -> f64 {
        self.0
    }
}

/// The plant: everything the original keeps in `COMMON` that survives between
/// calls and is not part of the integrated state.
///
/// Owned, `Send` and `Clone`, so a process can run as many of these as it
/// likes. The original supports exactly one, which is a large part of why this
/// port exists.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Plant {
    /// The latched valve commands, `VCV` (`teprob.f:799-804`).
    valve_command: [f64; 12],
}

impl Plant {
    /// A plant with nothing latched yet.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// The latched valve commands, `VCV`.
    #[must_use]
    pub const fn valve_command(&self) -> &[f64; 12] {
        &self.valve_command
    }

    /// **Impure.** Exactly once per outer step, *before* any derivative call.
    ///
    /// Advances the disturbance walks, draws the channel 10-12 spikes, and
    /// latches the valve commands. See the module documentation for why this
    /// cannot be merged with [`Plant::sample_measurements`].
    ///
    /// Stubbed: the walk machinery is Phase 3 and the latch is B-0025.
    pub fn advance_discrete(&mut self, _t: SimTime, _u: &Inputs) {}

    /// **Pure.** The derivative of the state, and the noise-free signals.
    ///
    /// Takes `&self`, so the compiler enforces that it cannot mutate the plant;
    /// the crate is `#![forbid(unsafe_code)]` and holds no statics, so there is
    /// no interior mutability or global state either. Purity is therefore a
    /// property of the signature, not a convention.
    ///
    /// # Errors
    ///
    /// [`PlantError::NotImplemented`] until B-0025. After that, a failed
    /// temperature solve (delta D-001).
    pub fn derivatives(
        &self,
        _t: SimTime,
        _y: &State,
        _u: &Inputs,
    ) -> Result<(Derivative, Signals), PlantError> {
        Err(PlantError::NotImplemented)
    }

    /// **Impure.** Exactly once per outer step, *after* the derivative call.
    ///
    /// Adds measurement noise and ticks the three sampled analysers.
    ///
    /// Stubbed: this is Phase 3.
    pub fn sample_measurements(&mut self, _t: SimTime, _signals: &Signals) -> Measurements {
        Measurements::default()
    }

    /// One explicit Euler step, sequencing the three phases in the only order
    /// that reproduces the original.
    ///
    /// This is the whole point of the split written down as code: the impure
    /// pre-phase runs once, the pure phase runs as many times as the integrator
    /// needs (once, here), and the impure post-phase runs once. An RK4 driver
    /// differs only in calling [`Plant::derivatives`] four times between the
    /// same two impure calls.
    ///
    /// # Errors
    ///
    /// Whatever [`Plant::derivatives`] returns.
    pub fn euler_step(
        &mut self,
        t: SimTime,
        y: &State,
        u: &Inputs,
        dt: f64,
    ) -> Result<(State, Measurements), PlantError> {
        self.advance_discrete(t, u);
        let (derivative, signals) = self.derivatives(t, y, u)?;
        let measurements = self.sample_measurements(t, &signals);
        Ok((y.step(dt, &derivative), measurements))
    }
}

/// Assert that [`Plant::derivatives`] is a pure function at this point.
///
/// Evaluates twice and requires bit-identical results and an unchanged plant.
/// This is the property the whole three-phase split exists to provide and the
/// one Tier 2 rests on, so every item that adds physics should call it.
///
/// It is deliberately in the library rather than in a test file: it is a
/// contract the plant offers, and it should be as easy to check from a
/// downstream crate as from here.
///
/// # Panics
///
/// If two evaluations differ in any bit, or if the plant changed.
pub fn assert_derivatives_are_pure(plant: &Plant, t: SimTime, y: &State, u: &Inputs) {
    let before = plant.clone();
    let first = plant.derivatives(t, y, u);
    let second = plant.derivatives(t, y, u);

    assert!(
        plant == &before,
        "evaluating the derivative changed the plant, so it is not pure"
    );

    match (first, second) {
        (Ok((d1, s1)), Ok((d2, s2))) => {
            for (slot, (a, b)) in d1.to_flat().iter().zip(d2.to_flat().iter()).enumerate() {
                assert_eq!(
                    a.to_bits(),
                    b.to_bits(),
                    "derivative slot {slot} differs between two evaluations at \
                     the same point: {a:?} then {b:?}"
                );
            }
            for (i, (a, b)) in s1
                .continuous
                .iter()
                .zip(s2.continuous.iter())
                .chain(s1.compositions.iter().zip(s2.compositions.iter()))
                .enumerate()
            {
                assert_eq!(
                    a.to_bits(),
                    b.to_bits(),
                    "signal {i} differs between two evaluations: {a:?} then {b:?}"
                );
            }
        }
        (Err(a), Err(b)) => assert_eq!(a, b, "the same point failed two different ways"),
        (a, b) => panic!("one evaluation succeeded and the other did not: {a:?} then {b:?}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_derivative_is_pure() {
        // Currently weak: there is no physics, so both evaluations fail the
        // same way. It strengthens by itself as B-0017 onward land, which is
        // why the assertion exists now rather than later.
        let plant = Plant::new();
        assert_derivatives_are_pure(&plant, SimTime(0.0), &State::default(), &Inputs::default());
    }

    #[test]
    fn the_unported_plant_says_so_rather_than_returning_zeros() {
        let plant = Plant::new();
        let outcome = plant.derivatives(SimTime(0.0), &State::default(), &Inputs::default());
        assert_eq!(outcome, Err(PlantError::NotImplemented));
    }

    #[test]
    fn a_euler_step_propagates_the_failure_rather_than_stepping() {
        let mut plant = Plant::new();
        let outcome = plant.euler_step(
            SimTime(0.0),
            &State::default(),
            &Inputs::default(),
            1.0 / 3600.0,
        );
        assert_eq!(outcome, Err(PlantError::NotImplemented));
    }

    #[test]
    fn a_plant_is_send_and_clonable_unlike_the_original() {
        fn assert_send<T: Send>() {}
        assert_send::<Plant>();
        let plant = Plant::new();
        assert_eq!(plant.clone(), plant);
    }

    #[test]
    fn measurements_are_indexed_by_name_not_by_number() {
        let measurements = Measurements::default();
        assert_eq!(measurements.as_array().len(), 41);
        // Indexing by `MeasIndex` must reach the same slot the array does.
        for one_based in 1..=MeasIndex::COUNT {
            let index = MeasIndex::new(one_based).expect("in range");
            assert_eq!(measurements[index].to_bits(), 0.0_f64.to_bits());
        }
    }
}
