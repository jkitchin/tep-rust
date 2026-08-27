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
//! The pure phase is complete as of B-0025: [`Plant::derivatives`] runs the
//! whole model, from unpacking the state to the fifty balances. The two impure
//! phases are still stubs and land with Phase 3, which brings the RNG and the
//! disturbance walks.

use crate::balances::{CoolantInlet, QuirkFixes, VALVE_STICTION, balances};
use crate::equilibrium::equilibrium;
use crate::flows::{FlowDrift, flows};
use crate::heat::{HeatDrift, heat_transfer};
use crate::kinetics::{ReactionDrift, kinetics};
use crate::measurements::{Shutdown, measurements};
use crate::rng::TepRng;
use crate::state::{Derivative, State};
use crate::streams::{FeedConditions, streams};
use crate::stripper::stripper;
use crate::thermo::TemperatureError;
use crate::variables::MeasIndex;
use crate::vessels::{TemperatureSeeds, unpack};
use crate::walk::{Walks, advance};

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

impl Inputs {
    /// The twenty disturbance flags as `teprob.f:341-344` leaves them.
    ///
    /// ```fortran
    ///       IF(IDV(I).GT.0)THEN
    ///       IDV(I)=1
    ///       ELSE
    ///       IDV(I)=0
    ///       ENDIF
    /// ```
    ///
    /// *Any* positive magnitude becomes exactly one, so this is a threshold
    /// and not a rounding: 0.4 gives 1, not 0. That is the Class C quirk this
    /// type's `disturbances` field exists to make fixable later, and until
    /// there is a sign-off the faithful path is the only path.
    #[must_use]
    pub fn clamped_disturbances(&self) -> [f64; 20] {
        let mut out = [0.0; 20];
        for (slot, raw) in out.iter_mut().zip(self.disturbances) {
            *slot = if raw > 0.0 { 1.0 } else { 0.0 };
        }
        out
    }
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
    /// Whether the plant is down, and why (`teprob.f:702-710`).
    ///
    /// Carried here rather than left to the derivative because the post-phase
    /// needs it: `teprob.f:711` skips the noise entirely when `ISD` is set, so
    /// [`Plant::sample_measurements`] cannot do its job without knowing.
    pub shutdown: Shutdown,
}

impl Default for Signals {
    fn default() -> Self {
        Self {
            continuous: [0.0; N_CONTINUOUS],
            compositions: [0.0; N_SAMPLED],
            shutdown: Shutdown::default(),
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
    /// A temperature solve failed. See [`TemperatureError`] and delta D-001.
    ///
    /// The only way an evaluation can fail. The original cannot fail at all:
    /// `TESUB2` returns its guess and reports success, which is what delta
    /// D-001 is about.
    Temperature(TemperatureError),
}

impl core::fmt::Display for PlantError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
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
#[derive(Clone, Debug, PartialEq)]
pub struct Plant {
    /// The latched valve commands, `VCV` (`teprob.f:799-804`).
    valve_command: [f64; 12],
    /// The four Newton warm-start temperatures.
    ///
    /// Genuine persistent state, not a cache: `TESUB2` takes its temperature
    /// as both guess and answer, so every evaluation starts from the previous
    /// one's result. See [`mod@crate::vessels`].
    seeds: TemperatureSeeds,
    /// Which Class C quirks are fixed rather than reproduced. All off by
    /// default; see [`QuirkFixes`].
    pub quirks: QuirkFixes,
    /// The walk-driven inputs the pure phase reads (`teprob.f:407-416`).
    ///
    /// Recomputed by [`Plant::advance_discrete`] on every step.
    walks: WalkInputs,
    /// The twelve disturbance channels.
    channels: Walks,
    /// The generator. The *only* place it lives, and it moves only in
    /// [`Plant::advance_discrete`].
    rng: TepRng,
}

/// Everything `teprob.f:407-416` reads out of the disturbance walks.
///
/// Grouped so that the pure phase takes one argument rather than five, and so
/// that Phase 3 has one place to fill in.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WalkInputs {
    /// `XST(1..3,4)`, `TST(1)` and `TST(4)`.
    pub feed: FeedConditions,
    /// `TESUB8(9)` and `TESUB8(12)`.
    pub flow: FlowDrift,
    /// `TESUB8(10)` and `TESUB8(11)`.
    pub heat: HeatDrift,
    /// `TESUB8(7)` and `TESUB8(8)`.
    pub reaction: ReactionDrift,
    /// `TCWR` and `TCWS`.
    pub coolant: CoolantInlet,
}

impl Default for Plant {
    fn default() -> Self {
        Self {
            valve_command: [0.0; 12],
            seeds: TemperatureSeeds::default(),
            quirks: QuirkFixes::default(),
            walks: WalkInputs::default(),
            channels: Walks::default(),
            rng: TepRng::with_default_seed(),
        }
    }
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
    /// Advances the walks, reads what they drive, and latches the valve
    /// commands.
    ///
    /// This is the only place the generator moves, which is the promise the
    /// three-phase split makes: an RK4 driver calls [`Plant::derivatives`]
    /// four times between two of these and the disturbances advance once.
    // @port teprob.f:340-416, 793-804
    pub fn advance_discrete(&mut self, t: SimTime, u: &Inputs) {
        let idv = u.clamped_disturbances();

        // teprob.f:347-406.
        advance(&mut self.channels, &mut self.rng, t.hours(), &idv);

        // teprob.f:407-416. Reading the channels is `TESUB8`.
        let channel = |n: usize| self.channels.channels[n - 1].at(t.hours());
        // teprob.f:407-408. Two step disturbances on the same component, and
        // the second of them also moves B on the next line, so IDV(2) shifts
        // A down and B up while IDV(1) only shifts A down.
        let a = channel(1) - idv[0] * 0.03 - idv[1] * 2.43719e-3;
        let b = channel(2) + idv[1] * 0.005;
        self.walks = WalkInputs {
            feed: FeedConditions {
                // teprob.f:410. C is the remainder of one, never a draw.
                ac_feed_light: [a, b, 1.0 - a - b],
                // teprob.f:411-412
                d_feed_celsius: channel(3) + idv[2] * 5.0,
                ac_feed_celsius: channel(4),
            },
            // teprob.f:572 and 583.
            flow: FlowDrift {
                steam_capacity: channel(9),
                reactor_outlet: channel(12),
            },
            // teprob.f:673 and 676.
            heat: HeatDrift {
                reactor_coolant: channel(10),
                condenser_coolant: channel(11),
            },
            // teprob.f:415-416.
            reaction: ReactionDrift {
                first: channel(7),
                second: channel(8),
            },
            // teprob.f:413-414.
            coolant: CoolantInlet {
                reactor: channel(5) + idv[3] * 5.0,
                condenser: channel(6) + idv[4] * 5.0,
            },
        };

        self.latch_valves(t, u, &idv);
    }

    /// The generator word. Moves only in [`Plant::advance_discrete`].
    #[must_use]
    pub const fn rng(&self) -> f64 {
        self.rng.state()
    }

    /// Set the generator word, for a harness placing the plant in a known
    /// condition, or a scenario choosing its seed.
    pub const fn set_rng(&mut self, g: f64) {
        self.rng = TepRng::new(g);
    }

    /// The twelve disturbance channels.
    #[must_use]
    pub const fn channels(&self) -> &Walks {
        &self.channels
    }

    /// Set the disturbance channels.
    pub const fn set_channels(&mut self, channels: Walks) {
        self.channels = channels;
    }

    // @port teprob.f:793-804
    fn latch_valves(&mut self, t: SimTime, u: &Inputs, idv: &[f64; 20]) {
        // teprob.f:793-798. Only six valves can stick, and only under three
        // specific disturbances. `IVST` is zero for the rest, which makes
        // their threshold zero and so makes them track exactly.
        let idv = |n: usize| idv[n - 1];
        let mut sticking = [0.0; 12];
        sticking[9] = idv(14);
        sticking[10] = idv(15);
        for valve in [4, 6, 7, 8] {
            sticking[valve] = idv(19);
        }

        // teprob.f:799-804. The command is latched unless it has moved
        // further than the valve's stiction band, and then clamped to travel.
        for ((command, stuck), wanted) in self
            .valve_command
            .iter_mut()
            .zip(sticking)
            .zip(u.manipulated)
        {
            let threshold = VALVE_STICTION * stuck;
            if t.hours() == 0.0 || (*command - wanted).abs() > threshold {
                *command = wanted;
            }
            // teprob.f:803-804, written as the two guards the listing writes
            // rather than as a `clamp`, so it can be checked a line at a time.
            // Identical for every finite input; `VCV` is never NaN here
            // because it is either a manipulated variable or a previous
            // clamped value.
            *command = command.clamp(0.0, 100.0);
        }
    }

    /// **Pure.** The derivative of the state, and the noise-free signals.
    ///
    /// Takes `&self`, so the compiler enforces that it cannot mutate the plant;
    /// the crate is `#![forbid(unsafe_code)]` and holds no statics, so there is
    /// no interior mutability or global state either. Purity is therefore a
    /// property of the signature, not a convention.
    ///
    /// # Errors
    ///
    /// [`PlantError::Temperature`] if any of the four Newton solves runs out
    /// of iterations. The original silently returns its guess instead; that is
    /// delta D-001, and it has never fired on the physical domain.
    ///
    /// # Note on the warm start
    ///
    /// This takes `&self`, so it cannot carry the converged temperatures
    /// forward: that would make it impure and defeat the whole split. The
    /// seeds it uses are the ones [`Plant::step_seeds`] last stored, which the
    /// integrator advances once per step. Evaluating twice at the same point
    /// therefore gives bit-identical answers, which is what
    /// [`assert_derivatives_are_pure`] checks.
    // @port teprob.f:407-710, 762-792, 805
    pub fn derivatives(
        &self,
        _t: SimTime,
        y: &State,
        u: &Inputs,
    ) -> Result<(Derivative, Signals), PlantError> {
        let w = &self.walks;
        let unpacked = unpack(y, self.seeds)?;
        let eq = equilibrium(&unpacked);
        let mut table = streams(&unpacked, &eq, &w.feed);

        // `teprob.f:341-344` clamps these before use. The pre-phase has
        // already done it for its own purposes; doing it again here keeps
        // `derivatives` a function of its arguments rather than of what the
        // pre-phase happened to store.
        let idv = u.clamped_disturbances();

        let mut flow = flows(y, &unpacked, &eq, &table, &idv, w.flow);
        let _ = stripper(&mut table, &mut flow, unpacked.stripper.celsius);
        let heat = heat_transfer(y, &unpacked, &table, &flow, w.heat);
        let kin = kinetics(&eq.reactor, unpacked.reactor.kelvin(), w.reaction);
        let pressures = (
            eq.reactor.pressure,
            eq.separator.pressure,
            eq.mixing_pressure,
        );
        let measured = measurements(y, &unpacked, &table, &flow, &heat, pressures);
        let assembled = balances(
            y,
            &table,
            &flow,
            &kin,
            &heat,
            measured.shutdown,
            w.coolant,
            &self.valve_command,
            self.quirks,
        );

        Ok((
            assembled.derivative,
            Signals {
                continuous: measured.continuous,
                // `XCMP` is B-0024b; it needs nothing this phase does not
                // already have, but it belongs with the analysers that read it.
                compositions: [0.0; N_SAMPLED],
                shutdown: measured.shutdown,
            },
        ))
    }

    /// Advance the Newton warm-start temperatures, once per outer step.
    ///
    /// The original carries `TCR`, `TCS`, `TCC` and `TCV` in `COMMON`, so
    /// every evaluation seeds itself from the previous one's answer. That
    /// makes `TEFUNC` path-dependent in the last bits, and reproducing it is
    /// what makes the port bit-exact; see [`mod@crate::vessels`].
    ///
    /// It cannot live inside [`Plant::derivatives`], which is `&self` on
    /// purpose. The integrator calls it once per step, alongside the other two
    /// impure phases.
    ///
    /// # Errors
    ///
    /// As [`Plant::derivatives`].
    pub fn step_seeds(&mut self, y: &State) -> Result<(), PlantError> {
        self.seeds = unpack(y, self.seeds)?.seeds;
        Ok(())
    }

    /// The warm-start temperatures this plant will use next.
    #[must_use]
    pub const fn seeds(&self) -> TemperatureSeeds {
        self.seeds
    }

    /// Set the warm-start temperatures, for a harness that needs to place the
    /// plant in a known condition.
    pub const fn set_seeds(&mut self, seeds: TemperatureSeeds) {
        self.seeds = seeds;
    }

    /// Set the latched valve commands directly, for the same reason.
    pub const fn set_valve_command(&mut self, command: [f64; 12]) {
        self.valve_command = command;
    }

    /// The walk-driven inputs the pure phase reads. Phase 3 fills these in.
    #[must_use]
    pub const fn walk_inputs(&self) -> &WalkInputs {
        &self.walks
    }

    /// Set the walk-driven inputs.
    pub const fn set_walk_inputs(&mut self, walks: WalkInputs) {
        self.walks = walks;
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
        // Once per step, after the derivative and never inside it: the warm
        // start is persistent state, and folding it into the pure phase would
        // make an RK4 driver advance it four times.
        self.step_seeds(y)?;
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

    /// The plant's nominal operating point, from `TEINIT`.
    ///
    /// Hand-built fixtures are not usable here. Several were tried and every
    /// one tripped the shutdown detector, most often on *level low*: the
    /// vessels need tens of lbmol of liquid before `VLR/35.3145` clears 2
    /// cubic metres, and a plausible-looking small state is a plant that has
    /// already emptied itself. A tripped plant returns fifty zeros, which is
    /// exactly the answer these tests exist to distinguish from a bug.
    fn running_state() -> State {
        State::from_flat(&crate::constants::NOMINAL_STATE)
    }

    /// The nominal state really is a healthy plant. Everything below depends
    /// on it, and a fixture that quietly tripped would make every one of these
    /// tests vacuous.
    #[test]
    fn the_nominal_state_does_not_trip() {
        let mut plant = Plant::new();
        let u = Inputs {
            manipulated: [50.0; 12],
            disturbances: [0.0; 20],
        };
        plant.advance_discrete(SimTime(0.0), &u);
        let (_, signals) = plant
            .derivatives(SimTime(0.1), &running_state(), &u)
            .expect("converges");
        assert!(
            !signals.shutdown.is_tripped(),
            "the nominal operating point trips: {:?}",
            signals.shutdown.first()
        );
    }

    /// The property the whole three-phase split exists to provide, now that
    /// there is physics behind it to test.
    #[test]
    fn the_derivative_is_pure() {
        let mut plant = Plant::new();
        let u = Inputs {
            manipulated: [50.0; 12],
            disturbances: [0.0; 20],
        };
        plant.advance_discrete(SimTime(0.0), &u);
        assert_derivatives_are_pure(&plant, SimTime(0.1), &running_state(), &u);
    }

    /// The plant produces a derivative that actually moves. An all-zero answer
    /// would satisfy every structural test and describe a frozen plant.
    #[test]
    fn the_derivative_is_not_all_zeros_on_a_running_plant() {
        let mut plant = Plant::new();
        let u = Inputs {
            manipulated: [50.0; 12],
            disturbances: [0.0; 20],
        };
        plant.advance_discrete(SimTime(0.0), &u);
        let (derivative, signals) = plant
            .derivatives(SimTime(0.1), &running_state(), &u)
            .expect("converges");
        let moving = derivative
            .to_flat()
            .iter()
            .filter(|v| v.to_bits() != 0.0_f64.to_bits())
            .count();
        assert!(moving > 30, "only {moving} of 50 derivatives are non-zero");
        assert_eq!(signals.continuous.len(), N_CONTINUOUS);
    }

    /// A euler step advances the state and the warm-start seeds together.
    #[test]
    fn a_euler_step_advances_both_the_state_and_the_seeds() {
        let mut plant = Plant::new();
        let u = Inputs {
            manipulated: [50.0; 12],
            disturbances: [0.0; 20],
        };
        let y = running_state();
        let before = plant.seeds();
        let (next, _) = plant
            .euler_step(SimTime(0.0), &y, &u, 1.0 / 3600.0)
            .expect("converges");
        assert!(
            next.to_flat()
                .iter()
                .zip(y.to_flat())
                .any(|(a, b)| a.to_bits() != b.to_bits()),
            "the state did not move"
        );
        assert!(
            plant.seeds() != before,
            "the warm-start seeds did not advance, so the next evaluation \
             would start from a stale guess and the port would not be \
             bit-exact"
        );
    }

    /// The latch tracks the command exactly when nothing is sticking, and
    /// clamps to the valve's travel.
    #[test]
    fn the_valve_latch_tracks_the_command_and_clamps_to_travel() {
        let mut plant = Plant::new();
        let mut u = Inputs {
            manipulated: [42.0; 12],
            ..Inputs::default()
        };
        plant.advance_discrete(SimTime(1.0), &u);
        // Exact: the latch is a copy, not an arithmetic result.
        for command in plant.valve_command() {
            assert_eq!(command.to_bits(), 42.0_f64.to_bits());
        }

        u.manipulated[0] = -5.0;
        u.manipulated[1] = 150.0;
        plant.advance_discrete(SimTime(1.0), &u);
        assert!((plant.valve_command()[0] - 0.0).abs() < f64::EPSILON);
        assert!((plant.valve_command()[1] - 100.0).abs() < f64::EPSILON);
    }

    /// `IDV(19)` makes four valves stick: the command has to move more than
    /// the stiction band before the valve follows it.
    #[test]
    fn a_sticking_valve_ignores_a_small_command_change() {
        let mut plant = Plant::new();
        let mut u = Inputs {
            manipulated: [50.0; 12],
            disturbances: [0.0; 20],
        };
        // Latch at 50 first, at a non-zero time so the `TIME = 0` branch does
        // not force it.
        plant.advance_discrete(SimTime(1.0), &u);
        u.disturbances[18] = 1.0; // IDV(19)

        // Valve 5 sticks; valve 1 does not.
        u.manipulated[4] = 51.0;
        u.manipulated[0] = 51.0;
        plant.advance_discrete(SimTime(1.0), &u);
        assert!(
            (plant.valve_command()[4] - 50.0).abs() < f64::EPSILON,
            "valve 5 should have stuck at 50"
        );
        assert!(
            (plant.valve_command()[0] - 51.0).abs() < f64::EPSILON,
            "valve 1 does not stick under IDV(19)"
        );

        // Past the band, it moves.
        u.manipulated[4] = 55.0;
        plant.advance_discrete(SimTime(1.0), &u);
        assert!((plant.valve_command()[4] - 55.0).abs() < f64::EPSILON);
    }

    /// The disturbance clamp is a threshold, not a rounding. `teprob.f:341`
    /// tests `.GT.0`, so any positive magnitude becomes exactly one.
    #[test]
    fn the_disturbance_clamp_is_a_threshold_not_a_rounding() {
        let mut u = Inputs::default();
        // A magnitude the eventual Class C fix would honour, which the
        // faithful path must still turn into a full-strength fault.
        u.disturbances[0] = 0.4;
        u.disturbances[1] = 1.0;
        u.disturbances[2] = -1.0;
        let clamped = u.clamped_disturbances();
        assert!(
            (clamped[0] - 1.0).abs() < f64::EPSILON,
            "0.4 becomes 1, not 0"
        );
        assert!((clamped[1] - 1.0).abs() < f64::EPSILON);
        assert!((clamped[2] - 0.0).abs() < f64::EPSILON);
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
