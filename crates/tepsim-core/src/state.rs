//! The integrable state, its derivative, and the vector-space operations an
//! integrator needs.
//!
//! # What the fifty slots are
//!
//! The original keeps the state in a bare `YY(50)` and unpacks it by index
//! arithmetic at `teprob.f:417-440`. That mapping is reproduced here once, in
//! [`State::from_flat`] and [`State::to_flat`], and pinned by a test that reads
//! the corresponding `COMMON/TEPROC/` variables back out of the Fortran rather
//! than trusting this comment.
//!
//! | `YY` (1-based) | Here | `teprob.f` |
//! |---|---|---|
//! | 1-8 | `reactor.moles` | `UCVR(1..3)`, `UCLR(4..8)` (417-425) |
//! | 9 | `reactor.energy` | `ETR` (431) |
//! | 10-17 | `separator.moles` | `UCVS(1..3)`, `UCLS(4..8)` (417-425) |
//! | 18 | `separator.energy` | `ETS` (432) |
//! | 19-26 | `stripper.moles` | `UCLC(1..8)` (427) |
//! | 27 | `stripper.energy` | `ETC` (433) |
//! | 28-35 | `mixing.moles` | `UCVV(1..8)` (428) |
//! | 36 | `mixing.energy` | `ETV` (434) |
//! | 37 | `reactor_cw_out_c` | `TWR` (435) |
//! | 38 | `condenser_cw_out_c` | `TWS` (436) |
//! | 39-50 | `valve_pos` | `VPOS(1..12)` (438-440) |
//!
//! # The eight slots do not mean the same thing in every vessel
//!
//! This is the part that is easy to get wrong. For the reactor and the
//! separator the original splits the array by phase: slots 1 to 3 are the
//! **vapour** holdups of A, B and C, and slots 4 to 8 are the **liquid**
//! holdups of D through H. `UCLR(1..3)` is set to zero (`teprob.f:420-421`)
//! because the non-condensibles never form a liquid, and `UCVR(4..8)` does not
//! come from the state at all: it is derived from the vapour-liquid
//! equilibrium later in the same call (`teprob.f:500-501`).
//!
//! For the stripper the eight slots are all liquid, and for the mixing zone
//! they are all vapour. The names here say `moles` rather than picking a phase,
//! because no single phase word is true of all four vessels.
//!
//! # Why there is a flat view at all
//!
//! Only to talk to the oracle and to the recorded traces, both of which speak
//! `[f64; 50]`. Nothing in the model should use it: index arithmetic over a
//! bare array is exactly the failure mode the typed state exists to remove.

use core::ops::{Index, IndexMut};

use crate::component::{ByComponent, Component};

/// How many slots the flat view has. The original integrates 50
/// (`teprob.f:24-26`).
pub const N_STATES: usize = 50;

/// One vessel's inventory: eight component holdups and an internal energy.
///
/// What the eight mean depends on the vessel; see the module documentation.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Holdup {
    /// Component holdups, in lbmol.
    pub moles: ByComponent<f64>,
    /// Internal energy of the vessel's contents.
    pub energy: f64,
}

impl Holdup {
    /// Total moles, summed in Fortran order.
    ///
    /// The order is load-bearing: `teprob.f:443-448` accumulates in this order
    /// and reassociating would change the last bits.
    #[must_use]
    pub fn total(&self) -> f64 {
        let mut total = 0.0;
        for value in self.moles.iter() {
            total += *value;
        }
        total
    }
}

/// The complete integrable state.
///
/// Named fields, no index arithmetic. See the module documentation for how this
/// maps onto the original's `YY(50)`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct State {
    /// Reactor: vapour A/B/C in slots 1-3, liquid D-H in slots 4-8.
    pub reactor: Holdup,
    /// Separator, split by phase the same way as the reactor.
    pub separator: Holdup,
    /// Stripper. All eight slots are liquid.
    pub stripper: Holdup,
    /// Mixing zone. All eight slots are vapour.
    pub mixing: Holdup,
    /// Reactor cooling water outlet temperature, `TWR`, in degrees Celsius.
    pub reactor_cw_out_c: f64,
    /// Condenser cooling water outlet temperature, `TWS`, in degrees Celsius.
    pub condenser_cw_out_c: f64,
    /// The twelve valve positions, as percentages of full travel.
    pub valve_pos: [f64; 12],
}

/// The time derivative of a [`State`].
///
/// A separate type, so that a derivative cannot be passed where a state is
/// wanted. It wraps `State` rather than repeating the field list, because two
/// hand-maintained copies of a fifty-slot layout would eventually disagree and
/// the disagreement would be silent.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(transparent)]
pub struct Derivative(State);

impl Derivative {
    /// Wrap a per-unit-time state.
    #[must_use]
    pub const fn new(rates: State) -> Self {
        Self(rates)
    }

    /// The rates, laid out like a state.
    #[must_use]
    pub const fn rates(&self) -> &State {
        &self.0
    }

    /// The rates, mutably.
    pub const fn rates_mut(&mut self) -> &mut State {
        &mut self.0
    }

    /// Unwrap.
    #[must_use]
    pub const fn into_inner(self) -> State {
        self.0
    }
}

impl From<State> for Derivative {
    fn from(rates: State) -> Self {
        Self(rates)
    }
}

/// Elementwise operations over the fifty slots.
///
/// Enough for Euler, RK4 and an embedded Dormand-Prince, and no more. Every
/// operation is elementwise, so none of them involves a reduction and none can
/// be reassociated: the determinism invariant costs nothing here.
pub trait VectorSpace: Sized {
    /// Combine two values slot by slot.
    fn zip_with(&self, other: &Self, f: impl Fn(f64, f64) -> f64) -> Self;

    /// Apply a function to every slot.
    fn map(&self, f: impl Fn(f64) -> f64) -> Self;

    /// Multiply every slot by `k`.
    fn scale(&self, k: f64) -> Self {
        self.map(|x| x * k)
    }

    /// `self + k * other`, slot by slot.
    fn add_scaled(&self, k: f64, other: &Self) -> Self {
        // Two roundings, matching what an explicit integrator written out by
        // hand would do. Not `mul_add`: see `thermo`'s module documentation for
        // why fusing is the wrong default in this crate.
        #[allow(clippy::suboptimal_flops, reason = "rounding is part of the contract")]
        self.zip_with(other, |a, b| a + k * b)
    }
}

impl VectorSpace for State {
    fn zip_with(&self, other: &Self, f: impl Fn(f64, f64) -> f64) -> Self {
        let mut out = *self;
        {
            let (mine, theirs) = (out.holdups_mut(), other.holdups());
            for (vessel, source) in mine.into_iter().zip(theirs) {
                for component in Component::ALL {
                    vessel.moles[component] = f(vessel.moles[component], source.moles[component]);
                }
                vessel.energy = f(vessel.energy, source.energy);
            }
        }
        out.reactor_cw_out_c = f(out.reactor_cw_out_c, other.reactor_cw_out_c);
        out.condenser_cw_out_c = f(out.condenser_cw_out_c, other.condenser_cw_out_c);
        for (mine, theirs) in out.valve_pos.iter_mut().zip(other.valve_pos) {
            *mine = f(*mine, theirs);
        }
        out
    }

    fn map(&self, f: impl Fn(f64) -> f64) -> Self {
        self.zip_with(self, |a, _| f(a))
    }
}

impl VectorSpace for Derivative {
    fn zip_with(&self, other: &Self, f: impl Fn(f64, f64) -> f64) -> Self {
        Self(self.0.zip_with(&other.0, f))
    }

    fn map(&self, f: impl Fn(f64) -> f64) -> Self {
        Self(self.0.map(f))
    }
}

impl State {
    /// The four vessels, in the order the flat layout stores them.
    const fn holdups(&self) -> [&Holdup; 4] {
        [&self.reactor, &self.separator, &self.stripper, &self.mixing]
    }

    /// The four vessels, mutably.
    const fn holdups_mut(&mut self) -> [&mut Holdup; 4] {
        [
            &mut self.reactor,
            &mut self.separator,
            &mut self.stripper,
            &mut self.mixing,
        ]
    }

    /// Advance by `dt` using `derivative`: the explicit Euler step.
    ///
    /// Written here rather than in an integrator so that the one place a state
    /// and a derivative are combined is next to the types themselves.
    #[must_use]
    pub fn step(&self, dt: f64, derivative: &Derivative) -> Self {
        self.add_scaled(dt, derivative.rates())
    }

    /// Read the flat `YY(50)` layout the oracle and the traces use.
    ///
    /// See the module documentation for the mapping, and
    /// `crates/tepsim-oracle/tests/state_layout.rs` for the test that pins it
    /// against the Fortran's own `COMMON` variables.
    // @port teprob.f:417-440
    #[must_use]
    pub fn from_flat(flat: &[f64; N_STATES]) -> Self {
        let mut state = Self::default();
        for (vessel, base) in state.holdups_mut().into_iter().zip([0, 9, 18, 27]) {
            for component in Component::ALL {
                vessel.moles[component] = flat[base + component.index()];
            }
            vessel.energy = flat[base + Component::COUNT];
        }
        state.reactor_cw_out_c = flat[36];
        state.condenser_cw_out_c = flat[37];
        state.valve_pos.copy_from_slice(&flat[38..N_STATES]);
        state
    }

    /// Write the flat `YY(50)` layout.
    #[must_use]
    pub fn to_flat(&self) -> [f64; N_STATES] {
        let mut flat = [0.0; N_STATES];
        for (vessel, base) in self.holdups().into_iter().zip([0, 9, 18, 27]) {
            for component in Component::ALL {
                flat[base + component.index()] = vessel.moles[component];
            }
            flat[base + Component::COUNT] = vessel.energy;
        }
        flat[36] = self.reactor_cw_out_c;
        flat[37] = self.condenser_cw_out_c;
        flat[38..N_STATES].copy_from_slice(&self.valve_pos);
        flat
    }
}

impl Derivative {
    /// Read the flat `YP(50)` layout.
    #[must_use]
    pub fn from_flat(flat: &[f64; N_STATES]) -> Self {
        Self(State::from_flat(flat))
    }

    /// Write the flat `YP(50)` layout.
    #[must_use]
    pub fn to_flat(&self) -> [f64; N_STATES] {
        self.0.to_flat()
    }
}

/// Which vessel, for code that has to name one.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Vessel {
    /// The reactor.
    Reactor,
    /// The condenser and separator.
    Separator,
    /// The stripper.
    Stripper,
    /// The mixing zone upstream of the reactor.
    Mixing,
}

impl Vessel {
    /// All four, in flat-layout order.
    pub const ALL: [Self; 4] = [Self::Reactor, Self::Separator, Self::Stripper, Self::Mixing];
}

impl Index<Vessel> for State {
    type Output = Holdup;

    fn index(&self, vessel: Vessel) -> &Holdup {
        match vessel {
            Vessel::Reactor => &self.reactor,
            Vessel::Separator => &self.separator,
            Vessel::Stripper => &self.stripper,
            Vessel::Mixing => &self.mixing,
        }
    }
}

impl IndexMut<Vessel> for State {
    fn index_mut(&mut self, vessel: Vessel) -> &mut Holdup {
        match vessel {
            Vessel::Reactor => &mut self.reactor,
            Vessel::Separator => &mut self.separator,
            Vessel::Stripper => &mut self.stripper,
            Vessel::Mixing => &mut self.mixing,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::assert_exact;

    /// A state whose every slot is distinguishable, so a mapping that swapped
    /// two of them could not pass.
    fn distinguishable() -> [f64; N_STATES] {
        let mut flat = [0.0; N_STATES];
        for (i, slot) in flat.iter_mut().enumerate() {
            *slot = (i as f64 + 1.0) * 1.5;
        }
        flat
    }

    #[test]
    fn the_flat_view_round_trips_every_slot() {
        let flat = distinguishable();
        let recovered = State::from_flat(&flat).to_flat();
        for (i, (a, b)) in flat.iter().zip(recovered.iter()).enumerate() {
            assert_exact(*b, *a, &alloc::format!("slot {i}"));
        }
    }

    /// The mapping must be a bijection: fifty distinct inputs, fifty distinct
    /// outputs. A layout that dropped a slot and duplicated another would still
    /// round trip if the duplicate happened to be read back from the same
    /// place, so count distinct values as well.
    #[test]
    fn the_flat_view_reaches_all_fifty_slots() {
        let state = State::from_flat(&distinguishable());
        let flat = state.to_flat();
        let mut seen: alloc::vec::Vec<u64> = flat.iter().map(|v| v.to_bits()).collect();
        seen.sort_unstable();
        seen.dedup();
        assert_eq!(seen.len(), N_STATES, "a slot was dropped or duplicated");
    }

    /// The index arithmetic in `teprob.f:417-440`, spelled out, so that the
    /// table in the module docs is checked and not merely written down.
    #[test]
    fn the_named_fields_sit_where_the_fortran_puts_them() {
        let flat = distinguishable();
        let state = State::from_flat(&flat);

        assert_exact(
            state.reactor.moles[Component::A],
            flat[0],
            "YY(1) = UCVR(1)",
        );
        assert_exact(
            state.reactor.moles[Component::H],
            flat[7],
            "YY(8) = UCLR(8)",
        );
        assert_exact(state.reactor.energy, flat[8], "YY(9) = ETR");
        assert_exact(
            state.separator.moles[Component::A],
            flat[9],
            "YY(10) = UCVS(1)",
        );
        assert_exact(state.separator.energy, flat[17], "YY(18) = ETS");
        assert_exact(
            state.stripper.moles[Component::A],
            flat[18],
            "YY(19) = UCLC(1)",
        );
        assert_exact(state.stripper.energy, flat[26], "YY(27) = ETC");
        assert_exact(
            state.mixing.moles[Component::A],
            flat[27],
            "YY(28) = UCVV(1)",
        );
        assert_exact(state.mixing.energy, flat[35], "YY(36) = ETV");
        assert_exact(state.reactor_cw_out_c, flat[36], "YY(37) = TWR");
        assert_exact(state.condenser_cw_out_c, flat[37], "YY(38) = TWS");
        assert_exact(state.valve_pos[0], flat[38], "YY(39) = VPOS(1)");
        assert_exact(state.valve_pos[11], flat[49], "YY(50) = VPOS(12)");
    }

    #[test]
    fn indexing_by_vessel_reaches_the_same_field() {
        let mut state = State::default();
        state[Vessel::Stripper].energy = 42.0;
        assert_exact(state.stripper.energy, 42.0, "stripper energy");
        for vessel in Vessel::ALL {
            state[vessel].moles[Component::G] = 1.0;
        }
        assert_exact(state.mixing.moles[Component::G], 1.0, "mixing G");
    }

    #[test]
    fn the_vector_space_operations_act_on_all_fifty_slots() {
        let state = State::from_flat(&distinguishable());
        let doubled = state.scale(2.0);
        for (i, (a, b)) in state
            .to_flat()
            .iter()
            .zip(doubled.to_flat().iter())
            .enumerate()
        {
            assert_exact(*b, a * 2.0, &alloc::format!("slot {i} doubled"));
        }

        let summed = state.zip_with(&state, |a, b| a + b);
        for (i, (a, b)) in state
            .to_flat()
            .iter()
            .zip(summed.to_flat().iter())
            .enumerate()
        {
            assert_exact(*b, a + a, &alloc::format!("slot {i} summed"));
        }
    }

    /// An Euler step is `y + dt * dy`, in every slot, with no slot forgotten.
    #[test]
    fn an_euler_step_advances_every_slot() {
        let state = State::from_flat(&distinguishable());
        let derivative = Derivative::new(State::from_flat(&[1.0; N_STATES]));
        let dt = 1.0 / 3600.0;
        let stepped = state.step(dt, &derivative);
        for (i, (before, after)) in state
            .to_flat()
            .iter()
            .zip(stepped.to_flat().iter())
            .enumerate()
        {
            assert_exact(*after, before + dt * 1.0, &alloc::format!("slot {i}"));
        }
    }

    #[test]
    fn a_holdup_totals_its_components_in_fortran_order() {
        let holdup = Holdup {
            moles: ByComponent::new([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]),
            energy: 0.0,
        };
        assert_exact(holdup.total(), 36.0, "total moles");
    }

    #[test]
    fn a_derivative_is_not_a_state_but_shares_its_layout() {
        let flat = distinguishable();
        let derivative = Derivative::from_flat(&flat);
        for (i, (a, b)) in flat.iter().zip(derivative.to_flat().iter()).enumerate() {
            assert_exact(*b, *a, &alloc::format!("derivative slot {i}"));
        }
        assert_exact(
            derivative.rates().reactor.energy,
            flat[8],
            "derivative slot 9",
        );
    }
}
