//! Integrators, and the reason the port could have them at all.
//!
//! # Why this is not trivial
//!
//! `TEFUNC` is not a pure function of `(t, y)`. Each call advances the
//! disturbance walks, draws from the generator for measurement noise, ticks the
//! sampled analysers and latches the valve positions through the sticking
//! logic. That is harmless for fixed-step Euler, which evaluates the right-hand
//! side exactly once per step, and it makes every multi-stage method *wrong*:
//! RK4 would advance the walks four times per step and draw four sets of noise,
//! so the disturbances would run at four times the intended rate and the
//! stochastic call order would no longer match the original's at all.
//!
//! This is why the port splits the right-hand side into three phases
//! (B-0012, and the decision entry of 2026-08-11):
//!
//! - **Before**, impure: `Plant::advance_discrete`, once per outer step.
//! - **Pure**: `Plant::derivatives`, as many times as the method needs.
//! - **After**, impure: `Plant::sample_measurements`, once per outer step.
//!
//! An integrator therefore calls `derivatives` `s` times for an `s`-stage
//! method and the discrete work still happens exactly once. That property is
//! the whole point of the split and it is asserted in
//! `the_walks_advance_once_per_step_whichever_method`.
//!
//! # Which to use
//!
//! [`Integrator::Euler`] is what the original does, and it is the only one that
//! reproduces the Fortran bit for bit. Everything the validation ladder claims
//! is a claim about Euler. The others exist because a one-second explicit Euler
//! step is a poor way to integrate a stiff-ish plant, and a user who wants
//! accuracy rather than fidelity should be able to ask for it.
//!
//! Choosing anything other than Euler is therefore a deliberate departure from
//! the original, and [`Integrator::is_faithful`] says so in code.

use tepsim_core::{Derivative, Inputs, Plant, PlantError, Signals, SimTime, State, VectorSpace};

/// How to advance the state.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Integrator {
    /// Fixed-step explicit Euler: `y + h f(t, y)`.
    ///
    /// What `temain_mod.f`'s `INTGTR` does, and the default. First order, one
    /// evaluation per step, and the only choice under which this port is
    /// bit-identical to the Fortran.
    #[default]
    Euler,
    /// Classical fixed-step fourth-order Runge-Kutta.
    ///
    /// Four evaluations per step. Fourth order in the step size, so halving the
    /// step divides the local error by sixteen rather than by two.
    Rk4,
    /// Fixed-step Dormand-Prince 5(4), using only its fifth-order solution.
    ///
    /// Seven stages, of which the last is the next step's first (the FSAL
    /// property), so six evaluations per step after the first. The embedded
    /// fourth-order solution is computed as well and reported through
    /// [`Step::error_estimate`], which is what an adaptive driver would need;
    /// the step size is not yet adapted, because a variable step changes when
    /// the discrete phases run and that is a decision about fidelity rather
    /// than about numerics.
    DormandPrince,
}

impl Integrator {
    /// Whether this method reproduces the original.
    ///
    /// Only Euler does. The validation ladder's every claim is a claim about
    /// Euler, so a run using anything else is a *better* integration of the
    /// same equations and not a reproduction of the same numbers.
    #[must_use]
    pub const fn is_faithful(self) -> bool {
        matches!(self, Self::Euler)
    }

    /// How many derivative evaluations one step costs.
    #[must_use]
    pub const fn stages(self) -> usize {
        match self {
            Self::Euler => 1,
            Self::Rk4 => 4,
            Self::DormandPrince => 7,
        }
    }

    /// The order of the method's local error.
    #[must_use]
    pub const fn order(self) -> u32 {
        match self {
            Self::Euler => 1,
            Self::Rk4 => 4,
            Self::DormandPrince => 5,
        }
    }

    /// A short name, for a report or a CLI flag.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Euler => "euler",
            Self::Rk4 => "rk4",
            Self::DormandPrince => "dopri5",
        }
    }

    /// Parse a name, for the CLI.
    #[must_use]
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "euler" => Some(Self::Euler),
            "rk4" => Some(Self::Rk4),
            "dopri5" | "dormand-prince" => Some(Self::DormandPrince),
            _ => None,
        }
    }
}

/// What the first stage produced, passed on to a method that needs it.
///
/// A struct rather than three arguments: `dormand_prince` already takes the
/// closure, the time, the state and the step, and four more would put it past
/// the point where a reader can hold the signature in their head.
struct FirstStage {
    /// `f(t, y)`, the first stage's slope.
    k1: State,
    /// The same, as the type the caller wants back.
    first: Derivative,
    /// The noise-free signals at the start of the step.
    signals: Signals,
}

/// The result of one integrator step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Step {
    /// The state at the end of the step.
    pub state: State,
    /// The derivative at the *start* of the step.
    pub first: Derivative,
    /// The noise-free signals at the start of the step.
    ///
    /// Returned rather than left to the caller because the measurement vector
    /// is a function of them and the first stage has already computed them.
    /// Making the caller re-evaluate would double Euler's cost for nothing.
    pub signals: Signals,
    /// An estimate of the local truncation error, where the method provides
    /// one.
    ///
    /// `None` for Euler and RK4, which have no embedded lower-order solution.
    /// For Dormand-Prince it is the difference between the fifth and fourth
    /// order solutions, in the state's own units.
    pub error_estimate: Option<State>,
}

impl Integrator {
    /// Advance one step.
    ///
    /// **Calls only the pure phase.** The caller is responsible for running
    /// `Plant::advance_discrete` once before this and
    /// `Plant::sample_measurements` once after, which is what keeps the
    /// disturbance walks and the generator advancing once per outer step
    /// however many stages the method has.
    ///
    /// # Errors
    ///
    /// Whatever `Plant::derivatives` returns: a temperature solve that failed
    /// to converge. A multi-stage method evaluates at intermediate states that
    /// are not physical states of the plant, so it can fail where Euler would
    /// not, and the caller is told rather than handed a guess.
    pub fn advance(
        self,
        plant: &Plant,
        t: SimTime,
        y: &State,
        u: &Inputs,
        h: f64,
    ) -> Result<Step, PlantError> {
        let f = |time: f64, state: &State| -> Result<State, PlantError> {
            plant
                .derivatives(SimTime(time), state, u)
                .map(|(derivative, _)| *derivative.rates())
        };
        let t0 = t.hours();

        // The first stage is evaluated through the full `derivatives` call
        // rather than the closure, because the caller wants its signals.
        let (first, signals) = plant.derivatives(t, y, u)?;
        let k1 = *first.rates();

        let state = match self {
            Self::Euler => y.add_scaled(h, &k1),
            Self::Rk4 => {
                let k2 = f(t0 + 0.5 * h, &y.add_scaled(0.5 * h, &k1))?;
                let k3 = f(t0 + 0.5 * h, &y.add_scaled(0.5 * h, &k2))?;
                let k4 = f(t0 + h, &y.add_scaled(h, &k3))?;
                // y + h/6 (k1 + 2 k2 + 2 k3 + k4), accumulated in stage order.
                // Written as four `add_scaled` calls rather than as one fused
                // expression so the summation order is explicit and fixed;
                // determinism is a hard invariant here.
                let sixth = h / 6.0;
                y.add_scaled(sixth, &k1)
                    .add_scaled(2.0 * sixth, &k2)
                    .add_scaled(2.0 * sixth, &k3)
                    .add_scaled(sixth, &k4)
            }
            Self::DormandPrince => {
                return self.dormand_prince(&f, t0, y, h, FirstStage { k1, first, signals });
            }
        };

        Ok(Step {
            state,
            first,
            signals,
            error_estimate: None,
        })
    }

    /// Dormand-Prince 5(4), kept separate because its tableau is long enough
    /// to drown the other two methods if inlined above.
    ///
    /// Takes the first stage's results as one value rather than as three
    /// arguments; splitting them out pushed the signature past the point where
    /// a reader can hold it.
    ///
    /// Coefficients from Dormand, J. R. and Prince, P. J. (1980), "A family of
    /// embedded Runge-Kutta formulae", *Journal of Computational and Applied
    /// Mathematics* 6(1), 19-26. Transcribed as exact rational quotients
    /// rather than as decimals, because several of them are not representable
    /// and a rounded tableau is a different method with a different order.
    fn dormand_prince(
        self,
        f: &impl Fn(f64, &State) -> Result<State, PlantError>,
        t0: f64,
        y: &State,
        h: f64,
        stage_one: FirstStage,
    ) -> Result<Step, PlantError> {
        let FirstStage { k1, first, signals } = stage_one;
        let k2 = f(t0 + h / 5.0, &y.add_scaled(h / 5.0, &k1))?;
        let k3 = f(
            t0 + 3.0 * h / 10.0,
            &y.add_scaled(3.0 * h / 40.0, &k1)
                .add_scaled(9.0 * h / 40.0, &k2),
        )?;
        let k4 = f(
            t0 + 4.0 * h / 5.0,
            &y.add_scaled(44.0 * h / 45.0, &k1)
                .add_scaled(-56.0 * h / 15.0, &k2)
                .add_scaled(32.0 * h / 9.0, &k3),
        )?;
        let k5 = f(
            t0 + 8.0 * h / 9.0,
            &y.add_scaled(19372.0 * h / 6561.0, &k1)
                .add_scaled(-25360.0 * h / 2187.0, &k2)
                .add_scaled(64448.0 * h / 6561.0, &k3)
                .add_scaled(-212.0 * h / 729.0, &k4),
        )?;
        let k6 = f(
            t0 + h,
            &y.add_scaled(9017.0 * h / 3168.0, &k1)
                .add_scaled(-355.0 * h / 33.0, &k2)
                .add_scaled(46732.0 * h / 5247.0, &k3)
                .add_scaled(49.0 * h / 176.0, &k4)
                .add_scaled(-5103.0 * h / 18656.0, &k5),
        )?;

        // The fifth-order solution. Note there is no k2 term: b2 is zero.
        let fifth = y
            .add_scaled(35.0 * h / 384.0, &k1)
            .add_scaled(500.0 * h / 1113.0, &k3)
            .add_scaled(125.0 * h / 192.0, &k4)
            .add_scaled(-2187.0 * h / 6784.0, &k5)
            .add_scaled(11.0 * h / 84.0, &k6);

        // The seventh stage, evaluated at the new point. It is the next step's
        // k1 (the FSAL property), which an adaptive driver would carry over;
        // here it is needed only for the embedded solution.
        let k7 = f(t0 + h, &fifth)?;

        let fourth = y
            .add_scaled(5179.0 * h / 57600.0, &k1)
            .add_scaled(7571.0 * h / 16695.0, &k3)
            .add_scaled(393.0 * h / 640.0, &k4)
            .add_scaled(-92097.0 * h / 339_200.0, &k5)
            .add_scaled(187.0 * h / 2100.0, &k6)
            .add_scaled(h / 40.0, &k7);

        Ok(Step {
            state: fifth,
            first,
            signals,
            error_estimate: Some(fifth.zip_with(&fourth, |a, b| a - b)),
        })
    }
}
