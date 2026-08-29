//! Tier 8: unguided differential search over the derivative.
//!
//! Tier 2 compares the fifty derivatives over three pools that somebody chose:
//! the nominal trajectory, perturbations of it, and a hand-built catalogue of
//! branch boundaries. Every one of those states exists because a human thought
//! of it. `PLAN.org` says of this tier that it "is how we find the branch we
//! did not think to test", and that sentence is the whole specification: Tier 8
//! samples the input space without a hypothesis about where the interesting
//! part is.
//!
//! # What a tuple is
//!
//! [`Tuple`] is the five things `PLAN.org` names: a state, the twelve
//! manipulated variables, the twenty `IDV` flags, the generator word, and the
//! time. Everything else a [`Scenario`] needs (the walk state, the four Newton
//! seed temperatures, the latched valve commands, the previous measurements)
//! is taken from the nominal starting condition, so that a tuple is a small,
//! writable, comparable object and two implementations forced into it differ
//! only in the model.
//!
//! # Deterministic from a seed, and reproducible one at a time
//!
//! A fuzz finding is worth nothing if it cannot be replayed. [`Generator`] is
//! therefore not a stream: `Generator::tuple(index)` is a pure function of the
//! seed and the index, built by re-seeding a [`Sampler`] from a mix of the two.
//! Tuple 91,432 of seed `0x7E2_0062` can be regenerated on its own, in a
//! debugger, without drawing the 91,431 before it, and the report says exactly
//! that pair.
//!
//! No new dependency and no new generator: the [`Sampler`] Tier 1 already uses
//! is SplitMix64, is already proved to be full-period, and is already the thing
//! every other pool in this repository draws from.
//!
//! # Sampling on a logarithmic scale, and not uniformly
//!
//! The fifty slots are not one population. The stripper's B holdup sits at
//! 8.0e-3 lbmol and the reactor's G holdup at 154; the four energies are order
//! one; the two cooling-water temperatures are degrees Celsius near 90, and the
//! twelve valve positions are percentages. A uniform draw over any single
//! interval would be a rounding error for one group and nonsense for another,
//! so [`Shape`] scales the extensive slots multiplicatively by `10^u` and
//! treats the temperatures and the percentages on their own terms.
//!
//! The four shapes exist for four different reasons, and [`Generator::tuple`]
//! mixes them deliberately rather than picking one:
//!
//! - [`Shape::Jitter`] stays within a tenth of a percent of nominal, where the
//!   plant is well conditioned and the comparison is at its sharpest.
//! - [`Shape::LogScale`] moves every slot over four orders of magnitude, which
//!   is the range the state itself spans.
//! - [`Shape::Sparse`] moves one to five slots and leaves the rest exactly
//!   nominal. A counterexample from this shape arrives nearly minimal already,
//!   and a mistake confined to one branch shows up here first.
//! - [`Shape::Wild`] goes to twelve orders of magnitude with zeros and sign
//!   flips. These states are not physical and mostly do not converge; they are
//!   run so that the harness's own behaviour on an unphysical state is under
//!   test rather than assumed.
//!
//! # A random state is not a physical state, and that is delta D-001
//!
//! `TESUB2` (`teprob.f:1415-1442`) runs Newton for a hundred iterations and, on
//! failure, executes `T=TIN` and returns *exactly as it does on success*. There
//! is no status flag. The port returns `Err(PlantError::Temperature)` instead,
//! which is delta D-001, and on a state drawn from twelve orders of magnitude
//! that path fires often.
//!
//! So a `PlantError` here is not a counterexample. It is counted as
//! [`Outcome::PortDidNotConverge`] and reported as its own line, because the
//! fraction of tuples that even reach a comparison is a number the next session
//! needs.
//!
//! The reverse case is caught too. If the Fortran fell through and the port
//! converged, the two would disagree wildly for a reason that is D-001 and not
//! a porting mistake. `TESUB2`'s fall-through writes back the caller's guess
//! *bit for bit*, so [`fortran_fell_through`] detects it by comparing the four
//! temperatures the call returned against the four it was given. That
//! misclassifies a solve that converged to exactly its own guess, which cannot
//! happen on a state drawn at random and is documented rather than guarded.
//!
//! # The gate is 1e-12 of the scale of the terms
//!
//! Not of the derivative. A balance is inflow minus outflow and those nearly
//! cancel, so the derivative's own magnitude says nothing about how accurately
//! it can be computed. See the decision of 2026-08-27 in `BACKLOG.org`,
//! `tepsim_core::balances::Balances::scale` for where the budget comes from,
//! and [`Comparison::observe_against`] for the accumulator. Tier 8 asks exactly
//! the question Tier 2 asks, of states nobody chose.
//!
//! # Non-finite results are scoped out, deliberately and narrowly
//!
//! On a state twelve orders of magnitude from anything physical, both
//! implementations can overflow. The Fortran has no error handling whatsoever,
//! so its behaviour there is not a specification and cannot be a reference.
//! [`slot_verdict`] therefore treats two `NaN`s as agreeing, because IEEE-754
//! does not define which payload an arithmetic `NaN` carries and gfortran and
//! LLVM need not choose the same one.
//!
//! What is *not* scoped out is disagreement about whether an answer exists at
//! all: a `NaN` against a number, or two infinities of opposite sign, is a
//! finding with an infinite error. That is the part of the non-finite domain
//! where the two really can be compared, and it is compared.
//!
//! # Shrinking is what makes a finding useful
//!
//! A tuple straight out of the generator has fifty perturbed slots, twelve
//! manipulated variables and twenty flags away from nominal, and says nothing
//! about *which* of them matters. [`shrink`] walks the eighty-four knobs
//! resetting each to its nominal value and keeping the reset whenever the
//! failure survives, repeats that to a fixpoint, then bisects each surviving
//! knob geometrically toward nominal. What comes out is a tuple that still
//! fails and differs from nominal in as few places, by as little, as the search
//! could manage.
//!
//! Shrinking preserves *a* failure, not necessarily the same one. That is the
//! standard contract, and it is the right one: the minimal tuple is a
//! reproducer, and the reproducer is what goes in the corpus.
//!
//! # Proving the harness has teeth
//!
//! "No counterexamples" is not a result unless the search could have produced
//! one. [`Mutation`] is a deliberate corruption applied to the port's answer
//! after the fact, standing in for the porting mistakes this tier exists to
//! catch: a term of the wrong size everywhere, and a term of the wrong size
//! only on one side of a branch. `tests/tier8_fuzz.rs` runs the whole search
//! against each mutant and requires it to find and shrink the fault.
//!
//! The mutation is expressed as a fraction of the balance's own scale rather
//! than of its value, for the same reason the gate is. `1e-11` of the scale is
//! ten times the gate and must be caught; `1e-13` is a tenth of it and must
//! not be. Both are asserted, so the harness's sensitivity is a measured
//! number and not a hope.
//!
//! # What it found
//!
//! Five million tuples, seed `0x7E2_0062`: 2,551,618 reached a comparison,
//! 2,263,855 froze, 184,526 lost `TESUB2`, and one disagreed. Shrunk to three
//! knobs: the nominal state, `IDV(13)` on, one generator word, `TIME = 973 h`.
//! At that point the kinetic drift has walked to 3e7 and the reaction rates to
//! 1.1e10, so the one ULP the vendored `exp` costs becomes 4.6e-12 of the
//! reactor energy balance's scale, past the 1e-12 gate. Bit-identical under
//! `libm-system`, so the algebra is exact and the amplification is the whole
//! story. No state in any Tier 2 pool has a drift factor above about one, which
//! is the case for building this tier in one sentence.
//!
//! `tests/tier8_fuzz.rs` carries the reproducer and the attribution.
//!
//! # Size
//!
//! [`Budget`] reads `TEP_TIER8`, in the idiom of `TEP_TIER5` and `TEP_TIER7`.
//! Unset is a few hundred tuples and runs in a second, which is what belongs in
//! a per-commit gate; `full` is millions and belongs in a nightly. A decimal
//! count is also accepted, for a bisection run that wants a specific size.

use core::fmt;

use tepsim_core::state::N_STATES;
use tepsim_core::{Inputs, Plant, SimTime, State, plant, vessels};

use crate::Oracle;
use crate::tier1::{Comparison, Sampler};
use crate::tier2::{Scenario, Snapshot};

/// The Tier 2 gate, which is also Tier 8's: 1e-12 of the scale of the terms.
///
/// `PLAN.org`, "Tier 2", as amended by the decision of 2026-08-27.
pub const TOLERANCE: f64 = 1e-12;

/// How many independently shrinkable knobs a [`Tuple`] has.
///
/// Fifty state slots, twelve manipulated variables, twenty disturbance flags,
/// the generator word and the time.
pub const N_KNOBS: usize = N_STATES + 12 + 20 + 2;

/// The knob index of the generator word.
const KNOB_RNG: usize = N_STATES + 12 + 20;
/// The knob index of the simulation time.
const KNOB_TIME: usize = KNOB_RNG + 1;

// ---------------------------------------------------------------------------
// The base scenario
// ---------------------------------------------------------------------------

/// The nominal starting condition every tuple is an overlay on.
///
/// # Why this is not `Pools::collect(oracle, 0, dt).nominal`
///
/// Because that one is not the same twice. `TCR`, `TCS`, `TCC` and `TCV` are
/// the Newton warm starts for the four vessel temperatures and `TEINIT` never
/// assigns them (see [`Oracle::init_cold`]); they are whatever the previous run
/// in this process left in `COMMON`. Tier 2 can live with that, because
/// [`Scenario::force`] restores the whole block and its scenarios are only ever
/// compared against themselves.
///
/// Tier 8 cannot. Two consequences, both measured rather than reasoned about:
///
/// - The census moved between tests in the same binary. The same 400 tuples of
///   the same seed reported 168 frozen and 24 fall-throughs in one test and 184
///   and 8 in another, purely because the seeds differed by a few ULP and that
///   moved sixteen states across `TESUB2`'s convergence boundary.
/// - Worse, a corpus entry would not be a fixed input. A tuple recorded today
///   would be evaluated against a different warm start tomorrow, depending on
///   which tests ran first, which is exactly the property a regression corpus
///   exists to have.
///
/// [`Oracle::init_cold`] zeroes the four first, which reproduces the
/// freshly-loaded-process result whatever the history. So the base here is a
/// function of the Fortran alone.
pub fn nominal_scenario(oracle: &mut Oracle) -> Scenario {
    let (time, state) = oracle.init_cold();
    Scenario {
        time,
        state,
        manipulated: oracle.manipulated(),
        disturbances: oracle.disturbances(),
        walk: oracle.wlk(),
        rng: oracle.rng(),
        measurements: oracle.measurements(),
        common: oracle.teproc(),
    }
}

// ---------------------------------------------------------------------------
// The tuple
// ---------------------------------------------------------------------------

/// One generated input: everything Tier 8 varies.
///
/// The five fields `PLAN.org` names. Everything else a [`Scenario`] carries is
/// inherited from the nominal starting condition by [`Tuple::scenario`], which
/// is what makes a tuple small enough to write down as a literal and put in a
/// regression corpus.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Tuple {
    /// Simulation time, in hours. Zero is special: `teprob.f:397-406` resets
    /// the whole walk state when `TIME.EQ.0.D0`, so it is drawn often enough to
    /// exercise that path and rarely enough not to dominate.
    pub time: f64,
    /// The fifty integrated states, `YY`.
    pub state: [f64; N_STATES],
    /// The twelve manipulated variables, `XMV`, as percentages.
    pub manipulated: [f64; 12],
    /// The twenty disturbance flags, `IDV`.
    ///
    /// Values outside `{0, 1}` are generated on purpose. `teprob.f:341-346`
    /// clamps anything positive to one and anything else to zero, and
    /// `Inputs::clamped_disturbances` is the port's copy of that; a fuzzer that
    /// only ever produced zeros and ones would never compare the two.
    pub disturbances: [i32; 20],
    /// The generator word, `COMMON/RANDSD/ G`.
    pub rng: f64,
}

impl Tuple {
    /// The tuple that reproduces the nominal starting condition exactly.
    ///
    /// This is what [`shrink`] moves toward, one knob at a time.
    #[must_use]
    pub fn nominal(base: &Scenario) -> Self {
        Self {
            time: base.time,
            state: base.state,
            manipulated: base.manipulated,
            disturbances: base.disturbances,
            rng: base.rng,
        }
    }

    /// Overlay this tuple on the nominal starting condition.
    ///
    /// The walk state, the four Newton seed temperatures, the latched valve
    /// commands and the previous measurements all come from `base`. Both
    /// implementations then start from the same seeds, so a difference in the
    /// answer is a difference in the model rather than in the warm start.
    #[must_use]
    pub fn scenario(&self, base: &Scenario) -> Scenario {
        Scenario {
            time: self.time,
            state: self.state,
            manipulated: self.manipulated,
            disturbances: self.disturbances,
            rng: self.rng,
            ..base.clone()
        }
    }

    /// Whether knob `k` holds the same bits as it does in `other`.
    ///
    /// Bits rather than values, so that a shrink step which changed nothing is
    /// recognised as having changed nothing even for a zero of either sign.
    #[must_use]
    pub fn knob_matches(&self, other: &Self, k: usize) -> bool {
        match Knob::of(k) {
            Knob::State(i) => self.state[i].to_bits() == other.state[i].to_bits(),
            Knob::Manipulated(i) => self.manipulated[i].to_bits() == other.manipulated[i].to_bits(),
            Knob::Disturbance(i) => self.disturbances[i] == other.disturbances[i],
            Knob::Rng => self.rng.to_bits() == other.rng.to_bits(),
            Knob::Time => self.time.to_bits() == other.time.to_bits(),
        }
    }

    /// Copy knob `k` from `other`, leaving every other knob alone.
    pub fn take_knob(&mut self, other: &Self, k: usize) {
        match Knob::of(k) {
            Knob::State(i) => self.state[i] = other.state[i],
            Knob::Manipulated(i) => self.manipulated[i] = other.manipulated[i],
            Knob::Disturbance(i) => self.disturbances[i] = other.disturbances[i],
            Knob::Rng => self.rng = other.rng,
            Knob::Time => self.time = other.time,
        }
    }

    /// Move knob `k` halfway toward `other`, geometrically where that is
    /// meaningful.
    ///
    /// A state slot that is a factor of a thousand away from nominal is not
    /// usefully halved arithmetically: thirty arithmetic halvings would still
    /// leave it a factor of one away, and the interesting structure in these
    /// slots is multiplicative. So a slot with the same sign as its nominal and
    /// a non-zero nominal is moved by the square root of its ratio. Everything
    /// else, the flags included, is moved arithmetically or reset outright.
    ///
    /// Returns whether anything moved.
    pub fn halve_knob(&mut self, other: &Self, k: usize) -> bool {
        let before = *self;
        match Knob::of(k) {
            Knob::State(i) => self.state[i] = halve_toward(self.state[i], other.state[i]),
            Knob::Manipulated(i) => {
                self.manipulated[i] = midpoint(self.manipulated[i], other.manipulated[i]);
            }
            // A flag has nothing between it and its target.
            Knob::Disturbance(i) => self.disturbances[i] = other.disturbances[i],
            Knob::Rng => self.rng = midpoint(self.rng, other.rng),
            Knob::Time => self.time = midpoint(self.time, other.time),
        }
        !before.knob_matches(self, k)
    }

    /// How many knobs differ from `nominal`, which is the size a shrink is
    /// trying to reduce.
    #[must_use]
    pub fn distance(&self, nominal: &Self) -> usize {
        (0..N_KNOBS)
            .filter(|k| !self.knob_matches(nominal, *k))
            .count()
    }

    /// The knobs that differ from `nominal`, named.
    #[must_use]
    pub fn differences(&self, nominal: &Self) -> Vec<String> {
        (0..N_KNOBS)
            .filter(|k| !self.knob_matches(nominal, *k))
            .map(|k| match Knob::of(k) {
                Knob::State(i) => {
                    format!(
                        "YY({}) {:?} <- {:?}",
                        i + 1,
                        self.state[i],
                        nominal.state[i]
                    )
                }
                Knob::Manipulated(i) => format!(
                    "XMV({}) {:?} <- {:?}",
                    i + 1,
                    self.manipulated[i],
                    nominal.manipulated[i]
                ),
                Knob::Disturbance(i) => format!(
                    "IDV({}) {} <- {}",
                    i + 1,
                    self.disturbances[i],
                    nominal.disturbances[i]
                ),
                Knob::Rng => format!("G {:?} <- {:?}", self.rng, nominal.rng),
                Knob::Time => format!("TIME {:?} <- {:?}", self.time, nominal.time),
            })
            .collect()
    }

    /// Emit this tuple as Rust source, ready to paste into the corpus.
    ///
    /// A counterexample is only permanently useful if it stops depending on the
    /// generator that produced it, so the harness prints the whole tuple as
    /// literals rather than as a seed and an index. `{:?}` on an `f64` is the
    /// shortest decimal that round-trips, so pasting this back gives the same
    /// bits.
    #[must_use]
    pub fn as_rust_literal(&self, name: &str, why: &str) -> String {
        let mut out = String::new();
        out.push_str(&format!("/// {why}\nconst {name}: Tuple = Tuple {{\n"));
        out.push_str(&format!("    time: {:?},\n    state: [\n", self.time));
        for chunk in self.state.chunks(4) {
            out.push_str("        ");
            for value in chunk {
                out.push_str(&format!("{value:?}, "));
            }
            out.push('\n');
        }
        out.push_str("    ],\n    manipulated: [\n");
        for chunk in self.manipulated.chunks(4) {
            out.push_str("        ");
            for value in chunk {
                out.push_str(&format!("{value:?}, "));
            }
            out.push('\n');
        }
        out.push_str("    ],\n    disturbances: [\n");
        for chunk in self.disturbances.chunks(10) {
            out.push_str("        ");
            for value in chunk {
                out.push_str(&format!("{value}, "));
            }
            out.push('\n');
        }
        out.push_str(&format!("    ],\n    rng: {:?},\n}};\n", self.rng));
        out
    }
}

/// One scalar degree of freedom of a [`Tuple`], as the shrinker sees it.
///
/// A flat index over the eighty-four, so the shrinker can loop rather than
/// repeat itself four times.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Knob {
    /// A slot of `YY`, zero-based.
    State(usize),
    /// A slot of `XMV`, zero-based.
    Manipulated(usize),
    /// A slot of `IDV`, zero-based.
    Disturbance(usize),
    /// The generator word.
    Rng,
    /// The simulation time.
    Time,
}

impl Knob {
    /// Decode a flat knob index.
    ///
    /// # Panics
    ///
    /// If `k` is not below [`N_KNOBS`].
    #[must_use]
    pub const fn of(k: usize) -> Self {
        if k < N_STATES {
            Self::State(k)
        } else if k < N_STATES + 12 {
            Self::Manipulated(k - N_STATES)
        } else if k < KNOB_RNG {
            Self::Disturbance(k - N_STATES - 12)
        } else if k == KNOB_RNG {
            Self::Rng
        } else if k == KNOB_TIME {
            Self::Time
        } else {
            panic!("knob index out of range")
        }
    }
}

/// Halfway from `value` to `target` on a logarithmic scale where that is
/// defined, arithmetically where it is not.
#[allow(
    clippy::suboptimal_flops,
    reason = "an explicit square root of the ratio, not a polynomial to fuse"
)]
fn halve_toward(value: f64, target: f64) -> f64 {
    if !value.is_finite() || !target.is_finite() {
        return target;
    }
    if target != 0.0 && value != 0.0 && value.signum() == target.signum() {
        let moved = target * (value / target).sqrt();
        // A ratio within one ULP of one has nowhere left to go; snap, so the
        // shrinker terminates instead of grinding on the last bit.
        if moved.to_bits() == value.to_bits() {
            return target;
        }
        return moved;
    }
    midpoint(value, target)
}

/// The arithmetic midpoint, snapping when there is nothing between the two.
fn midpoint(value: f64, target: f64) -> f64 {
    let moved = 0.5 * (value + target);
    if !moved.is_finite() || moved.to_bits() == value.to_bits() {
        target
    } else {
        moved
    }
}

// ---------------------------------------------------------------------------
// The generator
// ---------------------------------------------------------------------------

/// How one tuple's state is drawn. See the module documentation for why there
/// are four of these and not one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    /// Every slot nudged by a log-uniform relative amount between 1e-12 and
    /// 1e-3, with a random sign.
    Jitter,
    /// Every slot scaled by `10^u`, `u` uniform on `[-2, 2]`.
    LogScale,
    /// One to five slots scaled that way; the rest exactly nominal.
    Sparse,
    /// `10^u` on `[-6, 6]`, with zeros and sign flips. Mostly unphysical.
    Wild,
}

impl Shape {
    /// Pick a shape.
    ///
    /// Weighted toward [`Shape::Sparse`], because a sparse counterexample is
    /// already close to minimal and therefore the cheapest kind to act on, and
    /// away from [`Shape::Wild`], because most wild states never reach a
    /// comparison at all.
    fn draw(sampler: &mut Sampler) -> Self {
        let u = sampler.unit();
        if u < 0.40 {
            Self::Sparse
        } else if u < 0.65 {
            Self::Jitter
        } else if u < 0.90 {
            Self::LogScale
        } else {
            Self::Wild
        }
    }
}

/// What kind of quantity a state slot holds.
///
/// The fifty slots are three populations, and scaling a valve position by a
/// million is not a test of anything. `teprob.f:417-440` is where the split
/// comes from; `tepsim_core::state` documents it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum SlotKind {
    /// A component holdup in lbmol or a vessel's internal energy: slots 1-36.
    /// Positive, and spanning four orders of magnitude at nominal alone.
    Extensive,
    /// A cooling-water outlet temperature in degrees Celsius: slots 37-38.
    Celsius,
    /// A valve position as a percentage of travel: slots 39-50.
    Percent,
}

const fn slot_kind(slot: usize) -> SlotKind {
    if slot < 36 {
        SlotKind::Extensive
    } else if slot < 38 {
        SlotKind::Celsius
    } else {
        SlotKind::Percent
    }
}

/// Deterministic tuples from a seed.
///
/// Not a stream. [`Generator::tuple`] is a pure function of the seed and the
/// index, so a reported counterexample can be regenerated on its own; see the
/// module documentation.
#[derive(Clone, Debug)]
pub struct Generator {
    seed: u64,
    nominal: Tuple,
}

impl Generator {
    /// A generator seeded with `seed`, drawing around the nominal condition.
    #[must_use]
    pub fn new(seed: u64, base: &Scenario) -> Self {
        Self {
            seed,
            nominal: Tuple::nominal(base),
        }
    }

    /// The seed this generator was built with, for the report.
    #[must_use]
    pub const fn seed(&self) -> u64 {
        self.seed
    }

    /// The nominal tuple, which is what a shrink moves toward.
    #[must_use]
    pub const fn nominal(&self) -> &Tuple {
        &self.nominal
    }

    /// Tuple number `index`, independently of every other index.
    #[must_use]
    #[allow(
        clippy::suboptimal_flops,
        reason = "reproducibility across targets beats one rounding; see Sampler::range"
    )]
    pub fn tuple(&self, index: u64) -> Tuple {
        let mut s = Sampler::new(stream_seed(self.seed, index));
        let shape = Shape::draw(&mut s);
        let mut tuple = self.nominal;

        // Which slots move. Sparse picks a handful; the rest move everything,
        // so that a mistake needing two slots at once is reachable.
        let touched: Vec<usize> = if shape == Shape::Sparse {
            let count = 1 + (s.unit() * 5.0) as usize;
            (0..count)
                .map(|_| (s.unit() * N_STATES as f64) as usize % N_STATES)
                .collect()
        } else {
            (0..N_STATES).collect()
        };

        for slot in touched {
            let nominal = self.nominal.state[slot];
            tuple.state[slot] = match (shape, slot_kind(slot)) {
                (Shape::Jitter, _) => {
                    // Log-uniform between 1e-12 and 1e-3 of the slot's own
                    // magnitude, signed. Relative, because the slots do not
                    // share a scale.
                    let magnitude = 1e-12 * 1e9_f64.powf(s.unit());
                    let signed = if s.unit() < 0.5 {
                        -magnitude
                    } else {
                        magnitude
                    };
                    nominal * (1.0 + signed)
                }
                (Shape::LogScale | Shape::Sparse, SlotKind::Extensive) => {
                    nominal * 10.0_f64.powf(s.range(-2.0, 2.0))
                }
                (Shape::Wild, SlotKind::Extensive) => {
                    let u = s.unit();
                    if u < 0.05 {
                        0.0
                    } else if u < 0.10 {
                        -nominal
                    } else {
                        nominal * 10.0_f64.powf(s.range(-6.0, 6.0))
                    }
                }
                // Cooling-water outlet temperatures. Additive, over the range
                // the two heat balances can plausibly be asked about, and past
                // it under Wild.
                (Shape::LogScale | Shape::Sparse, SlotKind::Celsius) => s.range(0.0, 200.0),
                (Shape::Wild, SlotKind::Celsius) => s.range(-300.0, 800.0),
                // Valve positions. The plant does not clamp these: the
                // controller does, so out-of-range values are reachable in the
                // model and are generated.
                (Shape::LogScale | Shape::Sparse, SlotKind::Percent) => s.range(0.0, 100.0),
                (Shape::Wild, SlotKind::Percent) => s.range(-50.0, 150.0),
            };
        }

        for slot in &mut tuple.manipulated {
            *slot = if shape == Shape::Wild {
                s.range(-50.0, 150.0)
            } else {
                s.range(0.0, 100.0)
            };
        }

        // Half the tuples carry no fault at all, so that the fault-free
        // right-hand side keeps most of the sample. Of the rest, most flags are
        // binary and a few are out of range, to compare the two clamps.
        if s.unit() < 0.5 {
            for slot in &mut tuple.disturbances {
                *slot = match s.unit() {
                    u if u < 0.80 => 0,
                    u if u < 0.95 => 1,
                    u if u < 0.98 => 2,
                    _ => -1,
                };
            }
        }

        // `TESUB7` reduces the word modulo 2^32 on the first draw, so any
        // non-negative value is a legal starting point.
        tuple.rng = (s.unit() * 4_294_967_296.0).floor();

        // Zero often, because `teprob.f:397-406` resets the whole walk state
        // there and a pool that never visited it would leave that path
        // uncompared.
        tuple.time = match s.unit() {
            u if u < 0.25 => 0.0,
            u if u < 0.90 => s.range(0.0, 48.0),
            _ => s.range(0.0, 1000.0),
        };

        tuple
    }
}

/// SplitMix64's finalizer, so that adjacent indices give unrelated streams.
///
/// Without it, `Sampler::new(seed + index)` would hand consecutive indices
/// states one apart, and SplitMix64's first output for two nearby states is
/// well mixed but its *sequence* is the same shifted by one draw. Finalizing
/// first removes the question.
const fn stream_seed(seed: u64, index: u64) -> u64 {
    let mut z = seed ^ index.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

// ---------------------------------------------------------------------------
// The differential
// ---------------------------------------------------------------------------

/// Which compared number a report is talking about.
#[derive(Clone, Copy, Debug)]
pub struct FuzzCase {
    /// The tuple's index in the generator's stream.
    pub index: u64,
    /// The one-based `YP` subscript, matching the Fortran.
    pub component: usize,
}

impl fmt::Display for FuzzCase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "fuzz#{}[YP({})]", self.index, self.component)
    }
}

/// Why a tuple was rejected as a counterexample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Disagreement {
    /// The two disagree about whether the plant tripped
    /// (`teprob.f:702-710`), which decides whether all fifty derivatives are
    /// zeroed.
    Trip,
    /// One side produced a number and the other did not, or two infinities of
    /// opposite sign.
    Existence,
    /// The scaled error exceeded the tolerance.
    Magnitude,
}

/// A counterexample, with the number that made it one.
#[derive(Clone, Copy, Debug)]
pub struct Finding {
    /// What kind of disagreement it is.
    pub kind: Disagreement,
    /// The one-based `YP` subscript, or zero for a trip disagreement.
    pub component: usize,
    /// What the port produced.
    pub ours: f64,
    /// What the Fortran produced.
    pub theirs: f64,
    /// The balance's own error budget, from
    /// `tepsim_core::balances::Balances::scale`.
    pub scale: f64,
    /// The error over that budget.
    pub error: f64,
}

impl fmt::Display for Finding {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.kind {
            Disagreement::Trip => write!(
                f,
                "trip disagreement: the port says {}, the Fortran says {}",
                self.ours != 0.0,
                self.theirs != 0.0
            ),
            Disagreement::Existence | Disagreement::Magnitude => write!(
                f,
                "YP({}) ours {:?} theirs {:?} scale {:?} error {:e} of scale",
                self.component, self.ours, self.theirs, self.scale, self.error
            ),
        }
    }
}

/// What one tuple produced.
#[derive(Clone, Copy, Debug)]
pub enum Outcome {
    /// Both ran, both converged, and every component met the tolerance.
    Agreed {
        /// The worst scaled error over the fifty.
        worst: f64,
        /// The one-based subscript it came from.
        component: usize,
    },
    /// Both tripped, so `teprob.f:807-811` zeroed all fifty on both sides.
    /// Coverage of the freeze, and no evidence at all about the model.
    Frozen,
    /// The port's Newton solve gave up where the Fortran silently returned its
    /// guess. Delta D-001, not a counterexample.
    PortDidNotConverge,
    /// The Fortran's Newton solve gave up and the port's did not. The same
    /// delta, seen from the other side; see [`fortran_fell_through`].
    FortranFellThrough,
    /// A counterexample.
    Disagreed(Finding),
}

/// Whether `TESUB2` fell through on any of the four temperature solves.
///
/// `teprob.f:1439` executes `T=TIN` when the hundred iterations run out, so a
/// fall-through returns the caller's guess *bit for bit*. The four calls at
/// `teprob.f:460-465` take their guesses from `COMMON/TEPROC/`, which
/// [`Scenario::force`] has just written, so a returned temperature that still
/// equals the guess is the first sign of one.
///
/// It is only the first sign, and taking it for proof would throw away good
/// tuples. A solve that starts *at* the answer also returns its guess
/// unchanged, which is exactly what happens under [`Shape::Sparse`] whenever
/// the touched slots miss a vessel: the state is nominal, the seed is the
/// nominal temperature, and Newton's first step is smaller than an ULP.
///
/// So a bit-equal return is confirmed by evaluating Newton's step at the
/// returned temperature, with the Fortran's own `TESUB1` and `TESUB3`. If that
/// step already meets `teprob.f:1440`'s `1.D-12`, the solve converged and the
/// tuple is kept. If it does not, the hundred iterations really did run out.
#[must_use]
pub fn fortran_fell_through(oracle: &mut Oracle, scenario: &Scenario, snapshot: &Snapshot) -> bool {
    let before = &scenario.common;
    let after = &snapshot.common;
    // The four calls at `teprob.f:460-465`, with the composition, the target
    // specific energy and the basis each was given.
    let solves = [
        (before.tcr, after.tcr, after.xlr, after.esr, 0),
        (before.tcs, after.tcs, after.xls, after.ess, 0),
        (before.tcc, after.tcc, after.xlc, after.esc, 0),
        (before.tcv, after.tcv, after.xvv, after.esv, 2),
    ];
    for (guess, solved, z, h, ity) in solves {
        if guess.to_bits() != solved.to_bits() {
            continue;
        }
        // `teprob.f:1435-1438`, one iteration, at the temperature that came
        // back.
        let residual = oracle.tesub1(&z, solved, ity) - h;
        let slope = oracle.tesub3(&z, solved, ity);
        let step = -residual / slope;
        // Spelled through `partial_cmp` rather than as `>= 1e-12`, because a
        // NaN step is not "converged": it compares false against every bound,
        // and the reason this check exists is to catch the states where the
        // arithmetic came apart.
        let converged = step.abs().partial_cmp(&1e-12) == Some(core::cmp::Ordering::Less);
        if !converged {
            return true;
        }
    }
    false
}

/// Put a `Plant` into exactly the condition a scenario describes.
///
/// Mirrors `tests/tier2_balances.rs`, which is the canonical version and lives
/// in a test file rather than in the library. It is repeated here rather than
/// shared because Tier 8 is a library module and the Tier 2 harness is not, and
/// because copying eleven lines is cheaper than moving a file every other
/// session depends on.
///
/// Two things it deliberately does not do. It does not ask the plant for its
/// own walk inputs: those come from the oracle's `TESUB8` after the evaluation,
/// so that a bug in the port's walk cannot feed both sides of the comparison
/// and pass. And it reads the *clamped* `IDV` back out of the Fortran for the
/// two step disturbances at `teprob.f:407-408`, because that is what the
/// Fortran used, while the port keeps the raw values and clamps them itself.
///
/// # `VCV` has to come from after the call, not before it
///
/// `teprob.f:799-805` latches `VCV(I)=XMV(I)` and clamps it to `[0, 100]`
/// *immediately* before using it in `YP(I+38)`, so the value the twelve valve
/// lags were computed from is the one the call left behind. Tier 2 gets away
/// with the value from before, because on its pools `XMV` never moves and the
/// two are equal. Tier 8 draws `XMV` at random, and taking the value from
/// before made all twelve valve derivatives disagree on every tuple: 208
/// counterexamples out of 400, before a line of the model was in question.
///
/// The consequence is that Tier 8 compares the valve *lag* and not the valve
/// *latch*; the latch is `tests/hoist_valve_latch.rs`.
fn configure(oracle: &mut Oracle, scenario: &Scenario, snapshot: &Snapshot) -> (Plant, Inputs) {
    let t = scenario.time;
    // Post-call, so this is `IDV` as `teprob.f:341-346` left it.
    let clamped = oracle.disturbances();
    let raw = |n: usize| f64::from(clamped[n - 1]);

    // teprob.f:407-408 subtracts two terms from A, on two source lines.
    let a = oracle.tesub8(1, t) - raw(1) * 0.03 - raw(2) * 2.43719e-3;
    let b = oracle.tesub8(2, t) + raw(2) * 0.005;

    let walks = plant::WalkInputs {
        feed: tepsim_core::streams::FeedConditions {
            ac_feed_light: [a, b, 1.0 - a - b],
            d_feed_celsius: oracle.tesub8(3, t) + raw(3) * 5.0,
            ac_feed_celsius: oracle.tesub8(4, t),
        },
        flow: tepsim_core::flows::FlowDrift {
            steam_capacity: oracle.tesub8(9, t),
            reactor_outlet: oracle.tesub8(12, t),
        },
        heat: tepsim_core::heat::HeatDrift {
            reactor_coolant: oracle.tesub8(10, t),
            condenser_coolant: oracle.tesub8(11, t),
        },
        reaction: tepsim_core::kinetics::ReactionDrift {
            first: oracle.tesub8(7, t),
            second: oracle.tesub8(8, t),
        },
        // teprob.f:413-414.
        coolant: tepsim_core::balances::CoolantInlet {
            reactor: oracle.tesub8(5, t) + raw(4) * 5.0,
            condenser: oracle.tesub8(6, t) + raw(5) * 5.0,
        },
    };

    let mut p = Plant::new();
    p.set_walk_inputs(walks);
    p.set_seeds(vessels::TemperatureSeeds {
        reactor: scenario.common.tcr,
        separator: scenario.common.tcs,
        stripper: scenario.common.tcc,
        mixing: scenario.common.tcv,
    });
    p.set_valve_command(snapshot.common.vcv);

    let inputs = Inputs {
        manipulated: scenario.manipulated,
        disturbances: core::array::from_fn(|i| f64::from(scenario.disturbances[i])),
    };
    (p, inputs)
}

/// The error of one comparison against the balance's own scale.
///
/// A zero scale means the balance has no terms at all, so the only acceptable
/// answer is bit equality. Matches
/// [`Comparison::observe_against`](crate::tier1::Comparison::observe_against).
#[must_use]
pub fn scaled_error(ours: f64, theirs: f64, scale: f64) -> f64 {
    if scale == 0.0 {
        if ours.to_bits() == theirs.to_bits() {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        (ours - theirs).abs() / scale.abs()
    }
}

/// One component's verdict: `None` when the two agree, the error when they do
/// not.
///
/// See the module documentation for why two `NaN`s agree and an infinity
/// against a number does not.
#[must_use]
pub fn slot_verdict(
    ours: f64,
    theirs: f64,
    scale: f64,
    tolerance: f64,
) -> Option<(f64, Disagreement)> {
    match (ours.is_nan(), theirs.is_nan()) {
        (true, true) => None,
        (true, false) | (false, true) => Some((f64::INFINITY, Disagreement::Existence)),
        (false, false) => {
            if ours.is_infinite() || theirs.is_infinite() {
                return if ours.to_bits() == theirs.to_bits() {
                    None
                } else {
                    Some((f64::INFINITY, Disagreement::Existence))
                };
            }
            let error = scaled_error(ours, theirs, scale);
            (error > tolerance).then_some((error, Disagreement::Magnitude))
        }
    }
}

/// Force both implementations into `tuple` and compare the fifty derivatives.
///
/// `mutation` is normally [`Mutation::None`]; the other variants exist so the
/// harness can be shown to have teeth. See the module documentation.
///
/// `sink`, when given, receives every finite-against-finite pair so the caller
/// can accumulate a ULP histogram. Non-finite pairs are withheld from it on
/// purpose: `Comparison` counts two `NaN`s with different payloads as a
/// mismatch, and on a wild state that is not a claim about the model.
pub fn check(
    oracle: &mut Oracle,
    base: &Scenario,
    tuple: &Tuple,
    tolerance: f64,
    mutation: Mutation,
    index: u64,
    mut sink: Option<&mut Comparison<FuzzCase>>,
) -> Outcome {
    let scenario = tuple.scenario(base);
    let snapshot = scenario.force(oracle);
    if fortran_fell_through(oracle, &scenario, &snapshot) {
        return Outcome::FortranFellThrough;
    }
    let (plant, inputs) = configure(oracle, &scenario, &snapshot);
    let state = State::from_flat(&scenario.state);

    let Ok((derivative, scale, signals)) =
        plant.derivatives_with_scale(SimTime(scenario.time), &state, &inputs)
    else {
        return Outcome::PortDidNotConverge;
    };

    if signals.shutdown.is_tripped() != snapshot.tripped {
        return Outcome::Disagreed(Finding {
            kind: Disagreement::Trip,
            component: 0,
            ours: f64::from(u8::from(signals.shutdown.is_tripped())),
            theirs: f64::from(u8::from(snapshot.tripped)),
            scale: 1.0,
            error: f64::INFINITY,
        });
    }

    let mut ours = derivative.to_flat();
    let budgets = scale.to_flat();
    mutation.apply(&mut ours, &budgets, &snapshot);

    let mut worst = 0.0_f64;
    let mut worst_component = 0_usize;
    let mut finding: Option<Finding> = None;
    for (slot, ((ours, theirs), budget)) in ours
        .iter()
        .zip(snapshot.derivative.iter())
        .zip(budgets.iter())
        .enumerate()
    {
        let case = FuzzCase {
            index,
            component: slot + 1,
        };
        if let Some(comparison) = sink.as_deref_mut()
            && ours.is_finite()
            && theirs.is_finite()
        {
            comparison.observe_against(case, *ours, *theirs, *budget);
        }
        match slot_verdict(*ours, *theirs, *budget, tolerance) {
            Some((error, kind)) => {
                if finding.is_none() || error > worst {
                    finding = Some(Finding {
                        kind,
                        component: slot + 1,
                        ours: *ours,
                        theirs: *theirs,
                        scale: *budget,
                        error,
                    });
                    worst = error;
                    worst_component = slot + 1;
                }
            }
            None if finding.is_none() => {
                // Only meaningful while nothing has failed: once there is a
                // finding, `worst` belongs to it.
                let error = scaled_error(*ours, *theirs, *budget);
                if error.is_finite() && error > worst {
                    worst = error;
                    worst_component = slot + 1;
                }
            }
            None => {}
        }
    }

    if let Some(finding) = finding {
        return Outcome::Disagreed(finding);
    }
    if snapshot.tripped {
        // Fifty zeros on both sides. `Balances::scale` is zeroed too, so the
        // comparison above already required bit equality; what it did not do is
        // prove anything about the model.
        return Outcome::Frozen;
    }
    Outcome::Agreed {
        worst,
        component: worst_component,
    }
}

// ---------------------------------------------------------------------------
// Shrinking
// ---------------------------------------------------------------------------

/// The result of shrinking a counterexample.
#[derive(Clone, Copy, Debug)]
pub struct Shrunk {
    /// The minimal tuple that still fails.
    pub tuple: Tuple,
    /// How many knobs it still differs from nominal in.
    pub knobs: usize,
    /// How many knobs the tuple it started from differed in.
    pub knobs_before: usize,
    /// How many differential evaluations the search cost.
    pub evaluations: usize,
    /// What the minimal tuple fails with.
    pub finding: Finding,
}

/// Shrink a failing tuple toward nominal, one knob at a time.
///
/// Two passes, each to a fixpoint. The first resets each knob to its nominal
/// value and keeps the reset when the failure survives, which is what removes
/// the forty-odd irrelevant slots a generated tuple carries. The second bisects
/// each surviving knob toward nominal, geometrically for the state slots, which
/// is what turns "the reactor holdup is a thousand times nominal" into "it is
/// 1.6 times nominal, and here is the branch that crosses at 1.55".
///
/// The predicate is "still fails", not "still fails the same way". A shrink
/// that wandered from one fault to another would still hand back a reproducer,
/// which is what the corpus wants.
///
/// # Panics
///
/// If `failing` does not actually fail, which would make the whole search
/// meaningless.
pub fn shrink(
    oracle: &mut Oracle,
    base: &Scenario,
    nominal: &Tuple,
    failing: &Tuple,
    tolerance: f64,
    mutation: Mutation,
) -> Shrunk {
    let mut evaluations = 0;
    let fails = |oracle: &mut Oracle, t: &Tuple, evaluations: &mut usize| -> Option<Finding> {
        *evaluations += 1;
        match check(oracle, base, t, tolerance, mutation, 0, None) {
            Outcome::Disagreed(f) => Some(f),
            _ => None,
        }
    };

    let Some(mut finding) = fails(oracle, failing, &mut evaluations) else {
        panic!("shrink was handed a tuple that does not fail");
    };
    let mut best = *failing;
    let knobs_before = best.distance(nominal);

    // Pass one: reset knobs outright. Capped, because each round is eighty-four
    // evaluations and a pathological case could otherwise loop for as many
    // rounds as there are knobs.
    for _ in 0..8 {
        let mut improved = false;
        for k in 0..N_KNOBS {
            if best.knob_matches(nominal, k) {
                continue;
            }
            let mut candidate = best;
            candidate.take_knob(nominal, k);
            if let Some(f) = fails(oracle, &candidate, &mut evaluations) {
                best = candidate;
                finding = f;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }

    // Pass two: bisect what is left toward nominal.
    for _ in 0..40 {
        let mut improved = false;
        for k in 0..N_KNOBS {
            if best.knob_matches(nominal, k) {
                continue;
            }
            let mut candidate = best;
            if !candidate.halve_knob(nominal, k) {
                continue;
            }
            if let Some(f) = fails(oracle, &candidate, &mut evaluations) {
                best = candidate;
                finding = f;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }

    Shrunk {
        knobs: best.distance(nominal),
        knobs_before,
        tuple: best,
        evaluations,
        finding,
    }
}

// ---------------------------------------------------------------------------
// Teeth
// ---------------------------------------------------------------------------

/// A branch of `TEFUNC`, for a mutation that fires only on one side of it.
///
/// Read off the oracle's own snapshot, which is legitimate: the mutant stands
/// in for a port that mishandles the branch, and the Fortran is the authority
/// on which side of it a state sits.
/// Both of these are branches the *nominal* plant does not take. That is the
/// whole point: a mutation behind a branch the nominal plant already sits on
/// would fire on almost every tuple and would demonstrate nothing beyond
/// [`Mutation::WrongConstant`]. `VLR/7.8` is about 75 at nominal, on the far
/// side of `teprob.f:663`, which is how the first attempt at this measured 208
/// counterexamples out of 208 compared tuples.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Branch {
    /// `teprob.f:665`. The reactor heat-transfer coefficient pins to zero below
    /// `VLR/7.8 = 10`, which needs the reactor most of the way drained.
    ReactorBelowHeatTransferRamp,
    /// `teprob.f:615`. The stripping factor goes linear above 170 C.
    StripperAbove170C,
}

impl Branch {
    /// Whether the snapshot sits on the far side of this branch.
    #[must_use]
    pub fn taken(self, snapshot: &Snapshot) -> bool {
        match self {
            Self::ReactorBelowHeatTransferRamp => snapshot.common.vlr / 7.8 < 10.0,
            Self::StripperAbove170C => snapshot.common.tcc > 170.0,
        }
    }
}

/// A deliberate corruption of the port's answer, so that "no counterexamples"
/// can be shown to mean something.
///
/// The magnitude is a fraction of the balance's own scale rather than of its
/// value, because that is the quantity the gate is written against. A mutation
/// of `1e-11` is ten times the gate and must be caught; `1e-13` is a tenth of
/// it and must not be.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Mutation {
    /// The port exactly as it ships.
    None,
    /// A term of the wrong size in one balance, on every state. Stands in for a
    /// mistyped constant.
    WrongConstant {
        /// The one-based `YP` subscript to corrupt.
        component: usize,
        /// The size of the error, as a fraction of that balance's scale.
        relative: f64,
    },
    /// A term of the wrong size in one balance, but only when a branch is
    /// taken. Stands in for the case this tier exists for: a branch that the
    /// hand-built pools never enter.
    DroppedTermBehindBranch {
        /// The one-based `YP` subscript to corrupt.
        component: usize,
        /// The size of the error, as a fraction of that balance's scale.
        relative: f64,
        /// Which branch has to be taken for it to fire.
        branch: Branch,
    },
}

impl Mutation {
    /// Apply the corruption in place.
    #[allow(
        clippy::suboptimal_flops,
        reason = "the added term is the mutation; fusing it would change what is injected"
    )]
    fn apply(self, ours: &mut [f64; N_STATES], scale: &[f64; N_STATES], snapshot: &Snapshot) {
        let (component, relative) = match self {
            Self::None => return,
            Self::WrongConstant {
                component,
                relative,
            } => (component, relative),
            Self::DroppedTermBehindBranch {
                component,
                relative,
                branch,
            } => {
                if !branch.taken(snapshot) {
                    return;
                }
                (component, relative)
            }
        };
        let slot = component - 1;
        // A frozen plant has a zero scale (`teprob.f:807-811` zeroes the
        // derivative and `Balances::scale` follows), so this injects nothing
        // there. That is correct: the freeze is not what is being mutated.
        ours[slot] += relative * scale[slot].abs();
    }
}

// ---------------------------------------------------------------------------
// The run
// ---------------------------------------------------------------------------

/// How many tuples to run, from `TEP_TIER8`.
///
/// The same idiom as `TEP_TIER5` and `TEP_TIER7`. Unset is [`Budget::SMOKE`],
/// which is what belongs in a per-commit gate; `full` is [`Budget::FULL`],
/// which is a nightly. A decimal count is also accepted, for a bisection run
/// that wants a particular size.
///
/// This project has put an expensive test in the per-commit gate four times.
/// The default here is a few hundred tuples for that reason and no other.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Budget {
    /// How many tuples to generate.
    pub tuples: u64,
}

impl Budget {
    /// The environment variable that selects the size.
    pub const ENV: &'static str = "TEP_TIER8";

    /// What `cargo xtask ci` runs: enough to exercise every path in the
    /// harness, and a second of wall clock.
    pub const SMOKE: Self = Self { tuples: 400 };

    /// What a nightly runs. About six minutes at the measured 68 microseconds
    /// per tuple.
    pub const FULL: Self = Self { tuples: 5_000_000 };

    /// Which budget to run.
    #[must_use]
    pub fn selected() -> Self {
        match std::env::var(Self::ENV).as_deref() {
            Ok("full") => Self::FULL,
            Ok(count) => count.parse::<u64>().map_or(Self::SMOKE, |tuples| Self {
                tuples: tuples.max(1),
            }),
            Err(_) => Self::SMOKE,
        }
    }
}

/// The census of one fuzzing run.
#[derive(Clone, Debug)]
pub struct Report {
    /// The seed the tuples came from.
    pub seed: u64,
    /// How many tuples were generated.
    pub tuples: u64,
    /// How many produced a running, converged comparison of all fifty
    /// components. This is the number that carries the claim.
    pub compared: u64,
    /// How many tripped, so `teprob.f:807-811` zeroed both sides.
    pub frozen: u64,
    /// How many the port's Newton solve gave up on. Delta D-001.
    pub port_did_not_converge: u64,
    /// How many the Fortran's Newton solve gave up on. The same delta.
    pub fortran_fell_through: u64,
    /// How many counterexamples were found in total.
    ///
    /// Counted separately from [`Report::counterexamples`], which stops
    /// collecting so that a systematically broken port cannot fill memory with
    /// a million copies of one fault.
    pub disagreed: u64,
    /// The counterexamples that were kept, with the index that produced each.
    pub counterexamples: Vec<(u64, Tuple, Finding)>,
    /// The worst scaled error among the tuples that were compared.
    pub worst: f64,
    /// The index and component it came from.
    pub worst_at: (u64, usize),
    /// The ULP histogram over every finite-against-finite pair.
    pub comparison: Comparison<FuzzCase>,
}

impl Report {
    /// The fraction of generated tuples that reached a comparison of the model.
    ///
    /// Frozen states are excluded: fifty zeros on both sides is not evidence.
    #[must_use]
    pub fn physical_fraction(&self) -> f64 {
        if self.tuples == 0 {
            return 0.0;
        }
        self.compared as f64 / self.tuples as f64
    }
}

impl fmt::Display for Report {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "tier8 seed 0x{:X}", self.seed)?;
        writeln!(f, "  tuples             : {}", self.tuples)?;
        writeln!(
            f,
            "  compared           : {} ({:.1}% of tuples)",
            self.compared,
            100.0 * self.physical_fraction()
        )?;
        writeln!(f, "  frozen (tripped)   : {}", self.frozen)?;
        writeln!(
            f,
            "  port did not solve : {} (delta D-001)",
            self.port_did_not_converge
        )?;
        writeln!(
            f,
            "  Fortran fell through: {} (delta D-001)",
            self.fortran_fell_through
        )?;
        writeln!(
            f,
            "  counterexamples    : {} ({} kept)",
            self.disagreed,
            self.counterexamples.len()
        )?;
        writeln!(
            f,
            "  worst err/scale    : {:.3e} at fuzz#{}[YP({})]",
            self.worst, self.worst_at.0, self.worst_at.1
        )?;
        write!(f, "{}", self.comparison)
    }
}

/// Run `budget.tuples` tuples of `generator` against the oracle.
///
/// Stops collecting counterexamples at `keep`, so that a systematically broken
/// port does not fill memory with a million copies of the same fault, but keeps
/// counting them.
pub fn run(
    oracle: &mut Oracle,
    base: &Scenario,
    generator: &Generator,
    budget: Budget,
    tolerance: f64,
    mutation: Mutation,
    keep: usize,
) -> Report {
    let mut report = Report {
        seed: generator.seed(),
        tuples: budget.tuples,
        compared: 0,
        frozen: 0,
        port_did_not_converge: 0,
        fortran_fell_through: 0,
        disagreed: 0,
        counterexamples: Vec::new(),
        worst: 0.0,
        worst_at: (0, 0),
        comparison: Comparison::new("YP(1..50) over random tuples, against the scale of the terms"),
    };

    for index in 0..budget.tuples {
        let tuple = generator.tuple(index);
        let outcome = check(
            oracle,
            base,
            &tuple,
            tolerance,
            mutation,
            index,
            Some(&mut report.comparison),
        );
        match outcome {
            Outcome::Agreed { worst, component } => {
                report.compared += 1;
                if worst > report.worst {
                    report.worst = worst;
                    report.worst_at = (index, component);
                }
            }
            Outcome::Frozen => report.frozen += 1,
            Outcome::PortDidNotConverge => report.port_did_not_converge += 1,
            Outcome::FortranFellThrough => report.fortran_fell_through += 1,
            Outcome::Disagreed(finding) => {
                report.disagreed += 1;
                if report.counterexamples.len() < keep {
                    report.counterexamples.push((index, tuple, finding));
                }
            }
        }
    }
    report
}
