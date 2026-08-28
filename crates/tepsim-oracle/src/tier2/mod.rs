//! The Tier 2 measuring apparatus: forcing the Fortran into a chosen state,
//! and diffing any named `COMMON/TEPROC/` field against the port.
//!
//! Tier 1 proved the four utility routines. Tier 2 is the real evidence that
//! the plant is ported correctly: force both implementations into the identical
//! state, evaluate once, and compare. `PLAN.org` describes it as a comparison
//! of the fifty derivative components, and it ends there, but it does not have
//! to *start* there.
//!
//! # Every intermediate is observable, so no item has to be written blind
//!
//! `COMMON/TEPROC/` holds the plant's entire working set: `PPR`, `RR`, `CRXR`,
//! `FTM`, `FCM`, `HST`, `SFR`, `QUR` and the rest. The oracle already mirrors
//! all of it. So each Phase 2 item can diff the fields it computes as soon as
//! it is written, instead of nine items landing unverified and the first signal
//! arriving at the end. That is what [`Snapshot`] and [`compare_field`] are
//! for.
//!
//! # `TEFUNC` is not reproducible unless you make it so
//!
//! Evaluating the Fortran twice at the same state does *not* give the same
//! answer. `TEFUNC` advances the twelve disturbance walks, draws measurement
//! noise, and ticks the sampled analysers, all as side effects of what presents
//! itself as a derivative evaluation. A harness that ignored this would report
//! differences that were entirely its own doing, and Tier 2 would degenerate
//! into noise that someone would eventually "fix" by loosening a tolerance.
//!
//! [`Scenario::force`] therefore restores the walk state and the generator word
//! before every evaluation, and [`reproducible`] asserts that this is
//! sufficient rather than assuming it. It was not: see below.
//!
//! # The four solved temperatures are hidden state
//!
//! `TESUB2` takes its temperature argument as *both* the initial guess and the
//! result (`teprob.f:1432`, `1438`). The four call sites at `teprob.f:460-465`
//! pass `TCR`, `TCS`, `TCC` and `TCV`, which live in `COMMON/TEPROC/` and
//! therefore survive between calls. Every derivative evaluation warm-starts its
//! Newton solves from *the previous evaluation's answers*.
//!
//! Newton stops when the step falls below 1e-12, so the converged value depends
//! on where it started, in the last bits. That makes `TEFUNC` path-dependent at
//! the 1e-13 level even with the walks and the generator pinned. Measured here:
//! evaluating the same state twice with one unrelated evaluation in between
//! moved `TST(6)`, which is `TCV`, from 86.12011310855937 to 86.12011310855931,
//! and that propagated through `PTV` into `FTM(6)` and on into `YP(27)`.
//!
//! A [`Scenario`] therefore carries the entire `COMMON/TEPROC/` block, not a
//! chosen subset of it. Enumerating which fields carry over is exactly the kind
//! of judgement that was wrong the first time; restoring all of it cannot be.
//! It also picks up `VCV`, which the valve latch carries over, and the
//! analyser state `XDEL`, `TGAS` and `TPROD`.
//!
//! The consequence for the port is not confined to the harness. The four
//! temperatures are genuinely part of the plant's persistent state, alongside
//! the latched valve commands, and B-0017 must carry them.
//!
//! # A tripped plant is not evidence of anything
//!
//! `teprob.f:807-811` zeroes all fifty derivatives when any of the eight
//! shutdown conditions fires. A sampled state that trips therefore has an
//! all-zero derivative, which *any* port reproduces, correct or not. Such
//! states are counted and reported separately by [`Snapshot::tripped`]; they
//! are coverage, not evidence, and a Tier 2 result that did not say how many it
//! contained would be overstating itself.

pub mod adversarial;

use tepsim_core::state::N_STATES;

use crate::tier1::{Comparison, Sampler};
use crate::{Oracle, Teproc, Wlk};

/// A fully specified starting condition: everything `TEFUNC` reads.
///
/// If two implementations are put into this and evaluated once, any difference
/// in the answer is a difference in the model.
#[derive(Clone, Debug)]
pub struct Scenario {
    /// Simulation time, in hours.
    pub time: f64,
    /// The fifty integrated states, `YY`.
    pub state: [f64; N_STATES],
    /// The twelve manipulated variables, `XMV`.
    pub manipulated: [f64; 12],
    /// The twenty disturbance flags, `IDV`.
    pub disturbances: [i32; 20],
    /// The disturbance walk state, `COMMON/WLK/`.
    pub walk: Wlk,
    /// The generator word, `COMMON/RANDSD/ G`.
    pub rng: f64,
    /// `XMEAS(1..41)` as the previous evaluation left them.
    ///
    /// Only 23 through 41 matter, and only to the analysers: `teprob.f:744`
    /// and `755` write them inside a schedule check, so between samples the
    /// previous *reported* value persists and is read again. They do not reach
    /// the derivative at all, which is why Tier 2 got away without them; a
    /// whole-step comparison does not.
    pub measurements: [f64; 41],
    /// The whole plant working set carried in from the previous evaluation.
    ///
    /// Most of this is recomputed on every call and restoring it changes
    /// nothing. Four fields are not: `TCR`, `TCS`, `TCC` and `TCV` seed the
    /// Newton solves, `VCV` holds the latched valve commands, and `XDEL`,
    /// `TGAS` and `TPROD` hold the analyser state. Restoring the block whole
    /// is cheaper than being right about which.
    pub common: Teproc,
}

/// What one evaluation produced.
#[derive(Clone, Debug)]
pub struct Snapshot {
    /// The fifty derivatives, `YP`.
    pub derivative: [f64; N_STATES],
    /// The whole plant working set after the call.
    pub common: Teproc,
    /// The forty-one measurements, noise included.
    pub measurements: [f64; 41],
    /// Whether any of the eight shutdown conditions fired
    /// (`teprob.f:702-710`). If so, the derivative is all zeros and proves
    /// nothing.
    pub tripped: bool,
}

impl Scenario {
    /// Put the Fortran into this condition and evaluate once.
    ///
    /// Restores the walk state and the generator word first, so that repeated
    /// calls with the same scenario give the same answer. See the module
    /// documentation for why that is not automatic.
    pub fn force(&self, oracle: &mut Oracle) -> Snapshot {
        oracle.set_teproc(&self.common);
        oracle.set_wlk(&self.walk);
        oracle.set_rng(self.rng);
        oracle.set_measurements(&self.measurements);
        oracle.set_manipulated(&self.manipulated);
        oracle.set_disturbances(&self.disturbances);
        let derivative = oracle.derivatives(self.time, &self.state);
        Snapshot {
            derivative,
            common: oracle.teproc(),
            measurements: oracle.measurements(),
            tripped: oracle.shutdown_flag() != 0,
        }
    }
}

/// Where a sampled state came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pool {
    /// States visited by the nominal closed-loop trajectory.
    Nominal,
    /// Random perturbations of a nominal state, scaled across several orders of
    /// magnitude.
    Perturbed,
    /// States placed deliberately at a discontinuity or a clamp. B-0016.
    Adversarial,
}

impl core::fmt::Display for Pool {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Pool::Nominal => "nominal",
            Pool::Perturbed => "perturbed",
            Pool::Adversarial => "adversarial",
        })
    }
}

/// Identifies one compared number: which sampled state, and which component of
/// the field.
#[derive(Clone, Copy, Debug)]
pub struct Case {
    /// Which pool the state came from.
    pub pool: Pool,
    /// Its index within that pool.
    pub index: usize,
    /// The one-based component within the field, matching Fortran subscripts.
    pub component: usize,
}

impl core::fmt::Display for Case {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}#{}[{}]", self.pool, self.index, self.component)
    }
}

/// Record one field's worth of comparison, component by component.
///
/// The component index in the report is one-based, so it reads as the Fortran
/// subscript rather than as a Rust offset. `YP(27)` in a log entry should be
/// findable in `teprob.f` without arithmetic.
pub fn compare_field(
    comparison: &mut Comparison<Case>,
    pool: Pool,
    index: usize,
    ours: &[f64],
    theirs: &[f64],
) {
    assert_eq!(
        ours.len(),
        theirs.len(),
        "comparing fields of different lengths is a harness bug, not a finding"
    );
    for (offset, (a, b)) in ours.iter().zip(theirs.iter()).enumerate() {
        comparison.observe(
            Case {
                pool,
                index,
                component: offset + 1,
            },
            *a,
            *b,
        );
    }
}

/// Assert that forcing a scenario twice gives bit-identical results.
///
/// The property [`Scenario::force`] exists to provide. If it fails, every Tier
/// 2 number is meaningless, so it is checked directly rather than inferred.
///
/// # Panics
///
/// If any of the fifty derivatives, or the shutdown flag, differs between two
/// evaluations of the same scenario.
pub fn reproducible(oracle: &mut Oracle, scenario: &Scenario) {
    let first = scenario.force(oracle);
    // Evaluate something else in between, so a stale-state bug cannot pass by
    // the Fortran simply not having moved.
    let mut disturbed = scenario.clone();
    disturbed.state[0] *= 1.01;
    let _ = disturbed.force(oracle);
    let second = scenario.force(oracle);

    for (slot, (a, b)) in first
        .derivative
        .iter()
        .zip(second.derivative.iter())
        .enumerate()
    {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "YP({}) differs between two forced evaluations of the same \
             scenario: {a:?} then {b:?}. Something is carrying over between \
             evaluations that the scenario does not restore.",
            slot + 1
        );
    }
    assert_eq!(
        first.tripped, second.tripped,
        "the shutdown flag is not reproducible"
    );
}

/// Generators for the two sampling pools that do not need hand construction.
///
/// The third, the adversarial catalogue, is enumerable from the source and is
/// built by hand in B-0016, for the same reason Tier 1's boundary pool was
/// worth its weight: the interesting numerics live at the discontinuities.
#[derive(Clone, Debug)]
pub struct Pools {
    /// The nominal starting condition, from `TEINIT`.
    pub nominal: Scenario,
    /// States along the nominal trajectory, each with the time it was visited.
    ///
    /// The time matters and is not decoration. `teprob.f:397-406` resets the
    /// whole walk state whenever `TIME.EQ.0.D0`, so a pool that evaluated every
    /// state at t=0 would silently never exercise the disturbance walks at all.
    pub trajectory: Vec<(f64, [f64; N_STATES])>,
}

impl Pools {
    /// Run the Fortran forward from `TEINIT` with fixed-step Euler, recording
    /// the state before each step.
    ///
    /// This is the original's own integrator (`temain_mod.f`'s `INTGTR`), so
    /// the states visited are the ones a real run visits, which is the point of
    /// the pool.
    pub fn collect(oracle: &mut Oracle, steps: usize, dt: f64) -> Self {
        let (time, state) = oracle.init();
        let nominal = Scenario {
            time,
            state,
            manipulated: oracle.manipulated(),
            disturbances: oracle.disturbances(),
            walk: oracle.wlk(),
            rng: oracle.rng(),
            measurements: oracle.measurements(),
            common: oracle.teproc(),
        };

        let mut trajectory = Vec::with_capacity(steps);
        let mut y = state;
        let mut t = time;
        for _ in 0..steps {
            trajectory.push((t, y));
            let yp = oracle.derivatives(t, &y);
            for (state, rate) in y.iter_mut().zip(yp) {
                // Two roundings, not a fused multiply-add. This reproduces
                // `INTGTR`, and gfortran at the pinned `-O0` does not fuse.
                #[allow(clippy::suboptimal_flops, reason = "matches INTGTR's rounding")]
                {
                    *state += dt * rate;
                }
            }
            t += dt;
        }
        Self {
            nominal,
            trajectory,
        }
    }

    /// A scenario taken from the recorded trajectory.
    ///
    /// The walk state and generator word are reset to the nominal ones so that
    /// the scenario is self-contained: it can be evaluated in any order, any
    /// number of times, and give the same answer.
    #[must_use]
    pub fn nominal_case(&self, index: usize) -> Scenario {
        let (time, state) = self.trajectory[index % self.trajectory.len()];
        Scenario {
            time,
            state,
            ..self.nominal.clone()
        }
    }

    /// A nominal state with every component multiplied by `1 + e`, where `e` is
    /// drawn uniformly on a logarithmic scale between 1e-9 and 1e-1.
    ///
    /// Scaling rather than adding, so that a component's perturbation is
    /// meaningful relative to its own magnitude: the holdups run to hundreds
    /// and the valve positions to a hundred, but the specific energies are
    /// order one, and a fixed absolute perturbation would be a rounding error
    /// for some and a catastrophe for others.
    #[must_use]
    pub fn perturbed_case(&self, index: usize, sampler: &mut Sampler) -> Scenario {
        let mut scenario = self.nominal_case(index);
        for slot in &mut scenario.state {
            // 1e-9 to 1e-1, log-uniform, with a random sign.
            let magnitude = 1e-9 * 1e8_f64.powf(sampler.unit());
            let signed = if sampler.unit() < 0.5 {
                -magnitude
            } else {
                magnitude
            };
            *slot *= 1.0 + signed;
        }
        scenario
    }
}
