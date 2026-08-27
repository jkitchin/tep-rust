//! Tier 2 for the vapour-liquid equilibrium: `tepsim_core::equilibrium`
//! against `teprob.f:473-502`, field by field, over all three sampling pools.
//!
//! # This test is run twice, and the two runs answer different questions
//!
//! `teprob.f:485` and `teprob.f:488` are the model's first transcendental
//! calls. The port answers them from the vendored pure-Rust `libm` so that
//! x86-64, aarch64 and wasm32 agree bit for bit, and that library does not
//! round identically to the one gfortran links. So a difference here has two
//! possible causes, and separating them is the whole job of this file.
//!
//! - **Default build.** Vendored `libm`. Asserts the Tier 2 gate of 1e-12 and
//!   records the ULP histogram. This is the configuration that ships, so this
//!   is the number that describes the port.
//! - **`--features oracle,libm-system`.** Platform libm, which on the
//!   recording machine *is* what gfortran calls. Asserts bit equality. Any
//!   failure here is the algebra, because the transcendental has been taken
//!   out of the comparison.
//!
//! Without the second run, the bit-exactness assertion that held through
//! B-0017 would have to be dropped at the first call to `exp`, and every later
//! Phase 2 item would be checked against a 1e-12 tolerance with four orders of
//! magnitude of room to hide a real mistake in.

#![cfg(feature = "oracle")]

use tepsim_core::{State, equilibrium, math, vessels};
use tepsim_oracle::Oracle;
use tepsim_oracle::tier1::Comparison;
use tepsim_oracle::tier2::{Case, Pool, Pools, Scenario, Snapshot, adversarial, compare_field};

const DT: f64 = 1.0 / 3600.0;

/// `PLAN.org`, "Tier 2".
const TIER2_TOLERANCE: f64 = 1e-12;

/// Every field this item computes, paired with where to read it from each side.
struct Field {
    name: &'static str,
    ours: fn(&equilibrium::Equilibrium) -> Vec<f64>,
    theirs: fn(&Snapshot) -> Vec<f64>,
}

fn fields() -> Vec<Field> {
    vec![
        Field {
            name: "VVR/VVS",
            ours: |e| vec![e.reactor.volume, e.separator.volume],
            theirs: |s| vec![s.common.vvr, s.common.vvs],
        },
        Field {
            name: "PPR",
            ours: |e| e.reactor.partial.as_array().to_vec(),
            theirs: |s| s.common.ppr.to_vec(),
        },
        Field {
            name: "PPS",
            ours: |e| e.separator.partial.as_array().to_vec(),
            theirs: |s| s.common.pps.to_vec(),
        },
        Field {
            name: "PTR/PTS/PTV",
            ours: |e| vec![e.reactor.pressure, e.separator.pressure, e.mixing_pressure],
            theirs: |s| vec![s.common.ptr, s.common.pts, s.common.ptv],
        },
        Field {
            name: "XVR",
            ours: |e| e.reactor.fractions.fractions().as_array().to_vec(),
            theirs: |s| s.common.xvr.to_vec(),
        },
        Field {
            name: "XVS",
            ours: |e| e.separator.fractions.fractions().as_array().to_vec(),
            theirs: |s| s.common.xvs.to_vec(),
        },
        Field {
            name: "UTVR/UTVS",
            ours: |e| vec![e.reactor.total, e.separator.total],
            theirs: |s| vec![s.common.utvr, s.common.utvs],
        },
        Field {
            name: "UCVR",
            ours: |e| e.reactor.moles.as_array().to_vec(),
            theirs: |s| s.common.ucvr.to_vec(),
        },
        Field {
            name: "UCVS",
            ours: |e| e.separator.moles.as_array().to_vec(),
            theirs: |s| s.common.ucvs.to_vec(),
        },
    ]
}

/// Force the scenario, run the port through unpack and then equilibrium with
/// the same warm-start seeds, and record every field.
///
/// The seeds come out of the scenario's own `COMMON`, which is what the
/// Fortran will use. See `tier2_unpack.rs` for why that is not optional.
fn observe(
    oracle: &mut Oracle,
    scenario: &Scenario,
    pool: Pool,
    index: usize,
    comparisons: &mut [Comparison<Case>],
) -> bool {
    let snapshot = scenario.force(oracle);
    let seeds = vessels::TemperatureSeeds {
        reactor: scenario.common.tcr,
        separator: scenario.common.tcs,
        stripper: scenario.common.tcc,
        mixing: scenario.common.tcv,
    };
    let state = State::from_flat(&scenario.state);
    let Ok(unpacked) = vessels::unpack(&state, seeds) else {
        return false;
    };
    let solved = equilibrium::equilibrium(&unpacked);
    for (field, comparison) in fields().iter().zip(comparisons.iter_mut()) {
        compare_field(
            comparison,
            pool,
            index,
            &(field.ours)(&solved),
            &(field.theirs)(&snapshot),
        );
    }
    true
}

fn fresh() -> Vec<Comparison<Case>> {
    fields()
        .iter()
        .map(|f| Comparison::new(format!("equilibrium {}", f.name)))
        .collect()
}

/// Sweep all three pools into a fresh set of comparisons.
fn sweep(
    oracle: &mut Oracle,
    steps: usize,
    perturbations: usize,
    seed: u64,
) -> (Vec<Comparison<Case>>, usize) {
    let pools = Pools::collect(oracle, steps, DT);
    let mut comparisons = fresh();
    let mut skipped = 0;

    for index in 0..pools.trajectory.len() {
        if !observe(
            oracle,
            &pools.nominal_case(index),
            Pool::Nominal,
            index,
            &mut comparisons,
        ) {
            skipped += 1;
        }
    }

    let mut sampler = tepsim_oracle::tier1::Sampler::new(seed);
    for index in 0..perturbations {
        if !observe(
            oracle,
            &pools.perturbed_case(index, &mut sampler),
            Pool::Perturbed,
            index,
            &mut comparisons,
        ) {
            skipped += 1;
        }
    }

    let (boundaries, missed) = adversarial::build(oracle, &pools.nominal_case(0));
    assert!(missed.is_empty(), "the catalogue regressed: {missed:?}");
    for (index, boundary) in boundaries.iter().enumerate() {
        if !observe(
            oracle,
            &boundary.scenario,
            Pool::Adversarial,
            index,
            &mut comparisons,
        ) {
            skipped += 1;
        }
    }

    (comparisons, skipped)
}

#[test]
fn the_equilibrium_matches_the_fortran_over_all_three_pools() {
    let mut oracle = Oracle::lock();
    let (comparisons, skipped) = sweep(&mut oracle, 400, 2_000, 0x7E2_0018);

    println!(
        "exp comes from the {} libm",
        if math::USES_SYSTEM_LIBM {
            "platform"
        } else {
            "vendored"
        }
    );
    for comparison in &comparisons {
        println!("{comparison}");
    }
    println!("states skipped because the port would not converge: {skipped}");

    assert_eq!(
        skipped, 0,
        "{skipped} states failed to converge in the port. The Fortran would \
         have returned its guess instead (delta D-001), so these are not \
         disagreements, but they are also not evidence, and they need naming."
    );
    for comparison in &comparisons {
        comparison.assert_within(TIER2_TOLERANCE);
    }
}

/// With the transcendental taken out of the comparison, the algebra must be
/// bit-identical, exactly as it was through B-0017.
///
/// Run by `cargo xtask ci` as `--features oracle,libm-system`. It also carries
/// the premise the whole `libm-system` argument rests on: that gfortran's
/// `DEXP` and the platform `exp` are the same function. If they ever part, this
/// test is where it shows up, so its failure message must not claim the cause
/// is the algebra.
#[test]
#[cfg(feature = "libm-system")]
fn the_algebra_is_bit_identical_once_exp_agrees() {
    let mut oracle = Oracle::lock();
    let (comparisons, _) = sweep(&mut oracle, 200, 500, 1);

    let mut wrong = Vec::new();
    for comparison in &comparisons {
        println!("{comparison}");
        if comparison.max_ulp() != 0 {
            wrong.push(format!("{comparison}"));
        }
    }
    assert!(
        wrong.is_empty(),
        "{} field(s) are not bit-identical under the platform libm.\n\
         Two things can cause this, and they need different responses.\n\
         Either the algebra in `tepsim_core::equilibrium` does not associate \
         the way `teprob.f:473-502` does, which is a port bug; or gfortran no \
         longer resolves `DEXP` to the same code as the platform `exp`, in \
         which case this configuration no longer removes the transcendental \
         from the comparison and the claim it supports has to be restated. \
         `the_vendored_and_platform_exp_differ_only_by_rounding` distinguishes \
         them.\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}

/// Measure the price of the vendored `libm`, over the exact argument range
/// `teprob.f:485` and `teprob.f:488` reach.
///
/// This is the number every later Phase 2 item inherits, so it is measured on
/// whatever machine is running rather than quoted from a log entry. It
/// compares the vendored crate against the platform libm, and the sibling test
/// `the_algebra_is_bit_identical_once_exp_agrees` is what establishes that the
/// platform libm is what gfortran calls; together they pin the difference
/// against the Fortran.
///
/// It asserts only that the disagreement stays at one ULP. A larger
/// disagreement would not be wrong exactly, but it would invalidate the error
/// budget every downstream item is written against, and it should be a
/// decision rather than a surprise.
#[test]
fn the_vendored_and_platform_exp_differ_only_by_rounding() {
    use tepsim_core::constants::{AVP, BVP, CVP};
    use tepsim_core::{Component, math};

    let condensible = [
        Component::D,
        Component::E,
        Component::F,
        Component::G,
        Component::H,
    ];
    let (mut cases, mut differing, mut worst) = (0u64, 0u64, 0i64);
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);

    // Every tenth of a degree over a range far wider than the plant occupies.
    for tenths in -500..2500 {
        let celsius = f64::from(tenths) / 10.0;
        for component in condensible {
            let argument = AVP[component] + BVP[component] / (celsius + CVP[component]);
            if !argument.is_finite() {
                continue;
            }
            lo = lo.min(argument);
            hi = hi.max(argument);
            cases += 1;
            let delta = (math::exp(argument).to_bits() as i64) - (argument.exp().to_bits() as i64);
            if delta != 0 {
                differing += 1;
            }
            if delta.abs() > worst.abs() {
                worst = delta;
            }
        }
    }

    println!(
        "exp over the Antoine range [{lo:.4}, {hi:.4}]: {cases} arguments, \
         {differing} differ ({:.3}%), worst {worst} ulp",
        100.0 * differing as f64 / cases as f64
    );
    assert!(cases > 10_000, "the sweep collapsed to {cases} arguments");
    if math::USES_SYSTEM_LIBM {
        assert_eq!(differing, 0, "this build *is* the platform libm");
    } else {
        assert!(
            differing > 0,
            "the vendored libm agreed with the platform everywhere, so the \
             `libm-system` feature is measuring nothing and the bit-exactness \
             claim it supports is vacuous on this machine"
        );
        assert_eq!(
            worst.abs(),
            1,
            "the vendored libm is {worst} ulp from the platform, not 1. Every \
             Phase 2 error budget is written against a one-ulp difference."
        );
    }
}

/// The parts of this block that never touch `exp` must stay bit-identical
/// under the *default* build too.
///
/// Two of them: the A/B/C partial pressures, and the mixing zone pressure
/// `PTV`, whose vessel is single-phase and so has no Antoine term anywhere in
/// its derivation. `PTV` matters out of proportion to its one line, because it
/// drives the flow into the reactor and is measurement 13.
///
/// This keeps the ideal-gas half of the block held to zero ULP for the rest of
/// Phase 2 without needing the feature flag, and localises any future
/// regression to one side of the equilibrium immediately.
#[test]
#[cfg(not(feature = "libm-system"))]
fn the_ideal_gas_partial_pressures_are_bit_identical_under_the_vendored_libm() {
    use tepsim_core::Component;

    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 200, DT);
    let mut comparison: Comparison<Case> = Comparison::new("PPR/PPS for A, B and C");
    let mut mixing: Comparison<Case> = Comparison::new("PTV");

    for index in 0..pools.trajectory.len() {
        let scenario = pools.nominal_case(index);
        let snapshot = scenario.force(&mut oracle);
        let seeds = vessels::TemperatureSeeds {
            reactor: scenario.common.tcr,
            separator: scenario.common.tcs,
            stripper: scenario.common.tcc,
            mixing: scenario.common.tcv,
        };
        let state = State::from_flat(&scenario.state);
        let Ok(unpacked) = vessels::unpack(&state, seeds) else {
            continue;
        };
        let solved = equilibrium::equilibrium(&unpacked);
        mixing.observe(
            Case {
                pool: Pool::Nominal,
                index,
                component: 1,
            },
            solved.mixing_pressure,
            snapshot.common.ptv,
        );
        for component in [Component::A, Component::B, Component::C] {
            let case = Case {
                pool: Pool::Nominal,
                index,
                component: component.fortran_index(),
            };
            comparison.observe(
                case,
                solved.reactor.partial[component],
                snapshot.common.ppr[component.index()],
            );
            comparison.observe(
                case,
                solved.separator.partial[component],
                snapshot.common.pps[component.index()],
            );
        }
    }

    println!("{comparison}");
    println!("{mixing}");
    assert_eq!(
        comparison.max_ulp(),
        0,
        "the ideal-gas partial pressures involve no transcendental, so \
         nothing excuses a difference here"
    );
    assert_eq!(
        mixing.max_ulp(),
        0,
        "the mixing zone is single-phase and its pressure has no Antoine term \
         anywhere upstream, so nothing excuses a difference here either"
    );
}
