//! The Tier 1 harness measures the port, so something has to measure the
//! harness. These tests need no Fortran toolchain: they check that the pools
//! generate what they claim, that the ULP metric is exact on cases worked out
//! by hand, and that the gate fails when it should.
//!
//! The last part matters most. A comparison helper that reports "0 ULP, all
//! good" for a broken port is worse than no helper, so every failure mode it is
//! supposed to catch is provoked here on purpose.

use tepsim_core::{Component, Composition};
use tepsim_oracle::tier1::{
    BREAKPOINTS, Case, Comparison, EXACT_BUCKETS, Pool, Sampler, SimplexGrid, Sweep,
    TemperatureRange, relative_error, ulp_distance,
};

/// Cheap stand-in case label for the tests that do not care about compositions.
#[derive(Clone, Copy, Debug)]
struct Label(u32);

impl std::fmt::Display for Label {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "label#{}", self.0)
    }
}

fn smoke_case() -> Case {
    Case {
        pool: Pool::Grid,
        index: 0,
        composition: Composition::new([0.125; Component::COUNT]),
        celsius: 120.4,
    }
}

// ---------------------------------------------------------------- simplex grid

#[test]
fn the_grid_has_exactly_the_number_of_points_the_formula_predicts() {
    // C(divisions + 7, 7). Small cases are checked against an independent
    // count, so a wrong closed form cannot agree with a wrong enumeration.
    for divisions in 1_u32..=6 {
        let counted = SimplexGrid::new(divisions).count();
        assert_eq!(
            counted,
            SimplexGrid::len(divisions),
            "grid of {divisions} divisions"
        );
    }
    assert_eq!(
        SimplexGrid::len(1),
        8,
        "one division reaches the 8 vertices"
    );
    assert_eq!(SimplexGrid::len(8), 6435, "C(15,7)");
    assert_eq!(SimplexGrid::len(16), 245_157, "C(23,7)");
}

#[test]
fn every_grid_point_is_distinct() {
    let points: Vec<[u64; Component::COUNT]> = SimplexGrid::new(5)
        .map(|c| c.fractions().as_array().map(f64::to_bits))
        .collect();
    let mut sorted = points.clone();
    sorted.sort_unstable();
    sorted.dedup();
    assert_eq!(sorted.len(), points.len(), "the odometer repeated a point");
}

#[test]
fn a_power_of_two_grid_sums_to_exactly_one() {
    for point in SimplexGrid::new(8) {
        assert_eq!(
            point.sum().to_bits(),
            1.0_f64.to_bits(),
            "dyadic fractions must sum without rounding: {point:?}"
        );
    }
}

#[test]
fn the_grid_reaches_every_vertex_and_stays_in_the_simplex() {
    let mut vertices = 0;
    for point in SimplexGrid::new(4) {
        assert!(point.sums_to_one(), "left the simplex: {point:?}");
        for component in Component::ALL {
            assert!(point[component] >= 0.0, "negative fraction: {point:?}");
        }
        if point.fractions().iter().any(|f| *f > 0.999) {
            vertices += 1;
        }
    }
    assert_eq!(vertices, Component::COUNT, "one pure-species point each");
}

// -------------------------------------------------------------------- sampler

#[test]
fn the_sampler_is_reproducible_and_platform_independent() {
    // Integer arithmetic only, so these words are the same on every target.
    // If this test ever fails, every recorded Tier 1 number is invalidated,
    // because the sweep no longer visits the same cases.
    let mut sampler = Sampler::new(0x050E_15EE);
    let drawn: Vec<u64> = (0..4).map(|_| sampler.next_u64()).collect();
    let mut again = Sampler::new(0x050E_15EE);
    let repeat: Vec<u64> = (0..4).map(|_| again.next_u64()).collect();
    assert_eq!(drawn, repeat, "same seed, same stream");

    let mut different = Sampler::new(0x050E_15EF);
    assert_ne!(
        different.next_u64(),
        drawn[0],
        "adjacent seeds must not collide"
    );
}

#[test]
fn unit_draws_stay_in_the_half_open_interval() {
    let mut sampler = Sampler::new(1);
    let mut smallest = f64::INFINITY;
    let mut largest = f64::NEG_INFINITY;
    for _ in 0..100_000 {
        let u = sampler.unit();
        assert!((0.0..1.0).contains(&u), "unit draw out of range: {u}");
        smallest = smallest.min(u);
        largest = largest.max(u);
    }
    assert!(smallest < 0.001, "never approached zero: {smallest}");
    assert!(largest > 0.999, "never approached one: {largest}");
}

#[test]
fn dirichlet_samples_are_valid_compositions() {
    let mut sampler = Sampler::new(7);
    for _ in 0..20_000 {
        let c = sampler.dirichlet();
        assert!(c.sums_to_one(), "not normalised: {c:?} sums to {}", c.sum());
        for component in Component::ALL {
            assert!(c[component] >= 0.0, "negative fraction: {c:?}");
        }
    }
}

#[test]
fn dirichlet_is_uniform_on_the_simplex() {
    // Under Dirichlet(1,...,1) every marginal is Beta(1,7), so E[x] = 1/8 and
    // P(x > t) = (1-t)^7. Both are checked, because matching the mean alone
    // would also admit a fixed point at 1/8.
    const SAMPLES: usize = 200_000;
    let mut sampler = Sampler::new(11);
    let mut total = 0.0_f64;
    let mut above_quarter = 0_u32;
    for _ in 0..SAMPLES {
        let c = sampler.dirichlet();
        total += c[Component::A];
        if c[Component::A] > 0.25 {
            above_quarter += 1;
        }
    }
    let mean = total / SAMPLES as f64;
    assert!(
        (mean - 0.125).abs() < 0.002,
        "marginal mean {mean} is not 1/8"
    );
    let tail = f64::from(above_quarter) / SAMPLES as f64;
    let expected = 0.75_f64.powi(7);
    assert!(
        (tail - expected).abs() < 0.005,
        "P(x > 1/4) was {tail}, expected {expected}"
    );
}

#[test]
fn face_samples_land_on_faces_and_reach_every_support_size() {
    let mut sampler = Sampler::new(13);
    let mut seen_support = [0_u32; Component::COUNT + 1];
    for _ in 0..50_000 {
        let c = sampler.simplex_face();
        assert!(c.sums_to_one(), "not normalised: {c:?}");
        let support = c.fractions().iter().filter(|f| **f > 0.0).count();
        assert!(support >= 1, "empty support: {c:?}");
        for component in Component::ALL {
            assert!(c[component] >= 0.0, "negative fraction: {c:?}");
        }
        seen_support[support] += 1;
    }
    for (support, count) in seen_support.iter().enumerate().skip(1) {
        assert!(
            *count > 0,
            "no sample with support {support}; the boundary pool is not \
             covering what it claims"
        );
    }
    assert!(
        seen_support[1] > 0,
        "the vertices must be reachable exactly"
    );
}

// --------------------------------------------------------------- temperatures

#[test]
fn the_temperature_ladder_spans_the_range_and_includes_every_breakpoint() {
    let sweep = Sweep::SMOKE;
    let temperatures = sweep.temperatures();

    assert_eq!(
        temperatures[0].to_bits(),
        sweep.range.low.to_bits(),
        "the low bound must be hit exactly"
    );
    assert_eq!(
        temperatures[temperatures.len() - 1].to_bits(),
        sweep.range.high.to_bits(),
        "the high bound must be hit exactly, not approached"
    );

    for breakpoint in BREAKPOINTS {
        assert!(
            temperatures
                .iter()
                .any(|t| t.to_bits() == breakpoint.celsius.to_bits()),
            "{} is missing from the ladder, so nothing tests it: {}",
            breakpoint.celsius,
            breakpoint.why
        );
    }

    let mut sorted = temperatures.clone();
    sorted.sort_by(f64::total_cmp);
    assert_eq!(sorted, temperatures, "the ladder must be ascending");
    sorted.dedup_by(|a, b| a.to_bits() == b.to_bits());
    assert_eq!(sorted.len(), temperatures.len(), "duplicated temperature");
}

/// The domain claims to be the band the plant can reach. `TESUB4` divides by
/// `AD + (BD + CD*T)*T` (`teprob.f:1500-1501`), which has a root, and a sweep
/// that crossed it would be reporting on a singularity rather than on the
/// port. This computes each root from the constants and checks the ceiling
/// clears the nearest one.
#[test]
fn the_sweep_ceiling_stays_below_the_liquid_density_singularity() {
    use tepsim_core::constants::{AD, BD, CD};

    let mut nearest = f64::INFINITY;
    for component in Component::ALL {
        let (a, b, c) = (AD[component], BD[component], CD[component]);
        if c.abs() > 0.0 {
            // Positive root of c*T^2 + b*T + a = 0.
            let discriminant = (4.0 * c).mul_add(-a, b * b);
            assert!(discriminant > 0.0, "{component:?} has no real root");
            for root in [
                (-b + discriminant.sqrt()) / (2.0 * c),
                (-b - discriminant.sqrt()) / (2.0 * c),
            ] {
                if root > 0.0 {
                    nearest = nearest.min(root);
                }
            }
        }
    }

    assert!(
        nearest > TemperatureRange::PLANT.high,
        "the density correlation goes singular at {nearest} C, inside the \
         sweep range"
    );
    assert!(
        nearest < 250.0,
        "sanity: the nearest singularity should be a few tens of degrees \
         above the shutdown limit, got {nearest}"
    );
}

// ---------------------------------------------------------------------- sweep

#[test]
fn a_sweep_yields_exactly_as_many_cases_as_it_promises() {
    let sweep = Sweep::SMOKE;
    assert_eq!(sweep.cases().count(), sweep.len());
    assert!(!sweep.is_empty());
}

/// `Sweep::FULL` documents itself as "about ten million cases", which is the
/// volume `PLAN.org` asks for at Tier 1. Pin the claim so an innocent-looking
/// edit to one pool cannot quietly halve the coverage.
#[test]
fn the_full_sweep_is_the_volume_the_plan_asks_for() {
    let cases = Sweep::FULL.len();
    println!("Sweep::FULL: {cases} cases");
    assert!(
        (9_000_000..=11_000_000).contains(&cases),
        "Sweep::FULL is {cases} cases, not about ten million"
    );
}

#[test]
fn a_sweep_is_reproducible_case_for_case() {
    let sweep = Sweep::SMOKE;
    let first: Vec<_> = sweep
        .cases()
        .map(|c| {
            (
                c.pool,
                c.index,
                c.z().map(f64::to_bits),
                c.celsius.to_bits(),
            )
        })
        .collect();
    let second: Vec<_> = sweep
        .cases()
        .map(|c| {
            (
                c.pool,
                c.index,
                c.z().map(f64::to_bits),
                c.celsius.to_bits(),
            )
        })
        .collect();
    assert_eq!(
        first, second,
        "the same sweep must enumerate the same cases"
    );
}

#[test]
fn every_pool_is_represented_and_every_case_is_in_range() {
    let sweep = Sweep::SMOKE;
    let mut counts = [0_usize; 3];
    for case in sweep.cases() {
        assert!(
            sweep.range.contains(case.celsius),
            "{case} left the temperature range"
        );
        assert!(case.composition.sums_to_one(), "{case} is not normalised");
        counts[match case.pool {
            Pool::Grid => 0,
            Pool::Dirichlet => 1,
            Pool::Face => 2,
        }] += 1;
    }
    assert_eq!(counts[1], sweep.dirichlet_samples);
    assert_eq!(counts[2], sweep.face_samples);
    assert_eq!(
        counts[0],
        SimplexGrid::len(sweep.grid_divisions) * sweep.temperatures().len()
    );
}

/// Resizing one random pool must not disturb the other's cases, or every number
/// ever recorded for it would have to be re-measured.
#[test]
fn the_two_random_pools_have_independent_streams() {
    let small = Sweep {
        dirichlet_samples: 5,
        ..Sweep::SMOKE
    };
    let large = Sweep {
        dirichlet_samples: 500,
        ..Sweep::SMOKE
    };
    let faces = |sweep: Sweep| -> Vec<[u64; Component::COUNT]> {
        sweep
            .cases()
            .filter(|c| c.pool == Pool::Face)
            .take(20)
            .map(|c| c.z().map(f64::to_bits))
            .collect()
    };
    assert_eq!(faces(small), faces(large));
}

// ----------------------------------------------------------------- ulp metric

#[test]
fn ulp_distance_counts_representable_values() {
    assert_eq!(ulp_distance(1.0, 1.0), 0);
    assert_eq!(ulp_distance(1.0, f64::from_bits(1.0_f64.to_bits() + 1)), 1);
    assert_eq!(ulp_distance(1.0, f64::from_bits(1.0_f64.to_bits() + 7)), 7);
    assert_eq!(ulp_distance(-1.0, -1.0), 0);
    assert_eq!(ulp_distance(0.0, -0.0), 0, "the two zeros are one number");

    // Across an exponent boundary: the largest f64 below 1.0 and 1.0 itself.
    assert_eq!(ulp_distance(1.0, f64::from_bits(1.0_f64.to_bits() - 1)), 1);
    // The smallest subnormal is one step from zero, from either zero.
    assert_eq!(ulp_distance(0.0, f64::from_bits(1)), 1);
    assert_eq!(ulp_distance(-0.0, f64::from_bits(1)), 1);
    // A distance spanning zero is enormous, which is the point of reporting it.
    assert!(ulp_distance(f64::MIN_POSITIVE, -f64::MIN_POSITIVE) > (1 << 53));
    // Symmetric.
    assert_eq!(ulp_distance(1.0, 2.0), ulp_distance(2.0, 1.0));
}

#[test]
fn relative_error_is_what_it_says() {
    assert_eq!(relative_error(1.0, 1.0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(relative_error(0.0, 0.0).to_bits(), 0.0_f64.to_bits());
    assert_eq!(relative_error(-0.0, 0.0).to_bits(), 0.0_f64.to_bits());
    assert!((relative_error(1.000_001, 1.0) - 1e-6).abs() < 1e-12);
    assert!((relative_error(-2.0, -1.0) - 1.0).abs() < 1e-15);
    assert!(
        relative_error(1e-300, 0.0).is_infinite(),
        "a non-zero result where the Fortran returns exactly zero is a real \
         failure and must not be softened into a small number"
    );
    assert!(relative_error(f64::NAN, 1.0).is_nan());
}

// ------------------------------------------------------------ the report gate

#[test]
fn an_exact_match_reports_all_zeros() {
    let mut comparison = Comparison::new("exact");
    for i in 0..1000 {
        let value = f64::from(i) * 0.001;
        comparison.observe(Label(i as u32), value, value);
    }
    assert_eq!(comparison.cases(), 1000);
    assert_eq!(comparison.max_ulp(), 0);
    assert_eq!(comparison.max_relative_error().to_bits(), 0.0_f64.to_bits());
    assert_eq!(comparison.ulp_percentile(1.0), Some(0));
    comparison.assert_within(0.0);
}

/// The metric must report the size of a discrepancy, not merely its existence.
/// A known relative offset goes in; the same number has to come out.
#[test]
fn a_known_relative_offset_is_reported_at_its_true_size() {
    let mut comparison = Comparison::new("offset");
    for i in 1..1000 {
        let reference = f64::from(i);
        comparison.observe(Label(i as u32), reference * (1.0 + 1e-9), reference);
    }
    let reported = comparison.max_relative_error();
    assert!(
        (reported - 1e-9).abs() < 1e-15,
        "injected 1e-9, harness reported {reported:e}"
    );
    assert!(comparison.max_ulp() > 0, "and it is many ULP, not zero");
}

#[test]
fn the_histogram_and_percentiles_describe_the_distribution() {
    let mut comparison = Comparison::new("mixed");
    // 90 exact, 9 one ULP apart, 1 three ULP apart.
    for i in 0..90 {
        comparison.observe(Label(i), 1.0, 1.0);
    }
    for i in 0..9 {
        comparison.observe(Label(i), f64::from_bits(1.0_f64.to_bits() + 1), 1.0);
    }
    comparison.observe(Label(99), f64::from_bits(1.0_f64.to_bits() + 3), 1.0);

    assert_eq!(comparison.ulp_percentile(0.50), Some(0));
    assert_eq!(comparison.ulp_percentile(0.90), Some(0));
    assert_eq!(comparison.ulp_percentile(0.99), Some(1));
    assert_eq!(comparison.ulp_percentile(1.0), Some(3));
    assert_eq!(comparison.max_ulp(), 3);

    let report = comparison.to_string();
    assert!(
        report.contains("0:90"),
        "histogram missing exact bucket:\n{report}"
    );
    assert!(
        report.contains("1:9"),
        "histogram missing 1-ULP bucket:\n{report}"
    );
    assert!(
        report.contains("3:1"),
        "histogram missing 3-ULP bucket:\n{report}"
    );
    assert!(report.contains("p99=1"), "percentiles missing:\n{report}");
}

#[test]
fn distances_past_the_last_bucket_still_report_their_exact_worst() {
    let mut comparison = Comparison::new("far");
    let far = f64::from_bits(1.0_f64.to_bits() + 4096);
    comparison.observe(Label(0), far, 1.0);
    assert_eq!(comparison.max_ulp(), 4096);
    let report = comparison.to_string();
    assert!(
        report.contains(&format!(">={EXACT_BUCKETS}:1")),
        "the overflow bucket must be visible:\n{report}"
    );
    assert!(
        report.contains("p100>=16"),
        "and a percentile inside it must not be reported as a small number:\n{report}"
    );
}

#[test]
#[should_panic(expected = "exceeded its tolerance")]
fn the_gate_fails_when_the_error_is_too_large() {
    let mut comparison = Comparison::new("too big");
    comparison.observe(smoke_case(), 1.0 + 1e-9, 1.0);
    comparison.assert_within(1e-13);
}

#[test]
#[should_panic(expected = "no cases were compared")]
fn the_gate_fails_when_nothing_was_compared() {
    let comparison: Comparison<Label> = Comparison::new("empty");
    comparison.assert_within(1e-13);
}

#[test]
#[should_panic(expected = "non-finite")]
fn the_gate_fails_on_a_nan_that_the_fortran_did_not_produce() {
    let mut comparison = Comparison::new("nan");
    comparison.observe(smoke_case(), f64::NAN, 1.0);
    comparison.assert_within(1.0);
}

#[test]
fn two_matching_infinities_are_not_a_mismatch() {
    let mut comparison = Comparison::new("inf");
    comparison.observe(smoke_case(), f64::INFINITY, f64::INFINITY);
    assert_eq!(comparison.non_finite(), (1, 0));
    comparison.assert_within(0.0);
}

/// The report is pasted into `LOG.org` by hand, so its shape is part of the
/// contract with the next session.
#[test]
fn the_report_names_the_worst_case_well_enough_to_find_it_again() {
    let mut comparison = Comparison::new("TESUB1 ity=0");
    let case = Case {
        pool: Pool::Face,
        index: 4211,
        composition: Composition::new([0.125; Component::COUNT]),
        celsius: 120.4,
    };
    comparison.observe(smoke_case(), 1.0, 1.0);
    comparison.observe(case, 1.0 + 1e-12, 1.0);

    let report = comparison.to_string();
    assert!(
        report.starts_with("tier1 TESUB1 ity=0\n"),
        "header:\n{report}"
    );
    assert!(report.contains("face#4211"), "pool and index:\n{report}");
    assert!(
        report.contains("T=120.4"),
        "temperature, round-trippable:\n{report}"
    );
    assert!(report.contains("cases          : 2"), "count:\n{report}");
    assert_eq!(report.lines().count(), 7, "seven lines, stable:\n{report}");
}
