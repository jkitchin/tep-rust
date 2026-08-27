//! Builds the adversarial state catalogue and checks every entry landed.
//!
//! A state built for a boundary that misses it is worse than a missing state,
//! because it appears in the report as coverage while exercising the ordinary
//! branch. So every constructed state is verified against the observable it was
//! aimed at, and any boundary the bisection could not reach is named rather
//! than quietly dropped.

#![cfg(feature = "oracle")]

use tepsim_oracle::Oracle;
use tepsim_oracle::tier2::{Pools, adversarial};

const DT: f64 = 1.0 / 3600.0;

/// How close a constructed state has to sit to its boundary.
///
/// Sixty bisections exhaust an `f64` bracket, so the limit is not the search:
/// it is that some observables are step functions of the knob near the
/// boundary. One part in 1e-6 is far tighter than any branch is wide.
const PLACEMENT_TOLERANCE: f64 = 1e-6;

#[test]
fn every_boundary_in_the_catalogue_is_reached_and_verified() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 20, DT);
    let base = pools.nominal_case(0);

    let (built, missed) = adversarial::build(&mut oracle, &base);

    println!(
        "{} of {} boundaries constructed",
        built.len(),
        built.len() + missed.len()
    );
    for boundary in &built {
        println!(
            "  {:<52} knob {:<10.6} reached {:>14.6}  miss {:.1e}{}",
            boundary.target.name,
            boundary.setting,
            boundary.reached,
            boundary.miss(),
            if boundary.tripped { "  TRIPPED" } else { "" }
        );
    }
    for name in &missed {
        println!("  NOT REACHED: {name}");
    }

    for boundary in &built {
        boundary.verify(PLACEMENT_TOLERANCE);
    }
    assert!(
        missed.is_empty(),
        "{} boundaries could not be constructed: {missed:?}. Either the \
         bracket does not straddle the target or the observable is not \
         monotone in the knob. A named gap is acceptable only once it is \
         understood; an unnamed one is not.",
        missed.len()
    );
}

/// The catalogue has to cover the branches `PLAN.org` enumerates, and the count
/// is asserted so that deleting an entry is a test failure rather than a
/// quieter report.
#[test]
fn the_catalogue_covers_the_enumerated_branches() {
    let entries = adversarial::catalogue();
    // 17 from B-0016, plus three from B-0021 once `CPPRMX` existed to compare
    // against: the two pressure-ratio boundaries, and one state inside the
    // clamped region, because the boundary state itself does not clamp.
    assert_eq!(entries.len(), 20, "the catalogue changed size");

    let names: Vec<&str> = entries.iter().map(|(_, t, _)| t.name).collect();
    for expected in [
        "VLR at the 10% heat-transfer breakpoint",
        "VLR at the 50% heat-transfer breakpoint",
        "TCC at the lower stripping-factor branch",
        "TCC at the upper stripping-factor branch",
        "FTM(11) at the stripping-factor threshold",
        "PR = 1, the compressor reverse-flow clamp",
        "PR = CPPRMX, the compressor maximum-ratio clamp",
        "PR above CPPRMX, inside the clamped region",
    ] {
        assert!(names.contains(&expected), "{expected} is missing");
    }

    // Every entry must cite the line it sits on, or the catalogue decays into
    // a list of numbers nobody can check.
    for (_, target, _) in &entries {
        assert!(
            target.why.contains("teprob.f:"),
            "{} does not name its source line",
            target.name
        );
    }
}

/// Several boundaries trip the shutdown by construction. That is correct and
/// wanted, and it has to be *reported*, because a tripped plant has all fifty
/// derivatives zeroed (`teprob.f:807-811`) and any port reproduces that.
#[test]
fn the_shutdown_boundaries_actually_trip_and_are_counted() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 20, DT);
    let base = pools.nominal_case(0);
    let (built, _) = adversarial::build(&mut oracle, &base);

    let tripped: Vec<&str> = built
        .iter()
        .filter(|b| b.tripped)
        .map(|b| b.target.name)
        .collect();
    println!(
        "{} of {} boundaries trip: {tripped:?}",
        tripped.len(),
        built.len()
    );

    assert!(
        !tripped.is_empty(),
        "no boundary trips the shutdown, so the states built for the eight \
         shutdown conditions are not actually reaching them"
    );
    assert!(
        tripped.len() < built.len(),
        "every boundary trips, so the catalogue contains no state whose \
         derivative is evidence of anything"
    );
}

/// The four `DLP` clamps are reached, which the model never does on its own.
///
/// `teprob.f:576-598` clamps four pressure differences at zero. Each is a kink
/// in the flow as a function of the state, and the nominal plant sits well away
/// from all of them, so Tier 2 would never visit one without this pool.
#[test]
fn the_pressure_clamps_are_reached_from_the_correct_side() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 20, DT);
    let base = pools.nominal_case(0);

    let nominal = base.force(&mut oracle);
    println!(
        "nominal PTV {:.1}, PTR {:.1}, PTS {:.1}",
        nominal.common.ptv, nominal.common.ptr, nominal.common.pts
    );
    assert!(
        nominal.common.ptv > nominal.common.ptr
            && nominal.common.ptr > nominal.common.pts
            && nominal.common.pts > 760.0,
        "the nominal plant is already at one of the clamps, so placing a state \
         there proves nothing"
    );

    let (built, _) = adversarial::build(&mut oracle, &base);
    let clamps: Vec<&str> = built
        .iter()
        .filter(|b| b.target.why.contains("DLP is clamped"))
        .map(|b| b.target.name)
        .collect();
    println!("clamp boundaries constructed: {clamps:?}");
    assert_eq!(
        clamps.len(),
        3,
        "only {} of the three reachable DLP clamps were placed. The fourth, \
         the purge clamp, is unreachable; see the test below.",
        clamps.len()
    );
}

/// The purge-flow clamp cannot be reached from a nominal separator
/// composition, and this measures how far away it is.
///
/// `teprob.f:585-586` clamps `DLP = PTS - 760` at zero. Driving `PTS` down as
/// far as the state allows, by removing the whole separator vapour inventory
/// and cooling it to 0 C, leaves it at 811 mmHg: the separator liquid's own
/// vapour pressure. The ideal-gas term floors at `TKS = 273.15` rather than at
/// zero, so there is nothing further to remove.
///
/// The number is pinned because "we could not build that state" decays into
/// "nobody remembered to". If a future change to the pool or the constants
/// brings the floor below 760, this fails and the boundary joins the
/// catalogue.
#[test]
fn the_purge_clamp_is_unreachable_from_a_nominal_composition() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 20, DT);
    let base = pools.nominal_case(0);

    let mut floored = base.clone();
    (adversarial::Knob::SEPARATOR_DEPRESSURISE.apply)(&mut floored, 1e-9);
    let snapshot = floored.force(&mut oracle);

    let ideal: f64 = snapshot.common.pps[0..3].iter().sum();
    let antoine: f64 = snapshot.common.pps[3..8].iter().sum();
    println!(
        "PTS floor {:.2} mmHg at TCS {:.4} C: ideal-gas {ideal:.2}, Antoine {antoine:.2}",
        snapshot.common.pts, snapshot.common.tcs
    );

    assert!(
        snapshot.common.pts > 760.0,
        "PTS reached {:.2} mmHg, which is below the 760 clamp, so the purge \
         boundary IS reachable and belongs in the catalogue",
        snapshot.common.pts
    );
    assert!(
        (snapshot.common.pts - 811.49).abs() < 1.0,
        "the floor moved from the recorded 811.49 mmHg to {:.2}; something \
         about the separator changed",
        snapshot.common.pts
    );
    assert!(
        ideal < 1.0,
        "the ideal-gas contribution is {ideal:.2}, so the vapour inventory was \
         not actually removed and the floor is not a floor"
    );
}

/// A constructed boundary is still a valid, reproducible scenario.
#[test]
fn adversarial_scenarios_are_reproducible_like_any_other() {
    let mut oracle = Oracle::lock();
    let pools = Pools::collect(&mut oracle, 20, DT);
    let base = pools.nominal_case(0);
    let (built, _) = adversarial::build(&mut oracle, &base);

    for boundary in built.iter().take(5) {
        tepsim_oracle::tier2::reproducible(&mut oracle, &boundary.scenario);
    }
}

/// The knobs between them move most of the state vector, so the pool is not
/// concentrated in one corner.
#[test]
fn the_knobs_reach_most_of_the_state_vector() {
    let touched = adversarial::touched_slots();
    println!(
        "{} of 50 state slots are reachable by a knob",
        touched.len()
    );
    assert!(
        touched.len() >= 30,
        "only {} slots are reachable, so the adversarial pool leaves most of \
         the state at its nominal value",
        touched.len()
    );
}
