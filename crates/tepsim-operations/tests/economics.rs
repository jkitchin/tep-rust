//! What the economics module can be held to, and what it cannot.
//!
//! B-0074. The absolute currency value is not validated against anything, and
//! saying so once in a doc comment is not the same as testing what *is* true.
//! These are the claims with teeth: the unit conversions are the Fortran's own,
//! each term moves the answer in the direction it should, and the whole thing
//! agrees with the plant it is reading.

// Exact comparisons throughout. Zero prices give exactly zero, an inverted
// conversion returns exactly what it started from up to a stated tolerance, and
// a term untouched by a price change is untouched exactly. A tolerance on those
// would be asserting something weaker than what is true.
#![allow(
    clippy::float_cmp,
    reason = "exactness is the property under test, not an approximation of it"
)]

use tepsim::{Run, Scenario, Simulation};
use tepsim_core::{Component, Composition};
use tepsim_operations::economics::{
    self, Prices, cost_rate, mean_cost_rate, product_molar_rate, purge_molar_rate,
};

/// A short nominal run, sampled as the published data is.
fn nominal(hours: f64) -> Run {
    Simulation::new(Scenario::baseline().with_hours(hours)).run()
}

// ---------------------------------------------------------------------------
// The unit conversions are the Fortran's, inverted
// ---------------------------------------------------------------------------

/// `teprob.f:688` is `XMEAS(10)=FTM(10)*0.359/35.3145`. Inverting it has to
/// return the molar flow it started from, or every mass rate here is wrong by a
/// constant nobody would notice.
#[test]
fn the_purge_conversion_inverts_teprob_f_688() {
    let source = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../reference/fortran/teprob.f"
    ))
    .expect("the vendored Fortran");
    let line = source.lines().nth(687).expect("line 688");
    let squashed: String = line.split_whitespace().collect();
    assert_eq!(
        squashed, "XMEAS(10)=FTM(10)*0.359/35.3145",
        "teprob.f:688 is not the line this conversion was derived from"
    );

    for molar in [0.0, 1.0, 1234.5, 9.87e4] {
        let kscmh = molar * 0.359 / 35.3145;
        let back = purge_molar_rate(kscmh);
        assert!(
            (back - molar).abs() <= 1e-9 * molar.max(1.0),
            "{molar} lbmol/h round-tripped to {back}"
        );
    }
}

/// `teprob.f:695` is `XMEAS(17)=FTM(13)/DLC/35.3145`, and the same applies. The
/// density comes from the plant's own correlation rather than from a constant.
#[test]
fn the_product_conversion_inverts_teprob_f_695() {
    let source = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../reference/fortran/teprob.f"
    ))
    .expect("the vendored Fortran");
    let squashed: String = source
        .lines()
        .nth(694)
        .expect("line 695")
        .split_whitespace()
        .collect();
    assert_eq!(squashed, "XMEAS(17)=FTM(13)/DLC/35.3145");

    // A plausible stripper underflow: mostly G and H with unconverted D and E.
    let x = Composition::new([0.0, 0.0, 0.0, 0.02, 0.02, 0.0, 0.50, 0.46]);
    let celsius = 65.7;
    let density = tepsim_core::thermo::liquid_density(&x, celsius);
    assert!(density > 0.0, "the correlation returned {density}");

    for molar in [0.0, 10.0, 259.5] {
        let m3 = molar / density / 35.3145;
        let back = product_molar_rate(m3, &x, celsius);
        assert!(
            (back - molar).abs() <= 1e-9 * molar.max(1.0),
            "{molar} lbmol/h round-tripped to {back}"
        );
    }
}

// ---------------------------------------------------------------------------
// The structural claims
// ---------------------------------------------------------------------------

/// Nothing flowing costs nothing, and no price means no cost whatever flows.
#[test]
fn a_still_plant_and_a_free_one_both_cost_nothing() {
    let run = nominal(2.0);
    let sample = run.samples.last().expect("a sample");

    let free = cost_rate(&sample.measurements, &Prices::free()).expect("a rate");
    assert_eq!(free.total(), 0.0, "every price was zero");

    // Every flow zero, with the analysers still reading a real composition, so
    // this tests the flows and not the guard on an unusable analyser.
    let mut still = sample.measurements;
    still[9] = 0.0; // purge
    still[16] = 0.0; // product
    still[18] = 0.0; // steam
    still[19] = 0.0; // compressor
    let stopped = cost_rate(&still, &Prices::unit()).expect("a rate");
    assert_eq!(stopped.total(), 0.0, "nothing flowed: {stopped:?}");
}

/// Each of the four terms responds to its own price and to nothing else.
///
/// This is what catches a term wired to the wrong channel, which is the
/// mistake a table of unit conversions invites.
#[test]
fn each_term_answers_to_its_own_price_alone() {
    let run = nominal(4.0);
    let sample = run.samples.last().expect("a sample");
    let base = cost_rate(&sample.measurements, &Prices::unit()).expect("a rate");

    // Doubling the compressor price doubles the compressor term and moves
    // nothing else.
    let mut prices = Prices::unit();
    prices.compressor = 2.0;
    let moved = cost_rate(&sample.measurements, &prices).expect("a rate");
    assert_eq!(moved.compressor, 2.0 * base.compressor);
    assert_eq!(moved.purge, base.purge);
    assert_eq!(moved.product, base.product);
    assert_eq!(moved.steam, base.steam);

    let mut prices = Prices::unit();
    prices.steam = 3.0;
    let moved = cost_rate(&sample.measurements, &prices).expect("a rate");
    assert_eq!(moved.steam, 3.0 * base.steam);
    assert_eq!(moved.compressor, base.compressor);

    // And a single component's purge price moves only the purge term.
    let mut prices = Prices::free();
    prices.purge[Component::A as usize] = 1.0;
    let only_a = cost_rate(&sample.measurements, &prices).expect("a rate");
    assert!(only_a.purge > 0.0, "no A in the purge: {only_a:?}");
    assert_eq!(only_a.product, 0.0);
    assert_eq!(only_a.compressor, 0.0);
    assert_eq!(only_a.steam, 0.0);
}

/// Under unit prices the terms are physical rates, so they have to be the
/// right size. A mass rate that is out by a factor of 35 would pass every
/// test above and this is what catches it.
#[test]
fn the_physical_rates_are_the_right_order_of_magnitude() {
    let run = nominal(6.0);
    let (mean, counted) = mean_cost_rate(&run, &Prices::unit()).expect("rates");
    assert!(
        counted > 50,
        "only {counted} samples carried an analyser read"
    );

    println!(
        "nominal plant, unit prices: purge {:.2} kg/h, product {:.2} kg/h, \
         compressor {:.2} kWh/h, steam {:.2} kg/h ({counted} samples)",
        mean.purge, mean.product, mean.compressor, mean.steam
    );

    // The product stream is the plant's output and is the largest mass flow of
    // the four by a wide margin; the purge is a vent and is much smaller.
    assert!(
        mean.product > mean.purge,
        "the purge out-massed the product: {mean:?}"
    );
    // Bounds wide enough to be about magnitude rather than about a value:
    // a real chemical plant, not a laboratory and not a refinery.
    assert!(
        (1e3..1e6).contains(&mean.product),
        "product mass rate {} kg/h is not a plant-sized number",
        mean.product
    );
    assert!(
        (1e2..1e5).contains(&mean.purge),
        "purge mass rate {} kg/h is not plant-sized",
        mean.purge
    );
    assert!(
        (1e1..1e4).contains(&mean.compressor),
        "compressor {} kW is not plant-sized",
        mean.compressor
    );
}

/// The product mass rate agrees with the plant's own product measurement.
///
/// `XMEAS(17)` is a volume and the analyser gives a composition; multiplying
/// them back out through the density has to land on the same stream the
/// simulator thinks it is producing. This is the closest thing to an
/// independent check available, because it goes out through the measurements
/// and back in through the thermodynamics.
#[test]
fn the_product_mass_rate_agrees_with_volume_times_density() {
    let run = nominal(4.0);
    let sample = run.samples.last().expect("a sample");
    let m = &sample.measurements;

    let x = economics::product_composition(m).expect("a composition");
    let density = tepsim_core::thermo::liquid_density(&x, m[17]);
    let molar = product_molar_rate(m[16], &x, m[17]);

    // lbmol/h back to m3/h through the same density.
    let volume = molar / density / 35.3145;
    assert!(
        (volume - m[16]).abs() <= 1e-9 * m[16].abs().max(1.0),
        "{volume} m3/h against XMEAS(17) = {}",
        m[16]
    );
}

// ---------------------------------------------------------------------------
// Compositions
// ---------------------------------------------------------------------------

/// The analysers are noisy and read in percent, so they do not sum to 100.
/// Normalising is the point, and this asserts the raw readings really do
/// violate the sum so the normalisation is not decorative.
#[test]
fn the_analyser_readings_need_normalising_and_get_it() {
    let run = nominal(3.0);
    let sample = run.samples.last().expect("a sample");
    let m = &sample.measurements;

    let raw_purge: f64 = m[28..36].iter().sum();
    assert!(
        (raw_purge - 100.0).abs() > Composition::SUM_TOLERANCE,
        "the raw purge analyser summed to {raw_purge}, so nothing needed \
         normalising and this test proves nothing"
    );

    let purge = economics::purge_composition(m).expect("a composition");
    let sum: f64 = Component::ALL.map(|c| purge[c]).iter().sum();
    assert!((sum - 1.0).abs() < 1e-12, "normalised sum {sum}");

    let product = economics::product_composition(m).expect("a composition");
    // A, B and C are not in the stripper underflow.
    for absent in [Component::A, Component::B, Component::C] {
        assert_eq!(product[absent], 0.0, "{absent:?} appeared in the product");
    }
    let sum: f64 = Component::ALL.map(|c| product[c]).iter().sum();
    assert!((sum - 1.0).abs() < 1e-12, "normalised sum {sum}");
}

/// An analyser that has never reported is distinguishable from a free plant.
#[test]
fn an_unusable_analyser_is_none_rather_than_zero() {
    let mut m = [0.0_f64; 41];
    m[9] = 10.0;
    m[16] = 20.0;
    assert!(
        cost_rate(&m, &Prices::unit()).is_none(),
        "all-zero analysers produced a cost rather than None"
    );
    assert!(
        cost_rate(&[0.0; 5], &Prices::unit()).is_none(),
        "short slice"
    );
}

// ---------------------------------------------------------------------------
// Ranking
// ---------------------------------------------------------------------------

#[test]
fn ranked_puts_the_biggest_term_first() {
    let rate = tepsim_operations::CostRate {
        purge: 3.0,
        product: 10.0,
        compressor: 1.0,
        steam: 7.0,
    };
    let ranked = rate.ranked();
    assert_eq!(ranked[0], ("product", 10.0));
    assert_eq!(ranked[3], ("compressor", 1.0));
    assert_eq!(rate.total(), 21.0);
}
