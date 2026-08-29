//! What the nominal plant is doing, in physical rates, over 48 hours.
//!
//! Run under [`Prices::unit`], so each term reads as the rate itself: kg/h for
//! the two material streams, kWh/h for the compressor, kg/h for steam. That is
//! the part of this crate the repository can validate, so it is the part worth
//! recording as a baseline.

use tepsim::{Scenario, Simulation};
use tepsim_operations::economics::{Prices, cost_rate, mean_cost_rate};

fn main() {
    let scenario = Scenario::baseline().with_hours(48.0);
    let run = Simulation::new(scenario).run();
    let (mean, counted) = mean_cost_rate(&run, &Prices::unit()).expect("rates");

    println!(
        "nominal, 48 h, {} samples, {counted} with analyser reads",
        run.samples.len()
    );
    println!("  purge       {:>12.3} kg/h", mean.purge);
    println!("  product     {:>12.3} kg/h", mean.product);
    println!("  compressor  {:>12.3} kWh/h", mean.compressor);
    println!("  steam       {:>12.3} kg/h", mean.steam);

    let last = run.samples.last().expect("a sample");
    println!("\nraw channels at the last sample:");
    println!(
        "  XMEAS(10) purge      {:>10.4} kscmh",
        last.measurements[9]
    );
    println!(
        "  XMEAS(17) product    {:>10.4} m3/h",
        last.measurements[16]
    );
    println!(
        "  XMEAS(19) steam      {:>10.4} kg/h",
        last.measurements[18]
    );
    println!("  XMEAS(20) compressor {:>10.4} kW", last.measurements[19]);

    // A worked example of the thing this crate is actually for: the same run,
    // costed under a price vector, and the reason it is expensive.
    let mut prices = Prices::free();
    prices.purge = [1.0, 0.1, 1.0, 2.0, 2.0, 0.0, 0.0, 0.0];
    prices.product = [2.0, 2.0, 2.0, 2.0, 2.0, 0.0, -3.0, -3.0];
    prices.compressor = 0.06;
    prices.steam = 0.02;
    let (costed, _) = mean_cost_rate(&run, &prices).expect("rates");
    println!("\nunder an ILLUSTRATIVE price vector (not Downs and Vogel's):");
    for (name, value) in costed.ranked() {
        println!("  {name:<11} {value:>12.2} /h");
    }
    println!("  {:<11} {:>12.2} /h", "total", costed.total());

    let sample = cost_rate(&last.measurements, &prices).expect("a rate");
    println!(
        "  instantaneous at the last sample: {:.2} /h",
        sample.total()
    );
}
