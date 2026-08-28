//! The integrators, and the property the three-phase split exists to provide.
//!
//! Two things need proving. That each method has the order it claims, which is
//! a statement about arithmetic and is checked against a problem with a
//! closed-form solution. And that the discrete phases run *once per outer
//! step* whichever method is chosen, which is a statement about the port's
//! structure and is the entire reason `TEFUNC` was split into three.

#![allow(
    clippy::float_cmp,
    reason = "determinism means exact equality is the property under test"
)]
// The scalar harness transcribes the published tableaux. Rearranging a
// coefficient into `mul_add` would make it harder to check against the paper,
// which is the only thing that makes the harness trustworthy.
#![allow(
    clippy::suboptimal_flops,
    reason = "the tableaux are transcribed to be checkable against the source"
)]

use tepsim::{Integrator, Outcome, Scenario, Simulation};
use tepsim_core::{State, VectorSpace};

// ---------------------------------------------------------------------------
// Order of accuracy, against a closed-form solution
// ---------------------------------------------------------------------------

/// Integrate `y' = y` from `y(0) = 1` and compare against `exp(t)`.
///
/// The plant is not used here on purpose: a convergence-order test needs an
/// exact answer, and the TEP has none. The integrators' tableaux are exercised
/// through the same `VectorSpace` operations either way, so an error in a
/// coefficient shows up here.
fn exponential_error(method: Integrator, steps: usize) -> f64 {
    let h = 1.0 / steps as f64;
    // A `State` whose every slot holds the same value, so the scalar problem
    // rides on the real vector-space implementation.
    let mut y = State::from_flat(&[1.0; 50]);
    let mut t = 0.0;
    for _ in 0..steps {
        y = advance_scalar(method, &y, t, h);
        t += h;
    }
    let exact = core::f64::consts::E;
    (y.to_flat()[0] - exact).abs() / exact
}

/// One step of `y' = y`, using the method's own tableau.
///
/// Written out here rather than reached for through `Integrator::advance`,
/// which needs a `Plant`. The coefficients are the same ones; this is a
/// scalar harness for them.
fn advance_scalar(method: Integrator, y: &State, t: f64, h: f64) -> State {
    let f = |_t: f64, s: &State| *s;
    match method {
        Integrator::Euler => y.add_scaled(h, &f(t, y)),
        Integrator::Rk4 => {
            let k1 = f(t, y);
            let k2 = f(t + 0.5 * h, &y.add_scaled(0.5 * h, &k1));
            let k3 = f(t + 0.5 * h, &y.add_scaled(0.5 * h, &k2));
            let k4 = f(t + h, &y.add_scaled(h, &k3));
            let sixth = h / 6.0;
            y.add_scaled(sixth, &k1)
                .add_scaled(2.0 * sixth, &k2)
                .add_scaled(2.0 * sixth, &k3)
                .add_scaled(sixth, &k4)
        }
        Integrator::DormandPrince => {
            let k1 = f(t, y);
            let k2 = f(t + h / 5.0, &y.add_scaled(h / 5.0, &k1));
            let k3 = f(
                t + 3.0 * h / 10.0,
                &y.add_scaled(3.0 * h / 40.0, &k1)
                    .add_scaled(9.0 * h / 40.0, &k2),
            );
            let k4 = f(
                t + 4.0 * h / 5.0,
                &y.add_scaled(44.0 * h / 45.0, &k1)
                    .add_scaled(-56.0 * h / 15.0, &k2)
                    .add_scaled(32.0 * h / 9.0, &k3),
            );
            let k5 = f(
                t + 8.0 * h / 9.0,
                &y.add_scaled(19372.0 * h / 6561.0, &k1)
                    .add_scaled(-25360.0 * h / 2187.0, &k2)
                    .add_scaled(64448.0 * h / 6561.0, &k3)
                    .add_scaled(-212.0 * h / 729.0, &k4),
            );
            let k6 = f(
                t + h,
                &y.add_scaled(9017.0 * h / 3168.0, &k1)
                    .add_scaled(-355.0 * h / 33.0, &k2)
                    .add_scaled(46732.0 * h / 5247.0, &k3)
                    .add_scaled(49.0 * h / 176.0, &k4)
                    .add_scaled(-5103.0 * h / 18656.0, &k5),
            );
            y.add_scaled(35.0 * h / 384.0, &k1)
                .add_scaled(500.0 * h / 1113.0, &k3)
                .add_scaled(125.0 * h / 192.0, &k4)
                .add_scaled(-2187.0 * h / 6784.0, &k5)
                .add_scaled(11.0 * h / 84.0, &k6)
        }
    }
}

/// Halving the step divides the error by `2^order`.
///
/// This is what "fourth order" *means*, and it is the only test that can
/// distinguish a correct tableau from a plausible one: a wrong coefficient
/// usually still converges, just at a lower order.
///
/// # The step counts are chosen, not arbitrary
///
/// An order can only be measured where *truncation* error dominates *rounding*
/// error, and the higher the order the sooner rounding takes over. Measured on
/// `y' = y` over one unit of time:
///
/// | steps | euler   | rk4     | dopri5  |
/// |-------|---------|---------|---------|
/// | 5     | 8.46e-2 | 1.13e-5 | 6.21e-8 |
/// | 20    | 2.39e-2 | 5.00e-8 | 7.96e-11|
/// | 40    | 1.22e-2 | 3.19e-9 | 2.60e-12|
/// | 160   | 3.11e-3 | 1.27e-11| 4.57e-15|
/// | 320   | 1.56e-3 | 7.93e-13| 3.27e-15|
///
/// Dormand-Prince is at the `f64` floor by 160 steps: the last halving divides
/// its error by 1.4 rather than by 32, and measuring its order there would say
/// 0.5. Each method is therefore measured on the coarsest pair where its fine
/// error is still comfortably above the floor. Choosing one pair for all three
/// is what made the first version of this test report order 4.6 for a method
/// that is order 5.
#[test]
fn each_method_converges_at_the_order_it_claims() {
    // (method, coarse steps, fine steps), chosen from the table above.
    let cases = [
        (Integrator::Euler, 100, 200),
        (Integrator::Rk4, 20, 40),
        (Integrator::DormandPrince, 20, 40),
    ];
    for (method, coarse_steps, fine_steps) in cases {
        let coarse = exponential_error(method, coarse_steps);
        let fine = exponential_error(method, fine_steps);
        let observed = (coarse / fine).log2();
        println!(
            "{:<8} order {} claimed, {observed:.3} observed  ({coarse:.3e} at {coarse_steps} steps -> {fine:.3e} at {fine_steps})",
            method.name(),
            method.order()
        );
        assert!(
            fine > 1e3 * f64::EPSILON,
            "{} is at the rounding floor at {fine_steps} steps ({fine:.3e}), so \
             its order cannot be measured there",
            method.name()
        );
        assert!(
            (observed - f64::from(method.order())).abs() < 0.15,
            "{} claims order {} and converges at {observed:.3}",
            method.name(),
            method.order()
        );
    }
}

/// The higher-order methods are dramatically more accurate at the same step,
/// which is the reason to have them.
#[test]
fn the_higher_order_methods_are_more_accurate() {
    // Twenty steps, where all three are above the rounding floor and the
    // comparison is between methods rather than against machine epsilon.
    let euler = exponential_error(Integrator::Euler, 20);
    let rk4 = exponential_error(Integrator::Rk4, 20);
    let dopri = exponential_error(Integrator::DormandPrince, 20);
    println!("20 steps of y'=y: euler {euler:.3e}, rk4 {rk4:.3e}, dopri5 {dopri:.3e}");

    assert!(
        rk4 < euler / 1e5,
        "rk4 {rk4:.3e} is not much better than euler {euler:.3e}"
    );
    assert!(
        dopri < rk4 / 100.0,
        "dopri5 {dopri:.3e} is not much better than rk4 {rk4:.3e}"
    );
}

/// Dormand-Prince reports an error estimate and the others do not.
///
/// The estimate is the gap between the fifth and fourth order solutions, which
/// is what an adaptive driver would compare against a tolerance. The step size
/// is not yet adapted: a variable step changes *when the discrete phases run*,
/// and that is a decision about fidelity rather than about numerics.
#[test]
fn only_dormand_prince_reports_an_error_estimate() {
    use tepsim_core::{Inputs, Plant, SimTime, constants};

    let plant = Plant::new();
    let state = State::from_flat(&constants::NOMINAL_STATE);
    let inputs = Inputs {
        manipulated: core::array::from_fn(|i| constants::NOMINAL_STATE[38 + i]),
        disturbances: [0.0; 20],
    };
    let h = 1.0 / 3600.0;

    for method in [Integrator::Euler, Integrator::Rk4] {
        let step = method
            .advance(&plant, SimTime(0.0), &state, &inputs, h)
            .expect("converges");
        assert!(step.error_estimate.is_none(), "{}", method.name());
    }

    let step = Integrator::DormandPrince
        .advance(&plant, SimTime(0.0), &state, &inputs, h)
        .expect("converges");
    let estimate = step.error_estimate.expect("dopri5 estimates its error");

    // Measured as an adaptive driver would measure it: the mixed
    // absolute-and-relative norm
    //
    //     max_i |e_i| / (atol + rtol |y_i|)
    //
    // which is below one exactly when the step would be accepted. Two simpler
    // normalisations were tried first and both are wrong here. An absolute
    // bound is a statement about the plant's units, since the state's
    // components run from order one to order 1e5. Dividing by the step's own
    // change, `h y'`, is worse: at the nominal operating point most
    // derivatives are near zero by cancellation, which is the whole subject of
    // the Tier 2 decision of 2026-08-27, so the ratio reached 2e4 on a
    // component whose derivative is essentially nil.
    const ATOL: f64 = 1e-6;
    const RTOL: f64 = 1e-6;
    let flat = estimate.to_flat();
    let y = state.to_flat();
    let mut worst = (0.0_f64, 0usize);
    for (slot, (error, value)) in flat.iter().zip(y).enumerate() {
        let scaled = error.abs() / RTOL.mul_add(value.abs(), ATOL);
        if scaled > worst.0 {
            worst = (scaled, slot + 1);
        }
    }
    let largest = flat.into_iter().map(f64::abs).fold(0.0_f64, f64::max);
    println!(
        "dopri5 at h=1 s: largest estimate {largest:.3e} absolute; mixed norm \
         {:.4} at slot {} (a driver accepts the step below 1)",
        worst.0, worst.1
    );
    assert!(
        largest > 0.0 && largest.is_finite(),
        "the estimate is {largest}, so the embedded solution is not being computed"
    );
    assert!(
        worst.0 < 1.0,
        "the mixed error norm is {:.4} at slot {}, so an adaptive driver would \
         reject a one-second step on the nominal plant",
        worst.0,
        worst.1
    );
}

// ---------------------------------------------------------------------------
// The property the three-phase split exists for
// ---------------------------------------------------------------------------

/// The disturbance walks advance once per outer step, whatever the method.
///
/// **This is the point of the whole split.** `TEFUNC` advances the walks and
/// draws noise on every call, so a four-stage method built naively on it would
/// run the disturbances at four times their intended rate and consume four
/// times the generator stream. The port evaluates only the pure phase inside
/// the integrator, so the generator word after N steps must be identical
/// whichever method ran.
#[test]
fn the_walks_advance_once_per_step_whichever_method() {
    let hours = 0.2;
    let mut words = Vec::new();
    for method in [
        Integrator::Euler,
        Integrator::Rk4,
        Integrator::DormandPrince,
    ] {
        let scenario = Scenario::fault(8).with_hours(hours).with_integrator(method);
        let run = Simulation::new(scenario).run();
        // The analyser schedule and the noise both come off the generator, so
        // an identical *count* of draws shows up as an identical number of
        // samples and, more sharply, in the sampled compositions being drawn
        // at the same instants.
        words.push((method, run.samples.len()));
        assert_eq!(
            run.samples.len(),
            scenario.samples(),
            "{} produced a different number of samples",
            method.name()
        );
    }
    println!("sample counts by method: {words:?}");
    assert!(words.windows(2).all(|w| w[0].1 == w[1].1));
}

/// A stronger form: the generator is in the same state after the run.
///
/// If a multi-stage method drew noise per stage, the generator would have
/// advanced four or seven times as far and this would differ.
#[test]
fn the_generator_ends_in_the_same_state_whichever_method() {
    // An hour, not a few minutes. The gas analyser updates every 0.1 h, so a
    // shorter run sees zero updates and the test passes without measuring
    // anything. It did exactly that when first written.
    let hours = 1.0;
    let mut analyser_columns = Vec::new();
    for method in [
        Integrator::Euler,
        Integrator::Rk4,
        Integrator::DormandPrince,
    ] {
        let run = Simulation::new(
            Scenario::baseline()
                .with_hours(hours)
                .with_integrator(method)
                .sampling_every(30),
        )
        .run();
        // XMEAS(23) is a sampled analyser reading. Its *schedule* is driven by
        // the generator-independent clock, but its noise is not, so the
        // sequence of instants at which it changes is a fingerprint of how
        // many draws have happened.
        let changes = run
            .measurement(23)
            .windows(2)
            .filter(|w| w[0] != w[1])
            .count();
        analyser_columns.push((method.name(), changes));
    }
    println!("analyser update counts: {analyser_columns:?}");
    let first = analyser_columns[0].1;
    assert!(
        first > 5,
        "the analyser updated only {first} times, which is too few for this \
         test to distinguish one draw per step from four"
    );
    for (name, count) in &analyser_columns {
        assert_eq!(
            *count, first,
            "{name} updated the analyser {count} times against {first}"
        );
    }
}

// ---------------------------------------------------------------------------
// Fidelity
// ---------------------------------------------------------------------------

/// Euler is unchanged, and it is the only faithful method.
///
/// Everything the validation ladder proved is a claim about Euler. Adding
/// integrators must not have moved it by a single bit.
#[test]
fn euler_is_the_default_and_the_only_faithful_method() {
    assert_eq!(Scenario::baseline().integrator, Integrator::Euler);
    assert!(Integrator::Euler.is_faithful());
    assert!(!Integrator::Rk4.is_faithful());
    assert!(!Integrator::DormandPrince.is_faithful());

    // Explicitly asking for Euler is the same run as not asking.
    let implicit = Simulation::new(Scenario::baseline().with_hours(0.5)).run();
    let explicit = Simulation::new(
        Scenario::baseline()
            .with_hours(0.5)
            .with_integrator(Integrator::Euler),
    )
    .run();
    assert_eq!(implicit, explicit);
}

/// The methods differ on the plant, and the difference is *Euler's error*.
///
/// This is the most interesting number the integrators produced, and it is a
/// statement about the original rather than about the port.
///
/// Two independent high-order methods agree with each other to about 1.5e-6.
/// Both differ from Euler by about 1e-2. When two methods of different order
/// agree four orders of magnitude better with each other than either does with
/// a third, the third is the outlier, and here the third is the first-order
/// method with a one-second step.
///
/// Measured, worst relative difference over all 53 channels:
///
/// | horizon | euler vs rk4 | rk4 vs dopri5 | euler vs dopri5 |
/// |---------|--------------|---------------|-----------------|
/// | 0.25 h  | 8.579e-3     | 1.180e-6      | 8.578e-3        |
/// | 1 h     | 1.113e-2     | 1.496e-6      | 1.112e-2        |
/// | 4 h     | 1.176e-2     | 1.723e-6      | 1.176e-2        |
///
/// So the published Tennessee Eastman data carries roughly one percent of
/// integration error against an accurate solution of the same equations. That
/// is not a defect of this port: reproducing it is the point, and `Euler` is
/// the default for exactly that reason. It does mean that "the TEP" names a
/// particular discretisation and not only a set of differential equations.
///
/// Note also that the gap *saturates* rather than growing: 8.6e-3 at a quarter
/// hour and 1.18e-2 at four hours. The controllers hold the plant at its
/// setpoints, so the integration error does not accumulate the way it would in
/// an open-loop run.
#[test]
fn the_gap_between_the_methods_is_eulers_own_error() {
    let base = Scenario::baseline().with_hours(1.0);
    let euler = Simulation::new(base).run();
    let rk4 = Simulation::new(base.with_integrator(Integrator::Rk4)).run();
    let dopri = Simulation::new(base.with_integrator(Integrator::DormandPrince)).run();

    for run in [&euler, &rk4, &dopri] {
        assert_eq!(run.outcome, Outcome::Completed);
    }

    let worst = |a: &tepsim::Run, b: &tepsim::Run| {
        a.samples
            .iter()
            .zip(&b.samples)
            .flat_map(|(x, y)| x.row().into_iter().zip(y.row()))
            .filter(|(_, y)| *y != 0.0)
            .map(|(x, y)| (x - y).abs() / y.abs())
            .fold(0.0_f64, f64::max)
    };

    let euler_rk4 = worst(&euler, &rk4);
    let rk4_dopri = worst(&rk4, &dopri);
    let euler_dopri = worst(&euler, &dopri);
    println!(
        "1 h: euler-rk4 {euler_rk4:.3e}, rk4-dopri5 {rk4_dopri:.3e}, \
         euler-dopri5 {euler_dopri:.3e}"
    );

    // The two high-order methods agree with each other.
    assert!(
        rk4_dopri < 1e-4,
        "rk4 and dopri5 disagree by {rk4_dopri:.3e}, so they are not both \
         converging to the same solution and one of the tableaux is wrong"
    );
    // And they agree with each other far better than either does with Euler,
    // which is what identifies Euler as the outlier.
    assert!(
        rk4_dopri < euler_rk4 / 1000.0,
        "rk4-dopri5 {rk4_dopri:.3e} is not far below euler-rk4 {euler_rk4:.3e}, \
         so the gap cannot be attributed to Euler"
    );
    // Euler is about a percent out, and the run is still recognisably the
    // same plant.
    assert!(
        (1e-3..5e-2).contains(&euler_rk4),
        "Euler's integration error over an hour is {euler_rk4:.3e}, outside \
         the 1e-3 to 5e-2 band this test records"
    );
    // RK4 is not bit-identical to Euler, which would mean it never ran its
    // extra stages.
    assert!(euler_rk4 > 0.0);
}

#[test]
fn the_method_names_round_trip() {
    for method in [
        Integrator::Euler,
        Integrator::Rk4,
        Integrator::DormandPrince,
    ] {
        assert_eq!(Integrator::parse(method.name()), Some(method));
    }
    assert_eq!(
        Integrator::parse("dormand-prince"),
        Some(Integrator::DormandPrince)
    );
    assert_eq!(Integrator::parse("nonsense"), None);
    assert_eq!(Integrator::Rk4.stages(), 4);
    assert_eq!(Integrator::DormandPrince.stages(), 7);
}
