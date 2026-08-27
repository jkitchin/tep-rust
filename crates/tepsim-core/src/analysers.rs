//! Measurement noise and the three sampled composition analysers.
//!
//! Ported from `teprob.f:711-761`. This is the whole of
//! [`crate::plant::Plant::sample_measurements`], the impure post-phase of the
//! three-phase split, and the last piece of the model.
//!
//! # Two kinds of measurement
//!
//! `XMEAS(1..22)` are continuous instruments. They are read every step, and
//! noise is added to the value the model just computed.
//!
//! `XMEAS(23..41)` are *analysers*. They are read on a schedule, 0.1 hours for
//! the two gas analysers and 0.25 for the product one, and each reports the
//! composition from its *previous* sample rather than the current one. That
//! delay is the dead time, and it is the reason composition control in this
//! plant is hard.
//!
//! # The dead time is a latch, and the order of two lines makes it
//!
//! ```fortran
//!       XMEAS(I)=XDEL(I)
//!       CALL TESUB6(XNS(I),XMNS)
//!       XMEAS(I)=XMEAS(I)+XMNS
//!       XDEL(I)=XCMP(I)
//! ```
//!
//! The reported value is taken from the store *before* the store is updated.
//! Swapping those two lines gives an analyser with no dead time at all, which
//! produces entirely plausible numbers and a plant that is much easier to
//! control than the real one.
//!
//! # Three guards, and each one is about *when*
//!
//! **Noise is skipped at `TIME = 0` and on a tripped plant**
//! (`teprob.f:711`). Only the continuous noise, though: the analyser blocks at
//! `744-761` have no such guard, so a tripped plant still draws. B-0027
//! measured 258 draws in a tripped evaluation against 522 in a healthy one,
//! and a port that silenced everything on a trip would leave the generator 264
//! steps behind and desynchronise every later draw.
//!
//! **The schedules advance from their own previous value**, not from the
//! current time: `TGAS = TGAS + 0.1` at `teprob.f:751`. So a step arriving
//! late does not shift the schedule, exactly as a walk segment hands over at
//! its own `TNEXT` rather than at `t` (see [`mod@crate::walk`]).
//!
//! Both `0.1` literals, at `teprob.f:741` and `751`, are **single precision**.
//! The gas interval is therefore 0.10000000149011612, not 0.1, and a step
//! landing on exactly 0.1 does *not* sample. `0.25` is exactly representable
//! so the product analyser is unaffected, which is precisely the kind of
//! inconsistency that has to be read off the line rather than inferred from
//! its neighbour.
//!
//! **At `TIME = 0` the analysers are primed rather than sampled**
//! (`teprob.f:736-743`): the store and the reported value are both set to the
//! current composition, with no noise and no draw, and the two schedules are
//! set to their first due times.
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `XCMP(23..41)` | [`compositions`] | what the analysers would read now |
//! | `XDEL(23..41)` | [`Analysers::stored`] | what they will report next time |
//! | `XMEAS(23..41)` | [`Analysers::reported`] | what they reported last time |
//! | `TGAS`, `TPROD` | [`Analysers::next_gas`] etc. | when each is next due |

// Every float expression here reproduces `teprob.f`'s association and
// rounding exactly; see `crate::thermo`.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]

use crate::component::Component;
use crate::constants::{MEASUREMENT_NOISE, single};
use crate::disturbance::{Draws, noise};
use crate::measurements::Shutdown;
use crate::plant::{N_CONTINUOUS, N_SAMPLED};
use crate::stream::Stream;
use crate::streams::Streams;

/// The gas analysers' sampling interval, hours (`teprob.f:741`, `751`).
const GAS_INTERVAL: f64 = single(0.1);
/// The product analyser's, hours (`teprob.f:742`, `760`).
const PRODUCT_INTERVAL: f64 = single(0.25);

/// Mole fraction to percent (`teprob.f:717`). Single precision.
const PERCENT: f64 = single(100.0);

/// How many compositions the two gas analysers report, `XMEAS(23..36)`.
pub const GAS_COMPOSITIONS: usize = 14;

/// The three analysers' persistent state.
///
/// Genuine state, like [`crate::TemperatureSeeds`] and the latched valve
/// commands: what an analyser reports depends on when it last sampled.
///
/// # There are *two* stored values per analyser, not one
///
/// `XDEL(i)` holds the composition the analyser will report next time it is
/// due, noise-free. `XMEAS(i)` holds what it reported *last* time, noise
/// included, and the original simply leaves it alone between samples
/// (`teprob.f:744` writes it only inside the schedule check).
///
/// They are not the same number: they differ by one noise draw. Carrying only
/// `XDEL` and reporting it between samples gives a reading that is the right
/// composition with the noise stripped off, which looks entirely plausible,
/// stays within a fraction of a percent, and is wrong. It cost a differential
/// failure at step 362 of a 2,000-step run to find.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Analysers {
    /// `XDEL(23..41)`: the composition each analyser will report next time it
    /// is due, noise-free. Indexed from `XMEAS(23)`.
    pub stored: [f64; N_SAMPLED],
    /// `XMEAS(23..41)`: what each analyser reported last time, noise included.
    ///
    /// Held unchanged between samples. See the note above.
    pub reported: [f64; N_SAMPLED],
    /// `TGAS`: when the two gas analysers are next due.
    pub next_gas: f64,
    /// `TPROD`: when the product analyser is next due.
    pub next_product: f64,
}

impl Default for Analysers {
    /// Before the first sample. `teprob.f:741-742` sets the two schedules at
    /// `TIME = 0`, and the store is filled at the same moment.
    fn default() -> Self {
        Self {
            stored: [0.0; N_SAMPLED],
            reported: [0.0; N_SAMPLED],
            next_gas: GAS_INTERVAL,
            next_product: PRODUCT_INTERVAL,
        }
    }
}

/// The nineteen compositions the analysers would read right now, `XCMP`.
///
/// Reactor feed from stream 7, purge from stream 10, product from stream 13.
/// Note the reactor feed reports only A through F and the product only D
/// through H: the analysers do not measure what cannot be there in quantity.
// @port teprob.f:717-735
#[must_use]
pub fn compositions(stream_table: &Streams) -> [f64; N_SAMPLED] {
    let mut out = [0.0; N_SAMPLED];
    let fraction = |stream: Stream, c: Component| stream_table.composition[stream][c] * PERCENT;

    // teprob.f:717-722. XMEAS(23..28), the reactor feed, A through F.
    for (slot, c) in Component::ALL[..6].iter().enumerate() {
        out[slot] = fraction(Stream::ReactorInlet, *c);
    }
    // teprob.f:723-730. XMEAS(29..36), the purge, all eight.
    for (slot, c) in Component::ALL.iter().enumerate() {
        out[6 + slot] = fraction(Stream::Purge, *c);
    }
    // teprob.f:731-735. XMEAS(37..41), the product, D through H.
    for (slot, c) in Component::ALL[3..].iter().enumerate() {
        out[14 + slot] = fraction(Stream::Product, *c);
    }
    out
}

/// Add noise to the continuous measurements and tick the three analysers.
///
/// Returns all forty-one, in `XMEAS` order.
// @port teprob.f:711-761
pub fn sample<R: Draws + ?Sized>(
    state: &mut Analysers,
    rng: &mut R,
    t: f64,
    continuous: &[f64; N_CONTINUOUS],
    compositions: &[f64; N_SAMPLED],
    shutdown: Shutdown,
) -> [f64; N_CONTINUOUS + N_SAMPLED] {
    let mut out = [0.0; N_CONTINUOUS + N_SAMPLED];
    out[..N_CONTINUOUS].copy_from_slice(continuous);

    // teprob.f:711-716. Skipped at t=0 and on a trip; see the module docs.
    if t > 0.0 && !shutdown.is_tripped() {
        for (index, value) in out[..N_CONTINUOUS].iter_mut().enumerate() {
            *value += noise(rng, MEASUREMENT_NOISE[index]);
        }
    }

    // teprob.f:736-743. Prime rather than sample: no noise, no draw.
    if t == 0.0 {
        // teprob.f:738-739. Both stores take the current composition.
        state.stored = *compositions;
        state.reported = *compositions;
        out[N_CONTINUOUS..].copy_from_slice(compositions);
        state.next_gas = GAS_INTERVAL;
        state.next_product = PRODUCT_INTERVAL;
        return out;
    }

    // Between samples the original leaves `XMEAS(23..41)` exactly as it was,
    // so the *reported* value persists, noise and all. Not `XDEL`; see the
    // note on `Analysers`.
    out[N_CONTINUOUS..].copy_from_slice(&state.reported);

    // teprob.f:744-752. The two gas analysers, XMEAS(23..36).
    if t >= state.next_gas {
        for index in 0..GAS_COMPOSITIONS {
            // The store first, then the noise, then the update. Swapping the
            // last two lines removes the dead time; see the module docs.
            let mut value = state.stored[index];
            value += noise(rng, MEASUREMENT_NOISE[N_CONTINUOUS + index]);
            out[N_CONTINUOUS + index] = value;
            state.reported[index] = value;
            state.stored[index] = compositions[index];
        }
        state.next_gas += GAS_INTERVAL;
    }

    // teprob.f:753-761. The product analyser, XMEAS(37..41).
    if t >= state.next_product {
        for index in GAS_COMPOSITIONS..N_SAMPLED {
            let mut value = state.stored[index];
            value += noise(rng, MEASUREMENT_NOISE[N_CONTINUOUS + index]);
            out[N_CONTINUOUS + index] = value;
            state.reported[index] = value;
            state.stored[index] = compositions[index];
        }
        state.next_product += PRODUCT_INTERVAL;
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::TepRng;
    use crate::testing::assert_exact;

    fn draws(before: f64, after: f64) -> usize {
        let mut probe = TepRng::new(before);
        for step in 0..2_000 {
            if probe.state().to_bits() == after.to_bits() {
                return step;
            }
            let _ = probe.unit();
        }
        panic!("more than two thousand draws");
    }

    fn comps() -> [f64; N_SAMPLED] {
        core::array::from_fn(|i| 1.0 + i as f64)
    }

    /// At time zero the analysers are primed: the store and the report are set
    /// to the current composition, and nothing is drawn.
    #[test]
    fn time_zero_primes_the_analysers_without_drawing() {
        let mut state = Analysers::default();
        let mut rng = TepRng::with_default_seed();
        let before = rng.state();
        let compositions = comps();

        let out = sample(
            &mut state,
            &mut rng,
            0.0,
            &[7.0; N_CONTINUOUS],
            &compositions,
            Shutdown::default(),
        );
        assert_eq!(draws(before, rng.state()), 0, "priming must not draw");
        for (a, b) in state.stored.iter().zip(compositions) {
            assert_eq!(a.to_bits(), b.to_bits(), "the store takes XCMP");
        }
        for (a, b) in out[N_CONTINUOUS..].iter().zip(compositions) {
            assert_eq!(a.to_bits(), b.to_bits(), "and so does the report");
        }
        // And the continuous measurements pass through unchanged.
        for value in &out[..N_CONTINUOUS] {
            assert_exact(*value, 7.0, "no noise at t=0");
        }
    }

    /// The reported composition is the *previous* sample. That is the dead
    /// time, and it is the single most consequential line in the module.
    #[test]
    fn an_analyser_reports_the_previous_sample_not_the_current_one() {
        let mut state = Analysers::default();
        let mut rng = TepRng::with_default_seed();
        let first = comps();
        // Prime at t=0.
        let _ = sample(
            &mut state,
            &mut rng,
            0.0,
            &[0.0; N_CONTINUOUS],
            &first,
            Shutdown::default(),
        );

        // A completely different composition, sampled once the gas analysers
        // are due.
        let second: [f64; N_SAMPLED] = core::array::from_fn(|i| 100.0 + i as f64);
        // Just past the schedule, which is *not* 0.1: `teprob.f:741` writes
        // `TGAS=0.1` with no `D` suffix, so it is 0.10000000149011612 and a
        // step landing on exactly 0.1 does not sample.
        let out = sample(
            &mut state,
            &mut rng,
            GAS_INTERVAL,
            &[0.0; N_CONTINUOUS],
            &second,
            Shutdown::default(),
        );
        for index in 0..GAS_COMPOSITIONS {
            let reported = out[N_CONTINUOUS + index];
            assert!(
                (reported - first[index]).abs() < 5.0,
                "slot {index} reported {reported}, which is the *current* \
                 composition {} rather than the previous {}: the store is \
                 being updated before it is read, and the plant has no dead \
                 time",
                second[index],
                first[index]
            );
            // And the store now holds the current one, ready for next time.
            assert_exact(state.stored[index], second[index], "store updated");
        }
    }

    /// A healthy step past both schedules draws for all forty-one
    /// measurements, twelve each.
    ///
    /// The trip case is checked against the oracle in
    /// `tier3_analysers.rs`, where a genuinely tripped `Shutdown` is available
    /// from the measurement layer rather than constructed by hand.
    #[test]
    fn a_full_sample_draws_twelve_times_for_every_measurement() {
        let mut state = Analysers::default();
        let mut rng = TepRng::with_default_seed();
        let compositions = comps();
        let _ = sample(
            &mut state,
            &mut rng,
            0.0,
            &[0.0; N_CONTINUOUS],
            &compositions,
            Shutdown::default(),
        );

        let before = rng.state();
        let _ = sample(
            &mut state,
            &mut rng,
            0.3,
            &[0.0; N_CONTINUOUS],
            &compositions,
            Shutdown::default(),
        );
        // 22 continuous + 14 gas + 5 product, twelve draws each.
        assert_eq!(draws(before, rng.state()), (22 + 14 + 5) * 12);
    }

    /// The schedules advance from their own previous value, so a late step
    /// does not shift them.
    #[test]
    fn the_schedules_do_not_drift_when_a_step_arrives_late() {
        let mut state = Analysers::default();
        let mut rng = TepRng::with_default_seed();
        let compositions = comps();
        let _ = sample(
            &mut state,
            &mut rng,
            0.0,
            &[0.0; N_CONTINUOUS],
            &compositions,
            Shutdown::default(),
        );
        assert_exact(state.next_gas, GAS_INTERVAL, "first gas sample");

        // A step that arrives well after the due time.
        let _ = sample(
            &mut state,
            &mut rng,
            0.19,
            &[0.0; N_CONTINUOUS],
            &compositions,
            Shutdown::default(),
        );
        assert_exact(
            state.next_gas,
            GAS_INTERVAL + GAS_INTERVAL,
            "the schedule must advance by its own interval, not from t",
        );
    }

    /// Between samples an analyser holds its reading rather than reporting
    /// zero.
    #[test]
    fn an_analyser_holds_its_reading_between_samples() {
        let mut state = Analysers::default();
        let mut rng = TepRng::with_default_seed();
        let compositions = comps();
        let _ = sample(
            &mut state,
            &mut rng,
            0.0,
            &[0.0; N_CONTINUOUS],
            &compositions,
            Shutdown::default(),
        );
        let out = sample(
            &mut state,
            &mut rng,
            0.05,
            &[0.0; N_CONTINUOUS],
            &compositions,
            Shutdown::default(),
        );
        for index in 0..N_SAMPLED {
            assert_exact(
                out[N_CONTINUOUS + index],
                compositions[index],
                "a held reading",
            );
        }
        // Held from `reported`, not from `stored`: after a real sample the two
        // differ by a noise draw.
        for (a, b) in out[N_CONTINUOUS..].iter().zip(state.reported) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }

    /// The three analysers read three different streams, and two of them
    /// report only part of the composition.
    #[test]
    fn the_analysers_read_the_streams_they_are_named_for() {
        use crate::equilibrium::equilibrium;
        use crate::flows::{FlowDrift, flows};
        use crate::state::State;
        use crate::streams::{FeedConditions, streams};
        use crate::stripper::stripper;
        use crate::vessels::{TemperatureSeeds, unpack};

        let y = State::from_flat(&crate::constants::NOMINAL_STATE);
        let unpacked = unpack(&y, TemperatureSeeds::default()).expect("converges");
        let eq = equilibrium(&unpacked);
        let mut table = streams(&unpacked, &eq, &FeedConditions::default());
        let mut flow = flows(&y, &unpacked, &eq, &table, &[0.0; 20], FlowDrift::default());
        let _ = stripper(&mut table, &mut flow, unpacked.stripper.celsius);

        let xcmp = compositions(&table);
        // XMEAS(23..28): the reactor feed, A through F only.
        for (slot, c) in Component::ALL[..6].iter().enumerate() {
            assert_exact(
                xcmp[slot],
                table.composition[Stream::ReactorInlet][*c] * 100.0,
                "reactor feed",
            );
        }
        // XMEAS(37..41): the product, D through H only.
        for (slot, c) in Component::ALL[3..].iter().enumerate() {
            assert_exact(
                xcmp[14 + slot],
                table.composition[Stream::Product][*c] * 100.0,
                "product",
            );
        }
        // The purge is the only one reporting all eight.
        assert_exact(
            xcmp[6],
            table.composition[Stream::Purge][Component::A] * 100.0,
            "purge A",
        );
        // Percent, not fraction: the values should sum near 100 for the purge.
        let purge: f64 = xcmp[6..14].iter().sum();
        assert!(
            (purge - 100.0).abs() < 1.0,
            "the purge compositions sum to {purge}, not 100"
        );
    }
}
