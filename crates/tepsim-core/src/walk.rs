//! The twelve disturbance channels, and how they advance.
//!
//! Ported from `teprob.f:340-406`. Each channel is a chain of cubic segments;
//! this is what decides when a segment has run out and builds the next one.
//!
//! # Nine walks and three spike trains
//!
//! The twelve channels are not the same kind of object, which the shared array
//! names hide.
//!
//! Channels 1 to 9 are **random walks**. When one runs out, `teprob.f:359-371`
//! evaluates the old segment at its endpoint, takes value and slope there, and
//! hands both to [`crate::disturbance::walk_segment`] to build a segment that
//! continues smoothly from them.
//!
//! Channels 10 to 12 are **spike trains**, and `teprob.f:372-396` gives them
//! their own rule. They alternate between two states:
//!
//! - **Dwelling.** The channel sits at zero for a randomly drawn interval. Its
//!   segment is `CDIST = IDVWLK / h^2` with everything else zero, which is a
//!   parabola rising from zero, and the dwell ends when it reaches 0.1.
//! - **Spiking.** Once the value exceeds 0.1, the channel is given a cubic
//!   through the current value and slope that lasts exactly 0.1 hours, with
//!   coefficients that drive it hard and then back down.
//!
//! So a spike channel is off, off, off, then briefly on. That is what
//! `IDV(17)`, `IDV(18)` and `IDV(20)` are: intermittent faults rather than
//! sustained drifts, which is why they are the hardest of the twenty to detect
//! and why the literature reports them as such.
//!
//! # The flag scales the dwell, not the schedule
//!
//! `CDIST(I) = IDVWLK(I) / h^2` at `teprob.f:391` is the only place the flag
//! enters a spike channel. With the disturbance off it is zero, so the parabola
//! is flat at zero, never reaches 0.1, and the channel dwells forever, drawing
//! a fresh interval each time. With it on, the parabola climbs and the channel
//! eventually spikes.
//!
//! The channel therefore keeps drawing at the same rate either way, exactly as
//! [`crate::disturbance::walk_segment`] does. See that module: the same
//! discipline, for the same reason.
//!
//! # Two draws that look alike and are not
//!
//! `teprob.f:388` draws an interval directly rather than through `TESUB5`:
//!
//! ```fortran
//!       ISD=-1
//!       HWLK=HSPAN(I)*TESUB7(ISD)+HZERO(I)
//! ```
//!
//! It is the *same expression* as `TESUB5`'s first line, and it is one draw
//! rather than three. So a spike channel entering a dwell costs one draw and a
//! walk channel re-segmenting costs three, which is where B-0027's
//! `9 * 3 + 3 = 30` comes from.
//!
//! The `ISD` there is a *borrowed local*, reused as the sign flag, and has
//! nothing to do with the shutdown. The oracle's instrumentation renames it to
//! `IRAND` for exactly that reason. This port uses a plain argument.
//!
//! # Variables
//!
//! | Fortran | Here | Meaning |
//! |---|---|---|
//! | `ADIST`..`DDIST` | [`Channel::segment`] | the current cubic |
//! | `TLAST`, `TNEXT` | [`Channel::since`], [`Segment::until`] | its span |
//! | `IDVWLK` | [`Walks::flags`] | which disturbances are active |

// Every float expression here reproduces `teprob.f`'s association and
// rounding exactly; see `crate::thermo`.
#![allow(
    clippy::suboptimal_flops,
    reason = "expressions must round exactly as teprob.f does; see module docs"
)]

use crate::disturbance::{CHANNEL_SPANS, Draws, Segment, SegmentStart, walk_segment};

/// How many disturbance channels there are.
pub const CHANNELS: usize = 12;

/// The first channel driven by the spike rule rather than a walk, zero-based.
///
/// `teprob.f:372` starts its loop at 10.
pub const FIRST_SPIKE_CHANNEL: usize = 9;

/// The value a spike channel's dwell must exceed before it fires
/// (`teprob.f:380`).
const SPIKE_THRESHOLD: f64 = 0.1;

/// How long a spike lasts, in hours (`teprob.f:385`).
const SPIKE_DURATION: f64 = 0.1;

/// One channel's current segment and when it started.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Channel {
    /// The cubic currently in force, and when it ends.
    pub segment: Segment,
    /// `TLAST`: when it started. The cubic is evaluated in time since this.
    pub since: f64,
}

impl Channel {
    /// The channel's value at time `t`.
    // @port teprob.f:1556-1588
    #[must_use]
    pub fn at(&self, t: f64) -> f64 {
        self.segment.at(t, self.since)
    }
}

/// The twelve channels, and which disturbances drive them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Walks {
    /// Each channel's current segment.
    pub channels: [Channel; CHANNELS],
    /// `IDVWLK`: whether each channel's disturbance is active.
    pub flags: [i32; CHANNELS],
}

impl Default for Walks {
    /// The state `TEINIT` leaves at `teprob.f:1360-1367`: every channel flat at
    /// its own centre, with its first segment ending at 0.1 hours.
    ///
    /// That is what makes a `d00` run a nominal run rather than a run with
    /// quiet noise on it.
    fn default() -> Self {
        let mut channels = [Channel {
            segment: Segment::default(),
            since: 0.0,
        }; CHANNELS];
        for (channel, spans) in channels.iter_mut().zip(CHANNEL_SPANS) {
            channel.segment.constant = spans.value_centre;
            channel.segment.until = 0.1;
        }
        Self {
            channels,
            flags: [0; CHANNELS],
        }
    }
}

/// Map the twenty disturbance flags onto the twelve channels.
///
/// `teprob.f:347-358`. Not a bijection: `IDV(8)` drives channels 1 and 2, and
/// `IDV(13)` drives 7 and 8, so those two faults move a pair of channels
/// together. Eight of the twenty disturbances drive no channel at all; they are
/// the step faults and the valve-sticking flags.
// @port teprob.f:347-358
#[must_use]
pub fn channel_flags(disturbances: &[f64; 20]) -> [i32; CHANNELS] {
    // teprob.f:340-346 has already reduced each to 0 or 1 by this point; see
    // `crate::plant::Inputs::clamped_disturbances`.
    let idv = |n: usize| i32::from(disturbances[n - 1] > 0.0);
    [
        idv(8),  // IDVWLK(1): A/C feed composition
        idv(8),  // IDVWLK(2): the same fault, second channel
        idv(9),  // IDVWLK(3): D feed temperature
        idv(10), // IDVWLK(4): C feed temperature
        idv(11), // IDVWLK(5): reactor coolant inlet temperature
        idv(12), // IDVWLK(6): condenser coolant inlet temperature
        idv(13), // IDVWLK(7): reaction kinetics
        idv(13), // IDVWLK(8): the same fault, second channel
        idv(16), // IDVWLK(9): stripper steam valve
        idv(17), // IDVWLK(10): reactor coolant valve, a spike train
        idv(18), // IDVWLK(11): condenser coolant valve, a spike train
        idv(20), // IDVWLK(12): reactor outlet, a spike train
    ]
}

/// Advance every channel that has run out, drawing as the original draws.
///
/// Call once per outer step, before the derivative; see [`mod@crate::plant`].
// @port teprob.f:347-406
pub fn advance<R: Draws + ?Sized>(
    walks: &mut Walks,
    rng: &mut R,
    t: f64,
    disturbances: &[f64; 20],
) {
    let flags = channel_flags(disturbances);
    walks.flags = flags;

    // teprob.f:359-371. The nine true walks.
    for (index, (channel, spans)) in walks
        .channels
        .iter_mut()
        .zip(&CHANNEL_SPANS)
        .enumerate()
        .take(FIRST_SPIKE_CHANNEL)
    {
        if t < channel.segment.until {
            continue;
        }
        let start = end_of_segment(channel);
        *channel = Channel {
            segment: walk_segment(rng, start, spans, flags[index]),
            since: start.since,
        };
    }

    // teprob.f:372-396. The three spike trains.
    for (index, (channel, spans)) in walks
        .channels
        .iter_mut()
        .zip(&CHANNEL_SPANS)
        .enumerate()
        .skip(FIRST_SPIKE_CHANNEL)
    {
        if t < channel.segment.until {
            continue;
        }
        let start = end_of_segment(channel);
        let segment = if start.value > SPIKE_THRESHOLD {
            // teprob.f:381-385. Fire: a cubic through the current value and
            // slope, lasting exactly 0.1 hours. No draw.
            Segment {
                constant: start.value,
                linear: start.slope,
                quadratic: -(3.0 * start.value + 0.2 * start.slope) / 0.01,
                cubic: (2.0 * start.value + 0.1 * start.slope) / 0.001,
                until: start.since + SPIKE_DURATION,
            }
        } else {
            // teprob.f:387-393. Dwell: flat at zero with a curvature set by
            // the flag, for a freshly drawn interval. One draw, and it is the
            // same expression as `TESUB5`'s first line.
            let h = spans.duration_span * rng.signed() + spans.duration_centre;
            Segment {
                constant: 0.0,
                linear: 0.0,
                quadratic: f64::from(flags[index]) / (h * h),
                cubic: 0.0,
                until: start.since + h,
            }
        };
        *channel = Channel {
            segment,
            since: start.since,
        };
    }

    // teprob.f:397-406. At time zero everything is reset, whatever else
    // happened above. The advance still ran and still drew, which is why a
    // t=0 evaluation is not simply "no walk activity".
    if t == 0.0 {
        let reset = Walks::default();
        walks.channels = reset.channels;
    }
}

/// Where the current segment ends: its value and slope at `TNEXT`.
///
/// `teprob.f:361-366`. Note it evaluates at `TNEXT`, not at `t`: the segment
/// hands over at its own scheduled end, so a step that arrives late does not
/// shift the chain.
fn end_of_segment(channel: &Channel) -> SegmentStart {
    let h = channel.segment.until - channel.since;
    let s = channel.segment;
    SegmentStart {
        // teprob.f:362-363, Horner.
        value: s.constant + h * (s.linear + h * (s.quadratic + h * s.cubic)),
        // teprob.f:364-365, its derivative.
        slope: s.linear + h * (2.0 * s.quadratic + 3.0 * h * s.cubic),
        // teprob.f:366
        since: channel.segment.until,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::TepRng;
    use crate::testing::assert_exact;

    fn draws(before: f64, after: f64) -> usize {
        let mut probe = TepRng::new(before);
        for step in 0..200 {
            if probe.state().to_bits() == after.to_bits() {
                return step;
            }
            let _ = probe.unit();
        }
        panic!("more than two hundred draws");
    }

    /// The mapping is not a bijection: two faults each drive two channels, and
    /// eight of the twenty drive none.
    #[test]
    fn the_flag_mapping_is_many_to_many() {
        let mut idv = [0.0; 20];
        idv[7] = 1.0; // IDV(8)
        let flags = channel_flags(&idv);
        assert_eq!(flags[0], 1, "IDV(8) drives channel 1");
        assert_eq!(flags[1], 1, "and channel 2");
        assert_eq!(flags[2..].iter().sum::<i32>(), 0, "and nothing else");

        let mut idv = [0.0; 20];
        idv[12] = 1.0; // IDV(13)
        let flags = channel_flags(&idv);
        assert_eq!(flags[6], 1, "IDV(13) drives channel 7");
        assert_eq!(flags[7], 1, "and channel 8");

        // The eight that drive nothing.
        for fault in [1, 2, 3, 4, 5, 6, 7, 14, 15, 19] {
            let mut idv = [0.0; 20];
            idv[fault - 1] = 1.0;
            assert_eq!(
                channel_flags(&idv),
                [0; CHANNELS],
                "IDV({fault}) should drive no walk channel"
            );
        }
    }

    /// A channel that has not run out is not touched, and draws nothing.
    #[test]
    fn a_channel_with_time_left_does_not_advance() {
        let mut walks = Walks::default();
        let mut rng = TepRng::with_default_seed();
        let before = rng.state();
        let snapshot = walks;

        advance(&mut walks, &mut rng, 0.05, &[0.0; 20]);
        assert_eq!(draws(before, rng.state()), 0, "nothing should have drawn");
        assert_eq!(walks.channels, snapshot.channels);
    }

    /// Nine walks at three draws each, plus one each from the three spike
    /// channels dwelling: thirty, which is B-0027's measurement.
    #[test]
    fn a_full_advance_costs_thirty_draws_at_the_nominal_point() {
        let mut walks = Walks::default();
        let mut rng = TepRng::with_default_seed();
        let before = rng.state();

        advance(&mut walks, &mut rng, 0.1, &[0.0; 20]);
        assert_eq!(draws(before, rng.state()), 9 * 3 + 3);
    }

    /// A spike channel with its disturbance off dwells forever: the parabola
    /// is flat at zero and never reaches the threshold.
    #[test]
    fn an_inactive_spike_channel_never_fires() {
        let mut walks = Walks::default();
        let mut rng = TepRng::with_default_seed();
        let mut t = 0.1;

        for _ in 0..200 {
            advance(&mut walks, &mut rng, t, &[0.0; 20]);
            for index in FIRST_SPIKE_CHANNEL..CHANNELS {
                let channel = &walks.channels[index];
                assert_exact(channel.segment.quadratic, 0.0, "flat while inactive");
                assert!(
                    channel.at(t).abs() < 1e-12,
                    "an inactive spike channel reached {}",
                    channel.at(t)
                );
            }
            t += 0.1;
        }
    }

    /// With its disturbance on, a spike channel climbs, fires, and comes back.
    ///
    /// Both branches of `teprob.f:380` in one run, which is the coverage this
    /// item owes.
    #[test]
    fn an_active_spike_channel_dwells_then_fires() {
        let mut walks = Walks::default();
        let mut rng = TepRng::with_default_seed();
        let mut idv = [0.0; 20];
        idv[16] = 1.0; // IDV(17), channel 10

        let mut t = 0.0;
        let mut fired = 0;
        let mut dwelled = 0;
        let mut peak: f64 = 0.0;
        for _ in 0..4_000 {
            advance(&mut walks, &mut rng, t, &idv);
            let channel = &walks.channels[FIRST_SPIKE_CHANNEL];
            // A fired segment lasts exactly 0.1 hours and has a nonzero cubic;
            // a dwell is flat with a nonzero quadratic.
            let span = channel.segment.until - channel.since;
            if channel.segment.cubic != 0.0 {
                fired += 1;
                assert!(
                    (span - SPIKE_DURATION).abs() < 1e-12,
                    "a spike lasted {span} hours, not 0.1"
                );
            } else if channel.segment.quadratic != 0.0 {
                dwelled += 1;
            }
            peak = peak.max(channel.at(t));
            t += 0.01;
        }
        assert!(fired > 0, "the channel never fired in 40 hours");
        assert!(dwelled > 0, "the channel never dwelled");
        assert!(
            peak > SPIKE_THRESHOLD,
            "the channel peaked at {peak}, below the threshold it must cross"
        );
    }

    /// The handover is at the segment's own scheduled end, not at the current
    /// time. Evaluating at `t` instead would let a late step shift the chain
    /// and drift the whole schedule.
    #[test]
    fn a_segment_hands_over_at_its_own_end_not_at_the_current_time() {
        let mut early = Walks::default();
        let mut late = Walks::default();
        let mut rng_a = TepRng::with_default_seed();
        let mut rng_b = TepRng::with_default_seed();

        // Both channels end at 0.1; advance one at 0.1 and the other at 0.19.
        advance(&mut early, &mut rng_a, 0.1, &[0.0; 20]);
        advance(&mut late, &mut rng_b, 0.19, &[0.0; 20]);
        assert_eq!(
            early.channels, late.channels,
            "the new segment depends on when the step arrived, so a late step \
             shifts the chain"
        );
    }

    /// Time zero resets everything, and the reset is exactly `TEINIT`'s state.
    #[test]
    fn time_zero_resets_every_channel() {
        let mut walks = Walks::default();
        let mut rng = TepRng::with_default_seed();
        // Move the walks well away from their starting point.
        let mut t = 0.1;
        for _ in 0..20 {
            advance(&mut walks, &mut rng, t, &[0.0; 20]);
            t += 0.5;
        }
        assert_ne!(walks.channels, Walks::default().channels);

        advance(&mut walks, &mut rng, 0.0, &[0.0; 20]);
        assert_eq!(
            walks.channels,
            Walks::default().channels,
            "the t=0 reset does not reproduce TEINIT's state"
        );
    }

    /// The reset happens *after* the advance, so a t=0 call still draws.
    ///
    /// The obvious implementation returns early at t=0, which produces the
    /// right state and the wrong stream. Same shape as B-0028's bug.
    #[test]
    fn the_time_zero_reset_does_not_skip_the_draws() {
        let mut walks = Walks::default();
        let mut rng = TepRng::with_default_seed();
        let before = rng.state();
        // Every channel is due at 0.1, and t=0 is not past that, so nothing
        // advances and nothing draws.
        advance(&mut walks, &mut rng, 0.0, &[0.0; 20]);
        assert_eq!(draws(before, rng.state()), 0);

        // But with a channel already due, the advance runs and draws before
        // the reset discards its result.
        let mut walks = Walks::default();
        for channel in &mut walks.channels {
            channel.segment.until = 0.0;
        }
        let before = rng.state();
        advance(&mut walks, &mut rng, 0.0, &[0.0; 20]);
        assert!(
            draws(before, rng.state()) > 0,
            "the reset short-circuited the advance, so the state is right and \
             the generator is behind"
        );
        assert_eq!(walks.channels, Walks::default().channels);
    }
}
