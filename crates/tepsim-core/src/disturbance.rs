//! The two generator-consuming utilities: walk segments and measurement noise.
//!
//! Ported from `teprob.f:1506-1546`. `TESUB1` through `TESUB4` are pure
//! functions of their arguments and live in [`mod@crate::thermo`]; these two
//! are not, because they draw.
//!
//! # `TESUB5`: one segment of a random walk
//!
//! A disturbance channel is a chain of cubic segments. Each one is built from
//! where the previous ended and three fresh draws, and it is a *Hermite*
//! interpolant: it matches value and slope at both ends, which is what keeps
//! the walk continuously differentiable across knots rather than kinked.
//!
//! \\[
//!   h = H_\\text{span} u_1 + H_0, \\qquad
//!   s_1 = S_\\text{span} u_2 f + S_0, \\qquad
//!   s_1' = S'_\\text{span} u_3 f
//! \\]
//!
//! with \\(u_i\\) signed draws on [-1, 1) and \\(f\\) the channel's disturbance
//! flag. Then the cubic through \\((0, s, s')\\) and \\((h, s_1, s_1')\\) is
//!
//! \\[
//!   c = \\frac{3(s_1 - s) - h(s_1' + 2s')}{h^2}, \\qquad
//!   d = \\frac{2(s - s_1) + h(s_1' + s')}{h^3}
//! \\]
//!
//! # The flag multiplies the endpoints but not the duration
//!
//! `IDVFLAG` is zero unless the channel's disturbance is switched on, and it
//! multiplies the *endpoint* draws at `teprob.f:1529-1530` but not the
//! *duration* draw at `1528`. So an inactive channel still consumes all three
//! draws and still re-segments at random intervals; it simply lands on
//! \\(S_0\\) with zero slope every time, which is a constant.
//!
//! That matters twice over. The walk is a flat line when the disturbance is
//! off, as it should be. And the generator advances exactly the same amount
//! either way, so switching a disturbance on does not shift the stream for
//! every *other* channel. A port that skipped the draws when the flag was zero
//! would produce identical disturbance values and a completely different noise
//! sequence.
//!
//! # `TESUB6`: twelve uniforms, summed
//!
//! \\[ x = \\left(\\sum_{i=1}^{12} u_i - 6\\right)\\sigma \\]
//!
//! The Irwin-Hall construction: twelve uniforms on [0, 1) have mean 6 and
//! variance 1, so subtracting 6 gives an approximately standard normal with a
//! variance that is exactly 1 and tails that stop at ±6. It is not Gaussian,
//! and Tier 5 compares distributions rather than assuming one, so the
//! difference is measured rather than assumed away.
//!
//! Twelve draws per sample, and 22 measurements plus 19 compositions are
//! sampled per step, which is where the bulk of the generator traffic goes:
//! B-0027 measured 264 draws per evaluation for the continuous measurements
//! alone.
//!
//! ## The summation order is *not* load-bearing here, unusually
//!
//! `teprob.f:1542-1543` writes `X=X+TESUB7(I)` in a loop from zero, and the
//! port writes the same loop. The obvious assumption, from every balance in
//! [`mod@crate::balances`], is that the order is load-bearing and reassociating
//! it would change the last bits.
//!
//! Measured, it does not. Over 500,000 consecutive samples from the real
//! stream, summing the twelve draws forwards and backwards gives *bit-identical*
//! results every time. The reason is conditioning: all twelve addends lie in
//! [0, 1) and the running total never exceeds 12, so every addition is between
//! operands within a factor of about twelve, and nothing is lost to shift out.
//!
//! It is still written in call order, because that is what the listing says and
//! this module is meant to be checkable against it a line at a time. But the
//! *draw* order remains absolutely load-bearing: which draw feeds which
//! quantity is what Tier 3 exists to check. Only the arithmetic is forgiving.
//!
//! No transcendental is involved in either routine, so the `libm-system` run
//! offers no extra protection here: both must be bit-exact in *both*
//! configurations.

extern crate alloc;

use crate::rng::TepRng;

/// Where a walk segment starts.
///
/// `S`, `SP` and `TLAST` in the original's argument list.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SegmentStart {
    /// The value the previous segment ended at.
    pub value: f64,
    /// The slope it ended with.
    pub slope: f64,
    /// The time it ended.
    pub since: f64,
}

/// One channel's five span parameters and its disturbance flag.
///
/// Constant per channel; `TEINIT` sets them at `teprob.f:1297-1359` and
/// nothing writes them again. B-0030 transcribes the table.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChannelSpans {
    /// `HSPAN`, the half-range of the segment duration.
    pub duration_span: f64,
    /// `HZERO`, the centre of the segment duration.
    pub duration_centre: f64,
    /// `SSPAN`, the half-range of the endpoint value.
    pub value_span: f64,
    /// `SZERO`, the centre of the endpoint value, and the value an inactive
    /// channel holds exactly.
    pub value_centre: f64,
    /// `SPSPAN`, the half-range of the endpoint slope.
    pub slope_span: f64,
}

/// A cubic walk segment, and when it ends.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Segment {
    /// `ADIST`, the constant term.
    pub constant: f64,
    /// `BDIST`, the linear term.
    pub linear: f64,
    /// `CDIST`, the quadratic term.
    pub quadratic: f64,
    /// `DDIST`, the cubic term.
    pub cubic: f64,
    /// `TNEXT`, when this segment ends.
    pub until: f64,
}

impl Segment {
    /// Evaluate the segment at time `t`, as `TESUB8` does.
    ///
    /// Horner in the time since the segment started, matching
    /// `teprob.f:1585-1586`.
    #[must_use]
    pub fn at(&self, t: f64, since: f64) -> f64 {
        let h = t - since;
        self.constant + h * (self.linear + h * (self.quadratic + h * self.cubic))
    }
}

/// Build the next walk segment, consuming three draws.
///
/// `flag` is `IDVWLK(I)`: zero unless the channel's disturbance is active. It
/// scales the endpoint draws only; see the module documentation.
// @port teprob.f:1506-1537
#[must_use]
pub fn walk_segment(
    rng: &mut TepRng,
    start: SegmentStart,
    spans: &ChannelSpans,
    flag: i32,
) -> Segment {
    // teprob.f:1528-1530. Three draws, in this order, always.
    let h = spans.duration_span * rng.signed() + spans.duration_centre;
    let end_value = spans.value_span * rng.signed() * f64::from(flag) + spans.value_centre;
    let end_slope = spans.slope_span * rng.signed() * f64::from(flag);

    // teprob.f:1531-1535. `H**2` and `H**3` are integer powers, so they expand
    // to multiplication rather than `pow`; B-0023 measured that gfortran uses
    // binary exponentiation, which for these two exponents is unambiguous.
    Segment {
        constant: start.value,
        linear: start.slope,
        quadratic: (3.0 * (end_value - start.value) - h * (end_slope + 2.0 * start.slope))
            / (h * h),
        cubic: (2.0 * (start.value - end_value) + h * (end_slope + start.slope)) / (h * h * h),
        until: start.since + h,
    }
}

/// One measurement-noise sample of standard deviation `std`, consuming twelve
/// draws.
///
/// See the module documentation: the sum accumulates in call order and must
/// not be reassociated.
// @port teprob.f:1538-1546
#[must_use]
pub fn noise(rng: &mut TepRng, std: f64) -> f64 {
    // teprob.f:1541-1543
    let mut x = 0.0;
    for _ in 0..12 {
        x += rng.unit();
    }
    // teprob.f:1544
    (x - 6.0) * std
}

/// One draw from the generator, as Tier 3 records it.
///
/// The value alone is not enough. `TESUB7` returns two different scalings
/// depending on the sign of its argument (`teprob.f:1552-1553`), so a port
/// could call the wrong one the right number of times and produce a stream
/// that is wrong in a way no count would show.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Draw {
    /// What the draw returned.
    pub value: f64,
    /// Whether it was the signed form on [-1, 1) rather than [0, 1).
    ///
    /// `TESUB7(I)` with `I < 0` gives the signed form. `TESUB5` uses it for
    /// all three of its draws and `TESUB6` uses the unsigned form for all
    /// twelve, so a whole routine calling the wrong one is possible and is
    /// exactly what this field catches.
    pub signed: bool,
}

/// A generator that records what it hands out.
///
/// Not used by the model, and deliberately not a variant of [`TepRng`]: the
/// shipped plant holds a plain generator with no buffer behind it, and you
/// only get a trace by constructing one of these on purpose.
#[derive(Clone, Debug)]
pub struct TracingRng {
    rng: TepRng,
    draws: alloc::vec::Vec<Draw>,
}

impl TracingRng {
    /// Wrap a generator and start recording.
    #[must_use]
    pub fn new(rng: TepRng) -> Self {
        Self {
            rng,
            draws: alloc::vec::Vec::new(),
        }
    }

    /// One draw on [0, 1), recorded.
    pub fn unit(&mut self) -> f64 {
        let value = self.rng.unit();
        self.draws.push(Draw {
            value,
            signed: false,
        });
        value
    }

    /// One draw on [-1, 1), recorded.
    pub fn signed(&mut self) -> f64 {
        let value = self.rng.signed();
        self.draws.push(Draw {
            value,
            signed: true,
        });
        value
    }

    /// Everything drawn so far, in order.
    #[must_use]
    pub fn draws(&self) -> &[Draw] {
        &self.draws
    }

    /// Forget the recorded draws, keeping the generator where it is.
    ///
    /// Tier 3 compares one step at a time, so the buffer is cleared between
    /// steps rather than growing across a run: a 48-hour run makes tens of
    /// millions of draws, and comparing them one step at a time reports the
    /// first divergence at the step it happened rather than as an index into
    /// something enormous.
    pub fn clear(&mut self) {
        self.draws.clear();
    }

    /// The underlying generator.
    #[must_use]
    pub const fn state(&self) -> f64 {
        self.rng.state()
    }
}

/// How many draws [`noise`] consumes. Asserted rather than assumed, because
/// Tier 3 depends on it.
pub const NOISE_DRAWS: usize = 12;

/// How many draws [`walk_segment`] consumes, whatever the flag.
pub const SEGMENT_DRAWS: usize = 3;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::assert_exact;

    fn spans() -> ChannelSpans {
        // Channel 1's, from teprob.f:1297-1301.
        ChannelSpans {
            duration_span: 0.2,
            duration_centre: 0.5,
            value_span: 0.03,
            value_centre: 0.485,
            slope_span: 0.0,
        }
    }

    /// How many draws a call makes, recovered from the generator word.
    fn draws(before: f64, after: f64) -> usize {
        let mut probe = TepRng::new(before);
        for step in 0..100 {
            // Exact by construction: the generator word is a value the
            // sequence takes, not an arithmetic result to compare loosely.
            if probe.state().to_bits() == after.to_bits() {
                return step;
            }
            let _ = probe.unit();
        }
        panic!("more than a hundred draws");
    }

    /// Twelve draws per noise sample, whatever the standard deviation.
    #[test]
    fn a_noise_sample_consumes_exactly_twelve_draws() {
        for std in [0.0, 0.0012, 22.0] {
            let mut rng = TepRng::with_default_seed();
            let before = rng.state();
            let _ = noise(&mut rng, std);
            assert_eq!(draws(before, rng.state()), NOISE_DRAWS);
        }
    }

    /// Three draws per segment, *including* when the channel is inactive.
    /// Skipping them would keep the values right and desynchronise everything
    /// else that draws.
    #[test]
    fn a_segment_consumes_three_draws_whether_or_not_the_channel_is_active() {
        for flag in [0, 1] {
            let mut rng = TepRng::with_default_seed();
            let before = rng.state();
            let _ = walk_segment(
                &mut rng,
                SegmentStart {
                    value: 0.485,
                    slope: 0.0,
                    since: 0.0,
                },
                &spans(),
                flag,
            );
            assert_eq!(draws(before, rng.state()), SEGMENT_DRAWS, "flag {flag}");
        }
    }

    /// An inactive channel lands exactly on its centre with zero slope, so the
    /// walk is a constant. The duration still varies.
    #[test]
    fn an_inactive_channel_is_a_constant_at_its_centre() {
        let mut rng = TepRng::with_default_seed();
        let s = spans();
        let start = SegmentStart {
            value: s.value_centre,
            slope: 0.0,
            since: 0.0,
        };
        let segment = walk_segment(&mut rng, start, &s, 0);

        assert_exact(segment.constant, s.value_centre, "starts at the centre");
        assert_exact(segment.linear, 0.0, "with zero slope");
        // And it stays there for the whole segment.
        for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let t = fraction * (segment.until - start.since);
            let value = segment.at(t, start.since);
            assert!(
                (value - s.value_centre).abs() < 1e-12,
                "an inactive channel moved to {value} at t={t}"
            );
        }
        // The duration is still random, so two segments differ in length.
        let second = walk_segment(&mut rng, start, &s, 0);
        assert!(
            second.until.to_bits() != segment.until.to_bits(),
            "the duration draw was skipped, so the stream is wrong"
        );
    }

    /// The segment is a Hermite interpolant: it matches value and slope at
    /// both ends. That is what makes the walk differentiable across knots, and
    /// it is an independent statement of `teprob.f:1533-1534`.
    #[test]
    fn a_segment_matches_value_and_slope_at_both_ends() {
        let mut rng = TepRng::new(12345.0);
        let s = ChannelSpans {
            slope_span: 0.1,
            ..spans()
        };
        let start = SegmentStart {
            value: 0.4,
            slope: 0.02,
            since: 1.5,
        };
        let segment = walk_segment(&mut rng, start, &s, 1);

        assert_exact(
            segment.at(start.since, start.since),
            start.value,
            "value at t0",
        );
        // The slope at the start is the linear coefficient by construction.
        assert_exact(segment.linear, start.slope, "slope at t0");

        // At the far end, recompute what the draws must have produced and
        // check the cubic lands there.
        let h = segment.until - start.since;
        let end = segment.at(segment.until, start.since);
        let numeric_slope = (segment.at(segment.until + 1e-6, start.since) - end) / 1e-6;
        // Rebuild the endpoint from the same three draws.
        let mut replay = TepRng::new(12345.0);
        let h2 = s.duration_span * replay.signed() + s.duration_centre;
        let end_value = s.value_span * replay.signed() + s.value_centre;
        let end_slope = s.slope_span * replay.signed();
        assert!((h - h2).abs() < 1e-12, "the duration draw came first");
        assert!(
            (end - end_value).abs() < 1e-9,
            "the cubic ends at {end}, not at the drawn {end_value}"
        );
        assert!(
            (numeric_slope - end_slope).abs() < 1e-4,
            "the cubic ends with slope {numeric_slope}, not the drawn {end_slope}"
        );
    }

    /// Twelve uniforms minus six has mean zero and variance one, so the sample
    /// is a standard deviation in the units its argument names. Checked
    /// statistically, since that is the only way to check a distribution.
    #[test]
    fn the_noise_has_the_standard_deviation_it_is_given() {
        let mut rng = TepRng::with_default_seed();
        let target = 0.25;
        let n = 20_000;
        let (mut sum, mut sum_sq) = (0.0, 0.0);
        let mut extreme: f64 = 0.0;
        for _ in 0..n {
            let x = noise(&mut rng, target);
            sum += x;
            sum_sq += x * x;
            extreme = extreme.max(x.abs());
        }
        let mean = sum / f64::from(n);
        let variance = sum_sq / f64::from(n) - mean * mean;
        let sd = libm::sqrt(variance);
        assert!(mean.abs() < 0.01, "mean {mean}");
        assert!(
            (sd - target).abs() < 0.01,
            "standard deviation {sd}, wanted {target}"
        );
        // Irwin-Hall is bounded, unlike a Gaussian: nothing can exceed 6 sigma.
        assert!(
            extreme <= 6.0 * target,
            "a sample reached {extreme}, past the 6-sigma bound the \
             construction guarantees"
        );
    }

    /// The sum is the twelve draws in call order, and the *draw* order is what
    /// matters: sample 1 must be the first draw, not the twelfth.
    #[test]
    fn the_sample_is_the_twelve_draws_in_call_order() {
        let mut rng = TepRng::new(999.0);
        let ours = noise(&mut rng, 1.0);

        let mut replay = TepRng::new(999.0);
        let draws: [f64; 12] = core::array::from_fn(|_| replay.unit());
        let mut forward = 0.0;
        for value in draws {
            forward += value;
        }
        assert_exact(ours, (forward - 6.0) * 1.0, "forward order");
        assert_exact(rng.state(), replay.state(), "the same twelve draws");
    }

    /// The *summation* order, unlike the draw order, turns out not to matter.
    ///
    /// This is the opposite of every balance in `crate::balances`, and it was
    /// asserted the other way round first. Over 500,000 consecutive samples
    /// from the real stream, forward and reversed summation are bit-identical
    /// every time: the twelve addends are all in [0, 1) and the running total
    /// stays under 12, so nothing is lost to shift out.
    ///
    /// Recorded as a test rather than as a remark so that a future change to
    /// `unit()`'s range, which would break the conditioning argument, fails
    /// here instead of somewhere subtler.
    #[test]
    fn the_summation_order_does_not_change_the_result() {
        let mut rng = TepRng::with_default_seed();
        let mut differing = 0;
        let samples = 20_000;
        for _ in 0..samples {
            let draws: [f64; 12] = core::array::from_fn(|_| rng.unit());
            let mut forward = 0.0;
            for value in draws {
                forward += value;
            }
            let mut backward = 0.0;
            for value in draws.iter().rev() {
                backward += value;
            }
            if forward.to_bits() != backward.to_bits() {
                differing += 1;
            }
            // The conditioning argument, checked directly.
            assert!(
                draws.iter().all(|d| (0.0..1.0).contains(d)),
                "a draw left [0, 1), so the reassociation argument no longer \
                 holds and the summation order becomes load-bearing"
            );
        }
        assert_eq!(
            differing, 0,
            "{differing} of {samples} samples changed under reassociation, so \
             the summation order *is* load-bearing after all and the module \
             documentation is wrong"
        );
    }
}
