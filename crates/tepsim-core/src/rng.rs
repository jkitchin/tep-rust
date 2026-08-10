//! The Tennessee Eastman random number generator.
//!
//! Every stochastic thing in the simulator, measurement noise and disturbance
//! walks alike, comes from this one generator. Reproducing it exactly is what
//! makes the whole validation ladder possible, and it is where the existing
//! Python port of TEP diverges from the Fortran on its very first sample.
//!
//! # Do not tidy this arithmetic
//!
//! The original is one line, `teprob.f:1551`:
//!
//! ```fortran
//! G = DMOD(G*9228907.D0, 4294967296.D0)
//! ```
//!
//! It looks like a textbook Lehmer generator, `state = state * a mod m` with
//! `a = 9228907` and `m = 2^32`, and the obvious modernisation is to do it in
//! `u64` or `u128` integer arithmetic. **That produces a different sequence.**
//!
//! `G` is a `DOUBLE PRECISION`, so the multiply happens in floating point.
//! `G` can be just under `2^32` and the multiplier is about `2^23.14`, so the
//! product reaches about `2^55.1`, comfortably past the `2^53` where `f64`
//! stops representing consecutive integers. The product is therefore rounded to
//! even before the modulus is taken, and that rounding is not an artefact to be
//! corrected: it is part of the algorithm that generated thirty years of
//! published TEP results.
//!
//! Both operations are exactly specified by IEEE 754. Multiplication rounds to
//! nearest even, and remainder is exact. So the sequence *is* bit-reproducible
//! in any language with `f64`, provided it is written the same way. See
//! `exact_integer_arithmetic_diverges_almost_immediately` below, which measures
//! how quickly the "corrected" version parts company.
//!
//! # The seed is larger than the modulus
//!
//! `teprob.f:1187` compiles in `G = 4651207995`, which is greater than `2^32`.
//! The first step reduces it. That is harmless, and it is preserved rather than
//! normalised, because normalising would change the stream.

/// The Tennessee Eastman linear congruential generator.
///
/// State is a single `f64`. See the module documentation for why it must stay
/// that way.
///
/// The original selects between two output scalings using the *sign of an
/// integer argument*, which is why `TESUB7(-1)` and `TESUB7(1)` do different
/// things. Here they are two methods: [`TepRng::unit`] and [`TepRng::signed`].
//
// Claims the whole subroutine, declarations and RETURN/END included, because
// coverage is about whether a routine has been ported. The per-method claims
// below additionally pin which line each piece came from.
//
// @port teprob.f:1547-1555
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TepRng {
    g: f64,
}

/// The multiplier, `teprob.f:1551`.
const MULTIPLIER: f64 = 9_228_907.0;

/// The modulus, `2^32`, `teprob.f:1551`.
const MODULUS: f64 = 4_294_967_296.0;

impl TepRng {
    /// The seed compiled into the original at `teprob.f:1187`.
    ///
    /// Larger than the modulus `2^32`, deliberately; see the module docs.
    pub const DEFAULT_SEED: f64 = 4_651_207_995.0;

    /// Start from an explicit seed.
    ///
    /// `reference/README.org` lists the seeds the published `d00` through `d21`
    /// datasets were generated with.
    #[must_use]
    pub const fn new(seed: f64) -> Self {
        Self { g: seed }
    }

    /// Start from [`TepRng::DEFAULT_SEED`].
    #[must_use]
    pub const fn with_default_seed() -> Self {
        Self::new(Self::DEFAULT_SEED)
    }

    /// The raw generator state, `COMMON/RANDSD/ G`.
    ///
    /// Exposed because the Tier 3 call-order diff compares it against the
    /// Fortran after every step.
    #[must_use]
    pub const fn state(&self) -> f64 {
        self.g
    }

    /// Overwrite the raw state.
    #[inline]
    pub const fn set_state(&mut self, g: f64) {
        self.g = g;
    }

    /// Advance the state and return it.
    ///
    /// Written exactly as the Fortran writes it. Rust's `%` on `f64` is IEEE
    /// remainder-after-truncation, the same operation as Fortran's `DMOD`, and
    /// both it and the multiply are exactly specified, so this reproduces the
    /// original bit for bit.
    //
    // @port teprob.f:1551
    #[inline]
    fn advance(&mut self) -> f64 {
        self.g = (self.g * MULTIPLIER) % MODULUS;
        self.g
    }

    /// A draw in `[0, 1)`. The original's `TESUB7(i)` for `i >= 0`.
    //
    // @port teprob.f:1552
    #[inline]
    pub fn unit(&mut self) -> f64 {
        self.advance() / MODULUS
    }

    /// A draw in `[-1, 1)`. The original's `TESUB7(i)` for `i < 0`.
    ///
    /// The expression is `2*g/m - 1`, in that association, matching the
    /// Fortran. Writing it as `2*(g/m) - 1` or `2*(g/m - 0.5)` would be
    /// algebraically identical and numerically different.
    //
    // @port teprob.f:1553
    #[inline]
    pub fn signed(&mut self) -> f64 {
        2.0 * self.advance() / MODULUS - 1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The whole point of the module docs, made measurable.
    ///
    /// Redoing the recurrence in exact integer arithmetic is the obvious
    /// "cleanup", and it is wrong. This records how fast the two part company,
    /// so nobody has to take the warning on faith.
    #[test]
    fn exact_integer_arithmetic_diverges_almost_immediately() {
        let mut float = TepRng::with_default_seed();
        // The same recurrence with no rounding at all. u128 keeps the full
        // 2^55 product that f64 cannot hold.
        let mut exact: u128 = TepRng::DEFAULT_SEED as u128;

        let mut first_divergence = None;
        for draw in 1..=1000_u32 {
            float.advance();
            exact = (exact * 9_228_907) % 4_294_967_296;
            if first_divergence.is_none() && float.state().to_bits() != (exact as f64).to_bits() {
                first_divergence = Some(draw);
            }
        }

        let at = first_divergence.expect(
            "integer arithmetic must diverge from the original; if it does not, \
             the premise of the module docs is wrong and they need rewriting",
        );
        assert!(
            at <= 3,
            "expected divergence within the first few draws, saw it at {at}"
        );
    }

    /// Products above 2^53 are the reason for all of the above, so confirm the
    /// recurrence really gets there, and at the rate the arithmetic predicts.
    ///
    /// A product overflows exactly when the state exceeds `2^53 / MULTIPLIER`.
    /// With a state spread roughly uniformly over `[0, 2^32)`, that is about
    /// 77% of draws: lossy rounding is the normal case here, not a rare edge.
    #[test]
    fn the_product_exceeds_2_pow_53_at_the_predicted_rate() {
        const TWO_POW_53: f64 = 9_007_199_254_740_992.0;
        let threshold = TWO_POW_53 / MULTIPLIER;
        let predicted = 1.0 - threshold / MODULUS;
        assert!(
            (0.75..0.80).contains(&predicted),
            "the analytic overflow fraction should be about 0.77, got {predicted}"
        );

        const DRAWS: usize = 100_000;
        let mut rng = TepRng::with_default_seed();
        let mut exceeded = 0_usize;
        for _ in 0..DRAWS {
            if rng.state() * MULTIPLIER > TWO_POW_53 {
                exceeded += 1;
            }
            rng.advance();
        }

        let observed = exceeded as f64 / DRAWS as f64;
        assert!(
            (observed - predicted).abs() < 0.01,
            "observed overflow fraction {observed} does not match the predicted \
             {predicted}; either the generator is not spread as assumed or the \
             threshold arithmetic is wrong"
        );
    }

    #[test]
    fn unit_draws_stay_in_range() {
        let mut rng = TepRng::with_default_seed();
        for _ in 0..100_000 {
            let x = rng.unit();
            assert!((0.0..1.0).contains(&x), "unit draw {x} out of range");
        }
    }

    #[test]
    fn signed_draws_stay_in_range() {
        let mut rng = TepRng::with_default_seed();
        for _ in 0..100_000 {
            let x = rng.signed();
            assert!((-1.0..1.0).contains(&x), "signed draw {x} out of range");
        }
    }

    /// A zero state would collapse the generator to a constant stream. The
    /// original has no guard against it, so confirm it does not happen in
    /// practice for the seeds we care about.
    #[test]
    fn the_state_never_reaches_zero() {
        for seed in [
            TepRng::DEFAULT_SEED,
            1_431_655_765.0,
            4_243_534_565.0,
            5_687_912_315.0,
        ] {
            let mut rng = TepRng::new(seed);
            for draw in 0..1_000_000 {
                if rng.advance() == 0.0 {
                    panic!("state hit zero at draw {draw} from seed {seed}");
                }
            }
        }
    }

    #[test]
    fn state_round_trips() {
        let mut rng = TepRng::with_default_seed();
        rng.advance();
        let saved = rng.state();
        rng.advance();
        assert_ne!(rng.state().to_bits(), saved.to_bits());
        rng.set_state(saved);
        assert_eq!(rng.state().to_bits(), saved.to_bits());
    }

    /// Both output modes read the same underlying state, so interleaving them
    /// must not change the state sequence.
    #[test]
    fn both_modes_consume_exactly_one_draw() {
        let mut a = TepRng::with_default_seed();
        let mut b = TepRng::with_default_seed();
        for i in 0..1000 {
            if i % 2 == 0 {
                a.unit();
            } else {
                a.signed();
            }
            b.advance();
        }
        assert_eq!(a.state().to_bits(), b.state().to_bits());
    }
}
