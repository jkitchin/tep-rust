//! Accumulating a comparison, and reporting it as numbers rather than a
//! verdict.
//!
//! # Why a histogram and not a pass mark
//!
//! A threshold answers "is it broken today". It cannot answer "did it get
//! worse", which is the question that actually catches a regression, because
//! two runs that both pass at 1e-13 can differ by two orders of magnitude
//! between them. So the unit of output here is a distribution: how many cases
//! agreed bit for bit, how many were one ULP apart, where the worst one was and
//! what its inputs were. Those go in `LOG.org`, and the next session compares
//! against them. See `CLAUDE.md`, "Record numbers, not verdicts".

use core::fmt;

/// The ULP distance between two `f64`s, as a count of representable values
/// between them.
///
/// Uses the same monotonic key as [`f64::total_cmp`]: flipping the sign bit for
/// positives and every bit for negatives turns the IEEE-754 layout into a
/// two's-complement integer whose order matches the real line. Adjacent floats
/// then differ by exactly one, including across an exponent boundary and across
/// the subnormal-to-normal boundary.
///
/// Positive and negative zero are reported as zero apart, since they are the
/// same number even though they are different keys.
///
/// A distance that spans zero counts the two zeros as two steps rather than
/// one, so it is one larger than the true count. Such a distance exceeds 2^53
/// by construction, which already means "catastrophically different", so the
/// off-by-one buys nothing worth the special case.
#[must_use]
pub fn ulp_distance(a: f64, b: f64) -> u64 {
    if a.to_bits() == (-0.0_f64).to_bits() {
        return ulp_distance(0.0, b);
    }
    if b.to_bits() == (-0.0_f64).to_bits() {
        return ulp_distance(a, 0.0);
    }
    let distance = i128::from(total_order_key(a)) - i128::from(total_order_key(b));
    u64::try_from(distance.unsigned_abs()).unwrap_or(u64::MAX)
}

fn total_order_key(x: f64) -> i64 {
    let bits = x.to_bits() as i64;
    bits ^ ((((bits >> 63) as u64) >> 1) as i64)
}

/// The relative error of `actual` against `reference`.
///
/// Zero when the two are numerically equal, including when both are zero.
/// Infinite when `reference` is zero and `actual` is not, which is a real
/// failure rather than a division to be papered over: `TESUB1` returns exactly
/// zero at 0 degrees for `ITY` 0 and 1, so a port that returns anything else
/// there is wrong by an unboundedly large factor and should say so.
#[must_use]
pub fn relative_error(actual: f64, reference: f64) -> f64 {
    let difference = (actual - reference).abs();
    if difference.is_nan() {
        return f64::NAN;
    }
    if difference <= 0.0 {
        // An absolute value, so this is exactly zero: the two agree.
        return 0.0;
    }
    let scale = reference.abs();
    if scale > 0.0 {
        difference / scale
    } else {
        f64::INFINITY
    }
}

/// ULP distances counted individually. Anything past this is already a
/// catastrophe by Tier 1 standards, so it goes in one overflow bucket while
/// [`Comparison::max_ulp`] keeps the exact worst value.
pub const EXACT_BUCKETS: usize = 16;

/// The worst case seen for one metric.
#[derive(Clone, Copy, Debug)]
struct Worst<T, C> {
    value: T,
    case: Option<C>,
}

impl<T: Default, C> Default for Worst<T, C> {
    fn default() -> Self {
        Self {
            value: T::default(),
            case: None,
        }
    }
}

/// A running comparison of a Rust routine against the Fortran, over some pool
/// of cases.
///
/// Generic over the case label so the same accumulator serves Tier 2 later; it
/// only needs the label to be cheap to copy and, at report time, printable.
#[derive(Clone, Debug)]
pub struct Comparison<C> {
    what: String,
    cases: u64,
    histogram: [u64; EXACT_BUCKETS],
    beyond: u64,
    worst_ulp: Worst<u64, C>,
    worst_relative: Worst<f64, C>,
    non_finite: u64,
    mismatched_non_finite: u64,
}

impl<C: Copy> Comparison<C> {
    /// Start a comparison labelled with what is being compared, for example
    /// `"TESUB1 ity=0"`.
    #[must_use]
    pub fn new(what: impl Into<String>) -> Self {
        Self {
            what: what.into(),
            cases: 0,
            histogram: [0; EXACT_BUCKETS],
            beyond: 0,
            worst_ulp: Worst::default(),
            worst_relative: Worst::default(),
            non_finite: 0,
            mismatched_non_finite: 0,
        }
    }

    /// Record one evaluation: what the port produced, and what the Fortran did.
    pub fn observe(&mut self, case: C, actual: f64, reference: f64) {
        self.cases += 1;

        if !actual.is_finite() || !reference.is_finite() {
            self.non_finite += 1;
            // Two identical infinities agree; anything else, NaN included, does
            // not. Comparing bits rather than values is deliberate: NaN is not
            // equal to itself, and here that would silently pass.
            if actual.to_bits() != reference.to_bits() {
                self.mismatched_non_finite += 1;
                self.worst_relative = Worst {
                    value: f64::INFINITY,
                    case: Some(case),
                };
                self.worst_ulp = Worst {
                    value: u64::MAX,
                    case: Some(case),
                };
            }
            return;
        }

        let ulp = ulp_distance(actual, reference);
        if (ulp as usize) < EXACT_BUCKETS {
            self.histogram[ulp as usize] += 1;
        } else {
            self.beyond += 1;
        }
        if self.worst_ulp.case.is_none() || ulp > self.worst_ulp.value {
            self.worst_ulp = Worst {
                value: ulp,
                case: Some(case),
            };
        }

        let relative = relative_error(actual, reference);
        if self.worst_relative.case.is_none() || relative > self.worst_relative.value {
            self.worst_relative = Worst {
                value: relative,
                case: Some(case),
            };
        }
    }

    /// How many evaluations were recorded.
    #[must_use]
    pub fn cases(&self) -> u64 {
        self.cases
    }

    /// The largest relative error seen. Zero if nothing was recorded.
    #[must_use]
    pub fn max_relative_error(&self) -> f64 {
        self.worst_relative.value
    }

    /// The case that produced it.
    #[must_use]
    pub fn worst_relative_case(&self) -> Option<C> {
        self.worst_relative.case
    }

    /// The largest ULP distance seen, exactly, however far into the overflow
    /// bucket it fell. [`u64::MAX`] if it came from a non-finite mismatch.
    #[must_use]
    pub fn max_ulp(&self) -> u64 {
        self.worst_ulp.value
    }

    /// The case that produced it.
    #[must_use]
    pub fn worst_ulp_case(&self) -> Option<C> {
        self.worst_ulp.case
    }

    /// How many evaluations had a non-finite value on either side, and how many
    /// of those disagreed.
    #[must_use]
    pub fn non_finite(&self) -> (u64, u64) {
        (self.non_finite, self.mismatched_non_finite)
    }

    /// The smallest ULP distance at or below which at least `fraction` of the
    /// finite cases fall.
    ///
    /// `None` when that point lies in the overflow bucket, which means only
    /// "at least [`EXACT_BUCKETS`]"; the exact worst value is
    /// [`Comparison::max_ulp`] either way.
    #[must_use]
    pub fn ulp_percentile(&self, fraction: f64) -> Option<u64> {
        let finite: u64 = self.histogram.iter().sum::<u64>() + self.beyond;
        if finite == 0 {
            return None;
        }
        // Ceiling, so p100 means every case and p0 means at least one.
        let wanted = (fraction * finite as f64).ceil().max(1.0) as u64;
        let mut seen = 0;
        for (distance, count) in self.histogram.iter().enumerate() {
            seen += count;
            if seen >= wanted {
                return Some(distance as u64);
            }
        }
        None
    }

    /// Fail if the comparison does not meet `tolerance`, printing the whole
    /// report rather than one number.
    ///
    /// Also fails when nothing was recorded. A sweep that silently generated no
    /// cases would otherwise pass every tolerance ever set, which is the one
    /// way a validation tier can be worse than useless.
    ///
    /// # Panics
    ///
    /// If no cases were observed, if any non-finite values disagreed, or if the
    /// maximum relative error exceeds `tolerance`.
    pub fn assert_within(&self, tolerance: f64)
    where
        C: fmt::Display,
    {
        assert!(
            self.cases > 0,
            "{}: no cases were compared, so this proves nothing",
            self.what
        );
        assert!(
            self.mismatched_non_finite == 0,
            "{} disagreed on {} non-finite value(s)\n{self}",
            self.what,
            self.mismatched_non_finite
        );
        assert!(
            self.max_relative_error() <= tolerance,
            "{} exceeded its tolerance of {tolerance:e}\n{self}",
            self.what
        );
    }
}

/// The report format. Stable enough to paste into a `LOG.org` entry, and
/// stable enough that two entries can be diffed against each other.
impl<C: Copy + fmt::Display> fmt::Display for Comparison<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "tier1 {}", self.what)?;
        writeln!(f, "  cases          : {}", self.cases)?;

        match self.worst_relative.case {
            Some(case) => writeln!(
                f,
                "  max rel err    : {:.3e} at {case}",
                self.worst_relative.value
            )?,
            None => writeln!(f, "  max rel err    : (no cases)")?,
        }
        match self.worst_ulp.case {
            Some(case) => writeln!(f, "  max ulp        : {} at {case}", self.max_ulp())?,
            None => writeln!(f, "  max ulp        : (no cases)")?,
        }

        write!(f, "  ulp percentiles:")?;
        for (name, fraction) in [("p50", 0.50), ("p90", 0.90), ("p99", 0.99), ("p100", 1.0)] {
            match self.ulp_percentile(fraction) {
                Some(distance) => write!(f, " {name}={distance}")?,
                None => write!(f, " {name}>={}", EXACT_BUCKETS)?,
            }
        }
        writeln!(f)?;

        write!(f, "  ulp histogram  :")?;
        let mut printed = false;
        for (distance, count) in self.histogram.iter().enumerate() {
            if *count > 0 {
                write!(f, " {distance}:{count}")?;
                printed = true;
            }
        }
        if self.beyond > 0 {
            write!(f, " >={}:{}", EXACT_BUCKETS, self.beyond)?;
            printed = true;
        }
        if !printed {
            write!(f, " (empty)")?;
        }
        writeln!(f)?;

        write!(
            f,
            "  non-finite     : {} seen, {} mismatched",
            self.non_finite, self.mismatched_non_finite
        )
    }
}
