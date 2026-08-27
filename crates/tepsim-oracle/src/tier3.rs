//! Tier 3: the generator draw trace, and the differ that compares two of them.
//!
//! # What Tier 3 is for, and why nothing else covers it
//!
//! Every stochastic quantity in the model comes from one generator word. A
//! port can therefore consume the right *number* of draws in the wrong
//! *order*, or the right values through the wrong scaling, and be wrong in a
//! way that no other tier notices:
//!
//! - Tier 2 pins the generator before each evaluation, so it never sees the
//!   stream advance at all.
//! - Tier 5 compares distributions, and a permuted stream has the same
//!   distribution as the original by construction.
//! - Tier 4 would eventually diverge, but so does a correct port, for
//!   `libm` reasons, so a divergence there proves nothing.
//!
//! B-0028 demonstrated the failure concretely at Tier 1: skipping the two
//! endpoint draws of an inactive walk channel, which are multiplied by zero
//! anyway, produces *identical* segment values and leaves the generator two
//! steps behind. Every later draw in the run then differs.
//!
//! # The trace is per-step, not per-run
//!
//! `COMMON/RNGTRC/` holds `TRACE_CAPACITY` draws, and the harness
//! clears it before each evaluation. A 48-hour run makes tens of millions of
//! draws; comparing them one step at a time reports the first divergence *at
//! the step it happened*, which is what a person debugging needs, rather than
//! as an index into something enormous.
//!
//! B-0027 measured the worst case at 522 draws in one evaluation, so the
//! capacity is eight times what is needed. The Fortran counter keeps counting
//! past it, so an overflow is reported rather than silently truncating.

use core::fmt;

use tepsim_core::disturbance::Draw;

/// Where two traces first disagree, and how.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Divergence {
    /// The index of the first differing draw, zero-based.
    pub index: usize,
    /// What the port produced, if it produced anything at that index.
    pub ours: Option<Draw>,
    /// What the Fortran produced, likewise.
    pub theirs: Option<Draw>,
}

impl fmt::Display for Divergence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "draw {} differs: ", self.index)?;
        match (self.ours, self.theirs) {
            (Some(a), Some(b)) if a.signed != b.signed => write!(
                f,
                "the port used the {} form and the Fortran the {} form. \
                 TESUB7 returns [0,1) for a non-negative argument and [-1,1) \
                 for a negative one, so this is a call with the wrong sign \
                 flag, not a wrong value.",
                form(a.signed),
                form(b.signed)
            ),
            (Some(a), Some(b)) => write!(
                f,
                "port {:?} ({}), Fortran {:?}. Same scaling, different value, \
                 so the generator is already in a different place: look for \
                 the first *count* mismatch in an earlier step.",
                a.value,
                form(a.signed),
                b.value
            ),
            (Some(a), None) => write!(
                f,
                "the port drew {:?} ({}) and the Fortran had stopped: the port \
                 draws more than the original.",
                a.value,
                form(a.signed)
            ),
            (None, Some(b)) => write!(
                f,
                "the Fortran drew {:?} and the port had stopped: the port \
                 draws fewer than the original.",
                b.value
            ),
            (None, None) => write!(f, "neither side has a draw here, which is a differ bug"),
        }
    }
}

fn form(signed: bool) -> &'static str {
    if signed { "signed" } else { "unit" }
}

/// Compare two traces, returning where they first disagree.
///
/// `None` means they are identical, which is what Tier 3 requires.
///
/// Values are compared *by bits*, not within a tolerance. A draw is a value
/// the generator produced, not an arithmetic result, and two implementations
/// of the same sequence either agree exactly or are not the same sequence.
#[must_use]
pub fn diff(ours: &[Draw], theirs: &[Draw]) -> Option<Divergence> {
    for index in 0..ours.len().max(theirs.len()) {
        let a = ours.get(index).copied();
        let b = theirs.get(index).copied();
        let same = match (a, b) {
            (Some(x), Some(y)) => x.signed == y.signed && x.value.to_bits() == y.value.to_bits(),
            _ => false,
        };
        if !same {
            return Some(Divergence {
                index,
                ours: a,
                theirs: b,
            });
        }
    }
    None
}

#[cfg(feature = "oracle")]
pub use reader::{clear, trace};

#[cfg(feature = "oracle")]
mod reader {
    use tepsim_core::disturbance::Draw;

    use crate::TRACE_CAPACITY;
    use crate::ffi;

    /// Reset the Fortran's draw counter.
    ///
    /// Call before each evaluation. Takes `&mut Oracle` rather than being free
    /// so that the process-wide `COMMON` is only touched under the lock.
    pub fn clear(_oracle: &mut crate::Oracle) {
        // SAFETY: the caller holds the oracle lock, and `Rngtrc` mirrors the
        // layout `instrument.rs` declares. Writing the counter alone is enough
        // to reset the trace: the arrays are only read up to it.
        unsafe {
            (&raw mut ffi::rngtrc_.count).write(0);
        }
    }

    /// Everything the Fortran drew since the last [`clear`].
    ///
    /// # Panics
    ///
    /// If the evaluation made more draws than `COMMON/RNGTRC/` can hold. The
    /// counter keeps counting past the capacity precisely so that this is an
    /// error rather than a quietly truncated trace.
    pub fn trace(_oracle: &mut crate::Oracle) -> Vec<Draw> {
        // SAFETY: as above; the arrays are `[f64; N]` and `[i32; N]` written by
        // Fortran, read only up to the count.
        let count = unsafe { (&raw const ffi::rngtrc_.count).read() };
        let count = usize::try_from(count).expect("a non-negative draw count");
        assert!(
            count <= TRACE_CAPACITY,
            "the evaluation made {count} draws and COMMON/RNGTRC/ holds \
             {TRACE_CAPACITY}. The trace is truncated, so raise TRCCAP in \
             instrument.rs and TRACE_CAPACITY beside it; do not compare a \
             truncated trace."
        );

        let values = unsafe { (&raw const ffi::rngtrc_.value).read() };
        let signs = unsafe { (&raw const ffi::rngtrc_.sign).read() };
        (0..count)
            .map(|i| Draw {
                value: values[i],
                signed: signs[i] < 0,
            })
            .collect()
    }
}
