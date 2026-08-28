//! A content hash for a scenario.
//!
//! # Why FNV-1a and not BLAKE3
//!
//! This runs in `no_std`, compiles to wasm32, and must not spend the browser
//! bundle's size budget. What it has to do is distinguish two scenarios that
//! differ, so that a dataset carrying its scenario's digest cannot be confused
//! with one produced by a different experiment. It is not defending against an
//! adversary constructing a collision, and saying so plainly is better than
//! implying a guarantee that is not being provided.
//!
//! Tier 9 compares BLAKE3 digests of whole validation runs across
//! architectures, and that comparison lives in `xtask`, which can afford the
//! dependency.
//!
//! FNV-1a runs over `f64::to_bits`, so the digest is integer arithmetic end to
//! end. A hash that used floating-point arithmetic could itself vary between
//! architectures, which would make it useless for the one job it has.
//!
//! # Negative zero
//!
//! `0.0` and `-0.0` are equal as numbers and have different bit patterns.
//! [`Digest::push_f64`] normalises negative zero to positive, so two scenarios
//! that describe the same experiment hash the same even if one was built by
//! arithmetic that produced a signed zero. `NaN` is normalised to a single
//! pattern for the same reason.

/// A 64-bit FNV-1a hash over a stream of values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Digest(u64);

impl Default for Digest {
    fn default() -> Self {
        Self::new()
    }
}

impl Digest {
    /// The FNV-1a 64-bit offset basis.
    const OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    /// The FNV-1a 64-bit prime.
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    /// An empty digest.
    #[must_use]
    pub const fn new() -> Self {
        Self(Self::OFFSET_BASIS)
    }

    /// The accumulated value.
    #[must_use]
    pub const fn value(self) -> u64 {
        self.0
    }

    /// Absorb one byte.
    pub const fn push_byte(&mut self, byte: u8) {
        self.0 ^= byte as u64;
        self.0 = self.0.wrapping_mul(Self::PRIME);
    }

    /// Absorb a `u64`, most significant byte first.
    pub const fn push_u64(&mut self, value: u64) {
        let mut shift = 56_i32;
        while shift >= 0 {
            self.push_byte((value >> shift) as u8);
            shift -= 8;
        }
    }

    /// Absorb a `usize`.
    pub const fn push_usize(&mut self, value: usize) {
        self.push_u64(value as u64);
    }

    /// Absorb an `f64` by its bit pattern.
    ///
    /// Negative zero is normalised to positive and every `NaN` to one pattern;
    /// see the module documentation.
    pub fn push_f64(&mut self, value: f64) {
        let bits = if value == 0.0 {
            0
        } else if value.is_nan() {
            0x7ff8_0000_0000_0000
        } else {
            value.to_bits()
        };
        self.push_u64(bits);
    }

    /// Absorb a `bool`.
    pub const fn push_bool(&mut self, value: bool) {
        self.push_byte(value as u8);
    }

    /// Absorb a string's bytes.
    pub const fn push_str(&mut self, value: &str) {
        let bytes = value.as_bytes();
        let mut index = 0;
        while index < bytes.len() {
            self.push_byte(bytes[index]);
            index += 1;
        }
        // Length as well, so "ab" then "c" cannot collide with "a" then "bc".
        self.push_usize(bytes.len());
    }
}
