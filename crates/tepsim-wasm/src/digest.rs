//! A small, exact, dependency-free fingerprint over a stream of `f64`s.
//!
//! # Why not BLAKE3
//!
//! Tier 9 compares BLAKE3 digests of full validation runs across x86-64,
//! aarch64 and wasm32, and that comparison belongs to `xtask`, which can afford
//! the dependency. What the browser needs is smaller: something a page can
//! print next to a native run to show the two agree, cheap enough to compute on
//! every chunk, and small enough not to spend the wasm size budget on a hash.
//!
//! FNV-1a over `f64::to_bits` is integer arithmetic end to end, so the digest
//! cannot itself become a source of the cross-platform disagreement it exists
//! to detect. It is not cryptographic and not collision-resistant against an
//! adversary. It catches drift, which is the failure this project actually has.

/// The 64-bit FNV-1a hash of a stream of `f64` bit patterns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fnv1a64(u64);

impl Fnv1a64 {
    /// The FNV-1a 64-bit offset basis.
    const OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    /// The FNV-1a 64-bit prime.
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    /// An empty digest.
    #[must_use]
    pub const fn new() -> Self {
        Self(Self::OFFSET_BASIS)
    }

    /// Absorb one byte.
    #[inline]
    pub const fn write_u8(&mut self, byte: u8) {
        self.0 ^= byte as u64;
        self.0 = self.0.wrapping_mul(Self::PRIME);
    }

    /// Absorb the eight little-endian bytes of a value's bit pattern.
    ///
    /// Bits, not the value, so byte order is fixed by IEEE 754 rather than by
    /// the host, and so `-0.0` stays distinguishable from `0.0`.
    #[inline]
    pub const fn write_f64(&mut self, value: f64) {
        let bytes = value.to_bits().to_le_bytes();
        let mut i = 0;
        while i < bytes.len() {
            self.write_u8(bytes[i]);
            i += 1;
        }
    }

    /// Absorb every value of a slice, in order.
    #[inline]
    pub fn write_slice(&mut self, values: &[f64]) {
        for value in values {
            self.write_f64(*value);
        }
    }

    /// Absorb a boolean as one byte.
    #[inline]
    pub const fn write_bool(&mut self, value: bool) {
        self.write_u8(value as u8);
    }

    /// Absorb the eight little-endian bytes of an integer.
    #[inline]
    pub const fn write_u64(&mut self, value: u64) {
        let bytes = value.to_le_bytes();
        let mut i = 0;
        while i < bytes.len() {
            self.write_u8(bytes[i]);
            i += 1;
        }
    }

    /// The digest so far.
    #[must_use]
    pub const fn finish(self) -> u64 {
        self.0
    }
}

impl Default for Fnv1a64 {
    fn default() -> Self {
        Self::new()
    }
}

/// Lower-case hex, fixed width, so two digests line up when printed one above
/// the other. That is the only way anybody actually compares them.
#[must_use]
pub fn hex64(value: u64) -> String {
    format!("{value:016x}")
}
