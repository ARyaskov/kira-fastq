//! Stream positions.

/// A position in an input stream.
///
/// The meaning depends on the backend the offset came from, and offsets from different backends
/// are not interchangeable:
///
/// * **BGZF** — a htslib-style virtual offset: the compressed offset of the block in the high
///   48 bits, the offset inside the decompressed block in the low 16. Use
///   [`VirtualOffset::new`], [`VirtualOffset::compressed`] and [`VirtualOffset::uncompressed`]
///   rather than doing the arithmetic by hand.
/// * **Plain files and in-memory buffers** — a plain byte offset in the low 64 bits.
/// * **gzip and arbitrary streams** — the number of decoded bytes consumed so far. Reported by
///   `tell()` for progress purposes; `seek()` is not supported on those sources.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct VirtualOffset(pub u64);

impl VirtualOffset {
    /// Build a BGZF virtual offset from a block's compressed offset and an offset inside that
    /// block's decompressed data.
    #[inline]
    pub const fn new(compressed: u64, uncompressed: u16) -> Self {
        Self((compressed << 16) | (uncompressed as u64))
    }

    /// Compressed offset of the block, for BGZF offsets.
    #[inline]
    pub const fn compressed(&self) -> u64 {
        self.0 >> 16
    }

    /// Offset inside the decompressed block, for BGZF offsets.
    #[inline]
    pub const fn uncompressed(&self) -> u16 {
        (self.0 & 0xFFFF) as u16
    }

    /// Raw value: a byte offset for plain and stream sources, the packed pair for BGZF.
    #[inline]
    pub const fn get(&self) -> u64 {
        self.0
    }
}

impl From<u64> for VirtualOffset {
    #[inline]
    fn from(value: u64) -> Self {
        Self(value)
    }
}

impl From<VirtualOffset> for u64 {
    #[inline]
    fn from(value: VirtualOffset) -> Self {
        value.0
    }
}

#[cfg(feature = "noodles-bgzf")]
impl From<noodles_bgzf::VirtualPosition> for VirtualOffset {
    #[inline]
    fn from(value: noodles_bgzf::VirtualPosition) -> Self {
        Self(u64::from(value))
    }
}

#[cfg(feature = "noodles-bgzf")]
impl From<VirtualOffset> for noodles_bgzf::VirtualPosition {
    #[inline]
    fn from(value: VirtualOffset) -> Self {
        noodles_bgzf::VirtualPosition::from(value.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packs_and_unpacks() {
        let voff = VirtualOffset::new(4096, 1234);
        assert_eq!(voff.compressed(), 4096);
        assert_eq!(voff.uncompressed(), 1234);
        assert_eq!(voff, VirtualOffset((4096 << 16) | 1234));
    }
}
