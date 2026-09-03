//! Contiguous input: a memory-mapped file or an owned buffer.
//!
//! Both are parsed in place by the slice parser, so records borrow directly out of the source
//! with no copying.

use crate::backend::mmap::MmapBackend;

pub(crate) enum ContiguousBackend {
    Mmap(MmapBackend),
    /// Owned bytes: `FastqReader::from_vec`, and the buffered file reader that reads a whole
    /// file up front instead of mapping it.
    Owned(Vec<u8>),
}

impl ContiguousBackend {
    #[inline]
    pub(crate) fn as_slice(&self) -> &[u8] {
        match self {
            ContiguousBackend::Mmap(m) => m.as_slice(),
            ContiguousBackend::Owned(v) => v,
        }
    }

    #[inline]
    pub(crate) fn len(&self) -> usize {
        self.as_slice().len()
    }
}
