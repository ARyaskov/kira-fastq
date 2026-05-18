use std::fs::File;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};

use crate::error::FastqError;

/// Caveat: the mapping is a snapshot at `open()`. If the file is truncated by another
/// process while the mapping is live, accessing removed pages raises `SIGBUS` on Unix or
/// an access violation on Windows.
pub struct MmapBackend {
    mmap: Option<Mmap>,
}

impl MmapBackend {
    pub fn open(path: &Path) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        let metadata = file.metadata()?;
        if metadata.len() == 0 {
            return Ok(Self { mmap: None });
        }
        // SAFETY: file is kept alive for the duration of the mmap; mapping spans the whole file.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(Self { mmap: Some(mmap) })
    }

    pub fn advise_sequential(&self) {
        #[cfg(unix)]
        {
            if let Some(m) = &self.mmap {
                let _ = m.advise(memmap2::Advice::Sequential);
            }
        }
    }

    #[inline]
    pub fn slice_from(&self, pos: usize) -> &[u8] {
        match &self.mmap {
            Some(m) if pos < m.len() => &m[pos..],
            _ => &[],
        }
    }

    #[inline]
    pub fn as_slice(&self) -> &[u8] {
        match &self.mmap {
            Some(m) => &m[..],
            None => &[],
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.mmap.as_ref().map_or(0, |m| m.len())
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
