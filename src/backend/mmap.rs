use std::fs::File;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};

use crate::error::FastqError;

pub struct MmapBackend {
    mmap: Mmap,
}

impl MmapBackend {
    pub fn open(path: &Path) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        // SAFETY: file is kept alive for the duration of the mmap; mapping spans the whole file.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(Self { mmap })
    }

    #[inline]
    pub fn slice_from(&self, pos: usize) -> &[u8] {
        if pos >= self.mmap.len() {
            &[]
        } else {
            &self.mmap[pos..]
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.mmap.len()
    }
}
