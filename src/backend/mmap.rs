use std::fs::File;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};

use crate::error::FastqError;

/// A whole file mapped into the address space.
///
/// Caveat: the mapping is a snapshot taken at `open()`. If another process truncates the file
/// while the mapping is live, touching the removed pages raises `SIGBUS` on Unix or an access
/// violation on Windows. Read a file that is still being written with
/// [`crate::FastqReader::from_path_buffered`] instead.
pub(crate) struct MmapBackend {
    mmap: Option<Mmap>,
}

impl MmapBackend {
    pub(crate) fn open(path: &Path) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        let metadata = file.metadata()?;
        if metadata.len() == 0 {
            return Ok(Self { mmap: None });
        }
        // SAFETY: the mapping spans the whole file and the file handle stays alive for the
        // duration of the call; `memmap2` keeps the mapping valid afterwards.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(Self { mmap: Some(mmap) })
    }

    /// Hint the kernel that access will be sequential. No-op off Unix.
    pub(crate) fn advise_sequential(&self) {
        #[cfg(unix)]
        {
            if let Some(m) = &self.mmap {
                let _ = m.advise(memmap2::Advice::Sequential);
            }
        }
    }

    #[inline]
    pub(crate) fn as_slice(&self) -> &[u8] {
        match &self.mmap {
            Some(m) => &m[..],
            None => &[],
        }
    }
}
