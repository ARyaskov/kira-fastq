//! Optional `noodles-bgzf` adapter (feature `noodles-bgzf`).
//!
//! Wraps `noodles_bgzf::io::Reader<BufReader<File>>` behind kira's line reader. Use it when the
//! rest of a pipeline already shares noodles' BGZF semantics, e.g. when virtual offsets travel
//! between this reader and `noodles-bam` or `noodles-vcf`. Unlike kira's own BGZF backend it
//! reads through `Read` instead of a mapping, but it supports the same `tell`/`seek` pair and
//! its offsets convert to and from `noodles_bgzf::VirtualPosition`.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles_bgzf::io::Reader as BgzfReader;

use crate::backend::{LineStatus, read_line_bufread};
use crate::error::FastqError;
use crate::offset::VirtualOffset;

pub(crate) struct NoodlesBgzfBackend {
    inner: BgzfReader<BufReader<File>>,
    logical_offset: u64,
}

impl NoodlesBgzfBackend {
    pub(crate) fn open(path: &Path) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        let reader = BgzfReader::new(BufReader::with_capacity(64 * 1024, file));
        Ok(Self::from_reader(reader))
    }

    /// Construct from a pre-built `noodles_bgzf` reader.
    pub(crate) fn from_reader(reader: BgzfReader<BufReader<File>>) -> Self {
        Self {
            inner: reader,
            logical_offset: 0,
        }
    }

    #[inline]
    pub(crate) fn logical_offset(&self) -> u64 {
        self.logical_offset
    }

    /// Virtual offset of the next byte, in noodles' encoding (identical to htslib's).
    #[inline]
    pub(crate) fn tell(&self) -> VirtualOffset {
        VirtualOffset::from(self.inner.virtual_position())
    }

    pub(crate) fn seek(&mut self, voff: VirtualOffset) -> Result<(), FastqError> {
        self.inner
            .seek(voff.into())
            .map_err(FastqError::Io)
            .map(|_| ())
    }

    #[inline]
    pub(crate) fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        read_line_bufread(&mut self.inner, out, &mut self.logical_offset)
    }
}
