//! Gzip input: mmap over the compressed file, streaming inflate into a reusable window.
//!
//! Three properties this backend guarantees, all of which mainstream FASTQ tooling relies on:
//!
//! 1. **Multi-member streams are fully decoded.** `bgzip`, `pigz` and plain `cat a.gz b.gz`
//!    all produce concatenated members; stopping at the first one silently truncates the data.
//! 2. **CRC32 and ISIZE are always verified.** Corruption in a sequencing run should surface as
//!    an error, not as mangled bases. The check costs a few percent next to inflate.
//! 3. **A stream that ends mid-member is an error**, never a clean end of file.
//!
//! Inflate runs through `flate2` with the `zlib-rs` backend, which benchmarked ~1.7x faster than
//! the `miniz_oxide` implementation this backend used previously.

use std::fs::File;
use std::path::Path;

use flate2::{Decompress, FlushDecompress, Status};
use memmap2::{Mmap, MmapOptions};

use crate::backend::{ByteSource, LineStatus, read_line_from};
use crate::error::{FastqError, InvalidKind};

/// Decoded-byte window. Large enough to amortise inflate calls, small enough to stay in L2.
const DEFAULT_OUT_BUF_SIZE: usize = 512 * 1024;

pub(crate) struct GzipBackend {
    mmap: Mmap,
    src_pos: usize,
    buf: Vec<u8>,
    buf_pos: usize,
    buf_len: usize,
    /// Decoded bytes already dropped from the front of `buf`.
    dropped: u64,
    decomp: Decompress,
    need_header: bool,
    finished: bool,
    crc: crc32fast::Hasher,
    isize: u32,
}

impl GzipBackend {
    pub(crate) fn new(path: &Path) -> Result<Self, FastqError> {
        Self::new_with_buf_size(path, DEFAULT_OUT_BUF_SIZE)
    }

    pub(crate) fn new_with_buf_size(path: &Path, out_buf_size: usize) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        let metadata = file.metadata()?;
        if metadata.len() == 0 {
            return Err(FastqError::invalid(0, InvalidKind::GzipHeader));
        }
        // SAFETY: the mapping spans the whole file; see `MmapBackend::open`.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(Self {
            mmap,
            src_pos: 0,
            buf: vec![0u8; out_buf_size.max(64 * 1024)],
            buf_pos: 0,
            buf_len: 0,
            dropped: 0,
            decomp: Decompress::new(false),
            need_header: true,
            finished: false,
            crc: crc32fast::Hasher::new(),
            isize: 0,
        })
    }

    #[inline]
    pub(crate) fn logical_offset(&self) -> u64 {
        self.dropped + self.buf_pos as u64
    }

    #[inline]
    pub(crate) fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        read_line_from(self, out)
    }

    /// Drop the consumed prefix of the window so the next inflate call has room.
    #[inline]
    fn compact(&mut self) {
        if self.buf_pos == 0 {
            return;
        }
        let remaining = self.buf_len - self.buf_pos;
        if remaining > 0 {
            self.buf.copy_within(self.buf_pos..self.buf_len, 0);
        }
        self.dropped += self.buf_pos as u64;
        self.buf_len = remaining;
        self.buf_pos = 0;
    }

    /// Position the source at the next member, or finish the stream.
    ///
    /// Trailing NUL padding is skipped and any other trailing garbage ends the stream, which is
    /// what `gzip(1)` and `zcat` do.
    fn start_member(&mut self) -> Result<bool, FastqError> {
        while self.src_pos < self.mmap.len() && self.mmap[self.src_pos] == 0 {
            self.src_pos += 1;
        }
        if self.src_pos >= self.mmap.len() {
            return Ok(false);
        }
        if self.mmap.len() - self.src_pos < 2
            || self.mmap[self.src_pos] != 0x1f
            || self.mmap[self.src_pos + 1] != 0x8b
        {
            return Ok(false);
        }
        let header_end = self.parse_header_at(self.src_pos)?;
        self.src_pos = header_end;
        self.decomp.reset(false);
        self.crc = crc32fast::Hasher::new();
        self.isize = 0;
        Ok(true)
    }

    /// Verify the 8-byte member trailer: CRC32 then ISIZE, both little-endian.
    fn finish_member(&mut self) -> Result<(), FastqError> {
        let offset = self.logical_offset();
        let end = self
            .src_pos
            .checked_add(8)
            .ok_or_else(|| FastqError::eof(offset))?;
        if end > self.mmap.len() {
            return Err(FastqError::invalid(offset, InvalidKind::GzipTruncated));
        }
        let trailer = &self.mmap[self.src_pos..end];
        let expected_crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
        let expected_isize = u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
        let actual_crc = std::mem::replace(&mut self.crc, crc32fast::Hasher::new()).finalize();
        if actual_crc != expected_crc {
            return Err(FastqError::invalid(offset, InvalidKind::GzipTrailerCrc));
        }
        if self.isize != expected_isize {
            return Err(FastqError::invalid(offset, InvalidKind::GzipTrailerIsize));
        }
        self.src_pos = end;
        self.need_header = true;
        Ok(())
    }

    fn parse_header_at(&self, start: usize) -> Result<usize, FastqError> {
        let offset = self.logical_offset();
        let bytes = &self.mmap[start..];
        let bad = || FastqError::invalid(offset, InvalidKind::GzipHeader);
        if bytes.len() < 10 || bytes[0] != 0x1f || bytes[1] != 0x8b || bytes[2] != 8 {
            return Err(bad());
        }
        let flg = bytes[3];
        let mut i = 10usize;
        if flg & 0x04 != 0 {
            // FEXTRA
            let end = i.checked_add(2).ok_or_else(bad)?;
            if end > bytes.len() {
                return Err(bad());
            }
            let xlen = u16::from_le_bytes([bytes[i], bytes[i + 1]]) as usize;
            i = end.checked_add(xlen).ok_or_else(bad)?;
            if i > bytes.len() {
                return Err(bad());
            }
        }
        for flag in [0x08u8, 0x10] {
            // FNAME, FCOMMENT: NUL-terminated
            if flg & flag != 0 {
                while i < bytes.len() && bytes[i] != 0 {
                    i += 1;
                }
                if i >= bytes.len() {
                    return Err(bad());
                }
                i += 1;
            }
        }
        if flg & 0x02 != 0 {
            // FHCRC
            i = i.checked_add(2).ok_or_else(bad)?;
            if i > bytes.len() {
                return Err(bad());
            }
        }
        Ok(start.saturating_add(i))
    }
}

impl ByteSource for GzipBackend {
    #[inline]
    fn available(&self) -> &[u8] {
        if self.buf_pos >= self.buf_len {
            &[]
        } else {
            &self.buf[self.buf_pos..self.buf_len]
        }
    }

    #[inline]
    fn consume(&mut self, n: usize) {
        self.buf_pos = self.buf_pos.saturating_add(n).min(self.buf_len);
    }

    fn refill(&mut self) -> Result<bool, FastqError> {
        loop {
            if self.finished {
                return Ok(false);
            }
            if self.need_header {
                if !self.start_member()? {
                    self.finished = true;
                    return Ok(false);
                }
                self.need_header = false;
            }

            self.compact();
            if self.buf_len >= self.buf.len() {
                return Err(FastqError::invalid(
                    self.logical_offset(),
                    InvalidKind::BufferOverflow,
                ));
            }

            let before_in = self.decomp.total_in();
            let before_out = self.decomp.total_out();
            let status = self
                .decomp
                .decompress(
                    &self.mmap[self.src_pos..],
                    &mut self.buf[self.buf_len..],
                    FlushDecompress::None,
                )
                .map_err(|_| FastqError::invalid(self.logical_offset(), InvalidKind::GzipData))?;
            let consumed_in = (self.decomp.total_in() - before_in) as usize;
            let produced = (self.decomp.total_out() - before_out) as usize;
            self.src_pos += consumed_in;

            if produced > 0 {
                let fresh = &self.buf[self.buf_len..self.buf_len + produced];
                self.crc.update(fresh);
                self.isize = self.isize.wrapping_add(produced as u32);
                self.buf_len += produced;
            }

            match status {
                Status::StreamEnd => {
                    self.finish_member()?;
                    if self.buf_pos < self.buf_len {
                        return Ok(true);
                    }
                    continue;
                }
                Status::Ok | Status::BufError => {
                    if consumed_in == 0 && produced == 0 {
                        // No progress: either the member is cut short or the stream is corrupt.
                        let kind = if self.src_pos >= self.mmap.len() {
                            InvalidKind::GzipTruncated
                        } else {
                            InvalidKind::GzipData
                        };
                        return Err(FastqError::invalid(self.logical_offset(), kind));
                    }
                    if self.buf_pos < self.buf_len {
                        return Ok(true);
                    }
                    continue;
                }
            }
        }
    }
}
