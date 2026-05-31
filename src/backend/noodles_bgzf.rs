//! Optional `noodles-bgzf` adapter (enabled via the `noodles-bgzf` feature).
//!
//! Wraps `noodles_bgzf::Reader<File>` and exposes it through kira's [`LineSource`] trait.
//! Same FASTQ-over-BGZF semantics as [`crate::backend::bgzf::BgzfBackend`], but the inflate
//! is delegated to the noodles implementation — useful when the rest of your pipeline shares
//! noodles BGZF semantics (virtual offsets compatible with `noodles-bam` / `noodles-vcf`).
//!
//! Trade-off vs. kira's own BGZF backend: noodles uses streaming reads instead of mmap, so
//! the hot LF scan still benefits from kira's SIMD, but block I/O goes through `Read`.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use noodles_bgzf::io::Reader as BgzfReader;

use crate::backend::gzip::LineStatus;
use crate::error::FastqError;
use crate::simd::newline::find_lf;

/// Inflate buffer: large enough to amortize syscalls; small enough to stay in L2.
const REFILL_TARGET: usize = 256 * 1024;

pub struct NoodlesBgzfBackend {
    inner: BgzfReader<BufReader<File>>,
    buf: Vec<u8>,
    buf_pos: usize,
    buf_len: usize,
    logical_offset: u64,
    eof: bool,
}

impl NoodlesBgzfBackend {
    pub fn open(path: &Path) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        let reader = BgzfReader::new(BufReader::with_capacity(64 * 1024, file));
        Ok(Self {
            inner: reader,
            buf: vec![0u8; REFILL_TARGET],
            buf_pos: 0,
            buf_len: 0,
            logical_offset: 0,
            eof: false,
        })
    }

    /// Construct from a pre-built `noodles_bgzf::Reader`.
    pub fn from_reader(reader: BgzfReader<BufReader<File>>) -> Self {
        Self {
            inner: reader,
            buf: vec![0u8; REFILL_TARGET],
            buf_pos: 0,
            buf_len: 0,
            logical_offset: 0,
            eof: false,
        }
    }

    #[inline]
    pub fn logical_offset(&self) -> u64 {
        self.logical_offset + self.buf_pos as u64
    }

    pub fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        out.clear();
        loop {
            let slice = self.available_slice();
            if !slice.is_empty() {
                if let Some(lf) = find_lf(slice, 0) {
                    let mut end = lf;
                    if end > 0 && slice[end - 1] == b'\r' {
                        end -= 1;
                    }
                    out.extend_from_slice(&slice[..end]);
                    self.advance(lf + 1);
                    return Ok(LineStatus::Line);
                }
                let n = slice.len();
                out.extend_from_slice(slice);
                self.advance(n);
            }
            if !self.refill()? {
                if out.is_empty() {
                    return Ok(LineStatus::EofClean);
                }
                if out.last() == Some(&b'\r') {
                    out.pop();
                }
                return Ok(LineStatus::EofPartial);
            }
        }
    }

    #[inline]
    fn available_slice(&self) -> &[u8] {
        if self.buf_pos >= self.buf_len {
            &[]
        } else {
            &self.buf[self.buf_pos..self.buf_len]
        }
    }

    #[inline]
    fn advance(&mut self, n: usize) {
        self.buf_pos = self.buf_pos.saturating_add(n).min(self.buf_len);
    }

    fn refill(&mut self) -> Result<bool, FastqError> {
        if self.eof {
            return Ok(false);
        }
        // Account for bytes already yielded before resetting the scratch.
        self.logical_offset += self.buf_len as u64;
        self.buf_pos = 0;
        self.buf_len = 0;
        // Pull a chunk via BufRead::fill_buf — avoids double-buffering decoded bytes.
        let filled = self.inner.fill_buf().map_err(FastqError::Io)?;
        if filled.is_empty() {
            self.eof = true;
            return Ok(false);
        }
        let n = filled.len().min(self.buf.len());
        self.buf[..n].copy_from_slice(&filled[..n]);
        self.inner.consume(n);
        self.buf_len = n;
        Ok(true)
    }
}
