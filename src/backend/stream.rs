//! Streaming backend over a user-provided [`std::io::BufRead`].
//!
//! Unlike the mmap/gzip/bgzf backends which control their own input layout, this backend
//! adapts an arbitrary `BufRead` source — stdin, a TCP socket, a custom decoder, a noodles
//! reader, anything that yields bytes. There is **no mmap** on this path, so per-record
//! lines are copied into the reader's scratch buffer. SIMD-LF scan still applies.

use std::io::BufRead;

use crate::backend::gzip::LineStatus;
use crate::error::FastqError;
use crate::simd::newline::find_lf;

/// Adapter over `Box<dyn BufRead>`. Implements [`crate::multiline::LineSource`].
pub struct StreamBackend {
    inner: Box<dyn BufRead + Send>,
    logical_offset: u64,
}

impl StreamBackend {
    #[inline]
    pub fn new(inner: Box<dyn BufRead + Send>) -> Self {
        Self {
            inner,
            logical_offset: 0,
        }
    }

    #[inline]
    pub fn logical_offset(&self) -> u64 {
        self.logical_offset
    }

    /// Reads bytes up to and including the next LF, strips trailing CR/LF, and pushes the
    /// payload into `out`. Tracks bytes consumed (including line terminator) in
    /// `logical_offset`.
    pub fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        out.clear();
        loop {
            let available = match self.inner.fill_buf() {
                Ok(b) => b,
                Err(e) => return Err(FastqError::Io(e)),
            };

            if available.is_empty() {
                // EOF
                if out.is_empty() {
                    return Ok(LineStatus::EofClean);
                }
                if out.last() == Some(&b'\r') {
                    out.pop();
                }
                return Ok(LineStatus::EofPartial);
            }

            if let Some(lf) = find_lf(available, 0) {
                let mut end = lf;
                if end > 0 && available[end - 1] == b'\r' {
                    end -= 1;
                }
                out.extend_from_slice(&available[..end]);
                let consumed = lf + 1;
                self.inner.consume(consumed);
                self.logical_offset += consumed as u64;
                return Ok(LineStatus::Line);
            }

            let n = available.len();
            out.extend_from_slice(available);
            self.inner.consume(n);
            self.logical_offset += n as u64;
        }
    }

    /// Peek the next available byte without consuming it.
    pub fn peek_byte(&mut self) -> Result<Option<u8>, FastqError> {
        let available = self.inner.fill_buf().map_err(FastqError::Io)?;
        Ok(available.first().copied())
    }
}
