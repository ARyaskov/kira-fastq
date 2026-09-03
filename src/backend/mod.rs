//! Input backends and the line-reading machinery shared between them.

pub(crate) mod bgzf;
pub(crate) mod gzip;
pub(crate) mod memory;
pub(crate) mod mmap;
pub(crate) mod parallel;
pub(crate) mod stream;

#[cfg(feature = "noodles-bgzf")]
pub(crate) mod noodles_bgzf;

use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::GzipBackend;
use crate::backend::memory::ContiguousBackend;
use crate::backend::stream::StreamBackend;
use crate::error::FastqError;
use crate::simd::newline::find_lf;

#[cfg(feature = "noodles-bgzf")]
use crate::backend::noodles_bgzf::NoodlesBgzfBackend;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LineStatus {
    Line,
    EofClean,
    /// EOF reached with bytes still pending, i.e. the last line has no terminator. Readers treat
    /// this as a complete line: the missing final newline is the single most common deviation in
    /// FASTQ files in the wild, and every mainstream tool accepts it.
    EofPartial,
}

// One backend per reader; boxing would just add indirection on the hot path.
pub(crate) enum Backend {
    /// A whole file mapped into memory, or an in-memory buffer. Parsed in place: records borrow
    /// straight out of the mapping, no copying at all.
    Plain(ContiguousBackend),
    Gzip(GzipBackend),
    Bgzf(BgzfBackend),
    /// Arbitrary `BufRead` source. No random access; lines are copied into the reader's scratch.
    Stream(StreamBackend),
    /// Optional `noodles-bgzf` adapter, for pipelines that share virtual-offset semantics with
    /// the rest of the noodles ecosystem.
    #[cfg(feature = "noodles-bgzf")]
    NoodlesBgzf(NoodlesBgzfBackend),
}

/// A source that hands out contiguous runs of decoded bytes.
///
/// Implementing this gives a backend `read_line`/`peek_byte` for free, which is why the four
/// streaming backends no longer carry four copies of the same scan-copy-advance loop.
pub(crate) trait ByteSource {
    /// Decoded bytes ready to be consumed. Empty when the buffer is drained.
    fn available(&self) -> &[u8];
    /// Mark `n` bytes of [`ByteSource::available`] as consumed.
    fn consume(&mut self, n: usize);
    /// Pull more bytes. Returns `false` at end of input.
    fn refill(&mut self) -> Result<bool, FastqError>;
}

/// Read one line into `out`, stripping the `\n` and an optional preceding `\r`.
///
/// Inlined on purpose: the per-line call overhead is a measurable share of the streaming read
/// paths, which move one line at a time.
#[inline]
pub(crate) fn read_line_from<S: ByteSource>(
    src: &mut S,
    out: &mut Vec<u8>,
) -> Result<LineStatus, FastqError> {
    out.clear();
    loop {
        let slice = src.available();
        if !slice.is_empty() {
            if let Some(lf) = find_lf(slice, 0) {
                let mut end = lf;
                if end > 0 && slice[end - 1] == b'\r' {
                    end -= 1;
                }
                out.extend_from_slice(&slice[..end]);
                src.consume(lf + 1);
                return Ok(LineStatus::Line);
            }
            out.extend_from_slice(slice);
            let n = slice.len();
            src.consume(n);
        }
        if !src.refill()? {
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

/// Same as [`read_line_from`] for sources that are already a [`std::io::BufRead`], which cannot
/// implement [`ByteSource`] because `fill_buf` needs `&mut self`. `offset` tracks bytes consumed
/// including the terminator.
#[inline]
pub(crate) fn read_line_bufread<R: std::io::BufRead + ?Sized>(
    reader: &mut R,
    out: &mut Vec<u8>,
    offset: &mut u64,
) -> Result<LineStatus, FastqError> {
    out.clear();
    loop {
        let available = reader.fill_buf().map_err(FastqError::Io)?;
        if available.is_empty() {
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
            reader.consume(consumed);
            *offset += consumed as u64;
            return Ok(LineStatus::Line);
        }
        let n = available.len();
        out.extend_from_slice(available);
        reader.consume(n);
        *offset += n as u64;
    }
}
