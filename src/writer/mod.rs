//! FASTQ writer.
//!
//! Hot-path design choices, in priority order:
//!
//! 1. **One syscall per record.** Each record is assembled into a reusable scratch buffer
//!    (`@header\nseq\n+\nqual\n`) and written with a single [`Write::write_all`]. This cuts
//!    the per-record syscall (or buffer-fill-check) count from four to one. The scratch
//!    grows monotonically up to the largest record seen, then is reused without
//!    reallocation.
//!
//! 2. **Large output buffer.** The path-based constructors wrap the destination in a 1 MiB
//!    [`BufWriter`]. This amortises kernel write syscalls across many records.
//!
//! 3. **Opt-in SIMD pre-write validation.** [`WriteValidation`] piggybacks on the same SIMD
//!    base/quality validators used on the read side. Off by default — write is hot.
//!
//! 4. **mmap-write is not provided.** Streaming FASTQ output has unknown final size, and
//!    growing an mmap on every record extension is slower than buffered streaming on every
//!    OS we support. For known-size outputs (FAI-driven exports, fixed BAM→FASTQ dumps),
//!    use a custom `W: Write` over a pre-allocated mmap and pass it to
//!    [`FastqWriter::from_writer`].
//!
//! ## Output formats
//!
//! - Plain — direct `BufWriter<File>`.
//! - Gzip — `flate2::write::GzEncoder` (zlib-rs backend; SIMD-aware deflate).
//! - BGZF — feature-gated via `noodles-bgzf`; see [`FastqWriter::to_noodles_bgzf_path`].

mod record;

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;

use crate::error::FastqError;
use crate::record::{FastqRecord, FastqRecordOwned};
use crate::simd::bases::validate_bases_with;
use crate::simd::qual::validate_qual;
use crate::validation::Alphabet;

pub use self::record::assemble_record;

/// Default kernel-write buffer. 1 MiB is enough to absorb tens of thousands of typical
/// short-read records and stays in the L2 of any current CPU.
const DEFAULT_WRITE_BUF: usize = 1024 * 1024;

/// Default gzip compression level. `6` matches `gzip(1)` and trades CPU for size sanely.
const DEFAULT_GZIP_LEVEL: u32 = 6;

/// Optional pre-write validation. Off by default; turn on to fail fast on malformed inputs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum WriteValidation {
    #[default]
    None,
    /// Validate the base alphabet via SIMD.
    Bases,
    /// Validate Phred+33 qualities via SIMD.
    Qualities,
    /// Both.
    BasesAndQualities,
}

/// Streaming FASTQ writer over a generic [`Write`].
///
/// Constructed via [`FastqWriter::from_path`] (auto-detects gzip/BGZF by extension) or
/// [`FastqWriter::from_writer`] (generic over any `W: Write`).
pub struct FastqWriter<W: Write> {
    inner: W,
    scratch: Vec<u8>,
    validation: WriteValidation,
    alphabet: Alphabet,
}

impl<W: Write> FastqWriter<W> {
    /// Build a writer over an arbitrary `Write`. Does **not** add buffering — wrap your
    /// writer in `BufWriter` first for hot loops.
    #[inline]
    pub fn from_writer(inner: W) -> Self {
        Self {
            inner,
            scratch: Vec::with_capacity(8 * 1024),
            validation: WriteValidation::None,
            alphabet: Alphabet::default(),
        }
    }

    #[inline]
    pub fn with_validation(mut self, mode: WriteValidation) -> Self {
        self.validation = mode;
        self
    }

    #[inline]
    pub fn with_alphabet(mut self, alphabet: Alphabet) -> Self {
        self.alphabet = alphabet;
        self
    }

    /// Write a borrowed record. Single `write_all` syscall (plus whatever the underlying
    /// `W` chooses to do).
    #[inline]
    pub fn write_record(&mut self, rec: &FastqRecord<'_>) -> Result<(), FastqError> {
        self.validate(rec.seq(), rec.qual())?;
        assemble_record(&mut self.scratch, rec.header(), rec.seq(), rec.qual());
        self.inner.write_all(&self.scratch).map_err(FastqError::Io)
    }

    /// Write an owned record. Same code path as [`Self::write_record`].
    #[inline]
    pub fn write_record_owned(&mut self, rec: &FastqRecordOwned) -> Result<(), FastqError> {
        self.write_record(&rec.as_borrowed())
    }

    /// Write a record from raw parts. Useful when the upstream produced bytes directly
    /// (BAM→FASTQ extraction, format conversion) without first constructing a record type.
    #[inline]
    pub fn write_parts(
        &mut self,
        header: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), FastqError> {
        self.validate(seq, qual)?;
        assemble_record(&mut self.scratch, header, seq, qual);
        self.inner.write_all(&self.scratch).map_err(FastqError::Io)
    }

    /// Flush the underlying writer.
    #[inline]
    pub fn flush(&mut self) -> Result<(), FastqError> {
        self.inner.flush().map_err(FastqError::Io)
    }

    /// Consume the writer and return the inner sink.
    #[inline]
    pub fn into_inner(self) -> W {
        self.inner
    }

    #[inline]
    pub fn get_ref(&self) -> &W {
        &self.inner
    }

    #[inline]
    pub fn get_mut(&mut self) -> &mut W {
        &mut self.inner
    }

    #[inline]
    fn validate(&self, seq: &[u8], qual: &[u8]) -> Result<(), FastqError> {
        match self.validation {
            WriteValidation::None => Ok(()),
            WriteValidation::Bases => check_bases(seq, self.alphabet),
            WriteValidation::Qualities => check_qual(qual),
            WriteValidation::BasesAndQualities => {
                check_bases(seq, self.alphabet)?;
                check_qual(qual)
            }
        }
    }
}

// No `Drop` impl — that would conflict with `into_inner` (cannot move out of `Drop`).
// Inner writers (BufWriter, GzEncoder, noodles_bgzf::Writer) all auto-flush/finish on
// their own Drop. Call `flush()` explicitly if you need to surface errors.

/// Type-erased writer returned by path-based constructors. The variant is selected by
/// extension or explicit constructor; downstream code treats them uniformly.
pub type BoxedWriter = FastqWriter<Box<dyn Write + Send>>;

impl FastqWriter<Box<dyn Write + Send>> {
    /// Open an output path. Format is chosen by extension:
    ///
    /// - `.gz` → gzip (`flate2`/zlib-rs)
    /// - `.bgz` / `.bgzf` → BGZF (requires `noodles-bgzf` feature)
    /// - anything else → plain
    ///
    /// All paths get a 1 MiB [`BufWriter`].
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        if is_bgzf(path) {
            #[cfg(feature = "noodles-bgzf")]
            {
                return Self::to_noodles_bgzf_path(path);
            }
            #[cfg(not(feature = "noodles-bgzf"))]
            {
                return Err(FastqError::Io(std::io::Error::new(
                    std::io::ErrorKind::Unsupported,
                    "BGZF output requires the `noodles-bgzf` feature",
                )));
            }
        }
        if is_gz(path) {
            return Self::to_gz_path(path, DEFAULT_GZIP_LEVEL);
        }
        let file = File::create(path).map_err(FastqError::Io)?;
        let buf: Box<dyn Write + Send> = Box::new(BufWriter::with_capacity(DEFAULT_WRITE_BUF, file));
        Ok(Self::from_writer(buf))
    }

    /// Plain-text output to a path with explicit buffer size.
    pub fn to_plain_path<P: AsRef<Path>>(path: P, buf_size: usize) -> Result<Self, FastqError> {
        let file = File::create(path).map_err(FastqError::Io)?;
        let buf: Box<dyn Write + Send> = Box::new(BufWriter::with_capacity(buf_size, file));
        Ok(Self::from_writer(buf))
    }

    /// Gzip-compressed output. Level is a `flate2::Compression` raw value (0–9).
    pub fn to_gz_path<P: AsRef<Path>>(path: P, level: u32) -> Result<Self, FastqError> {
        let file = File::create(path).map_err(FastqError::Io)?;
        let buffered = BufWriter::with_capacity(DEFAULT_WRITE_BUF, file);
        let encoder = GzEncoder::new(buffered, Compression::new(level));
        let boxed: Box<dyn Write + Send> = Box::new(encoder);
        Ok(Self::from_writer(boxed))
    }

    /// BGZF output via the optional `noodles-bgzf` writer.
    #[cfg(feature = "noodles-bgzf")]
    pub fn to_noodles_bgzf_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let file = File::create(path.as_ref()).map_err(FastqError::Io)?;
        let buffered = BufWriter::with_capacity(DEFAULT_WRITE_BUF, file);
        let writer = noodles_bgzf::io::Writer::new(buffered);
        let boxed: Box<dyn Write + Send> = Box::new(writer);
        Ok(Self::from_writer(boxed))
    }
}

#[inline]
fn is_gz(path: &Path) -> bool {
    path.extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

#[inline]
fn is_bgzf(path: &Path) -> bool {
    path.extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case("bgz") || ext.eq_ignore_ascii_case("bgzf"))
}

#[inline]
fn check_bases(seq: &[u8], alphabet: Alphabet) -> Result<(), FastqError> {
    match validate_bases_with(seq, alphabet) {
        Ok(()) => Ok(()),
        Err(idx) => Err(FastqError::InvalidBase {
            offset: idx as u64,
            byte: seq[idx],
        }),
    }
}

#[inline]
fn check_qual(qual: &[u8]) -> Result<(), FastqError> {
    match validate_qual(qual) {
        Ok(()) => Ok(()),
        Err(idx) => Err(FastqError::InvalidQuality {
            offset: idx as u64,
            byte: qual[idx],
        }),
    }
}
