//! FASTQ writer.
//!
//! Design, in priority order:
//!
//! 1. **One `write_all` per record.** Each record is assembled into a reusable scratch buffer
//!    (`@header\nseq\n+\nqual\n`) and handed to the sink in one call. The scratch grows to the
//!    largest record seen and is then reused.
//!
//! 2. **Records that cannot round-trip are rejected.** Sequence and quality length must agree
//!    and the header must not contain a line break, always, on every write. Those checks are
//!    O(1) and O(header) respectively; without them the writer can emit a file that no reader,
//!    including this crate's own, can parse back. Scanning the sequence and quality for line
//!    breaks costs a pass over the payload and is opt-in through [`WriteValidation`].
//!
//! 3. **Opt-in SIMD content validation.** Alphabet and quality-range checks reuse the read
//!    side's vector kernels.
//!
//! 4. **Errors at close are reported.** Compressed formats write a trailer when they finish;
//!    [`FastqWriter::finish`] surfaces failures there instead of swallowing them in `Drop`.
//!
//! ## Output formats
//!
//! Chosen from the file extension: `.gz` is gzip, `.bgz`/`.bgzf` is BGZF, `.zst` is zstd with
//! the `zstd` feature, everything else is plain text. Output always uses LF line endings.

mod bgzf;
mod record;

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use flate2::write::GzEncoder;

#[cfg_attr(feature = "zstd", allow(unused_imports))]
use crate::error::UnsupportedOperation;
use crate::error::{FastqError, InvalidKind};
use crate::record::{FastqRecord, FastqRecordOwned};
use crate::simd::bases::validate_bases_with;
use crate::simd::newline::contains_line_break;
use crate::simd::qual::validate_qual_encoding;
use crate::validation::{Alphabet, QualityEncoding};

pub use self::bgzf::{BgzfWriter, ParallelBgzfWriter};
pub use self::record::assemble_record;

use self::bgzf::checked_level;

/// Default kernel-write buffer, big enough to absorb tens of thousands of short reads.
const DEFAULT_WRITE_BUF: usize = 1024 * 1024;

/// Default gzip level. `6` matches `gzip(1)` and trades CPU for size sanely.
const DEFAULT_GZIP_LEVEL: u32 = 6;

/// Default zstd level, matching the `zstd(1)` default.
#[cfg(feature = "zstd")]
const DEFAULT_ZSTD_LEVEL: i32 = 3;

/// Optional pre-write validation.
///
/// Length agreement and header line breaks are checked whatever this is set to; the variants
/// below add checks that cost a pass over the payload.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum WriteValidation {
    /// Structural checks only.
    #[default]
    None,
    /// Also reject line breaks inside the sequence or the qualities.
    LineBreaks,
    /// Line breaks plus the base alphabet.
    Bases,
    /// Line breaks plus the quality range.
    Qualities,
    /// Everything.
    BasesAndQualities,
}

/// Streaming FASTQ writer over a generic [`Write`].
pub struct FastqWriter<W: Write> {
    inner: W,
    scratch: Vec<u8>,
    validation: WriteValidation,
    alphabet: Alphabet,
    quality: QualityEncoding,
}

impl<W: Write> FastqWriter<W> {
    /// Build a writer over an arbitrary sink. Does **not** add buffering; wrap the sink in a
    /// [`BufWriter`] first for hot loops.
    #[inline]
    pub fn from_writer(inner: W) -> Self {
        Self {
            inner,
            scratch: Vec::with_capacity(8 * 1024),
            validation: WriteValidation::None,
            alphabet: Alphabet::default(),
            quality: QualityEncoding::default(),
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

    /// Quality encoding used by [`WriteValidation::Qualities`]. Defaults to Phred+33.
    #[inline]
    pub fn with_quality_encoding(mut self, encoding: QualityEncoding) -> Self {
        self.quality = encoding;
        self
    }

    /// Write a borrowed record.
    #[inline]
    pub fn write_record(&mut self, rec: &FastqRecord<'_>) -> Result<(), FastqError> {
        self.write_parts(rec.header(), rec.seq(), rec.qual())
    }

    /// Write an owned record. Same path as [`FastqWriter::write_record`].
    #[inline]
    pub fn write_record_owned(&mut self, rec: &FastqRecordOwned) -> Result<(), FastqError> {
        self.write_parts(rec.header(), rec.seq(), rec.qual())
    }

    /// Write a record from raw parts, for producers that already hold the bytes (BAM to FASTQ
    /// extraction, format conversion) and never build a record type.
    pub fn write_parts(
        &mut self,
        header: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), FastqError> {
        self.validate(header, seq, qual)?;
        assemble_record(&mut self.scratch, header, seq, qual);
        self.inner.write_all(&self.scratch).map_err(FastqError::Io)
    }

    /// Flush the underlying writer. For compressed output this does not write the trailer; see
    /// [`FastqWriter::finish`].
    #[inline]
    pub fn flush(&mut self) -> Result<(), FastqError> {
        self.inner.flush().map_err(FastqError::Io)
    }

    /// Consume the writer and return the sink, without finalising compressed output.
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

    fn validate(&self, header: &[u8], seq: &[u8], qual: &[u8]) -> Result<(), FastqError> {
        // Always: a record that breaks these cannot be read back at all.
        if seq.len() != qual.len() {
            return Err(FastqError::length_mismatch(0, seq.len(), qual.len()));
        }
        if contains_line_break(header) {
            return Err(FastqError::invalid(0, InvalidKind::HeaderContainsNewline));
        }
        if self.validation == WriteValidation::None {
            return Ok(());
        }
        if contains_line_break(seq) {
            return Err(FastqError::invalid(0, InvalidKind::SeqContainsNewline));
        }
        if contains_line_break(qual) {
            return Err(FastqError::invalid(0, InvalidKind::QualContainsNewline));
        }
        match self.validation {
            WriteValidation::None | WriteValidation::LineBreaks => Ok(()),
            WriteValidation::Bases => check_bases(seq, self.alphabet),
            WriteValidation::Qualities => check_qual(qual, self.quality),
            WriteValidation::BasesAndQualities => {
                check_bases(seq, self.alphabet)?;
                check_qual(qual, self.quality)
            }
        }
    }
}

/// Sink used by the path-based constructors.
///
/// An enum rather than `Box<dyn Write>` so that closing a compressed file can run the right
/// finaliser and report its errors: a gzip trailer that fails to write on a full disk produces
/// a file that decompresses to truncated data.
pub enum FastqSink {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
    Bgzf(BgzfWriter<BufWriter<File>>),
    BgzfParallel(ParallelBgzfWriter<BufWriter<File>>),
    #[cfg(feature = "zstd")]
    Zstd(zstd::stream::write::Encoder<'static, BufWriter<File>>),
    #[cfg(feature = "noodles-bgzf")]
    NoodlesBgzf(noodles_bgzf::io::Writer<BufWriter<File>>),
    /// Any other sink, for [`FastqWriter::from_writer`] users who still want the enum's type.
    Boxed(Box<dyn Write + Send>),
}

impl Write for FastqSink {
    #[inline]
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            FastqSink::Plain(w) => w.write(buf),
            FastqSink::Gzip(w) => w.write(buf),
            FastqSink::Bgzf(w) => w.write(buf),
            FastqSink::BgzfParallel(w) => w.write(buf),
            #[cfg(feature = "zstd")]
            FastqSink::Zstd(w) => w.write(buf),
            #[cfg(feature = "noodles-bgzf")]
            FastqSink::NoodlesBgzf(w) => w.write(buf),
            FastqSink::Boxed(w) => w.write(buf),
        }
    }

    #[inline]
    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            FastqSink::Plain(w) => w.write_all(buf),
            FastqSink::Gzip(w) => w.write_all(buf),
            FastqSink::Bgzf(w) => w.write_all(buf),
            FastqSink::BgzfParallel(w) => w.write_all(buf),
            #[cfg(feature = "zstd")]
            FastqSink::Zstd(w) => w.write_all(buf),
            #[cfg(feature = "noodles-bgzf")]
            FastqSink::NoodlesBgzf(w) => w.write_all(buf),
            FastqSink::Boxed(w) => w.write_all(buf),
        }
    }

    #[inline]
    fn flush(&mut self) -> io::Result<()> {
        match self {
            FastqSink::Plain(w) => w.flush(),
            FastqSink::Gzip(w) => w.flush(),
            FastqSink::Bgzf(w) => w.flush(),
            FastqSink::BgzfParallel(w) => w.flush(),
            #[cfg(feature = "zstd")]
            FastqSink::Zstd(w) => w.flush(),
            #[cfg(feature = "noodles-bgzf")]
            FastqSink::NoodlesBgzf(w) => w.flush(),
            FastqSink::Boxed(w) => w.flush(),
        }
    }
}

impl FastqSink {
    /// Finish the stream: write whatever trailer the format needs, flush, and close.
    pub fn finish(self) -> Result<(), FastqError> {
        match self {
            FastqSink::Plain(mut w) => w.flush().map_err(FastqError::Io),
            FastqSink::Gzip(w) => {
                let mut inner = w.finish().map_err(FastqError::Io)?;
                inner.flush().map_err(FastqError::Io)
            }
            FastqSink::Bgzf(w) => {
                let mut inner = w.finish()?;
                inner.flush().map_err(FastqError::Io)
            }
            FastqSink::BgzfParallel(w) => {
                let mut inner = w.finish()?;
                inner.flush().map_err(FastqError::Io)
            }
            #[cfg(feature = "zstd")]
            FastqSink::Zstd(w) => {
                let mut inner = w.finish().map_err(FastqError::Io)?;
                inner.flush().map_err(FastqError::Io)
            }
            #[cfg(feature = "noodles-bgzf")]
            FastqSink::NoodlesBgzf(w) => {
                let mut inner = w.finish().map_err(FastqError::Io)?;
                inner.flush().map_err(FastqError::Io)
            }
            FastqSink::Boxed(mut w) => w.flush().map_err(FastqError::Io),
        }
    }
}

/// Writer returned by the path-based constructors.
pub type PathWriter = FastqWriter<FastqSink>;

/// Former name of [`PathWriter`], kept so 0.3 code keeps compiling.
pub type BoxedWriter = PathWriter;

impl FastqWriter<FastqSink> {
    /// Open an output path. The format comes from the extension:
    ///
    /// - `.gz` gzip
    /// - `.bgz` / `.bgzf` BGZF
    /// - `.zst` zstd (feature `zstd`)
    /// - anything else plain text
    ///
    /// All variants sit on a 1 MiB [`BufWriter`]. Call [`FastqWriter::finish`] when done.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        if is_bgzf(path) {
            return Self::to_bgzf_path(path, DEFAULT_GZIP_LEVEL);
        }
        if is_gz(path) {
            return Self::to_gz_path(path, DEFAULT_GZIP_LEVEL);
        }
        if is_zst(path) {
            #[cfg(feature = "zstd")]
            {
                return Self::to_zstd_path(path, DEFAULT_ZSTD_LEVEL);
            }
            #[cfg(not(feature = "zstd"))]
            {
                return Err(FastqError::Unsupported(UnsupportedOperation::Zstd));
            }
        }
        Self::to_plain_path(path, DEFAULT_WRITE_BUF)
    }

    /// Plain output with an explicit buffer size.
    pub fn to_plain_path<P: AsRef<Path>>(path: P, buf_size: usize) -> Result<Self, FastqError> {
        let file = File::create(path).map_err(FastqError::Io)?;
        Ok(Self::from_writer(FastqSink::Plain(
            BufWriter::with_capacity(buf_size, file),
        )))
    }

    /// Gzip output. `level` is 0 to 9; anything higher is rejected rather than panicking inside
    /// the codec.
    pub fn to_gz_path<P: AsRef<Path>>(path: P, level: u32) -> Result<Self, FastqError> {
        let level = checked_level(level)?;
        let file = File::create(path).map_err(FastqError::Io)?;
        let buffered = BufWriter::with_capacity(DEFAULT_WRITE_BUF, file);
        Ok(Self::from_writer(FastqSink::Gzip(GzEncoder::new(
            buffered, level,
        ))))
    }

    /// BGZF output, block-compressed and indexable, readable by htslib and `bgzip -d`.
    pub fn to_bgzf_path<P: AsRef<Path>>(path: P, level: u32) -> Result<Self, FastqError> {
        let file = File::create(path).map_err(FastqError::Io)?;
        let buffered = BufWriter::with_capacity(DEFAULT_WRITE_BUF, file);
        Ok(Self::from_writer(FastqSink::Bgzf(BgzfWriter::new(
            buffered, level,
        )?)))
    }

    /// BGZF output compressed on `threads` worker threads. Pass `0` to size the pool from the
    /// machine's parallelism. Output is identical to [`FastqWriter::to_bgzf_path`].
    pub fn to_bgzf_path_parallel<P: AsRef<Path>>(
        path: P,
        level: u32,
        threads: usize,
    ) -> Result<Self, FastqError> {
        let threads = if threads == 0 {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        } else {
            threads
        };
        let file = File::create(path).map_err(FastqError::Io)?;
        let buffered = BufWriter::with_capacity(DEFAULT_WRITE_BUF, file);
        Ok(Self::from_writer(FastqSink::BgzfParallel(
            ParallelBgzfWriter::new(buffered, level, threads)?,
        )))
    }

    /// zstd output (feature `zstd`). `level` follows `zstd(1)`, 1 to 22.
    ///
    /// [`FastqWriter::finish`] is not optional here: unlike gzip and BGZF, the zstd encoder
    /// writes nothing on drop, so a writer that is only dropped leaves a truncated file.
    #[cfg(feature = "zstd")]
    pub fn to_zstd_path<P: AsRef<Path>>(path: P, level: i32) -> Result<Self, FastqError> {
        let file = File::create(path).map_err(FastqError::Io)?;
        let buffered = BufWriter::with_capacity(DEFAULT_WRITE_BUF, file);
        let encoder = zstd::stream::write::Encoder::new(buffered, level).map_err(FastqError::Io)?;
        Ok(Self::from_writer(FastqSink::Zstd(encoder)))
    }

    /// BGZF output through the `noodles-bgzf` writer (feature `noodles-bgzf`).
    #[cfg(feature = "noodles-bgzf")]
    pub fn to_noodles_bgzf_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let file = File::create(path.as_ref()).map_err(FastqError::Io)?;
        let buffered = BufWriter::with_capacity(DEFAULT_WRITE_BUF, file);
        Ok(Self::from_writer(FastqSink::NoodlesBgzf(
            noodles_bgzf::io::Writer::new(buffered),
        )))
    }

    /// Finish the output: write the format's trailer, flush, and close.
    ///
    /// Dropping the writer also finalises the stream, but only this reports errors. Data written
    /// to a compressed file is not complete until one of the two happens.
    pub fn finish(self) -> Result<(), FastqError> {
        self.into_inner().finish()
    }
}

#[inline]
fn is_gz(path: &Path) -> bool {
    has_extension(path, &["gz"])
}

#[inline]
fn is_bgzf(path: &Path) -> bool {
    has_extension(path, &["bgz", "bgzf"])
}

#[inline]
fn is_zst(path: &Path) -> bool {
    has_extension(path, &["zst", "zstd"])
}

#[inline]
fn has_extension(path: &Path, wanted: &[&str]) -> bool {
    path.extension().is_some_and(|ext| {
        wanted
            .iter()
            .any(|candidate| ext.eq_ignore_ascii_case(candidate))
    })
}

#[inline]
fn check_bases(seq: &[u8], alphabet: Alphabet) -> Result<(), FastqError> {
    match validate_bases_with(seq, alphabet) {
        Ok(()) => Ok(()),
        Err(idx) => Err(FastqError::invalid_base(idx as u64, seq[idx])),
    }
}

#[inline]
fn check_qual(qual: &[u8], encoding: QualityEncoding) -> Result<(), FastqError> {
    match validate_qual_encoding(qual, encoding) {
        Ok(()) => Ok(()),
        Err(idx) => Err(FastqError::invalid_quality(idx as u64, qual[idx])),
    }
}
