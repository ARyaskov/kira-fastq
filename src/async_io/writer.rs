use std::path::Path;
use std::pin::Pin;

use async_compression::tokio::write::GzipEncoder;
use tokio::fs::File;
use tokio::io::{AsyncWrite, AsyncWriteExt, BufWriter};

use crate::error::{FastqError, InvalidKind, UnsupportedOperation};
use crate::record::{FastqRecord, FastqRecordOwned};
use crate::simd::bases::validate_bases_with;
use crate::simd::newline::contains_line_break;
use crate::simd::qual::validate_qual_encoding;
use crate::validation::{Alphabet, QualityEncoding};
use crate::writer::{WriteValidation, assemble_record};

const DEFAULT_WRITE_BUF: usize = 1024 * 1024;
const DEFAULT_GZIP_LEVEL: u32 = 6;

/// Async FASTQ writer over any `AsyncWrite`.
///
/// Same discipline as the sync writer: one `write_all` per record, structural checks on every
/// record, opt-in content validation. Call [`AsyncFastqWriter::shutdown`] before dropping a
/// compressed writer, otherwise the trailer is never written.
pub struct AsyncFastqWriter<W: AsyncWrite + Unpin + Send> {
    inner: W,
    scratch: Vec<u8>,
    validation: WriteValidation,
    alphabet: Alphabet,
    quality: QualityEncoding,
}

impl<W: AsyncWrite + Unpin + Send> AsyncFastqWriter<W> {
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

    #[inline]
    pub fn with_quality_encoding(mut self, encoding: QualityEncoding) -> Self {
        self.quality = encoding;
        self
    }

    pub async fn write_record(&mut self, rec: &FastqRecord<'_>) -> Result<(), FastqError> {
        self.write_parts(rec.header(), rec.seq(), rec.qual()).await
    }

    pub async fn write_record_owned(&mut self, rec: &FastqRecordOwned) -> Result<(), FastqError> {
        self.write_parts(rec.header(), rec.seq(), rec.qual()).await
    }

    pub async fn write_parts(
        &mut self,
        header: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), FastqError> {
        self.validate(header, seq, qual)?;
        assemble_record(&mut self.scratch, header, seq, qual);
        self.inner
            .write_all(&self.scratch)
            .await
            .map_err(FastqError::Io)
    }

    /// Flush buffered bytes. For gzip output the trailer is still pending; see
    /// [`AsyncFastqWriter::shutdown`].
    pub async fn flush(&mut self) -> Result<(), FastqError> {
        self.inner.flush().await.map_err(FastqError::Io)
    }

    /// Final flush and close. Required for compressed output to be valid.
    pub async fn shutdown(&mut self) -> Result<(), FastqError> {
        self.inner.shutdown().await.map_err(FastqError::Io)
    }

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

/// Type-erased async writer returned by the path constructor.
pub type BoxedAsyncWriter = AsyncFastqWriter<Pin<Box<dyn AsyncWrite + Send + Unpin>>>;

impl AsyncFastqWriter<Pin<Box<dyn AsyncWrite + Send + Unpin>>> {
    /// Open an output path. `.gz` selects streaming gzip. BGZF and zstd are not available on
    /// the async path; use the sync writer, which supports both.
    pub async fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        if is_bgzf(path) {
            return Err(FastqError::Unsupported(UnsupportedOperation::AsyncBgzf));
        }
        if is_zst(path) {
            return Err(FastqError::Unsupported(UnsupportedOperation::Zstd));
        }
        let file = File::create(path).await.map_err(FastqError::Io)?;
        let buffered = BufWriter::with_capacity(DEFAULT_WRITE_BUF, file);
        if is_gz(path) {
            let encoder = GzipEncoder::with_quality(
                buffered,
                async_compression::Level::Precise(DEFAULT_GZIP_LEVEL as i32),
            );
            let boxed: Pin<Box<dyn AsyncWrite + Send + Unpin>> = Box::pin(encoder);
            Ok(Self::from_writer(boxed))
        } else {
            let boxed: Pin<Box<dyn AsyncWrite + Send + Unpin>> = Box::pin(buffered);
            Ok(Self::from_writer(boxed))
        }
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
