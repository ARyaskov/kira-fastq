use std::path::Path;
use std::pin::Pin;

use async_compression::tokio::write::GzipEncoder;
use tokio::fs::File;
use tokio::io::{AsyncWrite, AsyncWriteExt, BufWriter};

use crate::error::FastqError;
use crate::record::{FastqRecord, FastqRecordOwned};
use crate::simd::bases::validate_bases_with;
use crate::simd::qual::validate_qual;
use crate::validation::Alphabet;
use crate::writer::{WriteValidation, assemble_record};

const DEFAULT_WRITE_BUF: usize = 1024 * 1024;
const DEFAULT_GZIP_LEVEL: u32 = 6;

/// Async FASTQ writer over any `AsyncWrite`. Same single-syscall-per-record discipline
/// as the sync writer: each record is assembled into a scratch buffer and emitted with
/// one `write_all().await`.
pub struct AsyncFastqWriter<W: AsyncWrite + Unpin + Send> {
    inner: W,
    scratch: Vec<u8>,
    validation: WriteValidation,
    alphabet: Alphabet,
}

impl<W: AsyncWrite + Unpin + Send> AsyncFastqWriter<W> {
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

    pub async fn write_record(&mut self, rec: &FastqRecord<'_>) -> Result<(), FastqError> {
        self.validate(rec.seq(), rec.qual())?;
        assemble_record(&mut self.scratch, rec.header(), rec.seq(), rec.qual());
        self.inner
            .write_all(&self.scratch)
            .await
            .map_err(FastqError::Io)
    }

    pub async fn write_record_owned(
        &mut self,
        rec: &FastqRecordOwned,
    ) -> Result<(), FastqError> {
        self.write_record(&rec.as_borrowed()).await
    }

    pub async fn write_parts(
        &mut self,
        header: &[u8],
        seq: &[u8],
        qual: &[u8],
    ) -> Result<(), FastqError> {
        self.validate(seq, qual)?;
        assemble_record(&mut self.scratch, header, seq, qual);
        self.inner
            .write_all(&self.scratch)
            .await
            .map_err(FastqError::Io)
    }

    /// Flush the encoder/buffer. For gzip output you also want [`Self::shutdown`] at the
    /// end to flush the final deflate block + trailer.
    pub async fn flush(&mut self) -> Result<(), FastqError> {
        self.inner.flush().await.map_err(FastqError::Io)
    }

    /// Final flush + close. Required for gzip/BGZF output to be valid.
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

/// Type-erased async writer (returned by path-based constructors). Hides the gzip-vs-plain
/// generic split behind a `Pin<Box<dyn AsyncWrite>>`.
pub type BoxedAsyncWriter = AsyncFastqWriter<Pin<Box<dyn AsyncWrite + Send + Unpin>>>;

impl AsyncFastqWriter<Pin<Box<dyn AsyncWrite + Send + Unpin>>> {
    /// Open an output path. `.gz` triggers gzip via `async-compression`. BGZF is not
    /// supported on this async path (use the sync writer or wrap a `noodles_bgzf::AsyncWriter`).
    pub async fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        if is_bgzf(path) {
            return Err(FastqError::Io(std::io::Error::new(
                std::io::ErrorKind::Unsupported,
                "BGZF async output is not supported; use sync writer or wrap noodles_bgzf::AsyncWriter via from_writer",
            )));
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
