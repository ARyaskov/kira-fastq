use std::path::Path;
use std::pin::Pin;

use async_compression::tokio::bufread::GzipDecoder;
use futures_util::stream::{Stream, try_unfold};
use tokio::fs::File;
use tokio::io::{AsyncBufRead, AsyncBufReadExt, AsyncReadExt, BufReader};

use crate::backend::gzip::LineStatus;
use crate::error::{FastqError, InvalidKind};
use crate::format::FastqFormat;
use crate::parser::{RecordScratch, Segment};
use crate::record::{FastqRecord, FastqRecordOwned};
use crate::simd::bases::validate_bases_with;
use crate::simd::qual::validate_qual;
use crate::validation::{Alphabet, ValidationMode};

/// Async FASTQ reader over any `AsyncBufRead`. See [`crate::async_io`] module docs for the
/// performance contract.
pub struct AsyncFastqReader<R: AsyncBufRead + Unpin + Send> {
    inner: R,
    scratch: RecordScratch,
    seq_segs: Vec<Segment>,
    qual_segs: Vec<Segment>,
    validation: ValidationMode,
    alphabet: Alphabet,
    format: FastqFormat,
    logical_offset: u64,
}

impl<R: AsyncBufRead + Unpin + Send> AsyncFastqReader<R> {
    #[inline]
    pub fn from_reader(inner: R) -> Self {
        Self {
            inner,
            scratch: RecordScratch::new(),
            seq_segs: Vec::new(),
            qual_segs: Vec::new(),
            validation: ValidationMode::None,
            alphabet: Alphabet::default(),
            format: FastqFormat::default(),
            logical_offset: 0,
        }
    }

    #[inline]
    pub fn with_validation(mut self, mode: ValidationMode) -> Self {
        self.validation = mode;
        self
    }

    #[inline]
    pub fn with_alphabet(mut self, alphabet: Alphabet) -> Self {
        self.alphabet = alphabet;
        self
    }

    #[inline]
    pub fn with_format(mut self, format: FastqFormat) -> Self {
        self.format = format;
        self
    }

    /// Logical byte offset of the next-to-be-read line, counting decoded bytes.
    #[inline]
    pub fn logical_offset(&self) -> u64 {
        self.logical_offset
    }

    /// Read the next record. Borrowed into the reader's scratch buffer — the borrow is
    /// invalidated on the next call.
    pub async fn next(&mut self) -> Result<Option<FastqRecord<'_>>, FastqError> {
        match self.format {
            FastqFormat::SingleLine => self.next_single().await,
            FastqFormat::MultiLine => self.next_multi().await,
        }
    }

    async fn next_single(&mut self) -> Result<Option<FastqRecord<'_>>, FastqError> {
        let header_start = self.logical_offset;
        match read_line(&mut self.inner, &mut self.scratch.header, &mut self.logical_offset).await?
        {
            LineStatus::Line => {}
            LineStatus::EofClean => return Ok(None),
            LineStatus::EofPartial => {
                return Err(FastqError::UnexpectedEof {
                    offset: header_start,
                });
            }
        }
        if self.scratch.header.is_empty() || self.scratch.header[0] != b'@' {
            return Err(FastqError::InvalidFormat {
                offset: header_start,
                kind: InvalidKind::HeaderMissingAt,
            });
        }
        self.scratch.header.drain(..1);

        let seq_start = self.logical_offset;
        if read_line(&mut self.inner, &mut self.scratch.seq, &mut self.logical_offset).await?
            != LineStatus::Line
        {
            return Err(FastqError::UnexpectedEof { offset: seq_start });
        }
        if self.scratch.seq.is_empty() {
            return Err(FastqError::InvalidFormat {
                offset: seq_start,
                kind: InvalidKind::SeqLineEmpty,
            });
        }

        let plus_start = self.logical_offset;
        if read_line(&mut self.inner, &mut self.scratch.plus, &mut self.logical_offset).await?
            != LineStatus::Line
        {
            return Err(FastqError::UnexpectedEof { offset: plus_start });
        }
        if self.scratch.plus.is_empty() || self.scratch.plus[0] != b'+' {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::PlusMissing,
            });
        }

        let qual_start = self.logical_offset;
        if read_line(&mut self.inner, &mut self.scratch.qual, &mut self.logical_offset).await?
            != LineStatus::Line
        {
            return Err(FastqError::UnexpectedEof { offset: qual_start });
        }

        if self.scratch.seq.len() != self.scratch.qual.len() {
            return Err(FastqError::LengthMismatch {
                offset: qual_start,
                seq_len: self.scratch.seq.len(),
                qual_len: self.scratch.qual.len(),
            });
        }

        validate(
            self.validation,
            self.alphabet,
            &self.scratch.seq,
            &self.scratch.qual,
            seq_start,
            qual_start,
        )?;

        Ok(Some(FastqRecord::new(
            &self.scratch.header,
            &self.scratch.seq,
            &self.scratch.qual,
        )))
    }

    async fn next_multi(&mut self) -> Result<Option<FastqRecord<'_>>, FastqError> {
        let header_start = self.logical_offset;
        match read_line(&mut self.inner, &mut self.scratch.header, &mut self.logical_offset).await?
        {
            LineStatus::Line => {}
            LineStatus::EofClean => return Ok(None),
            LineStatus::EofPartial => {
                return Err(FastqError::UnexpectedEof {
                    offset: header_start,
                });
            }
        }
        if self.scratch.header.is_empty() || self.scratch.header[0] != b'@' {
            return Err(FastqError::InvalidFormat {
                offset: header_start,
                kind: InvalidKind::HeaderMissingAt,
            });
        }
        self.scratch.header.drain(..1);

        self.scratch.seq.clear();
        self.seq_segs.clear();
        let seq_start = self.logical_offset;
        let plus_start;
        let mut tmp = std::mem::take(&mut self.scratch.plus);
        loop {
            let line_start = self.logical_offset;
            tmp.clear();
            match read_line(&mut self.inner, &mut tmp, &mut self.logical_offset).await? {
                LineStatus::Line => {}
                LineStatus::EofClean | LineStatus::EofPartial => {
                    self.scratch.plus = tmp;
                    return Err(FastqError::UnexpectedEof { offset: line_start });
                }
            }
            if tmp.is_empty() {
                self.scratch.plus = tmp;
                return Err(FastqError::InvalidFormat {
                    offset: line_start,
                    kind: InvalidKind::SeqLineEmpty,
                });
            }
            if tmp[0] == b'+' {
                plus_start = line_start;
                break;
            }
            self.scratch.seq.extend_from_slice(&tmp);
            self.seq_segs.push(Segment {
                offset: line_start,
                len: tmp.len(),
            });
        }
        self.scratch.plus = tmp;

        if self.scratch.seq.is_empty() {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::SeqLineEmpty,
            });
        }

        let qual_start = self.logical_offset;
        self.scratch.qual.clear();
        self.qual_segs.clear();
        let mut remaining = self.scratch.seq.len();
        let mut tmp = std::mem::take(&mut self.scratch.plus);
        while remaining > 0 {
            let line_start = self.logical_offset;
            tmp.clear();
            match read_line(&mut self.inner, &mut tmp, &mut self.logical_offset).await? {
                LineStatus::Line => {}
                LineStatus::EofClean | LineStatus::EofPartial => {
                    self.scratch.plus = tmp;
                    return Err(FastqError::UnexpectedEof { offset: line_start });
                }
            }
            if tmp.is_empty() {
                self.scratch.plus = tmp;
                return Err(FastqError::InvalidFormat {
                    offset: line_start,
                    kind: InvalidKind::QualLineEmpty,
                });
            }
            if tmp.len() > remaining {
                let qlen = self.scratch.qual.len() + tmp.len();
                self.scratch.plus = tmp;
                return Err(FastqError::LengthMismatch {
                    offset: line_start,
                    seq_len: self.scratch.seq.len(),
                    qual_len: qlen,
                });
            }
            self.scratch.qual.extend_from_slice(&tmp);
            self.qual_segs.push(Segment {
                offset: line_start,
                len: tmp.len(),
            });
            remaining -= tmp.len();
        }
        self.scratch.plus = tmp;

        validate(
            self.validation,
            self.alphabet,
            &self.scratch.seq,
            &self.scratch.qual,
            seq_start,
            qual_start,
        )?;

        Ok(Some(FastqRecord::new(
            &self.scratch.header,
            &self.scratch.seq,
            &self.scratch.qual,
        )))
    }
}

impl<R: AsyncBufRead + Unpin + Send + 'static> AsyncFastqReader<R> {
    /// Convert the reader into a [`Stream`] of owned records. Use when records must
    /// outlive `next()` calls (channel sends, async pipelines, axum responses).
    /// Per-record allocation is the cost of moving from borrowed to owned.
    pub fn records(self) -> RecordStream<R> {
        RecordStream::new(self)
    }
}

/// Path-based async reader: an enum that dispatches between plain and gzip variants
/// at runtime. Keeps each variant statically typed (no `Box<dyn AsyncBufRead>`).
///
/// Construct via [`AnyAsyncReader::from_path`].
//
// The variants differ by ~250 B (gzip stack carries a `GzipDecoder` plus an extra
// `BufReader`); boxing would just add an indirection on the per-record hot path.
#[allow(clippy::large_enum_variant)]
pub enum AnyAsyncReader {
    Plain(AsyncFastqReader<BufReader<File>>),
    Gzip(AsyncFastqReader<BufReader<GzipDecoder<BufReader<File>>>>),
}

impl AnyAsyncReader {
    /// Open a path with magic-byte sniffing for gzip. BGZF is decoded as gzip on this path
    /// (no virtual-offset semantics).
    pub async fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        let mut file = File::open(path).await.map_err(FastqError::Io)?;
        let mut magic = [0u8; 2];
        let peeked = match file.read_exact(&mut magic).await {
            Ok(_) => true,
            Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => false,
            Err(e) => return Err(FastqError::Io(e)),
        };
        drop(file);
        // We re-open to reset the stream position rather than keeping a holding buffer —
        // a single extra `open(2)` is negligible compared to the per-record cost.
        let file = File::open(path).await.map_err(FastqError::Io)?;
        let buffered = BufReader::with_capacity(256 * 1024, file);
        if peeked && magic == [0x1f, 0x8b] {
            let decoder = GzipDecoder::new(buffered);
            let buffered = BufReader::with_capacity(256 * 1024, decoder);
            return Ok(AnyAsyncReader::Gzip(AsyncFastqReader::from_reader(
                buffered,
            )));
        }
        Ok(AnyAsyncReader::Plain(AsyncFastqReader::from_reader(
            buffered,
        )))
    }

    pub async fn next(&mut self) -> Result<Option<FastqRecord<'_>>, FastqError> {
        match self {
            Self::Plain(r) => r.next().await,
            Self::Gzip(r) => r.next().await,
        }
    }

    pub fn with_validation(self, mode: ValidationMode) -> Self {
        match self {
            Self::Plain(r) => Self::Plain(r.with_validation(mode)),
            Self::Gzip(r) => Self::Gzip(r.with_validation(mode)),
        }
    }

    pub fn with_alphabet(self, alphabet: Alphabet) -> Self {
        match self {
            Self::Plain(r) => Self::Plain(r.with_alphabet(alphabet)),
            Self::Gzip(r) => Self::Gzip(r.with_alphabet(alphabet)),
        }
    }

    pub fn with_format(self, format: FastqFormat) -> Self {
        match self {
            Self::Plain(r) => Self::Plain(r.with_format(format)),
            Self::Gzip(r) => Self::Gzip(r.with_format(format)),
        }
    }
}

/// Stream wrapper yielding owned records. Constructed via
/// [`AsyncFastqReader::records`].
pub struct RecordStream<R: AsyncBufRead + Unpin + Send + 'static> {
    #[allow(clippy::type_complexity)]
    inner: Pin<
        Box<
            dyn Stream<Item = Result<FastqRecordOwned, FastqError>>
                + Send
                + 'static,
        >,
    >,
    _marker: std::marker::PhantomData<R>,
}

impl<R: AsyncBufRead + Unpin + Send + 'static> RecordStream<R> {
    fn new(reader: AsyncFastqReader<R>) -> Self {
        let stream = try_unfold(reader, |mut reader| async move {
            match reader.next().await {
                Ok(Some(rec)) => {
                    let owned = rec.to_owned();
                    Ok(Some((owned, reader)))
                }
                Ok(None) => Ok(None),
                Err(e) => Err(e),
            }
        });
        Self {
            inner: Box::pin(stream),
            _marker: std::marker::PhantomData,
        }
    }
}

impl<R: AsyncBufRead + Unpin + Send + 'static> Stream for RecordStream<R> {
    type Item = Result<FastqRecordOwned, FastqError>;

    #[inline]
    fn poll_next(
        mut self: Pin<&mut Self>,
        cx: &mut std::task::Context<'_>,
    ) -> std::task::Poll<Option<Self::Item>> {
        self.inner.as_mut().poll_next(cx)
    }
}

/// Internal async line reader. Uses tokio's `read_until` (which uses `memchr::memchr`
/// under the hood — still SIMD-y, just not our AVX-512 path; the syscall cost dominates
/// async I/O anyway).
async fn read_line<R: AsyncBufRead + Unpin>(
    reader: &mut R,
    out: &mut Vec<u8>,
    logical_offset: &mut u64,
) -> Result<LineStatus, FastqError> {
    out.clear();
    let n = reader.read_until(b'\n', out).await.map_err(FastqError::Io)?;
    if n == 0 {
        return Ok(LineStatus::EofClean);
    }
    *logical_offset += n as u64;
    if out.last() == Some(&b'\n') {
        out.pop();
        if out.last() == Some(&b'\r') {
            out.pop();
        }
        return Ok(LineStatus::Line);
    }
    if out.last() == Some(&b'\r') {
        out.pop();
    }
    Ok(LineStatus::EofPartial)
}

#[inline]
fn validate(
    mode: ValidationMode,
    alphabet: Alphabet,
    seq: &[u8],
    qual: &[u8],
    seq_start: u64,
    qual_start: u64,
) -> Result<(), FastqError> {
    match mode {
        ValidationMode::None => Ok(()),
        ValidationMode::Bases => check_bases(seq, alphabet, seq_start),
        ValidationMode::Qualities => check_qual(qual, qual_start),
        ValidationMode::BasesAndQualities => {
            check_bases(seq, alphabet, seq_start)?;
            check_qual(qual, qual_start)
        }
    }
}

#[inline]
fn check_bases(seq: &[u8], alphabet: Alphabet, base: u64) -> Result<(), FastqError> {
    match validate_bases_with(seq, alphabet) {
        Ok(()) => Ok(()),
        Err(idx) => Err(FastqError::InvalidBase {
            offset: base + idx as u64,
            byte: seq[idx],
        }),
    }
}

#[inline]
fn check_qual(qual: &[u8], base: u64) -> Result<(), FastqError> {
    match validate_qual(qual) {
        Ok(()) => Ok(()),
        Err(idx) => Err(FastqError::InvalidQuality {
            offset: base + idx as u64,
            byte: qual[idx],
        }),
    }
}
