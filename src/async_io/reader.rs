use std::path::Path;
use std::pin::Pin;

use async_compression::tokio::bufread::GzipDecoder;
use futures_util::stream::{Stream, try_unfold};
use tokio::fs::File;
use tokio::io::{AsyncBufRead, AsyncBufReadExt, AsyncReadExt, BufReader};

use crate::backend::LineStatus;
use crate::error::{FastqError, InvalidKind, UnsupportedOperation};
use crate::format::FastqFormat;
use crate::parser::{RecordScratch, Segment};
use crate::record::{FastqRecord, FastqRecordOwned};
use crate::simd::bases::validate_bases_with;
use crate::simd::qual::validate_qual_encoding;
use crate::validation::{Alphabet, QualityEncoding, ValidationMode};

/// Async FASTQ reader over any `AsyncBufRead`.
///
/// Accepts the same deviations as the sync reader: a missing final newline, blank lines between
/// records, zero-length reads, and CRLF.
///
/// # Cancellation
///
/// [`AsyncFastqReader::next`] is **not** cancel-safe. A record spans four reads of the
/// underlying stream, so dropping the future partway leaves the reader between lines and the
/// next call will misparse. Do not call it in a `tokio::select!` branch that can lose the race;
/// drive the reader from one task, and use [`AsyncFastqReader::records`] with a channel if other
/// tasks need the records.
pub struct AsyncFastqReader<R: AsyncBufRead + Unpin + Send> {
    inner: R,
    scratch: RecordScratch,
    seq_segs: Vec<Segment>,
    qual_segs: Vec<Segment>,
    validation: ValidationMode,
    alphabet: Alphabet,
    quality: QualityEncoding,
    format: FastqFormat,
    logical_offset: u64,
    records: u64,
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
            quality: QualityEncoding::default(),
            format: FastqFormat::default(),
            logical_offset: 0,
            records: 0,
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
    pub fn with_quality_encoding(mut self, encoding: QualityEncoding) -> Self {
        self.quality = encoding;
        self
    }

    #[inline]
    pub fn with_format(mut self, format: FastqFormat) -> Self {
        self.format = format;
        self
    }

    /// Logical byte offset of the next line, counting decoded bytes.
    #[inline]
    pub fn logical_offset(&self) -> u64 {
        self.logical_offset
    }

    /// Number of records returned so far.
    #[inline]
    pub fn records_read(&self) -> u64 {
        self.records
    }

    /// Read the next record, borrowed from the reader's scratch. The borrow ends at the next
    /// call. Not cancel-safe; see the type's documentation.
    pub async fn next(&mut self) -> Result<Option<FastqRecord<'_>>, FastqError> {
        let record_index = self.records + 1;
        let found = match self.format {
            FastqFormat::SingleLine => self.next_single(record_index).await,
            FastqFormat::MultiLine => self.next_multi(record_index).await,
        }
        .map_err(|e| e.with_record(record_index))?;
        if !found {
            return Ok(None);
        }
        self.records = record_index;
        Ok(Some(FastqRecord::new(
            &self.scratch.header,
            &self.scratch.seq,
            &self.scratch.qual,
        )))
    }

    /// Read a header line, skipping blank lines. `false` means clean end of input.
    async fn read_header(&mut self) -> Result<Option<u64>, FastqError> {
        loop {
            let start = self.logical_offset;
            match read_line(
                &mut self.inner,
                &mut self.scratch.header,
                &mut self.logical_offset,
            )
            .await?
            {
                LineStatus::Line | LineStatus::EofPartial => {}
                LineStatus::EofClean => return Ok(None),
            }
            if self.scratch.header.is_empty() {
                continue;
            }
            if self.scratch.header[0] != b'@' {
                return Err(FastqError::invalid(start, InvalidKind::HeaderMissingAt));
            }
            self.scratch.header.drain(..1);
            return Ok(Some(start));
        }
    }

    async fn next_single(&mut self, _record: u64) -> Result<bool, FastqError> {
        if self.read_header().await?.is_none() {
            return Ok(false);
        }

        let seq_start = self.logical_offset;
        if read_line(
            &mut self.inner,
            &mut self.scratch.seq,
            &mut self.logical_offset,
        )
        .await?
            == LineStatus::EofClean
        {
            return Err(FastqError::eof(seq_start));
        }

        let plus_start = self.logical_offset;
        if read_line(
            &mut self.inner,
            &mut self.scratch.plus,
            &mut self.logical_offset,
        )
        .await?
            == LineStatus::EofClean
        {
            return Err(FastqError::eof(plus_start));
        }
        if self.scratch.plus.first() != Some(&b'+') {
            return Err(FastqError::invalid(plus_start, InvalidKind::PlusMissing));
        }

        let qual_start = self.logical_offset;
        if read_line(
            &mut self.inner,
            &mut self.scratch.qual,
            &mut self.logical_offset,
        )
        .await?
            == LineStatus::EofClean
        {
            return Err(FastqError::eof(qual_start));
        }

        if self.scratch.seq.len() != self.scratch.qual.len() {
            return Err(FastqError::length_mismatch(
                qual_start,
                self.scratch.seq.len(),
                self.scratch.qual.len(),
            ));
        }

        validate(
            self.validation,
            self.alphabet,
            self.quality,
            &self.scratch.seq,
            &self.scratch.qual,
            seq_start,
            qual_start,
        )?;
        Ok(true)
    }

    async fn next_multi(&mut self, _record: u64) -> Result<bool, FastqError> {
        if self.read_header().await?.is_none() {
            return Ok(false);
        }

        self.scratch.seq.clear();
        self.seq_segs.clear();
        let seq_start = self.logical_offset;
        let mut tmp = std::mem::take(&mut self.scratch.plus);
        let outcome = self.read_body(&mut tmp, seq_start).await;
        self.scratch.plus = tmp;
        let qual_start = outcome?;

        validate(
            self.validation,
            self.alphabet,
            self.quality,
            &self.scratch.seq,
            &self.scratch.qual,
            seq_start,
            qual_start,
        )?;
        Ok(true)
    }

    /// Sequence lines up to the `+` line, then quality lines. Returns the quality offset.
    async fn read_body(&mut self, tmp: &mut Vec<u8>, _seq_start: u64) -> Result<u64, FastqError> {
        loop {
            let line_start = self.logical_offset;
            if read_line(&mut self.inner, tmp, &mut self.logical_offset).await?
                == LineStatus::EofClean
            {
                return Err(FastqError::eof(line_start));
            }
            if tmp.first() == Some(&b'+') {
                break;
            }
            if tmp.is_empty() {
                return Err(FastqError::invalid(line_start, InvalidKind::SeqLineEmpty));
            }
            self.scratch.seq.extend_from_slice(tmp);
            self.seq_segs.push(Segment {
                offset: line_start,
                len: tmp.len(),
            });
        }

        let qual_start = self.logical_offset;
        self.scratch.qual.clear();
        self.qual_segs.clear();
        let mut remaining = self.scratch.seq.len();
        if remaining == 0 {
            if read_line(&mut self.inner, tmp, &mut self.logical_offset).await?
                == LineStatus::EofClean
            {
                return Err(FastqError::eof(qual_start));
            }
            if !tmp.is_empty() {
                return Err(FastqError::length_mismatch(qual_start, 0, tmp.len()));
            }
            return Ok(qual_start);
        }
        while remaining > 0 {
            let line_start = self.logical_offset;
            if read_line(&mut self.inner, tmp, &mut self.logical_offset).await?
                == LineStatus::EofClean
            {
                return Err(FastqError::eof(line_start));
            }
            if tmp.is_empty() {
                return Err(FastqError::invalid(line_start, InvalidKind::QualLineEmpty));
            }
            if tmp.len() > remaining {
                return Err(FastqError::length_mismatch(
                    line_start,
                    self.scratch.seq.len(),
                    self.scratch.qual.len() + tmp.len(),
                ));
            }
            self.scratch.qual.extend_from_slice(tmp);
            self.qual_segs.push(Segment {
                offset: line_start,
                len: tmp.len(),
            });
            remaining -= tmp.len();
        }
        Ok(qual_start)
    }
}

impl<R: AsyncBufRead + Unpin + Send + 'static> AsyncFastqReader<R> {
    /// Convert into a [`Stream`](futures_util::Stream) of owned records, for channel sends and
    /// async pipelines. One allocation per field per record is the price of owning them.
    pub fn records(self) -> RecordStream<R> {
        RecordStream::new(self)
    }
}

/// Path-based async reader: plain or gzip, chosen by magic bytes, each variant statically typed.
//
// The variants differ by a few hundred bytes; boxing would add an indirection per record.
#[allow(clippy::large_enum_variant)]
pub enum AnyAsyncReader {
    Plain(AsyncFastqReader<BufReader<File>>),
    Gzip(AsyncFastqReader<BufReader<GzipDecoder<BufReader<File>>>>),
}

impl AnyAsyncReader {
    /// Open a path, sniffing gzip from the magic bytes.
    ///
    /// BGZF is decoded here as ordinary gzip, which is valid but drops virtual-offset
    /// semantics: multi-member decoding is enabled, so every block is read. For true BGZF
    /// semantics use the sync reader or wrap a `noodles_bgzf` async reader with
    /// [`AsyncFastqReader::from_reader`].
    pub async fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        let mut file = File::open(path).await.map_err(FastqError::Io)?;
        let mut magic = [0u8; 6];
        let mut filled = 0usize;
        while filled < magic.len() {
            let n = file
                .read(&mut magic[filled..])
                .await
                .map_err(FastqError::Io)?;
            if n == 0 {
                break;
            }
            filled += n;
        }
        drop(file);

        let head = &magic[..filled];
        if head.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]) {
            return Err(FastqError::Unsupported(UnsupportedOperation::Zstd));
        }
        if head.starts_with(b"BZh") {
            return Err(FastqError::Unsupported(UnsupportedOperation::Bzip2));
        }
        if head.starts_with(&[0xfd, b'7', b'z', b'X', b'Z', 0x00]) {
            return Err(FastqError::Unsupported(UnsupportedOperation::Xz));
        }

        // Re-open to rewind: one extra `open` is nothing next to the per-record cost.
        let file = File::open(path).await.map_err(FastqError::Io)?;
        let buffered = BufReader::with_capacity(256 * 1024, file);
        if filled >= 2 && magic[0] == 0x1f && magic[1] == 0x8b {
            let mut decoder = GzipDecoder::new(buffered);
            // gzip files in bioinformatics are routinely concatenations of members: bgzip emits
            // one per 64 KiB block, pigz one per chunk. Without this the reader stops after the
            // first member and silently drops the rest of the file.
            decoder.multiple_members(true);
            let buffered = BufReader::with_capacity(256 * 1024, decoder);
            return Ok(AnyAsyncReader::Gzip(AsyncFastqReader::from_reader(
                buffered,
            )));
        }
        Ok(AnyAsyncReader::Plain(AsyncFastqReader::from_reader(
            buffered,
        )))
    }

    /// Read the next record. Not cancel-safe; see [`AsyncFastqReader`].
    pub async fn next(&mut self) -> Result<Option<FastqRecord<'_>>, FastqError> {
        match self {
            Self::Plain(r) => r.next().await,
            Self::Gzip(r) => r.next().await,
        }
    }

    /// Number of records returned so far.
    pub fn records_read(&self) -> u64 {
        match self {
            Self::Plain(r) => r.records_read(),
            Self::Gzip(r) => r.records_read(),
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

    pub fn with_quality_encoding(self, encoding: QualityEncoding) -> Self {
        match self {
            Self::Plain(r) => Self::Plain(r.with_quality_encoding(encoding)),
            Self::Gzip(r) => Self::Gzip(r.with_quality_encoding(encoding)),
        }
    }

    pub fn with_format(self, format: FastqFormat) -> Self {
        match self {
            Self::Plain(r) => Self::Plain(r.with_format(format)),
            Self::Gzip(r) => Self::Gzip(r.with_format(format)),
        }
    }
}

/// Stream of owned records, from [`AsyncFastqReader::records`].
pub struct RecordStream<R: AsyncBufRead + Unpin + Send + 'static> {
    #[allow(clippy::type_complexity)]
    inner: Pin<Box<dyn Stream<Item = Result<FastqRecordOwned, FastqError>> + Send + 'static>>,
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

/// Read one line, stripping `\n` and an optional preceding `\r`. A final line without a
/// terminator is reported as [`LineStatus::EofPartial`] and treated as a line by the callers.
async fn read_line<R: AsyncBufRead + Unpin>(
    reader: &mut R,
    out: &mut Vec<u8>,
    logical_offset: &mut u64,
) -> Result<LineStatus, FastqError> {
    out.clear();
    let n = reader
        .read_until(b'\n', out)
        .await
        .map_err(FastqError::Io)?;
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
    quality: QualityEncoding,
    seq: &[u8],
    qual: &[u8],
    seq_start: u64,
    qual_start: u64,
) -> Result<(), FastqError> {
    match mode {
        ValidationMode::None => Ok(()),
        ValidationMode::Bases => check_bases(seq, alphabet, seq_start),
        ValidationMode::Qualities => check_qual(qual, quality, qual_start),
        ValidationMode::BasesAndQualities => {
            check_bases(seq, alphabet, seq_start)?;
            check_qual(qual, quality, qual_start)
        }
    }
}

#[inline]
fn check_bases(seq: &[u8], alphabet: Alphabet, base: u64) -> Result<(), FastqError> {
    match validate_bases_with(seq, alphabet) {
        Ok(()) => Ok(()),
        Err(idx) => Err(FastqError::invalid_base(base + idx as u64, seq[idx])),
    }
}

#[inline]
fn check_qual(qual: &[u8], encoding: QualityEncoding, base: u64) -> Result<(), FastqError> {
    match validate_qual_encoding(qual, encoding) {
        Ok(()) => Ok(()),
        Err(idx) => Err(FastqError::invalid_quality(base + idx as u64, qual[idx])),
    }
}
