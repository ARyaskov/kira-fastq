use std::fs::File;
use std::io::{BufRead, BufReader, Cursor, Read};
use std::path::Path;

use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::GzipBackend;
use crate::backend::memory::ContiguousBackend;
use crate::backend::mmap::MmapBackend;
use crate::backend::parallel::ParallelBgzfReader;
use crate::backend::stream::StreamBackend;
use crate::backend::{Backend, LineStatus};
use crate::error::{FastqError, UnsupportedOperation};
use crate::format::FastqFormat;
use crate::multiline::MultiLineFastqParser;
use crate::offset::VirtualOffset;
use crate::parser::{FastqParser, ParsedRecord, RecordScratch, Segment};
use crate::record::FastqRecord;
use crate::simd::bases::validate_bases_with;
use crate::simd::newline::find_lf;
use crate::simd::qual::validate_qual_encoding;
use crate::validation::{Alphabet, QualityEncoding, ValidationMode};

#[cfg(feature = "noodles-bgzf")]
use crate::backend::noodles_bgzf::NoodlesBgzfBackend;

/// Buffer used by the constructors that read through `BufRead` instead of mapping.
const STREAM_BUF: usize = 256 * 1024;

/// How many candidate line starts `seek` will examine before giving up on resynchronising.
const RESYNC_LINE_BUDGET: usize = 4096;

/// Reads FASTQ records from a file, a buffer, or any stream.
///
/// Construct with [`FastqReader::from_path`], which picks the backend from the file's magic
/// bytes, or with one of the explicit constructors when the source is not a file.
pub struct FastqReader {
    backend: Backend,
    parser: FastqParser,
    parser_multi: MultiLineFastqParser,
    pos: usize,
    validation: ValidationMode,
    alphabet: Alphabet,
    quality: QualityEncoding,
    format: FastqFormat,
    scratch: RecordScratch,
    seq_segs: Vec<Segment>,
    qual_segs: Vec<Segment>,
    records: u64,
}

impl FastqReader {
    /// Open a path, choosing the backend from the file's magic bytes: plain, gzip, or BGZF.
    ///
    /// The file name is not consulted, because it routinely lies: `bgzip` writes BGZF to a plain
    /// `.gz` name, and a `.fastq` produced by a pipeline can be gzip. Compressed formats this
    /// crate does not decode (zstd without the `zstd` feature, bzip2, xz) are reported as
    /// [`FastqError::Unsupported`] rather than as a parse error.
    ///
    /// Plain files are memory-mapped and parsed in place. See
    /// [`FastqReader::from_path_buffered`] for the streaming alternative.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        let backend = match detect_compression(path)? {
            CompressionKind::Plain => {
                let mmap = MmapBackend::open(path)?;
                mmap.advise_sequential();
                Backend::Plain(ContiguousBackend::Mmap(mmap))
            }
            CompressionKind::Gzip => Backend::Gzip(GzipBackend::new(path)?),
            CompressionKind::Bgzf => Backend::Bgzf(BgzfBackend::new(path)?),
            CompressionKind::Zstd => return open_zstd(path).map(Self::with_backend),
            CompressionKind::Bzip2 => {
                return Err(FastqError::Unsupported(UnsupportedOperation::Bzip2));
            }
            CompressionKind::Xz => return Err(FastqError::Unsupported(UnsupportedOperation::Xz)),
        };
        Ok(Self::with_backend(backend))
    }

    /// Same as [`FastqReader::from_path`], kept for compatibility with 0.3.
    #[inline]
    pub fn from_path_auto<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        Self::from_path(path)
    }

    /// Open a path without memory-mapping it, reading through a buffered stream instead.
    ///
    /// Use this for files another process is still writing (a mapping would take a snapshot and
    /// can fault if the file shrinks), for network file systems, and on platforms where mapping
    /// measures slower than `read`. Compression is detected exactly as in
    /// [`FastqReader::from_path`]. Records are copied into the reader's scratch rather than
    /// borrowed from a mapping; `seek` is not available on this path.
    pub fn from_path_buffered<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        let kind = detect_compression(path)?;
        let file = File::open(path)?;
        let buffered = BufReader::with_capacity(STREAM_BUF, file);
        match kind {
            CompressionKind::Plain => Ok(Self::from_reader(buffered)),
            CompressionKind::Gzip | CompressionKind::Bgzf => Ok(Self::from_reader(
                BufReader::with_capacity(STREAM_BUF, flate2::read::MultiGzDecoder::new(buffered)),
            )),
            CompressionKind::Zstd => open_zstd(path).map(Self::with_backend),
            CompressionKind::Bzip2 => Err(FastqError::Unsupported(UnsupportedOperation::Bzip2)),
            CompressionKind::Xz => Err(FastqError::Unsupported(UnsupportedOperation::Xz)),
        }
    }

    /// Force the BGZF backend, which supports `tell`/`seek` by virtual offset.
    pub fn from_bgzf_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        Ok(Self::with_backend(Backend::Bgzf(BgzfBackend::new(
            path.as_ref(),
        )?)))
    }

    /// Open a BGZF file and inflate its blocks on `threads` worker threads.
    ///
    /// BGZF blocks are independent, so this scales inflate across cores the way `samtools -@`
    /// does. Decoded blocks are reassembled in file order, so records come out unchanged.
    /// `tell` and `seek` are not available on this path; pass `threads = 0` to let the crate
    /// pick a thread count from the machine's parallelism.
    pub fn from_bgzf_path_parallel<P: AsRef<Path>>(
        path: P,
        threads: usize,
    ) -> Result<Self, FastqError> {
        let threads = if threads == 0 {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        } else {
            threads
        };
        let reader = ParallelBgzfReader::open(path.as_ref(), threads, true)?;
        Ok(Self::from_reader(reader))
    }

    /// Parse an in-memory buffer. Records borrow out of it, so this path copies nothing.
    pub fn from_vec(data: Vec<u8>) -> Self {
        Self::with_backend(Backend::Plain(ContiguousBackend::Owned(data)))
    }

    /// Parse a slice by copying it into the reader. For a `&'static [u8]`, or to avoid the copy,
    /// use [`FastqReader::from_reader`] instead.
    pub fn from_bytes(data: &[u8]) -> Self {
        Self::from_vec(data.to_vec())
    }

    /// Read from any [`BufRead`] source: stdin, a socket, a decoder from another crate.
    ///
    /// There is no mapping and no random access on this path; every record's lines are copied
    /// into the reader's scratch buffer. Validation, multi-line parsing and paired reading all
    /// work the same.
    pub fn from_reader<R: BufRead + Send + 'static>(reader: R) -> Self {
        Self::with_backend(Backend::Stream(StreamBackend::new(Box::new(reader))))
    }

    /// Like [`FastqReader::from_reader`], but sniffs the stream and decompresses it when it
    /// turns out to be gzip or BGZF. Use it for stdin, where compressed input is the norm.
    pub fn from_reader_auto<R: BufRead + Send + 'static>(
        mut reader: R,
    ) -> Result<Self, FastqError> {
        let mut magic = [0u8; 6];
        let mut filled = 0usize;
        while filled < magic.len() {
            let n = reader.read(&mut magic[filled..]).map_err(FastqError::Io)?;
            if n == 0 {
                break;
            }
            filled += n;
        }
        let head = Cursor::new(magic[..filled].to_vec());
        let rest = head.chain(reader);
        match classify_magic(&magic[..filled]) {
            CompressionKind::Plain => Ok(Self::from_unbuffered(rest)),
            CompressionKind::Gzip | CompressionKind::Bgzf => Ok(Self::from_reader(
                BufReader::with_capacity(STREAM_BUF, flate2::read::MultiGzDecoder::new(rest)),
            )),
            CompressionKind::Zstd => {
                #[cfg(feature = "zstd")]
                {
                    let decoder = zstd::stream::read::Decoder::new(rest).map_err(FastqError::Io)?;
                    Ok(Self::from_unbuffered(decoder))
                }
                #[cfg(not(feature = "zstd"))]
                Err(FastqError::Unsupported(UnsupportedOperation::Zstd))
            }
            CompressionKind::Bzip2 => Err(FastqError::Unsupported(UnsupportedOperation::Bzip2)),
            CompressionKind::Xz => Err(FastqError::Unsupported(UnsupportedOperation::Xz)),
        }
    }

    /// Wrap a [`Read`] in a [`BufReader`] and delegate to [`FastqReader::from_reader`].
    pub fn from_unbuffered<R: Read + Send + 'static>(reader: R) -> Self {
        Self::from_reader(BufReader::with_capacity(STREAM_BUF, reader))
    }

    /// Open via the optional `noodles-bgzf` adapter, for pipelines that share virtual offsets
    /// with the rest of the noodles ecosystem. Supports `tell` and `seek`.
    #[cfg(feature = "noodles-bgzf")]
    pub fn from_noodles_bgzf_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        Ok(Self::with_backend(Backend::NoodlesBgzf(
            NoodlesBgzfBackend::open(path.as_ref())?,
        )))
    }

    fn with_backend(backend: Backend) -> Self {
        Self {
            backend,
            parser: FastqParser::new(),
            parser_multi: MultiLineFastqParser::new(),
            pos: 0,
            validation: ValidationMode::None,
            alphabet: Alphabet::default(),
            quality: QualityEncoding::default(),
            format: FastqFormat::default(),
            scratch: RecordScratch::new(),
            seq_segs: Vec::new(),
            qual_segs: Vec::new(),
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

    /// Set the quality encoding used when validating quality bytes. Defaults to Phred+33; pass
    /// [`QualityEncoding::PHRED64`] for Illumina 1.3 to 1.7 data.
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

    /// Turn off the check for the BGZF end-of-file marker.
    ///
    /// The marker is how a truncated BGZF file is detected, so leave it on unless you knowingly
    /// read a partial file, e.g. one still being written. No effect on other backends.
    #[inline]
    pub fn with_bgzf_eof_check(mut self, enabled: bool) -> Self {
        if let Backend::Bgzf(bgzf) = &mut self.backend {
            bgzf.set_eof_check(enabled);
        }
        self
    }

    /// Number of records returned so far.
    #[inline]
    pub fn records_read(&self) -> u64 {
        self.records
    }

    /// Next record, borrowed from the reader. The borrow ends at the next call.
    #[inline]
    pub fn next(&mut self) -> Result<Option<FastqRecord<'_>>, FastqError> {
        Ok(self.next_parsed()?.map(|p| p.record))
    }

    /// Run `f` over every remaining record, stopping at the first error from either side.
    pub fn try_for_each<E, F>(&mut self, mut f: F) -> Result<(), TryForEachError<E>>
    where
        F: FnMut(FastqRecord<'_>) -> Result<(), E>,
    {
        loop {
            match self.next() {
                Ok(Some(rec)) => {
                    if let Err(e) = f(rec) {
                        return Err(TryForEachError::User(e));
                    }
                }
                Ok(None) => return Ok(()),
                Err(e) => return Err(TryForEachError::Fastq(e)),
            }
        }
    }

    pub(crate) fn next_parsed<'a>(&'a mut self) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        let mode = self.validation;
        let format = self.format;
        let alphabet = self.alphabet;
        let quality = self.quality;
        let record_index = self.records + 1;

        let parsed = match (format, &mut self.backend) {
            (FastqFormat::SingleLine, Backend::Plain(src)) => self
                .parser
                .next_record_in_slice(src.as_slice(), &mut self.pos),
            (FastqFormat::SingleLine, Backend::Gzip(gzip)) => {
                self.parser.next_record_stream(gzip, &mut self.scratch)
            }
            (FastqFormat::SingleLine, Backend::Bgzf(bgzf)) => {
                self.parser.next_record_stream(bgzf, &mut self.scratch)
            }
            (FastqFormat::SingleLine, Backend::Stream(s)) => {
                self.parser.next_record_stream(s, &mut self.scratch)
            }
            #[cfg(feature = "noodles-bgzf")]
            (FastqFormat::SingleLine, Backend::NoodlesBgzf(b)) => {
                self.parser.next_record_stream(b, &mut self.scratch)
            }
            (FastqFormat::MultiLine, Backend::Plain(src)) => {
                self.parser_multi.next_record_in_slice(
                    src.as_slice(),
                    &mut self.pos,
                    &mut self.scratch,
                    &mut self.seq_segs,
                    &mut self.qual_segs,
                )
            }
            (FastqFormat::MultiLine, Backend::Gzip(gzip)) => self.parser_multi.next_record_stream(
                gzip,
                &mut self.scratch,
                &mut self.seq_segs,
                &mut self.qual_segs,
            ),
            (FastqFormat::MultiLine, Backend::Bgzf(bgzf)) => self.parser_multi.next_record_stream(
                bgzf,
                &mut self.scratch,
                &mut self.seq_segs,
                &mut self.qual_segs,
            ),
            (FastqFormat::MultiLine, Backend::Stream(s)) => self.parser_multi.next_record_stream(
                s,
                &mut self.scratch,
                &mut self.seq_segs,
                &mut self.qual_segs,
            ),
            #[cfg(feature = "noodles-bgzf")]
            (FastqFormat::MultiLine, Backend::NoodlesBgzf(b)) => {
                self.parser_multi.next_record_stream(
                    b,
                    &mut self.scratch,
                    &mut self.seq_segs,
                    &mut self.qual_segs,
                )
            }
        }
        .map_err(|e| e.with_record(record_index))?;

        let Some(parsed) = parsed else {
            return Ok(None);
        };

        if format == FastqFormat::SingleLine {
            validate_record_singleline(mode, alphabet, quality, &parsed)
        } else {
            validate_record_multiline(
                mode,
                alphabet,
                quality,
                &parsed,
                &self.seq_segs,
                &self.qual_segs,
            )
        }
        .map_err(|e| e.with_record(record_index))?;

        self.records = record_index;
        Ok(Some(parsed))
    }

    /// Current position: a byte offset for plain, in-memory and stream sources, a BGZF virtual
    /// offset for the BGZF backends. Valid as an argument to [`FastqReader::seek`] where that is
    /// supported.
    pub fn tell(&self) -> VirtualOffset {
        match &self.backend {
            Backend::Plain(_) => VirtualOffset(self.pos as u64),
            Backend::Gzip(gz) => VirtualOffset(gz.logical_offset()),
            Backend::Bgzf(bgzf) => bgzf.tell(),
            Backend::Stream(s) => VirtualOffset(s.logical_offset()),
            #[cfg(feature = "noodles-bgzf")]
            Backend::NoodlesBgzf(b) => b.tell(),
        }
    }

    /// Move to `voff` and resynchronise to the next record boundary.
    ///
    /// An offset taken from [`FastqReader::tell`] lands exactly where it was taken. An arbitrary
    /// offset scans forward for a record start, checking the whole record shape rather than just
    /// a leading `@`: `@` is a legal quality byte, so a bare `@`-at-line-start test lands inside
    /// records. After this call [`FastqReader::tell`] reports the resynchronised position.
    ///
    /// Supported on plain, in-memory and BGZF sources. gzip and arbitrary streams return
    /// [`FastqError::Unsupported`].
    pub fn seek(&mut self, voff: VirtualOffset) -> Result<(), FastqError> {
        let format = self.format;
        match &mut self.backend {
            Backend::Plain(src) => {
                self.pos = (voff.0 as usize).min(src.len());
                let slice = src.as_slice();
                resync_slice(slice, &mut self.pos, format);
                self.records = 0;
                Ok(())
            }
            Backend::Bgzf(bgzf) => {
                bgzf.seek(voff)?;
                resync_seekable(bgzf, format)?;
                self.records = 0;
                Ok(())
            }
            #[cfg(feature = "noodles-bgzf")]
            Backend::NoodlesBgzf(b) => {
                b.seek(voff)?;
                resync_seekable(b, format)?;
                self.records = 0;
                Ok(())
            }
            Backend::Gzip(_) | Backend::Stream(_) => {
                Err(FastqError::Unsupported(UnsupportedOperation::Seek))
            }
        }
    }
}

#[derive(Debug)]
pub enum TryForEachError<E> {
    Fastq(FastqError),
    User(E),
}

impl<E: std::fmt::Display> std::fmt::Display for TryForEachError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fastq(e) => write!(f, "FASTQ error: {e}"),
            Self::User(e) => write!(f, "user callback error: {e}"),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for TryForEachError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fastq(e) => Some(e),
            Self::User(e) => Some(e),
        }
    }
}

#[inline]
fn validate_record_singleline(
    mode: ValidationMode,
    alphabet: Alphabet,
    quality: QualityEncoding,
    parsed: &ParsedRecord<'_>,
) -> Result<(), FastqError> {
    match mode {
        ValidationMode::None => Ok(()),
        ValidationMode::Bases => check_bases(parsed, alphabet),
        ValidationMode::Qualities => check_qual(parsed, quality),
        ValidationMode::BasesAndQualities => {
            check_bases(parsed, alphabet)?;
            check_qual(parsed, quality)
        }
    }
}

#[inline]
fn check_bases(parsed: &ParsedRecord<'_>, alphabet: Alphabet) -> Result<(), FastqError> {
    if let Err(idx) = validate_bases_with(parsed.record.seq(), alphabet) {
        let b = parsed.record.seq()[idx];
        return Err(FastqError::invalid_base(parsed.seq_start + idx as u64, b));
    }
    Ok(())
}

#[inline]
fn check_qual(parsed: &ParsedRecord<'_>, quality: QualityEncoding) -> Result<(), FastqError> {
    if let Err(idx) = validate_qual_encoding(parsed.record.qual(), quality) {
        let b = parsed.record.qual()[idx];
        return Err(FastqError::invalid_quality(
            parsed.qual_start + idx as u64,
            b,
        ));
    }
    Ok(())
}

#[inline]
fn validate_record_multiline(
    mode: ValidationMode,
    alphabet: Alphabet,
    quality: QualityEncoding,
    parsed: &ParsedRecord<'_>,
    seq_segs: &[Segment],
    qual_segs: &[Segment],
) -> Result<(), FastqError> {
    match mode {
        ValidationMode::None => Ok(()),
        ValidationMode::Bases => check_bases_ml(parsed, alphabet, seq_segs),
        ValidationMode::Qualities => check_qual_ml(parsed, quality, qual_segs),
        ValidationMode::BasesAndQualities => {
            check_bases_ml(parsed, alphabet, seq_segs)?;
            check_qual_ml(parsed, quality, qual_segs)
        }
    }
}

#[inline]
fn check_bases_ml(
    parsed: &ParsedRecord<'_>,
    alphabet: Alphabet,
    segs: &[Segment],
) -> Result<(), FastqError> {
    if let Err(idx) = validate_bases_with(parsed.record.seq(), alphabet) {
        let b = parsed.record.seq()[idx];
        return Err(FastqError::invalid_base(
            map_offset(segs, idx, parsed.seq_start),
            b,
        ));
    }
    Ok(())
}

#[inline]
fn check_qual_ml(
    parsed: &ParsedRecord<'_>,
    quality: QualityEncoding,
    segs: &[Segment],
) -> Result<(), FastqError> {
    if let Err(idx) = validate_qual_encoding(parsed.record.qual(), quality) {
        let b = parsed.record.qual()[idx];
        return Err(FastqError::invalid_quality(
            map_offset(segs, idx, parsed.qual_start),
            b,
        ));
    }
    Ok(())
}

/// Map an index inside a joined multi-line field back to a byte offset in the file.
#[inline]
fn map_offset(segs: &[Segment], mut idx: usize, fallback: u64) -> u64 {
    for seg in segs {
        if idx < seg.len {
            return seg.offset + idx as u64;
        }
        idx -= seg.len;
    }
    fallback + idx as u64
}

/// Does a record start at these lines?
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Verdict {
    Yes,
    No,
    NeedMore,
}

/// Check the shape of a whole record, which is the only reliable way to find a record boundary:
/// `@` is a valid quality byte, so quality lines routinely start with it.
fn verify_record_start(lines: &[&[u8]], format: FastqFormat) -> Verdict {
    let Some(header) = lines.first() else {
        return Verdict::NeedMore;
    };
    if header.first() != Some(&b'@') {
        return Verdict::No;
    }
    match format {
        FastqFormat::SingleLine => {
            if lines.len() < 4 {
                return Verdict::NeedMore;
            }
            if lines[2].first() != Some(&b'+') {
                return Verdict::No;
            }
            if lines[1].len() != lines[3].len() {
                return Verdict::No;
            }
            Verdict::Yes
        }
        FastqFormat::MultiLine => {
            let mut seq_len = 0usize;
            let mut i = 1usize;
            loop {
                let Some(line) = lines.get(i) else {
                    return Verdict::NeedMore;
                };
                if line.first() == Some(&b'+') {
                    break;
                }
                if line.is_empty() {
                    return Verdict::No;
                }
                seq_len += line.len();
                i += 1;
            }
            if i == 1 {
                return Verdict::No;
            }
            let mut qual_len = 0usize;
            let mut j = i + 1;
            while qual_len < seq_len {
                let Some(line) = lines.get(j) else {
                    return Verdict::NeedMore;
                };
                if line.is_empty() {
                    return Verdict::No;
                }
                qual_len += line.len();
                j += 1;
            }
            if qual_len == seq_len {
                Verdict::Yes
            } else {
                Verdict::No
            }
        }
    }
}

/// Scan forward in a contiguous buffer for the next record start.
fn resync_slice(buf: &[u8], pos: &mut usize, format: FastqFormat) {
    let mut p = (*pos).min(buf.len());
    let mut budget = RESYNC_LINE_BUDGET;
    loop {
        if p >= buf.len() {
            *pos = buf.len();
            return;
        }
        let at_line_start = p == 0 || buf[p - 1] == b'\n';
        if at_line_start && buf[p] == b'@' && slice_record_starts_at(buf, p, format) {
            *pos = p;
            return;
        }
        match find_lf(buf, p) {
            Some(lf) => p = lf + 1,
            None => {
                *pos = buf.len();
                return;
            }
        }
        budget -= 1;
        if budget == 0 {
            *pos = p;
            return;
        }
    }
}

fn slice_record_starts_at(buf: &[u8], start: usize, format: FastqFormat) -> bool {
    let mut lines: Vec<&[u8]> = Vec::new();
    let mut p = start;
    loop {
        match verify_record_start(&lines, format) {
            Verdict::Yes => return true,
            Verdict::No => return false,
            Verdict::NeedMore => {}
        }
        if p >= buf.len() || lines.len() >= RESYNC_LINE_BUDGET {
            return false;
        }
        let line_end = find_lf(buf, p).unwrap_or(buf.len());
        let mut end = line_end;
        if end > p && buf[end - 1] == b'\r' {
            end -= 1;
        }
        lines.push(&buf[p..end]);
        p = line_end + 1;
    }
}

/// A source that can be repositioned by virtual offset, which is what lets `seek` resynchronise
/// on compressed input.
pub(crate) trait SeekableLines {
    fn tell_offset(&self) -> VirtualOffset;
    fn seek_offset(&mut self, voff: VirtualOffset) -> Result<(), FastqError>;
    fn next_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError>;
}

impl SeekableLines for BgzfBackend {
    #[inline]
    fn tell_offset(&self) -> VirtualOffset {
        self.tell()
    }
    #[inline]
    fn seek_offset(&mut self, voff: VirtualOffset) -> Result<(), FastqError> {
        self.seek(voff)
    }
    #[inline]
    fn next_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        self.read_line(out)
    }
}

#[cfg(feature = "noodles-bgzf")]
impl SeekableLines for NoodlesBgzfBackend {
    #[inline]
    fn tell_offset(&self) -> VirtualOffset {
        self.tell()
    }
    #[inline]
    fn seek_offset(&mut self, voff: VirtualOffset) -> Result<(), FastqError> {
        self.seek(voff)
    }
    #[inline]
    fn next_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        self.read_line(out)
    }
}

/// Scan forward from the current position for the next record start, leaving the source
/// positioned there. Candidates are re-examined by seeking back, so nothing is consumed.
fn resync_seekable<S: SeekableLines>(src: &mut S, format: FastqFormat) -> Result<(), FastqError> {
    let mut candidate = src.tell_offset();
    // One budget for the whole operation rather than one per candidate: a multi-line record can
    // legitimately run to hundreds of lines, but a file of noise must not be scanned forever.
    let mut budget = RESYNC_LINE_BUDGET;
    let mut owned: Vec<Vec<u8>> = Vec::new();
    let mut line = Vec::new();

    loop {
        src.seek_offset(candidate)?;
        owned.clear();
        let mut verdict = Verdict::NeedMore;
        while verdict == Verdict::NeedMore && budget > 0 {
            match src.next_line(&mut line)? {
                LineStatus::EofClean => break,
                LineStatus::Line | LineStatus::EofPartial => {
                    owned.push(std::mem::take(&mut line));
                    budget -= 1;
                }
            }
            let borrowed: Vec<&[u8]> = owned.iter().map(|l| l.as_slice()).collect();
            verdict = verify_record_start(&borrowed, format);
        }

        // Leave the source where the caller expects it, whatever the verdict.
        src.seek_offset(candidate)?;
        if verdict == Verdict::Yes || budget == 0 {
            return Ok(());
        }
        // Not a record start: try the next line as the candidate.
        match src.next_line(&mut line)? {
            LineStatus::EofClean => return Ok(()),
            _ => candidate = src.tell_offset(),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CompressionKind {
    Plain,
    Gzip,
    Bgzf,
    Zstd,
    Bzip2,
    Xz,
}

#[cfg(feature = "zstd")]
fn open_zstd(path: &Path) -> Result<Backend, FastqError> {
    let file = File::open(path)?;
    let decoder = zstd::stream::read::Decoder::new(BufReader::with_capacity(STREAM_BUF, file))
        .map_err(FastqError::Io)?;
    Ok(Backend::Stream(StreamBackend::new(Box::new(
        BufReader::with_capacity(STREAM_BUF, decoder),
    ))))
}

#[cfg(not(feature = "zstd"))]
fn open_zstd(_path: &Path) -> Result<Backend, FastqError> {
    Err(FastqError::Unsupported(UnsupportedOperation::Zstd))
}

/// Classify a file by its first bytes.
fn detect_compression(path: &Path) -> Result<CompressionKind, FastqError> {
    let mut file = File::open(path)?;
    let mut buf = [0u8; 128];
    let mut filled = 0usize;
    while filled < buf.len() {
        let n = file.read(&mut buf[filled..])?;
        if n == 0 {
            break;
        }
        filled += n;
    }
    Ok(classify_magic(&buf[..filled]))
}

fn classify_magic(buf: &[u8]) -> CompressionKind {
    if buf.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]) {
        return CompressionKind::Zstd;
    }
    if buf.starts_with(b"BZh") {
        return CompressionKind::Bzip2;
    }
    if buf.starts_with(&[0xfd, b'7', b'z', b'X', b'Z', 0x00]) {
        return CompressionKind::Xz;
    }
    if buf.len() < 2 || buf[0] != 0x1f || buf[1] != 0x8b {
        return CompressionKind::Plain;
    }
    if buf.len() < 12 {
        return CompressionKind::Gzip;
    }
    // BGZF is gzip carrying a `BC` extra subfield in the first member's header.
    if (buf[3] & 0x04) == 0 {
        return CompressionKind::Gzip;
    }
    let xlen = u16::from_le_bytes([buf[10], buf[11]]) as usize;
    let extra_end = 12usize.saturating_add(xlen).min(buf.len());
    let mut i = 12usize;
    while i + 4 <= extra_end {
        let si1 = buf[i];
        let si2 = buf[i + 1];
        let slen = u16::from_le_bytes([buf[i + 2], buf[i + 3]]) as usize;
        i += 4;
        let sub_end = i.saturating_add(slen);
        if sub_end > extra_end {
            break;
        }
        if si1 == b'B' && si2 == b'C' && slen == 2 {
            return CompressionKind::Bgzf;
        }
        i = sub_end;
    }
    CompressionKind::Gzip
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_magic_bytes() {
        assert_eq!(classify_magic(b"@r1\nAC\n"), CompressionKind::Plain);
        assert_eq!(
            classify_magic(&[0x1f, 0x8b, 0x08, 0x00]),
            CompressionKind::Gzip
        );
        assert_eq!(
            classify_magic(&[0x28, 0xb5, 0x2f, 0xfd, 0, 0]),
            CompressionKind::Zstd
        );
        assert_eq!(classify_magic(b"BZh9"), CompressionKind::Bzip2);
        assert_eq!(
            classify_magic(&[0xfd, b'7', b'z', b'X', b'Z', 0x00]),
            CompressionKind::Xz
        );
        let bgzf = [
            0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0x1b,
            0x00,
        ];
        assert_eq!(classify_magic(&bgzf), CompressionKind::Bgzf);
    }

    #[test]
    fn verifies_single_line_record_shape() {
        let lines: Vec<&[u8]> = vec![b"@r1", b"ACGT", b"+", b"!!!!"];
        assert_eq!(
            verify_record_start(&lines, FastqFormat::SingleLine),
            Verdict::Yes
        );
        let quality_line_starting_with_at: Vec<&[u8]> = vec![b"@AAA", b"@r2", b"TTTT", b"+"];
        assert_eq!(
            verify_record_start(&quality_line_starting_with_at, FastqFormat::SingleLine),
            Verdict::No
        );
    }

    #[test]
    fn verifies_multi_line_record_shape() {
        let lines: Vec<&[u8]> = vec![b"@r1", b"ACGT", b"AC", b"+", b"!!!!", b"!!"];
        assert_eq!(
            verify_record_start(&lines, FastqFormat::MultiLine),
            Verdict::Yes
        );
        let short: Vec<&[u8]> = vec![b"@r1", b"ACGT", b"+"];
        assert_eq!(
            verify_record_start(&short, FastqFormat::MultiLine),
            Verdict::NeedMore
        );
    }
}
