use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::backend::Backend;
use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::GzipBackend;
use crate::backend::mmap::MmapBackend;
use crate::error::FastqError;
use crate::error::UnsupportedOperation;
use crate::format::FastqFormat;
use crate::multiline::MultiLineFastqParser;
use crate::offset::VirtualOffset;
use crate::parser::{FastqParser, ParsedRecord, RecordScratch, Segment};
use crate::record::FastqRecord;
use crate::simd::bases::validate_bases_with;
use crate::simd::newline::find_lf;
use crate::simd::qual::validate_qual;
use crate::validation::{Alphabet, ValidationMode};

pub struct FastqReader {
    backend: Backend,
    parser: FastqParser,
    parser_multi: MultiLineFastqParser,
    pos: usize,
    validation: ValidationMode,
    alphabet: Alphabet,
    format: FastqFormat,
    scratch: RecordScratch,
    seq_segs: Vec<Segment>,
    qual_segs: Vec<Segment>,
    resync: bool,
}

impl FastqReader {
    /// Backend is chosen from the file extension.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        let backend = if is_bgzf(path) {
            Backend::Bgzf(BgzfBackend::new(path)?)
        } else if is_gz(path) {
            Backend::Gzip(GzipBackend::new(path)?)
        } else {
            let mmap = MmapBackend::open(path)?;
            mmap.advise_sequential();
            Backend::Plain(mmap)
        };
        Ok(Self::with_backend(backend))
    }

    /// Backend is chosen by sniffing the file's magic bytes.
    pub fn from_path_auto<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        let kind = detect_compression(path)?;
        let backend = match kind {
            CompressionKind::Plain => {
                let mmap = MmapBackend::open(path)?;
                mmap.advise_sequential();
                Backend::Plain(mmap)
            }
            CompressionKind::Gzip => Backend::Gzip(GzipBackend::new(path)?),
            CompressionKind::Bgzf => Backend::Bgzf(BgzfBackend::new(path)?),
        };
        Ok(Self::with_backend(backend))
    }

    pub fn from_bgzf_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        Ok(Self::with_backend(Backend::Bgzf(BgzfBackend::new(
            path.as_ref(),
        )?)))
    }

    fn with_backend(backend: Backend) -> Self {
        Self {
            backend,
            parser: FastqParser::new(),
            parser_multi: MultiLineFastqParser::new(),
            pos: 0,
            validation: ValidationMode::None,
            alphabet: Alphabet::default(),
            format: FastqFormat::default(),
            scratch: RecordScratch::new(),
            seq_segs: Vec::new(),
            qual_segs: Vec::new(),
            resync: false,
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

    #[inline]
    pub fn next(&mut self) -> Result<Option<FastqRecord<'_>>, FastqError> {
        Ok(self.next_parsed()?.map(|p| p.record))
    }

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

        if self.resync {
            self.resync = false;
            match &mut self.backend {
                Backend::Plain(mmap) => {
                    let buf = mmap.slice_from(0);
                    resync_plain(buf, &mut self.pos);
                }
                Backend::Bgzf(bgzf) => {
                    resync_streamed_bgzf(bgzf)?;
                }
                Backend::Gzip(_) => {}
            }
        }

        let parsed = match (format, &mut self.backend) {
            (FastqFormat::SingleLine, Backend::Plain(mmap)) => {
                let buf = mmap.slice_from(0);
                self.parser.next_record_in_slice(buf, &mut self.pos)
            }
            (FastqFormat::SingleLine, Backend::Gzip(gzip)) => {
                self.parser.next_record_gzip(gzip, &mut self.scratch)
            }
            (FastqFormat::SingleLine, Backend::Bgzf(bgzf)) => {
                self.parser.next_record_bgzf(bgzf, &mut self.scratch)
            }
            (FastqFormat::MultiLine, Backend::Plain(mmap)) => {
                let buf = mmap.slice_from(0);
                self.parser_multi.next_record_in_slice(
                    buf,
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
        }?;

        let Some(parsed) = parsed else {
            return Ok(None);
        };

        if format == FastqFormat::SingleLine {
            validate_record_singleline(mode, alphabet, &parsed)?;
        } else {
            validate_record_multiline(mode, alphabet, &parsed, &self.seq_segs, &self.qual_segs)?;
        }
        Ok(Some(parsed))
    }

    pub fn tell(&self) -> VirtualOffset {
        match &self.backend {
            Backend::Plain(_) => VirtualOffset(self.pos as u64),
            Backend::Gzip(gz) => VirtualOffset(gz.logical_offset()),
            Backend::Bgzf(bgzf) => bgzf.tell(),
        }
    }

    pub fn seek(&mut self, voff: VirtualOffset) -> Result<(), FastqError> {
        match &mut self.backend {
            Backend::Plain(_) => {
                self.pos = voff.0 as usize;
                self.resync = true;
                Ok(())
            }
            Backend::Bgzf(bgzf) => {
                bgzf.seek(voff)?;
                self.resync = true;
                Ok(())
            }
            Backend::Gzip(_) => Err(FastqError::Unsupported(UnsupportedOperation::Seek)),
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
    parsed: &ParsedRecord<'_>,
) -> Result<(), FastqError> {
    match mode {
        ValidationMode::None => Ok(()),
        ValidationMode::Bases => check_bases(parsed, alphabet),
        ValidationMode::Qualities => check_qual(parsed),
        ValidationMode::BasesAndQualities => {
            check_bases(parsed, alphabet)?;
            check_qual(parsed)
        }
    }
}

#[inline]
fn check_bases(parsed: &ParsedRecord<'_>, alphabet: Alphabet) -> Result<(), FastqError> {
    if let Err(idx) = validate_bases_with(parsed.record.seq(), alphabet) {
        let b = parsed.record.seq()[idx];
        return Err(FastqError::InvalidBase {
            offset: parsed.seq_start + idx as u64,
            byte: b,
        });
    }
    Ok(())
}

#[inline]
fn check_qual(parsed: &ParsedRecord<'_>) -> Result<(), FastqError> {
    if let Err(idx) = validate_qual(parsed.record.qual()) {
        let b = parsed.record.qual()[idx];
        return Err(FastqError::InvalidQuality {
            offset: parsed.qual_start + idx as u64,
            byte: b,
        });
    }
    Ok(())
}

#[inline]
fn validate_record_multiline(
    mode: ValidationMode,
    alphabet: Alphabet,
    parsed: &ParsedRecord<'_>,
    seq_segs: &[Segment],
    qual_segs: &[Segment],
) -> Result<(), FastqError> {
    match mode {
        ValidationMode::None => Ok(()),
        ValidationMode::Bases => check_bases_ml(parsed, alphabet, seq_segs),
        ValidationMode::Qualities => check_qual_ml(parsed, qual_segs),
        ValidationMode::BasesAndQualities => {
            check_bases_ml(parsed, alphabet, seq_segs)?;
            check_qual_ml(parsed, qual_segs)
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
        let offset = map_offset(segs, idx, parsed.seq_start);
        return Err(FastqError::InvalidBase { offset, byte: b });
    }
    Ok(())
}

#[inline]
fn check_qual_ml(parsed: &ParsedRecord<'_>, segs: &[Segment]) -> Result<(), FastqError> {
    if let Err(idx) = validate_qual(parsed.record.qual()) {
        let b = parsed.record.qual()[idx];
        let offset = map_offset(segs, idx, parsed.qual_start);
        return Err(FastqError::InvalidQuality { offset, byte: b });
    }
    Ok(())
}

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

// Quartet anchoring (`@…\n…\n+…\n…\n`) avoids the classic FASTQ resync pitfall
// where a quality line legitimately starts with `@`.
fn resync_plain(buf: &[u8], pos: &mut usize) {
    let mut p = *pos;
    if p > buf.len() {
        p = buf.len();
    }
    loop {
        if p >= buf.len() {
            *pos = p;
            return;
        }
        let at_line_start = p == 0 || buf[p - 1] == b'\n';
        if at_line_start && buf[p] == b'@' && looks_like_record_start(buf, p) {
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
    }
}

fn looks_like_record_start(buf: &[u8], p: usize) -> bool {
    let mut cursor = p;
    let lf1 = match find_lf(buf, cursor) {
        Some(v) => v,
        None => return false,
    };
    cursor = lf1 + 1;
    let lf2 = match find_lf(buf, cursor) {
        Some(v) => v,
        None => return false,
    };
    let seq_len = lf2.saturating_sub(cursor);
    let seq_len = if seq_len > 0 && buf[lf2 - 1] == b'\r' {
        seq_len - 1
    } else {
        seq_len
    };
    cursor = lf2 + 1;
    if cursor >= buf.len() || buf[cursor] != b'+' {
        return false;
    }
    let lf3 = match find_lf(buf, cursor) {
        Some(v) => v,
        None => return false,
    };
    cursor = lf3 + 1;
    let lf4 = match find_lf(buf, cursor) {
        Some(v) => v,
        None => return false,
    };
    let qual_len = lf4.saturating_sub(cursor);
    let qual_len = if qual_len > 0 && buf[lf4 - 1] == b'\r' {
        qual_len - 1
    } else {
        qual_len
    };
    seq_len == qual_len && seq_len > 0
}

// Streaming BGZF cannot afford the multi-line quartet lookahead; this matches htslib semantics.
fn resync_streamed_bgzf(bgzf: &mut BgzfBackend) -> Result<(), FastqError> {
    let mut buf = Vec::new();
    loop {
        match bgzf.peek_byte()? {
            Some(b'@') => return Ok(()),
            Some(_) => {
                buf.clear();
                let status = bgzf.read_line(&mut buf)?;
                if matches!(
                    status,
                    crate::backend::gzip::LineStatus::EofClean
                        | crate::backend::gzip::LineStatus::EofPartial
                ) {
                    return Ok(());
                }
            }
            None => return Ok(()),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CompressionKind {
    Plain,
    Gzip,
    Bgzf,
}

fn detect_compression(path: &Path) -> Result<CompressionKind, FastqError> {
    let mut file = File::open(path)?;
    let mut buf = [0u8; 128];
    let n = file.read(&mut buf)?;
    if n < 2 {
        return Ok(CompressionKind::Plain);
    }
    if buf[0] != 0x1f || buf[1] != 0x8b {
        return Ok(CompressionKind::Plain);
    }
    if n < 12 {
        return Ok(CompressionKind::Gzip);
    }
    let flg = buf[3];
    if (flg & 0x04) == 0 {
        return Ok(CompressionKind::Gzip);
    }
    let xlen = u16::from_le_bytes([buf[10], buf[11]]) as usize;
    let extra_end = 12usize.saturating_add(xlen);
    if extra_end > n {
        return Ok(CompressionKind::Gzip);
    }
    let mut i = 12;
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
            return Ok(CompressionKind::Bgzf);
        }
        i = sub_end;
    }
    Ok(CompressionKind::Gzip)
}
