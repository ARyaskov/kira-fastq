use std::fs::File;
use std::io::Read;
use std::marker::PhantomData;
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
use crate::parser::Segment;
use crate::parser::{FastqParser, ParsedRecord};
use crate::record::FastqRecord;
use crate::simd::bases::validate_bases;
use crate::simd::newline::find_lf;
use crate::simd::qual::validate_qual;
use crate::validation::ValidationMode;
use memchr::memchr_iter;

pub struct FastqReader {
    backend: Backend,
    parser: FastqParser,
    parser_multi: MultiLineFastqParser,
    pos: usize,
    validation: ValidationMode,
    format: FastqFormat,
    ml_seq_buf: Vec<u8>,
    ml_qual_buf: Vec<u8>,
    ml_seq_segs: Vec<Segment>,
    ml_qual_segs: Vec<Segment>,
    resync: bool,
}

impl FastqReader {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        let backend = if is_bgzf(path) {
            Backend::Bgzf(BgzfBackend::new(path)?)
        } else if is_gz(path) {
            Backend::Gzip(GzipBackend::new(path)?)
        } else {
            Backend::Plain(MmapBackend::open(path)?)
        };
        Ok(Self {
            backend,
            parser: FastqParser::new(),
            parser_multi: MultiLineFastqParser::new(),
            pos: 0,
            validation: ValidationMode::None,
            format: FastqFormat::SingleLine,
            ml_seq_buf: Vec::new(),
            ml_qual_buf: Vec::new(),
            ml_seq_segs: Vec::new(),
            ml_qual_segs: Vec::new(),
            resync: false,
        })
    }

    pub fn from_path_auto<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        let kind = detect_compression(path).map_err(FastqError::Io)?;
        let backend = match kind {
            CompressionKind::Plain => Backend::Plain(MmapBackend::open(path)?),
            CompressionKind::Gzip => Backend::Gzip(GzipBackend::new(path)?),
            CompressionKind::Bgzf => Backend::Bgzf(BgzfBackend::new(path)?),
        };
        Ok(Self {
            backend,
            parser: FastqParser::new(),
            parser_multi: MultiLineFastqParser::new(),
            pos: 0,
            validation: ValidationMode::None,
            format: FastqFormat::SingleLine,
            ml_seq_buf: Vec::new(),
            ml_qual_buf: Vec::new(),
            ml_seq_segs: Vec::new(),
            ml_qual_segs: Vec::new(),
            resync: false,
        })
    }

    #[inline]
    pub fn next(&mut self) -> Result<Option<FastqRecord<'_>>, FastqError> {
        let parsed = match self.next_parsed()? {
            Some(p) => p,
            None => return Ok(None),
        };
        Ok(Some(parsed.record))
    }

    #[inline]
    pub fn records(&mut self) -> RecordsIter<'_> {
        RecordsIter {
            reader: self as *mut FastqReader,
            _marker: PhantomData,
        }
    }

    #[inline]
    pub fn with_validation(mut self, mode: ValidationMode) -> Self {
        self.validation = mode;
        self
    }

    #[inline]
    pub fn with_format(mut self, format: FastqFormat) -> Self {
        self.format = format;
        self
    }
    pub(crate) fn next_parsed<'a>(&'a mut self) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        let mode = self.validation;
        let format = self.format;
        if self.resync {
            self.resync = false;
            match &mut self.backend {
                Backend::Plain(mmap) => {
                    let buf = mmap.slice_from(0);
                    resync_plain(buf, &mut self.pos);
                }
                Backend::Bgzf(bgzf) => {
                    resync_bgzf(bgzf)?;
                }
                Backend::Gzip(_) => {}
            }
        }
        let parsed = match (format, &mut self.backend) {
            (FastqFormat::SingleLine, Backend::Plain(mmap)) => {
                let buf = mmap.slice_from(0);
                self.parser.next_record(buf, &mut self.pos)
            }
            (FastqFormat::SingleLine, Backend::Gzip(gzip)) => self.parser.next_record_gzip(gzip),
            (FastqFormat::SingleLine, Backend::Bgzf(bgzf)) => self.parser.next_record_bgzf(bgzf),
            (FastqFormat::MultiLine, Backend::Plain(mmap)) => {
                let buf = mmap.slice_from(0);
                self.parser_multi.next_record(
                    buf,
                    &mut self.pos,
                    &mut self.ml_seq_buf,
                    &mut self.ml_qual_buf,
                    &mut self.ml_seq_segs,
                    &mut self.ml_qual_segs,
                )
            }
            (FastqFormat::MultiLine, Backend::Gzip(gzip)) => self.parser_multi.next_record_gzip(
                gzip,
                &mut self.ml_seq_buf,
                &mut self.ml_qual_buf,
                &mut self.ml_seq_segs,
                &mut self.ml_qual_segs,
            ),
            (FastqFormat::MultiLine, Backend::Bgzf(bgzf)) => self.parser_multi.next_record_bgzf(
                bgzf,
                &mut self.ml_seq_buf,
                &mut self.ml_qual_buf,
                &mut self.ml_seq_segs,
                &mut self.ml_qual_segs,
            ),
        }?;

        let parsed = match parsed {
            Some(p) => p,
            None => return Ok(None),
        };

        if format == FastqFormat::SingleLine {
            validate_record(mode, &parsed)?;
        } else {
            validate_record_multiline(mode, &parsed, &self.ml_seq_segs, &self.ml_qual_segs)?;
        }
        Ok(Some(parsed))
    }
}

#[inline]
fn validate_record(mode: ValidationMode, parsed: &ParsedRecord<'_>) -> Result<(), FastqError> {
    match mode {
        ValidationMode::None => Ok(()),
        ValidationMode::Bases => {
            if let Err(idx) = validate_bases(parsed.record.seq()) {
                let b = parsed.record.seq()[idx];
                return Err(FastqError::InvalidBase {
                    offset: parsed.seq_start + idx as u64,
                    byte: b,
                });
            }
            Ok(())
        }
        ValidationMode::Qualities => {
            if let Err(idx) = validate_qual(parsed.record.qual()) {
                let b = parsed.record.qual()[idx];
                return Err(FastqError::InvalidQuality {
                    offset: parsed.qual_start + idx as u64,
                    byte: b,
                });
            }
            Ok(())
        }
        ValidationMode::BasesAndQualities => {
            if let Err(idx) = validate_bases(parsed.record.seq()) {
                let b = parsed.record.seq()[idx];
                return Err(FastqError::InvalidBase {
                    offset: parsed.seq_start + idx as u64,
                    byte: b,
                });
            }
            if let Err(idx) = validate_qual(parsed.record.qual()) {
                let b = parsed.record.qual()[idx];
                return Err(FastqError::InvalidQuality {
                    offset: parsed.qual_start + idx as u64,
                    byte: b,
                });
            }
            Ok(())
        }
    }
}

#[inline]
fn validate_record_multiline(
    mode: ValidationMode,
    parsed: &ParsedRecord<'_>,
    seq_segs: &[Segment],
    qual_segs: &[Segment],
) -> Result<(), FastqError> {
    match mode {
        ValidationMode::None => Ok(()),
        ValidationMode::Bases => {
            if let Err(idx) = validate_bases(parsed.record.seq()) {
                let b = parsed.record.seq()[idx];
                let offset = map_offset(seq_segs, idx, parsed.seq_start);
                return Err(FastqError::InvalidBase { offset, byte: b });
            }
            Ok(())
        }
        ValidationMode::Qualities => {
            if let Err(idx) = validate_qual(parsed.record.qual()) {
                let b = parsed.record.qual()[idx];
                let offset = map_offset(qual_segs, idx, parsed.qual_start);
                return Err(FastqError::InvalidQuality { offset, byte: b });
            }
            Ok(())
        }
        ValidationMode::BasesAndQualities => {
            if let Err(idx) = validate_bases(parsed.record.seq()) {
                let b = parsed.record.seq()[idx];
                let offset = map_offset(seq_segs, idx, parsed.seq_start);
                return Err(FastqError::InvalidBase { offset, byte: b });
            }
            if let Err(idx) = validate_qual(parsed.record.qual()) {
                let b = parsed.record.qual()[idx];
                let offset = map_offset(qual_segs, idx, parsed.qual_start);
                return Err(FastqError::InvalidQuality { offset, byte: b });
            }
            Ok(())
        }
    }
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

pub struct RecordsIter<'a> {
    reader: *mut FastqReader,
    _marker: PhantomData<&'a mut FastqReader>,
}

impl<'a> Iterator for RecordsIter<'a> {
    type Item = Result<FastqRecord<'a>, FastqError>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        // SAFETY: RecordsIter owns the exclusive borrow of FastqReader for 'a.
        let reader = unsafe { &mut *self.reader };
        match reader.next() {
            Ok(Some(rec)) => Some(Ok(rec)),
            Ok(None) => None,
            Err(err) => Some(Err(err)),
        }
    }
}

#[inline]
fn is_gz(path: &Path) -> bool {
    match path.extension() {
        Some(ext) => ext.eq_ignore_ascii_case("gz"),
        None => false,
    }
}

#[inline]
fn is_bgzf(path: &Path) -> bool {
    match path.extension() {
        Some(ext) => ext.eq_ignore_ascii_case("bgz") || ext.eq_ignore_ascii_case("bgzf"),
        None => false,
    }
}

impl FastqReader {
    pub fn from_bgzf_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        let path = path.as_ref();
        let backend = Backend::Bgzf(BgzfBackend::new(path)?);
        Ok(Self {
            backend,
            parser: FastqParser::new(),
            parser_multi: MultiLineFastqParser::new(),
            pos: 0,
            validation: ValidationMode::None,
            format: FastqFormat::SingleLine,
            ml_seq_buf: Vec::new(),
            ml_qual_buf: Vec::new(),
            ml_seq_segs: Vec::new(),
            ml_qual_segs: Vec::new(),
            resync: false,
        })
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

#[inline]
fn resync_plain(buf: &[u8], pos: &mut usize) {
    let mut p = *pos;
    if p == 0 && p < buf.len() && buf[p] == b'@' {
        *pos = p;
        return;
    }
    loop {
        if p >= buf.len() {
            *pos = p;
            return;
        }
        if p == 0 || buf[p - 1] == b'\n' {
            if buf[p] == b'@' {
                *pos = p;
                return;
            }
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

fn resync_bgzf(bgzf: &mut BgzfBackend) -> Result<(), FastqError> {
    loop {
        let slice = bgzf.available_slice();
        if slice.is_empty() {
            if !bgzf.refill()? {
                return Ok(());
            }
            continue;
        }
        if (bgzf.logical_offset() == 0 || slice[0] == b'@') && slice[0] == b'@' {
            return Ok(());
        }
        if let Some(lf) = find_lf(slice, 0) {
            let next_pos = lf + 1;
            bgzf.advance(next_pos);
            let next = bgzf.available_slice();
            if !next.is_empty() && next[0] == b'@' {
                return Ok(());
            }
            continue;
        } else {
            let len = slice.len();
            bgzf.advance(len);
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum CompressionKind {
    Plain,
    Gzip,
    Bgzf,
}

fn detect_compression(path: &Path) -> Result<CompressionKind, std::io::Error> {
    let mut file = File::open(path)?;
    let mut buf = [0u8; 128];
    let n = file.read(&mut buf)?;
    if n < 2 {
        return Ok(CompressionKind::Plain);
    }
    if buf[0] != 0x1f || buf[1] != 0x8b {
        return Ok(CompressionKind::Plain);
    }
    if n < 10 {
        return Ok(CompressionKind::Gzip);
    }
    let flg = buf[3];
    if (flg & 0x04) == 0 {
        return Ok(CompressionKind::Gzip);
    }
    if n < 12 {
        return Ok(CompressionKind::Gzip);
    }
    let xlen = u16::from_le_bytes([buf[10], buf[11]]) as usize;
    let extra_end = 12 + xlen;
    if extra_end > n {
        return Ok(CompressionKind::Gzip);
    }
    let extra = &buf[12..extra_end];
    for i in memchr_iter(b'B', extra) {
        if i + 3 < extra.len() && extra[i + 1] == b'C' {
            let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
            if slen == 2 {
                return Ok(CompressionKind::Bgzf);
            }
        }
    }
    Ok(CompressionKind::Gzip)
}
