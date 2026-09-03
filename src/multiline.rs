//! Multi-line FASTQ parsing.
//!
//! Long-read basecallers and some older tools wrap sequence and quality over several lines. The
//! sequence runs until a line starting with `+`; the quality then runs until it is exactly as
//! long as the sequence, which is the only way to tell a wrapped quality line from the `@` of the
//! next record.
//!
//! Leniency matches [`crate::parser`]: blank lines between records are skipped, the final record
//! may lack its newline, and `\r\n` is accepted.

use crate::backend::LineStatus;
use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::GzipBackend;
use crate::backend::stream::StreamBackend;
use crate::error::{FastqError, InvalidKind};
use crate::parser::{
    ParsedRecord, RecordScratch, Segment, SliceLine, next_line_in_slice, skip_blank_lines,
};
use crate::record::FastqRecord;

#[cfg(feature = "noodles-bgzf")]
use crate::backend::noodles_bgzf::NoodlesBgzfBackend;

/// A source of lines: implemented by every streaming backend.
pub(crate) trait LineSource {
    fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError>;
    fn logical_offset(&self) -> u64;
}

impl LineSource for GzipBackend {
    #[inline]
    fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        GzipBackend::read_line(self, out)
    }
    #[inline]
    fn logical_offset(&self) -> u64 {
        GzipBackend::logical_offset(self)
    }
}

impl LineSource for BgzfBackend {
    #[inline]
    fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        BgzfBackend::read_line(self, out)
    }
    #[inline]
    fn logical_offset(&self) -> u64 {
        BgzfBackend::logical_offset(self)
    }
}

impl LineSource for StreamBackend {
    #[inline]
    fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        StreamBackend::read_line(self, out)
    }
    #[inline]
    fn logical_offset(&self) -> u64 {
        StreamBackend::logical_offset(self)
    }
}

#[cfg(feature = "noodles-bgzf")]
impl LineSource for NoodlesBgzfBackend {
    #[inline]
    fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        NoodlesBgzfBackend::read_line(self, out)
    }
    #[inline]
    fn logical_offset(&self) -> u64 {
        NoodlesBgzfBackend::logical_offset(self)
    }
}

pub(crate) struct MultiLineFastqParser;

impl MultiLineFastqParser {
    #[inline]
    pub(crate) fn new() -> Self {
        Self
    }

    pub(crate) fn next_record_in_slice<'a>(
        &mut self,
        buf: &'a [u8],
        pos: &mut usize,
        scratch: &'a mut RecordScratch,
        seq_segs: &mut Vec<Segment>,
        qual_segs: &mut Vec<Segment>,
    ) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        skip_blank_lines(buf, pos);

        let header_start = *pos as u64;
        let header = match next_line_in_slice(buf, pos) {
            SliceLine::Line(line) => line,
            SliceLine::Eof => return Ok(None),
        };
        if header.first() != Some(&b'@') {
            return Err(FastqError::invalid(
                header_start,
                InvalidKind::HeaderMissingAt,
            ));
        }
        scratch.header.clear();
        scratch.header.extend_from_slice(&header[1..]);

        scratch.seq.clear();
        seq_segs.clear();
        let seq_start = *pos as u64;
        loop {
            let line_start = *pos as u64;
            let line = match next_line_in_slice(buf, pos) {
                SliceLine::Line(line) => line,
                SliceLine::Eof => return Err(FastqError::eof(line_start)),
            };
            if line.first() == Some(&b'+') {
                break;
            }
            if line.is_empty() {
                return Err(FastqError::invalid(line_start, InvalidKind::SeqLineEmpty));
            }
            scratch.seq.extend_from_slice(line);
            seq_segs.push(Segment {
                offset: line_start,
                len: line.len(),
            });
        }

        let qual_start = *pos as u64;
        scratch.qual.clear();
        qual_segs.clear();
        let mut remaining = scratch.seq.len();
        if remaining == 0 {
            // Zero-length read: exactly one, empty, quality line.
            match next_line_in_slice(buf, pos) {
                SliceLine::Line(line) => {
                    if !line.is_empty() {
                        return Err(FastqError::length_mismatch(qual_start, 0, line.len()));
                    }
                }
                SliceLine::Eof => return Err(FastqError::eof(qual_start)),
            }
        }
        while remaining > 0 {
            let line_start = *pos as u64;
            let line = match next_line_in_slice(buf, pos) {
                SliceLine::Line(line) => line,
                SliceLine::Eof => return Err(FastqError::eof(line_start)),
            };
            if line.is_empty() {
                return Err(FastqError::invalid(line_start, InvalidKind::QualLineEmpty));
            }
            if line.len() > remaining {
                return Err(FastqError::length_mismatch(
                    line_start,
                    scratch.seq.len(),
                    scratch.qual.len() + line.len(),
                ));
            }
            scratch.qual.extend_from_slice(line);
            qual_segs.push(Segment {
                offset: line_start,
                len: line.len(),
            });
            remaining -= line.len();
        }

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(&scratch.header, &scratch.seq, &scratch.qual),
            header_start,
            seq_start,
            qual_start,
        }))
    }

    pub(crate) fn next_record_stream<'a, B: LineSource>(
        &mut self,
        backend: &mut B,
        scratch: &'a mut RecordScratch,
        seq_segs: &mut Vec<Segment>,
        qual_segs: &mut Vec<Segment>,
    ) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        let mut header_start = backend.logical_offset();
        loop {
            match backend.read_line(&mut scratch.header)? {
                LineStatus::Line | LineStatus::EofPartial => {}
                LineStatus::EofClean => return Ok(None),
            }
            if !scratch.header.is_empty() {
                break;
            }
            header_start = backend.logical_offset();
        }
        if scratch.header[0] != b'@' {
            return Err(FastqError::invalid(
                header_start,
                InvalidKind::HeaderMissingAt,
            ));
        }
        scratch.header.drain(..1);

        scratch.seq.clear();
        seq_segs.clear();
        let seq_start = backend.logical_offset();
        // `mem::take` keeps the borrows of the individual scratch fields disjoint.
        let mut tmp = std::mem::take(&mut scratch.plus);
        let outcome =
            read_multiline_body(backend, &mut tmp, scratch, seq_segs, qual_segs, seq_start);
        scratch.plus = tmp;
        let qual_start = outcome?;

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(&scratch.header, &scratch.seq, &scratch.qual),
            header_start,
            seq_start,
            qual_start,
        }))
    }
}

/// Read sequence lines up to the `+` line and then the quality lines, returning the offset the
/// quality field starts at. Split out so the caller can always restore its scratch buffer.
fn read_multiline_body<B: LineSource>(
    backend: &mut B,
    tmp: &mut Vec<u8>,
    scratch: &mut RecordScratch,
    seq_segs: &mut Vec<Segment>,
    qual_segs: &mut Vec<Segment>,
    _seq_start: u64,
) -> Result<u64, FastqError> {
    loop {
        let line_start = backend.logical_offset();
        if backend.read_line(tmp)? == LineStatus::EofClean {
            return Err(FastqError::eof(line_start));
        }
        if tmp.first() == Some(&b'+') {
            break;
        }
        if tmp.is_empty() {
            return Err(FastqError::invalid(line_start, InvalidKind::SeqLineEmpty));
        }
        scratch.seq.extend_from_slice(tmp);
        seq_segs.push(Segment {
            offset: line_start,
            len: tmp.len(),
        });
    }

    let qual_start = backend.logical_offset();
    scratch.qual.clear();
    qual_segs.clear();
    let mut remaining = scratch.seq.len();
    if remaining == 0 {
        if backend.read_line(tmp)? == LineStatus::EofClean {
            return Err(FastqError::eof(qual_start));
        }
        if !tmp.is_empty() {
            return Err(FastqError::length_mismatch(qual_start, 0, tmp.len()));
        }
        return Ok(qual_start);
    }
    while remaining > 0 {
        let line_start = backend.logical_offset();
        if backend.read_line(tmp)? == LineStatus::EofClean {
            return Err(FastqError::eof(line_start));
        }
        if tmp.is_empty() {
            return Err(FastqError::invalid(line_start, InvalidKind::QualLineEmpty));
        }
        if tmp.len() > remaining {
            return Err(FastqError::length_mismatch(
                line_start,
                scratch.seq.len(),
                scratch.qual.len() + tmp.len(),
            ));
        }
        scratch.qual.extend_from_slice(tmp);
        qual_segs.push(Segment {
            offset: line_start,
            len: tmp.len(),
        });
        remaining -= tmp.len();
    }
    Ok(qual_start)
}
