use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::{GzipBackend, LineStatus};
use crate::backend::stream::StreamBackend;
use crate::error::{FastqError, InvalidKind};
use crate::parser::{ParsedRecord, RecordScratch, Segment};
use crate::record::FastqRecord;
use crate::simd::newline::find_lf;

#[cfg(feature = "noodles-bgzf")]
use crate::backend::noodles_bgzf::NoodlesBgzfBackend;

pub trait LineSource {
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

pub struct MultiLineFastqParser;

impl MultiLineFastqParser {
    #[inline]
    pub fn new() -> Self {
        Self
    }

    pub fn next_record_in_slice<'a>(
        &mut self,
        buf: &'a [u8],
        pos: &mut usize,
        scratch: &'a mut RecordScratch,
        seq_segs: &mut Vec<Segment>,
        qual_segs: &mut Vec<Segment>,
    ) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        let header_start = *pos as u64;
        let header = match next_line_in_slice(buf, pos) {
            Line::Line(line) => line,
            Line::EofClean => return Ok(None),
            Line::EofPartial => {
                return Err(FastqError::UnexpectedEof {
                    offset: header_start,
                });
            }
        };
        if header.is_empty() || header[0] != b'@' {
            return Err(FastqError::InvalidFormat {
                offset: header_start,
                kind: InvalidKind::HeaderMissingAt,
            });
        }
        scratch.header.clear();
        scratch.header.extend_from_slice(&header[1..]);

        scratch.seq.clear();
        seq_segs.clear();
        let seq_start = *pos as u64;
        let plus_start;
        loop {
            let line_start = *pos as u64;
            let line = match next_line_in_slice(buf, pos) {
                Line::Line(line) => line,
                Line::EofClean | Line::EofPartial => {
                    return Err(FastqError::UnexpectedEof { offset: line_start });
                }
            };
            if line.is_empty() {
                return Err(FastqError::InvalidFormat {
                    offset: line_start,
                    kind: InvalidKind::SeqLineEmpty,
                });
            }
            if line[0] == b'+' {
                plus_start = line_start;
                break;
            }
            scratch.seq.extend_from_slice(line);
            seq_segs.push(Segment {
                offset: line_start,
                len: line.len(),
            });
        }

        if scratch.seq.is_empty() {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::SeqLineEmpty,
            });
        }

        let qual_start = *pos as u64;
        scratch.qual.clear();
        qual_segs.clear();
        let mut remaining = scratch.seq.len();
        while remaining > 0 {
            let line_start = *pos as u64;
            let line = match next_line_in_slice(buf, pos) {
                Line::Line(line) => line,
                Line::EofClean | Line::EofPartial => {
                    return Err(FastqError::UnexpectedEof { offset: line_start });
                }
            };
            if line.is_empty() {
                return Err(FastqError::InvalidFormat {
                    offset: line_start,
                    kind: InvalidKind::QualLineEmpty,
                });
            }
            if line.len() > remaining {
                return Err(FastqError::LengthMismatch {
                    offset: line_start,
                    seq_len: scratch.seq.len(),
                    qual_len: scratch.qual.len() + line.len(),
                });
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

    pub fn next_record_stream<'a, B: LineSource>(
        &mut self,
        backend: &mut B,
        scratch: &'a mut RecordScratch,
        seq_segs: &mut Vec<Segment>,
        qual_segs: &mut Vec<Segment>,
    ) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        let header_start = backend.logical_offset();
        match backend.read_line(&mut scratch.header)? {
            LineStatus::Line => {}
            LineStatus::EofClean => return Ok(None),
            LineStatus::EofPartial => {
                return Err(FastqError::UnexpectedEof {
                    offset: header_start,
                });
            }
        }
        if scratch.header.is_empty() || scratch.header[0] != b'@' {
            return Err(FastqError::InvalidFormat {
                offset: header_start,
                kind: InvalidKind::HeaderMissingAt,
            });
        }
        scratch.header.drain(..1);

        scratch.seq.clear();
        seq_segs.clear();
        let seq_start = backend.logical_offset();
        let plus_start;
        // `mem::take` keeps disjoint mutable borrows on scratch fields possible.
        let mut tmp = std::mem::take(&mut scratch.plus);
        loop {
            let line_start = backend.logical_offset();
            tmp.clear();
            match backend.read_line(&mut tmp)? {
                LineStatus::Line => {}
                LineStatus::EofClean | LineStatus::EofPartial => {
                    scratch.plus = tmp;
                    return Err(FastqError::UnexpectedEof { offset: line_start });
                }
            }
            if tmp.is_empty() {
                scratch.plus = tmp;
                return Err(FastqError::InvalidFormat {
                    offset: line_start,
                    kind: InvalidKind::SeqLineEmpty,
                });
            }
            if tmp[0] == b'+' {
                plus_start = line_start;
                break;
            }
            scratch.seq.extend_from_slice(&tmp);
            seq_segs.push(Segment {
                offset: line_start,
                len: tmp.len(),
            });
        }
        scratch.plus = tmp;

        if scratch.seq.is_empty() {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::SeqLineEmpty,
            });
        }

        let qual_start = backend.logical_offset();
        scratch.qual.clear();
        qual_segs.clear();
        let mut remaining = scratch.seq.len();
        let mut tmp = std::mem::take(&mut scratch.plus);
        while remaining > 0 {
            let line_start = backend.logical_offset();
            tmp.clear();
            match backend.read_line(&mut tmp)? {
                LineStatus::Line => {}
                LineStatus::EofClean | LineStatus::EofPartial => {
                    scratch.plus = tmp;
                    return Err(FastqError::UnexpectedEof { offset: line_start });
                }
            }
            if tmp.is_empty() {
                scratch.plus = tmp;
                return Err(FastqError::InvalidFormat {
                    offset: line_start,
                    kind: InvalidKind::QualLineEmpty,
                });
            }
            if tmp.len() > remaining {
                let qlen = scratch.qual.len() + tmp.len();
                scratch.plus = tmp;
                return Err(FastqError::LengthMismatch {
                    offset: line_start,
                    seq_len: scratch.seq.len(),
                    qual_len: qlen,
                });
            }
            scratch.qual.extend_from_slice(&tmp);
            qual_segs.push(Segment {
                offset: line_start,
                len: tmp.len(),
            });
            remaining -= tmp.len();
        }
        scratch.plus = tmp;

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(&scratch.header, &scratch.seq, &scratch.qual),
            header_start,
            seq_start,
            qual_start,
        }))
    }
}

#[allow(clippy::enum_variant_names)]
enum Line<'a> {
    Line(&'a [u8]),
    EofClean,
    EofPartial,
}

#[inline]
fn next_line_in_slice<'a>(buf: &'a [u8], pos: &mut usize) -> Line<'a> {
    let start = *pos;
    if start >= buf.len() {
        return Line::EofClean;
    }
    let lf = match find_lf(buf, start) {
        Some(idx) => idx,
        None => return Line::EofPartial,
    };
    let mut end = lf;
    if end > start && buf[end - 1] == b'\r' {
        end -= 1;
    }
    *pos = lf + 1;
    Line::Line(&buf[start..end])
}
