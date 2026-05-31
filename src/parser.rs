use crate::backend::gzip::LineStatus;
use crate::error::{FastqError, InvalidKind};
use crate::multiline::LineSource;
use crate::record::FastqRecord;
use crate::simd::newline::find_lf;

pub struct ParsedRecord<'a> {
    pub record: FastqRecord<'a>,
    pub header_start: u64,
    pub seq_start: u64,
    pub qual_start: u64,
}

#[derive(Debug, Clone, Copy)]
pub struct Segment {
    pub offset: u64,
    pub len: usize,
}

#[derive(Default)]
pub struct RecordScratch {
    pub header: Vec<u8>,
    pub seq: Vec<u8>,
    pub plus: Vec<u8>,
    pub qual: Vec<u8>,
}

impl RecordScratch {
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    #[inline]
    pub fn clear(&mut self) {
        self.header.clear();
        self.seq.clear();
        self.plus.clear();
        self.qual.clear();
    }
}

pub struct FastqParser;

impl FastqParser {
    #[inline]
    pub fn new() -> Self {
        Self
    }

    /// Zero-copy single-line parse over a contiguous slice (mmap, in-memory buffer).
    pub fn next_record_in_slice<'a>(
        &mut self,
        buf: &'a [u8],
        pos: &mut usize,
    ) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        let header_start = *pos as u64;
        let header = match next_line_in_slice(buf, pos) {
            LineInSlice::Line(line) => line,
            LineInSlice::EofClean => return Ok(None),
            LineInSlice::EofPartial => {
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
        let header = &header[1..];

        let seq_start = *pos as u64;
        let seq = match next_line_in_slice(buf, pos) {
            LineInSlice::Line(line) => line,
            LineInSlice::EofClean | LineInSlice::EofPartial => {
                return Err(FastqError::UnexpectedEof { offset: seq_start });
            }
        };
        if seq.is_empty() {
            return Err(FastqError::InvalidFormat {
                offset: seq_start,
                kind: InvalidKind::SeqLineEmpty,
            });
        }

        let plus_start = *pos as u64;
        let plus = match next_line_in_slice(buf, pos) {
            LineInSlice::Line(line) => line,
            LineInSlice::EofClean | LineInSlice::EofPartial => {
                return Err(FastqError::UnexpectedEof { offset: plus_start });
            }
        };
        if plus.is_empty() || plus[0] != b'+' {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::PlusMissing,
            });
        }

        let qual_start = *pos as u64;
        let qual = match next_line_in_slice(buf, pos) {
            LineInSlice::Line(line) => line,
            LineInSlice::EofClean | LineInSlice::EofPartial => {
                return Err(FastqError::UnexpectedEof { offset: qual_start });
            }
        };

        if seq.len() != qual.len() {
            return Err(FastqError::LengthMismatch {
                offset: qual_start,
                seq_len: seq.len(),
                qual_len: qual.len(),
            });
        }

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(header, seq, qual),
            header_start,
            seq_start,
            qual_start,
        }))
    }

    /// Generic single-line parser over any [`LineSource`] (Gzip, Bgzf, Stream, or the
    /// `noodles-bgzf` adapter when enabled). Records borrow into `scratch`.
    pub fn next_record_stream<'a, B: LineSource>(
        &mut self,
        backend: &mut B,
        scratch: &'a mut RecordScratch,
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

        let seq_start = backend.logical_offset();
        if backend.read_line(&mut scratch.seq)? != LineStatus::Line {
            return Err(FastqError::UnexpectedEof { offset: seq_start });
        }
        if scratch.seq.is_empty() {
            return Err(FastqError::InvalidFormat {
                offset: seq_start,
                kind: InvalidKind::SeqLineEmpty,
            });
        }

        let plus_start = backend.logical_offset();
        if backend.read_line(&mut scratch.plus)? != LineStatus::Line {
            return Err(FastqError::UnexpectedEof { offset: plus_start });
        }
        if scratch.plus.is_empty() || scratch.plus[0] != b'+' {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::PlusMissing,
            });
        }

        let qual_start = backend.logical_offset();
        if backend.read_line(&mut scratch.qual)? != LineStatus::Line {
            return Err(FastqError::UnexpectedEof { offset: qual_start });
        }

        if scratch.seq.len() != scratch.qual.len() {
            return Err(FastqError::LengthMismatch {
                offset: qual_start,
                seq_len: scratch.seq.len(),
                qual_len: scratch.qual.len(),
            });
        }

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(&scratch.header, &scratch.seq, &scratch.qual),
            header_start,
            seq_start,
            qual_start,
        }))
    }
}

#[allow(clippy::enum_variant_names)]
enum LineInSlice<'a> {
    Line(&'a [u8]),
    EofClean,
    EofPartial,
}

#[inline]
fn next_line_in_slice<'a>(buf: &'a [u8], pos: &mut usize) -> LineInSlice<'a> {
    let start = *pos;
    if start >= buf.len() {
        return LineInSlice::EofClean;
    }
    let lf = match find_lf(buf, start) {
        Some(idx) => idx,
        None => return LineInSlice::EofPartial,
    };
    let mut end = lf;
    if end > start && buf[end - 1] == b'\r' {
        end -= 1;
    }
    *pos = lf + 1;
    LineInSlice::Line(&buf[start..end])
}
