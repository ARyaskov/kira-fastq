use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::GzipBackend;
use crate::error::{FastqError, InvalidKind};
use crate::parser::{ParsedRecord, Segment};
use crate::record::FastqRecord;
use crate::simd::newline::find_lf;

pub struct MultiLineFastqParser;

impl MultiLineFastqParser {
    #[inline]
    pub fn new() -> Self {
        Self
    }

    pub fn next_record<'a>(
        &mut self,
        buf: &'a [u8],
        pos: &mut usize,
        seq_buf: &'a mut Vec<u8>,
        qual_buf: &'a mut Vec<u8>,
        seq_segs: &mut Vec<Segment>,
        qual_segs: &mut Vec<Segment>,
    ) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        let header_start = *pos;
        let header = match next_line(buf, pos) {
            LineResult::Line(line) => line,
            LineResult::EofNoBytes => return Ok(None),
            LineResult::EofPartial => {
                return Err(FastqError::UnexpectedEof {
                    offset: header_start as u64,
                });
            }
        };
        if header.is_empty() || header[0] != b'@' {
            return Err(FastqError::InvalidFormat {
                offset: header_start as u64,
                kind: InvalidKind::HeaderMissingAt,
            });
        }
        let header = &header[1..];

        seq_buf.clear();
        seq_segs.clear();
        let seq_start = *pos as u64;
        let plus_start = loop {
            let line_start = *pos;
            let line = match next_line(buf, pos) {
                LineResult::Line(line) => line,
                LineResult::EofNoBytes | LineResult::EofPartial => {
                    return Err(FastqError::UnexpectedEof {
                        offset: line_start as u64,
                    });
                }
            };
            if line.is_empty() {
                return Err(FastqError::InvalidFormat {
                    offset: line_start as u64,
                    kind: InvalidKind::SeqLineEmpty,
                });
            }
            if line[0] == b'+' {
                break line_start as u64;
            }
            seq_buf.extend_from_slice(line);
            seq_segs.push(Segment {
                offset: line_start as u64,
                len: line.len(),
            });
        };

        if seq_buf.is_empty() {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::SeqLineEmpty,
            });
        }

        let qual_start = *pos as u64;
        qual_buf.clear();
        qual_segs.clear();
        let mut remaining = seq_buf.len();
        while remaining > 0 {
            let line_start = *pos;
            let line = match next_line(buf, pos) {
                LineResult::Line(line) => line,
                LineResult::EofNoBytes | LineResult::EofPartial => {
                    return Err(FastqError::UnexpectedEof {
                        offset: line_start as u64,
                    });
                }
            };
            if line.is_empty() {
                return Err(FastqError::InvalidFormat {
                    offset: line_start as u64,
                    kind: InvalidKind::QualLineEmpty,
                });
            }
            if line.len() > remaining {
                return Err(FastqError::LengthMismatch {
                    offset: line_start as u64,
                    seq_len: seq_buf.len(),
                    qual_len: qual_buf.len() + line.len(),
                });
            }
            qual_buf.extend_from_slice(line);
            qual_segs.push(Segment {
                offset: line_start as u64,
                len: line.len(),
            });
            remaining -= line.len();
        }

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(header, seq_buf.as_slice(), qual_buf.as_slice()),
            header_start: header_start as u64,
            seq_start,
            qual_start,
        }))
    }

    pub fn next_record_gzip<'a>(
        &mut self,
        backend: &'a mut GzipBackend,
        seq_buf: &'a mut Vec<u8>,
        qual_buf: &'a mut Vec<u8>,
        seq_segs: &mut Vec<Segment>,
        qual_segs: &mut Vec<Segment>,
    ) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        let backend_ptr = backend as *mut GzipBackend;
        // SAFETY: backend_ptr is valid for the duration of this call.
        let header_start = unsafe { (*backend_ptr).logical_offset() };
        let header = match unsafe { next_line_gzip(backend_ptr)? } {
            LineResult::Line(line) => line,
            LineResult::EofNoBytes => return Ok(None),
            LineResult::EofPartial => {
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

        seq_buf.clear();
        seq_segs.clear();
        // SAFETY: backend_ptr is valid for the duration of this call.
        let seq_start = unsafe { (*backend_ptr).logical_offset() };
        let plus_start = loop {
            // SAFETY: backend_ptr is valid for the duration of this call.
            let line_start = unsafe { (*backend_ptr).logical_offset() };
            let line = match unsafe { next_line_gzip(backend_ptr)? } {
                LineResult::Line(line) => line,
                LineResult::EofNoBytes | LineResult::EofPartial => {
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
                break line_start;
            }
            seq_buf.extend_from_slice(line);
            seq_segs.push(Segment {
                offset: line_start,
                len: line.len(),
            });
        };

        if seq_buf.is_empty() {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::SeqLineEmpty,
            });
        }

        // SAFETY: backend_ptr is valid for the duration of this call.
        let qual_start = unsafe { (*backend_ptr).logical_offset() };
        qual_buf.clear();
        qual_segs.clear();
        let mut remaining = seq_buf.len();
        while remaining > 0 {
            // SAFETY: backend_ptr is valid for the duration of this call.
            let line_start = unsafe { (*backend_ptr).logical_offset() };
            let line = match unsafe { next_line_gzip(backend_ptr)? } {
                LineResult::Line(line) => line,
                LineResult::EofNoBytes | LineResult::EofPartial => {
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
                    seq_len: seq_buf.len(),
                    qual_len: qual_buf.len() + line.len(),
                });
            }
            qual_buf.extend_from_slice(line);
            qual_segs.push(Segment {
                offset: line_start,
                len: line.len(),
            });
            remaining -= line.len();
        }

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(header, seq_buf.as_slice(), qual_buf.as_slice()),
            header_start,
            seq_start,
            qual_start,
        }))
    }

    pub fn next_record_bgzf<'a>(
        &mut self,
        backend: &'a mut BgzfBackend,
        seq_buf: &'a mut Vec<u8>,
        qual_buf: &'a mut Vec<u8>,
        seq_segs: &mut Vec<Segment>,
        qual_segs: &mut Vec<Segment>,
    ) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        let backend_ptr = backend as *mut BgzfBackend;
        // SAFETY: backend_ptr is valid for the duration of this call.
        let header_start = unsafe { (*backend_ptr).logical_offset() };
        let header = match unsafe { next_line_bgzf(backend_ptr)? } {
            LineResult::Line(line) => line,
            LineResult::EofNoBytes => return Ok(None),
            LineResult::EofPartial => {
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

        seq_buf.clear();
        seq_segs.clear();
        // SAFETY: backend_ptr is valid for the duration of this call.
        let seq_start = unsafe { (*backend_ptr).logical_offset() };
        let plus_start = loop {
            // SAFETY: backend_ptr is valid for the duration of this call.
            let line_start = unsafe { (*backend_ptr).logical_offset() };
            let line = match unsafe { next_line_bgzf(backend_ptr)? } {
                LineResult::Line(line) => line,
                LineResult::EofNoBytes | LineResult::EofPartial => {
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
                break line_start;
            }
            seq_buf.extend_from_slice(line);
            seq_segs.push(Segment {
                offset: line_start,
                len: line.len(),
            });
        };

        if seq_buf.is_empty() {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::SeqLineEmpty,
            });
        }

        // SAFETY: backend_ptr is valid for the duration of this call.
        let qual_start = unsafe { (*backend_ptr).logical_offset() };
        qual_buf.clear();
        qual_segs.clear();
        let mut remaining = seq_buf.len();
        while remaining > 0 {
            // SAFETY: backend_ptr is valid for the duration of this call.
            let line_start = unsafe { (*backend_ptr).logical_offset() };
            let line = match unsafe { next_line_bgzf(backend_ptr)? } {
                LineResult::Line(line) => line,
                LineResult::EofNoBytes | LineResult::EofPartial => {
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
                    seq_len: seq_buf.len(),
                    qual_len: qual_buf.len() + line.len(),
                });
            }
            qual_buf.extend_from_slice(line);
            qual_segs.push(Segment {
                offset: line_start,
                len: line.len(),
            });
            remaining -= line.len();
        }

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(header, seq_buf.as_slice(), qual_buf.as_slice()),
            header_start,
            seq_start,
            qual_start,
        }))
    }
}

enum LineResult<'a> {
    Line(&'a [u8]),
    EofNoBytes,
    EofPartial,
}

#[inline]
fn next_line<'a>(buf: &'a [u8], pos: &mut usize) -> LineResult<'a> {
    let start = *pos;
    if start >= buf.len() {
        return LineResult::EofNoBytes;
    }
    let lf = match find_lf(buf, start) {
        Some(idx) => idx,
        None => return LineResult::EofPartial,
    };
    let mut end = lf;
    if end > start && buf[end - 1] == b'\r' {
        end -= 1;
    }
    *pos = lf + 1;
    LineResult::Line(&buf[start..end])
}

#[inline]
unsafe fn next_line_gzip<'a>(backend: *mut GzipBackend) -> Result<LineResult<'a>, FastqError> {
    loop {
        // SAFETY: backend is a valid pointer for the duration of this function.
        let slice = unsafe { (*backend).available_slice() };
        if !slice.is_empty() {
            if let Some(lf) = find_lf(slice, 0) {
                let mut end = lf;
                if end > 0 && slice[end - 1] == b'\r' {
                    end -= 1;
                }
                let ptr = slice.as_ptr();
                // SAFETY: advance does not invalidate `ptr` for the current slice.
                unsafe { (*backend).advance(lf + 1) };
                // SAFETY: `ptr` points to `slice` and `end` is within bounds.
                let line = unsafe { std::slice::from_raw_parts(ptr, end) };
                return Ok(LineResult::Line(line));
            }
            // For gzip, do not advance here; refill appends more output.
        }

        let had_bytes = !slice.is_empty();
        // SAFETY: backend is valid for refill.
        if !unsafe { (*backend).refill()? } {
            return Ok(if had_bytes {
                LineResult::EofPartial
            } else {
                LineResult::EofNoBytes
            });
        }
    }
}

#[inline]
unsafe fn next_line_bgzf<'a>(backend: *mut BgzfBackend) -> Result<LineResult<'a>, FastqError> {
    loop {
        // SAFETY: backend is a valid pointer for the duration of this function.
        let slice = unsafe { (*backend).available_slice() };
        if !slice.is_empty() {
            if let Some(lf) = find_lf(slice, 0) {
                let mut end = lf;
                if end > 0 && slice[end - 1] == b'\r' {
                    end -= 1;
                }
                let ptr = slice.as_ptr();
                // SAFETY: advance does not invalidate `ptr` for the current slice.
                unsafe { (*backend).advance(lf + 1) };
                // SAFETY: `ptr` points to `slice` and `end` is within bounds.
                let line = unsafe { std::slice::from_raw_parts(ptr, end) };
                return Ok(LineResult::Line(line));
            }
        }

        let had_bytes = !slice.is_empty();
        // SAFETY: backend is valid for refill.
        if !unsafe { (*backend).refill()? } {
            return Ok(if had_bytes {
                LineResult::EofPartial
            } else {
                LineResult::EofNoBytes
            });
        }
    }
}
