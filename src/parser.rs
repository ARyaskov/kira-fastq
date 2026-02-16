use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::GzipBackend;
use crate::error::{FastqError, InvalidKind};
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

pub struct FastqParser;

impl FastqParser {
    #[inline]
    pub fn new() -> Self {
        Self
    }

    pub fn next_record<'a>(
        &mut self,
        buf: &'a [u8],
        pos: &mut usize,
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

        let seq_start = *pos;
        let seq = match next_line(buf, pos) {
            LineResult::Line(line) => line,
            LineResult::EofNoBytes | LineResult::EofPartial => {
                return Err(FastqError::UnexpectedEof {
                    offset: seq_start as u64,
                });
            }
        };

        let plus_start = *pos;
        let plus = match next_line(buf, pos) {
            LineResult::Line(line) => line,
            LineResult::EofNoBytes | LineResult::EofPartial => {
                return Err(FastqError::UnexpectedEof {
                    offset: plus_start as u64,
                });
            }
        };
        if plus.is_empty() || plus[0] != b'+' {
            return Err(FastqError::InvalidFormat {
                offset: plus_start as u64,
                kind: InvalidKind::PlusMissing,
            });
        }

        let qual_start = *pos;
        let qual = match next_line(buf, pos) {
            LineResult::Line(line) => line,
            LineResult::EofNoBytes | LineResult::EofPartial => {
                return Err(FastqError::UnexpectedEof {
                    offset: qual_start as u64,
                });
            }
        };

        if seq.len() != qual.len() {
            return Err(FastqError::LengthMismatch {
                offset: qual_start as u64,
                seq_len: seq.len(),
                qual_len: qual.len(),
            });
        }

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(header, seq, qual),
            header_start: header_start as u64,
            seq_start: seq_start as u64,
            qual_start: qual_start as u64,
        }))
    }

    pub fn next_record_gzip<'a>(
        &mut self,
        backend: &'a mut GzipBackend,
    ) -> Result<Option<ParsedRecord<'a>>, FastqError> {
        // SAFETY: backend pointer is valid for this call and returned slices borrow backend data.
        unsafe { next_record_gzip_fast(backend as *mut GzipBackend) }
    }

    pub fn next_record_bgzf<'a>(
        &mut self,
        backend: &'a mut BgzfBackend,
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

        // SAFETY: backend_ptr is valid for the duration of this call.
        let seq_start = unsafe { (*backend_ptr).logical_offset() };
        let seq = match unsafe { next_line_bgzf(backend_ptr)? } {
            LineResult::Line(line) => line,
            LineResult::EofNoBytes | LineResult::EofPartial => {
                return Err(FastqError::UnexpectedEof { offset: seq_start });
            }
        };

        // SAFETY: backend_ptr is valid for the duration of this call.
        let plus_start = unsafe { (*backend_ptr).logical_offset() };
        let plus = match unsafe { next_line_bgzf(backend_ptr)? } {
            LineResult::Line(line) => line,
            LineResult::EofNoBytes | LineResult::EofPartial => {
                return Err(FastqError::UnexpectedEof { offset: plus_start });
            }
        };
        if plus.is_empty() || plus[0] != b'+' {
            return Err(FastqError::InvalidFormat {
                offset: plus_start,
                kind: InvalidKind::PlusMissing,
            });
        }

        // SAFETY: backend_ptr is valid for the duration of this call.
        let qual_start = unsafe { (*backend_ptr).logical_offset() };
        let qual = match unsafe { next_line_bgzf(backend_ptr)? } {
            LineResult::Line(line) => line,
            LineResult::EofNoBytes | LineResult::EofPartial => {
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
}

#[inline]
unsafe fn next_record_gzip_fast<'a>(
    backend: *mut GzipBackend,
) -> Result<Option<ParsedRecord<'a>>, FastqError> {
    loop {
        // SAFETY: backend pointer is valid and available_slice borrows from backend buffer.
        let slice = unsafe { (*backend).available_slice() };
        if !slice.is_empty() {
            // SAFETY: backend pointer is valid for logical offset read.
            let base_offset = unsafe { (*backend).logical_offset() };
            if let Some((record, consumed)) = parse_singleline_record_in_slice(slice, base_offset)?
            {
                // SAFETY: advancing after parsing keeps record references valid because data stays in backend buffer.
                unsafe { (*backend).advance(consumed) };
                return Ok(Some(record));
            }
            // Record crosses buffer boundary (or is malformed in a way fast path cannot decide): use robust slow path.
            return unsafe { next_record_gzip_slow(backend) };
        }

        // SAFETY: backend pointer is valid for refill.
        if !unsafe { (*backend).refill()? } {
            return Ok(None);
        }
    }
}

#[inline]
fn parse_singleline_record_in_slice<'a>(
    slice: &'a [u8],
    base_offset: u64,
) -> Result<Option<(ParsedRecord<'a>, usize)>, FastqError> {
    let lf1 = match find_lf(slice, 0) {
        Some(v) => v,
        None => return Ok(None),
    };
    let seq_start_idx = lf1 + 1;
    let lf2 = match find_lf(slice, seq_start_idx) {
        Some(v) => v,
        None => return Ok(None),
    };
    let plus_start_idx = lf2 + 1;
    let lf3 = match find_lf(slice, plus_start_idx) {
        Some(v) => v,
        None => return Ok(None),
    };
    let qual_start_idx = lf3 + 1;
    let lf4 = match find_lf(slice, qual_start_idx) {
        Some(v) => v,
        None => return Ok(None),
    };

    let header_start = base_offset;
    let seq_start = base_offset + seq_start_idx as u64;
    let plus_start = base_offset + plus_start_idx as u64;
    let qual_start = base_offset + qual_start_idx as u64;

    let header = trim_cr(&slice[0..lf1]);
    if header.is_empty() || header[0] != b'@' {
        return Err(FastqError::InvalidFormat {
            offset: header_start,
            kind: InvalidKind::HeaderMissingAt,
        });
    }
    let header = &header[1..];

    let seq = trim_cr(&slice[seq_start_idx..lf2]);
    let plus = trim_cr(&slice[plus_start_idx..lf3]);
    if plus.is_empty() || plus[0] != b'+' {
        return Err(FastqError::InvalidFormat {
            offset: plus_start,
            kind: InvalidKind::PlusMissing,
        });
    }
    let qual = trim_cr(&slice[qual_start_idx..lf4]);
    if seq.len() != qual.len() {
        return Err(FastqError::LengthMismatch {
            offset: qual_start,
            seq_len: seq.len(),
            qual_len: qual.len(),
        });
    }

    let consumed = lf4 + 1;
    Ok(Some((
        ParsedRecord {
            record: FastqRecord::new(header, seq, qual),
            header_start,
            seq_start,
            qual_start,
        },
        consumed,
    )))
}

#[inline(always)]
fn trim_cr(line: &[u8]) -> &[u8] {
    match line.last() {
        Some(b'\r') => &line[..line.len() - 1],
        _ => line,
    }
}

enum LineResult<'a> {
    Line(&'a [u8]),
    EofNoBytes,
    EofPartial,
}

#[inline]
unsafe fn next_record_gzip_slow<'a>(
    backend: *mut GzipBackend,
) -> Result<Option<ParsedRecord<'a>>, FastqError> {
    // SAFETY: backend pointer is valid for duration of this function.
    let header_start = unsafe { (*backend).logical_offset() };
    let header = match unsafe { next_line_gzip(backend)? } {
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

    // SAFETY: backend pointer is valid for duration of this function.
    let seq_start = unsafe { (*backend).logical_offset() };
    let seq = match unsafe { next_line_gzip(backend)? } {
        LineResult::Line(line) => line,
        LineResult::EofNoBytes | LineResult::EofPartial => {
            return Err(FastqError::UnexpectedEof { offset: seq_start });
        }
    };

    // SAFETY: backend pointer is valid for duration of this function.
    let plus_start = unsafe { (*backend).logical_offset() };
    let plus = match unsafe { next_line_gzip(backend)? } {
        LineResult::Line(line) => line,
        LineResult::EofNoBytes | LineResult::EofPartial => {
            return Err(FastqError::UnexpectedEof { offset: plus_start });
        }
    };
    if plus.is_empty() || plus[0] != b'+' {
        return Err(FastqError::InvalidFormat {
            offset: plus_start,
            kind: InvalidKind::PlusMissing,
        });
    }

    // SAFETY: backend pointer is valid for duration of this function.
    let qual_start = unsafe { (*backend).logical_offset() };
    let qual = match unsafe { next_line_gzip(backend)? } {
        LineResult::Line(line) => line,
        LineResult::EofNoBytes | LineResult::EofPartial => {
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
