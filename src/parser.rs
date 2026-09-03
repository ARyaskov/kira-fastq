//! Single-line FASTQ parsing.
//!
//! Two implementations share one set of rules: a zero-copy one over a contiguous slice (mmap or
//! an in-memory buffer) and a generic one over any streaming line source.
//!
//! ## What counts as a valid record
//!
//! The rules follow what `kseq.h`, seqtk, seqkit and BioPython accept, because a reader that
//! rejects files every other tool reads is not usable in a pipeline:
//!
//! * The final record may be missing its trailing newline.
//! * Blank lines between records are skipped.
//! * A record may have an empty sequence and an empty quality line. Adapter trimmers such as
//!   cutadapt, fastp and Trimmomatic emit zero-length reads, and aligners accept them.
//! * `\r\n` is accepted anywhere `\n` is.
//!
//! Sequence and quality length must still agree, the header must start with `@` and the third
//! line with `+`.

use crate::backend::LineStatus;
use crate::error::{FastqError, InvalidKind};
use crate::multiline::LineSource;
use crate::record::FastqRecord;
use crate::simd::newline::find_lf;

pub(crate) struct ParsedRecord<'a> {
    pub(crate) record: FastqRecord<'a>,
    pub(crate) header_start: u64,
    pub(crate) seq_start: u64,
    pub(crate) qual_start: u64,
}

/// Byte range of one line of a multi-line field, used to map an index inside the joined
/// sequence back to a file offset.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Segment {
    pub(crate) offset: u64,
    pub(crate) len: usize,
}

#[derive(Default)]
pub(crate) struct RecordScratch {
    pub(crate) header: Vec<u8>,
    pub(crate) seq: Vec<u8>,
    pub(crate) plus: Vec<u8>,
    pub(crate) qual: Vec<u8>,
}

impl RecordScratch {
    #[inline]
    pub(crate) fn new() -> Self {
        Self::default()
    }
}

pub(crate) struct FastqParser;

impl FastqParser {
    #[inline]
    pub(crate) fn new() -> Self {
        Self
    }

    /// Zero-copy parse over a contiguous slice. Records borrow out of `buf`.
    pub(crate) fn next_record_in_slice<'a>(
        &mut self,
        buf: &'a [u8],
        pos: &mut usize,
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
        let header = &header[1..];

        let seq_start = *pos as u64;
        let seq = match next_line_in_slice(buf, pos) {
            SliceLine::Line(line) => line,
            SliceLine::Eof => return Err(FastqError::eof(seq_start)),
        };

        let plus_start = *pos as u64;
        // Fast path: the separator is almost always a bare "+".
        match buf.get(*pos) {
            Some(b'+') => {
                if buf.get(*pos + 1) == Some(&b'\n') {
                    *pos += 2;
                } else if buf.get(*pos + 1) == Some(&b'\r') && buf.get(*pos + 2) == Some(&b'\n') {
                    *pos += 3;
                } else if next_line_in_slice(buf, pos) == SliceLine::Eof {
                    return Err(FastqError::eof(plus_start));
                }
            }
            Some(_) => return Err(FastqError::invalid(plus_start, InvalidKind::PlusMissing)),
            None => return Err(FastqError::eof(plus_start)),
        }

        // The quality line is as long as the sequence, which invites locating its end by
        // arithmetic instead of by scanning. That is not safe: on a malformed record the byte at
        // `qual_start + seq_len` can be a newline belonging to a *later* line, and the record
        // would silently absorb it instead of failing. The scan stays.
        let qual_start = *pos;
        let qual = match next_line_in_slice(buf, pos) {
            SliceLine::Line(line) => line,
            SliceLine::Eof => return Err(FastqError::eof(qual_start as u64)),
        };

        if seq.len() != qual.len() {
            return Err(FastqError::length_mismatch(
                qual_start as u64,
                seq.len(),
                qual.len(),
            ));
        }

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(header, seq, qual),
            header_start,
            seq_start,
            qual_start: qual_start as u64,
        }))
    }

    /// Parse from any [`LineSource`] (gzip, BGZF, a stream, or the noodles adapter). Records
    /// borrow out of `scratch`.
    pub(crate) fn next_record_stream<'a, B: LineSource>(
        &mut self,
        backend: &mut B,
        scratch: &'a mut RecordScratch,
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
            // Blank line between records.
            header_start = backend.logical_offset();
        }
        if scratch.header[0] != b'@' {
            return Err(FastqError::invalid(
                header_start,
                InvalidKind::HeaderMissingAt,
            ));
        }
        scratch.header.drain(..1);

        let seq_start = backend.logical_offset();
        if backend.read_line(&mut scratch.seq)? == LineStatus::EofClean {
            return Err(FastqError::eof(seq_start));
        }

        let plus_start = backend.logical_offset();
        if backend.read_line(&mut scratch.plus)? == LineStatus::EofClean {
            return Err(FastqError::eof(plus_start));
        }
        if scratch.plus.first() != Some(&b'+') {
            return Err(FastqError::invalid(plus_start, InvalidKind::PlusMissing));
        }

        let qual_start = backend.logical_offset();
        if backend.read_line(&mut scratch.qual)? == LineStatus::EofClean {
            return Err(FastqError::eof(qual_start));
        }

        if scratch.seq.len() != scratch.qual.len() {
            return Err(FastqError::length_mismatch(
                qual_start,
                scratch.seq.len(),
                scratch.qual.len(),
            ));
        }

        Ok(Some(ParsedRecord {
            record: FastqRecord::new(&scratch.header, &scratch.seq, &scratch.qual),
            header_start,
            seq_start,
            qual_start,
        }))
    }
}

#[derive(PartialEq, Eq)]
pub(crate) enum SliceLine<'a> {
    /// A complete line. The terminator is consumed; a missing terminator at end of input is
    /// accepted, since that is how most files end in practice.
    Line(&'a [u8]),
    Eof,
}

/// Advance `pos` past any empty lines.
#[inline]
pub(crate) fn skip_blank_lines(buf: &[u8], pos: &mut usize) {
    loop {
        match buf.get(*pos) {
            Some(b'\n') => *pos += 1,
            Some(b'\r') if buf.get(*pos + 1) == Some(&b'\n') => *pos += 2,
            _ => return,
        }
    }
}

#[inline]
pub(crate) fn next_line_in_slice<'a>(buf: &'a [u8], pos: &mut usize) -> SliceLine<'a> {
    let start = *pos;
    if start >= buf.len() {
        return SliceLine::Eof;
    }
    match find_lf(buf, start) {
        Some(lf) => {
            let mut end = lf;
            if end > start && buf[end - 1] == b'\r' {
                end -= 1;
            }
            *pos = lf + 1;
            SliceLine::Line(&buf[start..end])
        }
        None => {
            // Final line without a terminator.
            let mut end = buf.len();
            if end > start && buf[end - 1] == b'\r' {
                end -= 1;
            }
            *pos = buf.len();
            SliceLine::Line(&buf[start..end])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (header, sequence, quality) of one parsed record.
    type OwnedRecord = (Vec<u8>, Vec<u8>, Vec<u8>);

    fn parse_all(data: &[u8]) -> Result<Vec<OwnedRecord>, FastqError> {
        let mut parser = FastqParser::new();
        let mut pos = 0usize;
        let mut out = Vec::new();
        while let Some(parsed) = parser.next_record_in_slice(data, &mut pos)? {
            out.push((
                parsed.record.header().to_vec(),
                parsed.record.seq().to_vec(),
                parsed.record.qual().to_vec(),
            ));
        }
        Ok(out)
    }

    #[test]
    fn accepts_missing_final_newline() {
        let recs = parse_all(b"@a\nAC\n+\n!!\n@b\nGT\n+\n##").expect("parse");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[1].2, b"##");
    }

    #[test]
    fn accepts_empty_reads() {
        let recs = parse_all(b"@a\n\n+\n\n@b\nGT\n+\n##\n").expect("parse");
        assert_eq!(recs.len(), 2);
        assert!(recs[0].1.is_empty());
        assert!(recs[0].2.is_empty());
    }

    #[test]
    fn skips_blank_lines_between_and_after_records() {
        let recs = parse_all(b"\n@a\nAC\n+\n!!\n\n\n@b\nGT\n+\n##\n\n").expect("parse");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].0, b"a");
        assert_eq!(recs[1].0, b"b");
    }

    #[test]
    fn accepts_crlf_and_plus_with_id() {
        let recs = parse_all(b"@a\r\nAC\r\n+a\r\n!!\r\n").expect("parse");
        assert_eq!(recs, vec![(b"a".to_vec(), b"AC".to_vec(), b"!!".to_vec())]);
    }

    #[test]
    fn rejects_length_mismatch_even_when_next_line_aligns() {
        // The byte at qual_start + seq_len is a newline here, which a length-based fast path
        // would have mistaken for the end of the quality line.
        let err = parse_all(b"@a\nACGT\n+\n!!\n@\n").expect_err("must reject");
        match err {
            FastqError::LengthMismatch {
                seq_len, qual_len, ..
            } => {
                assert_eq!((seq_len, qual_len), (4, 2));
            }
            other => panic!("unexpected error: {other}"),
        }
    }

    #[test]
    fn rejects_truncated_record() {
        let err = parse_all(b"@a\nACGT\n+\n").expect_err("must reject");
        assert!(matches!(err, FastqError::UnexpectedEof { .. }));
    }

    #[test]
    fn rejects_missing_at_and_plus() {
        assert!(matches!(
            parse_all(b"a\nACGT\n+\n!!!!\n"),
            Err(FastqError::InvalidFormat {
                kind: InvalidKind::HeaderMissingAt,
                ..
            })
        ));
        assert!(matches!(
            parse_all(b"@a\nACGT\nx\n!!!!\n"),
            Err(FastqError::InvalidFormat {
                kind: InvalidKind::PlusMissing,
                ..
            })
        ));
    }
}
