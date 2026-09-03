use std::fmt;
use std::io;

/// Errors returned by every reader and writer in this crate.
///
/// Parse errors carry two positions: `offset` is the byte offset of the offending line in the
/// *decoded* stream, and `record` is the 1-based index of the record being parsed (`0` when the
/// error is raised outside record context, e.g. by a decompressor).
#[derive(Debug)]
#[non_exhaustive]
pub enum FastqError {
    Io(io::Error),
    InvalidFormat {
        offset: u64,
        record: u64,
        kind: InvalidKind,
    },
    UnexpectedEof {
        offset: u64,
        record: u64,
    },
    LengthMismatch {
        offset: u64,
        record: u64,
        seq_len: usize,
        qual_len: usize,
    },
    InvalidBase {
        offset: u64,
        record: u64,
        byte: u8,
    },
    InvalidQuality {
        offset: u64,
        record: u64,
        byte: u8,
    },
    /// One mate file ran out of records before the other. `which` names the file that was
    /// exhausted first.
    PairedCountMismatch {
        which: PairedWhich,
        record: u64,
    },
    PairedIdMismatch {
        offset_r1: u64,
        offset_r2: u64,
        record: u64,
        id_r1: Box<[u8]>,
        id_r2: Box<[u8]>,
    },
    Unsupported(UnsupportedOperation),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum InvalidKind {
    HeaderMissingAt,
    PlusMissing,
    GzipHeader,
    GzipData,
    GzipTrailerCrc,
    GzipTrailerIsize,
    /// Input ended in the middle of a deflate stream.
    GzipTruncated,
    SeqLineEmpty,
    QualLineEmpty,
    BgzfHeader,
    BgzfBlock,
    BgzfVirtualOffset,
    BgzfBlockTooLarge,
    BgzfBlockCrc,
    BgzfBlockIsize,
    /// The 28-byte BGZF end-of-file marker is missing: the file is truncated or was never
    /// finalised. Disable the check with `FastqReader::with_bgzf_eof_check(false)`.
    BgzfMissingEofMarker,
    BufferOverflow,
    /// A record handed to the writer has a line break inside its header.
    HeaderContainsNewline,
    /// A record handed to the writer has a line break inside its sequence.
    SeqContainsNewline,
    /// A record handed to the writer has a line break inside its qualities.
    QualContainsNewline,
}

impl fmt::Display for InvalidKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let msg = match self {
            InvalidKind::HeaderMissingAt => "header line does not start with '@'",
            InvalidKind::PlusMissing => "separator line does not start with '+'",
            InvalidKind::GzipHeader => "malformed gzip header",
            InvalidKind::GzipData => "corrupt deflate stream",
            InvalidKind::GzipTrailerCrc => "gzip CRC32 mismatch",
            InvalidKind::GzipTrailerIsize => "gzip ISIZE mismatch",
            InvalidKind::GzipTruncated => "gzip stream ends mid-member",
            InvalidKind::SeqLineEmpty => "empty sequence line",
            InvalidKind::QualLineEmpty => "empty quality line",
            InvalidKind::BgzfHeader => "malformed BGZF block header",
            InvalidKind::BgzfBlock => "malformed BGZF block",
            InvalidKind::BgzfVirtualOffset => "virtual offset does not address a BGZF block",
            InvalidKind::BgzfBlockTooLarge => "BGZF block exceeds 64 KiB of uncompressed data",
            InvalidKind::BgzfBlockCrc => "BGZF block CRC32 mismatch",
            InvalidKind::BgzfBlockIsize => "BGZF block ISIZE mismatch",
            InvalidKind::BgzfMissingEofMarker => "BGZF end-of-file marker is missing",
            InvalidKind::BufferOverflow => "decompression buffer overflow",
            InvalidKind::HeaderContainsNewline => "header contains a line break",
            InvalidKind::SeqContainsNewline => "sequence contains a line break",
            InvalidKind::QualContainsNewline => "qualities contain a line break",
        };
        f.write_str(msg)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum UnsupportedOperation {
    /// `seek` on a source that cannot be positioned (gzip, arbitrary streams).
    Seek,
    /// zstd input without the `zstd` feature.
    Zstd,
    /// bzip2 input: not handled by this crate.
    Bzip2,
    /// xz/LZMA input: not handled by this crate.
    Xz,
    /// Compression level outside the range accepted by the codec.
    CompressionLevel,
    /// BGZF output on the async path.
    AsyncBgzf,
}

impl fmt::Display for UnsupportedOperation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let msg = match self {
            UnsupportedOperation::Seek => "seek is not supported by this source",
            UnsupportedOperation::Zstd => "zstd input requires the `zstd` feature of kira-fastq",
            UnsupportedOperation::Bzip2 => "bzip2 input is not supported; decompress it first",
            UnsupportedOperation::Xz => "xz input is not supported; decompress it first",
            UnsupportedOperation::CompressionLevel => "compression level out of range",
            UnsupportedOperation::AsyncBgzf => {
                "BGZF output is not supported on the async path; use the sync writer"
            }
        };
        f.write_str(msg)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PairedWhich {
    R1,
    R2,
}

impl FastqError {
    #[inline]
    pub(crate) fn invalid(offset: u64, kind: InvalidKind) -> Self {
        FastqError::InvalidFormat {
            offset,
            record: 0,
            kind,
        }
    }

    #[inline]
    pub(crate) fn eof(offset: u64) -> Self {
        FastqError::UnexpectedEof { offset, record: 0 }
    }

    #[inline]
    pub(crate) fn length_mismatch(offset: u64, seq_len: usize, qual_len: usize) -> Self {
        FastqError::LengthMismatch {
            offset,
            record: 0,
            seq_len,
            qual_len,
        }
    }

    #[inline]
    pub(crate) fn invalid_base(offset: u64, byte: u8) -> Self {
        FastqError::InvalidBase {
            offset,
            record: 0,
            byte,
        }
    }

    #[inline]
    pub(crate) fn invalid_quality(offset: u64, byte: u8) -> Self {
        FastqError::InvalidQuality {
            offset,
            record: 0,
            byte,
        }
    }

    /// Attach the 1-based index of the record that was being parsed.
    #[inline]
    pub(crate) fn with_record(mut self, n: u64) -> Self {
        match &mut self {
            FastqError::InvalidFormat { record, .. }
            | FastqError::UnexpectedEof { record, .. }
            | FastqError::LengthMismatch { record, .. }
            | FastqError::InvalidBase { record, .. }
            | FastqError::InvalidQuality { record, .. }
            | FastqError::PairedCountMismatch { record, .. }
            | FastqError::PairedIdMismatch { record, .. } => *record = n,
            FastqError::Io(_) | FastqError::Unsupported(_) => {}
        }
        self
    }

    /// 1-based index of the record this error refers to, if known.
    pub fn record(&self) -> Option<u64> {
        let n = match self {
            FastqError::InvalidFormat { record, .. }
            | FastqError::UnexpectedEof { record, .. }
            | FastqError::LengthMismatch { record, .. }
            | FastqError::InvalidBase { record, .. }
            | FastqError::InvalidQuality { record, .. }
            | FastqError::PairedCountMismatch { record, .. }
            | FastqError::PairedIdMismatch { record, .. } => *record,
            FastqError::Io(_) | FastqError::Unsupported(_) => 0,
        };
        (n > 0).then_some(n)
    }

    /// Byte offset in the decoded stream this error refers to, if known.
    pub fn offset(&self) -> Option<u64> {
        match self {
            FastqError::InvalidFormat { offset, .. }
            | FastqError::UnexpectedEof { offset, .. }
            | FastqError::LengthMismatch { offset, .. }
            | FastqError::InvalidBase { offset, .. }
            | FastqError::InvalidQuality { offset, .. } => Some(*offset),
            FastqError::PairedIdMismatch { offset_r1, .. } => Some(*offset_r1),
            _ => None,
        }
    }
}

impl From<io::Error> for FastqError {
    #[inline]
    fn from(err: io::Error) -> Self {
        FastqError::Io(err)
    }
}

struct RecordSuffix(u64);

impl fmt::Display for RecordSuffix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0 == 0 {
            Ok(())
        } else {
            write!(f, " (record {})", self.0)
        }
    }
}

struct Printable(u8);

impl fmt::Display for Printable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_ascii_graphic() || self.0 == b' ' {
            write!(f, "{:?}", self.0 as char)
        } else {
            f.write_str("non-printable")
        }
    }
}

impl fmt::Display for FastqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FastqError::Io(e) => write!(f, "I/O error: {e}"),
            FastqError::InvalidFormat {
                offset,
                record,
                kind,
            } => write!(
                f,
                "invalid FASTQ at offset {offset}{}: {kind}",
                RecordSuffix(*record)
            ),
            FastqError::UnexpectedEof { offset, record } => write!(
                f,
                "unexpected end of input at offset {offset}{}",
                RecordSuffix(*record)
            ),
            FastqError::LengthMismatch {
                offset,
                record,
                seq_len,
                qual_len,
            } => write!(
                f,
                "sequence/quality length mismatch at offset {offset}{}: seq_len={seq_len}, qual_len={qual_len}",
                RecordSuffix(*record)
            ),
            FastqError::InvalidBase {
                offset,
                record,
                byte,
            } => write!(
                f,
                "invalid base byte 0x{byte:02x} ({}) at offset {offset}{}",
                Printable(*byte),
                RecordSuffix(*record)
            ),
            FastqError::InvalidQuality {
                offset,
                record,
                byte,
            } => write!(
                f,
                "invalid quality byte 0x{byte:02x} ({}) at offset {offset}{}",
                Printable(*byte),
                RecordSuffix(*record)
            ),
            FastqError::PairedCountMismatch { which, record } => write!(
                f,
                "paired-end files hold different numbers of records: {which:?} ended first{}",
                RecordSuffix(*record)
            ),
            FastqError::PairedIdMismatch {
                offset_r1,
                offset_r2,
                record,
                id_r1,
                id_r2,
            } => write!(
                f,
                "paired-end read ID mismatch{}: R1 {:?} at offset {offset_r1}, R2 {:?} at offset {offset_r2}",
                RecordSuffix(*record),
                String::from_utf8_lossy(id_r1),
                String::from_utf8_lossy(id_r2)
            ),
            FastqError::Unsupported(op) => write!(f, "unsupported operation: {op}"),
        }
    }
}

impl std::error::Error for FastqError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            FastqError::Io(e) => Some(e),
            _ => None,
        }
    }
}
