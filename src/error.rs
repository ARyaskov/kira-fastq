use std::fmt;
use std::io;

#[derive(Debug)]
#[non_exhaustive]
pub enum FastqError {
    Io(io::Error),
    InvalidFormat {
        offset: u64,
        kind: InvalidKind,
    },
    UnexpectedEof {
        offset: u64,
    },
    LengthMismatch {
        offset: u64,
        seq_len: usize,
        qual_len: usize,
    },
    InvalidBase {
        offset: u64,
        byte: u8,
    },
    InvalidQuality {
        offset: u64,
        byte: u8,
    },
    PairedLengthMismatch {
        which: PairedWhich,
    },
    PairedIdMismatch {
        offset_r1: u64,
        offset_r2: u64,
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
    SeqLineEmpty,
    QualLineEmpty,
    BgzfHeader,
    BgzfBlock,
    BgzfVirtualOffset,
    BgzfBlockTooLarge,
    BufferOverflow,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum UnsupportedOperation {
    Seek,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PairedWhich {
    R1,
    R2,
}

impl From<io::Error> for FastqError {
    #[inline]
    fn from(err: io::Error) -> Self {
        FastqError::Io(err)
    }
}

impl fmt::Display for FastqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FastqError::Io(e) => write!(f, "I/O error: {e}"),
            FastqError::InvalidFormat { offset, kind } => {
                write!(f, "invalid FASTQ format at offset {offset}: {kind:?}")
            }
            FastqError::UnexpectedEof { offset } => {
                write!(f, "unexpected end of input at offset {offset}")
            }
            FastqError::LengthMismatch {
                offset,
                seq_len,
                qual_len,
            } => write!(
                f,
                "sequence/quality length mismatch at offset {offset}: seq_len={seq_len}, qual_len={qual_len}"
            ),
            FastqError::InvalidBase { offset, byte } => write!(
                f,
                "invalid base byte 0x{byte:02x} ({:?}) at offset {offset}",
                *byte as char
            ),
            FastqError::InvalidQuality { offset, byte } => {
                write!(f, "invalid quality byte 0x{byte:02x} at offset {offset}")
            }
            FastqError::PairedLengthMismatch { which } => {
                write!(f, "paired-end length mismatch: {which:?} exhausted first")
            }
            FastqError::PairedIdMismatch {
                offset_r1,
                offset_r2,
            } => write!(
                f,
                "paired-end record ID mismatch (R1 offset {offset_r1}, R2 offset {offset_r2})"
            ),
            FastqError::Unsupported(op) => write!(f, "unsupported operation: {op:?}"),
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
