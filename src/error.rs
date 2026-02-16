use std::io;

#[derive(Debug)]
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
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnsupportedOperation {
    Seek,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
