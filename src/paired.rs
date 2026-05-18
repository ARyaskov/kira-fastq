use std::path::Path;

use crate::error::{FastqError, PairedWhich};
use crate::format::FastqFormat;
use crate::parser::ParsedRecord;
use crate::reader::FastqReader;
use crate::record::FastqRecord;
use crate::validation::{Alphabet, ValidationMode};

pub struct PairedFastqReader {
    r1: FastqReader,
    r2: FastqReader,
    id_check: bool,
}

impl PairedFastqReader {
    pub fn from_paths<P: AsRef<Path>, Q: AsRef<Path>>(
        r1_path: P,
        r2_path: Q,
    ) -> Result<Self, FastqError> {
        let r1 = FastqReader::from_path(r1_path)?;
        let r2 = FastqReader::from_path(r2_path)?;
        Ok(Self {
            r1,
            r2,
            id_check: false,
        })
    }

    #[inline]
    pub fn with_validation(self, mode: ValidationMode) -> Self {
        let Self { r1, r2, id_check } = self;
        Self {
            r1: r1.with_validation(mode),
            r2: r2.with_validation(mode),
            id_check,
        }
    }

    #[inline]
    pub fn with_alphabet(self, alphabet: Alphabet) -> Self {
        let Self { r1, r2, id_check } = self;
        Self {
            r1: r1.with_alphabet(alphabet),
            r2: r2.with_alphabet(alphabet),
            id_check,
        }
    }

    #[inline]
    pub fn with_id_check(mut self, enabled: bool) -> Self {
        self.id_check = enabled;
        self
    }

    #[inline]
    pub fn with_format(self, format: FastqFormat) -> Self {
        let Self { r1, r2, id_check } = self;
        Self {
            r1: r1.with_format(format),
            r2: r2.with_format(format),
            id_check,
        }
    }

    pub fn next(&mut self) -> Result<Option<(FastqRecord<'_>, FastqRecord<'_>)>, FastqError> {
        let Self { r1, r2, id_check } = self;
        let a = r1.next_parsed()?;
        let b = r2.next_parsed()?;

        match (a, b) {
            (None, None) => Ok(None),
            (Some(_), None) => Err(FastqError::PairedLengthMismatch {
                which: PairedWhich::R2,
            }),
            (None, Some(_)) => Err(FastqError::PairedLengthMismatch {
                which: PairedWhich::R1,
            }),
            (Some(a), Some(b)) => {
                if *id_check && !ids_match(&a, &b) {
                    return Err(FastqError::PairedIdMismatch {
                        offset_r1: a.header_start,
                        offset_r2: b.header_start,
                    });
                }
                Ok(Some((a.record, b.record)))
            }
        }
    }
}

#[inline]
fn ids_match(a: &ParsedRecord<'_>, b: &ParsedRecord<'_>) -> bool {
    canonical_read_id(a.record.header()) == canonical_read_id(b.record.header())
}

// Trim whitespace-suffix (Casava 1.8+ comments) then trailing `/N` (classic Illumina ≤1.7).
pub fn canonical_read_id(header: &[u8]) -> &[u8] {
    let mut end = header.len();
    for (i, &b) in header.iter().enumerate() {
        if b == b' ' || b == b'\t' {
            end = i;
            break;
        }
    }
    let prefix = &header[..end];
    if prefix.len() >= 2
        && prefix[prefix.len() - 2] == b'/'
        && prefix[prefix.len() - 1].is_ascii_digit()
    {
        return &prefix[..prefix.len() - 2];
    }
    prefix
}
