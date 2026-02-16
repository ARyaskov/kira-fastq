use std::path::Path;

use crate::error::{FastqError, PairedWhich};
use crate::format::FastqFormat;
use crate::parser::ParsedRecord;
use crate::reader::FastqReader;
use crate::record::FastqRecord;
use crate::validation::ValidationMode;

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
    pub fn with_validation(mut self, mode: ValidationMode) -> Self {
        self.r1 = self.r1.with_validation(mode);
        self.r2 = self.r2.with_validation(mode);
        self
    }

    #[inline]
    pub fn with_id_check(mut self, enabled: bool) -> Self {
        self.id_check = enabled;
        self
    }

    #[inline]
    pub fn with_format(mut self, format: FastqFormat) -> Self {
        self.r1 = self.r1.with_format(format);
        self.r2 = self.r2.with_format(format);
        self
    }

    pub fn next(&mut self) -> Result<Option<(FastqRecord<'_>, FastqRecord<'_>)>, FastqError> {
        let r1_ptr = &mut self.r1 as *mut FastqReader;
        let r2_ptr = &mut self.r2 as *mut FastqReader;

        // SAFETY: r1 and r2 are distinct fields; we never alias mutable borrows.
        let a = unsafe { (*r1_ptr).next_parsed()? };
        // SAFETY: r1 and r2 are distinct fields; we never alias mutable borrows.
        let b = unsafe { (*r2_ptr).next_parsed()? };

        match (a, b) {
            (None, None) => Ok(None),
            (Some(_), None) => Err(FastqError::PairedLengthMismatch {
                which: PairedWhich::R2,
            }),
            (None, Some(_)) => Err(FastqError::PairedLengthMismatch {
                which: PairedWhich::R1,
            }),
            (Some(a), Some(b)) => {
                if self.id_check && !ids_match(&a, &b) {
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
    id_prefix(a.record.header()) == id_prefix(b.record.header())
}

#[inline]
fn id_prefix(header: &[u8]) -> &[u8] {
    let mut i = 0usize;
    while i < header.len() {
        let b = header[i];
        if b == b' ' || b == b'\t' {
            break;
        }
        i += 1;
    }
    &header[..i]
}
