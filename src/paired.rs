//! Paired-end reading: two mate files read in lock step, or one interleaved file.

use std::path::Path;

use crate::error::{FastqError, PairedWhich};
use crate::format::FastqFormat;
use crate::parser::ParsedRecord;
use crate::reader::FastqReader;
use crate::record::{FastqRecord, FastqRecordOwned};
use crate::validation::{Alphabet, QualityEncoding, ValidationMode};

/// Reads R1 and R2 from two files, yielding one pair per call.
pub struct PairedFastqReader {
    r1: FastqReader,
    r2: FastqReader,
    id_check: bool,
    pairs: u64,
}

impl PairedFastqReader {
    /// Open both mate files. Each is opened with [`FastqReader::from_path`], so the two may use
    /// different compression.
    pub fn from_paths<P: AsRef<Path>, Q: AsRef<Path>>(
        r1_path: P,
        r2_path: Q,
    ) -> Result<Self, FastqError> {
        Ok(Self::new(
            FastqReader::from_path(r1_path)?,
            FastqReader::from_path(r2_path)?,
        ))
    }

    /// Pair up two readers built elsewhere, e.g. over stdin and a file.
    pub fn new(r1: FastqReader, r2: FastqReader) -> Self {
        Self {
            r1,
            r2,
            id_check: false,
            pairs: 0,
        }
    }

    #[inline]
    pub fn with_validation(mut self, mode: ValidationMode) -> Self {
        self.r1 = self.r1.with_validation(mode);
        self.r2 = self.r2.with_validation(mode);
        self
    }

    #[inline]
    pub fn with_alphabet(mut self, alphabet: Alphabet) -> Self {
        self.r1 = self.r1.with_alphabet(alphabet);
        self.r2 = self.r2.with_alphabet(alphabet);
        self
    }

    #[inline]
    pub fn with_quality_encoding(mut self, encoding: QualityEncoding) -> Self {
        self.r1 = self.r1.with_quality_encoding(encoding);
        self.r2 = self.r2.with_quality_encoding(encoding);
        self
    }

    /// Compare read IDs of every pair, after [`canonical_read_id`] canonicalisation.
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

    /// Number of pairs returned so far.
    #[inline]
    pub fn pairs_read(&self) -> u64 {
        self.pairs
    }

    pub fn next(&mut self) -> Result<Option<(FastqRecord<'_>, FastqRecord<'_>)>, FastqError> {
        let Self {
            r1,
            r2,
            id_check,
            pairs,
        } = self;
        let record = *pairs + 1;
        let a = r1.next_parsed()?;
        let b = r2.next_parsed()?;

        match (a, b) {
            (None, None) => Ok(None),
            (Some(_), None) => Err(FastqError::PairedCountMismatch {
                which: PairedWhich::R2,
                record,
            }),
            (None, Some(_)) => Err(FastqError::PairedCountMismatch {
                which: PairedWhich::R1,
                record,
            }),
            (Some(a), Some(b)) => {
                if *id_check {
                    check_ids(&a, &b, record)?;
                }
                *pairs = record;
                Ok(Some((a.record, b.record)))
            }
        }
    }
}

/// Reads an interleaved file, where R1 and R2 alternate in a single stream.
///
/// This is the layout `samtools fastq` writes without `-1/-2`, and what many aligners accept on
/// stdin. The first mate of each pair is copied into the reader so both records can be handed
/// out together.
pub struct InterleavedFastqReader {
    inner: FastqReader,
    first: FastqRecordOwned,
    id_check: bool,
    pairs: u64,
}

impl InterleavedFastqReader {
    /// Wrap an existing reader.
    pub fn new(inner: FastqReader) -> Self {
        Self {
            inner,
            first: FastqRecordOwned::default(),
            id_check: false,
            pairs: 0,
        }
    }

    /// Open an interleaved file, choosing the backend by magic bytes.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self, FastqError> {
        Ok(Self::new(FastqReader::from_path(path)?))
    }

    #[inline]
    pub fn with_validation(mut self, mode: ValidationMode) -> Self {
        self.inner = self.inner.with_validation(mode);
        self
    }

    #[inline]
    pub fn with_alphabet(mut self, alphabet: Alphabet) -> Self {
        self.inner = self.inner.with_alphabet(alphabet);
        self
    }

    #[inline]
    pub fn with_quality_encoding(mut self, encoding: QualityEncoding) -> Self {
        self.inner = self.inner.with_quality_encoding(encoding);
        self
    }

    #[inline]
    pub fn with_format(mut self, format: FastqFormat) -> Self {
        self.inner = self.inner.with_format(format);
        self
    }

    #[inline]
    pub fn with_id_check(mut self, enabled: bool) -> Self {
        self.id_check = enabled;
        self
    }

    #[inline]
    pub fn pairs_read(&self) -> u64 {
        self.pairs
    }

    pub fn next(&mut self) -> Result<Option<(FastqRecord<'_>, FastqRecord<'_>)>, FastqError> {
        let record = self.pairs + 1;
        match self.inner.next()? {
            Some(rec) => self.first.copy_from(&rec),
            None => return Ok(None),
        }
        let Some(second) = self.inner.next()? else {
            return Err(FastqError::PairedCountMismatch {
                which: PairedWhich::R2,
                record,
            });
        };
        if self.id_check {
            let id1 = canonical_read_id(self.first.header());
            let id2 = canonical_read_id(second.header());
            if id1 != id2 {
                return Err(FastqError::PairedIdMismatch {
                    offset_r1: 0,
                    offset_r2: 0,
                    record,
                    id_r1: id1.into(),
                    id_r2: id2.into(),
                });
            }
        }
        self.pairs = record;
        Ok(Some((self.first.as_borrowed(), second)))
    }
}

#[inline]
fn check_ids(a: &ParsedRecord<'_>, b: &ParsedRecord<'_>, record: u64) -> Result<(), FastqError> {
    let id_r1 = canonical_read_id(a.record.header());
    let id_r2 = canonical_read_id(b.record.header());
    if id_r1 == id_r2 {
        return Ok(());
    }
    Err(FastqError::PairedIdMismatch {
        offset_r1: a.header_start,
        offset_r2: b.header_start,
        record,
        id_r1: id_r1.into(),
        id_r2: id_r2.into(),
    })
}

/// Strip the parts of a header that differ between mates, so R1 and R2 IDs compare equal.
///
/// Handles the three conventions in circulation:
///
/// * Casava 1.8+ (`@INST:... 1:N:0:INDEX`): everything from the first space or tab is dropped.
/// * Classic Illumina (`@READ/1`, `@READ/2`): the trailing `/N` is dropped.
/// * SRA with `fastq-dump -I` (`@SRR000001.1.1`, `@SRR000001.1.2`): a trailing `.1` or `.2` is
///   dropped. Other trailing `.N` are kept, since in plain `fastq-dump` output the number after
///   the dot is the spot ID and identifies the read rather than the mate.
pub fn canonical_read_id(header: &[u8]) -> &[u8] {
    let end = header
        .iter()
        .position(|&b| b == b' ' || b == b'\t')
        .unwrap_or(header.len());
    let prefix = &header[..end];
    if prefix.len() >= 3 {
        let last = prefix[prefix.len() - 1];
        let sep = prefix[prefix.len() - 2];
        if (sep == b'/' && last.is_ascii_digit()) || (sep == b'.' && (last == b'1' || last == b'2'))
        {
            return &prefix[..prefix.len() - 2];
        }
    }
    prefix
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonicalises_casava_and_classic_and_sra() {
        assert_eq!(
            canonical_read_id(b"M01234:23:000:1:1101:1:1 1:N:0:ATCACG"),
            b"M01234:23:000:1:1101:1:1"
        );
        assert_eq!(
            canonical_read_id(b"HWUSI-EAS100R:6:73:941:1973#0/1"),
            b"HWUSI-EAS100R:6:73:941:1973#0"
        );
        assert_eq!(canonical_read_id(b"SRR000001.1.2"), b"SRR000001.1");
        // A bare spot ID must survive: it names the read, not the mate.
        assert_eq!(canonical_read_id(b"SRR000001.17"), b"SRR000001.17");
        assert_eq!(canonical_read_id(b"plain"), b"plain");
    }

    #[test]
    fn mates_match_after_canonicalisation() {
        assert_eq!(canonical_read_id(b"read/1"), canonical_read_id(b"read/2"));
        assert_eq!(
            canonical_read_id(b"read 1:N:0:AA"),
            canonical_read_id(b"read 2:N:0:AA")
        );
    }
}
