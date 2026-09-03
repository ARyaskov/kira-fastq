//! Validation knobs shared by the readers and the writers.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ValidationMode {
    #[default]
    None,
    Bases,
    Qualities,
    BasesAndQualities,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Alphabet {
    /// `ACGTN` only.
    AcgtnStrict,
    /// `ACGTN` plus the lower-case forms.
    AcgtnCase,
    /// IUPAC ambiguity codes (case-insensitive) plus the `.` and `-` gap characters.
    #[default]
    Iupac,
}

/// Accepted range of quality bytes.
///
/// FASTQ quality lines are ASCII; which sub-range is legal depends on the encoding the
/// producer used. Sanger/Illumina 1.8+ files are Phred+33, older Illumina pipelines
/// (1.3 to 1.7) emitted Phred+64.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QualityEncoding {
    offset: u8,
    min: u8,
    max: u8,
}

impl Default for QualityEncoding {
    #[inline]
    fn default() -> Self {
        Self::PHRED33
    }
}

impl QualityEncoding {
    /// Sanger / Illumina 1.8+ / everything modern: `!` (33) through `~` (126).
    pub const PHRED33: Self = Self {
        offset: 33,
        min: 33,
        max: 126,
    };

    /// Illumina 1.3 to 1.7: `@` (64) through `~` (126).
    pub const PHRED64: Self = Self {
        offset: 64,
        min: 64,
        max: 126,
    };

    /// A custom inclusive byte range. `offset` is the Phred zero point used by
    /// [`QualityEncoding::phred_score`].
    #[inline]
    pub const fn custom(offset: u8, min: u8, max: u8) -> Self {
        Self { offset, min, max }
    }

    #[inline]
    pub const fn offset(&self) -> u8 {
        self.offset
    }

    #[inline]
    pub const fn min(&self) -> u8 {
        self.min
    }

    #[inline]
    pub const fn max(&self) -> u8 {
        self.max
    }

    #[inline]
    pub const fn contains(&self, byte: u8) -> bool {
        byte >= self.min && byte <= self.max
    }

    /// Phred score of a quality byte under this encoding. Saturates at 0 for bytes below the
    /// zero point.
    #[inline]
    pub const fn phred_score(&self, byte: u8) -> u8 {
        byte.saturating_sub(self.offset)
    }
}

/// Best guess at the quality encoding of a sample of quality bytes.
///
/// Returns [`QualityEncoding::PHRED33`] as soon as a byte below 64 is seen, and
/// [`QualityEncoding::PHRED64`] when every byte is at or above 64 and at least one is above 74
/// (the range where the two encodings cannot both be plausible). `None` means the sample is
/// consistent with both encodings, which is normal for very short or very high-quality samples.
pub fn guess_quality_encoding(qual: &[u8]) -> Option<QualityEncoding> {
    let mut min = u8::MAX;
    let mut max = 0u8;
    for &b in qual {
        min = min.min(b);
        max = max.max(b);
    }
    if qual.is_empty() {
        return None;
    }
    if min < 64 {
        return Some(QualityEncoding::PHRED33);
    }
    if max > 74 {
        return Some(QualityEncoding::PHRED64);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn guesses_phred33_from_low_bytes() {
        assert_eq!(
            guess_quality_encoding(b"!!!!IIII"),
            Some(QualityEncoding::PHRED33)
        );
    }

    #[test]
    fn guesses_phred64_from_high_bytes() {
        assert_eq!(
            guess_quality_encoding(b"hhhhhfff"),
            Some(QualityEncoding::PHRED64)
        );
    }

    #[test]
    fn ambiguous_sample_is_none() {
        assert_eq!(guess_quality_encoding(b"@ABCDEFG"), None);
        assert_eq!(guess_quality_encoding(b""), None);
    }

    #[test]
    fn phred_score_uses_offset() {
        assert_eq!(QualityEncoding::PHRED33.phred_score(b'!'), 0);
        assert_eq!(QualityEncoding::PHRED33.phred_score(b'I'), 40);
        assert_eq!(QualityEncoding::PHRED64.phred_score(b'h'), 40);
    }
}
