mod common;

use kira_fastq::{Alphabet, FastqError, FastqReader, ValidationMode};

fn open_with(path: &std::path::PathBuf, mode: ValidationMode, alpha: Alphabet) -> FastqReader {
    FastqReader::from_path(path)
        .expect("open")
        .with_validation(mode)
        .with_alphabet(alpha)
}

#[test]
fn iupac_default_accepts_softmasked_lowercase() {
    let path = common::unique_path("softmask.fastq");
    common::write_plain(&path, b"@r1\nACGTacgtNn\n+\n!!!!!!!!!!\n");
    let mut reader = open_with(&path, ValidationMode::Bases, Alphabet::Iupac);
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGTacgtNn");
}

#[test]
fn iupac_accepts_ambiguity_codes() {
    let path = common::unique_path("iupac.fastq");
    common::write_plain(&path, b"@r1\nACGTRYSWKMBDHVN\n+\n!!!!!!!!!!!!!!!\n");
    let mut reader = open_with(&path, ValidationMode::Bases, Alphabet::Iupac);
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGTRYSWKMBDHVN");
}

#[test]
fn acgtn_strict_rejects_lowercase() {
    let path = common::unique_path("strict_lower.fastq");
    common::write_plain(&path, b"@r1\nACGTa\n+\n!!!!!\n");
    let mut reader = open_with(&path, ValidationMode::Bases, Alphabet::AcgtnStrict);
    let err = reader.next().expect_err("should reject lowercase");
    match err {
        FastqError::InvalidBase { byte, .. } => assert_eq!(byte, b'a'),
        other => panic!("unexpected: {other:?}"),
    }
}

#[test]
fn acgtn_case_accepts_lowercase() {
    let path = common::unique_path("case_lower.fastq");
    common::write_plain(&path, b"@r1\nACGTacgtn\n+\n!!!!!!!!!\n");
    let mut reader = open_with(&path, ValidationMode::Bases, Alphabet::AcgtnCase);
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGTacgtn");
}

#[test]
fn iupac_accepts_gap_chars() {
    let path = common::unique_path("gap.fastq");
    common::write_plain(&path, b"@r1\nAC-GT.N\n+\n!!!!!!!\n");
    let mut reader = open_with(&path, ValidationMode::Bases, Alphabet::Iupac);
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"AC-GT.N");
}
