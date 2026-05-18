mod common;

use kira_fastq::{Alphabet, FastqError, FastqReader, ValidationMode};

fn seq_start(data: &[u8]) -> u64 {
    let mut i = 0usize;
    while i < data.len() && data[i] != b'\n' {
        i += 1;
    }
    (i + 1) as u64
}

fn qual_start(data: &[u8]) -> u64 {
    let mut nl = 0usize;
    let mut i = 0usize;
    while i < data.len() {
        if data[i] == b'\n' {
            nl += 1;
            if nl == 3 {
                return (i + 1) as u64;
            }
        }
        i += 1;
    }
    data.len() as u64
}

#[test]
fn validation_disabled_accepts_invalid() {
    let data = b"@r1\nACGTX\n+\n!!!!?\n";
    let path = common::unique_path("val_disabled.fastq");
    common::write_plain(&path, data);
    let mut reader = FastqReader::from_path(&path).expect("open");
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGTX");
}

#[test]
fn bases_validation_plain() {
    let data = b"@r1\nACGTN\n+\n!!!!!\n";
    let path = common::unique_path("val_bases.fastq");
    common::write_plain(&path, data);
    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_validation(ValidationMode::Bases)
        .with_alphabet(Alphabet::AcgtnStrict);
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGTN");
}

#[test]
fn bases_validation_invalid_positions_plain() {
    let cases = [
        (b"@r1\nXCGTN\n+\n!!!!!\n".as_slice(), 0usize),
        (b"@r1\nACXTN\n+\n!!!!!\n".as_slice(), 2usize),
        (b"@r1\nACGTX\n+\n!!!!!\n".as_slice(), 4usize),
    ];
    for (data, bad_idx) in cases {
        let path = common::unique_path("val_bad_base.fastq");
        common::write_plain(&path, data);
        let mut reader = FastqReader::from_path(&path)
            .expect("open")
            .with_validation(ValidationMode::Bases)
            .with_alphabet(Alphabet::AcgtnStrict);
        let err = reader.next().expect_err("should error");
        match err {
            FastqError::InvalidBase { offset, byte } => {
                assert_eq!(byte, data[(seq_start(data) as usize) + bad_idx]);
                assert_eq!(offset, seq_start(data) + bad_idx as u64);
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }
}

#[test]
fn qualities_validation_invalid_plain() {
    let data = b"@r1\nACGT\n+\n!!!\x10\n";
    let path = common::unique_path("val_bad_qual.fastq");
    common::write_plain(&path, data);
    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_validation(ValidationMode::Qualities);
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::InvalidQuality { offset, byte } => {
            assert_eq!(byte, 0x10);
            assert_eq!(offset, qual_start(data) + 3);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn bases_and_qualities_gzip() {
    let data = b"@r1\nACGTN\n+\n!!!!!\n";
    let path = common::unique_path("val_both.fastq.gz");
    common::write_gzip(&path, data);
    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_validation(ValidationMode::BasesAndQualities)
        .with_alphabet(Alphabet::AcgtnStrict);
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGTN");
    assert_eq!(rec.qual(), b"!!!!!");
}

#[test]
fn bases_validation_gzip_invalid() {
    let data = b"@r1\nACGTX\n+\n!!!!!\n";
    let path = common::unique_path("val_bad_base.fastq.gz");
    common::write_gzip(&path, data);
    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_validation(ValidationMode::Bases)
        .with_alphabet(Alphabet::AcgtnStrict);
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::InvalidBase { offset, byte } => {
            assert_eq!(byte, b'X');
            assert_eq!(offset, seq_start(data) + 4);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn qualities_validation_accepts_full_phred_range() {
    // Bytes 33..=126 are all valid; check binary boundary.
    let mut data: Vec<u8> = Vec::from(b"@r1\nACGT\n+\n".as_slice());
    data.extend_from_slice(b"!~?@");
    data.push(b'\n');
    let path = common::unique_path("val_qual_full.fastq");
    common::write_plain(&path, &data);
    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_validation(ValidationMode::Qualities);
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.qual(), b"!~?@");
}
