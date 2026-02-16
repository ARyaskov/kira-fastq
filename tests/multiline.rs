use std::path::Path;

use kira_fastq::{FastqError, FastqFormat, FastqReader, PairedFastqReader, ValidationMode};

#[test]
fn multiline_plain_happy() {
    let path = Path::new("tests/data/multiline.fastq");
    let mut reader = FastqReader::from_path(path)
        .expect("open")
        .with_format(FastqFormat::MultiLine);

    let r1 = reader.next().expect("read").expect("rec1");
    assert_eq!(r1.header(), b"r1");
    assert_eq!(r1.seq(), b"ACGTN");
    assert_eq!(r1.qual(), b"!!!!!");

    let r2 = reader.next().expect("read").expect("rec2");
    assert_eq!(r2.header(), b"r2");
    assert_eq!(r2.seq(), b"TTGG");
    assert_eq!(r2.qual(), b"##!!");
}

#[test]
fn multiline_gzip_happy() {
    let path = Path::new("tests/data/multiline.fastq.gz");
    let mut reader = FastqReader::from_path(path)
        .expect("open")
        .with_format(FastqFormat::MultiLine);
    let r1 = reader.next().expect("read").expect("rec1");
    assert_eq!(r1.seq(), b"ACGTN");
}

#[test]
fn multiline_missing_plus() {
    let path = Path::new("tests/data/multiline_no_plus.fastq");
    let mut reader = FastqReader::from_path(path)
        .expect("open")
        .with_format(FastqFormat::MultiLine);
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::UnexpectedEof { .. } => {}
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn multiline_qual_short() {
    let path = Path::new("tests/data/multiline_qual_short.fastq");
    let mut reader = FastqReader::from_path(path)
        .expect("open")
        .with_format(FastqFormat::MultiLine);
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::UnexpectedEof { .. } => {}
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn multiline_qual_long() {
    let path = Path::new("tests/data/multiline_qual_long.fastq");
    let mut reader = FastqReader::from_path(path)
        .expect("open")
        .with_format(FastqFormat::MultiLine);
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGT");
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::InvalidFormat { .. } => {}
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn multiline_validation() {
    let dir = std::env::temp_dir();
    let path = dir.join("kira_multiline_invalid.fastq");
    std::fs::write(&path, b"@r1\nAC\nGTX\n+\n!!\n!!!\n").expect("write");

    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_format(FastqFormat::MultiLine)
        .with_validation(ValidationMode::Bases);
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::InvalidBase { byte, .. } => {
            assert_eq!(byte, b'X');
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn multiline_paired() {
    let dir = std::env::temp_dir();
    let r1 = dir.join("kira_multiline_r1.fastq");
    let r2 = dir.join("kira_multiline_r2.fastq");
    std::fs::write(&r1, b"@r1 1\nAC\nGT\n+\n!!!!\n!!\n").expect("write");
    std::fs::write(&r2, b"@r1 2\nAC\nGT\n+\n!!!!\n!!\n").expect("write");

    let mut pe = PairedFastqReader::from_paths(&r1, &r2)
        .expect("open")
        .with_format(FastqFormat::MultiLine)
        .with_id_check(true);
    let pair = pe.next().expect("read").expect("pair");
    assert_eq!(pair.0.seq(), b"ACGT");
    assert_eq!(pair.1.seq(), b"ACGT");
}
