mod common;

use std::path::Path;

use kira_fastq::{FastqError, FastqReader, InvalidKind};

#[test]
fn bad_header() {
    let path = Path::new("tests/data/bad_header.fastq");
    let mut reader = FastqReader::from_path(path).expect("open");
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::InvalidFormat { offset, kind, .. } => {
            assert_eq!(offset, 0);
            assert_eq!(kind, InvalidKind::HeaderMissingAt);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn bad_plus() {
    let path = Path::new("tests/data/bad_plus.fastq");
    let mut reader = FastqReader::from_path(path).expect("open");
    let err = reader.next().expect_err("should error");
    let data = std::fs::read(path).expect("read");
    let plus_offset = data.iter().position(|&b| b == b'x').unwrap() as u64;
    match err {
        FastqError::InvalidFormat { offset, kind, .. } => {
            assert_eq!(offset, plus_offset);
            assert_eq!(kind, InvalidKind::PlusMissing);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn length_mismatch() {
    let path = Path::new("tests/data/bad_len.fastq");
    let mut reader = FastqReader::from_path(path).expect("open");
    let err = reader.next().expect_err("should error");
    let data = std::fs::read(path).expect("read");
    let mut newline_count = 0;
    let mut qual_offset = 0u64;
    for (i, &b) in data.iter().enumerate() {
        if b == b'\n' {
            newline_count += 1;
            if newline_count == 3 {
                qual_offset = (i + 1) as u64;
                break;
            }
        }
    }
    match err {
        FastqError::LengthMismatch {
            offset,
            seq_len,
            qual_len,
            ..
        } => {
            assert_eq!(offset, qual_offset);
            assert_eq!(seq_len, 4);
            assert_eq!(qual_len, 3);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn bad_header_gzip() {
    let path = common::unique_path("bad_header.fastq.gz");
    common::write_gzip(&path, b"r1\nACGT\n+\n!!!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::InvalidFormat { offset, kind, .. } => {
            assert_eq!(offset, 0);
            assert_eq!(kind, InvalidKind::HeaderMissingAt);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn bad_plus_gzip() {
    let path = common::unique_path("bad_plus.fastq.gz");
    common::write_gzip(&path, b"@r1\nACGT\nx\n!!!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::InvalidFormat { offset, kind, .. } => {
            assert_eq!(offset, 9);
            assert_eq!(kind, InvalidKind::PlusMissing);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn length_mismatch_gzip() {
    let path = common::unique_path("bad_len.fastq.gz");
    common::write_gzip(&path, b"@r1\nACGT\n+\n!!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::LengthMismatch {
            offset,
            seq_len,
            qual_len,
            ..
        } => {
            assert_eq!(offset, 11);
            assert_eq!(seq_len, 4);
            assert_eq!(qual_len, 3);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}
