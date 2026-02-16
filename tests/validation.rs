use std::io::Write;
use std::path::PathBuf;

use flate2::Compression;
use flate2::write::GzEncoder;
use kira_fastq::{FastqError, FastqReader, ValidationMode};

fn write_plain(path: &PathBuf, data: &[u8]) {
    std::fs::write(path, data).expect("write");
}

fn write_gzip(path: &PathBuf, data: &[u8]) {
    let file = std::fs::File::create(path).expect("create");
    let mut enc = GzEncoder::new(file, Compression::default());
    enc.write_all(data).expect("write");
    enc.finish().expect("finish");
}

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
    let dir = std::env::temp_dir();
    let path = dir.join("kira_fastq_invalid_default.fastq");
    write_plain(&path, data);
    let mut reader = FastqReader::from_path(&path).expect("open");
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGTX");
}

#[test]
fn bases_validation_plain() {
    let data = b"@r1\nACGTN\n+\n!!!!!\n";
    let dir = std::env::temp_dir();
    let path = dir.join("kira_fastq_valid_bases.fastq");
    write_plain(&path, data);
    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_validation(ValidationMode::Bases);
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
        let dir = std::env::temp_dir();
        let path = dir.join("kira_fastq_invalid_base.fastq");
        write_plain(&path, data);
        let mut reader = FastqReader::from_path(&path)
            .expect("open")
            .with_validation(ValidationMode::Bases);
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
    let dir = std::env::temp_dir();
    let path = dir.join("kira_fastq_invalid_qual.fastq");
    write_plain(&path, data);
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
    let dir = std::env::temp_dir();
    let path = dir.join("kira_fastq_valid_both.fastq.gz");
    write_gzip(&path, data);
    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_validation(ValidationMode::BasesAndQualities);
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGTN");
    assert_eq!(rec.qual(), b"!!!!!");
}

#[test]
fn bases_validation_gzip_invalid() {
    let data = b"@r1\nACGTX\n+\n!!!!!\n";
    let dir = std::env::temp_dir();
    let path = dir.join("kira_fastq_invalid_base.fastq.gz");
    write_gzip(&path, data);
    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_validation(ValidationMode::Bases);
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::InvalidBase { offset, byte } => {
            assert_eq!(byte, b'X');
            assert_eq!(offset, seq_start(data) + 4);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}
