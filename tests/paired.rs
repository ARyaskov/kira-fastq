use std::io::Write;
use std::path::PathBuf;

use flate2::Compression;
use flate2::write::GzEncoder;
use kira_fastq::{FastqError, PairedFastqReader, PairedWhich, ValidationMode};

fn write_plain(path: &PathBuf, data: &[u8]) {
    std::fs::write(path, data).expect("write");
}

fn write_gzip(path: &PathBuf, data: &[u8]) {
    let file = std::fs::File::create(path).expect("create");
    let mut enc = GzEncoder::new(file, Compression::default());
    enc.write_all(data).expect("write");
    enc.finish().expect("finish");
}

#[test]
fn paired_plain_happy() {
    let dir = std::env::temp_dir();
    let r1 = dir.join("kira_pe_r1.fastq");
    let r2 = dir.join("kira_pe_r2.fastq");
    let d1 = b"@r1 1\nACGT\n+\n!!!!\n@r2 1\nTT\n+\n##\n";
    let d2 = b"@r1 2\nACGT\n+\n!!!!\n@r2 2\nTT\n+\n##\n";
    write_plain(&r1, d1);
    write_plain(&r2, d2);

    let mut pe = PairedFastqReader::from_paths(&r1, &r2)
        .expect("open")
        .with_id_check(true);
    let p1 = pe.next().expect("read").expect("pair1");
    assert_eq!(p1.0.header(), b"r1 1");
    assert_eq!(p1.1.header(), b"r1 2");
    let p2 = pe.next().expect("read").expect("pair2");
    assert_eq!(p2.0.header(), b"r2 1");
    assert_eq!(p2.1.header(), b"r2 2");
    assert!(pe.next().expect("read").is_none());
}

#[test]
fn paired_gzip_happy() {
    let dir = std::env::temp_dir();
    let r1 = dir.join("kira_pe_r1.fastq.gz");
    let r2 = dir.join("kira_pe_r2.fastq.gz");
    let d1 = b"@r1\nACGT\n+\n!!!!\n";
    let d2 = b"@r1\nACGT\n+\n!!!!\n";
    write_gzip(&r1, d1);
    write_gzip(&r2, d2);

    let mut pe = PairedFastqReader::from_paths(&r1, &r2)
        .expect("open")
        .with_id_check(true);
    let p = pe.next().expect("read").expect("pair");
    assert_eq!(p.0.header(), b"r1");
    assert_eq!(p.1.header(), b"r1");
    assert!(pe.next().expect("read").is_none());
}

#[test]
fn paired_length_mismatch_r2_shorter() {
    let dir = std::env::temp_dir();
    let r1 = dir.join("kira_pe_len_r1.fastq");
    let r2 = dir.join("kira_pe_len_r2.fastq");
    let d1 = b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n";
    let d2 = b"@r1\nACGT\n+\n!!!!\n";
    write_plain(&r1, d1);
    write_plain(&r2, d2);

    let mut pe = PairedFastqReader::from_paths(&r1, &r2).expect("open");
    let _ = pe.next().expect("read").expect("pair1");
    let err = pe.next().expect_err("should error");
    match err {
        FastqError::PairedLengthMismatch { which } => {
            assert_eq!(which, PairedWhich::R2);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn paired_length_mismatch_r1_shorter() {
    let dir = std::env::temp_dir();
    let r1 = dir.join("kira_pe_len_r1b.fastq");
    let r2 = dir.join("kira_pe_len_r2b.fastq");
    let d1 = b"@r1\nACGT\n+\n!!!!\n";
    let d2 = b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n";
    write_plain(&r1, d1);
    write_plain(&r2, d2);

    let mut pe = PairedFastqReader::from_paths(&r1, &r2).expect("open");
    let _ = pe.next().expect("read").expect("pair1");
    let err = pe.next().expect_err("should error");
    match err {
        FastqError::PairedLengthMismatch { which } => {
            assert_eq!(which, PairedWhich::R1);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn paired_id_mismatch_opt_in() {
    let dir = std::env::temp_dir();
    let r1 = dir.join("kira_pe_id_r1.fastq");
    let r2 = dir.join("kira_pe_id_r2.fastq");
    let d1 = b"@r1\nACGT\n+\n!!!!\n";
    let d2 = b"@x1\nACGT\n+\n!!!!\n";
    write_plain(&r1, d1);
    write_plain(&r2, d2);

    let mut pe = PairedFastqReader::from_paths(&r1, &r2)
        .expect("open")
        .with_id_check(true);
    let err = pe.next().expect_err("should error");
    match err {
        FastqError::PairedIdMismatch { .. } => {}
        other => panic!("unexpected error: {other:?}"),
    }

    let mut pe2 = PairedFastqReader::from_paths(&r1, &r2).expect("open");
    let pair = pe2.next().expect("read").expect("pair");
    assert_eq!(pair.0.header(), b"r1");
    assert_eq!(pair.1.header(), b"x1");
}

#[test]
fn paired_validation_applies_to_both() {
    let dir = std::env::temp_dir();
    let r1 = dir.join("kira_pe_val_r1.fastq");
    let r2 = dir.join("kira_pe_val_r2.fastq");
    let d1 = b"@r1\nACGT\n+\n!!!!\n";
    let d2 = b"@r1\nACGX\n+\n!!!!\n";
    write_plain(&r1, d1);
    write_plain(&r2, d2);

    let mut pe = PairedFastqReader::from_paths(&r1, &r2)
        .expect("open")
        .with_validation(ValidationMode::BasesAndQualities);
    let err = pe.next().expect_err("should error");
    match err {
        FastqError::InvalidBase { byte, .. } => {
            assert_eq!(byte, b'X');
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn paired_mixed_plain_gzip() {
    let dir = std::env::temp_dir();
    let r1 = dir.join("kira_pe_mixed_r1.fastq");
    let r2 = dir.join("kira_pe_mixed_r2.fastq.gz");
    let d1 = b"@r1\nACGT\n+\n!!!!\n";
    let d2 = b"@r1\nACGT\n+\n!!!!\n";
    write_plain(&r1, d1);
    write_gzip(&r2, d2);

    let mut pe = PairedFastqReader::from_paths(&r1, &r2).expect("open");
    let pair = pe.next().expect("read").expect("pair");
    assert_eq!(pair.0.header(), b"r1");
    assert_eq!(pair.1.header(), b"r1");
    assert!(pe.next().expect("read").is_none());
}
