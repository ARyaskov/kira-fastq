mod common;

use kira_fastq::{FastqError, PairedFastqReader, PairedWhich, ValidationMode};

#[test]
fn paired_plain_happy() {
    let r1 = common::unique_path("pe_r1.fastq");
    let r2 = common::unique_path("pe_r2.fastq");
    let d1 = b"@r1 1\nACGT\n+\n!!!!\n@r2 1\nTT\n+\n##\n";
    let d2 = b"@r1 2\nACGT\n+\n!!!!\n@r2 2\nTT\n+\n##\n";
    common::write_plain(&r1, d1);
    common::write_plain(&r2, d2);

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
    let r1 = common::unique_path("pe_r1.fastq.gz");
    let r2 = common::unique_path("pe_r2.fastq.gz");
    common::write_gzip(&r1, b"@r1\nACGT\n+\n!!!!\n");
    common::write_gzip(&r2, b"@r1\nACGT\n+\n!!!!\n");

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
    let r1 = common::unique_path("pe_len_r1.fastq");
    let r2 = common::unique_path("pe_len_r2.fastq");
    common::write_plain(&r1, b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n");
    common::write_plain(&r2, b"@r1\nACGT\n+\n!!!!\n");

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
    let r1 = common::unique_path("pe_len_r1b.fastq");
    let r2 = common::unique_path("pe_len_r2b.fastq");
    common::write_plain(&r1, b"@r1\nACGT\n+\n!!!!\n");
    common::write_plain(&r2, b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n");

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
    let r1 = common::unique_path("pe_id_r1.fastq");
    let r2 = common::unique_path("pe_id_r2.fastq");
    common::write_plain(&r1, b"@r1\nACGT\n+\n!!!!\n");
    common::write_plain(&r2, b"@x1\nACGT\n+\n!!!!\n");

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
    let r1 = common::unique_path("pe_val_r1.fastq");
    let r2 = common::unique_path("pe_val_r2.fastq");
    common::write_plain(&r1, b"@r1\nACGT\n+\n!!!!\n");
    common::write_plain(&r2, b"@r1\nACGX\n+\n!!!!\n");

    let mut pe = PairedFastqReader::from_paths(&r1, &r2)
        .expect("open")
        .with_validation(ValidationMode::BasesAndQualities)
        .with_alphabet(kira_fastq::Alphabet::AcgtnStrict);
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
    let r1 = common::unique_path("pe_mixed_r1.fastq");
    let r2 = common::unique_path("pe_mixed_r2.fastq.gz");
    common::write_plain(&r1, b"@r1\nACGT\n+\n!!!!\n");
    common::write_gzip(&r2, b"@r1\nACGT\n+\n!!!!\n");

    let mut pe = PairedFastqReader::from_paths(&r1, &r2).expect("open");
    let pair = pe.next().expect("read").expect("pair");
    assert_eq!(pair.0.header(), b"r1");
    assert_eq!(pair.1.header(), b"r1");
    assert!(pe.next().expect("read").is_none());
}

#[test]
fn paired_classic_illumina_slash_pair_matches() {
    let r1 = common::unique_path("pe_illumina_r1.fastq");
    let r2 = common::unique_path("pe_illumina_r2.fastq");
    common::write_plain(&r1, b"@HWUSI-EAS100R:6:73:941:1973#0/1\nACGT\n+\n!!!!\n");
    common::write_plain(&r2, b"@HWUSI-EAS100R:6:73:941:1973#0/2\nACGT\n+\n!!!!\n");

    let mut pe = PairedFastqReader::from_paths(&r1, &r2)
        .expect("open")
        .with_id_check(true);
    let pair = pe.next().expect("read").expect("pair");
    assert_eq!(pair.0.header(), b"HWUSI-EAS100R:6:73:941:1973#0/1");
    assert_eq!(pair.1.header(), b"HWUSI-EAS100R:6:73:941:1973#0/2");
}

#[test]
fn paired_casava_18_pair_matches() {
    let r1 = common::unique_path("pe_casava_r1.fastq");
    let r2 = common::unique_path("pe_casava_r2.fastq");
    common::write_plain(
        &r1,
        b"@M01234:23:000000000-A1BCD:1:1101:12345:6789 1:N:0:NNNN\nACGT\n+\n!!!!\n",
    );
    common::write_plain(
        &r2,
        b"@M01234:23:000000000-A1BCD:1:1101:12345:6789 2:N:0:NNNN\nACGT\n+\n!!!!\n",
    );

    let mut pe = PairedFastqReader::from_paths(&r1, &r2)
        .expect("open")
        .with_id_check(true);
    let pair = pe.next().expect("read").expect("pair");
    assert!(pair.0.header().starts_with(b"M01234"));
    assert!(pair.1.header().starts_with(b"M01234"));
}
