//! Interleaved paired-end input, the layout `samtools fastq` writes by default.

mod common;

use kira_fastq::{FastqError, InterleavedFastqReader, ValidationMode};

#[test]
fn yields_pairs_in_order() {
    let path = common::unique_path("interleaved.fastq");
    common::write_plain(
        &path,
        b"@r1/1\nACGT\n+\n!!!!\n@r1/2\nTTTT\n+\n####\n@r2/1\nAA\n+\n!!\n@r2/2\nCC\n+\n##\n",
    );
    let mut reader = InterleavedFastqReader::from_path(&path)
        .expect("open")
        .with_id_check(true);

    let (a, b) = reader.next().expect("read").expect("pair");
    assert_eq!(a.seq(), b"ACGT");
    assert_eq!(b.seq(), b"TTTT");
    let (a, b) = reader.next().expect("read").expect("pair");
    assert_eq!(a.id(), b"r2/1");
    assert_eq!(b.seq(), b"CC");
    assert!(reader.next().expect("read").is_none());
    assert_eq!(reader.pairs_read(), 2);
}

#[test]
fn odd_record_count_is_an_error() {
    let path = common::unique_path("interleaved_odd.fastq");
    common::write_plain(
        &path,
        b"@r1/1\nAC\n+\n!!\n@r1/2\nGT\n+\n##\n@r2/1\nAA\n+\n!!\n",
    );
    let mut reader = InterleavedFastqReader::from_path(&path).expect("open");
    reader.next().expect("read").expect("first pair");
    assert!(matches!(
        reader.next(),
        Err(FastqError::PairedCountMismatch { .. })
    ));
}

#[test]
fn id_check_catches_a_shuffled_file() {
    let path = common::unique_path("interleaved_bad.fastq");
    common::write_plain(&path, b"@r1/1\nAC\n+\n!!\n@r9/2\nGT\n+\n##\n");
    let mut reader = InterleavedFastqReader::from_path(&path)
        .expect("open")
        .with_id_check(true);
    match reader.next() {
        Err(FastqError::PairedIdMismatch { id_r1, id_r2, .. }) => {
            assert_eq!(&*id_r1, b"r1");
            assert_eq!(&*id_r2, b"r9");
        }
        other => panic!("unexpected: {other:?}"),
    }
}

#[test]
fn validation_applies_to_both_mates() {
    let path = common::unique_path("interleaved_validate.fastq");
    common::write_plain(&path, b"@r1/1\nAC\n+\n!!\n@r1/2\nGX\n+\n##\n");
    let mut reader = InterleavedFastqReader::from_path(&path)
        .expect("open")
        .with_validation(ValidationMode::Bases)
        .with_alphabet(kira_fastq::Alphabet::AcgtnStrict);
    assert!(matches!(
        reader.next(),
        Err(FastqError::InvalidBase { byte: b'X', .. })
    ));
}

#[test]
fn works_over_a_compressed_stream() {
    let path = common::unique_path("interleaved.fastq.gz");
    common::write_gzip(&path, b"@r1/1\nAC\n+\n!!\n@r1/2\nGT\n+\n##\n");
    let mut reader = InterleavedFastqReader::from_path(&path)
        .expect("open")
        .with_id_check(true);
    let (a, b) = reader.next().expect("read").expect("pair");
    assert_eq!(a.seq(), b"AC");
    assert_eq!(b.seq(), b"GT");
}
