//! The writer must not be able to emit a file that no reader can read back.

mod common;

use kira_fastq::{
    Alphabet, FastqError, FastqReader, FastqRecord, FastqWriter, InvalidKind, UnsupportedOperation,
    WriteValidation,
};

#[test]
fn rejects_length_mismatch_even_with_validation_off() {
    let path = common::unique_path("len_mismatch.fastq");
    let mut writer = FastqWriter::from_path(&path).expect("create");
    let err = writer
        .write_parts(b"r0", b"ACGT", b"!!!")
        .expect_err("seq and qual must agree");
    match err {
        FastqError::LengthMismatch {
            seq_len, qual_len, ..
        } => assert_eq!((seq_len, qual_len), (4, 3)),
        other => panic!("unexpected error: {other}"),
    }
    writer.finish().expect("finish");
    assert!(std::fs::read(&path).expect("read").is_empty());
}

#[test]
fn rejects_a_header_containing_a_line_break() {
    let path = common::unique_path("header_newline.fastq");
    let mut writer = FastqWriter::from_path(&path).expect("create");
    let err = writer
        .write_parts(b"r0\n@fake", b"AC", b"!!")
        .expect_err("header must be one line");
    assert!(matches!(
        err,
        FastqError::InvalidFormat {
            kind: InvalidKind::HeaderContainsNewline,
            ..
        }
    ));
}

#[test]
fn opt_in_check_rejects_line_breaks_in_sequence_and_quality() {
    let path = common::unique_path("seq_newline.fastq");
    let mut writer = FastqWriter::from_path(&path)
        .expect("create")
        .with_validation(WriteValidation::LineBreaks);
    assert!(matches!(
        writer.write_parts(b"r0", b"AC\nGT", b"!!!!!"),
        Err(FastqError::InvalidFormat {
            kind: InvalidKind::SeqContainsNewline,
            ..
        })
    ));
    assert!(matches!(
        writer.write_parts(b"r0", b"ACGTA", b"!!\n!!"),
        Err(FastqError::InvalidFormat {
            kind: InvalidKind::QualContainsNewline,
            ..
        })
    ));
}

#[test]
fn rejects_an_out_of_range_compression_level() {
    let path = common::unique_path("level.fastq.gz");
    match FastqWriter::to_gz_path(&path, 10) {
        Err(FastqError::Unsupported(UnsupportedOperation::CompressionLevel)) => {}
        Err(other) => panic!("unexpected error: {other}"),
        Ok(_) => panic!("level 10 is not a gzip level"),
    }
}

#[test]
fn still_catches_bad_bases_and_qualities() {
    let path = common::unique_path("content.fastq");
    let mut writer = FastqWriter::from_path(&path)
        .expect("create")
        .with_validation(WriteValidation::BasesAndQualities)
        .with_alphabet(Alphabet::AcgtnStrict);
    assert!(matches!(
        writer.write_record(&FastqRecord::new(b"r0", b"ACZT", b"!!!!")),
        Err(FastqError::InvalidBase { byte: b'Z', .. })
    ));
    assert!(matches!(
        writer.write_record(&FastqRecord::new(b"r0", b"ACGT", b"!!\x00!")),
        Err(FastqError::InvalidQuality { byte: 0x00, .. })
    ));
}

#[test]
fn phred64_output_can_be_validated() {
    let path = common::unique_path("phred64.fastq");
    let mut writer = FastqWriter::from_path(&path)
        .expect("create")
        .with_validation(WriteValidation::Qualities)
        .with_quality_encoding(kira_fastq::QualityEncoding::PHRED64);
    writer
        .write_parts(b"r0", b"ACGT", b"hhhh")
        .expect("Phred+64 bytes are in range");
    assert!(matches!(
        writer.write_parts(b"r0", b"ACGT", b"!!!!"),
        Err(FastqError::InvalidQuality { .. })
    ));
    writer.finish().expect("finish");
}

/// gzip and BGZF only become valid files once their trailer is written.
#[test]
fn finish_completes_compressed_output() {
    for name in ["finish.fastq.gz", "finish.fastq.bgz"] {
        let path = common::unique_path(name);
        let mut writer = FastqWriter::from_path(&path).expect("create");
        for i in 0..500 {
            writer
                .write_parts(format!("r{i}").as_bytes(), b"ACGTACGT", b"!!!!!!!!")
                .expect("write");
        }
        writer.finish().expect("finish");

        let mut reader = FastqReader::from_path(&path).expect("reopen");
        let mut n = 0;
        while reader.next().expect("read").is_some() {
            n += 1;
        }
        assert_eq!(n, 500, "{name}");
    }
}

/// BGZF output must be block-framed and end with the marker, or htslib calls it truncated.
#[test]
fn bgzf_output_is_seekable_and_terminated() {
    let path = common::unique_path("native.bgz");
    let mut writer = FastqWriter::to_bgzf_path(&path, 6).expect("create");
    for i in 0..20_000 {
        writer
            .write_parts(format!("r{i}").as_bytes(), b"ACGTACGTAC", b"!!!!!!!!!!")
            .expect("write");
    }
    writer.finish().expect("finish");

    let bytes = std::fs::read(&path).expect("read");
    let blocks = bytes
        .windows(4)
        .filter(|w| w == b"\x1f\x8b\x08\x04")
        .count();
    assert!(blocks > 3, "expected several blocks, found {blocks}");
    const BGZF_EOF: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
    assert!(
        bytes.ends_with(&BGZF_EOF),
        "file must end with the BGZF end-of-file marker"
    );

    let mut reader = FastqReader::from_bgzf_path(&path).expect("reopen");
    for _ in 0..10_000 {
        reader.next().expect("read").expect("record");
    }
    let checkpoint = reader.tell();
    let header = reader
        .next()
        .expect("read")
        .expect("record")
        .header()
        .to_vec();
    let mut resumed = FastqReader::from_bgzf_path(&path).expect("reopen");
    resumed.seek(checkpoint).expect("seek");
    assert_eq!(
        resumed.next().expect("read").expect("record").header(),
        header.as_slice()
    );
}

#[test]
fn parallel_bgzf_output_matches_single_threaded() {
    let sequential = common::unique_path("seq.bgz");
    let parallel = common::unique_path("par.bgz");
    for (path, threads) in [(&sequential, None), (&parallel, Some(4))] {
        let mut writer = match threads {
            None => FastqWriter::to_bgzf_path(path, 6).expect("create"),
            Some(t) => FastqWriter::to_bgzf_path_parallel(path, 6, t).expect("create"),
        };
        for i in 0..30_000 {
            writer
                .write_parts(format!("r{i}").as_bytes(), b"ACGTACGTAC", b"!!!!!!!!!!")
                .expect("write");
        }
        writer.finish().expect("finish");
    }
    assert_eq!(
        std::fs::read(&sequential).expect("read"),
        std::fs::read(&parallel).expect("read"),
        "parallel compression must not change the bytes"
    );
}

#[test]
fn writing_to_a_bgzf_path_needs_no_optional_feature() {
    let path = common::unique_path("default_feature.bgzf");
    let mut writer = FastqWriter::from_path(&path).expect("BGZF output in default features");
    writer.write_parts(b"r0", b"AC", b"!!").expect("write");
    writer.finish().expect("finish");
    let mut reader = FastqReader::from_path(&path).expect("reopen");
    assert_eq!(reader.next().expect("read").expect("record").seq(), b"AC");
}

#[test]
fn record_assembly_reserves_exactly() {
    // A wrong reserve means a reallocation on every record that sets a new size record.
    let mut buf = Vec::new();
    kira_fastq::writer::assemble_record(&mut buf, b"header", b"ACGT", b"!!!!");
    assert_eq!(buf, b"@header\nACGT\n+\n!!!!\n");
    assert_eq!(
        buf.len(),
        buf.capacity(),
        "reserve must match the written size"
    );
}
