mod common;

use kira_fastq::{
    Alphabet, FastqReader, FastqRecord, FastqWriter, ValidationMode, WriteValidation,
};

#[test]
fn plain_roundtrip_via_path() {
    let in_path = common::unique_path("in.fastq");
    common::write_plain(&in_path, b"@r1 lib=A\nACGT\n+\n!!!!\n@r2\nGGGG\n+\n@@@@\n");

    let out_path = common::unique_path("out.fastq");
    let mut writer = FastqWriter::from_path(&out_path).expect("open writer");
    let mut reader = FastqReader::from_path(&in_path).expect("open reader");
    while let Some(rec) = reader.next().expect("read") {
        writer.write_record(&rec).expect("write");
    }
    writer.flush().expect("flush");
    drop(writer);

    let written = std::fs::read(&out_path).expect("read written");
    assert_eq!(written, b"@r1 lib=A\nACGT\n+\n!!!!\n@r2\nGGGG\n+\n@@@@\n");
}

#[test]
fn gzip_roundtrip_via_path() {
    let in_path = common::unique_path("in.fastq");
    common::write_plain(&in_path, b"@a\nAC\n+\n!!\n@b\nGT\n+\n@@\n");

    let out_path = common::unique_path("out.fastq.gz");
    {
        let mut writer = FastqWriter::from_path(&out_path).expect("open gz writer");
        let mut reader = FastqReader::from_path(&in_path).expect("open reader");
        while let Some(rec) = reader.next().expect("read") {
            writer.write_record(&rec).expect("write");
        }
        writer.flush().expect("flush");
    }
    // Reopen via auto-detect — confirms valid gzip output.
    let mut reader = FastqReader::from_path_auto(&out_path).expect("reopen gz");
    let mut count = 0;
    while let Some(_rec) = reader.next().expect("read back") {
        count += 1;
    }
    assert_eq!(count, 2);
}

#[test]
fn write_parts_bypasses_record_type() {
    let out_path = common::unique_path("out_parts.fastq");
    let mut writer = FastqWriter::from_path(&out_path).expect("open");
    writer.write_parts(b"r0", b"AAAA", b"!!!!").expect("write");
    writer.flush().expect("flush");
    drop(writer);
    let data = std::fs::read(&out_path).expect("read");
    assert_eq!(data, b"@r0\nAAAA\n+\n!!!!\n");
}

#[test]
fn validation_catches_bad_base() {
    let out_path = common::unique_path("out_bad.fastq");
    let mut writer = FastqWriter::from_path(&out_path)
        .expect("open")
        .with_validation(WriteValidation::Bases)
        .with_alphabet(Alphabet::AcgtnStrict);
    let rec = FastqRecord::new(b"r0", b"ACZT", b"!!!!");
    let err = writer.write_record(&rec).expect_err("must reject Z");
    match err {
        kira_fastq::FastqError::InvalidBase { byte, .. } => assert_eq!(byte, b'Z'),
        other => panic!("expected InvalidBase, got {other:?}"),
    }
}

#[test]
fn validation_catches_bad_qual() {
    let out_path = common::unique_path("out_badq.fastq");
    let mut writer = FastqWriter::from_path(&out_path)
        .expect("open")
        .with_validation(WriteValidation::Qualities);
    // 0x00 is below Phred+33 minimum (33).
    let rec = FastqRecord::new(b"r0", b"AAAA", b"!!\x00!");
    let err = writer.write_record(&rec).expect_err("must reject low qual");
    match err {
        kira_fastq::FastqError::InvalidQuality { byte, .. } => assert_eq!(byte, 0x00),
        other => panic!("expected InvalidQuality, got {other:?}"),
    }
}

#[test]
fn write_to_vec_via_from_writer() {
    let mut buf: Vec<u8> = Vec::new();
    {
        let mut writer = FastqWriter::from_writer(&mut buf);
        writer
            .write_parts(b"r0", b"ACGT", b"!!!!")
            .expect("write parts");
        writer.flush().expect("flush");
    }
    assert_eq!(buf, b"@r0\nACGT\n+\n!!!!\n");
}

#[test]
fn read_from_in_memory_buffer() {
    let data: &[u8] = b"@r0\nACGT\n+\n!!!!\n@r1\nTT\n+\n##\n";
    let mut reader = FastqReader::from_reader(std::io::BufReader::new(data));
    let mut count = 0;
    while let Some(_rec) = reader.next().expect("read") {
        count += 1;
    }
    assert_eq!(count, 2);
}

#[test]
fn read_with_validation_from_buf_read() {
    let data: &[u8] = b"@r0\nACGT\n+\n!!!!\n";
    let mut reader = FastqReader::from_reader(std::io::BufReader::new(data))
        .with_validation(ValidationMode::BasesAndQualities)
        .with_alphabet(Alphabet::Iupac);
    let rec = reader.next().expect("read").expect("some");
    assert_eq!(rec.seq(), b"ACGT");
    assert!(reader.next().expect("eof").is_none());
}
