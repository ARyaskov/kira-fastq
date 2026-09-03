mod common;

use kira_fastq::{FastqError, FastqReader, InvalidKind};

/// (header, sequence) of one record.
type HeaderAndSeq = (Vec<u8>, Vec<u8>);

/// Read every record of a gzip file through the public API.
fn read_all(path: &std::path::PathBuf) -> Result<Vec<HeaderAndSeq>, FastqError> {
    let mut reader = FastqReader::from_path(path)?;
    let mut out = Vec::new();
    while let Some(rec) = reader.next()? {
        out.push((rec.header().to_vec(), rec.seq().to_vec()));
    }
    Ok(out)
}

#[test]
fn gzip_roundtrip() {
    let path = common::unique_path("rt.fastq.gz");
    common::write_gzip(&path, b"@a\nACGT\n+\n!!!!\n@b\nTT\n+\n##\n");
    let recs = read_all(&path).expect("read");
    assert_eq!(recs.len(), 2);
    assert_eq!(recs[0].1, b"ACGT");
    assert_eq!(recs[1].0, b"b");
}

#[test]
fn gzip_strips_crlf() {
    let path = common::unique_path("crlf.fastq.gz");
    common::write_gzip(&path, b"@a\r\nACGT\r\n+\r\n!!!!\r\n");
    let recs = read_all(&path).expect("read");
    assert_eq!(recs[0].1, b"ACGT");
}

#[test]
fn gzip_accepts_missing_final_newline() {
    let path = common::unique_path("partial.fastq.gz");
    common::write_gzip(&path, b"@a\nACGT\n+\n!!!!");
    let recs = read_all(&path).expect("read");
    assert_eq!(recs.len(), 1);
}

/// bgzip, pigz and `cat a.gz b.gz` all produce concatenated members. Reading only the first one
/// silently drops most of the file.
#[test]
fn gzip_reads_every_member() {
    let path = common::unique_path("multi.fastq.gz");
    common::write_multi_member_gzip(
        &path,
        &[b"@a\nAC\n+\n!!\n", b"@b\nGT\n+\n##\n", b"@c\nTT\n+\n@@\n"],
    );
    let recs = read_all(&path).expect("read");
    assert_eq!(recs.len(), 3);
    assert_eq!(recs[2].0, b"c");
}

/// A deflate stream that stops mid-member must be an error, not a short read.
#[test]
fn gzip_truncated_member_is_an_error() {
    let path = common::unique_path("trunc.fastq.gz");
    let mut data = Vec::new();
    {
        use std::io::Write;
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(6));
        let mut payload = Vec::new();
        for i in 0..2000 {
            payload.extend_from_slice(format!("@r{i}\nACGTACGTAC\n+\n!!!!!!!!!!\n").as_bytes());
        }
        enc.write_all(&payload).expect("write");
        data.extend_from_slice(&enc.finish().expect("finish"));
    }
    std::fs::write(&path, &data[..data.len() / 2]).expect("write");

    let err = read_all(&path).expect_err("truncated input must fail");
    match err {
        FastqError::InvalidFormat { kind, .. } => {
            assert_eq!(kind, InvalidKind::GzipTruncated);
        }
        other => panic!("unexpected error: {other}"),
    }
}

/// CRC is verified on every member, with no feature flag to turn it off.
#[test]
fn gzip_crc_mismatch_is_an_error() {
    use std::io::Write;

    let path = common::unique_path("bad_crc.fastq.gz");
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(6));
    enc.write_all(b"@r1\nACGT\n+\n!!!!\n").expect("write");
    let mut gz = enc.finish().expect("finish");
    let len = gz.len();
    gz[len - 8] ^= 0xFF;
    std::fs::write(&path, gz).expect("write");

    let err = read_all(&path).expect_err("corrupt CRC must fail");
    match err {
        FastqError::InvalidFormat { kind, .. } => {
            assert_eq!(kind, InvalidKind::GzipTrailerCrc);
        }
        other => panic!("unexpected error: {other}"),
    }
}

#[test]
fn gzip_isize_mismatch_is_an_error() {
    use std::io::Write;

    let path = common::unique_path("bad_isize.fastq.gz");
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::new(6));
    enc.write_all(b"@r1\nACGT\n+\n!!!!\n").expect("write");
    let mut gz = enc.finish().expect("finish");
    let len = gz.len();
    gz[len - 4] ^= 0x0F;
    std::fs::write(&path, gz).expect("write");

    let err = read_all(&path).expect_err("corrupt ISIZE must fail");
    match err {
        FastqError::InvalidFormat { kind, .. } => {
            assert_eq!(kind, InvalidKind::GzipTrailerIsize);
        }
        other => panic!("unexpected error: {other}"),
    }
}
