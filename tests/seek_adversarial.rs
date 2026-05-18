mod common;

use kira_fastq::{FastqError, FastqReader, VirtualOffset};

fn lcg(seed: &mut u64) -> u8 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    (*seed >> 32) as u8
}

fn write_bgzf_fixture() -> std::path::PathBuf {
    let path = common::unique_path("seek_fixture.bgz");
    let data = b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n";
    common::write_bgzf(&path, data, 20);
    path
}

#[test]
fn bgzf_sequential() {
    let path = write_bgzf_fixture();
    let mut reader = FastqReader::from_bgzf_path(&path).expect("open");
    let r1 = reader.next().expect("read").expect("rec1");
    assert_eq!(r1.header(), b"r1");
    let r2 = reader.next().expect("read").expect("rec2");
    assert_eq!(r2.header(), b"r2");
    assert!(reader.next().expect("read").is_none());
}

#[test]
fn bgzf_tell_seek_roundtrip() {
    let path = write_bgzf_fixture();
    let mut reader = FastqReader::from_bgzf_path(&path).expect("open");
    let r1 = reader.next().expect("read").expect("rec1");
    let _ = r1.header();
    let pos = reader.tell();
    let r2 = reader.next().expect("read").expect("rec2");
    let r2_header = r2.header().to_vec();
    let r2_seq = r2.seq().to_vec();
    reader.seek(pos).expect("seek");
    let r2b = reader.next().expect("read").expect("rec2b");
    assert_eq!(r2_header.as_slice(), r2b.header());
    assert_eq!(r2_seq.as_slice(), r2b.seq());
}

#[test]
fn bgzf_seek_middle_block() {
    let path = write_bgzf_fixture();
    let mut reader = FastqReader::from_bgzf_path(&path).expect("open");
    let pos = reader.tell();
    let mid = VirtualOffset(pos.0 + 5);
    reader.seek(mid).expect("seek");
    let rec = reader.next().expect("read").expect("rec");
    assert!(rec.header() == b"r1" || rec.header() == b"r2");
}

#[test]
fn bgzf_seek_invalid() {
    let path = write_bgzf_fixture();
    let mut reader = FastqReader::from_bgzf_path(&path).expect("open");
    let err = reader
        .seek(VirtualOffset(0x0001_0000))
        .expect_err("should error");
    match err {
        FastqError::InvalidFormat { .. } => {}
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn bgzf_truncated_block() {
    let base = b"@r1\nACGT\n+\n!!!!\n";
    let mut seed = 0x1234_5678_9abc_def0u64;
    let good = common::unique_path(&format!("trunc_good_{}.bgz", lcg(&mut seed)));
    common::write_bgzf(&good, base, 10);

    let data = std::fs::read(&good).expect("read");
    let bad = common::unique_path("trunc_bad.bgz");
    std::fs::write(&bad, &data[..data.len() / 2]).expect("write");
    let mut r = FastqReader::from_bgzf_path(&bad).expect("open");
    let _ = r.next();
}

#[test]
fn bgzf_bad_bsize() {
    let base = b"@r1\nACGT\n+\n!!!!\n";
    let mut seed = 0x1234_5678_9abc_def0u64;
    let good = common::unique_path(&format!("bsize_good_{}.bgz", lcg(&mut seed)));
    common::write_bgzf(&good, base, 10);

    let mut data = std::fs::read(&good).expect("read");
    data[16] = 0xFF;
    data[17] = 0xFF;
    let bad = common::unique_path("bsize_bad.bgz");
    std::fs::write(&bad, &data).expect("write");
    let mut r = FastqReader::from_bgzf_path(&bad).expect("open");
    let _ = r.next();
}

#[test]
fn bgzf_random_garbage_no_panic() {
    let path = common::unique_path("garbage.bgz");
    let mut seed = 0x9e37_79b9_7f4a_7c15u64;
    let mut data = vec![0u8; 1024];
    for b in data.iter_mut() {
        *b = lcg(&mut seed);
    }
    std::fs::write(&path, data).expect("write");
    let res = FastqReader::from_bgzf_path(&path);
    if let Ok(mut r) = res {
        let _ = r.next();
    }
}

#[test]
fn bgzf_error_variants_allowed() {
    let path = common::unique_path("bad4.bgz");
    std::fs::write(&path, b"bad").expect("write");
    let res = FastqReader::from_bgzf_path(&path);
    // Opening may fail OR opening succeeds then next() errors. Either is acceptable.
    if let Ok(mut r) = res {
        let _ = r.next();
    }
}

#[test]
fn seek_adversarial_plain() {
    let path = common::unique_path("seek_plain.fastq");
    common::write_plain(&path, b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    for off in [0u64, 1, 5, 10, 100] {
        let _ = reader.seek(VirtualOffset(off));
        let _ = reader.next();
    }
}

#[test]
fn seek_unsupported_gzip() {
    let path = common::unique_path("seek_gzip.fastq.gz");
    common::write_gzip(&path, b"@r1\nACGT\n+\n!!!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    let err = reader.seek(VirtualOffset(0)).expect_err("should error");
    match err {
        FastqError::Unsupported(_) => {}
        other => panic!("unexpected error: {other:?}"),
    }
}
