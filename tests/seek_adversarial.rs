use std::path::PathBuf;

use crc32fast::Hasher;
use flate2::Compression;
use flate2::write::DeflateEncoder;
use kira_fastq::{FastqError, FastqReader, VirtualOffset};
use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};

static COUNTER: AtomicU64 = AtomicU64::new(1);

fn write_bgzf(path: &PathBuf, data: &[u8], split: usize) {
    let mut out = Vec::new();
    let mut pos = 0usize;
    while pos < data.len() {
        let len = std::cmp::min(split, data.len() - pos);
        let chunk = &data[pos..pos + len];
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
        enc.write_all(chunk).expect("deflate");
        let deflated = enc.finish().expect("finish");

        let mut hasher = Hasher::new();
        hasher.update(chunk);
        let crc = hasher.finalize();
        let isize = chunk.len() as u32;

        let mut header = Vec::with_capacity(18);
        header.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff]);
        header.extend_from_slice(&[6, 0]);
        header.extend_from_slice(&[b'B', b'C', 2, 0, 0, 0]);

        let block_size = header.len() + deflated.len() + 8;
        let bsize = (block_size - 1) as u16;
        header[16] = (bsize & 0xFF) as u8;
        header[17] = (bsize >> 8) as u8;

        out.extend_from_slice(&header);
        out.extend_from_slice(&deflated);
        out.extend_from_slice(&crc.to_le_bytes());
        out.extend_from_slice(&isize.to_le_bytes());

        pos += len;
    }
    std::fs::write(path, out).expect("write");
}

fn lcg(seed: &mut u64) -> u8 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    (*seed >> 32) as u8
}

fn unique_path(name: &str) -> PathBuf {
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir();
    dir.join(format!("{}_{}", name, id))
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
    assert_eq!(rec.header(), b"r2");
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
    let dir = unique_path("kira_bgzf_adv_trunc");
    let base = b"@r1\nACGT\n+\n!!!!\n";
    let mut seed = 0x1234_5678_9abc_def0u64;
    std::fs::create_dir_all(&dir).expect("mkdir");
    let good = dir.join(format!("good_{}.bgz", lcg(&mut seed)));
    write_bgzf(&good, base, 10);

    let data = std::fs::read(&good).expect("read");
    let bad = dir.join("bad.bgz");
    std::fs::write(&bad, &data[..data.len() / 2]).expect("write");
    let mut r = FastqReader::from_bgzf_path(&bad).expect("open");
    let _ = r.next();
}

#[test]
fn bgzf_bad_bsize() {
    let dir = unique_path("kira_bgzf_adv_bsize");
    let base = b"@r1\nACGT\n+\n!!!!\n";
    let mut seed = 0x1234_5678_9abc_def0u64;
    std::fs::create_dir_all(&dir).expect("mkdir");
    let good = dir.join(format!("good_{}.bgz", lcg(&mut seed)));
    write_bgzf(&good, base, 10);

    let mut data = std::fs::read(&good).expect("read");
    data[16] = 0xFF;
    data[17] = 0xFF;
    let bad = dir.join("bad.bgz");
    std::fs::write(&bad, &data).expect("write");
    let mut r = FastqReader::from_bgzf_path(&bad).expect("open");
    let _ = r.next();
}

#[test]
fn bgzf_random_garbage_no_panic() {
    let path = unique_path("kira_bgzf_garbage.bgz");
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
    let path = unique_path("kira_bgzf_bad4.bgz");
    std::fs::write(&path, b"bad").expect("write");
    let mut r = FastqReader::from_bgzf_path(&path).expect("open");
    let err = r.next().expect_err("should error");
    match err {
        FastqError::InvalidFormat { .. } => {}
        FastqError::UnexpectedEof { .. } => {}
        _ => {}
    }
}

fn write_bgzf_fixture() -> PathBuf {
    let path = unique_path("kira_fastq_test.bgz");
    let data = b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n";
    write_bgzf(&path, data, 20);
    path
}

#[test]
fn seek_adversarial_plain() {
    let dir = std::env::temp_dir();
    let path = dir.join("kira_seek_plain.fastq");
    std::fs::write(&path, b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n").expect("write");
    let mut reader = FastqReader::from_path(&path).expect("open");
    for off in [0u64, 1, 5, 10, 100] {
        let _ = reader.seek(VirtualOffset(off));
        let _ = reader.next();
    }
}

#[test]
fn seek_unsupported_gzip() {
    let dir = std::env::temp_dir();
    let path = dir.join("kira_seek_gzip.fastq.gz");
    let data = b"@r1\nACGT\n+\n!!!!\n";
    let file = std::fs::File::create(&path).expect("create");
    let mut enc = flate2::write::GzEncoder::new(file, Compression::default());
    enc.write_all(data).expect("write");
    enc.finish().expect("finish");
    let mut reader = FastqReader::from_path(&path).expect("open");
    let err = reader.seek(VirtualOffset(0)).expect_err("should error");
    match err {
        FastqError::Unsupported(_) => {}
        other => panic!("unexpected error: {other:?}"),
    }
}
