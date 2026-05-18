#![allow(dead_code)]

use std::io::Write;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use crc32fast::Hasher;
use flate2::Compression;
use flate2::write::{DeflateEncoder, GzEncoder};

static COUNTER: AtomicU64 = AtomicU64::new(1);

pub fn unique_path(suffix: &str) -> PathBuf {
    let id = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    dir.join(format!("kira_fixture_{pid}_{id}_{suffix}"))
}

pub fn write_plain(path: &PathBuf, data: &[u8]) {
    std::fs::write(path, data).expect("write plain");
}

pub fn write_gzip(path: &PathBuf, data: &[u8]) {
    let file = std::fs::File::create(path).expect("create gzip");
    let mut enc = GzEncoder::new(file, Compression::default());
    enc.write_all(data).expect("gzip write");
    enc.finish().expect("gzip finish");
}

pub fn write_multi_member_gzip(path: &PathBuf, parts: &[&[u8]]) {
    let mut out = Vec::new();
    for p in parts {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(p).expect("write member");
        out.extend_from_slice(&enc.finish().expect("finish member"));
    }
    std::fs::write(path, out).expect("write multi-member gzip");
}

pub fn write_bgzf(path: &PathBuf, data: &[u8], split: usize) {
    let mut out = Vec::new();
    let mut pos = 0usize;
    let split = split.max(1);
    if data.is_empty() {
        out.extend_from_slice(&empty_bgzf_block());
    }
    while pos < data.len() {
        let len = std::cmp::min(split, data.len() - pos);
        let chunk = &data[pos..pos + len];

        let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
        enc.write_all(chunk).expect("deflate write");
        let deflated = enc.finish().expect("deflate finish");

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
    out.extend_from_slice(&empty_bgzf_block());
    std::fs::write(path, out).expect("write bgzf");
}

fn empty_bgzf_block() -> [u8; 28] {
    [
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ]
}
