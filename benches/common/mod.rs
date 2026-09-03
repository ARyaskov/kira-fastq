//! Shared corpus generation for the benchmarks.
//!
//! The old benches ran over a 213 KB file, which fits in L2 and mostly measured the cost of
//! opening a reader. These generate a corpus large enough to leave cache, write it once per
//! machine into the target directory, and reuse it across runs.

#![allow(dead_code)]

use std::io::Write;
use std::path::{Path, PathBuf};

/// 64 MiB of reads: big enough to fall out of any current L3.
pub const DEFAULT_BYTES: usize = 64 * 1024 * 1024;

pub struct Corpus {
    pub plain: PathBuf,
    pub gzip: PathBuf,
    pub bgzf: PathBuf,
    pub bytes: usize,
}

/// Generate (or reuse) the benchmark corpus. `read_len` is the read length in bases.
pub fn corpus(read_len: usize, bytes: usize) -> Corpus {
    let dir = std::env::temp_dir().join("kira-fastq-bench");
    std::fs::create_dir_all(&dir).expect("create bench dir");
    let plain = dir.join(format!("reads_{read_len}_{bytes}.fastq"));
    let gzip = dir.join(format!("reads_{read_len}_{bytes}.fastq.gz"));
    let bgzf = dir.join(format!("reads_{read_len}_{bytes}.fastq.bgz"));

    if !plain.exists() || !gzip.exists() || !bgzf.exists() {
        let data = generate(read_len, bytes);
        std::fs::write(&plain, &data).expect("write plain");
        write_gzip(&gzip, &data);
        write_bgzf(&bgzf, &data, 65280);
    }
    let size = std::fs::metadata(&plain).expect("stat").len() as usize;
    Corpus {
        plain,
        gzip,
        bgzf,
        bytes: size,
    }
}

/// Illumina-shaped records: Casava 1.8 headers, random bases, realistic quality spread.
pub fn generate(read_len: usize, bytes: usize) -> Vec<u8> {
    let mut state = 0x9e37_79b9_7f4a_7c15u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    let mut out = Vec::with_capacity(bytes + 1024);
    let mut i = 0u64;
    while out.len() < bytes {
        out.extend_from_slice(
            format!(
                "@A00123:45:HXXXXDRXX:1:{}:{}:{} 1:N:0:ATCACG\n",
                1101 + i % 4,
                1000 + (i * 7) % 30000,
                1000 + (i * 13) % 30000
            )
            .as_bytes(),
        );
        for _ in 0..read_len {
            out.push(b"ACGT"[(next() % 4) as usize]);
        }
        out.extend_from_slice(b"\n+\n");
        for _ in 0..read_len {
            // Skewed towards high quality, like a real run.
            let r = next() % 100;
            out.push(if r < 80 { b'F' } else { 33 + (r % 40) as u8 });
        }
        out.push(b'\n');
        i += 1;
    }
    out
}

pub fn write_gzip(path: &Path, data: &[u8]) {
    let file = std::fs::File::create(path).expect("create gzip");
    let mut enc = flate2::write::GzEncoder::new(file, flate2::Compression::new(6));
    enc.write_all(data).expect("write");
    enc.finish().expect("finish");
}

pub fn write_bgzf(path: &Path, data: &[u8], block: usize) {
    let mut out = Vec::new();
    for chunk in data.chunks(block) {
        let mut enc = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::new(6));
        enc.write_all(chunk).expect("deflate");
        let payload = enc.finish().expect("finish");
        let mut hasher = crc32fast::Hasher::new();
        hasher.update(chunk);
        let mut header = vec![
            0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00, 0, 0,
        ];
        let bsize = (header.len() + payload.len() + 8 - 1) as u16;
        header[16] = (bsize & 0xFF) as u8;
        header[17] = (bsize >> 8) as u8;
        out.extend_from_slice(&header);
        out.extend_from_slice(&payload);
        out.extend_from_slice(&hasher.finalize().to_le_bytes());
        out.extend_from_slice(&(chunk.len() as u32).to_le_bytes());
    }
    out.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02,
        0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ]);
    std::fs::write(path, out).expect("write bgzf");
}
