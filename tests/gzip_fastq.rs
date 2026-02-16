use std::path::Path;

use flate2::Compression;
use flate2::write::GzEncoder;
use kira_fastq::FastqReader;
use std::io::Write;

#[test]
fn reads_gzip_single() {
    let path = Path::new("tests/data/gzip_single.fastq.gz");
    let mut reader = FastqReader::from_path(path).expect("open");
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.header(), b"r1");
    assert_eq!(rec.seq(), b"ACGT");
    assert_eq!(rec.qual(), b"!!!!");
    assert!(reader.next().expect("read").is_none());
}

#[test]
fn reads_gzip_multi() {
    let path = Path::new("tests/data/gzip_multi.fastq.gz");
    let mut reader = FastqReader::from_path(path).expect("open");
    let rec1 = reader.next().expect("read").expect("rec1");
    assert_eq!(rec1.header(), b"r1");
    assert_eq!(rec1.seq(), b"ACGT");
    assert_eq!(rec1.qual(), b"!!!!");

    let rec2 = reader.next().expect("read").expect("rec2");
    assert_eq!(rec2.header(), b"r2");
    assert_eq!(rec2.seq(), b"TT");
    assert_eq!(rec2.qual(), b"##");

    assert!(reader.next().expect("read").is_none());
}

#[test]
fn reads_gzip_multi_member() {
    let dir = std::env::temp_dir();
    let path = dir.join("kira_fastq_multi_member.fastq.gz");

    let mut part1 = Vec::new();
    part1.extend_from_slice(b"@r1\nACGT\n+\n!!!!\n");
    let mut part2 = Vec::new();
    part2.extend_from_slice(b"@r2\nTT\n+\n##\n");

    let mut out = Vec::new();
    {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(&part1).expect("write");
        out.extend_from_slice(&enc.finish().expect("finish"));
    }
    {
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(&part2).expect("write");
        out.extend_from_slice(&enc.finish().expect("finish"));
    }
    std::fs::write(&path, out).expect("write");

    let mut reader = FastqReader::from_path(&path).expect("open");
    let rec1 = reader.next().expect("read").expect("rec1");
    assert_eq!(rec1.header(), b"r1");
    let rec2 = reader.next().expect("read").expect("rec2");
    assert_eq!(rec2.header(), b"r2");
    assert!(reader.next().expect("read").is_none());
}
