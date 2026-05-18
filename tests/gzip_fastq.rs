mod common;

use kira_fastq::FastqReader;

#[test]
fn reads_gzip_single() {
    let path = common::unique_path("single.fastq.gz");
    common::write_gzip(&path, b"@r1\nACGT\n+\n!!!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.header(), b"r1");
    assert_eq!(rec.seq(), b"ACGT");
    assert_eq!(rec.qual(), b"!!!!");
    assert!(reader.next().expect("read").is_none());
}

#[test]
fn reads_gzip_multi() {
    let path = common::unique_path("multi.fastq.gz");
    common::write_gzip(&path, b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n");
    let mut reader = FastqReader::from_path(&path).expect("open");

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
    let path = common::unique_path("multi_member.fastq.gz");
    common::write_multi_member_gzip(&path, &[b"@r1\nACGT\n+\n!!!!\n", b"@r2\nTT\n+\n##\n"]);

    let mut reader = FastqReader::from_path(&path).expect("open");
    let rec1 = reader.next().expect("read").expect("rec1");
    assert_eq!(rec1.header(), b"r1");
    let rec2 = reader.next().expect("read").expect("rec2");
    assert_eq!(rec2.header(), b"r2");
    assert!(reader.next().expect("read").is_none());
}

#[test]
fn auto_detect_gzip() {
    let path = common::unique_path("auto.fastq");
    common::write_gzip(&path, b"@r1\nACGT\n+\n!!!!\n");
    let mut reader = FastqReader::from_path_auto(&path).expect("open");
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.seq(), b"ACGT");
}
