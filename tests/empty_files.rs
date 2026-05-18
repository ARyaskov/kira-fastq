mod common;

use kira_fastq::FastqReader;

#[test]
fn empty_plain_file_is_clean_eof() {
    let path = common::unique_path("empty.fastq");
    common::write_plain(&path, b"");
    let mut reader = FastqReader::from_path(&path).expect("open empty");
    assert!(reader.next().expect("read").is_none());
}

#[test]
fn empty_gzip_member_is_clean_eof() {
    let path = common::unique_path("empty.fastq.gz");
    common::write_gzip(&path, b"");
    let mut reader = FastqReader::from_path(&path).expect("open empty gzip");
    assert!(reader.next().expect("read").is_none());
}
