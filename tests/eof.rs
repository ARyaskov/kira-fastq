use std::path::Path;

use kira_fastq::FastqError;
use kira_fastq::FastqReader;

#[test]
fn unexpected_eof() {
    let path = Path::new("tests/data/truncated.fastq");
    let mut reader = FastqReader::from_path(path).expect("open");
    let err = reader.next().expect_err("should error");
    let len = std::fs::read(path).expect("read").len() as u64;
    match err {
        FastqError::UnexpectedEof { offset } => {
            assert_eq!(offset, len);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn unexpected_eof_gzip() {
    let path = Path::new("tests/data/gzip_truncated.fastq.gz");
    let mut reader = FastqReader::from_path(path).expect("open");
    let err = reader.next().expect_err("should error");
    match err {
        FastqError::UnexpectedEof { offset } => {
            assert_eq!(offset, 11);
        }
        other => panic!("unexpected error: {other:?}"),
    }
}
