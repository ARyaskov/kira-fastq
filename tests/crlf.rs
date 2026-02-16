use std::path::Path;

use kira_fastq::FastqReader;

#[test]
fn reads_crlf() {
    let path = Path::new("tests/data/crlf.fastq");
    let mut reader = FastqReader::from_path(path).expect("open");
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.header(), b"r1");
    assert_eq!(rec.seq(), b"ACGTA");
    assert_eq!(rec.qual(), b"!!!!!");
    assert_eq!(rec.len(), 5);
    assert!(reader.next().expect("read").is_none());
}
