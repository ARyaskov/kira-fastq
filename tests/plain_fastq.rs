use std::path::Path;

use kira_fastq::FastqReader;

#[test]
fn reads_plain_fastq() {
    let path = Path::new("tests/data/plain.fastq");
    let mut reader = FastqReader::from_path(path).expect("open");

    let rec1 = reader.next().expect("read").expect("rec1");
    assert_eq!(rec1.header(), b"r1");
    assert_eq!(rec1.seq(), b"ACGT");
    assert_eq!(rec1.qual(), b"!!!!");
    assert_eq!(rec1.len(), 4);

    let rec2 = reader.next().expect("read").expect("rec2");
    assert_eq!(rec2.header(), b"r2");
    assert_eq!(rec2.seq(), b"TT");
    assert_eq!(rec2.qual(), b"##");
    assert_eq!(rec2.len(), 2);

    let end = reader.next().expect("read");
    assert!(end.is_none());
}
