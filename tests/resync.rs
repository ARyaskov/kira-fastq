mod common;

use kira_fastq::{FastqReader, VirtualOffset};

// `@` is a valid Phred+33 quality byte, so a quality line can legitimately start with it;
// naive `@`-at-line-start resync would misalign. Quartet-pattern resync handles this.
#[test]
fn resync_skips_at_in_quality_line() {
    let payload = b"@r1\nACGT\n+\n@AAA\n@r2\nTTTT\n+\nBBBB\n";
    let path = common::unique_path("resync_at_qual.fastq");
    common::write_plain(&path, payload);

    let mut reader = FastqReader::from_path(&path).expect("open");
    let _ = reader.next().expect("read").expect("rec1");
    reader.seek(VirtualOffset(5)).expect("seek");

    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.header(), b"r2");
}

#[test]
fn resync_at_start_of_file_works() {
    let payload = b"@r1\nACGT\n+\n!!!!\n";
    let path = common::unique_path("resync_start.fastq");
    common::write_plain(&path, payload);

    let mut reader = FastqReader::from_path(&path).expect("open");
    reader.seek(VirtualOffset(0)).expect("seek");
    let rec = reader.next().expect("read").expect("rec");
    assert_eq!(rec.header(), b"r1");
}
