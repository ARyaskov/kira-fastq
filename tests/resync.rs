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

/// `tell` must describe where the reader actually is, so a checkpoint taken right after a seek
/// is usable.
#[test]
fn tell_after_seek_reports_the_resynchronised_position() {
    let payload = b"@r1\nACGT\n+\n@AAA\n@r2\nTTTT\n+\nBBBB\n";
    let path = common::unique_path("tell_after_seek.fastq");
    common::write_plain(&path, payload);

    let mut reader = FastqReader::from_path(&path).expect("open");
    reader.seek(VirtualOffset(5)).expect("seek");
    let after_seek = reader.tell();
    // Offset 16 is the start of the second record; the '@' at offset 11 is a quality byte.
    assert_eq!(after_seek.0, 16);

    let mut again = FastqReader::from_path(&path).expect("reopen");
    again.seek(after_seek).expect("seek");
    assert_eq!(again.next().expect("read").expect("record").header(), b"r2");
}

#[test]
fn seek_past_the_end_yields_no_records() {
    let path = common::unique_path("seek_past_end.fastq");
    common::write_plain(&path, b"@r1\nACGT\n+\n!!!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    reader.seek(VirtualOffset(10_000)).expect("seek");
    assert!(reader.next().expect("read").is_none());
}
