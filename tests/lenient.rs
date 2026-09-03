//! Files that every mainstream FASTQ tool accepts and this crate used to reject.

mod common;

use kira_fastq::{FastqFormat, FastqReader};

/// Read the same bytes through every backend that can take them.
fn read_all_backends(data: &[u8]) -> Vec<(String, Vec<Vec<u8>>)> {
    let plain = common::unique_path("lenient.fastq");
    common::write_plain(&plain, data);
    let gz = common::unique_path("lenient.fastq.gz");
    common::write_gzip(&gz, data);
    let bgz = common::unique_path("lenient.fastq.bgz");
    common::write_bgzf(&bgz, data, 1000);

    vec![
        (
            "mmap".to_string(),
            drain(FastqReader::from_path(&plain).expect("open plain")),
        ),
        (
            "buffered".to_string(),
            drain(FastqReader::from_path_buffered(&plain).expect("open buffered")),
        ),
        (
            "memory".to_string(),
            drain(FastqReader::from_vec(data.to_vec())),
        ),
        (
            "stream".to_string(),
            drain(FastqReader::from_reader(std::io::BufReader::new(
                std::fs::File::open(&plain).expect("open"),
            ))),
        ),
        (
            "gzip".to_string(),
            drain(FastqReader::from_path(&gz).expect("open gzip")),
        ),
        (
            "bgzf".to_string(),
            drain(FastqReader::from_bgzf_path(&bgz).expect("open bgzf")),
        ),
    ]
}

fn drain(mut reader: FastqReader) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    while let Some(rec) = reader.next().expect("read") {
        out.push(rec.seq().to_vec());
    }
    out
}

#[test]
fn final_record_without_a_trailing_newline() {
    for (name, seqs) in read_all_backends(b"@a\nACGT\n+\n!!!!\n@b\nTT\n+\n##") {
        assert_eq!(seqs.len(), 2, "{name}");
        assert_eq!(seqs[1], b"TT", "{name}");
    }
}

/// cutadapt, fastp and Trimmomatic all emit zero-length reads; BWA, bowtie2 and STAR read them.
#[test]
fn zero_length_reads() {
    for (name, seqs) in read_all_backends(b"@a\n\n+\n\n@b\nGT\n+\n##\n") {
        assert_eq!(seqs.len(), 2, "{name}");
        assert!(seqs[0].is_empty(), "{name}");
    }
}

#[test]
fn blank_lines_between_and_after_records() {
    for (name, seqs) in read_all_backends(b"\n@a\nAC\n+\n!!\n\n\n@b\nGT\n+\n##\n\n") {
        assert_eq!(seqs.len(), 2, "{name}");
    }
}

#[test]
fn crlf_everywhere() {
    for (name, seqs) in read_all_backends(b"@a\r\nACGT\r\n+\r\n!!!!\r\n@b\r\nTT\r\n+\r\n##\r\n") {
        assert_eq!(seqs, vec![b"ACGT".to_vec(), b"TT".to_vec()], "{name}");
    }
}

#[test]
fn separator_line_may_repeat_the_id() {
    for (name, seqs) in read_all_backends(b"@a desc\nACGT\n+a desc\n!!!!\n") {
        assert_eq!(seqs.len(), 1, "{name}");
    }
}

#[test]
fn round_trips_through_the_writer() {
    // A file the reader accepts must be writable again, and reading that back must agree.
    let data = b"@a\n\n+\n\n@b\nGT\n+\n##";
    let in_path = common::unique_path("lenient_rt_in.fastq");
    common::write_plain(&in_path, data);
    let out_path = common::unique_path("lenient_rt_out.fastq");

    let mut reader = FastqReader::from_path(&in_path).expect("open");
    let mut writer = kira_fastq::FastqWriter::from_path(&out_path).expect("create");
    while let Some(rec) = reader.next().expect("read") {
        writer.write_record(&rec).expect("write");
    }
    writer.finish().expect("finish");

    assert_eq!(
        std::fs::read(&out_path).expect("read back"),
        b"@a\n\n+\n\n@b\nGT\n+\n##\n"
    );
}

#[test]
fn multi_line_records_stay_lenient() {
    let path = common::unique_path("lenient_multi.fastq");
    common::write_plain(&path, b"\n@a\nAC\nGT\n+\n!!\n!!\n\n@b\nTT\n+\n##");
    let mut reader = FastqReader::from_path(&path)
        .expect("open")
        .with_format(FastqFormat::MultiLine);
    let mut seqs = Vec::new();
    while let Some(rec) = reader.next().expect("read") {
        seqs.push(rec.seq().to_vec());
    }
    assert_eq!(seqs, vec![b"ACGT".to_vec(), b"TT".to_vec()]);
}
