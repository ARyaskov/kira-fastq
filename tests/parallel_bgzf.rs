//! Multi-threaded BGZF input must produce exactly what the sequential backend does.

mod common;

use kira_fastq::{FastqReader, ValidationMode};

fn corpus(records: usize) -> Vec<u8> {
    let mut seed = 0x2545_f491_4f6c_dd1du64;
    let mut rand = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    let mut out = Vec::new();
    for i in 0..records {
        out.extend_from_slice(format!("@read{i} len=120\n").as_bytes());
        let len = 40 + (rand() % 120) as usize;
        for _ in 0..len {
            out.push(b"ACGTN"[(rand() % 5) as usize]);
        }
        out.extend_from_slice(b"\n+\n");
        for _ in 0..len {
            out.push(33 + (rand() % 40) as u8);
        }
        out.push(b'\n');
    }
    out
}

fn collect(mut reader: FastqReader) -> Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> {
    let mut out = Vec::new();
    while let Some(rec) = reader.next().expect("read") {
        out.push((
            rec.header().to_vec(),
            rec.seq().to_vec(),
            rec.qual().to_vec(),
        ));
    }
    out
}

#[test]
fn matches_the_sequential_reader_across_many_blocks() {
    let data = corpus(4000);
    let path = common::unique_path("parallel.bgz");
    // Small blocks so records straddle block boundaries in every possible way.
    common::write_bgzf(&path, &data, 4096);

    let sequential = collect(FastqReader::from_bgzf_path(&path).expect("open"));
    assert!(sequential.len() == 4000);

    for threads in [1usize, 2, 4, 8] {
        let parallel =
            collect(FastqReader::from_bgzf_path_parallel(&path, threads).expect("open parallel"));
        assert_eq!(parallel, sequential, "threads={threads}");
    }
}

#[test]
fn thread_count_zero_picks_a_default() {
    let data = corpus(500);
    let path = common::unique_path("parallel_default.bgz");
    common::write_bgzf(&path, &data, 16384);
    let parallel = collect(FastqReader::from_bgzf_path_parallel(&path, 0).expect("open"));
    assert_eq!(parallel.len(), 500);
}

#[test]
fn validation_still_applies() {
    let path = common::unique_path("parallel_validate.bgz");
    common::write_bgzf(&path, b"@r1\nACGX\n+\n!!!!\n", 1000);
    let mut reader = FastqReader::from_bgzf_path_parallel(&path, 2)
        .expect("open")
        .with_validation(ValidationMode::Bases)
        .with_alphabet(kira_fastq::Alphabet::AcgtnStrict);
    let err = reader.next().expect_err("X is not a base");
    assert!(matches!(
        err,
        kira_fastq::FastqError::InvalidBase { byte: b'X', .. }
    ));
}

#[test]
fn corrupt_block_surfaces_as_an_error() {
    let path = common::unique_path("parallel_corrupt.bgz");
    common::write_bgzf(&path, &corpus(200), 4096);
    let mut bytes = std::fs::read(&path).expect("read");
    let len = bytes.len();
    bytes[len / 2] ^= 0xFF;
    std::fs::write(&path, bytes).expect("write");

    let mut reader = FastqReader::from_bgzf_path_parallel(&path, 4).expect("open");
    let mut failed = false;
    loop {
        match reader.next() {
            Ok(Some(_)) => {}
            Ok(None) => break,
            Err(_) => {
                failed = true;
                break;
            }
        }
    }
    assert!(failed, "a corrupted block must not read as clean data");
}

#[test]
fn dropping_early_does_not_hang() {
    let path = common::unique_path("parallel_drop.bgz");
    common::write_bgzf(&path, &corpus(5000), 4096);
    let mut reader = FastqReader::from_bgzf_path_parallel(&path, 4).expect("open");
    reader.next().expect("read").expect("first record");
    drop(reader);
}
