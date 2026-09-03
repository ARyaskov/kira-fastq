//! Every backend must return the same records, and agree on which inputs are errors.
//!
//! The backends share a parser but not their line plumbing, so this is where a bug that only
//! shows up on, say, a record straddling a BGZF block boundary would surface.

mod common;

use kira_fastq::{FastqError, FastqFormat, FastqReader};

type RecordTuple = (Vec<u8>, Vec<u8>, Vec<u8>);

fn xorshift(state: &mut u64) -> u64 {
    *state ^= *state << 13;
    *state ^= *state >> 7;
    *state ^= *state << 17;
    *state
}

/// Random FASTQ with variable read lengths and, when `messy`, the deviations real files carry.
fn generate(seed: u64, records: usize, messy: bool) -> Vec<u8> {
    let mut state = seed | 1;
    let mut out = Vec::new();
    for i in 0..records {
        let len = (xorshift(&mut state) % 200) as usize;
        let len = if messy { len } else { len.max(1) };
        out.extend_from_slice(format!("@read{i} len={len}\n").as_bytes());
        for _ in 0..len {
            out.push(b"ACGTN"[(xorshift(&mut state) % 5) as usize]);
        }
        out.push(b'\n');
        out.push(b'+');
        if messy && xorshift(&mut state).is_multiple_of(4) {
            out.extend_from_slice(format!("read{i} len={len}").as_bytes());
        }
        out.push(b'\n');
        for _ in 0..len {
            out.push(33 + (xorshift(&mut state) % 60) as u8);
        }
        out.push(b'\n');
        if messy && xorshift(&mut state).is_multiple_of(8) {
            out.push(b'\n');
        }
    }
    if messy && seed.is_multiple_of(2) {
        // No trailing newline on the last record.
        out.pop();
    }
    out
}

fn drain(mut reader: FastqReader) -> Result<Vec<RecordTuple>, FastqError> {
    let mut out = Vec::new();
    while let Some(rec) = reader.next()? {
        out.push((
            rec.header().to_vec(),
            rec.seq().to_vec(),
            rec.qual().to_vec(),
        ));
    }
    Ok(out)
}

#[test]
fn all_backends_agree_on_generated_files() {
    for seed in 1..12u64 {
        for messy in [false, true] {
            let data = generate(seed, 60, messy);
            let plain = common::unique_path("agree.fastq");
            common::write_plain(&plain, &data);
            let gz = common::unique_path("agree.fastq.gz");
            common::write_gzip(&gz, &data);
            let multi_gz = common::unique_path("agree_multi.fastq.gz");
            let half = data.len() / 2;
            common::write_multi_member_gzip(&multi_gz, &[&data[..half], &data[half..]]);
            let bgz = common::unique_path("agree.fastq.bgz");
            // A block size that is not a multiple of any record length, so records straddle
            // block boundaries in every possible way.
            common::write_bgzf(&bgz, &data, 997);

            let expected = drain(FastqReader::from_path(&plain).expect("open plain"));
            let cases: Vec<(&str, Result<Vec<RecordTuple>, FastqError>)> = vec![
                (
                    "buffered",
                    drain(FastqReader::from_path_buffered(&plain).expect("open")),
                ),
                ("memory", drain(FastqReader::from_vec(data.clone()))),
                (
                    "stream",
                    drain(FastqReader::from_reader(std::io::BufReader::new(
                        std::fs::File::open(&plain).expect("open"),
                    ))),
                ),
                ("gzip", drain(FastqReader::from_path(&gz).expect("open"))),
                (
                    "gzip_multi_member",
                    drain(FastqReader::from_path(&multi_gz).expect("open")),
                ),
                (
                    "bgzf",
                    drain(FastqReader::from_bgzf_path(&bgz).expect("open")),
                ),
                (
                    "bgzf_parallel",
                    drain(FastqReader::from_bgzf_path_parallel(&bgz, 4).expect("open")),
                ),
            ];

            for (name, got) in cases {
                match (&expected, &got) {
                    (Ok(a), Ok(b)) => assert_eq!(a, b, "seed={seed} messy={messy} {name}"),
                    (Err(a), Err(b)) => assert_eq!(
                        std::mem::discriminant(a),
                        std::mem::discriminant(b),
                        "seed={seed} messy={messy} {name}: {a} vs {b}"
                    ),
                    (a, b) => panic!("seed={seed} messy={messy} {name}: {a:?} vs {b:?}"),
                }
            }
        }
    }
}

#[test]
fn multi_line_and_single_line_agree_on_single_line_files() {
    for seed in 1..8u64 {
        let data = generate(seed, 40, false);
        let path = common::unique_path("agree_ml.fastq");
        common::write_plain(&path, &data);

        let single = drain(FastqReader::from_path(&path).expect("open")).expect("single-line");
        let multi = drain(
            FastqReader::from_path(&path)
                .expect("open")
                .with_format(FastqFormat::MultiLine),
        )
        .expect("multi-line");
        assert_eq!(single, multi, "seed={seed}");
    }
}

/// Whatever the reader accepts, the writer must be able to write back, and the result must read
/// identically.
#[test]
fn round_trip_through_every_output_format() {
    for seed in 1..6u64 {
        let data = generate(seed, 50, true);
        let src = common::unique_path("rt_src.fastq");
        common::write_plain(&src, &data);
        let Ok(expected) = drain(FastqReader::from_path(&src).expect("open")) else {
            continue; // A generated file that does not parse is not interesting here.
        };

        for suffix in ["rt_out.fastq", "rt_out.fastq.gz", "rt_out.fastq.bgz"] {
            let out = common::unique_path(suffix);
            let mut writer = kira_fastq::FastqWriter::from_path(&out).expect("create");
            let mut reader = FastqReader::from_path(&src).expect("open");
            while let Some(rec) = reader.next().expect("read") {
                writer.write_record(&rec).expect("write");
            }
            writer.finish().expect("finish");

            let back = drain(FastqReader::from_path(&out).expect("reopen")).expect("read back");
            assert_eq!(expected, back, "seed={seed} {suffix}");
        }
    }
}
