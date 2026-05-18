mod common;

use kira_fastq::{FastqError, FastqReader};

fn lcg(seed: &mut u64) -> u8 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    (*seed >> 32) as u8
}

fn gen_fastq(seed: &mut u64, records: usize, well_formed: bool) -> Vec<u8> {
    let mut out = Vec::new();
    for i in 0..records {
        let len = 1 + (lcg(seed) as usize % 50);
        out.extend_from_slice(b"@r");
        out.extend_from_slice(i.to_string().as_bytes());
        out.push(b'\n');
        for _ in 0..len {
            let r = lcg(seed) % 5;
            let b = match r {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            };
            out.push(b);
        }
        out.push(b'\n');
        out.push(b'+');
        out.push(b'\n');
        for _ in 0..len {
            out.push(33 + (lcg(seed) % 40));
        }
        if !well_formed && (lcg(seed) & 3) == 0 {
            out.pop();
        }
        if well_formed || (lcg(seed) & 7) != 0 {
            out.push(b'\n');
        }
        if !well_formed && (lcg(seed) & 7) == 0 {
            out.push(b'@');
        }
    }
    if !well_formed && (lcg(seed) & 1) == 0 {
        out.pop();
    }
    out
}

type RecordTuple = (Vec<u8>, Vec<u8>, Vec<u8>);

fn run_reader(path: &std::path::PathBuf) -> Result<Vec<RecordTuple>, FastqError> {
    let mut out = Vec::new();
    let mut reader = FastqReader::from_path(path).expect("open");
    while let Some(rec) = reader.next()? {
        out.push((
            rec.header().to_vec(),
            rec.seq().to_vec(),
            rec.qual().to_vec(),
        ));
        if out.len() > 10_000 {
            break;
        }
    }
    Ok(out)
}

#[test]
fn fuzz_fastq_no_panic() {
    let mut seed = 0x1234_5678_9abc_def0u64;
    for _ in 0..50 {
        let records = 1 + (lcg(&mut seed) as usize % 50);
        let data = gen_fastq(&mut seed, records, false);
        let plain = common::unique_path("fuzz_plain.fastq");
        let gzip = common::unique_path("fuzz_plain.fastq.gz");
        common::write_plain(&plain, &data);
        common::write_gzip(&gzip, &data);

        let _ = run_reader(&plain);
        let _ = run_reader(&gzip);
    }
}

#[test]
fn fuzz_plain_gzip_parity() {
    let mut seed = 0x0ddc_55aa_1122_3344u64;
    for _ in 0..20 {
        let records = 1 + (lcg(&mut seed) as usize % 50);
        let data = gen_fastq(&mut seed, records, true);
        let plain = common::unique_path("fuzz_parity.fastq");
        let gzip = common::unique_path("fuzz_parity.fastq.gz");
        common::write_plain(&plain, &data);
        common::write_gzip(&gzip, &data);
        let a = run_reader(&plain).expect("plain");
        let b = run_reader(&gzip).expect("gzip");
        assert_eq!(a, b);
    }
}

#[test]
fn fuzz_error_offsets_monotonic() {
    let mut seed = 0xabc0_1234_9876_5aa1u64;
    for _ in 0..20 {
        let records = 1 + (lcg(&mut seed) as usize % 50);
        let data = gen_fastq(&mut seed, records, false);
        let plain = common::unique_path("fuzz_offsets.fastq");
        common::write_plain(&plain, &data);
        let mut reader = FastqReader::from_path(&plain).expect("open");
        let mut last = 0u64;
        loop {
            match reader.next() {
                Ok(Some(_)) => {}
                Ok(None) => break,
                Err(err) => {
                    let off = match err {
                        FastqError::InvalidFormat { offset, .. } => offset,
                        FastqError::UnexpectedEof { offset } => offset,
                        FastqError::LengthMismatch { offset, .. } => offset,
                        FastqError::InvalidBase { offset, .. } => offset,
                        FastqError::InvalidQuality { offset, .. } => offset,
                        _ => 0,
                    };
                    assert!(off >= last);
                    assert!(off <= data.len() as u64);
                    break;
                }
            }
            last = reader.tell().0;
        }
    }
}
