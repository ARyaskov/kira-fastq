mod common;

use kira_fastq::FastqReader;

fn lcg(seed: &mut u64) -> u8 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    (*seed >> 32) as u8
}

fn make_records(count: usize, seed: &mut u64) -> Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> {
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let len = 20 + (lcg(seed) as usize % 80);
        let header = format!("r{}", i + 1).into_bytes();
        let mut seq = vec![0u8; len];
        let mut qual = vec![0u8; len];
        for j in 0..len {
            let r = lcg(seed) & 3;
            seq[j] = match r {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            qual[j] = 33 + (lcg(seed) % 40);
        }
        out.push((header, seq, qual));
    }
    out
}

fn write_fastq(path: &std::path::PathBuf, records: &[(Vec<u8>, Vec<u8>, Vec<u8>)]) {
    let mut out = Vec::new();
    for (h, s, q) in records {
        out.extend_from_slice(b"@");
        out.extend_from_slice(h);
        out.extend_from_slice(b"\n");
        out.extend_from_slice(s);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(q);
        out.extend_from_slice(b"\n");
    }
    std::fs::write(path, out).expect("write");
}

#[test]
fn plain_and_gzip_match() {
    let mut seed = 0x1357_9bdf_2468_ace0u64;
    let records = make_records(500, &mut seed);

    let plain_path = common::unique_path("fuzz_like_plain.fastq");
    let gzip_path = common::unique_path("fuzz_like_plain.fastq.gz");

    write_fastq(&plain_path, &records);
    let data = std::fs::read(&plain_path).expect("read");
    common::write_gzip(&gzip_path, &data);

    let mut plain = FastqReader::from_path(&plain_path).expect("open");
    let mut gzip = FastqReader::from_path(&gzip_path).expect("open");

    loop {
        let a = plain.next().expect("read");
        let b = gzip.next().expect("read");
        match (a, b) {
            (None, None) => break,
            (Some(x), Some(y)) => {
                assert_eq!(x.header(), y.header());
                assert_eq!(x.seq(), y.seq());
                assert_eq!(x.qual(), y.qual());
            }
            _ => panic!("length mismatch"),
        }
    }
}
