use kira_fastq::simd::newline::find_lf;

fn scalar_find_lf(buf: &[u8], start: usize) -> Option<usize> {
    let mut i = start;
    while i < buf.len() {
        if buf[i] == b'\n' {
            return Some(i);
        }
        i += 1;
    }
    None
}

fn make_buf(size: usize, newlines: &[usize]) -> Vec<u8> {
    let mut buf = vec![b'A'; size];
    for &pos in newlines {
        if pos < size {
            buf[pos] = b'\n';
        }
    }
    buf
}

fn lcg(seed: &mut u64) -> u8 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    (*seed >> 32) as u8
}

#[test]
fn newline_finds_match_scalar() {
    let sizes = [0usize, 1, 15, 16, 31, 32, 63, 64, 1024, 1024 * 1024];
    for &size in &sizes {
        let cases = [vec![], vec![0], vec![size.saturating_sub(1)], vec![1, 5, 9]];
        for newlines in cases {
            let buf = make_buf(size, &newlines);
            let starts = [0usize, 1, 7, 15, 31, 63, size.saturating_sub(1), size];
            for &start in &starts {
                let s = if start > size { size } else { start };
                let expected = scalar_find_lf(&buf, s);
                let got = find_lf(&buf, s);
                assert_eq!(expected, got);
            }
        }
    }
}

#[test]
fn newline_random_matches_scalar() {
    let mut seed = 0x1234_5678_9abc_def0u64;
    let mut buf = vec![0u8; 1 << 20];
    for b in buf.iter_mut() {
        *b = lcg(&mut seed);
    }
    for start in [
        0usize,
        1,
        7,
        15,
        31,
        63,
        127,
        255,
        1023,
        4096,
        buf.len() - 1,
    ] {
        let expected = scalar_find_lf(&buf, start);
        let got = find_lf(&buf, start);
        assert_eq!(expected, got);
    }
}
