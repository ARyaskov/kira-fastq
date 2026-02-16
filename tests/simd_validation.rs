use kira_fastq::simd::bases::validate_bases;
use kira_fastq::simd::qual::validate_qual;

fn scalar_bases(buf: &[u8]) -> Result<(), usize> {
    let mut i = 0usize;
    while i < buf.len() {
        let b = buf[i];
        if b != b'A' && b != b'C' && b != b'G' && b != b'T' && b != b'N' {
            return Err(i);
        }
        i += 1;
    }
    Ok(())
}

fn scalar_qual(buf: &[u8]) -> Result<(), usize> {
    let mut i = 0usize;
    while i < buf.len() {
        let b = buf[i];
        if b < 33 || b > 126 {
            return Err(i);
        }
        i += 1;
    }
    Ok(())
}

fn lcg(seed: &mut u64) -> u8 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    (*seed >> 32) as u8
}

#[test]
fn bases_simd_matches_scalar() {
    let mut seed = 0xabc0_1234_9876_5aa1u64;
    let mut buf = vec![0u8; 1 << 20];
    for b in buf.iter_mut() {
        *b = match lcg(&mut seed) % 6 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            4 => b'N',
            _ => b'X',
        };
    }
    let expected = scalar_bases(&buf);
    let got = validate_bases(&buf);
    assert_eq!(expected, got);
}

#[test]
fn qual_simd_matches_scalar() {
    let mut seed = 0x0ddc_55aa_1122_3344u64;
    let mut buf = vec![0u8; 1 << 20];
    for b in buf.iter_mut() {
        let r = lcg(&mut seed);
        *b = if (r & 1) == 0 { 33 + (r % 94) } else { r };
    }
    let expected = scalar_qual(&buf);
    let got = validate_qual(&buf);
    assert_eq!(expected, got);
}
