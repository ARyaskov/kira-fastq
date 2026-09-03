//! Quality-byte range validation.

use crate::simd::cpu_features;
use crate::validation::QualityEncoding;

/// Validate against Phred+33, the encoding of every current sequencer.
///
/// Returns the index of the first out-of-range byte on failure.
#[inline]
pub fn validate_qual(buf: &[u8]) -> Result<(), usize> {
    validate_qual_encoding(buf, QualityEncoding::PHRED33)
}

/// Validate against an explicit encoding, e.g. [`QualityEncoding::PHRED64`] for Illumina 1.3
/// to 1.7 data.
#[inline]
pub fn validate_qual_encoding(buf: &[u8], encoding: QualityEncoding) -> Result<(), usize> {
    validate_qual_range(buf, encoding.min(), encoding.max())
}

/// Validate that every byte falls in `min ..= max`.
#[inline]
pub fn validate_qual_range(buf: &[u8], min: u8, max: u8) -> Result<(), usize> {
    if min > max {
        return Ok(());
    }

    #[cfg(target_arch = "x86_64")]
    {
        if cpu_features().avx2 {
            // SAFETY: AVX2 confirmed at runtime.
            return unsafe { validate_qual_avx2(buf, min, max) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if cpu_features().neon {
            // SAFETY: NEON confirmed at runtime.
            return unsafe { validate_qual_neon(buf, min, max) };
        }
    }

    let _ = cpu_features;
    validate_qual_scalar(buf, min, max)
}

#[inline]
fn validate_qual_scalar(buf: &[u8], min: u8, max: u8) -> Result<(), usize> {
    for (i, &b) in buf.iter().enumerate() {
        if b < min || b > max {
            return Err(i);
        }
    }
    Ok(())
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn validate_qual_avx2(buf: &[u8], min: u8, max: u8) -> Result<(), usize> {
    use std::arch::x86_64::*;
    let len = buf.len();
    let ptr = buf.as_ptr();
    // Unsigned range via saturating subtraction: bad iff (min - v) | (v - max) != 0.
    let vmin = _mm256_set1_epi8(min as i8);
    let vmax = _mm256_set1_epi8(max as i8);
    let zero = _mm256_setzero_si256();
    let mut i = 0usize;
    while i + 32 <= len {
        // SAFETY: 32 bytes remain at `i`.
        let v = unsafe { _mm256_loadu_si256(ptr.add(i) as *const __m256i) };
        let lt = _mm256_subs_epu8(vmin, v);
        let gt = _mm256_subs_epu8(v, vmax);
        let bad = _mm256_or_si256(lt, gt);
        let eqz = _mm256_cmpeq_epi8(bad, zero);
        let mask = _mm256_movemask_epi8(eqz) as u32 ^ 0xFFFF_FFFFu32;
        if mask != 0 {
            return Err(i + mask.trailing_zeros() as usize);
        }
        i += 32;
    }
    validate_qual_scalar(&buf[i..], min, max).map_err(|e| i + e)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn validate_qual_neon(buf: &[u8], min: u8, max: u8) -> Result<(), usize> {
    use std::arch::aarch64::*;
    let len = buf.len();
    let ptr = buf.as_ptr();
    let vmin = vdupq_n_u8(min);
    let vmax = vdupq_n_u8(max);
    let mut i = 0usize;
    while i + 16 <= len {
        // SAFETY: 16 bytes remain at `i`.
        let v = unsafe { vld1q_u8(ptr.add(i)) };
        let lt = vcltq_u8(v, vmin);
        let gt = vcgtq_u8(v, vmax);
        let bad = vorrq_u8(lt, gt);
        let narrow = vshrn_n_u16(vreinterpretq_u16_u8(bad), 4);
        let bits = vget_lane_u64(vreinterpret_u64_u8(narrow), 0);
        if bits != 0 {
            return Err(i + (bits.trailing_zeros() as usize >> 2));
        }
        i += 16;
    }
    validate_qual_scalar(&buf[i..], min, max).map_err(|e| i + e)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_scalar_for_every_byte_value() {
        for b in 0u8..=255 {
            let buf = [b; 48];
            assert_eq!(
                validate_qual(&buf),
                validate_qual_scalar(&buf, 33, 126),
                "byte 0x{b:02x}"
            );
        }
    }

    #[test]
    fn reports_first_bad_index_at_every_position() {
        for len in [1usize, 15, 31, 32, 33, 64, 129] {
            for bad_at in 0..len {
                let mut buf = vec![b'I'; len];
                buf[bad_at] = 0x01;
                assert_eq!(validate_qual(&buf), Err(bad_at));
            }
        }
    }

    #[test]
    fn phred64_rejects_sanger_bytes() {
        assert!(validate_qual_encoding(b"hhhh", QualityEncoding::PHRED64).is_ok());
        assert_eq!(
            validate_qual_encoding(b"hh!h", QualityEncoding::PHRED64),
            Err(2)
        );
    }
}
