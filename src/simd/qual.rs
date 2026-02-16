#[inline(always)]
pub fn validate_qual(buf: &[u8]) -> Result<(), usize> {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: runtime feature check guarantees AVX2 support.
            return unsafe { validate_qual_avx2(buf) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            // SAFETY: runtime feature check guarantees NEON support.
            return unsafe { validate_qual_neon(buf) };
        }
    }

    validate_qual_scalar(buf)
}

#[inline(always)]
fn validate_qual_scalar(buf: &[u8]) -> Result<(), usize> {
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

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// SAFETY: caller must ensure AVX2 is enabled.
// We only read within `buf` and use unaligned loads.
unsafe fn validate_qual_avx2(buf: &[u8]) -> Result<(), usize> {
    use std::arch::x86_64::*;
    let mut i = 0usize;
    let len = buf.len();
    let ptr = buf.as_ptr();
    let min = _mm256_set1_epi8(33i8);
    let max = _mm256_set1_epi8(126i8);
    while i + 32 <= len {
        let v = unsafe { _mm256_loadu_si256(ptr.add(i) as *const __m256i) };
        let lt_min = _mm256_cmpgt_epi8(min, v);
        let gt_max = _mm256_cmpgt_epi8(v, max);
        let bad = _mm256_or_si256(lt_min, gt_max);
        let mask = _mm256_movemask_epi8(bad) as u32;
        if mask != 0 {
            return Err(i + mask.trailing_zeros() as usize);
        }
        i += 32;
    }
    while i < len {
        let b = unsafe { *ptr.add(i) };
        if b < 33 || b > 126 {
            return Err(i);
        }
        i += 1;
    }
    Ok(())
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
// SAFETY: caller must ensure NEON is enabled.
// We only read within `buf` and use unaligned loads.
unsafe fn validate_qual_neon(buf: &[u8]) -> Result<(), usize> {
    use std::arch::aarch64::*;
    let mut i = 0usize;
    let len = buf.len();
    let ptr = buf.as_ptr();
    let min = vdupq_n_u8(33);
    let max = vdupq_n_u8(126);
    while i + 16 <= len {
        let v = unsafe { vld1q_u8(ptr.add(i)) };
        let lt_min = vcltq_u8(v, min);
        let gt_max = vcgtq_u8(v, max);
        let bad = vorrq_u8(lt_min, gt_max);
        let mut tmp = [0u8; 16];
        vst1q_u8(tmp.as_mut_ptr(), bad);
        for j in 0..16 {
            if tmp[j] != 0 {
                return Err(i + j);
            }
        }
        i += 16;
    }
    while i < len {
        let b = unsafe { *ptr.add(i) };
        if b < 33 || b > 126 {
            return Err(i);
        }
        i += 1;
    }
    Ok(())
}
