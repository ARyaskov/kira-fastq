#[inline(always)]
pub fn validate_bases(buf: &[u8]) -> Result<(), usize> {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: runtime feature check guarantees AVX2 support.
            return unsafe { validate_bases_avx2(buf) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            // SAFETY: runtime feature check guarantees NEON support.
            return unsafe { validate_bases_neon(buf) };
        }
    }

    validate_bases_scalar(buf)
}

#[inline(always)]
fn validate_bases_scalar(buf: &[u8]) -> Result<(), usize> {
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

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// SAFETY: caller must ensure AVX2 is enabled.
// We only read within `buf` and use unaligned loads.
unsafe fn validate_bases_avx2(buf: &[u8]) -> Result<(), usize> {
    use std::arch::x86_64::*;
    let mut i = 0usize;
    let len = buf.len();
    let ptr = buf.as_ptr();
    let a = _mm256_set1_epi8(b'A' as i8);
    let c = _mm256_set1_epi8(b'C' as i8);
    let g = _mm256_set1_epi8(b'G' as i8);
    let t = _mm256_set1_epi8(b'T' as i8);
    let n = _mm256_set1_epi8(b'N' as i8);
    while i + 32 <= len {
        let v = unsafe { _mm256_loadu_si256(ptr.add(i) as *const __m256i) };
        let mut ok = _mm256_cmpeq_epi8(v, a);
        ok = _mm256_or_si256(ok, _mm256_cmpeq_epi8(v, c));
        ok = _mm256_or_si256(ok, _mm256_cmpeq_epi8(v, g));
        ok = _mm256_or_si256(ok, _mm256_cmpeq_epi8(v, t));
        ok = _mm256_or_si256(ok, _mm256_cmpeq_epi8(v, n));
        let bad = _mm256_movemask_epi8(ok) ^ 0xFFFF_FFFFu32 as i32;
        if bad != 0 {
            let idx = (bad as u32).trailing_zeros() as usize;
            return Err(i + idx);
        }
        i += 32;
    }
    while i < len {
        let b = unsafe { *ptr.add(i) };
        if b != b'A' && b != b'C' && b != b'G' && b != b'T' && b != b'N' {
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
unsafe fn validate_bases_neon(buf: &[u8]) -> Result<(), usize> {
    use std::arch::aarch64::*;
    let mut i = 0usize;
    let len = buf.len();
    let ptr = buf.as_ptr();
    let a = vdupq_n_u8(b'A');
    let c = vdupq_n_u8(b'C');
    let g = vdupq_n_u8(b'G');
    let t = vdupq_n_u8(b'T');
    let n = vdupq_n_u8(b'N');
    while i + 16 <= len {
        let v = unsafe { vld1q_u8(ptr.add(i)) };
        let mut ok = vceqq_u8(v, a);
        ok = vorrq_u8(ok, vceqq_u8(v, c));
        ok = vorrq_u8(ok, vceqq_u8(v, g));
        ok = vorrq_u8(ok, vceqq_u8(v, t));
        ok = vorrq_u8(ok, vceqq_u8(v, n));
        let bad = vmvnq_u8(ok);
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
        if b != b'A' && b != b'C' && b != b'G' && b != b'T' && b != b'N' {
            return Err(i);
        }
        i += 1;
    }
    Ok(())
}
