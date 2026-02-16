#[inline(always)]
pub fn find_lf(buf: &[u8], start: usize) -> Option<usize> {
    if start >= buf.len() {
        return None;
    }
    let tail = &buf[start..];

    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: runtime feature check guarantees AVX2 support.
            return unsafe { find_lf_avx2(tail).map(|idx| start + idx) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if std::arch::is_aarch64_feature_detected!("neon") {
            // SAFETY: runtime feature check guarantees NEON support.
            return unsafe { find_lf_neon(tail).map(|idx| start + idx) };
        }
    }

    memchr::memchr(b'\n', tail).map(|idx| start + idx)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
// SAFETY: caller must ensure AVX2 is enabled.
unsafe fn find_lf_avx2(buf: &[u8]) -> Option<usize> {
    use std::arch::x86_64::*;
    let mut i = 0usize;
    let len = buf.len();
    let needle = _mm256_set1_epi8(b'\n' as i8);
    let ptr = buf.as_ptr();
    while i + 32 <= len {
        let chunk = unsafe { _mm256_loadu_si256(ptr.add(i) as *const __m256i) };
        let cmp = _mm256_cmpeq_epi8(chunk, needle);
        let mask = _mm256_movemask_epi8(cmp) as u32;
        if mask != 0 {
            return Some(i + mask.trailing_zeros() as usize);
        }
        i += 32;
    }
    while i < len {
        if unsafe { *ptr.add(i) } == b'\n' {
            return Some(i);
        }
        i += 1;
    }
    None
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
// SAFETY: caller must ensure NEON is enabled.
unsafe fn find_lf_neon(buf: &[u8]) -> Option<usize> {
    use std::arch::aarch64::*;
    let mut i = 0usize;
    let len = buf.len();
    let needle = vdupq_n_u8(b'\n');
    let ptr = buf.as_ptr();
    while i + 16 <= len {
        let chunk = unsafe { vld1q_u8(ptr.add(i)) };
        let cmp = vceqq_u8(chunk, needle);
        let mut tmp = [0u8; 16];
        vst1q_u8(tmp.as_mut_ptr(), cmp);
        for j in 0..16 {
            if tmp[j] != 0 {
                return Some(i + j);
            }
        }
        i += 16;
    }
    while i < len {
        if unsafe { *ptr.add(i) } == b'\n' {
            return Some(i);
        }
        i += 1;
    }
    None
}
