use crate::simd::cpu_features;

#[inline]
pub fn find_lf(buf: &[u8], start: usize) -> Option<usize> {
    if start >= buf.len() {
        return None;
    }
    let tail = &buf[start..];

    #[cfg(target_arch = "x86_64")]
    {
        let f = cpu_features();
        if f.avx512bw {
            // SAFETY: AVX-512BW confirmed at runtime.
            return unsafe { find_lf_avx512(tail) }.map(|idx| start + idx);
        }
        if f.avx2 {
            // SAFETY: AVX2 confirmed at runtime.
            return unsafe { find_lf_avx2(tail) }.map(|idx| start + idx);
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if cpu_features().neon {
            // SAFETY: NEON confirmed at runtime.
            return unsafe { find_lf_neon(tail) }.map(|idx| start + idx);
        }
    }

    let _ = cpu_features;
    memchr::memchr(b'\n', tail).map(|idx| start + idx)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512bw,avx512f")]
unsafe fn find_lf_avx512(buf: &[u8]) -> Option<usize> {
    use std::arch::x86_64::*;
    let len = buf.len();
    let ptr = buf.as_ptr();
    let needle = _mm512_set1_epi8(b'\n' as i8);
    let mut i = 0usize;
    while i + 64 <= len {
        let chunk = unsafe { _mm512_loadu_si512(ptr.add(i) as *const __m512i) };
        let mask = _mm512_cmpeq_epi8_mask(chunk, needle);
        if mask != 0 {
            return Some(i + mask.trailing_zeros() as usize);
        }
        i += 64;
    }
    memchr::memchr(b'\n', &buf[i..]).map(|idx| i + idx)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn find_lf_avx2(buf: &[u8]) -> Option<usize> {
    use std::arch::x86_64::*;
    let len = buf.len();
    let ptr = buf.as_ptr();
    let needle = _mm256_set1_epi8(b'\n' as i8);
    let mut i = 0usize;
    while i + 32 <= len {
        let chunk = unsafe { _mm256_loadu_si256(ptr.add(i) as *const __m256i) };
        let cmp = _mm256_cmpeq_epi8(chunk, needle);
        let mask = _mm256_movemask_epi8(cmp) as u32;
        if mask != 0 {
            return Some(i + mask.trailing_zeros() as usize);
        }
        i += 32;
    }
    memchr::memchr(b'\n', &buf[i..]).map(|idx| i + idx)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn find_lf_neon(buf: &[u8]) -> Option<usize> {
    use std::arch::aarch64::*;
    let len = buf.len();
    let ptr = buf.as_ptr();
    let needle = vdupq_n_u8(b'\n');
    let mut i = 0usize;
    while i + 16 <= len {
        let chunk = unsafe { vld1q_u8(ptr.add(i)) };
        let cmp = vceqq_u8(chunk, needle);
        // vshrn nibble-mask: each set lane → 4 bits, then bit-index / 4 = lane index.
        let narrow = vshrn_n_u16(vreinterpretq_u16_u8(cmp), 4);
        let bits = vget_lane_u64(vreinterpret_u64_u8(narrow), 0);
        if bits != 0 {
            return Some(i + (bits.trailing_zeros() as usize >> 2));
        }
        i += 16;
    }
    memchr::memchr(b'\n', &buf[i..]).map(|idx| i + idx)
}
