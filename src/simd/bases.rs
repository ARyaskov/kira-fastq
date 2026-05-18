use crate::simd::cpu_features;
use crate::validation::Alphabet;

type AlphabetLut = [u8; 256];

#[inline]
const fn build_lut(allowed: &[u8]) -> AlphabetLut {
    let mut lut = [0u8; 256];
    let mut i = 0;
    while i < allowed.len() {
        lut[allowed[i] as usize] = 1;
        i += 1;
    }
    lut
}

static LUT_ACGTN_STRICT: AlphabetLut = build_lut(b"ACGTN");
static LUT_ACGTN_CASE: AlphabetLut = build_lut(b"ACGTNacgtn");
static LUT_IUPAC: AlphabetLut = build_lut(b"ACGTURYSWKMBDHVNacgturyswkmbdhvn.-");

#[inline]
fn lut_for(alphabet: Alphabet) -> &'static AlphabetLut {
    match alphabet {
        Alphabet::AcgtnStrict => &LUT_ACGTN_STRICT,
        Alphabet::AcgtnCase => &LUT_ACGTN_CASE,
        Alphabet::Iupac => &LUT_IUPAC,
    }
}

#[inline]
pub fn validate_bases(buf: &[u8]) -> Result<(), usize> {
    validate_bases_with(buf, Alphabet::default())
}

#[inline]
pub fn validate_bases_with(buf: &[u8], alphabet: Alphabet) -> Result<(), usize> {
    let lut = lut_for(alphabet);

    #[cfg(target_arch = "x86_64")]
    {
        if cpu_features().avx2 {
            // SAFETY: AVX2 confirmed at runtime.
            return unsafe { validate_bases_avx2(buf, lut) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if cpu_features().neon {
            // SAFETY: NEON confirmed at runtime.
            return unsafe { validate_bases_neon(buf, lut) };
        }
    }

    let _ = cpu_features;
    validate_bases_scalar(buf, lut)
}

#[inline]
fn validate_bases_scalar(buf: &[u8], lut: &AlphabetLut) -> Result<(), usize> {
    for (i, &b) in buf.iter().enumerate() {
        if lut[b as usize] == 0 {
            return Err(i);
        }
    }
    Ok(())
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn validate_bases_avx2(buf: &[u8], lut: &AlphabetLut) -> Result<(), usize> {
    let len = buf.len();
    let mut i = 0usize;
    while i + 32 <= len {
        let chunk = &buf[i..i + 32];
        let mut bad: u32 = 0;
        for (j, &b) in chunk.iter().enumerate() {
            if lut[b as usize] == 0 {
                bad |= 1u32 << j;
            }
        }
        if bad != 0 {
            return Err(i + bad.trailing_zeros() as usize);
        }
        i += 32;
    }
    validate_bases_scalar(&buf[i..], lut).map_err(|e| i + e)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn validate_bases_neon(buf: &[u8], lut: &AlphabetLut) -> Result<(), usize> {
    use std::arch::aarch64::*;
    let len = buf.len();
    let mut i = 0usize;
    while i + 16 <= len {
        let chunk_ptr = unsafe { buf.as_ptr().add(i) };
        let mut tmp = [0u8; 16];
        for j in 0..16 {
            tmp[j] = lut[unsafe { *chunk_ptr.add(j) } as usize];
        }
        let v = unsafe { vld1q_u8(tmp.as_ptr()) };
        let zero = vdupq_n_u8(0);
        let cmp = vceqq_u8(v, zero);
        let narrow = vshrn_n_u16(vreinterpretq_u16_u8(cmp), 4);
        let bits = vget_lane_u64(vreinterpret_u64_u8(narrow), 0);
        if bits != 0 {
            return Err(i + (bits.trailing_zeros() as usize >> 2));
        }
        i += 16;
    }
    validate_bases_scalar(&buf[i..], lut).map_err(|e| i + e)
}
