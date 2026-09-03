//! Base-alphabet validation.
//!
//! The vector kernels test set membership with the nibble-table trick: two 16-byte tables are
//! indexed by the low and high nibble of each byte with one shuffle each, and a byte is valid
//! when the two lookups share a bit. That keeps the whole check in registers, unlike a 256-entry
//! byte LUT which needs one dependent load per base. Bytes >= 0x80 index past the high table and
//! read zero, so non-ASCII input is rejected without a special case.

use crate::simd::cpu_features;
use crate::validation::Alphabet;

type AlphabetLut = [u8; 256];
/// `(low-nibble table, high-nibble table)`, both 16 bytes wide.
type NibbleTables = ([u8; 16], [u8; 16]);

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

#[inline]
const fn build_nibble_tables(allowed: &[u8]) -> NibbleTables {
    let mut lo = [0u8; 16];
    let mut hi = [0u8; 16];
    let mut i = 0;
    while i < allowed.len() {
        let b = allowed[i];
        // The bit index is the high nibble, so the scheme covers ASCII only. Every alphabet
        // this crate ships is ASCII; a non-ASCII entry is a compile-time error.
        assert!(b < 0x80, "alphabet must be ASCII");
        let low = (b & 0x0F) as usize;
        let high = (b >> 4) as usize;
        lo[low] |= 1 << high;
        hi[high] |= 1 << high;
        i += 1;
    }
    (lo, hi)
}

const ACGTN_STRICT: &[u8] = b"ACGTN";
const ACGTN_CASE: &[u8] = b"ACGTNacgtn";
const IUPAC: &[u8] = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn.-";

static LUT_ACGTN_STRICT: AlphabetLut = build_lut(ACGTN_STRICT);
static LUT_ACGTN_CASE: AlphabetLut = build_lut(ACGTN_CASE);
static LUT_IUPAC: AlphabetLut = build_lut(IUPAC);

static NIB_ACGTN_STRICT: NibbleTables = build_nibble_tables(ACGTN_STRICT);
static NIB_ACGTN_CASE: NibbleTables = build_nibble_tables(ACGTN_CASE);
static NIB_IUPAC: NibbleTables = build_nibble_tables(IUPAC);

#[inline]
fn lut_for(alphabet: Alphabet) -> &'static AlphabetLut {
    match alphabet {
        Alphabet::AcgtnStrict => &LUT_ACGTN_STRICT,
        Alphabet::AcgtnCase => &LUT_ACGTN_CASE,
        Alphabet::Iupac => &LUT_IUPAC,
    }
}

#[inline]
fn tables_for(alphabet: Alphabet) -> &'static NibbleTables {
    match alphabet {
        Alphabet::AcgtnStrict => &NIB_ACGTN_STRICT,
        Alphabet::AcgtnCase => &NIB_ACGTN_CASE,
        Alphabet::Iupac => &NIB_IUPAC,
    }
}

/// Validate against the default alphabet ([`Alphabet::Iupac`]).
///
/// On success returns `Ok(())`; on failure returns the index of the first invalid byte.
#[inline]
pub fn validate_bases(buf: &[u8]) -> Result<(), usize> {
    validate_bases_with(buf, Alphabet::default())
}

/// Validate against `alphabet`, returning the index of the first invalid byte on failure.
#[inline]
pub fn validate_bases_with(buf: &[u8], alphabet: Alphabet) -> Result<(), usize> {
    let lut = lut_for(alphabet);

    #[cfg(target_arch = "x86_64")]
    {
        if cpu_features().avx2 {
            // SAFETY: AVX2 confirmed at runtime.
            return unsafe { validate_bases_avx2(buf, tables_for(alphabet), lut) };
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        if cpu_features().neon {
            // SAFETY: NEON confirmed at runtime.
            return unsafe { validate_bases_neon(buf, tables_for(alphabet), lut) };
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
unsafe fn validate_bases_avx2(
    buf: &[u8],
    tables: &NibbleTables,
    lut: &AlphabetLut,
) -> Result<(), usize> {
    use std::arch::x86_64::*;

    let len = buf.len();
    let ptr = buf.as_ptr();
    // SAFETY: both tables are 16 bytes; the loads stay in bounds.
    let lo_tbl = unsafe {
        _mm256_broadcastsi128_si256(_mm_loadu_si128(tables.0.as_ptr() as *const __m128i))
    };
    let hi_tbl = unsafe {
        _mm256_broadcastsi128_si256(_mm_loadu_si128(tables.1.as_ptr() as *const __m128i))
    };
    let low_mask = _mm256_set1_epi8(0x0F);
    let zero = _mm256_setzero_si256();

    let mut i = 0usize;
    while i + 32 <= len {
        // SAFETY: 32 bytes remain at `i`.
        let v = unsafe { _mm256_loadu_si256(ptr.add(i) as *const __m256i) };
        let lo_idx = _mm256_and_si256(v, low_mask);
        // `srli_epi16` shifts in neighbouring bits; the mask drops them again.
        let hi_idx = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask);
        let lo = _mm256_shuffle_epi8(lo_tbl, lo_idx);
        let hi = _mm256_shuffle_epi8(hi_tbl, hi_idx);
        let ok = _mm256_and_si256(lo, hi);
        let bad = _mm256_cmpeq_epi8(ok, zero);
        let mask = _mm256_movemask_epi8(bad) as u32;
        if mask != 0 {
            return Err(i + mask.trailing_zeros() as usize);
        }
        i += 32;
    }
    validate_bases_scalar(&buf[i..], lut).map_err(|e| i + e)
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
unsafe fn validate_bases_neon(
    buf: &[u8],
    tables: &NibbleTables,
    lut: &AlphabetLut,
) -> Result<(), usize> {
    use std::arch::aarch64::*;

    let len = buf.len();
    let ptr = buf.as_ptr();
    // SAFETY: both tables are 16 bytes.
    let lo_tbl = unsafe { vld1q_u8(tables.0.as_ptr()) };
    let hi_tbl = unsafe { vld1q_u8(tables.1.as_ptr()) };
    let zero = vdupq_n_u8(0);
    let low_mask = vdupq_n_u8(0x0F);

    let mut i = 0usize;
    while i + 16 <= len {
        // SAFETY: 16 bytes remain at `i`.
        let v = unsafe { vld1q_u8(ptr.add(i)) };
        let lo_idx = vandq_u8(v, low_mask);
        // Indices >= 16 (bytes >= 0x80 have none after the shift) read as zero from `vqtbl1q`.
        let hi_idx = vshrq_n_u8(v, 4);
        let lo = vqtbl1q_u8(lo_tbl, lo_idx);
        let hi = vqtbl1q_u8(hi_tbl, hi_idx);
        let ok = vandq_u8(lo, hi);
        let bad = vceqq_u8(ok, zero);
        // vshrn nibble-mask: each set lane contributes 4 bits, so bit index / 4 = lane index.
        let narrow = vshrn_n_u16(vreinterpretq_u16_u8(bad), 4);
        let bits = vget_lane_u64(vreinterpret_u64_u8(narrow), 0);
        if bits != 0 {
            return Err(i + (bits.trailing_zeros() as usize >> 2));
        }
        i += 16;
    }
    validate_bases_scalar(&buf[i..], lut).map_err(|e| i + e)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scalar(buf: &[u8], alphabet: Alphabet) -> Result<(), usize> {
        validate_bases_scalar(buf, lut_for(alphabet))
    }

    const ALPHABETS: [Alphabet; 3] = [Alphabet::AcgtnStrict, Alphabet::AcgtnCase, Alphabet::Iupac];

    #[test]
    fn matches_scalar_for_every_byte_value() {
        for alphabet in ALPHABETS {
            for b in 0u8..=255 {
                let buf = [b; 64];
                assert_eq!(
                    validate_bases_with(&buf, alphabet),
                    scalar(&buf, alphabet),
                    "byte 0x{b:02x} under {alphabet:?}"
                );
            }
        }
    }

    #[test]
    fn reports_first_bad_index_at_every_position() {
        for alphabet in ALPHABETS {
            for len in [1usize, 15, 16, 31, 32, 33, 64, 129] {
                for bad_at in 0..len {
                    let mut buf = vec![b'A'; len];
                    buf[bad_at] = b'!';
                    assert_eq!(validate_bases_with(&buf, alphabet), Err(bad_at));
                }
            }
        }
    }

    #[test]
    fn accepts_alphabet_members() {
        assert!(validate_bases_with(b"ACGTN", Alphabet::AcgtnStrict).is_ok());
        assert!(validate_bases_with(b"acgtn", Alphabet::AcgtnStrict).is_err());
        assert!(validate_bases_with(b"acgtnACGTN", Alphabet::AcgtnCase).is_ok());
        assert!(validate_bases_with(b"RYSWKMBDHVN.-", Alphabet::Iupac).is_ok());
        assert!(validate_bases_with(b"ACGT*", Alphabet::Iupac).is_err());
    }

    #[test]
    fn rejects_high_bytes() {
        for alphabet in ALPHABETS {
            let mut buf = vec![b'A'; 40];
            buf[37] = 0xC3;
            assert_eq!(validate_bases_with(&buf, alphabet), Err(37));
        }
    }

    #[test]
    fn empty_is_valid() {
        assert!(validate_bases(b"").is_ok());
    }
}
