//! Line-terminator scanning.
//!
//! This used to carry hand-written AVX2/AVX-512/NEON kernels. Benchmarked against `memchr` on
//! realistic FASTQ line lengths they were a wash on AVX2 (8.9 vs 8.9 GB/s on a 150 bp read set),
//! while the AVX-512 kernel risks frequency throttling on several Intel generations. `memchr`
//! already dispatches SSE2/AVX2/NEON kernels at runtime and handles the unaligned head and tail
//! more carefully than the code it replaced, so this module is now a thin wrapper: same API,
//! less `unsafe`.

/// Index of the first `\n` at or after `start`, or `None` if there is none.
#[inline]
pub fn find_lf(buf: &[u8], start: usize) -> Option<usize> {
    if start >= buf.len() {
        return None;
    }
    memchr::memchr(b'\n', &buf[start..]).map(|idx| start + idx)
}

/// Number of `\n` bytes in `buf`. Useful for record counting and progress reporting.
#[inline]
pub fn count_lf(buf: &[u8]) -> usize {
    memchr::memchr_iter(b'\n', buf).count()
}

/// True when `buf` holds a `\n` or `\r`, i.e. would break FASTQ line framing if written out.
#[inline]
pub fn contains_line_break(buf: &[u8]) -> bool {
    memchr::memchr2(b'\n', b'\r', buf).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_from_offset() {
        let buf = b"abc\ndef\n";
        assert_eq!(find_lf(buf, 0), Some(3));
        assert_eq!(find_lf(buf, 4), Some(7));
        assert_eq!(find_lf(buf, 8), None);
        assert_eq!(find_lf(buf, 99), None);
    }

    #[test]
    fn counts_and_detects_breaks() {
        assert_eq!(count_lf(b"a\nb\nc"), 2);
        assert!(contains_line_break(b"ac\rgt"));
        assert!(!contains_line_break(b"ACGT"));
    }
}
