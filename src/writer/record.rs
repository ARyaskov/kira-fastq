//! Record assembly into a contiguous scratch buffer.
//!
//! The layout is:
//!
//! ```text
//! @ <header> \n
//! <seq>     \n
//! +         \n
//! <qual>    \n
//! ```
//!
//! No CRLF on output — kira always writes LF, matching every modern aligner's expectation.
//! Readers handle both; writers shouldn't propagate the ambiguity.

/// Pack a FASTQ record into `out`. `out` is cleared first and grown to exactly fit.
/// Uses [`Vec::extend_from_slice`] which lowers to `memcpy`; no per-byte writes.
#[inline]
pub fn assemble_record(out: &mut Vec<u8>, header: &[u8], seq: &[u8], qual: &[u8]) {
    // Layout: '@' + header + '\n' + seq + '\n' + '+' + '\n' + qual + '\n', so six bytes of
    // framing. Getting this count wrong costs a reallocation every time a record sets a new
    // high-water mark.
    let total = header.len() + seq.len() + qual.len() + 6;
    out.clear();
    out.reserve(total);
    out.push(b'@');
    out.extend_from_slice(header);
    out.push(b'\n');
    out.extend_from_slice(seq);
    out.push(b'\n');
    out.push(b'+');
    out.push(b'\n');
    out.extend_from_slice(qual);
    out.push(b'\n');
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_minimal() {
        let mut buf = Vec::new();
        assemble_record(&mut buf, b"r0", b"ACGT", b"!!!!");
        assert_eq!(buf, b"@r0\nACGT\n+\n!!!!\n");
    }

    #[test]
    fn empty_seq_qual() {
        let mut buf = Vec::new();
        assemble_record(&mut buf, b"empty", b"", b"");
        assert_eq!(buf, b"@empty\n\n+\n\n");
    }

    #[test]
    fn reuses_capacity() {
        let mut buf = Vec::with_capacity(0);
        assemble_record(&mut buf, b"a", b"AC", b"!!");
        let cap_first = buf.capacity();
        assemble_record(&mut buf, b"b", b"GT", b"@@");
        let cap_second = buf.capacity();
        assert_eq!(cap_first, cap_second, "no realloc on equal-sized record");
    }
}
