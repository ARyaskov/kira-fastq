/// Borrowed FASTQ record. Fields are slices into a backing buffer (mmap, an in-memory buffer,
/// or the reader's scratch). Cheap to copy by value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FastqRecord<'a> {
    header: &'a [u8],
    seq: &'a [u8],
    qual: &'a [u8],
}

impl<'a> FastqRecord<'a> {
    #[inline]
    pub fn new(header: &'a [u8], seq: &'a [u8], qual: &'a [u8]) -> Self {
        Self { header, seq, qual }
    }

    /// Full header line with the leading `@` stripped: read ID plus any description.
    #[inline]
    pub fn header(&self) -> &'a [u8] {
        self.header
    }

    /// Read ID: the header up to the first space or tab.
    #[inline]
    pub fn id(&self) -> &'a [u8] {
        split_header(self.header).0
    }

    /// Description: whatever follows the first run of spaces or tabs in the header. Empty when
    /// the header holds only an ID.
    #[inline]
    pub fn description(&self) -> &'a [u8] {
        split_header(self.header).1
    }

    #[inline]
    pub fn seq(&self) -> &'a [u8] {
        self.seq
    }

    /// Quality bytes, Phred+33 unless the producer said otherwise. See
    /// [`crate::validation::QualityEncoding`].
    #[inline]
    pub fn qual(&self) -> &'a [u8] {
        self.qual
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Copy fields into a heap-owned [`FastqRecordOwned`]. Used at the boundary of zero-copy
    /// pipelines, e.g. when yielding from an async stream or crossing thread boundaries.
    #[inline]
    #[allow(clippy::wrong_self_convention)]
    pub fn to_owned(&self) -> FastqRecordOwned {
        FastqRecordOwned {
            header: self.header.to_vec(),
            seq: self.seq.to_vec(),
            qual: self.qual.to_vec(),
        }
    }
}

/// Owned FASTQ record. Use this when the record must outlive the reader's scratch buffer
/// (async stream items, channel sends, and so on). Holds three [`Vec<u8>`]s; reuse one across
/// records with [`FastqRecordOwned::copy_from`] to keep the allocations.
#[derive(Debug, Default, Clone, PartialEq, Eq, Hash)]
pub struct FastqRecordOwned {
    header: Vec<u8>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl FastqRecordOwned {
    #[inline]
    pub fn new(header: Vec<u8>, seq: Vec<u8>, qual: Vec<u8>) -> Self {
        Self { header, seq, qual }
    }

    #[inline]
    pub fn header(&self) -> &[u8] {
        &self.header
    }

    /// Read ID: the header up to the first space or tab.
    #[inline]
    pub fn id(&self) -> &[u8] {
        split_header(&self.header).0
    }

    /// Description: whatever follows the first run of spaces or tabs in the header.
    #[inline]
    pub fn description(&self) -> &[u8] {
        split_header(&self.header).1
    }

    #[inline]
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    #[inline]
    pub fn qual(&self) -> &[u8] {
        &self.qual
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Borrow as a [`FastqRecord`]. Zero-cost.
    #[inline]
    pub fn as_borrowed(&self) -> FastqRecord<'_> {
        FastqRecord {
            header: &self.header,
            seq: &self.seq,
            qual: &self.qual,
        }
    }

    /// Overwrite this record with `rec`, reusing the existing allocations. Prefer this over
    /// [`FastqRecord::to_owned`] in loops: no allocation once the buffers are large enough.
    #[inline]
    pub fn copy_from(&mut self, rec: &FastqRecord<'_>) {
        self.header.clear();
        self.header.extend_from_slice(rec.header());
        self.seq.clear();
        self.seq.extend_from_slice(rec.seq());
        self.qual.clear();
        self.qual.extend_from_slice(rec.qual());
    }

    /// Empty all three buffers, keeping their capacity.
    #[inline]
    pub fn clear(&mut self) {
        self.header.clear();
        self.seq.clear();
        self.qual.clear();
    }

    /// Mutable header buffer. Useful for in-place transformations.
    #[inline]
    pub fn header_mut(&mut self) -> &mut Vec<u8> {
        &mut self.header
    }

    #[inline]
    pub fn seq_mut(&mut self) -> &mut Vec<u8> {
        &mut self.seq
    }

    #[inline]
    pub fn qual_mut(&mut self) -> &mut Vec<u8> {
        &mut self.qual
    }
}

impl From<&FastqRecord<'_>> for FastqRecordOwned {
    #[inline]
    fn from(rec: &FastqRecord<'_>) -> Self {
        rec.to_owned()
    }
}

/// Split a header into (id, description) at the first space or tab.
#[inline]
fn split_header(header: &[u8]) -> (&[u8], &[u8]) {
    match header.iter().position(|&b| b == b' ' || b == b'\t') {
        Some(i) => {
            let rest = &header[i..];
            let start = rest
                .iter()
                .position(|&b| b != b' ' && b != b'\t')
                .unwrap_or(rest.len());
            (&header[..i], &rest[start..])
        }
        None => (header, &[]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_id_and_description() {
        let rec = FastqRecord::new(b"SRR001 1:N:0:ATCACG", b"AC", b"!!");
        assert_eq!(rec.id(), b"SRR001");
        assert_eq!(rec.description(), b"1:N:0:ATCACG");
    }

    #[test]
    fn id_only_header_has_empty_description() {
        let rec = FastqRecord::new(b"SRR001", b"AC", b"!!");
        assert_eq!(rec.id(), b"SRR001");
        assert_eq!(rec.description(), b"");
    }

    #[test]
    fn copy_from_reuses_allocation() {
        let mut owned = FastqRecordOwned::default();
        owned.copy_from(&FastqRecord::new(b"a", b"ACGT", b"!!!!"));
        let cap = owned.seq_mut().capacity();
        owned.copy_from(&FastqRecord::new(b"b", b"TTTT", b"####"));
        assert_eq!(owned.header(), b"b");
        assert_eq!(owned.seq(), b"TTTT");
        assert_eq!(cap, owned.seq_mut().capacity(), "no realloc on equal sizes");
    }
}
