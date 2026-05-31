/// Borrowed FASTQ record. Fields are slices into a backing buffer (mmap, gzip-scratch,
/// bgzf-scratch, or a stream-scratch). Cheap to copy by value.
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

    /// Header bytes with the leading `@` stripped.
    #[inline]
    pub fn header(&self) -> &'a [u8] {
        self.header
    }

    #[inline]
    pub fn seq(&self) -> &'a [u8] {
        self.seq
    }

    /// Phred+33 by convention.
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
    /// pipelines (e.g. when yielding from an async [`Stream`] or crossing thread boundaries).
    #[inline]
    pub fn to_owned(&self) -> FastqRecordOwned {
        FastqRecordOwned {
            header: self.header.to_vec(),
            seq: self.seq.to_vec(),
            qual: self.qual.to_vec(),
        }
    }
}

/// Owned FASTQ record. Use this when the record must outlive the reader's scratch buffer
/// (async [`Stream`] items, channel sends, etc.). Allocates three [`Vec<u8>`]s per record.
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
