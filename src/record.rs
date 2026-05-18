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
}
