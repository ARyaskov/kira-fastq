#[derive(Debug, Clone, Copy)]
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

    #[inline]
    pub fn header(&self) -> &'a [u8] {
        self.header
    }

    #[inline]
    pub fn seq(&self) -> &'a [u8] {
        self.seq
    }

    #[inline]
    pub fn qual(&self) -> &'a [u8] {
        self.qual
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.seq.len()
    }
}
