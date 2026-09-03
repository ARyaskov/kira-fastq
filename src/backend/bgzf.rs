//! BGZF input: mmap over the compressed file, one block inflated at a time.
//!
//! BGZF is gzip with a per-block `BC` extra subfield giving the compressed block size, which is
//! what makes random access by virtual offset possible. This backend keeps htslib-compatible
//! semantics:
//!
//! * `tell()` returns `coffset << 16 | uoffset`. When a block is fully consumed it reports the
//!   *next* block at uoffset 0 rather than an unrepresentable uoffset of 65536, so the offset
//!   round-trips for blocks that hold the full 64 KiB the spec allows.
//! * Every block's CRC32 and ISIZE are verified.
//! * A missing 28-byte end-of-file marker is reported as an error, since that is how a truncated
//!   BGZF file presents itself.
//! * `seek` uses a lazily built block index, so repeated seeks do not rescan the file.

use std::fs::File;
use std::path::Path;

#[cfg(not(feature = "libdeflate"))]
use flate2::{Decompress, FlushDecompress, Status};
use memmap2::{Mmap, MmapOptions};

use crate::backend::{ByteSource, LineStatus, read_line_from};
use crate::error::{FastqError, InvalidKind};
use crate::offset::VirtualOffset;

/// Maximum uncompressed payload of a BGZF block, per the spec.
pub(crate) const MAX_BLOCK_UNCOMP: usize = 1 << 16;

/// The canonical 28-byte empty block every BGZF writer appends to mark the end of the file.
pub(crate) const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, b'B', b'C', 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// One entry of the lazily built block index: compressed offset and the number of decoded bytes
/// that precede the block.
#[derive(Clone, Copy)]
struct IndexEntry {
    coffset: usize,
    uncompressed_before: u64,
}

pub(crate) struct BgzfBackend {
    mmap: Mmap,
    /// Compressed offset of the next block to decode.
    src_pos: usize,
    /// Compressed offset of the block currently in `buf`.
    block_start: usize,
    buf: Vec<u8>,
    buf_pos: usize,
    buf_len: usize,
    /// Decoded bytes before the block currently in `buf`.
    total_uncomp_before: u64,
    finished: bool,
    eof_check: bool,
    has_eof_marker: bool,
    index: Vec<IndexEntry>,
    decoder: BlockDecoder,
}

impl BgzfBackend {
    pub(crate) fn new(path: &Path) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        let metadata = file.metadata()?;
        if metadata.len() == 0 {
            return Err(FastqError::invalid(0, InvalidKind::BgzfHeader));
        }
        // SAFETY: the mapping spans the whole file; see `MmapBackend::open`.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        let has_eof_marker =
            mmap.len() >= BGZF_EOF.len() && mmap[mmap.len() - BGZF_EOF.len()..] == BGZF_EOF[..];
        Ok(Self {
            mmap,
            src_pos: 0,
            block_start: 0,
            buf: vec![0u8; MAX_BLOCK_UNCOMP],
            buf_pos: 0,
            buf_len: 0,
            total_uncomp_before: 0,
            finished: false,
            eof_check: true,
            has_eof_marker,
            index: Vec::new(),
            decoder: BlockDecoder::new(),
        })
    }

    /// Turn the end-of-file marker check off, for readers of deliberately partial files.
    #[inline]
    pub(crate) fn set_eof_check(&mut self, enabled: bool) {
        self.eof_check = enabled;
    }

    #[inline]
    pub(crate) fn logical_offset(&self) -> u64 {
        self.total_uncomp_before + self.buf_pos as u64
    }

    /// Virtual offset of the next byte the parser will see.
    pub(crate) fn tell(&self) -> VirtualOffset {
        if self.buf_pos >= self.buf_len {
            // Current block is drained: address the start of the next one. A uoffset of 65536
            // does not fit the 16-bit field, so this is the only correct encoding.
            VirtualOffset::new(self.src_pos as u64, 0)
        } else {
            VirtualOffset::new(self.block_start as u64, self.buf_pos as u16)
        }
    }

    pub(crate) fn seek(&mut self, voff: VirtualOffset) -> Result<(), FastqError> {
        let block_off = voff.compressed() as usize;
        let uoff = voff.uncompressed() as usize;

        if block_off == self.mmap.len() && uoff == 0 {
            // The end-of-file virtual offset: a legitimate checkpoint at the end of the data.
            self.total_uncomp_before = self.uncompressed_before(block_off)?;
            self.src_pos = block_off;
            self.block_start = block_off;
            self.buf_pos = 0;
            self.buf_len = 0;
            self.finished = true;
            return Ok(());
        }
        if block_off >= self.mmap.len() {
            return Err(FastqError::invalid(0, InvalidKind::BgzfVirtualOffset));
        }

        let total = self.uncompressed_before(block_off)?;
        self.total_uncomp_before = total;
        self.src_pos = block_off;
        self.buf_pos = 0;
        self.buf_len = 0;
        self.finished = false;
        self.decode_block()?;
        if uoff > self.buf_len {
            return Err(FastqError::invalid(0, InvalidKind::BgzfVirtualOffset));
        }
        self.buf_pos = uoff;
        Ok(())
    }

    #[inline]
    pub(crate) fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        read_line_from(self, out)
    }

    /// Decoded bytes preceding the block at `target`, walking (and caching) block headers only.
    fn uncompressed_before(&mut self, target: usize) -> Result<u64, FastqError> {
        if target == 0 {
            return Ok(0);
        }
        let (mut pos, mut total) = match self
            .index
            .binary_search_by_key(&target, |entry| entry.coffset)
        {
            Ok(i) => return Ok(self.index[i].uncompressed_before),
            Err(0) => (0usize, 0u64),
            Err(i) => {
                let entry = self.index[i - 1];
                (entry.coffset, entry.uncompressed_before)
            }
        };
        while pos < target {
            let (block_size, isize) = parse_block_sizes(&self.mmap, pos)?;
            total = total
                .checked_add(u64::from(isize))
                .ok_or_else(|| FastqError::invalid(0, InvalidKind::BgzfVirtualOffset))?;
            pos = pos
                .checked_add(block_size)
                .ok_or_else(|| FastqError::invalid(0, InvalidKind::BgzfVirtualOffset))?;
            if pos > self.mmap.len() {
                return Err(FastqError::invalid(0, InvalidKind::BgzfVirtualOffset));
            }
            self.index.push(IndexEntry {
                coffset: pos,
                uncompressed_before: total,
            });
        }
        if pos != target {
            // The offset points into the middle of a block, not at a block boundary.
            return Err(FastqError::invalid(0, InvalidKind::BgzfVirtualOffset));
        }
        Ok(total)
    }

    fn decode_block(&mut self) -> Result<(), FastqError> {
        if self.src_pos >= self.mmap.len() {
            self.finished = true;
            self.block_start = self.src_pos;
            self.buf_pos = 0;
            self.buf_len = 0;
            return Ok(());
        }
        let offset = self.logical_offset();
        self.block_start = self.src_pos;
        let decoded = self
            .decoder
            .decode(&self.mmap, self.src_pos, &mut self.buf, offset)?;
        self.buf_pos = 0;
        self.buf_len = decoded.uncompressed_len;
        self.src_pos += decoded.block_size;
        Ok(())
    }
}

/// Result of decoding one BGZF block.
pub(crate) struct DecodedBlock {
    /// Size of the compressed block, i.e. how far to advance in the file.
    pub(crate) block_size: usize,
    /// Bytes written to the output buffer.
    pub(crate) uncompressed_len: usize,
}

/// Inflates single BGZF blocks and verifies their trailer.
///
/// Shared by the sequential backend and the worker threads of the parallel reader, so both
/// enforce the same integrity checks.
pub(crate) struct BlockDecoder {
    #[cfg(not(feature = "libdeflate"))]
    decomp: Decompress,
    #[cfg(feature = "libdeflate")]
    libdeflater: libdeflater::Decompressor,
}

impl BlockDecoder {
    pub(crate) fn new() -> Self {
        Self {
            #[cfg(not(feature = "libdeflate"))]
            decomp: Decompress::new(false),
            #[cfg(feature = "libdeflate")]
            libdeflater: libdeflater::Decompressor::new(),
        }
    }

    /// Decode the block starting at `start` into `out`, which must hold at least
    /// [`MAX_BLOCK_UNCOMP`] bytes. `offset` only labels errors.
    pub(crate) fn decode(
        &mut self,
        data: &[u8],
        start: usize,
        out: &mut [u8],
        offset: u64,
    ) -> Result<DecodedBlock, FastqError> {
        debug_assert!(out.len() >= MAX_BLOCK_UNCOMP);
        let (block_size, isize) = parse_block_sizes(data, start)?;
        let block_end = start
            .checked_add(block_size)
            .filter(|end| *end <= data.len())
            .ok_or_else(|| FastqError::invalid(offset, InvalidKind::BgzfBlock))?;
        let (deflate_start, deflate_end) = parse_deflate_range(data, start)?;
        if deflate_end > block_end || deflate_start > deflate_end {
            return Err(FastqError::invalid(offset, InvalidKind::BgzfBlock));
        }
        let input = &data[deflate_start..deflate_end];

        #[cfg(feature = "libdeflate")]
        let produced = if isize == 0 {
            // libdeflate rejects a zero-length output buffer; the EOF marker is exactly that.
            0usize
        } else {
            self.libdeflater
                .deflate_decompress(input, &mut out[..isize as usize])
                .map_err(|_| FastqError::invalid(offset, InvalidKind::BgzfBlock))?
        };

        #[cfg(not(feature = "libdeflate"))]
        let produced = {
            self.decomp.reset(false);
            let before_out = self.decomp.total_out();
            let status = self
                .decomp
                .decompress(input, out, FlushDecompress::Finish)
                .map_err(|_| FastqError::invalid(offset, InvalidKind::BgzfBlock))?;
            let produced = (self.decomp.total_out() - before_out) as usize;
            // A zero-length deflate payload never reports StreamEnd; that is the EOF marker.
            if status != Status::StreamEnd && !(produced == 0 && input.is_empty()) {
                return Err(FastqError::invalid(offset, InvalidKind::BgzfBlock));
            }
            produced
        };

        if produced > MAX_BLOCK_UNCOMP {
            return Err(FastqError::invalid(offset, InvalidKind::BgzfBlockTooLarge));
        }
        if produced != isize as usize {
            return Err(FastqError::invalid(offset, InvalidKind::BgzfBlockIsize));
        }
        let mut hasher = crc32fast::Hasher::new();
        hasher.update(&out[..produced]);
        if hasher.finalize() != block_crc(data, start, block_size)? {
            return Err(FastqError::invalid(offset, InvalidKind::BgzfBlockCrc));
        }

        Ok(DecodedBlock {
            block_size,
            uncompressed_len: produced,
        })
    }
}

impl ByteSource for BgzfBackend {
    #[inline]
    fn available(&self) -> &[u8] {
        if self.buf_pos >= self.buf_len {
            &[]
        } else {
            &self.buf[self.buf_pos..self.buf_len]
        }
    }

    #[inline]
    fn consume(&mut self, n: usize) {
        self.buf_pos = self.buf_pos.saturating_add(n).min(self.buf_len);
    }

    fn refill(&mut self) -> Result<bool, FastqError> {
        loop {
            if self.buf_pos < self.buf_len {
                return Ok(true);
            }
            if self.finished {
                return Ok(false);
            }
            if self.src_pos >= self.mmap.len() {
                self.finished = true;
                if self.eof_check && !self.has_eof_marker {
                    return Err(FastqError::invalid(
                        self.logical_offset(),
                        InvalidKind::BgzfMissingEofMarker,
                    ));
                }
                return Ok(false);
            }
            self.total_uncomp_before = self
                .total_uncomp_before
                .checked_add(self.buf_len as u64)
                .ok_or_else(|| {
                    FastqError::invalid(self.logical_offset(), InvalidKind::BgzfBlock)
                })?;
            self.decode_block()?;
        }
    }
}

/// `(block size in bytes, uncompressed size)` of the block starting at `start`.
pub(crate) fn parse_block_sizes(data: &[u8], start: usize) -> Result<(usize, u32), FastqError> {
    let (_header_end, bsize) = parse_bgzf_header(data, start)?;
    let block_size = (bsize as usize)
        .checked_add(1)
        .ok_or_else(|| FastqError::invalid(0, InvalidKind::BgzfBlock))?;
    // 12-byte fixed header + 6-byte BC subfield + >= 0 deflate bytes + 8-byte trailer.
    if block_size < 26 {
        return Err(FastqError::invalid(0, InvalidKind::BgzfBlock));
    }
    let trailer_pos = start
        .checked_add(block_size)
        .and_then(|v| v.checked_sub(8))
        .ok_or_else(|| FastqError::invalid(0, InvalidKind::BgzfBlock))?;
    if trailer_pos.saturating_add(8) > data.len() {
        return Err(FastqError::invalid(0, InvalidKind::BgzfBlock));
    }
    let isize = u32::from_le_bytes([
        data[trailer_pos + 4],
        data[trailer_pos + 5],
        data[trailer_pos + 6],
        data[trailer_pos + 7],
    ]);
    if (isize as usize) > MAX_BLOCK_UNCOMP {
        return Err(FastqError::invalid(0, InvalidKind::BgzfBlockTooLarge));
    }
    Ok((block_size, isize))
}

fn block_crc(data: &[u8], start: usize, block_size: usize) -> Result<u32, FastqError> {
    let trailer_pos = start
        .checked_add(block_size)
        .and_then(|v| v.checked_sub(8))
        .filter(|pos| pos.saturating_add(8) <= data.len())
        .ok_or_else(|| FastqError::invalid(0, InvalidKind::BgzfBlock))?;
    Ok(u32::from_le_bytes([
        data[trailer_pos],
        data[trailer_pos + 1],
        data[trailer_pos + 2],
        data[trailer_pos + 3],
    ]))
}

fn parse_deflate_range(data: &[u8], start: usize) -> Result<(usize, usize), FastqError> {
    let (header_end, bsize) = parse_bgzf_header(data, start)?;
    let block_size = (bsize as usize)
        .checked_add(1)
        .ok_or_else(|| FastqError::invalid(0, InvalidKind::BgzfBlock))?;
    let deflate_end = start
        .checked_add(block_size)
        .and_then(|v| v.checked_sub(8))
        .ok_or_else(|| FastqError::invalid(0, InvalidKind::BgzfBlock))?;
    Ok((header_end, deflate_end))
}

/// `(offset just past the extra field, BSIZE)` of the BGZF block starting at `start`.
fn parse_bgzf_header(data: &[u8], start: usize) -> Result<(usize, u16), FastqError> {
    let bad = || FastqError::invalid(0, InvalidKind::BgzfHeader);
    let after_fixed = start.checked_add(12).ok_or_else(bad)?;
    if after_fixed > data.len() {
        return Err(bad());
    }
    if data[start] != 0x1f || data[start + 1] != 0x8b || data[start + 2] != 8 {
        return Err(bad());
    }
    let flg = data[start + 3];
    if (flg & 0x04) == 0 {
        return Err(bad());
    }
    let mut i = start + 10;
    let xlen = u16::from_le_bytes([data[i], data[i + 1]]) as usize;
    i += 2;
    let extra_end = i.checked_add(xlen).ok_or_else(bad)?;
    if extra_end > data.len() {
        return Err(bad());
    }
    let mut bsize = None;
    while i.saturating_add(4) <= extra_end {
        let si1 = data[i];
        let si2 = data[i + 1];
        let slen = u16::from_le_bytes([data[i + 2], data[i + 3]]) as usize;
        i += 4;
        let sub_end = i.checked_add(slen).ok_or_else(bad)?;
        if sub_end > extra_end {
            return Err(bad());
        }
        if si1 == b'B' && si2 == b'C' && slen == 2 {
            bsize = Some(u16::from_le_bytes([data[i], data[i + 1]]));
        }
        i = sub_end;
    }
    let bsize = bsize.ok_or_else(bad)?;
    Ok((extra_end, bsize))
}
