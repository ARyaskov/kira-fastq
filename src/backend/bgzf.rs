use std::fs::File;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};
#[cfg(not(feature = "libdeflate"))]
use miniz_oxide::inflate::TINFLStatus;
#[cfg(not(feature = "libdeflate"))]
use miniz_oxide::inflate::core::{DecompressorOxide, decompress, inflate_flags};

use crate::backend::gzip::LineStatus;
use crate::error::{FastqError, InvalidKind};
use crate::offset::VirtualOffset;
use crate::simd::newline::find_lf;

const OUT_BUF_SIZE: usize = 256 * 1024;
// BGZF spec bound; required for `tell()`/`seek()` virtual offsets to round-trip in 16 bits.
const MAX_BLOCK_UNCOMP: usize = 1 << 16;

pub struct BgzfBackend {
    mmap: Mmap,
    src_pos: usize,
    buf: Vec<u8>,
    buf_pos: usize,
    buf_len: usize,
    block_start: usize,
    block_uncomp_len: usize,
    total_uncomp_before: u64,
    finished: bool,
    #[cfg(not(feature = "libdeflate"))]
    decomp: DecompressorOxide,
    #[cfg(feature = "libdeflate")]
    libdeflater: libdeflater::Decompressor,
}

impl BgzfBackend {
    pub fn new(path: &Path) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        let metadata = file.metadata()?;
        if metadata.len() == 0 {
            return Err(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfHeader,
            });
        }
        // SAFETY: file is kept alive for the duration of the mmap; mapping spans the whole file.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(Self {
            mmap,
            src_pos: 0,
            buf: vec![0u8; OUT_BUF_SIZE],
            buf_pos: 0,
            buf_len: 0,
            block_start: 0,
            block_uncomp_len: 0,
            total_uncomp_before: 0,
            finished: false,
            #[cfg(not(feature = "libdeflate"))]
            decomp: DecompressorOxide::new(),
            #[cfg(feature = "libdeflate")]
            libdeflater: libdeflater::Decompressor::new(),
        })
    }

    #[inline]
    pub fn logical_offset(&self) -> u64 {
        self.total_uncomp_before + self.buf_pos as u64
    }

    #[inline]
    pub fn tell(&self) -> VirtualOffset {
        debug_assert!(self.buf_pos < MAX_BLOCK_UNCOMP);
        VirtualOffset(((self.block_start as u64) << 16) | (self.buf_pos as u64 & 0xFFFF))
    }

    pub fn seek(&mut self, voff: VirtualOffset) -> Result<(), FastqError> {
        let block_off = (voff.0 >> 16) as usize;
        let uoff = (voff.0 & 0xFFFF) as usize;
        if block_off >= self.mmap.len() {
            return Err(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfVirtualOffset,
            });
        }
        let (bsize, isize) = parse_block_sizes(&self.mmap, block_off)?;
        let block_end = block_off
            .checked_add(bsize)
            .ok_or(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfVirtualOffset,
            })?;
        if block_end > self.mmap.len() {
            return Err(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfVirtualOffset,
            });
        }
        let total = sum_uncomp_before(&self.mmap, block_off)?;
        self.total_uncomp_before = total;
        self.block_start = block_off;
        self.src_pos = block_off;
        self.buf_pos = 0;
        self.buf_len = 0;
        self.block_uncomp_len = isize as usize;
        self.finished = false;
        self.decode_block()?;
        if uoff > self.buf_len {
            return Err(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfVirtualOffset,
            });
        }
        self.buf_pos = uoff;
        Ok(())
    }

    pub fn read_line(&mut self, out: &mut Vec<u8>) -> Result<LineStatus, FastqError> {
        out.clear();
        loop {
            let slice = self.available_slice();
            if !slice.is_empty() {
                if let Some(lf) = find_lf(slice, 0) {
                    let mut end = lf;
                    if end > 0 && slice[end - 1] == b'\r' {
                        end -= 1;
                    }
                    out.extend_from_slice(&slice[..end]);
                    self.advance(lf + 1);
                    return Ok(LineStatus::Line);
                }
                out.extend_from_slice(slice);
                let n = slice.len();
                self.advance(n);
            }
            if !self.refill()? {
                if out.is_empty() {
                    return Ok(LineStatus::EofClean);
                }
                if out.last() == Some(&b'\r') {
                    out.pop();
                }
                return Ok(LineStatus::EofPartial);
            }
        }
    }

    pub fn peek_byte(&mut self) -> Result<Option<u8>, FastqError> {
        loop {
            let slice = self.available_slice();
            if let Some(&b) = slice.first() {
                return Ok(Some(b));
            }
            if !self.refill()? {
                return Ok(None);
            }
        }
    }

    #[inline]
    fn available_slice(&self) -> &[u8] {
        if self.buf_pos >= self.buf_len {
            &[]
        } else {
            &self.buf[self.buf_pos..self.buf_len]
        }
    }

    #[inline]
    fn advance(&mut self, n: usize) {
        self.buf_pos = self.buf_pos.saturating_add(n).min(self.buf_len);
    }

    fn refill(&mut self) -> Result<bool, FastqError> {
        if self.buf_pos < self.buf_len {
            return Ok(true);
        }
        if self.finished {
            return Ok(false);
        }
        self.total_uncomp_before = self
            .total_uncomp_before
            .checked_add(self.block_uncomp_len as u64)
            .ok_or_else(|| FastqError::InvalidFormat {
                offset: self.logical_offset(),
                kind: InvalidKind::BgzfBlock,
            })?;
        self.decode_block()?;
        Ok(self.buf_len > self.buf_pos)
    }

    fn decode_block(&mut self) -> Result<(), FastqError> {
        if self.src_pos >= self.mmap.len() {
            self.finished = true;
            self.buf_pos = 0;
            self.buf_len = 0;
            self.block_uncomp_len = 0;
            return Ok(());
        }
        self.block_start = self.src_pos;
        let (block_size, _isize) = parse_block_sizes(&self.mmap, self.src_pos)?;
        let block_end =
            self.src_pos
                .checked_add(block_size)
                .ok_or_else(|| FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::BgzfBlock,
                })?;
        if block_end > self.mmap.len() {
            return Err(FastqError::InvalidFormat {
                offset: self.logical_offset(),
                kind: InvalidKind::BgzfBlock,
            });
        }
        let (deflate_start, deflate_end) = parse_deflate_range(&self.mmap, self.src_pos)?;
        if deflate_end > block_end || deflate_start > deflate_end {
            return Err(FastqError::InvalidFormat {
                offset: self.logical_offset(),
                kind: InvalidKind::BgzfBlock,
            });
        }
        let in_buf = &self.mmap[deflate_start..deflate_end];

        #[cfg(feature = "libdeflate")]
        let out_consumed = match self
            .libdeflater
            .deflate_decompress(in_buf, &mut self.buf[..MAX_BLOCK_UNCOMP])
        {
            Ok(n) => n,
            Err(_) => {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::BgzfBlock,
                });
            }
        };

        #[cfg(not(feature = "libdeflate"))]
        let out_consumed = {
            self.decomp = DecompressorOxide::new();
            let flags = inflate_flags::TINFL_FLAG_USING_NON_WRAPPING_OUTPUT_BUF;
            let (status, _in_consumed, out_consumed) =
                decompress(&mut self.decomp, in_buf, &mut self.buf, 0, flags);
            match status {
                TINFLStatus::Done => out_consumed,
                _ => {
                    return Err(FastqError::InvalidFormat {
                        offset: self.logical_offset(),
                        kind: InvalidKind::BgzfBlock,
                    });
                }
            }
        };

        if out_consumed > MAX_BLOCK_UNCOMP {
            return Err(FastqError::InvalidFormat {
                offset: self.logical_offset(),
                kind: InvalidKind::BgzfBlockTooLarge,
            });
        }

        self.buf_pos = 0;
        self.buf_len = out_consumed;
        self.block_uncomp_len = out_consumed;
        self.src_pos = block_end;
        if self.buf_len == 0 && self.src_pos >= self.mmap.len() {
            self.finished = true;
        }
        Ok(())
    }
}

fn parse_block_sizes(data: &[u8], start: usize) -> Result<(usize, u32), FastqError> {
    let (_header_end, bsize) = parse_bgzf_header(data, start)?;
    let block_size = (bsize as usize)
        .checked_add(1)
        .ok_or(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfBlock,
        })?;
    // 12-byte fixed header + 6-byte BC subfield + ≥0 deflate + 8-byte trailer.
    if block_size < 26 {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfBlock,
        });
    }
    let trailer_pos = start
        .checked_add(block_size)
        .and_then(|v| v.checked_sub(8))
        .ok_or(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfBlock,
        })?;
    let trailer_end = trailer_pos
        .checked_add(8)
        .ok_or(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfBlock,
        })?;
    if trailer_end > data.len() {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfBlock,
        });
    }
    let isize = u32::from_le_bytes([
        data[trailer_pos + 4],
        data[trailer_pos + 5],
        data[trailer_pos + 6],
        data[trailer_pos + 7],
    ]);
    if (isize as usize) > MAX_BLOCK_UNCOMP {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfBlockTooLarge,
        });
    }
    Ok((block_size, isize))
}

fn parse_deflate_range(data: &[u8], start: usize) -> Result<(usize, usize), FastqError> {
    let (header_end, bsize) = parse_bgzf_header(data, start)?;
    let block_size = (bsize as usize)
        .checked_add(1)
        .ok_or(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfBlock,
        })?;
    let deflate_start = header_end;
    let deflate_end = start
        .checked_add(block_size)
        .and_then(|v| v.checked_sub(8))
        .ok_or(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfBlock,
        })?;
    Ok((deflate_start, deflate_end))
}

fn parse_bgzf_header(data: &[u8], start: usize) -> Result<(usize, u16), FastqError> {
    let after_fixed = start.checked_add(10).ok_or(FastqError::InvalidFormat {
        offset: 0,
        kind: InvalidKind::BgzfHeader,
    })?;
    if after_fixed > data.len() {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfHeader,
        });
    }
    if data[start] != 0x1f || data[start + 1] != 0x8b || data[start + 2] != 8 {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfHeader,
        });
    }
    let flg = data[start + 3];
    if (flg & 0x04) == 0 {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfHeader,
        });
    }
    let mut i = after_fixed;
    if i.checked_add(2).is_none_or(|v| v > data.len()) {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfHeader,
        });
    }
    let xlen = u16::from_le_bytes([data[i], data[i + 1]]) as usize;
    i += 2;
    let extra_end = i.checked_add(xlen).ok_or(FastqError::InvalidFormat {
        offset: 0,
        kind: InvalidKind::BgzfHeader,
    })?;
    if extra_end > data.len() {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfHeader,
        });
    }
    let mut bsize = None;
    while i.saturating_add(4) <= extra_end {
        let si1 = data[i];
        let si2 = data[i + 1];
        let slen = u16::from_le_bytes([data[i + 2], data[i + 3]]) as usize;
        i += 4;
        let sub_end = i.checked_add(slen).ok_or(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfHeader,
        })?;
        if sub_end > extra_end {
            return Err(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfHeader,
            });
        }
        if si1 == b'B' && si2 == b'C' && slen == 2 {
            bsize = Some(u16::from_le_bytes([data[i], data[i + 1]]));
        }
        i = sub_end;
    }
    let bsize = bsize.ok_or(FastqError::InvalidFormat {
        offset: 0,
        kind: InvalidKind::BgzfHeader,
    })?;
    Ok((extra_end, bsize))
}

fn sum_uncomp_before(data: &[u8], target: usize) -> Result<u64, FastqError> {
    let mut pos = 0usize;
    let mut total = 0u64;
    while pos < target {
        let (block_size, isize) = parse_block_sizes(data, pos)?;
        total = total
            .checked_add(isize as u64)
            .ok_or(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfVirtualOffset,
            })?;
        pos = pos
            .checked_add(block_size)
            .ok_or(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfVirtualOffset,
            })?;
    }
    if pos != target {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfVirtualOffset,
        });
    }
    Ok(total)
}
