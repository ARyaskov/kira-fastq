use std::fs::File;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};
use miniz_oxide::inflate::TINFLStatus;
use miniz_oxide::inflate::core::{DecompressorOxide, decompress, inflate_flags};

use crate::error::{FastqError, InvalidKind};
use crate::offset::VirtualOffset;

const OUT_BUF_SIZE: usize = 256 * 1024;

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
}

impl BgzfBackend {
    pub fn new(path: &Path) -> Result<Self, FastqError> {
        let file = File::open(path)?;
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
        })
    }

    #[inline]
    pub fn available_slice(&self) -> &[u8] {
        if self.buf_pos >= self.buf_len {
            &[]
        } else {
            &self.buf[self.buf_pos..self.buf_len]
        }
    }

    #[inline]
    pub fn advance(&mut self, n: usize) {
        self.buf_pos += n;
    }

    #[inline]
    pub fn logical_offset(&self) -> u64 {
        self.total_uncomp_before + self.buf_pos as u64
    }

    #[inline]
    pub fn tell(&self) -> VirtualOffset {
        VirtualOffset(((self.block_start as u64) << 16) | (self.buf_pos as u64))
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
        let block_end = block_off + bsize;
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

    pub fn refill(&mut self) -> Result<bool, FastqError> {
        if self.buf_pos < self.buf_len {
            return Ok(true);
        }
        if self.finished {
            return Ok(false);
        }
        self.total_uncomp_before += self.block_uncomp_len as u64;
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
        let block_end = self.src_pos + block_size;
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
        let flags = inflate_flags::TINFL_FLAG_USING_NON_WRAPPING_OUTPUT_BUF;
        let mut decomp = DecompressorOxide::new();
        let (status, _in_consumed, out_consumed) =
            decompress(&mut decomp, in_buf, &mut self.buf, 0, flags);
        match status {
            TINFLStatus::Done => {}
            TINFLStatus::HasMoreOutput => {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::BgzfBlock,
                });
            }
            _ => {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::BgzfBlock,
                });
            }
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
    let (header_end, bsize) = parse_bgzf_header(data, start)?;
    let block_size = bsize as usize + 1;
    let trailer_pos = start + block_size - 8;
    if trailer_pos + 8 > data.len() {
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
    let _ = header_end;
    Ok((block_size, isize))
}

fn parse_deflate_range(data: &[u8], start: usize) -> Result<(usize, usize), FastqError> {
    let (header_end, bsize) = parse_bgzf_header(data, start)?;
    let block_size = bsize as usize + 1;
    let deflate_start = header_end;
    let deflate_end = start + block_size - 8;
    Ok((deflate_start, deflate_end))
}

fn parse_bgzf_header(data: &[u8], start: usize) -> Result<(usize, u16), FastqError> {
    if start + 10 > data.len() {
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
    let mut i = start + 10;
    if i + 2 > data.len() {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfHeader,
        });
    }
    let xlen = u16::from_le_bytes([data[i], data[i + 1]]) as usize;
    i += 2;
    if i + xlen > data.len() {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfHeader,
        });
    }
    let extra_end = i + xlen;
    let mut bsize = None;
    while i + 4 <= extra_end {
        let si1 = data[i];
        let si2 = data[i + 1];
        let slen = u16::from_le_bytes([data[i + 2], data[i + 3]]) as usize;
        i += 4;
        if i + slen > extra_end {
            return Err(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfHeader,
            });
        }
        if si1 == b'B' && si2 == b'C' && slen == 2 {
            bsize = Some(u16::from_le_bytes([data[i], data[i + 1]]));
        }
        i += slen;
    }
    let bsize = match bsize {
        Some(v) => v,
        None => {
            return Err(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::BgzfHeader,
            });
        }
    };
    Ok((extra_end, bsize))
}

fn sum_uncomp_before(data: &[u8], target: usize) -> Result<u64, FastqError> {
    let mut pos = 0usize;
    let mut total = 0u64;
    while pos < target {
        let (block_size, isize) = parse_block_sizes(data, pos)?;
        total = total.wrapping_add(isize as u64);
        pos += block_size;
    }
    if pos != target {
        return Err(FastqError::InvalidFormat {
            offset: 0,
            kind: InvalidKind::BgzfVirtualOffset,
        });
    }
    Ok(total)
}
