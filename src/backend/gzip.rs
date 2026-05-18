use std::fs::File;
use std::io;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};
use miniz_oxide::inflate::TINFLStatus;
use miniz_oxide::inflate::core::{DecompressorOxide, decompress, inflate_flags};

use crate::error::{FastqError, InvalidKind};
use crate::simd::newline::find_lf;

const DEFAULT_OUT_BUF_SIZE: usize = 4 * 1024 * 1024;
const HISTORY_KEEP: usize = 64 * 1024;
const VALIDATE_GZIP: bool = cfg!(feature = "gzip-validate");

const _: () = {
    assert!(DEFAULT_OUT_BUF_SIZE > HISTORY_KEEP);
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LineStatus {
    Line,
    EofClean,
    /// EOF with non-empty bytes in the destination (missing final `\n`).
    EofPartial,
}

pub struct GzipBackend {
    mmap: Mmap,
    src_pos: usize,
    buf: Vec<u8>,
    buf_pos: usize,
    buf_len: usize,
    logical_offset: u64,
    decomp: DecompressorOxide,
    need_header: bool,
    finished: bool,
    crc32: u32,
    isize: u32,
}

impl GzipBackend {
    pub fn new(path: &Path) -> Result<Self, FastqError> {
        Self::new_with_buf_size(path, DEFAULT_OUT_BUF_SIZE)
    }

    pub fn new_with_buf_size(path: &Path, out_buf_size: usize) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        let metadata = file.metadata()?;
        if metadata.len() == 0 {
            return Err(FastqError::InvalidFormat {
                offset: 0,
                kind: InvalidKind::GzipHeader,
            });
        }
        let out_buf_size = out_buf_size.max(HISTORY_KEEP * 2);
        // SAFETY: file is kept alive for the duration of the mmap; mapping spans the whole file.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(Self {
            mmap,
            src_pos: 0,
            buf: vec![0u8; out_buf_size],
            buf_pos: 0,
            buf_len: 0,
            logical_offset: 0,
            decomp: DecompressorOxide::new(),
            need_header: true,
            finished: false,
            crc32: 0,
            isize: 0,
        })
    }

    #[inline]
    pub fn logical_offset(&self) -> u64 {
        self.logical_offset + self.buf_pos as u64
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

    pub fn skip_line(&mut self) -> Result<(), FastqError> {
        loop {
            let slice = self.available_slice();
            if !slice.is_empty() {
                if let Some(lf) = find_lf(slice, 0) {
                    self.advance(lf + 1);
                    return Ok(());
                }
                let n = slice.len();
                self.advance(n);
            }
            if !self.refill()? {
                return Ok(());
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
        loop {
            if self.finished {
                return Ok(false);
            }
            if self.need_header {
                if self.src_pos >= self.mmap.len() {
                    self.finished = true;
                    return Ok(false);
                }
                self.start_member()?;
                self.need_header = false;
            }

            self.compact_for_inflate();

            if self.buf_len >= self.buf.len() {
                return Err(FastqError::Io(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "inflate buffer full",
                )));
            }

            let in_buf = &self.mmap[self.src_pos..];
            let flags = inflate_flags::TINFL_FLAG_USING_NON_WRAPPING_OUTPUT_BUF;
            let out_start = self.buf_len;
            let (status, in_consumed, out_consumed) =
                decompress(&mut self.decomp, in_buf, &mut self.buf, self.buf_len, flags);
            self.src_pos =
                self.src_pos
                    .checked_add(in_consumed)
                    .ok_or_else(|| FastqError::InvalidFormat {
                        offset: self.logical_offset(),
                        kind: InvalidKind::GzipData,
                    })?;
            self.buf_len = self.buf_len.checked_add(out_consumed).ok_or_else(|| {
                FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipData,
                }
            })?;

            if VALIDATE_GZIP && out_consumed > 0 {
                let chunk = &self.buf[out_start..out_start + out_consumed];
                let mut hasher = crc32fast::Hasher::new_with_initial(self.crc32);
                hasher.update(chunk);
                self.crc32 = hasher.finalize();
                self.isize = self.isize.wrapping_add(out_consumed as u32);
            }

            match status {
                TINFLStatus::Done => {
                    self.finish_member()?;
                    self.decomp = DecompressorOxide::new();
                    self.need_header = true;
                    if self.buf_pos < self.buf_len {
                        return Ok(true);
                    }
                    continue;
                }
                TINFLStatus::HasMoreOutput => return Ok(true),
                TINFLStatus::NeedsMoreInput | TINFLStatus::FailedCannotMakeProgress => {
                    if in_consumed == 0 && out_consumed == 0 {
                        return Err(FastqError::InvalidFormat {
                            offset: self.logical_offset(),
                            kind: InvalidKind::GzipData,
                        });
                    }
                    if self.src_pos >= self.mmap.len() {
                        return Ok(self.buf_len > self.buf_pos);
                    }
                    continue;
                }
                TINFLStatus::Failed | TINFLStatus::BadParam | TINFLStatus::Adler32Mismatch => {
                    return Err(FastqError::InvalidFormat {
                        offset: self.logical_offset(),
                        kind: InvalidKind::GzipData,
                    });
                }
                _ => {
                    return Err(FastqError::InvalidFormat {
                        offset: self.logical_offset(),
                        kind: InvalidKind::GzipData,
                    });
                }
            }
        }
    }

    #[inline]
    fn compact_for_inflate(&mut self) {
        if self.buf_pos <= HISTORY_KEEP {
            return;
        }
        let drop = self.buf_pos - HISTORY_KEEP;
        self.buf.copy_within(drop..self.buf_len, 0);
        self.buf_pos -= drop;
        self.buf_len -= drop;
        self.logical_offset += drop as u64;
    }

    fn start_member(&mut self) -> Result<(), FastqError> {
        let header_end = self.parse_header_at(self.src_pos)?;
        self.src_pos = header_end;
        self.crc32 = 0;
        self.isize = 0;
        Ok(())
    }

    fn finish_member(&mut self) -> Result<(), FastqError> {
        let trailer_end = self
            .src_pos
            .checked_add(8)
            .ok_or_else(|| FastqError::UnexpectedEof {
                offset: self.logical_offset(),
            })?;
        if trailer_end > self.mmap.len() {
            return Err(FastqError::UnexpectedEof {
                offset: self.logical_offset(),
            });
        }
        let trailer = &self.mmap[self.src_pos..trailer_end];
        if VALIDATE_GZIP {
            let expected_crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
            let expected_isize =
                u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
            if self.crc32 != expected_crc {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipTrailerCrc,
                });
            }
            if self.isize != expected_isize {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipTrailerIsize,
                });
            }
        }
        self.src_pos = trailer_end;
        Ok(())
    }

    fn parse_header_at(&self, start: usize) -> Result<usize, FastqError> {
        let bytes = &self.mmap[start..];
        if bytes.len() < 10 {
            return Err(FastqError::InvalidFormat {
                offset: self.logical_offset(),
                kind: InvalidKind::GzipHeader,
            });
        }
        if bytes[0] != 0x1f || bytes[1] != 0x8b || bytes[2] != 8 {
            return Err(FastqError::InvalidFormat {
                offset: self.logical_offset(),
                kind: InvalidKind::GzipHeader,
            });
        }
        let flg = bytes[3];
        let mut i = 10usize;
        if flg & 0x04 != 0 {
            let end = i.checked_add(2).ok_or_else(|| FastqError::InvalidFormat {
                offset: self.logical_offset(),
                kind: InvalidKind::GzipHeader,
            })?;
            if end > bytes.len() {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipHeader,
                });
            }
            let xlen = u16::from_le_bytes([bytes[i], bytes[i + 1]]) as usize;
            i = end;
            let after_extra = i
                .checked_add(xlen)
                .ok_or_else(|| FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipHeader,
                })?;
            if after_extra > bytes.len() {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipHeader,
                });
            }
            i = after_extra;
        }
        if flg & 0x08 != 0 {
            while i < bytes.len() && bytes[i] != 0 {
                i += 1;
            }
            if i >= bytes.len() {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipHeader,
                });
            }
            i += 1;
        }
        if flg & 0x10 != 0 {
            while i < bytes.len() && bytes[i] != 0 {
                i += 1;
            }
            if i >= bytes.len() {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipHeader,
                });
            }
            i += 1;
        }
        if flg & 0x02 != 0 {
            let end = i.checked_add(2).ok_or_else(|| FastqError::InvalidFormat {
                offset: self.logical_offset(),
                kind: InvalidKind::GzipHeader,
            })?;
            if end > bytes.len() {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipHeader,
                });
            }
            i = end;
        }
        Ok(start.saturating_add(i))
    }
}
