use std::fs::File;
use std::io;
use std::path::Path;

use memmap2::{Mmap, MmapOptions};
use miniz_oxide::inflate::TINFLStatus;
use miniz_oxide::inflate::core::{DecompressorOxide, decompress, inflate_flags};

use crate::error::{FastqError, InvalidKind};

const OUT_BUF_SIZE: usize = 4 * 1024 * 1024;
const HISTORY_KEEP: usize = 64 * 1024;
const VALIDATE_GZIP: bool = cfg!(feature = "gzip-validate");

const _: () = {
    assert!(OUT_BUF_SIZE > HISTORY_KEEP);
};

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
        let file = File::open(path)?;
        // SAFETY: file is kept alive for the duration of the mmap; mapping spans the whole file.
        let mmap = unsafe { MmapOptions::new().map(&file)? };
        Ok(Self {
            mmap,
            src_pos: 0,
            buf: vec![0u8; OUT_BUF_SIZE],
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
        self.logical_offset + self.buf_pos as u64
    }

    pub fn refill(&mut self) -> Result<bool, FastqError> {
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
            self.src_pos += in_consumed;
            self.buf_len += out_consumed;

            if VALIDATE_GZIP && out_consumed > 0 {
                let chunk = &self.buf[out_start..out_start + out_consumed];
                self.crc32 = crc32_update(self.crc32, chunk);
                self.isize = self.isize.wrapping_add(out_consumed as u32);
            }

            match status {
                TINFLStatus::Done => {
                    self.finish_member()?;
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
        self.decomp = DecompressorOxide::new();
        self.crc32 = 0;
        self.isize = 0;
        Ok(())
    }

    fn finish_member(&mut self) -> Result<(), FastqError> {
        if self.src_pos + 8 > self.mmap.len() {
            return Err(FastqError::UnexpectedEof {
                offset: self.logical_offset(),
            });
        }
        let trailer = &self.mmap[self.src_pos..self.src_pos + 8];
        if VALIDATE_GZIP {
            let expected_crc = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
            let expected_isize =
                u32::from_le_bytes([trailer[4], trailer[5], trailer[6], trailer[7]]);
            let actual_crc = crc32_finalize(self.crc32);
            if actual_crc != expected_crc {
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
        self.src_pos += 8;
        Ok(())
    }

    fn parse_header_at(&self, start: usize) -> Result<usize, FastqError> {
        let bytes = &self.mmap[start..];
        let mut i = 10usize;
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
        if flg & 0x04 != 0 {
            if i + 2 > bytes.len() {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipHeader,
                });
            }
            let xlen = u16::from_le_bytes([bytes[i], bytes[i + 1]]) as usize;
            i += 2;
            if i + xlen > bytes.len() {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipHeader,
                });
            }
            i += xlen;
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
            if i + 2 > bytes.len() {
                return Err(FastqError::InvalidFormat {
                    offset: self.logical_offset(),
                    kind: InvalidKind::GzipHeader,
                });
            }
            i += 2;
        }
        Ok(start + i)
    }
}

fn crc32_update(mut crc: u32, buf: &[u8]) -> u32 {
    let mut c = !crc;
    for &b in buf {
        let idx = ((c ^ b as u32) & 0xFF) as usize;
        c = CRC32_TABLE[idx] ^ (c >> 8);
    }
    crc = !c;
    crc
}

fn crc32_finalize(crc: u32) -> u32 {
    crc
}

const CRC32_TABLE: [u32; 256] = make_crc32_table();

const fn make_crc32_table() -> [u32; 256] {
    let mut table = [0u32; 256];
    let mut i = 0usize;
    while i < 256 {
        let mut c = i as u32;
        let mut k = 0;
        while k < 8 {
            if (c & 1) != 0 {
                c = 0xEDB88320 ^ (c >> 1);
            } else {
                c >>= 1;
            }
            k += 1;
        }
        table[i] = c;
        i += 1;
    }
    table
}
