//! BGZF output.
//!
//! Blocks are framed exactly as `bgzip` writes them: a gzip member per block carrying the `BC`
//! extra subfield with the block size, capped at 65280 bytes of input so the framed block always
//! fits the 16-bit size field, and the canonical empty block as the end-of-file marker. Files
//! written here are readable by htslib, `bgzip -d`, and this crate's own BGZF reader, and their
//! virtual offsets mean the same thing everywhere.
//!
//! Two implementations: [`BgzfWriter`] compresses on the calling thread, and
//! [`ParallelBgzfWriter`] hands blocks to a pool and writes the results in order, which is how
//! `samtools -@` scales compression.

use std::collections::VecDeque;
use std::io::{self, Write};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::thread::JoinHandle;

use flate2::{Compress, Compression, FlushCompress, Status};

use crate::backend::bgzf::BGZF_EOF;
use crate::error::{FastqError, UnsupportedOperation};

/// Uncompressed bytes per block. `bgzip` uses the same value; the slack under 64 KiB leaves room
/// for the header, trailer and worst-case deflate expansion.
pub(crate) const BLOCK_SIZE: usize = 65280;

/// Check a deflate level before handing it to a codec that would panic on an invalid one.
pub(crate) fn checked_level(level: u32) -> Result<Compression, FastqError> {
    if level > 9 {
        return Err(FastqError::Unsupported(
            UnsupportedOperation::CompressionLevel,
        ));
    }
    Ok(Compression::new(level))
}

/// Raw-deflate `input` into `out`, replacing its contents.
fn deflate_raw(comp: &mut Compress, input: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    comp.reset();
    out.clear();
    out.reserve(input.len() / 2 + 128);
    let mut consumed = 0usize;
    loop {
        let before_in = comp.total_in();
        let status = comp
            .compress_vec(&input[consumed..], out, FlushCompress::Finish)
            .map_err(io::Error::other)?;
        consumed += (comp.total_in() - before_in) as usize;
        match status {
            Status::StreamEnd => return Ok(()),
            // Ran out of spare capacity: grow and continue.
            Status::Ok | Status::BufError => out.reserve(out.len().max(1024)),
        }
    }
}

/// Frame one already-deflated payload as a BGZF block.
fn frame_block(
    payload: &[u8],
    crc: u32,
    uncompressed_len: usize,
    out: &mut Vec<u8>,
) -> io::Result<()> {
    let block_size = 18 + payload.len() + 8;
    if block_size > u16::MAX as usize + 1 {
        return Err(io::Error::other("BGZF block exceeds 64 KiB"));
    }
    let bsize = (block_size - 1) as u16;
    out.clear();
    out.reserve(block_size);
    out.extend_from_slice(&[0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff]);
    out.extend_from_slice(&[0x06, 0x00, b'B', b'C', 0x02, 0x00]);
    out.extend_from_slice(&bsize.to_le_bytes());
    out.extend_from_slice(payload);
    out.extend_from_slice(&crc.to_le_bytes());
    out.extend_from_slice(&(uncompressed_len as u32).to_le_bytes());
    Ok(())
}

/// Compress one block's worth of data into a complete framed BGZF block.
fn compress_block(comp: &mut Compress, data: &[u8]) -> io::Result<Vec<u8>> {
    let mut payload = Vec::new();
    deflate_raw(comp, data, &mut payload)?;
    let mut hasher = crc32fast::Hasher::new();
    hasher.update(data);
    let mut block = Vec::new();
    frame_block(&payload, hasher.finalize(), data.len(), &mut block)?;
    Ok(block)
}

/// Single-threaded BGZF writer.
pub struct BgzfWriter<W: Write> {
    inner: Option<W>,
    staging: Vec<u8>,
    payload: Vec<u8>,
    block: Vec<u8>,
    comp: Compress,
    finished: bool,
}

impl<W: Write> BgzfWriter<W> {
    pub fn new(inner: W, level: u32) -> Result<Self, FastqError> {
        let level = checked_level(level)?;
        Ok(Self {
            inner: Some(inner),
            staging: Vec::with_capacity(BLOCK_SIZE),
            payload: Vec::with_capacity(BLOCK_SIZE),
            block: Vec::with_capacity(BLOCK_SIZE + 64),
            comp: Compress::new(level, false),
            finished: false,
        })
    }

    fn write_staged(&mut self) -> io::Result<()> {
        if self.staging.is_empty() {
            return Ok(());
        }
        deflate_raw(&mut self.comp, &self.staging, &mut self.payload)?;
        let mut hasher = crc32fast::Hasher::new();
        hasher.update(&self.staging);
        frame_block(
            &self.payload,
            hasher.finalize(),
            self.staging.len(),
            &mut self.block,
        )?;
        if let Some(inner) = self.inner.as_mut() {
            inner.write_all(&self.block)?;
        }
        self.staging.clear();
        Ok(())
    }

    /// Flush the pending block, append the end-of-file marker, and return the sink.
    ///
    /// Call this rather than relying on `Drop` when write errors matter: the `Drop` path has
    /// nowhere to report them.
    pub fn finish(mut self) -> Result<W, FastqError> {
        self.finalize()?;
        Ok(self.inner.take().expect("writer already finished"))
    }

    fn finalize(&mut self) -> Result<(), FastqError> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;
        self.write_staged().map_err(FastqError::Io)?;
        if let Some(inner) = self.inner.as_mut() {
            inner.write_all(&BGZF_EOF).map_err(FastqError::Io)?;
            inner.flush().map_err(FastqError::Io)?;
        }
        Ok(())
    }
}

impl<W: Write> Write for BgzfWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let space = BLOCK_SIZE - self.staging.len();
        let n = space.min(buf.len());
        self.staging.extend_from_slice(&buf[..n]);
        if self.staging.len() == BLOCK_SIZE {
            self.write_staged()?;
        }
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.write_staged()?;
        match self.inner.as_mut() {
            Some(inner) => inner.flush(),
            None => Ok(()),
        }
    }
}

impl<W: Write> Drop for BgzfWriter<W> {
    fn drop(&mut self) {
        // A BGZF file without its end-of-file marker reads as truncated, so write it even on the
        // drop path; errors can only be reported by `finish`.
        let _ = self.finalize();
    }
}

type BlockJob = (Vec<u8>, SyncSender<io::Result<Vec<u8>>>);

/// BGZF writer that compresses blocks on a thread pool.
///
/// Blocks are handed out to workers as they fill up and collected in submission order, so the
/// output is byte-identical to what [`BgzfWriter`] would produce at the same level.
pub struct ParallelBgzfWriter<W: Write> {
    inner: Option<W>,
    staging: Vec<u8>,
    jobs: Option<SyncSender<BlockJob>>,
    pending: VecDeque<Receiver<io::Result<Vec<u8>>>>,
    workers: Vec<JoinHandle<()>>,
    max_pending: usize,
    finished: bool,
}

impl<W: Write> ParallelBgzfWriter<W> {
    pub fn new(inner: W, level: u32, threads: usize) -> Result<Self, FastqError> {
        let level = checked_level(level)?;
        let threads = threads.max(1);
        let (job_tx, job_rx) = sync_channel::<BlockJob>(threads * 2);
        let job_rx = std::sync::Arc::new(std::sync::Mutex::new(job_rx));
        let mut workers = Vec::with_capacity(threads);
        for id in 0..threads {
            let job_rx = std::sync::Arc::clone(&job_rx);
            let handle = std::thread::Builder::new()
                .name(format!("kira-bgzf-w{id}"))
                .spawn(move || {
                    let mut comp = Compress::new(level, false);
                    loop {
                        let job = {
                            let guard = match job_rx.lock() {
                                Ok(guard) => guard,
                                Err(_) => return,
                            };
                            guard.recv()
                        };
                        let Ok((data, reply)) = job else { return };
                        let _ = reply.send(compress_block(&mut comp, &data));
                    }
                })
                .map_err(FastqError::Io)?;
            workers.push(handle);
        }
        Ok(Self {
            inner: Some(inner),
            staging: Vec::with_capacity(BLOCK_SIZE),
            jobs: Some(job_tx),
            pending: VecDeque::new(),
            workers,
            max_pending: threads * 2,
            finished: false,
        })
    }

    fn submit_staged(&mut self) -> io::Result<()> {
        if self.staging.is_empty() {
            return Ok(());
        }
        let data = std::mem::replace(&mut self.staging, Vec::with_capacity(BLOCK_SIZE));
        let (tx, rx) = sync_channel(1);
        if let Some(jobs) = self.jobs.as_ref() {
            jobs.send((data, tx))
                .map_err(|_| io::Error::other("BGZF compression workers stopped"))?;
        }
        self.pending.push_back(rx);
        while self.pending.len() > self.max_pending {
            self.collect_one()?;
        }
        Ok(())
    }

    fn collect_one(&mut self) -> io::Result<()> {
        let Some(rx) = self.pending.pop_front() else {
            return Ok(());
        };
        let block = rx
            .recv()
            .map_err(|_| io::Error::other("BGZF compression worker died"))??;
        if let Some(inner) = self.inner.as_mut() {
            inner.write_all(&block)?;
        }
        Ok(())
    }

    fn drain(&mut self) -> io::Result<()> {
        while !self.pending.is_empty() {
            self.collect_one()?;
        }
        Ok(())
    }

    /// Flush every pending block, append the end-of-file marker, stop the pool, and return the
    /// sink.
    pub fn finish(mut self) -> Result<W, FastqError> {
        self.finalize()?;
        Ok(self.inner.take().expect("writer already finished"))
    }

    fn finalize(&mut self) -> Result<(), FastqError> {
        if self.finished {
            return Ok(());
        }
        self.finished = true;
        self.submit_staged().map_err(FastqError::Io)?;
        self.drain().map_err(FastqError::Io)?;
        self.jobs = None;
        for handle in self.workers.drain(..) {
            let _ = handle.join();
        }
        if let Some(inner) = self.inner.as_mut() {
            inner.write_all(&BGZF_EOF).map_err(FastqError::Io)?;
            inner.flush().map_err(FastqError::Io)?;
        }
        Ok(())
    }
}

impl<W: Write> Write for ParallelBgzfWriter<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        let space = BLOCK_SIZE - self.staging.len();
        let n = space.min(buf.len());
        self.staging.extend_from_slice(&buf[..n]);
        if self.staging.len() == BLOCK_SIZE {
            self.submit_staged()?;
        }
        Ok(n)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.submit_staged()?;
        self.drain()?;
        match self.inner.as_mut() {
            Some(inner) => inner.flush(),
            None => Ok(()),
        }
    }
}

impl<W: Write> Drop for ParallelBgzfWriter<W> {
    fn drop(&mut self) {
        let _ = self.finalize();
    }
}
