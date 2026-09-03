//! Multi-threaded BGZF input.
//!
//! BGZF blocks are independent deflate streams, so inflate parallelises the same way htslib's
//! `-@` option and `noodles_bgzf::MultithreadedReader` do it: one thread walks block headers
//! (which needs no decompression), worker threads inflate blocks, and the consumer collects them
//! in order. Jobs go out round-robin and each worker answers on its own channel, so the results
//! come back in file order without a reordering buffer.
//!
//! The result is a plain [`BufRead`], which is what lets the rest of the crate treat it like any
//! other stream: validation, multi-line parsing and paired reading all work unchanged. `tell`
//! and `seek` do not: the file is being consumed out of order behind the scenes.

use std::fs::File;
use std::io::{self, BufRead, Read};
use std::path::Path;
use std::sync::Arc;
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::thread::JoinHandle;

use memmap2::MmapOptions;

use crate::backend::bgzf::{BGZF_EOF, BlockDecoder, MAX_BLOCK_UNCOMP, parse_block_sizes};
use crate::error::{FastqError, InvalidKind};

/// Blocks in flight per worker. Two is enough to keep a worker busy while the consumer drains
/// the previous block, and caps memory at `threads * 2 * 64 KiB`.
const QUEUE_DEPTH: usize = 2;

type BlockResult = Result<Vec<u8>, FastqError>;

pub(crate) struct ParallelBgzfReader {
    results: Vec<Receiver<BlockResult>>,
    next: usize,
    current: Vec<u8>,
    pos: usize,
    done: bool,
    workers: Vec<JoinHandle<()>>,
    dispatcher: Option<JoinHandle<()>>,
}

impl ParallelBgzfReader {
    pub(crate) fn open(path: &Path, threads: usize, eof_check: bool) -> Result<Self, FastqError> {
        let file = File::open(path)?;
        if file.metadata()?.len() == 0 {
            return Err(FastqError::invalid(0, InvalidKind::BgzfHeader));
        }
        // SAFETY: the mapping spans the whole file; see `MmapBackend::open`.
        let mmap = Arc::new(unsafe { MmapOptions::new().map(&file)? });
        if eof_check
            && !(mmap.len() >= BGZF_EOF.len()
                && mmap[mmap.len() - BGZF_EOF.len()..] == BGZF_EOF[..])
        {
            return Err(FastqError::invalid(0, InvalidKind::BgzfMissingEofMarker));
        }

        let threads = threads.max(1);
        let mut job_txs: Vec<SyncSender<usize>> = Vec::with_capacity(threads);
        let mut results = Vec::with_capacity(threads);
        let mut workers = Vec::with_capacity(threads);

        for id in 0..threads {
            let (job_tx, job_rx) = sync_channel::<usize>(QUEUE_DEPTH);
            let (res_tx, res_rx) = sync_channel::<BlockResult>(QUEUE_DEPTH);
            let mmap = Arc::clone(&mmap);
            let handle = std::thread::Builder::new()
                .name(format!("kira-bgzf-{id}"))
                .spawn(move || {
                    let mut decoder = BlockDecoder::new();
                    let mut scratch = vec![0u8; MAX_BLOCK_UNCOMP];
                    while let Ok(start) = job_rx.recv() {
                        let outcome = match decoder.decode(&mmap, start, &mut scratch, 0) {
                            Ok(block) => Ok(scratch[..block.uncompressed_len].to_vec()),
                            Err(e) => Err(e),
                        };
                        let failed = outcome.is_err();
                        if res_tx.send(outcome).is_err() || failed {
                            break;
                        }
                    }
                })
                .map_err(FastqError::Io)?;
            job_txs.push(job_tx);
            results.push(res_rx);
            workers.push(handle);
        }

        let dispatcher = {
            let mmap = Arc::clone(&mmap);
            std::thread::Builder::new()
                .name("kira-bgzf-dispatch".to_string())
                .spawn(move || {
                    let mut pos = 0usize;
                    let mut seq = 0usize;
                    while pos < mmap.len() {
                        // Header parsing only touches 18 bytes per block, so one thread keeps
                        // every worker fed.
                        let size = match parse_block_sizes(&mmap, pos) {
                            Ok((size, _)) => size,
                            Err(_) => {
                                // Hand the bad offset to a worker so the consumer sees the same
                                // error it would get from the sequential backend, then stop.
                                let _ = job_txs[seq % job_txs.len()].send(pos);
                                return;
                            }
                        };
                        if job_txs[seq % job_txs.len()].send(pos).is_err() {
                            return;
                        }
                        pos += size;
                        seq += 1;
                    }
                })
                .map_err(FastqError::Io)?
        };

        Ok(Self {
            results,
            next: 0,
            current: Vec::new(),
            pos: 0,
            done: false,
            workers,
            dispatcher: Some(dispatcher),
        })
    }
}

impl Drop for ParallelBgzfReader {
    fn drop(&mut self) {
        // Dropping the receivers makes every blocked `send` fail, which unwinds the workers and
        // then the dispatcher.
        self.results.clear();
        if let Some(handle) = self.dispatcher.take() {
            let _ = handle.join();
        }
        for handle in self.workers.drain(..) {
            let _ = handle.join();
        }
    }
}

impl BufRead for ParallelBgzfReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        while self.pos >= self.current.len() {
            if self.done || self.results.is_empty() {
                return Ok(&[]);
            }
            match self.results[self.next].recv() {
                Ok(Ok(block)) => {
                    self.next = (self.next + 1) % self.results.len();
                    self.current = block;
                    self.pos = 0;
                }
                Ok(Err(e)) => {
                    self.done = true;
                    return Err(io::Error::other(e));
                }
                // The worker that owed us the next block is gone: end of file.
                Err(_) => {
                    self.done = true;
                    return Ok(&[]);
                }
            }
        }
        Ok(&self.current[self.pos..])
    }

    fn consume(&mut self, amt: usize) {
        self.pos = (self.pos + amt).min(self.current.len());
    }
}

impl Read for ParallelBgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf()?;
        let n = available.len().min(buf.len());
        buf[..n].copy_from_slice(&available[..n]);
        self.consume(n);
        Ok(n)
    }
}
