//! Async FASTQ reader and writer (`tokio`-based).
//!
//! **No mmap on this path.** mmap is inherently synchronous; async I/O implies streaming
//! reads/writes via `tokio::io`. The zero-copy borrowed-record API survives (records borrow
//! into the reader's scratch buffer), but the underlying source/sink goes through tokio's
//! syscall-backed I/O. Expect ~10–30 % lower throughput than the sync mmap path for plain
//! files; gzip/BGZF are typically I/O-bound and the difference shrinks to near zero.
//!
//! ## Reader
//!
//! - [`AsyncFastqReader::from_reader`] — any `R: AsyncBufRead + Unpin + Send`.
//! - [`AsyncFastqReader::from_path`] — `tokio::fs::File`, with auto-detect of gzip
//!   compression via the file's magic bytes (BGZF is treated as gzip in async mode; for
//!   true BGZF semantics use the sync path or a `noodles_bgzf::AsyncReader` wrapped via
//!   `from_reader`).
//! - [`AsyncFastqReader::next`] — borrowed record (`FastqRecord<'_>`).
//! - [`AsyncFastqReader::records`] — owned-record [`Stream`].
//!
//! ## Writer
//!
//! - [`AsyncFastqWriter::from_writer`] — any `W: AsyncWrite + Unpin + Send`.
//! - [`AsyncFastqWriter::from_path`] — `tokio::fs::File`, with `.gz` triggering streaming
//!   gzip via `async-compression`.
//! - Single `write_all().await` per record via the same assemble-into-scratch trick as
//!   the sync writer.

mod reader;
mod writer;

pub use self::reader::{AnyAsyncReader, AsyncFastqReader, RecordStream};
pub use self::writer::{AsyncFastqWriter, BoxedAsyncWriter};
