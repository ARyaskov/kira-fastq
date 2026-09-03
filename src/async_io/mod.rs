//! Async FASTQ reader and writer, on tokio.
//!
//! There is no memory mapping here: mapping is synchronous by nature, so async I/O goes through
//! `tokio::io`. Records still borrow out of the reader's scratch buffer, so per-record parsing
//! allocates nothing.
//!
//! **Reach for the sync API first.** It is faster and simpler for a tool that reads a file and
//! exits. Use this when a tokio runtime is already in the picture: long-lived servers,
//! backpressure-aware pipelines, integration with axum, tonic or sqlx.
//!
//! ## Reader
//!
//! - [`AsyncFastqReader::from_reader`] over any `AsyncBufRead`.
//! - [`AnyAsyncReader::from_path`] for a file, with gzip detected from the magic bytes and
//!   multi-member streams decoded in full.
//! - [`AsyncFastqReader::next`] yields a borrowed record; it is **not** cancel-safe, see the
//!   type's documentation.
//! - [`AsyncFastqReader::records`] turns the reader into a stream of owned records.
//!
//! ## Writer
//!
//! - [`AsyncFastqWriter::from_writer`] over any `AsyncWrite`, [`AsyncFastqWriter::from_path`]
//!   for a file with `.gz` support.
//! - One `write_all` per record, as on the sync side.
//! - [`AsyncFastqWriter::shutdown`] writes the gzip trailer; call it before dropping.
//!
//! BGZF and zstd output are sync-only. For BGZF input with real virtual offsets, use the sync
//! reader or wrap a `noodles_bgzf` async reader with [`AsyncFastqReader::from_reader`].

mod reader;
mod writer;

pub use self::reader::{AnyAsyncReader, AsyncFastqReader, RecordStream};
pub use self::writer::{AsyncFastqWriter, BoxedAsyncWriter};
