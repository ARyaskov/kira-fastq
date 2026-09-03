//! # kira-fastq
//!
//! FASTQ input and output for Rust: plain, gzip and BGZF, single- and multi-line records,
//! paired-end and interleaved reading, optional async, optional `noodles-bgzf` interop.
//!
//! ```no_run
//! use kira_fastq::FastqReader;
//!
//! # fn main() -> Result<(), kira_fastq::FastqError> {
//! let mut reader = FastqReader::from_path("reads.fastq.gz")?;
//! let mut bases = 0u64;
//! while let Some(rec) = reader.next()? {
//!     bases += rec.len() as u64;
//! }
//! println!("bases={bases}");
//! # Ok(()) }
//! ```
//!
//! ## What it accepts
//!
//! Reading follows what `kseq.h`, seqtk and seqkit accept: a missing newline on the final
//! record, blank lines between records, zero-length reads from adapter trimmers, and CRLF. The
//! compression backend comes from the file's magic bytes rather than its name, because `bgzip`
//! writes BGZF into `.gz` files. Corrupt input is reported, never silently truncated: gzip CRC
//! and ISIZE, BGZF per-block CRC, and the BGZF end-of-file marker are all checked.
//!
//! See `README.md` for the feature matrix and measured throughput.

/// The README's examples are compiled as doctests, so they cannot drift from the API.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeExamples;

pub mod error;
pub mod format;
pub mod offset;
pub mod paired;
pub mod reader;
pub mod record;
pub mod simd;
pub mod validation;
pub mod writer;

pub(crate) mod backend;
mod multiline;
mod parser;

#[cfg(feature = "async")]
#[path = "async_io/mod.rs"]
pub mod r#async;

pub use crate::error::{FastqError, InvalidKind, PairedWhich, UnsupportedOperation};
pub use crate::format::FastqFormat;
pub use crate::offset::VirtualOffset;
pub use crate::paired::{InterleavedFastqReader, PairedFastqReader, canonical_read_id};
pub use crate::reader::{FastqReader, TryForEachError};
pub use crate::record::{FastqRecord, FastqRecordOwned};
pub use crate::validation::{Alphabet, QualityEncoding, ValidationMode, guess_quality_encoding};
pub use crate::writer::{
    BgzfWriter, BoxedWriter, FastqSink, FastqWriter, ParallelBgzfWriter, PathWriter,
    WriteValidation,
};

#[cfg(feature = "async")]
pub use crate::r#async::{AnyAsyncReader, AsyncFastqReader, AsyncFastqWriter, RecordStream};
