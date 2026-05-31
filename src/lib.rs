//! # kira-fastq
//!
//! High-performance FASTQ I/O. mmap-first reading, single-syscall-per-record writing,
//! SIMD validation on every hot path, paired-end and multi-line support, optional async,
//! optional `noodles-bgzf` interop.
//!
//! See the crate `README.md` for the full feature matrix and a perf cheat-sheet.

pub mod backend;
pub mod error;
pub mod format;
mod multiline;
pub mod offset;
pub mod paired;
pub mod parser;
pub mod reader;
pub mod record;
pub mod simd;
pub mod validation;
pub mod writer;

#[cfg(feature = "async")]
#[path = "async_io/mod.rs"]
pub mod r#async;

pub use crate::error::{FastqError, InvalidKind, PairedWhich, UnsupportedOperation};
pub use crate::format::FastqFormat;
pub use crate::offset::VirtualOffset;
pub use crate::paired::{PairedFastqReader, canonical_read_id};
pub use crate::reader::{FastqReader, TryForEachError};
pub use crate::record::{FastqRecord, FastqRecordOwned};
pub use crate::validation::{Alphabet, ValidationMode};
pub use crate::writer::{BoxedWriter, FastqWriter, WriteValidation};

#[cfg(feature = "async")]
pub use crate::r#async::{AnyAsyncReader, AsyncFastqReader, AsyncFastqWriter, RecordStream};
