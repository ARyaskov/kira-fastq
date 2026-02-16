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

pub use crate::error::{FastqError, InvalidKind, PairedWhich, UnsupportedOperation};
pub use crate::format::FastqFormat;
pub use crate::offset::VirtualOffset;
pub use crate::paired::PairedFastqReader;
pub use crate::reader::{FastqReader, RecordsIter};
pub use crate::record::FastqRecord;
pub use crate::validation::ValidationMode;
