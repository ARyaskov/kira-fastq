pub mod bgzf;
pub mod gzip;
pub mod mmap;

use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::GzipBackend;
use crate::backend::mmap::MmapBackend;

// One backend per reader; boxing would just add indirection on the hot path.
#[allow(clippy::large_enum_variant)]
pub enum Backend {
    Plain(MmapBackend),
    Gzip(GzipBackend),
    Bgzf(BgzfBackend),
}
