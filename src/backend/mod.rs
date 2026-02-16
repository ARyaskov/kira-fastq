pub mod bgzf;
pub mod gzip;
pub mod mmap;

use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::GzipBackend;
use crate::backend::mmap::MmapBackend;

pub enum Backend {
    Plain(MmapBackend),
    Gzip(GzipBackend),
    Bgzf(BgzfBackend),
}
