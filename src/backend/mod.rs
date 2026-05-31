pub mod bgzf;
pub mod gzip;
pub mod mmap;
pub mod stream;

#[cfg(feature = "noodles-bgzf")]
pub mod noodles_bgzf;

use crate::backend::bgzf::BgzfBackend;
use crate::backend::gzip::GzipBackend;
use crate::backend::mmap::MmapBackend;
use crate::backend::stream::StreamBackend;

#[cfg(feature = "noodles-bgzf")]
use crate::backend::noodles_bgzf::NoodlesBgzfBackend;

// One backend per reader; boxing would just add indirection on the hot path.
#[allow(clippy::large_enum_variant)]
pub enum Backend {
    Plain(MmapBackend),
    Gzip(GzipBackend),
    Bgzf(BgzfBackend),
    /// Arbitrary `BufRead` source. No mmap, no random-access; lines are copied into the
    /// reader's scratch buffer.
    Stream(StreamBackend),
    /// Optional `noodles-bgzf` adapter. Same semantics as [`Backend::Bgzf`] but the inflate
    /// is performed by `noodles_bgzf::Reader` (which can use a thread-pool via
    /// `noodles_bgzf::MultithreadedReader` upstream).
    #[cfg(feature = "noodles-bgzf")]
    NoodlesBgzf(NoodlesBgzfBackend),
}
