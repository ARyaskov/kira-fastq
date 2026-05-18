pub mod bases;
pub mod newline;
pub mod qual;

use std::sync::OnceLock;

#[derive(Debug, Clone, Copy, Default)]
#[allow(dead_code)]
pub(crate) struct CpuFeatures {
    pub avx512bw: bool,
    pub avx2: bool,
    pub sse2: bool,
    pub neon: bool,
}

#[inline]
pub(crate) fn cpu_features() -> CpuFeatures {
    static CACHE: OnceLock<CpuFeatures> = OnceLock::new();
    *CACHE.get_or_init(detect_cpu_features)
}

#[inline]
fn detect_cpu_features() -> CpuFeatures {
    #[allow(unused_mut)]
    let mut f = CpuFeatures::default();

    #[cfg(target_arch = "x86_64")]
    {
        f.avx512bw = std::arch::is_x86_feature_detected!("avx512bw");
        f.avx2 = std::arch::is_x86_feature_detected!("avx2");
        f.sse2 = std::arch::is_x86_feature_detected!("sse2");
    }

    #[cfg(target_arch = "aarch64")]
    {
        f.neon = std::arch::is_aarch64_feature_detected!("neon");
    }

    f
}
