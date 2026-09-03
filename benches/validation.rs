//! SIMD validation kernels in isolation.

mod common;

use std::hint::black_box;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use kira_fastq::simd::bases::validate_bases_with;
use kira_fastq::simd::qual::{validate_qual, validate_qual_encoding};
use kira_fastq::{Alphabet, QualityEncoding};

fn bench_validation(c: &mut Criterion) {
    const LEN: usize = 32 * 1024 * 1024;
    let seq: Vec<u8> = (0..LEN).map(|i| b"ACGT"[i % 4]).collect();
    let qual: Vec<u8> = (0..LEN).map(|i| 33 + (i % 40) as u8).collect();
    // Phred+64 data has to be validated against Phred+64 bytes, or the kernel exits on the first
    // byte and the benchmark measures nothing.
    let qual64: Vec<u8> = (0..LEN).map(|i| 64 + (i % 40) as u8).collect();

    let mut group = c.benchmark_group("validation");
    group.throughput(Throughput::Bytes(LEN as u64));
    group.sample_size(20);

    for (name, alphabet) in [
        ("bases_acgtn", Alphabet::AcgtnStrict),
        ("bases_acgtn_case", Alphabet::AcgtnCase),
        ("bases_iupac", Alphabet::Iupac),
    ] {
        group.bench_function(name, |b| {
            b.iter(|| black_box(validate_bases_with(black_box(&seq), alphabet)))
        });
    }
    group.bench_function("qual_phred33", |b| {
        b.iter(|| black_box(validate_qual(black_box(&qual))))
    });
    group.bench_function("qual_phred64", |b| {
        b.iter(|| {
            black_box(validate_qual_encoding(
                black_box(&qual64),
                QualityEncoding::PHRED64,
            ))
        })
    });
    group.finish();
}

criterion_group!(benches, bench_validation);
criterion_main!(benches);
