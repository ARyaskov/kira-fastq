//! Mapped plain input against compressed input, at several read lengths.

mod common;

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use kira_fastq::FastqReader;

fn bench_mmap_vs_gzip(c: &mut Criterion) {
    let mut group = c.benchmark_group("mmap_vs_gzip");
    group.sample_size(10);
    // 150 bp is Illumina, 20 kb approximates a long-read run: per-record overhead dominates in
    // the first case, payload movement in the second.
    for read_len in [150usize, 20_000] {
        let corpus = common::corpus(read_len, common::DEFAULT_BYTES);
        group.throughput(Throughput::Bytes(corpus.bytes as u64));
        group.bench_with_input(
            BenchmarkId::new("plain_mmap", read_len),
            &corpus.plain,
            |b, path| b.iter(|| drain(FastqReader::from_path(path).expect("open"))),
        );
        group.bench_with_input(
            BenchmarkId::new("gzip", read_len),
            &corpus.gzip,
            |b, path| b.iter(|| drain(FastqReader::from_path(path).expect("open"))),
        );
        group.bench_with_input(
            BenchmarkId::new("bgzf", read_len),
            &corpus.bgzf,
            |b, path| b.iter(|| drain(FastqReader::from_path(path).expect("open"))),
        );
    }
    group.finish();
}

fn drain(mut reader: FastqReader) -> usize {
    let mut bases = 0usize;
    while let Some(rec) = reader.next().expect("read") {
        bases += rec.len();
        black_box(rec.seq());
    }
    black_box(bases)
}

criterion_group!(benches, bench_mmap_vs_gzip);
criterion_main!(benches);
