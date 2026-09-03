//! Every input path over the same data, so the numbers are comparable.

mod common;

use std::hint::black_box;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use kira_fastq::FastqReader;

fn bench_end_to_end(c: &mut Criterion) {
    let corpus = common::corpus(150, common::DEFAULT_BYTES);
    let mut group = c.benchmark_group("end_to_end");
    // Throughput is reported against the uncompressed size for every backend, so the numbers
    // say how fast records come out, not how fast bytes leave the disk.
    group.throughput(Throughput::Bytes(corpus.bytes as u64));
    group.sample_size(10);

    group.bench_function("plain_mmap", |b| {
        b.iter(|| drain(FastqReader::from_path(&corpus.plain).expect("open")))
    });
    group.bench_function("plain_buffered", |b| {
        b.iter(|| drain(FastqReader::from_path_buffered(&corpus.plain).expect("open")))
    });
    group.bench_function("gzip", |b| {
        b.iter(|| drain(FastqReader::from_path(&corpus.gzip).expect("open")))
    });
    group.bench_function("bgzf", |b| {
        b.iter(|| drain(FastqReader::from_path(&corpus.bgzf).expect("open")))
    });
    for threads in [2usize, 4, 8] {
        group.bench_function(format!("bgzf_parallel_{threads}"), |b| {
            b.iter(|| {
                drain(FastqReader::from_bgzf_path_parallel(&corpus.bgzf, threads).expect("open"))
            })
        });
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

criterion_group!(benches, bench_end_to_end);
criterion_main!(benches);
