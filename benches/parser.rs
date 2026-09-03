//! Record parsing throughput on a mapped plain file.

mod common;

use std::hint::black_box;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use kira_fastq::{Alphabet, FastqReader, ValidationMode};

fn bench_parser(c: &mut Criterion) {
    let corpus = common::corpus(150, common::DEFAULT_BYTES);
    let mut group = c.benchmark_group("parse_plain");
    group.throughput(Throughput::Bytes(corpus.bytes as u64));
    group.sample_size(20);

    group.bench_function("no_validation", |b| {
        b.iter(|| drain(FastqReader::from_path(&corpus.plain).expect("open")))
    });
    group.bench_function("validate_bases", |b| {
        b.iter(|| {
            drain(
                FastqReader::from_path(&corpus.plain)
                    .expect("open")
                    .with_validation(ValidationMode::Bases)
                    .with_alphabet(Alphabet::Iupac),
            )
        })
    });
    group.bench_function("validate_bases_and_qualities", |b| {
        b.iter(|| {
            drain(
                FastqReader::from_path(&corpus.plain)
                    .expect("open")
                    .with_validation(ValidationMode::BasesAndQualities),
            )
        })
    });
    group.finish();
}

fn drain(mut reader: FastqReader) -> (usize, usize) {
    let mut records = 0usize;
    let mut bases = 0usize;
    while let Some(rec) = reader.next().expect("read") {
        records += 1;
        bases += rec.len();
        black_box(rec.seq());
        black_box(rec.qual());
    }
    (black_box(records), black_box(bases))
}

criterion_group!(benches, bench_parser);
criterion_main!(benches);
