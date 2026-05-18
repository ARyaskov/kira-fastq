use std::hint::black_box;
use std::path::Path;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use kira_fastq::FastqReader;

fn bench_mmap_vs_gzip(c: &mut Criterion) {
    let plain = Path::new("benches/data/SRR014966.fastq");
    let gzip = Path::new("benches/data/SRR014966.fastq.gz");
    if !plain.exists() || !gzip.exists() {
        eprintln!("mmap_vs_gzip: skipping — supply benches/data/SRR014966.fastq{{,.gz}} to enable");
        return;
    }
    let plain_size = std::fs::metadata(plain).expect("stat").len();
    let mut group = c.benchmark_group("mmap_vs_gzip");

    group.throughput(Throughput::Bytes(plain_size));
    group.bench_function("plain_mmap", |b| {
        b.iter(|| {
            let mut reader = FastqReader::from_path(plain).expect("open");
            let mut records = 0usize;
            let mut bases = 0usize;
            while let Some(rec) = reader.next().expect("read") {
                records += 1;
                bases += rec.len();
                black_box(rec.seq());
            }
            black_box(records);
            black_box(bases);
        })
    });

    group.throughput(Throughput::Bytes(plain_size));
    group.bench_function("gzip_mmap_inflate", |b| {
        b.iter(|| {
            let mut reader = FastqReader::from_path(gzip).expect("open");
            let mut records = 0usize;
            let mut bases = 0usize;
            while let Some(rec) = reader.next().expect("read") {
                records += 1;
                bases += rec.len();
                black_box(rec.seq());
            }
            black_box(records);
            black_box(bases);
        })
    });

    group.finish();
}

criterion_group!(benches, bench_mmap_vs_gzip);
criterion_main!(benches);
