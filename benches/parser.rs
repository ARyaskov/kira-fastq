use std::path::Path;

use criterion::{Criterion, Throughput, black_box, criterion_group, criterion_main};
use kira_fastq::FastqReader;

fn bench_parser(c: &mut Criterion) {
    let path = Path::new("benches/data/sample.fastq");
    let size = std::fs::metadata(path).expect("stat").len();
    let mut group = c.benchmark_group("parser_plain");
    group.throughput(Throughput::Bytes(size));
    group.bench_function("parse_plain", |b| {
        b.iter(|| {
            let mut reader = FastqReader::from_path(path).expect("open");
            let mut records = 0usize;
            let mut bases = 0usize;
            while let Some(rec) = reader.next().expect("read") {
                records += 1;
                bases += rec.len();
                black_box(rec.seq());
                black_box(rec.qual());
            }
            black_box(records);
            black_box(bases);
        })
    });
    group.finish();
}

criterion_group!(benches, bench_parser);
criterion_main!(benches);
