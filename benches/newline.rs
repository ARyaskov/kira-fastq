use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};

use kira_fastq::simd::newline::find_lf;

fn make_buf(size: usize, step: usize) -> Vec<u8> {
    let mut buf = vec![b'A'; size];
    if step > 0 {
        let mut i = step;
        while i <= size {
            buf[i - 1] = b'\n';
            i += step;
        }
    }
    buf
}

fn scan_all(buf: &[u8]) -> usize {
    let mut count = 0usize;
    let mut start = 0usize;
    while let Some(pos) = find_lf(buf, start) {
        count += 1;
        start = pos + 1;
    }
    count
}

fn bench_newline(c: &mut Criterion) {
    let sizes = [1usize << 20, 4usize << 20, 16usize << 20];
    let densities = [
        ("sparse_4k", 4096usize),
        ("moderate_256", 256usize),
        ("dense_32", 32usize),
    ];

    let mut group = c.benchmark_group("newline");
    for &size in &sizes {
        for (label, step) in densities {
            let buf = make_buf(size, step);
            group.throughput(Throughput::Bytes(size as u64));
            group.bench_with_input(BenchmarkId::new(label, size), &buf, |b, data| {
                b.iter(|| {
                    let count = scan_all(black_box(data));
                    black_box(count);
                })
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_newline);
criterion_main!(benches);
