# kira-fastq

FASTQ input and output for Rust. Plain, gzip and BGZF; mapped or streaming input; single- and
multi-line records; paired-end and interleaved reading; multi-threaded BGZF; optional async and
`noodles-bgzf` interop.

```toml
[dependencies]
kira-fastq = "0.4"
```

## Quick start

```rust,no_run
use kira_fastq::FastqReader;

fn main() -> Result<(), kira_fastq::FastqError> {
    let mut reader = FastqReader::from_path("reads.fastq.gz")?;
    let mut bases = 0u64;
    while let Some(rec) = reader.next()? {
        bases += rec.len() as u64;
    }
    println!("bases={bases}");
    Ok(())
}
```

`from_path` picks the backend from the file's **magic bytes**, not from its name: `bgzip` writes
BGZF into files called `.gz`, and pipelines routinely hand you gzip data called `.fastq`. Formats
this crate does not decode (bzip2, xz, and zstd without the `zstd` feature) are reported as
`FastqError::Unsupported` rather than as a parse error.

## What it accepts, and what it refuses

Reading follows what `kseq.h`, seqtk, seqkit and BioPython accept, because a reader that rejects
files the rest of the ecosystem reads cannot sit in a pipeline:

- the final record may be missing its trailing newline
- blank lines between records are skipped
- zero-length reads, as emitted by cutadapt, fastp and Trimmomatic, are records like any other
- `\r\n` is accepted wherever `\n` is
- the `+` line may repeat the read ID

Corruption is never passed off as clean data:

- gzip CRC32 and ISIZE are verified on every member, with no way to turn the check off
- a gzip stream that ends mid-member is an error, not a short read
- every BGZF block's CRC32 and ISIZE are verified
- a missing BGZF end-of-file marker is an error, which is how a truncated BGZF file presents
  itself; use `with_bgzf_eof_check(false)` when you knowingly read a file still being written
- every gzip member is decoded, so concatenated output from `bgzip`, `pigz` or `cat a.gz b.gz` is
  read in full

Errors carry both the byte offset in the decoded stream and the 1-based record index:

```text
sequence/quality length mismatch at offset 4211 (record 137): seq_len=151, qual_len=150
```

## Choosing an input path

| Constructor | Backend | Zero copy | `tell`/`seek` |
|---|---|---|---|
| `from_path` | mmap for plain, streaming inflate for gzip/BGZF | plain only | plain, BGZF |
| `from_path_buffered` | buffered `read`, no mapping | no | no |
| `from_bgzf_path_parallel(path, threads)` | BGZF inflate on a thread pool | no | no |
| `from_vec` / `from_bytes` | in-memory buffer | `from_vec` only | byte offsets |
| `from_reader` / `from_reader_auto` | any `BufRead`, optionally sniffed for gzip | no | no |
| `from_noodles_bgzf_path` | `noodles-bgzf` (feature) | no | yes |

`from_path` maps plain files and parses them in place, so records borrow straight out of the
mapping. `from_path_buffered` avoids mapping: use it for files another process is still writing
(a mapping is a snapshot and can fault if the file shrinks), for network file systems, and where
mapping measures slower, which it does on Windows.

## Reading from stdin, sockets, other decoders

```rust
use std::io::BufReader;
use kira_fastq::FastqReader;

fn from_stdin() -> Result<(), kira_fastq::FastqError> {
    // The reader owns its source, so pass the handle itself rather than a lock borrowed from it.
    // `from_reader_auto` decompresses the stream if it turns out to be gzip or BGZF.
    let mut reader = FastqReader::from_reader_auto(BufReader::new(std::io::stdin()))?;
    while let Some(rec) = reader.next()? {
        let _ = (rec.id(), rec.description(), rec.seq(), rec.qual());
    }
    Ok(())
}
```

## Multi-threaded BGZF

BGZF blocks are independent deflate streams, so inflate scales across cores the way `samtools -@`
does. Blocks are reassembled in file order, so the records are identical to the sequential path.

```rust
use kira_fastq::FastqReader;

# fn run() -> Result<(), kira_fastq::FastqError> {
// 0 sizes the pool from the machine's parallelism.
let mut reader = FastqReader::from_bgzf_path_parallel("reads.fastq.bgz", 0)?;
while let Some(rec) = reader.next()? {
    let _ = rec.seq();
}
# Ok(()) }
```

## Writing

```rust
use kira_fastq::{FastqRecord, FastqWriter};

fn write_some() -> Result<(), kira_fastq::FastqError> {
    let mut writer = FastqWriter::from_path("out.fastq.gz")?;
    writer.write_record(&FastqRecord::new(b"r0 lib=A", b"ACGT", b"!!!!"))?;
    // `finish` writes the format's trailer and reports errors from it. Dropping the writer also
    // finalises the file, but then a failure has nowhere to go.
    writer.finish()
}
```

Each record is assembled in a reusable buffer and emitted with one `write_all`; path-based
constructors wrap the file in a 1 MiB `BufWriter`. Output is always LF-terminated.

The extension picks the format: `.gz` gzip, `.bgz`/`.bgzf` BGZF, `.zst` zstd (feature `zstd`),
anything else plain. BGZF output is native, needs no optional feature, and can compress on a
thread pool with `to_bgzf_path_parallel(path, level, threads)` for byte-identical output.

Every write checks that sequence and quality have the same length and that the header holds no
line break, because a record that fails those cannot be read back by anything. Content checks are
opt-in:

```rust
use kira_fastq::{Alphabet, FastqWriter, WriteValidation};

# fn run() -> Result<(), kira_fastq::FastqError> {
let mut writer = FastqWriter::from_path("out.fastq")?
    .with_validation(WriteValidation::BasesAndQualities)
    .with_alphabet(Alphabet::Iupac);
# writer.finish() }
```

| `WriteValidation` | Adds |
|---|---|
| `None` | nothing beyond the always-on structural checks |
| `LineBreaks` | no `\n` or `\r` inside sequence or qualities |
| `Bases` | line breaks plus the base alphabet |
| `Qualities` | line breaks plus the quality range |
| `BasesAndQualities` | both |

## Paired-end and interleaved

```rust
use kira_fastq::{PairedFastqReader, ValidationMode};

# fn run() -> Result<(), kira_fastq::FastqError> {
let mut reader = PairedFastqReader::from_paths("r1.fastq.gz", "r2.fastq.gz")?
    .with_id_check(true)
    .with_validation(ValidationMode::BasesAndQualities);
while let Some((r1, r2)) = reader.next()? {
    let _ = (r1.id(), r2.id());
}
# Ok(()) }
```

`with_id_check(true)` compares canonicalised IDs: the Casava 1.8+ whitespace comment, a classic
Illumina `/1`,`/2` suffix, and an SRA `fastq-dump -I` style `.1`,`.2` suffix are all stripped
before comparison. A mismatch names both IDs in the error.

`InterleavedFastqReader` reads R1 and R2 alternating in one stream, the layout `samtools fastq`
writes by default.

## Multi-line records

```rust
use kira_fastq::{FastqFormat, FastqReader};
# fn run() -> Result<(), kira_fastq::FastqError> {
let mut reader = FastqReader::from_path("longreads.fastq")?
    .with_format(FastqFormat::MultiLine);
# Ok(()) }
```

Sequence lines run until a line starting with `+`, then quality lines run until they match the
sequence length, which is the only way to tell a wrapped quality line from the next record's `@`.

## Validation

| Knob | Values |
|---|---|
| `with_validation` | `None`, `Bases`, `Qualities`, `BasesAndQualities` |
| `with_alphabet` | `AcgtnStrict`, `AcgtnCase`, `Iupac` |
| `with_quality_encoding` | `QualityEncoding::PHRED33` (default), `PHRED64`, `custom(offset, min, max)` |

`guess_quality_encoding` reports the encoding of a sample of quality bytes when the two are
distinguishable, for tools that must handle Illumina 1.3 to 1.7 files.

## Checkpoint and resume

```rust
use kira_fastq::{FastqReader, VirtualOffset};

fn checkpoint() -> Result<(), kira_fastq::FastqError> {
    let mut reader = FastqReader::from_bgzf_path("reads.fastq.bgz")?;
    while let Some(_rec) = reader.next()? {
        break;
    }
    let voff = reader.tell();
    reader.seek(voff)?;
    Ok(())
}
```

`tell` returns an htslib-compatible virtual offset for BGZF (`coffset << 16 | uoffset`, with
`VirtualOffset::new`, `compressed()` and `uncompressed()` to build and inspect it) and a byte
offset for plain, in-memory and stream sources. Offsets round-trip even across a block holding
the full 64 KiB the spec allows: a drained block is reported as the next block at offset 0, since
an in-block offset of 65536 does not fit the field.

`seek` accepts an arbitrary offset too and resynchronises to the next record boundary. It checks
the shape of a whole record rather than just a leading `@`, because `@` is a legal quality byte,
and it understands multi-line records. `seek` works on plain, in-memory and BGZF sources; gzip
and arbitrary streams return `FastqError::Unsupported`.

## Async (feature `async`)

```rust,ignore
use kira_fastq::AnyAsyncReader;

let mut reader = AnyAsyncReader::from_path("reads.fastq.gz").await?;
while let Some(rec) = reader.next().await? {
    let _ = rec.seq();
}
```

No mapping on this path: async I/O goes through `tokio::io`, and records borrow the reader's
scratch buffer. Multi-member gzip, and therefore BGZF read as gzip, is decoded in full; BGZF
virtual offsets are not available here. `AsyncFastqReader::next` is **not cancel-safe**: a record
spans four reads, so a dropped future leaves the reader between lines. Drive it from one task and
use `records()` with a channel if other tasks need the records. Call `shutdown()` before dropping
a compressed writer.

## `noodles-bgzf` interop (feature `noodles-bgzf`)

`FastqReader::from_noodles_bgzf_path` reads through `noodles_bgzf`, supports `tell` and `seek`,
and its offsets convert both ways with `noodles_bgzf::VirtualPosition`. Use it when offsets travel
between this crate and `noodles-bam` or `noodles-vcf`. `FastqWriter::to_noodles_bgzf_path` is the
matching writer. The crate's own BGZF backend is mmap-backed and needs no optional dependency.

## Feature flags

| Feature | Default | What it adds |
|---|---|---|
| `async` | no | `AsyncFastqReader`, `AsyncFastqWriter`, `RecordStream` (tokio, async-compression) |
| `noodles-bgzf` | no | `noodles-bgzf` reader/writer adapters and `VirtualOffset` conversions |
| `zstd` | no | zstd input and `.zst` output |
| `libdeflate` | no | libdeflate instead of zlib-rs for BGZF block inflate |

Plain, gzip and BGZF input, BGZF output, paired-end, interleaved, multi-line, multi-threaded
BGZF and SIMD validation are all in the default build.

## Performance

Measured on one machine: Windows 11, AMD x86-64 with AVX2 and no AVX-512, 200 MiB corpus of
150 bp Illumina-shaped reads, warm page cache, best of three runs. Numbers are throughput against
the **uncompressed** size, so the compressed rows say how fast records come out. Rerun them with
`cargo bench` on your own hardware before trusting any of it: the compressed paths track the CPU
rather than the disk, and the validation kernels run fast enough that their numbers move with
clock speed from run to run.

| Path | 0.3.0 | 0.4.0 |
|---|---|---|
| plain, mmap | 2.1 GB/s | 2.1 GB/s |
| plain, buffered | 3.6 GB/s | 3.6 GB/s |
| gzip | 0.5 GB/s | 0.9 GB/s |
| BGZF, 1 thread | 0.5 GB/s | 0.9 GB/s |
| BGZF, 8 threads | not available | 4.4 GB/s |
| base validation, IUPAC | 3.3 GB/s | 10-19 GB/s |
| quality validation | 13 GB/s | 13 GB/s |

Where the gains come from: inflate moved from `miniz_oxide` to zlib-rs (and the gzip and BGZF
figures now include CRC checking, which the 0.3 numbers did not do at all); base validation moved
from a scalar table lookup wearing a `#[target_feature]` hat to a real nibble-table shuffle
kernel; BGZF gained a thread pool.

Two claims from the 0.3 README that measurement did not support have been dropped rather than
defended. Mapping is not uniformly faster than buffered reading: on Windows it is slower, which is
why `from_path_buffered` exists. The hand-written AVX2/AVX-512 line scanner was no faster than
`memchr`, which dispatches its own kernels, so it is gone along with its `unsafe`.

## Public API

### `FastqReader`

`from_path`, `from_path_buffered`, `from_bgzf_path`, `from_bgzf_path_parallel`, `from_vec`,
`from_bytes`, `from_reader`, `from_reader_auto`, `from_unbuffered`, `from_noodles_bgzf_path`
(feature); builders `with_validation`, `with_alphabet`, `with_quality_encoding`, `with_format`,
`with_bgzf_eof_check`; `next`, `try_for_each`, `records_read`, `tell`, `seek`.

### `FastqWriter`

`from_path`, `to_plain_path`, `to_gz_path`, `to_bgzf_path`, `to_bgzf_path_parallel`,
`to_zstd_path` (feature), `to_noodles_bgzf_path` (feature), `from_writer`; builders
`with_validation`, `with_alphabet`, `with_quality_encoding`; `write_record`,
`write_record_owned`, `write_parts`, `flush`, `finish`, `into_inner`.

### Records

`FastqRecord` (borrowed) and `FastqRecordOwned` (owned) with `header`, `id`, `description`,
`seq`, `qual`, `len`; `to_owned` and `copy_from` cross between them, the latter reusing the
existing allocations.

### Errors

`FastqError` with `Io`, `InvalidFormat`, `UnexpectedEof`, `LengthMismatch`, `InvalidBase`,
`InvalidQuality`, `PairedCountMismatch`, `PairedIdMismatch`, `Unsupported`; `record()` and
`offset()` accessors. All enums are `#[non_exhaustive]`.

## Changes in 0.4

Bug fixes that change behaviour:

- the async reader no longer stops at the first gzip member, which silently dropped most of a
  BGZF or `pigz` file
- BGZF `tell` is correct at a full 64 KiB block boundary, so resuming from a checkpoint no longer
  replays a block
- truncated gzip and BGZF files are errors instead of clean end of input
- the writer refuses records whose sequence and quality lengths disagree, and headers containing
  line breaks
- files that every mainstream tool reads are accepted: missing final newline, zero-length reads,
  blank lines between records
- `from_path` selects the backend by content, not by file extension
- `seek` resynchronises correctly for multi-line records and on BGZF, and `tell` afterwards
  reports the resynchronised position

API changes: `FastqError` variants carry a `record` field, `PairedLengthMismatch` is now
`PairedCountMismatch`, `WriteValidation` gained `LineBreaks`, the `gzip-validate` feature is gone
because the checks it enabled are now always on, `BoxedWriter` is an alias for `PathWriter`, and
the `backend` and `parser` modules are no longer public.

## MSRV

Rust 1.89. Edition 2024 sets the floor at 1.85; the SIMD kernels call `core::arch`
intrinsics that only became callable without an `unsafe` block in 1.87, and the optional
`noodles-bgzf` dependency requires 1.89. MSRV changes only in a minor or major release.

## License

MIT.
