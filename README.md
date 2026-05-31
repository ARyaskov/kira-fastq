# kira-fastq

High-performance FASTQ I/O for Rust. mmap-first reading, single-syscall-per-record writing,
SIMD validation on every hot path, paired-end and multi-line parsing, optional async,
optional `noodles-bgzf` interop.

```toml
[dependencies]
kira-fastq = "0.3"
```

## Quick start

### Reading: plain, gzip, BGZF — one API, auto-detected by magic bytes

```rust
use kira_fastq::FastqReader;

fn main() -> Result<(), kira_fastq::FastqError> {
    let mut reader = FastqReader::from_path_auto("reads.fastq.gz")?;
    let mut bases = 0u64;
    while let Some(rec) = reader.next()? {
        bases += rec.len() as u64;
    }
    println!("bases={bases}");
    Ok(())
}
```

`from_path_auto` sniffs the file's magic bytes (gzip `1f 8b`, BGZF `BC` extra subfield) and
selects the right backend. `from_path` chooses by extension; `from_bgzf_path` forces BGZF.
**Plain files always go through mmap; gzip and BGZF use mmap for the compressed input plus
a streaming inflate buffer.**

### Reading from arbitrary streams (stdin, sockets, custom decoders)

```rust
use std::io::BufReader;
use kira_fastq::FastqReader;

fn from_stdin() -> Result<(), kira_fastq::FastqError> {
    let stdin = std::io::stdin().lock();
    let mut reader = FastqReader::from_reader(BufReader::new(stdin));
    while let Some(rec) = reader.next()? {
        // process rec.header() / rec.seq() / rec.qual()
        let _ = rec;
    }
    Ok(())
}
```

`from_reader` accepts any `BufRead + Send + 'static`. Use it for pipes, TCP, in-memory
buffers, or anything wrapping a non-file source. **No mmap on this path** — lines are
copied into the reader's scratch buffer (≈ 30 % slower than mmap-on-file in the common
case; the difference shrinks for short records and disappears for compressed input).

### Writing: plain, gzip, BGZF

```rust
use kira_fastq::{FastqWriter, FastqRecord};

fn write_some() -> Result<(), kira_fastq::FastqError> {
    let mut writer = FastqWriter::from_path("out.fastq.gz")?;
    let rec = FastqRecord::new(b"r0 lib=A", b"ACGT", b"!!!!");
    writer.write_record(&rec)?;
    writer.flush()?;
    Ok(())
}
```

Each record goes through a reusable scratch buffer and is emitted with **one
`write_all`** — four lines, one syscall (modulo the 1 MiB `BufWriter` that all path-based
constructors wrap around the file). Output is always LF-terminated.

Extension-driven format selection: `.gz` → gzip (zlib-rs), `.bgz` / `.bgzf` → BGZF (requires
the `noodles-bgzf` feature), anything else → plain.

For a custom sink (in-memory `Vec<u8>`, a `BufWriter` over a `TcpStream`, a pre-allocated
`mmap`), use `FastqWriter::from_writer(w)` with any `W: Write`. Buffering is your
responsibility on that path.

### Opt-in SIMD validation on the writer

```rust
use kira_fastq::{FastqWriter, WriteValidation, Alphabet};

fn validated() -> Result<(), kira_fastq::FastqError> {
    let mut writer = FastqWriter::from_path("out.fastq")?
        .with_validation(WriteValidation::BasesAndQualities)
        .with_alphabet(Alphabet::Iupac);
    let rec = kira_fastq::FastqRecord::new(b"r0", b"ACGT", b"!!!!");
    writer.write_record(&rec)?;
    Ok(())
}
```

Default is `WriteValidation::None` — the write path stays branch-light. Turn validation on
for forensic exports, public-data publishing, or in tests.

### Paired-end (R1/R2) with ID synchronisation

```rust
use kira_fastq::{PairedFastqReader, ValidationMode};

fn paired() -> Result<(), kira_fastq::FastqError> {
    let mut reader = PairedFastqReader::from_paths("r1.fastq.gz", "r2.fastq.gz")?
        .with_id_check(true)
        .with_validation(ValidationMode::BasesAndQualities);
    while let Some((r1, r2)) = reader.next()? {
        let _ = (r1, r2);
    }
    Ok(())
}
```

`with_id_check(true)` canonicalises read IDs (Casava 1.8+ whitespace suffix and classic
Illumina `/1`,`/2` are stripped before comparison).

### Multi-line FASTQ

```rust
use kira_fastq::{FastqReader, FastqFormat};
let mut reader = FastqReader::from_path("longreads.fastq")?
    .with_format(FastqFormat::MultiLine);
# Ok::<_, kira_fastq::FastqError>(())
```

### Checkpoint / resume on BGZF

```rust
use kira_fastq::{FastqReader, VirtualOffset};

fn checkpoint() -> Result<(), kira_fastq::FastqError> {
    let mut reader = FastqReader::from_bgzf_path("reads.fastq.bgz")?;
    // ...read some records...
    let voff = reader.tell();
    // later:
    reader.seek(voff)?;
    Ok(())
}
```

`tell()` is also defined for plain (returns byte offset) and stream backends; `seek()` is
supported for plain and BGZF, returns `Unsupported(Seek)` for gzip and stream sources.

## Async (feature = `async`)

```toml
[dependencies]
kira-fastq = { version = "0.3", features = ["async"] }
```

> ⚠️ **No mmap on the async path.** mmap is inherently synchronous; async I/O goes through
> `tokio::io`. The zero-copy borrowed-record API survives (records borrow into the reader's
> scratch buffer), but the underlying source/sink goes through tokio's syscall-backed I/O.
> Expect **~10–30 % lower throughput** than the sync mmap path for plain files; gzip is
> typically I/O-bound and the gap shrinks to near zero.
>
> **Reach for sync first.** Use async when you genuinely need it — long-lived servers,
> backpressure-aware pipelines, integration with axum / tonic / sqlx that already share a
> tokio runtime. For a CLI tool that reads a file and exits, sync is faster and simpler.

### Async reader — owned-record stream

```rust
use futures_util::StreamExt;
use kira_fastq::{AsyncFastqReader, ValidationMode};
use tokio::io::BufReader;
use tokio::fs::File;

#[tokio::main]
async fn main() -> Result<(), kira_fastq::FastqError> {
    let file = File::open("reads.fastq").await.map_err(kira_fastq::FastqError::Io)?;
    let reader = AsyncFastqReader::from_reader(BufReader::new(file))
        .with_validation(ValidationMode::BasesAndQualities);

    let mut stream = reader.records();
    let mut n = 0u64;
    while let Some(rec) = stream.next().await {
        let _rec = rec?;
        n += 1;
    }
    println!("n={n}");
    Ok(())
}
```

Or with borrowed records and explicit `next().await`:

```rust
use kira_fastq::AsyncFastqReader;
use tokio::fs::File;
use tokio::io::BufReader;

# async fn run() -> Result<(), kira_fastq::FastqError> {
let file = File::open("reads.fastq").await.map_err(kira_fastq::FastqError::Io)?;
let mut reader = AsyncFastqReader::from_reader(BufReader::new(file));
while let Some(rec) = reader.next().await? {
    let _ = rec.seq();
}
# Ok(()) }
```

For path-based input with gzip magic-byte sniffing, use [`AnyAsyncReader::from_path`]:

```rust
use kira_fastq::AnyAsyncReader;

# async fn run() -> Result<(), kira_fastq::FastqError> {
let mut reader = AnyAsyncReader::from_path("reads.fastq.gz").await?;
while let Some(rec) = reader.next().await? {
    let _ = rec.seq();
}
# Ok(()) }
```

**BGZF is decoded as gzip on the async path** — virtual-offset semantics are not preserved.
If you need true BGZF async, wrap a `noodles_bgzf::AsyncReader` via
`AsyncFastqReader::from_reader(...)`.

### Async writer

```rust
use kira_fastq::{AsyncFastqWriter, FastqRecord};

# async fn run() -> Result<(), kira_fastq::FastqError> {
let mut writer = AsyncFastqWriter::from_path("out.fastq.gz").await?;
let rec = FastqRecord::new(b"r0", b"ACGT", b"!!!!");
writer.write_record(&rec).await?;
writer.shutdown().await?;  // important for gzip — flushes deflate trailer
# Ok(()) }
```

Same single-`write_all`-per-record discipline as the sync writer.

## `noodles-bgzf` interop (feature = `noodles-bgzf`)

```toml
[dependencies]
kira-fastq = { version = "0.3", features = ["noodles-bgzf"] }
```

Adds a thin adapter that wraps `noodles_bgzf::Reader` / `noodles_bgzf::Writer` behind
kira's reader/writer API. Use this when downstream code expects BGZF virtual offsets
compatible with the rest of the noodles ecosystem (`noodles-bam`, `noodles-vcf`, etc.).

```rust
# #[cfg(feature = "noodles-bgzf")]
# fn run() -> Result<(), kira_fastq::FastqError> {
let mut reader = kira_fastq::FastqReader::from_noodles_bgzf_path("reads.fastq.bgz")?;
let mut writer = kira_fastq::FastqWriter::to_noodles_bgzf_path("out.fastq.bgz")?;
while let Some(rec) = reader.next()? {
    writer.write_record(&rec)?;
}
writer.flush()?;
# Ok(()) }
```

The default `from_bgzf_path` / `to_*` writers use kira's own BGZF implementation, which is
mmap-backed and a touch faster on single-threaded inflate.

## Feature flags

| Feature           | Default | What it enables |
|-------------------|---------|-----------------|
| `default`         | ✅      | sync reader + writer (plain, gzip via miniz_oxide, BGZF via kira's own decoder), SIMD validation, paired-end |
| `async`           | ❌      | `AsyncFastqReader`, `AsyncFastqWriter`, `RecordStream` (tokio + futures-util + async-compression for gzip) |
| `noodles-bgzf`    | ❌      | `FastqReader::from_noodles_bgzf_path`, `FastqWriter::to_noodles_bgzf_path` |
| `libdeflate`      | ❌      | use `libdeflater` instead of `miniz_oxide` for BGZF inflate |
| `gzip-validate`   | ❌      | CRC32 + ISIZE check on every gzip member (correctness over speed) |

## Performance contract

| Path                                | mmap | SIMD | Zero-copy |
|-------------------------------------|------|------|-----------|
| `FastqReader::from_path` (plain)    | ✅   | ✅   | ✅ |
| `FastqReader::from_path` (gzip)     | ✅\* | ✅   | ✅ via scratch |
| `FastqReader::from_path` (BGZF)     | ✅\* | ✅   | ✅ via scratch |
| `FastqReader::from_reader`          | ❌   | ✅   | ✅ via scratch |
| `FastqReader::from_noodles_bgzf_path` | ❌ | ✅   | ✅ via scratch |
| `AsyncFastqReader::*`               | ❌   | ✅   | ✅ via scratch (sync `next`); owned in `records()` |
| `FastqWriter::*`                    | n/a  | ✅\*\* | n/a |
| `AsyncFastqWriter::*`               | n/a  | ✅\*\* | n/a |

\* mmap covers the *compressed* input; decompressed bytes flow through an inflate buffer.
\*\* SIMD on writer side is only used when `WriteValidation` is non-`None`.

## Public API (0.3.0)

### `FastqReader`

| Method | Purpose |
|---|---|
| `from_path(path)` | Backend chosen by extension (`.gz`, `.bgz`, `.bgzf`, plain). |
| `from_path_auto(path)` | Backend chosen by sniffing the first 128 bytes of the file. |
| `from_bgzf_path(path)` | Force BGZF with virtual-offset semantics. |
| `from_reader<R: BufRead + Send + 'static>(r)` | Read from any `BufRead` source. No mmap. |
| `from_unbuffered<R: Read + Send + 'static>(r)` | Wraps `r` in a 256 KiB `BufReader`. |
| `from_noodles_bgzf_path(path)` | (feature `noodles-bgzf`) Use noodles BGZF semantics. |
| `with_validation(mode)` | `None` \| `Bases` \| `Qualities` \| `BasesAndQualities`. |
| `with_alphabet(alphabet)` | `AcgtnStrict` \| `AcgtnCase` \| `Iupac`. |
| `with_format(format)` | `SingleLine` (default) \| `MultiLine`. |
| `next()` | Borrowed `FastqRecord<'_>`. |
| `try_for_each(f)` | Callback iteration with user-error propagation. |
| `tell()` / `seek(voff)` | Checkpoint/resume; `seek` errors on gzip and stream. |

### `FastqWriter`

| Method | Purpose |
|---|---|
| `from_path(path)` | Auto-format from extension; wraps a 1 MiB `BufWriter`. |
| `to_plain_path(path, buf)` | Force plain output with explicit buffer size. |
| `to_gz_path(path, level)` | Force gzip; `level` is 0–9 (zlib-rs). |
| `to_noodles_bgzf_path(path)` | (feature `noodles-bgzf`) BGZF via noodles writer. |
| `from_writer<W: Write>(w)` | Wrap any `Write`. No automatic buffering. |
| `with_validation(mode)` | Opt-in pre-write SIMD validation. |
| `with_alphabet(alphabet)` | Alphabet for base validation. |
| `write_record(&rec)` | Single `write_all` per record. |
| `write_record_owned(&rec_owned)` | Same path, owned input. |
| `write_parts(header, seq, qual)` | Skip building a record type. |
| `flush()` / `into_inner()` | Standard. |

### `AsyncFastqReader<R: AsyncBufRead + Unpin + Send>` (feature `async`)

| Method | Purpose |
|---|---|
| `from_reader(r)` | Construct from any `AsyncBufRead`. |
| `with_validation` / `with_alphabet` / `with_format` | Same as sync. |
| `next()` | Borrowed record (`async`). |
| `records()` | Convert into a `RecordStream` of owned records (`futures::Stream`). |

### `AsyncFastqWriter<W: AsyncWrite + Unpin + Send>` (feature `async`)

| Method | Purpose |
|---|---|
| `from_writer(w)` | Construct from any `AsyncWrite`. |
| `from_path(path)` | `tokio::fs::File` with `.gz` → streaming gzip via `async-compression`. |
| `with_validation` / `with_alphabet` | Same as sync. |
| `write_record(&rec)` / `write_record_owned` / `write_parts` | All `async`. |
| `flush()` / `shutdown()` | `shutdown` flushes gzip trailer — call it before drop. |

### `FastqRecord<'a>` / `FastqRecordOwned`

Borrowed and owned variants. `FastqRecord::to_owned()` and
`FastqRecordOwned::as_borrowed()` cross the boundary. Use borrowed in tight read loops;
owned when the record must outlive `next()`-scope (channels, async streams, threads).

### `PairedFastqReader`

`from_paths(r1, r2)`, `with_id_check`, `with_validation`, `with_alphabet`, `with_format`,
`next()`.

### Errors

`FastqError` with structured variants:
`Io`, `InvalidFormat{offset, kind}`, `UnexpectedEof{offset}`, `LengthMismatch{...}`,
`InvalidBase{offset, byte}`, `InvalidQuality{offset, byte}`,
`PairedLengthMismatch{which}`, `PairedIdMismatch{...}`, `Unsupported(op)`. All
`#[non_exhaustive]`.

## MSRV

Rust **1.95+** (edition 2024). MSRV bumps only in minor or major releases.

## License

MIT.
