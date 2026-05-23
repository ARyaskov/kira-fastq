# kira-fastq

High-performance FASTQ reader with mmap-first design. Supports plain, gzip, and
BGZF inputs; optional multi-line parsing; and explicit paired-end reading.

## Usage

### Plain FASTQ: quick QC summary (reads, bases, GC%)

```rust
use kira_fastq::FastqReader;

fn main() -> Result<(), kira_fastq::FastqError> {
    let mut reader = FastqReader::from_path("reads.fastq")?;
    let mut reads: u64 = 0;
    let mut bases: u64 = 0;
    let mut gc: u64 = 0;

    while let Some(rec) = reader.next()? {
        reads += 1;
        bases += rec.len() as u64;
        gc += rec
            .seq()
            .iter()
            .filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
            .count() as u64;
    }

    let gc_pct = if bases > 0 {
        (gc as f64) * 100.0 / (bases as f64)
    } else {
        0.0
    };
    println!("reads={reads} bases={bases} gc_pct={gc_pct:.2}");
    Ok(())
}
```

### Gzipped FASTQ (`.fastq.gz`): filter reads for downstream analysis

```rust
use kira_fastq::FastqReader;

fn mean_q(qual: &[u8]) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    let sum: u64 = qual.iter().map(|&q| (q.saturating_sub(33)) as u64).sum();
    sum as f64 / qual.len() as f64
}

fn main() -> Result<(), kira_fastq::FastqError> {
    let mut reader = FastqReader::from_path_auto("reads.fastq.gz")?;

    let mut kept = 0u64;
    let mut dropped = 0u64;
    let min_len = 75usize;
    let min_mean_q = 25.0f64;

    while let Some(rec) = reader.next()? {
        let pass = rec.len() >= min_len && mean_q(rec.qual()) >= min_mean_q;
        if pass {
            kept += 1;
            // Example: push rec.seq() into k-mer counter / mapper / assembler queue.
        } else {
            dropped += 1;
        }
    }
    println!("kept={kept} dropped={dropped}");
    Ok(())
}
```

### BGZF FASTQ (`.fastq.bgz` / `.fastq.bgzf`): checkpoint and resume

```rust
use kira_fastq::{FastqReader, VirtualOffset};

fn main() -> Result<(), kira_fastq::FastqError> {
    let mut reader = FastqReader::from_bgzf_path("reads.fastq.bgz")?;

    // Save checkpoint every 100k reads for fault-tolerant pipelines.
    let mut checkpoint = VirtualOffset(0);
    for i in 0..100_000u64 {
        if reader.next()?.is_none() {
            break;
        }
        if i % 10_000 == 0 {
            checkpoint = reader.tell();
        }
    }

    // Later: resume from checkpoint and continue processing.
    reader.seek(checkpoint)?;
    let mut resumed = 0u64;
    while let Some(rec) = reader.next()? {
        resumed += 1;
        let _ = rec.seq(); // downstream processing
    }
    println!("resumed_reads={resumed}");
    Ok(())
}
```

### Paired FASTQ (R1/R2): synchronized processing with ID check

```rust
use kira_fastq::{PairedFastqReader, ValidationMode};

fn main() -> Result<(), kira_fastq::FastqError> {
    let mut reader = PairedFastqReader::from_paths("r1.fastq.gz", "r2.fastq.gz")?
        .with_id_check(true)
        .with_validation(ValidationMode::BasesAndQualities);

    let mut pairs = 0u64;
    let mut sum_bases_r1 = 0u64;
    let mut sum_bases_r2 = 0u64;

    while let Some((r1, r2)) = reader.next()? {
        pairs += 1;
        sum_bases_r1 += r1.len() as u64;
        sum_bases_r2 += r2.len() as u64;
        // Example: joint downstream step (paired mapper / UMI processing / consensus).
    }

    println!("pairs={pairs} bases_r1={sum_bases_r1} bases_r2={sum_bases_r2}");
    Ok(())
}
```

## Public API (basic reference)

### `FastqReader`

- `FastqReader::from_path(path)`  
  Open by extension (`.gz`, `.bgz`, `.bgzf`, plain).
- `FastqReader::from_path_auto(path)`  
  Detect compression by file signature (recommended for unknown input).
- `FastqReader::from_bgzf_path(path)`  
  Explicit BGZF reader with virtual offsets.
- `with_validation(mode)`  
  Enable validation: `None | Bases | Qualities | BasesAndQualities`.
- `with_format(format)`  
  Parsing mode: `FastqFormat::SingleLine` or `FastqFormat::MultiLine`.
- `next()`  
  Read next record: `Result<Option<FastqRecord>, FastqError>`.
- `records()`  
  Iterator API over records.
- `tell()`  
  Get current `VirtualOffset`.
- `seek(voff)`  
  Seek for plain and BGZF; returns `Unsupported(Seek)` for gzip.

### `PairedFastqReader`

- `PairedFastqReader::from_paths(r1, r2)`  
  Open synchronized paired-end readers.
- `with_validation(mode)`  
  Apply same validation mode to both ends.
- `with_id_check(true/false)`  
  Require matching read IDs (`header` prefix before whitespace).
- `with_format(format)`  
  Single-line or multi-line for both ends.
- `next()`  
  Read next pair: `Result<Option<(FastqRecord, FastqRecord)>, FastqError>`.

### `FastqRecord`

- `header()`  
  Read ID/header bytes (without `@`).
- `seq()`  
  Sequence bytes.
- `qual()`  
  Quality bytes (Phred+33).
- `len()`  
  Sequence length (same as quality length for valid records).

### Supporting enums/types

- `FastqFormat`  
  `SingleLine`, `MultiLine`.
- `ValidationMode`  
  `None`, `Bases`, `Qualities`, `BasesAndQualities`.
- `VirtualOffset(pub u64)`  
  BGZF virtual offset (`tell/seek`).
- `FastqError`  
  Main error enum (`InvalidFormat`, `UnexpectedEof`, `LengthMismatch`, etc.).
- `InvalidKind`  
  Detailed format failure kind (`HeaderMissingAt`, `PlusMissing`, `GzipData`, `BgzfBlock`, ...).

## MSRV

Rust **1.95+** (2024 edition). MSRV bumps only in minor or major releases.
