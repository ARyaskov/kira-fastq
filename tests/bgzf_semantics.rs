//! BGZF behaviour that has to match htslib: virtual offsets, integrity checks, resync.

mod common;

use kira_fastq::{FastqError, FastqFormat, FastqReader, InvalidKind, VirtualOffset};

/// 64-byte records, so a whole number of them fills a block exactly.
fn records(n: usize) -> Vec<u8> {
    let mut out = Vec::with_capacity(n * 64);
    for _ in 0..n {
        out.extend_from_slice(b"@rr\n");
        out.extend(std::iter::repeat_n(b'A', 28));
        out.extend_from_slice(b"\n+\n");
        out.extend(std::iter::repeat_n(b'!', 28));
        out.push(b'\n');
    }
    out
}

/// A block holding the full 64 KiB the spec allows leaves the in-block offset at 65536, which
/// does not fit the 16-bit field of a virtual offset. `tell` has to report the next block
/// instead, or resuming from the checkpoint replays the whole block.
#[test]
fn tell_and_seek_round_trip_across_a_full_64k_block() {
    let path = common::unique_path("full_block.bgz");
    let data = records(1024 + 64);
    common::write_bgzf(&path, &data, 65536);

    let mut reader = FastqReader::from_bgzf_path(&path).expect("open");
    // 1024 records * 64 bytes = exactly one full block.
    for _ in 0..1024 {
        reader.next().expect("read").expect("record");
    }
    let checkpoint = reader.tell();
    assert_eq!(
        checkpoint.uncompressed(),
        0,
        "a drained block must be reported as the next block at offset 0"
    );
    let after: Vec<Vec<u8>> = collect_seqs(&mut reader);

    let mut resumed = FastqReader::from_bgzf_path(&path).expect("reopen");
    resumed.seek(checkpoint).expect("seek");
    let replayed = collect_seqs(&mut resumed);
    assert_eq!(after.len(), 64);
    assert_eq!(
        after, replayed,
        "resuming must not duplicate or drop records"
    );
}

fn collect_seqs(reader: &mut FastqReader) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    while let Some(rec) = reader.next().expect("read") {
        out.push(rec.seq().to_vec());
    }
    out
}

#[test]
fn tell_after_every_record_round_trips() {
    let path = common::unique_path("checkpoints.bgz");
    let data = records(200);
    common::write_bgzf(&path, &data, 1000);

    let mut reader = FastqReader::from_bgzf_path(&path).expect("open");
    let mut checkpoints = vec![reader.tell()];
    while reader.next().expect("read").is_some() {
        checkpoints.push(reader.tell());
    }
    assert_eq!(checkpoints.len(), 201);

    for (i, voff) in checkpoints.iter().enumerate() {
        let mut r = FastqReader::from_bgzf_path(&path).expect("reopen");
        r.seek(*voff).expect("seek");
        let remaining = collect_seqs(&mut r).len();
        assert_eq!(remaining, 200 - i, "checkpoint {i}");
    }
}

/// `@` is a legal quality byte, so seeking into a quality line must not be mistaken for a
/// record start.
#[test]
fn seek_resyncs_past_an_at_sign_inside_a_quality_line() {
    let path = common::unique_path("at_in_qual.bgz");
    common::write_bgzf(&path, b"@r1\nACGT\n+\n!@!!\n@r2\nTTTT\n+\nBBBB\n", 1000);

    let mut reader = FastqReader::from_bgzf_path(&path).expect("open");
    // Offset 12 is the '@' inside the first record's quality line.
    reader.seek(VirtualOffset::new(0, 12)).expect("seek");
    let rec = reader.next().expect("read").expect("record");
    assert_eq!(rec.header(), b"r2");
}

#[test]
fn seek_resyncs_in_multi_line_mode() {
    let path = common::unique_path("multiline_resync.bgz");
    common::write_bgzf(&path, b"@r1\nACGT\nAC\n+\n!!!!\n!!\n@r2\nTT\n+\n##\n", 1000);

    let mut reader = FastqReader::from_bgzf_path(&path)
        .expect("open")
        .with_format(FastqFormat::MultiLine);
    reader.seek(VirtualOffset::new(0, 0)).expect("seek");
    let rec = reader.next().expect("read").expect("record");
    assert_eq!(
        rec.header(),
        b"r1",
        "seek(0) must not skip the first record"
    );
    assert_eq!(rec.seq(), b"ACGTAC");
}

/// The 28-byte marker is how a truncated BGZF file announces itself.
#[test]
fn missing_eof_marker_is_reported() {
    let path = common::unique_path("no_eof.bgz");
    common::write_bgzf_without_eof(&path, b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n", 1000);

    let mut reader = FastqReader::from_bgzf_path(&path).expect("open");
    let mut seen = 0;
    let err = loop {
        match reader.next() {
            Ok(Some(_)) => seen += 1,
            Ok(None) => panic!("clean EOF hides the truncation"),
            Err(e) => break e,
        }
    };
    assert_eq!(seen, 2, "records before the truncation still come through");
    match err {
        FastqError::InvalidFormat { kind, .. } => {
            assert_eq!(kind, InvalidKind::BgzfMissingEofMarker);
        }
        other => panic!("unexpected error: {other}"),
    }
}

#[test]
fn missing_eof_marker_check_can_be_turned_off() {
    let path = common::unique_path("no_eof_opt_out.bgz");
    common::write_bgzf_without_eof(&path, b"@r1\nACGT\n+\n!!!!\n", 1000);

    let mut reader = FastqReader::from_bgzf_path(&path)
        .expect("open")
        .with_bgzf_eof_check(false);
    assert_eq!(collect_seqs(&mut reader).len(), 1);
}

#[test]
fn corrupt_block_payload_is_reported() {
    let path = common::unique_path("bad_block.bgz");
    common::write_bgzf(&path, &records(40), 1000);
    let mut bytes = std::fs::read(&path).expect("read");
    // Flip a byte inside the first block's deflate payload.
    bytes[30] ^= 0xFF;
    std::fs::write(&path, bytes).expect("write");

    let mut reader = FastqReader::from_bgzf_path(&path).expect("open");
    let err = reader.next().expect_err("corrupt block must fail");
    match err {
        FastqError::InvalidFormat { kind, .. } => assert!(
            matches!(
                kind,
                InvalidKind::BgzfBlock | InvalidKind::BgzfBlockCrc | InvalidKind::BgzfBlockIsize
            ),
            "unexpected kind: {kind:?}"
        ),
        other => panic!("unexpected error: {other}"),
    }
}

/// Errors carry the record they refer to, which is what makes them actionable on a file with
/// millions of reads.
#[test]
fn errors_report_the_record_index() {
    let path = common::unique_path("record_index.fastq");
    common::write_plain(&path, b"@a\nAC\n+\n!!\n@b\nGT\n+\n##\n@c\nACGT\n+\n!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    reader.next().expect("read").expect("first");
    reader.next().expect("read").expect("second");
    let err = reader.next().expect_err("third record is malformed");
    assert_eq!(err.record(), Some(3));
    assert!(format!("{err}").contains("record 3"));
}

/// A long read wrapped over hundreds of lines still has to be recognised as a record start, so
/// resync cannot cap its lookahead at a handful of lines.
#[test]
fn seek_resyncs_across_a_long_wrapped_record() {
    let mut data = Vec::new();
    data.extend_from_slice(b"@long\n");
    for _ in 0..400 {
        data.extend_from_slice(b"ACGTACGTACGTACGTACGT\n");
    }
    data.extend_from_slice(b"+\n");
    for _ in 0..400 {
        data.extend_from_slice(b"!!!!!!!!!!!!!!!!!!!!\n");
    }
    data.extend_from_slice(b"@second\nACGT\n+\n!!!!\n");

    let path = common::unique_path("long_wrapped.bgz");
    common::write_bgzf(&path, &data, 4096);
    let mut reader = FastqReader::from_bgzf_path(&path)
        .expect("open")
        .with_format(FastqFormat::MultiLine);
    reader.seek(VirtualOffset::new(0, 0)).expect("seek");
    let rec = reader.next().expect("read").expect("record");
    assert_eq!(rec.header(), b"long");
    assert_eq!(rec.seq().len(), 400 * 20);
    assert_eq!(
        reader.next().expect("read").expect("second").header(),
        b"second"
    );
}
