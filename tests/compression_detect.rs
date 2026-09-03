//! Which backend a file gets, and what happens for formats this crate does not decode.

mod common;

use kira_fastq::{FastqError, FastqReader, UnsupportedOperation, VirtualOffset};

/// `bgzip` writes BGZF under a plain `.gz` name, and that file supports virtual offsets. Picking
/// the backend by extension would send it to the gzip path and lose `seek`.
#[test]
fn bgzf_named_gz_still_gets_the_bgzf_backend() {
    let path = common::unique_path("actually_bgzf.fastq.gz");
    common::write_bgzf(&path, b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n", 1000);

    let mut reader = FastqReader::from_path(&path).expect("open");
    let first = reader.next().expect("read").expect("record");
    assert_eq!(first.header(), b"r1");
    let checkpoint = reader.tell();
    reader.seek(checkpoint).expect("BGZF supports seek");
    assert_eq!(
        reader.next().expect("read").expect("record").header(),
        b"r2"
    );
}

/// The reverse case: plain text under a `.gz` name must not be fed to the inflater.
#[test]
fn plain_text_named_gz_is_read_as_text() {
    let path = common::unique_path("not_really.fastq.gz");
    common::write_plain(&path, b"@r1\nACGT\n+\n!!!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    assert_eq!(reader.next().expect("read").expect("record").seq(), b"ACGT");
}

#[test]
fn gzip_named_fastq_is_decompressed() {
    let path = common::unique_path("compressed.fastq");
    common::write_gzip(&path, b"@r1\nACGT\n+\n!!!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    assert_eq!(reader.next().expect("read").expect("record").seq(), b"ACGT");
}

#[test]
fn unsupported_compression_is_named() {
    let cases: [(&str, &[u8], UnsupportedOperation); 2] = [
        ("x.fastq.bz2", b"BZh91AY&SY", UnsupportedOperation::Bzip2),
        (
            "x.fastq.xz",
            &[0xfd, b'7', b'z', b'X', b'Z', 0x00, 0x00, 0x00],
            UnsupportedOperation::Xz,
        ),
    ];
    for (name, magic, expected) in cases {
        let path = common::unique_path(name);
        common::write_plain(&path, magic);
        match FastqReader::from_path(&path) {
            Err(FastqError::Unsupported(op)) => assert_eq!(op, expected, "{name}"),
            Err(other) => panic!("{name}: unexpected error: {other}"),
            Ok(_) => panic!("{name}: must not be read as FASTQ"),
        }
    }
}

#[cfg(not(feature = "zstd"))]
#[test]
fn zstd_without_the_feature_says_so() {
    let path = common::unique_path("x.fastq.zst");
    common::write_plain(&path, &[0x28, 0xb5, 0x2f, 0xfd, 0, 0, 0, 0]);
    match FastqReader::from_path(&path) {
        Err(FastqError::Unsupported(UnsupportedOperation::Zstd)) => {}
        Err(other) => panic!("unexpected error: {other}"),
        Ok(_) => panic!("zstd input needs the feature"),
    }
}

#[cfg(feature = "zstd")]
#[test]
fn zstd_round_trips_with_the_feature() {
    let path = common::unique_path("round.fastq.zst");
    let mut writer = kira_fastq::FastqWriter::from_path(&path).expect("create");
    for i in 0..100 {
        writer
            .write_parts(format!("r{i}").as_bytes(), b"ACGT", b"!!!!")
            .expect("write");
    }
    writer.finish().expect("finish");

    let mut reader = FastqReader::from_path(&path).expect("open");
    let mut n = 0;
    while reader.next().expect("read").is_some() {
        n += 1;
    }
    assert_eq!(n, 100);
}

/// Compressed input on stdin is the norm for FASTQ pipelines.
#[test]
fn from_reader_auto_detects_gzip() {
    let path = common::unique_path("stdin.fastq.gz");
    common::write_multi_member_gzip(&path, &[b"@a\nAC\n+\n!!\n", b"@b\nGT\n+\n##\n"]);
    let file = std::fs::File::open(&path).expect("open");
    let mut reader = FastqReader::from_reader_auto(std::io::BufReader::new(file)).expect("sniff");
    let mut n = 0;
    while reader.next().expect("read").is_some() {
        n += 1;
    }
    assert_eq!(n, 2, "both members must be read");
}

#[test]
fn from_reader_auto_passes_plain_text_through() {
    let data: &[u8] = b"@a\nACGT\n+\n!!!!\n";
    let mut reader = FastqReader::from_reader_auto(std::io::BufReader::new(data)).expect("sniff");
    assert_eq!(reader.next().expect("read").expect("record").seq(), b"ACGT");
}

#[test]
fn from_vec_parses_in_memory() {
    let mut reader = FastqReader::from_vec(b"@a\nACGT\n+\n!!!!\n@b\nTT\n+\n##\n".to_vec());
    let mut n = 0;
    while reader.next().expect("read").is_some() {
        n += 1;
    }
    assert_eq!(n, 2);
    assert_eq!(reader.records_read(), 2);
}

#[test]
fn buffered_path_reads_the_same_records_as_mmap() {
    let path = common::unique_path("buffered.fastq");
    common::write_plain(&path, b"@a\nACGT\n+\n!!!!\n@b\nTT\n+\n##\n");
    let mut mapped = FastqReader::from_path(&path).expect("open");
    let mut buffered = FastqReader::from_path_buffered(&path).expect("open buffered");
    loop {
        let a = mapped.next().expect("read").map(|r| r.seq().to_vec());
        let b = buffered.next().expect("read").map(|r| r.seq().to_vec());
        assert_eq!(a, b);
        if a.is_none() {
            break;
        }
    }
}

#[test]
fn seek_is_refused_on_sources_that_cannot_do_it() {
    let path = common::unique_path("plainseek.fastq.gz");
    common::write_gzip(&path, b"@r1\nACGT\n+\n!!!!\n");
    let mut reader = FastqReader::from_path(&path).expect("open");
    assert!(matches!(
        reader.seek(VirtualOffset(0)),
        Err(FastqError::Unsupported(UnsupportedOperation::Seek))
    ));

    let mut stream = FastqReader::from_reader(std::io::BufReader::new(&b"@a\nAC\n+\n!!\n"[..]));
    assert!(matches!(
        stream.seek(VirtualOffset(0)),
        Err(FastqError::Unsupported(UnsupportedOperation::Seek))
    ));
}
