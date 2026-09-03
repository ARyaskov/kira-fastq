#![cfg(feature = "noodles-bgzf")]

//! The `noodles-bgzf` adapter: same records, and virtual offsets that convert both ways.

mod common;

use kira_fastq::{FastqReader, VirtualOffset};

#[test]
fn reads_the_same_records_as_the_native_backend() {
    let path = common::unique_path("noodles.bgz");
    let mut payload = Vec::new();
    for i in 0..1000 {
        payload.extend_from_slice(format!("@r{i}\nACGTACGT\n+\n!!!!!!!!\n").as_bytes());
    }
    common::write_bgzf(&path, &payload, 4096);

    let mut native = FastqReader::from_bgzf_path(&path).expect("open native");
    let mut noodles = FastqReader::from_noodles_bgzf_path(&path).expect("open noodles");
    loop {
        let a = native.next().expect("read").map(|r| r.header().to_vec());
        let b = noodles.next().expect("read").map(|r| r.header().to_vec());
        assert_eq!(a, b);
        if a.is_none() {
            break;
        }
    }
}

/// The adapter exists so offsets can travel between this crate and the noodles ecosystem.
#[test]
fn virtual_offsets_round_trip_and_convert() {
    let path = common::unique_path("noodles_seek.bgz");
    common::write_bgzf(&path, b"@r1\nACGT\n+\n!!!!\n@r2\nTTTT\n+\n####\n", 20);

    let mut reader = FastqReader::from_noodles_bgzf_path(&path).expect("open");
    reader.next().expect("read").expect("first");
    let checkpoint = reader.tell();
    let second = reader
        .next()
        .expect("read")
        .expect("second")
        .header()
        .to_vec();

    reader
        .seek(checkpoint)
        .expect("noodles adapter supports seek");
    assert_eq!(
        reader.next().expect("read").expect("second again").header(),
        second.as_slice()
    );

    let as_noodles: noodles_bgzf::VirtualPosition = checkpoint.into();
    let back: VirtualOffset = as_noodles.into();
    assert_eq!(back, checkpoint);
}

/// An independent implementation must accept what this crate's native BGZF writer produces,
/// which is the real test of the block framing.
#[test]
fn noodles_reads_our_native_bgzf_output() {
    let path = common::unique_path("native_for_noodles.bgz");
    let mut writer = kira_fastq::FastqWriter::to_bgzf_path(&path, 6).expect("create");
    for i in 0..5_000 {
        writer
            .write_parts(format!("r{i}").as_bytes(), b"ACGTACGTAC", b"!!!!!!!!!!")
            .expect("write");
    }
    writer.finish().expect("finish");

    // Decode with noodles rather than with our own reader.
    let file = std::fs::File::open(&path).expect("open");
    let mut decoded = Vec::new();
    let mut reader = noodles_bgzf::io::Reader::new(std::io::BufReader::new(file));
    std::io::Read::read_to_end(&mut reader, &mut decoded).expect("noodles decode");
    let mut expected = Vec::new();
    for i in 0..5_000 {
        expected.extend_from_slice(format!("@r{i}\nACGTACGTAC\n+\n!!!!!!!!!!\n").as_bytes());
    }
    assert_eq!(
        decoded, expected,
        "noodles must decode our blocks byte for byte"
    );

    // And the same file read back through our reader agrees.
    let mut ours = FastqReader::from_bgzf_path(&path).expect("reopen");
    let mut n = 0;
    while ours.next().expect("read").is_some() {
        n += 1;
    }
    assert_eq!(n, 5_000);
}
