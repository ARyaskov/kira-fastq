#![cfg(feature = "async")]

mod common;

use futures_util::StreamExt;
use kira_fastq::{AnyAsyncReader, AsyncFastqReader, AsyncFastqWriter, FastqRecord};
use tokio::io::BufReader;

#[tokio::test]
async fn async_read_plain_from_path() {
    let path = common::unique_path("async_in.fastq");
    common::write_plain(&path, b"@a\nACGT\n+\n!!!!\n@b\nTT\n+\n##\n");
    let mut reader = AnyAsyncReader::from_path(&path).await.expect("open");
    let mut n = 0u32;
    while let Some(_rec) = reader.next().await.expect("read") {
        n += 1;
    }
    assert_eq!(n, 2);
}

#[tokio::test]
async fn async_read_gzip_from_path() {
    let path = common::unique_path("async_in.fastq.gz");
    common::write_gzip(&path, b"@a\nACGT\n+\n!!!!\n@b\nGT\n+\n@@\n");
    let mut reader = AnyAsyncReader::from_path(&path).await.expect("open gz");
    let mut n = 0u32;
    while let Some(_rec) = reader.next().await.expect("read") {
        n += 1;
    }
    assert_eq!(n, 2);
}

#[tokio::test]
async fn async_read_from_buf_read() {
    let data: &[u8] = b"@r0\nACGT\n+\n!!!!\n";
    let buf = BufReader::new(data);
    let mut reader = AsyncFastqReader::from_reader(buf);
    let rec = reader.next().await.expect("read").expect("some");
    assert_eq!(rec.seq(), b"ACGT");
    assert!(reader.next().await.expect("eof").is_none());
}

#[tokio::test]
async fn async_records_stream_yields_owned() {
    let data: &[u8] = b"@a\nAA\n+\n!!\n@b\nGG\n+\n@@\n";
    let reader = AsyncFastqReader::from_reader(BufReader::new(data));
    let mut stream = reader.records();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    while let Some(rec) = stream.next().await {
        seqs.push(rec.expect("ok").seq().to_vec());
    }
    assert_eq!(seqs, vec![b"AA".to_vec(), b"GG".to_vec()]);
}

#[tokio::test]
async fn async_write_plain_roundtrip() {
    let out_path = common::unique_path("async_out.fastq");
    {
        let mut writer = AsyncFastqWriter::from_path(&out_path).await.expect("open");
        let rec = FastqRecord::new(b"r0", b"ACGT", b"!!!!");
        writer.write_record(&rec).await.expect("write");
        writer.shutdown().await.expect("shutdown");
    }
    let data = std::fs::read(&out_path).expect("read");
    assert_eq!(data, b"@r0\nACGT\n+\n!!!!\n");
}

#[tokio::test]
async fn async_write_gzip_roundtrip() {
    let out_path = common::unique_path("async_out.fastq.gz");
    {
        let mut writer = AsyncFastqWriter::from_path(&out_path)
            .await
            .expect("open gz");
        let rec = FastqRecord::new(b"r0", b"ACGT", b"!!!!");
        writer.write_record(&rec).await.expect("write");
        writer.shutdown().await.expect("shutdown");
    }
    // Decode back synchronously via kira's reader to confirm the gzip stream is valid.
    let mut reader = kira_fastq::FastqReader::from_path_auto(&out_path).expect("reopen gz");
    let rec = reader.next().expect("read").expect("some");
    assert_eq!(rec.seq(), b"ACGT");
}

/// bgzip and pigz produce concatenated gzip members. Decoding only the first one drops most of
/// the file, silently: the async reader has to be told to keep going.
#[tokio::test]
async fn async_reads_every_gzip_member() {
    let path = common::unique_path("async_multi.fastq.gz");
    common::write_multi_member_gzip(
        &path,
        &[b"@a\nAC\n+\n!!\n", b"@b\nGT\n+\n@@\n", b"@c\nTT\n+\n##\n"],
    );
    let mut reader = AnyAsyncReader::from_path(&path).await.expect("open");
    let mut n = 0u32;
    while reader.next().await.expect("read").is_some() {
        n += 1;
    }
    assert_eq!(n, 3, "every member must be decoded");
}

/// BGZF is a multi-member gzip file, so the same fix makes BGZF readable on the async path.
#[tokio::test]
async fn async_reads_bgzf_as_gzip() {
    let path = common::unique_path("async.bgz");
    let mut payload = Vec::new();
    for i in 0..500 {
        payload.extend_from_slice(format!("@r{i}\nACGTACGT\n+\n!!!!!!!!\n").as_bytes());
    }
    common::write_bgzf(&path, &payload, 4096);
    let mut reader = AnyAsyncReader::from_path(&path).await.expect("open");
    let mut n = 0u32;
    while reader.next().await.expect("read").is_some() {
        n += 1;
    }
    assert_eq!(n, 500);
}

#[tokio::test]
async fn async_accepts_the_same_deviations_as_the_sync_reader() {
    let data: &[u8] = b"\n@a\n\n+\n\n\n@b\nGT\n+\n##";
    let mut reader = AsyncFastqReader::from_reader(BufReader::new(data));
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    while let Some(rec) = reader.next().await.expect("read") {
        seqs.push(rec.seq().to_vec());
    }
    assert_eq!(seqs, vec![Vec::new(), b"GT".to_vec()]);
}

#[tokio::test]
async fn async_writer_rejects_a_record_it_cannot_write_back() {
    let path = common::unique_path("async_bad.fastq");
    let mut writer = AsyncFastqWriter::from_path(&path).await.expect("create");
    assert!(writer.write_parts(b"r0", b"ACGT", b"!!!").await.is_err());
    assert!(writer.write_parts(b"r0\nx", b"AC", b"!!").await.is_err());
    writer.shutdown().await.expect("shutdown");
}

#[tokio::test]
async fn async_errors_carry_the_record_index() {
    let data: &[u8] = b"@a\nAC\n+\n!!\n@b\nACGT\n+\n!!\n";
    let mut reader = AsyncFastqReader::from_reader(BufReader::new(data));
    reader.next().await.expect("read").expect("first");
    let err = reader.next().await.expect_err("second record is malformed");
    assert_eq!(err.record(), Some(2));
}
