mod common;

use kira_fastq::backend::gzip::{GzipBackend, LineStatus};

#[test]
fn read_line_roundtrip() {
    let path = common::unique_path("rt.gz");
    common::write_gzip(&path, b"line one\nsecond line\n");
    let mut backend = GzipBackend::new(&path).expect("open");
    let mut out = Vec::new();
    assert_eq!(backend.read_line(&mut out).unwrap(), LineStatus::Line);
    assert_eq!(out, b"line one");
    assert_eq!(backend.read_line(&mut out).unwrap(), LineStatus::Line);
    assert_eq!(out, b"second line");
    assert_eq!(backend.read_line(&mut out).unwrap(), LineStatus::EofClean);
    assert!(out.is_empty());
}

#[test]
fn read_line_strips_cr() {
    let path = common::unique_path("crlf.gz");
    common::write_gzip(&path, b"crlf line\r\nlf line\n");
    let mut backend = GzipBackend::new(&path).expect("open");
    let mut out = Vec::new();
    assert_eq!(backend.read_line(&mut out).unwrap(), LineStatus::Line);
    assert_eq!(out, b"crlf line");
}

#[test]
fn read_line_eof_partial() {
    let path = common::unique_path("partial.gz");
    common::write_gzip(&path, b"unterminated");
    let mut backend = GzipBackend::new(&path).expect("open");
    let mut out = Vec::new();
    assert_eq!(backend.read_line(&mut out).unwrap(), LineStatus::EofPartial);
    assert_eq!(out, b"unterminated");
}

#[test]
fn multi_member() {
    let path = common::unique_path("multi.gz");
    common::write_multi_member_gzip(&path, &[b"first\n", b"second\n"]);
    let mut backend = GzipBackend::new(&path).expect("open");
    let mut out = Vec::new();
    assert_eq!(backend.read_line(&mut out).unwrap(), LineStatus::Line);
    assert_eq!(out, b"first");
    assert_eq!(backend.read_line(&mut out).unwrap(), LineStatus::Line);
    assert_eq!(out, b"second");
}
