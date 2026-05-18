mod common;

use kira_fastq::{FastqError, FastqReader, InvalidKind};

#[test]
fn fastq_error_implements_std_error() {
    let err: FastqError = FastqError::InvalidFormat {
        offset: 42,
        kind: InvalidKind::HeaderMissingAt,
    };
    let s = format!("{err}");
    assert!(s.contains("invalid FASTQ format"));
    assert!(s.contains("42"));
    let dyn_err: &dyn std::error::Error = &err;
    let _ = dyn_err.source();
}

#[test]
fn io_error_source_chain() {
    let inner = std::io::Error::new(std::io::ErrorKind::PermissionDenied, "nope");
    let err: FastqError = inner.into();
    let dyn_err: &dyn std::error::Error = &err;
    assert!(dyn_err.source().is_some());
}

#[test]
fn try_for_each_iteration() {
    let path = common::unique_path("for_each.fastq");
    common::write_plain(&path, b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n");
    let mut reader = FastqReader::from_path(&path).expect("open");

    let mut total_bases: u64 = 0;
    reader
        .try_for_each::<std::io::Error, _>(|rec| {
            total_bases += rec.len() as u64;
            Ok(())
        })
        .expect("iter");
    assert_eq!(total_bases, 4 + 2);
}

#[test]
fn try_for_each_user_error_propagates() {
    let path = common::unique_path("for_each_err.fastq");
    common::write_plain(&path, b"@r1\nACGT\n+\n!!!!\n@r2\nTT\n+\n##\n");
    let mut reader = FastqReader::from_path(&path).expect("open");

    let res = reader.try_for_each::<&str, _>(|_| Err("stop"));
    match res {
        Err(kira_fastq::TryForEachError::User("stop")) => {}
        other => panic!("unexpected: {other:?}"),
    }
}
