use kira_fastq::canonical_read_id;

#[test]
fn strips_mate_suffix() {
    assert_eq!(canonical_read_id(b"r1/1"), b"r1");
    assert_eq!(canonical_read_id(b"r1/2"), b"r1");
    assert_eq!(
        canonical_read_id(b"HWUSI-EAS100R:6:73:941:1973#0/1"),
        b"HWUSI-EAS100R:6:73:941:1973#0"
    );
}

#[test]
fn keeps_prefix_before_space() {
    assert_eq!(
        canonical_read_id(b"M01234:23:000000000-A1BCD:1:1101:12345:6789 1:N:0:NNNN"),
        b"M01234:23:000000000-A1BCD:1:1101:12345:6789"
    );
}

#[test]
fn casava_18_pair_matches() {
    let a = b"M01234:1:1:1101:1:1 1:N:0:NNNN";
    let b = b"M01234:1:1:1101:1:1 2:N:0:NNNN";
    assert_eq!(canonical_read_id(a), canonical_read_id(b));
}

#[test]
fn illumina_pair_matches() {
    let a = b"HWUSI:6:73:941:1973#0/1";
    let b = b"HWUSI:6:73:941:1973#0/2";
    assert_eq!(canonical_read_id(a), canonical_read_id(b));
}

#[test]
fn mismatched_ids() {
    assert_ne!(canonical_read_id(b"r1"), canonical_read_id(b"r2"));
}
