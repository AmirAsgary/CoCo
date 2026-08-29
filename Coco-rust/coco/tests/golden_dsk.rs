//! Pins the HDF5 reader against `dsk2ascii`, DSK's own text dump.
//!
//! This is the one place the port replaces a whole third-party library
//! (GATB-core's `Storage`/`Partition<Count>`), so it needs an independent ground
//! truth rather than a comparison against itself. `dsk2ascii` writes
//! `<kmer> <abundance>` in the same partition order the C++ iterator uses, which
//! pins the record decoding, the 2-bit base order and the iteration order all at
//! once.

mod common;

use coco::dsk::DskCountsFile;
use coco::types::INT2RES;

/// Decode a GATB `LargeInt` k-mer, most significant base first.
fn kmer_to_string(value: u128, k: usize) -> String {
    let mut out = vec![0u8; k];
    let mut v = value;
    for i in 0..k {
        out[k - 1 - i] = INT2RES[(v & 3) as usize];
        v >>= 2;
    }
    String::from_utf8(out).unwrap()
}

#[test]
fn reader_matches_dsk2ascii_exactly() {
    let counts = common::counts_file();
    let ascii = common::solid_ascii();

    let dsk = DskCountsFile::open(&counts).expect("open counts file");
    assert_eq!(dsk.kmer_size, 41, "the default pattern has span 41");

    let text = std::fs::read_to_string(&ascii).expect("read dsk2ascii output");
    let mut expected = text.lines();

    let mut n = 0usize;
    let mut mismatch: Option<String> = None;
    dsk.for_each(|c| {
        if mismatch.is_some() {
            return;
        }
        let Some(line) = expected.next() else {
            mismatch = Some(format!("reader produced more records than dsk2ascii ({n} so far)"));
            return;
        };
        let mut it = line.split_whitespace();
        let want_kmer = it.next().unwrap();
        let want_ab: u32 = it.next().unwrap().parse().unwrap();
        let got_kmer = kmer_to_string(c.value, 41);
        if got_kmer != want_kmer || c.abundance != want_ab {
            mismatch = Some(format!(
                "record {n}: got ({got_kmer}, {}), want ({want_kmer}, {want_ab})",
                c.abundance
            ));
        }
        n += 1;
    })
    .expect("iterate solid kmers");

    if let Some(m) = mismatch {
        panic!("{m}");
    }
    assert!(expected.next().is_none(), "reader produced fewer records than dsk2ascii");
    assert_eq!(n, 395_391, "solid k-mer count from the testdata README");
    eprintln!("verified {n} solid k-mers against dsk2ascii");
}

#[test]
fn item_count_matches_the_stream() {
    let dsk = DskCountsFile::open(&common::counts_file()).unwrap();
    let declared = dsk.num_items().unwrap();
    let mut seen = 0usize;
    dsk.for_each(|_| seen += 1).unwrap();
    assert_eq!(declared, seen, "num_items must agree with what for_each yields");
    assert_eq!(dsk.num_partitions(), 8);
}
