//! Pins the count table against an independently computed answer.
//!
//! `dsk2ascii` gives the 41-mers and their abundances. This test translates each
//! one to its spaced packed k-mer with an ordinary `HashMap`, accumulating exactly
//! the way `Options::COUNT_MODE_SUM` says to, and then requires the real table to
//! agree on every key. Nothing here reuses the table's own indexing, so a bug in
//! the grid arithmetic or the open addressing cannot hide.

mod common;

use coco::kmer::string_to_packed_kmer;
use coco::lookuptable::CountMode;
use coco::preprocessing::{build_lookuptable, TableKind};
use coco::translator::{KmerTranslator, DEFAULT_PATTERN};
use std::collections::HashMap;

/// Expected `spaced packed k-mer -> summed abundance`, computed from the text dump.
fn expected_counts(t: &KmerTranslator) -> HashMap<u64, u32> {
    let text = std::fs::read_to_string(common::solid_ascii()).unwrap();
    let mut m: HashMap<u64, u32> = HashMap::with_capacity(400_000);
    for line in text.lines() {
        let mut it = line.split_whitespace();
        let kmer = it.next().unwrap();
        let ab: u32 = it.next().unwrap().parse().unwrap();
        // A 41-mer packs into 82 bits, so build the span-wide value as u128.
        let mut v: u128 = 0;
        for &c in kmer.as_bytes() {
            v = (v << 2) | coco::types::res2int(c) as u128;
        }
        *m.entry(t.kmer2min_packed(v)).or_insert(0) += ab;
    }
    m
}

#[test]
fn compact_table_counts_match_an_independent_computation() {
    let t = KmerTranslator::new(DEFAULT_PATTERN).unwrap();
    let expect = expected_counts(&t);
    let table = build_lookuptable(
        &common::counts_file(),
        CountMode::Sum,
        &t,
        0,
        TableKind::Compact,
        4,
        &mut |_| {},
    )
    .expect("build lookup table");

    assert_eq!(table.len(), expect.len(), "number of distinct spaced k-mers");
    for (&k, &want) in &expect {
        assert_eq!(table.get_count(k), want, "count for packed kmer {k:#x}");
    }

    // Keys that are not present must read as zero, not as something stale.
    let mut absent = 0;
    for probe in 0..2000u64 {
        let k = probe.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 1;
        if !expect.contains_key(&k) {
            assert_eq!(table.get_count(k), 0, "absent kmer {k:#x} must read as 0");
            absent += 1;
        }
    }
    assert!(absent > 1000, "probe set was not actually absent");
    eprintln!("verified {} spaced k-mer counts", expect.len());
}

#[test]
fn spaced_collisions_really_happen_so_the_sum_path_is_exercised() {
    // If no two 41-mers ever collided under the spaced mask, COUNT_MODE_SUM would
    // be untested by the check above. Confirm the test set does exercise it.
    let t = KmerTranslator::new(DEFAULT_PATTERN).unwrap();
    let text = std::fs::read_to_string(common::solid_ascii()).unwrap();
    let mut seen: HashMap<u64, u32> = HashMap::with_capacity(400_000);
    for line in text.lines() {
        let kmer = line.split_whitespace().next().unwrap();
        let mut v: u128 = 0;
        for &c in kmer.as_bytes() {
            v = (v << 2) | coco::types::res2int(c) as u128;
        }
        *seen.entry(t.kmer2min_packed(v)).or_insert(0) += 1;
    }
    let collisions: usize = seen.values().filter(|&&n| n > 1).count();
    eprintln!("{collisions} spaced k-mers are shared by more than one 41-mer");
    assert!(collisions > 0, "no collisions: COUNT_MODE_SUM would be untested here");
}

#[test]
fn grid_and_compact_tables_agree_on_the_real_dataset() {
    let t = KmerTranslator::new(DEFAULT_PATTERN).unwrap();
    let counts = common::counts_file();
    let compact =
        build_lookuptable(&counts, CountMode::Sum, &t, 0, TableKind::Compact, 4, &mut |_| {}).unwrap();
    let grid =
        build_lookuptable(&counts, CountMode::Sum, &t, 0, TableKind::Grid, 1, &mut |_| {}).unwrap();
    assert_eq!(compact.len(), grid.len());

    let expect = expected_counts(&t);
    for &k in expect.keys() {
        assert_eq!(compact.get_count(k), grid.get_count(k), "kmer {k:#x}");
    }
}

#[test]
fn string_to_packed_kmer_is_the_inverse_of_the_decoder() {
    for s in [&b"ACGTACGTACGT"[..], b"AAAAAAAAAAAA", b"GTGTGTGTGTGT"] {
        let k = string_to_packed_kmer(s).unwrap();
        assert_eq!(coco::kmer::packed_kmer_to_string(k, s.len() as u16).as_bytes(), s);
    }
    assert!(string_to_packed_kmer(b"ACGN").is_none());
}
