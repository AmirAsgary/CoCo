//! Pins `KmerTranslator::kmer2packedKmer` against vectors captured from the C++.
//!
//! `testdata/kmer2packedkmer_golden.tsv` holds 1,206 `(pattern, kmer, packed)`
//! triples produced by the original binary. This is the first thing a port should
//! get right: the function is pure, sits in the hot path, and every later stage
//! is built on its output.

use coco::translator::KmerTranslator;

mod common;

fn testdata(name: &str) -> String {
    format!("{}/{}", common::testdata_dir(), name)
}

#[test]
fn kmer2packed_matches_the_cpp_golden_vectors() {
    let path = testdata("kmer2packedkmer_golden.tsv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let mut lines = text.lines();
    let header = lines.next().expect("golden file is empty");
    let cols: Vec<&str> = header.split('\t').collect();
    let ci = |name: &str| cols.iter().position(|c| *c == name).unwrap_or_else(|| panic!("no column {name}"));
    let (c_pattern, c_span, c_weight, c_kmer, c_packed) =
        (ci("pattern"), ci("span"), ci("weight"), ci("kmer_hex"), ci("packed_hex"));

    let mut checked = 0usize;
    let mut current: Option<(String, KmerTranslator)> = None;
    for (lineno, line) in lines.enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let pattern = f[c_pattern];
        let need_new = !matches!(&current, Some((p, _)) if p == pattern);
        if need_new {
            let t = KmerTranslator::new(pattern)
                .unwrap_or_else(|e| panic!("line {}: pattern {pattern} rejected: {e}", lineno + 2));
            current = Some((pattern.to_string(), t));
        }
        let (_, t) = current.as_ref().unwrap();

        assert_eq!(t.span() as usize, f[c_span].parse::<usize>().unwrap(), "span, line {}", lineno + 2);
        assert_eq!(t.weight() as usize, f[c_weight].parse::<usize>().unwrap(), "weight, line {}", lineno + 2);

        let kmer = u128::from_str_radix(f[c_kmer].trim_start_matches("0x"), 16).unwrap();
        let want = u64::from_str_radix(f[c_packed].trim_start_matches("0x"), 16).unwrap();
        let got = t.kmer2packed(kmer);
        assert_eq!(
            got, want,
            "line {}: pattern={pattern} kmer={kmer:#x}: got {got:#x}, want {want:#x}",
            lineno + 2
        );
        // The portable path must agree with whatever the dispatcher chose.
        assert_eq!(t.kmer2packed_scalar(kmer), want, "scalar path, line {}", lineno + 2);
        checked += 1;
    }
    assert!(checked > 1000, "expected the full golden set, only checked {checked}");
    eprintln!("checked {checked} golden kmer2packedKmer vectors");
}
