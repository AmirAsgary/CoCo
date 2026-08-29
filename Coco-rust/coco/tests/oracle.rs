//! Function-by-function comparison against the C++, via the oracle harness.
//!
//! End-to-end diffing catches whatever the pipeline happens to reach. Several
//! functions are not reachable that way at all: `calcNeighborhoodTolerance` is
//! file-static, `checkForSpuriousTransitionDrops` and
//! `checkForSpuriousTransitionDropsWithWindow` have no caller in the released
//! code, and the indel helpers only fire on ragged coverage. `Coco-rust/oracle`
//! is a small C++ binary that calls them directly and prints their results; these
//! tests feed it the same inputs as the Rust code and require identical output.

mod common;

use coco::countprofile::{CountProfile, TableAccess};
use coco::kmer::{min_index, packed_kmer_to_string, rev_complement};
use coco::lookuptable::{CompactLookupTable, CountMode};
use coco::seq::{get_avg_qual, SequenceInfo};
use coco::translator::KmerTranslator;
use std::io::Write;
use std::process::{Command, Stdio};

const DEFAULT_PATTERN: &str = coco::translator::DEFAULT_PATTERN;
/// A shorter symmetric pattern, so tests also cover a span the default never uses.
const SHORT_PATTERN: &str = "1111011011010101101101111";

fn oracle_bin() -> String {
    let p = std::env::var("COCO_ORACLE")
        .unwrap_or_else(|_| format!("{}/../oracle/build/oracle", common::crate_dir()));
    assert!(
        std::path::Path::new(&p).exists(),
        "missing oracle binary {p}\nBuild it with:\n  \
         cd Coco-rust/oracle && mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=RELEASE .. && make -j"
    );
    p
}

/// Run the oracle with `args`, feeding `input` on stdin, and return its stdout lines.
///
/// The input goes through a temporary file rather than a pipe. Writing it to a
/// pipe and only then reading stdout deadlocks as soon as the oracle's output
/// exceeds the 64 KB pipe buffer: it stops consuming stdin while blocked on its
/// own write, and this side blocks writing the rest of the input.
fn oracle(args: &[&str], input: &str) -> Vec<String> {
    static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("coco_oracle_in_{}_{n}", std::process::id()));
    {
        let mut f = std::fs::File::create(&path).expect("create oracle input");
        f.write_all(input.as_bytes()).unwrap();
        f.flush().unwrap();
    }
    let stdin = std::fs::File::open(&path).expect("open oracle input");
    let out = Command::new(oracle_bin())
        .args(args)
        .stdin(Stdio::from(stdin))
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .output()
        .expect("run oracle");
    let _ = std::fs::remove_file(&path);
    assert!(out.status.success(), "oracle {args:?} failed");
    String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .map(|s| s.to_string())
        .collect()
}

struct Rng(u64);
impl Rng {
    fn new(seed: u64) -> Self {
        Rng(seed | 1)
    }
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn dna(&mut self, len: usize) -> String {
        (0..len).map(|_| b"ACGT"[self.below(4) as usize] as char).collect()
    }
}

// ---------------------------------------------------------------- kmer.cpp ---

#[test]
fn rev_complement_and_min_index_match_the_cpp() {
    let mut rng = Rng::new(0xC0C0_1234);
    for k in [12u16, 20, 31, 32] {
        let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
        let values: Vec<u64> = (0..3000).map(|_| rng.next() & mask).collect();
        let input: String = values.iter().map(|v| format!("{v:016x}\n")).collect();
        let ks = k.to_string();

        let want = oracle(&["revcomp", &ks], &input);
        for (v, w) in values.iter().zip(&want) {
            let got = rev_complement(*v, k);
            assert_eq!(format!("{got:016x}"), *w, "revComplement k={k} v={v:#x}");
        }

        let want = oracle(&["minindex", &ks], &input);
        for (v, w) in values.iter().zip(&want) {
            let got = min_index(*v, k);
            assert_eq!(format!("{got:016x}"), *w, "minIndex k={k} v={v:#x}");
        }

        let want = oracle(&["packed2string", &ks], &input);
        for (v, w) in values.iter().zip(&want) {
            assert_eq!(packed_kmer_to_string(*v, k), *w, "packedKmer2String k={k}");
        }
    }
}

// ------------------------------------------------------- KmerTranslator.cpp ---

#[test]
fn translator_geometry_matches_the_cpp() {
    for pattern in [DEFAULT_PATTERN, SHORT_PATTERN, &"1".repeat(32), &"1".repeat(12)] {
        let t = KmerTranslator::new(pattern).unwrap_or_else(|e| panic!("{pattern}: {e}"));
        let lines = oracle(&["translator", pattern], "");
        let mut map = std::collections::HashMap::new();
        for l in &lines {
            let mut it = l.splitn(2, '\t');
            map.insert(it.next().unwrap().to_string(), it.next().unwrap_or("").to_string());
        }
        assert_eq!(map["span"], t.span().to_string(), "span for {pattern}");
        assert_eq!(map["weight"], t.weight().to_string(), "weight for {pattern}");
        assert_eq!(map["longestBlock"], t.longest_block().to_string(), "longestBlock for {pattern}");
        let (li, lo) = t.best_split();
        assert_eq!(map["logIndexSize"], li.to_string(), "logIndexSize for {pattern}");
        assert_eq!(map["logOffsetSize"], lo.to_string(), "logOffsetSize for {pattern}");
        let mask: Vec<String> = t.mask_array().iter().map(|v| v.to_string()).collect();
        assert_eq!(map["maskArray"], mask.join("\t"), "maskArray for {pattern}");
        let inv: Vec<String> = t.inverse_mask_array().iter().map(|v| v.to_string()).collect();
        assert_eq!(map["inverseMaskArray"], inv.join("\t"), "inverseMaskArray for {pattern}");
    }
}

// ------------------------------------------------------- SequenceInfo.h ------

#[test]
fn get_avg_qual_matches_the_cpp() {
    let mut rng = Rng::new(0xA11C);
    let mut input = String::new();
    let mut cases = Vec::new();
    for _ in 0..2000 {
        let len = 2 + rng.below(60) as usize;
        let q: String = (0..len).map(|_| (33 + rng.below(60) as u8) as char).collect();
        let pos = rng.below(len as u64) as usize;
        input.push_str(&format!("{q} {pos}\n"));
        cases.push((q, pos));
    }
    let want = oracle(&["avgqual"], &input);
    for ((q, pos), w) in cases.iter().zip(&want) {
        let got = get_avg_qual(q.as_bytes(), *pos);
        assert_eq!(got.to_string(), *w, "getAvgQual({q:?}, {pos})");
    }
}

// ------------------------------------------------------- filehandling.cpp ----

#[test]
fn path_helpers_match_the_cpp() {
    // getFilename strips the extension *before* taking the basename, which does
    // surprising things to paths whose directory contains a dot. Rather than
    // reason about basename(3)'s corner cases, ask the C++.
    let paths = [
        "/a/b/reads.fq",
        "reads.err0.1pct.fq",
        "noext",
        "/a/b/noext",
        "/a.b/c",
        "a.b/c.d",
        ".fq",
        "/",
        "x.",
        "./rel.fq",
        "dir/",
        "dir.d/",
    ];
    let input: String = paths.iter().map(|p| format!("{p}\n")).collect();

    let want = oracle(&["filename"], &input);
    for (p, w) in paths.iter().zip(&want) {
        assert_eq!(coco::filehandling::get_filename(p), *w, "getFilename({p:?})");
    }

    let want = oracle(&["fileext"], &input);
    for (p, w) in paths.iter().zip(&want) {
        assert_eq!(coco::filehandling::get_file_extension(p), *w, "getFileExtension({p:?})");
    }

    // endsWith takes the suffix first. "@" stands for the empty string, since the
    // protocol is whitespace separated.
    let cases: [(&str, &str); 7] = [
        (".gz", "reads.fq.gz"),
        (".gz", "reads.fq"),
        (".fastq", "a"),
        ("@", "anything"),
        (".fq", ".fq"),
        (".fq", "@"),
        ("a", "aa"),
    ];
    let input: String = cases.iter().map(|(s, t)| format!("{s} {t}\n")).collect();
    let want = oracle(&["endswith"], &input);
    for ((s, t), w) in cases.iter().zip(&want) {
        let (su, st) = (if *s == "@" { "" } else { *s }, if *t == "@" { "" } else { *t });
        assert_eq!(
            (coco::filehandling::ends_with(su, st) as u8).to_string(),
            *w,
            "endsWith({su:?}, {st:?})"
        );
    }
}

// ------------------------------------------------------- CountProfile.cpp ----

/// A profile line for the oracle: `-1` marks an invalid position.
fn profile_line(vals: &[i64]) -> String {
    vals.iter().map(|v| v.to_string()).collect::<Vec<_>>().join(" ") + "\n"
}

/// Install a profile directly, matching what the oracle's `setProfile` does.
fn make_profile<'a>(
    t: &'a KmerTranslator,
    table: &'a dyn coco::lookuptable::LookupTable,
    vals: &[i64],
) -> CountProfile<'a> {
    let mut cp = CountProfile::new(t, TableAccess::shared(table));
    cp.set_profile_for_test(vals);
    cp
}

fn empty_table() -> CompactLookupTable {
    CompactLookupTable::with_capacity(16, 30, 34, CountMode::Sum)
}

fn random_profiles(rng: &mut Rng, count: usize, min_len: usize, max_len: usize) -> Vec<Vec<i64>> {
    (0..count)
        .map(|_| {
            let len = min_len + rng.below((max_len - min_len + 1) as u64) as usize;
            // A plateau with dips, which is what a real count profile looks like.
            let base = 5 + rng.below(80) as i64;
            (0..len)
                .map(|_| match rng.below(10) {
                    0 => 0,
                    1 => rng.below(3) as i64,
                    2 if rng.below(4) == 0 => -1, // invalid position
                    _ => (base + rng.below(11) as i64 - 5).max(0),
                })
                .collect()
        })
        .collect()
}

#[test]
fn maximize_matches_the_cpp() {
    let mut rng = Rng::new(0x4242);
    for pattern in [DEFAULT_PATTERN, SHORT_PATTERN] {
        let t = KmerTranslator::new(pattern).unwrap();
        let table = empty_table();
        let profiles = random_profiles(&mut rng, 400, 1, 160);
        let input: String = profiles.iter().map(|p| profile_line(p)).collect();
        let want = oracle(&["maximize", pattern], &input);
        assert_eq!(want.len(), profiles.len());
        for (p, w) in profiles.iter().zip(&want) {
            let mut cp = make_profile(&t, &table, p);
            let got: Vec<String> = cp.maximize().iter().map(|v| v.to_string()).collect();
            assert_eq!(got.join(" "), *w, "maximize {pattern} on {p:?}");
        }
    }
}

#[test]
fn neighborhood_tolerance_matches_the_cpp() {
    let mut rng = Rng::new(0x9999);
    for pattern in [DEFAULT_PATTERN, SHORT_PATTERN] {
        let t = KmerTranslator::new(pattern).unwrap();
        let span = t.span() as usize;
        for &(threshold, pseudo, lower) in
            &[(0.01f64, 1u32, 5u32), (0.05, 3, 2), (0.002, 0, 20), (0.33, 1, 0)]
        {
            // The tolerance window reads maxProfile[0..=span], so the input must
            // be long enough for the C++ not to read past its own array. It also
            // takes a *maximized* profile, which is never negative -- feeding the
            // -1 invalid marker here would be read as UINT32_MAX on the C++ side.
            let profiles: Vec<Vec<i64>> = random_profiles(&mut rng, 120, span + 2, span + 140)
                .into_iter()
                .map(|p| p.into_iter().map(|v| v.max(0)).collect())
                .collect();
            let input: String = profiles.iter().map(|p| profile_line(p)).collect();
            let args = [
                "tolerance",
                pattern,
                &threshold.to_string(),
                &pseudo.to_string(),
                &lower.to_string(),
            ];
            let want = oracle(&args, &input);
            for (p, w) in profiles.iter().zip(&want) {
                let got = coco::countprofile::neighborhood_tolerance_for_test(
                    &t,
                    &p.iter().map(|&v| v.max(0) as u32).collect::<Vec<u32>>(),
                    threshold,
                    pseudo,
                    lower,
                );
                let got: Vec<String> = got.iter().map(|v| v.to_string()).collect();
                assert_eq!(
                    got.join(" "),
                    *w,
                    "calcNeighborhoodTolerance {pattern} t={threshold} p={pseudo} l={lower}"
                );
            }
        }
    }
}

#[test]
fn quantile_matches_the_cpp() {
    let mut rng = Rng::new(0x7777);
    let t = KmerTranslator::new(DEFAULT_PATTERN).unwrap();
    let table = empty_table();
    for q in [0.5f64, 0.67, 0.0, 0.9] {
        let profiles = random_profiles(&mut rng, 300, 1, 200);
        let input: String = profiles.iter().map(|p| profile_line(p)).collect();
        let want = oracle(&["quantile", &q.to_string()], &input);
        for (p, w) in profiles.iter().zip(&want) {
            let cp = make_profile(&t, &table, p);
            assert_eq!(cp.calc_x_quantile(q, &[]).to_string(), *w, "calcXquantile({q})");
        }
    }
}

#[test]
fn spurious_transition_drop_filter_matches_the_cpp() {
    let mut rng = Rng::new(0x1357);
    for pattern in [DEFAULT_PATTERN, SHORT_PATTERN] {
        let t = KmerTranslator::new(pattern).unwrap();
        let table = empty_table();
        let window = (t.longest_block() as usize + 1).max(5);
        for threshold in [0.1f64, 0.3, 0.5] {
            // Below 2*window the C++ reads past the end of its profile array, so
            // stay above it; that boundary is checked separately in the unit tests.
            let profiles = random_profiles(&mut rng, 250, 2 * window + 1, 2 * window + 150);
            let input: String = profiles.iter().map(|p| profile_line(p)).collect();
            let want = oracle(&["dropsnew", pattern, &threshold.to_string()], &input);
            for (p, w) in profiles.iter().zip(&want) {
                let cp = make_profile(&t, &table, p);
                let got = cp.check_for_spurious_transition_drops_with_window_new(threshold);
                assert_eq!(
                    (got as u8).to_string(),
                    *w,
                    "checkForSpuriousTransitionDropsWithWindowNew {pattern} thr={threshold} on {p:?}"
                );
            }
        }
    }
}

#[test]
fn legacy_drop_filters_match_the_cpp() {
    // Neither of these has a caller in the released C++; the oracle is the only
    // way to exercise them at all. checkForSpuriousTransitionDrops divides by the
    // previous position's count, so its inputs must stay strictly positive or the
    // C++ takes SIGFPE.
    let mut rng = Rng::new(0x8642);
    for pattern in [DEFAULT_PATTERN, SHORT_PATTERN] {
        let t = KmerTranslator::new(pattern).unwrap();
        let span = t.span() as usize;
        let table = empty_table();

        let positive: Vec<Vec<i64>> = random_profiles(&mut rng, 200, span + 10, span + 150)
            .into_iter()
            .map(|p| p.into_iter().map(|v| v.max(1)).collect())
            .collect();
        let input: String = positive.iter().map(|p| profile_line(p)).collect();

        for &(cov_est, local, global, mask_edges) in &[
            (40u32, 0.33f64, 0.33f64, 1i32),
            (40, 0.33, 0.33, 0),
            (10, 0.5, 0.2, 1),
            (100, 0.1, 0.5, 0),
        ] {
            let args = [
                "dropswindow",
                pattern,
                &cov_est.to_string(),
                &local.to_string(),
                &global.to_string(),
                &mask_edges.to_string(),
            ];
            let want = oracle(&args, &input);
            for (p, w) in positive.iter().zip(&want) {
                let mut cp = make_profile(&t, &table, p);
                let mx: Vec<u32> = cp.maximize().to_vec();
                let got = cp.check_for_spurious_transition_drops_with_window(
                    &mx, cov_est, local, global, mask_edges != 0,
                );
                assert_eq!(
                    (got as u8).to_string(), *w,
                    "checkForSpuriousTransitionDropsWithWindow {pattern} cov={cov_est} \
                     local={local} global={global} maskEdges={mask_edges} on {p:?}"
                );
            }
        }

        for &(drop_level, mask_edges) in &[(5u32, 1i32), (5, 0), (20, 1), (1, 0)] {
            let args = [
                "dropsold",
                pattern,
                &drop_level.to_string(),
                &mask_edges.to_string(),
            ];
            let want = oracle(&args, &input);
            for (p, w) in positive.iter().zip(&want) {
                let mut cp = make_profile(&t, &table, p);
                let mx: Vec<u32> = cp.maximize().to_vec();
                let got = cp.check_for_spurious_transition_drops(&mx, drop_level, mask_edges != 0);
                assert_eq!(
                    (got as u8).to_string(), *w,
                    "checkForSpuriousTransitionDrops {pattern} dropLevel={drop_level} \
                     maskEdges={mask_edges} on {p:?}"
                );
            }
        }
    }
}

// ---------------------------------------------- correction on a known table ---

/// Build the count table the oracle builds: every spaced k-mer of `reference`,
/// canonicalised, with weight `mult`.
fn table_from_reference(t: &KmerTranslator, reference: &str, mult: u32) -> CompactLookupTable {
    let mut table = CompactLookupTable::with_capacity(reference.len(), 30, 34, CountMode::Sum);
    let span = t.span() as usize;
    let bytes = reference.as_bytes();
    if bytes.len() < span {
        return table;
    }
    let mut kmer: u128 = 0;
    let mut n_store: u128 = 0;
    for idx in 0..bytes.len() {
        let code = coco::types::res2int(bytes[idx]);
        if code != -1 {
            kmer = (kmer << 2) | code as u128;
            n_store <<= 1;
        } else {
            kmer <<= 2;
            n_store = (n_store << 1) | 1;
        }
        if idx + 1 >= span {
            if (n_store & t.spaced_mask()) != 0 {
                continue;
            }
            table.add_element(t.kmer2min_packed(kmer), mult);
        }
    }
    table
}

/// `(reference, read, qual, multiplicity)` cases with substitutions and indels.
fn correction_cases(rng: &mut Rng, n: usize, span: usize, with_indels: bool) -> Vec<(String, String, String, u32)> {
    let mut out = Vec::new();
    while out.len() < n {
        let ref_len = span * 4 + rng.below(400) as usize;
        let reference = rng.dna(ref_len);
        let read_len = (span + 20 + rng.below(120) as usize).min(ref_len);
        let start = rng.below((ref_len - read_len + 1) as u64) as usize;
        let mut read: Vec<u8> = reference.as_bytes()[start..start + read_len].to_vec();

        // Inject errors.
        let n_err = 1 + rng.below(4) as usize;
        for _ in 0..n_err {
            if read.len() <= span + 2 {
                break;
            }
            let p = rng.below(read.len() as u64) as usize;
            match if with_indels { rng.below(3) } else { 0 } {
                0 => {
                    let cur = read[p];
                    let mut nb = b"ACGT"[rng.below(4) as usize];
                    while nb == cur {
                        nb = b"ACGT"[rng.below(4) as usize];
                    }
                    read[p] = nb;
                }
                1 => {
                    read.insert(p, b"ACGT"[rng.below(4) as usize]);
                }
                _ => {
                    read.remove(p);
                }
            }
        }
        if read.len() < span {
            continue;
        }
        let read = String::from_utf8(read).unwrap();
        let qual: String = (0..read.len()).map(|_| (33 + rng.below(41) as u8) as char).collect();
        let mult = 5 + rng.below(60) as u32;
        out.push((reference, read, qual, mult));
    }
    out
}

fn run_correction_oracle(
    mode: &str,
    pattern: &str,
    threshold: f64,
    pseudo: u32,
    lower: u32,
    flag: i32,
    cases: &[(String, String, String, u32)],
) -> Vec<String> {
    let input: String = cases
        .iter()
        .map(|(r, rd, q, m)| format!("{r} {rd} {q} {m}\n"))
        .collect();
    oracle(
        &[
            mode,
            pattern,
            &threshold.to_string(),
            &pseudo.to_string(),
            &lower.to_string(),
            &flag.to_string(),
        ],
        &input,
    )
}

#[test]
fn substitution_correction_matches_the_cpp() {
    let mut rng = Rng::new(0x2024);
    for pattern in [DEFAULT_PATTERN, SHORT_PATTERN] {
        let t = KmerTranslator::new(pattern).unwrap();
        let span = t.span() as usize;
        let cases = correction_cases(&mut rng, 300, span, false);
        for &(threshold, pseudo, lower) in &[(0.01f64, 1u32, 5u32), (0.05, 3, 2)] {
            for &multi in &[1i32, 0] {
                let want =
                    run_correction_oracle("correct", pattern, threshold, pseudo, lower, multi, &cases);
                for (i, (reference, read, qual, mult)) in cases.iter().enumerate() {
                    let table = table_from_reference(&t, reference, *mult);
                    let mut cp = CountProfile::new(&t, TableAccess::shared(&table));
                    let mut si = SequenceInfo {
                        name: b"t".to_vec(),
                        comment: Vec::new(),
                        seq: read.as_bytes().to_vec(),
                        qual: qual.as_bytes().to_vec(),
                        sep: b'@',
                    };
                    cp.fill(&si);
                    cp.maximize();
                    let mut sub = 0u32;
                    let status = cp.do_substitution_correction(
                        &mut si, threshold, pseudo, lower, multi != 0, false, &mut sub,
                    );
                    let got = format!(
                        "{}\t{}\t{}\t0\t0\t{}",
                        String::from_utf8_lossy(&si.seq),
                        if si.qual.is_empty() { "-".into() } else { String::from_utf8_lossy(&si.qual).to_string() },
                        sub,
                        status
                    );
                    assert_eq!(
                        got, want[i],
                        "doSubstitutionCorrection {pattern} thr={threshold} multi={multi} case {i}\n  ref={reference}\n  read={read}"
                    );
                }
            }
        }
    }
}

#[test]
fn indel_correction_matches_the_cpp() {
    let mut rng = Rng::new(0x2025);
    for pattern in [DEFAULT_PATTERN, SHORT_PATTERN] {
        let t = KmerTranslator::new(pattern).unwrap();
        let span = t.span() as usize;
        let cases = correction_cases(&mut rng, 400, span, true);
        for &(threshold, pseudo, lower) in &[(0.01f64, 1u32, 5u32), (0.05, 3, 2)] {
            for &trysub in &[1i32, 0] {
                let want =
                    run_correction_oracle("indel", pattern, threshold, pseudo, lower, trysub, &cases);
                for (i, (reference, read, qual, mult)) in cases.iter().enumerate() {
                    let table = table_from_reference(&t, reference, *mult);
                    let mut cp = CountProfile::new(&t, TableAccess::shared(&table));
                    let mut si = SequenceInfo {
                        name: b"t".to_vec(),
                        comment: Vec::new(),
                        seq: read.as_bytes().to_vec(),
                        qual: qual.as_bytes().to_vec(),
                        sep: b'@',
                    };
                    cp.fill(&si);
                    cp.maximize();
                    let (mut sub, mut ins, mut del) = (0u32, 0u32, 0u32);
                    cp.do_indel_correction(
                        &mut si, threshold, pseudo, lower, trysub != 0, false,
                        &mut sub, &mut ins, &mut del,
                    );
                    let got = format!(
                        "{}\t{}\t{sub}\t{ins}\t{del}",
                        String::from_utf8_lossy(&si.seq),
                        if si.qual.is_empty() { "-".into() } else { String::from_utf8_lossy(&si.qual).to_string() },
                    );
                    assert_eq!(
                        got, want[i],
                        "doIndelCorrection {pattern} thr={threshold} trysub={trysub} case {i}\n  ref={reference}\n  read={read}"
                    );
                }
            }
        }
    }
}

#[test]
fn trimming_matches_the_cpp() {
    let mut rng = Rng::new(0x2026);
    for pattern in [DEFAULT_PATTERN, SHORT_PATTERN] {
        let t = KmerTranslator::new(pattern).unwrap();
        let span = t.span() as usize;
        let cases = correction_cases(&mut rng, 250, span, true);
        for &max_trim in &[5i32, 20] {
            let want = run_correction_oracle("trim", pattern, 0.01, 1, 5, max_trim, &cases);
            for (i, (reference, read, qual, mult)) in cases.iter().enumerate() {
                let table = table_from_reference(&t, reference, *mult);
                let mut cp = CountProfile::new(&t, TableAccess::shared(&table));
                let mut si = SequenceInfo {
                    name: b"t".to_vec(),
                    comment: Vec::new(),
                    seq: read.as_bytes().to_vec(),
                    qual: qual.as_bytes().to_vec(),
                    sep: b'@',
                };
                cp.fill(&si);
                cp.maximize();
                let mut trimmed = 0u32;
                cp.do_trimming(&mut si, 0.01, 1, 5, max_trim as usize, false, &mut trimmed);
                let got = format!(
                    "{}\t{}\t{trimmed}",
                    String::from_utf8_lossy(&si.seq),
                    if si.qual.is_empty() { "-".into() } else { String::from_utf8_lossy(&si.qual).to_string() },
                );
                assert_eq!(got, want[i], "doTrimming {pattern} maxTrim={max_trim} case {i}");
            }
        }
    }
}
