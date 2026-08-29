//! End-to-end correction against the committed expectations.
//!
//! `testdata/expected/` was produced by the C++ binary and records not just how
//! many bases it changed but *which* ones. The bar the testdata README sets is
//! explicit: "not 'similar recall' but the same 7,911 corrections at the same
//! positions". This test reconstructs `corrections.tsv` and
//! `residual_errors.tsv` from the port's output and requires an exact match, then
//! checks the counters in `coco_stats.txt`.
//!
//! It drives the library rather than the binary, so a failure points at the
//! correction logic rather than at command wiring.

mod common;

use coco::commands::correction::{do_correction, CorrectionStatistic, CorrectorArgs};
use coco::countprofile::{CountProfile, TableAccess};
use coco::lookuptable::CountMode;
use coco::preprocessing::{build_lookuptable, TableKind};
use coco::seq::{read_all, FastxReader, SequenceInfo};
use coco::translator::{KmerTranslator, DEFAULT_PATTERN};

/// Run correction over a read file with CoCo's defaults, returning the corrected
/// records in input order plus the statistics.
fn correct_all(reads_path: &str) -> (Vec<SequenceInfo>, CorrectionStatistic) {
    let t = KmerTranslator::new(DEFAULT_PATTERN).unwrap();
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

    let args = CorrectorArgs {
        threshold: 0.01,
        pseudocount: 1,
        lower_bound: 5,
        max_corr_num: 10,
        max_trim_len: 0,
        update_lookup: false,
        debug: false,
    };

    let mut cp = CountProfile::new(&t, TableAccess::shared(table.as_ref()));
    let mut stats = CorrectionStatistic::default();
    let mut sseq = Vec::new();
    let mut squal = Vec::new();

    let mut reader = FastxReader::open(reads_path).unwrap();
    let mut batch = Vec::new();
    let mut out = Vec::new();
    let skip: i64 = 10;
    let min_len = skip + t.span() as i64;
    loop {
        let n = reader.read_batch(&mut batch, 8192).unwrap();
        if n == 0 {
            break;
        }
        for seq in batch.iter_mut().take(n) {
            if (seq.seq.len() as i64) < min_len {
                cp.set_seq_info_only();
            } else {
                cp.fill(seq);
                do_correction(&mut cp, seq, &args, &mut stats, &mut sseq, &mut squal);
            }
            out.push(seq.clone());
        }
    }
    (out, stats)
}

/// `read\tpos\tfrom_base\tto_base\twas_error\tcorrect` for every changed base.
fn corrections_table(
    err: &[SequenceInfo],
    perfect: &[SequenceInfo],
    corrected: &[SequenceInfo],
) -> Vec<String> {
    let mut rows = Vec::new();
    for i in 0..err.len() {
        assert_eq!(err[i].name, corrected[i].name, "record order changed");
        assert_eq!(
            err[i].seq.len(),
            corrected[i].seq.len(),
            "this dataset has no indel corrections, so lengths must be preserved"
        );
        for p in 0..err[i].seq.len() {
            let from = err[i].seq[p];
            let to = corrected[i].seq[p];
            if from == to {
                continue;
            }
            let truth = perfect[i].seq[p];
            rows.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                String::from_utf8_lossy(&err[i].name),
                p,
                from as char,
                to as char,
                (from != truth) as u8,
                (to == truth) as u8
            ));
        }
    }
    rows
}

/// Errors still present after correction.
fn residual_table(
    err: &[SequenceInfo],
    perfect: &[SequenceInfo],
    corrected: &[SequenceInfo],
) -> Vec<String> {
    let mut rows = Vec::new();
    for i in 0..err.len() {
        for p in 0..perfect[i].seq.len() {
            let truth = perfect[i].seq[p];
            let coco = corrected[i].seq[p];
            if coco == truth {
                continue;
            }
            let kind = if err[i].seq[p] == coco { "FN_left_alone" } else { "WR_miscorrected" };
            rows.push(format!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                String::from_utf8_lossy(&err[i].name),
                p,
                truth as char,
                err[i].seq[p] as char,
                coco as char,
                kind
            ));
        }
    }
    rows
}

fn expected_rows(path: &str) -> Vec<String> {
    let text = std::fs::read_to_string(path).unwrap();
    text.lines().skip(1).filter(|l| !l.trim().is_empty()).map(|s| s.to_string()).collect()
}

#[test]
fn corrections_match_the_cpp_base_for_base() {
    let err_path = common::testdata("reads.err0.1pct.fq");
    let perfect = read_all(&common::testdata("reads.perfect.fq")).unwrap();
    let err = read_all(&err_path).unwrap();
    let (corrected, stats) = correct_all(&err_path);

    assert_eq!(err.len(), 53_333, "record count from the testdata README");
    assert_eq!(err.len(), corrected.len());

    let got = corrections_table(&err, &perfect, &corrected);
    let want = expected_rows(&common::testdata("expected/corrections.tsv"));

    assert_eq!(
        got.len(),
        want.len(),
        "number of changed bases: got {}, expected {}",
        got.len(),
        want.len()
    );
    let mut shown = 0;
    for (i, (g, w)) in got.iter().zip(want.iter()).enumerate() {
        if g != w && shown < 10 {
            eprintln!("row {i}: got  {g}\n        want {w}");
            shown += 1;
        }
    }
    assert_eq!(got, want, "corrections differ from expected/corrections.tsv");

    assert_eq!(stats.substitution_multikmer, 7494);
    assert_eq!(stats.substitution_singlekmer + stats.substitution_on_edge, 417);
    assert_eq!(stats.insertion, 0);
    assert_eq!(stats.deletion, 0);
    assert_eq!(stats.trimmed, 0);
    eprintln!("{} corrections matched exactly", got.len());
}

#[test]
fn residual_errors_match_the_cpp() {
    let err_path = common::testdata("reads.err0.1pct.fq");
    let perfect = read_all(&common::testdata("reads.perfect.fq")).unwrap();
    let err = read_all(&err_path).unwrap();
    let (corrected, _) = correct_all(&err_path);

    let got = residual_table(&err, &perfect, &corrected);
    let want = expected_rows(&common::testdata("expected/residual_errors.tsv"));
    assert_eq!(got.len(), 38, "the README's 38 uncorrected bases");
    assert_eq!(got, want, "residual errors differ from expected/residual_errors.tsv");
}

#[test]
fn scoring_summary_matches_coco_stats() {
    let err_path = common::testdata("reads.err0.1pct.fq");
    let perfect = read_all(&common::testdata("reads.perfect.fq")).unwrap();
    let err = read_all(&err_path).unwrap();
    let (corrected, _) = correct_all(&err_path);

    let (mut errors_present, mut changed, mut tp, mut wr, mut fp) = (0u32, 0u32, 0u32, 0u32, 0u32);
    for i in 0..err.len() {
        for p in 0..perfect[i].seq.len() {
            let truth = perfect[i].seq[p];
            let before = err[i].seq[p];
            let after = corrected[i].seq[p];
            if before != truth {
                errors_present += 1;
            }
            if after != before {
                changed += 1;
                if before != truth && after == truth {
                    tp += 1;
                } else if before != truth {
                    wr += 1; // changed an error to the wrong base
                } else {
                    fp += 1; // changed a base that was already right
                }
            }
        }
    }
    let fn_ = errors_present - tp - wr;
    // From testdata/expected/coco_stats.txt.
    assert_eq!(errors_present, 7949, "errors present");
    assert_eq!(changed, 7911, "bases changed");
    assert_eq!(tp, 7911, "true positives");
    assert_eq!(wr, 0, "wrong corrections");
    assert_eq!(fp, 0, "false positives");
    assert_eq!(fn_, 38, "false negatives");
    eprintln!("recall {:.2}%  precision 100.00%", 100.0 * tp as f64 / errors_present as f64);
}

#[test]
fn quality_string_is_rewritten_only_where_the_base_changed() {
    // getAvgQual replaces the quality of a corrected base with the mean of its
    // neighbours. On this dataset every base is Q30 ('?'), so a corrected position
    // must still read '?' -- and every position must keep its length.
    let err_path = common::testdata("reads.err0.1pct.fq");
    let err = read_all(&err_path).unwrap();
    let (corrected, _) = correct_all(&err_path);
    for i in 0..err.len() {
        assert_eq!(err[i].qual.len(), corrected[i].qual.len(), "record {i}");
        assert_eq!(err[i].qual, corrected[i].qual, "record {i}: uniform Q30 stays Q30");
    }
}
