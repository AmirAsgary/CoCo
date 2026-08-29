//! Shared fixture locations for the golden tests.
//!
//! Paths are overridable so the suite can run against a different checkout, but
//! they are never silently skipped: a missing fixture fails with the command that
//! regenerates it, because a golden test that quietly does nothing is worse than
//! no test.

#![allow(dead_code)]

/// `Coco-rust/coco`, resolved at compile time.
pub fn crate_dir() -> &'static str {
    env!("CARGO_MANIFEST_DIR")
}

/// The CoCo checkout's `testdata/`, two levels up from the crate.
pub fn testdata_dir() -> String {
    std::env::var("COCO_TESTDATA").unwrap_or_else(|_| format!("{}/../../testdata", crate_dir()))
}

/// `Coco-rust/work`, where generated fixtures live.
pub fn work_dir() -> String {
    std::env::var("COCO_WORK").unwrap_or_else(|_| format!("{}/../work", crate_dir()))
}

pub fn testdata(name: &str) -> String {
    let p = format!("{}/{}", testdata_dir(), name);
    assert!(
        std::path::Path::new(&p).exists(),
        "missing fixture {p}\n(set COCO_TESTDATA to the checkout's testdata directory)"
    );
    p
}

/// The DSK counts file, which is not committed.
pub fn counts_file() -> String {
    let p = format!("{}/counts.err0.1pct.h5", work_dir());
    assert!(
        std::path::Path::new(&p).exists(),
        "missing {p}\nRegenerate with:\n  \
         dsk -file {}/reads.err0.1pct.fq -kmer-size 41 -abundance-min 2 \\\n    \
         -out counts.err0.1pct -out-dir {} -out-tmp $TMPDIR -max-memory 4000 -nb-cores 4",
        testdata_dir(),
        work_dir()
    );
    p
}

/// `dsk2ascii` output for the counts file: `<41-mer> <abundance>` per line.
pub fn solid_ascii() -> String {
    let p = format!("{}/solid_ascii.txt", work_dir());
    assert!(
        std::path::Path::new(&p).exists(),
        "missing {p}\nRegenerate with:\n  dsk2ascii -file {} -out {p}",
        counts_file()
    );
    p
}
