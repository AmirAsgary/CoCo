//! `coco abundance`, ported from `src/abundanceEstimator.cpp`.
//!
//! Reports the 67th percentile of a read's count profile as its abundance --
//! high enough to ignore error-depressed positions, low enough not to be pulled up
//! by repeats.

use crate::countprofile::CountProfile;
use crate::info::Info;
use crate::options::Options;
use crate::preprocessing::TableKind;
use crate::runner::{process_reads, ReadTask};
use crate::seq::SequenceInfo;
use crate::translator::KmerTranslator;

struct AbundanceTask;

impl ReadTask for AbundanceTask {
    type Stats = ();

    fn process(
        &self,
        cp: &mut CountProfile,
        seq: &mut SequenceInfo,
        skip: bool,
        out: &mut Vec<u8>,
        _stats: &mut (),
    ) {
        out.extend_from_slice(&seq.name);
        out.push(b'\t');
        if skip {
            out.push(b'-');
        } else {
            let est = cp.calc_x_quantile(0.67, &[]);
            out.extend_from_slice(est.to_string().as_bytes());
        }
        out.push(b'\n');
    }

    fn merge(_lhs: &mut (), _rhs: &()) {}
}

pub fn run(opt: &Options, info: &Info) -> Result<i32, String> {
    let translator = KmerTranslator::new(&opt.spaced_kmer_pattern).map_err(|e| format!("{e}\n"))?;
    let table = super::make_lookuptable(opt, &translator, info, TableKind::Compact)?;

    info.info("Step 2: Abundance estimation...\n");
    for (src, suffix) in [
        (&opt.reads, ".abundance.reads.tsv"),
        (&opt.forward_reads, ".abundance.1.tsv"),
        (&opt.reverse_reads, ".abundance.2.tsv"),
    ] {
        if src.is_empty() {
            continue;
        }
        let path = super::out_path(opt, src, suffix);
        let mut w = std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(&path).map_err(|e| e.to_string())?);
        info.info(&format!("...process reads from file {src}\n"));
        process_reads(src, table.as_ref(), &translator, &AbundanceTask, opt.skip as i64,
                      opt.threads as usize, &mut w).map_err(|e| format!("{e}\n"))?;
        use std::io::Write;
        w.flush().map_err(|e| e.to_string())?;
        info.info("...completed\n");
    }
    Ok(0)
}
