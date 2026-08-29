//! `coco filter`, ported from `src/filter.cpp`.
//!
//! Drops reads whose count profile has a drop narrower than one k-mer span -- the
//! signature of a chimeric junction, where each half of the read is well covered
//! but no k-mer bridges the join.

use crate::countprofile::CountProfile;
use crate::info::Info;
use crate::options::{CommandId, Options};
use crate::preprocessing::TableKind;
use crate::runner::{process_paired_reads, process_reads, PairedReadTask, ReadTask};
use crate::seq::{write_record, SeqInfoMode, SequenceInfo};
use crate::translator::KmerTranslator;

struct FilterTask {
    threshold: f64,
}

#[derive(Default, Clone, Copy)]
pub struct FilterStats {
    pub kept: u64,
    pub dropped: u64,
}

impl ReadTask for FilterTask {
    type Stats = FilterStats;

    fn process(
        &self,
        cp: &mut CountProfile,
        seq: &mut SequenceInfo,
        skip: bool,
        out: &mut Vec<u8>,
        stats: &mut FilterStats,
    ) {
        if skip {
            write_record(out, seq, SeqInfoMode::Auto);
            stats.kept += 1;
            return;
        }
        if cp.check_for_spurious_transition_drops_with_window_new(self.threshold) {
            stats.dropped += 1;
        } else {
            write_record(out, seq, SeqInfoMode::Auto);
            stats.kept += 1;
        }
    }

    fn merge(lhs: &mut FilterStats, rhs: &FilterStats) {
        lhs.kept += rhs.kept;
        lhs.dropped += rhs.dropped;
    }
}

impl PairedReadTask for FilterTask {
    type Stats = FilterStats;

    fn process(
        &self,
        cp1: &mut CountProfile,
        cp2: &mut CountProfile,
        r1: &mut SequenceInfo,
        r2: &mut SequenceInfo,
        skip: bool,
        out1: &mut Vec<u8>,
        out2: &mut Vec<u8>,
        stats: &mut FilterStats,
    ) {
        if skip {
            write_record(out1, r1, SeqInfoMode::Auto);
            write_record(out2, r2, SeqInfoMode::Auto);
            stats.kept += 1;
            return;
        }
        // The C++ short-circuits: mate 2 is only examined if mate 1 passed.
        let s1 = cp1.check_for_spurious_transition_drops_with_window_new(self.threshold);
        let s2 = if !s1 {
            cp2.check_for_spurious_transition_drops_with_window_new(self.threshold)
        } else {
            false
        };
        if !s1 && !s2 {
            write_record(out1, r1, SeqInfoMode::Auto);
            write_record(out2, r2, SeqInfoMode::Auto);
            stats.kept += 1;
        } else {
            stats.dropped += 1;
        }
    }

    fn merge(lhs: &mut FilterStats, rhs: &FilterStats) {
        lhs.kept += rhs.kept;
        lhs.dropped += rhs.dropped;
    }
}

pub fn run(opt: &Options, info: &Info) -> Result<i32, String> {
    let translator = KmerTranslator::new(&opt.spaced_kmer_pattern).map_err(|e| format!("{e}\n"))?;
    let ext = super::extension_for(opt, CommandId::Filter)?;
    let table = super::make_lookuptable(opt, &translator, info, TableKind::Compact)?;
    let task = FilterTask { threshold: opt.threshold };

    info.info("Step 2: Filter chimeric reads...\n");
    if !opt.reads.is_empty() {
        let path = super::out_path(opt, &opt.reads, &format!(".filter.reads{ext}"));
        let mut w = std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(&path).map_err(|e| e.to_string())?);
        info.info(&format!("...process reads from file {}\n", opt.reads));
        process_reads(&opt.reads, table.as_ref(), &translator, &task, opt.skip as i64,
                      opt.threads as usize, &mut w).map_err(|e| format!("{e}\n"))?;
        use std::io::Write;
        w.flush().map_err(|e| e.to_string())?;
        info.info("...completed\n");
    }
    if !opt.forward_reads.is_empty() && !opt.reverse_reads.is_empty() {
        let p1 = super::out_path(opt, &opt.forward_reads, &format!(".filter.1{ext}"));
        let p2 = super::out_path(opt, &opt.reverse_reads, &format!(".filter.2{ext}"));
        let mut w1 = std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(&p1).map_err(|e| e.to_string())?);
        let mut w2 = std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(&p2).map_err(|e| e.to_string())?);
        info.info(&format!("...process paired reads from files {} and {}\n", opt.forward_reads, opt.reverse_reads));
        process_paired_reads(&opt.forward_reads, &opt.reverse_reads, table.as_ref(), &translator,
                             &task, opt.skip as i64, opt.threads as usize, &mut w1, &mut w2)
            .map_err(|e| format!("{e}\n"))?;
        use std::io::Write;
        w1.flush().map_err(|e| e.to_string())?;
        w2.flush().map_err(|e| e.to_string())?;
        info.info("...completed\n");
    }
    Ok(0)
}
