//! `coco correction`, ported from `src/correction.cpp`.
//!
//! Per read: maximize the profile, run substitution correction until it stops
//! finding anything, then indel/edge correction, then a second substitution pass
//! that no longer insists a base be pinned by more than one k-mer, then optional
//! trimming. If the total number of edits exceeds `--max-corr-num` the read is
//! reverted -- that many changes usually means a strain difference, not errors.

use crate::countprofile::{CountProfile, ALL_CORRECTED, ERROR_FREE, SOME_CORRECTED};
use crate::info::Info;
use crate::options::{CommandId, OptId, Options};
use crate::preprocessing::TableKind;
use crate::runner::{process_paired_reads, process_reads, process_reads_updating, PairedReadTask, ReadTask};
use crate::seq::{write_record, SeqInfoMode, SequenceInfo};
use crate::translator::KmerTranslator;

#[derive(Debug, Default, Clone, Copy)]
pub struct CorrectionStatistic {
    pub substitution_multikmer: u32,
    pub substitution_singlekmer: u32,
    pub substitution_on_edge: u32,
    pub insertion: u32,
    pub deletion: u32,
    pub trimmed: u32,
}

impl CorrectionStatistic {
    fn total(&self) -> u32 {
        self.substitution_multikmer
            + self.substitution_singlekmer
            + self.substitution_on_edge
            + self.insertion
            + self.deletion
            + self.trimmed
    }
    fn add(&mut self, o: &CorrectionStatistic) {
        self.substitution_multikmer += o.substitution_multikmer;
        self.substitution_singlekmer += o.substitution_singlekmer;
        self.substitution_on_edge += o.substitution_on_edge;
        self.insertion += o.insertion;
        self.deletion += o.deletion;
        self.trimmed += o.trimmed;
    }
}

pub struct CorrectorArgs {
    pub threshold: f64,
    pub pseudocount: u32,
    pub lower_bound: u32,
    pub max_corr_num: i32,
    pub max_trim_len: i32,
    pub update_lookup: bool,
    pub debug: bool,
}

/// `doCorrection`.
pub fn do_correction(
    cp: &mut CountProfile,
    seq: &mut SequenceInfo,
    args: &CorrectorArgs,
    total: &mut CorrectionStatistic,
    saved_seq: &mut Vec<u8>,
    saved_qual: &mut Vec<u8>,
) {
    // Kept so the read can be restored if too many corrections accumulate.
    saved_seq.clear();
    saved_seq.extend_from_slice(&seq.seq);
    saved_qual.clear();
    saved_qual.extend_from_slice(&seq.qual);

    let mut stat = CorrectionStatistic::default();
    let update_lookup = args.update_lookup;
    let mut status;

    cp.maximize();
    loop {
        status = cp.do_substitution_correction(
            seq,
            args.threshold,
            args.pseudocount,
            args.lower_bound,
            true,
            update_lookup,
            &mut stat.substitution_multikmer,
        );
        if status == SOME_CORRECTED || status == ALL_CORRECTED {
            cp.maximize();
        }
        if status != SOME_CORRECTED {
            break;
        }
    }

    let changed = cp.do_indel_correction(
        seq,
        args.threshold,
        args.pseudocount,
        args.lower_bound,
        true,
        update_lookup,
        &mut stat.substitution_on_edge,
        &mut stat.insertion,
        &mut stat.deletion,
    );
    if changed {
        cp.maximize();
    }

    if status != ERROR_FREE {
        loop {
            status = cp.do_substitution_correction(
                seq,
                args.threshold,
                args.pseudocount,
                args.lower_bound,
                false,
                update_lookup,
                &mut stat.substitution_singlekmer,
            );
            if status == SOME_CORRECTED || status == ALL_CORRECTED {
                cp.maximize();
            }
            if status != SOME_CORRECTED {
                break;
            }
        }
    }

    if args.max_trim_len > 0 {
        cp.do_trimming(
            seq,
            args.threshold,
            args.pseudocount,
            args.lower_bound,
            args.max_trim_len as usize,
            update_lookup,
            &mut stat.trimmed,
        );
    }

    if args.max_corr_num > 0 && stat.total() as i32 > args.max_corr_num {
        // Too many edits for one read: more likely a different strain than errors.
        seq.seq.clear();
        seq.seq.extend_from_slice(saved_seq);
        seq.qual.clear();
        seq.qual.extend_from_slice(saved_qual);
        if update_lookup {
            cp.update(seq, true);
        }
    } else {
        total.add(&stat);
    }
}

struct CorrectionTask {
    args: CorrectorArgs,
}

/// Per-chunk state: the revert buffers, kept out of `do_correction` so they are
/// allocated once per worker rather than once per read.
#[derive(Default, Clone)]
pub struct CorrStats {
    pub stat: CorrectionStatistic,
    scratch_seq: Vec<u8>,
    scratch_qual: Vec<u8>,
}

impl ReadTask for CorrectionTask {
    type Stats = CorrStats;

    fn process(
        &self,
        cp: &mut CountProfile,
        seq: &mut SequenceInfo,
        skip: bool,
        out: &mut Vec<u8>,
        stats: &mut CorrStats,
    ) {
        if skip {
            if self.args.debug {
                eprint!(
                    "WARNING: sequence {} is too short, it'll be skipped\n",
                    String::from_utf8_lossy(&seq.name)
                );
            }
        } else {
            cp.debug = self.args.debug;
            let (mut s, mut q) = (
                std::mem::take(&mut stats.scratch_seq),
                std::mem::take(&mut stats.scratch_qual),
            );
            do_correction(cp, seq, &self.args, &mut stats.stat, &mut s, &mut q);
            stats.scratch_seq = s;
            stats.scratch_qual = q;
            if self.args.debug && !cp.debug_buf.is_empty() {
                use std::io::Write;
                let _ = std::io::stderr().write_all(&cp.debug_buf);
                cp.debug_buf.clear();
            }
        }
        write_record(out, seq, SeqInfoMode::Auto);
    }

    fn merge(lhs: &mut CorrStats, rhs: &CorrStats) {
        lhs.stat.add(&rhs.stat);
    }
}

impl PairedReadTask for CorrectionTask {
    type Stats = CorrStats;

    fn process(
        &self,
        cp1: &mut CountProfile,
        cp2: &mut CountProfile,
        r1: &mut SequenceInfo,
        r2: &mut SequenceInfo,
        skip: bool,
        out1: &mut Vec<u8>,
        out2: &mut Vec<u8>,
        stats: &mut CorrStats,
    ) {
        if !skip {
            // Mates are corrected independently; the C++ does the same.
            cp1.debug = self.args.debug;
            cp2.debug = self.args.debug;
            let (mut s, mut q) = (
                std::mem::take(&mut stats.scratch_seq),
                std::mem::take(&mut stats.scratch_qual),
            );
            do_correction(cp1, r1, &self.args, &mut stats.stat, &mut s, &mut q);
            do_correction(cp2, r2, &self.args, &mut stats.stat, &mut s, &mut q);
            stats.scratch_seq = s;
            stats.scratch_qual = q;
            if self.args.debug {
                use std::io::Write;
                let mut e = std::io::stderr();
                let _ = e.write_all(&cp1.debug_buf);
                let _ = e.write_all(&cp2.debug_buf);
                cp1.debug_buf.clear();
                cp2.debug_buf.clear();
            }
        }
        write_record(out1, r1, SeqInfoMode::Auto);
        write_record(out2, r2, SeqInfoMode::Auto);
    }

    fn merge(lhs: &mut CorrStats, rhs: &CorrStats) {
        lhs.stat.add(&rhs.stat);
    }
}

pub fn run(opt: &Options, info: &Info) -> Result<i32, String> {
    let translator = KmerTranslator::new(&opt.spaced_kmer_pattern).map_err(|e| format!("{e}\n"))?;
    let ext = super::extension_for(opt, CommandId::Correction)?;
    let mut table = super::make_lookuptable(opt, &translator, info, TableKind::Compact)?;

    let task = CorrectionTask {
        args: CorrectorArgs {
            threshold: opt.threshold,
            pseudocount: opt.pseudocount as u32,
            lower_bound: opt.lower_bound as u32,
            max_corr_num: opt.max_corr_num,
            max_trim_len: opt.max_trim_len,
            update_lookup: opt.update_lookup,
            debug: info.level() >= Info::DEBUG,
        },
    };

    info.info("Step 2: Sequencing error correction...\n");
    let mut total = CorrectionStatistic::default();

    if !opt.reads.is_empty() {
        let path = super::out_path(opt, &opt.reads, &format!(".corr.reads{ext}"));
        let mut w = std::io::BufWriter::with_capacity(
            8 << 20,
            std::fs::File::create(&path).map_err(|e| format!("ERROR: opening failed for file {path}: {e}\n"))?,
        );
        info.info(&format!("...process reads from file {}\n", opt.reads));
        let stats = if opt.update_lookup {
            process_reads_updating(&opt.reads, table.as_mut(), &translator, &task, opt.skip as i64, &mut w)
        } else {
            process_reads(
                &opt.reads,
                table.as_ref(),
                &translator,
                &task,
                opt.skip as i64,
                opt.threads as usize,
                &mut w,
            )
        }
        .map_err(|e| format!("{e}\n"))?;
        total.add(&stats.stat);
        use std::io::Write;
        w.flush().map_err(|e| e.to_string())?;
        info.info("...completed\n");
    }

    if !opt.forward_reads.is_empty() && !opt.reverse_reads.is_empty() {
        let p1 = super::out_path(opt, &opt.forward_reads, &format!(".corr.1{ext}"));
        let p2 = super::out_path(opt, &opt.reverse_reads, &format!(".corr.2{ext}"));
        let mut w1 = std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(&p1).map_err(|e| e.to_string())?);
        let mut w2 = std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(&p2).map_err(|e| e.to_string())?);
        info.info(&format!(
            "...process paired reads from files {} and {}\n",
            opt.forward_reads, opt.reverse_reads
        ));
        let stats = process_paired_reads(
            &opt.forward_reads,
            &opt.reverse_reads,
            table.as_ref(),
            &translator,
            &task,
            opt.skip as i64,
            if opt.update_lookup { 1 } else { opt.threads as usize },
            &mut w1,
            &mut w2,
        )
        .map_err(|e| format!("{e}\n"))?;
        total.add(&stats.stat);
        use std::io::Write;
        w1.flush().map_err(|e| e.to_string())?;
        w2.flush().map_err(|e| e.to_string())?;
        info.info("...completed\n");
    }

    let _ = &opt.is_set(OptId::Reads);
    info.info("### COCO ERROR CORRECTION STATISTIC ###\n");
    info.info(&format!(
        "substitution corrections (multi kmer step): {}\n",
        total.substitution_multikmer
    ));
    info.info(&format!(
        "substitution corrections (single kmer step): {}\n",
        total.substitution_singlekmer + total.substitution_on_edge
    ));
    info.info(&format!("insertion corrections: {}\n", total.insertion));
    info.info(&format!("deletion corrections: {}\n", total.deletion));
    info.info(&format!("trimmed nucleotides: {}\n", total.trimmed));
    Ok(0)
}
