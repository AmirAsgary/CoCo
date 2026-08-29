//! The five CoCo subcommands.

pub mod abundance;
pub mod correction;
pub mod counts2flat;
pub mod filter;
pub mod profile;

use crate::info::Info;
use crate::lookuptable::{CountMode, LookupTable};
use crate::options::{CommandId, OptId, Options};
use crate::preprocessing::{build_hash_table, build_lookuptable, TableKind};
use crate::translator::KmerTranslator;

/// The read files a command was pointed at, in the C++'s order.
pub fn read_filenames(opt: &Options) -> Vec<String> {
    let mut v = Vec::new();
    if opt.is_set(OptId::Reads) {
        v.push(opt.reads.clone());
    }
    if opt.is_set(OptId::ForwardReads) {
        v.push(opt.forward_reads.clone());
    }
    if opt.is_set(OptId::ReverseReads) {
        v.push(opt.reverse_reads.clone());
    }
    v
}

/// Step 1 of every command: get a count table, from `--counts` or from the reads.
pub fn make_lookuptable(
    opt: &Options,
    translator: &KmerTranslator,
    info: &Info,
    kind: TableKind,
) -> Result<Box<dyn LookupTable>, String> {
    info.info("Step 1: Generate lookuptable...\n");
    let mut log = |s: &str| info.info(s);
    if opt.is_set(OptId::CountFile) {
        build_lookuptable(
            &opt.count_file,
            CountMode::from_int(opt.count_mode),
            translator,
            0,
            kind,
            opt.threads.max(1) as usize,
            &mut log,
        )
    } else {
        build_hash_table(&read_filenames(opt), translator, &mut log)
    }
}

/// Output file path for a command, replicating the C++ naming.
pub fn out_path(opt: &Options, source: &str, suffix: &str) -> String {
    let prefix = if opt.is_set(OptId::Outprefix) {
        opt.outprefix.clone()
    } else {
        crate::options::get_filename(source)
    };
    format!("{}{}{}", opt.outdir, prefix, suffix)
}

/// `.fa` or `.fq`, decided from the input files as the C++ does.
pub fn extension_for(opt: &Options, cmd: CommandId) -> Result<&'static str, String> {
    let _ = cmd;
    let files = read_filenames(opt);
    let mut mode: Option<crate::seq::SeqInfoMode> = None;
    for f in &files {
        let m = crate::seq::get_seq_mode(f).map_err(|e| e.to_string())?;
        match mode {
            None => mode = Some(m),
            Some(prev) if prev != m => {
                return Err("ERROR: read files have inconsistent file formats\n".into())
            }
            _ => {}
        }
    }
    Ok(match mode {
        Some(crate::seq::SeqInfoMode::Fasta) => ".fa",
        _ => ".fq",
    })
}
