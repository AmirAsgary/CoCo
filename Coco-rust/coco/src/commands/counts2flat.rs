//! `coco counts2flat`, ported from `src/counts2flat.cpp`. Dumps the spaced k-mer
//! table as text.
//!
//! This is the one command that uses the grid table rather than the compact one:
//! its output order is the grid's own layout, so reproducing it byte for byte
//! means reproducing the grid. That costs the same `2^LOGINDEXSIZE`-word
//! allocation the C++ makes (8.6 GB for the default pattern).

use crate::info::Info;
use crate::options::Options;
use crate::preprocessing::TableKind;
use crate::translator::KmerTranslator;

pub fn run(opt: &Options, info: &Info) -> Result<i32, String> {
    let translator = KmerTranslator::new(&opt.spaced_kmer_pattern).map_err(|e| format!("{e}\n"))?;
    if opt.count_file.is_empty() {
        return Err("ERROR: Missing count file\n".into());
    }
    let table = super::make_lookuptable(opt, &translator, info, TableKind::Grid)?;
    let path = super::out_path(opt, &opt.count_file, ".counts2flat.tsv");
    let mut w = std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(&path).map_err(|e| e.to_string())?);
    table.iterate_over_all(&mut w).map_err(|e| e.to_string())?;
    use std::io::Write;
    w.flush().map_err(|e| e.to_string())?;
    Ok(0)
}
