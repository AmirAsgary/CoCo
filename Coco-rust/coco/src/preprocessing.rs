//! Building the count table, ported from `src/preprocessing.cpp`.

use crate::dsk::DskCountsFile;
use crate::lookuptable::{
    CompactLookupTable, CountMode, GridLookupTable, HashCountTable, LookupTable, ShardedLookupTable,
};
use crate::seq::FastxReader;
use crate::translator::KmerTranslator;
use crate::types::res2int;

/// Which table layout `build_lookuptable` should produce.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TableKind {
    /// Open addressing. Same answers, no `2^30`-word grid.
    Compact,
    /// The C++ grid, needed when `counts2flat` output must be byte-identical.
    Grid,
}

/// `buildLookuptable`: read DSK counts, translate to spaced k-mers, index them.
#[allow(clippy::too_many_arguments)]
pub fn build_lookuptable(
    count_file: &str,
    count_mode: CountMode,
    translator: &KmerTranslator,
    min_count: u32,
    kind: TableKind,
    threads: usize,
    log: &mut dyn FnMut(&str),
) -> Result<Box<dyn LookupTable>, String> {
    let dsk = DskCountsFile::open(count_file).map_err(|e| e.0)?;
    let kmer_span = translator.span() as u32;
    if dsk.kmer_size != kmer_span {
        return Err(format!(
            "ERROR: kmerSize {} used in hdf5 file {} is not supported.\n\
             Please pre-compute kmer counts with k={}\n",
            dsk.kmer_size, count_file, kmer_span
        ));
    }
    let (log_index, log_offset) = translator.best_split();
    let nb_items = dsk.num_items().map_err(|e| e.0)?;

    match kind {
        TableKind::Compact => {
            // Which table to build is a measured choice, not a guess. Best of
            // three, 2.9 M reads / 14.7 M solid k-mers, on one Raven GPU node:
            //
            //                     1 thread   8 threads   36 threads
            //   single table        93.9 s      13.3 s        4.59 s
            //   sharded table       93.6 s      12.5 s        4.10 s
            //
            // and on `abundance`, which does less work per read so the build
            // weighs more: 16.5 s vs 17.9 s at one thread, 2.25 s vs 1.68 s at 36.
            // Sharding costs a little single-threaded -- an extra indirection per
            // lookup and a bucketing pass -- and pays from two threads up.
            if threads <= 1 {
                let mut t =
                    CompactLookupTable::with_capacity(nb_items, log_index, log_offset, count_mode);
                log("...fill lookuptable with spaced k-mer count pairs\n");
                dsk.for_each(|c| t.add_element(translator.kmer2min_packed(c.value), c.abundance))
                    .map_err(|e| e.0)?;
                log("...final setup\n");
                t.final_setup_tables(min_count);
                log("...completed\n");
                return Ok(Box::new(t));
            }
            let mut t = ShardedLookupTable::with_capacity(
                nb_items, log_index, log_offset, count_mode, threads,
            );
            log("...fill lookuptable with spaced k-mer count pairs\n");
            dsk.for_each_partition(|recs| t.add_batch_parallel(recs, translator))
                .map_err(|e| e.0)?;
            log("...final setup\n");
            t.final_setup_tables(min_count);
            log("...completed\n");
            Ok(Box::new(t))
        }
        TableKind::Grid => {
            let mut t = GridLookupTable::new(nb_items, log_index, log_offset, count_mode);
            log("...construct grids\n");
            dsk.for_each(|c| {
                t.assign_kmer_to_grid(translator.kmer2min_packed(c.value));
            })
            .map_err(|e| e.0)?;
            t.setup_index_grid_table();
            log("...fill lookuptable with spaced k-mer count pairs\n");
            dsk.for_each(|c| {
                t.add_element(translator.kmer2min_packed(c.value), c.abundance);
            })
            .map_err(|e| e.0)?;
            log("...final setup\n");
            t.final_setup_tables(min_count);
            log("...completed\n");
            Ok(Box::new(t))
        }
    }
}

/// `buildHashTable`: count spaced k-mers straight from the reads.
pub fn build_hash_table(
    read_filenames: &[String],
    translator: &KmerTranslator,
    log: &mut dyn FnMut(&str),
) -> Result<Box<dyn LookupTable>, String> {
    log("WARNING: using internal hash table to count k-mers is not recommended for larger datasets: use --counts instead\n");
    let kmer_span = translator.span() as usize;
    let mut table = HashCountTable::new();
    log("count spaced k-mers...\n");

    let mask: u128 = if kmer_span * 2 >= 128 {
        u128::MAX
    } else {
        (1u128 << (kmer_span * 2)) - 1
    };

    for path in read_filenames {
        log(&format!("...for input file {path}\n"));
        let mut reader = FastxReader::open(path).map_err(|e| format!("{path}: {e}"))?;
        let mut batch = Vec::new();
        loop {
            let n = reader.read_batch(&mut batch, 8192).map_err(|e| e.to_string())?;
            if n == 0 {
                break;
            }
            for seq in batch.iter().take(n) {
                let bytes = &seq.seq;
                if bytes.len() < kmer_span {
                    continue;
                }
                let mut spaced: u128 = 0;
                let mut l = 0usize;
                for &b in bytes.iter() {
                    let c = res2int(b);
                    if c != -1 {
                        spaced = ((spaced << 2) | c as u128) & mask;
                        l += 1;
                        if l >= kmer_span {
                            table.increase_count(translator.kmer2min_packed(spaced));
                        }
                    } else {
                        l = 0;
                        spaced = 0;
                    }
                }
            }
        }
    }
    log("...completed\n");
    Ok(Box::new(table))
}
