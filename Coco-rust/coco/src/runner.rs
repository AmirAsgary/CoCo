//! Read iteration, ported from `src/runner.cpp` -- and the one place the port
//! deliberately does more than the original.
//!
//! The C++ `--threads` option was never implemented; `processReads` is a plain
//! loop. Here reads are pulled in batches and the per-read work is spread over a
//! thread pool, with each chunk writing into its own output buffer that is then
//! appended in input order. Output is therefore byte-identical to a sequential
//! run, which is what makes the parallel path testable against the C++ at all.
//!
//! `--update-lookup` is the exception: it feeds corrections back into the shared
//! count table, so results depend on the order reads are processed. That path stays
//! sequential and the command layer says so.

use crate::countprofile::{CountProfile, TableAccess};
use crate::lookuptable::LookupTable;
use crate::seq::{FastxReader, SequenceInfo};
use crate::translator::KmerTranslator;
use rayon::prelude::*;
use std::io::Write;

/// What a command does with one read.
pub trait ReadTask: Sync {
    type Stats: Default + Send + Clone;

    /// Handle one record. `skip` means the read was too short to profile, in which
    /// case `cp` holds no valid profile.
    fn process(
        &self,
        cp: &mut CountProfile,
        seq: &mut SequenceInfo,
        skip: bool,
        out: &mut Vec<u8>,
        stats: &mut Self::Stats,
    );

    fn merge(lhs: &mut Self::Stats, rhs: &Self::Stats);
}

/// Same, for paired reads processed together.
pub trait PairedReadTask: Sync {
    type Stats: Default + Send + Clone;

    #[allow(clippy::too_many_arguments)]
    fn process(
        &self,
        cp1: &mut CountProfile,
        cp2: &mut CountProfile,
        r1: &mut SequenceInfo,
        r2: &mut SequenceInfo,
        skip: bool,
        out1: &mut Vec<u8>,
        out2: &mut Vec<u8>,
        stats: &mut Self::Stats,
    );

    fn merge(lhs: &mut Self::Stats, rhs: &Self::Stats);
}

const BATCH_RECORDS: usize = 16384;

/// Minimum read length that gets a count profile: `skip + kmerSpan`.
///
/// The C++ writes `seq.sequence.l < skip + kmerSpan`, where `skip` is `int` and
/// `kmerSpan` is `unsigned int`, so the addition happens in 32-bit *unsigned*
/// arithmetic. A negative `--skip` therefore wraps to a huge threshold and skips
/// every read rather than skipping none. Computed the same way here.
#[inline]
fn skip_threshold(skip: i64, kmer_span: u16) -> u64 {
    (skip as u32).wrapping_add(kmer_span as u32) as u64
}

/// Records per parallel chunk. Small enough that a slow read cannot leave one
/// worker holding the whole batch, large enough that per-chunk setup disappears.
fn chunk_size(threads: usize) -> usize {
    (BATCH_RECORDS / (threads * 4).max(1)).clamp(32, 2048)
}

/// `processReads`.
pub fn process_reads<T: ReadTask>(
    reads_name: &str,
    table: &dyn LookupTable,
    translator: &KmerTranslator,
    task: &T,
    skip: i64,
    threads: usize,
    writer: &mut dyn Write,
) -> std::io::Result<T::Stats> {
    let mut reader = FastxReader::open(reads_name)?;
    let min_len = skip_threshold(skip, translator.span());

    if threads <= 1 {
        let mut batch: Vec<SequenceInfo> = Vec::new();
        let mut stats = T::Stats::default();
        let mut cp = CountProfile::new(translator, TableAccess::shared(table));
        let mut out = Vec::with_capacity(1 << 20);
        loop {
            let n = reader.read_batch(&mut batch, BATCH_RECORDS)?;
            if n == 0 {
                break;
            }
            out.clear();
            for seq in batch.iter_mut().take(n) {
                if seq.seq.is_empty() {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "ERROR: Invalid sequence record found",
                    ));
                }
                let too_short = (seq.seq.len() as u64) < min_len;
                if too_short {
                    cp.set_seq_info_only();
                } else {
                    cp.fill(seq);
                }
                task.process(&mut cp, seq, too_short, &mut out, &mut stats);
            }
            writer.write_all(&out)?;
        }
        return Ok(stats);
    }

    // Parallel path: a three-stage pipeline rather than a read/compute/write
    // cycle per batch. Reading 1.8 GB of FASTQ and writing 900 MB back are each
    // about a second of wall clock; done between parallel phases they serialise
    // and cap the speedup well below the core count, so the reader and writer run
    // as their own threads and overlap with the workers.
    //
    // Batches carry their sequence number and the writer reorders them, so the
    // output file is identical to a sequential run.
    let threads = threads.max(1);
    let batch_records = pipeline_batch_records(threads);
    let inflight = 2 * threads;

    let (tx_in, rx_in) = std::sync::mpsc::sync_channel::<(usize, Vec<SequenceInfo>, usize)>(inflight);
    let (tx_out, rx_out) =
        std::sync::mpsc::sync_channel::<(usize, Vec<u8>, T::Stats)>(inflight);
    // Batches are handed back to the reader so the per-record Vecs are allocated
    // once rather than once per batch.
    let (tx_recycle, rx_recycle) = std::sync::mpsc::sync_channel::<Vec<SequenceInfo>>(inflight + 4);
    let rx_in = std::sync::Mutex::new(rx_in);
    let read_error: std::sync::Mutex<Option<std::io::Error>> = std::sync::Mutex::new(None);

    let read_error_ref = &read_error;
    let stats = std::thread::scope(|scope| -> std::io::Result<T::Stats> {
        // Reader. Takes ownership of the input reader and both of its channel
        // endpoints; a Receiver is Send but not Sync, so it cannot be borrowed
        // into a scoped thread.
        scope.spawn(move || {
            let mut seq_no = 0usize;
            loop {
                let mut batch = rx_recycle.try_recv().unwrap_or_default();
                let n = match reader.read_batch(&mut batch, batch_records) {
                    Ok(n) => n,
                    Err(e) => {
                        *read_error_ref.lock().unwrap() = Some(e);
                        break;
                    }
                };
                if n == 0 {
                    break;
                }
                if batch[..n].iter().any(|s| s.seq.is_empty()) {
                    *read_error_ref.lock().unwrap() = Some(std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "ERROR: Invalid sequence record found",
                    ));
                    break;
                }
                if tx_in.send((seq_no, batch, n)).is_err() {
                    break;
                }
                seq_no += 1;
            }
            drop(tx_in);
        });

        // Workers.
        for _ in 0..threads {
            let tx_out = tx_out.clone();
            let tx_recycle = tx_recycle.clone();
            let rx_in = &rx_in;
            scope.spawn(move || {
                let mut cp = CountProfile::new(translator, TableAccess::shared(table));
                loop {
                    let item = { rx_in.lock().unwrap().recv() };
                    let Ok((seq_no, mut batch, n)) = item else { break };
                    let mut out = Vec::with_capacity(n * 256);
                    let mut st = T::Stats::default();
                    for seq in batch.iter_mut().take(n) {
                        let too_short = (seq.seq.len() as u64) < min_len;
                        if too_short {
                            cp.set_seq_info_only();
                        } else {
                            cp.fill(seq);
                        }
                        task.process(&mut cp, seq, too_short, &mut out, &mut st);
                    }
                    let _ = tx_recycle.send(batch);
                    if tx_out.send((seq_no, out, st)).is_err() {
                        break;
                    }
                }
            });
        }
        drop(tx_out);

        // Writer: this thread, emitting batches in input order.
        let mut pending: std::collections::HashMap<usize, (Vec<u8>, T::Stats)> =
            std::collections::HashMap::new();
        let mut next = 0usize;
        let mut stats = T::Stats::default();
        for (seq_no, out, st) in rx_out {
            pending.insert(seq_no, (out, st));
            while let Some((out, st)) = pending.remove(&next) {
                writer.write_all(&out)?;
                T::merge(&mut stats, &st);
                next += 1;
            }
        }
        Ok(stats)
    })?;

    if let Some(e) = read_error.lock().unwrap().take() {
        return Err(e);
    }
    Ok(stats)
}

/// Records per pipeline batch.
///
/// Small enough that every worker gets many batches -- one slow read must not
/// leave a thread holding a large share of the file -- and large enough that the
/// channel handoff and the writer's reordering stay negligible.
fn pipeline_batch_records(threads: usize) -> usize {
    let _ = threads;
    4096
}

/// `processReads` with a table the task may modify. Always sequential.
pub fn process_reads_updating<T: ReadTask>(
    reads_name: &str,
    table: &mut dyn LookupTable,
    translator: &KmerTranslator,
    task: &T,
    skip: i64,
    writer: &mut dyn Write,
) -> std::io::Result<T::Stats> {
    let mut reader = FastxReader::open(reads_name)?;
    let mut batch: Vec<SequenceInfo> = Vec::new();
    let mut stats = T::Stats::default();
    let min_len = skip_threshold(skip, translator.span());

    // Safety: this function is the only user of the table for its duration and
    // never spawns a thread.
    let access = unsafe { TableAccess::exclusive(table) };
    let mut cp = CountProfile::new(translator, access);
    let mut out = Vec::with_capacity(1 << 20);
    loop {
        let n = reader.read_batch(&mut batch, BATCH_RECORDS)?;
        if n == 0 {
            break;
        }
        out.clear();
        for seq in batch.iter_mut().take(n) {
            if seq.seq.is_empty() {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "ERROR: Invalid sequence record found",
                ));
            }
            let too_short = (seq.seq.len() as u64) < min_len;
            if too_short {
                cp.set_seq_info_only();
            } else {
                cp.fill(seq);
            }
            task.process(&mut cp, seq, too_short, &mut out, &mut stats);
        }
        writer.write_all(&out)?;
    }
    Ok(stats)
}

/// `processPairedReads`.
#[allow(clippy::too_many_arguments)]
pub fn process_paired_reads<T: PairedReadTask>(
    forward: &str,
    reverse: &str,
    table: &dyn LookupTable,
    translator: &KmerTranslator,
    task: &T,
    skip: i64,
    threads: usize,
    writer1: &mut dyn Write,
    writer2: &mut dyn Write,
) -> std::io::Result<T::Stats> {
    let mut r1 = FastxReader::open(forward)?;
    let mut r2 = FastxReader::open(reverse)?;
    let mut b1: Vec<SequenceInfo> = Vec::new();
    let mut b2: Vec<SequenceInfo> = Vec::new();
    let mut stats = T::Stats::default();
    let min_len = skip_threshold(skip, translator.span());
    let cs = chunk_size(threads.max(1));

    loop {
        let n1 = r1.read_batch(&mut b1, BATCH_RECORDS)?;
        let n2 = r2.read_batch(&mut b2, BATCH_RECORDS)?;
        let n = n1.min(n2);
        if n == 0 {
            break;
        }
        if b1[..n].iter().any(|s| s.seq.is_empty()) || b2[..n].iter().any(|s| s.seq.is_empty()) {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "ERROR: Invalid sequence record found when processing paired reads",
            ));
        }

        let pairs: Vec<(&mut SequenceInfo, &mut SequenceInfo)> = b1[..n]
            .iter_mut()
            .zip(b2[..n].iter_mut())
            .collect();

        let results: Vec<(Vec<u8>, Vec<u8>, T::Stats)> = if threads <= 1 {
            let mut cp1 = CountProfile::new(translator, TableAccess::shared(table));
            let mut cp2 = CountProfile::new(translator, TableAccess::shared(table));
            let mut o1 = Vec::new();
            let mut o2 = Vec::new();
            let mut st = T::Stats::default();
            for (s1, s2) in pairs {
                let too_short =
                    (s1.seq.len() as u64) < min_len || (s2.seq.len() as u64) < min_len;
                if too_short {
                    cp1.set_seq_info_only();
                    cp2.set_seq_info_only();
                } else {
                    cp1.fill(s1);
                    cp2.fill(s2);
                }
                task.process(&mut cp1, &mut cp2, s1, s2, too_short, &mut o1, &mut o2, &mut st);
            }
            vec![(o1, o2, st)]
        } else {
            let mut chunks: Vec<Vec<(&mut SequenceInfo, &mut SequenceInfo)>> = Vec::new();
            let mut cur = Vec::with_capacity(cs);
            for p in pairs {
                cur.push(p);
                if cur.len() == cs {
                    chunks.push(std::mem::replace(&mut cur, Vec::with_capacity(cs)));
                }
            }
            if !cur.is_empty() {
                chunks.push(cur);
            }
            chunks
                .into_par_iter()
                .map(|chunk| {
                    let mut cp1 = CountProfile::new(translator, TableAccess::shared(table));
                    let mut cp2 = CountProfile::new(translator, TableAccess::shared(table));
                    let mut o1 = Vec::with_capacity(chunk.len() * 256);
                    let mut o2 = Vec::with_capacity(chunk.len() * 256);
                    let mut st = T::Stats::default();
                    for (s1, s2) in chunk {
                        let too_short =
                            (s1.seq.len() as u64) < min_len || (s2.seq.len() as u64) < min_len;
                        if too_short {
                            cp1.set_seq_info_only();
                            cp2.set_seq_info_only();
                        } else {
                            cp1.fill(s1);
                            cp2.fill(s2);
                        }
                        task.process(
                            &mut cp1, &mut cp2, s1, s2, too_short, &mut o1, &mut o2, &mut st,
                        );
                    }
                    (o1, o2, st)
                })
                .collect()
        };

        for (o1, o2, st) in &results {
            writer1.write_all(o1)?;
            writer2.write_all(o2)?;
            T::merge(&mut stats, st);
        }
    }
    Ok(stats)
}
