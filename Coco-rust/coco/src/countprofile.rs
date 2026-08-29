//! Per-read spaced k-mer count profiles and the operations built on them.
//!
//! Port of `src/CountProfile.cpp`. The profile is one count per k-mer *start*
//! position in a read; a sequencing error depresses every k-mer that covers it, so
//! errors appear as dips. `maximize` spreads each k-mer's count back over the
//! positions it covers, turning the dip into a per-*base* signal.
//!
//! The layout differs from the C++ in one respect: the original stores an array of
//! 13-byte packed `CountProfileEntry` structs, while this keeps three parallel
//! arrays. `maximize`, the hottest loop here, then reads counts and validity
//! contiguously instead of striding over the k-mer field it never looks at.

use crate::lookuptable::LookupTable;
use crate::seq::{get_avg_qual, SequenceInfo};
use crate::sliding::{SlidingMax, SlidingMin};
use crate::translator::{KmerTranslator, NOT_INFORMATIVE};
use crate::types::{res2int, ALPHABET_SIZE, INT2RES};

pub const SIGNIFICANT_LEVEL_DIFF: u32 = 10;
const WINDOW_SIZE: usize = 5;

/// Return codes of `doSubstitutionCorrection`.
pub const ERROR_FREE: i32 = 0;
pub const NONE_CORRECTED: i32 = 1;
pub const SOME_CORRECTED: i32 = 2;
pub const ALL_CORRECTED: i32 = 3;
pub const TOO_MANY_ERRORS: i32 = 4;

/// Read access to the count table, plus optional write access.
///
/// `--update-lookup` makes correction feed its results back into the table, which
/// the C++ does through a plain pointer. That makes the result depend on read
/// order and rules out parallelism, so it is modelled explicitly here: `writable`
/// is set only on the sequential path, and the command layer refuses to combine it
/// with more than one thread.
pub struct TableAccess<'a> {
    ptr: *mut (dyn LookupTable + 'a),
    writable: bool,
}

impl<'a> TableAccess<'a> {
    /// Read-only handle. Safe to hand to any number of threads.
    pub fn shared(table: &'a dyn LookupTable) -> Self {
        TableAccess {
            ptr: table as *const (dyn LookupTable + 'a) as *mut (dyn LookupTable + 'a),
            writable: false,
        }
    }

    /// Read-write handle.
    ///
    /// # Safety
    /// The caller must guarantee that no other handle to the same table is in use
    /// concurrently -- in practice, that the run is single-threaded.
    pub unsafe fn exclusive(table: &'a mut (dyn LookupTable + 'a)) -> Self {
        TableAccess { ptr: table as *mut (dyn LookupTable + 'a), writable: true }
    }

    #[inline(always)]
    fn get(&self) -> &dyn LookupTable {
        // Safety: the pointer comes from a live borrow with lifetime 'a.
        unsafe { &*self.ptr }
    }

    #[inline(always)]
    pub fn count(&self, kmer: u64) -> u32 {
        self.get().get_count(kmer)
    }

    #[inline]
    fn increase(&self, kmer: u64) {
        if self.writable {
            // Safety: `writable` is only set through `exclusive`, whose contract is
            // that this handle is the only one in play.
            unsafe { (*self.ptr).increase_count(kmer) };
        }
    }

    #[inline]
    fn decrease(&self, kmer: u64) {
        if self.writable {
            unsafe { (*self.ptr).decrease_count(kmer) };
        }
    }
}

/// A read plus its count profile. Reused across reads to avoid reallocation.
pub struct CountProfile<'a> {
    translator: &'a KmerTranslator,
    table: TableAccess<'a>,
    /// Contiguous runs of informative positions in the spaced pattern, as
    /// `(offset, length)`. `maximize` walks these instead of the 32 individual
    /// offsets, which turns its inner loop into a contiguous max over a slice.
    runs: Vec<(usize, usize)>,

    counts: Vec<u32>,
    valid: Vec<u8>,
    kmers: Vec<u64>,
    profile_length: usize,

    /// Scratch reused by `maximize` and the correction passes.
    scratch_max: Vec<u32>,
    scratch_tol: Vec<u32>,
    scratch_affected: Vec<u64>,
    scratch_seq: Vec<u8>,
    scratch_qual: Vec<u8>,
    window_max: SlidingMax,

    pub ambig_corr: u64,

    /// When set, correction appends the C++'s `Info(Info::CDEBUG)` trace lines to
    /// [`Self::debug_buf`] instead of dropping them. The caller drains the buffer.
    pub debug: bool,
    pub debug_buf: Vec<u8>,
}

impl<'a> CountProfile<'a> {
    pub fn new(translator: &'a KmerTranslator, table: TableAccess<'a>) -> Self {
        let mut runs = Vec::new();
        let inv = translator.inverse_mask_array();
        let mut i = 0usize;
        while i < inv.len() {
            if inv[i] != NOT_INFORMATIVE {
                let start = i;
                while i < inv.len() && inv[i] != NOT_INFORMATIVE {
                    i += 1;
                }
                runs.push((start, i - start));
            } else {
                i += 1;
            }
        }
        CountProfile {
            translator,
            table,
            runs,
            counts: Vec::new(),
            valid: Vec::new(),
            kmers: Vec::new(),
            profile_length: 0,
            scratch_max: Vec::new(),
            scratch_tol: Vec::new(),
            scratch_affected: Vec::new(),
            scratch_seq: Vec::new(),
            scratch_qual: Vec::new(),
            window_max: SlidingMax::with_capacity(80),
            ambig_corr: 0,
            debug: false,
            debug_buf: Vec::new(),
        }
    }

    #[inline]
    pub fn profile_len(&self) -> usize {
        self.profile_length
    }

    #[inline]
    pub fn counts(&self) -> &[u32] {
        &self.counts[..self.profile_length]
    }

    #[inline]
    pub fn translator(&self) -> &KmerTranslator {
        self.translator
    }

    /// Mark the profile as not computed, matching `setSeqInfo`.
    pub fn set_seq_info_only(&mut self) {
        self.profile_length = 0;
    }

    /// Compute the profile for `seq`. `fill` in the C++.
    pub fn fill(&mut self, seq: &SequenceInfo) {
        self.update(seq, false);
    }

    /// Recompute the profile from the (possibly corrected) sequence.
    ///
    /// `update_lookuptable` first removes this read's old k-mers from the table and
    /// then adds the new ones, which is what makes `--update-lookup` order
    /// dependent.
    pub fn update(&mut self, seq: &SequenceInfo, update_lookuptable: bool) {
        if update_lookuptable {
            for i in 0..self.profile_length {
                if self.valid[i] != 0 {
                    self.table.decrease(self.kmers[i]);
                }
            }
        }

        let kmer_span = self.translator.span() as usize;
        let seq_len = seq.seq.len();
        // The C++ computes `len - span + 1` in unsigned arithmetic and relies on
        // the caller's `skip` check to keep it positive.
        self.profile_length = seq_len.saturating_sub(kmer_span - 1);
        let n = self.profile_length;
        if self.counts.len() < n {
            self.counts.resize(n, 0);
            self.valid.resize(n, 0);
            self.kmers.resize(n, 0);
        }
        if n == 0 {
            return;
        }

        let spaced_mask = self.translator.spaced_mask();
        let mut kmer: u128 = 0;
        let mut n_store: u128 = 0;
        let bytes = &seq.seq;

        for idx in 0..seq_len {
            let code = res2int(bytes[idx]);
            if code != -1 {
                kmer = (kmer << 2) | code as u128;
                n_store <<= 1;
            } else {
                kmer <<= 2;
                n_store = (n_store << 1) | 1;
            }

            if idx + 1 >= kmer_span {
                let p = idx + 1 - kmer_span;
                if (n_store & spaced_mask) != 0 {
                    self.counts[p] = 0;
                    self.valid[p] = 0;
                    self.kmers[p] = 0;
                    continue;
                }
                let packed = self.translator.kmer2min_packed(kmer);
                self.counts[p] = self.table.count(packed);
                self.valid[p] = 1;
                self.kmers[p] = packed;
                if update_lookuptable {
                    self.table.increase(packed);
                }
            }
        }
    }

    /// Spread each k-mer's count over the bases it covers. `maximize` in the C++.
    ///
    /// Returns a slice of `profile_length + span - 1` values, borrowed from
    /// internal scratch so repeated calls do not allocate.
    pub fn maximize(&mut self) -> &mut [u32] {
        let kmer_span = self.translator.span() as usize;
        let max_len = self.profile_length + kmer_span - 1;
        if self.scratch_max.len() < max_len {
            self.scratch_max.resize(max_len, 1);
        }
        let max_profile = &mut self.scratch_max[..max_len];
        max_profile.fill(1);

        for idx in 0..self.profile_length {
            if self.valid[idx] == 0 {
                continue;
            }
            let c = self.counts[idx];
            for &(off, len) in &self.runs {
                let s = idx + off;
                // Branchless: this lowers to cmov, or to vpmaxud for the longer
                // runs. The conditional form mispredicts on roughly every other
                // element, since counts along a read are noisy.
                for v in &mut max_profile[s..s + len] {
                    *v = (*v).max(c);
                }
            }
        }
        max_profile
    }

    /// `addCountPerPosition`.
    pub fn add_count_per_position(&self, summed: &mut Vec<u64>) {
        if summed.len() < self.profile_length {
            summed.resize(self.profile_length, 0);
        }
        for idx in 0..self.profile_length {
            summed[idx] += self.counts[idx] as u64;
        }
    }

    /// `calcXquantile`.
    ///
    /// Reproduces a bug in the original: when `positions_of_interest` is given, the
    /// C++ counts how many of those positions are in range but then copies the
    /// *first* that many entries of the list without re-filtering. Positions past
    /// the end are therefore still read whenever an earlier one was out of range.
    pub fn calc_x_quantile(&self, quantile: f64, positions_of_interest: &[u32]) -> u32 {
        let mut counts: Vec<u32>;
        if positions_of_interest.is_empty() {
            counts = self.counts[..self.profile_length].to_vec();
        } else {
            let total = positions_of_interest
                .iter()
                .filter(|&&p| (p as usize) < self.profile_length)
                .count();
            counts = Vec::with_capacity(total);
            for idx in 0..total {
                let p = positions_of_interest[idx] as usize;
                // Guard the out-of-range read the C++ performs unchecked.
                counts.push(if p < self.profile_length { self.counts[p] } else { 0 });
            }
        }
        if counts.is_empty() {
            return 0;
        }
        counts.sort_unstable();
        let idx = (quantile * counts.len() as f64) as usize;
        counts[idx.min(counts.len() - 1)]
    }

    pub fn calc_median(&self) -> u32 {
        self.calc_x_quantile(0.5, &[])
    }

    // -----------------------------------------------------------------------
    // Correction
    // -----------------------------------------------------------------------

    /// `calcNeighborhoodTolerance` over this profile's maximized counts.
    fn calc_neighborhood_tolerance(
        &mut self,
        max_profile_len: usize,
        threshold: f64,
        pseudocount: u32,
        lowerbound: u32,
    ) {
        let kmer_span = self.translator.span() as usize;
        calc_neighborhood_tolerance_into(
            &self.scratch_max[..max_profile_len.min(self.scratch_max.len())],
            max_profile_len,
            kmer_span,
            threshold,
            pseudocount,
            lowerbound,
            &mut self.scratch_tol,
            &mut self.window_max,
        );
    }

    /// `doSubstitutionCorrection`.
    ///
    /// Finds every base whose maximized count sits at or below the local tolerance,
    /// works out which k-mers each such base uniquely explains, and asks
    /// [`Self::first_last_unique_kmer_strategy`] which substitution repairs them.
    #[allow(clippy::too_many_arguments)]
    pub fn do_substitution_correction(
        &mut self,
        seq: &mut SequenceInfo,
        threshold: f64,
        pseudocount: u32,
        lowerbound: u32,
        need_multiple_kmers: bool,
        update_lookup: bool,
        corrected_substitutions: &mut u32,
    ) -> i32 {
        let kmer_span = self.translator.span() as usize;
        let kmer_weight = self.translator.weight() as usize;
        let max_profile_length = self.profile_length + kmer_span - 1;
        self.calc_neighborhood_tolerance(max_profile_length, threshold, pseudocount, lowerbound);

        let mut found_errors = 0usize;
        let mut error_positions = [0u32; 64];
        if self.scratch_affected.len() < self.profile_length {
            self.scratch_affected.resize(self.profile_length, 0);
        }
        self.scratch_affected[..self.profile_length].fill(0);

        let mask_array = self.translator.mask_array();
        for idx in 0..max_profile_length {
            if found_errors == 63 {
                return TOO_MANY_ERRORS;
            }
            if self.scratch_max[idx] <= self.scratch_tol[idx] {
                error_positions[found_errors] = idx as u32;
                for jdx in 0..kmer_weight {
                    let pos = idx as isize - mask_array[jdx] as isize;
                    if pos >= 0 && (pos as usize) < self.profile_length {
                        self.scratch_affected[pos as usize] |= 1u64 << found_errors;
                    }
                }
                found_errors += 1;
            }
        }
        if found_errors == 0 {
            return ERROR_FREE;
        }

        let mut corrected_errors = 0usize;
        for idx in 0..found_errors {
            let err_pos = error_positions[idx] as usize;
            let start = if err_pos >= kmer_span { err_pos - kmer_span + 1 } else { 0 };
            let end = if err_pos < self.profile_length { err_pos } else { self.profile_length - 1 };
            let mut first_unique = u32::MAX;
            let mut last_unique = 0u32;
            for jdx in start..=end {
                let a = self.scratch_affected[jdx];
                if a > 0 && a == (1u64 << idx) {
                    first_unique = first_unique.min(jdx as u32);
                    last_unique = last_unique.max(jdx as u32);
                }
            }

            if need_multiple_kmers && first_unique == last_unique {
                continue; // only one k-mer pins this base; not enough evidence
            }

            let target =
                self.first_last_unique_kmer_strategy(seq, err_pos, first_unique, last_unique);
            if target >= 0 {
                if self.debug {
                    // "corrected\t<name>\t<pos>\t<from>\t<to>", with <from> read
                    // before the substitution, as in the C++.
                    self.debug_buf.extend_from_slice(b"corrected\t");
                    self.debug_buf.extend_from_slice(&seq.name);
                    self.debug_buf.push(b'\t');
                    self.debug_buf.extend_from_slice(err_pos.to_string().as_bytes());
                    self.debug_buf.push(b'\t');
                    self.debug_buf.push(seq.seq[err_pos]);
                    self.debug_buf.push(b'\t');
                    self.debug_buf.push(INT2RES[target as usize]);
                    self.debug_buf.push(b'\n');
                }
                seq.seq[err_pos] = INT2RES[target as usize];
                corrected_errors += 1;
                if !seq.qual.is_empty() {
                    seq.qual[err_pos] = get_avg_qual(&seq.qual, err_pos);
                }
            }
        }

        *corrected_substitutions += corrected_errors as u32;

        if corrected_errors == 0 {
            NONE_CORRECTED
        } else if corrected_errors < found_errors {
            self.update(seq, update_lookup);
            SOME_CORRECTED
        } else {
            self.update(seq, update_lookup);
            ALL_CORRECTED
        }
    }

    /// `firstLastUniqueKmerCorrectionStrategy`.
    ///
    /// Returns the base code to substitute, `-1` if no candidate works, or `-2` if
    /// more than one does (ambiguous, so nothing is changed).
    fn first_last_unique_kmer_strategy(
        &mut self,
        seq: &SequenceInfo,
        error_pos: usize,
        first_unique_start: u32,
        last_unique_start: u32,
    ) -> i32 {
        if first_unique_start > last_unique_start {
            return -1;
        }
        let first_start = first_unique_start as usize;
        let last_start = last_unique_start as usize;
        let weight = self.translator.weight() as usize;
        let mask_array = self.translator.mask_array();
        let inv = self.translator.inverse_mask_array();
        let current_res = seq.seq[error_pos];

        let mut first_kmer: u64 = 0;
        let mut last_kmer: u64 = 0;
        for jdx in 0..weight {
            let m = mask_array[jdx] as usize;
            // `3 & res2int[..]` turns an invalid base into 3 ('G'), as in the C++.
            first_kmer = (first_kmer << 2) | (3 & res2int(seq.seq[first_start + m]) as u64);
            last_kmer = (last_kmer << 2) | (3 & res2int(seq.seq[last_start + m]) as u64);
        }

        let first_shift = mutation_shift(weight, inv[error_pos - first_start]);
        let last_shift = mutation_shift(weight, inv[error_pos - last_start]);
        let cur_code = res2int(current_res);

        let mut mutation_target: i32 = -1;
        for res_mutation in 0..ALPHABET_SIZE as u64 {
            if res_mutation as i8 == cur_code {
                continue;
            }
            first_kmer = (first_kmer & !(3u64 << first_shift)) | (res_mutation << first_shift);
            last_kmer = (last_kmer & !(3u64 << last_shift)) | (res_mutation << last_shift);

            let mut improvement = 0;
            if self.table.count(self.translator.packed2min_packed(first_kmer))
                > self.scratch_tol[error_pos]
            {
                improvement += 1;
            }
            if self.table.count(self.translator.packed2min_packed(last_kmer))
                > self.scratch_tol[error_pos]
            {
                improvement += 1;
            }

            if improvement == 2 {
                if mutation_target == -1 {
                    mutation_target = res_mutation as i32;
                } else {
                    mutation_target = -2; // ambiguous: leave the base alone
                    if first_unique_start != last_unique_start {
                        self.ambig_corr += 1;
                    }
                    break;
                }
            }
        }
        mutation_target
    }

    /// `edgeSubstitutionCorrection`: substitution near a read end, where the base
    /// is not uniquely pinned but is still covered by informative positions.
    fn edge_substitution_correction(&mut self, seq: &SequenceInfo, substitution_start: usize) -> i32 {
        let kmer_span = self.translator.span() as usize;
        let start = if substitution_start >= kmer_span {
            substitution_start - kmer_span + 1
        } else {
            0
        };
        let end = if substitution_start < self.profile_length {
            substitution_start
        } else {
            self.profile_length - 1
        };
        let inv = self.translator.inverse_mask_array();
        let mut first_unique = u32::MAX;
        let mut last_unique = 0u32;
        for jdx in start..=end {
            if inv[substitution_start - jdx] != NOT_INFORMATIVE {
                first_unique = first_unique.min(jdx as u32);
                last_unique = last_unique.max(jdx as u32);
            }
        }
        self.first_last_unique_kmer_strategy(seq, substitution_start, first_unique, last_unique)
    }

    /// `tryInsertionCorrection`: would deleting `insertion_len` bases at
    /// `insertion_start` make the three probe k-mers abundant again?
    fn try_insertion_correction(
        &self,
        seq: &SequenceInfo,
        insertion_start: usize,
        insertion_len: usize,
    ) -> bool {
        let kmer_span = self.translator.span() as usize;
        let first_kmer_start = if insertion_start >= kmer_span {
            insertion_start - kmer_span + 1
        } else {
            0
        };
        let second_last_kmer_start = if insertion_start < self.profile_length {
            insertion_start.wrapping_sub(1)
        } else {
            self.profile_length.wrapping_sub(insertion_len).wrapping_sub(1)
        };
        let mid_kmer_start = (first_kmer_start.wrapping_add(second_last_kmer_start)) / 2;

        if !(first_kmer_start < insertion_start
            && second_last_kmer_start.wrapping_add(kmer_span) >= insertion_start + insertion_len)
        {
            return false;
        }

        let weight = self.translator.weight() as usize;
        let mask_array = self.translator.mask_array();
        let mut first_kmer: u64 = 0;
        let mut mid_kmer: u64 = 0;
        let mut second_last_kmer: u64 = 0;
        for jdx in 0..weight {
            let m = mask_array[jdx] as usize;
            for (start, acc) in [
                (first_kmer_start, &mut first_kmer),
                (second_last_kmer_start, &mut second_last_kmer),
                (mid_kmer_start, &mut mid_kmer),
            ] {
                let mut p = start + m;
                if p >= insertion_start {
                    p += insertion_len;
                }
                if p >= seq.seq.len() {
                    return false;
                }
                let c = res2int(seq.seq[p]);
                if c == -1 {
                    return false;
                }
                *acc = (*acc << 2) | (3 & c as u64);
            }
        }

        let mut approved = 0;
        if self.table.count(self.translator.packed2min_packed(first_kmer))
            > self.tol_at(first_kmer_start)
        {
            approved += 1;
        }
        if self.table.count(self.translator.packed2min_packed(mid_kmer)) > self.tol_at(mid_kmer_start)
        {
            approved += 1;
        }
        if self.table.count(self.translator.packed2min_packed(second_last_kmer))
            > self.tol_at(second_last_kmer_start)
        {
            approved += 1;
        }
        approved == 3
    }

    #[inline]
    fn tol_at(&self, idx: usize) -> u32 {
        self.scratch_tol.get(idx).copied().unwrap_or(0)
    }

    /// `tryDeletionCorrection`: which base, inserted at `deletion_pos`, restores the
    /// two probe k-mers? `-1` for none, `-2` for ambiguous.
    fn try_deletion_correction(&self, seq: &SequenceInfo, deletion_pos: usize) -> i32 {
        let kmer_span = self.translator.span() as usize;
        let first_kmer_start = if deletion_pos >= kmer_span { deletion_pos - kmer_span + 1 } else { 0 };
        let last_kmer_start = if deletion_pos < self.profile_length {
            deletion_pos.wrapping_sub(1)
        } else {
            self.profile_length
        };

        let weight = self.translator.weight() as usize;
        let mask_array = self.translator.mask_array();
        let inv = self.translator.inverse_mask_array();
        let mut first_kmer: u64 = 0;
        let mut last_kmer: u64 = 0;
        for jdx in 0..weight {
            let m = mask_array[jdx] as usize;
            for (start, acc) in [(first_kmer_start, &mut first_kmer), (last_kmer_start, &mut last_kmer)] {
                let p = start.wrapping_add(m);
                if p < deletion_pos {
                    if p >= seq.seq.len() {
                        return -1;
                    }
                    let c = res2int(seq.seq[p]);
                    if c == -1 {
                        return -1;
                    }
                    *acc = (*acc << 2) | (3 & c as u64);
                } else if p == deletion_pos {
                    *acc <<= 2;
                } else {
                    if p == 0 || p - 1 >= seq.seq.len() {
                        return -1;
                    }
                    let c = res2int(seq.seq[p - 1]);
                    if c == -1 {
                        return -1;
                    }
                    *acc = (*acc << 2) | (3 & c as u64);
                }
            }
        }

        // Both indices are always within the span by construction, but the entry
        // they land on may be UCHAR_MAX; see `mutation_shift`.
        let fi = deletion_pos.wrapping_sub(first_kmer_start) % inv.len();
        let li = deletion_pos.wrapping_sub(last_kmer_start) % inv.len();
        let first_shift = mutation_shift(weight, inv[fi]);
        let last_shift = mutation_shift(weight, inv[li]);

        let mut mutation_target: i32 = -1;
        for res_mutation in 0..ALPHABET_SIZE as u64 {
            first_kmer = (first_kmer & !(3u64 << first_shift)) | (res_mutation << first_shift);
            last_kmer = (last_kmer & !(3u64 << last_shift)) | (res_mutation << last_shift);

            let mut improvement = 0;
            if self.table.count(self.translator.packed2min_packed(first_kmer)) > self.tol_at(deletion_pos)
            {
                improvement += 1;
            }
            if self.table.count(self.translator.packed2min_packed(last_kmer)) > self.tol_at(deletion_pos)
            {
                improvement += 1;
            }
            if improvement == 2 {
                if mutation_target == -1 {
                    mutation_target = res_mutation as i32;
                } else {
                    mutation_target = -2;
                    break;
                }
            }
        }
        mutation_target
    }

    /// `doIndelCorrection`: single insertions, single deletions, and substitutions
    /// too close to a read end for the substitution pass to pin down.
    #[allow(clippy::too_many_arguments)]
    pub fn do_indel_correction(
        &mut self,
        seq: &mut SequenceInfo,
        threshold: f64,
        pseudocount: u32,
        lowerbound: u32,
        try_substitution: bool,
        update_lookup: bool,
        corrected_substitutions: &mut u32,
        corrected_insertions: &mut u32,
        corrected_deletions: &mut u32,
    ) -> bool {
        let kmer_span = self.translator.span() as usize;
        let max_profile_length = self.profile_length + kmer_span - 1;
        self.calc_neighborhood_tolerance(max_profile_length, threshold, pseudocount, lowerbound);

        let mut changed = false;
        let mut drop_len: usize = 0;
        let mut offset: isize = 0;
        // The evaluation reads the *original* sequence while edits accumulate in
        // these copies, exactly as the C++ does.
        self.scratch_seq.clear();
        self.scratch_seq.extend_from_slice(&seq.seq);
        self.scratch_qual.clear();
        self.scratch_qual.extend_from_slice(&seq.qual);
        let mut sequence = std::mem::take(&mut self.scratch_seq);
        let mut qual = std::mem::take(&mut self.scratch_qual);

        // Mirrors the C++ loop shape exactly, including the `continue` on a
        // successful edge substitution that leaves `drop_len` *unreset*. That is
        // not an oversight to tidy up: it makes the same drop be reconsidered at
        // the next position for an insertion or deletion, and the counters differ
        // if the drop is cleared instead.
        for idx in 0..=self.profile_length {
            if idx < self.profile_length && self.counts[idx] <= self.scratch_tol[idx] {
                drop_len += 1;
            } else if drop_len > 0 {
                if drop_len <= kmer_span {
                    let mut insertion_approved;
                    let mut deletion_approved = false;
                    let mut deletion_pos = usize::MAX;
                    let mut insertion_pos = usize::MAX;
                    let mut substitution_pos = usize::MAX;

                    if idx == self.profile_length {
                        // Drop runs to the read end without recovering.
                        let p = idx - drop_len + kmer_span - 1;
                        deletion_pos = p;
                        insertion_pos = p;
                        substitution_pos = p;
                    } else if idx == drop_len {
                        // Drop starts at the read start.
                        deletion_pos = idx;
                        insertion_pos = idx - 1;
                        substitution_pos = idx - 1;
                    } else if drop_len >= kmer_span - 1 {
                        if drop_len == kmer_span - 1 {
                            deletion_pos = idx;
                        }
                        insertion_pos = idx - 1;
                    }

                    if try_substitution && substitution_pos != usize::MAX {
                        let res = self.edge_substitution_correction(seq, substitution_pos);
                        if res >= 0 {
                            let p = (substitution_pos as isize + offset) as usize;
                            if p < sequence.len() {
                                sequence[p] = INT2RES[res as usize];
                                if !qual.is_empty() && p < qual.len() {
                                    qual[p] = get_avg_qual(&qual, p);
                                }
                            }
                            changed = true;
                            *corrected_substitutions += 1;
                            continue; // deliberately does not reset drop_len
                        }
                    }

                    insertion_approved = false;
                    let mut res_to_add: i32 = -1;
                    if insertion_pos != usize::MAX && drop_len > 1 {
                        insertion_approved = self.try_insertion_correction(seq, insertion_pos, 1);
                    }
                    if deletion_pos != usize::MAX && drop_len != kmer_span {
                        res_to_add = self.try_deletion_correction(seq, deletion_pos);
                        if res_to_add >= 0 {
                            deletion_approved = true;
                        } else if res_to_add == -2 {
                            // Several bases would fit the deletion, which also
                            // invalidates the insertion hypothesis.
                            insertion_approved = false;
                        }
                    }

                    if insertion_approved && !deletion_approved {
                        let p = (insertion_pos as isize + offset) as usize;
                        if p < sequence.len() {
                            sequence.remove(p);
                            if !qual.is_empty() && p < qual.len() {
                                qual.remove(p);
                            }
                        }
                        changed = true;
                        offset -= 1;
                        *corrected_insertions += 1;
                    } else if deletion_approved && !insertion_approved && res_to_add >= 0 {
                        let p = (deletion_pos as isize + offset) as usize;
                        if p <= sequence.len() {
                            sequence.insert(p, INT2RES[res_to_add as usize]);
                            if !qual.is_empty() && p <= qual.len() {
                                qual.insert(p, 33);
                                qual[p] = get_avg_qual(&qual, p);
                            }
                        }
                        changed = true;
                        offset += 1;
                        *corrected_deletions += 1;
                    }
                    // Both approved means an ambiguous choice, so nothing changes.
                }
                drop_len = 0;
            }
        }

        seq.seq.clear();
        seq.seq.extend_from_slice(&sequence);
        seq.qual.clear();
        seq.qual.extend_from_slice(&qual);
        sequence.clear();
        qual.clear();
        self.scratch_seq = sequence;
        self.scratch_qual = qual;

        if changed {
            self.update(seq, update_lookup);
        }
        changed
    }

    /// `doTrimming`: cut leading/trailing runs that stayed below tolerance.
    #[allow(clippy::too_many_arguments)]
    pub fn do_trimming(
        &mut self,
        seq: &mut SequenceInfo,
        threshold: f64,
        pseudocount: u32,
        lowerbound: u32,
        max_trim_len: usize,
        update_lookup: bool,
        trimmed_counter: &mut u32,
    ) -> bool {
        let kmer_span = self.translator.span() as usize;
        let max_profile_length = self.profile_length + kmer_span - 1;
        self.calc_neighborhood_tolerance(max_profile_length, threshold, pseudocount, lowerbound);

        let mut changed = false;
        let mut drop_len = 0usize;
        let mut idx = 0usize;
        while idx < self.profile_length && self.counts[idx] <= self.scratch_tol[idx] {
            drop_len += 1;
            idx += 1;
        }
        if drop_len > 0 && drop_len <= max_trim_len && drop_len <= seq.seq.len() {
            seq.seq.drain(..drop_len);
            if !seq.qual.is_empty() {
                seq.qual.drain(..drop_len.min(seq.qual.len()));
            }
            changed = true;
            *trimmed_counter += drop_len as u32;
        }

        drop_len = 0;
        let mut idx = self.profile_length;
        while idx > 0 && self.counts[idx - 1] <= self.scratch_tol[idx - 1] {
            drop_len += 1;
            idx -= 1;
        }
        if drop_len > 0 && drop_len <= max_trim_len && drop_len <= seq.seq.len() {
            let n = seq.seq.len();
            seq.seq.truncate(n - drop_len);
            if !seq.qual.is_empty() {
                let q = seq.qual.len();
                seq.qual.truncate(q.saturating_sub(drop_len));
            }
            changed = true;
            *trimmed_counter += drop_len as u32;
        }

        if changed {
            self.update(seq, update_lookup);
        }
        changed
    }

    // -----------------------------------------------------------------------
    // Filtering
    // -----------------------------------------------------------------------

    /// `checkForSpuriousTransitionDropsWithWindowNew`: the chimera detector the
    /// `filter` command uses.
    ///
    /// A chimeric junction shows up as a count drop that is narrower than one
    /// k-mer span -- the two halves of the read are each abundant, but no k-mer
    /// spans the join.
    ///
    /// The C++ indexes `profile[idx + windowSize]` without checking, which reads
    /// past the end for reads shorter than `2 * windowSize + span - 1`. Those reads
    /// are reported as unfiltered here rather than reading out of bounds.
    pub fn check_for_spurious_transition_drops_with_window_new(&self, threshold: f64) -> bool {
        let longest_block = self.translator.longest_block() as usize;
        let window_size = (longest_block + 1).max(WINDOW_SIZE);
        let span = self.translator.span() as usize;

        if self.profile_length < 2 * window_size {
            return false;
        }

        let mut preceding = SlidingMinMax::new();
        let mut successive = SlidingMinMax::new();
        for idx in 0..window_size {
            preceding.push(self.counts[idx + window_size]);
            successive.push(self.counts[idx]);
        }

        let mut dropstart: isize = 0;
        for idx in (window_size - 1)..(self.profile_length - window_size) {
            successive.pop_front();
            successive.push(self.counts[idx]);
            preceding.pop_front();
            preceding.push(self.counts[idx + window_size]);

            let p_max = preceding.max();
            let p_min = preceding.min();
            let s_max = successive.max();
            let s_min = successive.min();

            if p_max == 0 || p_min == 0 || s_max == 0 || s_min == 0 {
                continue;
            }
            if (p_max as f64) / (s_min as f64) < threshold {
                dropstart = idx as isize + 1;
            }
            if (s_max as f64) / (p_min as f64) < threshold {
                let dropend = idx as isize;
                if dropstart >= 0 && dropend - dropstart + 1 < span as isize {
                    return true;
                }
                dropstart = -1;
            }
        }
        dropstart >= 0 && (self.profile_length as isize - dropstart) < span as isize
    }

    /// `checkForSpuriousTransitionDropsWithWindow`.
    ///
    /// **No command in the released C++ calls this.** It is the previous
    /// generation of the chimera filter, superseded by
    /// [`Self::check_for_spurious_transition_drops_with_window_new`], and is
    /// ported so the translation is complete. It is pinned by the oracle tests.
    ///
    /// Two passes: find positions whose count drops sharply below the surrounding
    /// window (candidates), then clear any candidate whose drop *survives*
    /// maximization -- a drop that survives is a genuinely low-coverage region,
    /// while one that disappears means no k-mer spans the position, which is what a
    /// chimeric junction looks like. Anything still flagged is a transition drop.
    ///
    /// The C++ subtracts `u32`s without checking order in two of the tests, so a
    /// larger subtrahend wraps to a huge value and the comparison succeeds.
    /// `wrapping_sub` keeps that.
    #[allow(unused_assignments)] // dropend mirrors the C++'s bookkeeping; some
                                 // stores are provably dead only given the
                                 // dropstart invariant, which the compiler cannot see
    pub fn check_for_spurious_transition_drops_with_window(
        &self,
        max_profile: &[u32],
        cov_est: u32,
        local_perc_drop: f64,
        global_perc_drop: f64,
        mask_only_drop_edges: bool,
    ) -> bool {
        let kmer_span = self.translator.span() as usize;
        let kmer_weight = self.translator.weight() as usize;
        let n = self.profile_length;
        if n == 0 {
            return false;
        }
        let mask_array = self.translator.mask_array();
        let half = kmer_span / 2;
        let corr_factor = 0.001f64;

        // Per-position allowance for seeing the same sequencing error repeatedly.
        let mut corr_values = vec![0u32; n];
        let mut corr_window = SlidingMax::new();
        for idx in 0..=half {
            corr_window.push(self.count_at(idx));
        }
        corr_values[0] = (corr_factor * corr_window.max() as f64 + 1.0) as u32;
        for idx in 1..n {
            if idx + half < n {
                corr_window.push(self.counts[idx + half]);
            }
            if idx > half {
                corr_window.pop_front();
            }
            corr_values[idx] = (corr_factor * corr_window.max() as f64 + 1.0) as u32;
        }

        let mut candidates = vec![0u32; n];
        let max_profile_length = n + kmer_span - 1;

        let mut window_back = SlidingMin::new();
        let mut window_front = SlidingMin::new();
        window_back.push(self.counts[0]);
        for idx in 1..WINDOW_SIZE + 1 {
            window_front.push(self.count_at(idx));
        }

        let mut dropstart_level = u32::MAX;
        let mut dropend_level;
        let mut dropstart = n;
        let mut dropend = n;

        if (self.counts[0] as f64) < cov_est as f64 * global_perc_drop {
            dropstart = 0;
        }

        for idx in 1..n {
            let back = window_back.min();
            let front = window_front.min();

            if (self.counts[idx] as f64) < local_perc_drop * (back as f64 - corr_values[idx - 1] as f64)
                && back.wrapping_sub(self.counts[idx]) >= SIGNIFICANT_LEVEL_DIFF
                && (self.counts[idx] as f64) < cov_est as f64 * global_perc_drop
            {
                dropstart = idx;
                dropend = n;
                dropstart_level = back;
            } else if window_front.len() == WINDOW_SIZE
                && ((local_perc_drop * (front as f64 - corr_values[idx] as f64) > back as f64
                    && front.wrapping_sub(back) >= SIGNIFICANT_LEVEL_DIFF)
                    || front as f64 > cov_est as f64 * global_perc_drop)
            {
                // A drop end only counts if a drop start was seen.
                if dropstart != n {
                    dropend = idx;
                    dropend_level = front;
                    let compare_level = dropstart_level.min(dropend_level);
                    for d in dropstart..dropend {
                        if (self.counts[d] as f64) < local_perc_drop * compare_level as f64 {
                            candidates[d] = compare_level;
                        }
                    }
                    for d in dropstart..dropend {
                        if (max_profile[d] as f64) < local_perc_drop * compare_level as f64
                            && (!mask_only_drop_edges
                                || ((d > 0
                                    && (max_profile[d] as f64)
                                        < local_perc_drop * max_profile[d - 1] as f64)
                                    || (max_profile[d] as f64)
                                        < local_perc_drop * max_profile[d + 1] as f64))
                        {
                            for jdx in 0..kmer_weight {
                                let pos = d as isize - mask_array[jdx] as isize;
                                if pos >= dropstart as isize && pos < n as isize {
                                    candidates[pos as usize] = 0;
                                }
                            }
                        }
                    }
                }
                dropstart = n;
                dropstart_level = u32::MAX;
                dropend = n;
            }

            if window_back.len() >= WINDOW_SIZE {
                window_back.pop_front();
            }
            window_back.push(self.counts[idx]);
            window_front.pop_front();
            if idx + WINDOW_SIZE < n {
                window_front.push(self.counts[idx + WINDOW_SIZE]);
            }
        }

        // A drop that runs to the end of the read never sees a drop end.
        if dropstart != 0 && dropstart != n {
            for d in dropstart..dropend {
                if (self.counts[d] as f64) < local_perc_drop * dropstart_level as f64 {
                    candidates[d] = dropstart_level;
                }
            }
            for d in dropstart..max_profile_length {
                if (max_profile[d] as f64) < local_perc_drop * dropstart_level as f64
                    && (!mask_only_drop_edges
                        || ((d > 0
                            && (max_profile[d] as f64) < local_perc_drop * max_profile[d - 1] as f64)
                            || (d + 1 < max_profile_length
                                && (max_profile[d] as f64)
                                    < local_perc_drop * max_profile[d + 1] as f64)))
                {
                    for jdx in 0..kmer_weight {
                        let pos = d as isize - mask_array[jdx] as isize;
                        if pos >= dropstart as isize && pos < n as isize {
                            candidates[pos as usize] = 0;
                        }
                    }
                }
            }
        }

        candidates.iter().any(|&c| c != 0)
    }

    /// `checkForSpuriousTransitionDrops`, marked "outdated" in the C++ header.
    ///
    /// **No command calls this either.** Ported for completeness and pinned by the
    /// oracle tests.
    ///
    /// Note the ratio tests are written as `a / b <= 0.5` on two `u32`s, so the
    /// division is *integer* division and the test really means `a < b`. That is
    /// reproduced. It also means the C++ divides by zero -- and takes SIGFPE --
    /// whenever the previous position has a count of 0. Since there is no result
    /// to reproduce there, this returns "not a drop" for a zero divisor.
    #[allow(unused_assignments)] // `dropend` mirrors the C++'s bookkeeping; the
                                 // stores the compiler flags are dead only given
                                 // the dropstart invariant, which it cannot see
    pub fn check_for_spurious_transition_drops(
        &self,
        max_profile: &[u32],
        drop_level_criterion: u32,
        mask_only_drop_edges: bool,
    ) -> bool {
        let n = self.profile_length;
        if n == 0 {
            return false;
        }
        let mut candidates = vec![0u8; n];
        let correction_factor = 0.001f64;
        let mut dropstart = 0usize;
        let mut dropend = n;

        // 1) find drops
        for idx in 1..n {
            let (cur, prev) = (self.counts[idx], self.counts[idx - 1]);
            let falls = prev != 0 && cur / prev == 0;
            let rises = cur != 0 && prev / cur == 0;
            if falls && cur <= drop_level_criterion && prev > drop_level_criterion {
                dropstart = idx;
            } else if rises && prev <= drop_level_criterion && cur > drop_level_criterion {
                if dropstart == n {
                    continue;
                }
                dropend = idx;
                for d in dropstart..dropend {
                    if self.counts[d] <= drop_level_criterion {
                        candidates[d] = 1;
                    }
                }
                dropstart = n;
                dropend = n;
            }
        }
        if dropstart > 0 && dropstart < n {
            for d in dropstart..n {
                if self.counts[d] <= drop_level_criterion {
                    candidates[d] = 1;
                }
            }
        }

        // 2) clear candidates whose drop disappears on the maximized profile
        let kmer_span = self.translator.span() as usize;
        let kmer_weight = self.translator.weight() as usize;
        let mask_array = self.translator.mask_array();
        let max_profile_len = n + kmer_span - 1;
        let mut dropstart = 0usize;
        let mut dropend = max_profile_len;

        let clear = |candidates: &mut [u8], from: usize| {
            for jdx in 0..kmer_weight {
                let pos = from as isize - mask_array[jdx] as isize;
                if pos >= 0 && (pos as usize) < n {
                    candidates[pos as usize] = 0;
                }
            }
        };

        for idx in 1..max_profile_len {
            let (cur, prev) = (max_profile[idx], max_profile[idx - 1]);
            // maxProfile is initialised to 1, so neither divisor can be zero here.
            let lo_prev = drop_level_criterion as f64 + correction_factor * prev as f64 + 1.0;
            let lo_cur = drop_level_criterion as f64 + correction_factor * cur as f64 + 1.0;
            if cur / prev == 0 && (cur as f64) <= lo_prev && (prev as f64) > lo_prev {
                dropstart = idx;
            } else if prev / cur == 0 && (prev as f64) <= lo_cur && (cur as f64) > lo_cur {
                if dropstart == max_profile_len {
                    continue;
                }
                dropend = idx;
                if mask_only_drop_edges {
                    clear(&mut candidates, dropstart);
                    if dropend > dropstart + 1 {
                        clear(&mut candidates, dropend - 1);
                    }
                } else {
                    for d in dropstart..dropend {
                        clear(&mut candidates, d);
                    }
                }
                dropstart = max_profile_len;
                dropend = max_profile_len;
            }
        }

        if dropstart > 0 && dropstart < max_profile_len {
            if mask_only_drop_edges {
                clear(&mut candidates, dropstart);
                if dropend > dropstart + 1 {
                    clear(&mut candidates, dropend - 1);
                }
            } else {
                for d in dropstart..max_profile_len {
                    clear(&mut candidates, d);
                }
            }
        }

        candidates.iter().any(|&c| c != 0)
    }

    /// Count at `idx`, or 0 past the end. The two legacy filters index a little
    /// past their own profile in the C++.
    #[inline]
    fn count_at(&self, idx: usize) -> u32 {
        self.counts.get(idx).copied().unwrap_or(0)
    }
}

/// A window that tracks both extremes, so the filter keeps one pass instead of
/// four `min_element`/`max_element` scans per position.
struct SlidingMinMax {
    lo: SlidingMin,
    hi: SlidingMax,
}

impl SlidingMinMax {
    fn new() -> Self {
        SlidingMinMax { lo: SlidingMin::new(), hi: SlidingMax::new() }
    }
    #[inline]
    fn push(&mut self, v: u32) {
        self.lo.push(v);
        self.hi.push(v);
    }
    #[inline]
    fn pop_front(&mut self) {
        self.lo.pop_front();
        self.hi.pop_front();
    }
    #[inline]
    fn min(&self) -> u32 {
        self.lo.min()
    }
    #[inline]
    fn max(&self) -> u32 {
        self.hi.max()
    }
}

/// `calcNeighborhoodTolerance` from `CountProfile.cpp`.
///
/// For each position this is `round(threshold * max(neighbourhood)) + pseudocount`,
/// where the neighbourhood is a roughly `span`-wide window over the maximized
/// profile -- unless the local maximum is below `lowerbound`, in which case
/// coverage is too thin to judge anything and the tolerance is 0.
///
/// The window is not a clean sliding one. It is seeded with `maxProfile[0..=span]`
/// and only starts advancing at `idx > span/2`, and the first value it appends is
/// one already in the window, so for a stretch of positions a duplicate sits
/// inside it. That is reproduced exactly; only the `max` is computed differently,
/// by monotonic deque instead of a linear scan.
///
/// The C++ is a file-static function with no external linkage, so this is exposed
/// for the oracle tests to reach it.
#[allow(clippy::too_many_arguments)]
pub fn calc_neighborhood_tolerance_into(
    max_profile: &[u32],
    max_profile_len: usize,
    kmer_span: usize,
    threshold: f64,
    pseudocount: u32,
    lowerbound: u32,
    out: &mut Vec<u32>,
    window: &mut SlidingMax,
) {
    if out.len() < max_profile_len {
        out.resize(max_profile_len, 0);
    }
    let at = |i: usize| max_profile.get(i).copied().unwrap_or(0);
    let half = kmer_span / 2;
    window.clear();
    for idx in 0..=kmer_span {
        window.push(at(idx));
    }
    for idx in 0..max_profile_len {
        if idx > half && idx + half < max_profile_len {
            window.push(at(idx + half));
            window.pop_front();
        }
        let max = window.max();
        out[idx] = if max >= lowerbound {
            ((threshold * max as f64).round() as i64 as u32).wrapping_add(pseudocount)
        } else {
            0
        };
    }
}

/// Bit offset used when substituting base `inv_index` of a packed k-mer.
///
/// `_inverse_mask_array` holds `UCHAR_MAX` for span positions the pattern does not
/// cover, and `tryDeletionCorrection` indexes it without checking that the
/// position is informative. The C++ then computes `2 * (weight - 255 - 1)` in
/// `int`, which is negative, converts it to `uint64_t` and passes it as a shift
/// count. x86-64 shift instructions take the count modulo 64, so instead of being
/// undefined in any visible way it lands on a different -- usually the lowest --
/// base, and the resulting k-mer counts change which corrections get approved.
///
/// Masking to 6 bits here reproduces that. Guarding it instead measurably changes
/// the insertion and deletion counters, because a `tryDeletionCorrection` that
/// returns "ambiguous" also cancels the competing insertion hypothesis.
#[inline]
fn mutation_shift(weight: usize, inv_index: u8) -> u32 {
    let s = 2i64 * (weight as i64 - inv_index as i64 - 1);
    (s as u64 & 63) as u32
}

// ---------------------------------------------------------------------------
// Hooks used by the oracle tests, which live in a separate crate and so cannot
// reach private state. Not used by the binary.
// ---------------------------------------------------------------------------

impl<'a> CountProfile<'a> {
    /// Install a profile directly, bypassing `fill`. `-1` marks an invalid
    /// position. Mirrors `setProfile` in the oracle harness.
    pub fn set_profile_for_test(&mut self, vals: &[i64]) {
        let n = vals.len();
        self.counts.clear();
        self.valid.clear();
        self.kmers.clear();
        self.counts.resize(n, 0);
        self.valid.resize(n, 0);
        self.kmers.resize(n, 0);
        for (i, &v) in vals.iter().enumerate() {
            self.valid[i] = if v < 0 { 0 } else { 1 };
            self.counts[i] = if v < 0 { 0 } else { v as u32 };
        }
        self.profile_length = n;
    }
}

/// Standalone `calcNeighborhoodTolerance`, for the oracle tests.
pub fn neighborhood_tolerance_for_test(
    translator: &KmerTranslator,
    max_profile: &[u32],
    threshold: f64,
    pseudocount: u32,
    lowerbound: u32,
) -> Vec<u32> {
    let mut out = Vec::new();
    let mut window = SlidingMax::new();
    calc_neighborhood_tolerance_into(
        max_profile,
        max_profile.len(),
        translator.span() as usize,
        threshold,
        pseudocount,
        lowerbound,
        &mut out,
        &mut window,
    );
    out.truncate(max_profile.len());
    out
}
