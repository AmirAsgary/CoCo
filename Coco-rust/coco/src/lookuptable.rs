//! Spaced k-mer count tables.
//!
//! Three implementations, all behind [`LookupTable`]:
//!
//! * [`GridLookupTable`] -- a faithful port of the templated C++ `Lookuptable`. It
//!   is the semantic reference and the oracle the fast table is tested against.
//! * [`CompactLookupTable`] -- the default. Same observable behaviour, but an open
//!   addressing hash table instead of a `2^30`-entry grid, which removes an 8.6 GB
//!   allocation and one dependent load per lookup.
//! * [`HashCountTable`] -- a port of C++ `HashTable`, used when no precomputed
//!   count file is supplied. Note its deliberately odd "absent means 1" contract.

use crate::kmer::{packed_kmer_to_string, PackedKmer};
use std::io::Write;

pub trait LookupTable: Send + Sync {
    fn get_count(&self, kmer: PackedKmer) -> u32;
    /// `+1` on an existing entry. Returns whether the k-mer was found.
    fn increase_count(&mut self, kmer: PackedKmer) -> bool;
    /// `-1` on an existing entry with a nonzero count.
    fn decrease_count(&mut self, kmer: PackedKmer) -> bool;
    /// Dump `<kmer>\t<count>` for every entry, in the C++ table's own order.
    fn iterate_over_all(&self, out: &mut dyn Write) -> std::io::Result<()>;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

// ---------------------------------------------------------------------------
// Grid table: faithful port of src/Lookuptable.cpp
// ---------------------------------------------------------------------------

/// Direct port of the C++ `Lookuptable<LOGINDEXSIZE, LOGOFFSETSIZE>`.
///
/// The C++ makes the two sizes template parameters and instantiates 21 variants;
/// here they are plain fields, which costs a couple of register reads and removes
/// the combinatorial explosion.
///
/// The C++ `IndexEntry` is a 9-byte packed bitfield (`indexOffset:LOGOFFSETSIZE`
/// plus a `u32` count). This uses two parallel arrays instead: 12 bytes per entry
/// rather than 9, but aligned loads and no bitfield extraction. The offset always
/// fits in `LOGOFFSETSIZE` bits by construction, so nothing is lost.
pub struct GridLookupTable {
    log_index_size: u32,
    log_offset_size: u32,
    index_grid_table: Vec<usize>,
    offsets: Vec<u64>,
    counts: Vec<u32>,
    /// Insertion order, needed to reproduce `iterateOverAll` ordering elsewhere.
    seq_no: Vec<u32>,
    number_items: usize,
    max_number_items: usize,
    offset_mask: u64,
    index_mask: u64,
    mode: CountMode,
    /// Count of `addElement` calls whose dedup scan reached back past the true
    /// start of the grid cell *and matched there*. See [`Self::cross_cell_merges`].
    cross_cell_merges: usize,
    /// Real start of each grid cell that has been written to, recorded on first
    /// insertion. Only populated when `track_cross_cell` is set, since a dense
    /// array would cost another `2^LOGINDEXSIZE` words.
    cell_starts: std::collections::HashMap<u32, usize>,
    track_cross_cell: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountMode {
    /// `Options::COUNT_MODE_SUM`
    Sum,
    /// `Options::COUNT_MODE_MAX`
    Max,
}

impl CountMode {
    pub fn from_int(v: i32) -> CountMode {
        if v == 1 {
            CountMode::Max
        } else {
            CountMode::Sum
        }
    }
}

impl GridLookupTable {
    pub fn new(nb_items: usize, log_index_size: u32, log_offset_size: u32, mode: CountMode) -> Self {
        let grid_size = (1usize << log_index_size) + 1;
        GridLookupTable {
            log_index_size,
            log_offset_size,
            index_grid_table: vec![0usize; grid_size],
            offsets: vec![0u64; nb_items],
            counts: vec![0u32; nb_items],
            seq_no: vec![0u32; nb_items],
            number_items: 0,
            max_number_items: nb_items,
            offset_mask: (1u64 << log_offset_size) - 1,
            index_mask: ((1u64 << log_index_size) - 1) << log_offset_size,
            mode,
            cross_cell_merges: 0,
            cell_starts: std::collections::HashMap::new(),
            track_cross_cell: false,
        }
    }

    /// Record cell starts so [`Self::cross_cell_merges`] is exact. Costs a hash
    /// map entry per occupied cell, so it is opt-in and used by tests.
    pub fn track_cross_cell_merges(&mut self, on: bool) {
        self.track_cross_cell = on;
    }

    #[inline(always)]
    fn grid_position(&self, kmer: PackedKmer) -> usize {
        ((kmer & self.index_mask) >> self.log_offset_size) as usize
    }

    #[inline(always)]
    fn offset(&self, kmer: PackedKmer) -> u64 {
        kmer & self.offset_mask
    }

    /// Pass 1: tally how many k-mers land in each grid cell.
    pub fn assign_kmer_to_grid(&mut self, kmer: PackedKmer) {
        let g = self.grid_position(kmer);
        self.index_grid_table[g] += 1;
    }

    /// Turn per-cell counts into cell start offsets (an exclusive prefix sum).
    pub fn setup_index_grid_table(&mut self) {
        let n = self.index_grid_table.len();
        let mut prev = self.index_grid_table[0];
        self.index_grid_table[0] = 0;
        for i in 0..n - 1 {
            let cur = self.index_grid_table[i + 1];
            self.index_grid_table[i + 1] = self.index_grid_table[i] + prev;
            prev = cur;
        }
    }

    /// Pass 2: insert, merging duplicate offsets within the cell.
    pub fn add_element(&mut self, kmer: PackedKmer, count: u32) -> usize {
        let grid_position = self.grid_position(kmer);
        let writing_position = self.index_grid_table[grid_position];
        let offset = self.offset(kmer);

        // The C++ takes the *previous* cell's running write pointer as the lower
        // bound of the dedup scan. Once that cell is full the pointer equals this
        // cell's start and the bound is correct, but while it is still filling the
        // scan reaches into the previous cell's live entries and can merge a k-mer
        // into the wrong cell. `cross_cell_merges` counts exactly that event so a
        // test can assert it never fires on real input.
        let prev_writing_position = if grid_position > 0 {
            self.index_grid_table[grid_position - 1]
        } else {
            0
        };

        assert!(
            writing_position < self.max_number_items,
            "ERROR: Lookuptable addElement overflows. Current writing position is {writing_position}"
        );

        let true_cell_start = if self.track_cross_cell {
            *self.cell_starts.entry(grid_position as u32).or_insert(writing_position)
        } else {
            prev_writing_position
        };
        for pos in prev_writing_position..writing_position {
            if self.counts[pos] != 0 && self.offsets[pos] == offset {
                if pos < true_cell_start {
                    self.cross_cell_merges += 1;
                }
                match self.mode {
                    CountMode::Sum => self.counts[pos] = self.counts[pos].wrapping_add(count),
                    CountMode::Max => self.counts[pos] = self.counts[pos].max(count),
                }
                return pos;
            }
        }

        self.offsets[writing_position] = offset;
        self.counts[writing_position] = count;
        self.seq_no[writing_position] = self.number_items as u32;
        self.index_grid_table[grid_position] += 1;
        self.number_items += 1;
        writing_position
    }

    /// Number of merges that crossed a grid-cell boundary during construction.
    pub fn cross_cell_merges(&self) -> usize {
        self.cross_cell_merges
    }

    /// Compact away entries at or below `count_threshold` and restore the grid
    /// table to cell-start form.
    pub fn final_setup_tables(&mut self, count_threshold: u32) {
        let mut prev = 0usize;
        let mut readpos = 0usize;
        let mut writepos = 0usize;
        for idx in 0..self.index_grid_table.len() {
            while readpos < self.index_grid_table[idx] {
                if self.counts[readpos] > count_threshold {
                    if readpos != writepos {
                        self.offsets[writepos] = self.offsets[readpos];
                        self.counts[writepos] = self.counts[readpos];
                        self.seq_no[writepos] = self.seq_no[readpos];
                    }
                    writepos += 1;
                }
                readpos += 1;
            }
            self.index_grid_table[idx] = prev;
            prev = writepos;
        }
        debug_assert!(writepos <= self.number_items);
        self.number_items = writepos;
        self.max_number_items = writepos;
        self.offsets.truncate(writepos);
        self.counts.truncate(writepos);
        self.seq_no.truncate(writepos);
        self.offsets.shrink_to_fit();
        self.counts.shrink_to_fit();
        self.seq_no.shrink_to_fit();
    }

    #[inline(always)]
    fn index_grid_range(&self, kmer: PackedKmer) -> (usize, usize) {
        let g = self.grid_position(kmer);
        (self.index_grid_table[g], self.index_grid_table[g + 1])
    }

    /// Entries as `(packed_kmer, count, insertion_seq)`, in table order.
    pub fn entries(&self) -> Vec<(PackedKmer, u32, u32)> {
        let mut out = Vec::with_capacity(self.number_items);
        let mut readpos = 0usize;
        for idx in 0..self.index_grid_table.len() - 1 {
            while readpos < self.index_grid_table[idx + 1] {
                let kmer = ((idx as u64) << self.log_offset_size) | self.offsets[readpos];
                out.push((kmer, self.counts[readpos], self.seq_no[readpos]));
                readpos += 1;
            }
        }
        out
    }
}

impl LookupTable for GridLookupTable {
    #[inline(always)]
    fn get_count(&self, kmer: PackedKmer) -> u32 {
        let (start, end) = self.index_grid_range(kmer);
        let offset = self.offset(kmer);
        for pos in start..end {
            if self.offsets[pos] == offset && self.counts[pos] != 0 {
                return self.counts[pos];
            }
        }
        0
    }

    fn increase_count(&mut self, kmer: PackedKmer) -> bool {
        let (start, end) = self.index_grid_range(kmer);
        let offset = self.offset(kmer);
        for pos in start..end {
            if self.offsets[pos] == offset && self.counts[pos] != 0 {
                self.counts[pos] += 1;
                return true;
            }
        }
        false
    }

    fn decrease_count(&mut self, kmer: PackedKmer) -> bool {
        let (start, end) = self.index_grid_range(kmer);
        let offset = self.offset(kmer);
        for pos in start..end {
            if self.offsets[pos] == offset && self.counts[pos] != 0 {
                self.counts[pos] -= 1;
                return true;
            }
        }
        false
    }

    fn iterate_over_all(&self, out: &mut dyn Write) -> std::io::Result<()> {
        let mut readpos = 0usize;
        let prefix_len = (self.log_index_size / 2) as u16;
        let suffix_len = (self.log_offset_size / 2) as u16;
        let mut buf = Vec::with_capacity(1 << 20);
        for idx in 0..self.index_grid_table.len() - 1 {
            if readpos >= self.index_grid_table[idx + 1] {
                continue;
            }
            let prefix = packed_kmer_to_string(idx as u64, prefix_len);
            while readpos < self.index_grid_table[idx + 1] {
                let suffix = packed_kmer_to_string(self.offsets[readpos], suffix_len);
                buf.extend_from_slice(prefix.as_bytes());
                buf.extend_from_slice(suffix.as_bytes());
                buf.push(b'\t');
                buf.extend_from_slice(self.counts[readpos].to_string().as_bytes());
                buf.push(b'\n');
                readpos += 1;
                if buf.len() > (1 << 20) - 128 {
                    out.write_all(&buf)?;
                    buf.clear();
                }
            }
        }
        out.write_all(&buf)
    }

    fn len(&self) -> usize {
        self.number_items
    }
}

// ---------------------------------------------------------------------------
// Compact table: same semantics, open addressing
// ---------------------------------------------------------------------------

/// Slot layout. `#[repr(packed)]` keeps this at 12 bytes so a slot never straddles
/// two cache lines more often than it has to; x86-64 handles the unaligned `u64`.
#[derive(Clone, Copy)]
#[repr(C, packed)]
struct Slot {
    key: u64,
    count: u32,
}

/// `u64::MAX` is `GGG...G`, whose reverse complement `CCC...C` is numerically
/// smaller, so a canonical k-mer can never equal it -- which makes it a free
/// "empty" sentinel with no extra occupancy bitmap.
const EMPTY: u64 = u64::MAX;

/// Open addressing replacement for [`GridLookupTable`] with identical observable
/// behaviour, minus the `2^LOGINDEXSIZE` grid array.
pub struct CompactLookupTable {
    slots: Vec<Slot>,
    mask: usize,
    len: usize,
    mode: CountMode,
    log_index_size: u32,
    log_offset_size: u32,
}

#[inline(always)]
fn mix64(mut x: u64) -> u64 {
    // murmur3 finalizer: good avalanche, 5 cheap ops.
    x ^= x >> 33;
    x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
    x ^= x >> 33;
    x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    x ^= x >> 33;
    x
}

impl CompactLookupTable {
    pub fn with_capacity(
        expected: usize,
        log_index_size: u32,
        log_offset_size: u32,
        mode: CountMode,
    ) -> Self {
        // Load factor 0.6: the probe chains stay short without doubling memory.
        let want = ((expected as f64 / 0.6) as usize).max(16);
        let cap = want.next_power_of_two();
        CompactLookupTable {
            slots: vec![Slot { key: EMPTY, count: 0 }; cap],
            mask: cap - 1,
            len: 0,
            mode,
            log_index_size,
            log_offset_size,
        }
    }

    #[inline(always)]
    fn probe(&self, key: u64) -> usize {
        let mut i = (mix64(key) as usize) & self.mask;
        loop {
            let k = self.slots[i].key;
            if k == key || k == EMPTY {
                return i;
            }
            i = (i + 1) & self.mask;
        }
    }

    /// Insert or merge, following `Lookuptable::addElement`'s count-mode rules.
    pub fn add_element(&mut self, kmer: PackedKmer, count: u32) {
        debug_assert_ne!(kmer, EMPTY, "a canonical k-mer can never be all-G");
        if self.len * 10 >= self.slots.len() * 7 {
            self.grow();
        }
        let i = self.probe(kmer);
        if self.slots[i].key == EMPTY {
            self.slots[i] = Slot { key: kmer, count };
            self.len += 1;
        } else {
            let cur = self.slots[i].count;
            self.slots[i].count = match self.mode {
                CountMode::Sum => cur.wrapping_add(count),
                CountMode::Max => cur.max(count),
            };
        }
    }

    fn grow(&mut self) {
        let cap = self.slots.len() * 2;
        let old = std::mem::replace(&mut self.slots, vec![Slot { key: EMPTY, count: 0 }; cap]);
        self.mask = cap - 1;
        for s in old {
            if s.key != EMPTY {
                let key = s.key;
                let i = self.probe(key);
                self.slots[i] = s;
            }
        }
    }

    /// Drop entries at or below `count_threshold`, matching `finalSetupTables`.
    pub fn final_setup_tables(&mut self, count_threshold: u32) {
        if !self.slots.iter().any(|s| s.key != EMPTY && s.count <= count_threshold) {
            return; // nothing to remove; keep the table (and its capacity) as is
        }
        let kept: Vec<Slot> = self
            .slots
            .iter()
            .copied()
            .filter(|s| s.key != EMPTY && s.count > count_threshold)
            .collect();
        let cap = ((kept.len() as f64 / 0.6) as usize).max(16).next_power_of_two();
        self.slots = vec![Slot { key: EMPTY, count: 0 }; cap];
        self.mask = cap - 1;
        self.len = 0;
        for s in kept {
            let key = s.key;
            let i = self.probe(key);
            self.slots[i] = s;
            self.len += 1;
        }
    }

    /// `(kmer, count)` pairs in grid order.
    ///
    /// The C++ orders entries *within* a grid cell by insertion, which an open
    /// addressing table does not preserve. Cells holding more than one k-mer are
    /// therefore ordered by offset here instead. With a `2^30` grid and a few
    /// hundred thousand k-mers such cells are rare but not impossible, so
    /// `counts2flat` uses [`GridLookupTable`] when byte-exact output is required.
    fn ordered_entries(&self) -> Vec<(u64, u32)> {
        let mut v: Vec<(u64, u32)> = self
            .slots
            .iter()
            .filter(|s| s.key != EMPTY)
            .map(|s| (s.key, s.count))
            .collect();
        let log_off = self.log_offset_size;
        v.sort_unstable_by_key(|&(k, _)| (k >> log_off, k & ((1u64 << log_off) - 1)));
        v
    }

    pub fn log_sizes(&self) -> (u32, u32) {
        (self.log_index_size, self.log_offset_size)
    }

    /// Every occupied `(kmer, count)` slot, in arbitrary order.
    pub fn entries(&self) -> impl Iterator<Item = (u64, u32)> + '_ {
        self.slots.iter().filter(|s| s.key != EMPTY).map(|s| {
            let (k, c) = (s.key, s.count);
            (k, c)
        })
    }
}

impl CompactLookupTable {
    /// Lookup reusing a hash the caller already computed, so a sharded table does
    /// not hash twice.
    #[inline(always)]
    pub fn get_count_hashed(&self, kmer: PackedKmer, hash: u64) -> u32 {
        let mut i = (hash as usize) & self.mask;
        loop {
            let s = unsafe { *self.slots.get_unchecked(i) };
            let k = s.key;
            if k == kmer {
                return s.count;
            }
            if k == EMPTY {
                return 0;
            }
            i = (i + 1) & self.mask;
        }
    }
}

impl LookupTable for CompactLookupTable {
    #[inline(always)]
    fn get_count(&self, kmer: PackedKmer) -> u32 {
        let mut i = (mix64(kmer) as usize) & self.mask;
        loop {
            let s = unsafe { *self.slots.get_unchecked(i) };
            let k = s.key;
            if k == kmer {
                return s.count;
            }
            if k == EMPTY {
                return 0;
            }
            i = (i + 1) & self.mask;
        }
    }

    fn increase_count(&mut self, kmer: PackedKmer) -> bool {
        let i = self.probe(kmer);
        if self.slots[i].key == kmer && self.slots[i].count != 0 {
            self.slots[i].count += 1;
            true
        } else {
            false
        }
    }

    fn decrease_count(&mut self, kmer: PackedKmer) -> bool {
        let i = self.probe(kmer);
        if self.slots[i].key == kmer && self.slots[i].count != 0 {
            self.slots[i].count -= 1;
            true
        } else {
            false
        }
    }

    fn iterate_over_all(&self, out: &mut dyn Write) -> std::io::Result<()> {
        let prefix_len = (self.log_index_size / 2) as u16;
        let suffix_len = (self.log_offset_size / 2) as u16;
        let mut buf = Vec::with_capacity(1 << 20);
        for (kmer, count) in self.ordered_entries() {
            let prefix = packed_kmer_to_string(kmer >> self.log_offset_size, prefix_len);
            let suffix = packed_kmer_to_string(kmer & ((1u64 << self.log_offset_size) - 1), suffix_len);
            buf.extend_from_slice(prefix.as_bytes());
            buf.extend_from_slice(suffix.as_bytes());
            buf.push(b'\t');
            buf.extend_from_slice(count.to_string().as_bytes());
            buf.push(b'\n');
            if buf.len() > (1 << 20) - 128 {
                out.write_all(&buf)?;
                buf.clear();
            }
        }
        out.write_all(&buf)
    }

    fn len(&self) -> usize {
        self.len
    }
}

// ---------------------------------------------------------------------------
// HashCountTable: port of src/HashTable.h
// ---------------------------------------------------------------------------

/// Port of the C++ `HashTable`, used when k-mers are counted from the reads rather
/// than read from a DSK file.
///
/// Keeps the C++ contract, oddities included: `get_count` returns **1** for a k-mer
/// that is not in the table, not 0. In the C++ this falls out of `kc_c1_put`
/// inserting the key and reporting `absent`, and it means an unseen k-mer looks
/// like a singleton rather than like an impossible k-mer. Correction thresholds
/// depend on it, so the port keeps it.
#[derive(Default)]
pub struct HashCountTable {
    map: std::collections::HashMap<u64, u32, BuildMix>,
}

#[derive(Default, Clone)]
pub struct BuildMix;
impl std::hash::BuildHasher for BuildMix {
    type Hasher = MixHasher;
    fn build_hasher(&self) -> MixHasher {
        MixHasher(0)
    }
}
pub struct MixHasher(u64);
impl std::hash::Hasher for MixHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = mix64(self.0 ^ b as u64);
        }
    }
    fn write_u64(&mut self, v: u64) {
        self.0 = mix64(v);
    }
}

impl HashCountTable {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn reserve(&mut self, n: usize) {
        self.map.reserve(n);
    }
}

impl LookupTable for HashCountTable {
    #[inline(always)]
    fn get_count(&self, kmer: PackedKmer) -> u32 {
        // See the type docs: absent means 1, matching kc_c1_put's `absent` path.
        *self.map.get(&kmer).unwrap_or(&1)
    }

    fn increase_count(&mut self, kmer: PackedKmer) -> bool {
        *self.map.entry(kmer).or_insert(0) += 1;
        true
    }

    fn decrease_count(&mut self, kmer: PackedKmer) -> bool {
        match self.map.get_mut(&kmer) {
            Some(v) => {
                *v = v.wrapping_sub(1);
                true
            }
            None => {
                // The C++ inserts the key via kc_c1_put and returns before writing
                // a value, leaving it uninitialised; a later getCount would read
                // garbage. Leaving the key out keeps that later read at the
                // well-defined "absent" answer of 1 instead.
                false
            }
        }
    }

    fn iterate_over_all(&self, _out: &mut dyn Write) -> std::io::Result<()> {
        Err(std::io::Error::new(
            std::io::ErrorKind::Unsupported,
            "iterateOverAll is not supported for the internal hash table",
        ))
    }

    fn len(&self) -> usize {
        self.map.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build both table kinds from the same `(kmer, count)` stream.
    fn build_both(
        items: &[(u64, u32)],
        log_index: u32,
        log_offset: u32,
        mode: CountMode,
    ) -> (GridLookupTable, CompactLookupTable) {
        let mut grid = GridLookupTable::new(items.len(), log_index, log_offset, mode);
        grid.track_cross_cell_merges(true);
        for &(k, _) in items {
            grid.assign_kmer_to_grid(k);
        }
        grid.setup_index_grid_table();
        for &(k, c) in items {
            grid.add_element(k, c);
        }
        grid.final_setup_tables(0);

        let mut compact = CompactLookupTable::with_capacity(items.len(), log_index, log_offset, mode);
        for &(k, c) in items {
            compact.add_element(k, c);
        }
        compact.final_setup_tables(0);
        (grid, compact)
    }

    fn sample_items(n: usize, key_bits: u32, seed: u64) -> Vec<(u64, u32)> {
        let mut state = seed | 1;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let key_mask = if key_bits >= 64 { u64::MAX } else { (1u64 << key_bits) - 1 };
        (0..n)
            .map(|_| ((next() & key_mask) & !(1u64 << 63), (next() % 50) as u32 + 1))
            .collect()
    }

    #[test]
    fn grid_and_compact_agree_on_every_key() {
        // weight 12 -> (22, 2), a 4 M-entry grid: small enough for a unit test.
        let items = sample_items(5000, 24, 0xDEAD_BEEF);
        for mode in [CountMode::Sum, CountMode::Max] {
            let (grid, compact) = build_both(&items, 22, 2, mode);
            assert_eq!(grid.len(), compact.len(), "entry counts differ for {mode:?}");
            for &(k, _) in &items {
                assert_eq!(grid.get_count(k), compact.get_count(k), "kmer {k:#x} {mode:?}");
            }
            // Absent keys must read as zero in both.
            for probe in [1u64 << 23, (1u64 << 23) | 1, 0xFFFF] {
                assert_eq!(grid.get_count(probe), compact.get_count(probe));
            }
        }
    }

    #[test]
    fn construction_never_merges_across_grid_cells() {
        // The C++ dedup scan can reach into the previous cell. Assert it does not
        // actually happen, which is what makes the compact table equivalent.
        let items = sample_items(5000, 24, 0x1234_5678);
        let (grid, _) = build_both(&items, 22, 2, CountMode::Sum);
        assert_eq!(grid.cross_cell_merges(), 0);
    }

    #[test]
    fn sum_mode_sums_duplicates() {
        let items = vec![(0x1000u64, 3u32), (0x1000, 4), (0x2000, 5)];
        let (grid, compact) = build_both(&items, 22, 2, CountMode::Sum);
        assert_eq!(grid.get_count(0x1000), 7);
        assert_eq!(compact.get_count(0x1000), 7);
        assert_eq!(grid.get_count(0x2000), 5);
        assert_eq!(grid.len(), 2);
    }

    #[test]
    fn max_mode_maximises_duplicates() {
        let items = vec![(0x1000u64, 3u32), (0x1000, 4), (0x1000, 2)];
        let (grid, compact) = build_both(&items, 22, 2, CountMode::Max);
        assert_eq!(grid.get_count(0x1000), 4);
        assert_eq!(compact.get_count(0x1000), 4);
    }

    #[test]
    fn increase_and_decrease_track_each_other() {
        let items = sample_items(500, 20, 0xABCD);
        let (mut grid, mut compact) = build_both(&items, 22, 2, CountMode::Sum);
        for &(k, _) in items.iter().take(100) {
            assert_eq!(grid.increase_count(k), compact.increase_count(k));
            assert_eq!(grid.get_count(k), compact.get_count(k));
            assert_eq!(grid.decrease_count(k), compact.decrease_count(k));
            assert_eq!(grid.get_count(k), compact.get_count(k));
        }
        // A k-mer that is not present is reported as not found by both.
        let absent = 1u64 << 21;
        assert_eq!(grid.increase_count(absent), compact.increase_count(absent));
    }

    #[test]
    fn iterate_over_all_matches_between_implementations() {
        let items = sample_items(3000, 24, 0x5555);
        let (grid, compact) = build_both(&items, 22, 2, CountMode::Sum);
        let mut a = Vec::new();
        let mut b = Vec::new();
        grid.iterate_over_all(&mut a).unwrap();
        compact.iterate_over_all(&mut b).unwrap();
        assert_eq!(
            String::from_utf8(a).unwrap(),
            String::from_utf8(b).unwrap(),
            "counts2flat output must be byte-identical"
        );
    }

    #[test]
    fn final_setup_drops_entries_at_or_below_threshold() {
        let items = vec![(0x10u64, 1u32), (0x20, 2), (0x30, 3)];
        let mut grid = GridLookupTable::new(items.len(), 22, 2, CountMode::Sum);
        for &(k, _) in &items {
            grid.assign_kmer_to_grid(k);
        }
        grid.setup_index_grid_table();
        for &(k, c) in &items {
            grid.add_element(k, c);
        }
        grid.final_setup_tables(1);
        assert_eq!(grid.len(), 2);
        assert_eq!(grid.get_count(0x10), 0);
        assert_eq!(grid.get_count(0x20), 2);
        assert_eq!(grid.get_count(0x30), 3);
    }

    #[test]
    fn sharded_table_agrees_with_the_single_table() {
        // Sharding must not change a single answer: same keys, same counts, same
        // zero for absent keys, whatever the shard count.
        let items = sample_items(20_000, 40, 0xFACE_B00C);
        for mode in [CountMode::Sum, CountMode::Max] {
            let mut single = CompactLookupTable::with_capacity(items.len(), 30, 34, mode);
            for &(k, c) in &items {
                single.add_element(k, c);
            }
            single.final_setup_tables(0);

            for threads in [1usize, 3, 8, 64] {
                let mut sharded =
                    ShardedLookupTable::with_capacity(items.len(), 30, 34, mode, threads);
                let recs: Vec<crate::dsk::DskCount> = items
                    .iter()
                    .map(|&(k, c)| crate::dsk::DskCount { value: k as u128, abundance: c })
                    .collect();
                // A contiguous pattern makes kmer2min_packed the identity on the
                // low 64 bits, so the two tables get exactly the same keys.
                let t = crate::translator::KmerTranslator::new(&"1".repeat(32)).unwrap();
                let keyed: Vec<(u64, u32)> =
                    recs.iter().map(|r| (t.kmer2min_packed(r.value), r.abundance)).collect();
                sharded.add_batch_parallel(&recs, &t);
                sharded.final_setup_tables(0);

                let mut expect = std::collections::HashMap::new();
                for &(k, c) in &keyed {
                    let e = expect.entry(k).or_insert(0u32);
                    *e = match mode {
                        CountMode::Sum => e.wrapping_add(c),
                        CountMode::Max => (*e).max(c),
                    };
                }
                assert_eq!(sharded.len(), expect.len(), "entry count, {threads} threads");
                for (&k, &want) in &expect {
                    assert_eq!(sharded.get_count(k), want, "kmer {k:#x}, {threads} threads");
                }
                for probe in [1u64 << 62, 3, 0x1234_5678_9ABC] {
                    if !expect.contains_key(&probe) {
                        assert_eq!(sharded.get_count(probe), 0);
                    }
                }
            }
        }
    }

    #[test]
    fn hash_count_table_reports_one_for_absent_keys() {
        // This is the C++ contract, not an accident -- see HashCountTable's docs.
        let mut h = HashCountTable::new();
        assert_eq!(h.get_count(42), 1);
        h.increase_count(42);
        assert_eq!(h.get_count(42), 1);
        h.increase_count(42);
        assert_eq!(h.get_count(42), 2);
    }
}

// ---------------------------------------------------------------------------
// Sharded table: same semantics again, built in parallel
// ---------------------------------------------------------------------------

/// [`CompactLookupTable`] split into independent shards so the table can be built
/// on all cores.
///
/// Building the table is the last part of a run that does not parallelise. On a
/// 1.8 GB input it is about 2.9 s of an 89 s single-threaded run -- invisible
/// there, but it sets the floor once the read loop is spread over 36 cores, and it
/// only grows with the size of the count file.
///
/// A k-mer's shard comes from the *top* bits of its hash and its slot from the
/// *low* bits, so the two are independent and each shard stays well distributed.
/// Every occurrence of a key lands in the same shard, so `COUNT_MODE_SUM` still
/// accumulates correctly across DSK partitions even though the shards are filled
/// concurrently.
///
/// Lookups cost one extra L1-resident indirection over the single table.
pub struct ShardedLookupTable {
    shards: Vec<CompactLookupTable>,
    shard_bits: u32,
    log_index_size: u32,
    log_offset_size: u32,
}

#[inline(always)]
fn shard_of(hash: u64, shard_bits: u32) -> usize {
    if shard_bits == 0 {
        0
    } else {
        (hash >> (64 - shard_bits)) as usize
    }
}

impl ShardedLookupTable {
    pub fn with_capacity(
        expected: usize,
        log_index_size: u32,
        log_offset_size: u32,
        mode: CountMode,
        threads: usize,
    ) -> Self {
        // Four shards per thread keeps the insert phase balanced when the hash
        // distribution is uneven, without making any shard too small to amortise
        // its own allocation.
        let want = (threads.max(1) * 4).next_power_of_two().clamp(1, 4096);
        let shard_bits = want.trailing_zeros();
        let per_shard = expected / want + 16;
        ShardedLookupTable {
            shards: (0..want)
                .map(|_| CompactLookupTable::with_capacity(per_shard, log_index_size, log_offset_size, mode))
                .collect(),
            shard_bits,
            log_index_size,
            log_offset_size,
        }
    }

    /// Translate and insert a whole DSK partition using every available core.
    pub fn add_batch_parallel(
        &mut self,
        records: &[crate::dsk::DskCount],
        translator: &crate::translator::KmerTranslator,
    ) {
        use rayon::prelude::*;
        let nshards = self.shards.len();
        let shard_bits = self.shard_bits;

        // Phase 1: translate to spaced k-mers and bucket by shard. Independent per
        // chunk, so it is a plain parallel map.
        let buckets: Vec<Vec<Vec<(u64, u32)>>> = records
            .par_chunks(64 * 1024)
            .map(|chunk| {
                let mut b: Vec<Vec<(u64, u32)>> =
                    (0..nshards).map(|_| Vec::with_capacity(chunk.len() / nshards + 16)).collect();
                for c in chunk {
                    let k = translator.kmer2min_packed(c.value);
                    b[shard_of(mix64(k), shard_bits)].push((k, c.abundance));
                }
                b
            })
            .collect();

        // Phase 2: insert. Each shard is touched by exactly one thread, so no
        // locking and no atomics.
        self.shards.par_iter_mut().enumerate().for_each(|(s, tbl)| {
            for b in &buckets {
                for &(k, ab) in &b[s] {
                    tbl.add_element(k, ab);
                }
            }
        });
    }

    pub fn final_setup_tables(&mut self, count_threshold: u32) {
        use rayon::prelude::*;
        self.shards.par_iter_mut().for_each(|t| t.final_setup_tables(count_threshold));
    }
}

impl LookupTable for ShardedLookupTable {
    #[inline(always)]
    fn get_count(&self, kmer: PackedKmer) -> u32 {
        let h = mix64(kmer);
        let shard = unsafe { self.shards.get_unchecked(shard_of(h, self.shard_bits)) };
        shard.get_count_hashed(kmer, h)
    }

    fn increase_count(&mut self, kmer: PackedKmer) -> bool {
        let s = shard_of(mix64(kmer), self.shard_bits);
        self.shards[s].increase_count(kmer)
    }

    fn decrease_count(&mut self, kmer: PackedKmer) -> bool {
        let s = shard_of(mix64(kmer), self.shard_bits);
        self.shards[s].decrease_count(kmer)
    }

    fn iterate_over_all(&self, out: &mut dyn Write) -> std::io::Result<()> {
        let log_off = self.log_offset_size;
        let mut v: Vec<(u64, u32)> = self.shards.iter().flat_map(|s| s.entries()).collect();
        v.sort_unstable_by_key(|&(k, _)| (k >> log_off, k & ((1u64 << log_off) - 1)));
        let prefix_len = (self.log_index_size / 2) as u16;
        let suffix_len = (log_off / 2) as u16;
        let mut buf = Vec::with_capacity(1 << 20);
        for (kmer, count) in v {
            buf.extend_from_slice(packed_kmer_to_string(kmer >> log_off, prefix_len).as_bytes());
            buf.extend_from_slice(
                packed_kmer_to_string(kmer & ((1u64 << log_off) - 1), suffix_len).as_bytes(),
            );
            buf.push(b'\t');
            buf.extend_from_slice(count.to_string().as_bytes());
            buf.push(b'\n');
            if buf.len() > (1 << 20) - 128 {
                out.write_all(&buf)?;
                buf.clear();
            }
        }
        out.write_all(&buf)
    }

    fn len(&self) -> usize {
        self.shards.iter().map(|s| s.len()).sum()
    }
}
