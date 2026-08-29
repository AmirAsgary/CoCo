# What differs from the C++, and how that was established

This is the companion to [README.md](README.md). It covers three things: how the
port was verified, where it deliberately departs from the original, and the three
bugs in the C++ that the verification turned up.

The short version: every command produces byte-identical output to the C++
binary, and the departures below are all either measured wins that preserve that
property, or places where the C++ has undefined behaviour and there is no output
to reproduce.

---

## Contents

- [Repository layout](#repository-layout)
- [Module map](#module-map)
- [How it was verified](#how-it-was-verified)
- [Setting up the reference build](#setting-up-the-reference-build)
- [Bugs found in the original C++](#bugs-found-in-the-original-c)
- [Deliberate departures](#deliberate-departures)
- [Performance](#performance)

---

## Repository layout

| path | what |
| --- | --- |
| `coco/` | the Rust crate — library plus the `coco` binary |
| `coco/src/` | the port, one module per C++ translation unit |
| `coco/tests/` | golden-fixture and oracle test suites |
| `oracle/` | a C++ harness exposing internal functions for per-function tests |
| `tools/` | build, test, comparison and benchmark scripts |
| `vendor/`, `.rust/`, `.conda/` | *generated* — isolated toolchain and dependencies |
| `work/` | *generated* — counts files, reference builds, outputs, logs |

Only the first four are in the repository; the rest are produced by
`tools/setup_reference.sh` and the install steps in the README.

There is **no copy of the C++ sources here**. This directory lives inside the CoCo
repository, so `../src` *is* the reference, and both the reference build and the
oracle compile against it directly.

## Module map

The layout mirrors the C++ closely enough to read both side by side.

| Rust | C++ |
| --- | --- |
| `types.rs` | `src/types.{h,cpp}` |
| `kmer.rs` | `src/kmer.{h,cpp}` |
| `translator.rs` | `src/KmerTranslator.{h,cpp}` |
| `lookuptable.rs` | `src/Lookuptable.{h,cpp}`, `src/HashTable.h` |
| `dsk.rs` | GATB-core's `Storage` / `Partition<Count>` |
| `seq.rs` | `lib/kseq/kseq.h`, `src/KSeqWrapper.cpp`, `src/SequenceInfo.h` |
| `countprofile.rs` | `src/CountProfile.{h,cpp}` |
| `sliding.rs` | *(new)* sliding-window extrema |
| `runner.rs` | `src/runner.{h,cpp}` |
| `options.rs` | `src/Options.{h,cpp}` |
| `filehandling.rs` | `src/filehandling.{h,cpp}` |
| `preprocessing.rs` | `src/preprocessing.{h,cpp}` |
| `commands/` | `src/{correction,filter,abundanceEstimator,profile,counts2flat}.cpp` |

Every function in `src/` is translated, including
`checkForSpuriousTransitionDrops` and `checkForSpuriousTransitionDropsWithWindow`,
which no command calls, and the `filehandling.cpp` helpers nothing references.
`preprocessing.cpp`'s `isValid` is the one exception: it is a debug-only
self-check, and `tests/golden_lookuptable.rs` verifies the same property more
strictly by recomputing every count independently.

---

## How it was verified

**329 automated checks, all passing.** Four layers, cheapest first, each
answering a question the previous one cannot.

### 1. Unit tests — 48

Pure functions against hand-derived values, and against slower reference
implementations of themselves:

- the `PEXT` spaced gather pinned to the scalar loop it replaces, over 20,000
  random values per pattern;
- the branch-free reverse complement pinned to the C++'s scalar loop, for every
  k from 1 to 32;
- the monotonic deque pinned to a naive `max_element` over random
  push/pop sequences;
- the open-addressing count table pinned to the grid table it replaces, and the
  sharded table pinned to both, at 1, 3, 8 and 64 shards.

### 2. Golden fixtures — 11

The committed expectations in `testdata/`:

- **`kmer2packedkmer_golden.tsv`** — 1,206 `(pattern, kmer, packed)` triples.
  This is the first thing a port should get right: pure, hot, and everything else
  is built on it.
- **`dsk2ascii` output** — all 395,391 solid k-mers, checked value and abundance,
  in the same partition order the C++ iterator uses. This pins the HDF5 record
  decoding, the 2-bit base order and the iteration order at once, against DSK's
  own dump rather than against the port itself.
- **`expected/corrections.tsv`** — every base the C++ changed. The testdata README
  sets the bar explicitly: *not "similar recall" but the same 7,911 corrections at
  the same positions*. Also checked are the 38 residual errors and the
  recall/precision summary in `coco_stats.txt`.

### 3. Oracle harness — 12

Several functions cannot be reached from the command line at all:
`calcNeighborhoodTolerance` is file-static, and the two drop-check functions above
have no caller anywhere in the released code. End-to-end diffing therefore cannot
pin them.

`oracle/oracle_main.cpp` `#include`s `CountProfile.cpp` rather than linking it, so
the file-static function is in scope, and redefines `private` before that include
so a profile can be installed without going through `fill()`. It prints results as
text; the Rust tests feed both implementations the same generated inputs and
require identical output.

Covered: `revComplement`, `minIndex`, `packedKmer2String`, the translator's
geometry, `getAvgQual`, `maximize`, `calcNeighborhoodTolerance`, `calcXquantile`,
all three drop checks, `doSubstitutionCorrection`, `doIndelCorrection`,
`doTrimming`, and the `filehandling.cpp` path helpers.

### 4. Differential testing — 245 + 13

Both binaries, every command, byte-diffed.

`tools/compare_all.sh` runs the committed test set: 13 comparisons covering
single- and multi-threaded correction, paired-end, FASTA input, all four other
commands, and the internal k-mer counting path.

`tools/differential_test.sh` generates what the fixtures cannot supply. **The
committed dataset applies zero indel corrections** — its coverage is a uniform
20×, so the count profile never dips and `doIndelCorrection` is never entered.
The generator cuts the reference into contigs with abundances spanning two orders
of magnitude, producing the ragged coverage that real metagenomes have, and injects
insertions and deletions as well as substitutions, plus `N`s, soft-masking, varied
read lengths and a realistic quality ramp.

Eight such datasets are swept across eleven option sets, four filter thresholds,
both compression formats, paired-end input and four thread counts: **245
comparisons, exercising 2.1 M insertion corrections, 1.45 M deletion corrections
and 178 k trimmed bases.** Two of the three bugs below were found here.

---

## Setting up the reference build

The comparison layers need a working C++ binary. The `lib/gatb-core` submodule is
empty in a fresh checkout, so the original does not build as-is.

```bash
tools/setup_reference.sh          # copies ../src etc. and clones gatb-core
tools/submit.sh coco-ref    18 120 build_reference.sh
tools/submit.sh coco-oracle 18  90 build_oracle.sh
```

`setup_reference.sh` checks out the `lib/gatb-core` submodule (empty in a fresh
checkout, which is why the original does not build as-is) and creates
`work/reference_fixed/` — the same sources with one function corrected, see the
`HashTable::getCount` bug below. That variant is what the internal-counting path
is compared against, since the unmodified binary's output there is not a
specification. Both binaries land under `work/`.

The counts file is not committed. Regenerate it with:

```bash
export PATH=$PWD/.conda/dsk/bin:$PATH
dsk -file ../testdata/reads.err0.1pct.fq -kmer-size 41 -abundance-min 2 \
    -out counts.err0.1pct -out-dir work -out-tmp $TMPDIR -max-memory 4000 -nb-cores 4
dsk2ascii -file work/counts.err0.1pct.h5 -out work/solid_ascii.txt
```

---

## Bugs found in the original C++

None of these is reachable from the committed test set, which is why the fixtures
alone were never going to be enough.

### A boolean flag in last position segfaults

`Options::parseOptions` checks that an option's argument exists —
`argIdx == argc-1` — only for non-boolean options. A boolean goes straight to
`argv[argIdx+1][0]`, so when the flag is last it reads `argv[argc]`, which C
guarantees to be `NULL`.

```
coco correction --reads r.fq --counts c.h5 --update-lookup   # segfault
coco correction --reads r.fq --update-lookup --counts c.h5   # fine
```

All three boolean options are affected: `--update-lookup`, `--soft`, `--aligned`.
A single 150 bp read reproduces it. The port treats a missing next argument as
"flag with no value", i.e. true, which is what the C++ intends when the flag is
not last.

This is why `--update-lookup` has no differential coverage against the *unmodified*
C++: it cannot be invoked that way without crashing. The suite passes
`--update-lookup true` instead.

### `HashTable::getCount` reads uninitialised memory

`getCount` calls `kc_c1_put`, which **inserts** the k-mer. khashl's map `put`
copies a bucket struct whose `val` field was never initialised, so the k-mer lands
in the table carrying stack garbage as its count. The first lookup of an absent
k-mer correctly returns 1; every later lookup of that same k-mer returns the
garbage.

It is deterministic in practice — three runs give identical output — but it
changes results: on a 10,000-read sample the C++ applies **seven insertion
corrections that have no basis in the data**. The path is reached whenever no
`--counts` file is given.

The port implements the contract the code plainly intends: absent means 1,
consistently. `reference_fixed/` is the C++ with exactly that one function
corrected, and the port is byte-identical to it for `correction`, `profile`,
`abundance` and `filter`. Against the unmodified C++ on that path, expect small
differences.

### A negative shift count that only works by accident

`tryDeletionCorrection` indexes `_inverse_mask_array` at a span position the
pattern may not cover, gets `UCHAR_MAX` back, and computes
`2 * (weight - 255 - 1)` — a negative shift count, converted to `uint64_t` and
handed to a shift instruction. On x86-64 the count is masked to six bits, so
instead of being visibly undefined it lands on a different base and quietly
changes which corrections get approved.

This is latent rather than live: it does not crash and it is deterministic on
x86-64, but it is undefined behaviour and would behave differently on an
architecture whose shifts saturate rather than wrap.

The port reproduces it, because guarding it instead measurably shifts the
insertion and deletion counters — a `tryDeletionCorrection` that reports
"ambiguous" also cancels the competing insertion hypothesis. See `mutation_shift`
in `countprofile.rs`.

### And one in the port

The same method caught a bug in the port: `filter()` overwrites `opt.threshold`
to `0.1` *before* calling `parseOptions`, so the filter command means something
different by `--threshold` than every other command does. The port had `0.01`.
Invisible on the committed fixtures, whose uniform coverage produces no spurious
drops either way; caught immediately by the differential sweep.

---

## Deliberate departures

Everything not listed here is intended to be bit-identical.

### 1. `--threads` actually works

The C++ has a `--threads` option that no command accepts and no code path uses;
`processReads` is a plain loop and `profile`/`abundance` `assert(threads == 1)`.

Here the per-read work runs as a three-stage pipeline — a reader thread, a worker
pool, a writer thread — over bounded channels. Batches carry a sequence number and
the writer restores input order, so parallel output is byte-identical to
sequential output. That property is what makes the parallel path testable against
the C++ at all, and the differential suite checks it explicitly at 1, 3, 8 and 32
threads.

`--threads` defaults to the machine's parallelism rather than to 1.
`--update-lookup` and `--verbose 4` both force it back to 1.

### 2. The count table is a hash table

The C++ `Lookuptable` indexes a `2^30`-entry grid — an 8.6 GB allocation for the
default pattern, regardless of how many k-mers it holds — and then scans a bucket.

The default here is open addressing over the actual k-mer count, which fits in
cache and removes a dependent load and most of the TLB pressure from the hottest
operation in the program. Above one thread it is split into shards by the *top*
bits of the hash (the slot index comes from the low bits, so the two are
independent) and built in parallel.

The grid is still implemented, and `counts2flat` uses it, because that command's
output order *is* the grid's layout. The three implementations are tested against
each other on the real dataset.

### 3. Bounds checks where the C++ has undefined behaviour

`checkForSpuriousTransitionDropsWithWindowNew` indexes `profile[idx + windowSize]`
without checking and reads past the end for any read shorter than
`2*windowSize + span - 1`. With the default `--skip 10` such reads do reach it.
The port reports them as unfiltered.

`tryDeletionCorrection` also indexes `seq[seqPos-1]` unchecked, and
`checkForSpuriousTransitionDrops` divides by the previous position's count — SIGFPE
whenever that count is zero. Both are guarded here.

### 4. `calcXquantile`'s position filter

Given an explicit position list, the C++ counts how many positions are in range
and then copies that many entries from the *front* of the list without
re-filtering, so an out-of-range position early in the list makes it read past the
profile. The port reproduces the selection and returns 0 for the out-of-range read
rather than reading out of bounds.

### 5. Startup and error output

The version banner is printed before options are parsed, as in the C++, so it
appears even at `--verbose 0`. Parse errors therefore follow the banner in both.

---

## Performance

2.94 M reads, 1.8 GB of FASTQ, 14.7 M solid k-mers. Best of three on one Raven
GPU node (Xeon Platinum 8360Y).

| | wall clock | peak RSS | vs. C++ |
| --- | --- | --- | --- |
| C++ (single-threaded by design) | 250.7 s | 8.15 GB | — |
| Rust `--threads 1` | 93.6 s | 0.52 GB | 2.7× |
| Rust `--threads 8` | 12.5 s | 0.62 GB | 20× |
| Rust `--threads 36` | 4.1 s | 0.79 GB | 61× |

On the committed 53,333-read set the single-threaded gap is wider — 6.44 s against
0.84 s — because the C++ pays for its fixed 8.6 GB allocation there regardless of
how few k-mers it holds. The 2.94 M-read figures are the ones to quote.

### Where it comes from

Not one trick. In rough order of contribution:

- **The read loop is a pipeline**, not a read → compute → write cycle per batch.
  The per-batch form left about 1.5 s of serialised I/O that no core count could
  hide; overlapping the reader and writer with the workers took 36 threads from
  8.1 s to 5.3 s. Batches are recycled back to the reader, so the per-record
  buffers are allocated once rather than once per batch.
- **The count table**, as above. This is most of the single-threaded win and all
  of the memory win.
- **Two `PEXT` instructions replace the spaced gather.** `kmer2packedKmer` walks
  32 mask positions per k-mer in the C++; `PEXT` compacts the selected bits in one
  operation, preserving their relative order, which is exactly what the scalar
  loop's shift-left accumulation produces. Detected at runtime, with the scalar
  path kept for machines without BMI2 and for AMD Zen 1/2, where `PEXT` is
  microcoded and slower than the loop it would replace.
- **O(1) sliding-window extrema.** `calcNeighborhoodTolerance` calls
  `max_element` over a 42-element window at every position and is rebuilt on every
  correction round; a monotonic deque makes it amortised constant. It has to
  support an arbitrary interleaving of push and pop, because the original's window
  is not a clean sliding one — for a stretch of positions it pushes a value already
  in the window, leaving a duplicate. That is reproduced, not tidied up.
- **`maximize` over runs.** The default pattern's 32 informative positions form 10
  contiguous runs, so the scatter becomes 10 contiguous max-updates that vectorise
  instead of 32 strided scalar ones.
- **Building the table in parallel**, above one thread — worth 4.59 s → 4.10 s at
  36 threads, and rather more on commands that do less work per read.
- **Fewer allocations.** The C++ allocates a `SequenceInfo` and two `std::string`s
  per read and a fresh `maxProfile` per correction round; all of that is reused.

### Measured and rejected

- **The run is compute-bound, not I/O-bound.** Reading the 1.8 GB input takes
  0.17 s and copying it 0.37 s, against 4.1 s of correction.
- **The parallel table build hurts at one thread**, so the single-table build is
  kept there. Cross-job single samples had suggested the parallel build was a
  regression everywhere; running both variants in one job with three repetitions
  showed that was node-to-node variance, not signal.
- **A GPU port was not attempted.** The inner loop is a latency-bound random
  gather into a multi-gigabyte table followed by data-dependent branching per read
  — the shape GPUs are worst at — and it would put byte-exactness at risk for no
  expected gain.

---

## Status

| command | byte-identical to the C++ | notes |
| --- | --- | --- |
| `correction` | yes | single- and multi-threaded, paired and unpaired |
| `filter` | yes | |
| `abundance` | yes | |
| `profile` | yes | |
| `counts2flat` | yes | uses the grid table to match output order |
| `consensus` | n/a | not implemented in the C++ either |
