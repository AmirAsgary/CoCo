# How CoCo works

CoCo corrects sequencing errors without an assembly and without an alignment. It
uses one observation: in data sequenced to reasonable depth, a k-mer that occurs
once or twice is almost certainly the product of a sequencing error, while a
k-mer that occurs twenty times is almost certainly real. Turn that into a per-base
signal and errors become visible as dips.

This document walks the whole pipeline, in the order the code runs it. Function
names in `monospace` refer to the C++; the Rust port keeps the same names in snake
case.

---

## Contents

1. [The alphabet](#1-the-alphabet)
2. [Spaced k-mers](#2-spaced-k-mers)
3. [The count table](#3-the-count-table)
4. [The count profile](#4-the-count-profile)
5. [Maximization: from k-mers to bases](#5-maximization-from-k-mers-to-bases)
6. [Neighbourhood tolerance](#6-neighbourhood-tolerance)
7. [Substitution correction](#7-substitution-correction)
8. [The correction pipeline](#8-the-correction-pipeline)
9. [Indel correction](#9-indel-correction)
10. [Trimming](#10-trimming)
11. [The other commands](#11-the-other-commands)
12. [Parameters, and how they interact](#12-parameters-and-how-they-interact)
13. [A worked example](#13-a-worked-example)

---

## 1. The alphabet

Bases are packed two bits each, in an order that is neither alphabetical nor the
usual `ACGT`:

```
A = 0    C = 1    T = 2    G = 3
```

The ordering is chosen so that **complementing a base is `x XOR 2`**: A↔T is 0↔2,
C↔G is 1↔3. Reverse-complementing a packed k-mer is then "reverse the 2-bit
groups, then flip bit 1 of each", with no lookup table.

Anything that is not `ACGTacgt` is invalid. (The C++ table has a quirk here that
the port preserves: characters in the range `[`–`z` that are not `acgt` decode as
`A` rather than as invalid, so lowercase `n` is silently read as `A` while
uppercase `N` is correctly rejected. See MODIFICATIONS.md.)

---

## 2. Spaced k-mers

### The pattern

A contiguous k-mer of length 41 would use all 41 bases. CoCo instead uses a
**spaced seed**: a fixed binary pattern over a window of 41 positions, of which
only 32 are read.

```
11110111111011011101010111011011111101111
^^^^ ^^^^^^ ^^ ^^^ ^ ^ ^^^ ^^ ^^^^^^ ^^^^
```

- **span** = 41, the width of the window
- **weight** = 32, the number of informative positions
- the gaps sit at offsets 4, 11, 14, 18, 20, 22, 26, 29, 36
- the pattern is **symmetric**, and starts and ends with `1`

Symmetry is required, and it is not cosmetic: it is what makes a spaced k-mer's
reverse complement a spaced k-mer of the same pattern, so a k-mer and its reverse
complement can be reduced to one canonical value.

### Why gaps

Consider a single wrong base at read position *p*. With a contiguous k-mer, every
one of the k k-mers covering *p* is destroyed, uniformly and indistinguishably.
With a spaced seed, a k-mer window that happens to place one of its **gaps** on
*p* reads the same as it would have without the error — it is still abundant.

That gives CoCo two things a contiguous k-mer cannot:

1. **A sharper signal.** Some k-mers overlapping the error stay high, so the
   depression is localised to the base rather than smeared over k positions.
2. **Attribution.** When a read has two errors, the set of k-mers depressed by
   each is different. CoCo exploits this directly (§7) to find k-mers that are
   explained by exactly one candidate error.

### Translation

Extracting a spaced k-mer means gathering the 32 informative bit-pairs out of the
41-position window and packing them contiguously:

```
window   A C T G ? A A C ...        ? = gap position, not read
         │ │ │ │   │ │ │
packed   A C T G   A A C ...        32 bases, 64 bits
```

`kmer2packedKmer` does this. The port does the same gather with two `PEXT`
instructions.

### Canonicalisation

A read and its reverse complement must give the same k-mer, or coverage would be
split between strands. `minIndex` takes the numerically smaller of the packed
value and its reverse complement.

Note the order: **the spaced positions are extracted first, then the 32-symbol
result is canonicalised.** This is well defined precisely because the pattern is
symmetric.

---

## 3. The count table

### Input

CoCo does not count k-mers itself by default — it reads a table produced by
[DSK](https://github.com/GATB/dsk). DSK's "solid" set is every canonical **41-mer**
occurring at least twice, with its abundance. The k-mer size DSK was run with must
equal the pattern's span.

### Building the table

Each 41-mer is put through the spaced translation, giving a 32-symbol packed
k-mer, and its abundance is added to that key's total.

Different 41-mers can produce the same spaced k-mer — they differ only at gap
positions. `--count-mode` decides what happens then: `0` sums the abundances
(the default; the spaced k-mer's true count is the sum over all 41-mers that
realise it), `1` takes the maximum. On the test dataset 395,391 solid 41-mers
collapse to 395,179 distinct spaced k-mers, so 212 collide.

### Storage

The C++ splits the 64-bit key into a 30-bit grid index and a 34-bit offset, and
allocates an array with one slot per grid index — `2^30` words, 8.6 GB, whatever
the k-mer count. The port uses an open-addressing hash table sized to the actual
number of k-mers instead, which is the same map with a much better constant. Both
are implemented and tested against each other.

### Counting internally

Without `--counts`, CoCo counts spaced k-mers from the reads directly into a hash
table. This path has a deliberate oddity: **a k-mer absent from the table reads
back as count 1, not 0.** An unseen k-mer looks like a singleton rather than like
something impossible, which keeps the correction thresholds meaningful.

---

## 4. The count profile

For a read of length *L*, there are `L - span + 1` k-mer start positions. The
**count profile** is the table count at each of them.

```
read     ......................................................
k-mers   [====== 41 ======]
          [====== 41 ======]
           [====== 41 ======]                 ...
profile   c₀  c₁  c₂  c₃  ...                 L-40 values
```

`fill` walks the read once, maintaining the packed window incrementally, and looks
up each position. A position is marked **invalid** if the window contains a
non-nucleotide at an *informative* offset — an `N` landing in one of the nine gaps
does not invalidate the k-mer, because it is never read.

A perfect read from a 20× region gives a flat profile near 20. One wrong base
gives a run of low counts.

---

## 5. Maximization: from k-mers to bases

The profile is indexed by k-mer *start*, but corrections happen at *bases*.
`maximize` converts between them:

> **maxProfile[p] = the largest count among all k-mers that read position p at an
> informative offset.**

The array has one entry per base of the read (`profile_length + span - 1`), floored
at 1.

This is where the spaced seed pays off. If position *p* is wrong, every k-mer that
*reads* it is depressed — and those are exactly the ones the maximum is taken over,
so `maxProfile[p]` collapses. Positions on either side are still covered by k-mers
that place a gap on *p*, or that do not overlap it at all, so they stay high.

A correct read's `maxProfile` is a plateau. An error is a narrow notch.

```
maxProfile   ▁▁▁████████████████▁████████████████▁▁▁
                                ^
                                one wrong base
```

---

## 6. Neighbourhood tolerance

"Low" has to be relative — coverage varies along a genome and between organisms.
`calcNeighborhoodTolerance` computes, for every base position, the count a k-mer
must beat to count as real:

```
tolerance[p] = round(threshold × max(neighbourhood of p)) + pseudocount
```

where the neighbourhood is a roughly span-wide window of `maxProfile` around *p*,
and `threshold` defaults to 0.01. In a region covered 60×, the tolerance is
`round(0.6) + 1 = 2`: a k-mer seen once or twice there is noise.

There is one guard. If the local maximum is below `--lowerbound` (default 5),
**the tolerance is 0** and nothing at that position can be flagged. Coverage that
thin cannot distinguish an error from a genuinely rare sequence, so CoCo declines
to guess.

The `+ pseudocount` term is what stops a repeated sequencing error — the same
mistake made independently in several reads, which does happen — from looking
real.

---

## 7. Substitution correction

`doSubstitutionCorrection` is the core of the tool. Four steps.

### Step 1 — find candidate positions

Every base position where `maxProfile[p] <= tolerance[p]` is an error candidate.
At most 63 are tracked per read; beyond that the read is abandoned as
`TOO_MANY_ERRORS`.

### Step 2 — attribute k-mers to candidates

For each candidate *i* at position *p*, and each informative offset *m* in the
pattern, the k-mer starting at `p - m` reads position *p*. That k-mer is recorded
as **affected by candidate i**, in a 64-bit mask per k-mer position — one bit per
candidate.

```
affected[j] = { i : k-mer starting at j reads the position of candidate i }
```

### Step 3 — find uniquely affected k-mers

A k-mer position *j* is **uniquely affected** by candidate *i* when
`affected[j] == {i}` exactly: this k-mer is depressed by that error and by nothing
else. Fixing the base should bring it back — and if it does not, the hypothesis was
wrong.

CoCo takes the **first** and **last** such k-mer. They overlap the error at
different offsets and therefore read different flanking sequence, so requiring both
to recover is close to two independent tests.

If a candidate has only one uniquely affected k-mer, the evidence is thin. The
first correction pass skips those (`needMultipleKmers`); a later pass allows them.

### Step 4 — try the three alternatives

`firstLastUniqueKmerCorrectionStrategy` rebuilds the two k-mers from the current
read, then for each of the three bases other than the one present:

- patch the two-bit field corresponding to position *p* in **both** k-mers;
- look both up in the table;
- count how many now exceed `tolerance[p]`.

A base that lifts **both** k-mers above tolerance is a candidate correction.

- exactly one such base → apply it
- **two or more → ambiguous, change nothing.** More than one alternative reading
  as real usually means a real variant or a repeat, not an error.
- none → leave the base alone

When a base is changed, its quality score is replaced by the mean of its two
neighbours' (`getAvgQual`) — the original call is no longer meaningful.

### Return value

The function reports whether it corrected everything it found, some of it, or
nothing, which is what drives the loop in the next section.

---

## 8. The correction pipeline

`doCorrection` runs per read:

```
maximize
┌── repeat ────────────────────────────────────────────────┐
│  substitution correction, requiring multiple k-mers      │
│  if anything changed: re-maximize                        │
└── until nothing more is corrected ───────────────────────┘

indel + edge-substitution correction
if anything changed: re-maximize

if the read was not already error-free:
┌── repeat ────────────────────────────────────────────────┐
│  substitution correction, single k-mer allowed           │
│  if anything changed: re-maximize                        │
└── until nothing more is corrected ───────────────────────┘

optional trimming

if total corrections > --max-corr-num: revert the entire read
```

Three things are worth drawing out.

**Why it iterates.** Correcting one error changes the counts of every k-mer
covering it, which can reveal a second error that the first was masking. The loop
runs until a pass finds nothing.

**Why the second pass is more permissive.** The first pass insists on two
uniquely-affected k-mers. That is the right default, but it cannot reach errors
near a read's ends, where fewer k-mers overlap. Once the confident corrections are
in and the profile has been rebuilt around them, the remaining candidates are
retried with a single k-mer allowed.

**Why a read can be reverted.** `--max-corr-num` (default 10) caps corrections per
read. A read needing more than that is usually not a badly-sequenced read — it is a
read from a different strain, and "correcting" it would erase a real biological
difference. The whole read is restored, corrections and all.

---

## 9. Indel correction

Insertions and deletions shift everything downstream, so they do not produce a
notch — they produce a **drop that runs to the end of the affected stretch**.
`doIndelCorrection` therefore works on the *unmaximized* profile and looks for runs
of consecutive k-mer positions below tolerance.

Only runs no longer than one span are considered — anything wider is a coverage
feature, not a single indel. Where the run sits determines which hypotheses are
worth testing:

| run position | delete a base at | insert a base at | substitute at |
| --- | --- | --- | --- |
| runs to the read's end | `idx - L + span - 1` | same | same |
| starts at the read's start | `idx` | `idx - 1` | `idx - 1` |
| inside the read, `L ≥ span-1` | `idx` if `L = span-1` | `idx - 1` | — |

### Edge substitution

Tried first, when a substitution position is defined. Near a read end a base is
not pinned by k-mers on both sides, so `doSubstitutionCorrection` cannot judge it.
`edgeSubstitutionCorrection` instead collects every k-mer that reads the position
at all, and hands the first and last to the same first/last strategy from §7.

### Insertion

Hypothesis: a base was inserted that should not be there. `tryInsertionCorrection`
builds three probe k-mers spanning the site **with that base removed** — the
leftmost, the rightmost, and one balanced between them — and requires **all three**
to exceed tolerance. Three probes rather than two because deleting a base shifts
every downstream position, so a coincidental match is easier to come by.

### Deletion

Hypothesis: a base is missing. `tryDeletionCorrection` builds two probe k-mers with
a gap opened at the site, then tries each of the four bases in it, exactly as in
§7. Exactly one base lifting both probes → insert it; more than one → ambiguous.

### Resolving the two

Insertion and deletion are tested independently and then compared:

- insertion approved, deletion not → delete the base
- deletion approved, insertion not → insert the base
- **both approved → ambiguous, change nothing**
- neither → leave it

An ambiguous *deletion* result also cancels the insertion hypothesis: if several
bases would fit a missing base, the region is not trustworthy enough to act on at
all.

> The committed test dataset exercises none of this. Its coverage is a uniform
> 20×, so the profile never drops, and CoCo applies exactly zero indel corrections
> to it. Verifying this path needed generated data with ragged coverage — see
> MODIFICATIONS.md.

---

## 10. Trimming

Off by default (`--max-trim-len 0`). When enabled, `doTrimming` measures the run of
below-tolerance k-mers at each end of the profile, and if a run is short enough,
cuts those bases off. It is the fallback for damage at a read's ends that neither
correction path could resolve — near an end there is not enough flanking sequence
to test a hypothesis, so removing the bases is better than keeping them or guessing.

---

## 11. The other commands

### `filter` — chimeric reads

A chimera is two unrelated fragments ligated together. Both halves are well
covered, but **no k-mer spans the junction**, so the profile shows a deep, narrow
hole exactly one k-mer span wide.

`checkForSpuriousTransitionDropsWithWindowNew` slides two windows — one ahead, one
behind — across the profile, both of width `max(longest run of 1s + 1, 5)`, which
is 7 for the default pattern. A drop starts where the level ahead falls to less
than `--threshold` (default **0.1** for this command) of the level behind, and ends
where it recovers by the same factor.

**A drop narrower than one span means the read is chimeric.** Real low-coverage
regions are wider than a single k-mer; a junction hole is not. A drop that starts
near the end and never recovers is judged the same way.

### `abundance`

Reports the **67th percentile** of the read's count profile. High enough to ignore
the positions depressed by errors, low enough not to be dragged up by a repeat.

### `profile` and `counts2flat`

Developer tools: the raw per-position counts for each read, and the whole spaced
k-mer table as text.

---

## 12. Parameters, and how they interact

| parameter | raising it | lowering it |
| --- | --- | --- |
| `--threshold` | corrects more aggressively; more false positives | corrects only obvious errors |
| `--pseudocount` | corrects less; guards against repeated sequencing errors | corrects more in low coverage |
| `--lowerbound` | refuses to judge more regions; safer on uneven coverage | attempts correction at lower depth |
| `--max-corr-num` | keeps more heavily-edited reads | discards more edits; protects strain variation |
| `--max-trim-len` | removes more unresolvable read ends | keeps reads intact |

The two that matter most are `--threshold` and `--lowerbound`, and they work at
different ends. `--threshold` sets how far below its neighbours a count must fall
to be suspect; `--lowerbound` sets the coverage below which CoCo declines to judge
at all. On uniform high-coverage data, raise `--threshold`. On uneven metagenomic
data, `--lowerbound` is what protects the low-abundance organisms from having their
real sequence "corrected" toward the abundant ones.

`--count-mode 1` (maximum rather than sum) is worth knowing about when different
organisms share sequence that differs only at the pattern's gap positions. Summing
inflates such k-mers; taking the maximum does not.

---

## 13. A worked example

From `testdata/`: 53,333 reads of 150 bp, 20× coverage of a 400 kb genome, with
0.0994 % substitution errors introduced.

**Setup.** k = 41 (the span). DSK finds 630,245 distinct 41-mers, of which 395,391
occur at least twice. Spaced translation collapses them to 395,179 keys.

**A read with one error.** 150 bases → 110 k-mer positions. The profile is near 20
everywhere except a run of low counts where the error sits. After maximization, one
base has `maxProfile` near 1 while its neighbours are near 20.

**Tolerance there.** `round(0.01 × 20) + 1 = 1`. The local maximum, 20, is above
`--lowerbound` 5, so the position is judged. `maxProfile[p] = 1 <= 1` → candidate.

**Attribution.** 32 k-mer positions read that base. Given only one error in the
read, all 32 are uniquely affected by it, so the first and last are span-wide apart.

**The three alternatives.** Two give k-mers absent from the table (count 0, below
tolerance 1). One gives k-mers with counts near 20 in both probes. Exactly one
candidate → apply it.

**The whole file.** 7,911 of the 7,949 introduced errors are corrected — 7,494 in
the multi-k-mer pass and 417 in the single-k-mer pass — with **zero** false
positives and zero miscorrections. The 38 that remain sit where the evidence was
genuinely ambiguous or too close to a read end. Recall 99.52 %, precision 100 %,
residual error rate 0.00048 % against the original 0.0994 %.

Zero indel corrections, because uniform coverage produces no drops for §9 to act on.
