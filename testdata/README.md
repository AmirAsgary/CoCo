# CoCo test dataset

A small, fully-labelled dataset for checking a reimplementation of CoCo against the
C++ original. Everything here is regenerable from `make_testdata.py` with a fixed
seed, so it is byte-reproducible.

Total committed size 33 MB. The DSK counts file is **not** committed (15 MB of
binary HDF5, exactly reproducible in one command — see below).

## What is in it

| file | what |
| --- | --- |
| `genome.fa` | the reference the reads come from: first 400,000 bp of `NC_002162.1` (*Ureaplasma parvum*), a complete closed genome, single sequence, pure ACGT, 25.5 % GC |
| `reads.perfect.fq` | 53,333 reads of 150 bp, **zero errors** — exact substrings of `genome.fa` (or its reverse complement) |
| `reads.err0.1pct.fq` | the **same reads, same order, same names, same quality strings**, with 0.0994 % substitution errors |
| `errors.tsv` | every introduced error: read, position in read, position in reference, strand, true base, observed base |
| `expected/corrections.tsv` | every base the C++ CoCo changed, and whether it was right |
| `expected/residual_errors.tsv` | the 38 bases still wrong after correction |
| `expected/coco_stats.txt` | CoCo's own counters plus the scoring summary |
| `kmer2packedkmer_golden.tsv` | 1,206 `(pattern, kmer, packed_kmer)` vectors for `KmerTranslator::kmer2packedKmer` |
| `make_testdata.py` | the generator |

## The property that makes scoring trivial

The two read sets differ **only in the sequence lines**. Same names, same order,
same quality strings, same lengths. So

```bash
diff <(awk 'NR%4==2' reads.perfect.fq) <(awk 'NR%4==2' reads.err0.1pct.fq)
```

is exactly the set of introduced errors, and a corrector's job is to turn the second
file back into the first. Scoring is a byte comparison at matching positions — no
aligner, no ambiguity about where an error "is".

Errors are **substitutions only**, deliberately. An indel changes read length, which
would break the positional identity between the two files and force every check
through an aligner. See "Known gap" below.

Quality is Phred **Q30 on every base**, and Q30 *is* p = 0.001, so the quality string
is honest rather than decorative. CoCo reads it (`getAvgQual`, used by
`doIndelCorrection`), so a port must too.

## Reproducing the counts file

CoCo takes k-mer counts from DSK. The k-mer size must equal the **span** of the
spaced pattern — 41 for CoCo's default pattern.

```bash
dsk -file reads.err0.1pct.fq -kmer-size 41 -abundance-min 2 \
    -out counts.err0.1pct -out-dir . -out-tmp /tmp -max-memory 4000 -nb-cores 4
```

`-abundance-min 2` is DSK's default and what CoCo expects: singleton k-mers are
overwhelmingly error k-mers. On this dataset it yields 395,391 solid and 234,854
weak k-mers out of 630,245 distinct.

## Reproducing the expected output

```bash
coco correction --reads reads.err0.1pct.fq --counts counts.err0.1pct.h5 \
     --outdir out --outprefix t
```

Expected counters:

```
substitution corrections (multi kmer step): 7494
substitution corrections (single kmer step): 417
insertion corrections: 0
deletion corrections: 0
```

and scored against `reads.perfect.fq`:

```
errors present 7949   changed 7911
TP 7911   WR 0   FN 38   FP 0
recall 99.52 %   precision 100.00 %   residual 0.00048 % (was 0.0994 %)
```

**Zero false positives** on this data. That is the bar for a port: not "similar
recall" but the same 7,911 corrections at the same positions, listed in
`expected/corrections.tsv`.

## Suggested order for a port

1. `kmer2packedKmer` — check against `kmer2packedkmer_golden.tsv` first. It is pure,
   has no dependencies, and is in the hot path.
2. Lookup table construction from the DSK counts, then spot-check counts against
   k-mer occurrences you can count directly in `genome.fa`.
3. `CountProfile::fill` / `maximize` — the count profile of a read.
4. `doSubstitutionCorrection` — then `expected/corrections.tsv` should match exactly.
5. `doIndelCorrection` — **not exercised by this dataset** (see below).

## Known gap: the indel path is not covered

CoCo applies **zero** indel corrections here, because coverage is uniform 20x and
`doIndelCorrection` triggers on count-profile drops. Real metagenomic data has
irregular coverage and CoCo then applies indel corrections at a high rate.

So this dataset exercises the substitution path only. A port should treat
`doIndelCorrection` as untested by these fixtures. Adding an indel-bearing variant
means giving up the equal-length property above, so it belongs in a separate
directory with its own alignment-based scoring.
