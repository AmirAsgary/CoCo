# CoCo-rust

A Rust port of [CoCo](https://github.com/soedinglab/CoCo) — Annika Jochheim's
suite for correcting sequencing errors in short reads using spaced k-mer count
profiles.

It reproduces the C++ binary's output **byte for byte** on every command, runs
2.7× faster on one core and 61× faster on 36, and uses about a tenth of the
memory. Same command-line interface, plus a `--threads` option that actually
does something.

- [ALGORITHM.md](ALGORITHM.md) — how CoCo works, end to end
- [MODIFICATIONS.md](MODIFICATIONS.md) — what differs from the C++, how the port
  was verified, and three bugs found in the original

---

## Install

You need `cargo` (Rust 1.75+) and a C compiler. Nothing else — HDF5 is built
from vendored source, so there is no system HDF5 to install.

```bash
git clone https://github.com/AmirAsgary/CoCo.git
cd CoCo/Coco-rust/coco
cargo build --release
```

The binary lands at `target/release/coco`. Put it on your `PATH`:

```bash
export PATH=$PWD/target/release:$PATH
coco --help
```

### Building without internet

By default `cargo build` fetches dependencies from crates.io. To build on a
machine with no network — an HPC compute node, for instance — vendor them once
from somewhere that does have network access:

```bash
cd Coco-rust/coco
cargo vendor ../vendor > .cargo/config.toml
```

That writes `vendor/` and a `.cargo/config.toml` pointing at it; neither is in the
repository, because the config file would break the build for anyone without the
directory. After that, `cargo build --release` works offline.

### On an HPC cluster

`Coco-rust/tools/submit.sh` submits a build or test script to a SLURM GPU
partition, running a snapshot of `tools/` so that editing a script cannot disturb
a job already in flight:

```bash
tools/submit.sh coco-verify 36 90 full_verify.sh
```

---

## Before you start: the k-mer counts

CoCo needs a table of k-mer counts. The usual source is
[DSK](https://github.com/GATB/dsk), and **the k-mer size must equal the span of
the spaced pattern** — 41 for the default pattern.

```bash
dsk -file reads.fq -kmer-size 41 -abundance-min 2 -out counts -out-dir .
# -> counts.h5
```

`-abundance-min 2` drops singleton k-mers, which are overwhelmingly errors.

You can skip this and let CoCo count k-mers itself by omitting `--counts`, but
that path builds an in-memory hash table of every k-mer in the input and is only
practical for small datasets.

---

## Usage

```
coco <command> [options]
```

| command | what it does |
| --- | --- |
| `correction` | find and fix sequencing errors |
| `filter` | drop chimeric reads |
| `abundance` | estimate each read's coverage |
| `profile` | write each read's raw count profile (developer tool) |
| `counts2flat` | dump the spaced k-mer table as text (developer tool) |

### Correcting reads

```bash
coco correction --reads reads.fq --counts counts.h5 --outdir out --threads 16
# -> out/reads.corr.reads.fq
```

Paired-end:

```bash
coco correction -1 reads_1.fq -2 reads_2.fq --counts counts.h5 --outdir out
# -> out/reads_1.corr.1.fq  and  out/reads_2.corr.2.fq
```

Mates are corrected independently; the pairing only keeps them in step, so a
read too short to profile causes both mates to be passed through untouched.

It prints a summary when it finishes:

```
### COCO ERROR CORRECTION STATISTIC ###
substitution corrections (multi kmer step): 7494
substitution corrections (single kmer step): 417
insertion corrections: 0
deletion corrections: 0
trimmed nucleotides: 0
```

### Filtering chimeras

```bash
coco filter --reads reads.fq --counts counts.h5 --outdir out --threads 16
# -> out/reads.filter.reads.fq   (chimeric reads are omitted)
```

### Estimating abundance

```bash
coco abundance --reads reads.fq --counts counts.h5 --outdir out
# -> out/reads.abundance.reads.tsv   ("<read name>\t<estimate>", or "-" if skipped)
```

### Developer tools

```bash
coco profile     --reads reads.fq --counts counts.h5 --outdir out
coco counts2flat --counts counts.h5 --outdir out
```

`profile` writes one block per read: a `#name` line, then `position<TAB>count`
for every k-mer start. `counts2flat` writes `<spaced k-mer><TAB><count>`.

---

## Options

Run `coco <command> -h` for the list a given command accepts.

### Input and output

| option | default | meaning |
| --- | --- | --- |
| `--reads <file>` | | unpaired, single or merged reads (FASTA/FASTQ, optionally `.gz` or `.bz2`) |
| `-1 <file>` `-2 <file>` | | paired-end reads; give either these or `--reads` |
| `--counts <file>` | | DSK counts file. Omit to count k-mers internally |
| `--outdir <dir>` | `coco_out/` | created if missing |
| `--outprefix <name>` | input file's basename | must not contain a path |
| `--threads <n>` | all cores | see the note below |
| `--verbose <0-4>` | `3` | 0 quiet, 1 errors, 2 +warnings, 3 +info, 4 +per-correction trace |

### Correction behaviour

| option | default | meaning |
| --- | --- | --- |
| `--threshold <float>` | `0.01` | a position is suspect when its count falls to this fraction of its neighbourhood's |
| `--pseudocount <int>` | `1` | added to every tolerance; raise it to correct less aggressively |
| `--lowerbound <int>` | `5` | neighbourhood counts below this are treated as too thin to judge, and nothing is corrected |
| `--max-corr-num <int>` | `10` | a read needing more corrections than this is left untouched — usually a different strain, not errors. `0` disables the cap |
| `--max-trim-len <int>` | `0` | trim up to this many bases from an end whose errors could not be fixed. `0` disables trimming |
| `--skip <int>` | `10` | reads shorter than `skip + span` bases are passed through unprocessed |
| `--update-lookup` | off | feed corrections back into the count table as reads are processed |

### Expert

| option | default | meaning |
| --- | --- | --- |
| `--spaced-pattern <bits>` | 41-position default | must be symmetric, start and end with `1`, span ≤ 64, weight 12–32. **The counts file's k must equal the span** |
| `--count-mode <0\|1>` | `0` | how to combine distinct k-mers that collide under the spaced mask: `0` sum, `1` maximum |

`filter` reads `--threshold` differently — there it is the permitted drop ratio
between two successive windows, and its default is `0.1`, not `0.01`.

### A note on `--threads`

Output is byte-identical at any thread count; results never depend on how many
cores you use.

The exception is `--update-lookup`, which feeds each read's corrections back into
the shared table and so genuinely depends on the order reads are processed. It
forces single-threaded operation. `--verbose 4` does too, so its per-correction
trace stays in input order.

---

## Requirements and limits

- **x86-64.** A `PEXT`-based fast path is used when the CPU supports BMI2 and is
  not Zen 1/2 (where `PEXT` is microcoded); otherwise a portable path runs. Other
  architectures work but have not been tested.
- **Reads must be at least `span` bases long** to be corrected — 41 by default.
  Shorter reads are copied to the output unchanged.
- **`consensus` is not implemented**, exactly as in the C++.
- Memory is roughly 12 bytes per distinct spaced k-mer, plus the input batches in
  flight. A 14.7 M k-mer table runs in well under 1 GB.

---

## Testing

```bash
tools/submit.sh coco-verify 36  90 full_verify.sh        # ~3 min
tools/submit.sh coco-diff   36 300 differential_test.sh  # ~45 min
```

Or directly, without a scheduler:

```bash
cd coco && cargo test --release
```

The full suite compares against a build of the original C++; see
[MODIFICATIONS.md](MODIFICATIONS.md) for what it covers and how to set the
reference up.

---

## Credits

CoCo was written by Annika Jochheim (annika.jochheim@mpinat.mpg.de) at the Max
Planck Institute for Multidisciplinary Sciences, with the `Lookuptable` design
by Martin Steinegger. This is a port of that work and carries the same GPL-3.0
licence.
