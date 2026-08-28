#!/usr/bin/env python3
"""Generate CoCo's test dataset: one reference, two read sets, exact error labels.

The two read sets are the SAME reads in the SAME order with the SAME quality
strings. They differ in the sequence lines and nowhere else, so

    diff <(awk 'NR%4==2' reads.perfect.fq) <(awk 'NR%4==2' reads.err0.1pct.fq)

is exactly the set of introduced errors, and a corrector's job is to turn the
second file back into the first. That makes scoring a byte comparison rather than
an alignment.

Design choices worth knowing when porting:

* Q30 everywhere, error rate 0.001. Phred Q = -10*log10(p), so Q30 IS p=0.001.
  The quality string is therefore honest rather than decorative, which matters
  because CoCo reads it (getAvgQual in CountProfile::doIndelCorrection).
* Errors are substitutions only, drawn uniformly from the three other bases.
  Indels are deliberately excluded: they change read length, so the positional
  identity between the two files would break and every downstream check would
  need an aligner. A separate indel fixture can be added later if wanted.
* Reads are drawn from both strands, uniformly over start positions.
* Fixed seed, so the output is byte-reproducible.
"""
import argparse, gzip, os, random

COMP = str.maketrans("ACGT", "TGCA")


def read_fasta(path):
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    break
                name = line[1:].split()[0]
            elif name:
                chunks.append(line)
    return name, "".join(chunks).upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, help="reference FASTA to draw from")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ref-len", type=int, default=400000)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--coverage", type=float, default=20.0)
    ap.add_argument("--error-rate", type=float, default=0.001)
    ap.add_argument("--qual", type=int, default=30, help="Phred; must satisfy 10^(-q/10)=error-rate")
    ap.add_argument("--seed", type=int, default=20260828)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    src_name, src = read_fasta(a.source)
    assert set(src) <= set("ACGT"), "reference must be pure ACGT"
    ref = src[:a.ref_len]
    assert len(ref) == a.ref_len, f"source only has {len(src)} bp"

    nreads = int(round(a.ref_len * a.coverage / a.read_len))
    rng = random.Random(a.seed)
    qline = chr(33 + a.qual) * a.read_len
    others = {b: [x for x in "ACGT" if x != b] for b in "ACGT"}

    gpath = os.path.join(a.outdir, "genome.fa")
    with open(gpath, "w") as fh:
        fh.write(f">{src_name} test reference, first {a.ref_len} bp\n")
        for i in range(0, len(ref), 70):
            fh.write(ref[i:i + 70] + "\n")

    p_path = os.path.join(a.outdir, "reads.perfect.fq")
    e_path = os.path.join(a.outdir, f"reads.err{a.error_rate:g}.fq")
    t_path = os.path.join(a.outdir, "errors.tsv")
    n_err = 0
    maxstart = a.ref_len - a.read_len
    with open(p_path, "w") as pf, open(e_path, "w") as ef, open(t_path, "w") as tf:
        tf.write("read\tread_pos\tref_pos\tstrand\ttrue_base\tobserved_base\n")
        for i in range(nreads):
            start = rng.randint(0, maxstart)
            rev = rng.random() < 0.5
            true = ref[start:start + a.read_len]
            if rev:
                true = true.translate(COMP)[::-1]
            obs = list(true)
            for j in range(a.read_len):
                if rng.random() < a.error_rate:
                    nb = rng.choice(others[true[j]])
                    obs[j] = nb
                    refpos = (start + a.read_len - 1 - j) if rev else (start + j)
                    tf.write(f"r{i:07d}\t{j}\t{refpos}\t{'-' if rev else '+'}\t{true[j]}\t{nb}\n")
                    n_err += 1
            name = f"r{i:07d}"
            pf.write(f"@{name}\n{true}\n+\n{qline}\n")
            ef.write(f"@{name}\n{''.join(obs)}\n+\n{qline}\n")

    tot = nreads * a.read_len
    print(f"reference   : {gpath}  ({src_name}, {a.ref_len:,} bp)")
    print(f"perfect     : {p_path}  ({nreads:,} reads)")
    print(f"errored     : {e_path}")
    print(f"error truth : {t_path}")
    print(f"reads {nreads:,}  bases {tot:,}  coverage {tot/a.ref_len:.1f}x")
    print(f"errors introduced {n_err:,}  = {100*n_err/tot:.4f} %  (target {100*a.error_rate:g} %)")


if __name__ == "__main__":
    main()
