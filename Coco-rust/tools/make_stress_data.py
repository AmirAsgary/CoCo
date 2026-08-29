#!/usr/bin/env python3
"""Generate read sets that exercise CoCo paths the committed fixtures do not.

The committed dataset is uniform 20x coverage of one genome with substitution
errors only, which by construction never triggers `doIndelCorrection` -- that code
only fires on count-profile drops, and uniform coverage does not produce them.
This builds metagenome-shaped input instead: the reference is cut into contigs
whose abundances span two orders of magnitude, so coverage is ragged and drops are
everywhere, and the reads carry insertions and deletions as well as substitutions.

Everything is seeded, so a failing differential run can be replayed exactly.
"""
import argparse, random, sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    seqs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--contigs", type=int, default=40,
                    help="pieces the reference is cut into, each with its own abundance")
    ap.add_argument("--min-depth", type=float, default=3.0)
    ap.add_argument("--max-depth", type=float, default=120.0)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--read-len-jitter", type=int, default=0,
                    help="vary read length by +/- this many bases")
    ap.add_argument("--sub-rate", type=float, default=0.005)
    ap.add_argument("--ins-rate", type=float, default=0.0015)
    ap.add_argument("--del-rate", type=float, default=0.0015)
    ap.add_argument("--n-rate", type=float, default=0.0,
                    help="fraction of bases replaced by 'N'")
    ap.add_argument("--lowercase-rate", type=float, default=0.0,
                    help="fraction of bases lowercased (soft-masking)")
    ap.add_argument("--qual", choices=["uniform", "varied"], default="varied")
    ap.add_argument("--format", choices=["fastq", "fasta"], default="fastq")
    ap.add_argument("--paired", action="store_true")
    ap.add_argument("--comments", action="store_true", help="add a description to each header")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    ref = "".join(s for _, s in read_fasta(args.genome)).upper()
    ref = "".join(c if c in "ACGT" else "A" for c in ref)

    # Cut into contigs and give each a depth, so coverage is ragged like a real
    # metagenome rather than flat.
    n = len(ref)
    cuts = sorted(rng.sample(range(1, n), args.contigs - 1))
    bounds = list(zip([0] + cuts, cuts + [n]))
    contigs = [(ref[a:b], rng.uniform(args.min_depth, args.max_depth)) for a, b in bounds]

    reads = []
    for seq, depth in contigs:
        if len(seq) < args.read_len + 2 * args.read_len_jitter + 5:
            continue
        count = max(1, int(depth * len(seq) / args.read_len))
        for _ in range(count):
            rl = args.read_len + rng.randint(-args.read_len_jitter, args.read_len_jitter)
            rl = max(50, min(rl, len(seq)))
            start = rng.randrange(0, len(seq) - rl + 1)
            frag = seq[start:start + rl]
            if rng.random() < 0.5:
                frag = revcomp(frag)
            reads.append(frag)
    rng.shuffle(reads)

    def mutate(s):
        out = []
        for base in s:
            r = rng.random()
            if r < args.del_rate:
                continue                                    # deletion
            if r < args.del_rate + args.ins_rate:
                out.append(rng.choice("ACGT"))              # insertion
                out.append(base)
                continue
            if r < args.del_rate + args.ins_rate + args.sub_rate:
                out.append(rng.choice([c for c in "ACGT" if c != base]))
                continue
            out.append(base)
        s2 = "".join(out)
        if args.n_rate > 0:
            s2 = "".join("N" if rng.random() < args.n_rate else c for c in s2)
        if args.lowercase_rate > 0:
            s2 = "".join(c.lower() if rng.random() < args.lowercase_rate else c for c in s2)
        return s2

    def qual_for(s):
        if args.qual == "uniform":
            return "?" * len(s)
        # A quality ramp that decays toward the 3' end, as real reads do.
        return "".join(chr(33 + max(2, min(40, int(40 - 25.0 * i / max(1, len(s)) + rng.randint(-3, 3)))))
                       for i in range(len(s)))

    def write(path, seqs, tag):
        with open(path, "w") as fh:
            for i, s in enumerate(seqs):
                m = mutate(s)
                if len(m) < 5:
                    m = s
                name = f"r{i:07d}"
                comment = f" depth_tag={tag} len={len(m)}" if args.comments else ""
                if args.format == "fasta":
                    fh.write(f">{name}{comment}\n{m}\n")
                else:
                    fh.write(f"@{name}{comment}\n{m}\n+\n{qual_for(m)}\n")

    ext = "fa" if args.format == "fasta" else "fq"
    if args.paired:
        write(f"{args.out_prefix}.1.{ext}", reads, "r1")
        write(f"{args.out_prefix}.2.{ext}", [revcomp(r) for r in reads], "r2")
        print(f"wrote {len(reads)} read pairs to {args.out_prefix}.{{1,2}}.{ext}", file=sys.stderr)
    else:
        write(f"{args.out_prefix}.{ext}", reads, "se")
        print(f"wrote {len(reads)} reads to {args.out_prefix}.{ext}", file=sys.stderr)


if __name__ == "__main__":
    main()
