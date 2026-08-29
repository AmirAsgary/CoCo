#!/bin/bash
# Run every command through both binaries and diff the outputs.
set -uo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
W=$ROOT/work/compare
COUNTS=$ROOT/work/counts.err0.1pct.h5
rm -rf "$W"; mkdir -p "$W"
fail=0

check () { # name file_cpp file_rs
  if cmp -s "$2" "$3"; then
    echo "PASS  $1"
  else
    echo "FAIL  $1"
    echo "      cpp=$2 rs=$3"
    diff <(head -40 "$2") <(head -40 "$3") | head -20
    fail=1
  fi
}

run () { # label binary outdir extra...
  local label=$1 bin=$2 out=$3; shift 3
  mkdir -p "$out"
  "$bin" "$@" --outdir "$out" --outprefix t > "$out/stdout.txt" 2> "$out/stderr.txt"
  echo "$?" > "$out/exit"
}

echo "### correction (unpaired) ###"
run cpp  "$CPP" "$W/corr_cpp" correction --reads "$TD/reads.err0.1pct.fq" --counts "$COUNTS"
run rs1  "$RS"  "$W/corr_rs1" correction --reads "$TD/reads.err0.1pct.fq" --counts "$COUNTS" --threads 1
run rsN  "$RS"  "$W/corr_rsN" correction --reads "$TD/reads.err0.1pct.fq" --counts "$COUNTS" --threads 16
check "correction 1 thread"   "$W/corr_cpp/t.corr.reads.fq" "$W/corr_rs1/t.corr.reads.fq"
check "correction 16 threads" "$W/corr_cpp/t.corr.reads.fq" "$W/corr_rsN/t.corr.reads.fq"

echo "### correction stats line ###"
for d in corr_rs1 corr_rsN; do
  grep -E "^(substitution|insertion|deletion|trimmed)" "$W/$d/stdout.txt" > "$W/$d/stats.txt"
done
grep -E "^(substitution|insertion|deletion|trimmed)" "$W/corr_cpp/stdout.txt" > "$W/corr_cpp/stats.txt"
check "correction stats 1t"  "$W/corr_cpp/stats.txt" "$W/corr_rs1/stats.txt"
check "correction stats 16t" "$W/corr_cpp/stats.txt" "$W/corr_rsN/stats.txt"

echo "### correction (paired: same file twice) ###"
run cpp "$W/x" "$W/pair_cpp" 2>/dev/null
run cppP "$CPP" "$W/pair_cpp" correction -1 "$TD/reads.err0.1pct.fq" -2 "$TD/reads.perfect.fq" --counts "$COUNTS"
run rsP  "$RS"  "$W/pair_rs"  correction -1 "$TD/reads.err0.1pct.fq" -2 "$TD/reads.perfect.fq" --counts "$COUNTS" --threads 16
check "correction paired r1" "$W/pair_cpp/t.corr.1.fq" "$W/pair_rs/t.corr.1.fq"
check "correction paired r2" "$W/pair_cpp/t.corr.2.fq" "$W/pair_rs/t.corr.2.fq"

echo "### profile ###"
run cpp "$CPP" "$W/prof_cpp" profile --reads "$TD/reads.err0.1pct.fq" --counts "$COUNTS"
run rs  "$RS"  "$W/prof_rs"  profile --reads "$TD/reads.err0.1pct.fq" --counts "$COUNTS" --threads 16
check "profile" "$W/prof_cpp/t.profile.reads.txt" "$W/prof_rs/t.profile.reads.txt"

echo "### abundance ###"
run cpp "$CPP" "$W/ab_cpp" abundance --reads "$TD/reads.err0.1pct.fq" --counts "$COUNTS"
run rs  "$RS"  "$W/ab_rs"  abundance --reads "$TD/reads.err0.1pct.fq" --counts "$COUNTS" --threads 16
check "abundance" "$W/ab_cpp/t.abundance.reads.tsv" "$W/ab_rs/t.abundance.reads.tsv"

echo "### filter ###"
run cpp "$CPP" "$W/filt_cpp" filter --reads "$TD/reads.err0.1pct.fq" --counts "$COUNTS"
run rs  "$RS"  "$W/filt_rs"  filter --reads "$TD/reads.err0.1pct.fq" --counts "$COUNTS" --threads 16
check "filter" "$W/filt_cpp/t.filter.reads.fq" "$W/filt_rs/t.filter.reads.fq"

echo "### counts2flat ###"
run cpp "$CPP" "$W/c2f_cpp" counts2flat --counts "$COUNTS"
run rs  "$RS"  "$W/c2f_rs"  counts2flat --counts "$COUNTS"
check "counts2flat" "$W/c2f_cpp/t.counts2flat.tsv" "$W/c2f_rs/t.counts2flat.tsv"

echo
echo "### FASTA input (genome.fa) ###"
run cpp "$CPP" "$W/fa_cpp" correction --reads "$TD/genome.fa" --counts "$COUNTS"
run rs  "$RS"  "$W/fa_rs"  correction --reads "$TD/genome.fa" --counts "$COUNTS" --threads 8
check "correction fasta" "$W/fa_cpp/t.corr.reads.fa" "$W/fa_rs/t.corr.reads.fa"

echo
echo "### internal k-mer counting (no --counts) ###"
# Compared against reference_fixed, not reference: the unmodified C++ reads
# uninitialised memory in HashTable::getCount on this path (see README), so its
# output is not a specification to match. reference_fixed is the same binary with
# that one function corrected to the contract the code clearly intends.
head -40000 "$TD/reads.err0.1pct.fq" > "$W/small.fq"
if [ -x "$FIXED" ]; then
  run cppfix "$FIXED" "$W/hash_fix" correction --reads "$W/small.fq"
  run rs     "$RS"    "$W/hash_rs"  correction --reads "$W/small.fq" --threads 8
  check "correction internal hashtable (vs reference_fixed)" \
        "$W/hash_fix/t.corr.reads.fq" "$W/hash_rs/t.corr.reads.fq"
else
  echo "SKIP  internal hashtable: build reference_fixed first"
fi

echo
if [ $fail -eq 0 ]; then echo "ALL COMPARISONS PASSED"; else echo "SOME COMPARISONS FAILED"; fi
exit $fail
