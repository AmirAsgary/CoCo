#!/bin/bash
# Large-input benchmark, and a breakdown of where the Rust binary spends its time.
# The committed test set is 53 k reads, which is too small to show which stage
# dominates on real data; this builds a multi-million-read set from the same genome.
set -uo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
W=${W:-$ROOT/work/scale}
DEPTH=${DEPTH:-1500}
mkdir -p "$W"
TMP=$(mktemp)

if [ ! -f "$W/big.fq" ]; then
  echo "generating reads (depth $DEPTH) ..."
  $PY "$TOOLS/make_stress_data.py" --genome "$TD/genome.fa" --out-prefix "$W/big" \
      --seed 42 --contigs 24 --min-depth $((DEPTH/2)) --max-depth "$DEPTH" \
      --sub-rate 0.004 --ins-rate 0.001 --del-rate 0.001
fi
NREADS=$(( $(wc -l < "$W/big.fq") / 4 ))
echo "reads: $NREADS  ($(du -h "$W/big.fq" | cut -f1))"

if [ ! -f "$W/big.counts.h5" ]; then
  echo "counting k-mers ..."
  dsk -file "$W/big.fq" -kmer-size 41 -abundance-min 2 -out "$W/big.counts" \
      -out-dir "$W" -out-tmp "${TMPDIR:-/tmp}" -max-memory 16000 -nb-cores 16 >"$W/dsk.log" 2>&1
fi
SOLID=$(grep -oP '(?<=kmers_nb_solid                          : )\d+' "$W/dsk.log" | head -1)
echo "solid k-mers: ${SOLID:-unknown}"
echo

timeit () { # label cmd...
  local label=$1; shift
  # CoCo creates its own --outdir, but /usr/bin/time needs the parent to exist.
  mkdir -p "$W"
  /usr/bin/time -f "%e %M" -o "$TMP" "$@" >/dev/null 2>/dev/null
  local sec kb; read -r sec kb < "$TMP"
  printf "%-40s %8ss   %7.2f GB\n" "$label" "$sec" "$(awk -v k="$kb" 'BEGIN{printf "%.2f", k/1048576}')"
}

mkdir -p "$W/c"
echo "=== correction ==="
timeit "C++" "$CPP" correction --reads "$W/big.fq" --counts "$W/big.counts.h5" --outdir "$W/c" --outprefix t
for t in 1 4 16 36; do
  mkdir -p "$W/r$t"
  timeit "Rust --threads $t" "$RS" correction --reads "$W/big.fq" --counts "$W/big.counts.h5" \
     --outdir "$W/r$t" --outprefix t --threads "$t"
done
echo
echo "=== byte-identical? ==="
for t in 1 16 36; do
  cmp -s "$W/c/t.corr.reads.fq" "$W/r$t/t.corr.reads.fq" && echo "  PASS threads=$t" || echo "  FAIL threads=$t"
done
rm -f "$TMP"
