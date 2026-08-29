#!/bin/bash
# Wall-clock and peak-RSS comparison, C++ vs Rust, on the committed test set.
set -uo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
COUNTS=${COUNTS:-$ROOT/work/counts.err0.1pct.h5}
W=$ROOT/work/bench
READS=${READS:-$TD/reads.err0.1pct.fq}
REPS=${REPS:-3}
rm -rf "$W"; mkdir -p "$W"
TMP=$(mktemp)

timeit () { # label outdir cmd...
  local label=$1 out=$2; shift 2
  mkdir -p "$out"
  local best="" rss=0
  for _ in $(seq "$REPS"); do
    /usr/bin/time -f "%e %M" -o "$TMP" "$@" >/dev/null 2>/dev/null
    local sec kb
    read -r sec kb < "$TMP"
    [ -z "$best" ] && best=$sec
    best=$(awk -v a="$sec" -v b="$best" 'BEGIN{print (a<b)?a:b}')
    [ "$kb" -gt "$rss" ] && rss=$kb
  done
  printf "%-40s %8ss   %7.2f GB\n" "$label" "$best" "$(awk -v k="$rss" 'BEGIN{printf "%.2f", k/1048576}')"
}

echo "reads: $READS  ($(( $(wc -l < "$READS") / 4 )) records)   best of $REPS   host $(hostname) $(nproc) cores"
echo
echo "=== correction ==="
timeit "C++ (single-threaded by design)" "$W/c1" \
  "$CPP" correction --reads "$READS" --counts "$COUNTS" --outdir "$W/c1" --outprefix t
for t in 1 2 4 8 16 32 72; do
  timeit "Rust --threads $t" "$W/r$t" \
    "$RS" correction --reads "$READS" --counts "$COUNTS" --outdir "$W/r$t" --outprefix t --threads "$t"
done

echo
echo "=== profile ==="
timeit "C++" "$W/pc" "$CPP" profile --reads "$READS" --counts "$COUNTS" --outdir "$W/pc" --outprefix t
timeit "Rust --threads 16" "$W/pr" "$RS" profile --reads "$READS" --counts "$COUNTS" --outdir "$W/pr" --outprefix t --threads 16

echo
echo "=== filter ==="
timeit "C++" "$W/lc" "$CPP" filter --reads "$READS" --counts "$COUNTS" --outdir "$W/lc" --outprefix t
timeit "Rust --threads 16" "$W/lr" "$RS" filter --reads "$READS" --counts "$COUNTS" --outdir "$W/lr" --outprefix t --threads 16

echo
echo "=== counts2flat (dominated by lookuptable build) ==="
timeit "C++" "$W/fc" "$CPP" counts2flat --counts "$COUNTS" --outdir "$W/fc" --outprefix t
timeit "Rust" "$W/fr" "$RS" counts2flat --counts "$COUNTS" --outdir "$W/fr" --outprefix t
rm -f "$TMP"
