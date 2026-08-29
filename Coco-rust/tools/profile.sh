#!/bin/bash
# Byte-identity at scale, plus a profile of where the single-threaded run spends
# its time. The committed test set is too small to show this: there, the C++'s
# fixed 8.6 GB table setup dominates everything else.
set -uo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
W=$ROOT/work/scale
mkdir -p "$W/r1" "$W/r36" "$W/prof"

echo "=== byte-identity at scale ($(( $(wc -l < "$W/big.fq") / 4 )) reads) ==="
$RS correction --reads "$W/big.fq" --counts "$W/big.counts.h5" --outdir "$W/r1"  --outprefix t --threads 1  >"$W/r1/log"  2>&1
$RS correction --reads "$W/big.fq" --counts "$W/big.counts.h5" --outdir "$W/r36" --outprefix t --threads 36 >"$W/r36/log" 2>&1
for d in r1 r36; do
  if cmp -s "$W/c/t.corr.reads.fq" "$W/$d/t.corr.reads.fq"; then echo "  PASS $d"; else echo "  FAIL $d"; fi
done
echo "--- statistics ---"
echo "C++ :"; grep -E "^(substitution|insertion|deletion|trimmed)" "$W/../scale_diag.log" 2>/dev/null | head -5
echo "Rust:"; grep -E "^(substitution|insertion|deletion|trimmed)" "$W/r36/log"

echo
echo "=== where the single-threaded run spends its time ==="
if command -v perf >/dev/null 2>&1; then
  perf record -q -g --call-graph=dwarf -F 199 -o "$W/prof/perf.data" -- \
     "$RS" correction --reads "$W/big.fq" --counts "$W/big.counts.h5" \
     --outdir "$W/prof" --outprefix t --threads 1 >/dev/null 2>&1
  perf report -i "$W/prof/perf.data" --no-children --percent-limit 1 --stdio 2>/dev/null | head -40
else
  echo "perf not available; timing the stages separately instead"
  echo "-- table build only (counts2flat writes after building, so subtract the write) --"
  /usr/bin/time -f "  counts2flat total %e s  %M KB" "$RS" counts2flat --counts "$W/big.counts.h5" \
     --outdir "$W/prof" --outprefix t 2>&1 | tail -1
  echo "-- profile command: fill + write, no correction --"
  /usr/bin/time -f "  profile  %e s" "$RS" profile --reads "$W/big.fq" --counts "$W/big.counts.h5" \
     --outdir "$W/prof" --outprefix t --threads 1 2>&1 | tail -1
  echo "-- abundance: fill + quantile, minimal output --"
  /usr/bin/time -f "  abundance %e s" "$RS" abundance --reads "$W/big.fq" --counts "$W/big.counts.h5" \
     --outdir "$W/prof" --outprefix t --threads 1 2>&1 | tail -1
  echo "-- correction: fill + all correction passes --"
  /usr/bin/time -f "  correction %e s" "$RS" correction --reads "$W/big.fq" --counts "$W/big.counts.h5" \
     --outdir "$W/prof" --outprefix t --threads 1 2>&1 | tail -1
fi
