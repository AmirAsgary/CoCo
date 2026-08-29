#!/bin/bash
# A/B the sharded parallel table build against the single-table build.
# Three reps each, best of, so the 5% differences seen in single-sample runs can
# be told apart from noise.
set -uo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
# Built into its own target directory so this does not swap the binary out from
# under a differential run that may be in flight.
export CARGO_TARGET_DIR=$ROOT/coco/target-ab
( cd "$ROOT/coco" && cargo build --release ) >/dev/null 2>&1 || { echo "build failed"; exit 1; }
W=$ROOT/work/scale
TMP=$(mktemp); mkdir -p "$W/ab"

run () { # label env threads cmd...
  local label=$1 tbl=$2 t=$3 cmd=$4
  local best=""
  for _ in 1 2 3; do
    COCO_TABLE=$tbl /usr/bin/time -f "%e" -o "$TMP" \
      "$RS" "$cmd" --reads "$W/big.fq" --counts "$W/big.counts.h5" \
      --outdir "$W/ab" --outprefix t --threads "$t" >/dev/null 2>&1
    local sec; sec=$(cat "$TMP")
    [ -z "$best" ] && best=$sec
    best=$(awk -v a="$sec" -v b="$best" 'BEGIN{print (a<b)?a:b}')
  done
  printf "%-34s %8ss\n" "$label" "$best"
}

echo "### correction (fill + correct + 924 MB of output) ###"
for t in 1 8 36; do
  run "single table, $t threads"  single   "$t" correction
  run "sharded table, $t threads" sharded  "$t" correction
done
echo
echo "### abundance (fill + quantile, ~60 MB of output) -- isolates output cost ###"
for t in 1 36; do
  run "single table, $t threads"  single  "$t" abundance
  run "sharded table, $t threads" sharded "$t" abundance
done
echo
echo "### how much of the run is just I/O? ###"
printf "  input  %s\n" "$(du -h "$W/big.fq" | cut -f1)"
printf "  correction output %s\n" "$(du -h "$W/ab/t.corr.reads.fq" 2>/dev/null | cut -f1)"
echo "  cat of the input to /dev/null:"
/usr/bin/time -f "    %e s" cat "$W/big.fq" > /dev/null
echo "  copy of the input (read + write):"
/usr/bin/time -f "    %e s" cp "$W/big.fq" "$W/ab/copy.fq"
rm -f "$W/ab/copy.fq" "$TMP"
