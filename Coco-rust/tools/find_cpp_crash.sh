#!/bin/bash
# Bisect a read file down to the smallest prefix that still crashes the C++.
set -uo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
W=$ROOT/work/crash
READS=${READS:-$ROOT/work/diff/uniform_qual.fq}
COUNTS=${COUNTS:-$ROOT/work/diff/uniform_qual.counts.h5}
OPTS=${OPTS:---update-lookup}
rm -rf "$W"; mkdir -p "$W"

crashes () { # nrecords
  head -n $(( $1 * 4 )) "$READS" > "$W/slice.fq"
  # shellcheck disable=SC2086
  $CPP correction --reads "$W/slice.fq" --counts "$COUNTS" --outdir "$W/o" --outprefix t $OPTS \
      >"$W/o.log" 2>&1
  local rc=$?
  [ $rc -ge 128 ] || [ $rc -eq 139 ]
}

total=$(( $(wc -l < "$READS") / 4 ))
echo "reads: $total"
if ! crashes "$total"; then echo "no crash on the full file"; exit 0; fi

lo=1; hi=$total
while [ $lo -lt $hi ]; do
  mid=$(( (lo + hi) / 2 ))
  if crashes "$mid"; then hi=$mid; else lo=$(( mid + 1 )); fi
done
echo "smallest crashing prefix: $lo reads"
head -n $(( lo * 4 )) "$READS" > "$W/minimal.fq"
tail -4 "$W/minimal.fq" | head -2
echo
echo "--- C++ on that prefix ---"
$CPP correction --reads "$W/minimal.fq" --counts "$COUNTS" --outdir "$W/oc" --outprefix t $OPTS 2>&1 | tail -6
echo "exit=$?"
echo
echo "--- Rust on the same prefix ---"
$RS correction --reads "$W/minimal.fq" --counts "$COUNTS" --outdir "$W/or" --outprefix t --threads 1 $OPTS 2>&1 | tail -8
echo "exit=$?"
echo
echo "--- length of the last read, and the span ---"
awk 'NR%4==2{print length($0)}' "$W/minimal.fq" | tail -3
