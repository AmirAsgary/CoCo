#!/bin/bash
# Build, then run every verification layer in order. Intended for tools/submit.sh.
set -uo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
export COCO_TESTDATA=$TD COCO_WORK=$ROOT/work
rc=0
banner () { echo; echo "=============== $1 ==============="; }

banner "build"
# Capture to a file first: piping cargo into grep hides its exit status behind
# pipefail's rightmost-failure rule, which silently let a broken build through.
if ! ( cd "$ROOT/coco" && cargo build --release ) > "$ROOT/work/build.log" 2>&1; then
  echo "BUILD FAILED"
  grep -E "^error" -A8 "$ROOT/work/build.log" | head -40
  exit 1
fi
grep -E "^warning: (unused|value|unreachable)" -A4 "$ROOT/work/build.log" | head -30
echo "build ok"

banner "unit tests"
( cd "$ROOT/coco" && cargo test --release --lib 2>&1 | grep -E "test result|^error|FAILED" ) || rc=1

banner "golden fixtures"
( cd "$ROOT/coco" && cargo test --release --test golden_kmer2packed --test golden_dsk \
    --test golden_lookuptable --test golden_correction -- --test-threads=2 2>&1 \
    | grep -E "test result|^error|FAILED|panicked" ) || rc=1

banner "oracle: per-function vs the C++"
( cd "$ROOT/coco" && cargo test --release --test oracle -- --test-threads=4 2>&1 \
    | grep -E "test result|^error|FAILED|panicked" ) || rc=1

banner "binary comparison, committed test set"
bash "$(dirname "$0")/compare_all.sh" || rc=1

banner "RESULT"
[ $rc -eq 0 ] && echo "ALL FAST SUITES PASSED" || echo "FAILURES ABOVE"
exit $rc
