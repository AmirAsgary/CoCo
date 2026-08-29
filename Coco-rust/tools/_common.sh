# Shared path resolution for the scripts in this directory. Source it, do not run it.
#
# tools/submit.sh runs a *snapshot* of tools/ so that editing a script cannot
# corrupt a job already reading it, which means $0's directory is not the
# repository. It exports COCO_RUST_ROOT for that case; direct invocation falls
# back to resolving the path from the sourcing script.
ROOT=${COCO_RUST_ROOT:-$(cd "$(dirname "${BASH_SOURCE[1]}")/.." && pwd)}
COCO=$(cd "$ROOT/.." && pwd)

RS=$ROOT/coco/target/release/coco                    # the Rust binary
CPP=$ROOT/work/reference-build/coco                  # the original C++
FIXED=$ROOT/work/reference_fixed-build/coco          # C++ with getCount corrected
ORACLE=$ROOT/oracle/build/oracle                     # per-function harness

TD=${COCO_TESTDATA:-$COCO/testdata}
WORK=${COCO_WORK:-$ROOT/work}
DSKBIN=$ROOT/.conda/dsk/bin
PY=$DSKBIN/python

export CARGO_HOME=$ROOT/.rust/cargo RUSTUP_HOME=$ROOT/.rust/rustup
export PATH=$ROOT/.rust/cargo/bin:$DSKBIN:$PATH

TOOLS=$(cd "$(dirname "${BASH_SOURCE[1]}")" && pwd)  # this script's own directory

need () { # <path> <what to run to create it>
  [ -e "$1" ] || { echo "missing: $1"; echo "create it with: $2"; exit 1; }
}
