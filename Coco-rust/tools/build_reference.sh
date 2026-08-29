#!/bin/bash
# Build the original C++ CoCo, and the one-function-corrected variant, out of
# source into work/. Run tools/setup_reference.sh first.
set -euo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
J=${SLURM_CPUS_PER_TASK:-8}

build () { # <source dir> <build dir> <label>
  echo "== $3 =="
  mkdir -p "$2"
  cmake -S "$1" -B "$2" -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. \
    > "$2/cmake.log" 2>&1 || { tail -30 "$2/cmake.log"; exit 1; }
  make -C "$2" -j "$J" coco > "$2/make.log" 2>&1 || { tail -40 "$2/make.log"; exit 1; }
  ls -la "$2/coco"
}

[ -f "$COCO/lib/gatb-core/gatb-core/CMakeLists.txt" ] \
  || { echo "gatb-core is missing -- run tools/setup_reference.sh first"; exit 1; }

build "$COCO"                       "$ROOT/work/reference-build"       "reference"
[ -d "$ROOT/work/reference_fixed" ] \
  && build "$ROOT/work/reference_fixed" "$ROOT/work/reference_fixed-build" "reference_fixed"
