#!/bin/bash
set -euo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
cd $ROOT/oracle
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE .. > cmake.log 2>&1 || { tail -30 cmake.log; exit 1; }
make -j "${SLURM_CPUS_PER_TASK:-8}" > make.log 2>&1 || { tail -40 make.log; exit 1; }
ls -la oracle
