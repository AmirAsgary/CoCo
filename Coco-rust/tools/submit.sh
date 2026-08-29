#!/bin/bash
# Submit a tools/ script to a SLURM GPU node, running a *snapshot* of tools/.
#
# Bash reads a script incrementally, so editing one while it runs shifts the byte
# offsets under the interpreter and produces a spurious syntax error partway
# through. Copying the directory first makes long runs immune to edits made while
# they are in flight. COCO_RUST_ROOT tells the snapshot where the real repository
# is, since $0 no longer points at it.
#
#   usage: tools/submit.sh <job-name> <cpus> <minutes> <script.sh> [args...]
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NAME=$1 CPUS=$2 MINUTES=$3 SCRIPT=$4; shift 4

# Raven GPU nodes have 72 cores and 4 A100s, and SLURM wants roughly one GPU per
# 18 cores requested.
GPUS=$(( (CPUS + 17) / 18 )); [ "$GPUS" -lt 1 ] && GPUS=1; [ "$GPUS" -gt 4 ] && GPUS=4

SNAP=$ROOT/work/snapshots/$NAME.$(date +%Y%m%d-%H%M%S)
mkdir -p "$SNAP" "$ROOT/work"
cp -r "$ROOT/tools" "$SNAP/"

sbatch --partition="${COCO_PARTITION:-gpu1}" --gres=gpu:a100:"$GPUS" \
       --cpus-per-task="$CPUS" --mem=$((GPUS*120))G --time="$MINUTES" \
       --job-name="$NAME" -o "$ROOT/work/$NAME.log" \
       --wrap="COCO_RUST_ROOT=$ROOT bash $SNAP/tools/$(basename "$SCRIPT") $*"
