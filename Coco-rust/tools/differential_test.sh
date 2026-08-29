#!/bin/bash
# Differential test: run the C++ and Rust binaries over many generated datasets and
# option combinations, and require byte-identical output every time.
#
# The committed fixtures cover the substitution path only. These datasets have
# ragged metagenome-like coverage plus real indels, N's, soft-masking, varied read
# lengths and varied quality, which is what reaches doIndelCorrection, doTrimming
# and the invalid-nucleotide branches.
set -uo pipefail
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
W=${W:-$ROOT/work/diff}
THREADS=${THREADS:-16}

mkdir -p "$W"
pass=0; fail=0

diffcheck () { # label f1 f2
  if cmp -s "$2" "$3"; then
    pass=$((pass+1))
  else
    fail=$((fail+1))
    echo "  FAIL  $1"
    if [ -f "$2" ] && [ -f "$3" ]; then
      echo "        first difference:"
      diff <(cat "$2") <(cat "$3") | head -6 | sed 's/^/          /'
    else
      echo "        missing output: $2 / $3"
    fi
  fi
}

# ---------------------------------------------------------------- datasets ----
# name : generator flags
declare -A DATASETS=(
  [indels]="--seed 1 --contigs 40 --sub-rate 0.004 --ins-rate 0.002 --del-rate 0.002"
  [heavy_indels]="--seed 2 --contigs 60 --min-depth 2 --max-depth 200 --sub-rate 0.008 --ins-rate 0.005 --del-rate 0.005"
  [with_n]="--seed 3 --contigs 30 --n-rate 0.004 --sub-rate 0.004 --ins-rate 0.002 --del-rate 0.002"
  [softmasked]="--seed 4 --contigs 30 --lowercase-rate 0.15 --sub-rate 0.004 --ins-rate 0.002 --del-rate 0.002"
  [varied_len]="--seed 5 --contigs 30 --read-len 160 --read-len-jitter 55 --sub-rate 0.005 --ins-rate 0.002 --del-rate 0.002"
  [uniform_qual]="--seed 6 --contigs 25 --qual uniform --sub-rate 0.006 --ins-rate 0.003 --del-rate 0.003"
  [comments]="--seed 7 --contigs 25 --comments --sub-rate 0.004 --ins-rate 0.002 --del-rate 0.002"
  [lowcov]="--seed 8 --contigs 80 --min-depth 1.5 --max-depth 12 --sub-rate 0.004 --ins-rate 0.002 --del-rate 0.002"
)

# ------------------------------------------------------------- option sets ----
declare -A OPTSETS=(
  [default]=""
  [trim]="--max-trim-len 5"
  [loose]="--threshold 0.05 --pseudocount 3 --lowerbound 2"
  [tight]="--threshold 0.002 --pseudocount 0 --lowerbound 20"
  [nocap]="--max-corr-num 0"
  [cap1]="--max-corr-num 1"
  [maxmode]="--count-mode 1"
  [noskip]="--skip 0"
  [bigskip]="--skip 60"
  # An explicit value, because a bare boolean flag as the *last* argument
  # dereferences argv[argc] (NULL) in the C++ and segfaults. See README.
  [updatelookup]="--update-lookup true"
  [trim_loose]="--max-trim-len 20 --threshold 0.05 --lowerbound 2"
)

echo "=== generating datasets ==="
for name in "${!DATASETS[@]}"; do
  # shellcheck disable=SC2086
  $PY "$TOOLS/make_stress_data.py" --genome "$TD/genome.fa" \
      --out-prefix "$W/$name" ${DATASETS[$name]} 2>>"$W/gen.log"
done
wc -l "$W"/*.fq | tail -1

echo "=== counting k-mers (dsk, k=41) ==="
for name in "${!DATASETS[@]}"; do
  dsk -file "$W/$name.fq" -kmer-size 41 -abundance-min 2 -out "$W/$name.counts" \
      -out-dir "$W" -out-tmp "${TMPDIR:-/tmp}" -max-memory 8000 -nb-cores 8 \
      >>"$W/dsk.log" 2>&1 || echo "  dsk failed for $name"
done

echo
echo "=== correction: C++ vs Rust over ${#DATASETS[@]} datasets x ${#OPTSETS[@]} option sets ==="
for name in "${!DATASETS[@]}"; do
  counts="$W/$name.counts.h5"
  [ -f "$counts" ] || { echo "  SKIP $name (no counts file)"; continue; }
  for oname in "${!OPTSETS[@]}"; do
    o=${OPTSETS[$oname]}
    dc="$W/out/${name}_${oname}_cpp"; dr="$W/out/${name}_${oname}_rs"
    mkdir -p "$dc" "$dr"
    # shellcheck disable=SC2086
    $CPP correction --reads "$W/$name.fq" --counts "$counts" --outdir "$dc" --outprefix t $o >"$dc/log" 2>&1
    # shellcheck disable=SC2086
    $RS  correction --reads "$W/$name.fq" --counts "$counts" --outdir "$dr" --outprefix t --threads $THREADS $o >"$dr/log" 2>&1
    diffcheck "correction $name/$oname reads" "$dc/t.corr.reads.fq" "$dr/t.corr.reads.fq"
    grep -E "^(substitution|insertion|deletion|trimmed)" "$dc/log" > "$dc/stats"
    grep -E "^(substitution|insertion|deletion|trimmed)" "$dr/log" > "$dr/stats"
    diffcheck "correction $name/$oname stats" "$dc/stats" "$dr/stats"
  done
done

echo
echo "=== other commands (default options) ==="
for name in "${!DATASETS[@]}"; do
  counts="$W/$name.counts.h5"
  [ -f "$counts" ] || continue
  for cmd in profile abundance filter; do
    dc="$W/out/${name}_${cmd}_cpp"; dr="$W/out/${name}_${cmd}_rs"
    mkdir -p "$dc" "$dr"
    $CPP $cmd --reads "$W/$name.fq" --counts "$counts" --outdir "$dc" --outprefix t >"$dc/log" 2>&1
    $RS  $cmd --reads "$W/$name.fq" --counts "$counts" --outdir "$dr" --outprefix t --threads $THREADS >"$dr/log" 2>&1
    f=$(ls "$dc" | grep -v '^log$' | head -1)
    diffcheck "$cmd $name" "$dc/$f" "$dr/$f"
  done
  # filter's --threshold means something different from correction's and has its
  # own default, so sweep it explicitly.
  for thr in 0.05 0.1 0.3 0.6; do
    dc="$W/out/${name}_filt${thr}_cpp"; dr="$W/out/${name}_filt${thr}_rs"
    mkdir -p "$dc" "$dr"
    $CPP filter --reads "$W/$name.fq" --counts "$counts" --threshold $thr --outdir "$dc" --outprefix t >"$dc/log" 2>&1
    $RS  filter --reads "$W/$name.fq" --counts "$counts" --threshold $thr --outdir "$dr" --outprefix t --threads $THREADS >"$dr/log" 2>&1
    diffcheck "filter $name thr=$thr" "$dc/t.filter.reads.fq" "$dr/t.filter.reads.fq"
  done
done

echo
echo "=== compressed input (.gz, .bz2) ==="
first=$(ls "$W"/*.fq 2>/dev/null | head -1)
if [ -n "$first" ]; then
  base=$(basename "$first" .fq)
  counts="$W/$base.counts.h5"
  [ -f "$W/$base.fq.gz" ]  || gzip  -kf "$first"
  [ -f "$W/$base.fq.bz2" ] || bzip2 -kf "$first"
  for z in gz bz2; do
    dc="$W/out/z${z}_cpp"; dr="$W/out/z${z}_rs"
    mkdir -p "$dc" "$dr"
    $CPP correction --reads "$W/$base.fq.$z" --counts "$counts" --outdir "$dc" --outprefix t >"$dc/log" 2>&1
    $RS  correction --reads "$W/$base.fq.$z" --counts "$counts" --outdir "$dr" --outprefix t --threads $THREADS >"$dr/log" 2>&1
    diffcheck "correction from .$z" "$dc/t.corr.reads.fq" "$dr/t.corr.reads.fq"
  done
fi

echo
echo "=== paired-end ==="
$PY "$TOOLS/make_stress_data.py" --genome "$TD/genome.fa" --out-prefix "$W/pe" \
    --seed 11 --contigs 30 --paired --sub-rate 0.004 --ins-rate 0.002 --del-rate 0.002 2>>"$W/gen.log"
cat "$W/pe.1.fq" "$W/pe.2.fq" > "$W/pe.both.fq"
[ -f "$W/pe.counts.h5" ] || dsk -file "$W/pe.both.fq" -kmer-size 41 -abundance-min 2 -out "$W/pe.counts" \
    -out-dir "$W" -out-tmp "${TMPDIR:-/tmp}" -max-memory 8000 -nb-cores 8 >>"$W/dsk.log" 2>&1
if [ -f "$W/pe.counts.h5" ]; then
  for oname in default trim nocap; do
    o=${OPTSETS[$oname]}
    dc="$W/out/pe_${oname}_cpp"; dr="$W/out/pe_${oname}_rs"
    mkdir -p "$dc" "$dr"
    # shellcheck disable=SC2086
    $CPP correction -1 "$W/pe.1.fq" -2 "$W/pe.2.fq" --counts "$W/pe.counts.h5" --outdir "$dc" --outprefix t $o >"$dc/log" 2>&1
    # shellcheck disable=SC2086
    $RS  correction -1 "$W/pe.1.fq" -2 "$W/pe.2.fq" --counts "$W/pe.counts.h5" --outdir "$dr" --outprefix t --threads $THREADS $o >"$dr/log" 2>&1
    diffcheck "paired correction $oname r1" "$dc/t.corr.1.fq" "$dr/t.corr.1.fq"
    diffcheck "paired correction $oname r2" "$dc/t.corr.2.fq" "$dr/t.corr.2.fq"
  done
  dc="$W/out/pe_filter_cpp"; dr="$W/out/pe_filter_rs"; mkdir -p "$dc" "$dr"
  $CPP filter -1 "$W/pe.1.fq" -2 "$W/pe.2.fq" --counts "$W/pe.counts.h5" --outdir "$dc" --outprefix t >"$dc/log" 2>&1
  $RS  filter -1 "$W/pe.1.fq" -2 "$W/pe.2.fq" --counts "$W/pe.counts.h5" --outdir "$dr" --outprefix t --threads $THREADS >"$dr/log" 2>&1
  diffcheck "paired filter r1" "$dc/t.filter.1.fq" "$dr/t.filter.1.fq"
  diffcheck "paired filter r2" "$dc/t.filter.2.fq" "$dr/t.filter.2.fq"
fi

echo
echo "=== thread-count invariance (Rust only) ==="
first=$(ls "$W"/*.fq 2>/dev/null | grep -v both | head -1)
base=$(basename "$first" .fq)
for t in 1 3 8 32; do
  d="$W/out/thr${t}"; mkdir -p "$d"
  $RS correction --reads "$first" --counts "$W/$base.counts.h5" --outdir "$d" --outprefix t --threads $t >"$d/log" 2>&1
done
for t in 3 8 32; do
  diffcheck "threads=$t matches threads=1" "$W/out/thr1/t.corr.reads.fq" "$W/out/thr${t}/t.corr.reads.fq"
done

echo
echo "=== indel coverage check ==="
tot_ins=0; tot_del=0; tot_trim=0
for f in "$W"/out/*_cpp/log; do
  i=$(grep -oP '(?<=^insertion corrections: )\d+' "$f" 2>/dev/null | head -1); i=${i:-0}
  d=$(grep -oP '(?<=^deletion corrections: )\d+' "$f" 2>/dev/null | head -1); d=${d:-0}
  t=$(grep -oP '(?<=^trimmed nucleotides: )\d+' "$f" 2>/dev/null | head -1); t=${t:-0}
  tot_ins=$((tot_ins+i)); tot_del=$((tot_del+d)); tot_trim=$((tot_trim+t))
done
echo "  insertion corrections exercised: $tot_ins"
echo "  deletion corrections exercised:  $tot_del"
echo "  trimmed nucleotides exercised:   $tot_trim"

echo
echo "=== RESULT: $pass passed, $fail failed ==="
[ "$fail" -eq 0 ]
