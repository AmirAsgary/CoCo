#!/bin/bash
# shellcheck source=tools/_common.sh
. "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
# Prepare everything needed to build the original C++ CoCo for comparison.
#
# There is no copy of the C++ sources in Coco-rust: this directory lives inside
# the CoCo repository, so `../src` *is* the reference. Only two things have to be
# created, and both are generated rather than committed:
#
#   1. lib/gatb-core -- the submodule is empty in a fresh checkout, which is why
#      the original does not build as-is.
#   2. work/reference_fixed -- the same tree with one function corrected, so the
#      internal-k-mer-counting path has a specification to compare against. See
#      MODIFICATIONS.md.
#
# Needs network access, so run it on a login node, not a compute node.
set -euo pipefail
WORK=$ROOT/work
mkdir -p "$WORK"

echo "== gatb-core =="
if [ -f "$COCO/lib/gatb-core/gatb-core/CMakeLists.txt" ]; then
  echo "   already present"
else
  git -C "$COCO" submodule update --init lib/gatb-core \
    || git clone --depth 1 --branch v1.4.2 https://github.com/GATB/gatb-core.git \
                 "$COCO/lib/gatb-core"
  echo "   ok"
fi

echo "== work/reference_fixed =="
FIX=$WORK/reference_fixed
rm -rf "$FIX"; mkdir -p "$FIX"
cp -r "$COCO/src" "$COCO/main.cpp" "$COCO/CMakeLists.txt" "$COCO/cmake" "$FIX/"
ln -s "$COCO/lib" "$FIX/lib"

python3 - "$FIX/src/HashTable.h" <<'PY'
import sys
p = sys.argv[1]
s = open(p).read()
old = """  unsigned int getCount(const packedKmerType kmer) const {
    khint_t itr;
    int absent;
    itr = kc_c1_put(this->hashTable, kmer, &absent);
    if (absent) return 1;
    return kh_val(this->hashTable, itr);
  }"""
new = """  unsigned int getCount(const packedKmerType kmer) const {
    /* PATCHED FOR COMPARISON WITH THE RUST PORT -- NOT THE ORIGINAL.
     * The original calls kc_c1_put here, which *inserts* the k-mer with an
     * uninitialised value taken from an uninitialised stack bucket struct. The
     * first lookup of an absent k-mer correctly returns 1, but every later lookup
     * of that same k-mer returns whatever was on the stack. This build uses the
     * non-inserting get so that "absent means 1" holds consistently, which is
     * what the code plainly intends and what the port implements. */
    khint_t itr = kc_c1_get(this->hashTable, kmer);
    if (itr == kh_end(this->hashTable)) return 1;
    return kh_val(this->hashTable, itr);
  }"""
if old not in s:
    sys.exit("HashTable.h does not contain the expected getCount; patch not applied")
open(p, "w").write(s.replace(old, new))
print("   patched HashTable::getCount")
PY

echo
echo "Now build them on a compute node:"
echo "  tools/submit.sh coco-ref    18 120 build_reference.sh"
echo "  tools/submit.sh coco-oracle 18  90 build_oracle.sh"
