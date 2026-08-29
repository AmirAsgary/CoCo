//! Packed k-mer primitives, ported from `src/kmer.cpp`.

use crate::types::INT2RES;

/// A k-mer over the full spaced *span*, two bits per position (`src/kmer.h`).
pub type SpacedKmer = u128;
/// A k-mer over only the *informative* positions, two bits each.
pub type PackedKmer = u64;

/// Reverse complement of a 2-bit packed k-mer.
///
/// The C++ has two implementations that must agree: a scalar loop and an SSSE3
/// `pshufb` version. This is a third, branch-free formulation that is faster than
/// both and provably identical: complementing is `^ 0b10` per base in CoCo's
/// `A=0,C=1,T=2,G=3` encoding, so the whole operation is "reverse the 2-bit groups,
/// then flip bit 1 of each". The group reversal is the standard swap-cascade
/// followed by `swap_bytes`.
#[inline(always)]
pub fn rev_complement(kmer: PackedKmer, kmer_size: u16) -> PackedKmer {
    debug_assert!(kmer_size >= 1 && kmer_size <= 32);
    let mut x = kmer;
    // Reverse the order of the 32 2-bit groups.
    x = ((x & 0x3333_3333_3333_3333) << 2) | ((x >> 2) & 0x3333_3333_3333_3333);
    x = ((x & 0x0F0F_0F0F_0F0F_0F0F) << 4) | ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F);
    x = x.swap_bytes();
    // Drop the groups beyond k, which the reversal pushed into the low bits.
    x >>= 64 - 2 * kmer_size as u32;
    // Complement every remaining base.
    let used = if kmer_size == 32 { u64::MAX } else { (1u64 << (2 * kmer_size)) - 1 };
    x ^ (0xAAAA_AAAA_AAAA_AAAA & used)
}

/// Scalar reference implementation, kept only so tests can pin `rev_complement`
/// against the exact loop the C++ falls back to when SSSE3 is unavailable.
#[cfg(test)]
pub fn rev_complement_scalar(kmer: PackedKmer, kmer_size: u16) -> PackedKmer {
    use crate::types::INT2REV;
    let mut rev: u64 = 0;
    let mut cp = kmer;
    for _ in 0..kmer_size {
        let nuc = (cp & 3) as usize;
        rev <<= 2;
        rev += INT2REV[nuc] as u64;
        cp >>= 2;
    }
    rev
}

/// Canonical form: the smaller of the k-mer and its reverse complement.
#[inline(always)]
pub fn min_index(kmer: PackedKmer, kmer_size: u16) -> PackedKmer {
    let rc = rev_complement(kmer, kmer_size);
    if kmer < rc {
        kmer
    } else {
        rc
    }
}

/// Decode a packed k-mer back to nucleotides, most significant base first.
pub fn packed_kmer_to_string(kmer: PackedKmer, kmer_size: u16) -> String {
    let n = kmer_size as usize;
    let mut out = vec![0u8; n];
    let mut k = kmer;
    for idx in 0..n {
        out[n - 1 - idx] = INT2RES[(k & 3) as usize];
        k >>= 2;
    }
    // Safe: every byte written comes from INT2RES, which is ASCII.
    unsafe { String::from_utf8_unchecked(out) }
}

/// Encode a nucleotide string into a packed k-mer. Inverse of
/// [`packed_kmer_to_string`]; returns `None` on a non-ACGT character.
pub fn string_to_packed_kmer(s: &[u8]) -> Option<PackedKmer> {
    let mut k: u64 = 0;
    for &c in s {
        let v = crate::types::res2int(c);
        if v < 0 {
            return None;
        }
        k = (k << 2) | v as u64;
    }
    Some(k)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rev_complement_matches_scalar_reference() {
        // Deterministic xorshift so the test is reproducible.
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for k in 1..=32u16 {
            let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
            for _ in 0..2000 {
                let v = next() & mask;
                assert_eq!(
                    rev_complement(v, k),
                    rev_complement_scalar(v, k),
                    "k={k} v={v:#x}"
                );
            }
            // Edge values as well.
            for &v in &[0u64, mask, 1, mask ^ 1] {
                assert_eq!(rev_complement(v, k), rev_complement_scalar(v, k), "k={k} v={v:#x}");
            }
        }
    }

    #[test]
    fn rev_complement_is_an_involution() {
        for k in 1..=32u16 {
            let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
            for v in [0u64, 1, 12345 & mask, mask] {
                assert_eq!(rev_complement(rev_complement(v, k), k), v);
            }
        }
    }

    #[test]
    fn rev_complement_known_values() {
        // "ACTG" -> codes 0,1,2,3 -> 0b00_01_10_11 = 0x1B
        let k = string_to_packed_kmer(b"ACTG").unwrap();
        assert_eq!(k, 0b00_01_10_11);
        // revcomp("ACTG") = "CAGT"
        assert_eq!(packed_kmer_to_string(rev_complement(k, 4), 4), "CAGT");
    }

    #[test]
    fn min_index_picks_the_smaller() {
        for k in 1..=32u16 {
            let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
            for v in [0u64, 1, 0x1234_5678 & mask, mask] {
                let m = min_index(v, k);
                assert!(m == v || m == rev_complement(v, k));
                assert!(m <= v && m <= rev_complement(v, k));
            }
        }
    }

    #[test]
    fn string_round_trip() {
        for s in [&b"AAAA"[..], b"ACTG", b"GTCA", b"TTTTTTTTTTTTTTTT"] {
            let k = string_to_packed_kmer(s).unwrap();
            assert_eq!(packed_kmer_to_string(k, s.len() as u16).as_bytes(), s);
        }
    }
}
