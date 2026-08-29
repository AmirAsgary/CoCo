//! Spaced-seed k-mer translation, ported from `src/KmerTranslator.cpp`.
//!
//! A *spaced* k-mer covers `span` consecutive bases but only keeps the `weight`
//! positions marked `1` in the pattern. Translation gathers those bit pairs out of
//! the span-wide value and packs them contiguously.

use crate::kmer::{min_index, PackedKmer, SpacedKmer};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranslatorError(pub String);

impl std::fmt::Display for TranslatorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for TranslatorError {}

/// Sentinel for "this span position is not informative", `UCHAR_MAX` in the C++.
pub const NOT_INFORMATIVE: u8 = u8::MAX;

/// CoCo's default spaced pattern: span 41, weight 32.
pub const DEFAULT_PATTERN: &str = "11110111111011011101010111011011111101111";

pub struct KmerTranslator {
    span: u16,
    weight: u16,
    /// One bit per span position, MSB = position 0. `_spaced_mask` in the C++.
    spaced_mask: SpacedKmer,
    /// Ascending list of the `weight` informative span offsets.
    mask_array: Vec<u8>,
    /// span offset -> index into `mask_array`, or [`NOT_INFORMATIVE`].
    inverse_mask_array: Vec<u8>,
    /// `spaced_mask` widened to two bits per position: the gather mask over a
    /// span-wide packed k-mer. Split into halves for 64-bit `PEXT`.
    pext_lo: u64,
    pext_hi: u64,
    /// Number of informative bits living in the low half, i.e. the shift that
    /// re-joins the two `PEXT` results.
    pext_lo_bits: u32,
    /// Whether BMI2 was detected at construction time.
    has_bmi2: bool,
}

impl KmerTranslator {
    pub fn new(pattern: &str) -> Result<Self, TranslatorError> {
        let bytes = pattern.as_bytes();
        let span = bytes.len();
        if span > 64 {
            return Err(TranslatorError("Error: Only kmerSpan <= 64 is supported".into()));
        }
        if span == 0 {
            return Err(TranslatorError(
                "Error: First and last position in spacedKmerPattern must 1".into(),
            ));
        }

        let mut weight: u16 = 0;
        for idx in 0..span {
            match bytes[idx] {
                b'1' => weight += 1,
                b'0' => continue,
                _ => {
                    return Err(TranslatorError(
                        "Error: Found invalid character in spacedKmerPattern".into(),
                    ))
                }
            }
            // Only reached for '1' positions -- the C++ `continue` above skips the
            // check for '0', which is equivalent to requiring full symmetry.
            if bytes[idx] != bytes[span - idx - 1] {
                return Err(TranslatorError(
                    "Error: spacedKmerPattern must be symmetric".into(),
                ));
            }
        }

        if bytes[0] != b'1' || bytes[span - 1] != b'1' {
            return Err(TranslatorError(
                "Error: First and last position in spacedKmerPattern must 1".into(),
            ));
        }
        if !(12..=32).contains(&weight) {
            return Err(TranslatorError(
                "Error: Only 12 <= weight <= 32 is supported".into(),
            ));
        }

        let mut mask_array = Vec::with_capacity(weight as usize);
        let mut inverse_mask_array = vec![NOT_INFORMATIVE; span];
        for idx in 0..span {
            if bytes[idx] == b'1' {
                inverse_mask_array[idx] = mask_array.len() as u8;
                mask_array.push(idx as u8);
            }
        }

        let mut spaced_mask: SpacedKmer = 0;
        for &b in bytes {
            spaced_mask = (spaced_mask << 1) | ((b == b'1') as SpacedKmer);
        }

        // Widen to two bits per position: span position p occupies bits
        // [2*(span-1-p), 2*(span-1-p)+1] of a span-wide packed k-mer.
        let mut gather: u128 = 0;
        for &p in &mask_array {
            gather |= 3u128 << (2 * (span - 1 - p as usize));
        }
        let pext_lo = gather as u64;
        let pext_hi = (gather >> 64) as u64;

        Ok(KmerTranslator {
            span: span as u16,
            weight,
            spaced_mask,
            mask_array,
            inverse_mask_array,
            pext_lo,
            pext_hi,
            pext_lo_bits: pext_lo.count_ones(),
            has_bmi2: bmi2_available(),
        })
    }

    #[inline(always)]
    pub fn span(&self) -> u16 {
        self.span
    }
    #[inline(always)]
    pub fn weight(&self) -> u16 {
        self.weight
    }
    #[inline(always)]
    pub fn spaced_mask(&self) -> SpacedKmer {
        self.spaced_mask
    }
    #[inline(always)]
    pub fn mask_array(&self) -> &[u8] {
        &self.mask_array
    }
    #[inline(always)]
    pub fn inverse_mask_array(&self) -> &[u8] {
        &self.inverse_mask_array
    }
    #[inline(always)]
    pub fn has_bmi2(&self) -> bool {
        self.has_bmi2
    }

    /// Longest run of consecutive informative positions.
    pub fn longest_block(&self) -> u16 {
        let mut cur = 0u16;
        let mut max = 0u16;
        for idx in 0..self.span as usize {
            if self.inverse_mask_array[idx] != NOT_INFORMATIVE {
                cur += 1;
            } else {
                cur = 0;
            }
            if cur > max {
                max = cur;
            }
        }
        max
    }

    /// How the packed k-mer is split into grid index and in-grid offset.
    pub fn best_split(&self) -> (u32, u32) {
        let w = self.weight as u32;
        if 2 * w > 30 {
            (30, 2 * w - 30)
        } else {
            (2 * w - 2, 2)
        }
    }

    /// Gather the informative bit pairs out of a span-wide k-mer.
    ///
    /// Portable formulation, bit-identical to the C++ loop.
    #[inline(always)]
    pub fn kmer2packed_scalar(&self, kmer: SpacedKmer) -> PackedKmer {
        let mut packed: u64 = 0;
        let span = self.span as usize;
        for &p in &self.mask_array {
            let shift = 2 * (span - 1 - p as usize);
            packed = (packed << 2) | ((kmer >> shift) & 3) as u64;
        }
        packed
    }

    /// `PEXT`-based gather. Two `PEXT`s replace the `weight`-iteration loop.
    ///
    /// `PEXT` compacts the selected bits toward the LSB while preserving their
    /// relative order. Span position `p` sits at bit `2*(span-1-p)`, so ascending
    /// `p` means descending bit index -- the same direction as the scalar loop's
    /// shift-left accumulation. The orders therefore agree and the results are
    /// bit-identical, which `pext_matches_scalar` pins down.
    ///
    /// # Safety
    /// Caller must have verified BMI2 support.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "bmi2")]
    #[inline]
    unsafe fn kmer2packed_pext(&self, kmer: SpacedKmer) -> PackedKmer {
        use std::arch::x86_64::_pext_u64;
        let lo = _pext_u64(kmer as u64, self.pext_lo);
        let hi = _pext_u64((kmer >> 64) as u64, self.pext_hi);
        (hi << self.pext_lo_bits) | lo
    }

    #[inline(always)]
    pub fn kmer2packed(&self, kmer: SpacedKmer) -> PackedKmer {
        #[cfg(target_arch = "x86_64")]
        {
            if self.has_bmi2 {
                // Safety: `has_bmi2` was set by runtime feature detection.
                return unsafe { self.kmer2packed_pext(kmer) };
            }
        }
        self.kmer2packed_scalar(kmer)
    }

    #[inline(always)]
    pub fn kmer2min_packed(&self, kmer: SpacedKmer) -> PackedKmer {
        min_index(self.kmer2packed(kmer), self.weight)
    }

    /// Canonicalise an already-packed k-mer (the C++ overload on `packedKmerType`).
    #[inline(always)]
    pub fn packed2min_packed(&self, kmer: PackedKmer) -> PackedKmer {
        min_index(kmer, self.weight)
    }
}

#[inline]
fn bmi2_available() -> bool {
    #[cfg(target_arch = "x86_64")]
    {
        // Zen 1/2 implement PEXT in microcode at ~18 cycles, far slower than the
        // scalar loop, so only take the fast path where PEXT is a single uop.
        std::is_x86_feature_detected!("bmi2") && !is_slow_pext()
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        false
    }
}

#[cfg(target_arch = "x86_64")]
fn is_slow_pext() -> bool {
    // AMD family 17h/19h before Zen 3 emulate PEXT. Distinguishing them needs
    // CPUID leaf 1; anything non-AMD, and Zen 3+, is fast.
    let vendor = {
        let r = std::arch::x86_64::__cpuid(0);
        let mut v = [0u8; 12];
        v[0..4].copy_from_slice(&r.ebx.to_le_bytes());
        v[4..8].copy_from_slice(&r.edx.to_le_bytes());
        v[8..12].copy_from_slice(&r.ecx.to_le_bytes());
        v
    };
    if &vendor != b"AuthenticAMD" {
        return false;
    }
    let (family, _model) = {
        let r = std::arch::x86_64::__cpuid(1);
        let base_family = (r.eax >> 8) & 0xF;
        let ext_family = (r.eax >> 20) & 0xFF;
        let family = if base_family == 0xF { base_family + ext_family } else { base_family };
        let base_model = (r.eax >> 4) & 0xF;
        let ext_model = (r.eax >> 16) & 0xF;
        (family, (ext_model << 4) | base_model)
    };
    // Zen 1/2 are family 0x17; Zen 3 (0x19) onwards has a hardware PEXT.
    family == 0x17
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_pattern_geometry() {
        let t = KmerTranslator::new(DEFAULT_PATTERN).unwrap();
        assert_eq!(t.span(), 41);
        assert_eq!(t.weight(), 32);
        assert_eq!(t.best_split(), (30, 34));
        assert_eq!(t.longest_block(), 6);
        assert_eq!(t.mask_array().len(), 32);
        // Pattern starts "1111 0 111111 0 ..." -> first four offsets are 0..3.
        assert_eq!(&t.mask_array()[..4], &[0, 1, 2, 3]);
        assert_eq!(t.inverse_mask_array()[4], NOT_INFORMATIVE);
        assert_eq!(t.inverse_mask_array()[3], 3);
    }

    #[test]
    fn spaced_mask_is_the_pattern_as_binary() {
        let t = KmerTranslator::new(DEFAULT_PATTERN).unwrap();
        let expect = u128::from_str_radix(DEFAULT_PATTERN, 2).unwrap();
        assert_eq!(t.spaced_mask(), expect);
    }

    #[test]
    fn rejects_invalid_patterns() {
        // Non-binary character.
        assert!(KmerTranslator::new("1112111").is_err());
        // Asymmetric.
        assert!(KmerTranslator::new("1101011").is_err());
        // Does not start/end with 1.
        assert!(KmerTranslator::new("0111110").is_err());
        // Weight below 12.
        assert!(KmerTranslator::new("1111111111").is_err());
        // Span above 64.
        let long = "1".repeat(65);
        assert!(KmerTranslator::new(&long).is_err());
        // Weight above 32.
        let heavy = "1".repeat(33);
        assert!(KmerTranslator::new(&heavy).is_err());
    }

    #[test]
    fn accepts_boundary_patterns() {
        assert!(KmerTranslator::new(&"1".repeat(12)).is_ok());
        assert!(KmerTranslator::new(&"1".repeat(32)).is_ok());
        let t = KmerTranslator::new(&"1".repeat(12)).unwrap();
        assert_eq!(t.best_split(), (22, 2));
    }

    #[test]
    fn pext_matches_scalar() {
        let all32 = "1".repeat(32);
        let all12 = "1".repeat(12);
        let patterns = [
            DEFAULT_PATTERN,
            all32.as_str(),
            all12.as_str(),
            "111010101010101010101010101010101010111",
            "11011011011011011011011011011011011011011",
        ];
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for p in patterns {
            let t = match KmerTranslator::new(p) {
                Ok(t) => t,
                Err(_) => continue,
            };
            let bits = 2 * t.span() as u32;
            let mask: u128 = if bits >= 128 { u128::MAX } else { (1u128 << bits) - 1 };
            for _ in 0..20_000 {
                let v = (((next() as u128) << 64) | next() as u128) & mask;
                assert_eq!(t.kmer2packed_scalar(v), t.kmer2packed(v), "pattern={p} v={v:#x}");
            }
            for v in [0u128, mask, 1, mask ^ 1] {
                assert_eq!(t.kmer2packed_scalar(v), t.kmer2packed(v), "pattern={p} v={v:#x}");
            }
        }
    }

    #[test]
    fn contiguous_pattern_is_the_identity() {
        // With every position informative, translation must be a no-op.
        let t = KmerTranslator::new(&"1".repeat(32)).unwrap();
        for v in [0u128, 1, 0xDEAD_BEEF, u64::MAX as u128] {
            assert_eq!(t.kmer2packed(v), v as u64);
        }
    }
}
