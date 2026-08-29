//! Nucleotide alphabet tables, ported verbatim from `src/types.cpp`.
//!
//! The encoding is CoCo's own: `A=0, C=1, T=2, G=3`. That ordering is not
//! alphabetical and is not the usual `ACGT`; it is chosen so that complementing a
//! base is `x ^ 2` (`A<->T`, `C<->G`), which the reverse-complement code relies on.

pub const ALPHABET_SIZE: usize = 4;

/// code -> character. `int2res` in the C++.
pub const INT2RES: [u8; ALPHABET_SIZE] = [b'A', b'C', b'T', b'G'];

/// code -> complementary code. `int2rev` in the C++.
pub const INT2REV: [u8; ALPHABET_SIZE] = [2, 3, 0, 1];

/// character -> code, `-1` for "not a nucleotide".
///
/// This reproduces a quirk of the original exactly. `types.cpp` declares
/// `int res2int['z' + 1]` at namespace scope (so zero-initialised), then fills only
/// `0..='Z'` with `-1` before writing the eight ACGT/acgt entries. Indices
/// `'['..='z'` that are not `acgt` are therefore left at **0**, which decodes as
/// `A` rather than as invalid. The visible consequence is that lowercase-masked
/// or ambiguity-coded bases in that range are silently read as `A` -- `'n'` (110)
/// becomes `A`, while uppercase `'N'` (78) is correctly invalid. Correction output
/// depends on this, so the port keeps it.
///
/// The C++ array stops at index `'z'`; indexing it with a byte above 122 is out of
/// bounds. Bytes `123..=255` are mapped to `-1` here, which is the only defined
/// choice available and matches the `'Z'`-and-below behaviour for non-nucleotides.
pub static RES2INT: [i8; 256] = build_res2int();

const fn build_res2int() -> [i8; 256] {
    let mut t = [-1i8; 256];
    // C++ leaves '[' ..= 'z' zero-initialised, i.e. code 0 == 'A'.
    let mut i = b'Z' as usize + 1;
    while i <= b'z' as usize {
        t[i] = 0;
        i += 1;
    }
    // Then the explicit ACGT/acgt assignments.
    let mut c = 0usize;
    while c < ALPHABET_SIZE {
        t[INT2RES[c] as usize] = c as i8;
        t[INT2RES[c] as usize + 32] = c as i8;
        c += 1;
    }
    t
}

/// `res2int[c]` with the C++ indexing behaviour.
#[inline(always)]
pub fn res2int(c: u8) -> i8 {
    // Safety-free: RES2INT covers the whole u8 range.
    RES2INT[c as usize]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alphabet_round_trips() {
        for c in 0..ALPHABET_SIZE {
            assert_eq!(res2int(INT2RES[c]), c as i8);
            assert_eq!(res2int(INT2RES[c] + 32), c as i8);
        }
    }

    #[test]
    fn complement_is_xor_two() {
        for c in 0..ALPHABET_SIZE {
            assert_eq!(INT2REV[c], (c ^ 2) as u8);
        }
    }

    #[test]
    fn matches_cpp_initialisation_quirks() {
        // Uppercase non-nucleotides are invalid.
        assert_eq!(res2int(b'N'), -1);
        assert_eq!(res2int(b'X'), -1);
        assert_eq!(res2int(b'-'), -1);
        // ... but the lowercase band silently decodes as 'A'.
        assert_eq!(res2int(b'n'), 0);
        assert_eq!(res2int(b'x'), 0);
        assert_eq!(res2int(b'['), 0);
        // Real lowercase nucleotides still decode correctly.
        assert_eq!(res2int(b'a'), 0);
        assert_eq!(res2int(b'c'), 1);
        assert_eq!(res2int(b't'), 2);
        assert_eq!(res2int(b'g'), 3);
    }
}
