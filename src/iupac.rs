//! Reverse-complement, ported from gffread/gclib's `gdna.cpp` IUPAC tables.
//!
//! Every ambiguity code complements to its own ambiguity code (`gdna.cpp`'s
//! `IUPAC_DEFS`/`IUPAC_COMP` pair), case is preserved, and anything outside
//! the IUPAC alphabet complements to `'N'` rather than erroring — this is
//! upstream's own fail-soft choice for stray bytes in a FASTA body, and
//! byte-identical output requires matching it rather than failing loud here.

fn complement_base(c: u8) -> u8 {
    match c {
        b'A' => b'T',
        b'a' => b't',
        b'C' => b'G',
        b'c' => b'g',
        b'T' => b'A',
        b't' => b'a',
        b'G' => b'C',
        b'g' => b'c',
        b'U' => b'A',
        b'u' => b'a',
        b'M' => b'K',
        b'm' => b'k',
        b'R' => b'Y',
        b'r' => b'y',
        b'W' => b'W',
        b'w' => b'w',
        b'S' => b'S',
        b's' => b's',
        b'Y' => b'R',
        b'y' => b'r',
        b'K' => b'M',
        b'k' => b'm',
        b'V' => b'B',
        b'v' => b'b',
        b'H' => b'D',
        b'h' => b'd',
        b'D' => b'H',
        b'd' => b'h',
        b'B' => b'V',
        b'b' => b'v',
        b'N' => b'N',
        b'n' => b'n',
        b'X' => b'X',
        b'x' => b'x',
        b'-' => b'-',
        b'*' => b'*',
        _ => b'N',
    }
}

/// Append the reverse-complement of `seq` to `out` (IUPAC-aware) without
/// allocating a temporary — the `-` strand splice hot path.
pub fn push_reverse_complement(out: &mut Vec<u8>, seq: &[u8]) {
    out.extend(seq.iter().rev().map(|&c| complement_base(c)));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rc(seq: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        push_reverse_complement(&mut out, seq);
        out
    }

    #[test]
    fn basic() {
        assert_eq!(rc(b"ATGC"), b"GCAT");
    }

    #[test]
    fn preserves_case() {
        assert_eq!(rc(b"atgC"), b"Gcat");
    }

    #[test]
    fn degenerate_codes_self_complement_pairs() {
        assert_eq!(rc(b"RYSWKMBDHVN"), b"NBDHVKMWSRY");
    }

    #[test]
    fn unknown_byte_becomes_n() {
        assert_eq!(rc(b"AZG"), b"CNT");
    }
}
