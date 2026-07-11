//! Codon-to-amino-acid table, ported verbatim from gffread/gclib's
//! `codons.cpp::codonData` (NCBI translation table 1, extended with the
//! standard IUPAC-degenerate entries gffread also resolves unambiguously).
//! Byte-identical translation requires this exact degenerate set, not just
//! the four unambiguous nucleotides.

/// Translate one codon (already uppercased) to its amino acid, or 'X' if
/// the triplet has no unambiguous translation.
pub fn translate_codon(a: u8, b: u8, c: u8) -> u8 {
    match (a, b, c) {
        (b'A', b'A', b'A') => b'K',
        (b'A', b'A', b'C') => b'N',
        (b'A', b'A', b'G') => b'K',
        (b'A', b'A', b'R') => b'K',
        (b'A', b'A', b'T') => b'N',
        (b'A', b'A', b'Y') => b'N',
        (b'A', b'C', b'A') => b'T',
        (b'A', b'C', b'B') => b'T',
        (b'A', b'C', b'C') => b'T',
        (b'A', b'C', b'D') => b'T',
        (b'A', b'C', b'G') => b'T',
        (b'A', b'C', b'H') => b'T',
        (b'A', b'C', b'K') => b'T',
        (b'A', b'C', b'M') => b'T',
        (b'A', b'C', b'N') => b'T',
        (b'A', b'C', b'R') => b'T',
        (b'A', b'C', b'S') => b'T',
        (b'A', b'C', b'T') => b'T',
        (b'A', b'C', b'V') => b'T',
        (b'A', b'C', b'W') => b'T',
        (b'A', b'C', b'Y') => b'T',
        (b'A', b'G', b'A') => b'R',
        (b'A', b'G', b'C') => b'S',
        (b'A', b'G', b'G') => b'R',
        (b'A', b'G', b'R') => b'R',
        (b'A', b'G', b'T') => b'S',
        (b'A', b'G', b'Y') => b'S',
        (b'A', b'T', b'A') => b'I',
        (b'A', b'T', b'C') => b'I',
        (b'A', b'T', b'G') => b'M',
        (b'A', b'T', b'H') => b'I',
        (b'A', b'T', b'M') => b'I',
        (b'A', b'T', b'T') => b'I',
        (b'A', b'T', b'W') => b'I',
        (b'A', b'T', b'Y') => b'I',
        (b'C', b'A', b'A') => b'Q',
        (b'C', b'A', b'C') => b'H',
        (b'C', b'A', b'G') => b'Q',
        (b'C', b'A', b'R') => b'Q',
        (b'C', b'A', b'T') => b'H',
        (b'C', b'A', b'Y') => b'H',
        (b'C', b'C', b'A') => b'P',
        (b'C', b'C', b'B') => b'P',
        (b'C', b'C', b'C') => b'P',
        (b'C', b'C', b'D') => b'P',
        (b'C', b'C', b'G') => b'P',
        (b'C', b'C', b'H') => b'P',
        (b'C', b'C', b'K') => b'P',
        (b'C', b'C', b'M') => b'P',
        (b'C', b'C', b'N') => b'P',
        (b'C', b'C', b'R') => b'P',
        (b'C', b'C', b'S') => b'P',
        (b'C', b'C', b'T') => b'P',
        (b'C', b'C', b'V') => b'P',
        (b'C', b'C', b'W') => b'P',
        (b'C', b'C', b'Y') => b'P',
        (b'C', b'G', b'A') => b'R',
        (b'C', b'G', b'B') => b'R',
        (b'C', b'G', b'C') => b'R',
        (b'C', b'G', b'D') => b'R',
        (b'C', b'G', b'G') => b'R',
        (b'C', b'G', b'H') => b'R',
        (b'C', b'G', b'K') => b'R',
        (b'C', b'G', b'M') => b'R',
        (b'C', b'G', b'N') => b'R',
        (b'C', b'G', b'R') => b'R',
        (b'C', b'G', b'S') => b'R',
        (b'C', b'G', b'T') => b'R',
        (b'C', b'G', b'V') => b'R',
        (b'C', b'G', b'W') => b'R',
        (b'C', b'G', b'Y') => b'R',
        (b'C', b'T', b'A') => b'L',
        (b'C', b'T', b'B') => b'L',
        (b'C', b'T', b'C') => b'L',
        (b'C', b'T', b'D') => b'L',
        (b'C', b'T', b'G') => b'L',
        (b'C', b'T', b'H') => b'L',
        (b'C', b'T', b'K') => b'L',
        (b'C', b'T', b'M') => b'L',
        (b'C', b'T', b'N') => b'L',
        (b'C', b'T', b'R') => b'L',
        (b'C', b'T', b'S') => b'L',
        (b'C', b'T', b'T') => b'L',
        (b'C', b'T', b'V') => b'L',
        (b'C', b'T', b'W') => b'L',
        (b'C', b'T', b'Y') => b'L',
        (b'G', b'A', b'A') => b'E',
        (b'G', b'A', b'C') => b'D',
        (b'G', b'A', b'G') => b'E',
        (b'G', b'A', b'R') => b'E',
        (b'G', b'A', b'T') => b'D',
        (b'G', b'A', b'Y') => b'D',
        (b'G', b'C', b'A') => b'A',
        (b'G', b'C', b'B') => b'A',
        (b'G', b'C', b'C') => b'A',
        (b'G', b'C', b'D') => b'A',
        (b'G', b'C', b'G') => b'A',
        (b'G', b'C', b'H') => b'A',
        (b'G', b'C', b'K') => b'A',
        (b'G', b'C', b'M') => b'A',
        (b'G', b'C', b'N') => b'A',
        (b'G', b'C', b'R') => b'A',
        (b'G', b'C', b'S') => b'A',
        (b'G', b'C', b'T') => b'A',
        (b'G', b'C', b'V') => b'A',
        (b'G', b'C', b'W') => b'A',
        (b'G', b'C', b'Y') => b'A',
        (b'G', b'G', b'A') => b'G',
        (b'G', b'G', b'B') => b'G',
        (b'G', b'G', b'C') => b'G',
        (b'G', b'G', b'D') => b'G',
        (b'G', b'G', b'G') => b'G',
        (b'G', b'G', b'H') => b'G',
        (b'G', b'G', b'K') => b'G',
        (b'G', b'G', b'M') => b'G',
        (b'G', b'G', b'N') => b'G',
        (b'G', b'G', b'R') => b'G',
        (b'G', b'G', b'S') => b'G',
        (b'G', b'G', b'T') => b'G',
        (b'G', b'G', b'V') => b'G',
        (b'G', b'G', b'W') => b'G',
        (b'G', b'G', b'Y') => b'G',
        (b'G', b'T', b'A') => b'V',
        (b'G', b'T', b'B') => b'V',
        (b'G', b'T', b'C') => b'V',
        (b'G', b'T', b'D') => b'V',
        (b'G', b'T', b'G') => b'V',
        (b'G', b'T', b'H') => b'V',
        (b'G', b'T', b'K') => b'V',
        (b'G', b'T', b'M') => b'V',
        (b'G', b'T', b'N') => b'V',
        (b'G', b'T', b'R') => b'V',
        (b'G', b'T', b'S') => b'V',
        (b'G', b'T', b'T') => b'V',
        (b'G', b'T', b'V') => b'V',
        (b'G', b'T', b'W') => b'V',
        (b'G', b'T', b'Y') => b'V',
        (b'M', b'G', b'A') => b'R',
        (b'M', b'G', b'G') => b'R',
        (b'M', b'G', b'R') => b'R',
        (b'N', b'N', b'N') => b'X',
        (b'R', b'A', b'Y') => b'B',
        (b'S', b'A', b'R') => b'Z',
        (b'T', b'A', b'A') => b'.',
        (b'T', b'A', b'C') => b'Y',
        (b'T', b'A', b'G') => b'.',
        (b'T', b'A', b'R') => b'.',
        (b'T', b'A', b'T') => b'Y',
        (b'T', b'A', b'Y') => b'Y',
        (b'T', b'C', b'A') => b'S',
        (b'T', b'C', b'B') => b'S',
        (b'T', b'C', b'C') => b'S',
        (b'T', b'C', b'D') => b'S',
        (b'T', b'C', b'G') => b'S',
        (b'T', b'C', b'H') => b'S',
        (b'T', b'C', b'K') => b'S',
        (b'T', b'C', b'M') => b'S',
        (b'T', b'C', b'N') => b'S',
        (b'T', b'C', b'R') => b'S',
        (b'T', b'C', b'S') => b'S',
        (b'T', b'C', b'T') => b'S',
        (b'T', b'C', b'V') => b'S',
        (b'T', b'C', b'W') => b'S',
        (b'T', b'C', b'Y') => b'S',
        (b'T', b'G', b'A') => b'.',
        (b'T', b'G', b'C') => b'C',
        (b'T', b'G', b'G') => b'W',
        (b'T', b'G', b'T') => b'C',
        (b'T', b'G', b'Y') => b'C',
        (b'T', b'R', b'A') => b'.',
        (b'T', b'T', b'A') => b'L',
        (b'T', b'T', b'C') => b'F',
        (b'T', b'T', b'G') => b'L',
        (b'T', b'T', b'R') => b'L',
        (b'T', b'T', b'T') => b'F',
        (b'T', b'T', b'Y') => b'F',
        (b'X', b'X', b'X') => b'X',
        (b'Y', b'T', b'A') => b'L',
        (b'Y', b'T', b'G') => b'L',
        (b'Y', b'T', b'R') => b'L',
        _ => b'X',
    }
}

/// Translate a spliced CDS nucleotide sequence into amino acids, matching
/// `gclib::translateDNA`: reads whole codons only, silently dropping any
/// trailing 1-2 leftover bases (no partial-codon error, no padding).
pub fn translate_dna(nt: &[u8]) -> Vec<u8> {
    let codon_count = nt.len() / 3;
    let mut aa = Vec::with_capacity(codon_count);
    for i in 0..codon_count {
        let base = i * 3;
        aa.push(translate_codon(
            nt[base].to_ascii_uppercase(),
            nt[base + 1].to_ascii_uppercase(),
            nt[base + 2].to_ascii_uppercase(),
        ));
    }
    aa
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_codons() {
        assert_eq!(translate_codon(b'A', b'T', b'G'), b'M');
        assert_eq!(translate_codon(b'T', b'A', b'A'), b'.');
        assert_eq!(translate_codon(b'T', b'A', b'G'), b'.');
        assert_eq!(translate_codon(b'T', b'G', b'A'), b'.');
    }

    #[test]
    fn degenerate_codons() {
        assert_eq!(translate_codon(b'A', b'C', b'N'), b'T');
        assert_eq!(translate_codon(b'C', b'G', b'N'), b'R');
        assert_eq!(translate_codon(b'N', b'N', b'N'), b'X');
        assert_eq!(translate_codon(b'A', b'T', b'N'), b'X');
    }

    #[test]
    fn case_insensitive_via_translate_dna() {
        assert_eq!(translate_dna(b"atgGCCaaa"), b"MAK");
    }

    #[test]
    fn drops_trailing_partial_codon() {
        assert_eq!(translate_dna(b"ATGGCCAA"), b"MA");
    }
}
