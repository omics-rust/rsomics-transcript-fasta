//! Exon splicing, ported from gffread/gclib's `GffObj::getSpliced` (`gff.cpp`).
//!
//! One walk over a transcript's (already merged, ascending) exon list
//! produces both plain outputs: the `-w` transcript is the walk over the
//! full exon span; the `-x`/`-y` CDS is the same walk narrowed to the
//! phase-adjusted CDS window. Passing the CDS window alongside a wider walk
//! window (as `-w` does) also recovers the *local*, spliced-sequence
//! position of that CDS — gffread's `CDS=<start>-<end>` defline tag.

use crate::iupac::push_reverse_complement;

/// Splice `exons` (genomic, 1-based inclusive, ascending, non-overlapping)
/// within `[walk_start, walk_end]`, in transcript order (reverse-complemented
/// per exon on `-`). `cds_bounds`, if given, is projected onto the output's
/// 1-based local coordinates — used for the `-w` `CDS=` tag; pass the same
/// window as `cds_bounds` and as the walk bounds to get the phase-adjusted
/// CDS/protein sequence directly (`-x`/`-y`).
pub fn splice(
    genome_seq: &[u8],
    exons: &[(u64, u64)],
    strand: u8,
    walk: (u64, u64),
    cds_bounds: Option<(u64, u64)>,
) -> (Vec<u8>, Option<(u64, u64)>) {
    let (walk_start, walk_end) = walk;
    let mut seq = Vec::with_capacity((walk_end - walk_start + 1) as usize);
    let mut cds_local_start: Option<u64> = None;
    let mut cds_local_end: Option<u64> = None;

    // Exons are ascending; `-` reads them in reverse transcript order. Index
    // both directions instead of cloning the exon list on every call.
    let minus = strand == b'-';
    let n = exons.len();
    for i in 0..n {
        let (seg_start, seg_end) = exons[if minus { n - 1 - i } else { i }];
        if walk_end < seg_start || walk_start > seg_end {
            continue;
        }
        let sg_start = seg_start.max(walk_start);
        let sg_end = seg_end.min(walk_end);
        let piece = &genome_seq[(sg_start - 1) as usize..sg_end as usize];
        if minus {
            push_reverse_complement(&mut seq, piece);
        } else {
            seq.extend_from_slice(piece);
        }
        let s = seq.len() as u64;

        if let Some((cds_start, cds_end)) = cds_bounds {
            if minus {
                if cds_end >= sg_start && cds_end <= sg_end {
                    cds_local_start = Some(s - (cds_end - sg_start));
                }
                if cds_start >= sg_start && cds_start <= sg_end {
                    cds_local_end = Some(s - (cds_start - sg_start));
                }
            } else {
                if cds_start >= sg_start && cds_start <= sg_end {
                    cds_local_start = Some(s - (sg_end - cds_start));
                }
                if cds_end >= sg_start && cds_end <= sg_end {
                    cds_local_end = Some(s - (sg_end - cds_end));
                }
            }
        }
    }

    let cds_local = cds_local_start.zip(cds_local_end);
    (seq, cds_local)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genome() -> Vec<u8> {
        // 1-based coordinates: position i lives at index i-1.
        b"NATGGCCAAACTGGTCAGCACCGGCGACGAGCGCTAAN".to_vec()
    }

    #[test]
    fn plus_strand_whole_span() {
        let g = genome();
        let (seq, cds) = splice(&g, &[(2, 37)], b'+', (2, 37), Some((2, 37)));
        assert_eq!(seq, b"ATGGCCAAACTGGTCAGCACCGGCGACGAGCGCTAA");
        assert_eq!(cds, Some((1, 36)));
    }

    #[test]
    fn minus_strand_revcomps_and_reorders() {
        // exon at genomic 2..=37 on minus strand: transcript-sense is the
        // reverse-complement of the plus-strand slice.
        let g = genome();
        let (seq, _) = splice(&g, &[(2, 37)], b'-', (2, 37), None);
        let mut expected = g[1..37].to_vec();
        expected.reverse();
        let expected: Vec<u8> = expected
            .iter()
            .map(|&c| match c {
                b'A' => b'T',
                b'T' => b'A',
                b'G' => b'C',
                b'C' => b'G',
                x => x,
            })
            .collect();
        assert_eq!(seq, expected);
    }

    #[test]
    fn narrowed_walk_window_is_just_the_cds() {
        let g = genome();
        // CDS occupies 2..=37 too in this fixture; narrow the walk to a sub-range.
        let (seq, cds) = splice(&g, &[(2, 37)], b'+', (5, 10), Some((5, 10)));
        assert_eq!(seq.len(), 6);
        assert_eq!(cds, Some((1, 6)));
    }
}
