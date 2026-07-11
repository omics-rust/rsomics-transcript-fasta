//! Spliced transcript/CDS/protein FASTA extraction from a genome + GFF3/GTF
//! gene model — a Rust port of `gffread`'s `-w`/`-x`/`-y` sequence-extraction
//! mode. See the crate README's `## Origin` section for how the exact
//! output semantics (defline format, wrap width, phase handling, ordering)
//! were nailed down from upstream source + black-box verification.

mod codon;
mod fasta_out;
mod genome;
mod gff;
mod iupac;
mod model;
mod splice;

use std::io::Write;
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

#[derive(Default)]
pub struct ExtractConfig {
    pub want_w: bool,
    pub want_x: bool,
    pub want_y: bool,
    /// `--w-nocds`: suppress the `CDS=<start>-<end>` tag on `-w` deflines.
    pub w_nocds: bool,
    /// `-S`: print `*` instead of `.` for any stop codon that survives
    /// translation (only non-terminal, premature stops ever do — see
    /// `protein_from_cds`).
    pub star_stop: bool,
}

#[derive(Default)]
pub struct ExtractStats {
    pub transcripts_seen: u64,
    pub w_records: u64,
    pub x_records: u64,
    pub y_records: u64,
}

/// Run extraction, writing to whichever of `w_out`/`x_out`/`y_out` are
/// `Some` (mirrors gffread accepting `-w`/`-x`/`-y` in any combination in a
/// single pass over the gene model).
pub fn extract(
    gff_input: &str,
    genome_path: &Path,
    cfg: &ExtractConfig,
    mut w_out: Option<Box<dyn Write>>,
    mut x_out: Option<Box<dyn Write>>,
    mut y_out: Option<Box<dyn Write>>,
) -> Result<ExtractStats> {
    let genome = genome::load(genome_path)?;
    let transcripts = model::load(gff_input)?;

    let mut stats = ExtractStats::default();
    let mut scratch = Vec::new();

    for t in &transcripts {
        stats.transcripts_seen += 1;
        let chrom_seq = genome.get(&t.seqid).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "transcript {}: chromosome {:?} not found in genome FASTA",
                t.id, t.seqid
            ))
        })?;

        if let Some(w) = w_out.as_deref_mut() {
            let full_span = (t.exons[0].0, t.exons[t.exons.len() - 1].1);
            let (spliced, cds_local) =
                splice::splice(chrom_seq, &t.exons, t.strand, full_span, t.cds_window);
            let defline = match (cds_local, cfg.w_nocds) {
                (Some((s, e)), false) => format!("{} CDS={s}-{e}", t.id),
                _ => t.id.clone(),
            };
            fasta_out::write_record(w, &mut scratch, &defline, &spliced, false)
                .map_err(RsomicsError::Io)?;
            if !spliced.is_empty() {
                stats.w_records += 1;
            }
        }

        if (cfg.want_x || cfg.want_y)
            && let Some(window) = t.cds_window
        {
            let (cds_nt, _) = splice::splice(chrom_seq, &t.exons, t.strand, window, None);

            if let Some(x) = x_out.as_deref_mut() {
                fasta_out::write_record(x, &mut scratch, &t.id, &cds_nt, false)
                    .map_err(RsomicsError::Io)?;
                if !cds_nt.is_empty() {
                    stats.x_records += 1;
                }
            }
            if let Some(y) = y_out.as_deref_mut() {
                let protein = protein_from_cds(&cds_nt);
                fasta_out::write_record(y, &mut scratch, &t.id, &protein, cfg.star_stop)
                    .map_err(RsomicsError::Io)?;
                if !protein.is_empty() {
                    stats.y_records += 1;
                }
            }
        }
    }

    if let Some(mut w) = w_out {
        w.flush().map_err(RsomicsError::Io)?;
    }
    if let Some(mut x) = x_out {
        x.flush().map_err(RsomicsError::Io)?;
    }
    if let Some(mut y) = y_out {
        y.flush().map_err(RsomicsError::Io)?;
    }
    Ok(stats)
}

/// Translate a phase-adjusted CDS nucleotide sequence, dropping the very
/// last amino acid only when it's itself a stop codon — matches gffread's
/// actual print-time check (`gff_utils.cpp`: `if (cdsaa[aalen-1]=='.') --aalen`),
/// which looks at the *final* codon specifically, not "the first stop found
/// anywhere in the string". A premature stop earlier in the CDS is left in
/// the printed string regardless of what the last codon turns out to be.
fn protein_from_cds(cds_nt: &[u8]) -> Vec<u8> {
    let mut aa = codon::translate_dna(cds_nt);
    if aa.last() == Some(&b'.') {
        aa.pop();
    }
    aa
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_drops_terminal_stop_only() {
        assert_eq!(protein_from_cds(b"ATGGCCTAA"), b"MA");
    }

    #[test]
    fn protein_keeps_premature_stop() {
        assert_eq!(protein_from_cds(b"ATGTAAGCC"), b"M.A");
    }

    #[test]
    fn protein_drops_terminal_stop_even_with_earlier_premature_stop() {
        assert_eq!(protein_from_cds(b"ATGTAAGCCTAA"), b"M.A");
    }

    #[test]
    fn protein_no_stop_translates_everything() {
        assert_eq!(protein_from_cds(b"ATGGCCAAA"), b"MAK");
    }
}
