//! Builds the transcript model (gene-level grouping is unnecessary here:
//! gffread's default FASTA deflines carry zero gene-derived attributes, so
//! this crate skips gene resolution entirely and groups straight to
//! transcripts, matching what's actually observable in `-w`/`-x`/`-y` output).

use std::collections::HashMap;
use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};

use crate::gff::{self, FeatureClass};

pub struct Transcript {
    pub id: String,
    pub seqid: String,
    pub strand: u8,
    /// Genomic start, used only to order output — the exon list (not this
    /// field) is what defines the transcript's actual splicing span.
    pub start: u64,
    /// Merged, ascending, gap-correct exon spans (falls back to the merged
    /// CDS spans when a transcript has CDS/start_codon/stop_codon lines but
    /// no explicit exon feature — a real minimal-GTF pattern).
    pub exons: Vec<(u64, u64)>,
    /// Phase-adjusted translatable window `[start, end]`, genomic, present
    /// only when the transcript has at least one CDS-classified segment.
    pub cds_window: Option<(u64, u64)>,
}

struct Building {
    seqid: String,
    strand: u8,
    start: u64,
    end: u64,
    exons: Vec<(u64, u64)>,
    cds_segments: Vec<(u64, u64, u8)>,
}

impl Building {
    fn new(seqid: &str, strand: u8, start: u64, end: u64) -> Self {
        Building {
            seqid: seqid.to_owned(),
            strand,
            start,
            end,
            exons: Vec::new(),
            cds_segments: Vec::new(),
        }
    }

    fn touch(&mut self, seqid: &str, strand: u8, start: u64, end: u64) {
        if self.strand == b'.' && strand != b'.' {
            self.strand = strand;
        }
        if self.seqid.is_empty() {
            self.seqid = seqid.to_owned();
        }
        self.start = self.start.min(start);
        self.end = self.end.max(end);
    }
}

fn phase_value(phase: u8) -> u8 {
    match phase {
        b'1' => 1,
        b'2' => 2,
        _ => 0,
    }
}

fn merge_intervals(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    v.sort_unstable_by_key(|&(s, _)| s);
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(v.len());
    for (start, end) in v {
        if let Some(last) = merged.last_mut()
            && start <= last.1
        {
            last.1 = last.1.max(end);
            continue;
        }
        merged.push((start, end));
    }
    merged
}

/// The phase attached to the CDS-type segment at the translation-start
/// boundary: the lowest-start segment on `+`, the highest-end segment on
/// `-` — internal segments' own phase values are never consulted (matches
/// gffread reading only `cdss->First()`/`cdss->Last()`, since a spliced,
/// concatenated CDS is frame-consistent from that one boundary onward).
fn boundary_phase(segments: &[(u64, u64, u8)], strand: u8) -> u8 {
    if strand == b'-' {
        segments
            .iter()
            .max_by_key(|&&(_, end, _)| end)
            .map(|&(_, _, p)| p)
            .unwrap_or(0)
    } else {
        segments
            .iter()
            .min_by_key(|&&(start, _, _)| start)
            .map(|&(_, _, p)| p)
            .unwrap_or(0)
    }
}

fn get_or_create<'b>(
    id: &str,
    seqid: &str,
    strand: u8,
    start: u64,
    end: u64,
    order: &mut Vec<String>,
    by_id: &'b mut HashMap<String, Building>,
) -> &'b mut Building {
    by_id
        .entry(id.to_owned())
        .and_modify(|b| b.touch(seqid, strand, start, end))
        .or_insert_with(|| {
            order.push(id.to_owned());
            Building::new(seqid, strand, start, end)
        })
}

pub fn load(input: &str) -> Result<Vec<Transcript>> {
    let reader: Box<dyn BufRead> = if input == "-" {
        Box::new(std::io::BufReader::new(std::io::stdin()))
    } else {
        let file = std::fs::File::open(input)
            .map_err(|e| RsomicsError::InvalidInput(format!("{input}: {e}")))?;
        Box::new(std::io::BufReader::new(file))
    };

    let mut order: Vec<String> = Vec::new();
    let mut by_id: HashMap<String, Building> = HashMap::new();
    let mut seqid_order: Vec<String> = Vec::new();
    let mut seqid_seen: std::collections::HashSet<String> = std::collections::HashSet::new();

    for raw_line in reader.lines() {
        let raw_line = raw_line.map_err(RsomicsError::Io)?;
        if raw_line.is_empty() || raw_line.starts_with('#') {
            continue;
        }
        let Some(line) = gff::parse_line(&raw_line) else {
            continue;
        };
        match gff::classify(line.ftype) {
            FeatureClass::Transcript => {
                let id = gff::get_attr(line.attrs, "ID")
                    .or_else(|| gff::get_attr(line.attrs, "transcript_id"));
                let Some(id) = id else { continue };
                get_or_create(
                    id,
                    line.seqid,
                    line.strand,
                    line.start,
                    line.end,
                    &mut order,
                    &mut by_id,
                );
                if seqid_seen.insert(line.seqid.to_owned()) {
                    seqid_order.push(line.seqid.to_owned());
                }
            }
            FeatureClass::Exon => {
                let owners = owning_transcripts(line.attrs);
                for owner in &owners {
                    let b = get_or_create(
                        owner,
                        line.seqid,
                        line.strand,
                        line.start,
                        line.end,
                        &mut order,
                        &mut by_id,
                    );
                    b.exons.push((line.start, line.end));
                }
                if !owners.is_empty() && seqid_seen.insert(line.seqid.to_owned()) {
                    seqid_order.push(line.seqid.to_owned());
                }
            }
            FeatureClass::Cds => {
                let owners = owning_transcripts(line.attrs);
                let phase = phase_value(line.phase);
                for owner in &owners {
                    let b = get_or_create(
                        owner,
                        line.seqid,
                        line.strand,
                        line.start,
                        line.end,
                        &mut order,
                        &mut by_id,
                    );
                    b.cds_segments.push((line.start, line.end, phase));
                }
                if !owners.is_empty() && seqid_seen.insert(line.seqid.to_owned()) {
                    seqid_order.push(line.seqid.to_owned());
                }
            }
            FeatureClass::Other => {}
        }
    }

    let seqid_rank: HashMap<&str, usize> = seqid_order
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let mut transcripts: Vec<Transcript> = order
        .into_iter()
        .filter_map(|id| {
            let b = by_id.remove(&id)?;
            let exons = if b.exons.is_empty() {
                merge_intervals(b.cds_segments.iter().map(|&(s, e, _)| (s, e)).collect())
            } else {
                merge_intervals(b.exons)
            };
            let cds_window = if b.cds_segments.is_empty() {
                None
            } else {
                let cds_start = b.cds_segments.iter().map(|&(s, _, _)| s).min().unwrap();
                let cds_end = b.cds_segments.iter().map(|&(_, e, _)| e).max().unwrap();
                let phase = boundary_phase(&b.cds_segments, b.strand) as u64;
                let window = if b.strand == b'-' {
                    (cds_start, cds_end.saturating_sub(phase))
                } else {
                    (cds_start + phase, cds_end)
                };
                Some(window)
            };
            Some(Transcript {
                id,
                seqid: b.seqid,
                strand: b.strand,
                start: b.start,
                exons,
                cds_window,
            })
        })
        .collect();

    transcripts.sort_by_key(|t| {
        let rank = seqid_rank
            .get(t.seqid.as_str())
            .copied()
            .unwrap_or(usize::MAX);
        (rank, t.start)
    });

    Ok(transcripts)
}

fn owning_transcripts(attrs: &str) -> Vec<&str> {
    let parents = gff::parent_ids(attrs);
    if !parents.is_empty() {
        return parents;
    }
    gff::get_attr(attrs, "transcript_id").into_iter().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_overlapping() {
        assert_eq!(
            merge_intervals(vec![(1, 10), (5, 15), (20, 30)]),
            vec![(1, 15), (20, 30)]
        );
    }

    #[test]
    fn merges_touching_zero_gap() {
        assert_eq!(merge_intervals(vec![(1, 10), (10, 20)]), vec![(1, 20)]);
    }
}
