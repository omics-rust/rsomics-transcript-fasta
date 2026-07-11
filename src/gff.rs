//! Minimal GFF3/GTF line + attribute parsing, in the same hand-rolled,
//! no-extra-dependency style as the rest of the `rsomics-gff-utils` family:
//! tab-split columns, and a single generic attribute reader that accepts
//! both GFF3's `key=value` and GTF's `key "value"` syntax on the same line
//! (matches gffread's own content-sniffed, not extension-gated, detection).

pub struct Line<'a> {
    pub seqid: &'a str,
    pub ftype: &'a str,
    pub start: u64,
    pub end: u64,
    pub strand: u8,
    pub phase: u8,
    pub attrs: &'a str,
}

/// Parse one non-comment, non-blank GFF/GTF line's 9 tab-separated columns.
/// Returns `None` for lines with too few columns (blank trailing lines,
/// stray whitespace) rather than failing the whole file over one blemish —
/// upstream is similarly tolerant of blank/short lines interleaved in
/// otherwise-valid GFF/GTF.
pub fn parse_line(line: &str) -> Option<Line<'_>> {
    let mut cols = line.splitn(9, '\t');
    let seqid = cols.next()?;
    let _source = cols.next()?;
    let ftype = cols.next()?;
    let start: u64 = cols.next()?.parse().ok()?;
    let end: u64 = cols.next()?.parse().ok()?;
    let _score = cols.next()?;
    let strand = cols.next()?.as_bytes().first().copied()?;
    let phase = cols.next()?.as_bytes().first().copied()?;
    let attrs = cols.next().unwrap_or("");
    let (start, end) = if start <= end {
        (start, end)
    } else {
        (end, start)
    };
    Some(Line {
        seqid,
        ftype,
        start,
        end,
        strand,
        phase,
        attrs,
    })
}

/// Iterate `(key, value)` pairs from a raw attribute column, transparently
/// handling `ID=foo;Parent=bar` and `gene_id "g1"; transcript_id "t1";`.
pub fn attr_pairs(attrs: &str) -> impl Iterator<Item = (&str, &str)> {
    attrs.split(';').filter_map(|chunk| {
        let chunk = chunk.trim();
        if chunk.is_empty() {
            return None;
        }
        if let Some(eq) = chunk.find('=') {
            let (k, v) = chunk.split_at(eq);
            Some((k.trim(), v[1..].trim()))
        } else if let Some(sp) = chunk.find(char::is_whitespace) {
            let (k, v) = chunk.split_at(sp);
            Some((k.trim(), v.trim().trim_matches('"')))
        } else {
            None
        }
    })
}

pub fn get_attr<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attr_pairs(attrs).find(|(k, _)| *k == key).map(|(_, v)| v)
}

/// GFF3's `Parent=t1,t2` lets one exon/CDS feed multiple transcripts.
pub fn parent_ids(attrs: &str) -> Vec<&str> {
    get_attr(attrs, "Parent")
        .into_iter()
        .flat_map(|v| v.split(','))
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .collect()
}

#[derive(PartialEq, Eq, Debug)]
pub enum FeatureClass {
    Transcript,
    Exon,
    Cds,
    Other,
}

/// Mirrors gffread's own `fnamelc`-based feature classification (`gff.cpp`):
/// `exon`-suffixed and `*utr*` types splice as exon content; a bare `CDS`,
/// or a `start_codon`/`stop_codon` line, extends the coding window (the
/// latter is how Ensembl-style GTF — CDS excluding the stop — still yields
/// a full-length CDS/protein); anything ending in `rna` or `transcript`
/// is a transcript record.
pub fn classify(ftype: &str) -> FeatureClass {
    let lc = ftype.to_ascii_lowercase();
    if lc == "cds" {
        return FeatureClass::Cds;
    }
    if lc.contains("start") && (lc.contains("codon") || lc.contains("cds")) {
        return FeatureClass::Cds;
    }
    if lc.contains("stop")
        && (lc.contains("codon") || lc.contains("cds"))
        && !lc.contains("redefined")
        && !lc.contains("selenocysteine")
    {
        return FeatureClass::Cds;
    }
    if lc.ends_with("exon") || lc.contains("utr") {
        return FeatureClass::Exon;
    }
    if lc.ends_with("rna") || lc.ends_with("transcript") {
        return FeatureClass::Transcript;
    }
    FeatureClass::Other
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_gff3_style() {
        let l = parse_line("chr1\tsrc\tCDS\t10\t20\t.\t+\t0\tID=c1;Parent=t1").unwrap();
        assert_eq!(l.seqid, "chr1");
        assert_eq!(l.ftype, "CDS");
        assert_eq!(l.start, 10);
        assert_eq!(l.end, 20);
        assert_eq!(l.strand, b'+');
        assert_eq!(l.phase, b'0');
        assert_eq!(get_attr(l.attrs, "Parent"), Some("t1"));
    }

    #[test]
    fn parses_gtf_style() {
        let l =
            parse_line("chr1\tsrc\texon\t10\t20\t.\t-\t.\tgene_id \"g1\"; transcript_id \"t1\";")
                .unwrap();
        assert_eq!(get_attr(l.attrs, "transcript_id"), Some("t1"));
        assert_eq!(get_attr(l.attrs, "gene_id"), Some("g1"));
    }

    #[test]
    fn swaps_inverted_coords() {
        let l = parse_line("chr1\tsrc\texon\t20\t10\t.\t+\t.\t.").unwrap();
        assert_eq!((l.start, l.end), (10, 20));
    }

    #[test]
    fn multi_parent_split() {
        assert_eq!(parent_ids("ID=e1;Parent=t1,t2"), vec!["t1", "t2"]);
    }

    #[test]
    fn classification() {
        assert_eq!(classify("CDS"), FeatureClass::Cds);
        assert_eq!(classify("cds"), FeatureClass::Cds);
        assert_eq!(classify("start_codon"), FeatureClass::Cds);
        assert_eq!(classify("stop_codon"), FeatureClass::Cds);
        assert_eq!(classify("exon"), FeatureClass::Exon);
        assert_eq!(classify("five_prime_UTR"), FeatureClass::Exon);
        assert_eq!(classify("mRNA"), FeatureClass::Transcript);
        assert_eq!(classify("transcript"), FeatureClass::Transcript);
        assert_eq!(classify("ncRNA"), FeatureClass::Transcript);
        assert_eq!(classify("gene"), FeatureClass::Other);
    }
}
