//! Byte-identity compat tests against frozen `gffread` v0.12.9 golden output.
//!
//! Runs the library directly (no live `gffread` binary at test time — the
//! goldens under `tests/golden/*.expected` were captured once and committed).
//! Parent re-verification re-diffs against the live binary separately.

use std::fs::File;
use std::path::Path;

use rsomics_transcript_fasta::{ExtractConfig, extract};

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

/// Scratch directory honouring the local KIOXIA build rule when that volume
/// is present, falling back to the system temp dir on CI runners.
fn scratch() -> tempfile::TempDir {
    let kioxia = Path::new("/Volumes/KIOXIA/tmp");
    if kioxia.is_dir() {
        tempfile::tempdir_in(kioxia).unwrap()
    } else {
        tempfile::tempdir().unwrap()
    }
}

fn run(gff: &str, genome: &str, star_stop: bool) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    let dir = scratch();
    let w_path = dir.path().join("out.w.fa");
    let x_path = dir.path().join("out.x.fa");
    let y_path = dir.path().join("out.y.fa");

    let cfg = ExtractConfig {
        want_w: true,
        want_x: true,
        want_y: true,
        w_nocds: false,
        star_stop,
    };
    let gff_path = golden(gff);
    let genome_path = golden(genome);
    extract(
        gff_path.to_str().unwrap(),
        &genome_path,
        &cfg,
        Some(Box::new(File::create(&w_path).unwrap()) as Box<dyn std::io::Write>),
        Some(Box::new(File::create(&x_path).unwrap()) as Box<dyn std::io::Write>),
        Some(Box::new(File::create(&y_path).unwrap()) as Box<dyn std::io::Write>),
    )
    .expect("extract should succeed");

    (
        std::fs::read(&w_path).unwrap(),
        std::fs::read(&x_path).unwrap(),
        std::fs::read(&y_path).unwrap(),
    )
}

fn run_w_nocds(gff: &str, genome: &str) -> Vec<u8> {
    let dir = scratch();
    let w_path = dir.path().join("out.w.fa");
    let cfg = ExtractConfig {
        want_w: true,
        want_x: false,
        want_y: false,
        w_nocds: true,
        star_stop: false,
    };
    let gff_path = golden(gff);
    let genome_path = golden(genome);
    extract(
        gff_path.to_str().unwrap(),
        &genome_path,
        &cfg,
        Some(Box::new(File::create(&w_path).unwrap()) as Box<dyn std::io::Write>),
        None,
        None,
    )
    .expect("extract should succeed");
    std::fs::read(&w_path).unwrap()
}

fn assert_matches_golden(actual: &[u8], expected_name: &str) {
    let expected = std::fs::read(golden(expected_name))
        .unwrap_or_else(|e| panic!("reading golden {expected_name}: {e}"));
    assert_eq!(
        actual,
        expected.as_slice(),
        "mismatch vs golden {expected_name}\n--- actual ---\n{}\n--- expected ---\n{}",
        String::from_utf8_lossy(actual),
        String::from_utf8_lossy(&expected),
    );
}

#[test]
fn gff3_main_fixture_byte_identical() {
    let (w, x, y) = run("genes.gff3", "genome.fa", false);
    assert_matches_golden(&w, "genes.w.expected");
    assert_matches_golden(&x, "genes.x.expected");
    assert_matches_golden(&y, "genes.y.expected");
}

#[test]
fn gtf_main_fixture_byte_identical_and_matches_gff3() {
    let (w, x, y) = run("genes.gtf", "genome.fa", false);
    assert_matches_golden(&w, "genes.w.expected");
    assert_matches_golden(&x, "genes.x.expected");
    assert_matches_golden(&y, "genes.y.expected");
}

#[test]
fn ensembl_style_split_stop_codon_merges_into_cds() {
    let (w, x, y) = run("ensembl_stopcodon.gtf", "small_genome.fa", false);
    assert_matches_golden(&w, "ensembl_stopcodon.w.expected");
    assert_matches_golden(&x, "ensembl_stopcodon.x.expected");
    assert_matches_golden(&y, "ensembl_stopcodon.y.expected");
}

#[test]
fn cds_only_transcript_with_no_exon_lines() {
    let (w, x, y) = run("no_exon_lines.gff3", "small_genome.fa", false);
    assert_matches_golden(&w, "no_exon_lines.w.expected");
    assert_matches_golden(&x, "no_exon_lines.x.expected");
    assert_matches_golden(&y, "no_exon_lines.y.expected");
}

#[test]
fn case_is_preserved_and_dot_phase_treated_as_zero() {
    let (w, x, y) = run("case_preserve.gff3", "small_genome.fa", false);
    assert_matches_golden(&w, "case_preserve.w.expected");
    assert_matches_golden(&x, "case_preserve.x.expected");
    assert_matches_golden(&y, "case_preserve.y.expected");
}

#[test]
fn partial_trailing_codon_kept_in_nt_dropped_in_protein() {
    let (_, x, y) = run("partial_codon.gff3", "partial_codon_genome.fa", false);
    assert_matches_golden(&x, "partial_codon.x.expected");
    assert_matches_golden(&y, "partial_codon.y.expected");
}

#[test]
fn comma_separated_parent_shares_exon_across_transcripts() {
    let (w, x, y) = run("multi_parent.gff3", "multi_parent_genome.fa", false);
    assert_matches_golden(&w, "multi_parent.w.expected");
    assert_matches_golden(&x, "multi_parent.x.expected");
    assert_matches_golden(&y, "multi_parent.y.expected");
}

/// Only the very last translated codon decides whether the trailing stop is
/// dropped — not "the first stop found anywhere in the CDS". A CDS with
/// several premature in-frame stops before a genuine terminal one separates
/// gffread's CDS-validity check (which scans for the *first* stop) from its
/// print-time truncation (`gff_utils.cpp`: `if (cdsaa[aalen-1]=='.') --aalen`),
/// which looks at the last codon only; matching the wrong one drops the wrong
/// residues.
#[test]
fn only_the_final_codon_decides_stop_truncation() {
    let (w, x, y) = run("multistop.gff3", "multistop_genome.fa", false);
    assert_matches_golden(&w, "multistop.w.expected");
    assert_matches_golden(&x, "multistop.x.expected");
    assert_matches_golden(&y, "multistop.y.expected");
}

#[test]
fn star_stop_flag_substitutes_every_remaining_dot() {
    let (_, _, y) = run("multistop.gff3", "multistop_genome.fa", true);
    assert_matches_golden(&y, "multistop.y.star.expected");
}

#[test]
fn w_nocds_flag_suppresses_the_cds_defline_tag() {
    let w = run_w_nocds("multistop.gff3", "multistop_genome.fa");
    assert_matches_golden(&w, "multistop.w.nocds.expected");
}

#[test]
fn childless_transcript_gets_synthesized_exon() {
    // gffread emits a transcript with no exon/CDS children as a single exon
    // spanning the transcript record itself.
    let dir = scratch();
    let g = dir.path().join("g.fa");
    std::fs::write(&g, ">chr1\nACGTACGTACGTACGTACGTACGTACGTAC\n").unwrap();
    let gff = dir.path().join("a.gff3");
    std::fs::write(&gff, "chr1\tsrc\tmRNA\t1\t30\t.\t+\t.\tID=t1\n").unwrap();
    let w = dir.path().join("w.fa");
    let cfg = ExtractConfig {
        want_w: true,
        ..Default::default()
    };
    extract(
        gff.to_str().unwrap(),
        &g,
        &cfg,
        Some(Box::new(File::create(&w).unwrap())),
        None,
        None,
    )
    .unwrap();
    assert_eq!(
        std::fs::read_to_string(&w).unwrap(),
        ">t1\nACGTACGTACGTACGTACGTACGTACGTAC\n"
    );
}

#[test]
fn feature_beyond_chromosome_length_is_a_hard_error() {
    let dir = scratch();
    let g = dir.path().join("g.fa");
    std::fs::write(&g, ">chr1\nACGTACGTAC\n").unwrap();
    let gff = dir.path().join("a.gff3");
    std::fs::write(
        &gff,
        "chr1\tsrc\tmRNA\t1\t50\t.\t+\t.\tID=t1\nchr1\tsrc\texon\t1\t50\t.\t+\t.\tParent=t1\n",
    )
    .unwrap();
    let cfg = ExtractConfig {
        want_w: true,
        ..Default::default()
    };
    let r = extract(
        gff.to_str().unwrap(),
        &g,
        &cfg,
        Some(Box::new(File::create(dir.path().join("w.fa")).unwrap())),
        None,
        None,
    );
    assert!(r.is_err());
}
