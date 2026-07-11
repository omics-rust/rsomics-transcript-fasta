use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use serde::Serialize;

use rsomics_transcript_fasta::{ExtractConfig, extract};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const TAGLINE: &str =
    "Extract spliced transcript/CDS/protein FASTA from a genome + GFF3/GTF — gffread port.";

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-transcript-fasta",
    version,
    about,
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Input GFF3/GTF gene model. `-` reads stdin.
    input: String,

    /// Genome FASTA (required — the source of every extracted base).
    #[arg(short = 'g', long = "genome")]
    genome: PathBuf,

    /// Write spliced exon (transcript) FASTA here.
    #[arg(short = 'w', long = "exons")]
    exons_out: Option<PathBuf>,

    /// Write spliced CDS nucleotide FASTA here.
    #[arg(short = 'x', long = "cds")]
    cds_out: Option<PathBuf>,

    /// Write translated CDS protein FASTA here.
    #[arg(short = 'y', long = "protein")]
    protein_out: Option<PathBuf>,

    /// For `-w`, omit the `CDS=<start>-<end>` defline tag.
    #[arg(long = "w-nocds")]
    w_nocds: bool,

    /// For `-y`, print `*` instead of `.` for any stop codon left in the
    /// translation (only a premature, non-terminal stop ever survives).
    #[arg(short = 'S', long = "star-stop")]
    star_stop: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

#[derive(Serialize)]
pub struct ExtractReport {
    pub input: String,
    pub genome: String,
    pub transcripts_seen: u64,
    pub exons_records: u64,
    pub cds_records: u64,
    pub protein_records: u64,
}

/// A whole-transcriptome FASTA run emits tens of megabytes; a 1 MiB block
/// buffer keeps the write-syscall count low (tens, not hundreds of thousands),
/// which matters most on high-per-syscall-latency filesystems and is free
/// everywhere else. Matches gffread's coarse buffered `FILE*` writes.
const OUT_BUF: usize = 1024 * 1024;

fn open_output(path: &std::path::Path, json_mode: bool) -> Result<Box<dyn std::io::Write>> {
    use std::io::BufWriter;
    if path.as_os_str() == "-" {
        if json_mode {
            Ok(Box::new(std::io::sink()))
        } else {
            Ok(Box::new(BufWriter::with_capacity(
                OUT_BUF,
                std::io::stdout().lock(),
            )))
        }
    } else {
        Ok(Box::new(BufWriter::with_capacity(
            OUT_BUF,
            std::fs::File::create(path).map_err(RsomicsError::Io)?,
        )))
    }
}

impl Cli {
    pub fn execute(&self) -> Result<ExtractReport> {
        if self.exons_out.is_none() && self.cds_out.is_none() && self.protein_out.is_none() {
            return Err(RsomicsError::InvalidInput(
                "at least one of -w/--exons, -x/--cds, -y/--protein is required".into(),
            ));
        }

        let cfg = ExtractConfig {
            want_w: self.exons_out.is_some(),
            want_x: self.cds_out.is_some(),
            want_y: self.protein_out.is_some(),
            w_nocds: self.w_nocds,
            star_stop: self.star_stop,
        };

        let w_writer = self
            .exons_out
            .as_deref()
            .map(|p| open_output(p, self.common.json))
            .transpose()?;
        let x_writer = self
            .cds_out
            .as_deref()
            .map(|p| open_output(p, self.common.json))
            .transpose()?;
        let y_writer = self
            .protein_out
            .as_deref()
            .map(|p| open_output(p, self.common.json))
            .transpose()?;

        let stats = extract(
            &self.input,
            &self.genome,
            &cfg,
            w_writer,
            x_writer,
            y_writer,
        )?;

        Ok(ExtractReport {
            input: self.input.clone(),
            genome: self.genome.display().to_string(),
            transcripts_seen: stats.transcripts_seen,
            exons_records: stats.w_records,
            cds_records: stats.x_records,
            protein_records: stats.y_records,
        })
    }
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        Cli::execute(&self)?;
        Ok(())
    }

    fn run(self) -> std::process::ExitCode {
        let common = self.common().clone();
        rsomics_common::run(&common, Self::meta(), move || Cli::execute(&self))
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: TAGLINE,
    origin: Some(Origin {
        upstream: "gffread",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.12688/f1000research.23297.2"),
    }),
    usage_lines: &[
        "-g <genome.fa> -w <exons.fa> [-x <cds.fa>] [-y <protein.fa>] <in.gff3|in.gtf|->",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('g'),
                long: "genome",
                aliases: &[],
                value: Some("<fasta>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Genome FASTA — the source of every extracted base.",
                why_default: None,
            },
            FlagSpec {
                short: Some('w'),
                long: "exons",
                aliases: &[],
                value: Some("<fasta>"),
                type_hint: Some("Path"),
                required: false,
                default: None,
                description: "Write spliced exon (transcript) FASTA here.",
                why_default: None,
            },
            FlagSpec {
                short: Some('x'),
                long: "cds",
                aliases: &[],
                value: Some("<fasta>"),
                type_hint: Some("Path"),
                required: false,
                default: None,
                description: "Write spliced CDS nucleotide FASTA here.",
                why_default: None,
            },
            FlagSpec {
                short: Some('y'),
                long: "protein",
                aliases: &[],
                value: Some("<fasta>"),
                type_hint: Some("Path"),
                required: false,
                default: None,
                description: "Write translated CDS protein FASTA here.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "w-nocds",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "For -w, omit the CDS=<start>-<end> defline tag.",
                why_default: None,
            },
            FlagSpec {
                short: Some('S'),
                long: "star-stop",
                aliases: &[],
                value: None,
                type_hint: None,
                required: false,
                default: None,
                description: "For -y, print '*' instead of '.' for a surviving stop codon.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Spliced transcript FASTA",
            command: "rsomics-transcript-fasta -g genome.fa -w transcripts.fa annotation.gff3",
        },
        Example {
            description: "CDS nucleotide + protein FASTA together",
            command: "rsomics-transcript-fasta -g genome.fa -x cds.fa -y protein.fa annotation.gtf",
        },
        Example {
            description: "Protein FASTA with '*' stop codons",
            command: "rsomics-transcript-fasta -g genome.fa -y protein.fa -S annotation.gff3",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
