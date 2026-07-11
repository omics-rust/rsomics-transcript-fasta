use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const N_GENES_PER_CHROM: usize = 400;
const N_CHROMS: usize = 3;
const SEED: u64 = 0x00FA_5704;

fn xorshift(x: &mut u64) -> u64 {
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    *x
}

fn rand_base(rng: &mut u64) -> u8 {
    b"ACGT"[((xorshift(rng) >> 33) & 3) as usize]
}

/// Multi-exon, both-strand gene models across several chromosomes — enough
/// bulk (~1,200 transcripts, a few MB genome) to keep wall-clock away from
/// pure process-startup noise while staying fast enough for routine `cargo
/// bench` runs; the campaign's release-gate numbers come from a much larger
/// fixture measured separately (see the crate's perf provenance notes).
fn synth_fixtures(fa: &Path, gff: &Path) {
    let mut rng = SEED;
    let f = File::create(fa).expect("create fa");
    let mut fa_w = BufWriter::new(f);
    let g = File::create(gff).expect("create gff");
    let mut gff_w = BufWriter::new(g);
    writeln!(gff_w, "##gff-version 3").unwrap();

    for c in 0..N_CHROMS {
        let chrom = format!("chr{}", c + 1);
        let mut seq = Vec::new();
        let mut records = Vec::new();

        for gidx in 0..N_GENES_PER_CHROM {
            let n_exons = 3 + (xorshift(&mut rng) % 4) as usize;
            let strand = if gidx % 2 == 0 { '+' } else { '-' };
            let mut exon_coords = Vec::new();
            for e in 0..n_exons {
                let exon_len = 80 + (xorshift(&mut rng) % 220) as usize;
                let start = seq.len() as u64 + 1;
                for _ in 0..exon_len {
                    seq.push(rand_base(&mut rng));
                }
                let end = seq.len() as u64;
                exon_coords.push((start, end));
                if e + 1 < n_exons {
                    let intron_len = 200 + (xorshift(&mut rng) % 800) as usize;
                    for _ in 0..intron_len {
                        seq.push(rand_base(&mut rng));
                    }
                }
            }
            let filler_len = 200 + (xorshift(&mut rng) % 800) as usize;
            for _ in 0..filler_len {
                seq.push(rand_base(&mut rng));
            }
            records.push((gidx, strand, exon_coords));
        }

        writeln!(fa_w, ">{chrom}").unwrap();
        for chunk in seq.chunks(70) {
            fa_w.write_all(chunk).unwrap();
            fa_w.write_all(b"\n").unwrap();
        }

        for (gidx, strand, exon_coords) in &records {
            let gene_id = format!("gene_{chrom}_{gidx}");
            let tx_id = format!("tx_{chrom}_{gidx}");
            let tx_start = exon_coords.first().unwrap().0;
            let tx_end = exon_coords.last().unwrap().1;
            writeln!(
                gff_w,
                "{chrom}\tbench\tgene\t{tx_start}\t{tx_end}\t.\t{strand}\t.\tID={gene_id}"
            )
            .unwrap();
            writeln!(
                gff_w,
                "{chrom}\tbench\tmRNA\t{tx_start}\t{tx_end}\t.\t{strand}\t.\tID={tx_id};Parent={gene_id}"
            )
            .unwrap();
            for (i, (s, e)) in exon_coords.iter().enumerate() {
                writeln!(gff_w, "{chrom}\tbench\texon\t{s}\t{e}\t.\t{strand}\t.\tID={tx_id}.exon{i};Parent={tx_id}")
                    .unwrap();
            }
            for (i, (s, e)) in exon_coords.iter().enumerate() {
                writeln!(gff_w, "{chrom}\tbench\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID={tx_id}.cds{i};Parent={tx_id}")
                    .unwrap();
            }
        }
    }
}

/// Scratch directory honouring the local KIOXIA build rule when that volume
/// is present, falling back to the system temp dir on CI runners.
fn scratch_dir() -> PathBuf {
    let kioxia = Path::new("/Volumes/KIOXIA/tmp");
    if kioxia.is_dir() {
        kioxia.to_path_buf()
    } else {
        std::env::temp_dir()
    }
}

fn ensure_fixtures() -> (PathBuf, PathBuf) {
    let dir = scratch_dir();
    let fa = dir.join(format!(
        "rsomics-transcript-fasta-bench-{N_CHROMS}x{N_GENES_PER_CHROM}.fa"
    ));
    let gff = fa.with_extension("gff3");
    if !fa.exists() || !gff.exists() {
        synth_fixtures(&fa, &gff);
    }
    (fa, gff)
}

fn bench(c: &mut Criterion) {
    let (fa, gff) = ensure_fixtures();
    let ours = env!("CARGO_BIN_EXE_rsomics-transcript-fasta");
    let devnull = if cfg!(windows) { "NUL" } else { "/dev/null" };

    let mut group = c.benchmark_group(format!(
        "transcript_fasta/{}x{}",
        N_CHROMS, N_GENES_PER_CHROM
    ));
    group.sample_size(10);

    group.bench_function("rsomics-transcript-fasta -w", |b| {
        b.iter(|| {
            let out = Command::new(ours)
                .args(["-g", fa.to_str().unwrap(), "-w", devnull])
                .arg(&gff)
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });

    group.bench_function("rsomics-transcript-fasta -y", |b| {
        b.iter(|| {
            let out = Command::new(ours)
                .args(["-g", fa.to_str().unwrap(), "-y", devnull])
                .arg(&gff)
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });

    if Command::new("gffread")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
    {
        group.bench_function("gffread -w", |b| {
            b.iter(|| {
                let out = Command::new("gffread")
                    .args(["-g", fa.to_str().unwrap(), "-w", devnull])
                    .arg(&gff)
                    .output()
                    .expect("gffread run");
                assert!(
                    out.status.success(),
                    "{}",
                    String::from_utf8_lossy(&out.stderr)
                );
            });
        });
        group.bench_function("gffread -y", |b| {
            b.iter(|| {
                let out = Command::new("gffread")
                    .args(["-g", fa.to_str().unwrap(), "-y", devnull])
                    .arg(&gff)
                    .output()
                    .expect("gffread run");
                assert!(
                    out.status.success(),
                    "{}",
                    String::from_utf8_lossy(&out.stderr)
                );
            });
        });
    } else {
        eprintln!("gffread not on PATH — skipping upstream comparison");
    }

    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
