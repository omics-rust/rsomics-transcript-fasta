# rsomics-transcript-fasta

Extract spliced transcript, CDS, and protein FASTA from a genome FASTA plus
a GFF3/GTF gene model — a Rust port of `gffread`'s sequence-extraction mode
(`-w`/`-x`/`-y`).

## Install

```
cargo install rsomics-transcript-fasta
```

## Usage

```
# Spliced transcript (exon-concatenated) FASTA
rsomics-transcript-fasta -g genome.fa -w transcripts.fa annotation.gff3

# CDS nucleotide + translated protein FASTA together, one pass
rsomics-transcript-fasta -g genome.fa -x cds.fa -y protein.fa annotation.gtf

# Protein FASTA with '*' instead of '.' for any surviving stop codon
rsomics-transcript-fasta -g genome.fa -y protein.fa -S annotation.gff3

# stdin
cat annotation.gff3 | rsomics-transcript-fasta -g genome.fa -w out.fa -
```

## Semantics (verified against `gffread` v0.12.9, both by reading
`gpertea/gffread`/`gpertea/gclib` source and by black-box comparison)

- **Defline**: `>id` for `-x`/`-y`; `-w` appends ` CDS=<start>-<end>` (1-based,
  positions within the *spliced* transcript sequence) whenever the
  transcript has a CDS — suppress with `--w-nocds`. No other attributes are
  ever added by default (gffread itself drops all non-essential GFF
  attributes from the FASTA defline unless its own `-F`/`-G` flags are
  given — see Deferred below).
- **Line wrap**: hardcoded 70 columns, exactly one trailing newline per
  record regardless of whether the sequence length is an exact multiple of
  the wrap width. Not configurable (gffread doesn't expose this either).
- **Record order**: transcripts are grouped by chromosome in the order that
  chromosome is first seen in the gene-model file (not alphabetical, not
  genome-FASTA order), then sorted by ascending start position within each
  chromosome; equal-start transcripts keep their relative file order (a
  stable sort, ties broken by input order — this matters for the common
  multi-isoform-same-start case).
- **Strand**: `-` strand transcripts are spliced 5'→3' in transcript order
  (highest-genomic-coordinate exon first) and each exon's bases are
  reverse-complemented; the full IUPAC ambiguity-code complement table is
  used (`M↔K`, `R↔Y`, `W↔W`, `S↔S`, `V↔B`, `H↔D`, `N↔N`, ...), matching
  `gclib`'s `gdna.cpp`. Case is preserved; a byte outside the IUPAC alphabet
  complements to `N` (upstream's own fail-soft choice for stray FASTA bytes,
  reproduced here for byte-identical output rather than failing loud).
- **Phase**: the GFF `phase` column on the *boundary* CDS segment — the
  lowest-start segment for `+`, the highest-end segment for `-` — is the
  only one that matters; phase values on internal CDS segments are read
  from the file (for GFF3 validity) but never consulted, exactly like
  gffread's `cdss->First()`/`cdss->Last()`. A phase of 1 or 2 skips that
  many bases from the translation-start boundary before the first codon —
  this is how a 5'-incomplete CDS annotation is handled.
- **Ensembl-style GTF** (`CDS` excluding the stop codon, with a separate
  `stop_codon` feature) is supported: `start_codon`/`stop_codon` lines
  extend the same coding window as `CDS` lines, so the full stop-codon-
  inclusive CDS/protein comes out whether the stop is annotated as part of
  the `CDS` feature or as a separate `stop_codon` feature — verified
  black-box against real gffread on both layouts.
- **No exon lines, CDS only**: a transcript whose only children are
  `CDS`/`start_codon`/`stop_codon` features (a real minimal-GTF pattern) has
  its exon list derived from the merged CDS spans, matching gffread.
- **`-x` (CDS nt)** always includes the stop codon's nucleotides when one is
  annotated; it is never truncated to a codon-count multiple — a trailing
  1-2 nt past the last full codon is emitted verbatim (only the *protein*
  translation drops a partial trailing codon, per `translateDNA`'s floor
  division).
- **`-y` (protein)**: standard genetic code (NCBI `transl_table=1`) plus the
  same IUPAC-degenerate codon resolutions gffread's `codons.cpp` table
  encodes (e.g. `ACN→T`, `CGN→R`) — an ambiguous triplet with no single
  unambiguous amino acid translates to `X`. The check for whether to drop
  the trailing stop looks at the **very last translated codon only**
  (gffread's `gff_utils.cpp`: `if (cdsaa[aalen-1]=='.') --aalen`) — it is
  *not* "the first stop found anywhere in the CDS". A CDS riddled with
  premature in-frame stops (common in synthetic/non-ORF test data, and in
  real pseudogenes) still has its trailing stop dropped whenever the *last*
  codon happens to be one, while every earlier stop stays in the string as
  `.` (or `*` with `-S`) regardless of how many there are. A CDS with no
  stop codon at all (3'-incomplete) translates every full codon with no
  truncation.
- **Case-insensitive translation**: `atg`/`ATG`/`Atg` all translate to `M`.
- Comma-separated `Parent=t1,t2` (one exon/CDS feeding multiple GFF3
  transcripts) is supported.
- **Not matched**: a single `CDS` feature whose genomic span crosses an
  intron into a second exon (i.e. one CDS line covering more than one exon)
  is not a well-formed GFF3/GTF gene model — every real annotation pipeline
  emits one CDS line per exon — and this crate's `-x`/`-y` clip such a
  feature to its per-exon overlap rather than reproducing gffread's raw,
  un-clipped internal segment list for that malformed shape. Standard,
  per-exon-split CDS annotations (verified extensively, including
  multi-parent shared exons) are unaffected and byte-identical.
- GFF3 and GTF are both accepted on the same code path — the attribute
  reader transparently parses `key=value` and `key "value"` syntax, so
  format is inferred from content rather than a file-extension check.

### Deferred (documented, not silently ignored)

The following gffread flags are out of scope for this crate and are
**rejected as unrecognized arguments** by clap rather than silently
ignored: `-F`/`-G`/`--keep-exon-attrs`/`-D` (attribute preservation —
default output already carries zero extra attributes, matching what
gffread itself does without these flags); `-W` (exon-coordinate annotation
in the defline); `--w-add`/`-u` (padding / unspliced span output); `-i`/`-l`
(intron/length filters); `-M`/`-K`/`-Q`/`-Y`/`--merge`/`--cluster-only`
(locus clustering); `-V`/`-H`/`-B`/`-J`/`-N`/`--adj-stop` (CDS validity
filtering/adjustment); `-r`/`-R`/`--jmatch` (region filters); `--ids`/
`--nids`; `--sort-alpha`/`--sort-by`; `--table`; `-C`/`--nc`/`-U`; `-m`;
`-t`; `--bed`/`--gtf`/`--tlf` output-format conversion (this crate's output
is always FASTA); `-g <dir>` (per-chromosome FASTA directory — only a
single multi-FASTA file is supported). These are gffread's filtering,
clustering, and format-conversion features; this crate's scope is
exclusively the `-g`+`-w`/`-x`/`-y` sequence-extraction path.

## Origin

This crate is an independent Rust reimplementation of `gffread`'s
sequence-extraction mode, informed by reading the upstream MIT-licensed
source (`gpertea/gffread`, using `gpertea/gclib`) — specifically
`gff_utils.cpp`'s `process_transcript`/`printFasta`, `gclib/gff.cpp`'s
`GffObj::getSpliced` and feature-type classification, and `gclib/codons.cpp`
/`gclib/gdna.cpp`'s translation and reverse-complement tables — cross-checked
against the installed `gffread` v0.12.9 binary's actual output on
purpose-built fixtures (multi-exon genes on both strands, non-zero CDS
phase, multi-transcript genes, GFF3 and GTF variants, Ensembl-style
split stop-codon CDS, CDS-only-no-exon-lines transcripts).

Reference: Pertea G, Pertea M. *GFF Utilities: GffRead and GffCompare.*
F1000Research 2020, 9:304 (doi: 10.12688/f1000research.23297.2).

License: MIT OR Apache-2.0.
Upstream credit: [gffread](https://github.com/gpertea/gffread) (MIT),
using [gclib](https://github.com/gpertea/gclib) (Artistic License 2.0).
