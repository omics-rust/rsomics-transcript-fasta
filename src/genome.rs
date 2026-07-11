//! Genome FASTA loading — one in-memory pass, `{seqid → sequence}`, the same
//! pattern as `rsomics-bed-getfasta`: CLI-scale genomes fit comfortably in
//! memory, and random access into an owned `Vec<u8>` beats gffread's own
//! on-disk `.fai`-seek approach for a batch, whole-file extraction run.

use std::collections::HashMap;
use std::path::Path;

use needletail::parse_fastx_file;
use rsomics_common::{Result, RsomicsError};

pub fn load(path: &Path) -> Result<HashMap<String, Vec<u8>>> {
    let mut map = HashMap::new();
    let mut reader = parse_fastx_file(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    while let Some(rec) = reader.next() {
        let rec = rec.map_err(|e| RsomicsError::InvalidInput(format!("FASTA parse error: {e}")))?;
        let id = String::from_utf8_lossy(rec.id())
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_owned();
        map.insert(id, rec.seq().into_owned());
    }
    Ok(map)
}
