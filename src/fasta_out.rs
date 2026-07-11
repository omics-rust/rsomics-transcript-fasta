//! FASTA record writer matching gffread's `printFasta` (`gff_utils.cpp`):
//! a hardcoded 70-column wrap (not configurable upstream, so not configurable
//! here either) and exactly one trailing newline regardless of whether the
//! last line landed on a wrap boundary.

use std::io::{self, Write};

const WRAP: usize = 70;

/// Write one `>defline` + wrapped-sequence record. A record with an empty
/// sequence is skipped entirely (no header, no body) — the CDS/protein
/// output modes never emit a record for a transcript whose translatable
/// window comes out empty.
///
/// The whole record is assembled in the caller-owned `scratch` buffer and
/// flushed with a single `write_all`, so a whole-transcriptome run costs one
/// (virtual) write per record rather than two per 70-column line.
pub fn write_record(
    out: &mut dyn Write,
    scratch: &mut Vec<u8>,
    defline: &str,
    seq: &[u8],
    star_stop: bool,
) -> io::Result<()> {
    if seq.is_empty() {
        return Ok(());
    }
    scratch.clear();
    scratch.push(b'>');
    scratch.extend_from_slice(defline.as_bytes());
    scratch.push(b'\n');
    for chunk in seq.chunks(WRAP) {
        if star_stop {
            scratch.extend(chunk.iter().map(|&b| if b == b'.' { b'*' } else { b }));
        } else {
            scratch.extend_from_slice(chunk);
        }
        scratch.push(b'\n');
    }
    out.write_all(scratch)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wraps_at_70() {
        let mut buf = Vec::new();
        let mut scratch = Vec::new();
        let seq = vec![b'A'; 75];
        write_record(&mut buf, &mut scratch, "id1", &seq, false).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], ">id1");
        assert_eq!(lines[1].len(), 70);
        assert_eq!(lines[2].len(), 5);
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn exact_multiple_of_70_has_no_blank_trailing_line() {
        let mut buf = Vec::new();
        let mut scratch = Vec::new();
        let seq = vec![b'A'; 70];
        write_record(&mut buf, &mut scratch, "id1", &seq, false).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert_eq!(text, format!(">id1\n{}\n", "A".repeat(70)));
    }

    #[test]
    fn empty_sequence_emits_nothing() {
        let mut buf = Vec::new();
        let mut scratch = Vec::new();
        write_record(&mut buf, &mut scratch, "id1", &[], false).unwrap();
        assert!(buf.is_empty());
    }

    #[test]
    fn star_stop_substitutes_dot() {
        let mut buf = Vec::new();
        let mut scratch = Vec::new();
        write_record(&mut buf, &mut scratch, "id1", b"MA.K", true).unwrap();
        assert_eq!(buf, b">id1\nMA*K\n");
    }
}
