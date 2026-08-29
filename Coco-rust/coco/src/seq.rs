//! FASTA/FASTQ reading and writing.
//!
//! Replaces `lib/kseq/kseq.h` plus `src/KSeqWrapper.cpp` and `src/SequenceInfo.h`.
//! The parsing rules are kseq's, reproduced deliberately:
//!
//! * the record name runs to the first *whitespace* character, not the first
//!   space -- so a name can be terminated by the newline itself;
//! * a comment exists only if the name was terminated by something other than a
//!   newline, and is the remainder of that line;
//! * sequence lines accumulate until a line begins with `>`, `@` or `+`, so a
//!   multi-line FASTA record is joined without separators;
//! * a trailing `\r` is stripped from every line, so CRLF input parses the same.
//!
//! Unlike kseq, this reads through one large buffer and hands out whole batches,
//! which is what lets the correction loop run in parallel.

use std::io::{BufRead, Read, Write};

/// The C++ `SequenceInfo`.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct SequenceInfo {
    pub name: Vec<u8>,
    pub comment: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// `'@'` for FASTQ records, `'>'` for FASTA.
    pub sep: u8,
}

impl SequenceInfo {
    pub fn clear(&mut self) {
        self.name.clear();
        self.comment.clear();
        self.seq.clear();
        self.qual.clear();
        self.sep = b'>';
    }
}

/// `SeqInfoMode` from `src/SequenceInfo.h`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqInfoMode {
    Auto,
    Fasta,
    Fastq,
}

/// Serialise one record, matching `sequenceInfo2FileEntry`.
pub fn write_record(out: &mut Vec<u8>, s: &SequenceInfo, mode: SeqInfoMode) {
    out.push(if mode == SeqInfoMode::Auto { s.sep } else { b'>' });
    out.extend_from_slice(&s.name);
    if !s.comment.is_empty() {
        out.push(b' ');
        out.extend_from_slice(&s.comment);
    }
    out.push(b'\n');
    out.extend_from_slice(&s.seq);
    out.push(b'\n');
    if mode == SeqInfoMode::Auto && s.sep == b'@' {
        out.extend_from_slice(b"+\n");
        out.extend_from_slice(&s.qual);
        out.push(b'\n');
    }
}

/// Average of the two neighbouring quality characters, from `getAvgQual`.
///
/// Note this averages the raw ASCII characters, not decoded Phred scores; since
/// Phred+33 is affine that is the same ranking, and it is what the C++ writes back
/// into the quality string.
#[inline]
pub fn get_avg_qual(qual: &[u8], pos: usize) -> u8 {
    if qual.is_empty() {
        return 33; // "error value, phred 33"
    }
    let mut qval: i32 = 0;
    let mut count: i32 = 0;
    if pos > 0 {
        qval += qual[pos - 1] as i32;
        count += 1;
    }
    if pos + 1 < qual.len() {
        qval += qual[pos + 1] as i32;
        count += 1;
    }
    if count == 0 {
        // Reachable only for a one-character quality string, where the C++ divides
        // by zero. Fall back to the same value it uses for "no quality at all".
        return 33;
    }
    (qval / count) as u8
}

// ---------------------------------------------------------------------------
// Reader
// ---------------------------------------------------------------------------

const DEFAULT_CAPACITY: usize = 8 << 20;

/// Buffered FASTA/FASTQ reader that yields batches of records.
pub struct FastxReader<R: Read> {
    inner: R,
    buf: Vec<u8>,
    /// Live region of `buf`.
    start: usize,
    end: usize,
    eof: bool,
}

impl FastxReader<Box<dyn Read + Send>> {
    /// Open a path, transparently decompressing `.gz`.
    pub fn open(path: &str) -> std::io::Result<Self> {
        let inner: Box<dyn Read + Send> = if path == "stdin" {
            Box::new(std::io::stdin())
        } else {
            let f = std::fs::File::open(path)?;
            if path.ends_with(".gz") {
                Box::new(flate2::read::MultiGzDecoder::new(std::io::BufReader::with_capacity(
                    1 << 20,
                    f,
                )))
            } else if path.ends_with(".bz2") {
                Box::new(bzip2::read::MultiBzDecoder::new(std::io::BufReader::with_capacity(
                    1 << 20,
                    f,
                )))
            } else {
                Box::new(f)
            }
        };
        Ok(FastxReader::new(inner))
    }
}

impl<R: Read> FastxReader<R> {
    pub fn new(inner: R) -> Self {
        FastxReader {
            inner,
            buf: vec![0u8; DEFAULT_CAPACITY],
            start: 0,
            end: 0,
            eof: false,
        }
    }

    /// Pull more bytes in, compacting the live region to the front first. Grows the
    /// buffer when a single record does not fit.
    fn refill(&mut self) -> std::io::Result<usize> {
        if self.start > 0 {
            self.buf.copy_within(self.start..self.end, 0);
            self.end -= self.start;
            self.start = 0;
        }
        if self.end == self.buf.len() {
            let n = self.buf.len();
            self.buf.resize(n * 2, 0);
        }
        let n = self.inner.read(&mut self.buf[self.end..])?;
        if n == 0 {
            self.eof = true;
        }
        self.end += n;
        Ok(n)
    }

    /// Fill `out` with up to `max` records. Existing `SequenceInfo`s are reused so
    /// a long run does not reallocate. Returns how many records were produced.
    pub fn read_batch(
        &mut self,
        out: &mut Vec<SequenceInfo>,
        max: usize,
    ) -> std::io::Result<usize> {
        let mut n = 0usize;
        while n < max {
            // Make sure the buffer holds at least one complete record.
            loop {
                if let Some(len) = record_len(&self.buf[self.start..self.end], self.eof) {
                    if out.len() <= n {
                        out.push(SequenceInfo::default());
                    }
                    let rec = &self.buf[self.start..self.start + len];
                    parse_record(rec, &mut out[n]);
                    self.start += len;
                    n += 1;
                    break;
                }
                if self.eof {
                    return Ok(n);
                }
                self.refill()?;
            }
        }
        Ok(n)
    }
}

/// Length of the first complete record in `buf`, or `None` if more input is needed.
///
/// A record ends where the next one begins: at a `>` or `@` that starts a line and
/// is not inside a quality block. FASTQ quality lines can themselves begin with
/// `@`, so the scan tracks how much quality has been seen rather than trusting the
/// first-character heuristic.
fn record_len(buf: &[u8], eof: bool) -> Option<usize> {
    if buf.is_empty() {
        return None;
    }
    let sep = buf[0];
    if sep != b'>' && sep != b'@' {
        // Leading junk: kseq skips to the first header character.
        let pos = buf.iter().position(|&c| c == b'>' || c == b'@')?;
        return record_len(&buf[pos..], eof).map(|l| l + pos);
    }

    let mut pos = line_end(buf, 0)?; // past the header line
    if sep == b'>' {
        // FASTA: sequence lines until the next '>' at the start of a line.
        loop {
            if pos >= buf.len() {
                return if eof { Some(buf.len()) } else { None };
            }
            if buf[pos] == b'>' {
                return Some(pos);
            }
            pos = match line_end(buf, pos) {
                Some(p) => p,
                None => return if eof { Some(buf.len()) } else { None },
            };
        }
    }

    // FASTQ: sequence lines until a '+' line, then as many quality lines as it
    // takes to match the sequence length.
    let mut seq_len = 0usize;
    loop {
        if pos >= buf.len() {
            return None;
        }
        if buf[pos] == b'+' {
            break;
        }
        let e = line_end(buf, pos)?;
        seq_len += trimmed_len(&buf[pos..e]);
        pos = e;
    }
    pos = line_end(buf, pos)?; // past the '+' line
    let mut qual_len = 0usize;
    while qual_len < seq_len {
        if pos >= buf.len() {
            return if eof && qual_len > 0 { Some(buf.len()) } else { None };
        }
        let e = match line_end(buf, pos) {
            Some(p) => p,
            None => return if eof { Some(buf.len()) } else { None },
        };
        qual_len += trimmed_len(&buf[pos..e]);
        pos = e;
    }
    Some(pos)
}

/// Index just past the newline ending the line that starts at `from`.
#[inline]
fn line_end(buf: &[u8], from: usize) -> Option<usize> {
    memchr::memchr(b'\n', &buf[from..]).map(|p| from + p + 1)
}

/// Length of a raw line slice (which still holds its `\n`) once `\r\n` is trimmed.
#[inline]
fn trimmed_len(line: &[u8]) -> usize {
    let mut l = line.len();
    if l > 0 && line[l - 1] == b'\n' {
        l -= 1;
    }
    if l > 0 && line[l - 1] == b'\r' {
        l -= 1;
    }
    l
}

#[inline]
fn trim(line: &[u8]) -> &[u8] {
    &line[..trimmed_len(line)]
}

/// Split one complete record into its fields.
fn parse_record(rec: &[u8], out: &mut SequenceInfo) {
    out.clear();
    out.sep = rec[0];

    let header_end = memchr::memchr(b'\n', rec).map(|p| p + 1).unwrap_or(rec.len());
    let header = trim(&rec[1..header_end]);
    // kseq's KS_SEP_SPACE stops at any whitespace.
    match header.iter().position(|c| c.is_ascii_whitespace()) {
        Some(i) => {
            out.name.extend_from_slice(&header[..i]);
            // The remainder after the single separator character is the comment.
            let rest = &header[i + 1..];
            out.comment.extend_from_slice(rest);
        }
        None => out.name.extend_from_slice(header),
    }

    let mut pos = header_end;
    if out.sep == b'>' {
        while pos < rec.len() {
            let e = memchr::memchr(b'\n', &rec[pos..]).map(|p| pos + p + 1).unwrap_or(rec.len());
            out.seq.extend_from_slice(trim(&rec[pos..e]));
            pos = e;
        }
        return;
    }

    while pos < rec.len() && rec[pos] != b'+' {
        let e = memchr::memchr(b'\n', &rec[pos..]).map(|p| pos + p + 1).unwrap_or(rec.len());
        out.seq.extend_from_slice(trim(&rec[pos..e]));
        pos = e;
    }
    if pos < rec.len() {
        // Skip the '+' line.
        pos = memchr::memchr(b'\n', &rec[pos..]).map(|p| pos + p + 1).unwrap_or(rec.len());
    }
    while pos < rec.len() && out.qual.len() < out.seq.len() {
        let e = memchr::memchr(b'\n', &rec[pos..]).map(|p| pos + p + 1).unwrap_or(rec.len());
        out.qual.extend_from_slice(trim(&rec[pos..e]));
        pos = e;
    }
}

/// `getSeqMode`: FASTA or FASTQ, decided from the first record.
pub fn get_seq_mode(path: &str) -> std::io::Result<SeqInfoMode> {
    let mut r = FastxReader::open(path)?;
    let mut batch = Vec::new();
    let n = r.read_batch(&mut batch, 1)?;
    if n == 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("ERROR: Can not read from {path}"),
        ));
    }
    let s = &batch[0];
    if s.qual.is_empty() {
        Ok(SeqInfoMode::Fasta)
    } else if s.qual.len() == s.seq.len() {
        Ok(SeqInfoMode::Fastq)
    } else {
        Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("ERROR: Invalid first entry found in {path}"),
        ))
    }
}

/// Large buffered writer; output is appended in input order by the caller.
pub struct RecordWriter {
    inner: std::io::BufWriter<std::fs::File>,
}

impl RecordWriter {
    pub fn create(path: &str) -> std::io::Result<Self> {
        let f = std::fs::File::create(path)?;
        Ok(RecordWriter { inner: std::io::BufWriter::with_capacity(4 << 20, f) })
    }
    pub fn write_all(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        self.inner.write_all(bytes)
    }
    pub fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

/// Read every record in a file. Convenience for tests and small inputs.
pub fn read_all(path: &str) -> std::io::Result<Vec<SequenceInfo>> {
    let mut r = FastxReader::open(path)?;
    let mut all = Vec::new();
    loop {
        let mut batch = Vec::new();
        let n = r.read_batch(&mut batch, 4096)?;
        if n == 0 {
            break;
        }
        batch.truncate(n);
        all.extend(batch);
    }
    Ok(all)
}

#[allow(dead_code)]
fn unused_bufread_marker<T: BufRead>(_: T) {}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_str(s: &str) -> Vec<SequenceInfo> {
        let mut r = FastxReader::new(std::io::Cursor::new(s.as_bytes().to_vec()));
        let mut out = Vec::new();
        let mut all = Vec::new();
        loop {
            let n = r.read_batch(&mut out, 8).unwrap();
            if n == 0 {
                break;
            }
            all.extend_from_slice(&out[..n]);
        }
        all
    }

    #[test]
    fn parses_a_plain_fastq_record() {
        let recs = parse_str("@r1\nACGT\n+\nIIII\n");
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].name, b"r1");
        assert_eq!(recs[0].seq, b"ACGT");
        assert_eq!(recs[0].qual, b"IIII");
        assert_eq!(recs[0].sep, b'@');
        assert!(recs[0].comment.is_empty());
    }

    #[test]
    fn splits_name_and_comment_at_the_first_whitespace() {
        let recs = parse_str("@r1 some comment here\nACGT\n+\nIIII\n");
        assert_eq!(recs[0].name, b"r1");
        assert_eq!(recs[0].comment, b"some comment here");
    }

    #[test]
    fn parses_multiline_fasta() {
        let recs = parse_str(">s1 desc\nACGT\nTTTT\n>s2\nGG\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].name, b"s1");
        assert_eq!(recs[0].comment, b"desc");
        assert_eq!(recs[0].seq, b"ACGTTTTT");
        assert_eq!(recs[0].sep, b'>');
        assert!(recs[0].qual.is_empty());
        assert_eq!(recs[1].name, b"s2");
        assert_eq!(recs[1].seq, b"GG");
    }

    #[test]
    fn handles_quality_lines_that_start_with_at() {
        // '@' is a legal Phred+33 quality character, so a naive line-oriented
        // parser would split this record in two.
        let recs = parse_str("@r1\nACGTACGT\n+\n@@@@@@@@\n@r2\nTTTT\n+\nIIII\n");
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].qual, b"@@@@@@@@");
        assert_eq!(recs[1].name, b"r2");
        assert_eq!(recs[1].seq, b"TTTT");
    }

    #[test]
    fn strips_carriage_returns() {
        let recs = parse_str("@r1 c\r\nACGT\r\n+\r\nIIII\r\n");
        assert_eq!(recs[0].name, b"r1");
        assert_eq!(recs[0].comment, b"c");
        assert_eq!(recs[0].seq, b"ACGT");
        assert_eq!(recs[0].qual, b"IIII");
    }

    #[test]
    fn handles_a_record_spanning_a_buffer_refill() {
        // Force many refills by using a tiny reader that dribbles bytes out.
        struct Dribble(Vec<u8>, usize);
        impl std::io::Read for Dribble {
            fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
                if self.1 >= self.0.len() {
                    return Ok(0);
                }
                let n = 3.min(self.0.len() - self.1).min(buf.len());
                buf[..n].copy_from_slice(&self.0[self.1..self.1 + n]);
                self.1 += n;
                Ok(n)
            }
        }
        let text = "@r1\nACGTACGTAC\n+\nIIIIIIIIII\n@r2\nTTTTTTTTTT\n+\nJJJJJJJJJJ\n";
        let mut r = FastxReader::new(Dribble(text.as_bytes().to_vec(), 0));
        let mut out = Vec::new();
        let mut all = Vec::new();
        loop {
            let n = r.read_batch(&mut out, 1).unwrap();
            if n == 0 {
                break;
            }
            all.extend_from_slice(&out[..n]);
        }
        assert_eq!(all.len(), 2);
        assert_eq!(all[0].seq, b"ACGTACGTAC");
        assert_eq!(all[1].qual, b"JJJJJJJJJJ");
    }

    #[test]
    fn round_trips_through_write_record() {
        let recs = parse_str("@r1 c\nACGT\n+\nIIII\n>s2\nGGGG\n");
        let mut out = Vec::new();
        write_record(&mut out, &recs[0], SeqInfoMode::Auto);
        assert_eq!(out, b"@r1 c\nACGT\n+\nIIII\n");
        out.clear();
        write_record(&mut out, &recs[1], SeqInfoMode::Auto);
        assert_eq!(out, b">s2\nGGGG\n");
    }

    #[test]
    fn avg_qual_averages_the_two_neighbours() {
        let q = b"ABCDE";
        // pos 2: neighbours 'B'(66) and 'D'(68) -> 67 = 'C'
        assert_eq!(get_avg_qual(q, 2), b'C');
        // pos 0: only the right neighbour counts
        assert_eq!(get_avg_qual(q, 0), b'B');
        // last position: only the left neighbour counts
        assert_eq!(get_avg_qual(q, 4), b'D');
        // empty quality -> phred 33
        assert_eq!(get_avg_qual(b"", 0), 33);
    }
}
