// Bio.SeqIO.IgIO - Parser and writer for the IntelliGenetics (MASE) sequence
// format, mirroring Biopython's Bio.SeqIO.IgIO module.
//
// The IntelliGenetics (a.k.a. MASE) file format is a multiple-sequence text
// format with the following structure:
//
//   1. Optional file-header lines, each starting with ";;" (two semicolons).
//      These are free-format and are ignored by the parser.
//   2. One or more records. Each record consists of:
//      a. One or more commentary lines starting with ";" (single semicolon).
//         The text following the ";" on each line is collected (stripped of
//         leading/trailing whitespace) and joined with newlines into a single
//         comment string.
//      b. A title line (the record's identifier / name).
//      c. One or more sequence lines. Whitespace within sequence lines is
//         removed. The sequence may optionally end with the digit "1" as a
//         terminator, which is stripped on parse.
//
// Reference:
//   https://biopython.org/docs/dev/api/Bio.SeqIO.IgIO.html

// ============================================================================
// IgRecord - a single IntelliGenetics record
// ============================================================================

///|
/// A single IntelliGenetics record: a title (id), the comment block, and the
/// sequence string.
pub struct IgRecord {
  title : String
  comment : String
  sequence : String
}

///|
/// Construct an `IgRecord`.
pub fn IgRecord::new(
  title : String,
  comment : String,
  sequence : String
) -> IgRecord {
  IgRecord::{ title, comment, sequence }
}

///|
/// The record title (also used as id and name).
pub fn IgRecord::title(self : IgRecord) -> String {
  self.title
}

///|
/// The record's free-form comment block (possibly empty).
pub fn IgRecord::comment(self : IgRecord) -> String {
  self.comment
}

///|
/// The record's sequence string (no whitespace, no terminator).
pub fn IgRecord::sequence(self : IgRecord) -> String {
  self.sequence
}

///|
/// The record id: the title string.
pub fn IgRecord::id(self : IgRecord) -> String {
  self.title
}

///|
/// The record name: the title string.
pub fn IgRecord::name(self : IgRecord) -> String {
  self.title
}

///|
/// The record description: the title string.
pub fn IgRecord::description(self : IgRecord) -> String {
  self.title
}

///|
/// Sequence length in nucleotides / amino acids.
pub fn IgRecord::length(self : IgRecord) -> Int {
  self.sequence.length()
}

///|
/// Human-readable representation of this record.
pub fn IgRecord::to_string(self : IgRecord) -> String {
  let buf = StringBuilder::new()
  buf.write_string("IgRecord(id=")
  buf.write_string(self.title)
  buf.write_string(", length=")
  buf.write_string(self.sequence.length().to_string())
  if self.comment.length() > 0 {
    buf.write_string(", comment=")
    buf.write_string(self.comment.length().to_string())
    buf.write_string(" chars")
  }
  buf.write_string(")")
  buf.to_string()
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Character code helpers for fast line classification.
fn ig_char_code(ch : Char) -> UInt16 {
  ch.to_int().to_uint16()
}

///|
/// True if a line starts with the given character.
fn ig_starts_with(line : String, ch : UInt16) -> Bool {
  if line.length() == 0 {
    return false
  }
  line.unsafe_get(0) == ch
}

///|
/// True if a line starts with the two-character prefix `c1` followed by `c2`.
fn ig_starts_with2(line : String, c1 : UInt16, c2 : UInt16) -> Bool {
  let n = line.length()
  if n < 2 {
    return false
  }
  line.unsafe_get(0) == c1 && line.unsafe_get(1) == c2
}

///|
/// Strip a trailing '\r' (CR) from a line (for CRLF input). Returns the
/// possibly-shortened string.
fn ig_strip_cr(line : String) -> String {
  let n = line.length()
  if n > 0 && line.unsafe_get(n - 1) == ig_char_code('\r') {
    line[0 : n - 1].to_owned()
  } else {
    line
  }
}

///|
/// Strip leading and trailing whitespace from a string.
fn ig_strip(s : String) -> String {
  let n = s.length()
  let mut start = 0
  let mut end = n
  let sp = ig_char_code(' ')
  let tab = ig_char_code('\t')
  let cr = ig_char_code('\r')
  let nl = ig_char_code('\n')
  while start < end {
    let c = s.unsafe_get(start)
    if c == sp || c == tab || c == cr || c == nl {
      start = start + 1
    } else {
      break
    }
  }
  while end > start {
    let c = s.unsafe_get(end - 1)
    if c == sp || c == tab || c == cr || c == nl {
      end = end - 1
    } else {
      break
    }
  }
  if start == 0 && end == n {
    s
  } else {
    s[start : end].to_owned()
  }
}

///|
/// Remove all whitespace from a string (used for sequence lines).
fn ig_remove_whitespace(s : String) -> String {
  let n = s.length()
  let buf = StringBuilder::new()
  let sp = ig_char_code(' ')
  let tab = ig_char_code('\t')
  let cr = ig_char_code('\r')
  let nl = ig_char_code('\n')
  for i = 0; i < n; i = i + 1 {
    let c = s.unsafe_get(i)
    if c != sp && c != tab && c != cr && c != nl {
      buf.write_char(c.unsafe_to_char())
    }
  }
  buf.to_string()
}

///|
/// Parse IntelliGenetics file content into an array of `IgRecord` objects.
///
/// - Lines starting with ";;" at the top of the file are file-header lines
///   and are skipped.
/// - Each record begins with one or more ";" commentary lines. The comment
///   text (after the ";") is collected, stripped, and joined with newlines.
/// - The first non-comment line is the record title.
/// - Subsequent lines until the next ";" (or EOF) are sequence lines.
///   Whitespace is removed. A trailing "1" terminator is stripped.
pub fn ig_parse(content : String) -> Array[IgRecord] {
  let records = Array::new()
  let raw_lines = content.split("\n").to_array()
  // Materialize lines as owned strings up-front for simpler downstream use.
  let lines : Array[String] = Array::new()
  for i in 0..<raw_lines.length() {
    lines.push(raw_lines[i].to_owned())
  }
  let n_lines = lines.length()
  let semi = ig_char_code(';')
  let one = ig_char_code('1')
  let mut idx = 0
  // Skip leading file-header lines (";;").
  while idx < n_lines {
    let line = ig_strip_cr(lines[idx])
    if ig_starts_with2(line, semi, semi) {
      idx = idx + 1
    } else {
      break
    }
  }
  // Parse records.
  while idx < n_lines {
    let line = ig_strip_cr(lines[idx])
    // A record must start with a ";" comment line. If not, skip the line
    // (tolerant of stray blank lines between records).
    if line.length() == 0 {
      idx = idx + 1
      continue
    }
    if !ig_starts_with(line, semi) {
      // Not a comment line; skip forward to the next ";". This is more
      // permissive than Biopython, which raises, but it keeps the parser
      // robust for hand-edited files.
      idx = idx + 1
      continue
    }
    // Collect comment lines.
    let comment_parts : Array[String] = Array::new()
    while idx < n_lines {
      let line2 = ig_strip_cr(lines[idx])
      if ig_starts_with(line2, semi) {
        // Strip the leading ";" (and one optional leading space) and trim.
        let body = if line2.length() >= 2 &&
          line2.unsafe_get(1) == ig_char_code(' ') {
          line2[2:].to_owned()
        } else {
          line2[1:].to_owned()
        }
        comment_parts.push(ig_strip(body))
        idx = idx + 1
      } else {
        break
      }
    }
    // The next non-empty line is the title.
    let mut title = ""
    while idx < n_lines {
      let line2 = ig_strip_cr(lines[idx])
      if line2.length() == 0 {
        idx = idx + 1
      } else {
        title = ig_strip(line2)
        idx = idx + 1
        break
      }
    }
    if title.length() == 0 {
      // No title found; stop parsing.
      break
    }
    // Collect sequence lines until the next ";" or EOF.
    let seq_buf = StringBuilder::new()
    let mut hit_next_record = false
    while idx < n_lines {
      let line2 = ig_strip_cr(lines[idx])
      if ig_starts_with(line2, semi) {
        hit_next_record = true
        break
      }
      if line2.length() > 0 {
        let cleaned = ig_remove_whitespace(line2)
        seq_buf.write_string(cleaned)
      }
      idx = idx + 1
    }
    let mut seq_str = seq_buf.to_string()
    // Strip an optional trailing "1" terminator.
    let sn = seq_str.length()
    if sn > 0 && seq_str.unsafe_get(sn - 1) == one {
      seq_str = seq_str[0 : sn - 1].to_owned()
    }
    let comment = ig_join_lines(comment_parts)
    records.push(IgRecord::new(title, comment, seq_str))
    if !hit_next_record {
      // Reached EOF.
      break
    }
  }
  records
}

///|
/// Join an array of comment lines with newlines (skipping trailing empties).
fn ig_join_lines(parts : Array[String]) -> String {
  let n = parts.length()
  // Trim trailing empty strings.
  let mut end = n
  while end > 0 && parts[end - 1].length() == 0 {
    end = end - 1
  }
  let buf = StringBuilder::new()
  for i = 0; i < end; i = i + 1 {
    if i > 0 {
      buf.write_char('\n')
    }
    buf.write_string(parts[i])
  }
  buf.to_string()
}

// ============================================================================
// Writing
// ============================================================================

///|
/// Default line width for wrapping IntelliGenetics sequence output.
let ig_line_width : Int = 60

///|
/// Render an `IgRecord` to IntelliGenetics text (without a trailing blank line).
///
/// The output is:
///   ;comment line 1
///   ;comment line 2
///   title
///   <wrapped sequence>
///   1
/// The trailing "1" terminator is emitted to match the canonical format.
pub fn ig_record_to_string(record : IgRecord) -> String {
  ig_record_to_string_width(record, ig_line_width)
}

///|
/// Render an `IgRecord` with a custom sequence line width. A width <= 0
/// disables wrapping (the sequence is emitted on a single line).
pub fn ig_record_to_string_width(
  record : IgRecord,
  width : Int
) -> String {
  let buf = StringBuilder::new()
  // Comment lines.
  if record.comment.length() > 0 {
    let parts = record.comment.split("\n").to_array()
    for i in 0..<parts.length() {
      buf.write_char(';')
      let p = parts[i].to_owned()
      if p.length() > 0 {
        buf.write_char(' ')
        buf.write_string(p)
      }
      buf.write_char('\n')
    }
  } else {
    // Emit a minimal ";" comment line so the record is identifiable as an
    // IntelliGenetics record (the format requires at least one comment line
    // before the title).
    buf.write_string(";\n")
  }
  // Title.
  buf.write_string(record.title)
  buf.write_char('\n')
  // Sequence, wrapped to `width` characters per line.
  let seq = record.sequence
  let n = seq.length()
  if width <= 0 || n <= width {
    buf.write_string(seq)
    buf.write_char('\n')
  } else {
    let mut i = 0
    while i < n {
      let end = if i + width < n { i + width } else { n }
      buf.write_string(seq[i : end].to_owned())
      buf.write_char('\n')
      i = end
    }
  }
  // Trailing terminator.
  buf.write_char('1')
  buf.write_char('\n')
  buf.to_string()
}

///|
/// Render an array of `IgRecord` objects as a full IntelliGenetics file. An
/// optional file header (the `;;` lines) is prepended when `file_header` is
/// non-empty; each header line is emitted verbatim with a leading ";;" if it
/// does not already start with one.
pub fn ig_records_to_string(
  records : Array[IgRecord],
  file_header : Array[String]?
) -> String {
  let buf = StringBuilder::new()
  // Optional file header.
  match file_header {
    Some(hs) =>
      for h in hs {
        let s = h
        if s.length() >= 2 &&
          s.unsafe_get(0) == ig_char_code(';') &&
          s.unsafe_get(1) == ig_char_code(';') {
          buf.write_string(s)
        } else {
          buf.write_string(";;")
          buf.write_string(s)
        }
        buf.write_char('\n')
      }
    None => ()
  }
  for i in 0..<records.length() {
    buf.write_string(ig_record_to_string(records[i]))
  }
  buf.to_string()
}

///|
/// Convenience: serialize a single record to a file-format string (with the
/// standard 60-char wrapping and a trailing terminator).
pub fn ig_to_string(record : IgRecord) -> String {
  ig_record_to_string(record)
}

// ============================================================================
// SeqRecord interop
// ============================================================================

///|
/// Convert an `IgRecord` to a `SeqRecord`. The record's title is used as the
/// id, name, and description. The comment is not preserved on SeqRecord
/// (which lacks a comment field); use `ig_from_record` to round-trip.
pub fn ig_to_seq_record(record : IgRecord) -> SeqRecord {
  SeqRecord::new(
    Seq::new(record.sequence),
    id=record.title,
    name=record.title,
    description=record.title,
  )
}

///|
/// Convert an array of `IgRecord` objects to an array of `SeqRecord` objects.
pub fn ig_to_seq_records(records : Array[IgRecord]) -> Array[SeqRecord] {
  let out = Array::new()
  for i in 0..<records.length() {
    out.push(ig_to_seq_record(records[i]))
  }
  out
}

///|
/// Convert a `SeqRecord` to an `IgRecord`. The record's id is used as the
/// title; the description is used as the comment (if non-empty and different
/// from the id), otherwise the comment is empty.
pub fn ig_from_seq_record(record : SeqRecord) -> IgRecord {
  let title = record.id
  let seq_str = record.seq.to_string()
  let comment = if record.description.length() > 0 &&
    record.description != title {
    record.description
  } else {
    ""
  }
  IgRecord::new(title, comment, seq_str)
}

///|
/// Convert an array of `SeqRecord` objects to an array of `IgRecord` objects.
pub fn ig_from_seq_records(records : Array[SeqRecord]) -> Array[IgRecord] {
  let out = Array::new()
  for i in 0..<records.length() {
    out.push(ig_from_seq_record(records[i]))
  }
  out
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// A small sample IntelliGenetics file with two records, used for demos and
/// tests. The sequences are short DNA fragments.
pub fn ig_sample_text() -> String {
  ";;\n;; Sample IntelliGenetics file\n;\nA_U455\nATGGCAGCTGTGATGAAGCAGAGACGGGTAAGAGCTC\n1\n;\nB_HXB2R\nATGGGTGCTGTGATGAAGCAGAGACGGGTAAGAGCTC\n1\n"
}

///|
/// Parse and return the sample IntelliGenetics records.
pub fn ig_sample_records() -> Array[IgRecord] {
  ig_parse(ig_sample_text())
}

///|
/// A multi-record IntelliGenetics file header (the `;;` lines) for use with
/// `ig_records_to_string`.
pub fn ig_sample_file_header() -> Array[String] {
  [";;", ";; Sample IntelliGenetics file"]
}

///|
/// Build a sample `IgRecord` for testing.
pub fn ig_sample_record() -> IgRecord {
  IgRecord::new(
    "A_U455",
    "HIV-1 group M subtype A",
    "ATGGCAGCTGTGATGAAGCAGAGACGGGTAAGAGCTC",
  )
}

///|
/// Build a second sample `IgRecord` for testing.
pub fn ig_sample_record_2() -> IgRecord {
  IgRecord::new(
    "B_HXB2R",
    "HIV-1 group M subtype B",
    "ATGGGTGCTGTGATGAAGCAGAGACGGGTAAGAGCTC",
  )
}
