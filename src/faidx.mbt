// pyfaidx-style FASTA fast index access.
//
// Mirrors pyfaidx (https://github.com/mdshw5/pyfaidx) which provides
// fast random access to FASTA files using a .fai index file.
//
// .fai index format (tab-separated, one line per sequence):
//   1. name      - sequence name
//   2. length    - total number of bases
//   3. offset    - byte offset of the first base in the file
//   4. line_len  - number of bases per line
//   5. line_blen - number of bytes per line (including line endings)

///|
/// Error raised by FastaIndex operations.
pub suberror FastaIndexError {
  FastaIndexError(String)
}

///|
/// A single index entry in the .fai file.
pub struct FaiRecord {
  name : String
  length : Int
  offset : Int
  line_len : Int
  line_blen : Int
} derive(Debug)

///|
/// Create a new FaiRecord.
pub fn FaiRecord::new(
  name : String,
  length : Int,
  offset : Int,
  line_len : Int,
  line_blen : Int,
) -> FaiRecord {
  FaiRecord::{ name, length, offset, line_len, line_blen }
}

///|
/// FastaIndex holds the parsed .fai index for fast random access.
pub struct FastaIndex {
  records : Array[FaiRecord]
  index : Map[String, Int]
} derive(Debug)

///|
/// Parse a .fai index file and return a FastaIndex.
pub fn parse_fai(content : String) -> FastaIndex raise FastaIndexError {
  let records : Array[FaiRecord] = Array::new()
  let index : Map[String, Int] = Map([], capacity=32)
  let lines = content.split("\n").to_array()
  let mut line_num = 0
  for line in lines {
    line_num = line_num + 1
    let mut trimmed = line
    if trimmed.length() > 0 &&
      trimmed.unsafe_get(trimmed.length() - 1) == '\r'.to_int().to_uint16() {
      trimmed = trimmed[0:trimmed.length() - 1]
    }
    if trimmed.length() == 0 {
      continue
    }
    let parts = trimmed.split("\t").to_array()
    if parts.length() < 5 {
      raise FastaIndexError(
        "Invalid .fai line " +
        line_num.to_string() +
        ": expected 5 fields, got " +
        parts.length().to_string(),
      )
    }
    let name = parts[0].to_owned()
    let length = parse_fai_int(parts[1].to_owned())
    let offset = parse_fai_int(parts[2].to_owned())
    let line_len = parse_fai_int(parts[3].to_owned())
    let line_blen = parse_fai_int(parts[4].to_owned())
    if length < 0 {
      raise FastaIndexError(
        "Invalid length for " + name + ": " + length.to_string(),
      )
    }
    if offset < 0 {
      raise FastaIndexError(
        "Invalid offset for " + name + ": " + offset.to_string(),
      )
    }
    if line_len <= 0 {
      raise FastaIndexError(
        "Invalid line_len for " + name + ": " + line_len.to_string(),
      )
    }
    if line_blen <= 0 {
      raise FastaIndexError(
        "Invalid line_blen for " + name + ": " + line_blen.to_string(),
      )
    }
    let idx = records.length()
    records.push(FaiRecord::new(name, length, offset, line_len, line_blen))
    index[name] = idx
  }
  FastaIndex::{ records, index }
}

///|
/// Get the number of sequences in the index.
pub fn FastaIndex::num_sequences(self : FastaIndex) -> Int {
  self.records.length()
}

///|
/// Get the names of all sequences in the index.
pub fn FastaIndex::get_names(self : FastaIndex) -> Array[String] {
  let names : Array[String] = Array::new()
  for rec in self.records {
    names.push(rec.name)
  }
  names
}

///|
/// Look up a FaiRecord by name.
pub fn FastaIndex::get_record(self : FastaIndex, name : String) -> FaiRecord? {
  match self.index.get(name) {
    Some(idx) => Some(self.records[idx])
    None => None
  }
}

///|
/// Check if a sequence exists in the index.
pub fn FastaIndex::contains(self : FastaIndex, name : String) -> Bool {
  self.index.contains(name)
}

///|
/// Build a .fai index from FASTA content.
///
/// This scans the FASTA content once and computes the offset, line_len,
/// and line_blen for each sequence, matching pyfaidx's index generation.
pub fn build_fai(content : String) -> FastaIndex {
  let records : Array[FaiRecord] = Array::new()
  let index : Map[String, Int] = Map([], capacity=32)
  let mut pos = 0
  let n = content.length()
  let gt : UInt16 = '>'.to_int().to_uint16()
  let nl : UInt16 = '\n'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()

  while pos < n {
    if content.unsafe_get(pos) != gt {
      pos = pos + 1
      continue
    }
    // Found header line '>'
    let header_start = pos
    // Find end of header line
    let mut header_end = pos + 1
    while header_end < n && content.unsafe_get(header_end) != nl {
      header_end = header_end + 1
    }
    // Parse name (first whitespace-delimited token after '>')
    let title_start = pos + 1
    let mut name_end = title_start
    while name_end < header_end {
      let c = content.unsafe_get(name_end)
      if c == ' '.to_int().to_uint16() || c == '\t'.to_int().to_uint16() {
        break
      }
      name_end = name_end + 1
    }
    let name = content[title_start:name_end].to_owned()
    // Sequence starts after header line
    let mut seq_start = header_end
    if seq_start < n && content.unsafe_get(seq_start) == nl {
      seq_start = seq_start + 1
    }
    let seq_offset = seq_start
    // Scan to find line_len and line_blen, and total length
    let mut line_len = 0
    let mut line_blen = 0
    let mut total_len = 0
    let mut found_line_len = false
    let mut scan_pos = seq_start
    while scan_pos < n && content.unsafe_get(scan_pos) != gt {
      // Find end of this line
      let mut line_end = scan_pos
      while line_end < n &&
            content.unsafe_get(line_end) != nl &&
            content.unsafe_get(line_end) != gt {
        line_end = line_end + 1
      }
      // Count bases (skip whitespace)
      let mut bases = 0
      for i = scan_pos; i < line_end; i = i + 1 {
        let c = content.unsafe_get(i)
        if c != ' '.to_int().to_uint16() &&
          c != '\t'.to_int().to_uint16() &&
          c != cr {
          bases = bases + 1
        }
      }
      if bases > 0 && !found_line_len {
        line_len = bases
        // line_blen = bytes per line (including line ending)
        let mut byte_count = line_end - scan_pos
        if line_end < n && content.unsafe_get(line_end) == nl {
          byte_count = byte_count + 1
        }
        line_blen = byte_count
        found_line_len = true
      }
      total_len = total_len + bases
      if line_end >= n || content.unsafe_get(line_end) == gt {
        break
      }
      scan_pos = line_end + 1
    }
    let idx = records.length()
    records.push(
      FaiRecord::new(name, total_len, seq_offset, line_len, line_blen),
    )
    index[name] = idx
    pos = scan_pos
  }
  FastaIndex::{ records, index }
}

///|
/// Write a FastaIndex to .fai format string.
pub fn write_fai(idx : FastaIndex) -> String {
  let buf = StringBuilder::new()
  for rec in idx.records {
    buf.write_string(rec.name)
    buf.write_char('\t')
    buf.write_string(rec.length.to_string())
    buf.write_char('\t')
    buf.write_string(rec.offset.to_string())
    buf.write_char('\t')
    buf.write_string(rec.line_len.to_string())
    buf.write_char('\t')
    buf.write_string(rec.line_blen.to_string())
    buf.write_char('\n')
  }
  buf.to_string()
}

///|
/// Fasta provides fast random access to sequences in a FASTA file,
/// using a pre-built index (like pyfaidx.Fasta).
///
/// Example:
///   let fasta = Fasta::new(fasta_content, fai_content)?
///   let seq = fasta.fetch("chr1", 1000, 2000)?
pub struct Fasta {
  content : String
  index : FastaIndex
} derive(Debug)

///|
/// Create a Fasta from FASTA content and .fai index content.
pub fn Fasta::new(
  fasta_content : String,
  fai_content : String,
) -> Fasta raise FastaIndexError {
  let index = parse_fai(fai_content)
  Fasta::{ content: fasta_content, index }
}

///|
/// Create a Fasta from FASTA content by building the index on-the-fly.
pub fn Fasta::from_content(content : String) -> Fasta {
  let index = build_fai(content)
  Fasta::{ content, index }
}

///|
/// Get the full sequence for a given name.
pub fn Fasta::get_seq(self : Fasta, name : String) -> Seq? {
  match self.index.get_record(name) {
    Some(rec) => Some(self.extract_sequence(rec, 0, rec.length))
    None => None
  }
}

///|
/// Fetch a subsequence (0-based, half-open: [start, end)).
/// Mirrors pyfaidx's Fasta[name][start:end].
pub fn Fasta::fetch(
  self : Fasta,
  name : String,
  start : Int,
  end : Int,
) -> Seq? raise FastaIndexError {
  match self.index.get_record(name) {
    Some(rec) => {
      if start < 0 {
        raise FastaIndexError("start must be >= 0, got " + start.to_string())
      }
      if end > rec.length {
        raise FastaIndexError(
          "end (" +
          end.to_string() +
          ") exceeds sequence length (" +
          rec.length.to_string() +
          ")",
        )
      }
      if start > end {
        raise FastaIndexError(
          "start (" + start.to_string() + ") > end (" + end.to_string() + ")",
        )
      }
      Some(self.extract_sequence(rec, start, end))
    }
    None => None
  }
}

///|
/// Get the length of a sequence.
pub fn Fasta::get_length(self : Fasta, name : String) -> Int? {
  match self.index.get_record(name) {
    Some(rec) => Some(rec.length)
    None => None
  }
}

///|
/// Check if a sequence exists.
pub fn Fasta::contains(self : Fasta, name : String) -> Bool {
  self.index.contains(name)
}

///|
/// Get all sequence names.
pub fn Fasta::get_names(self : Fasta) -> Array[String] {
  self.index.get_names()
}

///|
/// Extract a subsequence from the FASTA content using the index record.
///
/// Uses the formula from pyfaidx:
///   byte_offset = offset + (start / line_len) * line_blen + (start % line_len)
///
/// Then reads the required number of bases, skipping line endings.
fn Fasta::extract_sequence(
  self : Fasta,
  rec : FaiRecord,
  start : Int,
  end : Int,
) -> Seq {
  let total_bases = end - start
  if total_bases <= 0 {
    return Seq::new("")
  }
  // Calculate starting byte offset
  let start_line = start / rec.line_len
  let start_col = start % rec.line_len
  let byte_start_raw = rec.offset + start_line * rec.line_blen + start_col
  let content_len = self.content.length()
  let byte_start = if byte_start_raw > content_len {
    content_len
  } else {
    byte_start_raw
  }
  // Calculate ending byte offset (add extra for line endings)
  let end_line = end / rec.line_len
  let end_col = end % rec.line_len
  let byte_end_raw = rec.offset + end_line * rec.line_blen + end_col
  let byte_end = if byte_end_raw > content_len {
    content_len
  } else {
    byte_end_raw
  }
  // Extract raw bytes
  let raw = self.content[byte_start:byte_end].to_owned()
  // Filter out line endings and whitespace
  let buf = StringBuilder::new()
  let nl = '\n'.to_int().to_uint16()
  let cr = '\r'.to_int().to_uint16()
  let sp = ' '.to_int().to_uint16()
  let tab = '\t'.to_int().to_uint16()
  let mut count = 0
  let mut i = 0
  let n = raw.length()
  while i < n && count < total_bases {
    let c = raw.unsafe_get(i)
    if c != nl && c != cr && c != sp && c != tab {
      buf.write_char(c.unsafe_to_char())
      count = count + 1
    }
    i = i + 1
  }
  Seq::new(buf.to_string())
}

///|
/// Parse an integer from a string, raising FastaIndexError on failure.
fn parse_fai_int(s : String) -> Int raise FastaIndexError {
  let n = s.length()
  if n == 0 {
    raise FastaIndexError("Invalid integer: empty string")
  }
  let mut result = 0
  let mut i = 0
  let zero = '0'.to_int()
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 {
      raise FastaIndexError("Invalid integer: " + s)
    }
    result = result * 10 + (c - zero)
    i = i + 1
  }
  result
}
