///|
/// Bio.SeqIO.index - random access to sequence files by record key.
///
/// Port of Biopython's Bio.SeqIO.index, which scans a sequence file once,
/// records the (key, start offset, length) of each record, and then gives
/// dictionary-style random access to individual records without parsing
/// the whole file. This module is string-backed like the rest of BioSeqs:
/// `seqio_index(text, format)` maps keys to offsets in `text`, and
/// `index.get_raw(text, key)` returns the raw record substring while
/// `index.get(text, key)` parses it into a SeqRecord.
///
/// Supported formats: "fasta", "fastq"/"fastq-sanger"/"fastq-illumina"/
/// "fastq-solexa", "genbank"/"gb", "embl", "tab".

///|
/// One index entry: record key, start offset in the file text, and the
/// raw record length.
pub struct SeqIndexEntry {
  key : String
  offset : Int
  length : Int
} derive(Debug)
pub extend SeqIndexEntry with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn SeqIndexEntry::new(
  key : String,
  offset : Int,
  length : Int,
) -> SeqIndexEntry {
  SeqIndexEntry::{ key, offset, length, }
}

///|
/// Random-access index over a sequence file.
pub struct SeqIndex {
  format : String
  entries : Array[SeqIndexEntry]
  positions : Map[String, Int]
} derive(Debug)
pub extend SeqIndex with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Number of indexed records.
pub fn SeqIndex::length(self : SeqIndex) -> Int {
  self.entries.length()
}

///|
/// All record keys, in file order.
pub fn SeqIndex::keys(self : SeqIndex) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let mut i = 0
  while i < self.entries.length() {
    result.push(self.entries[i].key)
    i = i + 1
  }
  result
}

///|
/// Whether the index contains a record with the given key.
pub fn SeqIndex::contains(self : SeqIndex, key : String) -> Bool {
  self.positions.contains(key)
}

///|
/// Return the raw record text for a key.
pub fn SeqIndex::get_raw(
  self : SeqIndex,
  data : String,
  key : String,
) -> String raise SeqIOError {
  match self.positions.get(key) {
    Some(idx) => {
      let entry = self.entries[idx]
      data[entry.offset:entry.offset + entry.length].to_owned()
    }
    None => raise SeqIOError("Key " + key + " not found in index")
  }
}

///|
/// Parse and return the record for a key.
pub fn SeqIndex::get(
  self : SeqIndex,
  data : String,
  key : String,
) -> SeqRecord raise SeqIOError {
  let raw = self.get_raw(data, key)
  let records = seqio_parse(raw, self.format)
  if records.length() == 0 {
    raise SeqIOError("No record found for key " + key)
  }
  records[0]
}

///|
/// Index records in the file text for the given format.
pub fn seqio_index(data : String, format : String) -> SeqIndex raise SeqIOError {
  let lines = seq_index_split_lines(data)
  let entries = match format {
    "fasta" => seq_index_scan_fasta(data, lines)
    "fastq" | "fastq-sanger" | "fastq-illumina" | "fastq-solexa" =>
      seq_index_scan_fastq(data, lines)
    "genbank" | "gb" =>
      seq_index_scan_delimiter(data, lines, "LOCUS", "genbank")
    "embl" => seq_index_scan_delimiter(data, lines, "ID   ", "embl")
    "tab" => seq_index_scan_tab(data, lines)
    // qual files use ">" headers like FASTA; the quality numbers follow
    // until the next header (spec 5.4.1.2).
    "qual" => seq_index_scan_fasta(data, lines)
    // Swiss-Prot entries start at the "ID   " line; the key truncates at
    // the first semicolon (shared delimiter scanner, spec 5.4.1.2).
    "swiss" => seq_index_scan_delimiter(data, lines, "ID   ", "swiss")
    "seqxml" => seq_index_scan_seqxml(data, lines)
    _ => raise SeqIOError("Cannot index format: " + format)
  }
  let positions : Map[String, Int] = Map([], capacity=entries.length())
  let mut i = 0
  while i < entries.length() {
    positions.set(entries[i].key, i)
    i = i + 1
  }
  SeqIndex::{ format, entries, positions, }
}

// ============================================================================
// Line splitting with byte offsets
// ============================================================================

///|
struct IndexedLine {
  start : Int
  text : String
}

///|
/// Split text into lines, remembering the start offset of each line.
fn seq_index_split_lines(data : String) -> Array[IndexedLine] {
  let lines : Array[IndexedLine] = Array(capacity=0)
  let n = data.length()
  let mut line_start = 0
  let mut i = 0
  while i < n {
    let ch = data.unsafe_get(i).to_int()
    if ch == '\n'.to_int() {
      lines.push(IndexedLine::{
        start: line_start,
        text: data[line_start:i].to_owned(),
      })
      line_start = i + 1
    }
    i = i + 1
  }
  if line_start < n {
    lines.push(IndexedLine::{
      start: line_start,
      text: data[line_start:n].to_owned(),
    })
  }
  lines
}

///|
/// First whitespace-delimited token of a string, skipping leading spaces.
fn seq_index_first_token(text : String) -> String {
  let n = text.length()
  let mut start = 0
  while start < n {
    let ch = text.unsafe_get(start).to_int()
    if ch != ' '.to_int() && ch != '\t'.to_int() && ch != '\r'.to_int() {
      break
    }
    start = start + 1
  }
  let mut i = start
  while i < n {
    let ch = text.unsafe_get(i).to_int()
    if ch == ' '.to_int() || ch == '\t'.to_int() || ch == '\r'.to_int() {
      return text[start:i].to_owned()
    }
    i = i + 1
  }
  text[start:n].to_owned()
}

// ============================================================================
// Format scanners
// ============================================================================

///|
/// FASTA: each record starts at a line beginning with '>'; the key is the
/// first token of the header line.
fn seq_index_scan_fasta(
  data : String,
  lines : Array[IndexedLine],
) -> Array[SeqIndexEntry] raise SeqIOError {
  let entries : Array[SeqIndexEntry] = Array(capacity=0)
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.text.length() > 0 &&
      line.text.unsafe_get(0).to_int() == '>'.to_int() {
      let key = seq_index_first_token(line.text[1:].to_owned())
      if key.length() == 0 {
        raise SeqIOError(
          "FASTA record without an identifier at offset " +
          line.start.to_string(),
        )
      }
      let start = line.start
      let mut j = i + 1
      while j < lines.length() {
        if lines[j].text.length() > 0 &&
          lines[j].text.unsafe_get(0).to_int() == '>'.to_int() {
          break
        }
        j = j + 1
      }
      let end = if j < lines.length() { lines[j].start } else { data.length() }
      entries.push(SeqIndexEntry::new(key, start, end - start))
      i = j
    } else {
      i = i + 1
    }
  }
  entries
}

///|
/// Standard 4-line FASTQ: header, sequence, '+', quality.
fn seq_index_scan_fastq(
  data : String,
  lines : Array[IndexedLine],
) -> Array[SeqIndexEntry] raise SeqIOError {
  let entries : Array[SeqIndexEntry] = Array(capacity=0)
  let n = lines.length()
  let mut i = 0
  while i < n {
    let line = lines[i]
    if line.text.length() == 0 {
      i = i + 1
      continue
    }
    if line.text.unsafe_get(0).to_int() != '@'.to_int() {
      raise SeqIOError(
        "Expected FASTQ '@' header at line offset " + line.start.to_string(),
      )
    }
    if i + 3 >= n {
      raise SeqIOError(
        "Truncated FASTQ record starting at offset " + line.start.to_string(),
      )
    }
    let plus_line = lines[i + 2].text
    if plus_line.length() == 0 ||
      plus_line.unsafe_get(0).to_int() != '+'.to_int() {
      raise SeqIOError(
        "Malformed FASTQ record: expected '+' separator line at offset " +
        lines[i + 2].start.to_string(),
      )
    }
    let key = seq_index_first_token(line.text[1:].to_owned())
    let start = line.start
    let end = if i + 4 < n { lines[i + 4].start } else { data.length() }
    entries.push(SeqIndexEntry::new(key, start, end - start))
    i = i + 4
  }
  entries
}

///|
/// Formats whose records start at a fixed prefix line (LOCUS / ID) and
/// continue through the "//" terminator. `kind` selects key parsing.
fn seq_index_scan_delimiter(
  data : String,
  lines : Array[IndexedLine],
  prefix : String,
  kind : String,
) -> Array[SeqIndexEntry] raise SeqIOError {
  let entries : Array[SeqIndexEntry] = Array(capacity=0)
  let n = lines.length()
  let mut i = 0
  while i < n {
    let line = lines[i]
    if line.text.length() >= prefix.length() &&
      line.text[0:prefix.length()] == prefix {
      let start = line.start
      let key = if kind == "genbank" {
        // "LOCUS <name> length bp ..."
        let rest = line.text[prefix.length():].to_owned()
        seq_index_first_token(rest)
      } else {
        // EMBL: "ID   XX00001; ..." / Swiss: "ID   NAME   STANDARD; ..."
        let rest = line.text[prefix.length():].to_owned()
        let token = seq_index_first_token(rest)
        let semicolon = token.find(";")
        match semicolon {
          Some(idx) => token[0:idx].to_owned()
          None => token
        }
      }
      if key.length() == 0 {
        raise SeqIOError(
          "Record without an identifier at offset " + line.start.to_string(),
        )
      }
      // Find the start of the next record (or EOF).
      let mut j = i + 1
      while j < n {
        let candidate = lines[j].text
        if candidate.length() >= prefix.length() &&
          candidate[0:prefix.length()] == prefix {
          break
        }
        j = j + 1
      }
      let end = if j < n { lines[j].start } else { data.length() }
      entries.push(SeqIndexEntry::new(key, start, end - start))
      i = j
    } else {
      i = i + 1
    }
  }
  entries
}

///|
/// Tab format: one record per line, "id\tsequence".
fn seq_index_scan_tab(
  data : String,
  lines : Array[IndexedLine],
) -> Array[SeqIndexEntry] raise SeqIOError {
  let entries : Array[SeqIndexEntry] = Array(capacity=0)
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.text.length() == 0 {
      i = i + 1
      continue
    }
    let tab = line.text.find("\t")
    match tab {
      Some(pos) => {
        if pos == 0 {
          raise SeqIOError(
            "Tab record with empty identifier at offset " +
            line.start.to_string(),
          )
        }
        let key = line.text[0:pos].to_owned()
        let end = if i + 1 < lines.length() {
          lines[i + 1].start
        } else {
          data.length()
        }
        entries.push(SeqIndexEntry::new(key, line.start, end - line.start))
      }
      None =>
        raise SeqIOError(
          "Malformed tab record (no tab separator) at offset " +
          line.start.to_string(),
        )
    }
    i = i + 1
  }
  entries
}
