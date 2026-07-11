///|
/// FASTA index entry (FAI format), matching samtools faidx specification.
/// Each entry contains 5 fields:
/// - name: sequence identifier
/// - length: total length in bases (excluding newlines)
/// - offset: byte offset in FASTA content of first base
/// - line_bases: number of bases per line (excluding newline)
/// - line_width: number of bytes per line (including newline)
pub struct FaiEntry {
  name : String
  length : Int
  offset : Int
  line_bases : Int
  line_width : Int
} derive(Eq, Debug)

///|
/// FASTA index for fast random access.
/// Uses a .fai index format for O(1) sequence lookup
/// and O(m) subsequence retrieval where m is the length of the region.
pub struct FastaIndex {
  entries : @hashmap.HashMap[String, FaiEntry]
  fasta_content : String
  entry_list : Array[FaiEntry]
}

///|
/// Create a new FastaIndex from a map of entries.
pub fn FastaIndex::new(entries : @hashmap.HashMap[String, FaiEntry], fasta_content : String) -> FastaIndex {
  let entry_list : Array[FaiEntry] = Array::new()
  for entry in entries.values() {
    entry_list.push(entry)
  }
  FastaIndex::{ entries, fasta_content, entry_list }
}

///|
/// Parse a string to Int, returning 0 on error.
fn fai_parse_int(s : StringView) -> Int {
  let n = s.length()
  if n == 0 {
    return 0
  }
  let mut result = 0
  let mut i = 0
  let zero = '0'.to_int()
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 {
      break
    }
    result = result * 10 + (c - zero)
    i = i + 1
  }
  result
}

///|
/// Build a FASTA index from FASTA content string.
/// Returns the FastaIndex.
pub fn FastaIndex::build(fasta_content : String) -> FastaIndex {
  let entries = @hashmap.HashMap([])
  let lines = fasta_content.split("\n").to_array()
  
  let mut i = 0
  let n = lines.length()
  let mut total_offset = 0
  
  while i < n {
    let line = lines[i]
    let line_len = line.length()
    if line_len == 0 || line.unsafe_get(0).to_int() != '>'.to_int() {
      total_offset = total_offset + line_len + 1
      i = i + 1
      continue
    }
    
    let name = if line_len > 1 {
      let title = line[1:line_len].to_owned()
      match title.find(" ") {
        Some(pos) => title[0:pos].to_owned()
        None => title
      }
    } else {
      ""
    }
    
    let header_offset = total_offset + line_len + 1
    let mut sequence_len = 0
    let mut line_bases = 0
    let mut line_width = 0
    let mut first_seq_line = true
    
    total_offset = total_offset + line_len + 1
    i = i + 1
    while i < n {
      let seq_line = lines[i]
      let seq_len = seq_line.length()
      if seq_len > 0 && seq_line.unsafe_get(0).to_int() == '>'.to_int() {
        break
      }
      let mut bases_in_line = 0
      let mut j = 0
      while j < seq_len {
        let c = seq_line.unsafe_get(j)
        if c.to_int() != ' '.to_int() && c.to_int() != '\t'.to_int() {
          bases_in_line = bases_in_line + 1
        }
        j = j + 1
      }
      if bases_in_line > 0 {
        sequence_len = sequence_len + bases_in_line
        if first_seq_line {
          line_bases = bases_in_line
          line_width = seq_len + 1
          first_seq_line = false
        }
      }
      total_offset = total_offset + seq_len + 1
      i = i + 1
    }
    
    if line_width == 0 {
      line_width = line_bases + 1
    }
    
    let entry = FaiEntry::{
      name,
      length: sequence_len,
      offset: header_offset,
      line_bases,
      line_width,
    }
    entries.set(name, entry)
  }
  
  FastaIndex::new(entries, fasta_content)
}

///|
/// Get all sequence names (keys) in the index.
pub fn FastaIndex::keys(self : FastaIndex) -> Array[String] {
  self.entries.keys().to_array()
}

///|
/// Get the number of sequences in the index.
pub fn FastaIndex::len(self : FastaIndex) -> Int {
  self.entry_list.length()
}

///|
/// Check if a sequence name exists in the index.
pub fn FastaIndex::contains(self : FastaIndex, name : String) -> Bool {
  self.entries.contains(name)
}

///|
/// Get the FaiEntry for a sequence name, or None if not found.
pub fn FastaIndex::get_entry(self : FastaIndex, name : String) -> FaiEntry? {
  self.entries.get(name)
}

///|
/// Get the length of a sequence.
pub fn FastaIndex::get_length(self : FastaIndex, name : String) -> Int? {
  match self.entries.get(name) {
    Some(entry) => Some(entry.length)
    None => None
  }
}

///|
/// Calculate the byte offset for a given 1-based position in the sequence.
fn calc_byte_offset(entry : FaiEntry, pos : Int) -> Int {
  let base_pos = pos - 1
  let lines_before = base_pos / entry.line_bases
  entry.offset + base_pos + lines_before * (entry.line_width - entry.line_bases)
}

///|
/// Fetch a subsequence from the FASTA content.
/// start and end are 1-based coordinates, matching pyfaidx convention.
/// Returns None if the sequence name is not found.
/// The result includes both start and end positions.
pub fn FastaIndex::fetch(self : FastaIndex, name : String, start : Int, end : Int) -> String? {
  match self.entries.get(name) {
    None => None
    Some(entry) => {
      let s = if start < 1 { 1 } else { start }
      let e = if end > entry.length { entry.length } else { end }
      if s > e {
        return Some("")
      }
      
      let start_offset = calc_byte_offset(entry, s)
      let end_offset = calc_byte_offset(entry, e)
      let read_length = end_offset - start_offset + entry.line_width
      
      let content = self.fasta_content
      let raw = if start_offset + read_length <= content.length() {
        content[start_offset:(start_offset + read_length)].to_owned()
      } else {
        content[start_offset:content.length()].to_owned()
      }
      
      let result = StringBuilder::new()
      let mut i = 0
      while i < raw.length() {
        let c = raw.unsafe_get(i)
        let ci = c.to_int()
        if ci != '\n'.to_int() && ci != '\r'.to_int() && 
           ci != ' '.to_int() && ci != '\t'.to_int() {
          result.write_char(c.unsafe_to_char())
        }
        i = i + 1
      }
      
      let full_seq = result.to_string()
      Some(full_seq[0:(e - s + 1)].to_owned())
    }
  }
}

///|
/// Fetch the entire sequence for a given name.
pub fn FastaIndex::fetch_full(self : FastaIndex, name : String) -> String? {
  match self.entries.get(name) {
    None => None
    Some(entry) => self.fetch(name, 1, entry.length)
  }
}

///|
/// Generate the .fai index content string.
pub fn FastaIndex::to_fai_string(self : FastaIndex) -> String {
  let buf = StringBuilder::new()
  let keys = self.keys()
  let n = keys.length()
  let sorted = Array::make(n, "")
  let mut i = 0
  while i < n {
    sorted[i] = keys[i]
    i = i + 1
  }
  i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      if sorted[i] > sorted[j] {
        let tmp = sorted[i]
        sorted[i] = sorted[j]
        sorted[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  i = 0
  while i < n {
    let key = sorted[i]
    match self.entries.get(key) {
      None => ignore(self.entries)
      Some(entry) => {
        buf.write_string(entry.name)
        buf.write_char('\t')
        buf.write_string(entry.length.to_string())
        buf.write_char('\t')
        buf.write_string(entry.offset.to_string())
        buf.write_char('\t')
        buf.write_string(entry.line_bases.to_string())
        buf.write_char('\t')
        buf.write_string(entry.line_width.to_string())
        buf.write_char('\n')
      }
    }
    i = i + 1
  }
  buf.to_string()
}
