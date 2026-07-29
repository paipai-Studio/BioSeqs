///|
/// ACE contig assembly format parser.
/// Inspired by Biopython's Bio.Sequencing.Ace module.
///
/// The ACE format is produced by the ACE assembly program and describes
/// assembled contigs with their constituent reads, base qualities, and
/// alignment positions.

// ===== Data structures =====

///|
/// A single read from an ACE assembly.
pub struct AceRead {
  read_id : String
  sequence : String
  quality : Array[Int]
  clip_start : Int
  clip_end : Int
  sequence_length : Int
  strand : String
  chromat_file : String
  chemistry : String
  dye : String
}

///|
/// Create a new AceRead.
pub fn AceRead::new(
  read_id : String,
  sequence : String,
  quality : Array[Int],
  clip_start : Int,
  clip_end : Int,
  strand : String,
  chromat_file : String,
  chemistry : String,
  dye : String,
) -> AceRead {
  AceRead::{
    read_id,
    sequence,
    quality,
    clip_start,
    clip_end,
    sequence_length: sequence.length(),
    strand,
    chromat_file,
    chemistry,
    dye,
  }
}

///|
/// A contig from an ACE assembly.
pub struct AceContig {
  contig_name : String
  sequence : String
  length : Int
  reads : Array[AceRead]
  base_qualities : Array[Int]
  alignment_complement : String
  reverse_complemented : Bool
}

///|
/// Create a new AceContig.
pub fn AceContig::new(
  contig_name : String,
  sequence : String,
  reads : Array[AceRead],
  base_qualities : Array[Int],
  alignment_complement : String,
  reverse_complemented : Bool,
) -> AceContig {
  AceContig::{
    contig_name,
    sequence,
    length: sequence.length(),
    reads,
    base_qualities,
    alignment_complement,
    reverse_complemented,
  }
}

///|
/// An alignment entry within a contig.
pub struct AceAlignment {
  contig_id : String
  start_position : Int
  sequence : String
  quality : Array[Int]
}

///|
/// Create a new AceAlignment.
pub fn AceAlignment::new(
  contig_id : String,
  start_position : Int,
  sequence : String,
  quality : Array[Int],
) -> AceAlignment {
  AceAlignment::{ contig_id, start_position, sequence, quality }
}

///|
/// Top-level container for ACE assembly data.
pub struct AceData {
  contigs : Array[AceContig]
  total_reads : Int
  total_contigs : Int
  base_qual_threshold : Int
}

///|
/// Create a new AceData.
pub fn AceData::new(
  contigs : Array[AceContig],
  base_qual_threshold : Int,
) -> AceData {
  let mut sum = 0
  let mut i = 0
  while i < contigs.length() {
    sum = sum + contigs[i].reads.length()
    i = i + 1
  }
  AceData::{
    contigs,
    total_reads: sum,
    total_contigs: contigs.length(),
    base_qual_threshold,
  }
}

// ===== Internal parsing helpers =====

///|
/// Split a line into tokens by whitespace.
fn ace_tokenize(line : String) -> Array[String] {
  let tokens : Array[String] = Array::new()
  let n = line.length()
  let mut i = 0
  while i < n {
    let mut j = i
    while j < n && line.unsafe_get(j).to_int() == ' '.to_int() {
      j = j + 1
    }
    let start = j
    while j < n && line.unsafe_get(j).to_int() != ' '.to_int() {
      j = j + 1
    }
    if start < j {
      tokens.push(line[start:j].to_owned())
    }
    i = j
  }
  tokens
}

///|
/// Parse quality values from quality lines following sequence data.
fn ace_parse_quality_lines(
  lines : Array[String],
  start_idx : Int,
  expected_length : Int,
) -> (Array[Int], Int) {
  let qual : Array[Int] = Array::new()
  let mut idx = start_idx

  while idx < lines.length() && qual.length() < expected_length {
    let line = lines[idx].to_string()
    let trimmed = trim(line)
    if trimmed.length() == 0 {
      idx = idx + 1
      continue
    }

    let first_char = trimmed.unsafe_get(0).to_int()
    if first_char == 'q'.to_int() ||
      (first_char >= '0'.to_int() && first_char <= '9'.to_int()) {
      let rest = if first_char == 'q'.to_int() {
        trim(trimmed[1:].to_owned())
      } else {
        trimmed
      }
      let qtokens = rest.split(" ").to_array()
      let mut qi = 0
      while qi < qtokens.length() && qual.length() < expected_length {
        let qt = trim(qtokens[qi].to_owned())
        if qt.length() > 0 {
          qual.push(parse_int(qt))
        }
        qi = qi + 1
      }
      idx = idx + 1
      continue
    }

    break
  }

  (qual, idx)
}

///|
/// Write a string to StringBuilder, wrapping at given width.
fn ace_write_wrapped(buf : StringBuilder, s : String, line_width : Int) -> Unit {
  let n = s.length()
  let mut i = 0
  while i < n {
    let end = if i + line_width > n { n } else { i + line_width }
    buf.write_string(s[i:end].to_owned())
    if end < n {
      buf.write_string("\n")
    }
    i = end
  }
}

///|
/// Convert quality array to space-separated string via StringBuilder.
fn ace_quality_to_string(qualities : Array[Int]) -> String {
  let sb = StringBuilder::new()
  let mut i = 0
  while i < qualities.length() {
    if i > 0 {
      sb.write_string(" ")
    }
    sb.write_string(qualities[i].to_string())
    i = i + 1
  }
  sb.to_string()
}

// ===== Main parsing functions =====

///|
/// Parse ACE format content into AceData.
pub fn ace_parse(content : String) -> AceData {
  let contigs : Array[AceContig] = Array::new()

  let lines = split_by_char(content, '\n'.to_int())
  let n = lines.length()

  let reads_map : Map[String, AceRead] = Map([], capacity=64)

  let mut i = 0
  while i < n {
    let line = lines[i].to_string()
    let trimmed = trim(line)
    if trimmed.length() == 0 {
      i = i + 1
      continue
    }

    if trimmed.has_prefix("RD ") {
      let read = ace_parse_read_block(lines, i)
      match read {
        Some(r) => {
          reads_map[r.read_id] = r
          i = ace_skip_read_block(lines, i)
        }
        None => {
          i = i + 1
        }
      }
      continue
    }

    if trimmed.has_prefix("CT ") {
      let contig_and_new_idx = ace_parse_contig_block(lines, i, reads_map)
      match contig_and_new_idx {
        Some((contig, new_idx)) => {
          contigs.push(contig)
          i = new_idx
        }
        None => {
          i = i + 1
        }
      }
      continue
    }

    i = i + 1
  }

  AceData::new(contigs, 0)
}

///|
/// Parse only the reads section from ACE content.
pub fn ace_parse_reads(content : String) -> Array[AceRead] {
  let reads : Array[AceRead] = Array::new()
  let lines = split_by_char(content, '\n'.to_int())
  let n = lines.length()
  let mut i = 0
  while i < n {
    let line = lines[i].to_string()
    let trimmed = trim(line)
    if trimmed.has_prefix("RD ") {
      let read = ace_parse_read_block(lines, i)
      match read {
        Some(r) => reads.push(r)
        None => { () }
      }
      i = ace_skip_read_block(lines, i)
    } else {
      i = i + 1
    }
  }
  reads
}

///|
/// Parse only the contigs section from ACE content.
pub fn ace_parse_contigs(content : String) -> Array[AceContig] {
  let contigs : Array[AceContig] = Array::new()
  let lines = split_by_char(content, '\n'.to_int())
  let reads_map : Map[String, AceRead] = Map([], capacity=64)
  let n = lines.length()

  let mut i = 0
  while i < n {
    let line = lines[i].to_string()
    let trimmed = trim(line)
    if trimmed.has_prefix("RD ") {
      let read = ace_parse_read_block(lines, i)
      match read {
        Some(r) => {
          reads_map[r.read_id] = r
          i = ace_skip_read_block(lines, i)
        }
        None => {
          i = i + 1
        }
      }
    } else if trimmed.has_prefix("CT ") {
      let contig_and_new_idx = ace_parse_contig_block(lines, i, reads_map)
      match contig_and_new_idx {
        Some((contig, new_idx)) => {
          contigs.push(contig)
          i = new_idx
        }
        None => {
          i = i + 1
        }
      }
    } else {
      i = i + 1
    }
  }
  contigs
}

///|
/// Parse a single RD block starting at the given line index.
fn ace_parse_read_block(
  lines : Array[String],
  start_idx : Int,
) -> AceRead? {
  if start_idx >= lines.length() {
    return None
  }
  let header = trim(lines[start_idx].to_string())
  if !header.has_prefix("RD ") {
    return None
  }

  let tokens = ace_tokenize(header)
  if tokens.length() < 4 {
    return None
  }

  let read_id = tokens[1]
  let read_len = parse_int(tokens[2])
  let chemistry = if tokens.length() >= 5 { tokens[3] } else { "" }
  let dye = if tokens.length() >= 6 { tokens[4] } else { "" }

  let mut idx = start_idx + 1

  // Parse sequence
  let mut sequence = ""
  while idx < lines.length() && sequence.length() < read_len {
    let line = lines[idx].to_string()
    let trimmed = trim(line)
    if trimmed.length() == 0 {
      idx = idx + 1
      continue
    }
    let first_char = trimmed.unsafe_get(0).to_int()
    if first_char == 'q'.to_int() ||
      (first_char >= '0'.to_int() && first_char <= '9'.to_int()) {
      break
    }
    if (first_char >= 'A'.to_int() && first_char <= 'Z'.to_int()) ||
      (first_char >= 'a'.to_int() && first_char <= 'z'.to_int()) ||
      first_char == '-'.to_int() || first_char == '.'.to_int() {
      let remaining = read_len - sequence.length()
      let chunk = if trimmed.length() > remaining {
        trimmed[0:remaining].to_owned()
      } else {
        trimmed
      }
      sequence = sequence + chunk
      idx = idx + 1
    } else {
      break
    }
  }

  // Parse quality
  let (quality, new_idx) = ace_parse_quality_lines(lines, idx, read_len)
  idx = new_idx

  Some(AceRead::new(
    read_id,
    sequence,
    quality,
    1,
    read_len,
    "+",
    "",
    chemistry,
    dye,
  ))
}

///|
/// Skip past an RD block, returns next index.
fn ace_skip_read_block(lines : Array[String], start_idx : Int) -> Int {
  if start_idx >= lines.length() {
    return start_idx + 1
  }
  let header = trim(lines[start_idx].to_string())
  let tokens = ace_tokenize(header)
  if tokens.length() < 3 {
    return start_idx + 1
  }
  let read_len = parse_int(tokens[2])

  let mut idx = start_idx + 1
  let mut seq_chars = 0

  // Skip sequence lines
  while idx < lines.length() && seq_chars < read_len {
    let line = lines[idx].to_string()
    let trimmed = trim(line)
    if trimmed.length() == 0 {
      idx = idx + 1
      continue
    }
    let first_char = trimmed.unsafe_get(0).to_int()
    if first_char == 'q'.to_int() {
      idx = idx + 1
      break
    }
    if (first_char >= '0'.to_int() && first_char <= '9'.to_int()) {
      break
    }
    if (first_char >= 'A'.to_int() && first_char <= 'Z'.to_int()) ||
      (first_char >= 'a'.to_int() && first_char <= 'z'.to_int()) ||
      first_char == '-'.to_int() || first_char == '.'.to_int() {
      seq_chars = seq_chars + trimmed.length()
      idx = idx + 1
    } else {
      break
    }
  }

  // Skip quality lines
  while idx < lines.length() {
    let line = lines[idx].to_string()
    let trimmed = trim(line)
    if trimmed.length() == 0 {
      idx = idx + 1
      continue
    }
    let first_char = trimmed.unsafe_get(0).to_int()
    if first_char == 'q'.to_int() ||
      (first_char >= '0'.to_int() && first_char <= '9'.to_int()) {
      idx = idx + 1
    } else {
      break
    }
  }

  idx
}

///|
/// Parse a CT block (contig) starting at the given line index.
fn ace_parse_contig_block(
  lines : Array[String],
  start_idx : Int,
  reads_map : Map[String, AceRead],
) -> (AceContig, Int)? {
  if start_idx >= lines.length() {
    return None
  }
  let header = trim(lines[start_idx].to_string())
  if !header.has_prefix("CT ") {
    return None
  }

  let tokens = ace_tokenize(header)
  if tokens.length() < 5 {
    return None
  }

  let contig_name = tokens[1]
  let contig_len = parse_int(tokens[2])

  let mut idx = start_idx + 1

  // Parse consensus sequence
  let mut sequence = ""
  while idx < lines.length() && sequence.length() < contig_len {
    let line = lines[idx].to_string()
    let trimmed = trim(line)
    if trimmed.length() == 0 {
      idx = idx + 1
      continue
    }
    let first_char = trimmed.unsafe_get(0).to_int()
    if first_char == 'q'.to_int() ||
      (first_char >= '0'.to_int() && first_char <= '9'.to_int()) {
      break
    }
    if (first_char >= 'A'.to_int() && first_char <= 'Z'.to_int()) ||
      (first_char >= 'a'.to_int() && first_char <= 'z'.to_int()) ||
      first_char == '-'.to_int() || first_char == '.'.to_int() {
      let remaining = contig_len - sequence.length()
      let chunk = if trimmed.length() > remaining {
        trimmed[0:remaining].to_owned()
      } else {
        trimmed
      }
      sequence = sequence + chunk
      idx = idx + 1
    } else {
      break
    }
  }

  // Parse base qualities
  let (base_qualities, new_idx) = ace_parse_quality_lines(lines, idx, contig_len)
  idx = new_idx

  // Parse AF alignment lines
  let reads : Array[AceRead] = Array::new()
  let alignment_positions : Map[String, (Int, Int, Int, Int, Int, Int, String)] = Map([], capacity=32)

  while idx < lines.length() {
    let line = lines[idx].to_string()
    let trimmed = trim(line)
    if trimmed.length() == 0 {
      idx = idx + 1
      continue
    }

    if trimmed.has_prefix("AF ") {
      let af_tokens = ace_tokenize(trimmed)
      if af_tokens.length() >= 8 {
        let read_name = af_tokens[1]
        let contig_start = parse_int(af_tokens[2])
        let contig_end = parse_int(af_tokens[3])
        let read_start = parse_int(af_tokens[4])
        let read_end = parse_int(af_tokens[5])
        let qual_start = parse_int(af_tokens[6])
        let qual_end = parse_int(af_tokens[7])
        let strand = if af_tokens.length() >= 9 { af_tokens[8] } else { "+" }

        alignment_positions[read_name] = (
          contig_start,
          contig_end,
          read_start,
          read_end,
          qual_start,
          qual_end,
          strand,
        )
      }
      idx = idx + 1
      continue
    }

    if trimmed.has_prefix("MQ ") {
      idx = idx + 1
      continue
    }

    if trimmed.has_prefix("CT ") ||
      trimmed.has_prefix("AF ") ||
      trimmed.has_prefix("RD ") ||
      trimmed.has_prefix("AQ ") ||
      trimmed.has_prefix("CO ") {
      break
    }

    idx = idx + 1
  }

  // Build read list from alignment positions
  let read_names = alignment_positions.keys().collect()
  let mut ri = 0
  while ri < read_names.length() {
    let read_name = read_names[ri]
    let pos_info = alignment_positions[read_name]
    match reads_map.get(read_name) {
      Some(read) => {
        let strand = pos_info.6
        let new_read = AceRead::new(
          read.read_id,
          read.sequence,
          read.quality,
          pos_info.4,
          pos_info.5,
          strand,
          read.chromat_file,
          read.chemistry,
          read.dye,
        )
        reads.push(new_read)
      }
      None => {
        let placeholder = AceRead::new(
          read_name,
          "",
          [],
          pos_info.4,
          pos_info.5,
          pos_info.6,
          "",
          "",
          "",
        )
        reads.push(placeholder)
      }
    }
    ri = ri + 1
  }

  Some((
    AceContig::new(
      contig_name,
      sequence,
      reads,
      base_qualities,
      "",
      false,
    ),
    idx,
  ))
}

// ===== Serialization =====

///|
/// Serialize AceData back to ACE format string.
pub fn ace_to_string(data : AceData) -> String {
  let buf = StringBuilder::new()

  let mut ci = 0
  while ci < data.contigs.length() {
    let contig = data.contigs[ci]

    buf.write_string("AF ")
    buf.write_string(contig.contig_name)
    buf.write_string(" ")
    buf.write_string(contig.length.to_string())
    buf.write_string("\n")

    let mut ri = 0
    while ri < contig.reads.length() {
      let read = contig.reads[ri]
      buf.write_string("RD ")
      buf.write_string(read.read_id)
      buf.write_string(" ")
      buf.write_string(read.sequence_length.to_string())
      buf.write_string(" ")
      buf.write_string(read.chemistry)
      buf.write_string(" ")
      buf.write_string(read.dye)
      buf.write_string("\n")

      ace_write_wrapped(buf, read.sequence, 60)
      buf.write_string("\n")

      let qual_str = ace_quality_to_string(read.quality)
      ace_write_wrapped(buf, qual_str, 60)
      buf.write_string("\n")

      ri = ri + 1
    }

    buf.write_string("CT ")
    buf.write_string(contig.contig_name)
    buf.write_string(" ")
    buf.write_string(contig.length.to_string())
    buf.write_string(" 0 ")
    buf.write_string(contig.length.to_string())
    buf.write_string("\n")

    ace_write_wrapped(buf, contig.sequence, 60)
    buf.write_string("\n")

    let bq_str = ace_quality_to_string(contig.base_qualities)
    ace_write_wrapped(buf, bq_str, 60)
    buf.write_string("\n")

    let mut ri2 = 0
    while ri2 < contig.reads.length() {
      let read = contig.reads[ri2]
      buf.write_string("AF ")
      buf.write_string(read.read_id)
      buf.write_string(" 0 ")
      buf.write_string(contig.length.to_string())
      buf.write_string(" ")
      buf.write_string(read.clip_start.to_string())
      buf.write_string(" ")
      buf.write_string(read.clip_end.to_string())
      buf.write_string(" ")
      buf.write_string(read.clip_start.to_string())
      buf.write_string(" ")
      buf.write_string(read.clip_end.to_string())
      buf.write_string(" ")
      buf.write_string(read.strand)
      buf.write_string("\n")
      ri2 = ri2 + 1
    }

    ci = ci + 1
  }

  buf.to_string()
}

// ===== Query/utility functions =====

///|
/// Get a contig by name from AceData.
pub fn ace_get_contig(data : AceData, name : String) -> AceContig? {
  let mut i = 0
  while i < data.contigs.length() {
    if data.contigs[i].contig_name == name {
      return Some(data.contigs[i])
    }
    i = i + 1
  }
  None
}

///|
/// Get the length of a contig by name.
pub fn ace_contig_length(data : AceData, name : String) -> Int {
  match ace_get_contig(data, name) {
    Some(contig) => contig.length
    None => 0
  }
}

///|
/// Get all reads in a contig by name.
pub fn ace_contig_reads(data : AceData, name : String) -> Array[AceRead] {
  match ace_get_contig(data, name) {
    Some(contig) => contig.reads.copy()
    None => []
  }
}

///|
/// Calculate read coverage per position along a contig.
pub fn ace_read_coverage(contig : AceContig) -> Map[Int, Int] {
  let coverage : Map[Int, Int] = Map([], capacity=contig.length)

  let mut pos = 1
  while pos <= contig.length {
    coverage[pos] = 0
    pos = pos + 1
  }

  let mut ri = 0
  while ri < contig.reads.length() {
    let read = contig.reads[ri]
    let start = if read.clip_start > 0 { read.clip_start } else { 1 }
    let end = if read.clip_end > 0 {
      if read.clip_end > contig.length { contig.length } else { read.clip_end }
    } else {
      contig.length
    }

    let mut p = start
    while p <= end && p <= contig.length {
      match coverage.get(p) {
        Some(val) => coverage[p] = val + 1
        None => coverage[p] = 1
      }
      p = p + 1
    }
    ri = ri + 1
  }

  coverage
}

///|
/// Get the consensus sequence of a contig.
pub fn ace_consensus_sequence(contig : AceContig) -> String {
  contig.sequence
}

///|
/// Calculate GC content of a contig.
pub fn ace_contig_gc_content(contig : AceContig) -> Double {
  let seq = contig.sequence
  if seq.length() == 0 {
    return 0.0
  }

  let mut gc_count = 0
  let mut i = 0
  while i < seq.length() {
    let c = seq.unsafe_get(i).to_int()
    if c == 'G'.to_int() ||
      c == 'C'.to_int() ||
      c == 'g'.to_int() ||
      c == 'c'.to_int() {
      gc_count = gc_count + 1
    }
    i = i + 1
  }

  gc_count.to_double() / seq.length().to_double()
}

// ===== Sample data generators =====

///|
/// Create sample ACE content for testing.
pub fn sample_ace_content() -> String {
  "AF contig1 60\n" +
    "RD read1 60 chemistry1 dye1\n" +
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n" +
    "q 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30 30\n" +
    "RD read2 60 chemistry2 dye2\n" +
    "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n" +
    "q 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25 25\n" +
    "CT contig1 60 0 60\n" +
    "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n" +
    "q 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40 40\n" +
    "AF read1 1 60 1 60 1 60 +\n" +
    "AF read2 1 60 1 60 1 60 -\n"
}

pub fn sample_ace_data() -> AceData {
  ace_parse(sample_ace_content())
}