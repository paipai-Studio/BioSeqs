///|
/// PHD file parser - Phred base-calling output format.
///
/// References:
///   - Biopython Bio.Sequencing.Phd module
///   - Ewing B, Hillier L, Wendl MC, Green P (1998) Base-calling of
///     automated sequencer traces using phred. I. Accuracy assessment.
///     Genome Res 8(3):175-85.
///
/// The PHD format is produced by Phred and stores base calls with
/// per-base quality scores and peak positions (trace indices).
/// A PHD file contains one or more reads, each delimited by
/// BEGIN_SEQUENCE / END_SEQUENCE blocks.
///
/// All public identifiers are prefixed with `phd_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// Comment block metadata from a PHD read.
pub struct PhdComment {
  chromat_file : String
  phred_version : String
  phred_parameter_file : String
  trace_index_min : Int
  trace_index_max : Int
  trim_start : Int
  trim_end : Int
  trim_prob : Double
  chemistry : String
  dye : String
}

///|
/// Create a PhdComment.
pub fn PhdComment::new(
  chromat_file : String,
  phred_version : String,
  phred_parameter_file : String,
  trace_index_min : Int,
  trace_index_max : Int,
  trim_start : Int,
  trim_end : Int,
  trim_prob : Double,
  chemistry : String,
  dye : String,
) -> PhdComment {
  {
    chromat_file,
    phred_version,
    phred_parameter_file,
    trace_index_min,
    trace_index_max,
    trim_start,
    trim_end,
    trim_prob,
    chemistry,
    dye,
  }
}

///|
/// Get chromatogram file name.
pub fn PhdComment::chromat_file(self : PhdComment) -> String {
  self.chromat_file
}

///|
/// Get phred version.
pub fn PhdComment::phred_version(self : PhdComment) -> String {
  self.phred_version
}

///|
/// Get phred parameter file path.
pub fn PhdComment::phred_parameter_file(self : PhdComment) -> String {
  self.phred_parameter_file
}

///|
/// Get trace index minimum.
pub fn PhdComment::trace_index_min(self : PhdComment) -> Int {
  self.trace_index_min
}

///|
/// Get trace index maximum.
pub fn PhdComment::trace_index_max(self : PhdComment) -> Int {
  self.trace_index_max
}

///|
/// Get trim start position.
pub fn PhdComment::trim_start(self : PhdComment) -> Int {
  self.trim_start
}

///|
/// Get trim end position.
pub fn PhdComment::trim_end(self : PhdComment) -> Int {
  self.trim_end
}

///|
/// Get trim probability.
pub fn PhdComment::trim_prob(self : PhdComment) -> Double {
  self.trim_prob
}

///|
/// Get chemistry type.
pub fn PhdComment::chemistry(self : PhdComment) -> String {
  self.chemistry
}

///|
/// Get dye type.
pub fn PhdComment::dye(self : PhdComment) -> String {
  self.dye
}

///|
/// A single base call with quality and peak position.
pub struct PhdBase {
  base : String
  quality : Int
  peak_position : Int
}

///|
/// Create a PhdBase.
pub fn PhdBase::new(base : String, quality : Int, peak_position : Int) -> PhdBase {
  { base, quality, peak_position }
}

///|
/// Get base character.
pub fn PhdBase::base(self : PhdBase) -> String {
  self.base
}

///|
/// Get quality score.
pub fn PhdBase::quality(self : PhdBase) -> Int {
  self.quality
}

///|
/// Get peak position (trace index).
pub fn PhdBase::peak_position(self : PhdBase) -> Int {
  self.peak_position
}

///|
/// A complete PHD read: sequence name, base calls, and comment metadata.
pub struct PhdRead {
  name : String
  bases : Array[PhdBase]
  comment : PhdComment
}

///|
/// Create a PhdRead.
pub fn PhdRead::new(
  name : String,
  bases : Array[PhdBase],
  comment : PhdComment,
) -> PhdRead {
  { name, bases, comment }
}

///|
/// Get read name.
pub fn PhdRead::name(self : PhdRead) -> String {
  self.name
}

///|
/// Get base calls.
pub fn PhdRead::bases(self : PhdRead) -> Array[PhdBase] {
  self.bases
}

///|
/// Get comment block.
pub fn PhdRead::comment(self : PhdRead) -> PhdComment {
  self.comment
}

///|
/// Get the sequence string (concatenated bases).
pub fn PhdRead::sequence(self : PhdRead) -> String {
  let mut s = ""
  for b in self.bases {
    s = s + b.base
  }
  s
}

///|
/// Get quality scores as an array of integers.
pub fn PhdRead::quality(self : PhdRead) -> Array[Int] {
  let quals : Array[Int] = []
  for b in self.bases {
    quals.push(b.quality)
  }
  quals
}

///|
/// Get peak positions as an array of integers.
pub fn PhdRead::peak_positions(self : PhdRead) -> Array[Int] {
  let peaks : Array[Int] = []
  for b in self.bases {
    peaks.push(b.peak_position)
  }
  peaks
}

///|
/// Number of base calls.
pub fn PhdRead::length(self : PhdRead) -> Int {
  self.bases.length()
}

///|
/// Top-level container for a parsed PHD file (may contain multiple reads).
pub struct PhdFile {
  reads : Array[PhdRead]
}

///|
/// Create a PhdFile.
pub fn PhdFile::new(reads : Array[PhdRead]) -> PhdFile {
  { reads }
}

///|
/// Get all reads.
pub fn PhdFile::reads(self : PhdFile) -> Array[PhdRead] {
  self.reads
}

///|
/// Number of reads.
pub fn PhdFile::num_reads(self : PhdFile) -> Int {
  self.reads.length()
}

///|
/// Get read by index.
pub fn PhdFile::get(self : PhdFile, index : Int) -> PhdRead? {
  if index < 0 || index >= self.reads.length() {
    None
  } else {
    Some(self.reads[index])
  }
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Tokenize a line by whitespace (space or tab), skipping empty tokens.
fn phd_tokenize(line : String) -> Array[String] {
  let tokens : Array[String] = []
  let n = line.length()
  let mut i = 0
  while i < n {
    // Skip whitespace
    while i < n {
      let c = line.unsafe_get(i).to_int()
      if c == ' '.to_int() || c == '\t'.to_int() {
        i = i + 1
      } else {
        break
      }
    }
    let start = i
    // Read token
    while i < n {
      let c = line.unsafe_get(i).to_int()
      if c == ' '.to_int() || c == '\t'.to_int() {
        break
      }
      i = i + 1
    }
    if start < i {
      tokens.push(line[start:i].to_owned())
    }
  }
  tokens
}

///|
/// Extract the value part after "KEY:" in a comment line.
fn phd_extract_value(line : String) -> String {
  let n = line.length()
  let mut i = 0
  // Find the colon
  while i < n {
    if line.unsafe_get(i).to_int() == ':'.to_int() {
      i = i + 1
      // Skip whitespace after colon
      while i < n {
        let c = line.unsafe_get(i).to_int()
        if c == ' '.to_int() || c == '\t'.to_int() {
          i = i + 1
        } else {
          break
        }
      }
      return line[i:n].to_string()
    }
    i = i + 1
  }
  ""
}

///|
/// Parse a PHD comment block from lines between BEGIN_COMMENT and END_COMMENT.
fn phd_parse_comment(lines : Array[String], start : Int) -> (PhdComment, Int) {
  let mut chromat_file = ""
  let mut phred_version = ""
  let mut phred_parameter_file = ""
  let mut trace_index_min = 0
  let mut trace_index_max = 0
  let mut trim_start = 0
  let mut trim_end = 0
  let mut trim_prob = 0.0
  let mut chemistry = ""
  let mut dye = ""
  let mut idx = start
  while idx < lines.length() {
    let line = trim(lines[idx].to_string())
    if line == "END_COMMENT" {
      idx = idx + 1
      break
    }
    if starts_with(line, 0, "CHROMAT_FILE") {
      chromat_file = phd_extract_value(line)
    } else if starts_with(line, 0, "PHRED_VERSION") {
      phred_version = phd_extract_value(line)
    } else if starts_with(line, 0, "PHRED_PARAMETER_FILE") {
      phred_parameter_file = phd_extract_value(line)
    } else if starts_with(line, 0, "TRACE_INDEX_MIN") {
      trace_index_min = parse_int(phd_extract_value(line))
    } else if starts_with(line, 0, "TRACE_INDEX_MAX") {
      trace_index_max = parse_int(phd_extract_value(line))
    } else if starts_with(line, 0, "TRIM") {
      let parts = phd_tokenize(phd_extract_value(line))
      if parts.length() >= 3 {
        trim_start = parse_int(parts[0])
        trim_end = parse_int(parts[1])
        match parse_double(parts[2]) {
          Some(v) => trim_prob = v
          None => trim_prob = 0.0
        }
      }
    } else if starts_with(line, 0, "CHEMISTRY") {
      chemistry = phd_extract_value(line)
    } else if starts_with(line, 0, "DYE") {
      dye = phd_extract_value(line)
    }
    idx = idx + 1
  }
  (
    PhdComment::new(
      chromat_file,
      phred_version,
      phred_parameter_file,
      trace_index_min,
      trace_index_max,
      trim_start,
      trim_end,
      trim_prob,
      chemistry,
      dye,
    ),
    idx,
  )
}

///|
/// Parse a PHD DNA block from lines between BEGIN_DNA and END_DNA.
fn phd_parse_dna(lines : Array[String], start : Int) -> (Array[PhdBase], Int) {
  let bases : Array[PhdBase] = []
  let mut idx = start
  while idx < lines.length() {
    let line = trim(lines[idx].to_string())
    if line == "END_DNA" {
      idx = idx + 1
      break
    }
    if line.length() == 0 {
      idx = idx + 1
      continue
    }
    let tokens = phd_tokenize(line)
    if tokens.length() >= 3 {
      let base = tokens[0]
      let qual = parse_int(tokens[1])
      let peak = parse_int(tokens[2])
      bases.push(PhdBase::new(base, qual, peak))
    }
    idx = idx + 1
  }
  (bases, idx)
}

///|
/// Parse a single PHD sequence block (BEGIN_SEQUENCE ... END_SEQUENCE).
fn phd_parse_sequence(
  lines : Array[String],
  start : Int,
) -> (PhdRead?, Int) {
  let first_line = trim(lines[start].to_string())
  // Extract read name from "BEGIN_SEQUENCE name"
  let tokens = phd_tokenize(first_line)
  if tokens.length() < 2 {
    return (None, start + 1)
  }
  let read_name = tokens[1]
  let mut comment = PhdComment::new(
    "", "", "", 0, 0, 0, 0, 0.0, "", "",
  )
  let mut bases : Array[PhdBase] = []
  let mut idx = start + 1
  while idx < lines.length() {
    let line = trim(lines[idx].to_string())
    if line == "END_SEQUENCE" {
      idx = idx + 1
      break
    }
    if line == "BEGIN_COMMENT" {
      let (c, next) = phd_parse_comment(lines, idx + 1)
      comment = c
      idx = next
    } else if line == "BEGIN_DNA" {
      let (b, next) = phd_parse_dna(lines, idx + 1)
      bases = b
      idx = next
    } else {
      idx = idx + 1
    }
  }
  (Some(PhdRead::new(read_name, bases, comment)), idx)
}

///|
/// Parse a complete PHD file from text input.
///
/// Returns a PhdFile containing all reads found in the input.
pub fn phd_parse(text : String) -> PhdFile {
  let lines = split_by_char(text, '\n'.to_int())
  let reads : Array[PhdRead] = []
  let mut idx = 0
  while idx < lines.length() {
    let line = trim(lines[idx].to_string())
    if starts_with(line, 0, "BEGIN_SEQUENCE") {
      let (read_opt, next) = phd_parse_sequence(lines, idx)
      match read_opt {
        Some(r) => reads.push(r)
        None => ()
      }
      idx = next
    } else {
      idx = idx + 1
    }
  }
  PhdFile::new(reads)
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate sample PHD file text with one read.
/// The read "sample_read" has 10 bases with quality scores and peak positions.
pub fn phd_sample_text() -> String {
  let mut s = "BEGIN_SEQUENCE sample_read\n\n"
  s = s + "BEGIN_COMMENT\n"
  s = s + "CHROMAT_FILE: sample_read.ab1\n"
  s = s + "PHRED_VERSION: 0.020425.c\n"
  s = s + "PHRED_PARAMETER_FILE: /usr/local/etc/phred.dat\n"
  s = s + "TRACE_INDEX_MIN: 0\n"
  s = s + "TRACE_INDEX_MAX: 1000\n"
  s = s + "TRIM: 0 10 0.05\n"
  s = s + "CHEMISTRY: term\n"
  s = s + "DYE: big\n"
  s = s + "END_COMMENT\n\n"
  s = s + "BEGIN_DNA\n"
  let bases = ["A", "C", "G", "T", "A", "C", "G", "T", "A", "C"]
  let quals = [35, 40, 45, 50, 38, 42, 47, 33, 39, 44]
  let peaks = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
  for i in 0..<10 {
    s = s + bases[i] + " " + quals[i].to_string() + " " + peaks[i].to_string() + "\n"
  }
  s = s + "END_DNA\n\n"
  s = s + "END_SEQUENCE\n"
  s
}

// ============================================================================
// Utilities
// ============================================================================

///|
/// Format a PhdRead as a readable string.
pub fn PhdRead::to_string(self : PhdRead) -> String {
  let mut s = "PhdRead(name=" + self.name + ", length=" + self.length().to_string() + ")\n"
  s = s + "  chromat_file: " + self.comment.chromat_file + "\n"
  s = s + "  phred_version: " + self.comment.phred_version + "\n"
  s = s + "  sequence: " + self.sequence() + "\n"
  s = s + "  quality: "
  for i in 0..<self.bases.length() {
    if i > 0 {
      s = s + " "
    }
    s = s + self.bases[i].quality.to_string()
  }
  s
}

///|
/// Format a PhdFile as a readable summary.
pub fn PhdFile::to_string(self : PhdFile) -> String {
  let mut s = "PhdFile(" + self.num_reads().to_string() + " reads)\n"
  for r in self.reads {
    s = s + r.to_string() + "\n"
  }
  s
}
