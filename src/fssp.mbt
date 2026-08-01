///|
/// FSSP file parser - Fold classification based on Structure-Structure
/// alignment of Proteins database.
///
/// References:
///   - Biopython Bio.FSSP module
///   - Holm L, Sander C (1996) Mapping the protein universe. Science
///     273(5275):595-603.
///   - Holm L, Sander C (1994) The FSSP database of structurally aligned
///     protein fold families. Nucleic Acids Res 22(17):3600-9.
///
/// The FSSP database contains structural alignments of proteins in the
/// PDB. Each FSSP file has a reference structure and a set of aligned
/// structures with sequence alignments and Z-scores.
///
/// All public identifiers are prefixed with `fssp_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// FSSP file header metadata.
pub struct FsspHeader {
  pdbid : String
  pdbid_date : String
  title : String
  compound : String
  source : String
  author : String
  seq_length : Int
  n_alignments : Int
  threshold : Double
}

///|
/// Create an FsspHeader.
pub fn FsspHeader::new(
  pdbid : String,
  pdbid_date : String,
  title : String,
  compound : String,
  source : String,
  author : String,
  seq_length : Int,
  n_alignments : Int,
  threshold : Double,
) -> FsspHeader {
  {
    pdbid,
    pdbid_date,
    title,
    compound,
    source,
    author,
    seq_length,
    n_alignments,
    threshold,
  }
}

///|
/// Get PDB ID of the reference structure.
pub fn FsspHeader::pdbid(self : FsspHeader) -> String {
  self.pdbid
}

///|
/// Get PDB ID date string.
pub fn FsspHeader::pdbid_date(self : FsspHeader) -> String {
  self.pdbid_date
}

///|
/// Get title.
pub fn FsspHeader::title(self : FsspHeader) -> String {
  self.title
}

///|
/// Get compound description.
pub fn FsspHeader::compound(self : FsspHeader) -> String {
  self.compound
}

///|
/// Get source organism.
pub fn FsspHeader::source(self : FsspHeader) -> String {
  self.source
}

///|
/// Get author.
pub fn FsspHeader::author(self : FsspHeader) -> String {
  self.author
}

///|
/// Get reference sequence length.
pub fn FsspHeader::seq_length(self : FsspHeader) -> Int {
  self.seq_length
}

///|
/// Get number of alignments.
pub fn FsspHeader::n_alignments(self : FsspHeader) -> Int {
  self.n_alignments
}

///|
/// Get Z-score threshold.
pub fn FsspHeader::threshold(self : FsspHeader) -> String {
  self.threshold.to_string()
}

///|
/// A single aligned sequence in an FSSP alignment block.
pub struct FsspAlignment {
  pdbid : String
  alignment_id : String
  start_position : Int
  aligned_sequence : String
  z_score : Double
  rmsd : Double
  lali : Int
  pid : Double
}

///|
/// Create an FsspAlignment.
pub fn FsspAlignment::new(
  pdbid : String,
  alignment_id : String,
  start_position : Int,
  aligned_sequence : String,
  z_score : Double,
  rmsd : Double,
  lali : Int,
  pid : Double,
) -> FsspAlignment {
  {
    pdbid,
    alignment_id,
    start_position,
    aligned_sequence,
    z_score,
    rmsd,
    lali,
    pid,
  }
}

///|
/// Get PDB ID.
pub fn FsspAlignment::pdbid(self : FsspAlignment) -> String {
  self.pdbid
}

///|
/// Get alignment ID.
pub fn FsspAlignment::alignment_id(self : FsspAlignment) -> String {
  self.alignment_id
}

///|
/// Get start position in the reference.
pub fn FsspAlignment::start_position(self : FsspAlignment) -> Int {
  self.start_position
}

///|
/// Get aligned sequence (with gaps and insertions).
pub fn FsspAlignment::aligned_sequence(self : FsspAlignment) -> String {
  self.aligned_sequence
}

///|
/// Get Z-score.
pub fn FsspAlignment::z_score(self : FsspAlignment) -> Double {
  self.z_score
}

///|
/// Get RMSD (root mean square deviation).
pub fn FsspAlignment::rmsd(self : FsspAlignment) -> Double {
  self.rmsd
}

///|
/// Get length of alignment.
pub fn FsspAlignment::lali(self : FsspAlignment) -> Int {
  self.lali
}

///|
/// Get percent identity.
pub fn FsspAlignment::pid(self : FsspAlignment) -> Double {
  self.pid
}

///|
/// Top-level container for parsed FSSP data.
pub struct FsspData {
  header : FsspHeader
  alignments : Array[FsspAlignment]
}

///|
/// Create an FsspData.
pub fn FsspData::new(
  header : FsspHeader,
  alignments : Array[FsspAlignment],
) -> FsspData {
  { header, alignments }
}

///|
/// Get header.
pub fn FsspData::header(self : FsspData) -> FsspHeader {
  self.header
}

///|
/// Get all alignments.
pub fn FsspData::alignments(self : FsspData) -> Array[FsspAlignment] {
  self.alignments
}

///|
/// Number of alignments.
pub fn FsspData::n_alignments(self : FsspData) -> Int {
  self.alignments.length()
}

///|
/// Get the reference (first) alignment.
pub fn FsspData::reference(self : FsspData) -> FsspAlignment? {
  if self.alignments.length() == 0 {
    None
  } else {
    Some(self.alignments[0])
  }
}

///|
/// Get alignments with Z-score above a threshold.
pub fn FsspData::filter_by_zscore(
  self : FsspData,
  min_zscore : Double,
) -> Array[FsspAlignment] {
  let result : Array[FsspAlignment] = []
  for a in self.alignments {
    if a.z_score >= min_zscore {
      result.push(a)
    }
  }
  result
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Tokenize a line by whitespace (space or tab).
fn fssp_tokenize(line : String) -> Array[String] {
  let tokens : Array[String] = []
  let n = line.length()
  let mut i = 0
  while i < n {
    while i < n {
      let c = line.unsafe_get(i).to_int()
      if c == ' '.to_int() || c == '\t'.to_int() {
        i = i + 1
      } else {
        break
      }
    }
    let start = i
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
/// Extract text after a keyword tag (first token).
fn fssp_extract_after_keyword(line : String, keyword : String) -> String {
  let trimmed = trim(line)
  if starts_with(trimmed, 0, keyword) {
    let rest = trimmed[keyword.length():trimmed.length()].to_string()
    trim(rest)
  } else {
    ""
  }
}

///|
/// Parse the FSSP header from the first few lines of the file.
fn fssp_parse_header(lines : Array[String]) -> (FsspHeader, Int) {
  let mut pdbid = ""
  let mut pdbid_date = ""
  let mut title = ""
  let mut compound = ""
  let mut source = ""
  let mut author = ""
  let mut seq_length = 0
  let mut n_alignments = 0
  let mut threshold = 0.0
  let mut idx = 0
  while idx < lines.length() {
    let line = trim(lines[idx].to_string())
    if line.length() == 0 {
      idx = idx + 1
      continue
    }
    if starts_with(line, 0, "##") {
      break
    }
    if starts_with(line, 0, "HEADER") {
      // HEADER \_1dfa_A 1 30-jul-1998
      let rest = fssp_extract_after_keyword(line, "HEADER")
      let tokens = fssp_tokenize(rest)
      if tokens.length() >= 3 {
        // Strip \_ prefix from PDB ID if present
        let raw_pdbid = tokens[0]
        if starts_with(raw_pdbid, 0, "\\_") {
          pdbid = raw_pdbid[2:raw_pdbid.length()].to_string()
        } else {
          pdbid = raw_pdbid
        }
        pdbid_date = tokens[2]
      }
    } else if starts_with(line, 0, "TITLE") {
      title = fssp_extract_after_keyword(line, "TITLE")
    } else if starts_with(line, 0, "COMPND") {
      compound = fssp_extract_after_keyword(line, "COMPND")
    } else if starts_with(line, 0, "SOURCE") {
      source = fssp_extract_after_keyword(line, "SOURCE")
    } else if starts_with(line, 0, "AUTHOR") {
      author = fssp_extract_after_keyword(line, "AUTHOR")
    } else if starts_with(line, 0, "SEQLENGTH") {
      let rest = fssp_extract_after_keyword(line, "SEQLENGTH")
      seq_length = parse_int(rest)
    } else if starts_with(line, 0, "ALIGNMENTS") {
      let rest = fssp_extract_after_keyword(line, "ALIGNMENTS")
      let tokens = fssp_tokenize(rest)
      if tokens.length() >= 1 {
        n_alignments = parse_int(tokens[0])
      }
    } else if starts_with(line, 0, "THRESHOLD") {
      let rest = fssp_extract_after_keyword(line, "THRESHOLD")
      match parse_double(rest) {
        Some(v) => threshold = v
        None => threshold = 0.0
      }
    }
    idx = idx + 1
  }
  (
    FsspHeader::new(
      pdbid,
      pdbid_date,
      title,
      compound,
      source,
      author,
      seq_length,
      n_alignments,
      threshold,
    ),
    idx,
  )
}

///|
/// Parse an alignment record line.
/// Format: pdbid pos seq [zscore rmsd lali pid]
fn fssp_parse_alignment_line(line : String) -> FsspAlignment? {
  let tokens = fssp_tokenize(line)
  if tokens.length() < 3 {
    return None
  }
  let pdbid = tokens[0]
  let start_pos = parse_int(tokens[1])
  let aligned_seq = tokens[2]
  let mut zscore = 0.0
  let mut rmsd = 0.0
  let mut lali = 0
  let mut pid = 0.0
  if tokens.length() >= 7 {
    match parse_double(tokens[3]) {
      Some(v) => zscore = v
      None => zscore = 0.0
    }
    match parse_double(tokens[4]) {
      Some(v) => rmsd = v
      None => rmsd = 0.0
    }
    lali = parse_int(tokens[5])
    match parse_double(tokens[6]) {
      Some(v) => pid = v
      None => pid = 0.0
    }
  }
  Some(
    FsspAlignment::new(
      pdbid,
      pdbid,
      start_pos,
      aligned_seq,
      zscore,
      rmsd,
      lali,
      pid,
    ),
  )
}

///|
/// Parse a complete FSSP file from text input.
pub fn fssp_parse(text : String) -> FsspData {
  let lines = split_by_char(text, '\n'.to_int())
  let (header, header_end) = fssp_parse_header(lines)
  let mut idx = header_end
  let alignments : Array[FsspAlignment] = []
  let mut in_alignment_section = false
  while idx < lines.length() {
    let line = trim(lines[idx].to_string())
    if line.length() == 0 {
      idx = idx + 1
      continue
    }
    // Check for section markers
    if starts_with(line, 0, "## ALIGNMENTS") || starts_with(line, 0, "##ALIGNMENTS") {
      in_alignment_section = true
      idx = idx + 1
      continue
    }
    if starts_with(line, 0, "## SEQUENCE") || starts_with(line, 0, "##SEQUENCE") {
      in_alignment_section = false
      idx = idx + 1
      continue
    }
    if starts_with(line, 0, "##") {
      in_alignment_section = false
      idx = idx + 1
      continue
    }
    if in_alignment_section && !starts_with(line, 0, "#") {
      match fssp_parse_alignment_line(line) {
        Some(a) => alignments.push(a)
        None => ()
      }
    }
    idx = idx + 1
  }
  FsspData::new(header, alignments)
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate sample FSSP file text with a reference and two aligned structures.
pub fn fssp_sample_text() -> String {
  let mut s = ""
  s = s + "HEADER \\_1dfa_A 1 30-jul-1998\n"
  s = s + "TITLE  crystal structure of human type i 17beta-hydroxysteroid dehydrogenase\n"
  s = s + "COMPND  17beta-hsd\n"
  s = s + "SOURCE  human\n"
  s = s + "AUTHOR  holm\n"
  s = s + "SEQLENGTH 25\n"
  s = s + "ALIGNMENTS 3 3\n"
  s = s + "THRESHOLD 2.0\n"
  s = s + "## ALIGNMENTS\n"
  s = s + "1dfa_A 1 MNIFVHEKDLFRTIVS 15.3 1.5 16 75.0\n"
  s = s + "1csa_A 1 MNIFVHEKDLFRTIVS 12.1 2.0 16 68.8\n"
  s = s + "1hsd_A 1 MNIFVHEKDLFRTIVS 8.5 3.1 15 62.5\n"
  s = s + "## SEQUENCE\n"
  s = s + "1dfa_A 1 MNIFVHEKDLFRTIVS\n"
  s = s + "## END\n"
  s
}

// ============================================================================
// Utilities
// ============================================================================

///|
/// Format an FsspAlignment as a readable string.
pub fn FsspAlignment::to_string(self : FsspAlignment) -> String {
  "FsspAlignment(pdbid=" +
  self.pdbid +
  ", pos=" +
  self.start_position.to_string() +
  ", seq=" +
  self.aligned_sequence +
  ", z=" +
  self.z_score.to_string() +
  ", rmsd=" +
  self.rmsd.to_string() +
  ", pid=" +
  self.pid.to_string() +
  ")"
}

///|
/// Format an FsspHeader as a readable string.
pub fn FsspHeader::to_string(self : FsspHeader) -> String {
  "FsspHeader(pdbid=" +
  self.pdbid +
  ", date=" +
  self.pdbid_date +
  ", seq_length=" +
  self.seq_length.to_string() +
  ", n_alignments=" +
  self.n_alignments.to_string() +
  ", threshold=" +
  self.threshold.to_string() +
  ")"
}

///|
/// Format an FsspData as a readable summary.
pub fn FsspData::to_string(self : FsspData) -> String {
  let mut s = "FSSP File Summary\n"
  s = s + self.header.to_string() + "\n"
  s = s + "Alignments (" + self.alignments.length().to_string() + "):\n"
  for a in self.alignments {
    s = s + "  " + a.to_string() + "\n"
  }
  s
}
