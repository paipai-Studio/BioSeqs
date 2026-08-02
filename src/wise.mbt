// Bio.Wise - Parser for Wise2 (GeneWise/GenomeWise) output, mirroring
// BioPython's Bio.Wise.
//
// Wise2 is a package for comparing DNA sequences to protein sequences.
// The main tools are:
//   genewise   - compares a protein sequence to a genomic DNA sequence
//   genomewise - combines multiple genewise analyses across gene segments
//
// GeneWise output format (text):
//   1. Header with command line and parameters
//   2. Alignment between protein and gene (with introns marked)
//   3. Gene structure: exons, introns, splice sites
//   4. Score and statistics
//
// The output shows:
//   - Protein sequence aligned to translated gene
//   - Match/mismatch indicators
//   - Intron regions with phase information
//   - Exon boundaries and splice site scores
//   - Overall alignment score

// ============================================================================
// Data structures
// ============================================================================

///|
/// Type of a Wise2 alignment block.
pub enum WiseBlockType {
  WiseExon // exon region
  WiseIntron // intron region
  WiseMatch // match block
  WiseMismatch // mismatch block
  WiseInsertion // insertion (gap in protein)
  WiseDeletion // deletion (gap in gene)
} derive(Eq)

///|
/// Convert block type to string.
pub fn WiseBlockType::to_string(self : WiseBlockType) -> String {
  match self {
    WiseExon => "exon"
    WiseIntron => "intron"
    WiseMatch => "match"
    WiseMismatch => "mismatch"
    WiseInsertion => "insertion"
    WiseDeletion => "deletion"
  }
}

///|
/// Construct a WiseBlockType value by name.
/// Since enum variants cannot be constructed cross-package, this factory
/// allows external callers to obtain block type values for testing/inspection.
pub fn wise_make_block_type(name : String) -> WiseBlockType {
  match name {
    "exon" => WiseExon
    "intron" => WiseIntron
    "match" => WiseMatch
    "mismatch" => WiseMismatch
    "insertion" => WiseInsertion
    "deletion" => WiseDeletion
    _ => WiseMatch
  }
}

///|
/// An exon in a gene prediction.
pub struct WiseExon {
  pub start : Int // gene start position (0-based)
  pub end : Int // gene end position (0-based, inclusive)
  pub phase : Int // reading frame phase (0, 1, or 2)
  pub end_phase : Int // end phase (0, 1, or 2)
  pub score : Double // exon score
  pub protein_start : Int // start position in protein
  pub protein_end : Int // end position in protein
}

///|
/// Construct a WiseExon.
pub fn WiseExon::new(
  start~ : Int = 0,
  end~ : Int = 0,
  phase~ : Int = 0,
  end_phase~ : Int = 0,
  score~ : Double = 0.0,
  protein_start~ : Int = 0,
  protein_end~ : Int = 0,
) -> WiseExon {
  WiseExon::{
    start,
    end,
    phase,
    end_phase,
    score,
    protein_start,
    protein_end,
  }
}

///|
/// Get the length of the exon in nucleotides.
pub fn WiseExon::length(self : WiseExon) -> Int {
  self.end - self.start + 1
}

///|
/// An intron between exons.
pub struct WiseIntron {
  pub start : Int // intron start position
  pub end : Int // intron end position
  pub donor_score : Double // 5' splice site score
  pub acceptor_score : Double // 3' splice site score
  pub length : Int
}

///|
/// Construct a WiseIntron.
pub fn WiseIntron::new(
  start~ : Int = 0,
  end~ : Int = 0,
  donor_score~ : Double = 0.0,
  acceptor_score~ : Double = 0.0,
) -> WiseIntron {
  let len = if end >= start { end - start + 1 } else { 0 }
  WiseIntron::{
    start,
    end,
    donor_score,
    acceptor_score,
    length: len,
  }
}

///|
/// A single position in a GeneWise alignment.
pub struct WiseAlignmentColumn {
  pub protein_char : Char // amino acid in protein
  pub gene_codon : String // 3-nucleotide codon
  pub match_type : String // " " = match, " " = mismatch, "*" = stop
  pub gene_position : Int
  pub protein_position : Int
}

///|
/// Construct a WiseAlignmentColumn.
pub fn WiseAlignmentColumn::new(
  protein_char~ : Char,
  gene_codon~ : String,
  match_type~ : String = " ",
  gene_position~ : Int = 0,
  protein_position~ : Int = 0,
) -> WiseAlignmentColumn {
  WiseAlignmentColumn::{
    protein_char,
    gene_codon,
    match_type,
    gene_position,
    protein_position,
  }
}

///|
/// A complete GeneWise result.
pub struct WiseResult {
  pub mut protein_id : String
  pub mut gene_id : String
  pub mut score : Double // overall alignment score
  pub mut bits_score : Double // score in bits
  pub probability : Double // probability of alignment
  pub protein_sequence : String
  pub gene_sequence : String
  exons : Array[WiseExon]
  introns : Array[WiseIntron]
  alignment : Array[WiseAlignmentColumn]
  pub parameters : Map[String, String]
}

///|
/// Construct an empty WiseResult.
pub fn WiseResult::new() -> WiseResult {
  WiseResult::{
    protein_id: "",
    gene_id: "",
    score: 0.0,
    bits_score: 0.0,
    probability: 0.0,
    protein_sequence: "",
    gene_sequence: "",
    exons: [],
    introns: [],
    alignment: [],
    parameters: Map([], capacity=0),
  }
}

///|
/// Add an exon.
pub fn WiseResult::add_exon(self : WiseResult, e : WiseExon) -> Unit {
  self.exons.push(e)
}

///|
/// Add an intron.
pub fn WiseResult::add_intron(self : WiseResult, i : WiseIntron) -> Unit {
  self.introns.push(i)
}

///|
/// Add an alignment column.
pub fn WiseResult::add_alignment_column(
  self : WiseResult,
  c : WiseAlignmentColumn,
) -> Unit {
  self.alignment.push(c)
}

///|
/// Get all exons.
pub fn WiseResult::get_exons(self : WiseResult) -> Array[WiseExon] {
  self.exons
}

///|
/// Get all introns.
pub fn WiseResult::get_introns(self : WiseResult) -> Array[WiseIntron] {
  self.introns
}

///|
/// Get the alignment.
pub fn WiseResult::get_alignment(self : WiseResult) -> Array[WiseAlignmentColumn] {
  self.alignment
}

///|
/// Get the number of exons.
pub fn WiseResult::get_num_exons(self : WiseResult) -> Int {
  self.exons.length()
}

///|
/// Get the number of introns.
pub fn WiseResult::get_num_introns(self : WiseResult) -> Int {
  self.introns.length()
}

///|
/// Get total gene length covered by exons.
pub fn WiseResult::get_total_exon_length(self : WiseResult) -> Int {
  let mut total = 0
  for e in self.exons {
    total = total + e.length()
  }
  total
}

///|
/// Set the protein ID.
pub fn WiseResult::set_protein_id(self : WiseResult, id : String) -> Unit {
  self.protein_id = id
}

///|
/// Set the gene ID.
pub fn WiseResult::set_gene_id(self : WiseResult, id : String) -> Unit {
  self.gene_id = id
}

///|
/// Set the score.
pub fn WiseResult::set_score(self : WiseResult, s : Double) -> Unit {
  self.score = s
}

///|
/// Set the bits score.
pub fn WiseResult::set_bits_score(self : WiseResult, s : Double) -> Unit {
  self.bits_score = s
}

// ============================================================================
// GeneWise output parsing
// ============================================================================

///|
/// Parse GeneWise text output into a WiseResult.
///
/// Expected format:
///   Lines starting with "Gene:" contain gene info
///   Lines starting with "Protein:" contain protein info
///   Lines starting with "Score:" contain alignment score
///   Lines starting with "Exon" contain exon information
///   Lines starting with "Intron" contain intron information
///   Lines with "Protein" and "Gene" tags contain alignment blocks
pub fn wise_parse(content : String) -> WiseResult {
  let result = WiseResult::new()
  let lines = content.split("\n").to_array()
  for line in lines {
    let trimmed = line.to_string()
    if trimmed.length() == 0 {
      continue
    }
    // Parse score line: "Score: 123.45"
    if trimmed.starts_with("Score:") {
      let parts = trimmed.split(":").to_array()
      if parts.length() >= 2 {
        result.score = parse_double(parts[1].to_string().trim().to_string()).unwrap_or(0.0)
      }
    }
    // Parse bits score: "Bits: 45.67"
    if trimmed.starts_with("Bits:") {
      let parts = trimmed.split(":").to_array()
      if parts.length() >= 2 {
        result.bits_score = parse_double(parts[1].to_string().trim().to_string()).unwrap_or(0.0)
      }
    }
    // Parse protein ID: "Protein: P12345"
    if trimmed.starts_with("Protein:") {
      let parts = trimmed.split(":").to_array()
      if parts.length() >= 2 {
        result.protein_id = parts[1].to_string().trim().to_string()
      }
    }
    // Parse gene ID: "Gene: G67890"
    if trimmed.starts_with("Gene:") {
      let parts = trimmed.split(":").to_array()
      if parts.length() >= 2 {
        result.gene_id = parts[1].to_string().trim().to_string()
      }
    }
    // Parse exon line:
    // "Exon 1: 100-200 (phase 0) score 45.5 protein 1-33"
    if trimmed.starts_with("Exon") {
      let exon = wise_parse_exon_line(trimmed)
      match exon {
        Some(e) => result.add_exon(e)
        None => ()
      }
    }
    // Parse intron line:
    // "Intron 1: 201-300 donor 12.5 acceptor 8.3"
    if trimmed.starts_with("Intron") {
      let intron = wise_parse_intron_line(trimmed)
      match intron {
        Some(i) => result.add_intron(i)
        None => ()
      }
    }
    // Parse parameter lines: "Parameter: value"
    if trimmed.contains(":") && !trimmed.starts_with("Exon") && !trimmed.starts_with("Intron") && !trimmed.starts_with("Score") && !trimmed.starts_with("Bits") && !trimmed.starts_with("Protein") && !trimmed.starts_with("Gene") {
      let parts = trimmed.split(":").to_array()
      if parts.length() >= 2 {
        let key = parts[0].to_string().trim().to_string()
        let value = parts[1].to_string().trim().to_string()
        if key.length() > 0 && value.length() > 0 && !key.contains(" ") {
          result.parameters[key] = value
        }
      }
    }
  }
  result
}

///|
/// Parse an exon line from GeneWise output.
/// Format: "Exon <num>: <start>-<end> (phase <p>) score <s> protein <ps>-<pe>"
fn wise_parse_exon_line(line : String) -> WiseExon? {
  // Extract start-end
  let dash_pos = line.find("-")
  let mut start = 0
  let mut end = 0
  match dash_pos {
    Some(dp) => {
      // Find the start of the number before dash
      let before = line[0:dp].to_string()
      let after = line[dp + 1:].to_string()
      // Parse start number from before (last number)
      start = wise_extract_last_number(before)
      // Parse end number from after (first number)
      end = wise_extract_first_number(after)
    }
    None => return None
  }
  // Extract phase
  let mut phase = 0
  if line.contains("phase ") {
    let phase_str = wise_extract_after(line, "phase ")
    phase = wise_extract_first_number(phase_str)
  }
  // Extract score
  let mut score = 0.0
  if line.contains("score ") {
    let score_str = wise_extract_after(line, "score ")
    score = wise_extract_first_double(score_str)
  }
  // Extract protein start-end
  let mut protein_start = 0
  let mut protein_end = 0
  if line.contains("protein ") {
    let protein_part = wise_extract_after(line, "protein ")
    let dash = protein_part.find("-")
    match dash {
      Some(dp) => {
        protein_start = wise_extract_first_number(protein_part[0:dp].to_string())
        protein_end = wise_extract_first_number(protein_part[dp + 1:].to_string())
      }
      None => ()
    }
  }
  Some(
    WiseExon::new(
      start=start,
      end=end,
      phase=phase,
      score=score,
      protein_start=protein_start,
      protein_end=protein_end,
    ),
  )
}

///|
/// Parse an intron line from GeneWise output.
/// Format: "Intron <num>: <start>-<end> donor <d> acceptor <a>"
fn wise_parse_intron_line(line : String) -> WiseIntron? {
  let dash_pos = line.find("-")
  let mut start = 0
  let mut end = 0
  match dash_pos {
    Some(dp) => {
      let before = line[0:dp].to_string()
      let after = line[dp + 1:].to_string()
      start = wise_extract_last_number(before)
      end = wise_extract_first_number(after)
    }
    None => return None
  }
  let mut donor = 0.0
  let mut acceptor = 0.0
  if line.contains("donor ") {
    let donor_str = wise_extract_after(line, "donor ")
    donor = wise_extract_first_double(donor_str)
  }
  if line.contains("acceptor ") {
    let acc_str = wise_extract_after(line, "acceptor ")
    acceptor = wise_extract_first_double(acc_str)
  }
  Some(
    WiseIntron::new(
      start=start,
      end=end,
      donor_score=donor,
      acceptor_score=acceptor,
    ),
  )
}

///|
/// Extract the last integer from a string.
fn wise_extract_last_number(s : String) -> Int {
  let mut result = 0
  let mut current = ""
  let chars = s.to_array()
  for i in 0..<chars.length() {
    let c = chars[i]
    if c >= '0' && c <= '9' {
      current = current + c.to_string()
    } else {
      if current.length() > 0 {
        result = parse_int(current)
        current = ""
      }
    }
  }
  if current.length() > 0 {
    result = parse_int(current)
  }
  result
}

///|
/// Extract the first integer from a string.
fn wise_extract_first_number(s : String) -> Int {
  let chars = s.to_array()
  let mut current = ""
  for i in 0..<chars.length() {
    let c = chars[i]
    if c >= '0' && c <= '9' {
      current = current + c.to_string()
    } else {
      if current.length() > 0 {
        return parse_int(current)
      }
    }
  }
  if current.length() > 0 {
    return parse_int(current)
  }
  0
}

///|
/// Extract the first double from a string.
fn wise_extract_first_double(s : String) -> Double {
  let chars = s.to_array()
  let mut current = ""
  let mut found_digit = false
  for i in 0..<chars.length() {
    let c = chars[i]
    if (c >= '0' && c <= '9') || c == '.' || c == '-' {
      current = current + c.to_string()
      if c >= '0' && c <= '9' {
        found_digit = true
      }
    } else {
      if found_digit && current.length() > 0 {
        return parse_double(current).unwrap_or(0.0)
      }
      current = ""
    }
  }
  if found_digit && current.length() > 0 {
    return parse_double(current).unwrap_or(0.0)
  }
  0.0
}

///|
/// Extract the substring after a pattern.
fn wise_extract_after(s : String, pattern : String) -> String {
  match s.find(pattern) {
    Some(pos) => s[pos + pattern.length():].to_string()
    None => ""
  }
}

// ============================================================================
// Serialization
// ============================================================================

///|
/// Serialize a WiseResult to GeneWise-like text format.
pub fn wise_write(result : WiseResult) -> String {
  let sb = StringBuilder::new()
  sb.write_string("GeneWise Result\n")
  sb.write_string("Protein: " + result.protein_id + "\n")
  sb.write_string("Gene: " + result.gene_id + "\n")
  sb.write_string("Score: " + result.score.to_string() + "\n")
  sb.write_string("Bits: " + result.bits_score.to_string() + "\n")
  // Parameters
  for key in result.parameters.keys() {
    match result.parameters.get(key) {
      Some(v) => sb.write_string(key + ": " + v + "\n")
      None => ()
    }
  }
  // Exons
  for i in 0..<result.exons.length() {
    let e = result.exons[i]
    sb.write_string(
      "Exon " +
      (i + 1).to_string() +
      ": " +
      e.start.to_string() +
      "-" +
      e.end.to_string() +
      " (phase " +
      e.phase.to_string() +
      ") score " +
      e.score.to_string() +
      " protein " +
      e.protein_start.to_string() +
      "-" +
      e.protein_end.to_string() +
      "\n",
    )
  }
  // Introns
  for i in 0..<result.introns.length() {
    let intron = result.introns[i]
    sb.write_string(
      "Intron " +
      (i + 1).to_string() +
      ": " +
      intron.start.to_string() +
      "-" +
      intron.end.to_string() +
      " donor " +
      intron.donor_score.to_string() +
      " acceptor " +
      intron.acceptor_score.to_string() +
      "\n",
    )
  }
  sb.to_string()
}

// ============================================================================
// Analysis functions
// ============================================================================

///|
/// Get the gene structure as a list of exons and introns.
pub fn wise_gene_structure(result : WiseResult) -> String {
  let sb = StringBuilder::new()
  for i in 0..<result.exons.length() {
    let e = result.exons[i]
    sb.write_string("Exon " + (i + 1).to_string() + ": " + e.start.to_string() + "-" + e.end.to_string() + " (length " + e.length().to_string() + " nt, phase " + e.phase.to_string() + ")\n")
  }
  for i in 0..<result.introns.length() {
    let intron = result.introns[i]
    sb.write_string("Intron " + (i + 1).to_string() + ": " + intron.start.to_string() + "-" + intron.end.to_string() + " (length " + intron.length.to_string() + " nt)\n")
  }
  sb.to_string()
}

///|
/// Get the translated protein from the gene sequence.
pub fn wise_translate_gene(
  result : WiseResult,
  gene_sequence : String,
) -> String {
  let sb = StringBuilder::new()
  for exon in result.exons {
    let exon_seq = gene_sequence.substring(start=exon.start, end=exon.end + 1)
    // Translate codon by codon
    let mut i = 0
    while i + 2 < exon_seq.length() {
      let codon = exon_seq[i:i + 3].to_string()
      let aa = wise_translate_codon(codon)
      sb.write_char(aa)
      i = i + 3
    }
  }
  sb.to_string()
}

///|
/// Translate a single codon to an amino acid.
fn wise_translate_codon(codon : String) -> Char {
  match codon {
    "TTT" | "TTC" => 'F'
    "TTA" | "TTG" | "CTT" | "CTC" | "CTA" | "CTG" => 'L'
    "ATT" | "ATC" | "ATA" => 'I'
    "ATG" => 'M'
    "GTT" | "GTC" | "GTA" | "GTG" => 'V'
    "TCT" | "TCC" | "TCA" | "TCG" => 'S'
    "CCT" | "CCC" | "CCA" | "CCG" => 'P'
    "ACT" | "ACC" | "ACA" | "ACG" => 'T'
    "GCT" | "GCC" | "GCA" | "GCG" => 'A'
    "TAT" | "TAC" => 'Y'
    "TAA" | "TAG" | "TGA" => '*'
    "CAT" | "CAC" => 'H'
    "CAA" | "CAG" => 'Q'
    "AAT" | "AAC" => 'N'
    "AAA" | "AAG" => 'K'
    "GAT" | "GAC" => 'D'
    "GAA" | "GAG" => 'E'
    "TGT" | "TGC" => 'C'
    "TGG" => 'W'
    "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => 'R'
    "AGT" | "AGC" => 'S'
    "GGT" | "GGC" | "GGA" | "GGG" => 'G'
    _ => 'X'
  }
}

///|
/// Calculate the percent identity between protein and translated gene.
pub fn wise_percent_identity(
  result : WiseResult,
  protein : String,
  gene_sequence : String,
) -> Double {
  let translated = wise_translate_gene(result, gene_sequence)
  if translated.length() == 0 || protein.length() == 0 {
    return 0.0
  }
  let min_len = if translated.length() < protein.length() {
    translated.length()
  } else {
    protein.length()
  }
  let mut matches = 0
  for i in 0..<min_len {
    if translated.unsafe_get(i) == protein.unsafe_get(i) {
      matches = matches + 1
    }
  }
  (matches.to_double() / min_len.to_double()) * 100.0
}

// ============================================================================
// GenomeWise parsing
// ============================================================================

///|
/// A GenomeWise segment result.
pub struct GenomeWiseSegment {
  pub segment_id : String
  pub start : Int
  pub end : Int
  pub score : Double
  pub wise_result : WiseResult
}

///|
/// Construct a GenomeWiseSegment.
pub fn GenomeWiseSegment::new(
  segment_id~ : String = "",
  start~ : Int = 0,
  end~ : Int = 0,
  score~ : Double = 0.0,
) -> GenomeWiseSegment {
  GenomeWiseSegment::{
    segment_id,
    start,
    end,
    score,
    wise_result: WiseResult::new(),
  }
}

///|
/// A complete GenomeWise result containing multiple segments.
pub struct GenomeWiseResult {
  pub gene_id : String
  pub total_score : Double
  segments : Array[GenomeWiseSegment]
}

///|
/// Construct an empty GenomeWiseResult.
pub fn GenomeWiseResult::new() -> GenomeWiseResult {
  GenomeWiseResult::{
    gene_id: "",
    total_score: 0.0,
    segments: [],
  }
}

///|
/// Add a segment.
pub fn GenomeWiseResult::add_segment(
  self : GenomeWiseResult,
  s : GenomeWiseSegment,
) -> Unit {
  self.segments.push(s)
}

///|
/// Get all segments.
pub fn GenomeWiseResult::get_segments(self : GenomeWiseResult) -> Array[GenomeWiseSegment] {
  self.segments
}

///|
/// Get the number of segments.
pub fn GenomeWiseResult::get_num_segments(self : GenomeWiseResult) -> Int {
  self.segments.length()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Create a sample GeneWise output for testing.
pub fn wise_sample_output() -> String {
  "GeneWise output\nProtein: P12345\nGene: G67890\nScore: 145.32\nBits: 52.18\nExon 1: 0-99 (phase 0) score 45.5 protein 1-33\nIntron 1: 100-199 donor 12.5 acceptor 8.3\nExon 2: 200-299 (phase 0) score 52.3 protein 34-66\nIntron 2: 300-399 donor 10.2 acceptor 9.1\nExon 3: 400-479 (phase 0) score 47.5 protein 67-93\n"
}

///|
/// Create a sample WiseResult for testing.
pub fn wise_sample() -> WiseResult {
  wise_parse(wise_sample_output())
}

///|
/// Generate a summary string.
pub fn wise_summary(result : WiseResult) -> String {
  let sb = StringBuilder::new()
  sb.write_string("GeneWise Result Summary:\n")
  sb.write_string("  Protein ID: " + result.protein_id + "\n")
  sb.write_string("  Gene ID: " + result.gene_id + "\n")
  sb.write_string("  Score: " + result.score.to_string() + "\n")
  sb.write_string("  Bits: " + result.bits_score.to_string() + "\n")
  sb.write_string("  Exons: " + result.get_num_exons().to_string() + "\n")
  sb.write_string("  Introns: " + result.get_num_introns().to_string() + "\n")
  sb.write_string("  Total exon length: " + result.get_total_exon_length().to_string() + " nt\n")
  sb.write_string("  Parameters: " + result.parameters.size().to_string() + "\n")
  sb.to_string()
}
