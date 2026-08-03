// Bio.Motifs.AlignAce - Parser for AlignACE motif discovery output,
// mirroring BioPython's Bio.Motifs.AlignAce.
//
// AlignACE (Align Nucleic Acid Conserved Elements) is a Gibbs sampler-based
// motif discovery tool. Its output format contains:
//   - Header with version and command-line information
//   - Parameters used for the run
//   - One or more motifs, each represented as a position frequency matrix
//   - Optional site information (sequences and positions matching the motif)
//
// Format example:
//   # AlignACE version 4.0
//   # Program: AlignACE
//   # Command: alignace -i=input.fa -numcols=10
//   #
//   # Parameters used:
//   # numcols=10
//   # expect=10
//   #
//   i    T   G   A   C   T   C   G   A   T
//        0   1   2   3   4   5   6   7   8   9
//   A    1   0   8   0   0   0   0   0   0   0
//   C    0   0   0   0   1   0   8   0   0   1
//   G    0   8   0   0   0   0   0   8   0   0
//   T    7   0   0   8   7   8   0   0   8   7
//   #
//   # Sites:
//   # seq1   10   TGACTCGAT
//   #

// ============================================================================
// Data structures
// ============================================================================

///|
/// A motif discovered by AlignACE, represented as a position frequency matrix.
pub struct AlignAceMotif {
  count_matrix : Array[Array[Int]] // 4 rows (A, C, G, T) × width columns
  width : Int
  num_sites : Int // total number of sites contributing
  consensus : String // consensus sequence derived from the matrix
  sites : Array[AlignAceSite] // individual sites matching this motif
}

///|
/// Construct an AlignAceMotif from a count matrix.
pub fn AlignAceMotif::new(count_matrix : Array[Array[Int]]) -> AlignAceMotif {
  let width = if count_matrix.length() > 0 { count_matrix[0].length() } else { 0 }
  let mut total = 0
  if width > 0 {
    for base_idx in 0..<4 {
      if base_idx < count_matrix.length() {
        total = total + count_matrix[base_idx][0]
      }
    }
  }
  let motif = AlignAceMotif::{
    count_matrix,
    width,
    num_sites: total,
    consensus: alignace_consensus(count_matrix),
    sites: [],
  }
  motif
}

///|
/// A single site matching an AlignACE motif.
pub struct AlignAceSite {
  sequence_id : String
  position : Int // 0-based start position in the original sequence
  strand : String // "+" or "-"
  sequence : String // the actual subsequence matching the motif
}

///|
/// Construct an AlignAceSite.
pub fn AlignAceSite::new(
  sequence_id~ : String,
  position~ : Int,
  strand~ : String,
  sequence~ : String,
) -> AlignAceSite {
  AlignAceSite::{ sequence_id, position, strand, sequence }
}

///|
/// A complete AlignACE output record, containing version, parameters, and motifs.
pub struct AlignAceRecord {
  pub mut version : String
  pub mut command : String
  parameters : Map[String, String]
  motifs : Array[AlignAceMotif]
}

///|
/// Construct an empty AlignAceRecord.
pub fn AlignAceRecord::new() -> AlignAceRecord {
  AlignAceRecord::{
    version: "",
    command: "",
    parameters: Map([], capacity=0),
    motifs: [],
  }
}

///|
/// Set the version string.
pub fn AlignAceRecord::set_version(self : AlignAceRecord, s : String) -> Unit {
  self.version = s
}

///|
/// Set the command string.
pub fn AlignAceRecord::set_command(self : AlignAceRecord, s : String) -> Unit {
  self.command = s
}

///|
/// Add a site to the motif.
pub fn AlignAceMotif::add_site(self : AlignAceMotif, s : AlignAceSite) -> Unit {
  self.sites.push(s)
}

///|
/// Add a motif to the record.
pub fn AlignAceRecord::add_motif(self : AlignAceRecord, m : AlignAceMotif) -> Unit {
  self.motifs.push(m)
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Parse AlignACE output text into an AlignAceRecord.
///
/// Mirrors BioPython's Bio.Motifs.AlignAce.read(handle).
pub fn alignace_parse(text : String) -> AlignAceRecord {
  let record = AlignAceRecord::new()
  // Convert StringView array to String array for indexed access
  let lines_raw = text.split("\n").to_array()
  let lines : Array[String] = Array::new()
  for l in lines_raw {
    lines.push(l.to_string())
  }
  let n = lines.length()
  let mut i = 0
  let mut in_parameters = false
  while i < n {
    let line = lines[i]
    let trimmed = line.trim().to_string()
    // Header and comment lines
    if trimmed.starts_with("#") {
      let content = trimmed[1:].to_string().trim().to_string()
      if content.starts_with("AlignACE version") {
        record.version = content[16:].to_string().trim().to_string()
      } else if content.starts_with("Command:") {
        record.command = content[8:].to_string().trim().to_string()
      } else if content == "Parameters used:" {
        in_parameters = true
      } else if content == "" {
        in_parameters = false
      } else if in_parameters && content.contains("=") {
        let parts = content.split("=").to_array()
        if parts.length() == 2 {
          let key = parts[0].trim().to_string()
          let val = parts[1].trim().to_string()
          record.parameters[key] = val
        }
      } else if content.starts_with("Sites:") {
        // Next lines will be sites for the current motif
        i = i + 1
        if record.motifs.length() > 0 {
          let last_idx = record.motifs.length() - 1
          while i < n {
            let sline = lines[i].trim().to_string()
            if sline.starts_with("#") {
              let scontent = sline[1:].to_string().trim().to_string()
              if scontent == "" || scontent.starts_with("Sites:") {
                break
              }
              // Parse site line: "seq_id position strand sequence"
              let site = alignace_parse_site(scontent)
              match site {
                Some(s) => record.motifs[last_idx].sites.push(s)
                None => ()
              }
            } else {
              break
            }
            i = i + 1
          }
        }
        continue
      }
      i = i + 1
      continue
    }
    // Motif matrix starts with a line beginning with "i"
    if trimmed.starts_with("i ") || trimmed.starts_with("i\t") {
      let motif = alignace_parse_motif(lines, i)
      match motif {
        Some(m) => {
          record.motifs.push(m)
          // Advance past the matrix (4 data rows + header rows + "#")
          i = i + 7 // "i" line, position line, 4 base rows, "#"
        }
        None => i = i + 1
      }
      continue
    }
    i = i + 1
  }
  record
}

///|
/// Parse a motif matrix starting at the given line index.
/// Returns Some(motif) if successful, None otherwise.
fn alignace_parse_motif(
  lines : Array[String],
  start : Int,
) -> AlignAceMotif? {
  // Line 0: "i   T   G   A   C   T   C   G   A   T" (consensus letters)
  // Line 1: "    0   1   2   3   4   5   6   7   8   9" (column indices)
  // Line 2: "A   1   0   8   0   0   0   0   0   0   0"
  // Line 3: "C   0   0   0   0   1   0   8   0   0   1"
  // Line 4: "G   0   8   0   0   0   0   0   8   0   0"
  // Line 5: "T   7   0   0   8   7   8   0   0   8   7"
  let base_rows = ["A", "C", "G", "T"]
  let count_matrix : Array[Array[Int]] = Array::new()
  for row_idx in 0..<4 {
    let line_idx = start + 2 + row_idx
    if line_idx >= lines.length() {
      return None
    }
    let line = lines[line_idx].trim().to_string()
    // Parse "A 1 0 8 0 ..."
    let parts = alignace_split_whitespace(line)
    if parts.length() < 2 {
      return None
    }
    // First part should be the base letter
    if parts[0].to_upper() != base_rows[row_idx] {
      return None
    }
    let counts = Array::new()
    for j in 1..<parts.length() {
      counts.push(parse_int(parts[j]))
    }
    count_matrix.push(counts)
  }
  Some(AlignAceMotif::new(count_matrix))
}

///|
/// Parse a site line: "seq_id position strand sequence"
/// or "seq_id position sequence"
fn alignace_parse_site(text : String) -> AlignAceSite? {
  let parts = alignace_split_whitespace(text)
  if parts.length() < 3 {
    return None
  }
  let seq_id = parts[0]
  let pos = parse_int(parts[1])
  // Check if there's a strand field
  if parts.length() >= 4 && (parts[2] == "+" || parts[2] == "-") {
    let strand = parts[2]
    let sequence = parts[3]
    return Some(
      AlignAceSite::new(
        sequence_id=seq_id,
        position=pos,
        strand=strand,
        sequence=sequence,
      ),
    )
  }
  // No strand field, assume "+"
  let sequence = parts[2]
  Some(
    AlignAceSite::new(
      sequence_id=seq_id,
      position=pos,
      strand="+",
      sequence=sequence,
    ),
  )
}

///|
/// Split a string by whitespace, returning non-empty parts.
fn alignace_split_whitespace(s : String) -> Array[String] {
  let parts = Array::new()
  let sb = StringBuilder::new()
  for c in s.to_array() {
    if c == ' ' || c == '\t' {
      let built = sb.to_string()
      if built != "" {
        parts.push(built)
      }
      sb.reset()
    } else {
      sb.write_char(c)
    }
  }
  let last = sb.to_string()
  if last != "" {
    parts.push(last)
  }
  parts
}

// ============================================================================
// Consensus and matrix operations
// ============================================================================

///|
/// Compute the consensus sequence from a count matrix.
/// Uses the base with the highest count at each position.
fn alignace_consensus(count_matrix : Array[Array[Int]]) -> String {
  if count_matrix.length() < 4 {
    return ""
  }
  let width = count_matrix[0].length()
  let bases = ["A", "C", "G", "T"]
  let sb = StringBuilder::new()
  for col in 0..<width {
    let mut best_base = "A"
    let mut best_count = count_matrix[0][col]
    for row in 1..<4 {
      if count_matrix[row][col] > best_count {
        best_count = count_matrix[row][col]
        best_base = bases[row]
      }
    }
    sb.write_string(best_base)
  }
  sb.to_string()
}

///|
/// Get the consensus sequence for an AlignAceMotif.
pub fn alignace_motif_consensus(motif : AlignAceMotif) -> String {
  motif.consensus
}

///|
/// Convert a count matrix to a PWM (position weight matrix) by normalizing
/// counts to frequencies.
pub fn alignace_to_pwm(motif : AlignAceMotif) -> Array[Array[Double]] {
  let width = motif.width
  let pwm : Array[Array[Double]] = Array::new()
  for row in 0..<4 {
    let freq_row = Array::make(width, 0.0)
    for col in 0..<width {
      let total = motif.count_matrix[0][col] +
        motif.count_matrix[1][col] +
        motif.count_matrix[2][col] +
        motif.count_matrix[3][col]
      if total > 0 {
        freq_row[col] = motif.count_matrix[row][col].to_double() / total.to_double()
      } else {
        freq_row[col] = 0.25
      }
    }
    pwm.push(freq_row)
  }
  pwm
}

///|
/// Compute the information content (bits) at each position of the motif.
/// IC(i) = 2 + Σ_b f(b,i) * log2(f(b,i))
pub fn alignace_information_content(motif : AlignAceMotif) -> Array[Double] {
  let pwm = alignace_to_pwm(motif)
  let width = motif.width
  let ic = Array::make(width, 0.0)
  for col in 0..<width {
    let mut sum = 0.0
    for row in 0..<4 {
      let f = pwm[row][col]
      if f > 0.0 {
        sum = sum + f * (@math.ln(f) / 0.6931471805599453)
      }
    }
    ic[col] = 2.0 + sum
  }
  ic
}

///|
/// Total information content of the motif (sum of per-position IC).
pub fn alignace_total_ic(motif : AlignAceMotif) -> Double {
  let ic = alignace_information_content(motif)
  let mut total = 0.0
  for v in ic {
    total = total + v
  }
  total
}

// ============================================================================
// Serialization
// ============================================================================

///|
/// Serialize an AlignAceRecord to AlignACE output format text.
pub fn alignace_write(record : AlignAceRecord) -> String {
  let sb = StringBuilder::new()
  // Header
  sb.write_string("# AlignACE version " + record.version + "\n")
  sb.write_string("# Program: AlignACE\n")
  sb.write_string("# Command: " + record.command + "\n")
  sb.write_string("#\n")
  // Parameters
  sb.write_string("# Parameters used:\n")
  let keys = record.parameters.keys()
  for k in keys {
    sb.write_string("# " + k + "=" + record.parameters[k] + "\n")
  }
  sb.write_string("#\n")
  // Motifs
  for motif in record.motifs {
    sb.write_string(alignace_write_motif(motif))
    sb.write_string("\n")
  }
  sb.to_string()
}

///|
/// Serialize a single motif to AlignACE matrix format.
fn alignace_write_motif(motif : AlignAceMotif) -> String {
  let sb = StringBuilder::new()
  let width = motif.width
  // Consensus line
  sb.write_string("i")
  for col in 0..<width {
    let mut best_base = 'A'
    let mut best_count = motif.count_matrix[0][col]
    for row in 1..<4 {
      if motif.count_matrix[row][col] > best_count {
        best_count = motif.count_matrix[row][col]
        best_base = match row {
          0 => 'A'
          1 => 'C'
          2 => 'G'
          _ => 'T'
        }
      }
    }
    sb.write_string("\t" + best_base.to_string())
  }
  sb.write_string("\n")
  // Position indices
  sb.write_string(" ")
  for col in 0..<width {
    sb.write_string("\t" + col.to_string())
  }
  sb.write_string("\n")
  // Base rows
  let bases = ["A", "C", "G", "T"]
  for row in 0..<4 {
    sb.write_string(bases[row])
    for col in 0..<width {
      sb.write_string("\t" + motif.count_matrix[row][col].to_string())
    }
    sb.write_string("\n")
  }
  sb.write_string("#\n")
  // Sites
  if motif.sites.length() > 0 {
    sb.write_string("# Sites:\n")
    for site in motif.sites {
      sb.write_string(
        "# " + site.sequence_id + "\t" + site.position.to_string() + "\t" +
        site.strand + "\t" + site.sequence + "\n",
      )
    }
    sb.write_string("#\n")
  }
  sb.to_string()
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Count the number of motifs in a record.
pub fn alignace_num_motifs(record : AlignAceRecord) -> Int {
  record.motifs.length()
}

///|
/// Get a motif by index (0-based). Returns None if out of range.
pub fn alignace_get_motif(record : AlignAceRecord, index : Int) -> AlignAceMotif? {
  if index >= 0 && index < record.motifs.length() {
    Some(record.motifs[index])
  } else {
    None
  }
}

///|
/// Generate a summary string for an AlignAceRecord.
pub fn alignace_summary(record : AlignAceRecord) -> String {
  let sb = StringBuilder::new()
  sb.write_string("AlignACE Record Summary:\n")
  sb.write_string("  Version: " + record.version + "\n")
  sb.write_string("  Command: " + record.command + "\n")
  sb.write_string("  Parameters: " + record.parameters.size().to_string() + "\n")
  let keys = record.parameters.keys()
  for k in keys {
    sb.write_string("    " + k + " = " + record.parameters[k] + "\n")
  }
  sb.write_string("  Motifs: " + record.motifs.length().to_string() + "\n")
  for i in 0..<record.motifs.length() {
    let m = record.motifs[i]
    sb.write_string(
      "    Motif " + (i + 1).to_string() + ": width=" + m.width.to_string() +
      " sites=" + m.num_sites.to_string() + " IC=" +
      alignace_total_ic(m).to_string() + "\n",
    )
    sb.write_string("      Consensus: " + m.consensus + "\n")
  }
  sb.to_string()
}

///|
/// Create a sample AlignAceRecord for testing/demonstration.
pub fn alignace_sample() -> AlignAceRecord {
  let record = AlignAceRecord::new()
  record.version = "4.0"
  record.command = "alignace -i=input.fa -numcols=10 -expect=10"
  record.parameters["numcols"] = "10"
  record.parameters["expect"] = "10"
  // Motif 1: TGACTCGAT
  let m1 = AlignAceMotif::new(
    [
      [1, 0, 8, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 1, 0, 8, 0, 0, 1],
      [0, 8, 0, 0, 0, 0, 0, 8, 0, 0],
      [7, 0, 0, 8, 7, 8, 0, 0, 8, 7],
    ],
  )
  m1.sites.push(
    AlignAceSite::new(
      sequence_id="seq1",
      position=10,
      strand="+",
      sequence="TGACTCGATA",
    ),
  )
  m1.sites.push(
    AlignAceSite::new(
      sequence_id="seq2",
      position=25,
      strand="+",
      sequence="TGACTCGATA",
    ),
  )
  record.motifs.push(m1)
  // Motif 2: AATAAACAAA
  let m2 = AlignAceMotif::new(
    [
      [8, 8, 1, 8, 8, 8, 0, 8, 8, 8],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
      [1, 1, 8, 1, 1, 1, 9, 1, 1, 1],
    ],
  )
  record.motifs.push(m2)
  record
}
