// Bio.motifs.mast - MAST (Motif Alignment and Search Tool) output parsing.
// Inspired by Biopython's Bio.motifs.mast module.
//
// Parses MAST text output files containing motif search results. MAST searches
// sequences for matches to a set of motifs and reports the matches with
// p-values and E-values. The text output includes:
//   - Database and motif information
//   - Sequence matches with combined p-values and E-values
//   - Motif diagrams showing match positions
//   - Individual motif match details
//
// Reference:
//   Bailey TL, Gribskov M. (1998) Combining evidence using p-values:
//   application to sequence homology searches. Bioinformatics 14:48-54.

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// A single motif match at a specific position in a sequence.
pub struct MastMotifMatch {
  pub motif_id : String // motif identifier
  pub position : Int // start position (1-based)
  pub length : Int // motif length
  pub pvalue : Double // p-value of the match
  pub score : Double // match score
  pub strand : String // strand ("+" or "-")
  pub sequence : String // matched sequence substring
}

///|
/// Create a new MastMotifMatch.
pub fn MastMotifMatch::new(
  motif_id : String,
  position : Int,
  pvalue : Double,
) -> MastMotifMatch {
  MastMotifMatch::{
    motif_id,
    position,
    length: 0,
    pvalue,
    score: 0.0,
    strand: "+",
    sequence: "",
  }
}

///|
/// A sequence with its motif match results from MAST.
pub struct MastMotifSeq {
  pub name : String // sequence name
  pub mut description : String // sequence description
  pub mut length : Int // sequence length
  pub mut combined_pvalue : Double // combined p-value
  pub mut combined_evalue : Double // combined E-value
  pub matches : Array[MastMotifMatch] // individual motif matches
  pub diagram : String // motif diagram string
}

///|
/// Create a new MastMotifSeq.
pub fn MastMotifSeq::new(name : String) -> MastMotifSeq {
  MastMotifSeq::{
    name,
    description: "",
    length: 0,
    combined_pvalue: 1.0,
    combined_evalue: 0.0,
    matches: Array::new(),
    diagram: "",
  }
}

///|
/// Add a motif match to the sequence.
pub fn MastMotifSeq::add_match(
  self : MastMotifSeq,
  match_ : MastMotifMatch,
) -> Unit {
  self.matches.push(match_)
}

///|
/// Get the number of motif matches.
pub fn MastMotifSeq::num_matches(self : MastMotifSeq) -> Int {
  self.matches.length()
}

///|
/// Get the best (lowest p-value) motif match.
pub fn MastMotifSeq::best_match(self : MastMotifSeq) -> MastMotifMatch? {
  if self.matches.length() == 0 {
    return None
  }
  let mut best = self.matches[0]
  for m in self.matches {
    if m.pvalue < best.pvalue {
      best = m
    }
  }
  Some(best)
}

///|
/// A motif used in the MAST search.
pub struct MastMotifInfo {
  pub id : String // motif identifier
  pub name : String // motif name
  pub width : Int // motif width
  pub alphabet : String // alphabet used
}

///|
/// Create a new MastMotifInfo.
pub fn MastMotifInfo::new(id : String, width : Int) -> MastMotifInfo {
  MastMotifInfo::{ id, name: "", width, alphabet: "ACGT" }
}

///|
/// A MAST record holding all search results.
pub struct MastMotifRecord {
  pub mut version : String // MAST version
  pub mut database : String // database name
  pub motifs : Array[MastMotifInfo] // motifs used in search
  pub sequences : Array[MastMotifSeq] // sequences with matches
  pub mut motif_count : Int // number of motifs
  pub mut sequence_count : Int // number of sequences searched
  pub mut hit_count : Int // number of sequences with hits
}

///|
/// Create a new empty MastMotifRecord.
pub fn MastMotifRecord::new() -> MastMotifRecord {
  MastMotifRecord::{
    version: "",
    database: "",
    motifs: Array::new(),
    sequences: Array::new(),
    motif_count: 0,
    sequence_count: 0,
    hit_count: 0,
  }
}

///|
/// Get the number of sequences with significant hits.
pub fn MastMotifRecord::num_hits(self : MastMotifRecord) -> Int {
  self.sequences.length()
}

///|
/// Find a sequence by name.
pub fn MastMotifRecord::find_sequence(
  self : MastMotifRecord,
  name : String,
) -> MastMotifSeq? {
  for seq in self.sequences {
    if seq.name == name {
      return Some(seq)
    }
  }
  None
}

///|
/// Get all matches for a specific motif across all sequences.
pub fn MastMotifRecord::motif_matches(
  self : MastMotifRecord,
  motif_id : String,
) -> Array[MastMotifMatch] {
  let result = Array::new()
  for seq in self.sequences {
    for m in seq.matches {
      if m.motif_id == motif_id {
        result.push(m)
      }
    }
  }
  result
}

///|
/// Get the sequence with the best (lowest) E-value.
pub fn MastMotifRecord::best_sequence(self : MastMotifRecord) -> MastMotifSeq? {
  if self.sequences.length() == 0 {
    return None
  }
  let mut best = self.sequences[0]
  for seq in self.sequences {
    if seq.combined_evalue < best.combined_evalue {
      best = seq
    }
  }
  Some(best)
}

// ===========================================================================
// Parsing
// ===========================================================================

///|
/// Parse MAST text output into a MastMotifRecord.
///
/// Expected format sections:
///   * MAST version line
///   * DATABASE section
///   * MOTIF section with motif IDs and widths
///   * SECTION III (detailed motif matches per sequence)
pub fn mast_motif_parse(text : String) -> MastMotifRecord {
  let record = MastMotifRecord::new()
  let lines = text.split("\n").to_array()
  let mut i = 0
  let n = lines.length()
  let mut section = "" // current section name
  while i < n {
    let line = lines[i]
    let trimmed = line.trim().to_string()
    if trimmed.length() == 0 {
      i = i + 1
      continue
    }
    // MAST version
    if trimmed.starts_with("MAST") && trimmed.contains("version") {
      let parts = trimmed.split(" ").to_array()
      for p in parts {
        let s = p.trim().to_string()
        if s.length() > 0 {
          if s.starts_with("v") || (s[0:1].to_string() >= "0" && s[0:1].to_string() <= "9") {
            record.version = s
          }
        }
      }
      i = i + 1
      continue
    }
    // Database line
    if trimmed.starts_with("DATABASE") || trimmed.contains("database") {
      if trimmed.contains("=") {
        let parts = trimmed.split("=").to_array()
        if parts.length() >= 2 {
          record.database = parts[1].trim().to_string()
        }
      } else {
        let parts = trimmed.split(" ").to_array()
        if parts.length() >= 2 {
          record.database = parts[1].trim().to_string()
        }
      }
      i = i + 1
      continue
    }
    // Section header
    if trimmed.starts_with("SECTION") {
      section = trimmed
      i = i + 1
      continue
    }
    // Motif info lines (in MOTIF section)
    if section.contains("MOTIF") || section.contains("I") {
      if trimmed.starts_with("Motif") || trimmed.contains("width=") {
        let parts = trimmed.split(" ").to_array()
        let mut motif_id = ""
        let mut width = 0
        for p in parts {
          let s = p.trim().to_string()
          if s.starts_with("width=") {
            let w_str = s.substring(start=6)
            width = parse_int(w_str)
          } else if s.starts_with("Motif") {
            motif_id = s.substring(start=5)
          }
        }
        if motif_id.length() > 0 {
          record.motifs.push(MastMotifInfo::new(motif_id, width))
        }
      }
    }
    // Sequence match lines (in SECTION III or similar)
    if section.contains("III") || section.contains("RESULT") || section.contains(
        "DETAILED",
      ) {
      // Parse sequence lines with p-values
      // Format: "sequence_name description length p-value e-value"
      let parts = trimmed.split(" ").to_array()
      let filtered : Array[String] = Array::new()
      for p in parts {
        let s = p.trim().to_string()
        if s.length() > 0 {
          filtered.push(s)
        }
      }
      if filtered.length() >= 4 {
        let name = filtered[0]
        // Try to parse p-value and e-value from the end
        let last_idx = filtered.length() - 1
        let evalue_str = filtered[last_idx]
        let pvalue_str = filtered[last_idx - 1]
        let evalue_opt = parse_double(evalue_str)
        let pvalue_opt = parse_double(pvalue_str)
        let pvalue = match pvalue_opt {
          Some(v) => v
          None => 0.0
        }
        let evalue = match evalue_opt {
          Some(v) => v
          None => 0.0
        }
        if pvalue > 0.0 && pvalue <= 1.0 {
          let seq = MastMotifSeq::new(name)
          seq.combined_pvalue = pvalue
          seq.combined_evalue = evalue
          if filtered.length() >= 3 {
            let len_str = filtered[last_idx - 2]
            seq.length = parse_int(len_str)
          }
          if filtered.length() >= 3 {
            let mut desc = ""
            let mut j = 1
            while j < last_idx - 2 {
              if desc.length() > 0 {
                desc = desc + " "
              }
              desc = desc + filtered[j]
              j = j + 1
            }
            seq.description = desc
          }
          record.sequences.push(seq)
        }
      }
    }
    i = i + 1
  }
  record.motif_count = record.motifs.length()
  record.hit_count = record.sequences.length()
  record
}

// ===========================================================================
// Serialization
// ===========================================================================

///|
/// Serialize a MastMotifMatch to string.
pub fn mast_motif_match_to_string(m : MastMotifMatch) -> String {
  "MotifMatch(motif=" +
  m.motif_id +
  ", pos=" +
  m.position.to_string() +
  ", pvalue=" +
  m.pvalue.to_string() +
  ", strand=" +
  m.strand +
  ")"
}

///|
/// Serialize a MastMotifSeq to string.
pub fn mast_motif_seq_to_string(seq : MastMotifSeq) -> String {
  let sb = StringBuilder::new()
  sb.write_string("Sequence: " + seq.name)
  if seq.description.length() > 0 {
    sb.write_string(" (" + seq.description + ")")
  }
  sb.write_string(
    "\n  Length: " +
    seq.length.to_string() +
    ", p-value: " +
    seq.combined_pvalue.to_string() +
    ", E-value: " +
    seq.combined_evalue.to_string(),
  )
  sb.write_string("\n  Matches: " + seq.matches.length().to_string())
  for m in seq.matches {
    sb.write_string("\n    " + mast_motif_match_to_string(m))
  }
  if seq.diagram.length() > 0 {
    sb.write_string("\n  Diagram: " + seq.diagram)
  }
  sb.to_string()
}

///|
/// Serialize a MastMotifRecord to string summary.
pub fn mast_motif_record_to_string(record : MastMotifRecord) -> String {
  let sb = StringBuilder::new()
  sb.write_string("MAST Record\n")
  sb.write_string("  Version: " + record.version + "\n")
  sb.write_string("  Database: " + record.database + "\n")
  sb.write_string("  Motifs: " + record.motif_count.to_string() + "\n")
  sb.write_string("  Sequences with hits: " + record.hit_count.to_string() + "\n")
  sb.write_string("\nMotif summary:\n")
  for m in record.motifs {
    sb.write_string("  " + m.id + " (width=" + m.width.to_string() + ")\n")
  }
  sb.write_string("\nSequence results:\n")
  for seq in record.sequences {
    sb.write_string(mast_motif_seq_to_string(seq) + "\n")
  }
  sb.to_string()
}

// ===========================================================================
// Sample data
// ===========================================================================

///|
/// Create a sample MAST text output for testing.
pub fn mast_motif_sample_text() -> String {
  "MAST version 5.0.3\n\nDATABASE AND MOTifs\nDATABASE= sample_db.fasta\nMotifs= 3\n\nMOTIF SUMMARY\n\n  Motif 1 width= 10\n  Motif 2 width= 8\n  Motif 3 width= 12\n\nSECTION III: DETAILED RESULTS\n\nSEQUENCE NAME         LENGTH  COMBINED P-VALUE  COMBINED E-VALUE\nseq1 description1      500     1.2e-05           3.4e-03\nseq2 description2      750     5.6e-04           1.1e-01\nseq3 description3      300     2.3e-08           4.5e-06\n"
}

///|
/// Create a sample MastMotifRecord for testing.
pub fn mast_motif_sample_record() -> MastMotifRecord {
  let record = MastMotifRecord::new()
  record.version = "5.0.3"
  record.database = "sample_db.fasta"
  record.motifs.push(MastMotifInfo::new("1", 10))
  record.motifs.push(MastMotifInfo::new("2", 8))
  record.motifs.push(MastMotifInfo::new("3", 12))
  let seq1 = MastMotifSeq::new("seq1")
  seq1.description = "test sequence 1"
  seq1.length = 500
  seq1.combined_pvalue = 1.2e-05
  seq1.combined_evalue = 3.4e-03
  seq1.add_match(MastMotifMatch::new("1", 50, 0.0001))
  seq1.add_match(MastMotifMatch::new("2", 200, 0.0005))
  record.sequences.push(seq1)
  let seq2 = MastMotifSeq::new("seq2")
  seq2.description = "test sequence 2"
  seq2.length = 750
  seq2.combined_pvalue = 5.6e-04
  seq2.combined_evalue = 1.1e-01
  seq2.add_match(MastMotifMatch::new("3", 300, 0.001))
  record.sequences.push(seq2)
  record.motif_count = 3
  record.sequence_count = 2
  record.hit_count = 2
  record
}
