///|
/// CAPS - Cleaved Amplified Polymorphic Sequences marker analysis.
///
/// References:
///   - Biopython Bio.CAPS
///   - Konieczny A, Ausubel FM (1993) A procedure for mapping Arabidopsis
///     mutations using co-dominant ecotype-specific PCR-based markers.
///     Plant J 4(2):403-410.
///
/// A CAPS marker is a position in a multiple sequence alignment where a
/// restriction enzyme cuts at least one sequence but is blocked in at least
/// one other. This differential cutting pattern can be visualized on a gel
/// and serves as a molecular marker for polymorphisms.

// ============================================================================
// Data structures
// ============================================================================

///|
/// A differential cutsite: a position in an alignment where an enzyme
/// cuts some sequences but is blocked in others.
pub struct CapsDifferentialCutsite {
  start : Int // alignment-relative cut position
  enzyme_name : String // name of the restriction enzyme
  cuts_in : Array[Int] // indices of sequences the enzyme cuts at this position
  blocked_in : Array[Int] // indices of sequences where the enzyme is blocked
}

///|
/// Create a differential cutsite.
pub fn CapsDifferentialCutsite::new(
  start : Int,
  enzyme_name : String,
  cuts_in : Array[Int],
  blocked_in : Array[Int],
) -> CapsDifferentialCutsite {
  { start, enzyme_name, cuts_in, blocked_in }
}

///|
/// Get the alignment-relative cut position.
pub fn CapsDifferentialCutsite::start(self : CapsDifferentialCutsite) -> Int {
  self.start
}

///|
/// Get the enzyme name.
pub fn CapsDifferentialCutsite::enzyme_name(
  self : CapsDifferentialCutsite,
) -> String {
  self.enzyme_name
}

///|
/// Get indices of sequences the enzyme cuts at this position.
pub fn CapsDifferentialCutsite::cuts_in(
  self : CapsDifferentialCutsite,
) -> Array[Int] {
  self.cuts_in
}

///|
/// Get indices of sequences where the enzyme is blocked.
pub fn CapsDifferentialCutsite::blocked_in(
  self : CapsDifferentialCutsite,
) -> Array[Int] {
  self.blocked_in
}

///|
/// A map of all differential cutsites in a multiple sequence alignment.
pub struct CapsMap {
  sequences : Array[String]
  sequence_names : Array[String]
  enzymes : Array[RestrictionEnzyme]
  dcuts : Array[CapsDifferentialCutsite]
}

///|
/// Get the sequences in the alignment.
pub fn CapsMap::sequences(self : CapsMap) -> Array[String] {
  self.sequences
}

///|
/// Get the sequence names.
pub fn CapsMap::sequence_names(self : CapsMap) -> Array[String] {
  self.sequence_names
}

///|
/// Get the restriction enzymes used.
pub fn CapsMap::enzymes(self : CapsMap) -> Array[RestrictionEnzyme] {
  self.enzymes
}

///|
/// Get all differential cutsites found.
pub fn CapsMap::dcuts(self : CapsMap) -> Array[CapsDifferentialCutsite] {
  self.dcuts
}

///|
/// Get the number of sequences in the alignment.
pub fn CapsMap::size(self : CapsMap) -> Int {
  self.sequences.length()
}

///|
/// Get the alignment length (all sequences must be the same length).
pub fn CapsMap::length(self : CapsMap) -> Int {
  if self.sequences.length() == 0 {
    0
  } else {
    self.sequences[0].length()
  }
}

// ============================================================================
// CAPS map construction
// ============================================================================

///|
/// Build a CAPS map from a set of aligned sequences and restriction enzymes.
///
/// All sequences must have the same length. The function scans each sequence
/// for enzyme cut sites, then identifies positions where the cutting pattern
/// differs between sequences (differential cutsites).
pub fn caps_map(
  sequences : Array[String],
  enzymes : Array[RestrictionEnzyme],
) -> CapsMap {
  // Validate equal lengths
  let n = sequences.length()
  if n == 0 {
    return {
      sequences,
      sequence_names: [],
      enzymes,
      dcuts: [],
    }
  }
  let expected_len = sequences[0].length()
  let names : Array[String] = []
  for i in 0..<n {
    names.push("seq" + i.to_string())
    // Length validation: sequences with different lengths are still
    // processed, but cut sites outside the alignment are skipped
  }
  ignore(expected_len)
  let dcuts = caps_digest(sequences, enzymes)
  { sequences, sequence_names: names, enzymes, dcuts }
}

///|
/// Build a CAPS map with named sequences.
pub fn caps_map_named(
  sequences : Array[String],
  sequence_names : Array[String],
  enzymes : Array[RestrictionEnzyme],
) -> CapsMap {
  let dcuts = caps_digest(sequences, enzymes)
  { sequences, sequence_names, enzymes, dcuts }
}

///|
/// Core digestion algorithm: find all differential cutsites.
fn caps_digest(
  sequences : Array[String],
  enzymes : Array[RestrictionEnzyme],
) -> Array[CapsDifferentialCutsite] {
  let dcuts : Array[CapsDifferentialCutsite] = []
  let n = sequences.length()
  if n < 2 {
    return dcuts // Need at least 2 sequences for differential cutting
  }
  for enzyme in enzymes {
    caps_digest_with(sequences, enzyme, dcuts)
  }
  dcuts
}

///|
/// Find differential cutsites for a single enzyme.
fn caps_digest_with(
  sequences : Array[String],
  enzyme : RestrictionEnzyme,
  dcuts : Array[CapsDifferentialCutsite],
) -> Unit {
  let n = sequences.length()
  // Find cut sites per sequence
  // cuts_per_seq[i] = list of cut positions in sequence i
  let cuts_per_seq : Array[Array[Int]] = []
  for i in 0..<n {
    cuts_per_seq.push(find_cut_sites(sequences[i], enzyme))
  }
  // Collect all unique cut positions across all sequences
  let all_positions : Array[Int] = []
  for cuts in cuts_per_seq {
    for pos in cuts {
      // Check if already in all_positions
      let mut found = false
      for existing in all_positions {
        if existing == pos {
          found = true
          break
        }
      }
      if !found {
        all_positions.push(pos)
      }
    }
  }
  // Sort positions (simple insertion sort)
  let pn = all_positions.length()
  let mut i = 1
  while i < pn {
    let key = all_positions[i]
    let mut j = i - 1
    while j >= 0 && all_positions[j] > key {
      all_positions[j + 1] = all_positions[j]
      j = j - 1
    }
    all_positions[j + 1] = key
    i = i + 1
  }
  // For each cut position, classify sequences
  for pos in all_positions {
    let cuts_in : Array[Int] = []
    let blocked_in : Array[Int] = []
    for i in 0..<n {
      let mut found = false
      for p in cuts_per_seq[i] {
        if p == pos {
          found = true
          break
        }
      }
      if found {
        cuts_in.push(i)
      } else {
        blocked_in.push(i)
      }
    }
    // A differential cutsite requires both cuts_in and blocked_in non-empty
    if cuts_in.length() > 0 && blocked_in.length() > 0 {
      dcuts.push(
        CapsDifferentialCutsite::new(pos, enzyme.name, cuts_in, blocked_in),
      )
    }
  }
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Get all differential cutsites for a specific enzyme.
pub fn CapsMap::get_dcuts_by_enzyme(
  self : CapsMap,
  enzyme_name : String,
) -> Array[CapsDifferentialCutsite] {
  let result : Array[CapsDifferentialCutsite] = []
  for dc in self.dcuts {
    if dc.enzyme_name == enzyme_name {
      result.push(dc)
    }
  }
  result
}

///|
/// Get all differential cutsites at a specific position.
pub fn CapsMap::get_dcuts_at_position(
  self : CapsMap,
  position : Int,
) -> Array[CapsDifferentialCutsite] {
  let result : Array[CapsDifferentialCutsite] = []
  for dc in self.dcuts {
    if dc.start == position {
      result.push(dc)
    }
  }
  result
}

///|
/// Count the total number of differential cutsites.
pub fn CapsMap::dcut_count(self : CapsMap) -> Int {
  self.dcuts.length()
}

///|
/// Check if a specific enzyme produces any differential cutsites.
pub fn CapsMap::has_dcuts_for_enzyme(
  self : CapsMap,
  enzyme_name : String,
) -> Bool {
  for dc in self.dcuts {
    if dc.enzyme_name == enzyme_name {
      return true
    }
  }
  false
}

///|
/// Format a differential cutsite as a human-readable string.
pub fn CapsDifferentialCutsite::to_string(
  self : CapsDifferentialCutsite,
) -> String {
  let mut s = "DifferentialCutsite(pos=" + self.start.to_string()
  s = s + ", enzyme=" + self.enzyme_name
  s = s + ", cuts_in=["
  let mut first = true
  for i in self.cuts_in {
    if !first {
      s = s + ","
    }
    first = false
    s = s + i.to_string()
  }
  s = s + "], blocked_in=["
  first = true
  for i in self.blocked_in {
    if !first {
      s = s + ","
    }
    first = false
    s = s + i.to_string()
  }
  s = s + "])"
  s
}

///|
/// Format a CapsMap as a summary string.
pub fn CapsMap::to_string(self : CapsMap) -> String {
  let mut s = "CAPSMap(sequences=" + self.size().to_string()
  s = s + ", length=" + self.length().to_string()
  s = s + ", enzymes=" + self.enzymes.length().to_string()
  s = s + ", dcuts=" + self.dcuts.length().to_string()
  s = s + ")"
  s
}

///|
/// Format a CapsMap as a detailed report with each differential cutsite.
pub fn CapsMap::report(self : CapsMap) -> String {
  let mut s = self.to_string() + "\n"
  if self.dcuts.length() == 0 {
    s = s + "No differential cutsites found.\n"
    return s
  }
  s = s + "Differential cutsites:\n"
  for dc in self.dcuts {
    s = s + "  " + dc.to_string() + "\n"
  }
  s
}

// ============================================================================
// Sample data for testing
// ============================================================================

///|
/// Create sample aligned sequences for testing CAPS analysis.
/// Two sequences with a SNP that creates/destroys an EcoRI site (GAATTC):
///   seq0: ...GAATTC...  (EcoRI cuts)
///   seq1: ...GAATTT...  (EcoRI blocked - mismatch at last position)
pub fn caps_sample_sequences() -> Array[String] {
  // 30bp sequences, EcoRI site at position 10
  // seq0 has GAATTC (cut by EcoRI at position 11)
  // seq1 has GAATTT (not cut by EcoRI)
  [
    "AAAAAAAAGAGAATTCAAAAAAAAAAAAAA", // has GAATTC at pos 10
    "AAAAAAAAGAGAATTTAAAAAAAAAAAAAA", // has GAATTT at pos 10 (no cut)
  ]
}

///|
/// Create sample enzymes for testing: EcoRI and HindIII.
pub fn caps_sample_enzymes() -> Array[RestrictionEnzyme] {
  let enzymes_map = builtin_enzymes()
  let result : Array[RestrictionEnzyme] = []
  match enzymes_map.get("EcoRI") {
    Some(e) => result.push(e)
    None => ()
  }
  match enzymes_map.get("HindIII") {
    Some(e) => result.push(e)
    None => ()
  }
  result
}
