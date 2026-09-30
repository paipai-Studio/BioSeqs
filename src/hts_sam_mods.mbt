///|
/// htslib base modification parsing (sam_mods.c port).
/// Parses MM/ML tags from BAM records to extract modified base information
/// (5mC, 6mA, etc.) used in long-read methylation analysis.
///
/// MM tag format: "C+m,5,10;C+h,3,8;" means
///   - Canonical base C, + strand, modification code 'm', at 0-based delta positions 5,10
///   - Canonical base C, + strand, modification code 'h', at delta positions 3,8
/// Deltas are relative to the count of the canonical base in the sequence.
///
/// ML tag: array of quality scores (0-255), one per modification site.

///|
/// A single base modification call.
pub struct BaseMod {
  /// Position in the query sequence (0-based)
  pos : Int
  /// Canonical base (A/C/G/T/N)
  canonical : Char
  /// Modification code (e.g. 'm' for 5mC, 'h' for 5hmC, 'a' for 6mA)
  code : Char
  /// Strand: +1 for forward, -1 for reverse
  strand : Int
  /// Quality score (0-255), -1 if unavailable
  qual : Int
}

///|
/// Parsed modification type from MM tag.
pub struct ModType {
  /// Canonical base
  canonical : Char
  /// Modification code
  code : Char
  /// Strand: +1 or -1
  strand : Int
  /// Whether this is implicit (unlisted positions = unmodified)
  implicit : Bool
  /// Positions (0-based query positions) where this modification occurs
  positions : Array[Int]
  /// Quality scores for each position
  quals : Array[Int]
}

///|
/// Base modification state, parsed from MM/ML tags.
pub struct BaseModState {
  /// All modification types found in the MM tag
  types : Array[ModType]
}

///|
/// Parse the MM tag string to extract modification types and delta positions.
/// Returns an array of (canonical, code, strand, implicit, deltas) tuples.
fn parse_mm_tag(mm : String) -> Array[(Char, Char, Int, Bool, Array[Int])] {
  let result : Array[(Char, Char, Int, Bool, Array[Int])] = []
  let chars = mm.to_array()
  let n = chars.length()
  let mut i = 0
  while i < n {
    // Skip whitespace
    if chars[i] == ' ' || chars[i] == ';' {
      i = i + 1
      continue
    }
    // Parse canonical base
    let canonical = chars[i]
    if canonical != 'A' && canonical != 'C' && canonical != 'G' && canonical != 'T' && canonical != 'U' && canonical != 'N' {
      i = i + 1
      continue
    }
    let canonical = if canonical == 'U' { 'T' } else { canonical }
    i = i + 1
    if i >= n {
      break
    }
    // Parse strand
    let strand = if chars[i] == '+' {
      i = i + 1
      1
    } else if chars[i] == '-' {
      i = i + 1
      -1
    } else {
      1
    }
    // Check for explicit/implicit marker
    let mut implicit = true
    if i < n && (chars[i] == '.' || chars[i] == '?') {
      implicit = chars[i] == '.'
      i = i + 1
    }
    // Parse modification codes (one or more letters, or ChEBI number)
    // For simplicity, parse one code at a time
    while i < n && chars[i] != ',' && chars[i] != ';' {
      let code = chars[i]
      i = i + 1
      // Skip to delta list
      // Now parse delta list: ,d1,d2,d3...;
      let deltas : Array[Int] = []
      while i < n && chars[i] == ',' {
        i = i + 1
        let mut num = 0
        let mut has_digit = false
        while i < n && chars[i] >= '0' && chars[i] <= '9' {
          num = num * 10 + (chars[i].to_int() - '0'.to_int())
          has_digit = true
          i = i + 1
        }
        if has_digit {
          deltas.push(num)
        }
      }
      // Skip to ; or next code
      while i < n && chars[i] != ';' && chars[i] != ',' && chars[i] != 'A' && chars[i] != 'C' && chars[i] != 'G' && chars[i] != 'T' && chars[i] != 'N' {
        i = i + 1
      }
      result.push((canonical, code, strand, implicit, deltas))
      // If next char is a letter (another mod code), continue inner loop
      if i < n && chars[i] >= 'a' && chars[i] <= 'z' {
        continue
      }
      break
    }
    // Skip to next chain (after ;)
    while i < n && chars[i] != ';' {
      i = i + 1
    }
    if i < n {
      i = i + 1 // skip ;
    }
  }
  result
}

///|
/// Convert delta positions to absolute query positions.
/// Deltas are 0-based offsets into the positions of the canonical base.
/// E.g. for sequence "ACGTCGT" and canonical 'C', C appears at positions 1,5.
/// Delta 0 → position 1, delta 1 → position 5.
fn deltas_to_positions(
  deltas : Array[Int],
  canonical : Char,
  sequence : String,
) -> Array[Int] {
  // Find all positions of canonical base in sequence
  let base_positions : Array[Int] = []
  for i = 0; i < sequence.length(); i = i + 1 {
    if sequence.unsafe_get(i).to_int() == canonical.to_int() {
      base_positions.push(i)
    }
  }
  // Convert deltas to absolute positions
  let result : Array[Int] = []
  let mut prev_delta = -1
  for d in deltas {
    let curr = if prev_delta == -1 { d } else { prev_delta + d + 1 }
    prev_delta = curr
    if curr >= 0 && curr < base_positions.length() {
      result.push(base_positions[curr])
    }
  }
  result
}

///|
/// Parse base modifications from MM and ML tags.
/// MM tag: string like "C+m,5,10;C+h,3,8;"
/// ML tag: integer array of quality scores
/// sequence: the query sequence (for converting deltas to positions)
pub fn parse_base_mods(
  mm_tag : String,
  ml_tag : Array[Int]?,
  sequence : String,
) -> BaseModState {
  let parsed = parse_mm_tag(mm_tag)
  let types : Array[ModType] = []
  let mut ml_idx = 0
  for (canonical, code, strand, implicit, deltas) in parsed {
    let positions = deltas_to_positions(deltas, canonical, sequence)
    let quals : Array[Int] = []
    for _ in positions {
      match ml_tag {
        Some(arr) => {
          if ml_idx < arr.length() {
            quals.push(arr[ml_idx])
            ml_idx = ml_idx + 1
          } else {
            quals.push(-1)
          }
        }
        None => quals.push(-1)
      }
    }
    types.push(ModType::{ canonical, code, strand, implicit, positions, quals })
  }
  BaseModState::{ types }
}

///|
/// Get all base modification calls from the state.
pub fn BaseModState::all_mods(self : BaseModState) -> Array[BaseMod] {
  let result : Array[BaseMod] = []
  for t in self.types {
    for i = 0; i < t.positions.length(); i = i + 1 {
      result.push(BaseMod::{
        pos: t.positions[i],
        canonical: t.canonical,
        code: t.code,
        strand: t.strand,
        qual: if i < t.quals.length() { t.quals[i] } else { -1 },
      })
    }
  }
  result
}

///|
/// Get modifications at a specific query position.
pub fn BaseModState::mods_at(
  self : BaseModState,
  pos : Int,
) -> Array[BaseMod] {
  let result : Array[BaseMod] = []
  for t in self.types {
    for i = 0; i < t.positions.length(); i = i + 1 {
      if t.positions[i] == pos {
        result.push(BaseMod::{
          pos,
          canonical: t.canonical,
          code: t.code,
          strand: t.strand,
          qual: if i < t.quals.length() { t.quals[i] } else { -1 },
        })
      }
    }
  }
  result
}

///|
/// Get all modification types recorded (codes).
pub fn BaseModState::mod_codes(self : BaseModState) -> Array[Char] {
  let result : Array[Char] = []
  for t in self.types {
    let mut found = false
    for c in result {
      if c == t.code {
        found = true
        break
      }
    }
    if not(found) {
      result.push(t.code)
    }
  }
  result
}

///|
/// Get modification count by code.
pub fn BaseModState::count_by_code(
  self : BaseModState,
  code : Char,
) -> Int {
  let mut total = 0
  for t in self.types {
    if t.code == code {
      total = total + t.positions.length()
    }
  }
  total
}

///|
/// Methylation summary: counts of modified vs unmodified for each canonical base.
pub struct MethylationSummary {
  /// Canonical base
  canonical : Char
  /// Modification code
  code : Char
  /// Number of modified bases
  modified : Int
  /// Number of unmodified bases (total canonical - modified)
  unmodified : Int
  /// Methylation fraction (modified / total)
  fraction : Double
}

///|
/// Compute methylation summary from base mod state and sequence.
pub fn BaseModState::methylation_summary(
  self : BaseModState,
  sequence : String,
) -> Array[MethylationSummary] {
  let result : Array[MethylationSummary] = []
  for t in self.types {
    // Count total occurrences of canonical base
    let mut total = 0
    for i = 0; i < sequence.length(); i = i + 1 {
      if sequence.unsafe_get(i).to_int() == t.canonical.to_int() {
        total = total + 1
      }
    }
    let modified = t.positions.length()
    let unmodified = total - modified
    let fraction = if total > 0 {
      modified.to_double() / total.to_double()
    } else {
      0.0
    }
    result.push(MethylationSummary::{
      canonical: t.canonical,
      code: t.code,
      modified,
      unmodified,
      fraction,
    })
  }
  result
}
