///|
/// Secondary structure content analysis — fraction of helix, sheet, coil.
///
/// Biotite provides secondary structure fraction analysis. This module
/// provides the same for BioSeqs' DSSP output.

///|
/// Secondary structure content summary.
pub struct SSContent {
  helix_fraction : Double
  sheet_fraction : Double
  coil_fraction : Double
  turn_fraction : Double
  helix_count : Int32
  sheet_count : Int32
  coil_count : Int32
  turn_count : Int32
  total_residues : Int32
}

///|
/// Get helix fraction.
pub fn SSContent::get_helix_fraction(self : SSContent) -> Double {
  self.helix_fraction
}

///|
/// Get sheet fraction.
pub fn SSContent::get_sheet_fraction(self : SSContent) -> Double {
  self.sheet_fraction
}

///|
/// Get coil fraction.
pub fn SSContent::get_coil_fraction(self : SSContent) -> Double {
  self.coil_fraction
}

///|
/// DSSP single-letter secondary structure codes:
/// H = alpha helix, G = 3-10 helix, I = pi helix
/// E = beta strand, B = beta bridge
/// T = turn, S = bend, C = coil (or space)
fn is_helix(ss : Char) -> Bool {
  ss == 'H' || ss == 'G' || ss == 'I'
}

///|
fn is_sheet(ss : Char) -> Bool {
  ss == 'E' || ss == 'B'
}

///|
fn is_turn(ss : Char) -> Bool {
  ss == 'T' || ss == 'S'
}

///|
fn is_coil(ss : Char) -> Bool {
  ss == 'C' || ss == ' ' || (!is_helix(ss) && !is_sheet(ss) && !is_turn(ss))
}

///|
/// Calculate secondary structure content from DSSP assignments.
///
/// Takes an array of single-character DSSP codes and returns
/// the fraction of helix, sheet, turn, and coil residues.
pub fn calculate_ss_content(ss_codes : Array[Char]) -> SSContent {
  let total = ss_codes.length()
  let mut h = 0
  let mut e = 0
  let mut t = 0
  let mut c = 0
  for ss in ss_codes {
    if is_helix(ss) {
      h = h + 1
    } else if is_sheet(ss) {
      e = e + 1
    } else if is_turn(ss) {
      t = t + 1
    } else {
      c = c + 1
    }
  }
  let tot = if total > 0 { total.to_double() } else { 1.0 }
  SSContent::{
    helix_fraction: h.to_double() / tot,
    sheet_fraction: e.to_double() / tot,
    coil_fraction: c.to_double() / tot,
    turn_fraction: t.to_double() / tot,
    helix_count: h,
    sheet_count: e,
    coil_count: c,
    turn_count: t,
    total_residues: total,
  }
}

///|
/// Count number of helical segments (contiguous helix runs).
pub fn count_helix_segments(ss_codes : Array[Char]) -> Int32 {
  let mut count = 0
  let mut in_helix = false
  for ss in ss_codes {
    if is_helix(ss) {
      if !in_helix {
        count = count + 1
        in_helix = true
      }
    } else {
      in_helix = false
    }
  }
  count
}

///|
/// Count number of beta strand segments.
pub fn count_strand_segments(ss_codes : Array[Char]) -> Int32 {
  let mut count = 0
  let mut in_strand = false
  for ss in ss_codes {
    if is_sheet(ss) {
      if !in_strand {
        count = count + 1
        in_strand = true
      }
    } else {
      in_strand = false
    }
  }
  count
}

///|
/// Get helix segment boundaries (start, end) pairs.
pub fn helix_boundaries(ss_codes : Array[Char]) -> Array[(Int32, Int32)] {
  let result : Array[(Int32, Int32)] = []
  let mut start = -1
  for i in 0..<ss_codes.length() {
    if is_helix(ss_codes[i]) {
      if start < 0 {
        start = i
      }
    } else if start >= 0 {
      result.push((start, i - 1))
      start = -1
    }
  }
  if start >= 0 {
    result.push((start, ss_codes.length() - 1))
  }
  result
}

///|
/// Get strand segment boundaries (start, end) pairs.
pub fn strand_boundaries(ss_codes : Array[Char]) -> Array[(Int32, Int32)] {
  let result : Array[(Int32, Int32)] = []
  let mut start = -1
  for i in 0..<ss_codes.length() {
    if is_sheet(ss_codes[i]) {
      if start < 0 {
        start = i
      }
    } else if start >= 0 {
      result.push((start, i - 1))
      start = -1
    }
  }
  if start >= 0 {
    result.push((start, ss_codes.length() - 1))
  }
  result
}
