///|
/// Bio.SeqLocation module - Sequence location types for representing
/// feature annotations on sequences.
/// Inspired by Biopython Bio.SeqLocation.
///
/// Provides position types (ExactPosition, BeforePosition, AfterPosition,
/// OneOfPosition, WithinPosition) and location types (SimpleLocation,
/// CompoundLocation) used to represent the location of sequence features.
///
/// References:
///   - https://biopython.org/wiki/SeqLocation

// ============================================================================
// Position Types
// ============================================================================

///|
/// Position representation enum.
/// Each variant represents a different type of sequence position.
pub enum Pos {
  ExactPosition(Int)
  BeforePosition(Int)
  AfterPosition(Int)
  OneOfPosition(Int, Array[Pos])
  WithinPosition(Int, Pos, Pos)
}

///|
/// Get the numeric value of the position.
pub fn pos_get(pos : Pos) -> Int {
  match pos {
    ExactPosition(p) => p
    BeforePosition(p) => p
    AfterPosition(p) => p
    OneOfPosition(p, _) => p
    WithinPosition(p, _, _) => p
  }
}

///|
/// Get the string representation of a position.
pub fn pos_to_string(pos : Pos) -> String {
  match pos {
    ExactPosition(p) => p.to_string()
    BeforePosition(p) => "<" + p.to_string()
    AfterPosition(p) => ">" + p.to_string()
    OneOfPosition(_, opts) => {
      let mut s = "{"
      let mut i = 0
      while i < opts.length() {
        if i > 0 { s = s + "," }
        s = s + pos_to_string(opts[i])
        i = i + 1
      }
      s = s + "}"
      s
    }
    WithinPosition(_, left, right) => {
      let mut s = pos_to_string(left) + ".." + pos_to_string(right)
      s = "(" + s + ")"
      s
    }
  }
}

///|
/// Check if a position is exact.
pub fn pos_is_exact(pos : Pos) -> Bool {
  match pos {
    ExactPosition(_) => true
    _ => false
  }
}

///|
/// Create an exact position.
pub fn exact_position(pos : Int) -> Pos {
  ExactPosition(pos)
}

///|
/// Create a before position.
pub fn before_position(pos : Int) -> Pos {
  BeforePosition(pos)
}

///|
/// Create an after position.
pub fn after_position(pos : Int) -> Pos {
  AfterPosition(pos)
}

///|
/// Create a one-of position.
pub fn one_of_position(pos : Int, options : Array[Pos]) -> Pos {
  OneOfPosition(pos, options)
}

///|
/// Create a within position.
pub fn within_position(pos : Int, left : Pos, right : Pos) -> Pos {
  WithinPosition(pos, left, right)
}

// ============================================================================
// Location Types
// ============================================================================

///|
/// Location representation enum.
pub enum Loc {
  SimpleLocation(Pos, Pos, String, String)
  CompoundLocation(Array[Loc], String, String)
}

///|
/// Create a simple location.
pub fn simple_location(
  start_pos : Pos,
  end_pos : Pos,
  strand? : String = "+",
  ref_id? : String = "",
) -> Loc {
  SimpleLocation(start_pos, end_pos, strand, ref_id)
}

///|
/// Create a compound location.
pub fn compound_location(
  locs : Array[Loc],
  strand? : String = "+",
  ref_id? : String = "",
) -> Loc {
  CompoundLocation(locs, strand, ref_id)
}

///|
/// Create a compound location from arrays of start/end positions.
pub fn compound_from_simple(
  starts : Array[Int],
  ends : Array[Int],
  strand? : String = "+",
) -> Loc {
  let locs : Array[Loc] = Array::new()
  let mut i = 0
  while i < starts.length() {
    let sp = exact_position(starts[i])
    let ep = exact_position(ends[i])
    locs.push(simple_location(sp, ep, strand=strand))
    i = i + 1
  }
  CompoundLocation(locs, strand, "")
}

///|
/// Get the start position (0-based).
pub fn loc_start(loc : Loc) -> Int {
  match loc {
    SimpleLocation(sp, _, _, _) => pos_get(sp)
    CompoundLocation(locs, _, _) => {
      if locs.length() == 0 { return 0 }
      locs[0].start()
    }
  }
}

///|
/// Get the end position (0-based, exclusive).
pub fn loc_end(loc : Loc) -> Int {
  match loc {
    SimpleLocation(_, ep, _, _) => pos_get(ep)
    CompoundLocation(locs, _, _) => {
      if locs.length() == 0 { return 0 }
      let last = locs[locs.length() - 1]
      last.end()
    }
  }
}

///|
/// Get the strand.
pub fn loc_strand(loc : Loc) -> String {
  match loc {
    SimpleLocation(_, _, strand, _) => strand
    CompoundLocation(_, strand, _) => strand
  }
}

///|
/// Get the length of the location.
pub fn loc_len(loc : Loc) -> Int {
  match loc {
    SimpleLocation(sp, ep, _, _) => {
      let e = pos_get(ep)
      let s = pos_get(sp)
      if e > s { e - s } else { s - e }
    }
    CompoundLocation(locs, _, _) => {
      let mut total = 0
      let mut i = 0
      while i < locs.length() {
        total = total + locs[i].len()
        i = i + 1
      }
      total
    }
  }
}

///|
/// Check if this is a compound location.
pub fn loc_is_compound(loc : Loc) -> Bool {
  match loc {
    CompoundLocation(_, _, _) => true
    _ => false
  }
}

///|
/// Get the sub-locations.
pub fn loc_parts(loc : Loc) -> Array[Loc] {
  match loc {
    SimpleLocation(_, _, _, _) => {
      let result : Array[Loc] = Array::new()
      result.push(loc)
      result
    }
    CompoundLocation(locs, _, _) => locs.copy()
  }
}

///|
/// Check if two locations overlap.
pub fn loc_overlaps(loc : Loc, other : Loc) -> Bool {
  match loc {
    SimpleLocation(_, _, _, _) => {
      let s1 = loc.start()
      let e1 = loc.end()
      let s2 = other.start()
      let e2 = other.end()
      let (a1, b1) = if s1 <= e1 { (s1, e1) } else { (e1, s1) }
      let (a2, b2) = if s2 <= e2 { (s2, e2) } else { (e2, s2) }
      a1 < b2 && a2 < b1
    }
    CompoundLocation(locs, _, _) => {
      let mut i = 0
      while i < locs.length() {
        if locs[i].overlaps(other) {
          return true
        }
        i = i + 1
      }
      false
    }
  }
}

///|
/// Check if this location contains a position.
pub fn loc_contains(loc : Loc, pos : Int) -> Bool {
  match loc {
    SimpleLocation(sp, ep, _, _) => {
      let s = pos_get(sp)
      let e = pos_get(ep)
      if s <= e { pos >= s && pos < e } else { pos > e && pos <= s }
    }
    CompoundLocation(locs, _, _) => {
      let mut i = 0
      while i < locs.length() {
        if locs[i].contains(pos) {
          return true
        }
        i = i + 1
      }
      false
    }
  }
}

///|
/// Get location as string.
pub fn loc_to_string(loc : Loc) -> String {
  match loc {
    SimpleLocation(sp, ep, _, _) => pos_to_string(sp) + ".." + pos_to_string(ep)
    CompoundLocation(locs, _, _) => {
      let mut s = ""
      let mut i = 0
      while i < locs.length() {
        if i > 0 { s = s + "," }
        s = s + locs[i].to_string()
        i = i + 1
      }
      s
    }
  }
}

///|
/// Return location in 1-based coordinates (GFF/GenBank style).
pub fn loc_to_one_based(loc : Loc) -> String {
  match loc {
    SimpleLocation(sp, ep, _, _) => {
      let s = pos_get(sp) + 1
      let e = pos_get(ep)
      "(" + s.to_string() + ", " + e.to_string() + ")"
    }
    CompoundLocation(locs, _, _) => {
      let mut s = "join("
      let mut i = 0
      while i < locs.length() {
        if i > 0 { s = s + "," }
        s = s + locs[i].to_one_based()
        i = i + 1
      }
      s = s + ")"
      s
    }
  }
}

// ============================================================================
// Pos methods
// ============================================================================

/// Get the numeric value of the position.
pub fn Pos::get(self : Pos) -> Int { pos_get(self) }
/// Get the string representation of a position.
pub fn Pos::to_string(self : Pos) -> String { pos_to_string(self) }
/// Check if a position is exact.
pub fn Pos::is_exact(self : Pos) -> Bool { pos_is_exact(self) }

// ============================================================================
// Loc methods
// ============================================================================

/// Get the start position (0-based).
pub fn Loc::start(self : Loc) -> Int { loc_start(self) }
/// Get the end position (0-based, exclusive).
pub fn Loc::end(self : Loc) -> Int { loc_end(self) }
/// Get the strand.
pub fn Loc::strand(self : Loc) -> String { loc_strand(self) }
/// Get the length of the location.
pub fn Loc::len(self : Loc) -> Int { loc_len(self) }
/// Check if this is a compound location.
pub fn Loc::is_compound(self : Loc) -> Bool { loc_is_compound(self) }
/// Get the sub-locations.
pub fn Loc::parts(self : Loc) -> Array[Loc] { loc_parts(self) }
/// Check if two locations overlap.
pub fn Loc::overlaps(self : Loc, other : Loc) -> Bool { loc_overlaps(self, other) }
/// Check if this location contains a position.
pub fn Loc::contains(self : Loc, pos : Int) -> Bool { loc_contains(self, pos) }
/// Get location as string.
pub fn Loc::to_string(self : Loc) -> String { loc_to_string(self) }
/// Return location in 1-based coordinates.
pub fn Loc::to_one_based(self : Loc) -> String { loc_to_one_based(self) }

// ============================================================================
// Helper Functions
// ============================================================================

///|
/// Create a location string in GenBank format.
pub fn location_to_genbank(loc : Loc) -> String {
  loc.to_one_based()
}

///|
/// Create a simple location from 1-based GenBank coordinates.
pub fn from_genbank_coords(start : Int, end_ : Int, strand? : String = "+") -> Loc {
  let sp = exact_position(start - 1)
  let ep = exact_position(end_)
  simple_location(sp, ep, strand=strand)
}

///|
/// Parse a GenBank-style location string (simplified).
pub fn parse_genbank_location(loc_str : String) -> Loc {
  let s = loc_str
  if s.contains("..") {
    let parts = s.split("..").to_array()
    if parts.length() >= 2 {
      let start = parse_int_trimmed(parts[0]) - 1
      let end_ = parse_int_trimmed(parts[1])
      simple_location(exact_position(start), exact_position(end_))
    } else {
      simple_location(exact_position(0), exact_position(0))
    }
  } else {
    simple_location(exact_position(0), exact_position(0))
  }
}

/// Parse an integer from a trimmed string view, returning 1 on failure.
fn parse_int_trimmed(sv : StringView) -> Int {
  let s = sv.trim().to_owned()
  let n = s.length()
  if n == 0 { return 1 }
  let mut result = 0
  let mut i = 0
  let zero = '0'.to_int()
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c < zero || c > zero + 9 { break }
    result = result * 10 + (c - zero)
    i = i + 1
  }
  if result == 0 && n > 0 { 1 } else { result }
}

///|
/// Create a sample location for testing.
pub fn seq_location_sample() -> Array[Loc] {
  let locs : Array[Loc] = Array::new()
  locs.push(simple_location(exact_position(10), exact_position(50), strand="+"))
  locs.push(simple_location(exact_position(100), exact_position(200), strand="-"))
  let exon1 = simple_location(exact_position(0), exact_position(100), strand="+")
  let exon2 = simple_location(exact_position(200), exact_position(350), strand="+")
  let exon3 = simple_location(exact_position(500), exact_position(600), strand="+")
  let exon_locs : Array[Loc] = Array::new()
  exon_locs.push(exon1)
  exon_locs.push(exon2)
  exon_locs.push(exon3)
  locs.push(compound_location(exon_locs, strand="+"))
  locs.push(simple_location(before_position(50), after_position(100), strand="+"))
  let opts : Array[Pos] = Array::new()
  opts.push(exact_position(5))
  opts.push(exact_position(7))
  opts.push(exact_position(9))
  locs.push(simple_location(one_of_position(7, opts), exact_position(50), strand="+"))
  locs
}
