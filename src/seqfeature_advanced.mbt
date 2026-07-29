///|
/// Bio.SeqFeature advanced features: CompoundLocation extensions,
/// LocationParser, SeqFeatureExtended, and location utilities.
///
/// This module extends the base seqfeature module with:
///   - CompoundLocation helper constructors and methods
///   - GenBank/EMBL location string parser
///   - Extended sequence features with qualifiers
///   - Location manipulation utilities

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

///|
/// Trim leading and trailing whitespace (spaces, tabs, newlines) from a string.
fn trim_whitespace(s : String) -> String {
  let n = s.length()
  let mut start = 0
  while start < n {
    let c = s.unsafe_get(start)
    if c == ' '.to_int().to_uint16() ||
      c == '\t'.to_int().to_uint16() ||
      c == '\n'.to_int().to_uint16() ||
      c == '\r'.to_int().to_uint16() {
      start = start + 1
    } else {
      break
    }
  }
  let mut end = n
  while end > start {
    let c = s.unsafe_get(end - 1)
    if c == ' '.to_int().to_uint16() ||
      c == '\t'.to_int().to_uint16() ||
      c == '\n'.to_int().to_uint16() ||
      c == '\r'.to_int().to_uint16() {
      end = end - 1
    } else {
      break
    }
  }
  if start >= end {
    ""
  } else {
    s[start:end].to_owned()
  }
}

///|
/// Find the index of the matching closing parenthesis starting from `start`.
/// The string at `start` should be an opening '('.
/// Returns None if no matching paren is found.
fn find_matching_paren(s : String, start : Int) -> Int? {
  let n = s.length()
  let mut depth = 0
  let mut i = start
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c == '(' {
      depth = depth + 1
    } else if c == ')' {
      depth = depth - 1
      if depth == 0 {
        return Some(i)
      }
    }
    i = i + 1
  }
  None
}

///|
/// Parse a non-negative integer from a string starting at `pos`.
/// Returns the parsed integer and the position after the last digit.
fn parse_int_at(s : String, pos : Int) -> (Int, Int) raise LocationParserError {
  let n = s.length()
  let mut i = pos
  let mut result = 0
  let mut found = false
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c >= '0'.to_int() && c <= '9'.to_int() {
      result = result * 10 + (c - '0'.to_int())
      found = true
      i = i + 1
    } else {
      break
    }
  }
  if !found {
    raise LocationParserError(
      "Expected integer at position " +
      pos.to_string() +
      " in '" +
      s +
      "'",
    )
  }
  (result, i)
}

///|
/// Split a string by commas that are at the top parenthesis nesting level.
/// Returns an array of trimmed substring contents.
fn split_top_level_commas(s : String) -> Array[String] {
  let parts : Array[String] = Array::new()
  let n = s.length()
  let mut depth = 0
  let mut start = 0
  let mut i = 0
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c == '(' {
      depth = depth + 1
    } else if c == ')' {
      depth = depth - 1
    } else if c == ',' && depth == 0 {
      let part = trim_whitespace(s[start:i].to_owned())
      parts.push(part)
      start = i + 1
    }
    i = i + 1
  }
  let last = trim_whitespace(s[start:n].to_owned())
  if last.length() > 0 {
    parts.push(last)
  }
  parts
}

// ---------------------------------------------------------------------------
// ParsedLocation - union type for parse results
// ---------------------------------------------------------------------------

///|
/// Result of parsing a location string: either a simple FeatureLocation
/// or a CompoundLocation.
pub enum ParsedLocation {
  ParsedFeatureLocation(FeatureLocation)
  ParsedCompoundLocation(CompoundLocation)
}

///|
/// Extract the FeatureLocation from a ParsedLocation, if it is one.
pub fn ParsedLocation::to_feature_location(
  self : ParsedLocation,
) -> FeatureLocation? {
  match self {
    ParsedFeatureLocation(loc) => Some(loc)
    ParsedCompoundLocation(_) => None
  }
}

///|
/// Extract the CompoundLocation from a ParsedLocation, if it is one.
pub fn ParsedLocation::to_compound_location(
  self : ParsedLocation,
) -> CompoundLocation? {
  match self {
    ParsedFeatureLocation(_) => None
    ParsedCompoundLocation(comp) => Some(comp)
  }
}

///|
/// Check if this ParsedLocation is a FeatureLocation.
pub fn ParsedLocation::is_feature_location(self : ParsedLocation) -> Bool {
  match self {
    ParsedFeatureLocation(_) => true
    ParsedCompoundLocation(_) => false
  }
}

///|
/// Check if this ParsedLocation is a CompoundLocation.
pub fn ParsedLocation::is_compound_location(self : ParsedLocation) -> Bool {
  match self {
    ParsedFeatureLocation(_) => false
    ParsedCompoundLocation(_) => true
  }
}

// ---------------------------------------------------------------------------
// 1. CompoundLocation extensions
// ---------------------------------------------------------------------------

///|
/// Create a CompoundLocation with "join" operator.
/// `join` merges locations that are on the same strand and should be
/// combined (e.g., exons of a gene).
pub fn CompoundLocation::join(
  locations : Array[FeatureLocation],
) -> CompoundLocation {
  CompoundLocation::{
    parts: locations,
    operator: "join",
    ref_seq: None,
    ref_db: None,
  }
}

///|
/// Create a CompoundLocation with "order" operator.
/// `order` indicates locations that should be taken in order but
/// are not necessarily contiguous on the same strand.
pub fn CompoundLocation::order(
  locations : Array[FeatureLocation],
) -> CompoundLocation {
  CompoundLocation::{
    parts: locations,
    operator: "order",
    ref_seq: None,
    ref_db: None,
  }
}

///|
/// Get the start position of the first part in this compound location.
pub fn compound_start(self : CompoundLocation) -> Position {
  if self.parts.length() == 0 {
    return UnknownPosition
  }
  self.parts[0].start
}

///|
/// Get the end position of the last part in this compound location.
pub fn compound_end(self : CompoundLocation) -> Position {
  if self.parts.length() == 0 {
    return UnknownPosition
  }
  self.parts[self.parts.length() - 1].end
}

///|
/// Extract the compound feature sequence from a parent sequence.
/// Each part's sequence is extracted and concatenated in order.
/// If any part has strand -1, that part is reverse-complemented.
pub fn compound_extract(
  self : CompoundLocation,
  parent_sequence : Seq,
) -> String raise LocationParserError {
  if self.parts.length() == 0 {
    return ""
  }
  let buf = StringBuilder::new()
  for part in self.parts {
    let extracted = part.extract(parent_sequence)
    buf.write_string(extracted.to_string())
  }
  buf.to_string()
}

///|
/// Generate a GenBank-format string representation of this compound location.
/// Example: "join(1..100,complement(200..300))"
pub fn compound_to_string(self : CompoundLocation) -> String {
  let buf = StringBuilder::new()
  let strand = self._compound_strand()
  if strand == Some(-1) {
    buf.write_string("complement(")
  }
  buf.write_string(self.operator)
  buf.write_char('(')
  for i, part in self.parts {
    if i > 0 {
      buf.write_string(",")
    }
    buf.write_string(_location_gb_string(part))
  }
  buf.write_char(')')
  if strand == Some(-1) {
    buf.write_char(')')
  }
  buf.to_string()
}

///|
/// Determine the compound strand by checking all parts.
/// If all parts share the same strand, return that strand.
/// Otherwise return None.
fn _compound_strand(self : CompoundLocation) -> Int? {
  if self.parts.length() == 0 {
    return None
  }
  let first_strand = self.parts[0].strand
  if first_strand == None {
    return None
  }
  let expected = first_strand.unwrap()
  for i in 1..<self.parts.length() {
    if self.parts[i].strand != first_strand {
      return None
    }
  }
  Some(expected)
}

///|
/// Convert a FeatureLocation to its GenBank string representation.
/// Handles strand wrapping for complement locations.
fn _location_gb_string(loc : FeatureLocation) -> String {
  let inner = _location_inner_gb_string(loc)
  if loc.strand == Some(-1) {
    "complement(" + inner + ")"
  } else {
    inner
  }
}

///|
/// Inner GenBank string for a FeatureLocation (without complement wrapping).
fn _location_inner_gb_string(loc : FeatureLocation) -> String {
  _position_gb_string(loc.start) + ".." + _position_gb_string(loc.end)
}

///|
/// Convert a Position to its GenBank string representation.
fn _position_gb_string(pos : Position) -> String {
  match pos {
    ExactPosition(p) => p.to_string()
    BeforePosition(p) => "<" + p.to_string()
    AfterPosition(p) => ">" + p.to_string()
    WithinPosition(p, fs, fe) =>
      "<" +
      fs.to_string() +
      ".." +
      p.to_string() +
      ".." +
      fe.to_string() +
      ">"
    UncertainPosition(p) => "?" + p.to_string() + "?"
    UnknownPosition => "?"
  }
}

// ---------------------------------------------------------------------------
// 2. LocationParser
// ---------------------------------------------------------------------------

///|
/// Parse a position string into a Position enum value.
///
/// Supported formats:
///   - "100"       → ExactPosition(100)
///   - "<100"      → BeforePosition(100)
///   - ">100"      → AfterPosition(100)
///   - "?100?"     → UncertainPosition(100)
///   - "100^200"   → WithinPosition(100, 100, 200)
///   - "?"         → UnknownPosition
pub fn parse_position_string(
  pos_string : String,
) -> Position raise LocationParserError {
  let trimmed = trim_whitespace(pos_string)
  if trimmed.length() == 0 {
    raise LocationParserError("Empty position string")
  }
  if trimmed == "?" {
    return UnknownPosition
  }
  let first = trimmed.unsafe_get(0).to_int()
  let last = trimmed.unsafe_get(trimmed.length() - 1).to_int()

  if first == '<'.to_int() {
    let num_str = trim_whitespace(trimmed[1:trimmed.length()].to_owned())
    let (val, _) = parse_int_at(num_str, 0)
    return BeforePosition(val)
  }
  if first == '>'.to_int() {
    let num_str = trim_whitespace(trimmed[1:trimmed.length()].to_owned())
    let (val, _) = parse_int_at(num_str, 0)
    return AfterPosition(val)
  }
  if first == '?'.to_int() && last == '?'.to_int() && trimmed.length() > 2 {
    let num_str = trim_whitespace(trimmed[1:trimmed.length() - 1].to_owned())
    let (val, _) = parse_int_at(num_str, 0)
    return UncertainPosition(val)
  }

  let caret_pos = trimmed.find("^")
  if caret_pos != None {
    let cp = caret_pos.unwrap()
    let left = trim_whitespace(trimmed[0:cp].to_owned())
    let right = trim_whitespace(trimmed[cp + 1:trimmed.length()].to_owned())
    let (left_val, _) = parse_int_at(left, 0)
    let (right_val, _) = parse_int_at(right, 0)
    return WithinPosition(left_val, left_val, right_val)
  }

  let (val, _) = parse_int_at(trimmed, 0)
  ExactPosition(val)
}

///|
/// Parse a simple location string like "1..100" or "<1..>100".
///
/// Supported formats:
///   - "1..100"     → FeatureLocation(ExactPosition(1), ExactPosition(100))
///   - "<1..100"    → FeatureLocation(BeforePosition(1), ExactPosition(100))
///   - "1..>100"    → FeatureLocation(ExactPosition(1), AfterPosition(100))
///   - "<1..>100"   → FeatureLocation(BeforePosition(1), AfterPosition(100))
///   - "100^200..300" → Fuzzy start position supported
pub fn parse_simple_location(
  loc_string : String,
) -> FeatureLocation raise LocationParserError {
  let trimmed = trim_whitespace(loc_string)
  let dots = ".."
  let dot_pos = match trimmed.find(dots) {
    Some(p) => p
    None =>
      raise LocationParserError(
        "Expected '..' in location string: " + trimmed,
      )
  }
  let start_str = trim_whitespace(trimmed[0:dot_pos].to_owned())
  let end_str = trim_whitespace(
    trimmed[dot_pos + 2:trimmed.length()].to_owned(),
  )
  let start_pos = parse_position_string(start_str)
  let end_pos = parse_position_string(end_str)
  FeatureLocation::new_safe(start_pos, end_pos)
}

///|
/// Parse a complement location string like "complement(1..100)".
/// Returns the parsed ParsedLocation and its strand.
///
/// Supported formats:
///   - "complement(1..100)"          → (ParsedFeatureLocation, -1)
///   - "complement(join(1..100,200..300))" → (ParsedCompoundLocation, -1)
pub fn parse_complement_location(
  loc_string : String,
) -> (ParsedLocation, Int) raise LocationParserError {
  let trimmed = trim_whitespace(loc_string)
  let prefix = "complement("
  if !trimmed.starts_with(prefix) {
    raise LocationParserError(
      "Expected 'complement(' prefix in: " + trimmed,
    )
  }
  let paren_start = prefix.length() - 1
  let after_prefix = trim_whitespace(
    trimmed[paren_start:trimmed.length()].to_owned(),
  )
  let close_paren = match find_matching_paren(after_prefix, 0) {
    Some(p) => p
    None =>
      raise LocationParserError(
        "Missing closing parenthesis in complement: " + trimmed,
      )
  }
  let inner = trim_whitespace(after_prefix[1:close_paren].to_owned())
  let parsed = parse_location_string(inner)
  (parsed, -1)
}

///|
/// Parse a join() or order() compound location string.
///
/// Supported formats:
///   - "join(1..100,200..300)"
///   - "order(1..100,complement(200..300))"
pub fn parse_join_location(
  loc_string : String,
) -> CompoundLocation raise LocationParserError {
  let trimmed = trim_whitespace(loc_string)
  let (operator, after_op) = if trimmed.starts_with("join(") {
    ("join", trimmed[4:trimmed.length()].to_owned())
  } else if trimmed.starts_with("order(") {
    ("order", trimmed[5:trimmed.length()].to_owned())
  } else {
    raise LocationParserError(
      "Expected 'join(' or 'order(' in: " + trimmed,
    )
  }
  let close_paren = match find_matching_paren(after_op, 0) {
    Some(p) => p
    None =>
      raise LocationParserError(
        "Missing closing parenthesis in " +
        operator +
        ": " +
        trimmed,
      )
  }
  let inner = trim_whitespace(after_op[1:close_paren].to_owned())
  let part_strings = split_top_level_commas(inner)
  let parts : Array[FeatureLocation] = Array::new()
  for part_str in part_strings {
    if part_str.length() == 0 {
      continue
    }
    let parsed = parse_location_string(part_str)
    match parsed {
      ParsedFeatureLocation(loc) => {
        parts.push(loc)
      }
      ParsedCompoundLocation(_) =>
        raise LocationParserError(
          "Nested compound locations are not supported in " +
          operator +
          ": " +
          part_str,
        )
    }
  }
  CompoundLocation::{
    parts,
    operator,
    ref_seq: None,
    ref_db: None,
  }
}

///|
/// Parse a GenBank/EMBL feature location string into a FeatureLocation
/// or CompoundLocation.
///
/// Supported formats:
///   - Simple:     "1..100"                  → FeatureLocation
///   - Fuzzy:      "<1..100"                 → FeatureLocation with BeforePosition
///   - Complement: "complement(1..100)"      → FeatureLocation strand -1
///   - Compound:   "join(1..100,200..300)"  → CompoundLocation
///   - Order:      "order(1..100,200..300)" → CompoundLocation
///   - Nested:     "complement(join(1..100,200..300))" → CompoundLocation strand -1
pub fn parse_location_string(
  loc_string : String,
) -> ParsedLocation raise LocationParserError {
  let trimmed = trim_whitespace(loc_string)
  if trimmed.length() == 0 {
    raise LocationParserError("Empty location string")
  }

  if trimmed.starts_with("complement(") {
    let (parsed, strand) = parse_complement_location(trimmed)
    match parsed {
      ParsedFeatureLocation(loc) => {
        let new_loc = _set_strand(loc, strand)
        return ParsedFeatureLocation(new_loc)
      }
      ParsedCompoundLocation(comp) => {
        let new_parts : Array[FeatureLocation] = Array::new()
        for part in comp.parts {
          new_parts.push(_set_strand(part, strand))
        }
        return ParsedCompoundLocation(CompoundLocation::{
          parts: new_parts,
          operator: comp.operator,
          ref_seq: comp.ref_seq,
          ref_db: comp.ref_db,
        })
      }
    }
  }

  if trimmed.starts_with("join(") || trimmed.starts_with("order(") {
    let compound = parse_join_location(trimmed)
    return ParsedCompoundLocation(compound)
  }

  let loc = parse_simple_location(trimmed)
  ParsedFeatureLocation(loc)
}

///|
/// Set the strand on a FeatureLocation, creating a new copy.
fn _set_strand(loc : FeatureLocation, strand : Int) -> FeatureLocation {
  let new_strand : Int? = if strand == 1 {
    Some(1)
  } else if strand == -1 {
    Some(-1)
  } else {
    Some(0)
  }
  FeatureLocation::{
    start: loc.start,
    end: loc.end,
    strand: new_strand,
    ref_seq: loc.ref_seq,
    ref_db: loc.ref_db,
  }
}

// ---------------------------------------------------------------------------
// 3. SeqFeatureExtended
// ---------------------------------------------------------------------------

///|
/// Extended sequence feature with qualifiers (key-value pairs).
///
/// This provides a richer feature representation than SeqFeature,
/// with explicit qualifiers and a GenBank-compatible string representation.
pub struct SeqFeatureExtended {
  location : FeatureLocation
  type_ : String
  qualifiers : Map[String, Array[String]]
  id : String
  description : String
} derive(Eq)

///|
/// Create a new SeqFeatureExtended.
///
/// Parameters:
///   - location:    FeatureLocation for this feature
///   - type_:       Feature type (e.g., "gene", "CDS", "exon")
///   - id?:         Optional feature identifier
///   - description?: Optional feature description
pub fn SeqFeatureExtended::new(
  location : FeatureLocation,
  type_ : String,
  id? : String = "",
  description? : String = "",
) -> SeqFeatureExtended {
  SeqFeatureExtended::{
    location,
    type_,
    qualifiers: Map([], capacity=0),
    id,
    description,
  }
}

///|
/// Get the feature location.
pub fn SeqFeatureExtended::location(self : SeqFeatureExtended) -> FeatureLocation {
  self.location
}

///|
/// Get the feature type (e.g., "gene", "CDS", "exon").
pub fn SeqFeatureExtended::type_(self : SeqFeatureExtended) -> String {
  self.type_
}

///|
/// Get the feature id.
pub fn SeqFeatureExtended::id(self : SeqFeatureExtended) -> String {
  self.id
}

///|
/// Get the feature description.
pub fn SeqFeatureExtended::description(self : SeqFeatureExtended) -> String {
  self.description
}

///|
/// Get all qualifiers as a map.
pub fn SeqFeatureExtended::qualifiers(
  self : SeqFeatureExtended,
) -> Map[String, Array[String]] {
  self.qualifiers
}

///|
/// Add a qualifier key-value pair to this feature.
/// If the key already exists, the value is appended to the existing array.
pub fn add_qualifier(
  self : SeqFeatureExtended,
  key : String,
  value : String,
) -> SeqFeatureExtended {
  let new_quals : Map[String, Array[String]] = Map(
    [],
    capacity=self.qualifiers.length() + 1,
  )
  for k, v in self.qualifiers {
    new_quals[k] = v
  }
  match new_quals.get(key) {
    Some(arr) => {
      let new_arr = arr.copy()
      new_arr.push(value)
      new_quals[key] = new_arr
    }
    None => {
      new_quals[key] = [value]
    }
  }
  SeqFeatureExtended::{ ..self, qualifiers: new_quals }
}

///|
/// Get qualifier values for a given key.
/// Returns an empty array if the key does not exist.
pub fn get_qualifier(
  self : SeqFeatureExtended,
  key : String,
) -> Array[String] {
  match self.qualifiers.get(key) {
    Some(arr) => arr
    None => []
  }
}

///|
/// Generate a GenBank-format string representation of this feature.
///
/// Format:
///   type  location
///   /key="value"
///   /key="value"
pub fn seq_feature_to_string(self : SeqFeatureExtended) -> String {
  let buf = StringBuilder::new()
  buf.write_string(self.type_)
  buf.write_char(' ')
  buf.write_string(location_to_string(self.location))
  for key, values in self.qualifiers {
    for val in values {
      buf.write_string("\n/")
      buf.write_string(key)
      buf.write_string("=\"")
      buf.write_string(val)
      buf.write_char('"')
    }
  }
  buf.to_string()
}

// ---------------------------------------------------------------------------
// 4. Location Utilities
// ---------------------------------------------------------------------------

///|
/// Convert a FeatureLocation to its GenBank-format string representation.
///
/// Examples:
///   - FeatureLocation(ExactPosition(1), ExactPosition(100))
///     → "1..100"
///   - With strand -1
///     → "complement(1..100)"
///   - With fuzzy start
///     → "<1..100"
pub fn location_to_string(loc : FeatureLocation) -> String {
  _location_gb_string(loc)
}

///|
/// Create a new FeatureLocation with the strand complemented.
/// - strand +1 becomes -1
/// - strand -1 becomes +1
/// - strand 0 or None stays the same
pub fn feature_location_complement(
  loc : FeatureLocation,
) -> FeatureLocation {
  let new_strand = match loc.strand {
    Some(s) =>
      if s == 1 {
        Some(-1)
      } else if s == -1 {
        Some(1)
      } else {
        Some(0)
      }
    None => None
  }
  FeatureLocation::{
    start: loc.start,
    end: loc.end,
    strand: new_strand,
    ref_seq: loc.ref_seq,
    ref_db: loc.ref_db,
  }
}

///|
/// Shift all positions in a FeatureLocation by the given offset.
///
/// Positive offset moves positions forward, negative moves backward.
/// Fuzzy positions (BeforePosition, AfterPosition) are shifted by
/// modifying their internal value. WithinPosition is shifted
/// by updating both the main position and the fuzzy boundaries.
pub fn feature_location_shift(
  loc : FeatureLocation,
  offset : Int,
) -> FeatureLocation {
  let new_start = _shift_position(loc.start, offset)
  let new_end = _shift_position(loc.end, offset)
  FeatureLocation::{
    start: new_start,
    end: new_end,
    strand: loc.strand,
    ref_seq: loc.ref_seq,
    ref_db: loc.ref_db,
  }
}

///|
/// Shift a single Position by the given offset.
fn _shift_position(pos : Position, offset : Int) -> Position {
  match pos {
    ExactPosition(p) => ExactPosition(p + offset)
    WithinPosition(p, fs, fe) =>
      WithinPosition(p + offset, fs + offset, fe + offset)
    BeforePosition(p) => BeforePosition(p + offset)
    AfterPosition(p) => AfterPosition(p + offset)
    UncertainPosition(p) => UncertainPosition(p + offset)
    UnknownPosition => UnknownPosition
  }
}

///|
/// Merge two overlapping or adjacent FeatureLocations into one.
///
/// The merged location spans from the minimum start to the maximum end.
/// Both locations must have compatible strands (same or one None).
/// Returns a new FeatureLocation with the combined span.
///
/// Raises LocationParserError if the strands are incompatible.
pub fn merge_locations(
  loc1 : FeatureLocation,
  loc2 : FeatureLocation,
) -> FeatureLocation raise LocationParserError {
  let new_strand = match (loc1.strand, loc2.strand) {
    (Some(s1), Some(s2)) =>
      if s1 == s2 {
        Some(s1)
      } else {
        raise LocationParserError(
          "Cannot merge locations with different strands: " +
          s1.to_string() +
          " and " +
          s2.to_string(),
        )
      }
    (Some(s1), None) => Some(s1)
    (None, Some(s2)) => Some(s2)
    (None, None) => None
  }
  let start1 = loc1.start.to_int()
  let start2 = loc2.start.to_int()
  let end1 = loc1.end.to_int()
  let end2 = loc2.end.to_int()
  let new_start = if start1 <= start2 {
    loc1.start
  } else {
    loc2.start
  }
  let new_end = if end1 >= end2 {
    loc1.end
  } else {
    loc2.end
  }
  FeatureLocation::{
    start: new_start,
    end: new_end,
    strand: new_strand,
    ref_seq: loc1.ref_seq,
    ref_db: loc1.ref_db,
  }
}