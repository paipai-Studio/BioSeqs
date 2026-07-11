///|
/// Error raised when parsing location strings fails.
pub suberror LocationParserError {
  LocationParserError(String)
}

///|
/// Position types for feature locations.
/// These represent different types of fuzzy positions.
pub enum Position {
  ExactPosition(Int)
  WithinPosition(Int, Int, Int)
  BeforePosition(Int)
  AfterPosition(Int)
  UncertainPosition(Int)
  UnknownPosition
} derive(Eq)

///|
pub fn Position::new_exact(pos : Int) -> Position {
  ExactPosition(pos)
}

///|
pub fn Position::new_within(
  pos : Int,
  fuzzy_start : Int,
  fuzzy_end : Int,
) -> Position raise LocationParserError {
  if pos != fuzzy_start && pos != fuzzy_end {
    raise LocationParserError(
      "WithinPosition: " +
      pos.to_string() +
      " should match left " +
      fuzzy_start.to_string() +
      " or right " +
      fuzzy_end.to_string(),
    )
  }
  WithinPosition(pos, fuzzy_start, fuzzy_end)
}

///|
pub fn Position::new_before(pos : Int) -> Position {
  BeforePosition(pos)
}

///|
pub fn Position::new_after(pos : Int) -> Position {
  AfterPosition(pos)
}

///|
pub fn Position::new_uncertain(pos : Int) -> Position {
  UncertainPosition(pos)
}

///|
pub fn Position::new_unknown() -> Position {
  UnknownPosition
}

///|
pub fn Position::to_int(self : Position) -> Int {
  match self {
    ExactPosition(p) => p
    WithinPosition(p, _, _) => p
    BeforePosition(p) => p
    AfterPosition(p) => p
    UncertainPosition(p) => p
    UnknownPosition => -1
  }
}

///|
pub fn Position::repr(self : Position) -> String {
  match self {
    ExactPosition(p) => "ExactPosition(" + p.to_string() + ")"
    WithinPosition(p, fs, fe) =>
      "WithinPosition(" +
      p.to_string() +
      ", " +
      fs.to_string() +
      ", " +
      fe.to_string() +
      ")"
    BeforePosition(p) => "BeforePosition(" + p.to_string() + ")"
    AfterPosition(p) => "AfterPosition(" + p.to_string() + ")"
    UncertainPosition(p) => "UncertainPosition(" + p.to_string() + ")"
    UnknownPosition => "UnknownPosition()"
  }
}

///|
pub fn Position::str(self : Position) -> String {
  match self {
    ExactPosition(p) => p.to_string()
    WithinPosition(p, fs, fe) =>
      "<" + fs.to_string() + ".." + p.to_string() + ".." + fe.to_string() + ">"
    BeforePosition(p) => "<" + p.to_string()
    AfterPosition(p) => ">" + p.to_string()
    UncertainPosition(p) => "?" + p.to_string() + "?"
    UnknownPosition => "?"
  }
}

///|
/// FeatureLocation: specifies the start and end position of a feature.
/// Supports strand (+1, -1, 0, or None) and reference information.
pub struct FeatureLocation {
  start : Position
  end : Position
  strand : Int?
  ref_seq : String?
  ref_db : String?
} derive(Eq)

///|
pub fn FeatureLocation::new(
  start : Position,
  end : Position,
  strand? : Int? = None,
  ref_seq? : String? = None,
  ref_db? : String? = None,
) -> FeatureLocation raise LocationParserError {
  if strand != None {
    let s = strand.unwrap()
    if s != 1 && s != -1 && s != 0 {
      raise LocationParserError(
        "Strand should be +1, -1, 0 or None, not " + s.to_string(),
      )
    }
  }
  FeatureLocation::{ start, end, strand, ref_seq, ref_db }
}

///|
pub fn FeatureLocation::new_safe(
  start : Position,
  end : Position,
  strand? : Int? = None,
  ref_seq? : String? = None,
  ref_db? : String? = None,
) -> FeatureLocation {
  FeatureLocation::{ start, end, strand, ref_seq, ref_db }
}

///|
pub fn FeatureLocation::new_exact(
  start : Int,
  end : Int,
  strand? : Int? = None,
) -> FeatureLocation raise LocationParserError {
  if strand != None {
    let s = strand.unwrap()
    if s != 1 && s != -1 && s != 0 {
      raise LocationParserError(
        "Strand should be +1, -1, 0 or None, not " + s.to_string(),
      )
    }
  }
  FeatureLocation::{
    start: ExactPosition(start),
    end: ExactPosition(end),
    strand,
    ref_seq: None,
    ref_db: None,
  }
}

///|
pub fn FeatureLocation::start(self : FeatureLocation) -> Position {
  self.start
}

///|
pub fn FeatureLocation::end(self : FeatureLocation) -> Position {
  self.end
}

///|
pub fn FeatureLocation::strand(self : FeatureLocation) -> Int? {
  self.strand
}

///|
pub fn FeatureLocation::ref_seq(self : FeatureLocation) -> String? {
  self.ref_seq
}

///|
pub fn FeatureLocation::ref_db(self : FeatureLocation) -> String? {
  self.ref_db
}

///|
pub fn FeatureLocation::start_int(self : FeatureLocation) -> Int {
  self.start.to_int()
}

///|
pub fn FeatureLocation::end_int(self : FeatureLocation) -> Int {
  self.end.to_int()
}

///|
pub fn FeatureLocation::repr(self : FeatureLocation) -> String {
  let buf = StringBuilder::new()
  buf.write_string("FeatureLocation(")
  buf.write_string(self.start.repr())
  buf.write_string(", ")
  buf.write_string(self.end.repr())
  if self.strand != None {
    buf.write_string(", strand=")
    buf.write_string(self.strand.unwrap().to_string())
  }
  buf.write_char(')')
  buf.to_string()
}

///|
pub fn FeatureLocation::str(self : FeatureLocation) -> String {
  let buf = StringBuilder::new()
  buf.write_char('[')
  buf.write_string(self.start.str())
  buf.write_char(':')
  buf.write_string(self.end.str())
  buf.write_char(']')
  if self.strand != None {
    let s = self.strand.unwrap()
    buf.write_char('(')
    if s == 1 {
      buf.write_char('+')
    } else if s == -1 {
      buf.write_char('-')
    } else {
      buf.write_string(s.to_string())
    }
    buf.write_char(')')
  }
  buf.to_string()
}

///|
pub fn FeatureLocation::extract(
  self : FeatureLocation,
  parent_sequence : Seq,
) -> Seq raise LocationParserError {
  let start = self.start_int()
  let end = self.end_int()
  if start < 0 || end < 0 {
    raise LocationParserError("Cannot extract with negative position")
  }
  let seq_str = parent_sequence.to_string()
  let seq_len = seq_str.length()
  if start > seq_len || end > seq_len {
    raise LocationParserError("Position out of bounds")
  }
  let sliced = Seq::new(parent_sequence.to_string()[start:end].to_owned())
  if self.strand == Some(-1) {
    sliced.reverse_complement()
  } else {
    sliced
  }
}

///|
/// CompoundLocation: collection of FeatureLocation objects (for joins etc).
pub struct CompoundLocation {
  parts : Array[FeatureLocation]
  operator : String
  ref_seq : String?
  ref_db : String?
} derive(Eq)

///|
pub fn CompoundLocation::new(
  parts : Array[FeatureLocation],
  operator? : String = "join",
  ref_seq? : String? = None,
  ref_db? : String? = None,
) -> CompoundLocation {
  CompoundLocation::{ parts, operator, ref_seq, ref_db }
}

///|
pub fn CompoundLocation::parts(
  self : CompoundLocation,
) -> Array[FeatureLocation] {
  self.parts
}

///|
pub fn CompoundLocation::operator(self : CompoundLocation) -> String {
  self.operator
}

///|
pub fn CompoundLocation::repr(self : CompoundLocation) -> String {
  let buf = StringBuilder::new()
  buf.write_string("CompoundLocation(")
  buf.write_string(self.operator)
  buf.write_char(',')
  for i, part in self.parts {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(part.repr())
  }
  buf.write_char(')')
  buf.to_string()
}

///|
pub fn CompoundLocation::str(self : CompoundLocation) -> String {
  let buf = StringBuilder::new()
  buf.write_string(self.operator)
  buf.write_char('(')
  for i, part in self.parts {
    if i > 0 {
      buf.write_string(",")
    }
    buf.write_string(part.str())
  }
  buf.write_char(')')
  buf.to_string()
}

///|
/// SeqFeature: represents a sequence feature on an object.
pub struct SeqFeature {
  location : FeatureLocation?
  feature_type : String
  id : String
  qualifiers : Map[String, Array[String]]
  sub_features : Array[SeqFeature]
  ref_seq : String?
  ref_db : String?
  strand : Int?
} derive(Eq)

///|
pub fn SeqFeature::new(
  location? : FeatureLocation? = None,
  feature_type? : String = "",
  id? : String = "<unknown id>",
  qualifiers? : Map[String, Array[String]] = Map([], capacity=0),
  sub_features? : Array[SeqFeature] = [],
  strand? : Int? = None,
) -> SeqFeature {
  let loc_ref = match location {
    Some(l) => l.ref_seq()
    None => None
  }
  let loc_db = match location {
    Some(l) => l.ref_db()
    None => None
  }
  let strand_val = match strand {
    Some(s) => strand
    None =>
      match location {
        Some(l) => l.strand()
        None => None
      }
  }
  SeqFeature::{
    location,
    feature_type,
    id,
    qualifiers,
    sub_features,
    ref_seq: loc_ref,
    ref_db: loc_db,
    strand: strand_val,
  }
}

///|
pub fn SeqFeature::location(self : SeqFeature) -> FeatureLocation? {
  self.location
}

///|
pub fn SeqFeature::feature_type(self : SeqFeature) -> String {
  self.feature_type
}

///|
pub fn SeqFeature::id(self : SeqFeature) -> String {
  self.id
}

///|
pub fn SeqFeature::qualifiers(self : SeqFeature) -> Map[String, Array[String]] {
  self.qualifiers
}

///|
pub fn SeqFeature::sub_features(self : SeqFeature) -> Array[SeqFeature] {
  self.sub_features
}

///|
pub fn SeqFeature::strand(self : SeqFeature) -> Int? {
  self.strand
}

///|
pub fn SeqFeature::ref_seq(self : SeqFeature) -> String? {
  match self.location {
    Some(l) => l.ref_seq()
    None => self.ref_seq
  }
}

///|
pub fn SeqFeature::ref_db(self : SeqFeature) -> String? {
  match self.location {
    Some(l) => l.ref_db()
    None => self.ref_db
  }
}

///|
pub fn SeqFeature::repr(self : SeqFeature) -> String {
  let buf = StringBuilder::new()
  buf.write_string("SeqFeature(")
  match self.location {
    Some(l) => buf.write_string(l.repr())
    None => buf.write_string("None")
  }
  buf.write_string(", type='")
  buf.write_string(self.feature_type)
  buf.write_char('\'')
  buf.write_char(')')
  buf.to_string()
}

///|
pub fn SeqFeature::str(self : SeqFeature) -> String {
  let buf = StringBuilder::new()
  buf.write_string(self.feature_type)
  match self.location {
    Some(l) => {
      buf.write_char(' ')
      buf.write_string(l.str())
    }
    None => ()
  }
  for name, values in self.qualifiers {
    for val in values {
      buf.write_char('\n')
      buf.write_char('/')
      buf.write_string(name)
      buf.write_string("=\"")
      buf.write_string(val)
      buf.write_char('"')
    }
  }
  buf.to_string()
}

///|
pub fn SeqFeature::extract(
  self : SeqFeature,
  parent_sequence : Seq,
) -> Seq raise LocationParserError {
  match self.location {
    Some(l) => {
      let extracted = l.extract(parent_sequence)
      if self.strand == Some(-1) {
        extracted.reverse_complement()
      } else {
        extracted
      }
    }
    None => raise LocationParserError("The feature's .location is None")
  }
}

///|
pub fn SeqFeature::add_qualifier(
  self : SeqFeature,
  name : String,
  value : String,
) -> SeqFeature {
  let new_quals = Map([], capacity=self.qualifiers.length() + 1)
  for k, v in self.qualifiers {
    new_quals[k] = v
  }
  let existing = new_quals.get(name)
  match existing {
    Some(arr) => {
      let new_arr = arr.copy()
      new_arr.push(value)
      new_quals[name] = new_arr
    }
    None => new_quals[name] = [value]
  }
  SeqFeature::{ ..self, qualifiers: new_quals }
}

///|
/// Reference: holds information about a sequence reference.
pub struct Reference {
  authors : String
  title : String
  journal : String
  pubmed_id : String
  comment : String
} derive(Eq)

///|
pub fn Reference::new(
  authors? : String = "",
  title? : String = "",
  journal? : String = "",
  pubmed_id? : String = "",
  comment? : String = "",
) -> Reference {
  Reference::{ authors, title, journal, pubmed_id, comment }
}

///|
pub fn Reference::repr(self : Reference) -> String {
  "Reference(authors='" + self.authors + "', title='" + self.title + "')"
}
