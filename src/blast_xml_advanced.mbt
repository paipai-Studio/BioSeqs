// Biopython-compatible modern Bio.Blast XML support.
//
// This module parses and writes NCBI BLAST XML1 and XML2 artifacts. It does
// not submit remote BLAST jobs or launch local BLAST executables.

///|
/// Error raised for malformed or inconsistent BLAST XML.
pub suberror BlastXmlError {
  BlastXmlError(String)
}

///|
/// NCBI BLAST XML dialect.
pub(all) enum BlastXmlFormat {
  Xml1
  Xml2
} derive(Eq, Debug)
pub extend BlastXmlFormat with Eq::{equal, not_equal}
pub extend BlastXmlFormat with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One masked interval on the query, using zero-based half-open coordinates.
pub struct BlastXmlMask {
  start : Int
  end : Int
} derive(Eq, Debug)
pub extend BlastXmlMask with Eq::{equal, not_equal}
pub extend BlastXmlMask with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Query metadata for one BLAST record.
pub struct BlastXmlQuery {
  id : String
  description : String
  length : Int
  sequence : String?
  masks : Array[BlastXmlMask]
} derive(Eq, Debug)
pub extend BlastXmlQuery with Eq::{equal, not_equal}
pub extend BlastXmlQuery with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Search parameters stored in a BLAST XML header.
pub struct BlastXmlParameters {
  matrix : String?
  expect : Double
  inclusion_threshold : Double?
  score_match : Int?
  score_mismatch : Int?
  gap_open : Int?
  gap_extend : Int?
  filter : String?
  pattern : String?
  entrez_query : String?
  composition_based_statistics : Int?
  query_genetic_code : Int?
  database_genetic_code : Int?
  bl2seq_mode : Int?
} derive(Eq, Debug)
pub extend BlastXmlParameters with Eq::{equal, not_equal}
pub extend BlastXmlParameters with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Database and Karlin-Altschul statistics for one query.
///
/// `database_letters` is a `Double` because BLAST databases can exceed the
/// platform `Int` range. Integer values below 2^53 remain exactly represented.
pub struct BlastXmlStatistics {
  database_sequences : Int
  database_letters : Double
  effective_hsp_length : Int
  effective_search_space : Double
  kappa : Double
  lambda : Double
  entropy : Double
} derive(Eq, Debug)
pub extend BlastXmlStatistics with Eq::{equal, not_equal}
pub extend BlastXmlStatistics with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One target description. XML2 can attach several descriptions to one hit.
pub struct BlastXmlDescription {
  id : String
  accession : String
  title : String
  taxid : Int?
  scientific_name : String?
} derive(Eq, Debug)
pub extend BlastXmlDescription with Eq::{equal, not_equal}
pub extend BlastXmlDescription with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One breakpoint in an HSP coordinate path.
///
/// Coordinates are zero-based boundaries. Reverse-strand axes decrease.
/// Translated nucleotide axes move in steps of three.
pub struct BlastXmlCoordinate {
  target : Int
  query : Int
} derive(Eq, Debug)
pub extend BlastXmlCoordinate with Eq::{equal, not_equal}
pub extend BlastXmlCoordinate with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One high-scoring segment pair.
pub struct BlastXmlHsp {
  number : Int
  bit_score : Double
  score : Double
  evalue : Double
  query_from : Int
  query_to : Int
  target_from : Int
  target_to : Int
  query_frame : Int?
  target_frame : Int?
  query_strand : String?
  target_strand : String?
  identity : Int
  positive : Int?
  gaps : Int?
  alignment_length : Int?
  density : Int?
  pattern_from : Int?
  pattern_to : Int?
  query_sequence : String
  target_sequence : String
  midline : String?
  coordinates : Array[BlastXmlCoordinate]
  query_coded_by : String?
  target_coded_by : String?
} derive(Eq, Debug)
pub extend BlastXmlHsp with Eq::{equal, not_equal}
pub extend BlastXmlHsp with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One BLAST database hit and its HSPs.
pub struct BlastXmlHit {
  number : Int
  descriptions : Array[BlastXmlDescription]
  length : Int
  hsps : Array[BlastXmlHsp]
} derive(Eq, Debug)
pub extend BlastXmlHit with Eq::{equal, not_equal}
pub extend BlastXmlHit with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Results for one query.
pub struct BlastXmlRecord {
  number : Int
  query : BlastXmlQuery
  hits : Array[BlastXmlHit]
  statistics : BlastXmlStatistics?
  message : String?
} derive(Eq, Debug)
pub extend BlastXmlRecord with Eq::{equal, not_equal}
pub extend BlastXmlRecord with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// A complete BLAST XML1 or XML2 document.
pub struct BlastXmlDocument {
  format : BlastXmlFormat
  program : String
  version : String
  reference : String
  database : String
  global_query : BlastXmlQuery?
  parameters : BlastXmlParameters
  records : Array[BlastXmlRecord]
  megablast_statistics : BlastXmlStatistics?
} derive(Eq, Debug)
pub extend BlastXmlDocument with Eq::{equal, not_equal}
pub extend BlastXmlDocument with @moonbitlang/core/debug.Debug::{to_repr}

///|
priv struct BlastXmlNode {
  name : String
  text : String
  children : Array[BlastXmlNode]
}

///|
fn blast_xml_fail(message : String) -> Unit raise BlastXmlError {
  raise BlastXmlError(message)
}

///|
fn blast_xml_copy_masks(values : Array[BlastXmlMask]) -> Array[BlastXmlMask] {
  let copied : Array[BlastXmlMask] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn blast_xml_copy_descriptions(
  values : Array[BlastXmlDescription],
) -> Array[BlastXmlDescription] {
  let copied : Array[BlastXmlDescription] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn blast_xml_copy_coordinates(
  values : Array[BlastXmlCoordinate],
) -> Array[BlastXmlCoordinate] {
  let copied : Array[BlastXmlCoordinate] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn blast_xml_copy_hsps(values : Array[BlastXmlHsp]) -> Array[BlastXmlHsp] {
  let copied : Array[BlastXmlHsp] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn blast_xml_copy_hits(values : Array[BlastXmlHit]) -> Array[BlastXmlHit] {
  let copied : Array[BlastXmlHit] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn blast_xml_copy_records(
  values : Array[BlastXmlRecord],
) -> Array[BlastXmlRecord] {
  let copied : Array[BlastXmlRecord] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn blast_xml_copy_query(query : BlastXmlQuery) -> BlastXmlQuery {
  BlastXmlQuery::{
    id: query.id,
    description: query.description,
    length: query.length,
    sequence: query.sequence,
    masks: blast_xml_copy_masks(query.masks),
  }
}

///|
pub fn BlastXmlMask::start(self : BlastXmlMask) -> Int {
  self.start
}

///|
pub fn BlastXmlMask::end(self : BlastXmlMask) -> Int {
  self.end
}

///|
pub fn BlastXmlQuery::id(self : BlastXmlQuery) -> String {
  self.id
}

///|
pub fn BlastXmlQuery::description(self : BlastXmlQuery) -> String {
  self.description
}

///|
pub fn BlastXmlQuery::length(self : BlastXmlQuery) -> Int {
  self.length
}

///|
pub fn BlastXmlQuery::sequence(self : BlastXmlQuery) -> String? {
  self.sequence
}

///|
pub fn BlastXmlQuery::masks(self : BlastXmlQuery) -> Array[BlastXmlMask] {
  blast_xml_copy_masks(self.masks)
}

///|
pub fn BlastXmlParameters::matrix(self : BlastXmlParameters) -> String? {
  self.matrix
}

///|
pub fn BlastXmlParameters::expect(self : BlastXmlParameters) -> Double {
  self.expect
}

///|
pub fn BlastXmlParameters::inclusion_threshold(
  self : BlastXmlParameters,
) -> Double? {
  self.inclusion_threshold
}

///|
pub fn BlastXmlParameters::score_match(self : BlastXmlParameters) -> Int? {
  self.score_match
}

///|
pub fn BlastXmlParameters::score_mismatch(self : BlastXmlParameters) -> Int? {
  self.score_mismatch
}

///|
pub fn BlastXmlParameters::gap_open(self : BlastXmlParameters) -> Int? {
  self.gap_open
}

///|
pub fn BlastXmlParameters::gap_extend(self : BlastXmlParameters) -> Int? {
  self.gap_extend
}

///|
pub fn BlastXmlParameters::filter(self : BlastXmlParameters) -> String? {
  self.filter
}

///|
pub fn BlastXmlParameters::pattern(self : BlastXmlParameters) -> String? {
  self.pattern
}

///|
pub fn BlastXmlParameters::entrez_query(self : BlastXmlParameters) -> String? {
  self.entrez_query
}

///|
pub fn BlastXmlParameters::composition_based_statistics(
  self : BlastXmlParameters,
) -> Int? {
  self.composition_based_statistics
}

///|
pub fn BlastXmlParameters::query_genetic_code(
  self : BlastXmlParameters,
) -> Int? {
  self.query_genetic_code
}

///|
pub fn BlastXmlParameters::database_genetic_code(
  self : BlastXmlParameters,
) -> Int? {
  self.database_genetic_code
}

///|
pub fn BlastXmlParameters::bl2seq_mode(self : BlastXmlParameters) -> Int? {
  self.bl2seq_mode
}

///|
pub fn BlastXmlStatistics::database_sequences(self : BlastXmlStatistics) -> Int {
  self.database_sequences
}

///|
pub fn BlastXmlStatistics::database_letters(
  self : BlastXmlStatistics,
) -> Double {
  self.database_letters
}

///|
pub fn BlastXmlStatistics::effective_hsp_length(
  self : BlastXmlStatistics,
) -> Int {
  self.effective_hsp_length
}

///|
pub fn BlastXmlStatistics::effective_search_space(
  self : BlastXmlStatistics,
) -> Double {
  self.effective_search_space
}

///|
pub fn BlastXmlStatistics::kappa(self : BlastXmlStatistics) -> Double {
  self.kappa
}

///|
pub fn BlastXmlStatistics::lambda(self : BlastXmlStatistics) -> Double {
  self.lambda
}

///|
pub fn BlastXmlStatistics::entropy(self : BlastXmlStatistics) -> Double {
  self.entropy
}

///|
pub fn BlastXmlDescription::id(self : BlastXmlDescription) -> String {
  self.id
}

///|
pub fn BlastXmlDescription::accession(self : BlastXmlDescription) -> String {
  self.accession
}

///|
pub fn BlastXmlDescription::title(self : BlastXmlDescription) -> String {
  self.title
}

///|
pub fn BlastXmlDescription::taxid(self : BlastXmlDescription) -> Int? {
  self.taxid
}

///|
pub fn BlastXmlDescription::scientific_name(
  self : BlastXmlDescription,
) -> String? {
  self.scientific_name
}

///|
pub fn BlastXmlCoordinate::target(self : BlastXmlCoordinate) -> Int {
  self.target
}

///|
pub fn BlastXmlCoordinate::query(self : BlastXmlCoordinate) -> Int {
  self.query
}

///|
pub fn BlastXmlHsp::number(self : BlastXmlHsp) -> Int {
  self.number
}

///|
pub fn BlastXmlHsp::bit_score(self : BlastXmlHsp) -> Double {
  self.bit_score
}

///|
pub fn BlastXmlHsp::score(self : BlastXmlHsp) -> Double {
  self.score
}

///|
pub fn BlastXmlHsp::evalue(self : BlastXmlHsp) -> Double {
  self.evalue
}

///|
pub fn BlastXmlHsp::query_from(self : BlastXmlHsp) -> Int {
  self.query_from
}

///|
pub fn BlastXmlHsp::query_to(self : BlastXmlHsp) -> Int {
  self.query_to
}

///|
pub fn BlastXmlHsp::target_from(self : BlastXmlHsp) -> Int {
  self.target_from
}

///|
pub fn BlastXmlHsp::target_to(self : BlastXmlHsp) -> Int {
  self.target_to
}

///|
pub fn BlastXmlHsp::query_frame(self : BlastXmlHsp) -> Int? {
  self.query_frame
}

///|
pub fn BlastXmlHsp::target_frame(self : BlastXmlHsp) -> Int? {
  self.target_frame
}

///|
pub fn BlastXmlHsp::query_strand(self : BlastXmlHsp) -> String? {
  self.query_strand
}

///|
pub fn BlastXmlHsp::target_strand(self : BlastXmlHsp) -> String? {
  self.target_strand
}

///|
pub fn BlastXmlHsp::identity(self : BlastXmlHsp) -> Int {
  self.identity
}

///|
pub fn BlastXmlHsp::positive(self : BlastXmlHsp) -> Int? {
  self.positive
}

///|
pub fn BlastXmlHsp::gaps(self : BlastXmlHsp) -> Int? {
  self.gaps
}

///|
pub fn BlastXmlHsp::alignment_length(self : BlastXmlHsp) -> Int? {
  self.alignment_length
}

///|
pub fn BlastXmlHsp::density(self : BlastXmlHsp) -> Int? {
  self.density
}

///|
pub fn BlastXmlHsp::query_sequence(self : BlastXmlHsp) -> String {
  self.query_sequence
}

///|
pub fn BlastXmlHsp::target_sequence(self : BlastXmlHsp) -> String {
  self.target_sequence
}

///|
pub fn BlastXmlHsp::midline(self : BlastXmlHsp) -> String? {
  self.midline
}

///|
pub fn BlastXmlHsp::coordinates(
  self : BlastXmlHsp,
) -> Array[BlastXmlCoordinate] {
  blast_xml_copy_coordinates(self.coordinates)
}

///|
pub fn BlastXmlHsp::query_coded_by(self : BlastXmlHsp) -> String? {
  self.query_coded_by
}

///|
pub fn BlastXmlHsp::target_coded_by(self : BlastXmlHsp) -> String? {
  self.target_coded_by
}

///|
pub fn BlastXmlHsp::identity_fraction(self : BlastXmlHsp) -> Double {
  match self.alignment_length {
    Some(length) =>
      if length == 0 {
        0.0
      } else {
        self.identity.to_double() / length.to_double()
      }
    None => 0.0
  }
}

///|
pub fn BlastXmlHit::number(self : BlastXmlHit) -> Int {
  self.number
}

///|
pub fn BlastXmlHit::descriptions(
  self : BlastXmlHit,
) -> Array[BlastXmlDescription] {
  blast_xml_copy_descriptions(self.descriptions)
}

///|
pub fn BlastXmlHit::primary_description(
  self : BlastXmlHit,
) -> BlastXmlDescription {
  self.descriptions[0]
}

///|
pub fn BlastXmlHit::length(self : BlastXmlHit) -> Int {
  self.length
}

///|
pub fn BlastXmlHit::hsps(self : BlastXmlHit) -> Array[BlastXmlHsp] {
  blast_xml_copy_hsps(self.hsps)
}

///|
pub fn BlastXmlHit::best_hsp(self : BlastXmlHit) -> BlastXmlHsp? {
  if self.hsps.length() == 0 {
    return None
  }
  let mut best = self.hsps[0]
  for i = 1; i < self.hsps.length(); i = i + 1 {
    let candidate = self.hsps[i]
    if candidate.evalue < best.evalue ||
      (candidate.evalue == best.evalue && candidate.bit_score > best.bit_score) {
      best = candidate
    }
  }
  Some(best)
}

///|
pub fn BlastXmlRecord::number(self : BlastXmlRecord) -> Int {
  self.number
}

///|
pub fn BlastXmlRecord::query(self : BlastXmlRecord) -> BlastXmlQuery {
  blast_xml_copy_query(self.query)
}

///|
pub fn BlastXmlRecord::hits(self : BlastXmlRecord) -> Array[BlastXmlHit] {
  blast_xml_copy_hits(self.hits)
}

///|
pub fn BlastXmlRecord::statistics(self : BlastXmlRecord) -> BlastXmlStatistics? {
  self.statistics
}

///|
pub fn BlastXmlRecord::message(self : BlastXmlRecord) -> String? {
  self.message
}

///|
pub fn BlastXmlRecord::hit(
  self : BlastXmlRecord,
  identifier : String,
) -> BlastXmlHit? {
  for hit in self.hits {
    for description in hit.descriptions {
      if description.id == identifier || description.accession == identifier {
        return Some(hit)
      }
    }
  }
  None
}

///|
pub fn BlastXmlRecord::best_hit(self : BlastXmlRecord) -> BlastXmlHit? {
  if self.hits.length() == 0 {
    return None
  }
  let mut best_hit = self.hits[0]
  let mut best_hsp = best_hit.best_hsp()
  for i = 1; i < self.hits.length(); i = i + 1 {
    let candidate_hsp = self.hits[i].best_hsp()
    match (best_hsp, candidate_hsp) {
      (None, Some(value)) => {
        best_hit = self.hits[i]
        best_hsp = Some(value)
      }
      (Some(current), Some(candidate)) =>
        if candidate.evalue < current.evalue ||
          (
            candidate.evalue == current.evalue &&
            candidate.bit_score > current.bit_score
          ) {
          best_hit = self.hits[i]
          best_hsp = Some(candidate)
        }
      _ => ()
    }
  }
  Some(best_hit)
}

///|
pub fn BlastXmlRecord::filter_evalue(
  self : BlastXmlRecord,
  maximum : Double,
) -> BlastXmlRecord raise BlastXmlError {
  if maximum < 0.0 || maximum.abs() > 1.0e300 {
    blast_xml_fail("BLAST XML E-value threshold must be finite and nonnegative")
  }
  let hits : Array[BlastXmlHit] = []
  for hit in self.hits {
    let hsps : Array[BlastXmlHsp] = []
    for hsp in hit.hsps {
      if hsp.evalue <= maximum {
        hsps.push(hsp)
      }
    }
    if hsps.length() > 0 {
      hits.push(BlastXmlHit::{
        number: hit.number,
        descriptions: blast_xml_copy_descriptions(hit.descriptions),
        length: hit.length,
        hsps,
      })
    }
  }
  BlastXmlRecord::{
    number: self.number,
    query: blast_xml_copy_query(self.query),
    hits,
    statistics: self.statistics,
    message: self.message,
  }
}

///|
pub fn BlastXmlDocument::format(self : BlastXmlDocument) -> BlastXmlFormat {
  self.format
}

///|
pub fn BlastXmlDocument::program(self : BlastXmlDocument) -> String {
  self.program
}

///|
pub fn BlastXmlDocument::version(self : BlastXmlDocument) -> String {
  self.version
}

///|
pub fn BlastXmlDocument::reference(self : BlastXmlDocument) -> String {
  self.reference
}

///|
pub fn BlastXmlDocument::database(self : BlastXmlDocument) -> String {
  self.database
}

///|
pub fn BlastXmlDocument::global_query(
  self : BlastXmlDocument,
) -> BlastXmlQuery? {
  match self.global_query {
    Some(query) => Some(blast_xml_copy_query(query))
    None => None
  }
}

///|
pub fn BlastXmlDocument::parameters(
  self : BlastXmlDocument,
) -> BlastXmlParameters {
  self.parameters
}

///|
pub fn BlastXmlDocument::records(
  self : BlastXmlDocument,
) -> Array[BlastXmlRecord] {
  blast_xml_copy_records(self.records)
}

///|
pub fn BlastXmlDocument::megablast_statistics(
  self : BlastXmlDocument,
) -> BlastXmlStatistics? {
  self.megablast_statistics
}

///|
pub fn BlastXmlDocument::record(
  self : BlastXmlDocument,
  query_id : String,
) -> BlastXmlRecord? {
  for record in self.records {
    if record.query.id == query_id {
      return Some(record)
    }
  }
  None
}

///|
pub fn BlastXmlDocument::total_hits(self : BlastXmlDocument) -> Int {
  let mut total = 0
  for record in self.records {
    total = total + record.hits.length()
  }
  total
}

///|
pub fn BlastXmlDocument::total_hsps(self : BlastXmlDocument) -> Int {
  let mut total = 0
  for record in self.records {
    for hit in record.hits {
      total = total + hit.hsps.length()
    }
  }
  total
}

///|
fn blast_xml_is_space(code : Int) -> Bool {
  code == 32 || code == 9 || code == 10 || code == 13
}

///|
fn blast_xml_trim(value : String) -> String {
  let mut start = 0
  let mut end = value.length()
  while start < end && blast_xml_is_space(value.unsafe_get(start).to_int()) {
    start = start + 1
  }
  while end > start && blast_xml_is_space(value.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  value[start:end].to_owned()
}

///|
fn blast_xml_starts_at(text : String, index : Int, token : String) -> Bool {
  if index < 0 || index + token.length() > text.length() {
    return false
  }
  for i = 0; i < token.length(); i = i + 1 {
    if text.unsafe_get(index + i) != token.unsafe_get(i) {
      return false
    }
  }
  true
}

///|
fn blast_xml_find(text : String, token : String, start : Int) -> Int {
  if token.length() == 0 {
    return start
  }
  let mut index = start
  while index + token.length() <= text.length() {
    if blast_xml_starts_at(text, index, token) {
      return index
    }
    index = index + 1
  }
  -1
}

///|
fn blast_xml_parse_radix_digits(value : String, radix : Int) -> Int? {
  if value.length() == 0 {
    return None
  }
  let mut result = 0
  for i = 0; i < value.length(); i = i + 1 {
    let code = value.unsafe_get(i).to_int()
    let digit = if code >= 48 && code <= 57 {
      code - 48
    } else if code >= 65 && code <= 70 {
      code - 55
    } else if code >= 97 && code <= 102 {
      code - 87
    } else {
      return None
    }
    if digit >= radix {
      return None
    }
    result = result * radix + digit
  }
  Some(result)
}

///|
fn blast_xml_character(code : Int) -> Char {
  code.unsafe_to_char()
}

///|
fn blast_xml_entity_character(entity : String) -> Char? {
  match entity {
    "lt" => Some('<')
    "gt" => Some('>')
    "amp" => Some('&')
    "quot" => Some('"')
    "apos" => Some('\'')
    "nbsp" => Some(blast_xml_character(160))
    "Auml" => Some(blast_xml_character(196))
    "auml" => Some(blast_xml_character(228))
    "Ouml" => Some(blast_xml_character(214))
    "ouml" => Some(blast_xml_character(246))
    "Uuml" => Some(blast_xml_character(220))
    "uuml" => Some(blast_xml_character(252))
    "szlig" => Some(blast_xml_character(223))
    "eacute" => Some(blast_xml_character(233))
    _ =>
      if entity.length() > 1 && entity.unsafe_get(0).to_int() == 35 {
        let digits = entity[1:].to_owned()
        let parsed = if digits.length() > 1 &&
          (
            digits.unsafe_get(0).to_int() == 120 ||
            digits.unsafe_get(0).to_int() == 88
          ) {
          blast_xml_parse_radix_digits(digits[1:].to_owned(), 16)
        } else {
          blast_xml_parse_radix_digits(digits, 10)
        }
        match parsed {
          Some(code) =>
            if code > 0 && code <= 0x10ffff {
              Some(code.unsafe_to_char())
            } else {
              None
            }
          None => None
        }
      } else {
        None
      }
  }
}

///|
fn blast_xml_unescape_once(value : String) -> String {
  let output = StringBuilder::new()
  let mut index = 0
  while index < value.length() {
    if value.unsafe_get(index).to_int() == 38 {
      let semicolon = blast_xml_find(value, ";", index + 1)
      if semicolon >= 0 {
        let entity = value[index + 1:semicolon].to_owned()
        match blast_xml_entity_character(entity) {
          Some(character) => {
            output.write_char(character)
            index = semicolon + 1
            continue
          }
          None => ()
        }
      }
    }
    output.write_char(value.unsafe_get(index).unsafe_to_char())
    index = index + 1
  }
  output.to_string()
}

///|
fn blast_xml_unescape(value : String) -> String {
  blast_xml_unescape_once(blast_xml_unescape_once(value))
}

///|
fn blast_xml_normalize_tag(raw : String) -> String {
  let mut value = raw
  let colon = blast_xml_find(value, ":", 0)
  if colon >= 0 {
    value = value[colon + 1:].to_owned()
  }
  let mut last_underscore = -1
  for i = 0; i < value.length(); i = i + 1 {
    if value.unsafe_get(i).to_int() == 95 {
      last_underscore = i
    }
  }
  if last_underscore >= 0 {
    value = value[last_underscore + 1:].to_owned()
  }
  value.to_lower()
}

///|
fn blast_xml_open_end(text : String, start : Int) -> Int raise BlastXmlError {
  let mut quote = 0
  let mut index = start
  while index < text.length() {
    let code = text.unsafe_get(index).to_int()
    if quote == 0 && (code == 34 || code == 39) {
      quote = code
    } else if quote == code {
      quote = 0
    } else if quote == 0 && code == 62 {
      return index
    }
    index = index + 1
  }
  blast_xml_fail("Unterminated BLAST XML opening tag")
  0
}

///|
fn blast_xml_tag_name(opening : String) -> String raise BlastXmlError {
  let trimmed = blast_xml_trim(opening)
  let mut end = 0
  while end < trimmed.length() &&
        !blast_xml_is_space(trimmed.unsafe_get(end).to_int()) &&
        trimmed.unsafe_get(end).to_int() != 47 {
    end = end + 1
  }
  if end == 0 {
    blast_xml_fail("BLAST XML contains an empty tag name")
  }
  blast_xml_normalize_tag(trimmed[0:end].to_owned())
}

///|
fn blast_xml_skip_declaration(
  text : String,
  start : Int,
) -> Int raise BlastXmlError {
  if blast_xml_starts_at(text, start, "<!--") {
    let end = blast_xml_find(text, "-->", start + 4)
    if end < 0 {
      blast_xml_fail("Unterminated BLAST XML comment")
    }
    return end + 3
  }
  if blast_xml_starts_at(text, start, "<?") {
    let end = blast_xml_find(text, "?>", start + 2)
    if end < 0 {
      blast_xml_fail("Unterminated BLAST XML processing instruction")
    }
    return end + 2
  }
  if blast_xml_starts_at(text, start, "<!DOCTYPE") {
    let mut index = start + 9
    let mut quote = 0
    let mut brackets = 0
    while index < text.length() {
      let code = text.unsafe_get(index).to_int()
      if quote == 0 && (code == 34 || code == 39) {
        quote = code
      } else if quote == code {
        quote = 0
      } else if quote == 0 && code == 91 {
        brackets = brackets + 1
      } else if quote == 0 && code == 93 {
        brackets = brackets - 1
      } else if quote == 0 && brackets == 0 && code == 62 {
        return index + 1
      }
      index = index + 1
    }
    blast_xml_fail("Unterminated BLAST XML DOCTYPE")
  }
  blast_xml_fail("Unsupported BLAST XML declaration")
  0
}

///|
fn blast_xml_parse_node(
  text : String,
  start : Int,
) -> (BlastXmlNode, Int) raise BlastXmlError {
  if start >= text.length() || text.unsafe_get(start).to_int() != 60 {
    blast_xml_fail("Expected a BLAST XML opening tag")
  }
  if blast_xml_starts_at(text, start, "</") {
    blast_xml_fail("Unexpected BLAST XML closing tag")
  }
  let open_end = blast_xml_open_end(text, start + 1)
  let opening = text[start + 1:open_end].to_owned()
  let name = blast_xml_tag_name(opening)
  let trimmed_opening = blast_xml_trim(opening)
  let self_closing = trimmed_opening.length() > 0 &&
    trimmed_opening.unsafe_get(trimmed_opening.length() - 1).to_int() == 47
  if self_closing {
    return (BlastXmlNode::{ name, text: "", children: [], }, open_end + 1)
  }
  let children : Array[BlastXmlNode] = []
  let content = StringBuilder::new()
  let mut index = open_end + 1
  while index < text.length() {
    if blast_xml_starts_at(text, index, "</") {
      let close_end = blast_xml_open_end(text, index + 2)
      let close_name = blast_xml_tag_name(text[index + 2:close_end].to_owned())
      if close_name != name {
        blast_xml_fail(
          "Mismatched BLAST XML closing tag: expected " +
          name +
          ", found " +
          close_name,
        )
      }
      return (
        BlastXmlNode::{
          name,
          text: blast_xml_unescape_once(content.to_string()),
          children,
        },
        close_end + 1,
      )
    }
    if blast_xml_starts_at(text, index, "<![CDATA[") {
      let end = blast_xml_find(text, "]]>", index + 9)
      if end < 0 {
        blast_xml_fail("Unterminated BLAST XML CDATA section")
      }
      content.write_string(text[index + 9:end].to_owned())
      index = end + 3
      continue
    }
    if blast_xml_starts_at(text, index, "<!--") ||
      blast_xml_starts_at(text, index, "<?") ||
      blast_xml_starts_at(text, index, "<!DOCTYPE") {
      index = blast_xml_skip_declaration(text, index)
      continue
    }
    if text.unsafe_get(index).to_int() == 60 {
      let (child, next) = blast_xml_parse_node(text, index)
      children.push(child)
      index = next
      continue
    }
    content.write_char(text.unsafe_get(index).unsafe_to_char())
    index = index + 1
  }
  blast_xml_fail("Premature end of BLAST XML document")
  (BlastXmlNode::{ name: "", text: "", children: [], }, 0)
}

///|
fn blast_xml_parse_tree(text : String) -> BlastXmlNode raise BlastXmlError {
  let mut index = 0
  while index < text.length() {
    while index < text.length() &&
          blast_xml_is_space(text.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index >= text.length() {
      blast_xml_fail("BLAST XML input is empty")
    }
    if blast_xml_starts_at(text, index, "<?") ||
      blast_xml_starts_at(text, index, "<!--") ||
      blast_xml_starts_at(text, index, "<!DOCTYPE") {
      index = blast_xml_skip_declaration(text, index)
    } else {
      break
    }
  }
  let (root, next) = blast_xml_parse_node(text, index)
  index = next
  while index < text.length() {
    while index < text.length() &&
          blast_xml_is_space(text.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index < text.length() &&
      (
        blast_xml_starts_at(text, index, "<!--") ||
        blast_xml_starts_at(text, index, "<?")
      ) {
      index = blast_xml_skip_declaration(text, index)
    } else {
      break
    }
  }
  if index != text.length() {
    blast_xml_fail("Unexpected data after BLAST XML root element")
  }
  root
}

///|
fn blast_xml_children_named(
  node : BlastXmlNode,
  name : String,
) -> Array[BlastXmlNode] {
  let matches : Array[BlastXmlNode] = []
  for child in node.children {
    if child.name == name {
      matches.push(child)
    }
  }
  matches
}

///|
fn blast_xml_optional_child(
  node : BlastXmlNode,
  name : String,
) -> BlastXmlNode? raise BlastXmlError {
  let matches = blast_xml_children_named(node, name)
  if matches.length() > 1 {
    blast_xml_fail("Duplicate BLAST XML element: " + name)
  }
  if matches.length() == 0 {
    None
  } else {
    Some(matches[0])
  }
}

///|
fn blast_xml_required_child(
  node : BlastXmlNode,
  name : String,
) -> BlastXmlNode raise BlastXmlError {
  match blast_xml_optional_child(node, name) {
    Some(child) => child
    None => {
      blast_xml_fail("Missing BLAST XML element: " + name)
      BlastXmlNode::{ name: "", text: "", children: [], }
    }
  }
}

///|
fn blast_xml_node_text(
  node : BlastXmlNode,
  label : String,
) -> String raise BlastXmlError {
  if node.children.length() > 0 {
    blast_xml_fail("BLAST XML " + label + " must contain text only")
  }
  blast_xml_trim(node.text)
}

///|
fn blast_xml_required_text(
  node : BlastXmlNode,
  name : String,
) -> String raise BlastXmlError {
  let value = blast_xml_node_text(blast_xml_required_child(node, name), name)
  if value.length() == 0 {
    blast_xml_fail("BLAST XML " + name + " cannot be empty")
  }
  value
}

///|
fn blast_xml_optional_text(
  node : BlastXmlNode,
  name : String,
) -> String? raise BlastXmlError {
  match blast_xml_optional_child(node, name) {
    Some(child) => Some(blast_xml_node_text(child, name))
    None => None
  }
}

///|
fn blast_xml_parse_number(
  value : String,
  label : String,
) -> Double raise BlastXmlError {
  let number : Double = @string.parse_double(blast_xml_trim(value)) catch {
    _ => {
      blast_xml_fail("BLAST XML " + label + " must be numeric")
      0.0
    }
  }
  if number.is_nan() || number.is_inf() {
    blast_xml_fail("BLAST XML " + label + " must be finite")
  }
  number
}

///|
fn blast_xml_parse_integer(
  value : String,
  label : String,
) -> Int raise BlastXmlError {
  let number = blast_xml_parse_number(value, label)
  let integer = number.to_int()
  if number != integer.to_double() {
    blast_xml_fail("BLAST XML " + label + " must be an integer")
  }
  integer
}

///|
fn blast_xml_required_int(
  node : BlastXmlNode,
  name : String,
) -> Int raise BlastXmlError {
  blast_xml_parse_integer(blast_xml_required_text(node, name), name)
}

///|
fn blast_xml_optional_int(
  node : BlastXmlNode,
  name : String,
) -> Int? raise BlastXmlError {
  match blast_xml_optional_text(node, name) {
    Some(value) => Some(blast_xml_parse_integer(value, name))
    None => None
  }
}

///|
fn blast_xml_required_double(
  node : BlastXmlNode,
  name : String,
) -> Double raise BlastXmlError {
  blast_xml_parse_number(blast_xml_required_text(node, name), name)
}

///|
fn blast_xml_optional_double(
  node : BlastXmlNode,
  name : String,
) -> Double? raise BlastXmlError {
  match blast_xml_optional_text(node, name) {
    Some(value) => Some(blast_xml_parse_number(value, name))
    None => None
  }
}

///|
fn blast_xml_descendants(
  node : BlastXmlNode,
  name : String,
  output : Array[BlastXmlNode],
) -> Unit {
  for child in node.children {
    if child.name == name {
      output.push(child)
    }
    blast_xml_descendants(child, name, output)
  }
}

///|
fn blast_xml_first_descendant(
  node : BlastXmlNode,
  name : String,
) -> BlastXmlNode? {
  for child in node.children {
    if child.name == name {
      return Some(child)
    }
    match blast_xml_first_descendant(child, name) {
      Some(found) => return Some(found)
      None => ()
    }
  }
  None
}

///|
fn blast_xml_unwrap(node : BlastXmlNode, name : String) -> BlastXmlNode {
  let mut current = node
  while current.children.length() == 1 && current.children[0].name == name {
    current = current.children[0]
  }
  current
}

///|
fn blast_xml_parse_query(
  node : BlastXmlNode,
  fallback : BlastXmlQuery?,
) -> BlastXmlQuery raise BlastXmlError {
  let id = match blast_xml_optional_text(node, "query-id") {
    Some(value) =>
      if value.length() == 0 {
        blast_xml_fail("BLAST XML query-id cannot be empty")
        ""
      } else {
        value
      }
    None =>
      match fallback {
        Some(query) => query.id
        None => {
          blast_xml_fail("BLAST XML record is missing query-id")
          ""
        }
      }
  }
  let description = match blast_xml_optional_text(node, "query-def") {
    Some(value) => value
    None =>
      match blast_xml_optional_text(node, "query-title") {
        Some(value) => value
        None =>
          match fallback {
            Some(query) => query.description
            None => ""
          }
      }
  }
  let length = match blast_xml_optional_int(node, "query-len") {
    Some(value) => value
    None =>
      match fallback {
        Some(query) => query.length
        None => {
          blast_xml_fail("BLAST XML record is missing query-len")
          0
        }
      }
  }
  if length <= 0 {
    blast_xml_fail("BLAST XML query-len must be positive")
  }
  let sequence : String? = match blast_xml_optional_text(node, "query-seq") {
    Some(value) =>
      if value.length() != length {
        blast_xml_fail("BLAST XML query-seq length does not match query-len")
        None
      } else {
        Some(value)
      }
    None =>
      match fallback {
        Some(query) => if query.id == id { query.sequence } else { None }
        None => None
      }
  }
  let masks : Array[BlastXmlMask] = []
  let masking_nodes : Array[BlastXmlNode] = []
  blast_xml_descendants(node, "query-masking", masking_nodes)
  for masking in masking_nodes {
    let range_nodes : Array[BlastXmlNode] = []
    blast_xml_descendants(masking, "range", range_nodes)
    if range_nodes.length() == 0 {
      blast_xml_fail("BLAST XML query-masking is missing Range")
    }
    for range in range_nodes {
      let start = blast_xml_required_int(range, "from") - 1
      let end = blast_xml_required_int(range, "to")
      if start < 0 || end <= start || end > length {
        blast_xml_fail("BLAST XML query-masking range is out of bounds")
      }
      masks.push(BlastXmlMask::{ start, end, })
    }
  }
  BlastXmlQuery::{ id, description, length, sequence, masks, }
}

///|
fn blast_xml_parse_parameters(
  container : BlastXmlNode,
) -> BlastXmlParameters raise BlastXmlError {
  let node = match blast_xml_first_descendant(container, "parameters") {
    Some(value) => value
    None =>
      if container.name == "parameters" {
        container
      } else {
        blast_xml_fail("BLAST XML header is missing Parameters")
        container
      }
  }
  let expect = blast_xml_required_double(node, "expect")
  if expect < 0.0 {
    blast_xml_fail("BLAST XML expect must be nonnegative")
  }
  let inclusion_threshold = blast_xml_optional_double(node, "include")
  match inclusion_threshold {
    Some(value) =>
      if value < 0.0 {
        blast_xml_fail("BLAST XML include must be nonnegative")
      }
    None => ()
  }
  let gap_open = blast_xml_optional_int(node, "gap-open")
  let gap_extend = blast_xml_optional_int(node, "gap-extend")
  match gap_open {
    Some(value) =>
      if value < 0 {
        blast_xml_fail("BLAST XML gap-open must be nonnegative")
      }
    None => ()
  }
  match gap_extend {
    Some(value) =>
      if value < 0 {
        blast_xml_fail("BLAST XML gap-extend must be nonnegative")
      }
    None => ()
  }
  BlastXmlParameters::{
    matrix: blast_xml_optional_text(node, "matrix"),
    expect,
    inclusion_threshold,
    score_match: blast_xml_optional_int(node, "sc-match"),
    score_mismatch: blast_xml_optional_int(node, "sc-mismatch"),
    gap_open,
    gap_extend,
    filter: blast_xml_optional_text(node, "filter"),
    pattern: blast_xml_optional_text(node, "pattern"),
    entrez_query: blast_xml_optional_text(node, "entrez-query"),
    composition_based_statistics: blast_xml_optional_int(node, "cbs"),
    query_genetic_code: blast_xml_optional_int(node, "query-gencode"),
    database_genetic_code: blast_xml_optional_int(node, "db-gencode"),
    bl2seq_mode: blast_xml_optional_int(node, "bl2seq-mode"),
  }
}

///|
fn blast_xml_parse_statistics(
  node : BlastXmlNode,
) -> BlastXmlStatistics raise BlastXmlError {
  let statistics = if node.name == "statistics" {
    node
  } else {
    match blast_xml_first_descendant(node, "statistics") {
      Some(value) => value
      None => {
        blast_xml_fail("BLAST XML statistics container is empty")
        node
      }
    }
  }
  let database_sequences = blast_xml_required_int(statistics, "db-num")
  let database_letters = blast_xml_required_double(statistics, "db-len")
  let effective_hsp_length = blast_xml_required_int(statistics, "hsp-len")
  let effective_search_space = blast_xml_required_double(
    statistics, "eff-space",
  )
  let kappa = blast_xml_required_double(statistics, "kappa")
  let lambda = blast_xml_required_double(statistics, "lambda")
  let entropy = blast_xml_required_double(statistics, "entropy")
  if database_sequences < 0 ||
    database_letters < 0.0 ||
    effective_hsp_length < 0 ||
    effective_search_space < 0.0 ||
    kappa < 0.0 ||
    lambda < 0.0 ||
    entropy < 0.0 {
    blast_xml_fail("BLAST XML statistics cannot be negative")
  }
  BlastXmlStatistics::{
    database_sequences,
    database_letters,
    effective_hsp_length,
    effective_search_space,
    kappa,
    lambda,
    entropy,
  }
}

///|
fn blast_xml_program_supported(program : String) -> Bool {
  match program {
    "blastn"
    | "megablast"
    | "blastp"
    | "blastx"
    | "tblastn"
    | "tblastx"
    | "rpsblast"
    | "psiblast" => true
    _ => false
  }
}

///|
fn blast_xml_query_translated(program : String) -> Bool {
  program == "blastx" || program == "tblastx"
}

///|
fn blast_xml_target_translated(program : String) -> Bool {
  program == "tblastn" || program == "tblastx"
}

///|
fn blast_xml_validate_frame(
  program : String,
  frame : Int?,
  query_axis : Bool,
) -> Unit raise BlastXmlError {
  match frame {
    None => ()
    Some(value) => {
      let valid = if query_axis {
        match program {
          "blastn" | "megablast" => value == -1 || value == 1
          "blastx" | "tblastx" => value >= -3 && value <= 3 && value != 0
          "blastp" | "tblastn" | "rpsblast" | "psiblast" => value == 0
          _ => false
        }
      } else {
        match program {
          "blastn" | "megablast" => value == -1 || value == 1
          "tblastn" | "tblastx" => value >= -3 && value <= 3 && value != 0
          "blastp" | "blastx" | "rpsblast" | "psiblast" => value == 0
          _ => false
        }
      }
      if !valid {
        blast_xml_fail(
          "BLAST XML contains an invalid frame for program " + program,
        )
      }
    }
  }
}

///|
fn blast_xml_validate_strand(
  strand : String?,
  label : String,
) -> Unit raise BlastXmlError {
  match strand {
    Some(value) =>
      if value != "Plus" && value != "Minus" {
        blast_xml_fail("BLAST XML " + label + " must be Plus or Minus")
      }
    None => ()
  }
}

///|
fn blast_xml_axis_direction(
  translated : Bool,
  from : Int,
  to : Int,
  frame : Int?,
  strand : String?,
) -> Int {
  if translated {
    match frame {
      Some(value) => if value < 0 { -1 } else { 1 }
      None => if from > to { -1 } else { 1 }
    }
  } else {
    match strand {
      Some("Minus") => -1
      Some(_) => 1
      None =>
        match frame {
          Some(value) => if value < 0 { -1 } else if from > to { -1 } else { 1 }
          None => if from > to { -1 } else { 1 }
        }
    }
  }
}

///|
fn blast_xml_axis_start(from : Int, to : Int, direction : Int) -> Int {
  if direction > 0 {
    if from < to {
      from - 1
    } else {
      to - 1
    }
  } else if from > to {
    from
  } else {
    to
  }
}

///|
fn blast_xml_axis_expected_end(from : Int, to : Int, direction : Int) -> Int {
  if direction > 0 {
    if from > to {
      from
    } else {
      to
    }
  } else if from < to {
    from - 1
  } else {
    to - 1
  }
}

///|
fn blast_xml_ungapped_length(value : String) -> Int {
  let mut length = 0
  for i = 0; i < value.length(); i = i + 1 {
    if value.unsafe_get(i).to_int() != 45 {
      length = length + 1
    }
  }
  length
}

///|
fn blast_xml_validate_aligned_sequence(
  value : String,
  label : String,
) -> Unit raise BlastXmlError {
  for i = 0; i < value.length(); i = i + 1 {
    let code = value.unsafe_get(i).to_int()
    if blast_xml_is_space(code) {
      blast_xml_fail("BLAST XML " + label + " contains whitespace")
    }
  }
}

///|
fn blast_xml_build_coordinates(
  query_sequence : String,
  target_sequence : String,
  query_from : Int,
  query_to : Int,
  target_from : Int,
  target_to : Int,
  query_step : Int,
  target_step : Int,
  query_direction : Int,
  target_direction : Int,
) -> Array[BlastXmlCoordinate] raise BlastXmlError {
  let coordinates : Array[BlastXmlCoordinate] = []
  let mut query_coordinate = blast_xml_axis_start(
    query_from, query_to, query_direction,
  )
  let mut target_coordinate = blast_xml_axis_start(
    target_from, target_to, target_direction,
  )
  coordinates.push(BlastXmlCoordinate::{
    target: target_coordinate,
    query: query_coordinate,
  })
  let mut previous_query_consumed = false
  let mut previous_target_consumed = false
  for column = 0; column < query_sequence.length(); column = column + 1 {
    let query_consumed = query_sequence.unsafe_get(column).to_int() != 45
    let target_consumed = target_sequence.unsafe_get(column).to_int() != 45
    if !query_consumed && !target_consumed {
      blast_xml_fail("BLAST XML HSP contains a double-gap column")
    }
    if column > 0 &&
      (
        query_consumed != previous_query_consumed ||
        target_consumed != previous_target_consumed
      ) {
      coordinates.push(BlastXmlCoordinate::{
        target: target_coordinate,
        query: query_coordinate,
      })
    }
    if query_consumed {
      query_coordinate = query_coordinate + query_direction * query_step
    }
    if target_consumed {
      target_coordinate = target_coordinate + target_direction * target_step
    }
    previous_query_consumed = query_consumed
    previous_target_consumed = target_consumed
  }
  coordinates.push(BlastXmlCoordinate::{
    target: target_coordinate,
    query: query_coordinate,
  })
  let expected_query = blast_xml_axis_expected_end(
    query_from, query_to, query_direction,
  )
  let expected_target = blast_xml_axis_expected_end(
    target_from, target_to, target_direction,
  )
  if query_coordinate != expected_query {
    blast_xml_fail(
      "BLAST XML query coordinates do not match aligned sequence consumption",
    )
  }
  if target_coordinate != expected_target {
    blast_xml_fail(
      "BLAST XML target coordinates do not match aligned sequence consumption",
    )
  }
  coordinates
}

///|
fn blast_xml_coded_by(
  identifier : String,
  from : Int,
  to : Int,
  direction : Int,
) -> String {
  let lower = if from < to { from } else { to }
  let upper = if from < to { to } else { from }
  let location = identifier + ":" + lower.to_string() + ".." + upper.to_string()
  if direction < 0 {
    "complement(" + location + ")"
  } else {
    location
  }
}

///|
fn blast_xml_parse_hsp(
  node : BlastXmlNode,
  program : String,
  query : BlastXmlQuery,
  target : BlastXmlDescription,
  target_length : Int,
) -> BlastXmlHsp raise BlastXmlError {
  let number = blast_xml_required_int(node, "num")
  let bit_score = blast_xml_required_double(node, "bit-score")
  let score = blast_xml_required_double(node, "score")
  let evalue = blast_xml_required_double(node, "evalue")
  if number <= 0 || bit_score < 0.0 || score < 0.0 || evalue < 0.0 {
    blast_xml_fail("BLAST XML HSP number and scores must be nonnegative")
  }
  let query_from = blast_xml_required_int(node, "query-from")
  let query_to = blast_xml_required_int(node, "query-to")
  let target_from = blast_xml_required_int(node, "hit-from")
  let target_to = blast_xml_required_int(node, "hit-to")
  let query_frame = blast_xml_optional_int(node, "query-frame")
  let target_frame = blast_xml_optional_int(node, "hit-frame")
  let query_strand = blast_xml_optional_text(node, "query-strand")
  let target_strand = blast_xml_optional_text(node, "hit-strand")
  blast_xml_validate_frame(program, query_frame, true)
  blast_xml_validate_frame(program, target_frame, false)
  blast_xml_validate_strand(query_strand, "query-strand")
  blast_xml_validate_strand(target_strand, "hit-strand")
  let identity = blast_xml_required_int(node, "identity")
  let positive = blast_xml_optional_int(node, "positive")
  let gaps = blast_xml_optional_int(node, "gaps")
  let alignment_length = blast_xml_optional_int(node, "align-len")
  let density = blast_xml_optional_int(node, "density")
  let pattern_from = blast_xml_optional_int(node, "pattern-from")
  let pattern_to = blast_xml_optional_int(node, "pattern-to")
  if identity < 0 {
    blast_xml_fail("BLAST XML HSP identity cannot be negative")
  }
  match positive {
    Some(value) =>
      if value < identity {
        blast_xml_fail("BLAST XML HSP positive count is below identity")
      }
    None => ()
  }
  match gaps {
    Some(value) =>
      if value < 0 {
        blast_xml_fail("BLAST XML HSP gaps cannot be negative")
      }
    None => ()
  }
  let query_sequence = match blast_xml_optional_text(node, "qseq") {
    Some(value) => value
    None => ""
  }
  let target_sequence = match blast_xml_optional_text(node, "hseq") {
    Some(value) => value
    None => ""
  }
  let midline = blast_xml_optional_text(node, "midline")
  blast_xml_validate_aligned_sequence(query_sequence, "qseq")
  blast_xml_validate_aligned_sequence(target_sequence, "hseq")
  let coordinates : Array[BlastXmlCoordinate] = []
  let mut query_coded_by : String? = None
  let mut target_coded_by : String? = None
  match alignment_length {
    Some(length) => {
      if length <= 0 ||
        query_sequence.length() != length ||
        target_sequence.length() != length {
        blast_xml_fail("BLAST XML HSP sequence lengths must match align-len")
      }
      match midline {
        Some(value) =>
          if value.length() != length {
            blast_xml_fail("BLAST XML HSP midline length must match align-len")
          }
        None => ()
      }
      if identity > length {
        blast_xml_fail("BLAST XML HSP identity exceeds align-len")
      }
      match positive {
        Some(value) =>
          if value > length {
            blast_xml_fail("BLAST XML HSP positive count exceeds align-len")
          }
        None => ()
      }
      let mut observed_gaps = 0
      for i = 0; i < length; i = i + 1 {
        if query_sequence.unsafe_get(i).to_int() == 45 ||
          target_sequence.unsafe_get(i).to_int() == 45 {
          observed_gaps = observed_gaps + 1
        }
      }
      match gaps {
        Some(value) =>
          if value != observed_gaps {
            blast_xml_fail(
              "BLAST XML HSP gaps does not match aligned sequences",
            )
          }
        None => ()
      }
      if query_from <= 0 ||
        query_to <= 0 ||
        target_from <= 0 ||
        target_to <= 0 ||
        query_from > query.length ||
        query_to > query.length ||
        target_from > target_length ||
        target_to > target_length {
        blast_xml_fail("BLAST XML HSP coordinates are out of bounds")
      }
      let query_translated = blast_xml_query_translated(program)
      let target_translated = blast_xml_target_translated(program)
      if query_translated && query_frame is None {
        blast_xml_fail("BLAST XML translated query is missing query-frame")
      }
      if target_translated && target_frame is None {
        blast_xml_fail("BLAST XML translated target is missing hit-frame")
      }
      if (program == "blastn" || program == "megablast") &&
        (
          (query_frame is None && query_strand is None) ||
          (target_frame is None && target_strand is None)
        ) {
        blast_xml_fail("BLAST XML nucleotide HSP is missing strand or frame")
      }
      if !query_translated &&
        program != "blastn" &&
        program != "megablast" &&
        query_from > query_to {
        blast_xml_fail("BLAST XML protein query coordinates must increase")
      }
      if !target_translated &&
        program != "blastn" &&
        program != "megablast" &&
        target_from > target_to {
        blast_xml_fail("BLAST XML protein target coordinates must increase")
      }
      let query_step = if query_translated { 3 } else { 1 }
      let target_step = if target_translated { 3 } else { 1 }
      let query_direction = blast_xml_axis_direction(
        query_translated, query_from, query_to, query_frame, query_strand,
      )
      let target_direction = blast_xml_axis_direction(
        target_translated, target_from, target_to, target_frame, target_strand,
      )
      if query_from != query_to &&
        (
          (query_direction > 0 && query_from > query_to) ||
          (query_direction < 0 && query_from < query_to)
        ) {
        blast_xml_fail("BLAST XML query direction conflicts with coordinates")
      }
      if target_from != target_to &&
        (
          (target_direction > 0 && target_from > target_to) ||
          (target_direction < 0 && target_from < target_to)
        ) {
        blast_xml_fail("BLAST XML target direction conflicts with coordinates")
      }
      let query_span = (query_to - query_from).abs() + 1
      let target_span = (target_to - target_from).abs() + 1
      if query_span != blast_xml_ungapped_length(query_sequence) * query_step {
        blast_xml_fail("BLAST XML query span does not match qseq residue count")
      }
      if target_span != blast_xml_ungapped_length(target_sequence) * target_step {
        blast_xml_fail(
          "BLAST XML target span does not match hseq residue count",
        )
      }
      let built = blast_xml_build_coordinates(
        query_sequence, target_sequence, query_from, query_to, target_from, target_to,
        query_step, target_step, query_direction, target_direction,
      )
      for coordinate in built {
        coordinates.push(coordinate)
      }
      if query_translated {
        query_coded_by = Some(
          blast_xml_coded_by(query.id, query_from, query_to, query_direction),
        )
      }
      if target_translated {
        target_coded_by = Some(
          blast_xml_coded_by(
            target.id,
            target_from,
            target_to,
            target_direction,
          ),
        )
      }
    }
    None =>
      if query_sequence.length() > 0 ||
        target_sequence.length() > 0 ||
        query_from != 0 ||
        query_to != 0 ||
        target_from != 0 ||
        target_to != 0 {
        blast_xml_fail(
          "BLAST XML HSP without align-len must have empty alignment data",
        )
      }
  }
  BlastXmlHsp::{
    number,
    bit_score,
    score,
    evalue,
    query_from,
    query_to,
    target_from,
    target_to,
    query_frame,
    target_frame,
    query_strand,
    target_strand,
    identity,
    positive,
    gaps,
    alignment_length,
    density,
    pattern_from,
    pattern_to,
    query_sequence,
    target_sequence,
    midline,
    coordinates,
    query_coded_by,
    target_coded_by,
  }
}

///|
fn blast_xml_parse_description(
  node : BlastXmlNode,
  xml2 : Bool,
) -> BlastXmlDescription raise BlastXmlError {
  let id = blast_xml_required_text(node, "id")
  let accession = blast_xml_required_text(node, "accession")
  let title = if xml2 {
    blast_xml_required_text(node, "title")
  } else {
    blast_xml_required_text(node, "def")
  }
  let taxid = blast_xml_optional_int(node, "taxid")
  match taxid {
    Some(value) =>
      if value <= 0 {
        blast_xml_fail("BLAST XML taxid must be positive")
      }
    None => ()
  }
  BlastXmlDescription::{
    id,
    accession,
    title,
    taxid,
    scientific_name: blast_xml_optional_text(node, "sciname"),
  }
}

///|
fn blast_xml_parse_hit(
  node : BlastXmlNode,
  program : String,
  query : BlastXmlQuery,
  xml2 : Bool,
) -> BlastXmlHit raise BlastXmlError {
  let number = blast_xml_required_int(node, "num")
  let length = blast_xml_required_int(node, "len")
  if number <= 0 || length <= 0 {
    blast_xml_fail("BLAST XML hit number and length must be positive")
  }
  let descriptions : Array[BlastXmlDescription] = []
  if xml2 {
    let container = blast_xml_required_child(node, "description")
    let description_nodes = blast_xml_children_named(container, "hitdescr")
    if description_nodes.length() == 0 {
      blast_xml_fail("BLAST XML2 hit has no HitDescr entries")
    }
    for description_node in description_nodes {
      descriptions.push(blast_xml_parse_description(description_node, true))
    }
  } else {
    descriptions.push(blast_xml_parse_description(node, false))
  }
  for i = 0; i < descriptions.length(); i = i + 1 {
    for j = i + 1; j < descriptions.length(); j = j + 1 {
      if descriptions[i].id == descriptions[j].id {
        blast_xml_fail("BLAST XML hit contains duplicate target identifiers")
      }
    }
  }
  let hsps : Array[BlastXmlHsp] = []
  let container = blast_xml_required_child(node, "hsps")
  let hsp_nodes = blast_xml_children_named(container, "hsp")
  for hsp_node in hsp_nodes {
    let hsp = blast_xml_parse_hsp(
      hsp_node,
      program,
      query,
      descriptions[0],
      length,
    )
    for previous in hsps {
      if previous.number == hsp.number {
        blast_xml_fail("BLAST XML hit contains duplicate HSP numbers")
      }
    }
    hsps.push(hsp)
  }
  BlastXmlHit::{ number, descriptions, length, hsps, }
}

///|
fn blast_xml_find_statistics(
  node : BlastXmlNode,
) -> BlastXmlStatistics? raise BlastXmlError {
  match blast_xml_optional_child(node, "stat") {
    Some(stat) => Some(blast_xml_parse_statistics(stat))
    None =>
      match blast_xml_optional_child(node, "statistics") {
        Some(statistics) => Some(blast_xml_parse_statistics(statistics))
        None => None
      }
  }
}

///|
fn blast_xml_parse_record(
  node : BlastXmlNode,
  program : String,
  xml2 : Bool,
  fallback_query : BlastXmlQuery?,
  default_number : Int,
) -> BlastXmlRecord raise BlastXmlError {
  let number = match blast_xml_optional_int(node, "iter-num") {
    Some(value) => value
    None =>
      match blast_xml_optional_int(node, "num") {
        Some(value) => value
        None => default_number
      }
  }
  if number <= 0 {
    blast_xml_fail("BLAST XML record number must be positive")
  }
  let query = blast_xml_parse_query(node, fallback_query)
  let hits : Array[BlastXmlHit] = []
  match blast_xml_optional_child(node, "hits") {
    Some(container) =>
      for hit_node in blast_xml_children_named(container, "hit") {
        let hit = blast_xml_parse_hit(hit_node, program, query, xml2)
        for previous in hits {
          if previous.number == hit.number {
            blast_xml_fail("BLAST XML record contains duplicate hit numbers")
          }
        }
        hits.push(hit)
      }
    None => ()
  }
  BlastXmlRecord::{
    number,
    query,
    hits,
    statistics: blast_xml_find_statistics(node),
    message: blast_xml_optional_text(node, "message"),
  }
}

///|
fn blast_xml_parse_xml1(
  root : BlastXmlNode,
) -> BlastXmlDocument raise BlastXmlError {
  let program = blast_xml_required_text(root, "program").to_lower()
  if !blast_xml_program_supported(program) {
    blast_xml_fail("Unsupported BLAST program: " + program)
  }
  let global_query = blast_xml_parse_query(root, None)
  let parameters = blast_xml_parse_parameters(
    blast_xml_required_child(root, "param"),
  )
  let records : Array[BlastXmlRecord] = []
  let iterations = blast_xml_required_child(root, "iterations")
  let iteration_nodes = blast_xml_children_named(iterations, "iteration")
  if iteration_nodes.length() == 0 {
    blast_xml_fail("BLAST XML1 document has no Iteration records")
  }
  for i = 0; i < iteration_nodes.length(); i = i + 1 {
    let record = blast_xml_parse_record(
      iteration_nodes[i],
      program,
      false,
      Some(global_query),
      i + 1,
    )
    for previous in records {
      if previous.number == record.number {
        blast_xml_fail("BLAST XML document contains duplicate record numbers")
      }
    }
    records.push(record)
  }
  let megablast_statistics = match blast_xml_optional_child(root, "mbstat") {
    Some(node) => Some(blast_xml_parse_statistics(node))
    None => None
  }
  BlastXmlDocument::{
    format: BlastXmlFormat::Xml1,
    program,
    version: blast_xml_required_text(root, "version"),
    reference: blast_xml_unescape(blast_xml_required_text(root, "reference")),
    database: blast_xml_required_text(root, "db"),
    global_query: Some(global_query),
    parameters,
    records,
    megablast_statistics,
  }
}

///|
fn blast_xml_record_candidates(
  node : BlastXmlNode,
  output : Array[BlastXmlNode],
) -> Unit {
  let has_query = blast_xml_children_named(node, "query-id").length() > 0
  let has_hits = blast_xml_children_named(node, "hits").length() > 0
  if node.name == "search" && (has_query || has_hits) {
    output.push(node)
    return
  }
  for child in node.children {
    blast_xml_record_candidates(child, output)
  }
}

///|
fn blast_xml_parse_xml2(
  root : BlastXmlNode,
) -> BlastXmlDocument raise BlastXmlError {
  let outputs : Array[BlastXmlNode] = []
  if root.name == "blastoutput2" {
    outputs.push(root)
  } else {
    blast_xml_descendants(root, "blastoutput2", outputs)
  }
  if outputs.length() == 0 {
    blast_xml_fail("BLAST XML2 document has no BlastOutput2 reports")
  }
  let records : Array[BlastXmlRecord] = []
  let mut program = ""
  let mut version = ""
  let mut reference = ""
  let mut database = ""
  let mut parameters : BlastXmlParameters? = None
  for output in outputs {
    let report_wrapper = blast_xml_required_child(output, "report")
    let report = blast_xml_unwrap(report_wrapper, "report")
    let current_program = blast_xml_required_text(report, "program").to_lower()
    if !blast_xml_program_supported(current_program) {
      blast_xml_fail("Unsupported BLAST program: " + current_program)
    }
    let current_version = blast_xml_required_text(report, "version")
    let current_reference = blast_xml_unescape(
      blast_xml_required_text(report, "reference"),
    )
    let target = match blast_xml_first_descendant(report, "target") {
      Some(value) => value
      None => {
        blast_xml_fail("BLAST XML2 report is missing search target")
        report
      }
    }
    let current_database = blast_xml_required_text(target, "db")
    let params_container = blast_xml_required_child(report, "params")
    let current_parameters = blast_xml_parse_parameters(params_container)
    if program.length() == 0 {
      program = current_program
      version = current_version
      reference = current_reference
      database = current_database
      parameters = Some(current_parameters)
    } else if program != current_program ||
      version != current_version ||
      reference != current_reference ||
      database != current_database ||
      parameters != Some(current_parameters) {
      blast_xml_fail("BLAST XML2 reports have inconsistent headers")
    }
    let candidates : Array[BlastXmlNode] = []
    blast_xml_record_candidates(report, candidates)
    if candidates.length() == 0 {
      blast_xml_fail("BLAST XML2 report has no Search record")
    }
    for candidate in candidates {
      let record = blast_xml_parse_record(
        candidate,
        current_program,
        true,
        None,
        records.length() + 1,
      )
      for previous in records {
        if previous.number == record.number {
          blast_xml_fail("BLAST XML document contains duplicate record numbers")
        }
      }
      records.push(record)
    }
  }
  BlastXmlDocument::{
    format: BlastXmlFormat::Xml2,
    program,
    version,
    reference,
    database,
    global_query: None,
    parameters: parameters.unwrap(),
    records,
    megablast_statistics: None,
  }
}

///|
/// Parse and strictly validate an NCBI BLAST XML1 or XML2 document.
pub fn blast_xml_parse(text : String) -> BlastXmlDocument raise BlastXmlError {
  let root = blast_xml_parse_tree(text)
  match root.name {
    "blastoutput" => blast_xml_parse_xml1(root)
    "blastxml2" | "blastoutput2" => blast_xml_parse_xml2(root)
    _ => {
      blast_xml_fail(
        "Expected BlastOutput, BlastXML2, or BlastOutput2 root element",
      )
      blast_xml_parse_xml1(root)
    }
  }
}

///|
fn blast_xml_escape(value : String) -> String raise BlastXmlError {
  let output = StringBuilder::new()
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    match code {
      38 => output.write_string("&amp;")
      60 => output.write_string("&lt;")
      62 => output.write_string("&gt;")
      34 => output.write_string("&quot;")
      39 => output.write_string("&apos;")
      9 | 10 | 13 => output.write_char(code.unsafe_to_char())
      _ =>
        if code < 32 {
          blast_xml_fail("BLAST XML text contains an illegal control character")
        } else {
          output.write_char(value.unsafe_get(index).unsafe_to_char())
        }
    }
  }
  output.to_string()
}

///|
fn blast_xml_write_element(
  output : StringBuilder,
  indent : String,
  name : String,
  value : String,
) -> Unit raise BlastXmlError {
  output.write_string(indent)
  output.write_string("<")
  output.write_string(name)
  output.write_string(">")
  output.write_string(blast_xml_escape(value))
  output.write_string("</")
  output.write_string(name)
  output.write_string(">\n")
}

///|
fn blast_xml_write_optional_string(
  output : StringBuilder,
  indent : String,
  name : String,
  value : String?,
) -> Unit raise BlastXmlError {
  match value {
    Some(text) => blast_xml_write_element(output, indent, name, text)
    None => ()
  }
}

///|
fn blast_xml_write_optional_int(
  output : StringBuilder,
  indent : String,
  name : String,
  value : Int?,
) -> Unit raise BlastXmlError {
  match value {
    Some(number) =>
      blast_xml_write_element(output, indent, name, number.to_string())
    None => ()
  }
}

///|
fn blast_xml_write_optional_double(
  output : StringBuilder,
  indent : String,
  name : String,
  value : Double?,
) -> Unit raise BlastXmlError {
  match value {
    Some(number) =>
      blast_xml_write_element(output, indent, name, number.to_string())
    None => ()
  }
}

///|
fn blast_xml_validate_write(
  document : BlastXmlDocument,
  format : BlastXmlFormat,
) -> Unit raise BlastXmlError {
  if !blast_xml_program_supported(document.program) {
    blast_xml_fail("Unsupported BLAST program: " + document.program)
  }
  if document.version.length() == 0 ||
    document.reference.length() == 0 ||
    document.database.length() == 0 {
    blast_xml_fail("BLAST XML header fields cannot be empty")
  }
  if document.records.length() == 0 {
    blast_xml_fail("BLAST XML document has no records")
  }
  match format {
    BlastXmlFormat::Xml1 => {
      if document.parameters.gap_open is None ||
        document.parameters.gap_extend is None {
        blast_xml_fail("XML1 requires gap-open and gap-extend parameters")
      }
      if document.parameters.composition_based_statistics is Some(_) ||
        document.parameters.query_genetic_code is Some(_) ||
        document.parameters.database_genetic_code is Some(_) ||
        document.parameters.bl2seq_mode is Some(_) {
        blast_xml_fail(
          "XML1 cannot represent XML2-only composition or genetic-code parameters",
        )
      }
      let global_query = match document.global_query {
        Some(query) => query
        None => document.records[0].query
      }
      if global_query.masks.length() > 0 {
        blast_xml_fail("XML1 cannot represent query masking")
      }
      for record in document.records {
        if record.query.masks.length() > 0 {
          blast_xml_fail("XML1 cannot represent query masking")
        }
        match record.query.sequence {
          Some(sequence) =>
            if global_query.id != record.query.id ||
              global_query.sequence != Some(sequence) {
              blast_xml_fail("XML1 cannot represent per-record query sequences")
            }
          None => ()
        }
      }
    }
    BlastXmlFormat::Xml2 => {
      if document.megablast_statistics is Some(_) {
        blast_xml_fail("XML2 cannot represent BlastOutput_mbstat")
      }
      for index = 0; index < document.records.length(); index = index + 1 {
        let record = document.records[index]
        if record.number != index + 1 {
          blast_xml_fail("XML2 record numbers must be consecutive")
        }
        if record.query.sequence is Some(_) {
          blast_xml_fail("XML2 cannot represent query sequences")
        }
      }
    }
  }
}

///|
fn blast_xml_write_parameters(
  output : StringBuilder,
  parameters : BlastXmlParameters,
  xml2 : Bool,
) -> Unit raise BlastXmlError {
  let indent = if xml2 { "          " } else { "      " }
  let prefix = if xml2 { "" } else { "Parameters_" }
  blast_xml_write_optional_string(
    output,
    indent,
    prefix + "matrix",
    parameters.matrix,
  )
  blast_xml_write_element(
    output,
    indent,
    prefix + "expect",
    parameters.expect.to_string(),
  )
  blast_xml_write_optional_double(
    output,
    indent,
    prefix + "include",
    parameters.inclusion_threshold,
  )
  blast_xml_write_optional_int(
    output,
    indent,
    prefix + "sc-match",
    parameters.score_match,
  )
  blast_xml_write_optional_int(
    output,
    indent,
    prefix + "sc-mismatch",
    parameters.score_mismatch,
  )
  blast_xml_write_optional_int(
    output,
    indent,
    prefix + "gap-open",
    parameters.gap_open,
  )
  blast_xml_write_optional_int(
    output,
    indent,
    prefix + "gap-extend",
    parameters.gap_extend,
  )
  blast_xml_write_optional_string(
    output,
    indent,
    prefix + "filter",
    parameters.filter,
  )
  blast_xml_write_optional_string(
    output,
    indent,
    prefix + "pattern",
    parameters.pattern,
  )
  blast_xml_write_optional_string(
    output,
    indent,
    prefix + "entrez-query",
    parameters.entrez_query,
  )
  if xml2 {
    blast_xml_write_optional_int(
      output,
      indent,
      "cbs",
      parameters.composition_based_statistics,
    )
    blast_xml_write_optional_int(
      output,
      indent,
      "query-gencode",
      parameters.query_genetic_code,
    )
    blast_xml_write_optional_int(
      output,
      indent,
      "db-gencode",
      parameters.database_genetic_code,
    )
    blast_xml_write_optional_int(
      output,
      indent,
      "bl2seq-mode",
      parameters.bl2seq_mode,
    )
  }
}

///|
fn blast_xml_write_statistics(
  output : StringBuilder,
  statistics : BlastXmlStatistics,
  xml2 : Bool,
  indent : String,
) -> Unit raise BlastXmlError {
  output.write_string(indent)
  output.write_string("<Statistics>\n")
  let child_indent = indent + "  "
  let prefix = if xml2 { "" } else { "Statistics_" }
  blast_xml_write_element(
    output,
    child_indent,
    prefix + "db-num",
    statistics.database_sequences.to_string(),
  )
  blast_xml_write_element(
    output,
    child_indent,
    prefix + "db-len",
    statistics.database_letters.to_string(),
  )
  blast_xml_write_element(
    output,
    child_indent,
    prefix + "hsp-len",
    statistics.effective_hsp_length.to_string(),
  )
  blast_xml_write_element(
    output,
    child_indent,
    prefix + "eff-space",
    statistics.effective_search_space.to_string(),
  )
  blast_xml_write_element(
    output,
    child_indent,
    prefix + "kappa",
    statistics.kappa.to_string(),
  )
  blast_xml_write_element(
    output,
    child_indent,
    prefix + "lambda",
    statistics.lambda.to_string(),
  )
  blast_xml_write_element(
    output,
    child_indent,
    prefix + "entropy",
    statistics.entropy.to_string(),
  )
  output.write_string(indent)
  output.write_string("</Statistics>\n")
}

///|
fn blast_xml_legacy_description(hit : BlastXmlHit) -> String {
  let output = StringBuilder::new()
  output.write_string(hit.descriptions[0].title)
  for index = 1; index < hit.descriptions.length(); index = index + 1 {
    let description = hit.descriptions[index]
    output.write_string(" >")
    output.write_string(description.id)
    output.write_string(" ")
    output.write_string(description.title)
  }
  output.to_string()
}

///|
fn blast_xml_xml1_frame(frame : Int?, strand : String?) -> Int? {
  match frame {
    Some(value) => Some(value)
    None =>
      match strand {
        Some("Plus") => Some(1)
        Some("Minus") => Some(-1)
        _ => None
      }
  }
}

///|
fn blast_xml_xml2_strand(strand : String?, frame : Int?) -> String? {
  match strand {
    Some(value) => Some(value)
    None =>
      match frame {
        Some(value) =>
          if value < 0 {
            Some("Minus")
          } else if value > 0 {
            Some("Plus")
          } else {
            None
          }
        None => None
      }
  }
}

///|
fn blast_xml_write_hsp_xml1(
  output : StringBuilder,
  hsp : BlastXmlHsp,
) -> Unit raise BlastXmlError {
  output.write_string("    <Hsp>\n")
  blast_xml_write_element(output, "      ", "Hsp_num", hsp.number.to_string())
  blast_xml_write_element(
    output,
    "      ",
    "Hsp_bit-score",
    hsp.bit_score.to_string(),
  )
  blast_xml_write_element(output, "      ", "Hsp_score", hsp.score.to_string())
  blast_xml_write_element(
    output,
    "      ",
    "Hsp_evalue",
    hsp.evalue.to_string(),
  )
  blast_xml_write_element(
    output,
    "      ",
    "Hsp_query-from",
    hsp.query_from.to_string(),
  )
  blast_xml_write_element(
    output,
    "      ",
    "Hsp_query-to",
    hsp.query_to.to_string(),
  )
  blast_xml_write_element(
    output,
    "      ",
    "Hsp_hit-from",
    hsp.target_from.to_string(),
  )
  blast_xml_write_element(
    output,
    "      ",
    "Hsp_hit-to",
    hsp.target_to.to_string(),
  )
  blast_xml_write_optional_int(
    output,
    "      ",
    "Hsp_pattern-from",
    hsp.pattern_from,
  )
  blast_xml_write_optional_int(
    output,
    "      ",
    "Hsp_pattern-to",
    hsp.pattern_to,
  )
  blast_xml_write_optional_int(
    output,
    "      ",
    "Hsp_query-frame",
    blast_xml_xml1_frame(hsp.query_frame, hsp.query_strand),
  )
  blast_xml_write_optional_int(
    output,
    "      ",
    "Hsp_hit-frame",
    blast_xml_xml1_frame(hsp.target_frame, hsp.target_strand),
  )
  blast_xml_write_element(
    output,
    "      ",
    "Hsp_identity",
    hsp.identity.to_string(),
  )
  blast_xml_write_optional_int(output, "      ", "Hsp_positive", hsp.positive)
  blast_xml_write_optional_int(output, "      ", "Hsp_gaps", hsp.gaps)
  blast_xml_write_optional_int(
    output,
    "      ",
    "Hsp_align-len",
    hsp.alignment_length,
  )
  blast_xml_write_optional_int(output, "      ", "Hsp_density", hsp.density)
  blast_xml_write_element(output, "      ", "Hsp_qseq", hsp.query_sequence)
  blast_xml_write_element(output, "      ", "Hsp_hseq", hsp.target_sequence)
  blast_xml_write_optional_string(output, "      ", "Hsp_midline", hsp.midline)
  output.write_string("    </Hsp>\n")
}

///|
fn blast_xml_write_hit_xml1(
  output : StringBuilder,
  hit : BlastXmlHit,
) -> Unit raise BlastXmlError {
  output.write_string("<Hit>\n")
  blast_xml_write_element(output, "  ", "Hit_num", hit.number.to_string())
  let primary = hit.descriptions[0]
  blast_xml_write_element(output, "  ", "Hit_id", primary.id)
  blast_xml_write_element(
    output,
    "  ",
    "Hit_def",
    blast_xml_legacy_description(hit),
  )
  blast_xml_write_element(output, "  ", "Hit_accession", primary.accession)
  blast_xml_write_element(output, "  ", "Hit_len", hit.length.to_string())
  output.write_string("  <Hit_hsps>\n")
  for hsp in hit.hsps {
    blast_xml_write_hsp_xml1(output, hsp)
  }
  output.write_string("  </Hit_hsps>\n")
  output.write_string("</Hit>\n")
}

///|
fn blast_xml_write_query_xml1(
  output : StringBuilder,
  query : BlastXmlQuery,
  prefix : String,
) -> Unit raise BlastXmlError {
  blast_xml_write_element(output, "", prefix + "query-ID", query.id)
  blast_xml_write_element(output, "", prefix + "query-def", query.description)
  blast_xml_write_element(
    output,
    "",
    prefix + "query-len",
    query.length.to_string(),
  )
}

///|
fn blast_xml_write_record_xml1(
  output : StringBuilder,
  record : BlastXmlRecord,
) -> Unit raise BlastXmlError {
  output.write_string("<Iteration>\n")
  blast_xml_write_element(
    output,
    "  ",
    "Iteration_iter-num",
    record.number.to_string(),
  )
  blast_xml_write_query_xml1(output, record.query, "Iteration_")
  output.write_string("<Iteration_hits>\n")
  for hit in record.hits {
    blast_xml_write_hit_xml1(output, hit)
  }
  output.write_string("</Iteration_hits>\n")
  match record.statistics {
    Some(statistics) => {
      output.write_string("  <Iteration_stat>\n")
      blast_xml_write_statistics(output, statistics, false, "    ")
      output.write_string("  </Iteration_stat>\n")
    }
    None => ()
  }
  blast_xml_write_optional_string(
    output,
    "  ",
    "Iteration_message",
    record.message,
  )
  output.write_string("</Iteration>\n")
}

///|
fn blast_xml_write_xml1(
  document : BlastXmlDocument,
) -> String raise BlastXmlError {
  let output = StringBuilder::new()
  output.write_string("<?xml version=\"1.0\"?>\n")
  output.write_string(
    "<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd\">\n",
  )
  output.write_string("<BlastOutput>\n")
  blast_xml_write_element(output, "  ", "BlastOutput_program", document.program)
  blast_xml_write_element(output, "  ", "BlastOutput_version", document.version)
  blast_xml_write_element(
    output,
    "  ",
    "BlastOutput_reference",
    document.reference,
  )
  blast_xml_write_element(output, "  ", "BlastOutput_db", document.database)
  let query = match document.global_query {
    Some(value) => value
    None => document.records[0].query
  }
  blast_xml_write_element(output, "  ", "BlastOutput_query-ID", query.id)
  blast_xml_write_element(
    output,
    "  ",
    "BlastOutput_query-def",
    query.description,
  )
  blast_xml_write_element(
    output,
    "  ",
    "BlastOutput_query-len",
    query.length.to_string(),
  )
  blast_xml_write_optional_string(
    output,
    "  ",
    "BlastOutput_query-seq",
    query.sequence,
  )
  output.write_string("  <BlastOutput_param>\n")
  output.write_string("    <Parameters>\n")
  blast_xml_write_parameters(output, document.parameters, false)
  output.write_string("    </Parameters>\n")
  output.write_string("  </BlastOutput_param>\n")
  output.write_string("<BlastOutput_iterations>\n")
  for record in document.records {
    blast_xml_write_record_xml1(output, record)
  }
  output.write_string("</BlastOutput_iterations>\n")
  match document.megablast_statistics {
    Some(statistics) => {
      output.write_string("  <BlastOutput_mbstat>\n")
      blast_xml_write_statistics(output, statistics, false, "    ")
      output.write_string("  </BlastOutput_mbstat>\n")
    }
    None => ()
  }
  output.write_string("</BlastOutput>\n")
  output.to_string()
}

///|
fn blast_xml_write_hsp_xml2(
  output : StringBuilder,
  hsp : BlastXmlHsp,
  program : String,
  extra_indent : String,
) -> Unit raise BlastXmlError {
  output.write_string(extra_indent + "                    <Hsp>\n")
  let indent = extra_indent + "                      "
  blast_xml_write_element(output, indent, "num", hsp.number.to_string())
  blast_xml_write_element(
    output,
    indent,
    "bit-score",
    hsp.bit_score.to_string(),
  )
  blast_xml_write_element(output, indent, "score", hsp.score.to_string())
  blast_xml_write_element(output, indent, "evalue", hsp.evalue.to_string())
  blast_xml_write_element(output, indent, "identity", hsp.identity.to_string())
  blast_xml_write_optional_int(output, indent, "positive", hsp.positive)
  blast_xml_write_optional_int(output, indent, "density", hsp.density)
  blast_xml_write_optional_int(output, indent, "pattern-from", hsp.pattern_from)
  blast_xml_write_optional_int(output, indent, "pattern-to", hsp.pattern_to)
  blast_xml_write_element(
    output,
    indent,
    "query-from",
    hsp.query_from.to_string(),
  )
  blast_xml_write_element(output, indent, "query-to", hsp.query_to.to_string())
  if program == "blastn" || program == "megablast" {
    blast_xml_write_optional_string(
      output,
      indent,
      "query-strand",
      blast_xml_xml2_strand(hsp.query_strand, hsp.query_frame),
    )
  } else if blast_xml_query_translated(program) {
    blast_xml_write_optional_int(output, indent, "query-frame", hsp.query_frame)
  }
  blast_xml_write_element(
    output,
    indent,
    "hit-from",
    hsp.target_from.to_string(),
  )
  blast_xml_write_element(output, indent, "hit-to", hsp.target_to.to_string())
  if program == "blastn" || program == "megablast" {
    blast_xml_write_optional_string(
      output,
      indent,
      "hit-strand",
      blast_xml_xml2_strand(hsp.target_strand, hsp.target_frame),
    )
  } else if blast_xml_target_translated(program) {
    blast_xml_write_optional_int(output, indent, "hit-frame", hsp.target_frame)
  }
  blast_xml_write_optional_int(
    output,
    indent,
    "align-len",
    hsp.alignment_length,
  )
  blast_xml_write_optional_int(output, indent, "gaps", hsp.gaps)
  blast_xml_write_element(output, indent, "qseq", hsp.query_sequence)
  blast_xml_write_element(output, indent, "hseq", hsp.target_sequence)
  blast_xml_write_optional_string(output, indent, "midline", hsp.midline)
  output.write_string(extra_indent + "                    </Hsp>\n")
}

///|
fn blast_xml_write_hit_xml2(
  output : StringBuilder,
  hit : BlastXmlHit,
  program : String,
  extra_indent : String,
) -> Unit raise BlastXmlError {
  output.write_string(extra_indent + "                <Hit>\n")
  blast_xml_write_element(
    output,
    extra_indent + "                  ",
    "num",
    hit.number.to_string(),
  )
  output.write_string(extra_indent + "                  <description>\n")
  for description in hit.descriptions {
    output.write_string(extra_indent + "                    <HitDescr>\n")
    blast_xml_write_element(
      output,
      extra_indent + "                      ",
      "id",
      description.id,
    )
    blast_xml_write_element(
      output,
      extra_indent + "                      ",
      "accession",
      description.accession,
    )
    blast_xml_write_element(
      output,
      extra_indent + "                      ",
      "title",
      description.title,
    )
    blast_xml_write_optional_int(
      output,
      extra_indent + "                      ",
      "taxid",
      description.taxid,
    )
    blast_xml_write_optional_string(
      output,
      extra_indent + "                      ",
      "sciname",
      description.scientific_name,
    )
    output.write_string(extra_indent + "                    </HitDescr>\n")
  }
  output.write_string(extra_indent + "                  </description>\n")
  blast_xml_write_element(
    output,
    extra_indent + "                  ",
    "len",
    hit.length.to_string(),
  )
  output.write_string(extra_indent + "                  <hsps>\n")
  for hsp in hit.hsps {
    blast_xml_write_hsp_xml2(output, hsp, program, extra_indent)
  }
  output.write_string(extra_indent + "                  </hsps>\n")
  output.write_string(extra_indent + "                </Hit>\n")
}

///|
fn blast_xml_write_record_xml2(
  output : StringBuilder,
  record : BlastXmlRecord,
  program : String,
) -> Unit raise BlastXmlError {
  let content_indent = if program == "psiblast" { "    " } else { "" }
  if program == "psiblast" {
    output.write_string("            <Iteration>\n")
    blast_xml_write_element(
      output,
      "              ",
      "iter-num",
      record.number.to_string(),
    )
  }
  output.write_string(content_indent + "          <search>\n")
  output.write_string(content_indent + "            <Search>\n")
  blast_xml_write_element(
    output,
    content_indent + "              ",
    "query-id",
    record.query.id,
  )
  blast_xml_write_element(
    output,
    content_indent + "              ",
    "query-title",
    record.query.description,
  )
  blast_xml_write_element(
    output,
    content_indent + "              ",
    "query-len",
    record.query.length.to_string(),
  )
  if record.query.masks.length() > 0 {
    output.write_string(content_indent + "              <query-masking>\n")
    for mask in record.query.masks {
      output.write_string(content_indent + "                <Range>\n")
      blast_xml_write_element(
        output,
        content_indent + "                  ",
        "from",
        (mask.start + 1).to_string(),
      )
      blast_xml_write_element(
        output,
        content_indent + "                  ",
        "to",
        mask.end.to_string(),
      )
      output.write_string(content_indent + "                </Range>\n")
    }
    output.write_string(content_indent + "              </query-masking>\n")
  }
  output.write_string(content_indent + "              <hits>\n")
  for hit in record.hits {
    blast_xml_write_hit_xml2(output, hit, program, content_indent)
  }
  output.write_string(content_indent + "              </hits>\n")
  match record.statistics {
    Some(statistics) => {
      output.write_string(content_indent + "              <stat>\n")
      blast_xml_write_statistics(
        output,
        statistics,
        true,
        content_indent + "                ",
      )
      output.write_string(content_indent + "              </stat>\n")
    }
    None => ()
  }
  blast_xml_write_optional_string(
    output,
    content_indent + "              ",
    "message",
    record.message,
  )
  output.write_string(content_indent + "            </Search>\n")
  output.write_string(content_indent + "          </search>\n")
  if program == "psiblast" {
    output.write_string("            </Iteration>\n")
  }
}

///|
fn blast_xml_write_output_xml2(
  output : StringBuilder,
  document : BlastXmlDocument,
  start : Int,
  end : Int,
) -> Unit raise BlastXmlError {
  output.write_string("<BlastOutput2>\n")
  output.write_string("  <report>\n")
  output.write_string("    <Report>\n")
  blast_xml_write_element(output, "      ", "program", document.program)
  blast_xml_write_element(output, "      ", "version", document.version)
  blast_xml_write_element(output, "      ", "reference", document.reference)
  output.write_string("      <search-target>\n")
  output.write_string("        <Target>\n")
  blast_xml_write_element(output, "          ", "db", document.database)
  output.write_string("        </Target>\n")
  output.write_string("      </search-target>\n")
  output.write_string("      <params>\n")
  output.write_string("        <Parameters>\n")
  blast_xml_write_parameters(output, document.parameters, true)
  output.write_string("        </Parameters>\n")
  output.write_string("      </params>\n")
  output.write_string("      <results>\n")
  output.write_string("        <Results>\n")
  if document.program == "psiblast" {
    output.write_string("          <iterations>\n")
  }
  for index = start; index < end; index = index + 1 {
    blast_xml_write_record_xml2(
      output,
      document.records[index],
      document.program,
    )
  }
  if document.program == "psiblast" {
    output.write_string("          </iterations>\n")
  }
  output.write_string("        </Results>\n")
  output.write_string("      </results>\n")
  output.write_string("    </Report>\n")
  output.write_string("  </report>\n")
  output.write_string("</BlastOutput2>\n")
}

///|
fn blast_xml_write_xml2(
  document : BlastXmlDocument,
) -> String raise BlastXmlError {
  let output = StringBuilder::new()
  output.write_string("<?xml version=\"1.0\"?>\n")
  output.write_string(
    "<BlastXML2 xmlns=\"http://www.ncbi.nlm.nih.gov\" xmlns:xs=\"http://www.w3.org/2001/XMLSchema-instance\" xs:schemaLocation=\"http://www.ncbi.nlm.nih.gov http://www.ncbi.nlm.nih.gov/data_specs/schema_alt/NCBI_BlastOutput2.xsd\">\n",
  )
  if document.program == "psiblast" {
    blast_xml_write_output_xml2(output, document, 0, document.records.length())
  } else {
    for index = 0; index < document.records.length(); index = index + 1 {
      blast_xml_write_output_xml2(output, document, index, index + 1)
    }
  }
  output.write_string("</BlastXML2>\n")
  output.to_string()
}

///|
/// Write a canonical BLAST XML document using an explicit XML dialect.
///
/// XML1 collapses multiple XML2 descriptions into the legacy `Hit_def` text.
/// Fields that the selected schema cannot represent are rejected.
pub fn blast_xml_write_as(
  document : BlastXmlDocument,
  format : BlastXmlFormat,
) -> String raise BlastXmlError {
  blast_xml_validate_write(document, format)
  let result = match format {
    BlastXmlFormat::Xml1 => blast_xml_write_xml1(document)
    BlastXmlFormat::Xml2 => blast_xml_write_xml2(document)
  }
  ignore(blast_xml_parse(result))
  result
}

///|
/// Write a canonical BLAST XML document in its original dialect.
pub fn blast_xml_write(
  document : BlastXmlDocument,
) -> String raise BlastXmlError {
  blast_xml_write_as(document, document.format)
}
