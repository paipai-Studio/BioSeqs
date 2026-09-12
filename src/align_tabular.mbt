// Parse alignment-aware BLAST outfmt 7 and FASTA -m 8CB/-m 8CC output.
//
// This follows the field and traceback semantics of Biopython 1.86
// Bio.Align.tabular. Unlike hit-table parsers, BTOP and aln_code values are
// expanded into explicit pairwise coordinate paths.

///|
/// Error raised by alignment-aware tabular parsing and conversion.
pub suberror AlignTabularError {
  AlignTabularError(String)
}

///|
/// Source used to construct an alignment path.
pub enum AlignTabularTrace {
  NoTrace
  Btop
  Cigar
} derive(Eq, Debug)

///|
/// One run in a tabular alignment path.
///
/// `Aligned` consumes both rows, `QueryGap` consumes only the target row, and
/// `TargetGap` consumes only the query row.
pub enum AlignTabularOperation {
  Aligned
  QueryGap
  TargetGap
} derive(Eq, Debug)

///|
/// One original field/value pair from a tabular row.
pub struct AlignTabularField {
  name : String
  value : String
} derive(Eq, Debug)

///|
/// Metadata attached to one BLAST or FASTA query block.
pub struct AlignTabularMetadata {
  program : String
  version : String
  command_line : String
  database : String
  rid : String
} derive(Eq, Debug)

///|
/// One pairwise alignment parsed from a tabular row.
///
/// Intervals are zero-based and half-open. Coordinate paths preserve the
/// reported strand direction. For translated searches, nucleotide axes use
/// nucleotide coordinates and have `units_per_residue == 3`.
pub struct AlignTabularAlignment {
  metadata : AlignTabularMetadata
  query_id : String
  query_description : String
  target_id : String
  query_length : Int?
  target_length : Int?
  query_start : Int?
  query_end : Int?
  target_start : Int?
  target_end : Int?
  query_strand : Strand
  target_strand : Strand
  query_sequence : String
  target_sequence : String
  target_coordinates : Array[Int]
  query_coordinates : Array[Int]
  operations : Array[AlignTabularOperation]
  trace : AlignTabularTrace
  trace_text : String
  path_is_absolute : Bool
  target_units_per_residue : Int
  query_units_per_residue : Int
  fields : Array[AlignTabularField]
} derive(Debug)

///|
/// One query block, including zero-hit queries.
pub struct AlignTabularQueryResult {
  metadata : AlignTabularMetadata
  query_id : String
  query_description : String
  query_length : Int?
  query_unit : String
  field_names : Array[String]
  declared_hits : Int
  alignments : Array[AlignTabularAlignment]
} derive(Debug)

///|
/// Complete parsed document.
pub struct AlignTabularDocument {
  queries : Array[AlignTabularQueryResult]
  processed_queries : Int?
} derive(Debug)

///|
priv struct AlignTabularRelativePath {
  target_coordinates : Array[Int]
  query_coordinates : Array[Int]
  operations : Array[AlignTabularOperation]
  columns : Int
}

///|
fn align_tabular_starts_with(value : String, prefix : String) -> Bool {
  if prefix.length() > value.length() {
    return false
  }
  for index in 0..<prefix.length() {
    if value.unsafe_get(index).to_int() != prefix.unsafe_get(index).to_int() {
      return false
    }
  }
  true
}

///|
fn align_tabular_ends_with(value : String, suffix : String) -> Bool {
  if suffix.length() > value.length() {
    return false
  }
  let offset = value.length() - suffix.length()
  for index in 0..<suffix.length() {
    if value.unsafe_get(offset + index).to_int() !=
      suffix.unsafe_get(index).to_int() {
      return false
    }
  }
  true
}

///|
fn align_tabular_trim(value : String) -> String {
  value.trim().to_owned()
}

///|
fn align_tabular_split_char(value : String, separator : Char) -> Array[String] {
  let result : Array[String] = []
  let mut start = 0
  for index in 0..<value.length() {
    if value.unsafe_get(index).to_int() == separator.to_int() {
      result.push(value[start:index].to_owned())
      start = index + 1
    }
  }
  result.push(value[start:value.length()].to_owned())
  result
}

///|
fn align_tabular_split_whitespace(value : String) -> Array[String] {
  let result : Array[String] = []
  let mut start = 0
  let mut in_token = false
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    let whitespace = code == ' '.to_int() ||
      code == '\t'.to_int() ||
      code == '\r'.to_int()
    if whitespace {
      if in_token {
        result.push(value[start:index].to_owned())
        in_token = false
      }
    } else if !in_token {
      start = index
      in_token = true
    }
  }
  if in_token {
    result.push(value[start:value.length()].to_owned())
  }
  result
}

///|
fn align_tabular_join(values : Array[String], start : Int) -> String {
  let mut result = ""
  for index in start..<values.length() {
    if index > start {
      result = result + " "
    }
    result = result + values[index]
  }
  result
}

///|
fn align_tabular_last_delimiter(value : String, delimiter : String) -> Int? {
  if delimiter.length() == 0 || delimiter.length() > value.length() {
    return None
  }
  let mut found : Int? = None
  for start in 0..<=(value.length() - delimiter.length()) {
    let mut matches = true
    for offset in 0..<delimiter.length() {
      if value.unsafe_get(start + offset).to_int() !=
        delimiter.unsafe_get(offset).to_int() {
        matches = false
        break
      }
    }
    if matches {
      found = Some(start)
    }
  }
  found
}

///|
fn align_tabular_parse_int(
  value : String,
  label : String,
) -> Int raise AlignTabularError {
  if value.length() == 0 {
    raise AlignTabularError("missing integer field " + label)
  }
  let mut start = 0
  let negative = value.unsafe_get(0).to_int() == '-'.to_int()
  if negative || value.unsafe_get(0).to_int() == '+'.to_int() {
    start = 1
  }
  if start == value.length() {
    raise AlignTabularError("invalid integer in " + label + ": " + value)
  }
  for index in start..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise AlignTabularError("invalid integer in " + label + ": " + value)
    }
  }
  let mut first_digit = start
  while first_digit < value.length() &&
        value.unsafe_get(first_digit).to_int() == '0'.to_int() {
    first_digit = first_digit + 1
  }
  if first_digit == value.length() {
    return 0
  }
  let significant_length = value.length() - first_digit
  let limit = if negative { "2147483648" } else { "2147483647" }
  let mut exceeds_limit = significant_length > limit.length()
  let mut equals_limit = significant_length == limit.length()
  if significant_length == limit.length() {
    for offset in 0..<limit.length() {
      let digit = value.unsafe_get(first_digit + offset).to_int()
      let bound = limit.unsafe_get(offset).to_int()
      if digit > bound {
        exceeds_limit = true
        equals_limit = false
        break
      } else if digit < bound {
        equals_limit = false
        break
      }
    }
  }
  if exceeds_limit {
    raise AlignTabularError("integer overflow in " + label + ": " + value)
  }
  if negative && equals_limit {
    return -2147483647 - 1
  }
  let mut result = 0
  for index in first_digit..<value.length() {
    result = result * 10 + (value.unsafe_get(index).to_int() - '0'.to_int())
  }
  if negative {
    -result
  } else {
    result
  }
}

///|
fn align_tabular_parse_double(
  value : String,
  label : String,
) -> Double raise AlignTabularError {
  if value.length() == 0 {
    raise AlignTabularError("missing number field " + label)
  }
  let mut index = 0
  if value.unsafe_get(0).to_int() == '-'.to_int() ||
    value.unsafe_get(0).to_int() == '+'.to_int() {
    index = 1
  }
  let mut mantissa_digits = 0
  while index < value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      break
    }
    mantissa_digits = mantissa_digits + 1
    index = index + 1
  }
  if index < value.length() && value.unsafe_get(index).to_int() == '.'.to_int() {
    index = index + 1
    while index < value.length() {
      let code = value.unsafe_get(index).to_int()
      if code < '0'.to_int() || code > '9'.to_int() {
        break
      }
      mantissa_digits = mantissa_digits + 1
      index = index + 1
    }
  }
  if mantissa_digits == 0 {
    raise AlignTabularError("invalid number in " + label + ": " + value)
  }
  if index < value.length() &&
    (
      value.unsafe_get(index).to_int() == 'e'.to_int() ||
      value.unsafe_get(index).to_int() == 'E'.to_int()
    ) {
    index = index + 1
    let exponent_value_start = index
    if index < value.length() &&
      (
        value.unsafe_get(index).to_int() == '-'.to_int() ||
        value.unsafe_get(index).to_int() == '+'.to_int()
      ) {
      index = index + 1
    }
    let exponent_start = index
    while index < value.length() {
      let code = value.unsafe_get(index).to_int()
      if code < '0'.to_int() || code > '9'.to_int() {
        break
      }
      index = index + 1
    }
    if index == exponent_start {
      raise AlignTabularError("invalid number in " + label + ": " + value)
    }
    ignore(
      align_tabular_parse_int(
        value[exponent_value_start:index].to_owned(),
        label + " exponent",
      ),
    )
  }
  if index != value.length() {
    raise AlignTabularError("invalid number in " + label + ": " + value)
  }
  let parse_value = if value.unsafe_get(0).to_int() == '+'.to_int() {
    value[1:value.length()].to_owned()
  } else {
    value
  }
  match parse_double(parse_value) {
    Some(number) =>
      if number.is_nan() || number.abs() > 1.0e300 {
        raise AlignTabularError("non-finite number in " + label + ": " + value)
      } else {
        number
      }
    None => raise AlignTabularError("invalid number in " + label + ": " + value)
  }
}

///|
fn align_tabular_abs(value : Int) -> Int {
  if value < 0 {
    -value
  } else {
    value
  }
}

///|
fn align_tabular_copy_ints(values : Array[Int]) -> Array[Int] {
  let result : Array[Int] = []
  for value in values {
    result.push(value)
  }
  result
}

///|
fn align_tabular_copy_strings(values : Array[String]) -> Array[String] {
  let result : Array[String] = []
  for value in values {
    result.push(value)
  }
  result
}

///|
fn align_tabular_copy_operations(
  values : Array[AlignTabularOperation],
) -> Array[AlignTabularOperation] {
  let result : Array[AlignTabularOperation] = []
  for value in values {
    result.push(value)
  }
  result
}

///|
fn align_tabular_known_blast_program(program : String) -> Bool {
  match program {
    "BLASTN"
    | "BLASTP"
    | "BLASTX"
    | "TBLASTN"
    | "TBLASTX"
    | "DELTABLAST"
    | "PSIBLAST"
    | "RPSBLAST"
    | "RPSTBLASTN" => true
    _ => false
  }
}

///|
fn align_tabular_known_fasta_program(program : String) -> Bool {
  match program {
    "FASTA"
    | "SSEARCH"
    | "GGSEARCH"
    | "GLSEARCH"
    | "FASTX"
    | "FASTY"
    | "TFASTX"
    | "TFASTY" => true
    _ => false
  }
}

///|
fn align_tabular_program_header(value : String) -> (String, String)? {
  let parts = align_tabular_split_whitespace(value)
  if parts.length() < 2 {
    return None
  }
  let program = parts[0].to_upper()
  if !align_tabular_known_blast_program(program) &&
    !align_tabular_known_fasta_program(program) {
    return None
  }
  Some((program, align_tabular_join(parts, 1)))
}

///|
fn align_tabular_field_value(
  fields : Array[AlignTabularField],
  name : String,
) -> String? {
  for field in fields {
    if field.name == name {
      return Some(field.value)
    }
  }
  None
}

///|
fn align_tabular_optional_int(
  fields : Array[AlignTabularField],
  name : String,
) -> Int? raise AlignTabularError {
  match align_tabular_field_value(fields, name) {
    Some(value) => Some(align_tabular_parse_int(value, name))
    None => None
  }
}

///|
fn align_tabular_optional_double(
  fields : Array[AlignTabularField],
  name : String,
) -> Double? raise AlignTabularError {
  match align_tabular_field_value(fields, name) {
    Some(value) => Some(align_tabular_parse_double(value, name))
    None => None
  }
}

///|
fn align_tabular_is_known_field(name : String) -> Bool {
  match name {
    "query id"
    | "subject id"
    | "% identity"
    | "alignment length"
    | "mismatches"
    | "gap opens"
    | "q. start"
    | "q. end"
    | "s. start"
    | "s. end"
    | "evalue"
    | "bit score"
    | "BTOP"
    | "aln_code"
    | "query gi"
    | "query acc."
    | "query acc.ver"
    | "query length"
    | "subject ids"
    | "subject gi"
    | "subject gis"
    | "subject acc."
    | "subject acc.ver"
    | "subject accs."
    | "subject length"
    | "query seq"
    | "subject seq"
    | "score"
    | "identical"
    | "positives"
    | "gaps"
    | "% positives"
    | "% hsp coverage"
    | "query/sbjct frames"
    | "query frame"
    | "sbjct frame"
    | "subject tax ids"
    | "subject sci names"
    | "subject com names"
    | "subject blast names"
    | "subject super kingdoms"
    | "subject title"
    | "subject titles"
    | "subject strand"
    | "% subject coverage" => true
    _ => false
  }
}

///|
fn align_tabular_validate_field_names(
  names : Array[String],
) -> Unit raise AlignTabularError {
  if names.length() == 0 {
    raise AlignTabularError("Fields header must not be empty")
  }
  let seen : Map[String, Bool] = Map([], capacity=names.length())
  for name in names {
    if name.length() == 0 {
      raise AlignTabularError("Fields header contains an empty field")
    }
    if !align_tabular_is_known_field(name) {
      raise AlignTabularError("unexpected tabular field: " + name)
    }
    if seen.contains(name) {
      raise AlignTabularError("duplicate tabular field: " + name)
    }
    seen.set(name, true)
  }
}

///|
fn align_tabular_push_step(
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
  operations : Array[AlignTabularOperation],
  operation : AlignTabularOperation,
  target_count : Int,
  query_count : Int,
) -> Unit {
  let last = operations.length() - 1
  if last >= 0 && operations[last] == operation {
    let point = target_coordinates.length() - 1
    target_coordinates[point] = target_coordinates[point] + target_count
    query_coordinates[point] = query_coordinates[point] + query_count
  } else {
    let point = target_coordinates.length() - 1
    target_coordinates.push(target_coordinates[point] + target_count)
    query_coordinates.push(query_coordinates[point] + query_count)
    operations.push(operation)
  }
}

///|
/// Parse a BLAST traceback operations (BTOP) string into a relative path.
pub fn align_tabular_parse_btop(
  btop : String,
) -> (Array[Int], Array[Int], Array[AlignTabularOperation]) raise AlignTabularError {
  if btop.length() == 0 {
    raise AlignTabularError("BTOP value must not be empty")
  }
  let target_coordinates = [0]
  let query_coordinates = [0]
  let operations : Array[AlignTabularOperation] = []
  let mut index = 0
  while index < btop.length() {
    let code = btop.unsafe_get(index).to_int()
    if code >= '0'.to_int() && code <= '9'.to_int() {
      let start = index
      while index < btop.length() {
        let digit = btop.unsafe_get(index).to_int()
        if digit >= '0'.to_int() && digit <= '9'.to_int() {
          index = index + 1
        } else {
          break
        }
      }
      let count = align_tabular_parse_int(
        btop[start:index].to_owned(),
        "BTOP match length",
      )
      if count <= 0 {
        raise AlignTabularError("BTOP match lengths must be positive")
      }
      align_tabular_push_step(
        target_coordinates,
        query_coordinates,
        operations,
        AlignTabularOperation::Aligned,
        count,
        count,
      )
    } else {
      if index + 1 >= btop.length() {
        raise AlignTabularError("truncated BTOP residue pair")
      }
      let first = btop.unsafe_get(index).unsafe_to_char()
      let second = btop.unsafe_get(index + 1).unsafe_to_char()
      if first == '-' && second == '-' {
        raise AlignTabularError("BTOP residue pair cannot contain two gaps")
      }
      if first == '-' {
        align_tabular_push_step(
          target_coordinates,
          query_coordinates,
          operations,
          AlignTabularOperation::QueryGap,
          1,
          0,
        )
      } else if second == '-' {
        align_tabular_push_step(
          target_coordinates,
          query_coordinates,
          operations,
          AlignTabularOperation::TargetGap,
          0,
          1,
        )
      } else {
        let first_valid = first.is_ascii_alphabetic() || first == '*'
        let second_valid = second.is_ascii_alphabetic() || second == '*'
        if !first_valid || !second_valid {
          raise AlignTabularError("invalid BTOP residue pair")
        }
        align_tabular_push_step(
          target_coordinates,
          query_coordinates,
          operations,
          AlignTabularOperation::Aligned,
          1,
          1,
        )
      }
      index = index + 2
    }
  }
  (target_coordinates, query_coordinates, operations)
}

///|
/// Parse FASTA `aln_code` CIGAR syntax into a relative path.
///
/// FASTA uses `I` for residues present only in the subject/target and `D` for
/// residues present only in the query, matching Biopython's parser.
pub fn align_tabular_parse_cigar(
  cigar : String,
) -> (Array[Int], Array[Int], Array[AlignTabularOperation]) raise AlignTabularError {
  if cigar.length() == 0 {
    raise AlignTabularError("aln_code value must not be empty")
  }
  let target_coordinates = [0]
  let query_coordinates = [0]
  let operations : Array[AlignTabularOperation] = []
  let mut index = 0
  while index < cigar.length() {
    let start = index
    while index < cigar.length() {
      let code = cigar.unsafe_get(index).to_int()
      if code >= '0'.to_int() && code <= '9'.to_int() {
        index = index + 1
      } else {
        break
      }
    }
    if start == index || index >= cigar.length() {
      raise AlignTabularError("invalid aln_code operation")
    }
    let count = align_tabular_parse_int(
      cigar[start:index].to_owned(),
      "aln_code length",
    )
    if count <= 0 {
      raise AlignTabularError("aln_code lengths must be positive")
    }
    let operation = cigar.unsafe_get(index).unsafe_to_char()
    match operation {
      'M' =>
        align_tabular_push_step(
          target_coordinates,
          query_coordinates,
          operations,
          AlignTabularOperation::Aligned,
          count,
          count,
        )
      'I' =>
        align_tabular_push_step(
          target_coordinates,
          query_coordinates,
          operations,
          AlignTabularOperation::QueryGap,
          count,
          0,
        )
      'D' =>
        align_tabular_push_step(
          target_coordinates,
          query_coordinates,
          operations,
          AlignTabularOperation::TargetGap,
          0,
          count,
        )
      _ => raise AlignTabularError("unsupported aln_code operation")
    }
    index = index + 1
  }
  (target_coordinates, query_coordinates, operations)
}

///|
fn align_tabular_path(
  fields : Array[AlignTabularField],
) -> AlignTabularRelativePath? raise AlignTabularError {
  let btop = align_tabular_field_value(fields, "BTOP")
  let cigar = align_tabular_field_value(fields, "aln_code")
  if btop is Some(_) && cigar is Some(_) {
    raise AlignTabularError("a row cannot contain both BTOP and aln_code")
  }
  let parsed = match (btop, cigar) {
    (Some(value), None) => {
      let (target, query, operations) = align_tabular_parse_btop(value)
      Some((target, query, operations))
    }
    (None, Some(value)) => {
      let (target, query, operations) = align_tabular_parse_cigar(value)
      Some((target, query, operations))
    }
    _ => None
  }
  match parsed {
    Some((target, query, operations)) => {
      let mut columns = 0
      for index in 0..<operations.length() {
        let target_count = target[index + 1] - target[index]
        let query_count = query[index + 1] - query[index]
        columns = columns +
          (if target_count > query_count { target_count } else { query_count })
      }
      Some(AlignTabularRelativePath::{
        target_coordinates: target,
        query_coordinates: query,
        operations,
        columns,
      })
    }
    None => None
  }
}

///|
fn align_tabular_axis_factor(program : String, query : Bool) -> Int {
  if query {
    match program {
      "BLASTX" | "TBLASTX" | "FASTX" | "FASTY" => 3
      _ => 1
    }
  } else {
    match program {
      "TBLASTN" | "TBLASTX" | "RPSTBLASTN" | "TFASTX" | "TFASTY" => 3
      _ => 1
    }
  }
}

///|
fn align_tabular_interval(
  reported_start : Int?,
  reported_end : Int?,
  label : String,
) -> (Int?, Int?, Strand) raise AlignTabularError {
  match (reported_start, reported_end) {
    (None, None) => (None, None, Strand::Star)
    (Some(first), Some(last)) => {
      if first <= 0 || last <= 0 {
        raise AlignTabularError(label + " coordinates must be positive")
      }
      if first <= last {
        (Some(first - 1), Some(last), Strand::Plus)
      } else {
        (Some(last - 1), Some(first), Strand::Minus)
      }
    }
    _ => raise AlignTabularError(label + " start and end must appear together")
  }
}

///|
fn align_tabular_absolute_axis(
  relative : Array[Int],
  reported_start : Int,
  reported_end : Int,
  factor : Int,
  label : String,
) -> Array[Int] raise AlignTabularError {
  let consumed = relative[relative.length() - 1] - relative[0]
  let expected = align_tabular_abs(reported_end - reported_start) + 1
  if consumed * factor != expected {
    raise AlignTabularError(
      label +
      " traceback span " +
      (consumed * factor).to_string() +
      " does not match reported span " +
      expected.to_string(),
    )
  }
  let direction = if reported_start <= reported_end { 1 } else { -1 }
  let origin = if direction > 0 { reported_start - 1 } else { reported_start }
  let result : Array[Int] = []
  for coordinate in relative {
    result.push(origin + direction * coordinate * factor)
  }
  result
}

///|
fn align_tabular_scale_axis(relative : Array[Int], factor : Int) -> Array[Int] {
  let result : Array[Int] = []
  for coordinate in relative {
    result.push(coordinate * factor)
  }
  result
}

///|
fn align_tabular_ungapped_length(sequence : String) -> Int {
  let mut length = 0
  for index in 0..<sequence.length() {
    if sequence.unsafe_get(index).to_int() != '-'.to_int() {
      length = length + 1
    }
  }
  length
}

///|
fn align_tabular_validate_sequence(
  sequence : String,
  label : String,
) -> Unit raise AlignTabularError {
  for index in 0..<sequence.length() {
    let code = sequence.unsafe_get(index).to_int()
    if code == ' '.to_int() ||
      code == '\t'.to_int() ||
      code == '\r'.to_int() ||
      code == '\n'.to_int() {
      raise AlignTabularError(label + " contains whitespace")
    }
  }
}

///|
fn align_tabular_parse_query_header(
  value : String,
  fasta_style : Bool,
) -> (String, String, Int?, String) raise AlignTabularError {
  let mut query_text = align_tabular_trim(value)
  let mut length : Int? = None
  let mut unit = ""
  if fasta_style {
    match align_tabular_last_delimiter(query_text, " - ") {
      Some(delimiter) => {
        let suffix = align_tabular_split_whitespace(
          query_text[delimiter + 3:query_text.length()].to_owned(),
        )
        if suffix.length() == 2 && (suffix[1] == "nt" || suffix[1] == "aa") {
          let parsed_length = align_tabular_parse_int(suffix[0], "query size")
          if parsed_length <= 0 {
            raise AlignTabularError("query size must be positive")
          }
          length = Some(parsed_length)
          unit = suffix[1]
          query_text = align_tabular_trim(query_text[0:delimiter].to_owned())
        }
      }
      None => ()
    }
  }
  let parts = align_tabular_split_whitespace(query_text)
  if parts.length() == 0 {
    raise AlignTabularError("Query header must contain an identifier")
  }
  let identifier = parts[0]
  let description = align_tabular_join(parts, 1)
  (identifier, description, length, unit)
}

///|
fn align_tabular_fields_header(value : String) -> Array[String] {
  let raw = align_tabular_split_char(value, ',')
  let result : Array[String] = []
  for field in raw {
    result.push(align_tabular_trim(field))
  }
  result
}

///|
fn align_tabular_row_columns(line : String, expected : Int) -> Array[String] {
  let tabular = align_tabular_split_char(line, '\t')
  if tabular.length() == expected {
    let result : Array[String] = []
    for value in tabular {
      result.push(align_tabular_trim(value))
    }
    return result
  }
  align_tabular_split_whitespace(line)
}

///|
fn align_tabular_fields(
  names : Array[String],
  columns : Array[String],
) -> Array[AlignTabularField] raise AlignTabularError {
  if columns.length() != names.length() {
    raise AlignTabularError(
      "tabular row has " +
      columns.length().to_string() +
      " columns but Fields declares " +
      names.length().to_string(),
    )
  }
  let result : Array[AlignTabularField] = []
  for index in 0..<names.length() {
    result.push(AlignTabularField::{
      name: names[index],
      value: columns[index],
    })
  }
  result
}

///|
fn align_tabular_validate_numeric_fields(
  fields : Array[AlignTabularField],
) -> Unit raise AlignTabularError {
  for field in fields {
    match field.name {
      "alignment length"
      | "mismatches"
      | "gap opens"
      | "q. start"
      | "q. end"
      | "s. start"
      | "s. end"
      | "query length"
      | "subject length"
      | "score"
      | "identical"
      | "positives"
      | "gaps" => {
        let value = align_tabular_parse_int(field.value, field.name)
        if value < 0 {
          raise AlignTabularError(field.name + " must be non-negative")
        }
        if (
            field.name == "alignment length" ||
            field.name == "q. start" ||
            field.name == "q. end" ||
            field.name == "s. start" ||
            field.name == "s. end" ||
            field.name == "query length" ||
            field.name == "subject length"
          ) &&
          value == 0 {
          raise AlignTabularError(field.name + " must be positive")
        }
      }
      "% identity"
      | "evalue"
      | "bit score"
      | "% positives"
      | "% hsp coverage"
      | "% subject coverage" => {
        let value = align_tabular_parse_double(field.value, field.name)
        if value < 0.0 {
          raise AlignTabularError(field.name + " must be non-negative")
        }
        if (
            field.name == "% identity" ||
            field.name == "% positives" ||
            field.name == "% hsp coverage" ||
            field.name == "% subject coverage"
          ) &&
          value > 100.0 {
          raise AlignTabularError(field.name + " must not exceed 100")
        }
      }
      _ => ()
    }
  }
}

///|
fn align_tabular_alignment_from_row(
  metadata : AlignTabularMetadata,
  header_query_id : String,
  header_description : String,
  header_query_length : Int?,
  names : Array[String],
  line : String,
) -> AlignTabularAlignment raise AlignTabularError {
  let columns = align_tabular_row_columns(line, names.length())
  let fields = align_tabular_fields(names, columns)
  align_tabular_validate_numeric_fields(fields)
  let query_id = match align_tabular_field_value(fields, "query id") {
    Some(value) => value
    None =>
      match align_tabular_field_value(fields, "query acc.ver") {
        Some(value) => value
        None => header_query_id
      }
  }
  if query_id.length() == 0 {
    raise AlignTabularError("tabular row has no query identifier")
  }
  if header_query_id.length() > 0 && query_id != header_query_id {
    raise AlignTabularError(
      "row query identifier " +
      query_id +
      " does not match Query header " +
      header_query_id,
    )
  }
  let target_id = match align_tabular_field_value(fields, "subject id") {
    Some(value) => value
    None =>
      match align_tabular_field_value(fields, "subject acc.ver") {
        Some(value) => value
        None => ""
      }
  }
  if target_id.length() == 0 {
    raise AlignTabularError("tabular row has no subject identifier")
  }
  let row_query_length = align_tabular_optional_int(fields, "query length")
  let query_length = match (header_query_length, row_query_length) {
    (Some(header), Some(row)) =>
      if header != row {
        raise AlignTabularError("query length disagrees with the Query header")
      } else {
        Some(header)
      }
    (Some(header), None) => Some(header)
    (None, Some(row)) => Some(row)
    (None, None) => None
  }
  match query_length {
    Some(value) =>
      if value <= 0 {
        raise AlignTabularError("query length must be positive")
      }
    None => ()
  }
  let target_length = align_tabular_optional_int(fields, "subject length")
  match target_length {
    Some(value) =>
      if value <= 0 {
        raise AlignTabularError("subject length must be positive")
      }
    None => ()
  }
  let query_reported_start = align_tabular_optional_int(fields, "q. start")
  let query_reported_end = align_tabular_optional_int(fields, "q. end")
  let target_reported_start = align_tabular_optional_int(fields, "s. start")
  let target_reported_end = align_tabular_optional_int(fields, "s. end")
  let (query_start, query_end, query_strand) = align_tabular_interval(
    query_reported_start, query_reported_end, "query",
  )
  let (target_start, target_end, target_strand) = align_tabular_interval(
    target_reported_start, target_reported_end, "subject",
  )
  match (query_length, query_end) {
    (Some(length), Some(end)) =>
      if end > length {
        raise AlignTabularError("query coordinates exceed query length")
      }
    _ => ()
  }
  match (target_length, target_end) {
    (Some(length), Some(end)) =>
      if end > length {
        raise AlignTabularError("subject coordinates exceed subject length")
      }
    _ => ()
  }
  let query_sequence = align_tabular_field_value(fields, "query seq").unwrap_or(
    "",
  )
  let target_sequence = align_tabular_field_value(fields, "subject seq").unwrap_or(
    "",
  )
  align_tabular_validate_sequence(query_sequence, "query sequence")
  align_tabular_validate_sequence(target_sequence, "subject sequence")
  if query_sequence.length() > 0 &&
    target_sequence.length() > 0 &&
    query_sequence.length() != target_sequence.length() {
    raise AlignTabularError(
      "query and subject aligned sequences have different lengths",
    )
  }
  let path = align_tabular_path(fields)
  let query_factor = align_tabular_axis_factor(metadata.program, true)
  let target_factor = align_tabular_axis_factor(metadata.program, false)
  let mut target_coordinates : Array[Int] = []
  let mut query_coordinates : Array[Int] = []
  let mut operations : Array[AlignTabularOperation] = []
  let mut path_is_absolute = false
  let trace = if align_tabular_field_value(fields, "BTOP") is Some(_) {
    AlignTabularTrace::Btop
  } else if align_tabular_field_value(fields, "aln_code") is Some(_) {
    AlignTabularTrace::Cigar
  } else {
    AlignTabularTrace::NoTrace
  }
  let trace_text = match trace {
    Btop => align_tabular_field_value(fields, "BTOP").unwrap()
    Cigar => align_tabular_field_value(fields, "aln_code").unwrap()
    NoTrace => ""
  }
  match path {
    Some(relative) => {
      let alignment_length = align_tabular_optional_int(
        fields, "alignment length",
      )
      match alignment_length {
        Some(length) =>
          if length != relative.columns {
            raise AlignTabularError(
              "traceback columns do not match alignment length",
            )
          }
        None => ()
      }
      if query_sequence.length() > 0 {
        if query_sequence.length() != relative.columns {
          raise AlignTabularError(
            "query sequence length does not match traceback columns",
          )
        }
        let consumed = relative.query_coordinates[relative.query_coordinates.length() -
          1]
        if align_tabular_ungapped_length(query_sequence) != consumed {
          raise AlignTabularError(
            "query sequence consumption does not match traceback",
          )
        }
      }
      if target_sequence.length() > 0 {
        if target_sequence.length() != relative.columns {
          raise AlignTabularError(
            "subject sequence length does not match traceback columns",
          )
        }
        let consumed = relative.target_coordinates[relative.target_coordinates.length() -
          1]
        if align_tabular_ungapped_length(target_sequence) != consumed {
          raise AlignTabularError(
            "subject sequence consumption does not match traceback",
          )
        }
      }
      match
        (
          query_reported_start, query_reported_end, target_reported_start, target_reported_end,
        ) {
        (Some(qs), Some(qe), Some(ts), Some(te)) => {
          query_coordinates = align_tabular_absolute_axis(
            relative.query_coordinates,
            qs,
            qe,
            query_factor,
            "query",
          )
          target_coordinates = align_tabular_absolute_axis(
            relative.target_coordinates,
            ts,
            te,
            target_factor,
            "subject",
          )
          path_is_absolute = true
        }
        (None, None, None, None) => {
          query_coordinates = align_tabular_scale_axis(
            relative.query_coordinates,
            query_factor,
          )
          target_coordinates = align_tabular_scale_axis(
            relative.target_coordinates,
            target_factor,
          )
        }
        _ =>
          raise AlignTabularError(
            "traceback requires all four query and subject coordinates",
          )
      }
      operations = align_tabular_copy_operations(relative.operations)
    }
    None => ()
  }
  AlignTabularAlignment::{
    metadata,
    query_id,
    query_description: header_description,
    target_id,
    query_length,
    target_length,
    query_start,
    query_end,
    target_start,
    target_end,
    query_strand,
    target_strand,
    query_sequence,
    target_sequence,
    target_coordinates,
    query_coordinates,
    operations,
    trace,
    trace_text,
    path_is_absolute,
    target_units_per_residue: target_factor,
    query_units_per_residue: query_factor,
    fields,
  }
}

///|
fn align_tabular_parse_processed(
  value : String,
) -> Int? raise AlignTabularError {
  let suffix = " queries"
  if !align_tabular_ends_with(value, suffix) {
    return None
  }
  let prefixes = ["BLAST processed ", "FASTA processed "]
  for prefix in prefixes {
    if align_tabular_starts_with(value, prefix) {
      let number = value[prefix.length():value.length() - suffix.length()].to_owned()
      return Some(align_tabular_parse_int(number, "processed query count"))
    }
  }
  None
}

///|
fn align_tabular_parse_hits(value : String) -> Int? raise AlignTabularError {
  let suffix = " hits found"
  if !align_tabular_ends_with(value, suffix) {
    return None
  }
  let number = align_tabular_trim(
    value[0:value.length() - suffix.length()].to_owned(),
  )
  Some(align_tabular_parse_int(number, "hit count"))
}

///|
/// Parse BLAST outfmt 7 or FASTA `-m 8CB`/`-m 8CC` text.
pub fn align_tabular_parse(
  content : String,
) -> AlignTabularDocument raise AlignTabularError {
  if align_tabular_trim(content).length() == 0 {
    raise AlignTabularError("empty alignment tabular document")
  }
  let lines = align_tabular_split_char(content, '\n')
  let queries : Array[AlignTabularQueryResult] = []
  let mut processed_queries : Int? = None
  let mut index = 0
  while index < lines.length() {
    let line = align_tabular_trim(lines[index])
    if line.length() == 0 {
      index = index + 1
      continue
    }
    if !align_tabular_starts_with(line, "# ") {
      raise AlignTabularError("missing alignment tabular header")
    }
    let first_value = align_tabular_trim(line[2:line.length()].to_owned())
    match align_tabular_parse_processed(first_value) {
      Some(count) => {
        if count < 0 {
          raise AlignTabularError("processed query count must be non-negative")
        }
        if processed_queries is Some(_) {
          raise AlignTabularError("duplicate processed query summary")
        }
        processed_queries = Some(count)
        index = index + 1
        continue
      }
      None => ()
    }
    let mut command_line = ""
    let mut program = ""
    let mut version = ""
    let mut fasta_style = false
    match align_tabular_program_header(first_value) {
      Some((parsed_program, parsed_version)) => {
        program = parsed_program
        version = parsed_version
        fasta_style = align_tabular_known_fasta_program(program)
        index = index + 1
      }
      None => {
        command_line = first_value
        fasta_style = true
        index = index + 1
        if index >= lines.length() {
          raise AlignTabularError("missing FASTA program header")
        }
        let program_line = align_tabular_trim(lines[index])
        if !align_tabular_starts_with(program_line, "# ") {
          raise AlignTabularError("missing FASTA program header")
        }
        let program_value = align_tabular_trim(
          program_line[2:program_line.length()].to_owned(),
        )
        match align_tabular_program_header(program_value) {
          Some((parsed_program, parsed_version)) =>
            if !align_tabular_known_fasta_program(parsed_program) {
              raise AlignTabularError(
                "command line must be followed by a FASTA-suite program header",
              )
            } else {
              program = parsed_program
              version = parsed_version
            }
          None => raise AlignTabularError("invalid FASTA program header")
        }
        index = index + 1
      }
    }
    let mut query_id = ""
    let mut query_description = ""
    let mut query_length : Int? = None
    let mut query_unit = ""
    let mut database = ""
    let mut rid = ""
    let mut field_names : Array[String] = []
    let mut declared_hits : Int? = None
    let mut has_query = false
    let mut has_database = false
    let mut has_rid = false
    let mut has_fields = false
    while index < lines.length() {
      let header = align_tabular_trim(lines[index])
      if !align_tabular_starts_with(header, "# ") {
        break
      }
      let value = align_tabular_trim(header[2:header.length()].to_owned())
      match align_tabular_parse_hits(value) {
        Some(count) => {
          if count < 0 {
            raise AlignTabularError("hit count must be non-negative")
          }
          declared_hits = Some(count)
          index = index + 1
          break
        }
        None => ()
      }
      if align_tabular_starts_with(value, "Query:") {
        if has_query {
          raise AlignTabularError("duplicate Query header")
        }
        has_query = true
        let parsed = align_tabular_parse_query_header(
          value[6:value.length()].to_owned(),
          fasta_style,
        )
        query_id = parsed.0
        query_description = parsed.1
        query_length = parsed.2
        query_unit = parsed.3
      } else if align_tabular_starts_with(value, "Database:") {
        if has_database {
          raise AlignTabularError("duplicate Database header")
        }
        has_database = true
        database = align_tabular_trim(value[9:value.length()].to_owned())
      } else if align_tabular_starts_with(value, "RID:") {
        if has_rid {
          raise AlignTabularError("duplicate RID header")
        }
        has_rid = true
        rid = align_tabular_trim(value[4:value.length()].to_owned())
      } else if align_tabular_starts_with(value, "Fields:") {
        if has_fields {
          raise AlignTabularError("duplicate Fields header")
        }
        has_fields = true
        field_names = align_tabular_fields_header(
          value[7:value.length()].to_owned(),
        )
        align_tabular_validate_field_names(field_names)
      } else {
        raise AlignTabularError("unexpected header line: " + value)
      }
      index = index + 1
    }
    let hit_count = match declared_hits {
      Some(value) => value
      None => raise AlignTabularError("query block is missing a hit count")
    }
    if query_id.length() == 0 {
      raise AlignTabularError("query block is missing a Query header")
    }
    if database.length() == 0 {
      raise AlignTabularError("query block is missing a Database header")
    }
    if hit_count > 0 && field_names.length() == 0 {
      raise AlignTabularError("non-empty query block is missing Fields")
    }
    let metadata = AlignTabularMetadata::{
      program,
      version,
      command_line,
      database,
      rid,
    }
    let alignments : Array[AlignTabularAlignment] = []
    while index < lines.length() {
      let row = align_tabular_trim(lines[index])
      if row.length() == 0 {
        index = index + 1
        continue
      }
      if align_tabular_starts_with(row, "# ") {
        break
      }
      alignments.push(
        align_tabular_alignment_from_row(
          metadata, query_id, query_description, query_length, field_names, row,
        ),
      )
      index = index + 1
    }
    if alignments.length() != hit_count {
      raise AlignTabularError(
        "declared " +
        hit_count.to_string() +
        " hits but parsed " +
        alignments.length().to_string() +
        " rows",
      )
    }
    let mut inferred_query_length : Int? = None
    for alignment in alignments {
      if inferred_query_length is None && alignment.query_length is Some(_) {
        inferred_query_length = alignment.query_length
      }
    }
    let final_query_length = if query_length is Some(_) {
      query_length
    } else {
      inferred_query_length
    }
    for alignment in alignments {
      if final_query_length is Some(_) &&
        alignment.query_length is Some(_) &&
        final_query_length != alignment.query_length {
        raise AlignTabularError("query length changes within a query block")
      }
    }
    queries.push(AlignTabularQueryResult::{
      metadata,
      query_id,
      query_description,
      query_length: final_query_length,
      query_unit,
      field_names: align_tabular_copy_strings(field_names),
      declared_hits: hit_count,
      alignments,
    })
  }
  match processed_queries {
    Some(count) =>
      if count != queries.length() {
        raise AlignTabularError(
          "processed query count does not match parsed query blocks",
        )
      }
    None => ()
  }
  AlignTabularDocument::{ queries, processed_queries, }
}

///|
/// Return an original row field.
pub fn AlignTabularAlignment::field(
  self : AlignTabularAlignment,
  name : String,
) -> String? {
  align_tabular_field_value(self.fields, name)
}

///|
/// Return an integer row field.
pub fn AlignTabularAlignment::integer_field(
  self : AlignTabularAlignment,
  name : String,
) -> Int? raise AlignTabularError {
  align_tabular_optional_int(self.fields, name)
}

///|
/// Return a floating-point row field.
pub fn AlignTabularAlignment::number_field(
  self : AlignTabularAlignment,
  name : String,
) -> Double? raise AlignTabularError {
  align_tabular_optional_double(self.fields, name)
}

///|
/// Number of printed alignment columns represented by the path.
pub fn AlignTabularAlignment::alignment_columns(
  self : AlignTabularAlignment,
) -> Int {
  let mut result = 0
  for index in 0..<self.operations.length() {
    let target = align_tabular_abs(
        self.target_coordinates[index + 1] - self.target_coordinates[index],
      ) /
      self.target_units_per_residue
    let query = align_tabular_abs(
        self.query_coordinates[index + 1] - self.query_coordinates[index],
      ) /
      self.query_units_per_residue
    result = result + (if target > query { target } else { query })
  }
  result
}

///|
/// Number of aligned residue pairs represented by the path.
pub fn AlignTabularAlignment::aligned_residues(
  self : AlignTabularAlignment,
) -> Int {
  let mut result = 0
  for index in 0..<self.operations.length() {
    if self.operations[index] == AlignTabularOperation::Aligned {
      result = result +
        align_tabular_abs(
          self.query_coordinates[index + 1] - self.query_coordinates[index],
        ) /
        self.query_units_per_residue
    }
  }
  result
}

///|
/// Number of query-gap residues, i.e. residues present only in the target.
pub fn AlignTabularAlignment::query_gap_residues(
  self : AlignTabularAlignment,
) -> Int {
  let mut result = 0
  for index in 0..<self.operations.length() {
    if self.operations[index] == AlignTabularOperation::QueryGap {
      result = result +
        align_tabular_abs(
          self.target_coordinates[index + 1] - self.target_coordinates[index],
        ) /
        self.target_units_per_residue
    }
  }
  result
}

///|
/// Number of target-gap residues, i.e. residues present only in the query.
pub fn AlignTabularAlignment::target_gap_residues(
  self : AlignTabularAlignment,
) -> Int {
  let mut result = 0
  for index in 0..<self.operations.length() {
    if self.operations[index] == AlignTabularOperation::TargetGap {
      result = result +
        align_tabular_abs(
          self.query_coordinates[index + 1] - self.query_coordinates[index],
        ) /
        self.query_units_per_residue
    }
  }
  result
}

///|
/// Number of distinct gap runs.
pub fn AlignTabularAlignment::gap_events(self : AlignTabularAlignment) -> Int {
  let mut result = 0
  for operation in self.operations {
    if operation != AlignTabularOperation::Aligned {
      result = result + 1
    }
  }
  result
}

///|
/// True if the query was reported on the reverse strand.
pub fn AlignTabularAlignment::query_is_reverse(
  self : AlignTabularAlignment,
) -> Bool {
  self.query_strand == Strand::Minus
}

///|
/// True if the subject was reported on the reverse strand.
pub fn AlignTabularAlignment::target_is_reverse(
  self : AlignTabularAlignment,
) -> Bool {
  self.target_strand == Strand::Minus
}

///|
/// Convert the path to the repository's coordinate alignment abstraction.
///
/// The target axis is normalized to increasing order; the query axis then
/// represents relative strand orientation.
pub fn AlignTabularAlignment::coordinate_alignment(
  self : AlignTabularAlignment,
) -> CoordinatePairwiseAlignment raise AlignTabularError {
  if self.target_coordinates.length() == 0 || !self.path_is_absolute {
    raise AlignTabularError(
      "an absolute traceback path is required for coordinate conversion",
    )
  }
  let mut target = align_tabular_copy_ints(self.target_coordinates)
  let mut query = align_tabular_copy_ints(self.query_coordinates)
  if target[0] > target[target.length() - 1] {
    let reversed_target : Array[Int] = []
    let reversed_query : Array[Int] = []
    let mut index = target.length()
    while index > 0 {
      index = index - 1
      reversed_target.push(target[index])
      reversed_query.push(query[index])
    }
    target = reversed_target
    query = reversed_query
  }
  let mut inferred_target_length = 0
  let mut inferred_query_length = 0
  for coordinate in target {
    if coordinate > inferred_target_length {
      inferred_target_length = coordinate
    }
  }
  for coordinate in query {
    if coordinate > inferred_query_length {
      inferred_query_length = coordinate
    }
  }
  let target_length = self.target_length.unwrap_or(inferred_target_length)
  let query_length = self.query_length.unwrap_or(inferred_query_length)
  coordinate_pairwise_alignment_with_lengths(
    self.target_id,
    target_length,
    self.query_id,
    query_length,
    target,
    query,
  ) catch {
    AlignmentMapError(message) =>
      raise AlignTabularError("coordinate conversion failed: " + message)
  }
}

///|
/// Compact human-readable alignment summary.
pub fn AlignTabularAlignment::summary(self : AlignTabularAlignment) -> String {
  let mut result = self.query_id + " -> " + self.target_id
  match self.field("evalue") {
    Some(value) => result = result + ", E=" + value
    None => ()
  }
  match self.field("bit score") {
    Some(value) => result = result + ", bits=" + value
    None => ()
  }
  if self.operations.length() > 0 {
    result = result +
      ", columns=" +
      self.alignment_columns().to_string() +
      ", gaps=" +
      (self.query_gap_residues() + self.target_gap_residues()).to_string()
  }
  result
}

///|
/// Flatten all query blocks into input-order alignments.
pub fn AlignTabularDocument::alignments(
  self : AlignTabularDocument,
) -> Array[AlignTabularAlignment] {
  let result : Array[AlignTabularAlignment] = []
  for query in self.queries {
    for alignment in query.alignments {
      result.push(alignment)
    }
  }
  result
}

///|
/// Find a query block by identifier.
pub fn AlignTabularDocument::find_query(
  self : AlignTabularDocument,
  query_id : String,
) -> AlignTabularQueryResult? {
  for query in self.queries {
    if query.query_id == query_id {
      return Some(query)
    }
  }
  None
}

///|
/// Return the lowest-E-value alignment in a query block.
pub fn AlignTabularQueryResult::best_by_evalue(
  self : AlignTabularQueryResult,
) -> AlignTabularAlignment? raise AlignTabularError {
  let mut best : AlignTabularAlignment? = None
  let mut best_value = 0.0
  for alignment in self.alignments {
    match alignment.number_field("evalue") {
      Some(value) =>
        if best is None || value < best_value {
          best = Some(alignment)
          best_value = value
        }
      None => ()
    }
  }
  best
}

///|
/// Return alignments with E-value at or below a threshold.
pub fn AlignTabularQueryResult::filter_evalue(
  self : AlignTabularQueryResult,
  threshold : Double,
) -> Array[AlignTabularAlignment] raise AlignTabularError {
  if threshold.is_nan() || threshold < 0.0 || threshold.abs() > 1.0e300 {
    raise AlignTabularError("E-value threshold must be finite and non-negative")
  }
  let result : Array[AlignTabularAlignment] = []
  for alignment in self.alignments {
    match alignment.number_field("evalue") {
      Some(value) => if value <= threshold { result.push(alignment) }
      None => ()
    }
  }
  result
}
