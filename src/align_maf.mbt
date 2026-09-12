// Biopython-compatible Bio.Align.maf support.
//
// The legacy maf.mbt API remains available for permissive block analysis.
// This module models modern MAF documents, strict a/s/i/e/q records,
// absolute coordinate paths, and in-memory MafIndex-style interval queries.

///|
pub suberror AlignMafError {
  AlignMafError(String)
}

///|
pub struct AlignMafTrack {
  name : String?
  description : String?
  frames : String?
  maf_dot : String?
  visibility : String?
  species_order : Array[String]
} derive(Eq, Debug)

///|
pub struct AlignMafDocument {
  version : String
  scoring : String?
  program : String?
  comments : Array[String]
  track : AlignMafTrack?
  blocks : Array[BigMafBlock]
} derive(Eq, Debug)

///|
pub struct AlignMafCoordinate {
  column : Int
  positions : Array[Int]
} derive(Eq, Debug)

///|
pub struct AlignMafIndexEntry {
  block_index : Int
  start : Int
  end : Int
} derive(Eq, Debug)

///|
pub struct AlignMafIndex {
  reference : String
  source_size : Int
  reference_strand : String
  entries : Array[AlignMafIndexEntry]
  blocks : Array[BigMafBlock]
} derive(Eq, Debug)

///|
pub struct AlignMafSplicedSequence {
  source : String
  text : String
} derive(Eq, Debug)

///|
pub struct AlignMafSplicedAlignment {
  reference : String
  strand : String
  starts : Array[Int]
  ends : Array[Int]
  sequences : Array[AlignMafSplicedSequence]
  columns : Int
} derive(Eq, Debug)

///|
pub struct AlignMafSummary {
  block_count : Int
  component_count : Int
  empty_component_count : Int
  aligned_columns : Int
  reference_bases : Int
  source_count : Int
} derive(Eq, Debug)

///|
fn align_maf_fail(message : String) -> Unit raise AlignMafError {
  raise AlignMafError(message)
}

///|
fn align_maf_is_space(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\r'.to_int() ||
  code == '\n'.to_int()
}

///|
fn align_maf_trim(value : String) -> String {
  let mut start = 0
  while start < value.length() &&
        align_maf_is_space(value.unsafe_get(start).to_int()) {
    start = start + 1
  }
  let mut end = value.length()
  while end > start && align_maf_is_space(value.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  value[start:end].to_owned()
}

///|
fn align_maf_strip_cr(value : String) -> String {
  if value.length() > 0 &&
    value.unsafe_get(value.length() - 1).to_int() == '\r'.to_int() {
    value[0:value.length() - 1].to_owned()
  } else {
    value
  }
}

///|
fn align_maf_lines(value : String) -> Array[String] {
  let raw = value.split("\n").to_array()
  let lines : Array[String] = []
  for line in raw {
    lines.push(align_maf_strip_cr(line.to_owned()))
  }
  lines
}

///|
fn align_maf_words(value : String) -> Array[String] {
  let words : Array[String] = []
  let mut index = 0
  while index < value.length() {
    while index < value.length() &&
          align_maf_is_space(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index >= value.length() {
      break
    }
    let start = index
    while index < value.length() &&
          !align_maf_is_space(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    words.push(value[start:index].to_owned())
  }
  words
}

///|
fn align_maf_shell_words(value : String) -> Array[String] raise AlignMafError {
  let words : Array[String] = []
  let mut current = ""
  let mut quote = 0
  let mut escaped = false
  let mut has_value = false
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if escaped {
      current = current + value[index:index + 1].to_owned()
      escaped = false
      has_value = true
    } else if code == '\\'.to_int() {
      escaped = true
      has_value = true
    } else if quote != 0 {
      if code == quote {
        quote = 0
      } else {
        current = current + value[index:index + 1].to_owned()
      }
      has_value = true
    } else if code == '"'.to_int() || code == '\''.to_int() {
      quote = code
      has_value = true
    } else if align_maf_is_space(code) {
      if has_value {
        words.push(current)
        current = ""
        has_value = false
      }
    } else {
      current = current + value[index:index + 1].to_owned()
      has_value = true
    }
  }
  if escaped || quote != 0 {
    align_maf_fail("MAF track line contains an unterminated quote or escape")
  }
  if has_value {
    words.push(current)
  }
  words
}

///|
fn align_maf_key_value(
  field : String,
  label : String,
) -> (String, String) raise AlignMafError {
  let mut equal = -1
  for index in 0..<field.length() {
    if field.unsafe_get(index).to_int() == '='.to_int() {
      equal = index
      break
    }
  }
  if equal <= 0 || equal + 1 >= field.length() {
    align_maf_fail("Malformed " + label + " field '" + field + "'")
  }
  (field[0:equal].to_owned(), field[equal + 1:field.length()].to_owned())
}

///|
fn align_maf_validate_scalar(
  value : String,
  label : String,
) -> Unit raise AlignMafError {
  if value.length() == 0 {
    align_maf_fail(label + " must not be empty")
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code == 0 || code == '\r'.to_int() || code == '\n'.to_int() {
      align_maf_fail(label + " contains a control character")
    }
  }
}

///|
fn align_maf_validate_token(
  value : String,
  label : String,
) -> Unit raise AlignMafError {
  align_maf_validate_scalar(value, label)
  for index in 0..<value.length() {
    if align_maf_is_space(value.unsafe_get(index).to_int()) {
      align_maf_fail(label + " must not contain whitespace")
    }
  }
}

///|
fn align_maf_parse_uint(
  value : String,
  label : String,
) -> Int raise AlignMafError {
  if value.length() == 0 {
    align_maf_fail(label + " must be an unsigned integer")
  }
  let mut result = 0
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      align_maf_fail(label + " must be an unsigned integer")
    }
    let digit = code - '0'.to_int()
    if result > (2147483647 - digit) / 10 {
      align_maf_fail(label + " exceeds the supported integer range")
    }
    result = result * 10 + digit
  }
  result
}

///|
fn align_maf_is_number(value : String) -> Bool {
  if value.length() == 0 {
    return false
  }
  let mut index = 0
  if value.unsafe_get(index).to_int() == '-'.to_int() ||
    value.unsafe_get(index).to_int() == '+'.to_int() {
    index = index + 1
  }
  let mut digits = 0
  while index < value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      break
    }
    digits = digits + 1
    index = index + 1
  }
  if index < value.length() && value.unsafe_get(index).to_int() == '.'.to_int() {
    index = index + 1
    while index < value.length() {
      let code = value.unsafe_get(index).to_int()
      if code < '0'.to_int() || code > '9'.to_int() {
        break
      }
      digits = digits + 1
      index = index + 1
    }
  }
  if digits == 0 {
    return false
  }
  if index < value.length() &&
    (
      value.unsafe_get(index).to_int() == 'e'.to_int() ||
      value.unsafe_get(index).to_int() == 'E'.to_int()
    ) {
    index = index + 1
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
      return false
    }
  }
  index == value.length()
}

///|
fn align_maf_parse_double(
  value : String,
  label : String,
) -> Double raise AlignMafError {
  if !align_maf_is_number(value) {
    align_maf_fail(label + " must be a number")
  }
  let parsed_value = if value.unsafe_get(0).to_int() == '+'.to_int() {
    value[1:value.length()].to_owned()
  } else {
    value
  }
  let number = match parse_double(parsed_value) {
    Some(number) => number
    None => {
      align_maf_fail(label + " must be a number")
      0.0
    }
  }
  if number.is_nan() || number.abs() > 1.0e300 {
    align_maf_fail(label + " must be finite")
  }
  number
}

///|
fn align_maf_copy_strings(values : Array[String]) -> Array[String] {
  let copied : Array[String] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn align_maf_copy_blocks(values : Array[BigMafBlock]) -> Array[BigMafBlock] {
  let copied : Array[BigMafBlock] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn align_maf_copy_ints(values : Array[Int]) -> Array[Int] {
  let copied : Array[Int] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
pub fn AlignMafTrack::create(
  name? : String? = None,
  description? : String? = None,
  frames? : String? = None,
  maf_dot? : String? = None,
  visibility? : String? = None,
  species_order? : Array[String] = [],
) -> AlignMafTrack raise AlignMafError {
  match name {
    Some(value) => align_maf_validate_scalar(value, "MAF track name")
    None => ()
  }
  match description {
    Some(value) => align_maf_validate_scalar(value, "MAF track description")
    None => ()
  }
  match frames {
    Some(value) => align_maf_validate_scalar(value, "MAF track frames")
    None => ()
  }
  match maf_dot {
    Some(value) =>
      if value != "on" && value != "off" {
        align_maf_fail("MAF track mafDot must be 'on' or 'off'")
      }
    None => ()
  }
  match visibility {
    Some(value) =>
      if value != "dense" && value != "pack" && value != "full" {
        align_maf_fail(
          "MAF track visibility must be 'dense', 'pack', or 'full'",
        )
      }
    None => ()
  }
  for source in species_order {
    align_maf_validate_token(source, "MAF track species")
  }
  AlignMafTrack::{
    name,
    description,
    frames,
    maf_dot,
    visibility,
    species_order: align_maf_copy_strings(species_order),
  }
}

///|
pub fn AlignMafDocument::create(
  blocks : Array[BigMafBlock],
  version? : String = "1",
  scoring? : String? = None,
  program? : String? = None,
  comments? : Array[String] = [],
  track? : AlignMafTrack? = None,
) -> AlignMafDocument raise AlignMafError {
  if version != "1" {
    align_maf_fail("MAF version must be 1")
  }
  match scoring {
    Some(value) => align_maf_validate_token(value, "MAF scoring metadata")
    None => ()
  }
  match program {
    Some(value) => align_maf_validate_token(value, "MAF program metadata")
    None => ()
  }
  for comment in comments {
    align_maf_validate_scalar(comment, "MAF comment")
  }
  AlignMafDocument::{
    version,
    scoring,
    program,
    comments: align_maf_copy_strings(comments),
    track,
    blocks: align_maf_copy_blocks(blocks),
  }
}

///|
fn align_maf_parse_track(line : String) -> AlignMafTrack raise AlignMafError {
  let words = align_maf_shell_words(line)
  if words.length() < 2 || words[0] != "track" {
    align_maf_fail("Malformed MAF track line")
  }
  let mut name : String? = None
  let mut description : String? = None
  let mut frames : String? = None
  let mut maf_dot : String? = None
  let mut visibility : String? = None
  let species_order : Array[String] = []
  for index in 1..<words.length() {
    let (key, value) = align_maf_key_value(words[index], "MAF track")
    if key == "name" {
      if name is Some(_) {
        align_maf_fail("Duplicate MAF track name")
      }
      name = Some(value)
    } else if key == "description" {
      if description is Some(_) {
        align_maf_fail("Duplicate MAF track description")
      }
      description = Some(value)
    } else if key == "frames" {
      if frames is Some(_) {
        align_maf_fail("Duplicate MAF track frames")
      }
      frames = Some(value)
    } else if key == "mafDot" {
      if maf_dot is Some(_) {
        align_maf_fail("Duplicate MAF track mafDot")
      }
      maf_dot = Some(value)
    } else if key == "visibility" {
      if visibility is Some(_) {
        align_maf_fail("Duplicate MAF track visibility")
      }
      visibility = Some(value)
    } else if key == "speciesOrder" {
      if species_order.length() > 0 {
        align_maf_fail("Duplicate MAF track speciesOrder")
      }
      for source in align_maf_words(value) {
        species_order.push(source)
      }
      if species_order.length() == 0 {
        align_maf_fail("MAF track speciesOrder must not be empty")
      }
    } else {
      align_maf_fail("Unexpected MAF track variable '" + key + "'")
    }
  }
  AlignMafTrack::create(
    name~,
    description~,
    frames~,
    maf_dot~,
    visibility~,
    species_order~,
  )
}

///|
fn align_maf_normalize_text(value : String) -> String {
  let output = StringBuilder::new(size_hint=value.length())
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code == '.'.to_int() || code == '='.to_int() || code == '_'.to_int() {
      output.write_char('-')
    } else {
      output.write_string(value[index:index + 1].to_owned())
    }
  }
  output.to_string()
}

///|
fn align_maf_parse_block(
  lines : Array[String],
) -> BigMafBlock raise AlignMafError {
  if lines.length() == 0 {
    align_maf_fail("MAF alignment block must not be empty")
  }
  let header = align_maf_words(align_maf_trim(lines[0]))
  if header.length() == 0 || header[0] != "a" {
    align_maf_fail("MAF alignment block must start with an 'a' line")
  }
  let mut score : Double? = None
  let mut pass_number : Int? = None
  for index in 1..<header.length() {
    let (key, value) = align_maf_key_value(header[index], "MAF alignment")
    if key == "score" {
      if score is Some(_) {
        align_maf_fail("Duplicate MAF score annotation")
      }
      score = Some(align_maf_parse_double(value, "MAF score"))
    } else if key == "pass" {
      if pass_number is Some(_) {
        align_maf_fail("Duplicate MAF pass annotation")
      }
      let pass = align_maf_parse_uint(value, "MAF pass")
      if pass <= 0 {
        align_maf_fail("MAF pass annotation must be positive")
      }
      pass_number = Some(pass)
    } else {
      align_maf_fail("Unknown MAF alignment annotation '" + key + "'")
    }
  }
  let components : Array[BigMafComponent] = []
  let empty_components : Array[BigMafEmptyComponent] = []
  let mut last_component = -1
  for line_index in 1..<lines.length() {
    let line = align_maf_trim(lines[line_index])
    if line.length() == 0 || line.unsafe_get(0).to_int() == '#'.to_int() {
      continue
    }
    let words = align_maf_words(line)
    if words.length() == 0 {
      continue
    }
    if words[0] == "s" {
      if words.length() != 7 {
        align_maf_fail("MAF sequence lines must contain seven fields")
      }
      let component_start = align_maf_parse_uint(
        words[2],
        "MAF component start",
      )
      let component_size = align_maf_parse_uint(words[3], "MAF component size")
      let source_size = align_maf_parse_uint(words[5], "MAF source size")
      if source_size <= 0 ||
        component_size <= 0 ||
        component_start > source_size ||
        component_size > source_size - component_start {
        align_maf_fail("MAF component coordinates exceed the source sequence")
      }
      let component = BigMafComponent::create(
        words[1],
        component_start,
        component_size,
        words[4],
        source_size,
        align_maf_normalize_text(words[6]),
      ) catch {
        BigMafError(message) => {
          align_maf_fail(message)
          abort("")
        }
      }
      components.push(component)
      last_component = components.length() - 1
    } else if words[0] == "i" {
      if last_component < 0 || words.length() != 6 {
        align_maf_fail(
          "MAF insertion lines must follow an s line and have six fields",
        )
      }
      if words[1] != components[last_component].source {
        align_maf_fail(
          "MAF insertion source does not match the preceding s line",
        )
      }
      if components[last_component].insertion is Some(_) {
        align_maf_fail("Duplicate MAF insertion line for one component")
      }
      let left_count = align_maf_parse_uint(
        words[3],
        "MAF left insertion count",
      )
      let right_count = align_maf_parse_uint(
        words[5],
        "MAF right insertion count",
      )
      let insertion = BigMafInsertion::create(
        words[2],
        left_count,
        words[4],
        right_count,
      ) catch {
        BigMafError(message) => {
          align_maf_fail(message)
          abort("")
        }
      }
      components[last_component] = components[last_component].with_insertion(
        insertion,
      ) catch {
        BigMafError(message) => {
          align_maf_fail(message)
          abort("")
        }
      }
    } else if words[0] == "q" {
      if last_component < 0 || words.length() != 3 {
        align_maf_fail(
          "MAF quality lines must follow an s line and have three fields",
        )
      }
      if words[1] != components[last_component].source {
        align_maf_fail("MAF quality source does not match the preceding s line")
      }
      if components[last_component].quality is Some(_) {
        align_maf_fail("Duplicate MAF quality line for one component")
      }
      components[last_component] = components[last_component].with_quality(
        words[2],
      ) catch {
        BigMafError(message) => {
          align_maf_fail(message)
          abort("")
        }
      }
    } else if words[0] == "e" {
      if words.length() != 7 {
        align_maf_fail("MAF empty component lines must contain seven fields")
      }
      let empty_start = align_maf_parse_uint(
        words[2],
        "MAF empty component start",
      )
      let empty_size = align_maf_parse_uint(
        words[3],
        "MAF empty component size",
      )
      let empty_source_size = align_maf_parse_uint(
        words[5],
        "MAF empty source size",
      )
      if empty_source_size <= 0 ||
        empty_start > empty_source_size ||
        empty_size > empty_source_size - empty_start {
        align_maf_fail(
          "MAF empty component coordinates exceed the source sequence",
        )
      }
      let empty = BigMafEmptyComponent::create(
        words[1],
        empty_start,
        empty_size,
        words[4],
        empty_source_size,
        words[6],
      ) catch {
        BigMafError(message) => {
          align_maf_fail(message)
          abort("")
        }
      }
      empty_components.push(empty)
      last_component = -1
    } else {
      align_maf_fail("Unexpected MAF line type '" + words[0] + "'")
    }
  }
  BigMafBlock::create(components, score~, pass_number~, empty_components~) catch {
    BigMafError(message) => {
      align_maf_fail(message)
      abort("")
    }
  }
}

///|
pub fn align_maf_parse(text : String) -> AlignMafDocument raise AlignMafError {
  let lines = align_maf_lines(text)
  let mut index = 0
  while index < lines.length() && align_maf_trim(lines[index]).length() == 0 {
    index = index + 1
  }
  if index >= lines.length() {
    align_maf_fail("MAF input is empty")
  }
  let mut track : AlignMafTrack? = None
  let first = align_maf_trim(lines[index])
  if first == "track" || first.has_prefix("track ") {
    track = Some(align_maf_parse_track(first))
    index = index + 1
  }
  while index < lines.length() && align_maf_trim(lines[index]).length() == 0 {
    index = index + 1
  }
  if index >= lines.length() {
    align_maf_fail("MAF input is missing its header")
  }
  let header = align_maf_words(align_maf_trim(lines[index]))
  if header.length() == 0 || header[0] != "##maf" {
    align_maf_fail("MAF header line must start with ##maf")
  }
  let mut version : String? = None
  let mut scoring : String? = None
  let mut program : String? = None
  for field_index in 1..<header.length() {
    let (key, value) = align_maf_key_value(header[field_index], "MAF header")
    if key == "version" {
      if version is Some(_) {
        align_maf_fail("Duplicate MAF version metadata")
      }
      version = Some(value)
    } else if key == "scoring" {
      if scoring is Some(_) {
        align_maf_fail("Duplicate MAF scoring metadata")
      }
      scoring = Some(value)
    } else if key == "program" {
      if program is Some(_) {
        align_maf_fail("Duplicate MAF program metadata")
      }
      program = Some(value)
    } else {
      align_maf_fail("Unexpected MAF header variable '" + key + "'")
    }
  }
  let found_version = match version {
    Some(value) => value
    None => {
      align_maf_fail("MAF header must declare version=1")
      ""
    }
  }
  if found_version != "1" {
    align_maf_fail("MAF version must be 1")
  }
  index = index + 1
  let comments : Array[String] = []
  let blocks : Array[BigMafBlock] = []
  let mut current : Array[String] = []
  let mut saw_block = false
  while index < lines.length() {
    let line = align_maf_trim(lines[index])
    if line.length() == 0 {
      index = index + 1
      continue
    }
    if line.unsafe_get(0).to_int() == '#'.to_int() &&
      current.length() == 0 &&
      !saw_block {
      comments.push(align_maf_trim(line[1:line.length()].to_owned()))
    } else if line == "a" || line.has_prefix("a ") || line.has_prefix("a\t") {
      if current.length() > 0 {
        blocks.push(align_maf_parse_block(current))
      }
      current = [line]
      saw_block = true
    } else {
      if current.length() == 0 {
        align_maf_fail("MAF data line appears before an alignment line")
      }
      current.push(line)
    }
    index = index + 1
  }
  if current.length() > 0 {
    blocks.push(align_maf_parse_block(current))
  }
  AlignMafDocument::create(
    blocks,
    version=found_version,
    scoring~,
    program~,
    comments~,
    track~,
  )
}

///|
fn align_maf_quote_track(value : String) -> String {
  let mut needs_quotes = false
  let output = StringBuilder::new(size_hint=value.length())
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if align_maf_is_space(code) {
      needs_quotes = true
    }
    if code == '"'.to_int() || code == '\\'.to_int() {
      output.write_char('\\')
      needs_quotes = true
    }
    output.write_string(value[index:index + 1].to_owned())
  }
  let escaped = output.to_string()
  if needs_quotes {
    "\"" + escaped + "\""
  } else {
    escaped
  }
}

///|
fn align_maf_write_track(output : StringBuilder, track : AlignMafTrack) -> Unit {
  output.write_string("track")
  match track.name {
    Some(value) => output.write_string(" name=" + align_maf_quote_track(value))
    None => ()
  }
  match track.description {
    Some(value) =>
      output.write_string(" description=" + align_maf_quote_track(value))
    None => ()
  }
  match track.frames {
    Some(value) =>
      output.write_string(" frames=" + align_maf_quote_track(value))
    None => ()
  }
  match track.maf_dot {
    Some(value) => output.write_string(" mafDot=" + value)
    None => ()
  }
  match track.visibility {
    Some(value) => output.write_string(" visibility=" + value)
    None => ()
  }
  if track.species_order.length() > 0 {
    let species = StringBuilder::new()
    for index in 0..<track.species_order.length() {
      if index > 0 {
        species.write_char(' ')
      }
      species.write_string(track.species_order[index])
    }
    output.write_string(
      " speciesOrder=" + align_maf_quote_track(species.to_string()),
    )
  }
  output.write_char('\n')
}

///|
pub fn align_maf_write(
  document : AlignMafDocument,
) -> String raise AlignMafError {
  let validated = AlignMafDocument::create(
    document.blocks,
    version=document.version,
    scoring=document.scoring,
    program=document.program,
    comments=document.comments,
    track=document.track,
  )
  let output = StringBuilder::new()
  match validated.track {
    Some(track) => align_maf_write_track(output, track)
    None => ()
  }
  output.write_string("##maf version=1")
  match validated.scoring {
    Some(value) => output.write_string(" scoring=" + value)
    None => ()
  }
  match validated.program {
    Some(value) => output.write_string(" program=" + value)
    None => ()
  }
  output.write_char('\n')
  for comment in validated.comments {
    output.write_string("# " + comment + "\n")
  }
  output.write_char('\n')
  for block in validated.blocks {
    output.write_string(block.to_maf())
  }
  output.to_string()
}

///|
fn align_maf_pattern_equal(left : Array[Bool], right : Array[Bool]) -> Bool {
  if left.length() != right.length() {
    return false
  }
  for index in 0..<left.length() {
    if left[index] != right[index] {
      return false
    }
  }
  true
}

///|
fn align_maf_column_pattern(block : BigMafBlock, column : Int) -> Array[Bool] {
  let pattern : Array[Bool] = []
  for component in block.components {
    pattern.push(component.text.unsafe_get(column).to_int() != '-'.to_int())
  }
  pattern
}

///|
fn align_maf_copy_positions(positions : Array[Int]) -> Array[Int] {
  let copied : Array[Int] = []
  for position in positions {
    copied.push(position)
  }
  copied
}

///|
pub fn align_maf_coordinate_path(
  block : BigMafBlock,
) -> Array[AlignMafCoordinate] {
  let width = block.aligned_columns()
  let positions : Array[Int] = []
  for component in block.components {
    positions.push(
      if component.strand == "+" {
        component.start
      } else {
        component.source_size - component.start
      },
    )
  }
  let path : Array[AlignMafCoordinate] = [
    AlignMafCoordinate::{
      column: 0,
      positions: align_maf_copy_positions(positions),
    },
  ]
  if width == 0 {
    return path
  }
  let mut previous = align_maf_column_pattern(block, 0)
  for column in 0..<width {
    let pattern = align_maf_column_pattern(block, column)
    if column > 0 && !align_maf_pattern_equal(previous, pattern) {
      path.push(AlignMafCoordinate::{
        column,
        positions: align_maf_copy_positions(positions),
      })
    }
    for row in 0..<block.components.length() {
      if pattern[row] {
        positions[row] = positions[row] +
          (if block.components[row].strand == "+" { 1 } else { -1 })
      }
    }
    previous = pattern
  }
  path.push(AlignMafCoordinate::{
    column: width,
    positions: align_maf_copy_positions(positions),
  })
  path
}

///|
pub fn align_maf_map_position(
  block : BigMafBlock,
  from_source : String,
  position : Int,
  to_source : String,
) -> Int? raise AlignMafError {
  let from = match block.component(from_source) {
    Some(component) => component
    None => return None
  }
  let to = match block.component(to_source) {
    Some(component) => component
    None => return None
  }
  match from.source_to_column(position) {
    Some(column) =>
      to.column_to_source(column) catch {
        BigMafError(message) => {
          align_maf_fail(message)
          None
        }
      }
    None => None
  }
}

///|
pub fn AlignMafDocument::sources(self : AlignMafDocument) -> Array[String] {
  let sources : Array[String] = []
  for block in self.blocks {
    for component in block.components {
      if !sources.contains(component.source) {
        sources.push(component.source)
      }
    }
    for empty in block.empty_components {
      if !sources.contains(empty.source) {
        sources.push(empty.source)
      }
    }
  }
  sources
}

///|
pub fn AlignMafDocument::summary(self : AlignMafDocument) -> AlignMafSummary {
  let mut component_count = 0
  let mut empty_component_count = 0
  let mut aligned_columns = 0
  let mut reference_bases = 0
  for block in self.blocks {
    component_count = component_count + block.components.length()
    empty_component_count = empty_component_count +
      block.empty_components.length()
    aligned_columns = aligned_columns + block.aligned_columns()
    if block.components.length() > 0 {
      reference_bases = reference_bases + block.components[0].size
    }
  }
  AlignMafSummary::{
    block_count: self.blocks.length(),
    component_count,
    empty_component_count,
    aligned_columns,
    reference_bases,
    source_count: self.sources().length(),
  }
}

///|
fn align_maf_entry_sort(entries : Array[AlignMafIndexEntry]) -> Unit {
  for index in 1..<entries.length() {
    let value = entries[index]
    let mut position = index
    while position > 0 &&
          (
            entries[position - 1].start > value.start ||
            (
              entries[position - 1].start == value.start &&
              entries[position - 1].end > value.end
            )
          ) {
      entries[position] = entries[position - 1]
      position = position - 1
    }
    entries[position] = value
  }
}

///|
pub fn AlignMafIndex::create(
  document : AlignMafDocument,
  reference : String,
) -> AlignMafIndex raise AlignMafError {
  align_maf_validate_token(reference, "MAF index reference")
  if document.blocks.length() == 0 {
    align_maf_fail("MAF index requires at least one alignment block")
  }
  let entries : Array[AlignMafIndexEntry] = []
  let mut source_size = -1
  let mut reference_strand = ""
  for block_index in 0..<document.blocks.length() {
    let block = document.blocks[block_index]
    let mut found : BigMafComponent? = None
    for component in block.components {
      if component.source == reference {
        if found is Some(_) {
          align_maf_fail(
            "MAF index reference occurs more than once in one block",
          )
        }
        found = Some(component)
      }
    }
    let component = match found {
      Some(value) => value
      None => {
        align_maf_fail("MAF index reference is missing from an alignment block")
        block.components[0]
      }
    }
    if source_size < 0 {
      source_size = component.source_size
      reference_strand = component.strand
    } else {
      if component.source_size != source_size {
        align_maf_fail("MAF index reference source sizes are inconsistent")
      }
      if component.strand != reference_strand {
        align_maf_fail("MAF index reference strands are inconsistent")
      }
    }
    let (start, end) = component.forward_interval()
    entries.push(AlignMafIndexEntry::{ block_index, start, end, })
  }
  align_maf_entry_sort(entries)
  AlignMafIndex::{
    reference,
    source_size,
    reference_strand,
    entries,
    blocks: align_maf_copy_blocks(document.blocks),
  }
}

///|
fn align_maf_validate_range(
  start : Int,
  end : Int,
  source_size : Int,
) -> Unit raise AlignMafError {
  if start < 0 || end <= start || end > source_size {
    align_maf_fail(
      "MAF query interval must satisfy 0 <= start < end <= source size",
    )
  }
}

///|
pub fn AlignMafIndex::search(
  self : AlignMafIndex,
  start : Int,
  end : Int,
) -> Array[BigMafBlock] raise AlignMafError {
  align_maf_validate_range(start, end, self.source_size)
  let blocks : Array[BigMafBlock] = []
  for entry in self.entries {
    if entry.start < end && entry.end > start {
      blocks.push(self.blocks[entry.block_index])
    }
  }
  blocks
}

///|
pub fn AlignMafIndex::search_ranges(
  self : AlignMafIndex,
  starts : Array[Int],
  ends : Array[Int],
) -> Array[BigMafBlock] raise AlignMafError {
  if starts.length() == 0 || starts.length() != ends.length() {
    align_maf_fail("MAF range arrays must be non-empty and have equal lengths")
  }
  let selected = Array::make(self.blocks.length(), false)
  for index in 0..<starts.length() {
    align_maf_validate_range(starts[index], ends[index], self.source_size)
    for entry in self.entries {
      if entry.start < ends[index] && entry.end > starts[index] {
        selected[entry.block_index] = true
      }
    }
  }
  let blocks : Array[BigMafBlock] = []
  for index in 0..<self.blocks.length() {
    if selected[index] {
      blocks.push(self.blocks[index])
    }
  }
  blocks
}

///|
pub fn AlignMafIndex::entry_count(self : AlignMafIndex) -> Int {
  self.entries.length()
}

///|
fn align_maf_repeat(value : String, count : Int) -> String {
  let output = StringBuilder::new(size_hint=count)
  for _ in 0..<count {
    output.write_string(value)
  }
  output.to_string()
}

///|
fn align_maf_splice_sources(index : AlignMafIndex) -> Array[String] {
  let sources : Array[String] = [index.reference]
  for block in index.blocks {
    for component in block.components {
      if component.source != index.reference &&
        !sources.contains(component.source) {
        sources.push(component.source)
      }
    }
  }
  sources
}

///|
fn align_maf_append_missing(texts : Array[String], count : Int) -> Unit {
  if count <= 0 {
    return
  }
  texts[0] = texts[0] + align_maf_repeat("N", count)
  for index in 1..<texts.length() {
    texts[index] = texts[index] + align_maf_repeat("-", count)
  }
}

///|
pub fn AlignMafIndex::get_spliced(
  self : AlignMafIndex,
  starts : Array[Int],
  ends : Array[Int],
  strand? : String = "+",
) -> AlignMafSplicedAlignment raise AlignMafError {
  if starts.length() == 0 || starts.length() != ends.length() {
    align_maf_fail("MAF exon arrays must be non-empty and have equal lengths")
  }
  if strand != "+" && strand != "-" {
    align_maf_fail("MAF spliced strand must be '+' or '-'")
  }
  if self.reference_strand != "+" {
    align_maf_fail(
      "MAF splicing currently requires a plus-strand reference component",
    )
  }
  for index in 0..<starts.length() {
    align_maf_validate_range(starts[index], ends[index], self.source_size)
    if index > 0 && starts[index] < ends[index - 1] {
      align_maf_fail("MAF exon intervals must be sorted and non-overlapping")
    }
  }
  let sources = align_maf_splice_sources(self)
  let texts = Array::make(sources.length(), "")
  for exon_index in 0..<starts.length() {
    let exon_start = starts[exon_index]
    let exon_end = ends[exon_index]
    let mut current = exon_start
    for entry in self.entries {
      if entry.end <= exon_start || entry.start >= exon_end {
        continue
      }
      let overlap_start = if entry.start > exon_start {
        entry.start
      } else {
        exon_start
      }
      let overlap_end = if entry.end < exon_end { entry.end } else { exon_end }
      if overlap_start < current {
        align_maf_fail(
          "Overlapping MAF reference blocks make splicing ambiguous",
        )
      }
      align_maf_append_missing(texts, overlap_start - current)
      let block = self.blocks[entry.block_index]
      let reference = match block.component(self.reference) {
        Some(component) => component
        None => {
          align_maf_fail("Indexed MAF block lost its reference component")
          block.components[0]
        }
      }
      let first_column = match reference.source_to_column(overlap_start) {
        Some(column) => column
        None => {
          align_maf_fail("MAF reference interval does not map to a column")
          0
        }
      }
      let last_column = match reference.source_to_column(overlap_end - 1) {
        Some(column) => column
        None => {
          align_maf_fail("MAF reference interval does not map to a column")
          0
        }
      }
      let column_start = if first_column < last_column {
        first_column
      } else {
        last_column
      }
      let column_end = if first_column > last_column {
        first_column + 1
      } else {
        last_column + 1
      }
      let width = column_end - column_start
      for source_index in 0..<sources.length() {
        match block.component(sources[source_index]) {
          Some(component) =>
            texts[source_index] = texts[source_index] +
              component.text[column_start:column_end].to_owned()
          None =>
            texts[source_index] = texts[source_index] +
              align_maf_repeat("-", width)
        }
      }
      current = overlap_end
    }
    align_maf_append_missing(texts, exon_end - current)
  }
  if strand == "-" {
    for index in 0..<texts.length() {
      texts[index] = reverse_complement(texts[index])
    }
  }
  let sequences : Array[AlignMafSplicedSequence] = []
  for index in 0..<sources.length() {
    sequences.push(AlignMafSplicedSequence::{
      source: sources[index],
      text: texts[index],
    })
  }
  AlignMafSplicedAlignment::{
    reference: self.reference,
    strand,
    starts: align_maf_copy_ints(starts),
    ends: align_maf_copy_ints(ends),
    sequences,
    columns: if texts.length() == 0 {
      0
    } else {
      texts[0].length()
    },
  }
}

///|
pub fn AlignMafSplicedAlignment::sequence(
  self : AlignMafSplicedAlignment,
  source : String,
) -> String? {
  for sequence in self.sequences {
    if sequence.source == source {
      return Some(sequence.text)
    }
  }
  None
}

///|
pub fn AlignMafSplicedAlignment::ungapped_length(
  self : AlignMafSplicedAlignment,
  source : String,
) -> Int? {
  match self.sequence(source) {
    Some(text) => {
      let mut count = 0
      for index in 0..<text.length() {
        let code = text.unsafe_get(index).to_int()
        if code != '-'.to_int() && code != 'N'.to_int() {
          count = count + 1
        }
      }
      Some(count)
    }
    None => None
  }
}

///|
pub fn align_maf_example_text() -> String {
  "track name=\"Multiz demo\" description=\"MAF coordinate example\" " +
  "mafDot=on visibility=pack speciesOrder=\"hg38 mm39 rn7\"\n" +
  "##maf version=1 scoring=blastz program=multiz\n" +
  "# generated fixture\n\n" +
  "a score=23262.0 pass=1\n" +
  "s hg38.chr7 100 8 + 1000 AC-TGCAAT\n" +
  "i hg38.chr7 C 0 I 2\n" +
  "q hg38.chr7 98-765432\n" +
  "s mm39.chr5 50 8 - 500 ACCTG-AAT\n" +
  "e rn7.chr1 80 6 + 700 I\n\n" +
  "a score=125.5\n" +
  "s hg38.chr7 110 4 + 1000 GGTT\n" +
  "s mm39.chr5 70 3 - 500 G-TT\n\n" +
  "a score=80\n" +
  "s hg38.chr7 120 4 + 1000 AAAA\n" +
  "s rn7.chr1 90 4 + 700 AATA\n"
}
