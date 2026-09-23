// Alignment-aware SAM parsing, coordinate paths, typed tags, and writing.
//
// This module follows Biopython 1.86 Bio.Align.sam semantics. SAM POS and
// PNEXT fields are converted from one-based coordinates to zero-based values.

///|
pub suberror AlignSamError {
  AlignSamError(String)
}

///|
pub(all) enum AlignSamOperation {
  SamAligned
  SamInsertion
  SamDeletion
  SamSkipped
  SamEqual
  SamMismatch
} derive(Eq, Debug)
pub extend AlignSamOperation with Eq::{equal, not_equal}
pub extend AlignSamOperation with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub(all) enum AlignSamTagValue {
  SamTagInteger(Int)
  SamTagFloat(Double)
  SamTagCharacter(String)
  SamTagString(String)
  SamTagHex(String)
  SamTagIntegerArray(String, Array[Int])
  SamTagFloatArray(Array[Double])
} derive(Debug)
pub extend AlignSamTagValue with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct AlignSamHeaderField {
  key : String
  value : String
} derive(Eq, Debug)
pub extend AlignSamHeaderField with Eq::{equal, not_equal}
pub extend AlignSamHeaderField with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct AlignSamHeader {
  record_type : String
  fields : Array[AlignSamHeaderField]
  comment : String
} derive(Debug)
pub extend AlignSamHeader with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct AlignSamReference {
  name : String
  length : Int
  annotations : Array[AlignSamHeaderField]
} derive(Debug)
pub extend AlignSamReference with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct AlignSamCigarElement {
  operation : String
  length : Int
} derive(Eq, Debug)
pub extend AlignSamCigarElement with Eq::{equal, not_equal}
pub extend AlignSamCigarElement with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct AlignSamTag {
  name : String
  value : AlignSamTagValue
} derive(Debug)
pub extend AlignSamTag with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct AlignSamStats {
  reference_consumed : Int
  query_consumed : Int
  aligned : Int
  inserted : Int
  deleted : Int
  skipped : Int
  soft_clipped : Int
  hard_clipped : Int
} derive(Eq, Debug)
pub extend AlignSamStats with Eq::{equal, not_equal}
pub extend AlignSamStats with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct AlignSamAlignment {
  query_name : String
  flag : Int
  reference_name : String
  reference_length : Int?
  target_coordinates : Array[Int]
  query_coordinates : Array[Int]
  operations : Array[AlignSamOperation]
  cigar : Array[AlignSamCigarElement]
  mapq : Int?
  mate_reference : String?
  mate_position : Int?
  template_length : Int?
  query_sequence : String
  sequence_known : Bool
  qualities : Array[Int]
  qualities_known : Bool
  hard_clip_left : Int?
  hard_clip_right : Int?
  tags : Array[AlignSamTag]
}

///|
pub struct AlignSamDocument {
  headers : Array[AlignSamHeader]
  references : Array[AlignSamReference]
  alignments : Array[AlignSamAlignment]
}

///|
fn align_sam_fail(message : String) -> Unit raise AlignSamError {
  raise AlignSamError(message)
}

///|
fn align_sam_copy_ints(values : Array[Int]) -> Array[Int] {
  let copy : Array[Int] = []
  for value in values {
    copy.push(value)
  }
  copy
}

///|
fn align_sam_copy_doubles(values : Array[Double]) -> Array[Double] {
  let copy : Array[Double] = []
  for value in values {
    copy.push(value)
  }
  copy
}

///|
fn align_sam_copy_tag_value(value : AlignSamTagValue) -> AlignSamTagValue {
  match value {
    SamTagInteger(item) => SamTagInteger(item)
    SamTagFloat(item) => SamTagFloat(item)
    SamTagCharacter(item) => SamTagCharacter(item)
    SamTagString(item) => SamTagString(item)
    SamTagHex(item) => SamTagHex(item)
    SamTagIntegerArray(kind, items) =>
      SamTagIntegerArray(kind, align_sam_copy_ints(items))
    SamTagFloatArray(items) => SamTagFloatArray(align_sam_copy_doubles(items))
  }
}

///|
fn align_sam_copy_cigar(
  values : Array[AlignSamCigarElement],
) -> Array[AlignSamCigarElement] {
  let copy : Array[AlignSamCigarElement] = []
  for value in values {
    copy.push(AlignSamCigarElement::{
      operation: value.operation,
      length: value.length,
    })
  }
  copy
}

///|
fn align_sam_copy_fields(
  values : Array[AlignSamHeaderField],
) -> Array[AlignSamHeaderField] {
  let copy : Array[AlignSamHeaderField] = []
  for value in values {
    copy.push(AlignSamHeaderField::{ key: value.key, value: value.value, })
  }
  copy
}

///|
fn align_sam_strip_cr(value : String) -> String {
  if value.length() > 0 &&
    value.unsafe_get(value.length() - 1).to_int() == '\r'.to_int() {
    value[0:value.length() - 1].to_owned()
  } else {
    value
  }
}

///|
fn align_sam_split_lines(value : String) -> Array[String] {
  let lines : Array[String] = []
  let mut start = 0
  for index in 0..<value.length() {
    if value.unsafe_get(index).to_int() == '\n'.to_int() {
      lines.push(align_sam_strip_cr(value[start:index].to_owned()))
      start = index + 1
    }
  }
  if start < value.length() {
    lines.push(align_sam_strip_cr(value[start:value.length()].to_owned()))
  }
  lines
}

///|
fn align_sam_split_char(value : String, separator : Char) -> Array[String] {
  let parts : Array[String] = []
  let mut start = 0
  for index in 0..<value.length() {
    if value.unsafe_get(index).to_int() == separator.to_int() {
      parts.push(value[start:index].to_owned())
      start = index + 1
    }
  }
  parts.push(value[start:value.length()].to_owned())
  parts
}

///|
fn align_sam_find_char(value : String, separator : Char, start : Int) -> Int? {
  for index in start..<value.length() {
    if value.unsafe_get(index).to_int() == separator.to_int() {
      return Some(index)
    }
  }
  None
}

///|
fn align_sam_parse_int(
  value : String,
  label : String,
) -> Int raise AlignSamError {
  if value.length() == 0 {
    align_sam_fail(label + " is missing")
  }
  let mut index = 0
  let negative = value.unsafe_get(0).to_int() == '-'.to_int()
  if negative || value.unsafe_get(0).to_int() == '+'.to_int() {
    index = 1
  }
  if index == value.length() {
    align_sam_fail(label + " is not an integer")
  }
  if value == "-2147483648" {
    return -2147483647 - 1
  }
  let mut result = 0
  while index < value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      align_sam_fail(label + " is not an integer: " + value)
    }
    let digit = code - '0'.to_int()
    if result > (2147483647 - digit) / 10 {
      align_sam_fail(label + " exceeds the MoonBit Int range")
    }
    result = result * 10 + digit
    index = index + 1
  }
  if negative {
    -result
  } else {
    result
  }
}

///|
fn align_sam_parse_double(
  value : String,
  label : String,
) -> Double raise AlignSamError {
  match parse_double(value) {
    Some(result) => {
      if result.abs() > 1.0e300 {
        align_sam_fail(label + " must be finite")
      }
      result
    }
    None => {
      align_sam_fail(label + " is not a floating-point number: " + value)
      0.0
    }
  }
}

///|
fn align_sam_validate_token(
  value : String,
  label : String,
  allow_star : Bool,
) -> Unit raise AlignSamError {
  if value.length() == 0 {
    align_sam_fail(label + " must not be empty")
  }
  if value == "*" && allow_star {
    return
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code <= 32 || code == 127 {
      align_sam_fail(label + " contains whitespace or a control character")
    }
  }
}

///|
fn align_sam_validate_tag_name(name : String) -> Unit raise AlignSamError {
  if name.length() != 2 {
    align_sam_fail("SAM tag names must contain exactly two characters")
  }
  let first = name.unsafe_get(0).to_int()
  let second = name.unsafe_get(1).to_int()
  let first_ok = (first >= 'A'.to_int() && first <= 'Z'.to_int()) ||
    (first >= 'a'.to_int() && first <= 'z'.to_int())
  let second_ok = (second >= 'A'.to_int() && second <= 'Z'.to_int()) ||
    (second >= 'a'.to_int() && second <= 'z'.to_int()) ||
    (second >= '0'.to_int() && second <= '9'.to_int())
  if !first_ok || !second_ok {
    align_sam_fail("invalid SAM tag name " + name)
  }
}

///|
fn align_sam_is_hex(code : Int) -> Bool {
  (code >= '0'.to_int() && code <= '9'.to_int()) ||
  (code >= 'A'.to_int() && code <= 'F'.to_int()) ||
  (code >= 'a'.to_int() && code <= 'f'.to_int())
}

///|
fn align_sam_is_md_base(code : Int) -> Bool {
  code == 'A'.to_int() ||
  code == 'C'.to_int() ||
  code == 'G'.to_int() ||
  code == 'T'.to_int() ||
  code == 'N'.to_int() ||
  code == 'a'.to_int() ||
  code == 'c'.to_int() ||
  code == 'g'.to_int() ||
  code == 't'.to_int() ||
  code == 'n'.to_int()
}

///|
fn align_sam_operation_from_text(
  operation : String,
) -> AlignSamOperation raise AlignSamError {
  if operation == "M" {
    SamAligned
  } else if operation == "I" {
    SamInsertion
  } else if operation == "D" {
    SamDeletion
  } else if operation == "N" {
    SamSkipped
  } else if operation == "=" {
    SamEqual
  } else if operation == "X" {
    SamMismatch
  } else {
    align_sam_fail("CIGAR operation " + operation + " has no alignment path")
    SamAligned
  }
}

///|
fn align_sam_is_clip(operation : String) -> Bool {
  operation == "S" || operation == "H"
}

///|
fn align_sam_is_core(operation : String) -> Bool {
  operation == "M" ||
  operation == "I" ||
  operation == "D" ||
  operation == "N" ||
  operation == "=" ||
  operation == "X"
}

///|
pub fn align_sam_parse_cigar(
  text : String,
) -> Array[AlignSamCigarElement] raise AlignSamError {
  if text == "*" {
    return []
  }
  if text.length() == 0 {
    align_sam_fail("CIGAR must not be empty")
  }
  let elements : Array[AlignSamCigarElement] = []
  let mut start = 0
  let mut index = 0
  while index < text.length() {
    let code = text.unsafe_get(index).to_int()
    if code >= '0'.to_int() && code <= '9'.to_int() {
      index = index + 1
      continue
    }
    if index == start {
      align_sam_fail("CIGAR operation is missing a length")
    }
    let length = align_sam_parse_int(
      text[start:index].to_owned(),
      "CIGAR length",
    )
    if length <= 0 {
      align_sam_fail("CIGAR lengths must be positive")
    }
    let operation = text[index:index + 1].to_owned()
    if operation == "P" {
      align_sam_fail(
        "CIGAR padding operation P is not supported by Bio.Align.sam",
      )
    }
    if !align_sam_is_clip(operation) && !align_sam_is_core(operation) {
      align_sam_fail("unknown CIGAR operation " + operation)
    }
    elements.push(AlignSamCigarElement::{ operation, length, })
    index = index + 1
    start = index
  }
  if start != text.length() {
    align_sam_fail("CIGAR ends with a length but no operation")
  }
  let mut first_core = elements.length()
  let mut last_core = -1
  for element_index in 0..<elements.length() {
    if align_sam_is_core(elements[element_index].operation) {
      if first_core == elements.length() {
        first_core = element_index
      }
      last_core = element_index
    }
  }
  if last_core < 0 {
    align_sam_fail("mapped CIGAR must contain an alignment operation")
  }
  for element_index in first_core..<=last_core {
    if align_sam_is_clip(elements[element_index].operation) {
      align_sam_fail("soft and hard clipping are only valid at CIGAR ends")
    }
  }
  let mut seen_soft = false
  for element_index in 0..<first_core {
    let operation = elements[element_index].operation
    if operation == "S" {
      seen_soft = true
    } else if operation == "H" && seen_soft {
      align_sam_fail("left hard clipping must precede left soft clipping")
    }
  }
  seen_soft = false
  let mut element_index = elements.length()
  while element_index > last_core + 1 {
    element_index = element_index - 1
    let operation = elements[element_index].operation
    if operation == "S" {
      seen_soft = true
    } else if operation == "H" && !seen_soft {
      // Right hard clipping is outermost and therefore visited first.
    } else if operation == "H" {
      align_sam_fail("right hard clipping must follow right soft clipping")
    }
  }
  align_sam_copy_cigar(elements)
}

///|
pub fn align_sam_cigar_string(elements : Array[AlignSamCigarElement]) -> String {
  if elements.length() == 0 {
    return "*"
  }
  let output = StringBuilder::new()
  for element in elements {
    output.write_string(element.length.to_string())
    output.write_string(element.operation)
  }
  output.to_string()
}

///|
fn align_sam_cigar_query_consumed(
  elements : Array[AlignSamCigarElement],
) -> Int {
  let mut total = 0
  for element in elements {
    if element.operation == "M" ||
      element.operation == "I" ||
      element.operation == "S" ||
      element.operation == "=" ||
      element.operation == "X" {
      total = total + element.length
    }
  }
  total
}

///|
fn align_sam_left_soft(elements : Array[AlignSamCigarElement]) -> Int {
  let mut total = 0
  for element in elements {
    if align_sam_is_core(element.operation) {
      break
    }
    if element.operation == "S" {
      total = total + element.length
    }
  }
  total
}

///|
fn align_sam_raw_hard_clips(
  elements : Array[AlignSamCigarElement],
) -> (Int?, Int?) {
  let mut left = 0
  let mut right = 0
  for element in elements {
    if align_sam_is_core(element.operation) {
      break
    }
    if element.operation == "H" {
      left = left + element.length
    }
  }
  let mut index = elements.length()
  while index > 0 {
    index = index - 1
    let element = elements[index]
    if align_sam_is_core(element.operation) {
      break
    }
    if element.operation == "H" {
      right = right + element.length
    }
  }
  (
    if left == 0 {
      None
    } else {
      Some(left)
    },
    if right == 0 {
      None
    } else {
      Some(right)
    },
  )
}

///|
fn align_sam_build_path(
  target_start : Int,
  cigar : Array[AlignSamCigarElement],
  query_length : Int,
  reverse : Bool,
) -> (Array[Int], Array[Int], Array[AlignSamOperation]) raise AlignSamError {
  let target_coordinates : Array[Int] = [target_start]
  let query_coordinates : Array[Int] = [align_sam_left_soft(cigar)]
  let operations : Array[AlignSamOperation] = []
  let mut target = target_start
  let mut query = query_coordinates[0]
  for element in cigar {
    let operation = element.operation
    if operation == "S" || operation == "H" {
      continue
    }
    if operation == "M" || operation == "=" || operation == "X" {
      target = target + element.length
      query = query + element.length
    } else if operation == "I" {
      query = query + element.length
    } else if operation == "D" || operation == "N" {
      target = target + element.length
    }
    target_coordinates.push(target)
    query_coordinates.push(query)
    operations.push(align_sam_operation_from_text(operation))
  }
  if query_length != align_sam_cigar_query_consumed(cigar) {
    align_sam_fail(
      "query length does not match query-consuming CIGAR operations",
    )
  }
  if reverse {
    for index in 0..<query_coordinates.length() {
      query_coordinates[index] = query_length - query_coordinates[index]
    }
  }
  (target_coordinates, query_coordinates, operations)
}

///|
fn align_sam_parse_tag(field : String) -> AlignSamTag raise AlignSamError {
  let first = match align_sam_find_char(field, ':', 0) {
    Some(value) => value
    None => {
      align_sam_fail("SAM optional tag is missing its datatype")
      0
    }
  }
  let second = match align_sam_find_char(field, ':', first + 1) {
    Some(value) => value
    None => {
      align_sam_fail("SAM optional tag is missing its value")
      0
    }
  }
  let name = field[0:first].to_owned()
  let datatype = field[first + 1:second].to_owned()
  let text = field[second + 1:field.length()].to_owned()
  align_sam_validate_tag_name(name)
  let value = if datatype == "i" {
    SamTagInteger(align_sam_parse_int(text, "SAM integer tag " + name))
  } else if datatype == "f" {
    SamTagFloat(align_sam_parse_double(text, "SAM float tag " + name))
  } else if datatype == "A" {
    if text.length() != 1 {
      align_sam_fail("SAM A tag " + name + " must contain one character")
    }
    let code = text.unsafe_get(0).to_int()
    if code < 33 || code > 126 {
      align_sam_fail("SAM A tag " + name + " must be printable ASCII")
    }
    SamTagCharacter(text)
  } else if datatype == "Z" {
    for index in 0..<text.length() {
      let code = text.unsafe_get(index).to_int()
      if code < 32 || code > 126 {
        align_sam_fail("SAM Z tag " + name + " contains invalid text")
      }
    }
    SamTagString(text)
  } else if datatype == "H" {
    if text.length() % 2 != 0 {
      align_sam_fail("SAM H tag " + name + " must contain byte pairs")
    }
    for index in 0..<text.length() {
      if !align_sam_is_hex(text.unsafe_get(index).to_int()) {
        align_sam_fail("SAM H tag " + name + " contains non-hexadecimal text")
      }
    }
    SamTagHex(text)
  } else if datatype == "B" {
    let values = align_sam_split_char(text, ',')
    if values.length() < 2 || values[0].length() != 1 {
      align_sam_fail("SAM B tag " + name + " is missing its array subtype")
    }
    let subtype = values[0]
    if subtype == "f" {
      let items : Array[Double] = []
      for index in 1..<values.length() {
        items.push(align_sam_parse_double(values[index], "SAM B:f tag " + name))
      }
      SamTagFloatArray(items)
    } else if subtype == "c" ||
      subtype == "C" ||
      subtype == "s" ||
      subtype == "S" ||
      subtype == "i" ||
      subtype == "I" {
      let items : Array[Int] = []
      for index in 1..<values.length() {
        let item = align_sam_parse_int(values[index], "SAM B tag " + name)
        if subtype == "c" && (item < -128 || item > 127) {
          align_sam_fail("SAM B:c value is outside the signed byte range")
        }
        if subtype == "C" && (item < 0 || item > 255) {
          align_sam_fail("SAM B:C value is outside the unsigned byte range")
        }
        if subtype == "s" && (item < -32768 || item > 32767) {
          align_sam_fail("SAM B:s value is outside the signed short range")
        }
        if subtype == "S" && (item < 0 || item > 65535) {
          align_sam_fail("SAM B:S value is outside the unsigned short range")
        }
        if subtype == "I" && item < 0 {
          align_sam_fail("SAM B:I values must be non-negative")
        }
        items.push(item)
      }
      SamTagIntegerArray(subtype, items)
    } else {
      align_sam_fail("unknown SAM B tag subtype " + subtype)
      SamTagIntegerArray("i", [])
    }
  } else {
    align_sam_fail("unknown SAM optional tag datatype " + datatype)
    SamTagString("")
  }
  AlignSamTag::{ name, value, }
}

///|
fn align_sam_tag_text(tag : AlignSamTag) -> String {
  let prefix = tag.name + ":"
  match tag.value {
    SamTagInteger(value) => prefix + "i:" + value.to_string()
    SamTagFloat(value) => prefix + "f:" + value.to_string()
    SamTagCharacter(value) => prefix + "A:" + value
    SamTagString(value) => prefix + "Z:" + value
    SamTagHex(value) => prefix + "H:" + value
    SamTagIntegerArray(kind, values) => {
      let output = StringBuilder::new()
      output.write_string(prefix)
      output.write_string("B:")
      output.write_string(kind)
      for value in values {
        output.write_char(',')
        output.write_string(value.to_string())
      }
      output.to_string()
    }
    SamTagFloatArray(values) => {
      let output = StringBuilder::new()
      output.write_string(prefix)
      output.write_string("B:f")
      for value in values {
        output.write_char(',')
        output.write_string(value.to_string())
      }
      output.to_string()
    }
  }
}

///|
pub fn align_sam_integer_tag(
  name : String,
  value : Int,
) -> AlignSamTag raise AlignSamError {
  align_sam_validate_tag_name(name)
  AlignSamTag::{ name, value: SamTagInteger(value), }
}

///|
pub fn align_sam_float_tag(
  name : String,
  value : Double,
) -> AlignSamTag raise AlignSamError {
  align_sam_validate_tag_name(name)
  if value.abs() > 1.0e300 {
    align_sam_fail("SAM floating-point tag value must be finite")
  }
  AlignSamTag::{ name, value: SamTagFloat(value), }
}

///|
pub fn align_sam_string_tag(
  name : String,
  value : String,
) -> AlignSamTag raise AlignSamError {
  align_sam_parse_tag(name + ":Z:" + value)
}

///|
pub fn align_sam_character_tag(
  name : String,
  value : String,
) -> AlignSamTag raise AlignSamError {
  align_sam_parse_tag(name + ":A:" + value)
}

///|
pub fn align_sam_hex_tag(
  name : String,
  value : String,
) -> AlignSamTag raise AlignSamError {
  align_sam_parse_tag(name + ":H:" + value)
}

///|
pub fn align_sam_integer_array_tag(
  name : String,
  subtype : String,
  values : Array[Int],
) -> AlignSamTag raise AlignSamError {
  let output = StringBuilder::new()
  output.write_string(name)
  output.write_string(":B:")
  output.write_string(subtype)
  for value in values {
    output.write_char(',')
    output.write_string(value.to_string())
  }
  align_sam_parse_tag(output.to_string())
}

///|
pub fn align_sam_float_array_tag(
  name : String,
  values : Array[Double],
) -> AlignSamTag raise AlignSamError {
  let output = StringBuilder::new()
  output.write_string(name)
  output.write_string(":B:f")
  for value in values {
    output.write_char(',')
    output.write_string(value.to_string())
  }
  align_sam_parse_tag(output.to_string())
}

///|
fn align_sam_parse_header(line : String) -> AlignSamHeader raise AlignSamError {
  if line.length() < 3 || line.unsafe_get(0).to_int() != '@'.to_int() {
    align_sam_fail("invalid SAM header line")
  }
  let parts = align_sam_split_char(line, '\t')
  let record_type = parts[0][1:parts[0].length()].to_owned()
  if record_type.length() != 2 {
    align_sam_fail("SAM header record type must contain two characters")
  }
  if record_type == "CO" {
    let mut comment = ""
    for index in 1..<parts.length() {
      if index > 1 {
        comment = comment + "\t"
      }
      comment = comment + parts[index]
    }
    return AlignSamHeader::{ record_type, fields: [], comment, }
  }
  let fields : Array[AlignSamHeaderField] = []
  let seen : Map[String, Bool] = Map([])
  for index in 1..<parts.length() {
    let part = parts[index]
    let separator = match align_sam_find_char(part, ':', 0) {
      Some(value) => value
      None => {
        align_sam_fail("SAM header field is missing ':'")
        0
      }
    }
    let key = part[0:separator].to_owned()
    let value = part[separator + 1:part.length()].to_owned()
    if key.length() != 2 || value.length() == 0 {
      align_sam_fail("invalid SAM header field " + part)
    }
    if seen.contains(key) {
      align_sam_fail("duplicate SAM header field " + key)
    }
    seen[key] = true
    fields.push(AlignSamHeaderField::{ key, value, })
  }
  AlignSamHeader::{ record_type, fields, comment: "", }
}

///|
fn align_sam_header_value(header : AlignSamHeader, key : String) -> String? {
  for field in header.fields {
    if field.key == key {
      return Some(field.value)
    }
  }
  None
}

///|
fn align_sam_parse_reference(
  header : AlignSamHeader,
) -> AlignSamReference raise AlignSamError {
  let name = match align_sam_header_value(header, "SN") {
    Some(value) => value
    None => {
      align_sam_fail("@SQ header is missing SN")
      ""
    }
  }
  align_sam_validate_token(name, "@SQ SN", false)
  let length = match align_sam_header_value(header, "LN") {
    Some(value) => align_sam_parse_int(value, "@SQ LN")
    None => {
      align_sam_fail("@SQ header is missing LN")
      0
    }
  }
  if length <= 0 {
    align_sam_fail("@SQ LN must be positive")
  }
  let annotations : Array[AlignSamHeaderField] = []
  for field in header.fields {
    if field.key != "SN" && field.key != "LN" {
      if field.key == "TP" &&
        field.value != "linear" &&
        field.value != "circular" {
        align_sam_fail("@SQ TP must be linear or circular")
      }
      annotations.push(AlignSamHeaderField::{
        key: field.key,
        value: field.value,
      })
    }
  }
  AlignSamReference::{ name, length, annotations, }
}

///|
fn align_sam_validate_sequence(sequence : String) -> Unit raise AlignSamError {
  for index in 0..<sequence.length() {
    let code = sequence.unsafe_get(index).to_int()
    let valid = (code >= 'A'.to_int() && code <= 'Z'.to_int()) ||
      (code >= 'a'.to_int() && code <= 'z'.to_int()) ||
      code == '='.to_int() ||
      code == '.'.to_int()
    if !valid {
      align_sam_fail("SAM SEQ contains an invalid character")
    }
  }
}

///|
fn align_sam_reverse_ints(values : Array[Int]) -> Array[Int] {
  let reversed : Array[Int] = []
  let mut index = values.length()
  while index > 0 {
    index = index - 1
    reversed.push(values[index])
  }
  reversed
}

///|
fn align_sam_reverse_complement(sequence : String) -> String {
  Seq::new(sequence).reverse_complement().to_string()
}

///|
fn align_sam_reference_length(
  references : Array[AlignSamReference],
  name : String,
) -> Int? {
  for reference in references {
    if reference.name == name {
      return Some(reference.length)
    }
  }
  None
}

///|
fn align_sam_validate_md(
  md : String,
  cigar : Array[AlignSamCigarElement],
) -> Unit raise AlignSamError {
  if md.length() == 0 {
    align_sam_fail("MD tag must not be empty")
  }
  let expected : Array[String] = []
  for element in cigar {
    if element.operation == "M" ||
      element.operation == "=" ||
      element.operation == "X" ||
      element.operation == "D" {
      for _ in 0..<element.length {
        expected.push(element.operation)
      }
    }
  }
  let mut pointer = 0
  let mut index = 0
  let mut expect_count = true
  while index < md.length() {
    let code = md.unsafe_get(index).to_int()
    if expect_count {
      if code < '0'.to_int() || code > '9'.to_int() {
        align_sam_fail(
          "MD tag must start with and alternate through match counts",
        )
      }
      let start = index
      while index < md.length() {
        let digit = md.unsafe_get(index).to_int()
        if digit < '0'.to_int() || digit > '9'.to_int() {
          break
        }
        index = index + 1
      }
      let matches = align_sam_parse_int(
        md[start:index].to_owned(),
        "MD match count",
      )
      for _ in 0..<matches {
        if pointer >= expected.length() {
          align_sam_fail("MD match count exceeds the CIGAR reference path")
        }
        if expected[pointer] == "D" || expected[pointer] == "X" {
          align_sam_fail("MD match count conflicts with the CIGAR operation")
        }
        pointer = pointer + 1
      }
      expect_count = false
      continue
    }
    if code == '^'.to_int() {
      index = index + 1
      let start = index
      while index < md.length() &&
            align_sam_is_md_base(md.unsafe_get(index).to_int()) {
        if pointer >= expected.length() || expected[pointer] != "D" {
          align_sam_fail("MD deletion does not match a CIGAR D operation")
        }
        pointer = pointer + 1
        index = index + 1
      }
      if index == start {
        align_sam_fail("MD deletion marker must be followed by bases")
      }
      expect_count = true
      continue
    }
    if align_sam_is_md_base(code) {
      if pointer >= expected.length() || expected[pointer] == "D" {
        align_sam_fail("MD mismatch does not match an aligned CIGAR operation")
      }
      if expected[pointer] == "=" {
        align_sam_fail("MD mismatch conflicts with a CIGAR = operation")
      }
      pointer = pointer + 1
      index = index + 1
      expect_count = true
      continue
    }
    align_sam_fail("MD tag contains an invalid character")
  }
  if expect_count || pointer != expected.length() {
    align_sam_fail(
      "MD tag consumes " +
      pointer.to_string() +
      " reference bases but CIGAR requires " +
      expected.length().to_string(),
    )
  }
}

///|
fn align_sam_parse_alignment(
  line : String,
  references : Array[AlignSamReference],
) -> AlignSamAlignment raise AlignSamError {
  let fields = align_sam_split_char(line, '\t')
  if fields.length() < 11 {
    align_sam_fail(
      "SAM alignment has " +
      fields.length().to_string() +
      " fields; expected at least 11",
    )
  }
  let query_name = fields[0]
  align_sam_validate_token(query_name, "QNAME", false)
  let flag = align_sam_parse_int(fields[1], "FLAG")
  if flag < 0 || flag > 65535 {
    align_sam_fail("FLAG must be between 0 and 65535")
  }
  let unmapped = (flag & 0x4) != 0
  let reverse = (flag & 0x10) != 0
  let reference_name = fields[2]
  align_sam_validate_token(reference_name, "RNAME", true)
  let one_based_position = align_sam_parse_int(fields[3], "POS")
  let raw_mapq = align_sam_parse_int(fields[4], "MAPQ")
  if raw_mapq < 0 || raw_mapq > 255 {
    align_sam_fail("MAPQ must be between 0 and 255")
  }
  let cigar = align_sam_parse_cigar(fields[5]) catch {
    AlignSamError(message) =>
      if unmapped && fields[5] == "*" {
        []
      } else {
        raise AlignSamError(message)
      }
  }
  let raw_mate_reference = fields[6]
  align_sam_validate_token(raw_mate_reference, "RNEXT", true)
  let raw_mate_position = align_sam_parse_int(fields[7], "PNEXT")
  let raw_template_length = align_sam_parse_int(fields[8], "TLEN")
  let raw_sequence = fields[9]
  let sequence_known = raw_sequence != "*"
  if sequence_known {
    align_sam_validate_sequence(raw_sequence)
  }
  let quality_text = fields[10]
  let qualities_known = quality_text != "*"
  let raw_qualities : Array[Int] = []
  if qualities_known {
    if !sequence_known || quality_text.length() != raw_sequence.length() {
      align_sam_fail("QUAL length must equal SEQ length")
    }
    for index in 0..<quality_text.length() {
      let code = quality_text.unsafe_get(index).to_int()
      if code < 33 || code > 126 {
        align_sam_fail("QUAL contains a character outside PHRED+33 range")
      }
      raw_qualities.push(code - 33)
    }
  }
  let tags : Array[AlignSamTag] = []
  let seen_tags : Map[String, Bool] = Map([])
  let mut md : String? = None
  for index in 11..<fields.length() {
    let tag = align_sam_parse_tag(fields[index])
    if seen_tags.contains(tag.name) {
      align_sam_fail("duplicate SAM optional tag " + tag.name)
    }
    seen_tags[tag.name] = true
    if tag.name == "MD" {
      match tag.value {
        SamTagString(value) => md = Some(value)
        _ => align_sam_fail("MD tag must use Z datatype")
      }
    }
    if tag.name == "AS" {
      match tag.value {
        SamTagInteger(_) => ()
        _ => align_sam_fail("AS tag must use i datatype")
      }
    }
    tags.push(tag)
  }
  if unmapped {
    if reference_name != "*" || one_based_position != 0 || fields[5] != "*" {
      align_sam_fail(
        "unmapped SAM records must use RNAME *, POS 0, and CIGAR *",
      )
    }
  } else {
    if reference_name == "*" {
      align_sam_fail("mapped SAM records require a reference name")
    }
    if one_based_position <= 0 {
      align_sam_fail("mapped SAM POS must be positive")
    }
    if fields[5] == "*" {
      align_sam_fail("mapped SAM records require a CIGAR")
    }
  }
  if raw_mate_position < 0 {
    align_sam_fail("PNEXT must be non-negative")
  }
  let query_length = if sequence_known {
    raw_sequence.length()
  } else if unmapped {
    0
  } else {
    align_sam_cigar_query_consumed(cigar)
  }
  let reference_length : Int? = if unmapped {
    None
  } else {
    let length = align_sam_reference_length(references, reference_name)
    if references.length() > 0 && length is None {
      align_sam_fail(
        "mapped reference " + reference_name + " is absent from @SQ",
      )
    }
    length
  }
  let (target_coordinates, query_coordinates, operations) = if unmapped {
    ([], [], [])
  } else {
    let target_start = one_based_position - 1
    let path = align_sam_build_path(target_start, cigar, query_length, reverse)
    match reference_length {
      Some(length) =>
        if path.0[path.0.length() - 1] > length {
          align_sam_fail("alignment extends beyond its @SQ reference length")
        }
      None => ()
    }
    path
  }
  match md {
    Some(value) => {
      if unmapped {
        align_sam_fail("unmapped SAM records cannot carry an MD tag")
      }
      align_sam_validate_md(value, cigar)
    }
    None => ()
  }
  let raw_hard_clips : (Int?, Int?) = if unmapped {
    (None, None)
  } else {
    align_sam_raw_hard_clips(cigar)
  }
  let raw_hard_left = raw_hard_clips.0
  let raw_hard_right = raw_hard_clips.1
  let hard_clip_left = if reverse { raw_hard_right } else { raw_hard_left }
  let hard_clip_right = if reverse { raw_hard_left } else { raw_hard_right }
  let query_sequence = if !sequence_known {
    ""
  } else if reverse {
    align_sam_reverse_complement(raw_sequence)
  } else {
    raw_sequence
  }
  let qualities = if !qualities_known {
    []
  } else if reverse {
    align_sam_reverse_ints(raw_qualities)
  } else {
    raw_qualities
  }
  let mate_reference : String? = if raw_mate_reference == "*" {
    None
  } else if raw_mate_reference == "=" {
    if reference_name == "*" {
      align_sam_fail("RNEXT = requires a mapped reference")
    }
    Some(reference_name)
  } else {
    Some(raw_mate_reference)
  }
  AlignSamAlignment::{
    query_name,
    flag,
    reference_name,
    reference_length,
    target_coordinates,
    query_coordinates,
    operations,
    cigar,
    mapq: if raw_mapq == 255 {
      None
    } else {
      Some(raw_mapq)
    },
    mate_reference,
    mate_position: if raw_mate_position == 0 {
      None
    } else {
      Some(raw_mate_position - 1)
    },
    template_length: if raw_template_length == 0 {
      None
    } else {
      Some(raw_template_length)
    },
    query_sequence,
    sequence_known,
    qualities,
    qualities_known,
    hard_clip_left,
    hard_clip_right,
    tags,
  }
}

///|
pub fn align_sam_parse(text : String) -> AlignSamDocument raise AlignSamError {
  let headers : Array[AlignSamHeader] = []
  let references : Array[AlignSamReference] = []
  let alignments : Array[AlignSamAlignment] = []
  let reference_names : Map[String, Bool] = Map([])
  let mut saw_alignment = false
  let mut saw_hd = false
  for line in align_sam_split_lines(text) {
    if line.length() == 0 {
      continue
    }
    if line.unsafe_get(0).to_int() == '@'.to_int() {
      if saw_alignment {
        align_sam_fail("SAM headers must precede alignment records")
      }
      let header = align_sam_parse_header(line)
      if header.record_type == "HD" {
        if saw_hd || headers.length() > 0 {
          align_sam_fail("@HD must be the first and only @HD record")
        }
        saw_hd = true
        if align_sam_header_value(header, "VN") is None {
          align_sam_fail("@HD header is missing VN")
        }
      }
      if header.record_type == "SQ" {
        let reference = align_sam_parse_reference(header)
        if reference_names.contains(reference.name) {
          align_sam_fail("duplicate @SQ reference " + reference.name)
        }
        reference_names[reference.name] = true
        references.push(reference)
      }
      headers.push(header)
    } else {
      saw_alignment = true
      alignments.push(align_sam_parse_alignment(line, references))
    }
  }
  AlignSamDocument::{ headers, references, alignments, }
}

///|
pub fn AlignSamDocument::reference(
  self : AlignSamDocument,
  name : String,
) -> AlignSamReference? {
  for reference in self.references {
    if reference.name == name {
      return Some(AlignSamReference::{
        name: reference.name,
        length: reference.length,
        annotations: align_sam_copy_fields(reference.annotations),
      })
    }
  }
  None
}

///|
pub fn AlignSamAlignment::is_unmapped(self : AlignSamAlignment) -> Bool {
  (self.flag & 0x4) != 0
}

///|
pub fn AlignSamAlignment::is_reverse(self : AlignSamAlignment) -> Bool {
  (self.flag & 0x10) != 0
}

///|
pub fn AlignSamAlignment::is_paired(self : AlignSamAlignment) -> Bool {
  (self.flag & 0x1) != 0
}

///|
pub fn AlignSamAlignment::is_secondary(self : AlignSamAlignment) -> Bool {
  (self.flag & 0x100) != 0
}

///|
pub fn AlignSamAlignment::is_supplementary(self : AlignSamAlignment) -> Bool {
  (self.flag & 0x800) != 0
}

///|
pub fn AlignSamAlignment::query_length(self : AlignSamAlignment) -> Int {
  if self.sequence_known {
    self.query_sequence.length()
  } else {
    align_sam_cigar_query_consumed(self.cigar)
  }
}

///|
pub fn AlignSamAlignment::reference_start(self : AlignSamAlignment) -> Int? {
  if self.target_coordinates.length() == 0 {
    None
  } else {
    Some(self.target_coordinates[0])
  }
}

///|
pub fn AlignSamAlignment::reference_end(self : AlignSamAlignment) -> Int? {
  if self.target_coordinates.length() == 0 {
    None
  } else {
    Some(self.target_coordinates[self.target_coordinates.length() - 1])
  }
}

///|
pub fn AlignSamAlignment::tag(
  self : AlignSamAlignment,
  name : String,
) -> AlignSamTag? {
  for tag in self.tags {
    if tag.name == name {
      return Some(AlignSamTag::{
        name: tag.name,
        value: align_sam_copy_tag_value(tag.value),
      })
    }
  }
  None
}

///|
pub fn AlignSamAlignment::integer_tag(
  self : AlignSamAlignment,
  name : String,
) -> Int? {
  match self.tag(name) {
    Some(tag) =>
      match tag.value {
        SamTagInteger(value) => Some(value)
        _ => None
      }
    None => None
  }
}

///|
pub fn AlignSamAlignment::string_tag(
  self : AlignSamAlignment,
  name : String,
) -> String? {
  match self.tag(name) {
    Some(tag) =>
      match tag.value {
        SamTagCharacter(value) => Some(value)
        SamTagString(value) => Some(value)
        SamTagHex(value) => Some(value)
        _ => None
      }
    None => None
  }
}

///|
pub fn AlignSamAlignment::stats(self : AlignSamAlignment) -> AlignSamStats {
  let mut reference_consumed = 0
  let mut query_consumed = 0
  let mut aligned = 0
  let mut inserted = 0
  let mut deleted = 0
  let mut skipped = 0
  let mut soft_clipped = 0
  let mut hard_clipped = 0
  for element in self.cigar {
    let length = element.length
    if element.operation == "M" ||
      element.operation == "=" ||
      element.operation == "X" {
      reference_consumed = reference_consumed + length
      query_consumed = query_consumed + length
      aligned = aligned + length
    } else if element.operation == "I" {
      query_consumed = query_consumed + length
      inserted = inserted + length
    } else if element.operation == "D" {
      reference_consumed = reference_consumed + length
      deleted = deleted + length
    } else if element.operation == "N" {
      reference_consumed = reference_consumed + length
      skipped = skipped + length
    } else if element.operation == "S" {
      query_consumed = query_consumed + length
      soft_clipped = soft_clipped + length
    } else if element.operation == "H" {
      hard_clipped = hard_clipped + length
    }
  }
  AlignSamStats::{
    reference_consumed,
    query_consumed,
    aligned,
    inserted,
    deleted,
    skipped,
    soft_clipped,
    hard_clipped,
  }
}

///|
fn align_sam_is_aligned_operation(operation : AlignSamOperation) -> Bool {
  operation == SamAligned || operation == SamEqual || operation == SamMismatch
}

///|
pub fn AlignSamAlignment::target_to_query(
  self : AlignSamAlignment,
  position : Int,
) -> Int? {
  for index in 0..<self.operations.length() {
    if !align_sam_is_aligned_operation(self.operations[index]) {
      continue
    }
    let target_start = self.target_coordinates[index]
    let target_end = self.target_coordinates[index + 1]
    if position >= target_start && position < target_end {
      let query_start = self.query_coordinates[index]
      let query_end = self.query_coordinates[index + 1]
      let offset = position - target_start
      if query_end >= query_start {
        return Some(query_start + offset)
      } else {
        return Some(query_start - 1 - offset)
      }
    }
  }
  None
}

///|
pub fn AlignSamAlignment::query_to_target(
  self : AlignSamAlignment,
  position : Int,
) -> Int? {
  for index in 0..<self.operations.length() {
    if !align_sam_is_aligned_operation(self.operations[index]) {
      continue
    }
    let query_start = self.query_coordinates[index]
    let query_end = self.query_coordinates[index + 1]
    if query_end >= query_start {
      if position >= query_start && position < query_end {
        return Some(self.target_coordinates[index] + position - query_start)
      }
    } else if position >= query_end && position < query_start {
      return Some(self.target_coordinates[index] + query_start - 1 - position)
    }
  }
  None
}

///|
fn align_sam_raw_sequence(alignment : AlignSamAlignment) -> String {
  if !alignment.sequence_known {
    "*"
  } else if alignment.is_reverse() {
    align_sam_reverse_complement(alignment.query_sequence)
  } else {
    alignment.query_sequence
  }
}

///|
fn align_sam_md_positions(
  alignment : AlignSamAlignment,
) -> (Array[Int], Array[String]) {
  let positions : Array[Int] = []
  let defaults : Array[String] = []
  let raw_sequence = align_sam_raw_sequence(alignment)
  let mut target = match alignment.reference_start() {
    Some(value) => value
    None => 0
  }
  let mut query = 0
  for element in alignment.cigar {
    let operation = element.operation
    if operation == "S" {
      query = query + element.length
    } else if operation == "M" || operation == "=" || operation == "X" {
      for offset in 0..<element.length {
        positions.push(target + offset)
        if alignment.sequence_known {
          defaults.push(
            raw_sequence[query + offset:query + offset + 1].to_owned(),
          )
        } else {
          defaults.push("?")
        }
      }
      target = target + element.length
      query = query + element.length
    } else if operation == "I" {
      query = query + element.length
    } else if operation == "D" {
      for offset in 0..<element.length {
        positions.push(target + offset)
        defaults.push("?")
      }
      target = target + element.length
    } else if operation == "N" {
      target = target + element.length
    }
  }
  (positions, defaults)
}

///|
fn align_sam_reference_from_md(
  alignment : AlignSamAlignment,
) -> (Array[Int], Array[String]) raise AlignSamError {
  let (positions, values) = align_sam_md_positions(alignment)
  let md = match alignment.string_tag("MD") {
    Some(value) => value
    None => return (positions, values)
  }
  let mut pointer = 0
  let mut index = 0
  while index < md.length() {
    let code = md.unsafe_get(index).to_int()
    if code >= '0'.to_int() && code <= '9'.to_int() {
      let start = index
      while index < md.length() {
        let digit = md.unsafe_get(index).to_int()
        if digit < '0'.to_int() || digit > '9'.to_int() {
          break
        }
        index = index + 1
      }
      pointer = pointer +
        align_sam_parse_int(md[start:index].to_owned(), "MD match count")
    } else if code == '^'.to_int() {
      index = index + 1
      while index < md.length() &&
            align_sam_is_md_base(md.unsafe_get(index).to_int()) {
        if pointer >= values.length() {
          align_sam_fail("MD deletion exceeds the CIGAR path")
        }
        values[pointer] = md[index:index + 1].to_owned()
        pointer = pointer + 1
        index = index + 1
      }
    } else {
      if pointer >= values.length() {
        align_sam_fail("MD mismatch exceeds the CIGAR path")
      }
      values[pointer] = md[index:index + 1].to_owned()
      pointer = pointer + 1
      index = index + 1
    }
  }
  if pointer != values.length() {
    align_sam_fail("MD tag does not cover the complete CIGAR reference path")
  }
  (positions, values)
}

///|
fn align_sam_reference_char(
  positions : Array[Int],
  values : Array[String],
  position : Int,
) -> String {
  for index in 0..<positions.length() {
    if positions[index] == position {
      return values[index]
    }
  }
  "?"
}

///|
pub fn AlignSamAlignment::aligned_rows(
  self : AlignSamAlignment,
  reference_sequence? : String = "",
) -> (String, String) raise AlignSamError {
  if self.is_unmapped() {
    align_sam_fail("unmapped SAM records do not have aligned rows")
  }
  if reference_sequence.length() > 0 {
    match self.reference_length {
      Some(length) =>
        if reference_sequence.length() != length {
          align_sam_fail("reference sequence length does not match @SQ LN")
        }
      None => ()
    }
  }
  let (md_positions, md_values) = align_sam_reference_from_md(self)
  let raw_sequence = align_sam_raw_sequence(self)
  let target_row = StringBuilder::new()
  let query_row = StringBuilder::new()
  let mut target = match self.reference_start() {
    Some(value) => value
    None => 0
  }
  let mut query = 0
  for element in self.cigar {
    let operation = element.operation
    if operation == "S" {
      query = query + element.length
    } else if operation == "M" || operation == "=" || operation == "X" {
      for offset in 0..<element.length {
        let target_position = target + offset
        if reference_sequence.length() > 0 {
          target_row.write_string(
            reference_sequence[target_position:target_position + 1].to_owned(),
          )
        } else {
          target_row.write_string(
            align_sam_reference_char(md_positions, md_values, target_position),
          )
        }
        if self.sequence_known {
          query_row.write_string(
            raw_sequence[query + offset:query + offset + 1].to_owned(),
          )
        } else {
          query_row.write_char('?')
        }
      }
      target = target + element.length
      query = query + element.length
    } else if operation == "I" {
      target_row.write_string("-".repeat(element.length))
      if self.sequence_known {
        query_row.write_string(
          raw_sequence[query:query + element.length].to_owned(),
        )
      } else {
        query_row.write_string("?".repeat(element.length))
      }
      query = query + element.length
    } else if operation == "D" || operation == "N" {
      for offset in 0..<element.length {
        let target_position = target + offset
        if reference_sequence.length() > 0 {
          target_row.write_string(
            reference_sequence[target_position:target_position + 1].to_owned(),
          )
        } else {
          target_row.write_string(
            align_sam_reference_char(md_positions, md_values, target_position),
          )
        }
      }
      query_row.write_string("-".repeat(element.length))
      target = target + element.length
    }
  }
  (target_row.to_string(), query_row.to_string())
}

///|
pub fn AlignSamAlignment::computed_nm(self : AlignSamAlignment) -> Int? {
  if self.is_unmapped() {
    return None
  }
  let mut count = 0
  for element in self.cigar {
    if element.operation == "I" ||
      element.operation == "D" ||
      element.operation == "X" {
      count = count + element.length
    }
  }
  match self.string_tag("MD") {
    Some(md) => {
      let mut index = 0
      let mut in_deletion = false
      while index < md.length() {
        let code = md.unsafe_get(index).to_int()
        if code == '^'.to_int() {
          in_deletion = true
        } else if code >= '0'.to_int() && code <= '9'.to_int() {
          in_deletion = false
        } else if align_sam_is_md_base(code) && !in_deletion {
          count = count + 1
        }
        index = index + 1
      }
      // X operations and MD mismatch letters describe the same substitutions.
      let mut explicit_mismatches = 0
      for element in self.cigar {
        if element.operation == "X" {
          explicit_mismatches = explicit_mismatches + element.length
        }
      }
      Some(count - explicit_mismatches)
    }
    None => {
      for element in self.cigar {
        if element.operation == "M" {
          return None
        }
      }
      Some(count)
    }
  }
}

///|
pub fn align_sam_calculate_md(
  alignment : AlignSamAlignment,
  reference_sequence : String,
) -> String raise AlignSamError {
  if alignment.is_unmapped() {
    align_sam_fail("cannot calculate MD for an unmapped record")
  }
  match alignment.reference_length {
    Some(length) =>
      if reference_sequence.length() != length {
        align_sam_fail("reference sequence length does not match @SQ LN")
      }
    None => ()
  }
  let raw_sequence = align_sam_raw_sequence(alignment)
  if !alignment.sequence_known {
    align_sam_fail("cannot calculate MD without a query sequence")
  }
  let mut target = match alignment.reference_start() {
    Some(value) => value
    None => 0
  }
  let mut query = 0
  let mut matches = 0
  let output = StringBuilder::new()
  for element in alignment.cigar {
    let operation = element.operation
    if operation == "S" {
      query = query + element.length
    } else if operation == "M" || operation == "=" || operation == "X" {
      for offset in 0..<element.length {
        let target_base = reference_sequence[target + offset:target + offset + 1].to_owned()
        let query_base = raw_sequence[query + offset:query + offset + 1].to_owned()
        if target_base.to_upper() == query_base.to_upper() {
          matches = matches + 1
        } else {
          output.write_string(matches.to_string())
          output.write_string(target_base)
          matches = 0
        }
      }
      target = target + element.length
      query = query + element.length
    } else if operation == "I" {
      query = query + element.length
    } else if operation == "D" {
      output.write_string(matches.to_string())
      matches = 0
      output.write_char('^')
      output.write_string(
        reference_sequence[target:target + element.length].to_owned(),
      )
      target = target + element.length
    } else if operation == "N" {
      target = target + element.length
    }
  }
  output.write_string(matches.to_string())
  output.to_string()
}

///|
fn align_sam_quality_text(
  alignment : AlignSamAlignment,
) -> String raise AlignSamError {
  if !alignment.qualities_known {
    return "*"
  }
  let values = if alignment.is_reverse() {
    align_sam_reverse_ints(alignment.qualities)
  } else {
    align_sam_copy_ints(alignment.qualities)
  }
  let output = StringBuilder::new()
  for value in values {
    if value < 0 || value > 93 {
      align_sam_fail("PHRED quality must be between 0 and 93")
    }
    output.write_char((value + 33).unsafe_to_char())
  }
  output.to_string()
}

///|
pub fn align_sam_format(
  alignment : AlignSamAlignment,
  calculate_md? : Bool = false,
  reference_sequence? : String = "",
) -> String raise AlignSamError {
  let output = StringBuilder::new()
  output.write_string(alignment.query_name)
  output.write_char('\t')
  output.write_string(alignment.flag.to_string())
  output.write_char('\t')
  output.write_string(alignment.reference_name)
  output.write_char('\t')
  let position = match alignment.reference_start() {
    Some(value) => value + 1
    None => 0
  }
  output.write_string(position.to_string())
  output.write_char('\t')
  output.write_string(
    match alignment.mapq {
      Some(value) => value.to_string()
      None => "255"
    },
  )
  output.write_char('\t')
  output.write_string(align_sam_cigar_string(alignment.cigar))
  output.write_char('\t')
  let mate_reference = match alignment.mate_reference {
    Some(value) =>
      if value == alignment.reference_name && alignment.reference_name != "*" {
        "="
      } else {
        value
      }
    None => "*"
  }
  output.write_string(mate_reference)
  output.write_char('\t')
  output.write_string(
    match alignment.mate_position {
      Some(value) => (value + 1).to_string()
      None => "0"
    },
  )
  output.write_char('\t')
  output.write_string(
    match alignment.template_length {
      Some(value) => value.to_string()
      None => "0"
    },
  )
  output.write_char('\t')
  output.write_string(align_sam_raw_sequence(alignment))
  output.write_char('\t')
  output.write_string(align_sam_quality_text(alignment))
  for tag in alignment.tags {
    if calculate_md && tag.name == "MD" {
      continue
    }
    output.write_char('\t')
    output.write_string(align_sam_tag_text(tag))
  }
  if calculate_md {
    if reference_sequence.length() == 0 {
      align_sam_fail("reference_sequence is required when calculate_md is true")
    }
    output.write_string("\tMD:Z:")
    output.write_string(align_sam_calculate_md(alignment, reference_sequence))
  }
  output.write_char('\n')
  output.to_string()
}

///|
fn align_sam_header_text(header : AlignSamHeader) -> String {
  let output = StringBuilder::new()
  output.write_char('@')
  output.write_string(header.record_type)
  if header.record_type == "CO" {
    if header.comment.length() > 0 {
      output.write_char('\t')
      output.write_string(header.comment)
    }
  } else {
    for field in header.fields {
      output.write_char('\t')
      output.write_string(field.key)
      output.write_char(':')
      output.write_string(field.value)
    }
  }
  output.write_char('\n')
  output.to_string()
}

///|
pub fn align_sam_write(
  document : AlignSamDocument,
) -> String raise AlignSamError {
  let output = StringBuilder::new()
  for header in document.headers {
    output.write_string(align_sam_header_text(header))
  }
  for alignment in document.alignments {
    output.write_string(align_sam_format(alignment))
  }
  output.to_string()
}

///|
fn align_sam_format_tags(tags : Array[AlignSamTag]) -> String {
  let output = StringBuilder::new()
  for tag in tags {
    output.write_char('\t')
    output.write_string(align_sam_tag_text(tag))
  }
  output.to_string()
}

///|
fn align_sam_quality_from_values(
  qualities : Array[Int],
) -> String raise AlignSamError {
  let output = StringBuilder::new()
  for value in qualities {
    if value < 0 || value > 93 {
      align_sam_fail("PHRED quality must be between 0 and 93")
    }
    output.write_char((value + 33).unsafe_to_char())
  }
  output.to_string()
}

///|
pub fn align_sam_create(
  query_name : String,
  reference_name : String,
  reference_length : Int,
  target_start : Int,
  query_sequence : String,
  cigar : String,
  reverse? : Bool = false,
  mapq? : Int = 255,
  qualities? : Array[Int] = [],
  tags? : Array[AlignSamTag] = [],
  flag? : Int = 0,
) -> AlignSamAlignment raise AlignSamError {
  if reference_length <= 0 {
    align_sam_fail("reference_length must be positive")
  }
  if target_start < 0 {
    align_sam_fail("target_start must be non-negative")
  }
  if (flag & 0x4) != 0 {
    align_sam_fail("align_sam_create constructs mapped alignments")
  }
  let final_flag = if reverse { flag | 0x10 } else { flag & 0xffef }
  let stored_sequence = if reverse {
    align_sam_reverse_complement(query_sequence)
  } else {
    query_sequence
  }
  let stored_qualities = if reverse {
    align_sam_reverse_ints(qualities)
  } else {
    align_sam_copy_ints(qualities)
  }
  if qualities.length() > 0 && qualities.length() != query_sequence.length() {
    align_sam_fail("qualities length must equal query_sequence length")
  }
  let quality_text = if qualities.length() == 0 {
    "*"
  } else {
    align_sam_quality_from_values(stored_qualities)
  }
  let text = "@SQ\tSN:" +
    reference_name +
    "\tLN:" +
    reference_length.to_string() +
    "\n" +
    query_name +
    "\t" +
    final_flag.to_string() +
    "\t" +
    reference_name +
    "\t" +
    (target_start + 1).to_string() +
    "\t" +
    mapq.to_string() +
    "\t" +
    cigar +
    "\t*\t0\t0\t" +
    stored_sequence +
    "\t" +
    quality_text +
    align_sam_format_tags(tags) +
    "\n"
  let document = align_sam_parse(text)
  document.alignments[0]
}

///|
pub fn AlignSamAlignment::summary(self : AlignSamAlignment) -> String {
  if self.is_unmapped() {
    return "SAM alignment " + self.query_name + ": unmapped"
  }
  let stats = self.stats()
  let start = match self.reference_start() {
    Some(value) => value.to_string()
    None => "?"
  }
  let end = match self.reference_end() {
    Some(value) => value.to_string()
    None => "?"
  }
  "SAM alignment " +
  self.query_name +
  " -> " +
  self.reference_name +
  ":" +
  start +
  "-" +
  end +
  ", aligned=" +
  stats.aligned.to_string() +
  ", inserted=" +
  stats.inserted.to_string() +
  ", deleted=" +
  stats.deleted.to_string() +
  ", skipped=" +
  stats.skipped.to_string()
}

///|
pub fn align_sam_example_data() -> String {
  "@HD\tVN:1.6\tSO:coordinate\n" +
  "@SQ\tSN:chr1\tLN:40\tAS:demo\n" +
  "@RG\tID:rg1\tSM:sample1\n" +
  "read_forward\t0\tchr1\t5\t60\t3S4M1I2M1D3M2S\t*\t0\t0\tGGGACGTATCGTAAA\tIIIIIIIIIIIIIII\tNM:i:3\tMD:Z:6^G1A1\tAS:i:8\tRG:Z:rg1\n" +
  "read_reverse\t16\tchr1\t20\t40\t2H2S5M3N3M1S1H\t*\t0\t0\tTACGTACGTAA\tJJJJJJJJJJJ\tNM:i:0\tMD:Z:8\tAS:i:8\n" +
  "read_unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t!!!!\tRG:Z:rg1\n"
}
