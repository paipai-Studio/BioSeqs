// Biopython-compatible Bio.Align.chain support.
//
// This module is intentionally separate from chain_liftover.mbt. The older
// module provides a permissive rtracklayer-style liftOver API, while this file
// models modern Bio.Align pairwise alignments as absolute coordinate paths.

///|
/// Error raised for malformed chain data or invalid coordinate operations.
pub suberror AlignChainError {
  AlignChainError(String)
}

///|
/// One point in a two-row absolute coordinate path.
pub struct AlignChainCoordinate {
  target : Int
  query : Int
} derive(Eq, Debug)
pub extend AlignChainCoordinate with Eq::{equal, not_equal}
pub extend AlignChainCoordinate with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One canonical UCSC chain block.
///
/// A non-final row contains all three fields. The final row uses `size` only,
/// and therefore has zero target and query gaps.
pub struct AlignChainBlock {
  size : Int
  target_gap : Int
  query_gap : Int
} derive(Eq, Debug)
pub extend AlignChainBlock with Eq::{equal, not_equal}
pub extend AlignChainBlock with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One aligned part of a requested genomic range.
///
/// Coordinates are always zero-based, forward-axis, half-open intervals.
pub struct AlignChainRangeMapping {
  target_start : Int
  target_end : Int
  query_start : Int
  query_end : Int
  target_strand : String
  query_strand : String
} derive(Eq, Debug)
pub extend AlignChainRangeMapping with Eq::{equal, not_equal}
pub extend AlignChainRangeMapping with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One aligned residue pair in forward genomic coordinates.
pub struct AlignChainPair {
  target : Int
  query : Int
} derive(Eq, Debug)
pub extend AlignChainPair with Eq::{equal, not_equal}
pub extend AlignChainPair with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Coordinate-only chain alignment statistics.
pub struct AlignChainCounts {
  aligned : Int
  target_gap_bases : Int
  query_gap_bases : Int
  gap_columns : Int
  target_gap_opens : Int
  query_gap_opens : Int
  aligned_blocks : Int
  path_segments : Int
} derive(Eq, Debug)
pub extend AlignChainCounts with Eq::{equal, not_equal}
pub extend AlignChainCounts with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// A coordinate-aware pairwise chain alignment.
///
/// `coordinates` follows Biopython's `Alignment.coordinates` convention.
/// Reverse strands are represented by decreasing coordinates. Chain files do
/// not store sequence letters, only sequence lengths and the coordinate path.
pub struct AlignChainAlignment {
  target_id : String
  target_size : Int
  query_id : String
  query_size : Int
  score : Double
  chain_id : String
  coordinates : Array[AlignChainCoordinate]
} derive(Eq, Debug)
pub extend AlignChainAlignment with Eq::{equal, not_equal}
pub extend AlignChainAlignment with @moonbitlang/core/debug.Debug::{to_repr}

///|
fn align_chain_fail(message : String) -> Unit raise AlignChainError {
  raise AlignChainError(message)
}

///|
fn align_chain_abs(value : Int) -> Int {
  if value < 0 {
    -value
  } else {
    value
  }
}

///|
fn align_chain_min(left : Int, right : Int) -> Int {
  if left < right {
    left
  } else {
    right
  }
}

///|
fn align_chain_max(left : Int, right : Int) -> Int {
  if left > right {
    left
  } else {
    right
  }
}

///|
fn align_chain_sign(value : Int) -> Int {
  if value < 0 {
    -1
  } else if value > 0 {
    1
  } else {
    0
  }
}

///|
fn align_chain_copy_coordinates(
  coordinates : Array[AlignChainCoordinate],
) -> Array[AlignChainCoordinate] {
  let copied : Array[AlignChainCoordinate] = []
  for coordinate in coordinates {
    copied.push(coordinate)
  }
  copied
}

///|
fn align_chain_copy_alignments(
  alignments : Array[AlignChainAlignment],
) -> Array[AlignChainAlignment] {
  let copied : Array[AlignChainAlignment] = []
  for alignment in alignments {
    copied.push(alignment)
  }
  copied
}

///|
fn align_chain_validate_token(
  value : String,
  label : String,
  allow_empty : Bool,
) -> Unit raise AlignChainError {
  if !allow_empty && value.length() == 0 {
    align_chain_fail(label + " must not be empty")
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code <= 32 || code == 127 {
      align_chain_fail(
        label + " must not contain whitespace or control characters",
      )
    }
  }
}

///|
fn align_chain_validate_strand(
  strand : String,
  label : String,
) -> Unit raise AlignChainError {
  if strand != "+" && strand != "-" {
    align_chain_fail(label + " strand must be '+' or '-'")
  }
}

///|
fn align_chain_validate_interval(
  start : Int,
  end : Int,
  size : Int,
  label : String,
) -> Unit raise AlignChainError {
  if size < 0 {
    align_chain_fail(label + " size must be non-negative")
  }
  if start < 0 || end < start || end > size {
    align_chain_fail(label + " interval must satisfy 0 <= start <= end <= size")
  }
}

///|
fn align_chain_validate_alignment(
  alignment : AlignChainAlignment,
) -> Unit raise AlignChainError {
  align_chain_validate_token(alignment.target_id, "target identifier", false)
  align_chain_validate_token(alignment.query_id, "query identifier", false)
  align_chain_validate_token(alignment.chain_id, "chain identifier", true)
  if alignment.target_size < 0 {
    align_chain_fail("target size must be non-negative")
  }
  if alignment.query_size < 0 {
    align_chain_fail("query size must be non-negative")
  }
  if alignment.score.is_nan() || alignment.score.abs() > 1.0e300 {
    align_chain_fail("chain score must be finite")
  }
  if alignment.coordinates.length() < 2 {
    align_chain_fail("chain coordinate path must contain at least two points")
  }
  let mut target_direction = 0
  let mut query_direction = 0
  let mut aligned = 0
  for index in 0..<alignment.coordinates.length() {
    let point = alignment.coordinates[index]
    if point.target < 0 || point.target > alignment.target_size {
      align_chain_fail("target coordinate is outside the declared sequence")
    }
    if point.query < 0 || point.query > alignment.query_size {
      align_chain_fail("query coordinate is outside the declared sequence")
    }
    if index + 1 < alignment.coordinates.length() {
      let next = alignment.coordinates[index + 1]
      let target_step = next.target - point.target
      let query_step = next.query - point.query
      if target_step == 0 && query_step == 0 {
        align_chain_fail("chain coordinate path contains a zero-length segment")
      }
      if target_step != 0 {
        let direction = align_chain_sign(target_step)
        if target_direction == 0 {
          target_direction = direction
        } else if target_direction != direction {
          align_chain_fail("target coordinates must be monotonic")
        }
      }
      if query_step != 0 {
        let direction = align_chain_sign(query_step)
        if query_direction == 0 {
          query_direction = direction
        } else if query_direction != direction {
          align_chain_fail("query coordinates must be monotonic")
        }
      }
      if target_step != 0 && query_step != 0 {
        if align_chain_abs(target_step) != align_chain_abs(query_step) {
          align_chain_fail(
            "aligned target and query steps must have equal lengths",
          )
        }
        aligned = aligned + align_chain_abs(target_step)
      }
    }
  }
  if target_direction == 0 || query_direction == 0 || aligned == 0 {
    align_chain_fail("chain path must contain at least one aligned segment")
  }
}

///|
/// Construct one absolute coordinate.
pub fn AlignChainCoordinate::create(
  target : Int,
  query : Int,
) -> AlignChainCoordinate raise AlignChainError {
  if target < 0 || query < 0 {
    align_chain_fail("chain coordinates must be non-negative")
  }
  AlignChainCoordinate::{ target, query, }
}

///|
/// Construct one chain block.
pub fn AlignChainBlock::create(
  size : Int,
  target_gap? : Int = 0,
  query_gap? : Int = 0,
) -> AlignChainBlock raise AlignChainError {
  if size < 0 || target_gap < 0 || query_gap < 0 {
    align_chain_fail("chain block values must be non-negative")
  }
  AlignChainBlock::{ size, target_gap, query_gap, }
}

///|
/// Construct and validate a coordinate-aware chain alignment.
pub fn AlignChainAlignment::create(
  target_id : String,
  target_size : Int,
  query_id : String,
  query_size : Int,
  coordinates : Array[AlignChainCoordinate],
  score? : Double = 0.0,
  chain_id? : String = "",
) -> AlignChainAlignment raise AlignChainError {
  let alignment = AlignChainAlignment::{
    target_id,
    target_size,
    query_id,
    query_size,
    score,
    chain_id,
    coordinates: align_chain_copy_coordinates(coordinates),
  }
  align_chain_validate_alignment(alignment)
  alignment
}

///|
fn align_chain_oriented_coordinate(
  coordinate : Int,
  size : Int,
  strand : String,
) -> Int {
  if strand == "+" {
    coordinate
  } else {
    size - coordinate
  }
}

///|
/// Build an alignment from normalized UCSC header coordinates and blocks.
pub fn align_chain_from_blocks(
  target_id : String,
  target_size : Int,
  target_strand : String,
  target_start : Int,
  target_end : Int,
  query_id : String,
  query_size : Int,
  query_strand : String,
  query_start : Int,
  query_end : Int,
  blocks : Array[AlignChainBlock],
  score? : Double = 0.0,
  chain_id? : String = "",
) -> AlignChainAlignment raise AlignChainError {
  align_chain_validate_token(target_id, "target identifier", false)
  align_chain_validate_token(query_id, "query identifier", false)
  align_chain_validate_token(chain_id, "chain identifier", true)
  align_chain_validate_strand(target_strand, "target")
  align_chain_validate_strand(query_strand, "query")
  align_chain_validate_interval(target_start, target_end, target_size, "target")
  align_chain_validate_interval(query_start, query_end, query_size, "query")
  if score.is_nan() || score.abs() > 1.0e300 {
    align_chain_fail("chain score must be finite")
  }
  if blocks.length() == 0 {
    align_chain_fail("chain record must contain at least one block row")
  }
  let relative : Array[AlignChainCoordinate] = [
    AlignChainCoordinate::{ target: 0, query: 0, },
  ]
  let mut target_position = 0
  let mut query_position = 0
  let mut aligned = 0
  for index in 0..<blocks.length() {
    let block = blocks[index]
    if block.size < 0 || block.target_gap < 0 || block.query_gap < 0 {
      align_chain_fail("chain block values must be non-negative")
    }
    let is_final = index + 1 == blocks.length()
    if is_final {
      if block.target_gap != 0 || block.query_gap != 0 {
        align_chain_fail("final chain block must contain size only")
      }
    } else if block.target_gap == 0 && block.query_gap == 0 {
      align_chain_fail(
        "non-final chain block must contain a target or query gap",
      )
    }
    if block.size > 0 {
      target_position = target_position + block.size
      query_position = query_position + block.size
      relative.push(AlignChainCoordinate::{
        target: target_position,
        query: query_position,
      })
      aligned = aligned + block.size
    }
    if block.target_gap > 0 {
      target_position = target_position + block.target_gap
      relative.push(AlignChainCoordinate::{
        target: target_position,
        query: query_position,
      })
    }
    if block.query_gap > 0 {
      query_position = query_position + block.query_gap
      relative.push(AlignChainCoordinate::{
        target: target_position,
        query: query_position,
      })
    }
    if target_position > target_end - target_start {
      align_chain_fail("chain blocks exceed the declared target span")
    }
    if query_position > query_end - query_start {
      align_chain_fail("chain blocks exceed the declared query span")
    }
  }
  if aligned == 0 {
    align_chain_fail("chain record must contain aligned bases")
  }
  if target_position != target_end - target_start {
    align_chain_fail("chain blocks do not match the declared target span")
  }
  if query_position != query_end - query_start {
    align_chain_fail("chain blocks do not match the declared query span")
  }
  let coordinates : Array[AlignChainCoordinate] = []
  for point in relative {
    let target_oriented = target_start + point.target
    let query_oriented = query_start + point.query
    coordinates.push(AlignChainCoordinate::{
      target: align_chain_oriented_coordinate(
        target_oriented, target_size, target_strand,
      ),
      query: align_chain_oriented_coordinate(
        query_oriented, query_size, query_strand,
      ),
    })
  }
  AlignChainAlignment::create(
    target_id,
    target_size,
    query_id,
    query_size,
    coordinates,
    score~,
    chain_id~,
  )
}

///|
fn align_chain_strip_cr(value : String) -> String {
  if value.length() > 0 &&
    value.unsafe_get(value.length() - 1).to_int() == '\r'.to_int() {
    value[0:value.length() - 1].to_owned()
  } else {
    value
  }
}

///|
fn align_chain_split_lines(value : String) -> Array[String] {
  let lines : Array[String] = []
  let mut start = 0
  for index in 0..<value.length() {
    if value.unsafe_get(index).to_int() == '\n'.to_int() {
      lines.push(align_chain_strip_cr(value[start:index].to_owned()))
      start = index + 1
    }
  }
  if start < value.length() {
    lines.push(align_chain_strip_cr(value[start:value.length()].to_owned()))
  } else if value.length() > 0 {
    lines.push("")
  }
  lines
}

///|
fn align_chain_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\r'.to_int() ||
  code == '\n'.to_int()
}

///|
fn align_chain_trim(value : String) -> String {
  let mut start = 0
  while start < value.length() &&
        align_chain_is_whitespace(value.unsafe_get(start).to_int()) {
    start = start + 1
  }
  let mut end = value.length()
  while end > start &&
        align_chain_is_whitespace(value.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  value[start:end].to_owned()
}

///|
fn align_chain_split_whitespace(value : String) -> Array[String] {
  let words : Array[String] = []
  let mut index = 0
  while index < value.length() {
    while index < value.length() &&
          align_chain_is_whitespace(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index >= value.length() {
      break
    }
    let start = index
    while index < value.length() &&
          !align_chain_is_whitespace(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    words.push(value[start:index].to_owned())
  }
  words
}

///|
fn align_chain_parse_uint(
  value : String,
  label : String,
) -> Int raise AlignChainError {
  if value.length() == 0 {
    align_chain_fail(label + " must be an unsigned integer")
  }
  let mut result = 0
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      align_chain_fail(label + " must be an unsigned integer")
    }
    let digit = code - '0'.to_int()
    if result > (2147483647 - digit) / 10 {
      align_chain_fail(label + " exceeds the supported integer range")
    }
    result = result * 10 + digit
  }
  result
}

///|
fn align_chain_is_score(value : String) -> Bool {
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
fn align_chain_parse_score(value : String) -> Double raise AlignChainError {
  if !align_chain_is_score(value) {
    align_chain_fail("chain score must be a number")
  }
  let parse_value = if value.length() > 0 &&
    value.unsafe_get(0).to_int() == '+'.to_int() {
    value[1:value.length()].to_owned()
  } else {
    value
  }
  let number = match parse_double(parse_value) {
    Some(number) => number
    None => {
      align_chain_fail("chain score must be a number")
      0.0
    }
  }
  if number.is_nan() || number.abs() > 1.0e300 {
    align_chain_fail("chain score must be finite")
  }
  number
}

///|
fn align_chain_parse_record(
  lines : Array[String],
  start : Int,
) -> (AlignChainAlignment, Int) raise AlignChainError {
  let header = align_chain_split_whitespace(align_chain_trim(lines[start]))
  if header.length() != 12 && header.length() != 13 {
    align_chain_fail(
      "chain header must contain exactly 12 or 13 whitespace-delimited fields",
    )
  }
  if header[0] != "chain" {
    align_chain_fail("chain record must start with the 'chain' keyword")
  }
  let score = align_chain_parse_score(header[1])
  let target_id = header[2]
  let target_size = align_chain_parse_uint(header[3], "target size")
  let target_strand = header[4]
  let target_start = align_chain_parse_uint(header[5], "target start")
  let target_end = align_chain_parse_uint(header[6], "target end")
  let query_id = header[7]
  let query_size = align_chain_parse_uint(header[8], "query size")
  let query_strand = header[9]
  let query_start = align_chain_parse_uint(header[10], "query start")
  let query_end = align_chain_parse_uint(header[11], "query end")
  let chain_id = if header.length() == 13 { header[12] } else { "" }
  align_chain_validate_strand(target_strand, "target")
  align_chain_validate_strand(query_strand, "query")
  align_chain_validate_interval(target_start, target_end, target_size, "target")
  align_chain_validate_interval(query_start, query_end, query_size, "query")
  let blocks : Array[AlignChainBlock] = []
  let mut index = start + 1
  let mut terminated = false
  while index < lines.length() {
    let line = align_chain_trim(lines[index])
    if line.length() == 0 {
      align_chain_fail("blank line found before the final chain block")
    }
    let words = align_chain_split_whitespace(line)
    if words.length() > 0 && words[0] == "chain" {
      align_chain_fail("chain record is missing its final size-only block")
    }
    if words.length() != 1 && words.length() != 3 {
      align_chain_fail("chain block must contain one or three integer fields")
    }
    let size = align_chain_parse_uint(words[0], "chain block size")
    if words.length() == 1 {
      blocks.push(AlignChainBlock::{ size, target_gap: 0, query_gap: 0, })
      index = index + 1
      terminated = true
      break
    }
    let target_gap = align_chain_parse_uint(words[1], "target gap")
    let query_gap = align_chain_parse_uint(words[2], "query gap")
    if target_gap == 0 && query_gap == 0 {
      align_chain_fail("non-final chain block must contain a non-zero gap")
    }
    blocks.push(AlignChainBlock::{ size, target_gap, query_gap, })
    index = index + 1
  }
  if !terminated {
    align_chain_fail("chain record is missing its final size-only block")
  }
  let alignment = align_chain_from_blocks(
    target_id,
    target_size,
    target_strand,
    target_start,
    target_end,
    query_id,
    query_size,
    query_strand,
    query_start,
    query_end,
    blocks,
    score~,
    chain_id~,
  )
  (alignment, index)
}

///|
/// Parse all chain records in a string.
pub fn align_chain_parse_all(
  text : String,
) -> Array[AlignChainAlignment] raise AlignChainError {
  let lines = align_chain_split_lines(text)
  let alignments : Array[AlignChainAlignment] = []
  let mut index = 0
  while index < lines.length() {
    let line = align_chain_trim(lines[index])
    if line.length() == 0 {
      index = index + 1
      continue
    }
    if align_chain_split_whitespace(line)[0] != "chain" {
      align_chain_fail(
        "unexpected content before chain record at line " +
        (index + 1).to_string(),
      )
    }
    let (alignment, next) = align_chain_parse_record(lines, index)
    alignments.push(alignment)
    index = next
  }
  alignments
}

///|
/// Parse exactly one chain record.
pub fn align_chain_parse(
  text : String,
) -> AlignChainAlignment raise AlignChainError {
  let alignments = align_chain_parse_all(text)
  if alignments.length() == 0 {
    align_chain_fail("chain input contains no records")
  }
  if alignments.length() != 1 {
    align_chain_fail("expected exactly one chain record")
  }
  alignments[0]
}

///|
/// Return the target strand encoded by the coordinate path.
pub fn AlignChainAlignment::target_strand(self : AlignChainAlignment) -> String {
  if self.coordinates[self.coordinates.length() - 1].target >
    self.coordinates[0].target {
    "+"
  } else {
    "-"
  }
}

///|
/// Return the query strand encoded by the coordinate path.
pub fn AlignChainAlignment::query_strand(self : AlignChainAlignment) -> String {
  if self.coordinates[self.coordinates.length() - 1].query >
    self.coordinates[0].query {
    "+"
  } else {
    "-"
  }
}

///|
/// Return the normalized target start.
pub fn AlignChainAlignment::target_start(self : AlignChainAlignment) -> Int {
  align_chain_min(
    self.coordinates[0].target,
    self.coordinates[self.coordinates.length() - 1].target,
  )
}

///|
/// Return the normalized target end.
pub fn AlignChainAlignment::target_end(self : AlignChainAlignment) -> Int {
  align_chain_max(
    self.coordinates[0].target,
    self.coordinates[self.coordinates.length() - 1].target,
  )
}

///|
/// Return the normalized query start.
pub fn AlignChainAlignment::query_start(self : AlignChainAlignment) -> Int {
  align_chain_min(
    self.coordinates[0].query,
    self.coordinates[self.coordinates.length() - 1].query,
  )
}

///|
/// Return the normalized query end.
pub fn AlignChainAlignment::query_end(self : AlignChainAlignment) -> Int {
  align_chain_max(
    self.coordinates[0].query,
    self.coordinates[self.coordinates.length() - 1].query,
  )
}

///|
fn align_chain_target_oriented(
  alignment : AlignChainAlignment,
  coordinate : Int,
) -> Int {
  if alignment.target_strand() == "+" {
    coordinate
  } else {
    alignment.target_size - coordinate
  }
}

///|
fn align_chain_query_oriented(
  alignment : AlignChainAlignment,
  coordinate : Int,
) -> Int {
  if alignment.query_strand() == "+" {
    coordinate
  } else {
    alignment.query_size - coordinate
  }
}

///|
/// Convert the coordinate path to canonical chain block rows.
pub fn AlignChainAlignment::blocks(
  self : AlignChainAlignment,
) -> Array[AlignChainBlock] raise AlignChainError {
  align_chain_validate_alignment(self)
  let blocks : Array[AlignChainBlock] = []
  let mut size = 0
  let mut target_gap = 0
  let mut query_gap = 0
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    let target_step = align_chain_target_oriented(self, next.target) -
      align_chain_target_oriented(self, current.target)
    let query_step = align_chain_query_oriented(self, next.query) -
      align_chain_query_oriented(self, current.query)
    if target_step < 0 || query_step < 0 {
      align_chain_fail("chain coordinate path is not monotonic on its strand")
    }
    if target_step > 0 && query_step > 0 {
      if target_step != query_step {
        align_chain_fail(
          "aligned target and query steps must have equal lengths",
        )
      }
      if target_gap > 0 || query_gap > 0 {
        blocks.push(AlignChainBlock::{ size, target_gap, query_gap, })
        size = target_step
        target_gap = 0
        query_gap = 0
      } else {
        size = size + target_step
      }
    } else if target_step > 0 {
      target_gap = target_gap + target_step
    } else if query_step > 0 {
      query_gap = query_gap + query_step
    }
  }
  if target_gap > 0 || query_gap > 0 {
    blocks.push(AlignChainBlock::{ size, target_gap, query_gap, })
    size = 0
  }
  blocks.push(AlignChainBlock::{ size, target_gap: 0, query_gap: 0, })
  blocks
}

///|
fn align_chain_header_target_start(alignment : AlignChainAlignment) -> Int {
  align_chain_target_oriented(alignment, alignment.coordinates[0].target)
}

///|
fn align_chain_header_target_end(alignment : AlignChainAlignment) -> Int {
  align_chain_target_oriented(
    alignment,
    alignment.coordinates[alignment.coordinates.length() - 1].target,
  )
}

///|
fn align_chain_header_query_start(alignment : AlignChainAlignment) -> Int {
  align_chain_query_oriented(alignment, alignment.coordinates[0].query)
}

///|
fn align_chain_header_query_end(alignment : AlignChainAlignment) -> Int {
  align_chain_query_oriented(
    alignment,
    alignment.coordinates[alignment.coordinates.length() - 1].query,
  )
}

///|
/// Serialize one alignment as a canonical UCSC chain record.
pub fn align_chain_write(
  alignment : AlignChainAlignment,
) -> String raise AlignChainError {
  align_chain_validate_alignment(alignment)
  let output = StringBuilder::new()
  output.write_string("chain ")
  output.write_string(alignment.score.to_string())
  output.write_char(' ')
  output.write_string(alignment.target_id)
  output.write_char(' ')
  output.write_string(alignment.target_size.to_string())
  output.write_char(' ')
  output.write_string(alignment.target_strand())
  output.write_char(' ')
  output.write_string(align_chain_header_target_start(alignment).to_string())
  output.write_char(' ')
  output.write_string(align_chain_header_target_end(alignment).to_string())
  output.write_char(' ')
  output.write_string(alignment.query_id)
  output.write_char(' ')
  output.write_string(alignment.query_size.to_string())
  output.write_char(' ')
  output.write_string(alignment.query_strand())
  output.write_char(' ')
  output.write_string(align_chain_header_query_start(alignment).to_string())
  output.write_char(' ')
  output.write_string(align_chain_header_query_end(alignment).to_string())
  if alignment.chain_id.length() > 0 {
    output.write_char(' ')
    output.write_string(alignment.chain_id)
  }
  output.write_char('\n')
  let blocks = alignment.blocks()
  for index in 0..<blocks.length() {
    let block = blocks[index]
    output.write_string(block.size.to_string())
    if index + 1 < blocks.length() {
      output.write_char('\t')
      output.write_string(block.target_gap.to_string())
      output.write_char('\t')
      output.write_string(block.query_gap.to_string())
    }
    output.write_char('\n')
  }
  output.write_char('\n')
  output.to_string()
}

///|
/// Serialize multiple alignments as consecutive chain records.
pub fn align_chain_write_all(
  alignments : Array[AlignChainAlignment],
) -> String raise AlignChainError {
  let output = StringBuilder::new()
  for alignment in alignments {
    output.write_string(align_chain_write(alignment))
  }
  output.to_string()
}

///|
/// Return one operation per coordinate-path segment.
///
/// `M` consumes both rows, `D` consumes target only, and `I` consumes query
/// only.
pub fn AlignChainAlignment::operation_path(
  self : AlignChainAlignment,
) -> String {
  let output = StringBuilder::new()
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    let target_step = next.target - current.target
    let query_step = next.query - current.query
    if target_step != 0 && query_step != 0 {
      output.write_char('M')
    } else if target_step != 0 {
      output.write_char('D')
    } else {
      output.write_char('I')
    }
  }
  output.to_string()
}

///|
/// Count the number of alignment columns represented by the path.
pub fn AlignChainAlignment::alignment_length(self : AlignChainAlignment) -> Int {
  let mut length = 0
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    length = length +
      align_chain_max(
        align_chain_abs(next.target - current.target),
        align_chain_abs(next.query - current.query),
      )
  }
  length
}

///|
/// Count aligned bases.
pub fn AlignChainAlignment::aligned_bases(self : AlignChainAlignment) -> Int {
  let mut aligned = 0
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    let target_step = align_chain_abs(next.target - current.target)
    let query_step = align_chain_abs(next.query - current.query)
    if target_step > 0 && query_step > 0 {
      aligned = aligned + target_step
    }
  }
  aligned
}

///|
/// Compute coordinate-only alignment counts.
pub fn AlignChainAlignment::counts(
  self : AlignChainAlignment,
) -> AlignChainCounts {
  let mut aligned = 0
  let mut target_gap_bases = 0
  let mut query_gap_bases = 0
  let mut target_gap_opens = 0
  let mut query_gap_opens = 0
  let mut aligned_blocks = 0
  let mut previous = 'X'
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    let target_step = align_chain_abs(next.target - current.target)
    let query_step = align_chain_abs(next.query - current.query)
    if target_step > 0 && query_step > 0 {
      aligned = aligned + target_step
      if previous != 'M' {
        aligned_blocks = aligned_blocks + 1
      }
      previous = 'M'
    } else if target_step > 0 {
      target_gap_bases = target_gap_bases + target_step
      if previous != 'D' {
        target_gap_opens = target_gap_opens + 1
      }
      previous = 'D'
    } else {
      query_gap_bases = query_gap_bases + query_step
      if previous != 'I' {
        query_gap_opens = query_gap_opens + 1
      }
      previous = 'I'
    }
  }
  AlignChainCounts::{
    aligned,
    target_gap_bases,
    query_gap_bases,
    gap_columns: target_gap_bases + query_gap_bases,
    target_gap_opens,
    query_gap_opens,
    aligned_blocks,
    path_segments: self.coordinates.length() - 1,
  }
}

///|
fn align_chain_map_on_segment(
  source : Int,
  source_start : Int,
  source_end : Int,
  destination_start : Int,
  destination_end : Int,
) -> Int? {
  let source_step = source_end - source_start
  let destination_step = destination_end - destination_start
  if source_step == 0 ||
    destination_step == 0 ||
    align_chain_abs(source_step) != align_chain_abs(destination_step) {
    return None
  }
  let offset = if source_step > 0 {
    if source < source_start || source >= source_end {
      return None
    }
    source - source_start
  } else {
    if source < source_end || source >= source_start {
      return None
    }
    source_start - 1 - source
  }
  if destination_step > 0 {
    Some(destination_start + offset)
  } else {
    Some(destination_start - 1 - offset)
  }
}

///|
/// Map one target residue to the query, returning `None` in an unaligned gap.
pub fn AlignChainAlignment::map_target_position(
  self : AlignChainAlignment,
  position : Int,
) -> Int? {
  if position < 0 || position >= self.target_size {
    return None
  }
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    match
      align_chain_map_on_segment(
        position,
        current.target,
        next.target,
        current.query,
        next.query,
      ) {
      Some(mapped) => return Some(mapped)
      None => ()
    }
  }
  None
}

///|
/// Map one query residue to the target, returning `None` in an unaligned gap.
pub fn AlignChainAlignment::map_query_position(
  self : AlignChainAlignment,
  position : Int,
) -> Int? {
  if position < 0 || position >= self.query_size {
    return None
  }
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    match
      align_chain_map_on_segment(
        position,
        current.query,
        next.query,
        current.target,
        next.target,
      ) {
      Some(mapped) => return Some(mapped)
      None => ()
    }
  }
  None
}

///|
fn align_chain_range_mapping(
  alignment : AlignChainAlignment,
  target_start : Int,
  target_end : Int,
  query_start : Int,
  query_end : Int,
) -> AlignChainRangeMapping {
  AlignChainRangeMapping::{
    target_start,
    target_end,
    query_start,
    query_end,
    target_strand: alignment.target_strand(),
    query_strand: alignment.query_strand(),
  }
}

///|
/// Map a target interval to all aligned query pieces.
///
/// Target-only gap portions are omitted, so one input interval may yield
/// multiple output pieces.
pub fn AlignChainAlignment::map_target_range(
  self : AlignChainAlignment,
  start : Int,
  end : Int,
) -> Array[AlignChainRangeMapping] raise AlignChainError {
  if start < 0 || end <= start || end > self.target_size {
    align_chain_fail("target range must satisfy 0 <= start < end <= size")
  }
  let mappings : Array[AlignChainRangeMapping] = []
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    let target_low = align_chain_min(current.target, next.target)
    let target_high = align_chain_max(current.target, next.target)
    if current.target != next.target && current.query != next.query {
      let piece_start = align_chain_max(start, target_low)
      let piece_end = align_chain_min(end, target_high)
      if piece_start < piece_end {
        let query_first = match self.map_target_position(piece_start) {
          Some(value) => value
          None => {
            align_chain_fail("internal target range mapping failure")
            0
          }
        }
        let query_last = match self.map_target_position(piece_end - 1) {
          Some(value) => value
          None => {
            align_chain_fail("internal target range mapping failure")
            0
          }
        }
        mappings.push(
          align_chain_range_mapping(
            self,
            piece_start,
            piece_end,
            align_chain_min(query_first, query_last),
            align_chain_max(query_first, query_last) + 1,
          ),
        )
      }
    }
  }
  mappings
}

///|
/// Map a query interval to all aligned target pieces.
pub fn AlignChainAlignment::map_query_range(
  self : AlignChainAlignment,
  start : Int,
  end : Int,
) -> Array[AlignChainRangeMapping] raise AlignChainError {
  if start < 0 || end <= start || end > self.query_size {
    align_chain_fail("query range must satisfy 0 <= start < end <= size")
  }
  let mappings : Array[AlignChainRangeMapping] = []
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    let query_low = align_chain_min(current.query, next.query)
    let query_high = align_chain_max(current.query, next.query)
    if current.target != next.target && current.query != next.query {
      let piece_start = align_chain_max(start, query_low)
      let piece_end = align_chain_min(end, query_high)
      if piece_start < piece_end {
        let target_first = match self.map_query_position(piece_start) {
          Some(value) => value
          None => {
            align_chain_fail("internal query range mapping failure")
            0
          }
        }
        let target_last = match self.map_query_position(piece_end - 1) {
          Some(value) => value
          None => {
            align_chain_fail("internal query range mapping failure")
            0
          }
        }
        mappings.push(
          align_chain_range_mapping(
            self,
            align_chain_min(target_first, target_last),
            align_chain_max(target_first, target_last) + 1,
            piece_start,
            piece_end,
          ),
        )
      }
    }
  }
  mappings
}

///|
/// Materialize aligned residue pairs, subject to a caller-controlled limit.
pub fn AlignChainAlignment::aligned_pairs(
  self : AlignChainAlignment,
  limit? : Int = 1000000,
) -> Array[AlignChainPair] raise AlignChainError {
  if limit < 0 {
    align_chain_fail("aligned-pair limit must be non-negative")
  }
  let aligned = self.aligned_bases()
  if aligned > limit {
    align_chain_fail("aligned-pair count exceeds the requested limit")
  }
  let pairs : Array[AlignChainPair] = []
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    let target_step = next.target - current.target
    let query_step = next.query - current.query
    if target_step != 0 && query_step != 0 {
      let length = align_chain_abs(target_step)
      let target_direction = align_chain_sign(target_step)
      let query_direction = align_chain_sign(query_step)
      for offset in 0..<length {
        pairs.push(AlignChainPair::{
          target: if target_direction > 0 {
            current.target + offset
          } else {
            current.target - 1 - offset
          },
          query: if query_direction > 0 {
            current.query + offset
          } else {
            current.query - 1 - offset
          },
        })
      }
    }
  }
  pairs
}

///|
/// Swap target and query while preserving score, ID, and path orientation.
pub fn AlignChainAlignment::invert(
  self : AlignChainAlignment,
) -> AlignChainAlignment raise AlignChainError {
  let coordinates : Array[AlignChainCoordinate] = []
  for point in self.coordinates {
    coordinates.push(AlignChainCoordinate::{
      target: point.query,
      query: point.target,
    })
  }
  AlignChainAlignment::create(
    self.query_id,
    self.query_size,
    self.target_id,
    self.target_size,
    coordinates,
    score=self.score,
    chain_id=self.chain_id,
  )
}

///|
/// Find the first alignment with the requested chain ID.
pub fn align_chain_find_by_id(
  alignments : Array[AlignChainAlignment],
  chain_id : String,
) -> AlignChainAlignment? {
  for alignment in alignments {
    if alignment.chain_id == chain_id {
      return Some(alignment)
    }
  }
  None
}

///|
/// Return alignments whose target span overlaps a half-open interval.
pub fn align_chain_target_overlaps(
  alignments : Array[AlignChainAlignment],
  target_id : String,
  start : Int,
  end : Int,
) -> Array[AlignChainAlignment] raise AlignChainError {
  if start < 0 || end <= start {
    align_chain_fail("overlap query must satisfy 0 <= start < end")
  }
  let matches : Array[AlignChainAlignment] = []
  for alignment in alignments {
    if alignment.target_id == target_id &&
      alignment.target_start() < end &&
      alignment.target_end() > start {
      matches.push(alignment)
    }
  }
  align_chain_copy_alignments(matches)
}

///|
/// Return a concise alignment summary.
pub fn AlignChainAlignment::summary(self : AlignChainAlignment) -> String {
  let counts = self.counts()
  "AlignChainAlignment(target=" +
  self.target_id +
  ":" +
  self.target_start().to_string() +
  "-" +
  self.target_end().to_string() +
  self.target_strand() +
  ", query=" +
  self.query_id +
  ":" +
  self.query_start().to_string() +
  "-" +
  self.query_end().to_string() +
  self.query_strand() +
  ", aligned=" +
  counts.aligned.to_string() +
  ", gaps=" +
  counts.gap_columns.to_string() +
  ")"
}

///|
/// Return a compact two-record fixture based on Biopython's chain examples.
pub fn align_chain_example_text() -> String {
  "chain 176 chr3 198295559 + 42530895 42532606 NR_046654.1 181 - 0 181 1\n" +
  "63\t1062\t0\n" +
  "75\t468\t0\n" +
  "43\n\n" +
  "chain 500 chr1 1000 - 100 200 queryA 500 + 20 130 reverse-target\n" +
  "40\t10\t20\n" +
  "50\n"
}
