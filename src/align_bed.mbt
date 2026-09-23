// Biopython-compatible Bio.Align.bed support.
//
// This module models standalone BED records as pairwise alignment coordinate
// paths. BigBed binary storage remains implemented separately in bigbed.mbt.

///|
/// Error raised for malformed BED data or invalid coordinate operations.
pub suberror AlignBedError {
  AlignBedError(String)
}

///|
/// A BED score. Biopython preserves non-numeric score tokens as text.
pub enum AlignBedScore {
  AlignBedNumeric(Double)
  AlignBedText(String)
} derive(Eq, Debug)
pub extend AlignBedScore with Eq::{equal, not_equal}
pub extend AlignBedScore with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One point in the target/query absolute coordinate path.
pub struct AlignBedCoordinate {
  target : Int
  query : Int
} derive(Eq, Debug)
pub extend AlignBedCoordinate with Eq::{equal, not_equal}
pub extend AlignBedCoordinate with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One aligned BED block in absolute target and query coordinates.
pub struct AlignBedBlock {
  target_start : Int
  target_end : Int
  query_start : Int
  query_end : Int
} derive(Eq, Debug)
pub extend AlignBedBlock with Eq::{equal, not_equal}
pub extend AlignBedBlock with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Coordinate-only alignment counts.
pub struct AlignBedCounts {
  aligned : Int
  target_skip_bases : Int
  query_skip_bases : Int
  target_skip_opens : Int
  query_skip_opens : Int
  blocks : Int
  columns : Int
} derive(Eq, Debug)
pub extend AlignBedCounts with Eq::{equal, not_equal}
pub extend AlignBedCounts with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One BED pairwise alignment.
///
/// Target coordinates may be stored in either orientation by callers; BED
/// serialization normalizes them to increasing genomic coordinates. Query
/// coordinates decrease for reverse-strand alignments.
pub struct AlignBedAlignment {
  target_id : String
  query_id : String?
  coordinates : Array[AlignBedCoordinate]
  score : AlignBedScore?
  thick_start : Int?
  thick_end : Int?
  item_rgb : String?
  source_columns : Int
} derive(Eq, Debug)
pub extend AlignBedAlignment with Eq::{equal, not_equal}
pub extend AlignBedAlignment with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Summary across all records in a BED document.
pub struct AlignBedSummary {
  alignment_count : Int
  target_count : Int
  query_count : Int
  plus_strand_count : Int
  minus_strand_count : Int
  aligned_bases : Int
  target_skip_bases : Int
} derive(Eq, Debug)
pub extend AlignBedSummary with Eq::{equal, not_equal}
pub extend AlignBedSummary with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// A standalone BED document containing zero or more pairwise alignments.
pub struct AlignBedDocument {
  alignments : Array[AlignBedAlignment]
} derive(Eq, Debug)
pub extend AlignBedDocument with Eq::{equal, not_equal}
pub extend AlignBedDocument with @moonbitlang/core/debug.Debug::{to_repr}

///|
fn align_bed_fail(message : String) -> Unit raise AlignBedError {
  raise AlignBedError(message)
}

///|
fn align_bed_abs(value : Int) -> Int {
  if value < 0 {
    -value
  } else {
    value
  }
}

///|
fn align_bed_min(left : Int, right : Int) -> Int {
  if left < right {
    left
  } else {
    right
  }
}

///|
fn align_bed_max(left : Int, right : Int) -> Int {
  if left > right {
    left
  } else {
    right
  }
}

///|
fn align_bed_sign(value : Int) -> Int {
  if value < 0 {
    -1
  } else if value > 0 {
    1
  } else {
    0
  }
}

///|
fn align_bed_copy_coordinates(
  coordinates : Array[AlignBedCoordinate],
) -> Array[AlignBedCoordinate] {
  let copied : Array[AlignBedCoordinate] = []
  for coordinate in coordinates {
    copied.push(coordinate)
  }
  copied
}

///|
fn align_bed_copy_alignments(
  alignments : Array[AlignBedAlignment],
) -> Array[AlignBedAlignment] {
  let copied : Array[AlignBedAlignment] = []
  for alignment in alignments {
    copied.push(alignment)
  }
  copied
}

///|
fn align_bed_validate_token(
  value : String,
  label : String,
  allow_empty : Bool,
) -> Unit raise AlignBedError {
  if !allow_empty && value.length() == 0 {
    align_bed_fail(label + " must not be empty")
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code <= 32 || code == 127 {
      align_bed_fail(
        label + " must not contain whitespace or control characters",
      )
    }
  }
}

///|
fn align_bed_validate_score(score : AlignBedScore) -> Unit raise AlignBedError {
  match score {
    AlignBedNumeric(value) =>
      if value.is_nan() || value.abs() > 1.0e300 {
        align_bed_fail("BED score must be finite")
      }
    AlignBedText(value) => align_bed_validate_token(value, "BED score", false)
  }
}

///|
fn align_bed_reversed_coordinates(
  coordinates : Array[AlignBedCoordinate],
) -> Array[AlignBedCoordinate] {
  let reversed : Array[AlignBedCoordinate] = []
  for index in 0..<coordinates.length() {
    reversed.push(coordinates[coordinates.length() - 1 - index])
  }
  reversed
}

///|
fn align_bed_forward_coordinates(
  coordinates : Array[AlignBedCoordinate],
) -> Array[AlignBedCoordinate] {
  if coordinates.length() > 1 &&
    coordinates[0].target > coordinates[coordinates.length() - 1].target {
    align_bed_reversed_coordinates(coordinates)
  } else {
    align_bed_copy_coordinates(coordinates)
  }
}

///|
fn align_bed_blocks_from_coordinates(
  coordinates : Array[AlignBedCoordinate],
) -> Array[AlignBedBlock] {
  let forward = align_bed_forward_coordinates(coordinates)
  let blocks : Array[AlignBedBlock] = []
  if forward.length() < 2 {
    return blocks
  }
  let mut target_start = forward[0].target
  let mut query_start = forward[0].query
  for index in 1..<forward.length() {
    let target_end = forward[index].target
    let query_end = forward[index].query
    if target_start == target_end {
      query_start = query_end
    } else if query_start == query_end {
      target_start = target_end
    } else {
      blocks.push(AlignBedBlock::{
        target_start,
        target_end,
        query_start: align_bed_min(query_start, query_end),
        query_end: align_bed_max(query_start, query_end),
      })
      target_start = target_end
      query_start = query_end
    }
  }
  blocks
}

///|
fn align_bed_validate_alignment(
  alignment : AlignBedAlignment,
) -> Unit raise AlignBedError {
  align_bed_validate_token(alignment.target_id, "BED target identifier", false)
  match alignment.query_id {
    Some(value) =>
      align_bed_validate_token(value, "BED query identifier", false)
    None => ()
  }
  if alignment.source_columns < 3 || alignment.source_columns > 12 {
    align_bed_fail("BED column count must be between 3 and 12")
  }
  match alignment.score {
    Some(value) => align_bed_validate_score(value)
    None => ()
  }
  match alignment.item_rgb {
    Some(value) => align_bed_validate_token(value, "BED itemRgb", false)
    None => ()
  }
  if alignment.coordinates.length() < 2 {
    align_bed_fail("BED coordinate path must contain at least two points")
  }
  let mut target_direction = 0
  let mut query_direction = 0
  let mut aligned_segments = 0
  for index in 0..<alignment.coordinates.length() {
    let point = alignment.coordinates[index]
    if point.target < 0 || point.query < 0 {
      align_bed_fail("BED coordinates must be non-negative")
    }
    if index + 1 < alignment.coordinates.length() {
      let next = alignment.coordinates[index + 1]
      let target_step = next.target - point.target
      let query_step = next.query - point.query
      if target_step == 0 && query_step == 0 {
        align_bed_fail("BED coordinate path contains a zero-length segment")
      }
      if target_step != 0 {
        let direction = align_bed_sign(target_step)
        if target_direction == 0 {
          target_direction = direction
        } else if target_direction != direction {
          align_bed_fail("BED target coordinates must be monotonic")
        }
      }
      if query_step != 0 {
        let direction = align_bed_sign(query_step)
        if query_direction == 0 {
          query_direction = direction
        } else if query_direction != direction {
          align_bed_fail("BED query coordinates must be monotonic")
        }
      }
      if target_step != 0 && query_step != 0 {
        if align_bed_abs(target_step) != align_bed_abs(query_step) {
          align_bed_fail(
            "BED aligned coordinate steps must have equal target and query sizes",
          )
        }
        aligned_segments = aligned_segments + 1
      }
    }
  }
  if aligned_segments == 0 {
    align_bed_fail("BED alignment must contain at least one aligned block")
  }
  let blocks = align_bed_blocks_from_coordinates(alignment.coordinates)
  if blocks.length() == 0 {
    align_bed_fail("BED alignment has no serializable aligned blocks")
  }
  let interval_start = blocks[0].target_start
  let interval_end = blocks[blocks.length() - 1].target_end
  match alignment.thick_start {
    Some(value) =>
      if value < interval_start || value > interval_end {
        align_bed_fail("BED thickStart must lie inside the target interval")
      }
    None => ()
  }
  match alignment.thick_end {
    Some(value) =>
      if value < interval_start || value > interval_end {
        align_bed_fail("BED thickEnd must lie inside the target interval")
      }
    None => ()
  }
  match (alignment.thick_start, alignment.thick_end) {
    (Some(start), Some(end)) =>
      if start > end {
        align_bed_fail("BED thickStart must not exceed thickEnd")
      }
    _ => ()
  }
}

///|
/// Construct one coordinate point.
pub fn AlignBedCoordinate::create(
  target : Int,
  query : Int,
) -> AlignBedCoordinate raise AlignBedError {
  if target < 0 || query < 0 {
    align_bed_fail("BED coordinates must be non-negative")
  }
  AlignBedCoordinate::{ target, query, }
}

///|
/// Construct a finite numeric BED score.
pub fn AlignBedScore::numeric(
  value : Double,
) -> AlignBedScore raise AlignBedError {
  let score = AlignBedNumeric(value)
  align_bed_validate_score(score)
  score
}

///|
/// Construct a non-empty textual BED score token.
pub fn AlignBedScore::text(value : String) -> AlignBedScore raise AlignBedError {
  let score = AlignBedText(value)
  align_bed_validate_score(score)
  score
}

///|
/// Return the numeric score, if this score was parsed as a number.
pub fn AlignBedScore::number(self : AlignBedScore) -> Double? {
  match self {
    AlignBedNumeric(value) => Some(value)
    AlignBedText(_) => None
  }
}

///|
/// Return the original textual score, if this is a non-numeric token.
pub fn AlignBedScore::text_value(self : AlignBedScore) -> String? {
  match self {
    AlignBedNumeric(_) => None
    AlignBedText(value) => Some(value)
  }
}

///|
/// Construct and validate a pairwise BED alignment.
pub fn AlignBedAlignment::create(
  target_id : String,
  query_id : String?,
  coordinates : Array[AlignBedCoordinate],
  score? : AlignBedScore? = None,
  thick_start? : Int? = None,
  thick_end? : Int? = None,
  item_rgb? : String? = None,
  source_columns? : Int = 12,
) -> AlignBedAlignment raise AlignBedError {
  let alignment = AlignBedAlignment::{
    target_id,
    query_id,
    coordinates: align_bed_copy_coordinates(coordinates),
    score,
    thick_start,
    thick_end,
    item_rgb,
    source_columns,
  }
  align_bed_validate_alignment(alignment)
  alignment
}

///|
/// Construct a BED document and copy the input alignment array.
pub fn AlignBedDocument::create(
  alignments : Array[AlignBedAlignment],
) -> AlignBedDocument raise AlignBedError {
  let copied = align_bed_copy_alignments(alignments)
  for alignment in copied {
    align_bed_validate_alignment(alignment)
  }
  AlignBedDocument::{ alignments: copied, }
}

///|
fn align_bed_strip_cr(value : String) -> String {
  if value.length() > 0 &&
    value.unsafe_get(value.length() - 1).to_int() == '\r'.to_int() {
    value[0:value.length() - 1].to_owned()
  } else {
    value
  }
}

///|
fn align_bed_split_lines(value : String) -> Array[String] {
  let lines : Array[String] = []
  let mut start = 0
  for index in 0..<value.length() {
    if value.unsafe_get(index).to_int() == '\n'.to_int() {
      lines.push(align_bed_strip_cr(value[start:index].to_owned()))
      start = index + 1
    }
  }
  if start < value.length() {
    lines.push(align_bed_strip_cr(value[start:value.length()].to_owned()))
  } else if value.length() > 0 {
    lines.push("")
  }
  lines
}

///|
fn align_bed_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\r'.to_int() ||
  code == '\n'.to_int()
}

///|
fn align_bed_trim(value : String) -> String {
  let mut start = 0
  while start < value.length() &&
        align_bed_is_whitespace(value.unsafe_get(start).to_int()) {
    start = start + 1
  }
  let mut end = value.length()
  while end > start &&
        align_bed_is_whitespace(value.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  value[start:end].to_owned()
}

///|
fn align_bed_split_whitespace(value : String) -> Array[String] {
  let words : Array[String] = []
  let mut index = 0
  while index < value.length() {
    while index < value.length() &&
          align_bed_is_whitespace(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index >= value.length() {
      break
    }
    let start = index
    while index < value.length() &&
          !align_bed_is_whitespace(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    words.push(value[start:index].to_owned())
  }
  words
}

///|
fn align_bed_parse_uint(
  value : String,
  label : String,
) -> Int raise AlignBedError {
  if value.length() == 0 {
    align_bed_fail(label + " must be an unsigned integer")
  }
  let mut result = 0
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      align_bed_fail(label + " must be an unsigned integer")
    }
    let digit = code - '0'.to_int()
    if result > (2147483647 - digit) / 10 {
      align_bed_fail(label + " exceeds the supported integer range")
    }
    result = result * 10 + digit
  }
  result
}

///|
fn align_bed_is_number(value : String) -> Bool {
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
fn align_bed_parse_score(value : String) -> AlignBedScore raise AlignBedError {
  if !align_bed_is_number(value) {
    return AlignBedScore::text(value)
  }
  let parsed_value = if value.unsafe_get(0).to_int() == '+'.to_int() {
    value[1:value.length()].to_owned()
  } else {
    value
  }
  let number = match parse_double(parsed_value) {
    Some(result) => result
    None => {
      align_bed_fail("BED numeric score could not be parsed")
      0.0
    }
  }
  AlignBedScore::numeric(number)
}

///|
fn align_bed_parse_uint_list(
  value : String,
  label : String,
) -> Array[Int] raise AlignBedError {
  let values : Array[Int] = []
  let mut start = 0
  for index in 0..<value.length() {
    if value.unsafe_get(index).to_int() == ','.to_int() {
      if index == start {
        if index + 1 == value.length() {
          start = index + 1
          continue
        }
        align_bed_fail(label + " contains an empty item")
      }
      values.push(
        align_bed_parse_uint(value[start:index].to_owned(), label + " item"),
      )
      start = index + 1
    }
  }
  if start < value.length() {
    values.push(
      align_bed_parse_uint(
        value[start:value.length()].to_owned(),
        label + " item",
      ),
    )
  }
  if values.length() == 0 {
    align_bed_fail(label + " must contain at least one integer")
  }
  values
}

///|
fn align_bed_validate_blocks(
  span : Int,
  count : Int,
  sizes : Array[Int],
  starts : Array[Int],
) -> Unit raise AlignBedError {
  if count <= 0 {
    align_bed_fail("BED blockCount must be positive")
  }
  if sizes.length() != count {
    align_bed_fail("BED blockCount does not match blockSizes")
  }
  if starts.length() != count {
    align_bed_fail("BED blockCount does not match blockStarts")
  }
  if starts[0] != 0 {
    align_bed_fail("the first BED block must start at chromStart")
  }
  let mut previous_end = 0
  for index in 0..<count {
    let size = sizes[index]
    let start = starts[index]
    if size <= 0 {
      align_bed_fail("BED block sizes must be positive")
    }
    if start < previous_end {
      align_bed_fail("BED blocks must be sorted and non-overlapping")
    }
    if start > span || size > span - start {
      align_bed_fail("BED block extends beyond chromEnd")
    }
    previous_end = start + size
  }
  if previous_end != span {
    align_bed_fail("the last BED block must end at chromEnd")
  }
}

///|
fn align_bed_infer_starts(
  span : Int,
  sizes : Array[Int],
) -> Array[Int] raise AlignBedError {
  let starts : Array[Int] = []
  let mut position = 0
  for size in sizes {
    if size <= 0 || position > span || size > span - position {
      align_bed_fail("BED11 blockSizes cannot be placed inside chromEnd")
    }
    starts.push(position)
    position = position + size
  }
  if position != span {
    align_bed_fail(
      "BED11 blockStarts are ambiguous unless blockSizes cover the full span",
    )
  }
  starts
}

///|
fn align_bed_coordinates_from_blocks(
  chrom_start : Int,
  strand : String,
  sizes : Array[Int],
  starts : Array[Int],
) -> Array[AlignBedCoordinate] {
  let coordinates : Array[AlignBedCoordinate] = []
  let mut query_position = 0
  let mut target_position = chrom_start
  coordinates.push(AlignBedCoordinate::{ target: chrom_start, query: 0, })
  for index in 0..<sizes.length() {
    let block_start = chrom_start + starts[index]
    if block_start != target_position {
      coordinates.push(AlignBedCoordinate::{
        target: block_start,
        query: query_position,
      })
    }
    target_position = block_start + sizes[index]
    query_position = query_position + sizes[index]
    coordinates.push(AlignBedCoordinate::{
      target: target_position,
      query: query_position,
    })
  }
  if strand == "-" {
    for index in 0..<coordinates.length() {
      coordinates[index] = AlignBedCoordinate::{
        target: coordinates[index].target,
        query: query_position - coordinates[index].query,
      }
    }
  }
  coordinates
}

///|
fn align_bed_parse_line(
  line : String,
  line_number : Int,
) -> AlignBedAlignment raise AlignBedError {
  let words = align_bed_split_whitespace(line)
  let bed_columns = words.length()
  if bed_columns < 3 || bed_columns > 12 {
    align_bed_fail(
      "expected between 3 and 12 BED columns at line " +
      line_number.to_string() +
      ", found " +
      bed_columns.to_string(),
    )
  }
  let target_id = words[0]
  align_bed_validate_token(target_id, "BED target identifier", false)
  let chrom_start = align_bed_parse_uint(words[1], "BED chromStart")
  let chrom_end = align_bed_parse_uint(words[2], "BED chromEnd")
  if chrom_end <= chrom_start {
    align_bed_fail("BED coordinates must define a positive half-open interval")
  }
  let span = chrom_end - chrom_start
  let query_id = if bed_columns >= 4 {
    align_bed_validate_token(words[3], "BED query identifier", false)
    Some(words[3])
  } else {
    None
  }
  let score = if bed_columns >= 5 {
    Some(align_bed_parse_score(words[4]))
  } else {
    None
  }
  let strand = if bed_columns >= 6 {
    if words[5] != "+" && words[5] != "-" && words[5] != "." {
      align_bed_fail("BED strand must be '+', '-', or '.'")
    }
    if words[5] == "-" {
      "-"
    } else {
      "+"
    }
  } else {
    "+"
  }
  let thick_start = if bed_columns >= 7 {
    Some(align_bed_parse_uint(words[6], "BED thickStart"))
  } else {
    None
  }
  let thick_end = if bed_columns >= 8 {
    Some(align_bed_parse_uint(words[7], "BED thickEnd"))
  } else {
    None
  }
  let item_rgb = if bed_columns >= 9 {
    align_bed_validate_token(words[8], "BED itemRgb", false)
    Some(words[8])
  } else {
    None
  }
  let sizes : Array[Int] = []
  let starts : Array[Int] = []
  if bed_columns <= 9 {
    sizes.push(span)
    starts.push(0)
  } else {
    let block_count = align_bed_parse_uint(words[9], "BED blockCount")
    if bed_columns == 10 {
      if block_count != 1 {
        align_bed_fail(
          "BED10 records with multiple blocks lack blockSizes and blockStarts",
        )
      }
      sizes.push(span)
      starts.push(0)
    } else {
      let parsed_sizes = align_bed_parse_uint_list(words[10], "BED blockSizes")
      if parsed_sizes.length() != block_count {
        align_bed_fail("BED blockCount does not match blockSizes")
      }
      for size in parsed_sizes {
        sizes.push(size)
      }
      if bed_columns == 11 {
        let inferred = align_bed_infer_starts(span, sizes)
        for start in inferred {
          starts.push(start)
        }
      } else {
        let parsed_starts = align_bed_parse_uint_list(
          words[11],
          "BED blockStarts",
        )
        for start in parsed_starts {
          starts.push(start)
        }
      }
      align_bed_validate_blocks(span, block_count, sizes, starts)
    }
  }
  let coordinates = align_bed_coordinates_from_blocks(
    chrom_start, strand, sizes, starts,
  )
  AlignBedAlignment::create(
    target_id,
    query_id,
    coordinates,
    score~,
    thick_start~,
    thick_end~,
    item_rgb~,
    source_columns=bed_columns,
  )
}

///|
/// Parse all BED records from a string.
pub fn align_bed_parse(text : String) -> AlignBedDocument raise AlignBedError {
  if text.length() == 0 {
    return AlignBedDocument::create([])
  }
  let lines = align_bed_split_lines(text)
  let alignments : Array[AlignBedAlignment] = []
  for index in 0..<lines.length() {
    let line = align_bed_trim(lines[index])
    if line.length() == 0 {
      if index + 1 == lines.length() {
        continue
      }
      align_bed_fail("blank BED record at line " + (index + 1).to_string())
    }
    alignments.push(align_bed_parse_line(line, index + 1))
  }
  AlignBedDocument::create(alignments)
}

///|
/// Return aligned blocks after normalizing the target axis to increasing order.
pub fn AlignBedAlignment::blocks(
  self : AlignBedAlignment,
) -> Array[AlignBedBlock] {
  align_bed_blocks_from_coordinates(self.coordinates)
}

///|
/// Return the relative query strand encoded by this coordinate path.
pub fn AlignBedAlignment::strand(self : AlignBedAlignment) -> String {
  let forward = align_bed_forward_coordinates(self.coordinates)
  if forward[0].query > forward[forward.length() - 1].query {
    "-"
  } else {
    "+"
  }
}

///|
/// Return the first aligned target coordinate.
pub fn AlignBedAlignment::target_start(self : AlignBedAlignment) -> Int {
  self.blocks()[0].target_start
}

///|
/// Return the exclusive final aligned target coordinate.
pub fn AlignBedAlignment::target_end(self : AlignBedAlignment) -> Int {
  let blocks = self.blocks()
  blocks[blocks.length() - 1].target_end
}

///|
/// Return the query length represented by BED blocks.
pub fn AlignBedAlignment::query_size(self : AlignBedAlignment) -> Int {
  let mut size = 0
  for block in self.blocks() {
    size = size + (block.target_end - block.target_start)
  }
  size
}

///|
/// Return the number of aligned target/query residues.
pub fn AlignBedAlignment::aligned_bases(self : AlignBedAlignment) -> Int {
  self.query_size()
}

///|
/// Return the target span between the first and final aligned blocks.
pub fn AlignBedAlignment::target_span(self : AlignBedAlignment) -> Int {
  self.target_end() - self.target_start()
}

///|
/// Return coordinate-only alignment statistics.
pub fn AlignBedAlignment::counts(self : AlignBedAlignment) -> AlignBedCounts {
  let mut aligned = 0
  let mut target_skips = 0
  let mut query_skips = 0
  let mut target_opens = 0
  let mut query_opens = 0
  let mut columns = 0
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    let target_step = align_bed_abs(next.target - current.target)
    let query_step = align_bed_abs(next.query - current.query)
    if target_step > 0 && query_step > 0 {
      aligned = aligned + target_step
    } else if target_step > 0 {
      target_skips = target_skips + target_step
      target_opens = target_opens + 1
    } else {
      query_skips = query_skips + query_step
      query_opens = query_opens + 1
    }
    columns = columns + align_bed_max(target_step, query_step)
  }
  AlignBedCounts::{
    aligned,
    target_skip_bases: target_skips,
    query_skip_bases: query_skips,
    target_skip_opens: target_opens,
    query_skip_opens: query_opens,
    blocks: self.blocks().length(),
    columns,
  }
}

///|
fn align_bed_map_on_segment(
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
    align_bed_abs(source_step) != align_bed_abs(destination_step) {
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
/// Map one target residue to the query, returning `None` in an intron/gap.
pub fn AlignBedAlignment::map_target_position(
  self : AlignBedAlignment,
  position : Int,
) -> Int? {
  if position < 0 {
    return None
  }
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    match
      align_bed_map_on_segment(
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
pub fn AlignBedAlignment::map_query_position(
  self : AlignBedAlignment,
  position : Int,
) -> Int? {
  if position < 0 {
    return None
  }
  for index in 0..<(self.coordinates.length() - 1) {
    let current = self.coordinates[index]
    let next = self.coordinates[index + 1]
    match
      align_bed_map_on_segment(
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
/// Return true if the aligned target interval overlaps a half-open query.
pub fn AlignBedAlignment::overlaps(
  self : AlignBedAlignment,
  target_id : String,
  start : Int,
  end : Int,
) -> Bool {
  end > start &&
  self.target_id == target_id &&
  self.target_start() < end &&
  start < self.target_end()
}

///|
fn align_bed_join_tabs(fields : Array[String]) -> String {
  let output = StringBuilder()
  for index in 0..<fields.length() {
    if index > 0 {
      output.write_char('\t')
    }
    output.write_string(fields[index])
  }
  output.to_string()
}

///|
fn align_bed_integer_list(values : Array[Int]) -> String {
  let output = StringBuilder()
  for value in values {
    output.write_string(value.to_string())
    output.write_char(',')
  }
  output.to_string()
}

///|
fn align_bed_score_text(score : AlignBedScore) -> String {
  match score {
    AlignBedNumeric(value) => value.to_string()
    AlignBedText(value) => value
  }
}

///|
/// Format one alignment using BED3 through BED12.
pub fn align_bed_format(
  alignment : AlignBedAlignment,
  bed_columns? : Int = 12,
) -> String raise AlignBedError {
  if bed_columns < 3 || bed_columns > 12 {
    align_bed_fail("BED column count must be between 3 and 12")
  }
  align_bed_validate_alignment(alignment)
  let blocks = alignment.blocks()
  if blocks.length() == 0 {
    return ""
  }
  let chrom_start = blocks[0].target_start
  let chrom_end = blocks[blocks.length() - 1].target_end
  let fields : Array[String] = [
    alignment.target_id,
    chrom_start.to_string(),
    chrom_end.to_string(),
  ]
  if bed_columns == 3 {
    return align_bed_join_tabs(fields) + "\n"
  }
  fields.push(alignment.query_id.unwrap_or("query"))
  if bed_columns == 4 {
    return align_bed_join_tabs(fields) + "\n"
  }
  fields.push(
    match alignment.score {
      Some(score) => align_bed_score_text(score)
      None => "0"
    },
  )
  if bed_columns == 5 {
    return align_bed_join_tabs(fields) + "\n"
  }
  fields.push(alignment.strand())
  if bed_columns == 6 {
    return align_bed_join_tabs(fields) + "\n"
  }
  fields.push(alignment.thick_start.unwrap_or(chrom_start).to_string())
  if bed_columns == 7 {
    return align_bed_join_tabs(fields) + "\n"
  }
  fields.push(alignment.thick_end.unwrap_or(chrom_end).to_string())
  if bed_columns == 8 {
    return align_bed_join_tabs(fields) + "\n"
  }
  fields.push(alignment.item_rgb.unwrap_or("0"))
  if bed_columns == 9 {
    return align_bed_join_tabs(fields) + "\n"
  }
  fields.push(blocks.length().to_string())
  if bed_columns == 10 {
    return align_bed_join_tabs(fields) + "\n"
  }
  let block_sizes : Array[Int] = []
  let block_starts : Array[Int] = []
  for block in blocks {
    block_sizes.push(block.target_end - block.target_start)
    block_starts.push(block.target_start - chrom_start)
  }
  fields.push(align_bed_integer_list(block_sizes))
  if bed_columns == 11 {
    return align_bed_join_tabs(fields) + "\n"
  }
  fields.push(align_bed_integer_list(block_starts))
  align_bed_join_tabs(fields) + "\n"
}

///|
/// Serialize all records in a BED document.
pub fn align_bed_write(
  document : AlignBedDocument,
  bed_columns? : Int = 12,
) -> String raise AlignBedError {
  let validated = AlignBedDocument::create(document.alignments)
  let output = StringBuilder()
  for alignment in validated.alignments {
    output.write_string(align_bed_format(alignment, bed_columns~))
  }
  output.to_string()
}

///|
/// Find all records overlapping one target half-open interval.
pub fn AlignBedDocument::search(
  self : AlignBedDocument,
  target_id : String,
  start : Int,
  end : Int,
) -> Array[AlignBedAlignment] raise AlignBedError {
  if start < 0 || end <= start {
    align_bed_fail("BED search interval must satisfy 0 <= start < end")
  }
  let results : Array[AlignBedAlignment] = []
  for alignment in self.alignments {
    if alignment.overlaps(target_id, start, end) {
      results.push(alignment)
    }
  }
  results
}

///|
/// Return distinct target identifiers in first-seen order.
pub fn AlignBedDocument::targets(self : AlignBedDocument) -> Array[String] {
  let targets : Array[String] = []
  let seen : Map[String, Bool] = Map([])
  for alignment in self.alignments {
    if !seen.contains(alignment.target_id) {
      seen[alignment.target_id] = true
      targets.push(alignment.target_id)
    }
  }
  targets
}

///|
/// Return distinct non-missing query identifiers in first-seen order.
pub fn AlignBedDocument::queries(self : AlignBedDocument) -> Array[String] {
  let queries : Array[String] = []
  let seen : Map[String, Bool] = Map([])
  for alignment in self.alignments {
    match alignment.query_id {
      Some(query) =>
        if !seen.contains(query) {
          seen[query] = true
          queries.push(query)
        }
      None => ()
    }
  }
  queries
}

///|
/// Summarize all records in the document.
pub fn AlignBedDocument::summary(self : AlignBedDocument) -> AlignBedSummary {
  let mut plus = 0
  let mut minus = 0
  let mut aligned = 0
  let mut target_skips = 0
  for alignment in self.alignments {
    if alignment.strand() == "-" {
      minus = minus + 1
    } else {
      plus = plus + 1
    }
    let counts = alignment.counts()
    aligned = aligned + counts.aligned
    target_skips = target_skips + counts.target_skip_bases
  }
  AlignBedSummary::{
    alignment_count: self.alignments.length(),
    target_count: self.targets().length(),
    query_count: self.queries().length(),
    plus_strand_count: plus,
    minus_strand_count: minus,
    aligned_bases: aligned,
    target_skip_bases: target_skips,
  }
}

///|
/// A compact Biopython-style BED3/BED12 fixture.
pub fn align_bed_example_text() -> String {
  "chr22\t1000\t5000\tmRNA1\t960\t+\t1200\t4900\t255,0,0\t2\t567,488,\t0,3512,\n" +
  "chr22\t2000\t6000\tmRNA2\t900\t-\t2300\t5960\t0,255,0\t2\t433,399,\t0,3601,\n" +
  "chr7\t100\t180"
}
