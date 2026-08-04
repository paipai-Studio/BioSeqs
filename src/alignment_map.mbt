// Compose coordinate-based alignments without realigning their sequences.
//
// This implements the core semantics of Biopython Bio.Align.Alignment.map
// and mapall. Coordinates are zero-based, half-open sequence boundaries.

///|
/// Error raised by coordinate alignment construction or projection.
pub suberror AlignmentMapError {
  AlignmentMapError(String)
}

///|
/// One aligned block in a coordinate-based pairwise alignment.
pub struct CoordinateAlignmentBlock {
  target_start : Int
  target_end : Int
  query_start : Int
  query_end : Int
  size : Int
} derive(Eq, Debug)

///|
/// Coordinate representation of a pairwise alignment.
///
/// Sequence text may be empty when only lengths and coordinates are known.
/// Target coordinates are non-decreasing. Query coordinates may increase or
/// decrease, allowing reverse-strand alignments.
pub struct CoordinatePairwiseAlignment {
  target_name : String
  query_name : String
  target_sequence : String
  query_sequence : String
  target_length : Int
  query_length : Int
  target_coordinates : Array[Int]
  query_coordinates : Array[Int]
} derive(Eq, Debug)

///|
/// Coordinate representation of a multiple sequence alignment.
pub struct CoordinateMultipleAlignment {
  names : Array[String]
  sequences : Array[String]
  coordinates : Array[Array[Int]]
} derive(Eq, Debug)

///|
/// Basic alignment statistics derived from a coordinate path.
pub struct CoordinateAlignmentCounts {
  aligned : Int
  target_gap_bases : Int
  query_gap_bases : Int
  target_gap_events : Int
  query_gap_events : Int
  blocks : Int
} derive(Eq, Debug)

///|
fn alignment_map_copy_ints(values : Array[Int]) -> Array[Int] {
  let copy : Array[Int] = []
  for value in values {
    copy.push(value)
  }
  copy
}

///|
fn alignment_map_copy_strings(values : Array[String]) -> Array[String] {
  let copy : Array[String] = []
  for value in values {
    copy.push(value)
  }
  copy
}

///|
fn alignment_map_abs(value : Int) -> Int {
  if value < 0 {
    -value
  } else {
    value
  }
}

///|
fn alignment_map_min(left : Int, right : Int) -> Int {
  if left < right {
    left
  } else {
    right
  }
}

///|
fn alignment_map_max(left : Int, right : Int) -> Int {
  if left > right {
    left
  } else {
    right
  }
}

///|
fn alignment_map_validate_name(
  name : String,
  label : String,
) -> Unit raise AlignmentMapError {
  if name.length() == 0 {
    raise AlignmentMapError(label + " must not be empty")
  }
}

///|
fn alignment_map_validate_sequence(
  sequence : String,
  expected_length : Int,
  label : String,
) -> Unit raise AlignmentMapError {
  if expected_length < 0 {
    raise AlignmentMapError(label + " length must be non-negative")
  }
  if sequence.length() > 0 && sequence.length() != expected_length {
    raise AlignmentMapError(label + " length does not match its declaration")
  }
  let mut index = 0
  while index < sequence.length() {
    let residue = sequence.unsafe_get(index).unsafe_to_char()
    if residue == '-' {
      raise AlignmentMapError(label + " must not contain gap characters")
    }
    if residue == ' ' || residue == '\t' || residue == '\n' || residue == '\r' {
      raise AlignmentMapError(label + " must not contain whitespace")
    }
    index = index + 1
  }
}

///|
fn alignment_map_validate_coordinates(
  target_length : Int,
  query_length : Int,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> Unit raise AlignmentMapError {
  if target_coordinates.length() != query_coordinates.length() {
    raise AlignmentMapError(
      "target and query coordinates must have equal length",
    )
  }
  if target_coordinates.length() == 1 {
    raise AlignmentMapError(
      "an alignment path must be empty or contain at least two points",
    )
  }
  let mut index = 0
  while index < target_coordinates.length() {
    let target = target_coordinates[index]
    let query = query_coordinates[index]
    if target < 0 || target > target_length {
      raise AlignmentMapError("target coordinate is out of bounds")
    }
    if query < 0 || query > query_length {
      raise AlignmentMapError("query coordinate is out of bounds")
    }
    index = index + 1
  }
  let mut query_direction = 0
  index = 0
  while index + 1 < target_coordinates.length() {
    let target_step = target_coordinates[index + 1] - target_coordinates[index]
    let query_step = query_coordinates[index + 1] - query_coordinates[index]
    if target_step < 0 {
      raise AlignmentMapError("target coordinates must be non-decreasing")
    }
    if target_step == 0 && query_step == 0 {
      raise AlignmentMapError("consecutive alignment coordinates must differ")
    }
    if query_step > 0 {
      if query_direction < 0 {
        raise AlignmentMapError("query coordinates change direction")
      }
      query_direction = 1
    } else if query_step < 0 {
      if query_direction > 0 {
        raise AlignmentMapError("query coordinates change direction")
      }
      query_direction = -1
    }
    index = index + 1
  }
}

///|
fn alignment_map_create(
  target_name : String,
  target_sequence : String,
  target_length : Int,
  query_name : String,
  query_sequence : String,
  query_length : Int,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> CoordinatePairwiseAlignment raise AlignmentMapError {
  alignment_map_validate_name(target_name, "target name")
  alignment_map_validate_name(query_name, "query name")
  alignment_map_validate_sequence(
    target_sequence, target_length, "target sequence",
  )
  alignment_map_validate_sequence(
    query_sequence, query_length, "query sequence",
  )
  alignment_map_validate_coordinates(
    target_length, query_length, target_coordinates, query_coordinates,
  )
  CoordinatePairwiseAlignment::{
    target_name,
    query_name,
    target_sequence,
    query_sequence,
    target_length,
    query_length,
    target_coordinates: alignment_map_copy_ints(target_coordinates),
    query_coordinates: alignment_map_copy_ints(query_coordinates),
  }
}

///|
/// Construct a coordinate pairwise alignment from concrete sequences.
pub fn coordinate_pairwise_alignment(
  target_name : String,
  target_sequence : String,
  query_name : String,
  query_sequence : String,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> CoordinatePairwiseAlignment raise AlignmentMapError {
  alignment_map_create(
    target_name,
    target_sequence,
    target_sequence.length(),
    query_name,
    query_sequence,
    query_sequence.length(),
    target_coordinates,
    query_coordinates,
  )
}

///|
/// Construct an alignment when sequence text is unavailable.
pub fn coordinate_pairwise_alignment_with_lengths(
  target_name : String,
  target_length : Int,
  query_name : String,
  query_length : Int,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
  target_sequence? : String = "",
  query_sequence? : String = "",
) -> CoordinatePairwiseAlignment raise AlignmentMapError {
  alignment_map_create(
    target_name, target_sequence, target_length, query_name, query_sequence, query_length,
    target_coordinates, query_coordinates,
  )
}

///|
/// Adapt the repository's existing PairwiseAlignment representation.
pub fn coordinate_alignment_from_pairwise(
  alignment : PairwiseAlignment,
  target_name? : String = "target",
  query_name? : String = "query",
) -> CoordinatePairwiseAlignment raise AlignmentMapError {
  if alignment.aligned_target.length() != alignment.aligned_query.length() {
    raise AlignmentMapError("aligned pairwise rows must have equal length")
  }
  let target_coordinates = [alignment.target_start]
  let query_coordinates = [alignment.query_start]
  let mut target = alignment.target_start
  let mut query = alignment.query_start
  let mut column = 0
  while column < alignment.aligned_target.length() {
    if alignment.aligned_target.unsafe_get(column).unsafe_to_char() != '-' {
      target = target + 1
    }
    if alignment.aligned_query.unsafe_get(column).unsafe_to_char() != '-' {
      query = query + 1
    }
    target_coordinates.push(target)
    query_coordinates.push(query)
    column = column + 1
  }
  if target != alignment.target_end || query != alignment.query_end {
    raise AlignmentMapError("pairwise coordinates disagree with aligned rows")
  }
  coordinate_pairwise_alignment(
    target_name,
    alignment.target,
    query_name,
    alignment.query,
    target_coordinates,
    query_coordinates,
  )
}

///|
fn alignment_map_reverse(values : Array[Int]) -> Array[Int] {
  let result : Array[Int] = []
  let mut index = values.length()
  while index > 0 {
    index = index - 1
    result.push(values[index])
  }
  result
}

///|
fn alignment_map_relationship(
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
  label : String,
) -> Int raise AlignmentMapError {
  let mut relationship = 0
  let mut index = 0
  while index + 1 < target_coordinates.length() {
    let target_step = target_coordinates[index + 1] - target_coordinates[index]
    let query_step = query_coordinates[index + 1] - query_coordinates[index]
    if target_step != 0 && query_step != 0 {
      let current = if target_step * query_step > 0 { 1 } else { -1 }
      if relationship != 0 && relationship != current {
        raise AlignmentMapError("inconsistent steps in " + label)
      }
      relationship = current
    }
    index = index + 1
  }
  if relationship == 0 {
    1
  } else {
    relationship
  }
}

///|
fn alignment_map_validate_equal_steps(
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
  label : String,
) -> Unit raise AlignmentMapError {
  let mut index = 0
  while index + 1 < target_coordinates.length() {
    let target_step = target_coordinates[index + 1] - target_coordinates[index]
    let query_step = query_coordinates[index + 1] - query_coordinates[index]
    if target_step < 0 || query_step < 0 {
      raise AlignmentMapError("normalized coordinates decrease in " + label)
    }
    if target_step > 0 && query_step > 0 && target_step != query_step {
      raise AlignmentMapError("unequal aligned step sizes in " + label)
    }
    index = index + 1
  }
}

///|
fn alignment_map_transform_against_reverse_middle(
  alignment : CoordinatePairwiseAlignment,
) -> (Array[Int], Array[Int]) {
  let target_reversed = alignment_map_reverse(alignment.target_coordinates)
  let query_reversed = alignment_map_reverse(alignment.query_coordinates)
  let targets : Array[Int] = []
  let queries : Array[Int] = []
  for value in target_reversed {
    targets.push(alignment.target_length - value)
  }
  for value in query_reversed {
    queries.push(alignment.query_length - value)
  }
  (targets, queries)
}

///|
/// Compose `self` (outer target -> shared middle) with `next`
/// (shared middle -> final query).
pub fn CoordinatePairwiseAlignment::map(
  self : CoordinatePairwiseAlignment,
  next : CoordinatePairwiseAlignment,
) -> CoordinatePairwiseAlignment raise AlignmentMapError {
  if self.query_length != next.target_length {
    raise AlignmentMapError(
      "first query length must equal second target length",
    )
  }
  if self.target_coordinates.length() == 0 ||
    next.target_coordinates.length() == 0 {
    return alignment_map_create(
      self.target_name,
      self.target_sequence,
      self.target_length,
      next.query_name,
      next.query_sequence,
      next.query_length,
      [],
      [],
    )
  }
  let relationship1 = alignment_map_relationship(
    self.target_coordinates,
    self.query_coordinates,
    "first alignment",
  )
  let relationship2 = alignment_map_relationship(
    next.target_coordinates,
    next.query_coordinates,
    "second alignment",
  )
  let coordinates1_target = alignment_map_copy_ints(self.target_coordinates)
  let coordinates1_query = alignment_map_copy_ints(self.query_coordinates)
  let mut coordinates2_target = alignment_map_copy_ints(next.target_coordinates)
  let mut coordinates2_query = alignment_map_copy_ints(next.query_coordinates)
  if relationship1 > 0 {
    if relationship2 < 0 {
      let transformed : Array[Int] = []
      for value in coordinates2_query {
        transformed.push(next.query_length - value)
      }
      coordinates2_query = transformed
    }
  } else {
    let transformed_middle : Array[Int] = []
    for value in coordinates1_query {
      transformed_middle.push(self.query_length - value)
    }
    coordinates1_query.clear()
    for value in transformed_middle {
      coordinates1_query.push(value)
    }
    let (reversed_target, reversed_query) = if relationship2 > 0 {
      alignment_map_transform_against_reverse_middle(next)
    } else {
        let reversed_middle = alignment_map_reverse(next.target_coordinates)
        let normalized_middle : Array[Int] = []
        for value in reversed_middle {
          normalized_middle.push(next.target_length - value)
        }
      (
          normalized_middle,
        alignment_map_reverse(next.query_coordinates),
      )
    }
    coordinates2_target = reversed_target
    coordinates2_query = reversed_query
  }
  alignment_map_validate_equal_steps(
    coordinates1_target, coordinates1_query, "first alignment",
  )
  alignment_map_validate_equal_steps(
    coordinates2_target, coordinates2_query, "second alignment",
  )
  let mut first_segment = 0
  while first_segment + 1 < coordinates1_target.length() {
    if coordinates1_target[first_segment] <
      coordinates1_target[first_segment + 1] &&
      coordinates1_query[first_segment] < coordinates1_query[first_segment + 1] {
      break
    }
    first_segment = first_segment + 1
  }
  if first_segment + 1 >= coordinates1_target.length() {
    return alignment_map_create(
      self.target_name,
      self.target_sequence,
      self.target_length,
      next.query_name,
      next.query_sequence,
      next.query_length,
      [],
      [],
    )
  }
  let mut outer_start = coordinates1_target[first_segment]
  let mut middle_start = coordinates1_query[first_segment]
  let mut outer_end = coordinates1_target[first_segment + 1]
  let mut middle_end = coordinates1_query[first_segment + 1]
  let path_target : Array[Int] = []
  let path_query : Array[Int] = []
  let mut output_target_end = 2147483647
  let mut output_query_end = 2147483647
  let mut second_target_start = 2147483647
  let mut second_query_start = 2147483647
  let mut point2 = 0
  while point2 < coordinates2_target.length() {
    let second_target_end = coordinates2_target[point2]
    let second_query_end = coordinates2_query[point2]
    while second_query_start < second_query_end &&
          second_target_start < second_target_end {
      let mut size = 0
      let mut handled = false
      while !handled {
        if second_target_start < middle_start {
          size = alignment_map_min(second_target_end, middle_start) -
            second_target_start
          handled = true
        } else if second_target_start < middle_end {
          let offset = second_target_start - middle_start
          size = alignment_map_min(second_target_end, middle_end) -
            second_target_start
          let query_start = second_query_start
          let target_start = outer_start + offset
          if target_start != output_target_end ||
            query_start != output_query_end {
            if target_start > output_target_end &&
              query_start > output_query_end {
              path_target.push(target_start)
              path_query.push(output_query_end)
            }
            path_target.push(target_start)
            path_query.push(query_start)
          }
          output_query_end = query_start + size
          output_target_end = target_start + size
          path_target.push(output_target_end)
          path_query.push(output_query_end)
          handled = true
        } else {
          first_segment = first_segment + 1
          let mut found = false
          while first_segment + 1 < coordinates1_target.length() {
            outer_start = coordinates1_target[first_segment]
            middle_start = coordinates1_query[first_segment]
            outer_end = coordinates1_target[first_segment + 1]
            middle_end = coordinates1_query[first_segment + 1]
            if outer_start < outer_end && middle_start < middle_end {
              found = true
              break
            }
            first_segment = first_segment + 1
          }
          if !found {
            size = second_query_end - second_query_start
            handled = true
          }
        }
      }
      if size <= 0 {
        break
      }
      second_query_start = second_query_start + size
      second_target_start = second_target_start + size
    }
    second_target_start = second_target_end
    second_query_start = second_query_end
    point2 = point2 + 1
  }
  if relationship1 != relationship2 {
    let transformed : Array[Int] = []
    for value in path_query {
      transformed.push(next.query_length - value)
    }
    path_query.clear()
    for value in transformed {
      path_query.push(value)
    }
  }
  alignment_map_create(
    self.target_name,
    self.target_sequence,
    self.target_length,
    next.query_name,
    next.query_sequence,
    next.query_length,
    path_target,
    path_query,
  )
}

///|
/// Map several alignments through the same outer alignment.
pub fn CoordinatePairwiseAlignment::map_many(
  self : CoordinatePairwiseAlignment,
  alignments : Array[CoordinatePairwiseAlignment],
) -> Array[CoordinatePairwiseAlignment] raise AlignmentMapError {
  let results : Array[CoordinatePairwiseAlignment] = []
  for alignment in alignments {
    results.push(self.map(alignment))
  }
  results
}

///|
pub fn CoordinatePairwiseAlignment::is_empty(
  self : CoordinatePairwiseAlignment,
) -> Bool {
  self.target_coordinates.length() == 0
}

///|
pub fn CoordinatePairwiseAlignment::is_reverse(
  self : CoordinatePairwiseAlignment,
) -> Bool {
  if self.query_coordinates.length() < 2 {
    return false
  }
  self.query_coordinates[self.query_coordinates.length() - 1] <
  self.query_coordinates[0]
}

///|
pub fn CoordinatePairwiseAlignment::blocks(
  self : CoordinatePairwiseAlignment,
) -> Array[CoordinateAlignmentBlock] {
  let blocks : Array[CoordinateAlignmentBlock] = []
  let mut index = 0
  while index + 1 < self.target_coordinates.length() {
    let target_start = self.target_coordinates[index]
    let target_end = self.target_coordinates[index + 1]
    let query_start = self.query_coordinates[index]
    let query_end = self.query_coordinates[index + 1]
    let target_size = target_end - target_start
    let query_size = alignment_map_abs(query_end - query_start)
    if target_size > 0 && query_size > 0 {
      blocks.push(CoordinateAlignmentBlock::{
        target_start,
        target_end,
        query_start,
        query_end,
        size: alignment_map_min(target_size, query_size),
      })
    }
    index = index + 1
  }
  blocks
}

///|
pub fn CoordinatePairwiseAlignment::counts(
  self : CoordinatePairwiseAlignment,
) -> CoordinateAlignmentCounts {
  let mut aligned = 0
  let mut target_gap_bases = 0
  let mut query_gap_bases = 0
  let mut target_gap_events = 0
  let mut query_gap_events = 0
  let mut blocks = 0
  let mut index = 0
  while index + 1 < self.target_coordinates.length() {
    let target_step = self.target_coordinates[index + 1] -
      self.target_coordinates[index]
    let query_step = alignment_map_abs(
      self.query_coordinates[index + 1] - self.query_coordinates[index],
    )
    if target_step > 0 && query_step > 0 {
      aligned = aligned + alignment_map_min(target_step, query_step)
      blocks = blocks + 1
    } else if target_step == 0 && query_step > 0 {
      target_gap_bases = target_gap_bases + query_step
      target_gap_events = target_gap_events + 1
    } else if target_step > 0 && query_step == 0 {
      query_gap_bases = query_gap_bases + target_step
      query_gap_events = query_gap_events + 1
    }
    index = index + 1
  }
  CoordinateAlignmentCounts::{
    aligned,
    target_gap_bases,
    query_gap_bases,
    target_gap_events,
    query_gap_events,
    blocks,
  }
}

///|
pub fn CoordinatePairwiseAlignment::target_to_query(
  self : CoordinatePairwiseAlignment,
  position : Int,
) -> Int? {
  let mut index = 0
  while index + 1 < self.target_coordinates.length() {
    let target_start = self.target_coordinates[index]
    let target_end = self.target_coordinates[index + 1]
    let query_start = self.query_coordinates[index]
    let query_end = self.query_coordinates[index + 1]
    if target_start <= position &&
      position < target_end &&
      query_start != query_end {
      let offset = position - target_start
      if query_end > query_start {
        return Some(query_start + offset)
      }
      return Some(query_start - offset - 1)
    }
    index = index + 1
  }
  None
}

///|
pub fn CoordinatePairwiseAlignment::query_to_target(
  self : CoordinatePairwiseAlignment,
  position : Int,
) -> Int? {
  let mut index = 0
  while index + 1 < self.query_coordinates.length() {
    let target_start = self.target_coordinates[index]
    let target_end = self.target_coordinates[index + 1]
    let query_start = self.query_coordinates[index]
    let query_end = self.query_coordinates[index + 1]
    if target_start != target_end {
      if query_end > query_start &&
        query_start <= position &&
        position < query_end {
        return Some(target_start + position - query_start)
      }
      if query_end < query_start &&
        query_end <= position &&
        position < query_start {
        return Some(target_start + query_start - position - 1)
      }
    }
    index = index + 1
  }
  None
}

///|
fn alignment_map_repeat(
  builder : StringBuilder,
  value : Char,
  count : Int,
) -> Unit {
  let mut index = 0
  while index < count {
    builder.write_char(value)
    index = index + 1
  }
}

///|
fn alignment_map_complement(value : Char) -> Char {
  match value {
    'A' => 'T'
    'a' => 't'
    'C' => 'G'
    'c' => 'g'
    'G' => 'C'
    'g' => 'c'
    'T' => 'A'
    't' => 'a'
    'U' => 'A'
    'u' => 'a'
    _ => value
  }
}

///|
fn alignment_map_write_segment(
  builder : StringBuilder,
  sequence : String,
  start : Int,
  end : Int,
  complement_reverse : Bool,
) -> Unit {
  if end >= start {
    let mut index = start
    while index < end {
      builder.write_char(sequence.unsafe_get(index).unsafe_to_char())
      index = index + 1
    }
  } else {
    let mut index = start
    while index > end {
      index = index - 1
      let residue = sequence.unsafe_get(index).unsafe_to_char()
      builder.write_char(
        if complement_reverse {
          alignment_map_complement(residue)
        } else {
          residue
        },
      )
    }
  }
}

///|
/// Render gapped rows when both sequence contents are available.
pub fn CoordinatePairwiseAlignment::aligned_rows(
  self : CoordinatePairwiseAlignment,
  complement_reverse? : Bool = true,
) -> (String, String)? {
  if self.target_sequence.length() == 0 || self.query_sequence.length() == 0 {
    return None
  }
  let target = StringBuilder::new()
  let query = StringBuilder::new()
  let mut index = 0
  while index + 1 < self.target_coordinates.length() {
    let target_start = self.target_coordinates[index]
    let target_end = self.target_coordinates[index + 1]
    let query_start = self.query_coordinates[index]
    let query_end = self.query_coordinates[index + 1]
    let target_step = target_end - target_start
    let query_step = alignment_map_abs(query_end - query_start)
    if target_step > 0 && query_step > 0 {
      alignment_map_write_segment(
        target,
        self.target_sequence,
        target_start,
        target_end,
        false,
      )
      alignment_map_write_segment(
        query,
        self.query_sequence,
        query_start,
        query_end,
        complement_reverse,
      )
    } else if target_step > 0 {
      alignment_map_write_segment(
        target,
        self.target_sequence,
        target_start,
        target_end,
        false,
      )
      alignment_map_repeat(query, '-', target_step)
    } else if query_step > 0 {
      alignment_map_repeat(target, '-', query_step)
      alignment_map_write_segment(
        query,
        self.query_sequence,
        query_start,
        query_end,
        complement_reverse,
      )
    }
    index = index + 1
  }
  Some((target.to_string(), query.to_string()))
}

///|
fn alignment_map_join_ints(values : Array[Int]) -> String {
  let builder = StringBuilder::new()
  let mut index = 0
  while index < values.length() {
    if index > 0 {
      builder.write_char(',')
    }
    builder.write_string(values[index].to_string())
    index = index + 1
  }
  if values.length() > 0 {
    builder.write_char(',')
  }
  builder.to_string()
}

///|
/// Serialize the coordinate blocks as a PSL record.
///
/// Aligned positions are reported as matches because composition does not
/// require sequence contents; gap event/base counts remain exact.
pub fn CoordinatePairwiseAlignment::to_psl(
  self : CoordinatePairwiseAlignment,
) -> String raise AlignmentMapError {
  let blocks = self.blocks()
  if blocks.length() == 0 {
    raise AlignmentMapError("cannot serialize an empty alignment as PSL")
  }
  let counts = self.counts()
  let block_sizes : Array[Int] = []
  let query_starts : Array[Int] = []
  let target_starts : Array[Int] = []
  let reverse = self.is_reverse()
  let mut query_min = self.query_length
  let mut query_max = 0
  let mut target_min = self.target_length
  let mut target_max = 0
  for block in blocks {
    block_sizes.push(block.size)
    target_starts.push(block.target_start)
    let block_query_min = alignment_map_min(block.query_start, block.query_end)
    let block_query_max = alignment_map_max(block.query_start, block.query_end)
    query_min = alignment_map_min(query_min, block_query_min)
    query_max = alignment_map_max(query_max, block_query_max)
    target_min = alignment_map_min(target_min, block.target_start)
    target_max = alignment_map_max(target_max, block.target_end)
    query_starts.push(
      if reverse {
        self.query_length - block_query_max
      } else {
        block_query_min
      },
    )
  }
  let strand = if reverse { "-" } else { "+" }
  counts.aligned.to_string() +
  "\t0\t0\t0\t" +
  counts.target_gap_events.to_string() +
  "\t" +
  counts.target_gap_bases.to_string() +
  "\t" +
  counts.query_gap_events.to_string() +
  "\t" +
  counts.query_gap_bases.to_string() +
  "\t" +
  strand +
  "\t" +
  self.query_name +
  "\t" +
  self.query_length.to_string() +
  "\t" +
  query_min.to_string() +
  "\t" +
  query_max.to_string() +
  "\t" +
  self.target_name +
  "\t" +
  self.target_length.to_string() +
  "\t" +
  target_min.to_string() +
  "\t" +
  target_max.to_string() +
  "\t" +
  blocks.length().to_string() +
  "\t" +
  alignment_map_join_ints(block_sizes) +
  "\t" +
  alignment_map_join_ints(query_starts) +
  "\t" +
  alignment_map_join_ints(target_starts) +
  "\n"
}

///|
pub fn CoordinatePairwiseAlignment::summary(
  self : CoordinatePairwiseAlignment,
) -> String {
  let counts = self.counts()
  "CoordinatePairwiseAlignment(" +
  self.target_name +
  " <- " +
  self.query_name +
  ", blocks=" +
  counts.blocks.to_string() +
  ", aligned=" +
  counts.aligned.to_string() +
  ", strand=" +
  (if self.is_reverse() { "-" } else { "+" }) +
  ")"
}

///|
fn alignment_map_ungap(sequence : String) -> String {
  let builder = StringBuilder::new()
  let mut index = 0
  while index < sequence.length() {
    let residue = sequence.unsafe_get(index).unsafe_to_char()
    if residue != '-' {
      builder.write_char(residue)
    }
    index = index + 1
  }
  builder.to_string()
}

///|
/// Build a coordinate MSA from equal-width gapped rows.
pub fn coordinate_multiple_alignment(
  names : Array[String],
  aligned_sequences : Array[String],
) -> CoordinateMultipleAlignment raise AlignmentMapError {
  if names.length() == 0 {
    raise AlignmentMapError(
      "a multiple alignment must contain at least one row",
    )
  }
  if names.length() != aligned_sequences.length() {
    raise AlignmentMapError("MSA names and rows must have equal length")
  }
  let width = aligned_sequences[0].length()
  let sequences : Array[String] = []
  let coordinates : Array[Array[Int]] = []
  let mut row_index = 0
  while row_index < names.length() {
    alignment_map_validate_name(names[row_index], "MSA row name")
    if aligned_sequences[row_index].length() != width {
      raise AlignmentMapError("all MSA rows must have equal width")
    }
    let raw = alignment_map_ungap(aligned_sequences[row_index])
    alignment_map_validate_sequence(raw, raw.length(), "MSA sequence")
    sequences.push(raw)
    let row = [0]
    let mut coordinate = 0
    let mut column = 0
    while column < width {
      if aligned_sequences[row_index].unsafe_get(column).unsafe_to_char() != '-' {
        coordinate = coordinate + 1
      }
      row.push(coordinate)
      column = column + 1
    }
    coordinates.push(row)
    row_index = row_index + 1
  }
  CoordinateMultipleAlignment::{
    names: alignment_map_copy_strings(names),
    sequences,
    coordinates,
  }
}

///|
fn alignment_map_factor(
  alignment : CoordinatePairwiseAlignment,
) -> Int raise AlignmentMapError {
  let mut factor = 0
  let mut index = 0
  while index + 1 < alignment.target_coordinates.length() {
    let target_step = alignment_map_abs(
      alignment.target_coordinates[index + 1] -
      alignment.target_coordinates[index],
    )
    let query_step = alignment_map_abs(
      alignment.query_coordinates[index + 1] -
      alignment.query_coordinates[index],
    )
    if target_step > 0 && query_step > 0 {
      let current = if query_step == target_step {
        1
      } else if query_step == 3 * target_step {
        3
      } else {
        raise AlignmentMapError(
          "mapall supports only 1:1 or protein-to-codon mappings",
        )
      }
      if factor != 0 && factor != current {
        raise AlignmentMapError("mapping contains inconsistent step factors")
      }
      factor = current
    }
    index = index + 1
  }
  if factor == 0 {
    raise AlignmentMapError("mapping has no aligned blocks")
  }
  factor
}

///|
fn alignment_map_validate_msa(
  alignment : CoordinateMultipleAlignment,
) -> Unit raise AlignmentMapError {
  if alignment.names.length() == 0 ||
    alignment.names.length() != alignment.sequences.length() ||
    alignment.names.length() != alignment.coordinates.length() {
    raise AlignmentMapError("invalid multiple alignment dimensions")
  }
  let points = alignment.coordinates[0].length()
  if points < 2 {
    raise AlignmentMapError(
      "multiple alignment must contain at least one column",
    )
  }
  let mut row_index = 0
  while row_index < alignment.names.length() {
    alignment_map_validate_name(alignment.names[row_index], "MSA row name")
    if alignment.coordinates[row_index].length() != points {
      raise AlignmentMapError("MSA coordinate rows must have equal length")
    }
    let mut point = 0
    while point < points {
      let coordinate = alignment.coordinates[row_index][point]
      if coordinate < 0 || coordinate > alignment.sequences[row_index].length() {
        raise AlignmentMapError("MSA coordinate is out of bounds")
      }
      if point > 0 && coordinate < alignment.coordinates[row_index][point - 1] {
        raise AlignmentMapError("MSA coordinates must be non-decreasing")
      }
      point = point + 1
    }
    row_index = row_index + 1
  }
}

///|
/// Project every row of an MSA through its corresponding sequence mapping.
///
/// A 1:1 mapping produces a mapped sequence MSA. A 1:3 protein-to-nucleotide
/// mapping produces a codon-aware nucleotide MSA without translating or
/// realigning the nucleotide sequences.
pub fn CoordinateMultipleAlignment::mapall(
  self : CoordinateMultipleAlignment,
  mappings : Array[CoordinatePairwiseAlignment],
) -> CoordinateMultipleAlignment raise AlignmentMapError {
  alignment_map_validate_msa(self)
  if mappings.length() != self.sequences.length() {
    raise AlignmentMapError("mapall requires one mapping per MSA row")
  }
  let mut factor = 0
  let mut row_index = 0
  while row_index < mappings.length() {
    if mappings[row_index].target_length != self.sequences[row_index].length() {
      raise AlignmentMapError(
        "mapping target length does not match its MSA row",
      )
    }
    let current_factor = alignment_map_factor(mappings[row_index])
    if factor == 0 {
      factor = current_factor
    } else if factor != current_factor {
      raise AlignmentMapError("mapall mappings use inconsistent step factors")
    }
    row_index = row_index + 1
  }
  let points = self.coordinates[0].length()
  let master_coordinates = [0]
  let mut master = 0
  let mut point = 0
  while point + 1 < points {
    let mut step = 0
    row_index = 0
    while row_index < self.coordinates.length() {
      let current = alignment_map_abs(
        self.coordinates[row_index][point + 1] -
        self.coordinates[row_index][point],
      )
      step = alignment_map_max(step, current)
      row_index = row_index + 1
    }
    master = master + factor * step
    master_coordinates.push(master)
    point = point + 1
  }
  let projected : Array[CoordinatePairwiseAlignment] = []
  row_index = 0
  while row_index < mappings.length() {
    let scaled_row : Array[Int] = []
    for coordinate in self.coordinates[row_index] {
      scaled_row.push(factor * coordinate)
    }
    let outer = alignment_map_create(
      "alignment",
      "",
      master,
      self.names[row_index],
      "",
      factor * self.sequences[row_index].length(),
      master_coordinates,
      scaled_row,
    )
    let scaled_target : Array[Int] = []
    for coordinate in mappings[row_index].target_coordinates {
      scaled_target.push(factor * coordinate)
    }
    let inner = alignment_map_create(
      mappings[row_index].target_name,
      "",
      factor * mappings[row_index].target_length,
      mappings[row_index].query_name,
      mappings[row_index].query_sequence,
      mappings[row_index].query_length,
      scaled_target,
      mappings[row_index].query_coordinates,
    )
    projected.push(outer.map(inner))
    row_index = row_index + 1
  }
  let output_coordinates : Array[Array[Int]] = []
  let indices = Array::make(projected.length(), 0)
  row_index = 0
  while row_index < projected.length() {
    output_coordinates.push([])
    row_index = row_index + 1
  }
  let mut previous = 0
  while true {
    let mut found = false
    let mut position = 2147483647
    row_index = 0
    while row_index < projected.length() {
      let index = indices[row_index]
      if index < projected[row_index].target_coordinates.length() {
        position = alignment_map_min(
          position,
          projected[row_index].target_coordinates[index],
        )
        found = true
      }
      row_index = row_index + 1
    }
    if !found {
      break
    }
    row_index = 0
    while row_index < projected.length() {
      let index = indices[row_index]
      let output = output_coordinates[row_index]
      if index >= projected[row_index].target_coordinates.length() {
        output.push(
          if output.length() == 0 {
            0
          } else {
            output[output.length() - 1]
          },
        )
      } else {
        let current_target = projected[row_index].target_coordinates[index]
        let current_query = projected[row_index].query_coordinates[index]
        if current_target == position {
          output.push(current_query)
          indices[row_index] = index + 1
        } else if current_target > position {
          if output.length() == 0 {
            output.push(current_query)
          } else {
            let last = output[output.length() - 1]
            let step = if current_query > last {
              position - previous
            } else {
              0
            }
            output.push(last + step)
          }
        } else {
          raise AlignmentMapError("mapall coordinate merge became inconsistent")
        }
      }
      row_index = row_index + 1
    }
    previous = position
  }
  let names : Array[String] = []
  let sequences : Array[String] = []
  for alignment in projected {
    names.push(alignment.query_name)
    sequences.push(alignment.query_sequence)
  }
  CoordinateMultipleAlignment::{
    names,
    sequences,
    coordinates: output_coordinates,
  }
}

///|
pub fn CoordinateMultipleAlignment::num_sequences(
  self : CoordinateMultipleAlignment,
) -> Int {
  self.sequences.length()
}

///|
pub fn CoordinateMultipleAlignment::alignment_length(
  self : CoordinateMultipleAlignment,
) -> Int {
  if self.coordinates.length() == 0 || self.coordinates[0].length() < 2 {
    return 0
  }
  let mut length = 0
  let mut point = 0
  while point + 1 < self.coordinates[0].length() {
    let mut step = 0
    for row in self.coordinates {
      step = alignment_map_max(
        step,
        alignment_map_abs(row[point + 1] - row[point]),
      )
    }
    length = length + step
    point = point + 1
  }
  length
}

///|
pub fn CoordinateMultipleAlignment::row(
  self : CoordinateMultipleAlignment,
  index : Int,
  complement_reverse? : Bool = true,
) -> String? {
  if index < 0 ||
    index >= self.sequences.length() ||
    self.coordinates.length() == 0 {
    return None
  }
  let sequence = self.sequences[index]
  if sequence.length() == 0 {
    return None
  }
  let builder = StringBuilder::new()
  let mut point = 0
  while point + 1 < self.coordinates[index].length() {
    let start = self.coordinates[index][point]
    let end = self.coordinates[index][point + 1]
    let row_step = alignment_map_abs(end - start)
    let mut width = 0
    for row_coordinates in self.coordinates {
      width = alignment_map_max(
        width,
        alignment_map_abs(row_coordinates[point + 1] - row_coordinates[point]),
      )
    }
    if row_step == 0 {
      alignment_map_repeat(builder, '-', width)
    } else {
      alignment_map_write_segment(
        builder, sequence, start, end, complement_reverse,
      )
      alignment_map_repeat(builder, '-', width - row_step)
    }
    point = point + 1
  }
  Some(builder.to_string())
}

///|
pub fn CoordinateMultipleAlignment::to_fasta(
  self : CoordinateMultipleAlignment,
) -> String {
  let builder = StringBuilder::new()
  let mut index = 0
  while index < self.names.length() {
    builder.write_char('>')
    builder.write_string(self.names[index])
    builder.write_char('\n')
    match self.row(index) {
      Some(row) => builder.write_string(row)
      None => ()
    }
    builder.write_char('\n')
    index = index + 1
  }
  builder.to_string()
}

///|
pub fn CoordinateMultipleAlignment::summary(
  self : CoordinateMultipleAlignment,
) -> String {
  "CoordinateMultipleAlignment(rows=" +
  self.num_sequences().to_string() +
  ", columns=" +
  self.alignment_length().to_string() +
  ")"
}

///|
/// Small chromosome/transcript/read example from the Biopython tutorial.
pub fn alignment_map_example() -> CoordinatePairwiseAlignment raise AlignmentMapError {
  let chromosome = coordinate_pairwise_alignment(
    "chromosome",
    "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA",
    "transcript",
    "CCCCCCCGGGGGG",
    [8, 15, 26, 32],
    [0, 7, 7, 13],
  )
  let read = coordinate_pairwise_alignment(
    "transcript",
    "CCCCCCCGGGGGG",
    "read",
    "CCCCGGGG",
    [3, 11],
    [0, 8],
  )
  chromosome.map(read)
}
