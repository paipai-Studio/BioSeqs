// Alignment-aware UCSC PSL/PSLX support.
//
// This module follows Biopython Bio.Align.psl coordinate semantics. PSL uses
// zero-based half-open intervals. Nucleotide alignments may reverse the query
// axis; translated DNA-to-protein alignments may reverse the target axis.

///|
pub suberror PslError {
  PslError(String)
}

///|
pub(all) enum PslSequenceType {
  PslNucleotide
  PslProtein
} derive(Eq, Debug)
pub extend PslSequenceType with Eq::{equal, not_equal}
pub extend PslSequenceType with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub(all) enum PslMaskMode {
  PslNoMask
  PslMaskLower
  PslMaskUpper
} derive(Eq, Debug)
pub extend PslMaskMode with Eq::{equal, not_equal}
pub extend PslMaskMode with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct PslBlock {
  target_start : Int
  target_end : Int
  query_start : Int
  query_end : Int
  target_size : Int
  query_size : Int
  query_sequence : String
  target_sequence : String
} derive(Eq, Debug)
pub extend PslBlock with Eq::{equal, not_equal}
pub extend PslBlock with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct PslCounts {
  aligned_query_units : Int
  query_insert_count : Int
  query_insert_bases : Int
  target_insert_count : Int
  target_insert_bases : Int
  block_count : Int
} derive(Eq, Debug)
pub extend PslCounts with Eq::{equal, not_equal}
pub extend PslCounts with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct PslAlignment {
  target_name : String
  target_size : Int
  query_name : String
  query_size : Int
  target_coordinates : Array[Int]
  query_coordinates : Array[Int]
  sequence_type : PslSequenceType
  matches : Int
  mismatches : Int
  repeat_matches : Int
  n_count : Int
  target_sequence : String
  query_sequence : String
  query_block_sequences : Array[String]
  target_block_sequences : Array[String]
}

///|
pub struct PslDocument {
  version : String
  has_header : Bool
  alignments : Array[PslAlignment]
}

///|
pub struct PslWriteConfig {
  header : Bool
  pslx : Bool
  recount : Bool
  mask : PslMaskMode
  wildcard : Char
}

///|
pub struct PslSummary {
  alignment_count : Int
  pslx_count : Int
  nucleotide_count : Int
  protein_count : Int
  aligned_query_units : Int
  query_insert_bases : Int
  target_insert_bases : Int
}

///|
priv struct PslStorage {
  strand : String
  query_start : Int
  query_end : Int
  target_start : Int
  target_end : Int
  block_sizes : Array[Int]
  query_starts : Array[Int]
  target_starts : Array[Int]
}

///|
fn psl_fail(message : String) -> Unit raise PslError {
  raise PslError(message)
}

///|
fn psl_abs(value : Int) -> Int {
  if value < 0 {
    -value
  } else {
    value
  }
}

///|
fn psl_min(left : Int, right : Int) -> Int {
  if left < right {
    left
  } else {
    right
  }
}

///|
fn psl_max(left : Int, right : Int) -> Int {
  if left > right {
    left
  } else {
    right
  }
}

///|
fn psl_copy_ints(values : Array[Int]) -> Array[Int] {
  let copy : Array[Int] = []
  for value in values {
    copy.push(value)
  }
  copy
}

///|
fn psl_copy_strings(values : Array[String]) -> Array[String] {
  let copy : Array[String] = []
  for value in values {
    copy.push(value)
  }
  copy
}

///|
fn psl_reverse_ints(values : Array[Int]) -> Array[Int] {
  let reversed : Array[Int] = []
  let mut index = values.length()
  while index > 0 {
    index = index - 1
    reversed.push(values[index])
  }
  reversed
}

///|
fn psl_strip_cr(value : String) -> String {
  if value.length() > 0 &&
    value.unsafe_get(value.length() - 1).to_int() == '\r'.to_int() {
    value[0:value.length() - 1].to_owned()
  } else {
    value
  }
}

///|
fn psl_validate_text(
  value : String,
  label : String,
  allow_empty : Bool,
) -> Unit raise PslError {
  if !allow_empty && value.length() == 0 {
    psl_fail(label + " must not be empty")
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code == 0 ||
      code == '\t'.to_int() ||
      code == '\n'.to_int() ||
      code == '\r'.to_int() {
      psl_fail(label + " contains a forbidden control character")
    }
  }
}

///|
fn psl_validate_sequence(
  sequence : String,
  expected_size : Int,
  label : String,
) -> Unit raise PslError {
  if sequence.length() == 0 {
    return
  }
  if sequence.length() != expected_size {
    psl_fail(label + " length does not match its declaration")
  }
  for index in 0..<sequence.length() {
    let code = sequence.unsafe_get(index).to_int()
    if code == '-'.to_int() ||
      code == ' '.to_int() ||
      code == '\t'.to_int() ||
      code == '\n'.to_int() ||
      code == '\r'.to_int() {
      psl_fail(label + " must be ungapped and contain no whitespace")
    }
  }
}

///|
fn psl_validate_fragment(
  sequence : String,
  expected_size : Int,
  label : String,
) -> Unit raise PslError {
  if sequence.length() != expected_size {
    psl_fail(label + " length does not match its block size")
  }
  for index in 0..<sequence.length() {
    let code = sequence.unsafe_get(index).to_int()
    if code == '-'.to_int() ||
      code == ','.to_int() ||
      code == ' '.to_int() ||
      code == '\t'.to_int() ||
      code == '\n'.to_int() ||
      code == '\r'.to_int() {
      psl_fail(label + " contains an invalid PSLX residue")
    }
  }
}

///|
fn psl_axis_direction(
  coordinates : Array[Int],
  label : String,
) -> Int raise PslError {
  let mut direction = 0
  for index in 0..<(coordinates.length() - 1) {
    let step = coordinates[index + 1] - coordinates[index]
    if step > 0 {
      if direction < 0 {
        psl_fail(label + " coordinates change direction")
      }
      direction = 1
    } else if step < 0 {
      if direction > 0 {
        psl_fail(label + " coordinates change direction")
      }
      direction = -1
    }
  }
  direction
}

///|
fn psl_raw_blocks(
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
  query_sequences : Array[String],
  target_sequences : Array[String],
) -> Array[PslBlock] {
  let blocks : Array[PslBlock] = []
  let mut block_index = 0
  for index in 0..<(target_coordinates.length() - 1) {
    let target_start = target_coordinates[index]
    let target_end = target_coordinates[index + 1]
    let query_start = query_coordinates[index]
    let query_end = query_coordinates[index + 1]
    let target_size = psl_abs(target_end - target_start)
    let query_size = psl_abs(query_end - query_start)
    if target_size > 0 && query_size > 0 {
      blocks.push(PslBlock::{
        target_start,
        target_end,
        query_start,
        query_end,
        target_size,
        query_size,
        query_sequence: if query_sequences.length() > 0 {
          query_sequences[block_index]
        } else {
          ""
        },
        target_sequence: if target_sequences.length() > 0 {
          target_sequences[block_index]
        } else {
          ""
        },
      })
      block_index = block_index + 1
    }
  }
  blocks
}

///|
fn psl_path_counts(
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> PslCounts {
  let mut aligned_query_units = 0
  let mut query_insert_count = 0
  let mut query_insert_bases = 0
  let mut target_insert_count = 0
  let mut target_insert_bases = 0
  let mut block_count = 0
  for index in 0..<(target_coordinates.length() - 1) {
    let target_step = psl_abs(
      target_coordinates[index + 1] - target_coordinates[index],
    )
    let query_step = psl_abs(
      query_coordinates[index + 1] - query_coordinates[index],
    )
    if target_step > 0 && query_step > 0 {
      aligned_query_units = aligned_query_units + query_step
      block_count = block_count + 1
    } else if target_step == 0 && query_step > 0 {
      query_insert_count = query_insert_count + 1
      query_insert_bases = query_insert_bases + query_step
    } else if target_step > 0 && query_step == 0 {
      target_insert_count = target_insert_count + 1
      target_insert_bases = target_insert_bases + target_step
    }
  }
  PslCounts::{
    aligned_query_units,
    query_insert_count,
    query_insert_bases,
    target_insert_count,
    target_insert_bases,
    block_count,
  }
}

///|
pub fn PslAlignment::create(
  target_name : String,
  target_size : Int,
  query_name : String,
  query_size : Int,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
  sequence_type? : PslSequenceType = PslNucleotide,
  matches? : Int = -1,
  mismatches? : Int = -1,
  repeat_matches? : Int = -1,
  n_count? : Int = -1,
  target_sequence? : String = "",
  query_sequence? : String = "",
  query_block_sequences? : Array[String] = [],
  target_block_sequences? : Array[String] = [],
) -> PslAlignment raise PslError {
  psl_validate_text(target_name, "PSL target name", false)
  psl_validate_text(query_name, "PSL query name", false)
  if target_size <= 0 || query_size <= 0 {
    psl_fail("PSL target and query sizes must be positive")
  }
  if target_coordinates.length() != query_coordinates.length() {
    psl_fail("PSL target and query coordinate arrays must have equal length")
  }
  if target_coordinates.length() < 2 {
    psl_fail("PSL alignment path must contain at least two points")
  }
  let mut targets = psl_copy_ints(target_coordinates)
  let mut queries = psl_copy_ints(query_coordinates)
  let initial_target_direction = psl_axis_direction(targets, "PSL target")
  let initial_query_direction = psl_axis_direction(queries, "PSL query")
  match sequence_type {
    PslNucleotide =>
      if initial_target_direction < 0 {
        targets = psl_reverse_ints(targets)
        queries = psl_reverse_ints(queries)
      }
    PslProtein =>
      if initial_query_direction < 0 {
        targets = psl_reverse_ints(targets)
        queries = psl_reverse_ints(queries)
      }
  }
  let target_direction = psl_axis_direction(targets, "PSL target")
  let query_direction = psl_axis_direction(queries, "PSL query")
  match sequence_type {
    PslNucleotide =>
      if target_direction <= 0 || query_direction == 0 {
        psl_fail(
          "nucleotide PSL requires an increasing target and an oriented query",
        )
      }
    PslProtein =>
      if query_direction <= 0 || target_direction == 0 {
        psl_fail(
          "translated PSL requires an increasing query and an oriented target",
        )
      }
  }
  for index in 0..<targets.length() {
    if targets[index] < 0 || targets[index] > target_size {
      psl_fail("PSL target coordinate is out of bounds")
    }
    if queries[index] < 0 || queries[index] > query_size {
      psl_fail("PSL query coordinate is out of bounds")
    }
  }
  let mut block_count = 0
  for index in 0..<(targets.length() - 1) {
    let target_step = psl_abs(targets[index + 1] - targets[index])
    let query_step = psl_abs(queries[index + 1] - queries[index])
    if target_step == 0 && query_step == 0 {
      psl_fail("PSL path contains a zero-length step")
    }
    if target_step > 0 && query_step > 0 {
      match sequence_type {
        PslNucleotide =>
          if target_step != query_step {
            psl_fail("nucleotide PSL aligned steps must have equal lengths")
          }
        PslProtein =>
          if target_step != 3 * query_step {
            psl_fail("translated PSL aligned steps must have a 3:1 ratio")
          }
      }
      block_count = block_count + 1
    }
  }
  if block_count == 0 {
    psl_fail("PSL alignment must contain at least one aligned block")
  }
  let first_target_step = psl_abs(targets[1] - targets[0])
  let first_query_step = psl_abs(queries[1] - queries[0])
  let last = targets.length() - 1
  let last_target_step = psl_abs(targets[last] - targets[last - 1])
  let last_query_step = psl_abs(queries[last] - queries[last - 1])
  if first_target_step == 0 ||
    first_query_step == 0 ||
    last_target_step == 0 ||
    last_query_step == 0 {
    psl_fail("PSL path cannot begin or end with an unaligned gap")
  }
  psl_validate_sequence(target_sequence, target_size, "PSL target sequence")
  psl_validate_sequence(query_sequence, query_size, "PSL query sequence")
  let has_query_blocks = query_block_sequences.length() > 0
  let has_target_blocks = target_block_sequences.length() > 0
  if has_query_blocks != has_target_blocks {
    psl_fail("PSLX requires both query and target block sequences")
  }
  if has_query_blocks {
    if query_block_sequences.length() != block_count ||
      target_block_sequences.length() != block_count {
      psl_fail("PSLX block sequence count does not match blockCount")
    }
    let blocks = psl_raw_blocks(targets, queries, [], [])
    for index in 0..<blocks.length() {
      psl_validate_fragment(
        query_block_sequences[index],
        blocks[index].query_size,
        "PSLX query block",
      )
      psl_validate_fragment(
        target_block_sequences[index],
        blocks[index].query_size,
        "PSLX target block",
      )
    }
  }
  let counts = psl_path_counts(targets, queries)
  let all_missing = matches < 0 &&
    mismatches < 0 &&
    repeat_matches < 0 &&
    n_count < 0
  let normalized_matches = if all_missing {
    counts.aligned_query_units
  } else {
    matches
  }
  let normalized_mismatches = if all_missing { 0 } else { mismatches }
  let normalized_repeat = if all_missing { 0 } else { repeat_matches }
  let normalized_n = if all_missing { 0 } else { n_count }
  if normalized_matches < 0 ||
    normalized_mismatches < 0 ||
    normalized_repeat < 0 ||
    normalized_n < 0 {
    psl_fail("PSL match categories must all be supplied or all be omitted")
  }
  if normalized_matches +
    normalized_mismatches +
    normalized_repeat +
    normalized_n !=
    counts.aligned_query_units {
    psl_fail("PSL match categories must sum to aligned query units")
  }
  PslAlignment::{
    target_name,
    target_size,
    query_name,
    query_size,
    target_coordinates: targets,
    query_coordinates: queries,
    sequence_type,
    matches: normalized_matches,
    mismatches: normalized_mismatches,
    repeat_matches: normalized_repeat,
    n_count: normalized_n,
    target_sequence,
    query_sequence,
    query_block_sequences: psl_copy_strings(query_block_sequences),
    target_block_sequences: psl_copy_strings(target_block_sequences),
  }
}

///|
pub fn PslWriteConfig::create(
  header? : Bool = true,
  pslx? : Bool = false,
  recount? : Bool = false,
  mask? : PslMaskMode = PslNoMask,
  wildcard? : Char = 'N',
) -> PslWriteConfig raise PslError {
  let code = wildcard.to_int()
  if code == 0 ||
    code == ','.to_int() ||
    code == '\t'.to_int() ||
    code == '\n'.to_int() ||
    code == '\r'.to_int() {
    psl_fail("PSL wildcard must be a printable non-delimiter character")
  }
  PslWriteConfig::{ header, pslx, recount, mask, wildcard, }
}

///|
pub fn PslWriteConfig::default() -> PslWriteConfig {
  PslWriteConfig::{
    header: true,
    pslx: false,
    recount: false,
    mask: PslNoMask,
    wildcard: 'N',
  }
}

///|
pub fn PslAlignment::blocks(self : PslAlignment) -> Array[PslBlock] {
  psl_raw_blocks(
    self.target_coordinates,
    self.query_coordinates,
    self.query_block_sequences,
    self.target_block_sequences,
  )
}

///|
pub fn PslAlignment::counts(self : PslAlignment) -> PslCounts {
  psl_path_counts(self.target_coordinates, self.query_coordinates)
}

///|
pub fn PslAlignment::is_reverse(self : PslAlignment) -> Bool {
  match self.sequence_type {
    PslNucleotide =>
      self.query_coordinates[0] >
      self.query_coordinates[self.query_coordinates.length() - 1]
    PslProtein =>
      self.target_coordinates[0] >
      self.target_coordinates[self.target_coordinates.length() - 1]
  }
}

///|
pub fn PslAlignment::is_pslx(self : PslAlignment) -> Bool {
  self.query_block_sequences.length() > 0
}

///|
pub fn PslAlignment::identity(self : PslAlignment) -> Double {
  let total = self.matches +
    self.mismatches +
    self.repeat_matches +
    self.n_count
  if total == 0 {
    0.0
  } else {
    (self.matches + self.repeat_matches).to_double() / total.to_double()
  }
}

///|
pub fn PslAlignment::score(self : PslAlignment) -> Int {
  let counts = self.counts()
  self.matches +
  self.repeat_matches / 2 -
  self.mismatches -
  counts.query_insert_count -
  counts.target_insert_count
}

///|
pub fn PslAlignment::target_to_query(
  self : PslAlignment,
  position : Int,
) -> Int? {
  if position < 0 || position >= self.target_size {
    return None
  }
  for block in self.blocks() {
    let target_minimum = psl_min(block.target_start, block.target_end)
    let target_maximum = psl_max(block.target_start, block.target_end)
    if position >= target_minimum && position < target_maximum {
      let target_offset = if block.target_end > block.target_start {
        position - block.target_start
      } else {
        block.target_start - 1 - position
      }
      let query_offset = match self.sequence_type {
        PslNucleotide => target_offset
        PslProtein => target_offset / 3
      }
      return Some(
        if block.query_end > block.query_start {
          block.query_start + query_offset
        } else {
          block.query_start - 1 - query_offset
        },
      )
    }
  }
  None
}

///|
pub fn PslAlignment::query_to_target_interval(
  self : PslAlignment,
  position : Int,
) -> (Int, Int)? {
  if position < 0 || position >= self.query_size {
    return None
  }
  for block in self.blocks() {
    let query_minimum = psl_min(block.query_start, block.query_end)
    let query_maximum = psl_max(block.query_start, block.query_end)
    if position >= query_minimum && position < query_maximum {
      let query_offset = if block.query_end > block.query_start {
        position - block.query_start
      } else {
        block.query_start - 1 - position
      }
      let scale = match self.sequence_type {
        PslNucleotide => 1
        PslProtein => 3
      }
      let oriented = query_offset * scale
      let first = if block.target_end > block.target_start {
        block.target_start + oriented
      } else {
        block.target_start - oriented - scale
      }
      return Some(
        (psl_min(first, first + scale), psl_max(first, first + scale)),
      )
    }
  }
  None
}

///|
pub fn PslAlignment::query_to_target(
  self : PslAlignment,
  position : Int,
) -> Int? {
  match self.query_to_target_interval(position) {
    Some((start, _)) => Some(start)
    None => None
  }
}

///|
fn psl_upper_code(code : Int) -> Int {
  if code >= 'a'.to_int() && code <= 'z'.to_int() {
    code - 32
  } else {
    code
  }
}

///|
fn psl_is_lower(code : Int) -> Bool {
  code >= 'a'.to_int() && code <= 'z'.to_int()
}

///|
fn psl_is_upper(code : Int) -> Bool {
  code >= 'A'.to_int() && code <= 'Z'.to_int()
}

///|
fn psl_complement_code(code : Int) -> Int {
  match psl_upper_code(code) {
    65 => 84
    67 => 71
    71 => 67
    84 => 65
    85 => 65
    82 => 89
    89 => 82
    77 => 75
    75 => 77
    66 => 86
    86 => 66
    68 => 72
    72 => 68
    value => value
  }
}

///|
fn psl_translated_block(
  sequence : String,
  start : Int,
  end : Int,
) -> String raise PslError {
  let minimum = psl_min(start, end)
  let maximum = psl_max(start, end)
  let dna = sequence[minimum:maximum].to_owned()
  let oriented = if end < start {
    Seq::new(dna).reverse_complement()
  } else {
    Seq::new(dna)
  }
  let protein = oriented.translate() catch {
    SeqError(message) => {
      psl_fail("cannot translate PSL target block: " + message)
      Seq::new("")
    }
  }
  protein.to_string()
}

///|
pub fn PslAlignment::recount(
  self : PslAlignment,
  target_sequence? : String = "",
  query_sequence? : String = "",
  mask? : PslMaskMode = PslNoMask,
  wildcard? : Char = 'N',
) -> PslAlignment raise PslError {
  let target = if target_sequence.length() > 0 {
    target_sequence
  } else {
    self.target_sequence
  }
  let query = if query_sequence.length() > 0 {
    query_sequence
  } else {
    self.query_sequence
  }
  psl_validate_sequence(target, self.target_size, "PSL target sequence")
  psl_validate_sequence(query, self.query_size, "PSL query sequence")
  if target.length() == 0 || query.length() == 0 {
    psl_fail("PSL recount requires concrete target and query sequences")
  }
  if self.sequence_type == PslProtein && mask != PslNoMask {
    psl_fail("repeat masking is only defined for nucleotide PSL")
  }
  let wildcard_code = psl_upper_code(wildcard.to_int())
  let mut matches = 0
  let mut mismatches = 0
  let mut repeat_matches = 0
  let mut n_count = 0
  for block in self.blocks() {
    match self.sequence_type {
      PslProtein => {
        let translated = psl_translated_block(
          target,
          block.target_start,
          block.target_end,
        )
        if translated.length() != block.query_size {
          psl_fail("translated PSL target block has an unexpected length")
        }
        for offset in 0..<block.query_size {
          let query_index = block.query_start + offset
          let target_code = psl_upper_code(
            translated.unsafe_get(offset).to_int(),
          )
          let query_code = psl_upper_code(
            query.unsafe_get(query_index).to_int(),
          )
          if target_code == wildcard_code || query_code == wildcard_code {
            n_count = n_count + 1
          } else if target_code == query_code {
            matches = matches + 1
          } else {
            mismatches = mismatches + 1
          }
        }
      }
      PslNucleotide =>
        for offset in 0..<block.query_size {
          let target_index = block.target_start + offset
          let query_index = if block.query_end > block.query_start {
            block.query_start + offset
          } else {
            block.query_start - 1 - offset
          }
          let raw_target = target.unsafe_get(target_index).to_int()
          let target_code = psl_upper_code(raw_target)
          let raw_query = query.unsafe_get(query_index).to_int()
          let query_code = if block.query_end < block.query_start {
            psl_complement_code(raw_query)
          } else {
            psl_upper_code(raw_query)
          }
          if target_code == wildcard_code || query_code == wildcard_code {
            n_count = n_count + 1
          } else if target_code != query_code {
            mismatches = mismatches + 1
          } else {
            let masked = match mask {
              PslMaskLower => psl_is_lower(raw_target)
              PslMaskUpper => psl_is_upper(raw_target)
              PslNoMask => false
            }
            if masked {
              repeat_matches = repeat_matches + 1
            } else {
              matches = matches + 1
            }
          }
        }
    }
  }
  PslAlignment::create(
    self.target_name,
    self.target_size,
    self.query_name,
    self.query_size,
    self.target_coordinates,
    self.query_coordinates,
    sequence_type=self.sequence_type,
    matches~,
    mismatches~,
    repeat_matches~,
    n_count~,
    target_sequence=target,
    query_sequence=query,
  )
}

///|
fn psl_storage(alignment : PslAlignment) -> PslStorage {
  let blocks = alignment.blocks()
  let reverse = alignment.is_reverse()
  let block_sizes : Array[Int] = []
  let query_starts : Array[Int] = []
  let target_starts : Array[Int] = []
  let mut query_start = alignment.query_size
  let mut query_end = 0
  let mut target_start = alignment.target_size
  let mut target_end = 0
  for block in blocks {
    let query_minimum = psl_min(block.query_start, block.query_end)
    let query_maximum = psl_max(block.query_start, block.query_end)
    let target_minimum = psl_min(block.target_start, block.target_end)
    let target_maximum = psl_max(block.target_start, block.target_end)
    query_start = psl_min(query_start, query_minimum)
    query_end = psl_max(query_end, query_maximum)
    target_start = psl_min(target_start, target_minimum)
    target_end = psl_max(target_end, target_maximum)
    block_sizes.push(block.query_size)
    query_starts.push(
      if alignment.sequence_type == PslNucleotide && reverse {
        alignment.query_size - query_maximum
      } else {
        query_minimum
      },
    )
    target_starts.push(
      if alignment.sequence_type == PslProtein && reverse {
        alignment.target_size - target_maximum
      } else {
        target_minimum
      },
    )
  }
  let strand = match alignment.sequence_type {
    PslNucleotide => if reverse { "-" } else { "+" }
    PslProtein => if reverse { "+-" } else { "++" }
  }
  PslStorage::{
    strand,
    query_start,
    query_end,
    target_start,
    target_end,
    block_sizes,
    query_starts,
    target_starts,
  }
}

///|
fn psl_join_ints(values : Array[Int]) -> String {
  let output = StringBuilder::new()
  for index in 0..<values.length() {
    if index > 0 {
      output.write_char(',')
    }
    output.write_string(values[index].to_string())
  }
  if values.length() > 0 {
    output.write_char(',')
  }
  output.to_string()
}

///|
fn psl_join_strings(values : Array[String]) -> String {
  let output = StringBuilder::new()
  for index in 0..<values.length() {
    if index > 0 {
      output.write_char(',')
    }
    output.write_string(values[index])
  }
  if values.length() > 0 {
    output.write_char(',')
  }
  output.to_string()
}

///|
fn psl_block_sequences(
  alignment : PslAlignment,
) -> (Array[String], Array[String]) raise PslError {
  if alignment.query_block_sequences.length() > 0 {
    return (
      psl_copy_strings(alignment.query_block_sequences),
      psl_copy_strings(alignment.target_block_sequences),
    )
  }
  if alignment.target_sequence.length() == 0 ||
    alignment.query_sequence.length() == 0 {
    psl_fail(
      "PSLX output requires block sequences or concrete target/query sequences",
    )
  }
  let query_sequences : Array[String] = []
  let target_sequences : Array[String] = []
  for block in alignment.blocks() {
    match alignment.sequence_type {
      PslNucleotide => {
        let target = alignment.target_sequence[block.target_start:block.target_end].to_owned()
        let query_minimum = psl_min(block.query_start, block.query_end)
        let query_maximum = psl_max(block.query_start, block.query_end)
        let raw_query = alignment.query_sequence[query_minimum:query_maximum].to_owned()
        let query = if block.query_end < block.query_start {
          Seq::new(raw_query).reverse_complement().to_string()
        } else {
          raw_query
        }
        target_sequences.push(target)
        query_sequences.push(query)
      }
      PslProtein => {
        let target = psl_translated_block(
          alignment.target_sequence,
          block.target_start,
          block.target_end,
        )
        let query = alignment.query_sequence[block.query_start:block.query_end].to_owned()
        target_sequences.push(target)
        query_sequences.push(query)
      }
    }
  }
  (query_sequences, target_sequences)
}

///|
pub fn psl_header(version? : String = "3") -> String raise PslError {
  psl_validate_text(version, "PSL version", false)
  for index in 0..<version.length() {
    if version.unsafe_get(index).to_int() == ' '.to_int() {
      psl_fail("PSL version must not contain whitespace")
    }
  }
  "psLayout version " +
  version +
  "\n\n" +
  "match\tmis-\trep.\tN's\tQ gap\tQ gap\tT gap\tT gap\tstrand\tQ name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\n" +
  "\tmatch\tmatch\t\tcount\tbases\tcount\tbases\n" +
  "---------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
}

///|
pub fn PslAlignment::format(
  self : PslAlignment,
  config? : PslWriteConfig = PslWriteConfig::default(),
) -> String raise PslError {
  let alignment = if config.recount {
    self.recount(mask=config.mask, wildcard=config.wildcard)
  } else {
    self
  }
  let storage = psl_storage(alignment)
  let counts = alignment.counts()
  let output = StringBuilder::new()
  let fields = [
    alignment.matches.to_string(),
    alignment.mismatches.to_string(),
    alignment.repeat_matches.to_string(),
    alignment.n_count.to_string(),
    counts.query_insert_count.to_string(),
    counts.query_insert_bases.to_string(),
    counts.target_insert_count.to_string(),
    counts.target_insert_bases.to_string(),
    storage.strand,
    alignment.query_name,
    alignment.query_size.to_string(),
    storage.query_start.to_string(),
    storage.query_end.to_string(),
    alignment.target_name,
    alignment.target_size.to_string(),
    storage.target_start.to_string(),
    storage.target_end.to_string(),
    counts.block_count.to_string(),
    psl_join_ints(storage.block_sizes),
    psl_join_ints(storage.query_starts),
    psl_join_ints(storage.target_starts),
  ]
  for index in 0..<fields.length() {
    if index > 0 {
      output.write_char('\t')
    }
    output.write_string(fields[index])
  }
  if config.pslx {
    let (query_sequences, target_sequences) = psl_block_sequences(alignment)
    output.write_char('\t')
    output.write_string(psl_join_strings(query_sequences))
    output.write_char('\t')
    output.write_string(psl_join_strings(target_sequences))
  }
  output.write_char('\n')
  output.to_string()
}

///|
pub fn psl_write(
  alignments : Array[PslAlignment],
  config? : PslWriteConfig = PslWriteConfig::default(),
  version? : String = "3",
) -> String raise PslError {
  let output = StringBuilder::new()
  if config.header {
    output.write_string(psl_header(version~))
  }
  for alignment in alignments {
    output.write_string(alignment.format(config~))
  }
  output.to_string()
}

///|
fn psl_parse_int(value : String, label : String) -> Int raise PslError {
  if value.length() == 0 {
    psl_fail(label + " must not be empty")
  }
  let mut result = 0
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      psl_fail(label + " must be a non-negative integer")
    }
    let digit = code - '0'.to_int()
    if result > (2147483647 - digit) / 10 {
      psl_fail(label + " exceeds the supported integer range")
    }
    result = result * 10 + digit
  }
  result
}

///|
fn psl_parse_csv_ints(
  value : String,
  label : String,
) -> Array[Int] raise PslError {
  let parts = split_by_char(value, ','.to_int())
  let values : Array[Int] = []
  for index in 0..<parts.length() {
    if parts[index].length() == 0 {
      if index + 1 == parts.length() {
        continue
      }
      psl_fail(label + " contains an empty element")
    }
    values.push(psl_parse_int(parts[index], label + " element"))
  }
  if values.length() == 0 {
    psl_fail(label + " must contain at least one value")
  }
  values
}

///|
fn psl_parse_csv_strings(
  value : String,
  label : String,
) -> Array[String] raise PslError {
  let parts = split_by_char(value, ','.to_int())
  let values : Array[String] = []
  for index in 0..<parts.length() {
    if parts[index].length() == 0 {
      if index + 1 == parts.length() {
        continue
      }
      psl_fail(label + " contains an empty element")
    }
    values.push(parts[index])
  }
  if values.length() == 0 {
    psl_fail(label + " must contain at least one value")
  }
  values
}

///|
fn psl_parse_record(
  line : String,
  line_number : Int,
) -> PslAlignment raise PslError {
  let words = split_by_char(line, '\t'.to_int())
  if words.length() != 21 && words.length() != 23 {
    psl_fail(
      "PSL line " +
      line_number.to_string() +
      " has " +
      words.length().to_string() +
      " columns; expected 21 or 23",
    )
  }
  let matches = psl_parse_int(words[0], "PSL matches")
  let mismatches = psl_parse_int(words[1], "PSL misMatches")
  let repeat_matches = psl_parse_int(words[2], "PSL repMatches")
  let n_count = psl_parse_int(words[3], "PSL nCount")
  let declared_query_insert_count = psl_parse_int(words[4], "PSL qNumInsert")
  let declared_query_insert_bases = psl_parse_int(words[5], "PSL qBaseInsert")
  let declared_target_insert_count = psl_parse_int(words[6], "PSL tNumInsert")
  let declared_target_insert_bases = psl_parse_int(words[7], "PSL tBaseInsert")
  let strand = words[8]
  let sequence_type = match strand {
    "+" | "-" => PslNucleotide
    "++" | "+-" => PslProtein
    _ => {
      psl_fail("PSL strand must be '+', '-', '++', or '+-'")
      PslNucleotide
    }
  }
  let query_name = words[9]
  let query_size = psl_parse_int(words[10], "PSL qSize")
  let declared_query_start = psl_parse_int(words[11], "PSL qStart")
  let declared_query_end = psl_parse_int(words[12], "PSL qEnd")
  let target_name = words[13]
  let target_size = psl_parse_int(words[14], "PSL tSize")
  let declared_target_start = psl_parse_int(words[15], "PSL tStart")
  let declared_target_end = psl_parse_int(words[16], "PSL tEnd")
  let block_count = psl_parse_int(words[17], "PSL blockCount")
  if block_count <= 0 {
    psl_fail("PSL blockCount must be positive")
  }
  let block_sizes = psl_parse_csv_ints(words[18], "PSL blockSizes")
  let query_starts = psl_parse_csv_ints(words[19], "PSL qStarts")
  let target_starts = psl_parse_csv_ints(words[20], "PSL tStarts")
  if block_sizes.length() != block_count ||
    query_starts.length() != block_count ||
    target_starts.length() != block_count {
    psl_fail("PSL block arrays do not match blockCount")
  }
  let query_block_sequences : Array[String] = if words.length() == 23 {
    psl_parse_csv_strings(words[21], "PSLX query sequences")
  } else {
    []
  }
  let target_block_sequences : Array[String] = if words.length() == 23 {
    psl_parse_csv_strings(words[22], "PSLX target sequences")
  } else {
    []
  }
  if words.length() == 23 &&
    (
      query_block_sequences.length() != block_count ||
      target_block_sequences.length() != block_count
    ) {
    psl_fail("PSLX sequence arrays do not match blockCount")
  }
  let target_coordinates : Array[Int] = [target_starts[0]]
  let query_coordinates : Array[Int] = [query_starts[0]]
  let mut target_position = target_starts[0]
  let mut query_position = query_starts[0]
  for index in 0..<block_count {
    let block_size = block_sizes[index]
    if block_size <= 0 {
      psl_fail("PSL block sizes must be positive")
    }
    let target_block_size = match sequence_type {
      PslNucleotide => block_size
      PslProtein => 3 * block_size
    }
    let target_start = target_starts[index]
    let query_start = query_starts[index]
    if target_start < target_position || query_start < query_position {
      psl_fail("PSL blocks overlap or are not sorted")
    }
    if target_start + target_block_size > target_size ||
      query_start + block_size > query_size {
      psl_fail("PSL block extends beyond a declared sequence size")
    }
    if target_start != target_position {
      target_coordinates.push(target_start)
      query_coordinates.push(query_position)
      target_position = target_start
    }
    if query_start != query_position {
      target_coordinates.push(target_position)
      query_coordinates.push(query_start)
      query_position = query_start
    }
    target_position = target_position + target_block_size
    query_position = query_position + block_size
    target_coordinates.push(target_position)
    query_coordinates.push(query_position)
  }
  if strand == "-" {
    for index in 0..<query_coordinates.length() {
      query_coordinates[index] = query_size - query_coordinates[index]
    }
  } else if strand == "+-" {
    for index in 0..<target_coordinates.length() {
      target_coordinates[index] = target_size - target_coordinates[index]
    }
  }
  let alignment = PslAlignment::create(
    target_name,
    target_size,
    query_name,
    query_size,
    target_coordinates,
    query_coordinates,
    sequence_type~,
    matches~,
    mismatches~,
    repeat_matches~,
    n_count~,
    query_block_sequences~,
    target_block_sequences~,
  )
  let counts = alignment.counts()
  if counts.query_insert_count != declared_query_insert_count ||
    counts.query_insert_bases != declared_query_insert_bases {
    psl_fail("PSL query insertion statistics are inconsistent with blocks")
  }
  if counts.target_insert_count != declared_target_insert_count ||
    counts.target_insert_bases != declared_target_insert_bases {
    psl_fail("PSL target insertion statistics are inconsistent with blocks")
  }
  let storage = psl_storage(alignment)
  if storage.strand != strand {
    psl_fail("PSL strand is inconsistent with its coordinate path")
  }
  if storage.query_start != declared_query_start ||
    storage.query_end != declared_query_end {
    psl_fail("PSL qStart/qEnd are inconsistent with blocks")
  }
  if storage.target_start != declared_target_start ||
    storage.target_end != declared_target_end {
    psl_fail("PSL tStart/tEnd are inconsistent with blocks")
  }
  alignment
}

///|
fn psl_is_separator(line : String) -> Bool {
  if line.length() == 0 {
    return false
  }
  for index in 0..<line.length() {
    if line.unsafe_get(index).to_int() != '-'.to_int() {
      return false
    }
  }
  true
}

///|
pub fn psl_parse(content : String) -> PslDocument raise PslError {
  let raw_lines = split_by_char(content, '\n'.to_int())
  let lines : Array[String] = []
  for line in raw_lines {
    lines.push(psl_strip_cr(line))
  }
  let mut index = 0
  while index < lines.length() && lines[index].length() == 0 {
    index = index + 1
  }
  if index >= lines.length() {
    psl_fail("PSL input is empty")
  }
  let mut version = ""
  let mut has_header = false
  if starts_with(lines[index], 0, "psLayout ") {
    let words = split_by_whitespace(lines[index])
    if words.length() != 3 || words[1] != "version" {
      psl_fail("malformed PSL header")
    }
    version = words[2]
    psl_validate_text(version, "PSL version", false)
    has_header = true
    index = index + 1
    let mut found_separator = false
    while index < lines.length() {
      if psl_is_separator(lines[index]) {
        found_separator = true
        index = index + 1
        break
      }
      index = index + 1
    }
    if !found_separator {
      psl_fail("PSL header separator was not found")
    }
  }
  let alignments : Array[PslAlignment] = []
  while index < lines.length() {
    if lines[index].length() > 0 {
      alignments.push(psl_parse_record(lines[index], index + 1))
    }
    index = index + 1
  }
  PslDocument::{ version, has_header, alignments, }
}

///|
pub fn PslDocument::query(
  self : PslDocument,
  name : String,
) -> Array[PslAlignment] {
  let results : Array[PslAlignment] = []
  for alignment in self.alignments {
    if alignment.query_name == name {
      results.push(alignment)
    }
  }
  results
}

///|
pub fn PslDocument::target(
  self : PslDocument,
  name : String,
) -> Array[PslAlignment] {
  let results : Array[PslAlignment] = []
  for alignment in self.alignments {
    if alignment.target_name == name {
      results.push(alignment)
    }
  }
  results
}

///|
pub fn PslDocument::summary(self : PslDocument) -> PslSummary {
  let mut pslx_count = 0
  let mut nucleotide_count = 0
  let mut protein_count = 0
  let mut aligned_query_units = 0
  let mut query_insert_bases = 0
  let mut target_insert_bases = 0
  for alignment in self.alignments {
    if alignment.is_pslx() {
      pslx_count = pslx_count + 1
    }
    match alignment.sequence_type {
      PslNucleotide => nucleotide_count = nucleotide_count + 1
      PslProtein => protein_count = protein_count + 1
    }
    let counts = alignment.counts()
    aligned_query_units = aligned_query_units + counts.aligned_query_units
    query_insert_bases = query_insert_bases + counts.query_insert_bases
    target_insert_bases = target_insert_bases + counts.target_insert_bases
  }
  PslSummary::{
    alignment_count: self.alignments.length(),
    pslx_count,
    nucleotide_count,
    protein_count,
    aligned_query_units,
    query_insert_bases,
    target_insert_bases,
  }
}

///|
pub fn PslDocument::write(
  self : PslDocument,
  config? : PslWriteConfig = PslWriteConfig::default(),
) -> String raise PslError {
  let version = if self.version.length() > 0 { self.version } else { "3" }
  psl_write(self.alignments, config~, version~)
}

///|
pub fn PslAlignment::to_coordinate_alignment(
  self : PslAlignment,
) -> CoordinatePairwiseAlignment raise PslError {
  if self.sequence_type == PslProtein && self.is_reverse() {
    psl_fail(
      "reverse-target translated PSL cannot use the target-increasing coordinate adapter",
    )
  }
  coordinate_pairwise_alignment_with_lengths(
    self.target_name,
    self.target_size,
    self.query_name,
    self.query_size,
    self.target_coordinates,
    self.query_coordinates,
    target_sequence=self.target_sequence,
    query_sequence=self.query_sequence,
  ) catch {
    AlignmentMapError(message) => {
      psl_fail("cannot convert PSL coordinate path: " + message)
      coordinate_pairwise_alignment_with_lengths(
        "target",
        1,
        "query",
        1,
        [0, 1],
        [0, 1],
      ) catch {
        _ => abort("unreachable PSL coordinate fallback")
      }
    }
  }
}

///|
pub fn psl_from_coordinate_alignment(
  alignment : CoordinatePairwiseAlignment,
  sequence_type? : PslSequenceType = PslNucleotide,
) -> PslAlignment raise PslError {
  PslAlignment::create(
    alignment.target_name,
    alignment.target_length,
    alignment.query_name,
    alignment.query_length,
    alignment.target_coordinates,
    alignment.query_coordinates,
    sequence_type~,
    target_sequence=alignment.target_sequence,
    query_sequence=alignment.query_sequence,
  )
}

///|
pub fn psl_example_data() -> Array[PslAlignment] raise PslError {
  let nucleotide = PslAlignment::create(
    "chrExample",
    24,
    "readForward",
    12,
    [2, 6, 9, 9, 13],
    [0, 4, 4, 6, 10],
    target_sequence="TTACGTGGGACGTCCCAAATTTGG",
    query_sequence="ACGTGGACGTAA",
    matches=6,
    mismatches=1,
    repeat_matches=0,
    n_count=1,
  )
  let reverse = PslAlignment::create(
    "chrExample",
    24,
    "readReverse",
    8,
    [4, 8, 10, 14],
    [8, 4, 4, 0],
  )
  let protein = PslAlignment::create(
    "codingDna",
    12,
    "peptide",
    2,
    [0, 6],
    [0, 2],
    sequence_type=PslProtein,
    target_sequence="ATGGCTTAATAG",
    query_sequence="MA",
  ).recount()
  [nucleotide, reverse, protein]
}
