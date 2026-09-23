// Bio.Align.bigpsl - indexed UCSC BigPsl pairwise alignments.
//
// BigPsl is a bed12+13 BigBed specialization. The BigBed container owns the
// binary layout and indices; this module owns pairwise-coordinate semantics,
// translated DNA-to-protein blocks, match categories, and strict schema
// validation.

///|
pub suberror BigPslError {
  BigPslError(String)
}

///|
pub(all) enum BigPslSequenceType {
  BigPslUnknown
  BigPslNucleotide
  BigPslProtein
} derive(Eq, Debug)
pub extend BigPslSequenceType with Eq::{equal, not_equal}
pub extend BigPslSequenceType with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub(all) enum BigPslMaskMode {
  BigPslNoMask
  BigPslMaskLower
  BigPslMaskUpper
} derive(Eq, Debug)
pub extend BigPslMaskMode with Eq::{equal, not_equal}
pub extend BigPslMaskMode with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigPslBlock {
  target_start : Int
  target_end : Int
  query_start : Int
  query_end : Int
  target_size : Int
  query_size : Int
} derive(Eq, Debug)
pub extend BigPslBlock with Eq::{equal, not_equal}
pub extend BigPslBlock with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigPslAlignment {
  target_name : String
  target_size : Int
  query_name : String
  query_size : Int
  target_coordinates : Array[Int]
  query_coordinates : Array[Int]
  query_sequence : String
  cds : String
  score : Int
  thick_start : Int
  thick_end : Int
  item_rgb : String
  matches : Int
  mismatches : Int
  repeat_matches : Int
  n_count : Int
  sequence_type : BigPslSequenceType
} derive(Eq, Debug)
pub extend BigPslAlignment with Eq::{equal, not_equal}
pub extend BigPslAlignment with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigPslWriteConfig {
  compress : Bool
  block_size : Int
  items_per_slot : Int
  store_query_sequence : Bool
  store_cds : Bool
} derive(Eq, Debug)
pub extend BigPslWriteConfig with Eq::{equal, not_equal}
pub extend BigPslWriteConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigPslAlignmentCounts {
  aligned_query_units : Int
  query_insert_units : Int
  target_insert_bases : Int
  query_insert_events : Int
  target_insert_events : Int
  block_count : Int
} derive(Eq, Debug)
pub extend BigPslAlignmentCounts with Eq::{equal, not_equal}
pub extend BigPslAlignmentCounts with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigPslSummary {
  alignment_count : Int
  target_count : Int
  nucleotide_count : Int
  protein_count : Int
  unknown_count : Int
  aligned_query_units : Int
  covered_target_bases : Int
  compressed : Bool
} derive(Eq, Debug)
pub extend BigPslSummary with Eq::{equal, not_equal}
pub extend BigPslSummary with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct BigPslFile {
  bed : BigBedFile
}

///|
priv struct BigPslStorage {
  chrom_start : Int
  chrom_end : Int
  strand : String
  other_start : Int
  other_end : Int
  other_strand : String
  target_block_sizes : Array[Int]
  target_block_starts : Array[Int]
  other_block_starts : Array[Int]
}

///|
fn bigpsl_fail(message : String) -> Unit raise BigPslError {
  raise BigPslError(message)
}

///|
fn bigpsl_abs(value : Int) -> Int {
  if value < 0 {
    -value
  } else {
    value
  }
}

///|
fn bigpsl_min(left : Int, right : Int) -> Int {
  if left < right {
    left
  } else {
    right
  }
}

///|
fn bigpsl_max(left : Int, right : Int) -> Int {
  if left > right {
    left
  } else {
    right
  }
}

///|
fn bigpsl_copy_ints(values : Array[Int]) -> Array[Int] {
  let copy : Array[Int] = []
  for value in values {
    copy.push(value)
  }
  copy
}

///|
fn bigpsl_reverse_ints(values : Array[Int]) -> Array[Int] {
  let reversed : Array[Int] = []
  let mut index = values.length()
  while index > 0 {
    index = index - 1
    reversed.push(values[index])
  }
  reversed
}

///|
fn bigpsl_validate_text(
  value : String,
  label : String,
  allow_empty : Bool,
) -> Unit raise BigPslError {
  if !allow_empty && value.length() == 0 {
    bigpsl_fail(label + " must not be empty")
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code == 0 ||
      code == '\t'.to_int() ||
      code == '\n'.to_int() ||
      code == '\r'.to_int() {
      bigpsl_fail(label + " contains a forbidden control character")
    }
  }
}

///|
fn bigpsl_validate_sequence(
  sequence : String,
  expected_size : Int,
  label : String,
) -> Unit raise BigPslError {
  if sequence.length() == 0 {
    return
  }
  if sequence.length() != expected_size {
    bigpsl_fail(label + " length does not match its declaration")
  }
  for index in 0..<sequence.length() {
    let code = sequence.unsafe_get(index).to_int()
    if code == '-'.to_int() ||
      code == ' '.to_int() ||
      code == '\t'.to_int() ||
      code == '\n'.to_int() ||
      code == '\r'.to_int() {
      bigpsl_fail(label + " must be an ungapped sequence without whitespace")
    }
  }
}

///|
fn bigpsl_type_code(sequence_type : BigPslSequenceType) -> Int {
  match sequence_type {
    BigPslUnknown => 0
    BigPslNucleotide => 1
    BigPslProtein => 2
  }
}

///|
fn bigpsl_type_from_code(code : Int) -> BigPslSequenceType raise BigPslError {
  match code {
    0 => BigPslUnknown
    1 => BigPslNucleotide
    2 => BigPslProtein
    _ => {
      bigpsl_fail("BigPsl seqType must be 0, 1, or 2")
      BigPslUnknown
    }
  }
}

///|
fn bigpsl_axis_direction(
  coordinates : Array[Int],
  label : String,
) -> Int raise BigPslError {
  let mut direction = 0
  for index in 0..<(coordinates.length() - 1) {
    let step = coordinates[index + 1] - coordinates[index]
    if step > 0 {
      if direction < 0 {
        bigpsl_fail(label + " coordinates change direction")
      }
      direction = 1
    } else if step < 0 {
      if direction > 0 {
        bigpsl_fail(label + " coordinates change direction")
      }
      direction = -1
    }
  }
  direction
}

///|
fn bigpsl_validate_coordinates(
  target_size : Int,
  query_size : Int,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
  sequence_type : BigPslSequenceType,
) -> Unit raise BigPslError {
  if target_size <= 0 || query_size <= 0 {
    bigpsl_fail("BigPsl target and query sizes must be positive")
  }
  if target_coordinates.length() != query_coordinates.length() {
    bigpsl_fail(
      "BigPsl target and query coordinate arrays must have equal length",
    )
  }
  if target_coordinates.length() < 2 {
    bigpsl_fail("BigPsl alignment paths require at least two coordinate points")
  }
  for index in 0..<target_coordinates.length() {
    if target_coordinates[index] < 0 || target_coordinates[index] > target_size {
      bigpsl_fail("BigPsl target coordinate is out of bounds")
    }
    if query_coordinates[index] < 0 || query_coordinates[index] > query_size {
      bigpsl_fail("BigPsl query coordinate is out of bounds")
    }
  }
  let target_direction = bigpsl_axis_direction(
    target_coordinates, "BigPsl target",
  )
  let query_direction = bigpsl_axis_direction(query_coordinates, "BigPsl query")
  let mut aligned_blocks = 0
  for index in 0..<(target_coordinates.length() - 1) {
    let target_step = target_coordinates[index + 1] - target_coordinates[index]
    let query_step = query_coordinates[index + 1] - query_coordinates[index]
    if target_step == 0 && query_step == 0 {
      bigpsl_fail("Consecutive BigPsl coordinate points must differ")
    }
    if target_step != 0 && query_step != 0 {
      let target_count = bigpsl_abs(target_step)
      let query_count = bigpsl_abs(query_step)
      match sequence_type {
        BigPslProtein =>
          if target_count != 3 * query_count {
            bigpsl_fail(
              "Protein BigPsl aligned target blocks must be three times their query blocks",
            )
          }
        _ =>
          if target_count != query_count {
            bigpsl_fail(
              "Nucleotide BigPsl aligned target and query blocks must have equal size",
            )
          }
      }
      aligned_blocks = aligned_blocks + 1
    }
  }
  if aligned_blocks == 0 {
    bigpsl_fail("BigPsl alignment must contain at least one aligned block")
  }
  match sequence_type {
    BigPslProtein => {
      if query_direction != 1 {
        bigpsl_fail("Protein BigPsl query coordinates must increase")
      }
      if target_direction == 0 {
        bigpsl_fail("Protein BigPsl target coordinates must have a direction")
      }
    }
    _ => {
      if target_direction != 1 {
        bigpsl_fail("Nucleotide BigPsl target coordinates must increase")
      }
      if query_direction == 0 {
        bigpsl_fail("Nucleotide BigPsl query coordinates must have a direction")
      }
    }
  }
}

///|
fn bigpsl_raw_blocks(
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> Array[BigPslBlock] {
  let blocks : Array[BigPslBlock] = []
  for index in 0..<(target_coordinates.length() - 1) {
    let target_start = target_coordinates[index]
    let target_end = target_coordinates[index + 1]
    let query_start = query_coordinates[index]
    let query_end = query_coordinates[index + 1]
    let target_count = bigpsl_abs(target_end - target_start)
    let query_count = bigpsl_abs(query_end - query_start)
    if target_count > 0 && query_count > 0 {
      blocks.push(BigPslBlock::{
        target_start,
        target_end,
        query_start,
        query_end,
        target_size: target_count,
        query_size: query_count,
      })
    }
  }
  blocks
}

///|
fn bigpsl_aligned_query_units(
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> Int {
  let mut total = 0
  for block in bigpsl_raw_blocks(target_coordinates, query_coordinates) {
    total = total + block.query_size
  }
  total
}

///|
fn bigpsl_coordinate_interval(coordinates : Array[Int]) -> (Int, Int) {
  let mut start = coordinates[0]
  let mut end = coordinates[0]
  for coordinate in coordinates {
    start = bigpsl_min(start, coordinate)
    end = bigpsl_max(end, coordinate)
  }
  (start, end)
}

///|
pub fn BigPslAlignment::create(
  target_name : String,
  target_size : Int,
  query_name : String,
  query_size : Int,
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
  sequence_type? : BigPslSequenceType = BigPslUnknown,
  query_sequence? : String = "",
  cds? : String = "",
  score? : Int = 0,
  thick_start? : Int = -1,
  thick_end? : Int = -1,
  item_rgb? : String = "0",
  matches? : Int = -1,
  mismatches? : Int = -1,
  repeat_matches? : Int = -1,
  n_count? : Int = -1,
) -> BigPslAlignment raise BigPslError {
  bigpsl_validate_text(target_name, "BigPsl target name", false)
  bigpsl_validate_text(query_name, "BigPsl query name", false)
  bigpsl_validate_text(item_rgb, "BigPsl itemRgb", false)
  bigpsl_validate_text(cds, "BigPsl CDS", true)
  bigpsl_validate_sequence(query_sequence, query_size, "BigPsl query sequence")
  bigpsl_validate_coordinates(
    target_size, query_size, target_coordinates, query_coordinates, sequence_type,
  )
  if score < 0 || score > 1000 {
    bigpsl_fail("BigPsl score must be between 0 and 1000")
  }
  let (target_minimum, target_maximum) = bigpsl_coordinate_interval(
    target_coordinates,
  )
  let normalized_thick_start = if thick_start < 0 {
    target_minimum
  } else {
    thick_start
  }
  let normalized_thick_end = if thick_end < 0 {
    target_maximum
  } else {
    thick_end
  }
  if normalized_thick_start < target_minimum ||
    normalized_thick_start > normalized_thick_end ||
    normalized_thick_end > target_maximum {
    bigpsl_fail("BigPsl thick interval must lie inside its target interval")
  }
  let aligned = bigpsl_aligned_query_units(
    target_coordinates, query_coordinates,
  )
  let all_unspecified = matches < 0 &&
    mismatches < 0 &&
    repeat_matches < 0 &&
    n_count < 0
  let normalized_matches = if all_unspecified { aligned } else { matches }
  let normalized_mismatches = if all_unspecified { 0 } else { mismatches }
  let normalized_repeat = if all_unspecified { 0 } else { repeat_matches }
  let normalized_n = if all_unspecified { 0 } else { n_count }
  if normalized_matches < 0 ||
    normalized_mismatches < 0 ||
    normalized_repeat < 0 ||
    normalized_n < 0 {
    bigpsl_fail(
      "BigPsl match categories must all be supplied or all be omitted",
    )
  }
  if normalized_matches +
    normalized_mismatches +
    normalized_repeat +
    normalized_n !=
    aligned {
    bigpsl_fail("BigPsl match categories must sum to the aligned query length")
  }
  BigPslAlignment::{
    target_name,
    target_size,
    query_name,
    query_size,
    target_coordinates: bigpsl_copy_ints(target_coordinates),
    query_coordinates: bigpsl_copy_ints(query_coordinates),
    query_sequence,
    cds,
    score,
    thick_start: normalized_thick_start,
    thick_end: normalized_thick_end,
    item_rgb,
    matches: normalized_matches,
    mismatches: normalized_mismatches,
    repeat_matches: normalized_repeat,
    n_count: normalized_n,
    sequence_type,
  }
}

///|
pub fn BigPslWriteConfig::create(
  compress? : Bool = true,
  block_size? : Int = 256,
  items_per_slot? : Int = 512,
  store_query_sequence? : Bool = false,
  store_cds? : Bool = false,
) -> BigPslWriteConfig raise BigPslError {
  if block_size < 2 || block_size > 65535 {
    bigpsl_fail("BigPsl blockSize must be between 2 and 65535")
  }
  if items_per_slot <= 0 || items_per_slot > 65535 {
    bigpsl_fail("BigPsl itemsPerSlot must be between 1 and 65535")
  }
  BigPslWriteConfig::{
    compress,
    block_size,
    items_per_slot,
    store_query_sequence,
    store_cds,
  }
}

///|
pub fn BigPslWriteConfig::default() -> BigPslWriteConfig {
  BigPslWriteConfig::{
    compress: true,
    block_size: 256,
    items_per_slot: 512,
    store_query_sequence: false,
    store_cds: false,
  }
}

///|
pub fn BigPslAlignment::blocks(self : BigPslAlignment) -> Array[BigPslBlock] {
  bigpsl_raw_blocks(self.target_coordinates, self.query_coordinates)
}

///|
pub fn BigPslAlignment::target_interval(self : BigPslAlignment) -> (Int, Int) {
  bigpsl_coordinate_interval(self.target_coordinates)
}

///|
pub fn BigPslAlignment::query_interval(self : BigPslAlignment) -> (Int, Int) {
  bigpsl_coordinate_interval(self.query_coordinates)
}

///|
pub fn BigPslAlignment::is_reverse(self : BigPslAlignment) -> Bool {
  let last = self.target_coordinates.length() - 1
  self.target_coordinates[0] > self.target_coordinates[last] ||
  self.query_coordinates[0] > self.query_coordinates[last]
}

///|
pub fn BigPslAlignment::counts(self : BigPslAlignment) -> BigPslAlignmentCounts {
  let mut aligned_query_units = 0
  let mut query_insert_units = 0
  let mut target_insert_bases = 0
  let mut query_insert_events = 0
  let mut target_insert_events = 0
  let mut block_count = 0
  for index in 0..<(self.target_coordinates.length() - 1) {
    let target_step = bigpsl_abs(
      self.target_coordinates[index + 1] - self.target_coordinates[index],
    )
    let query_step = bigpsl_abs(
      self.query_coordinates[index + 1] - self.query_coordinates[index],
    )
    if target_step > 0 && query_step > 0 {
      aligned_query_units = aligned_query_units + query_step
      block_count = block_count + 1
    } else if target_step == 0 {
      query_insert_units = query_insert_units + query_step
      query_insert_events = query_insert_events + 1
    } else {
      target_insert_bases = target_insert_bases + target_step
      target_insert_events = target_insert_events + 1
    }
  }
  BigPslAlignmentCounts::{
    aligned_query_units,
    query_insert_units,
    target_insert_bases,
    query_insert_events,
    target_insert_events,
    block_count,
  }
}

///|
pub fn BigPslAlignment::target_to_query(
  self : BigPslAlignment,
  position : Int,
) -> Int? {
  if position < 0 || position >= self.target_size {
    return None
  }
  for block in self.blocks() {
    let minimum = bigpsl_min(block.target_start, block.target_end)
    let maximum = bigpsl_max(block.target_start, block.target_end)
    if position >= minimum && position < maximum {
      let offset = if block.target_end > block.target_start {
        position - block.target_start
      } else {
        block.target_start - 1 - position
      }
      let query_offset = match self.sequence_type {
        BigPslProtein => offset / 3
        _ => offset
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
pub fn BigPslAlignment::query_to_target_interval(
  self : BigPslAlignment,
  position : Int,
) -> (Int, Int)? {
  if position < 0 || position >= self.query_size {
    return None
  }
  for block in self.blocks() {
    let minimum = bigpsl_min(block.query_start, block.query_end)
    let maximum = bigpsl_max(block.query_start, block.query_end)
    if position >= minimum && position < maximum {
      let offset = if block.query_end > block.query_start {
        position - block.query_start
      } else {
        block.query_start - 1 - position
      }
      let scale = match self.sequence_type {
        BigPslProtein => 3
        _ => 1
      }
      let oriented = offset * scale
      let first = if block.target_end > block.target_start {
        block.target_start + oriented
      } else {
        block.target_start - oriented - scale
      }
      return Some(
        (bigpsl_min(first, first + scale), bigpsl_max(first, first + scale)),
      )
    }
  }
  None
}

///|
fn bigpsl_upper_code(code : Int) -> Int {
  if code >= 'a'.to_int() && code <= 'z'.to_int() {
    code - 32
  } else {
    code
  }
}

///|
fn bigpsl_is_lower(code : Int) -> Bool {
  code >= 'a'.to_int() && code <= 'z'.to_int()
}

///|
fn bigpsl_is_upper(code : Int) -> Bool {
  code >= 'A'.to_int() && code <= 'Z'.to_int()
}

///|
fn bigpsl_complement_code(code : Int) -> Int {
  match bigpsl_upper_code(code) {
    65 => 84 // A -> T
    67 => 71 // C -> G
    71 => 67 // G -> C
    84 => 65 // T -> A
    85 => 65 // U -> A
    82 => 89 // R -> Y
    89 => 82 // Y -> R
    77 => 75 // M -> K
    75 => 77 // K -> M
    66 => 86 // B -> V
    86 => 66 // V -> B
    68 => 72 // D -> H
    72 => 68 // H -> D
    value => value
  }
}

///|
pub fn BigPslAlignment::recount(
  self : BigPslAlignment,
  target_sequence : String,
  query_sequence : String,
  mask? : BigPslMaskMode = BigPslNoMask,
  wildcard? : Char = 'N',
) -> BigPslAlignment raise BigPslError {
  bigpsl_validate_sequence(
    target_sequence,
    self.target_size,
    "BigPsl target sequence",
  )
  bigpsl_validate_sequence(
    query_sequence,
    self.query_size,
    "BigPsl query sequence",
  )
  if target_sequence.length() == 0 || query_sequence.length() == 0 {
    bigpsl_fail("BigPsl recount requires concrete target and query sequences")
  }
  if self.sequence_type == BigPslProtein && mask != BigPslNoMask {
    bigpsl_fail("Repeat masking is only defined for nucleotide BigPsl records")
  }
  let wildcard_code = bigpsl_upper_code(wildcard.to_int())
  let mut matches = 0
  let mut mismatches = 0
  let mut repeat_matches = 0
  let mut n_count = 0
  for block in self.blocks() {
    match self.sequence_type {
      BigPslProtein => {
        let start = bigpsl_min(block.target_start, block.target_end)
        let end = bigpsl_max(block.target_start, block.target_end)
        let dna = target_sequence[start:end].to_owned()
        let oriented = if block.target_end < block.target_start {
          Seq::new(dna).reverse_complement()
        } else {
          Seq::new(dna)
        }
        let translated = oriented.translate() catch {
          SeqError(message) => {
            bigpsl_fail("Cannot translate BigPsl target block: " + message)
            Seq::new("")
          }
        }
        let amino = translated.to_string()
        if amino.length() != block.query_size {
          bigpsl_fail("Translated BigPsl target block has an unexpected length")
        }
        for offset in 0..<block.query_size {
          let query_index = block.query_start + offset
          let target_code = bigpsl_upper_code(amino.unsafe_get(offset).to_int())
          let query_code = bigpsl_upper_code(
            query_sequence.unsafe_get(query_index).to_int(),
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
      _ =>
        for offset in 0..<block.target_size {
          let target_index = block.target_start + offset
          let query_index = if block.query_end > block.query_start {
            block.query_start + offset
          } else {
            block.query_start - 1 - offset
          }
          let raw_target = target_sequence.unsafe_get(target_index).to_int()
          let target_code = bigpsl_upper_code(raw_target)
          let raw_query = query_sequence.unsafe_get(query_index).to_int()
          let query_code = if block.query_end < block.query_start {
            bigpsl_complement_code(raw_query)
          } else {
            bigpsl_upper_code(raw_query)
          }
          if target_code == wildcard_code || query_code == wildcard_code {
            n_count = n_count + 1
          } else if target_code != query_code {
            mismatches = mismatches + 1
          } else {
            let masked = match mask {
              BigPslMaskLower => bigpsl_is_lower(raw_target)
              BigPslMaskUpper => bigpsl_is_upper(raw_target)
              BigPslNoMask => false
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
  BigPslAlignment::create(
    self.target_name,
    self.target_size,
    self.query_name,
    self.query_size,
    self.target_coordinates,
    self.query_coordinates,
    sequence_type=self.sequence_type,
    query_sequence~,
    cds=self.cds,
    score=self.score,
    thick_start=self.thick_start,
    thick_end=self.thick_end,
    item_rgb=self.item_rgb,
    matches~,
    mismatches~,
    repeat_matches~,
    n_count~,
  )
}

///|
fn bigpsl_storage(alignment : BigPslAlignment) -> BigPslStorage {
  let blocks = alignment.blocks()
  let target_reverse = alignment.target_coordinates[0] >
    alignment.target_coordinates[alignment.target_coordinates.length() - 1]
  let query_reverse = alignment.query_coordinates[0] >
    alignment.query_coordinates[alignment.query_coordinates.length() - 1]
  let strand = if target_reverse || query_reverse { "-" } else { "+" }
  let other_strand = if target_reverse { "-" } else { "+" }
  let ordered : Array[BigPslBlock] = []
  if target_reverse {
    let mut index = blocks.length()
    while index > 0 {
      index = index - 1
      ordered.push(blocks[index])
    }
  } else {
    for block in blocks {
      ordered.push(block)
    }
  }
  let mut chrom_start = alignment.target_size
  let mut chrom_end = 0
  let mut other_start = alignment.query_size
  let mut other_end = 0
  for block in ordered {
    chrom_start = bigpsl_min(
      chrom_start,
      bigpsl_min(block.target_start, block.target_end),
    )
    chrom_end = bigpsl_max(
      chrom_end,
      bigpsl_max(block.target_start, block.target_end),
    )
    other_start = bigpsl_min(
      other_start,
      bigpsl_min(block.query_start, block.query_end),
    )
    other_end = bigpsl_max(
      other_end,
      bigpsl_max(block.query_start, block.query_end),
    )
  }
  let target_block_sizes : Array[Int] = []
  let target_block_starts : Array[Int] = []
  let other_block_starts : Array[Int] = []
  for block in ordered {
    let target_minimum = bigpsl_min(block.target_start, block.target_end)
    let query_minimum = bigpsl_min(block.query_start, block.query_end)
    let query_maximum = bigpsl_max(block.query_start, block.query_end)
    target_block_sizes.push(block.target_size)
    target_block_starts.push(target_minimum - chrom_start)
    other_block_starts.push(
      if strand == "-" {
        alignment.query_size - query_maximum
      } else {
        query_minimum
      },
    )
  }
  BigPslStorage::{
    chrom_start,
    chrom_end,
    strand,
    other_start,
    other_end,
    other_strand,
    target_block_sizes,
    target_block_starts,
    other_block_starts,
  }
}

///|
fn bigpsl_join_ints(values : Array[Int], trailing : Bool) -> String {
  let output = StringBuilder()
  for index in 0..<values.length() {
    if index > 0 {
      output.write_char(',')
    }
    output.write_string(values[index].to_string())
  }
  if trailing && values.length() > 0 {
    output.write_char(',')
  }
  output.to_string()
}

///|
pub fn bigpsl_schema() -> BigBedSchema {
  let standard_names = [
    "chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd",
    "reserved", "blockCount", "blockSizes", "chromStarts",
  ]
  let standard_types = [
    "string", "uint", "uint", "string", "uint", "char[1]", "uint", "uint", "uint",
    "int", "int[blockCount]", "int[blockCount]",
  ]
  let standard_comments = [
    "Reference sequence chromosome or scaffold", "Start position in chromosome",
    "End position in chromosome", "Name or ID of item, ideally both human readable and unique",
    "Score (0-1000)", "+ or - indicates whether the query aligns to the + or - strand on the reference",
    "Start of where display should be thick (start codon)", "End of where display should be thick (stop codon)",
    "RGB value (use R,G,B string in input file)", "Number of blocks", "Comma separated list of block sizes",
    "Start positions relative to chromStart",
  ]
  let custom_names = [
    "oChromStart", "oChromEnd", "oStrand", "oChromSize", "oChromStarts", "oSequence",
    "oCDS", "chromSize", "match", "misMatch", "repMatch", "nCount", "seqType",
  ]
  let custom_types = [
    "uint", "uint", "char[1]", "uint", "int[blockCount]", "lstring", "string", "uint",
    "uint", "uint", "uint", "uint", "uint",
  ]
  let custom_comments = [
    "Start position in other chromosome", "End position in other chromosome", "+ or -, - means that psl was reversed into BED-compatible coordinates",
    "Size of other chromosome.", "Start positions relative to oChromStart or from oChromStart+oChromSize depending on strand",
    "Sequence on other chrom (or edit list, or empty)", "CDS in NCBI format", "Size of target chromosome",
    "Number of bases matched.", "Number of bases that don't match", "Number of bases that match but are part of repeats",
    "Number of 'N' bases", "0=empty, 1=nucleotide, 2=amino_acid",
  ]
  let fields : Array[BigBedField] = []
  for index in 0..<standard_names.length() {
    fields.push(BigBedField::{
      as_type: standard_types[index],
      name: standard_names[index],
      comment: standard_comments[index],
    })
  }
  for index in 0..<custom_names.length() {
    fields.push(BigBedField::{
      as_type: custom_types[index],
      name: custom_names[index],
      comment: custom_comments[index],
    })
  }
  BigBedSchema::{
    name: "bigPsl",
    comment: "bigPsl pairwise alignment",
    fields,
  }
}

///|
fn bigpsl_custom_fields() -> Array[BigBedField] {
  let schema = bigpsl_schema()
  let fields : Array[BigBedField] = []
  for index in 12..<schema.fields.length() {
    fields.push(schema.fields[index])
  }
  fields
}

///|
fn bigpsl_target_index(targets : Array[BigBedTarget], name : String) -> Int {
  for index in 0..<targets.length() {
    if targets[index].name == name {
      return index
    }
  }
  -1
}

///|
fn bigpsl_validated_copy(
  alignment : BigPslAlignment,
) -> BigPslAlignment raise BigPslError {
  BigPslAlignment::create(
    alignment.target_name,
    alignment.target_size,
    alignment.query_name,
    alignment.query_size,
    alignment.target_coordinates,
    alignment.query_coordinates,
    sequence_type=alignment.sequence_type,
    query_sequence=alignment.query_sequence,
    cds=alignment.cds,
    score=alignment.score,
    thick_start=alignment.thick_start,
    thick_end=alignment.thick_end,
    item_rgb=alignment.item_rgb,
    matches=alignment.matches,
    mismatches=alignment.mismatches,
    repeat_matches=alignment.repeat_matches,
    n_count=alignment.n_count,
  )
}

///|
pub fn bigpsl_write(
  targets : Array[BigBedTarget],
  alignments : Array[BigPslAlignment],
  config? : BigPslWriteConfig = BigPslWriteConfig::default(),
) -> Array[Int] raise BigPslError {
  if targets.length() == 0 {
    bigpsl_fail("BigPsl requires at least one target")
  }
  for index in 0..<targets.length() {
    ignore(BigBedTarget::create(targets[index].name, targets[index].length)) catch {
      BigBedError(message) => bigpsl_fail("BigBed: " + message)
    }
    for previous in 0..<index {
      if targets[previous].name == targets[index].name {
        bigpsl_fail("Duplicate BigPsl target '" + targets[index].name + "'")
      }
    }
  }
  let records : Array[BigBedRecord] = []
  let target_order : Array[Int] = []
  for source in alignments {
    let alignment = bigpsl_validated_copy(source)
    let target_index = bigpsl_target_index(targets, alignment.target_name)
    if target_index < 0 {
      bigpsl_fail(
        "BigPsl alignment refers to unknown target '" +
        alignment.target_name +
        "'",
      )
    }
    if targets[target_index].length != alignment.target_size {
      bigpsl_fail("BigPsl alignment target size disagrees with its target")
    }
    let storage = bigpsl_storage(alignment)
    let query_sequence = if config.store_query_sequence {
      alignment.query_sequence
    } else {
      ""
    }
    if config.store_query_sequence && query_sequence.length() == 0 {
      bigpsl_fail("BigPsl store_query_sequence requires every query sequence")
    }
    let cds = if config.store_cds {
      if alignment.cds.length() == 0 {
        "n/a"
      } else {
        alignment.cds
      }
    } else {
      ""
    }
    let extra = [
      storage.other_start.to_string(),
      storage.other_end.to_string(),
      storage.other_strand,
      alignment.query_size.to_string(),
      bigpsl_join_ints(storage.other_block_starts, false),
      query_sequence,
      cds,
      alignment.target_size.to_string(),
      alignment.matches.to_string(),
      alignment.mismatches.to_string(),
      alignment.repeat_matches.to_string(),
      alignment.n_count.to_string(),
      bigpsl_type_code(alignment.sequence_type).to_string(),
    ]
    let record = BigBedRecord::create(
      alignment.target_name,
      storage.chrom_start,
      storage.chrom_end,
      name=alignment.query_name,
      score=alignment.score,
      strand=storage.strand,
      thick_start=alignment.thick_start,
      thick_end=alignment.thick_end,
      item_rgb=alignment.item_rgb,
      block_count=storage.target_block_sizes.length(),
      block_sizes=storage.target_block_sizes,
      block_starts=storage.target_block_starts,
      extra_fields=extra,
    ) catch {
      BigBedError(message) => {
        bigpsl_fail("BigBed: " + message)
        BigBedRecord::{
          chrom: "",
          chrom_start: 0,
          chrom_end: 1,
          name: ".",
          score: 0,
          strand: ".",
          thick_start: 0,
          thick_end: 1,
          item_rgb: "0",
          declared_block_count: 0,
          block_sizes: [],
          block_starts: [],
          extra_fields: [],
        }
      }
    }
    records.push(record)
    target_order.push(target_index)
  }
  let order : Array[Int] = []
  for index in 0..<records.length() {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if target_order[left] < target_order[right] {
      -1
    } else if target_order[left] > target_order[right] {
      1
    } else if records[left].chrom_start < records[right].chrom_start {
      -1
    } else if records[left].chrom_start > records[right].chrom_start {
      1
    } else if records[left].chrom_end < records[right].chrom_end {
      -1
    } else if records[left].chrom_end > records[right].chrom_end {
      1
    } else {
      0
    }
  })
  let sorted : Array[BigBedRecord] = []
  for index in order {
    sorted.push(records[index])
  }
  let bed_config = BigBedWriteConfig::create(
    bed_columns=12,
    items_per_slot=config.items_per_slot,
    block_size=config.block_size,
    compress=config.compress,
    schema_name="bigPsl",
    schema_comment="bigPsl pairwise alignment",
    custom_fields=bigpsl_custom_fields(),
  ) catch {
    BigBedError(message) => {
      bigpsl_fail("BigBed: " + message)
      BigBedWriteConfig::default()
    }
  }
  bigbed_write(targets, sorted, config=bed_config) catch {
    BigBedError(message) => {
      bigpsl_fail("BigBed: " + message)
      []
    }
  }
}

///|
fn bigpsl_parse_int(value : String, label : String) -> Int raise BigPslError {
  if value.length() == 0 {
    bigpsl_fail("Missing integer for " + label)
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      bigpsl_fail("Invalid integer for " + label + ": " + value)
    }
  }
  parse_int(value)
}

///|
fn bigpsl_split_csv(
  value : String,
  label : String,
) -> Array[Int] raise BigPslError {
  if value.length() == 0 {
    return []
  }
  let values : Array[Int] = []
  let mut start = 0
  for index in 0..<=value.length() {
    if index == value.length() ||
      value.unsafe_get(index).to_int() == ','.to_int() {
      if index == start {
        if index == value.length() && values.length() > 0 {
          return values
        }
        bigpsl_fail("Empty value in " + label)
      }
      values.push(bigpsl_parse_int(value[start:index].to_owned(), label))
      start = index + 1
    }
  }
  values
}

///|
fn bigpsl_validate_schema(bed : BigBedFile) -> Unit raise BigPslError {
  let expected = bigpsl_schema()
  if bed.header.defined_field_count != 12 ||
    bed.header.field_count != 25 ||
    bed.schema.name != "bigPsl" ||
    bed.schema.fields.length() != 25 {
    bigpsl_fail("BigPsl requires the standard bigPsl bed12+13 declaration")
  }
  for index in 0..<expected.fields.length() {
    if bed.schema.fields[index].name != expected.fields[index].name ||
      bed.schema.fields[index].as_type != expected.fields[index].as_type {
      bigpsl_fail(
        "BigPsl AutoSQL field " +
        index.to_string() +
        " does not match the standard declaration",
      )
    }
  }
}

///|
fn bigpsl_decode_record(
  bed : BigBedFile,
  record : BigBedRecord,
) -> BigPslAlignment raise BigPslError {
  if record.extra_fields.length() != 13 {
    bigpsl_fail("BigPsl binary records require thirteen custom fields")
  }
  if record.strand != "+" && record.strand != "-" {
    bigpsl_fail("BigPsl BED strand must be '+' or '-'")
  }
  let other_start = bigpsl_parse_int(
    record.extra_fields[0],
    "BigPsl oChromStart",
  )
  let other_end = bigpsl_parse_int(record.extra_fields[1], "BigPsl oChromEnd")
  let other_strand = record.extra_fields[2]
  if other_strand != "+" && other_strand != "-" {
    bigpsl_fail("BigPsl oStrand must be '+' or '-'")
  }
  if other_strand == "-" && record.strand != "-" {
    bigpsl_fail("BigPsl oStrand '-' requires BED strand '-'")
  }
  let query_size = bigpsl_parse_int(record.extra_fields[3], "BigPsl oChromSize")
  let stored_query_starts = bigpsl_split_csv(
    record.extra_fields[4],
    "BigPsl oChromStarts",
  )
  let query_sequence = record.extra_fields[5]
  let cds = if record.extra_fields[6] == "n/a" {
    ""
  } else {
    record.extra_fields[6]
  }
  let target_size = bigpsl_parse_int(record.extra_fields[7], "BigPsl chromSize")
  let matches = bigpsl_parse_int(record.extra_fields[8], "BigPsl match")
  let mismatches = bigpsl_parse_int(record.extra_fields[9], "BigPsl misMatch")
  let repeat_matches = bigpsl_parse_int(
    record.extra_fields[10],
    "BigPsl repMatch",
  )
  let n_count = bigpsl_parse_int(record.extra_fields[11], "BigPsl nCount")
  let sequence_type = bigpsl_type_from_code(
    bigpsl_parse_int(record.extra_fields[12], "BigPsl seqType"),
  )
  match bed.target(record.chrom) {
    Some(target) =>
      if target.length != target_size {
        bigpsl_fail("BigPsl chromSize disagrees with its target declaration")
      }
    None => bigpsl_fail("BigPsl record refers to an unknown target")
  }
  let block_count = record.block_sizes.length()
  if block_count == 0 ||
    record.block_starts.length() != block_count ||
    stored_query_starts.length() != block_count {
    bigpsl_fail("BigPsl block arrays disagree with blockCount")
  }
  let target_starts : Array[Int] = []
  let target_block_sizes : Array[Int] = []
  let query_block_sizes : Array[Int] = []
  let query_starts : Array[Int] = []
  for index in 0..<block_count {
    let target_block_size = record.block_sizes[index]
    let query_block_size = match sequence_type {
      BigPslProtein => {
        if target_block_size % 3 != 0 {
          bigpsl_fail(
            "Protein BigPsl target block size must be divisible by three",
          )
        }
        target_block_size / 3
      }
      _ => target_block_size
    }
    target_starts.push(record.chrom_start + record.block_starts[index])
    target_block_sizes.push(target_block_size)
    query_block_sizes.push(query_block_size)
    query_starts.push(stored_query_starts[index])
    if stored_query_starts[index] < 0 ||
      stored_query_starts[index] + query_block_size > query_size {
      bigpsl_fail("BigPsl query block extends beyond oChromSize")
    }
  }
  let (
    normalized_target_starts,
    normalized_query_starts,
    normalized_target_sizes,
    normalized_query_sizes,
  ) = if other_strand == "-" && record.strand == "-" {
    let transformed_targets : Array[Int] = []
    let transformed_queries : Array[Int] = []
    for index in 0..<block_count {
      transformed_targets.push(
        target_size - target_starts[index] - target_block_sizes[index],
      )
      transformed_queries.push(
        query_size - query_starts[index] - query_block_sizes[index],
      )
    }
    (
      bigpsl_reverse_ints(transformed_targets),
      bigpsl_reverse_ints(transformed_queries),
      bigpsl_reverse_ints(target_block_sizes),
      bigpsl_reverse_ints(query_block_sizes),
    )
  } else {
    (target_starts, query_starts, target_block_sizes, query_block_sizes)
  }
  let target_coordinates : Array[Int] = []
  let query_coordinates : Array[Int] = []
  let mut target_position = normalized_target_starts[0]
  let mut query_position = normalized_query_starts[0]
  target_coordinates.push(target_position)
  query_coordinates.push(query_position)
  for index in 0..<block_count {
    let target_start = normalized_target_starts[index]
    let query_start = normalized_query_starts[index]
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
    target_position = target_position + normalized_target_sizes[index]
    query_position = query_position + normalized_query_sizes[index]
    target_coordinates.push(target_position)
    query_coordinates.push(query_position)
  }
  if record.strand == "-" {
    if other_strand == "-" {
      for index in 0..<target_coordinates.length() {
        target_coordinates[index] = target_size - target_coordinates[index]
      }
    } else {
      for index in 0..<query_coordinates.length() {
        query_coordinates[index] = query_size - query_coordinates[index]
      }
    }
  }
  let (target_minimum, target_maximum) = bigpsl_coordinate_interval(
    target_coordinates,
  )
  let (query_minimum, query_maximum) = bigpsl_coordinate_interval(
    query_coordinates,
  )
  if target_minimum != record.chrom_start || target_maximum != record.chrom_end {
    bigpsl_fail("BigPsl BED interval disagrees with its alignment blocks")
  }
  if query_minimum != other_start || query_maximum != other_end {
    bigpsl_fail("BigPsl other interval disagrees with its alignment blocks")
  }
  BigPslAlignment::create(
    record.chrom,
    target_size,
    record.name,
    query_size,
    target_coordinates,
    query_coordinates,
    sequence_type~,
    query_sequence~,
    cds~,
    score=record.score,
    thick_start=record.thick_start,
    thick_end=record.thick_end,
    item_rgb=record.item_rgb,
    matches~,
    mismatches~,
    repeat_matches~,
    n_count~,
  )
}

///|
pub fn bigpsl_parse(data : Array[Int]) -> BigPslFile raise BigPslError {
  let bed = bigbed_parse(data) catch {
    BigBedError(message) => {
      bigpsl_fail("BigBed: " + message)
      BigBedFile::{
        header: BigBedHeader::{
          little_endian: true,
          version: 4,
          zoom_levels: 0,
          chromosome_tree_offset: 0,
          full_data_offset: 0,
          full_index_offset: 0,
          field_count: 0,
          defined_field_count: 0,
          auto_sql_offset: 0,
          total_summary_offset: 0,
          uncompress_buffer_size: 0,
          extra_indices_offset: 0,
        },
        schema: bigpsl_schema(),
        targets: [],
        source: [],
        blocks: [],
        record_count: 0,
      }
    }
  }
  bigpsl_validate_schema(bed)
  let records = bed.records() catch {
    BigBedError(message) => {
      bigpsl_fail("BigBed: " + message)
      []
    }
  }
  for record in records {
    ignore(bigpsl_decode_record(bed, record))
  }
  BigPslFile::{ bed, }
}

///|
pub fn BigPslFile::alignments(
  self : BigPslFile,
) -> Array[BigPslAlignment] raise BigPslError {
  let records = self.bed.records() catch {
    BigBedError(message) => {
      bigpsl_fail("BigBed: " + message)
      []
    }
  }
  let alignments : Array[BigPslAlignment] = []
  for record in records {
    alignments.push(bigpsl_decode_record(self.bed, record))
  }
  alignments
}

///|
pub fn BigPslFile::search(
  self : BigPslFile,
  target : String,
  start? : Int = 0,
  end? : Int = -1,
) -> Array[BigPslAlignment] raise BigPslError {
  let records = self.bed.search(target, start~, end~) catch {
    BigBedError(message) => {
      bigpsl_fail("BigBed: " + message)
      []
    }
  }
  let alignments : Array[BigPslAlignment] = []
  for record in records {
    alignments.push(bigpsl_decode_record(self.bed, record))
  }
  alignments
}

///|
pub fn BigPslFile::find_by_query(
  self : BigPslFile,
  query_name : String,
) -> Array[BigPslAlignment] raise BigPslError {
  let records = self.bed.find_by_name(query_name) catch {
    BigBedError(message) => {
      bigpsl_fail("BigBed: " + message)
      []
    }
  }
  let alignments : Array[BigPslAlignment] = []
  for record in records {
    alignments.push(bigpsl_decode_record(self.bed, record))
  }
  alignments
}

///|
pub fn BigPslFile::target(self : BigPslFile, name : String) -> BigBedTarget? {
  self.bed.target(name)
}

///|
pub fn BigPslFile::is_compressed(self : BigPslFile) -> Bool {
  self.bed.is_compressed()
}

///|
pub fn BigPslAlignment::to_psl(self : BigPslAlignment) -> String {
  let blocks = self.blocks()
  let counts = self.counts()
  let target_reverse = self.target_coordinates[0] >
    self.target_coordinates[self.target_coordinates.length() - 1]
  let query_reverse = self.query_coordinates[0] >
    self.query_coordinates[self.query_coordinates.length() - 1]
  let strand = if target_reverse {
    "+-"
  } else if query_reverse {
    "-"
  } else {
    "+"
  }
  let block_sizes : Array[Int] = []
  let query_starts : Array[Int] = []
  let target_starts : Array[Int] = []
  for block in blocks {
    block_sizes.push(block.query_size)
    let query_minimum = bigpsl_min(block.query_start, block.query_end)
    let query_maximum = bigpsl_max(block.query_start, block.query_end)
    let target_minimum = bigpsl_min(block.target_start, block.target_end)
    let target_maximum = bigpsl_max(block.target_start, block.target_end)
    query_starts.push(
      if query_reverse {
        self.query_size - query_maximum
      } else {
        query_minimum
      },
    )
    target_starts.push(
      if target_reverse {
        self.target_size - target_maximum
      } else {
        target_minimum
      },
    )
  }
  let (query_start, query_end) = self.query_interval()
  let (target_start, target_end) = self.target_interval()
  self.matches.to_string() +
  "\t" +
  self.mismatches.to_string() +
  "\t" +
  self.repeat_matches.to_string() +
  "\t" +
  self.n_count.to_string() +
  "\t" +
  counts.query_insert_events.to_string() +
  "\t" +
  counts.query_insert_units.to_string() +
  "\t" +
  counts.target_insert_events.to_string() +
  "\t" +
  counts.target_insert_bases.to_string() +
  "\t" +
  strand +
  "\t" +
  self.query_name +
  "\t" +
  self.query_size.to_string() +
  "\t" +
  query_start.to_string() +
  "\t" +
  query_end.to_string() +
  "\t" +
  self.target_name +
  "\t" +
  self.target_size.to_string() +
  "\t" +
  target_start.to_string() +
  "\t" +
  target_end.to_string() +
  "\t" +
  blocks.length().to_string() +
  "\t" +
  bigpsl_join_ints(block_sizes, true) +
  "\t" +
  bigpsl_join_ints(query_starts, true) +
  "\t" +
  bigpsl_join_ints(target_starts, true) +
  "\n"
}

///|
pub fn BigPslFile::to_psl(self : BigPslFile) -> String raise BigPslError {
  let output = StringBuilder()
  for alignment in self.alignments() {
    output.write_string(alignment.to_psl())
  }
  output.to_string()
}

///|
pub fn BigPslAlignment::summary(self : BigPslAlignment) -> String {
  let (start, end) = self.target_interval()
  "BigPslAlignment(" +
  self.target_name +
  ":" +
  start.to_string() +
  "-" +
  end.to_string() +
  " <- " +
  self.query_name +
  ", blocks=" +
  self.blocks().length().to_string() +
  ", type=" +
  (match self.sequence_type {
    BigPslUnknown => "unknown"
    BigPslNucleotide => "nucleotide"
    BigPslProtein => "protein"
  }) +
  ")"
}

///|
pub fn BigPslFile::summary(
  self : BigPslFile,
) -> BigPslSummary raise BigPslError {
  let alignments = self.alignments()
  let mut nucleotide_count = 0
  let mut protein_count = 0
  let mut unknown_count = 0
  let mut aligned_query_units = 0
  let mut covered_target_bases = 0
  for alignment in alignments {
    match alignment.sequence_type {
      BigPslNucleotide => nucleotide_count = nucleotide_count + 1
      BigPslProtein => protein_count = protein_count + 1
      BigPslUnknown => unknown_count = unknown_count + 1
    }
    aligned_query_units = aligned_query_units +
      alignment.counts().aligned_query_units
    let (start, end) = alignment.target_interval()
    covered_target_bases = covered_target_bases + end - start
  }
  BigPslSummary::{
    alignment_count: alignments.length(),
    target_count: self.bed.targets.length(),
    nucleotide_count,
    protein_count,
    unknown_count,
    aligned_query_units,
    covered_target_bases,
    compressed: self.is_compressed(),
  }
}

///|
pub fn bigpsl_example_data() -> (Array[BigBedTarget], Array[BigPslAlignment]) raise BigPslError {
  let targets = [
    BigBedTarget::create("chr1", 1000) catch {
      BigBedError(message) => {
        bigpsl_fail("BigBed: " + message)
        BigBedTarget::{ name: "chr1", length: 1000, }
      }
    },
    BigBedTarget::create("chr2", 800) catch {
      BigBedError(message) => {
        bigpsl_fail("BigBed: " + message)
        BigBedTarget::{ name: "chr2", length: 800, }
      }
    },
  ]
  let forward = BigPslAlignment::create(
    "chr1",
    1000,
    "rna-forward",
    30,
    [100, 110, 120, 130],
    [0, 10, 10, 20],
    sequence_type=BigPslNucleotide,
    score=900,
    thick_start=103,
    thick_end=127,
  )
  let reverse = BigPslAlignment::create(
    "chr1",
    1000,
    "rna-reverse",
    40,
    [200, 208, 208, 215],
    [30, 22, 20, 13],
    sequence_type=BigPslNucleotide,
    score=850,
  )
  let protein = BigPslAlignment::create(
    "chr2",
    800,
    "protein-reverse",
    20,
    [500, 488, 470, 470, 455],
    [0, 4, 4, 6, 11],
    sequence_type=BigPslProtein,
    score=1000,
  )
  (targets, [forward, reverse, protein])
}
