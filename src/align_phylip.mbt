// Biopython-compatible Bio.Align.phylip support.
//
// This module is separate from phylip_io.mbt. The older module exposes the
// historical AlignIO-style MultipleSeqAlignment API, while this file models
// the modern coordinate-bearing Alignment representation and strict PHYLIP
// sequential/interleaved parsing semantics.

///|
/// Error raised for malformed PHYLIP data or invalid alignment operations.
pub suberror AlignPhylipError {
  AlignPhylipError(String)
}

///|
/// Physical layout detected in or requested for a PHYLIP document.
pub(all) enum AlignPhylipLayout {
  PhylipSequential
  PhylipInterleaved
} derive(Eq, Debug)

///|
/// One coordinate-bearing PHYLIP row.
pub struct AlignPhylipSequence {
  id : String
  sequence : String
  aligned_sequence : String
} derive(Eq, Debug)

///|
/// A modern PHYLIP multiple sequence alignment.
///
/// PHYLIP identifiers read from files are the trimmed fixed-width 10-column
/// fields. Programmatically constructed rows may retain longer identifiers;
/// the writer applies Biopython's PHYLIP normalization when serializing them.
pub struct AlignPhylipAlignment {
  sequences : Array[AlignPhylipSequence]
  source_layout : AlignPhylipLayout
} derive(Eq, Debug)

///|
/// Pairwise or all-pairs alignment statistics.
pub struct AlignPhylipCounts {
  pairs : Int
  aligned : Int
  identities : Int
  mismatches : Int
  gap_columns : Int
  double_gap_columns : Int
  gap_opens : Int
} derive(Eq, Debug)

///|
/// Construct one validated PHYLIP row from its printed aligned sequence.
pub fn AlignPhylipSequence::create(
  id : String,
  aligned_sequence : String,
) -> AlignPhylipSequence raise AlignPhylipError {
  align_phylip_validate_id(id)
  if aligned_sequence.length() == 0 {
    raise AlignPhylipError("PHYLIP aligned sequence must not be empty")
  }
  align_phylip_validate_segment(aligned_sequence)
  AlignPhylipSequence::{
    id,
    sequence: align_phylip_remove_gaps(aligned_sequence),
    aligned_sequence,
  }
}

///|
/// Construct a validated modern PHYLIP alignment.
pub fn AlignPhylipAlignment::create(
  sequences : Array[AlignPhylipSequence],
  source_layout? : AlignPhylipLayout = PhylipSequential,
) -> AlignPhylipAlignment raise AlignPhylipError {
  let copied : Array[AlignPhylipSequence] = []
  for sequence in sequences {
    copied.push(sequence)
  }
  let alignment = AlignPhylipAlignment::{ sequences: copied, source_layout }
  align_phylip_validate_alignment(alignment)
  alignment
}

///|
/// Construct a coordinate-aware PHYLIP alignment from printed rows.
pub fn align_phylip_from_aligned(
  ids : Array[String],
  aligned_sequences : Array[String],
  source_layout? : AlignPhylipLayout = PhylipSequential,
) -> AlignPhylipAlignment raise AlignPhylipError {
  if ids.length() == 0 {
    raise AlignPhylipError(
      "PHYLIP alignment must contain at least one sequence",
    )
  }
  if ids.length() != aligned_sequences.length() {
    raise AlignPhylipError(
      "PHYLIP identifiers and aligned rows must have equal lengths",
    )
  }
  let sequences : Array[AlignPhylipSequence] = []
  for index = 0; index < ids.length(); index = index + 1 {
    sequences.push(
      AlignPhylipSequence::create(ids[index], aligned_sequences[index]),
    )
  }
  AlignPhylipAlignment::create(sequences, source_layout~)
}

///|
/// Parse one strict PHYLIP alignment.
///
/// The first line must contain the declared row and column counts. Identifiers
/// occupy exactly the first 10 columns of each first-block or first-row line.
/// As in Biopython 1.86, a blank line after the first `n` physical rows selects
/// interleaved parsing; otherwise wrapped sequential parsing is used. A file
/// containing exactly `n` complete rows is treated as sequential.
pub fn align_phylip_parse(
  text : String,
) -> AlignPhylipAlignment raise AlignPhylipError {
  let lines = align_phylip_lines(text)
  if lines.length() == 0 || lines[0].trim().length() == 0 {
    raise AlignPhylipError("Empty PHYLIP input")
  }
  let header = align_phylip_split_whitespace(lines[0])
  if header.length() != 2 {
    raise AlignPhylipError("PHYLIP header must contain exactly two integers")
  }
  let number_of_sequences = align_phylip_parse_positive_int(
    header[0],
    "PHYLIP sequence count",
  )
  let number_of_columns = align_phylip_parse_positive_int(
    header[1],
    "PHYLIP column count",
  )
  let mut last = lines.length() - 1
  while last > 0 && lines[last].trim().length() == 0 {
    last = last - 1
  }
  if last < number_of_sequences {
    raise AlignPhylipError(
      "PHYLIP input has fewer rows than declared in its header",
    )
  }
  for index = 1; index <= number_of_sequences; index = index + 1 {
    if lines[index].trim().length() == 0 {
      raise AlignPhylipError(
        "PHYLIP first block contains an empty sequence row",
      )
    }
  }
  let marker = number_of_sequences + 1
  let layout = if marker > last {
    PhylipSequential
  } else if lines[marker].trim().length() == 0 {
    PhylipInterleaved
  } else {
    PhylipSequential
  }
  let sequences = match layout {
    PhylipSequential =>
      align_phylip_parse_sequential(
        lines, last, number_of_sequences, number_of_columns,
      )
    PhylipInterleaved =>
      align_phylip_parse_interleaved(
        lines, last, number_of_sequences, number_of_columns,
      )
  }
  AlignPhylipAlignment::create(sequences, source_layout=layout)
}

///|
/// Serialize a PHYLIP alignment.
///
/// Defaults match Biopython's modern writer: sequential rows, no wrapping,
/// fixed 10-column normalized identifiers, and no grouping whitespace.
/// `block_width` and `group_width` enable canonical wrapped or interleaved
/// presentation without changing alignment coordinates.
pub fn align_phylip_write(
  alignment : AlignPhylipAlignment,
  layout? : AlignPhylipLayout = PhylipSequential,
  block_width? : Int = 0,
  group_width? : Int = 0,
) -> String raise AlignPhylipError {
  align_phylip_validate_alignment(alignment)
  if block_width < 0 {
    raise AlignPhylipError("PHYLIP block width must be non-negative")
  }
  if group_width < 0 {
    raise AlignPhylipError("PHYLIP group width must be non-negative")
  }
  let width = alignment.alignment_length()
  let effective_block_width = if block_width == 0 {
    if layout == PhylipInterleaved {
      60
    } else {
      width
    }
  } else {
    block_width
  }
  if effective_block_width <= 0 {
    raise AlignPhylipError("PHYLIP block width must be positive")
  }
  let names : Array[String] = []
  for sequence in alignment.sequences {
    names.push(align_phylip_normalize_id(sequence.id))
  }
  let output = StringBuilder::new()
  output.write_string(alignment.num_sequences().to_string())
  output.write_char(' ')
  output.write_string(width.to_string())
  output.write_char('\n')
  match layout {
    PhylipSequential =>
      for row = 0; row < alignment.num_sequences(); row = row + 1 {
        let aligned = alignment.sequences[row].aligned_sequence
        let mut start = 0
        while start < width {
          let stop = if start + effective_block_width < width {
            start + effective_block_width
          } else {
            width
          }
          if start == 0 {
            align_phylip_write_name(output, names[row])
          } else {
            output.write_string(" ".repeat(10))
          }
          align_phylip_write_grouped(
            output,
            aligned[start:stop].to_owned(),
            group_width,
          )
          output.write_char('\n')
          start = stop
        }
      }
    PhylipInterleaved => {
      let mut start = 0
      while start < width {
        let stop = if start + effective_block_width < width {
          start + effective_block_width
        } else {
          width
        }
        if start > 0 {
          output.write_char('\n')
        }
        for row = 0; row < alignment.num_sequences(); row = row + 1 {
          if start == 0 {
            align_phylip_write_name(output, names[row])
          } else {
            output.write_string(" ".repeat(10))
          }
          align_phylip_write_grouped(
            output,
            alignment.sequences[row].aligned_sequence[start:stop].to_owned(),
            group_width,
          )
          output.write_char('\n')
        }
        start = stop
      }
    }
  }
  let result = output.to_string()
  let reparsed = align_phylip_parse(result)
  if reparsed.num_sequences() != alignment.num_sequences() ||
    reparsed.alignment_length() != alignment.alignment_length() {
    raise AlignPhylipError("PHYLIP writer self-validation failed")
  }
  for row = 0; row < alignment.num_sequences(); row = row + 1 {
    if reparsed.sequences[row].id != names[row] ||
      reparsed.sequences[row].aligned_sequence !=
      alignment.sequences[row].aligned_sequence {
      raise AlignPhylipError("PHYLIP writer self-validation failed")
    }
  }
  result
}

///|
/// Apply Biopython's strict PHYLIP identifier normalization.
///
/// Outer whitespace is removed, `[](),` are deleted, `:;` become `|`, and
/// the result is truncated to the fixed 10-column identifier field.
pub fn align_phylip_normalize_id(id : String) -> String raise AlignPhylipError {
  align_phylip_validate_id(id)
  let trimmed = id.trim().to_owned()
  let output = StringBuilder::new(size_hint=trimmed.length())
  for index = 0; index < trimmed.length(); index = index + 1 {
    let code = trimmed.unsafe_get(index).to_int()
    if code == '['.to_int() ||
      code == ']'.to_int() ||
      code == '('.to_int() ||
      code == ')'.to_int() ||
      code == ','.to_int() {
      continue
    }
    if code == ':'.to_int() || code == ';'.to_int() {
      output.write_char('|')
    } else {
      output.write_char(code.unsafe_to_char())
    }
  }
  let normalized = output.to_string()
  if normalized.length() > 10 {
    normalized[0:10].to_owned()
  } else {
    normalized
  }
}

///|
/// Return serialized identifiers in row order.
pub fn AlignPhylipAlignment::normalized_ids(
  self : AlignPhylipAlignment,
) -> Array[String] {
  let result : Array[String] = []
  for sequence in self.sequences {
    result.push(
      align_phylip_normalize_id(sequence.id) catch {
        AlignPhylipError(_) => ""
      },
    )
  }
  result
}

///|
/// Return whether PHYLIP normalization maps two rows to the same identifier.
pub fn AlignPhylipAlignment::has_normalized_id_collisions(
  self : AlignPhylipAlignment,
) -> Bool {
  let names = self.normalized_ids()
  for first = 0; first < names.length(); first = first + 1 {
    for second = first + 1; second < names.length(); second = second + 1 {
      if names[first] == names[second] {
        return true
      }
    }
  }
  false
}

///|
/// Return the number of sequence rows.
pub fn AlignPhylipAlignment::num_sequences(self : AlignPhylipAlignment) -> Int {
  self.sequences.length()
}

///|
/// Return the alignment width including gaps.
pub fn AlignPhylipAlignment::alignment_length(
  self : AlignPhylipAlignment,
) -> Int {
  if self.sequences.length() == 0 {
    0
  } else {
    self.sequences[0].aligned_sequence.length()
  }
}

///|
/// Locate every row with an exact identifier. Duplicate strict identifiers
/// are legal in PHYLIP and therefore all matches are returned.
pub fn AlignPhylipAlignment::find_sequences(
  self : AlignPhylipAlignment,
  id : String,
) -> Array[Int] {
  let result : Array[Int] = []
  for index = 0; index < self.sequences.length(); index = index + 1 {
    if self.sequences[index].id == id {
      result.push(index)
    }
  }
  result
}

///|
/// Return one printed alignment column.
pub fn AlignPhylipAlignment::column(
  self : AlignPhylipAlignment,
  column : Int,
) -> String? {
  if column < 0 || column >= self.alignment_length() {
    return None
  }
  let result = StringBuilder::new(size_hint=self.sequences.length())
  for sequence in self.sequences {
    result.write_char(
      sequence.aligned_sequence.unsafe_get(column).unsafe_to_char(),
    )
  }
  Some(result.to_string())
}

///|
/// Map a zero-based ungapped row position to an alignment column.
pub fn AlignPhylipAlignment::sequence_position_to_column(
  self : AlignPhylipAlignment,
  row : Int,
  position : Int,
) -> Int? {
  if row < 0 ||
    row >= self.sequences.length() ||
    position < 0 ||
    position >= self.sequences[row].sequence.length() {
    return None
  }
  let aligned = self.sequences[row].aligned_sequence
  let mut coordinate = 0
  for column = 0; column < aligned.length(); column = column + 1 {
    if aligned.unsafe_get(column).to_int() != '-'.to_int() {
      if coordinate == position {
        return Some(column)
      }
      coordinate = coordinate + 1
    }
  }
  None
}

///|
/// Map one alignment column to a zero-based ungapped row position.
pub fn AlignPhylipAlignment::column_to_sequence_position(
  self : AlignPhylipAlignment,
  row : Int,
  column : Int,
) -> Int? {
  if row < 0 ||
    row >= self.sequences.length() ||
    column < 0 ||
    column >= self.alignment_length() {
    return None
  }
  let aligned = self.sequences[row].aligned_sequence
  if aligned.unsafe_get(column).to_int() == '-'.to_int() {
    return None
  }
  let mut coordinate = 0
  for index = 0; index < column; index = index + 1 {
    if aligned.unsafe_get(index).to_int() != '-'.to_int() {
      coordinate = coordinate + 1
    }
  }
  Some(coordinate)
}

///|
/// Project one row position through the alignment to another row.
pub fn AlignPhylipAlignment::map_position(
  self : AlignPhylipAlignment,
  source_row : Int,
  target_row : Int,
  source_position : Int,
) -> Int? {
  if target_row < 0 || target_row >= self.sequences.length() {
    return None
  }
  match self.sequence_position_to_column(source_row, source_position) {
    Some(column) => self.column_to_sequence_position(target_row, column)
    None => None
  }
}

///|
/// Return per-column residue positions for two rows.
pub fn AlignPhylipAlignment::aligned_pairs(
  self : AlignPhylipAlignment,
  first_row : Int,
  second_row : Int,
) -> Array[(Int?, Int?)] raise AlignPhylipError {
  align_phylip_validate_row(self, first_row)
  align_phylip_validate_row(self, second_row)
  let result : Array[(Int?, Int?)] = []
  let mut first_position = 0
  let mut second_position = 0
  for column = 0; column < self.alignment_length(); column = column + 1 {
    let first_gap = self.sequences[first_row].aligned_sequence
      .unsafe_get(column)
      .to_int() ==
      '-'.to_int()
    let second_gap = self.sequences[second_row].aligned_sequence
      .unsafe_get(column)
      .to_int() ==
      '-'.to_int()
    result.push(
      (
        if first_gap {
          None
        } else {
          Some(first_position)
        },
        if second_gap {
          None
        } else {
          Some(second_position)
        },
      ),
    )
    if !first_gap {
      first_position = first_position + 1
    }
    if !second_gap {
      second_position = second_position + 1
    }
  }
  result
}

///|
/// Return a compact Biopython-style coordinate path for all rows.
pub fn AlignPhylipAlignment::coordinate_path(
  self : AlignPhylipAlignment,
) -> Array[Array[Int]] {
  let paths : Array[Array[Int]] = []
  let coordinates = Array::make(self.sequences.length(), 0)
  for _row = 0; _row < self.sequences.length(); _row = _row + 1 {
    paths.push([0])
  }
  let width = self.alignment_length()
  for column = 0; column < width; column = column + 1 {
    for row = 0; row < self.sequences.length(); row = row + 1 {
      if self.sequences[row].aligned_sequence.unsafe_get(column).to_int() !=
        '-'.to_int() {
        coordinates[row] = coordinates[row] + 1
      }
    }
    let boundary = if column + 1 == width {
      true
    } else {
      align_phylip_movement_changes(self, column, column + 1)
    }
    if boundary {
      for row = 0; row < self.sequences.length(); row = row + 1 {
        paths[row].push(coordinates[row])
      }
    }
  }
  paths
}

///|
/// Compute statistics for one pair of rows.
pub fn AlignPhylipAlignment::pair_counts(
  self : AlignPhylipAlignment,
  first_row : Int,
  second_row : Int,
) -> AlignPhylipCounts raise AlignPhylipError {
  align_phylip_validate_row(self, first_row)
  align_phylip_validate_row(self, second_row)
  align_phylip_count_pair(self, first_row, second_row)
}

///|
/// Aggregate statistics across every unordered pair.
pub fn AlignPhylipAlignment::counts(
  self : AlignPhylipAlignment,
) -> AlignPhylipCounts {
  let mut pairs = 0
  let mut aligned = 0
  let mut identities = 0
  let mut mismatches = 0
  let mut gap_columns = 0
  let mut double_gap_columns = 0
  let mut gap_opens = 0
  for first = 0; first < self.sequences.length(); first = first + 1 {
    for second = first + 1
        second < self.sequences.length()
        second = second + 1 {
      let counts = align_phylip_count_pair(self, first, second)
      pairs = pairs + 1
      aligned = aligned + counts.aligned
      identities = identities + counts.identities
      mismatches = mismatches + counts.mismatches
      gap_columns = gap_columns + counts.gap_columns
      double_gap_columns = double_gap_columns + counts.double_gap_columns
      gap_opens = gap_opens + counts.gap_opens
    }
  }
  AlignPhylipCounts::{
    pairs,
    aligned,
    identities,
    mismatches,
    gap_columns,
    double_gap_columns,
    gap_opens,
  }
}

///|
/// Return exact identity over columns containing residues in both rows.
pub fn AlignPhylipCounts::identity(self : AlignPhylipCounts) -> Double {
  if self.aligned == 0 {
    0.0
  } else {
    self.identities.to_double() / self.aligned.to_double()
  }
}

///|
/// Return per-column non-gap occupancy.
pub fn AlignPhylipAlignment::occupancy(
  self : AlignPhylipAlignment,
) -> Array[Double] {
  let result : Array[Double] = []
  for column = 0; column < self.alignment_length(); column = column + 1 {
    let mut residues = 0
    for sequence in self.sequences {
      if sequence.aligned_sequence.unsafe_get(column).to_int() != '-'.to_int() {
        residues = residues + 1
      }
    }
    result.push(residues.to_double() / self.sequences.length().to_double())
  }
  result
}

///|
/// Calculate a majority residue consensus. Gaps do not vote.
pub fn AlignPhylipAlignment::majority_consensus(
  self : AlignPhylipAlignment,
  minimum_fraction? : Double = 0.0,
) -> String raise AlignPhylipError {
  if minimum_fraction != minimum_fraction ||
    minimum_fraction < 0.0 ||
    minimum_fraction > 1.0 {
    raise AlignPhylipError(
      "PHYLIP consensus minimum fraction must be between 0 and 1",
    )
  }
  let output = StringBuilder::new(size_hint=self.alignment_length())
  for column = 0; column < self.alignment_length(); column = column + 1 {
    let frequencies = Array::make(128, 0)
    let mut residues = 0
    for sequence in self.sequences {
      let code = sequence.aligned_sequence.unsafe_get(column).to_int()
      if code != '-'.to_int() {
        let upper = align_phylip_upper_code(code)
        if upper >= 0 && upper < frequencies.length() {
          frequencies[upper] = frequencies[upper] + 1
        }
        residues = residues + 1
      }
    }
    if residues == 0 {
      output.write_char('-')
      continue
    }
    let mut best_code = 0
    let mut best_count = -1
    for code = 0; code < frequencies.length(); code = code + 1 {
      if frequencies[code] > best_count {
        best_code = code
        best_count = frequencies[code]
      }
    }
    if best_count.to_double() / residues.to_double() < minimum_fraction {
      output.write_char('X')
    } else {
      output.write_char(best_code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
/// Return a compact alignment summary.
pub fn AlignPhylipAlignment::summary(self : AlignPhylipAlignment) -> String {
  "AlignPhylipAlignment(layout=" +
  (match self.source_layout {
    PhylipSequential => "sequential"
    PhylipInterleaved => "interleaved"
  }) +
  ", sequences=" +
  self.num_sequences().to_string() +
  ", columns=" +
  self.alignment_length().to_string() +
  ", normalized_collisions=" +
  self.has_normalized_id_collisions().to_string() +
  ")"
}

///|
/// Return a compact official-style interleaved PHYLIP fixture.
pub fn align_phylip_example_text() -> String {
  "3 12\n" +
  "Reference MKT--AIL\n" +
  "Query_one M-TQQAIL\n" +
  "Query_two MKT--AVL\n" +
  "\n" +
  "          GHAA\n" +
  "          G-AA\n" +
  "          GHAA\n"
}

///|
fn align_phylip_parse_sequential(
  lines : Array[String],
  last : Int,
  number_of_sequences : Int,
  number_of_columns : Int,
) -> Array[AlignPhylipSequence] raise AlignPhylipError {
  let result : Array[AlignPhylipSequence] = []
  let mut line_index = 1
  for _row = 0; _row < number_of_sequences; _row = _row + 1 {
    if line_index > last || lines[line_index].trim().length() == 0 {
      raise AlignPhylipError(
        "PHYLIP sequential data ended before all rows were read",
      )
    }
    let (id, first_segment) = align_phylip_parse_named_row(lines[line_index])
    line_index = line_index + 1
    let output = StringBuilder::new(size_hint=number_of_columns)
    output.write_string(first_segment)
    let mut length = first_segment.length()
    if length > number_of_columns {
      raise AlignPhylipError(
        "PHYLIP sequential row exceeds the declared column count",
      )
    }
    while length < number_of_columns {
      if line_index > last || lines[line_index].trim().length() == 0 {
        raise AlignPhylipError(
          "PHYLIP sequential row ended before its declared length",
        )
      }
      let segment = align_phylip_compact_segment(lines[line_index])
      align_phylip_validate_segment(segment)
      if length + segment.length() > number_of_columns {
        raise AlignPhylipError(
          "PHYLIP sequential row exceeds the declared column count",
        )
      }
      output.write_string(segment)
      length = length + segment.length()
      line_index = line_index + 1
    }
    result.push(AlignPhylipSequence::create(id, output.to_string()))
  }
  while line_index <= last {
    if lines[line_index].trim().length() != 0 {
      raise AlignPhylipError(
        "PHYLIP input contains trailing non-whitespace data",
      )
    }
    line_index = line_index + 1
  }
  result
}

///|
fn align_phylip_parse_interleaved(
  lines : Array[String],
  last : Int,
  number_of_sequences : Int,
  number_of_columns : Int,
) -> Array[AlignPhylipSequence] raise AlignPhylipError {
  let ids : Array[String] = []
  let builders : Array[StringBuilder] = []
  let lengths = Array::make(number_of_sequences, 0)
  let mut first_width = -1
  for row = 0; row < number_of_sequences; row = row + 1 {
    let (id, segment) = align_phylip_parse_named_row(lines[row + 1])
    if first_width < 0 {
      first_width = segment.length()
    } else if segment.length() != first_width {
      raise AlignPhylipError(
        "PHYLIP rows in an interleaved block must have equal widths",
      )
    }
    if segment.length() > number_of_columns {
      raise AlignPhylipError(
        "PHYLIP interleaved row exceeds the declared column count",
      )
    }
    ids.push(id)
    let output = StringBuilder::new(size_hint=number_of_columns)
    output.write_string(segment)
    builders.push(output)
    lengths[row] = segment.length()
  }
  let mut line_index = number_of_sequences + 1
  while lengths[0] < number_of_columns {
    while line_index <= last && lines[line_index].trim().length() == 0 {
      line_index = line_index + 1
    }
    if line_index > last {
      raise AlignPhylipError(
        "PHYLIP interleaved data ended before its declared length",
      )
    }
    let segments : Array[String] = []
    let mut block_width = -1
    for _row = 0; _row < number_of_sequences; _row = _row + 1 {
      if line_index > last || lines[line_index].trim().length() == 0 {
        raise AlignPhylipError(
          "PHYLIP interleaved block has fewer rows than declared",
        )
      }
      let segment = align_phylip_compact_segment(lines[line_index])
      align_phylip_validate_segment(segment)
      if block_width < 0 {
        block_width = segment.length()
      } else if segment.length() != block_width {
        raise AlignPhylipError(
          "PHYLIP rows in an interleaved block must have equal widths",
        )
      }
      segments.push(segment)
      line_index = line_index + 1
    }
    for row = 0; row < number_of_sequences; row = row + 1 {
      if lengths[row] + segments[row].length() > number_of_columns {
        raise AlignPhylipError(
          "PHYLIP interleaved row exceeds the declared column count",
        )
      }
      builders[row].write_string(segments[row])
      lengths[row] = lengths[row] + segments[row].length()
    }
  }
  for row = 0; row < number_of_sequences; row = row + 1 {
    if lengths[row] != number_of_columns {
      raise AlignPhylipError(
        "PHYLIP row length differs from the declared column count",
      )
    }
  }
  while line_index <= last {
    if lines[line_index].trim().length() != 0 {
      raise AlignPhylipError(
        "PHYLIP input contains trailing non-whitespace data",
      )
    }
    line_index = line_index + 1
  }
  let result : Array[AlignPhylipSequence] = []
  for row = 0; row < number_of_sequences; row = row + 1 {
    result.push(
      AlignPhylipSequence::create(ids[row], builders[row].to_string()),
    )
  }
  result
}

///|
fn align_phylip_parse_named_row(
  line : String,
) -> (String, String) raise AlignPhylipError {
  if line.length() < 10 {
    raise AlignPhylipError(
      "PHYLIP sequence row is shorter than the 10-column identifier field",
    )
  }
  let id = line[0:10].trim().to_owned()
  align_phylip_validate_id(id)
  let segment = align_phylip_compact_segment(line[10:line.length()].to_owned())
  align_phylip_validate_segment(segment)
  (id, segment)
}

///|
fn align_phylip_validate_alignment(
  alignment : AlignPhylipAlignment,
) -> Unit raise AlignPhylipError {
  if alignment.sequences.length() == 0 {
    raise AlignPhylipError(
      "PHYLIP alignment must contain at least one sequence",
    )
  }
  let width = alignment.sequences[0].aligned_sequence.length()
  if width == 0 {
    raise AlignPhylipError("PHYLIP alignment must contain at least one column")
  }
  for sequence in alignment.sequences {
    align_phylip_validate_id(sequence.id)
    if sequence.aligned_sequence.length() != width {
      raise AlignPhylipError("PHYLIP aligned rows must have equal widths")
    }
    align_phylip_validate_segment(sequence.aligned_sequence)
    if align_phylip_remove_gaps(sequence.aligned_sequence) != sequence.sequence {
      raise AlignPhylipError(
        "PHYLIP sequence row has inconsistent derived fields",
      )
    }
  }
  for column = 0; column < width; column = column + 1 {
    let mut residues = 0
    for sequence in alignment.sequences {
      if sequence.aligned_sequence.unsafe_get(column).to_int() != '-'.to_int() {
        residues = residues + 1
      }
    }
    if residues == 0 {
      raise AlignPhylipError(
        "PHYLIP alignment contains an all-gap column at " + column.to_string(),
      )
    }
  }
}

///|
fn align_phylip_count_pair(
  alignment : AlignPhylipAlignment,
  first_row : Int,
  second_row : Int,
) -> AlignPhylipCounts {
  let first = alignment.sequences[first_row].aligned_sequence
  let second = alignment.sequences[second_row].aligned_sequence
  let mut aligned = 0
  let mut identities = 0
  let mut mismatches = 0
  let mut gap_columns = 0
  let mut double_gap_columns = 0
  let mut gap_opens = 0
  let mut in_gap = false
  for column = 0; column < alignment.alignment_length(); column = column + 1 {
    let first_code = first.unsafe_get(column).to_int()
    let second_code = second.unsafe_get(column).to_int()
    let first_gap = first_code == '-'.to_int()
    let second_gap = second_code == '-'.to_int()
    if first_gap && second_gap {
      double_gap_columns = double_gap_columns + 1
      in_gap = false
    } else if first_gap || second_gap {
      gap_columns = gap_columns + 1
      if !in_gap {
        gap_opens = gap_opens + 1
      }
      in_gap = true
    } else {
      aligned = aligned + 1
      if first_code == second_code {
        identities = identities + 1
      } else {
        mismatches = mismatches + 1
      }
      in_gap = false
    }
  }
  AlignPhylipCounts::{
    pairs: 1,
    aligned,
    identities,
    mismatches,
    gap_columns,
    double_gap_columns,
    gap_opens,
  }
}

///|
fn align_phylip_movement_changes(
  alignment : AlignPhylipAlignment,
  first_column : Int,
  second_column : Int,
) -> Bool {
  for sequence in alignment.sequences {
    let first_moves = sequence.aligned_sequence
      .unsafe_get(first_column)
      .to_int() !=
      '-'.to_int()
    let second_moves = sequence.aligned_sequence
      .unsafe_get(second_column)
      .to_int() !=
      '-'.to_int()
    if first_moves != second_moves {
      return true
    }
  }
  false
}

///|
fn align_phylip_write_name(output : StringBuilder, name : String) -> Unit {
  output.write_string(name)
  if name.length() < 10 {
    output.write_string(" ".repeat(10 - name.length()))
  }
}

///|
fn align_phylip_write_grouped(
  output : StringBuilder,
  segment : String,
  group_width : Int,
) -> Unit {
  if group_width == 0 || group_width >= segment.length() {
    output.write_string(segment)
    return
  }
  let mut start = 0
  while start < segment.length() {
    let stop = if start + group_width < segment.length() {
      start + group_width
    } else {
      segment.length()
    }
    if start > 0 {
      output.write_char(' ')
    }
    output.write_string(segment[start:stop].to_owned())
    start = stop
  }
}

///|
fn align_phylip_validate_id(id : String) -> Unit raise AlignPhylipError {
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    if code < 32 || code > 126 || code == '\n'.to_int() || code == '\r'.to_int() {
      raise AlignPhylipError(
        "PHYLIP sequence identifier contains a control or non-ASCII character",
      )
    }
  }
}

///|
fn align_phylip_validate_segment(
  segment : String,
) -> Unit raise AlignPhylipError {
  if segment.length() == 0 {
    raise AlignPhylipError("PHYLIP sequence segment must not be empty")
  }
  for index = 0; index < segment.length(); index = index + 1 {
    let code = segment.unsafe_get(index).to_int()
    if code == '.'.to_int() {
      raise AlignPhylipError("PHYLIP format no longer allows dots in sequences")
    }
    if code <= 32 || code > 126 {
      raise AlignPhylipError(
        "PHYLIP sequence contains whitespace or a non-ASCII character",
      )
    }
  }
}

///|
fn align_phylip_validate_row(
  alignment : AlignPhylipAlignment,
  row : Int,
) -> Unit raise AlignPhylipError {
  if row < 0 || row >= alignment.sequences.length() {
    raise AlignPhylipError("PHYLIP row index is out of bounds")
  }
}

///|
fn align_phylip_compact_segment(value : String) -> String {
  let output = StringBuilder::new(size_hint=value.length())
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    if code != ' '.to_int() && code != '\t'.to_int() {
      output.write_char(code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn align_phylip_remove_gaps(value : String) -> String {
  let output = StringBuilder::new(size_hint=value.length())
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    if code != '-'.to_int() {
      output.write_char(code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn align_phylip_upper_code(code : Int) -> Int {
  if code >= 'a'.to_int() && code <= 'z'.to_int() {
    code - 'a'.to_int() + 'A'.to_int()
  } else {
    code
  }
}

///|
fn align_phylip_lines(text : String) -> Array[String] {
  let raw = text.split("\n").to_array()
  let lines : Array[String] = []
  for line in raw {
    if line.length() > 0 &&
      line.unsafe_get(line.length() - 1).to_int() == '\r'.to_int() {
      lines.push(line[0:line.length() - 1].to_owned())
    } else {
      lines.push(line.to_owned())
    }
  }
  lines
}

///|
fn align_phylip_split_whitespace(value : String) -> Array[String] {
  let fields : Array[String] = []
  let mut start = 0
  let mut in_field = false
  for index = 0; index < value.length(); index = index + 1 {
    let whitespace = align_phylip_is_whitespace(
      value.unsafe_get(index).to_int(),
    )
    if whitespace {
      if in_field {
        fields.push(value[start:index].to_owned())
        in_field = false
      }
    } else if !in_field {
      start = index
      in_field = true
    }
  }
  if in_field {
    fields.push(value[start:value.length()].to_owned())
  }
  fields
}

///|
fn align_phylip_parse_positive_int(
  text : String,
  label : String,
) -> Int raise AlignPhylipError {
  if text.length() == 0 {
    raise AlignPhylipError(label + " is empty")
  }
  let mut value = 0
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise AlignPhylipError(label + " must be a positive integer")
    }
    let digit = code - '0'.to_int()
    if value > 214748364 || (value == 214748364 && digit > 7) {
      raise AlignPhylipError(label + " is outside the supported integer range")
    }
    value = value * 10 + digit
  }
  if value <= 0 {
    raise AlignPhylipError(label + " must be positive")
  }
  value
}

///|
fn align_phylip_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\r'.to_int() ||
  code == '\n'.to_int()
}
