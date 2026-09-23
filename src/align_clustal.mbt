// Biopython-compatible Bio.Align.clustal support.
//
// This module is separate from clustal_io.mbt. The older module exposes the
// historical AlignIO-style MultipleSeqAlignment API, while this file retains
// modern file metadata, column annotations, and coordinate-aware rows.

///|
/// Error raised for malformed CLUSTAL data or invalid alignment operations.
pub suberror AlignClustalError {
  AlignClustalError(String)
}

///|
/// File-level metadata from a CLUSTAL alignment.
///
/// Known generators follow Biopython 1.86: CLUSTAL, PROBCONS, MUSCLE,
/// MSAPROBS, Kalign, and Biopython.
pub struct AlignClustalMetadata {
  program : String
  version : String
  header : String
} derive(Eq, Debug)
pub extend AlignClustalMetadata with Eq::{equal, not_equal}
pub extend AlignClustalMetadata with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One coordinate-bearing row in a CLUSTAL alignment.
pub struct AlignClustalSequence {
  id : String
  sequence : String
  aligned_sequence : String
} derive(Eq, Debug)
pub extend AlignClustalSequence with Eq::{equal, not_equal}
pub extend AlignClustalSequence with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// A modern CLUSTAL multiple sequence alignment.
///
/// `consensus` retains the `*`, `:`, `.`, and space column annotation when it
/// is present. Sparse source annotations are padded with spaces.
pub struct AlignClustalAlignment {
  metadata : AlignClustalMetadata
  sequences : Array[AlignClustalSequence]
  consensus : String?
} derive(Eq, Debug)
pub extend AlignClustalAlignment with Eq::{equal, not_equal}
pub extend AlignClustalAlignment with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Pairwise or all-pairs alignment statistics.
pub struct AlignClustalCounts {
  pairs : Int
  aligned : Int
  identities : Int
  mismatches : Int
  gap_columns : Int
  double_gap_columns : Int
  gap_opens : Int
} derive(Eq, Debug)
pub extend AlignClustalCounts with Eq::{equal, not_equal}
pub extend AlignClustalCounts with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One parsed interleaved block.
priv struct AlignClustalBlock {
  ids : Array[String]
  segments : Array[String]
  residue_counts : Array[Int?]
  consensus : String?
  width : Int
}

///|
/// Construct validated CLUSTAL file metadata.
pub fn AlignClustalMetadata::create(
  program : String,
  version? : String = "",
  header? : String = "",
) -> AlignClustalMetadata raise AlignClustalError {
  align_clustal_validate_program(program)
  align_clustal_validate_version(version)
  let canonical_header = if header.length() == 0 {
    align_clustal_canonical_header(program, version)
  } else {
    align_clustal_validate_plain_text(header, "CLUSTAL header")
    header
  }
  AlignClustalMetadata::{ program, version, header: canonical_header, }
}

///|
/// Construct one validated CLUSTAL row.
pub fn AlignClustalSequence::create(
  id : String,
  aligned_sequence : String,
) -> AlignClustalSequence raise AlignClustalError {
  align_clustal_validate_id(id)
  if aligned_sequence.length() == 0 {
    raise AlignClustalError("CLUSTAL aligned sequence must not be empty")
  }
  align_clustal_validate_segment(aligned_sequence)
  let sequence = align_clustal_remove_gaps(aligned_sequence)
  if sequence.length() == 0 {
    raise AlignClustalError(
      "CLUSTAL sequence row must contain at least one residue",
    )
  }
  AlignClustalSequence::{ id, sequence, aligned_sequence, }
}

///|
/// Construct a validated modern CLUSTAL alignment.
pub fn AlignClustalAlignment::create(
  metadata : AlignClustalMetadata,
  sequences : Array[AlignClustalSequence],
  consensus? : String? = None,
) -> AlignClustalAlignment raise AlignClustalError {
  let copied : Array[AlignClustalSequence] = []
  for sequence in sequences {
    copied.push(sequence)
  }
  let alignment = AlignClustalAlignment::{
    metadata,
    sequences: copied,
    consensus,
  }
  align_clustal_validate_alignment(alignment)
  alignment
}

///|
/// Construct a coordinate-aware alignment from printed rows.
pub fn align_clustal_from_aligned(
  ids : Array[String],
  aligned_sequences : Array[String],
  program? : String = "Biopython",
  version? : String = "",
  consensus? : String? = None,
) -> AlignClustalAlignment raise AlignClustalError {
  if ids.length() == 0 {
    raise AlignClustalError(
      "CLUSTAL alignment must contain at least one sequence",
    )
  }
  if ids.length() != aligned_sequences.length() {
    raise AlignClustalError(
      "CLUSTAL identifiers and aligned rows must have equal lengths",
    )
  }
  let sequences : Array[AlignClustalSequence] = []
  for index = 0; index < ids.length(); index = index + 1 {
    sequences.push(
      AlignClustalSequence::create(ids[index], aligned_sequences[index]),
    )
  }
  let metadata = AlignClustalMetadata::create(program, version~)
  AlignClustalAlignment::create(metadata, sequences, consensus~)
}

///|
/// Parse one CLUSTAL alignment.
///
/// Blocks must repeat the first block's unique identifiers in the same order.
/// Optional right-hand residue counts are validated as cumulative ungapped
/// counts. Consensus lines may be absent from some blocks; missing portions
/// are represented by spaces in the retained column annotation.
pub fn align_clustal_parse(
  text : String,
) -> AlignClustalAlignment raise AlignClustalError {
  let lines = align_clustal_lines(text)
  if lines.length() == 0 || lines[0].trim().length() == 0 {
    raise AlignClustalError("Empty CLUSTAL input")
  }
  let metadata = align_clustal_parse_header(lines[0])
  let blocks : Array[AlignClustalBlock] = []
  let mut block_lines : Array[String] = []
  for line_index = 1; line_index < lines.length(); line_index = line_index + 1 {
    let line = lines[line_index]
    if line.trim().length() == 0 {
      if line.length() >= 10 && block_lines.length() > 0 {
        block_lines.push(line)
      } else if block_lines.length() > 0 {
        blocks.push(align_clustal_parse_block(block_lines))
        block_lines = []
      }
    } else {
      block_lines.push(line)
    }
  }
  if block_lines.length() > 0 {
    blocks.push(align_clustal_parse_block(block_lines))
  }
  if blocks.length() == 0 {
    raise AlignClustalError("CLUSTAL input contains no alignment rows")
  }

  let ids = blocks[0].ids
  let builders : Array[StringBuilder] = []
  let cumulative = Array::make(ids.length(), 0)
  let seen_ids : Array[String] = []
  for id in ids {
    align_clustal_validate_id(id)
    if align_clustal_find_string(seen_ids, id) >= 0 {
      raise AlignClustalError(
        "Duplicated CLUSTAL sequence identifier '" + id + "'",
      )
    }
    seen_ids.push(id)
    builders.push(StringBuilder::new())
  }

  let consensus_builder = StringBuilder::new()
  let mut has_consensus = false
  for block_index = 0
      block_index < blocks.length()
      block_index = block_index + 1 {
    let block = blocks[block_index]
    if block.ids.length() != ids.length() {
      raise AlignClustalError(
        "CLUSTAL block sequence count differs from the first block",
      )
    }
    for row = 0; row < ids.length(); row = row + 1 {
      if block.ids[row] != ids[row] {
        raise AlignClustalError(
          "CLUSTAL block identifiers must repeat in the same order",
        )
      }
      builders[row].write_string(block.segments[row])
      cumulative[row] = cumulative[row] +
        align_clustal_residue_count(block.segments[row])
      match block.residue_counts[row] {
        Some(declared) =>
          if declared != cumulative[row] {
            raise AlignClustalError(
              "CLUSTAL cumulative residue count mismatch for '" + ids[row] + "'",
            )
          }
        None => ()
      }
    }
    match block.consensus {
      Some(annotation) => {
        has_consensus = true
        consensus_builder.write_string(annotation)
      }
      None => consensus_builder.write_string(" ".repeat(block.width))
    }
  }

  let sequences : Array[AlignClustalSequence] = []
  for row = 0; row < ids.length(); row = row + 1 {
    sequences.push(
      AlignClustalSequence::create(ids[row], builders[row].to_string()),
    )
  }
  let consensus : String? = if has_consensus {
    Some(consensus_builder.to_string())
  } else {
    None
  }
  AlignClustalAlignment::create(metadata, sequences, consensus~)
}

///|
/// Write canonical CLUSTAL text.
///
/// Biopython-compatible defaults use 50 alignment columns and a 36-character
/// name field. Identifiers are never silently truncated; names longer than 30
/// characters require a larger `name_width`.
pub fn align_clustal_write(
  alignment : AlignClustalAlignment,
  block_width? : Int = 50,
  name_width? : Int = 36,
  include_counts? : Bool = false,
  include_consensus? : Bool = true,
) -> String raise AlignClustalError {
  align_clustal_validate_alignment(alignment)
  if block_width <= 0 {
    raise AlignClustalError("CLUSTAL block width must be positive")
  }
  if name_width <= 1 {
    raise AlignClustalError("CLUSTAL name width must be greater than one")
  }
  for sequence in alignment.sequences {
    if sequence.id.length() >= name_width {
      raise AlignClustalError(
        "CLUSTAL name width must exceed every identifier length",
      )
    }
    if name_width == 36 && sequence.id.length() > 30 {
      raise AlignClustalError(
        "CLUSTAL canonical writer does not silently truncate identifiers " +
        "longer than 30 characters",
      )
    }
  }
  let output = StringBuilder::new()
  output.write_string(
    align_clustal_canonical_header(
      alignment.metadata.program,
      alignment.metadata.version,
    ),
  )
  output.write_string("\n\n\n")
  let cumulative = Array::make(alignment.num_sequences(), 0)
  let width = alignment.alignment_length()
  let mut start = 0
  while start < width {
    let stop = if start + block_width < width {
      start + block_width
    } else {
      width
    }
    for row = 0; row < alignment.num_sequences(); row = row + 1 {
      let sequence = alignment.sequences[row]
      let segment = sequence.aligned_sequence[start:stop].to_owned()
      output.write_string(sequence.id)
      output.write_string(" ".repeat(name_width - sequence.id.length()))
      output.write_string(segment)
      cumulative[row] = cumulative[row] + align_clustal_residue_count(segment)
      if include_counts {
        output.write_char(' ')
        output.write_string(cumulative[row].to_string())
      }
      output.write_char('\n')
    }
    if include_consensus {
      match alignment.consensus {
        Some(consensus) => {
          output.write_string(" ".repeat(name_width))
          output.write_string(consensus[start:stop].to_owned())
          output.write_char('\n')
        }
        None => ()
      }
    }
    output.write_char('\n')
    start = stop
  }
  output.write_char('\n')
  let result = output.to_string()
  ignore(align_clustal_parse(result))
  result
}

///|
/// Return the number of sequence rows.
pub fn AlignClustalAlignment::num_sequences(
  self : AlignClustalAlignment,
) -> Int {
  self.sequences.length()
}

///|
/// Return the alignment width including gaps.
pub fn AlignClustalAlignment::alignment_length(
  self : AlignClustalAlignment,
) -> Int {
  if self.sequences.length() == 0 {
    0
  } else {
    self.sequences[0].aligned_sequence.length()
  }
}

///|
/// Locate a row by exact identifier.
pub fn AlignClustalAlignment::find_sequence(
  self : AlignClustalAlignment,
  id : String,
) -> Int? {
  for index = 0; index < self.sequences.length(); index = index + 1 {
    if self.sequences[index].id == id {
      return Some(index)
    }
  }
  None
}

///|
/// Return one printed alignment column.
pub fn AlignClustalAlignment::column(
  self : AlignClustalAlignment,
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
pub fn AlignClustalAlignment::sequence_position_to_column(
  self : AlignClustalAlignment,
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
pub fn AlignClustalAlignment::column_to_sequence_position(
  self : AlignClustalAlignment,
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
pub fn AlignClustalAlignment::map_position(
  self : AlignClustalAlignment,
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
pub fn AlignClustalAlignment::aligned_pairs(
  self : AlignClustalAlignment,
  first_row : Int,
  second_row : Int,
) -> Array[(Int?, Int?)] raise AlignClustalError {
  align_clustal_validate_row(self, first_row)
  align_clustal_validate_row(self, second_row)
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
pub fn AlignClustalAlignment::coordinate_path(
  self : AlignClustalAlignment,
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
      align_clustal_movement_changes(self, column, column + 1)
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
pub fn AlignClustalAlignment::pair_counts(
  self : AlignClustalAlignment,
  first_row : Int,
  second_row : Int,
) -> AlignClustalCounts raise AlignClustalError {
  align_clustal_validate_row(self, first_row)
  align_clustal_validate_row(self, second_row)
  align_clustal_count_pair(self, first_row, second_row)
}

///|
/// Aggregate statistics across every unordered pair.
pub fn AlignClustalAlignment::counts(
  self : AlignClustalAlignment,
) -> AlignClustalCounts {
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
      let counts = align_clustal_count_pair(self, first, second)
      pairs = pairs + 1
      aligned = aligned + counts.aligned
      identities = identities + counts.identities
      mismatches = mismatches + counts.mismatches
      gap_columns = gap_columns + counts.gap_columns
      double_gap_columns = double_gap_columns + counts.double_gap_columns
      gap_opens = gap_opens + counts.gap_opens
    }
  }
  AlignClustalCounts::{
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
pub fn AlignClustalCounts::identity(self : AlignClustalCounts) -> Double {
  if self.aligned == 0 {
    0.0
  } else {
    self.identities.to_double() / self.aligned.to_double()
  }
}

///|
/// Return per-column non-gap occupancy.
pub fn AlignClustalAlignment::occupancy(
  self : AlignClustalAlignment,
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
pub fn AlignClustalAlignment::majority_consensus(
  self : AlignClustalAlignment,
  minimum_fraction? : Double = 0.0,
) -> String raise AlignClustalError {
  if minimum_fraction != minimum_fraction ||
    minimum_fraction < 0.0 ||
    minimum_fraction > 1.0 {
    raise AlignClustalError(
      "CLUSTAL consensus minimum fraction must be between 0 and 1",
    )
  }
  let output = StringBuilder::new(size_hint=self.alignment_length())
  for column = 0; column < self.alignment_length(); column = column + 1 {
    let frequencies = Array::make(128, 0)
    let mut residues = 0
    for sequence in self.sequences {
      let code = sequence.aligned_sequence.unsafe_get(column).to_int()
      if code != '-'.to_int() {
        let upper = align_clustal_upper_code(code)
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
/// Calculate CLUSTAL `*`, `:`, `.`, and space conservation symbols.
///
/// Nucleotide alignments receive `*` only for exact fully occupied columns.
/// Protein alignments additionally use the standard strong and weak ClustalW
/// residue groups.
pub fn AlignClustalAlignment::calculated_consensus(
  self : AlignClustalAlignment,
) -> String {
  let output = StringBuilder::new(size_hint=self.alignment_length())
  let nucleotide = align_clustal_is_nucleotide(self)
  let strong_groups = [
    "STA", "NEQK", "NHQK", "NDEQ", "QHRK", "MILV", "MILF", "HY", "FYW",
  ]
  let weak_groups = [
    "CSA", "ATV", "SAG", "STNK", "STPA", "SGND", "SNDEQK", "NDEQHK", "NEQHRK", "FVLIM",
    "HFY",
  ]
  for column = 0; column < self.alignment_length(); column = column + 1 {
    if align_clustal_column_has_gap(self, column) {
      output.write_char(' ')
    } else if align_clustal_column_identical(self, column) {
      output.write_char('*')
    } else if !nucleotide &&
      align_clustal_column_in_any_group(self, column, strong_groups) {
      output.write_char(':')
    } else if !nucleotide &&
      align_clustal_column_in_any_group(self, column, weak_groups) {
      output.write_char('.')
    } else {
      output.write_char(' ')
    }
  }
  output.to_string()
}

///|
/// Return a compact alignment summary.
pub fn AlignClustalAlignment::summary(self : AlignClustalAlignment) -> String {
  "AlignClustalAlignment(program=" +
  self.metadata.program +
  ", version=" +
  self.metadata.version +
  ", sequences=" +
  self.num_sequences().to_string() +
  ", columns=" +
  self.alignment_length().to_string() +
  ", consensus=" +
  (match self.consensus {
    Some(_) => "yes"
    None => "no"
  }) +
  ")"
}

///|
/// Return a compact interleaved CLUSTAL example with counts and consensus.
pub fn align_clustal_example_text() -> String {
  "CLUSTAL W (2.1) multiple sequence alignment\n" +
  "\n" +
  "\n" +
  "reference      MKT--AI 5\n" +
  "query_one     M-TQQAI 6\n" +
  "query_two     MKT--AV 5\n" +
  "               * *  *:\n" +
  "\n" +
  "reference      LGH 8\n" +
  "query_one      LG- 8\n" +
  "query_two      LGH 8\n" +
  "               ** \n"
}

///|
fn align_clustal_validate_alignment(
  alignment : AlignClustalAlignment,
) -> Unit raise AlignClustalError {
  align_clustal_validate_program(alignment.metadata.program)
  align_clustal_validate_version(alignment.metadata.version)
  align_clustal_validate_plain_text(alignment.metadata.header, "CLUSTAL header")
  if alignment.sequences.length() == 0 {
    raise AlignClustalError(
      "CLUSTAL alignment must contain at least one sequence",
    )
  }
  let width = alignment.sequences[0].aligned_sequence.length()
  if width == 0 {
    raise AlignClustalError(
      "CLUSTAL alignment must contain at least one column",
    )
  }
  let ids : Array[String] = []
  for sequence in alignment.sequences {
    align_clustal_validate_id(sequence.id)
    if align_clustal_find_string(ids, sequence.id) >= 0 {
      raise AlignClustalError(
        "Duplicated CLUSTAL sequence identifier '" + sequence.id + "'",
      )
    }
    ids.push(sequence.id)
    if sequence.aligned_sequence.length() != width {
      raise AlignClustalError("CLUSTAL aligned rows must have equal widths")
    }
    align_clustal_validate_segment(sequence.aligned_sequence)
    if align_clustal_remove_gaps(sequence.aligned_sequence) != sequence.sequence {
      raise AlignClustalError(
        "CLUSTAL sequence row has inconsistent derived fields",
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
      raise AlignClustalError(
        "CLUSTAL alignment contains an all-gap column at " + column.to_string(),
      )
    }
  }
  match alignment.consensus {
    Some(consensus) => {
      if consensus.length() != width {
        raise AlignClustalError(
          "CLUSTAL consensus length must equal the alignment width",
        )
      }
      align_clustal_validate_consensus(consensus)
    }
    None => ()
  }
}

///|
fn align_clustal_parse_header(
  line : String,
) -> AlignClustalMetadata raise AlignClustalError {
  let words = align_clustal_split_whitespace(line)
  if words.length() == 0 {
    raise AlignClustalError("Empty CLUSTAL header")
  }
  let program = words[0]
  align_clustal_validate_program(program)
  let mut version = ""
  for index = 1; index < words.length(); index = index + 1 {
    let word = align_clustal_strip_parentheses(words[index])
    if word.length() > 0 {
      let first = word.unsafe_get(0).to_int()
      if first >= '0'.to_int() && first <= '9'.to_int() {
        version = word
        break
      }
    }
  }
  AlignClustalMetadata::create(program, version~, header=line)
}

///|
fn align_clustal_parse_block(
  lines : Array[String],
) -> AlignClustalBlock raise AlignClustalError {
  let ids : Array[String] = []
  let segments : Array[String] = []
  let residue_counts : Array[Int?] = []
  let mut segment_start = -1
  let mut width = -1
  let mut consensus : String? = None
  let mut saw_consensus = false
  for line in lines {
    if line.length() == 0 {
      continue
    }
    let first = line.unsafe_get(0).to_int()
    if align_clustal_is_whitespace(first) {
      if ids.length() == 0 {
        raise AlignClustalError(
          "CLUSTAL consensus line appears before sequence rows",
        )
      }
      if saw_consensus {
        raise AlignClustalError(
          "CLUSTAL block contains multiple consensus lines",
        )
      }
      saw_consensus = true
      let annotation = align_clustal_extract_consensus(
        line, segment_start, width,
      )
      consensus = Some(annotation)
    } else {
      if saw_consensus {
        raise AlignClustalError(
          "CLUSTAL sequence row appears after a consensus line",
        )
      }
      let fields = align_clustal_split_whitespace(line)
      if fields.length() < 2 || fields.length() > 3 {
        raise AlignClustalError(
          "CLUSTAL sequence line must contain two or three fields",
        )
      }
      let id = fields[0]
      let segment = fields[1]
      align_clustal_validate_id(id)
      align_clustal_validate_segment(segment)
      if segment_start < 0 {
        segment_start = align_clustal_segment_start(line)
      }
      if width < 0 {
        width = segment.length()
      } else if segment.length() != width {
        raise AlignClustalError(
          "CLUSTAL rows in one block must have equal widths",
        )
      }
      ids.push(id)
      segments.push(segment)
      if fields.length() == 3 {
        residue_counts.push(
          Some(
            align_clustal_parse_nonnegative_int(
              fields[2],
              "CLUSTAL residue count",
            ),
          ),
        )
      } else {
        residue_counts.push(None)
      }
    }
  }
  if ids.length() == 0 || width <= 0 {
    raise AlignClustalError("CLUSTAL block contains no sequence rows")
  }
  AlignClustalBlock::{ ids, segments, residue_counts, consensus, width, }
}

///|
fn align_clustal_extract_consensus(
  line : String,
  start : Int,
  width : Int,
) -> String raise AlignClustalError {
  if start < 0 || width <= 0 {
    raise AlignClustalError("Invalid CLUSTAL consensus position")
  }
  for index = 0; index < start && index < line.length(); index = index + 1 {
    if !align_clustal_is_whitespace(line.unsafe_get(index).to_int()) {
      raise AlignClustalError(
        "CLUSTAL consensus line is not aligned with sequence data",
      )
    }
  }
  let output = StringBuilder::new(size_hint=width)
  for offset = 0; offset < width; offset = offset + 1 {
    let index = start + offset
    if index >= line.length() {
      output.write_char(' ')
    } else {
      let code = line.unsafe_get(index).to_int()
      if code == ' '.to_int() {
        output.write_char(' ')
      } else if code == '*'.to_int() ||
        code == ':'.to_int() ||
        code == '.'.to_int() {
        output.write_char(code.unsafe_to_char())
      } else {
        raise AlignClustalError("CLUSTAL consensus contains an invalid symbol")
      }
    }
  }
  for index = start + width; index < line.length(); index = index + 1 {
    if !align_clustal_is_whitespace(line.unsafe_get(index).to_int()) {
      raise AlignClustalError(
        "CLUSTAL consensus has non-whitespace trailing data",
      )
    }
  }
  output.to_string()
}

///|
fn align_clustal_count_pair(
  alignment : AlignClustalAlignment,
  first_row : Int,
  second_row : Int,
) -> AlignClustalCounts {
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
  AlignClustalCounts::{
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
fn align_clustal_movement_changes(
  alignment : AlignClustalAlignment,
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
fn align_clustal_column_has_gap(
  alignment : AlignClustalAlignment,
  column : Int,
) -> Bool {
  for sequence in alignment.sequences {
    if sequence.aligned_sequence.unsafe_get(column).to_int() == '-'.to_int() {
      return true
    }
  }
  false
}

///|
fn align_clustal_column_identical(
  alignment : AlignClustalAlignment,
  column : Int,
) -> Bool {
  let first = align_clustal_upper_code(
    alignment.sequences[0].aligned_sequence.unsafe_get(column).to_int(),
  )
  for index = 1; index < alignment.sequences.length(); index = index + 1 {
    let code = align_clustal_upper_code(
      alignment.sequences[index].aligned_sequence.unsafe_get(column).to_int(),
    )
    if code != first {
      return false
    }
  }
  true
}

///|
fn align_clustal_column_in_any_group(
  alignment : AlignClustalAlignment,
  column : Int,
  groups : Array[String],
) -> Bool {
  for group in groups {
    let mut all_present = true
    for sequence in alignment.sequences {
      let code = align_clustal_upper_code(
        sequence.aligned_sequence.unsafe_get(column).to_int(),
      )
      if !align_clustal_string_has_code(group, code) {
        all_present = false
        break
      }
    }
    if all_present {
      return true
    }
  }
  false
}

///|
fn align_clustal_is_nucleotide(alignment : AlignClustalAlignment) -> Bool {
  let symbols = "ACGTUNRYKMSWBDHV"
  for sequence in alignment.sequences {
    for index = 0; index < sequence.sequence.length(); index = index + 1 {
      let code = align_clustal_upper_code(
        sequence.sequence.unsafe_get(index).to_int(),
      )
      if code != '?'.to_int() && !align_clustal_string_has_code(symbols, code) {
        return false
      }
    }
  }
  true
}

///|
fn align_clustal_validate_program(
  program : String,
) -> Unit raise AlignClustalError {
  if program != "CLUSTAL" &&
    program != "PROBCONS" &&
    program != "MUSCLE" &&
    program != "MSAPROBS" &&
    program != "Kalign" &&
    program != "Biopython" {
    raise AlignClustalError(
      "Unknown CLUSTAL generator '" +
      program +
      "'; expected CLUSTAL, PROBCONS, MUSCLE, MSAPROBS, Kalign, or Biopython",
    )
  }
}

///|
fn align_clustal_validate_version(
  version : String,
) -> Unit raise AlignClustalError {
  for index = 0; index < version.length(); index = index + 1 {
    let code = version.unsafe_get(index).to_int()
    if align_clustal_is_whitespace(code) || code < 33 || code > 126 {
      raise AlignClustalError("CLUSTAL version must be one printable token")
    }
  }
}

///|
fn align_clustal_validate_id(id : String) -> Unit raise AlignClustalError {
  if id.length() == 0 {
    raise AlignClustalError("CLUSTAL sequence identifier must not be empty")
  }
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    if align_clustal_is_whitespace(code) || code < 33 || code > 126 {
      raise AlignClustalError(
        "CLUSTAL sequence identifier contains whitespace or a control character",
      )
    }
  }
}

///|
fn align_clustal_validate_segment(
  segment : String,
) -> Unit raise AlignClustalError {
  if segment.length() == 0 {
    raise AlignClustalError("CLUSTAL sequence segment must not be empty")
  }
  for index = 0; index < segment.length(); index = index + 1 {
    let code = segment.unsafe_get(index).to_int()
    let letter = (code >= 'A'.to_int() && code <= 'Z'.to_int()) ||
      (code >= 'a'.to_int() && code <= 'z'.to_int())
    if !letter &&
      code != '-'.to_int() &&
      code != '?'.to_int() &&
      code != '*'.to_int() {
      raise AlignClustalError(
        "CLUSTAL sequence contains an invalid residue or gap symbol",
      )
    }
  }
}

///|
fn align_clustal_validate_consensus(
  consensus : String,
) -> Unit raise AlignClustalError {
  for index = 0; index < consensus.length(); index = index + 1 {
    let code = consensus.unsafe_get(index).to_int()
    if code != ' '.to_int() &&
      code != '*'.to_int() &&
      code != ':'.to_int() &&
      code != '.'.to_int() {
      raise AlignClustalError("CLUSTAL consensus contains an invalid symbol")
    }
  }
}

///|
fn align_clustal_validate_plain_text(
  value : String,
  label : String,
) -> Unit raise AlignClustalError {
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    if code == '\n'.to_int() || code == '\r'.to_int() || code < 32 {
      raise AlignClustalError(
        label + " contains a line break or control character",
      )
    }
  }
}

///|
fn align_clustal_validate_row(
  alignment : AlignClustalAlignment,
  row : Int,
) -> Unit raise AlignClustalError {
  if row < 0 || row >= alignment.sequences.length() {
    raise AlignClustalError("CLUSTAL row index is out of bounds")
  }
}

///|
fn align_clustal_canonical_header(program : String, version : String) -> String {
  if version.length() == 0 {
    program + " multiple sequence alignment"
  } else {
    program + " " + version + " multiple sequence alignment"
  }
}

///|
fn align_clustal_lines(text : String) -> Array[String] {
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
fn align_clustal_split_whitespace(value : String) -> Array[String] {
  let fields : Array[String] = []
  let mut start = 0
  let mut in_field = false
  for index = 0; index < value.length(); index = index + 1 {
    let whitespace = align_clustal_is_whitespace(
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
fn align_clustal_segment_start(line : String) -> Int raise AlignClustalError {
  let mut index = 0
  while index < line.length() &&
        !align_clustal_is_whitespace(line.unsafe_get(index).to_int()) {
    index = index + 1
  }
  while index < line.length() &&
        align_clustal_is_whitespace(line.unsafe_get(index).to_int()) {
    index = index + 1
  }
  if index >= line.length() {
    raise AlignClustalError("CLUSTAL sequence line has no sequence segment")
  }
  index
}

///|
fn align_clustal_parse_nonnegative_int(
  text : String,
  label : String,
) -> Int raise AlignClustalError {
  if text.length() == 0 {
    raise AlignClustalError(label + " is empty")
  }
  let mut value = 0
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise AlignClustalError(label + " must be a non-negative integer")
    }
    let digit = code - '0'.to_int()
    if value > 214748364 || (value == 214748364 && digit > 7) {
      raise AlignClustalError(label + " is outside the supported integer range")
    }
    value = value * 10 + digit
  }
  value
}

///|
fn align_clustal_residue_count(segment : String) -> Int {
  let mut count = 0
  for index = 0; index < segment.length(); index = index + 1 {
    if segment.unsafe_get(index).to_int() != '-'.to_int() {
      count = count + 1
    }
  }
  count
}

///|
fn align_clustal_remove_gaps(value : String) -> String {
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
fn align_clustal_strip_parentheses(value : String) -> String {
  if value.length() >= 2 &&
    value.unsafe_get(0).to_int() == '('.to_int() &&
    value.unsafe_get(value.length() - 1).to_int() == ')'.to_int() {
    value[1:value.length() - 1].to_owned()
  } else {
    value
  }
}

///|
fn align_clustal_upper_code(code : Int) -> Int {
  if code >= 'a'.to_int() && code <= 'z'.to_int() {
    code - 'a'.to_int() + 'A'.to_int()
  } else {
    code
  }
}

///|
fn align_clustal_string_has_code(value : String, code : Int) -> Bool {
  for index = 0; index < value.length(); index = index + 1 {
    if value.unsafe_get(index).to_int() == code {
      return true
    }
  }
  false
}

///|
fn align_clustal_find_string(values : Array[String], target : String) -> Int {
  align_clustal_find_string_from(values, target, 0)
}

///|
fn align_clustal_find_string_from(
  values : Array[String],
  target : String,
  start : Int,
) -> Int {
  for index = start; index < values.length(); index = index + 1 {
    if values[index] == target {
      return index
    }
  }
  -1
}

///|
fn align_clustal_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\r'.to_int() ||
  code == '\n'.to_int()
}
