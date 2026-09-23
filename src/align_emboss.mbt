// Alignment-aware EMBOSS srspair/pair/simple support.
//
// This module follows Biopython 1.86 Bio.Align.emboss semantics. EMBOSS
// sequence lines use a fixed 21-character prefix and one-based inclusive
// printed coordinates. Public coordinates are normalized to zero-based
// residue positions and boundary coordinates.

///|
/// Error raised for malformed EMBOSS alignment output or invalid operations.
pub suberror AlignEmbossError {
  AlignEmbossError(String)
}

///|
/// Orientation of a sequence row in an EMBOSS alignment.
pub(all) enum AlignEmbossStrand {
  AlignEmbossForward
  AlignEmbossReverse
} derive(Eq, Debug)
pub extend AlignEmbossStrand with Eq::{equal, not_equal}
pub extend AlignEmbossStrand with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// File-level metadata from the EMBOSS report header.
pub struct AlignEmbossMetadata {
  program : String
  rundate : String
  report_file : String
  align_format : String
  command_line : String
} derive(Eq, Debug)
pub extend AlignEmbossMetadata with Eq::{equal, not_equal}
pub extend AlignEmbossMetadata with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Per-alignment annotations reported by EMBOSS.
pub struct AlignEmbossAnnotations {
  matrix : String
  gap_penalty : Double?
  extend_penalty : Double?
  length : Int
  identity : Int?
  similarity : Int?
  gaps : Int?
  score : Double?
  longest_identity : String
  longest_similarity : String
  shortest_identity : String
  shortest_similarity : String
} derive(Eq, Debug)
pub extend AlignEmbossAnnotations with Eq::{equal, not_equal}
pub extend AlignEmbossAnnotations with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One sequence row in an EMBOSS alignment.
///
/// `start` and `end` are zero-based boundary coordinates. Forward rows have
/// `start <= end`; reverse rows have `start > end`. `aligned_sequence` is in
/// displayed alignment orientation and `sequence` is the same row without
/// gap characters.
pub struct AlignEmbossSequence {
  id : String
  sequence : String
  aligned_sequence : String
  start : Int
  end : Int
  strand : AlignEmbossStrand
} derive(Eq, Debug)
pub extend AlignEmbossSequence with Eq::{equal, not_equal}
pub extend AlignEmbossSequence with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One parsed EMBOSS pairwise or multiple sequence alignment.
pub struct AlignEmbossAlignment {
  sequences : Array[AlignEmbossSequence]
  consensus : String
  annotations : AlignEmbossAnnotations
} derive(Eq, Debug)
pub extend AlignEmbossAlignment with Eq::{equal, not_equal}
pub extend AlignEmbossAlignment with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// A complete EMBOSS report containing one or more alignments.
pub struct AlignEmbossDocument {
  metadata : AlignEmbossMetadata
  alignments : Array[AlignEmbossAlignment]
} derive(Eq, Debug)
pub extend AlignEmbossDocument with Eq::{equal, not_equal}
pub extend AlignEmbossDocument with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Pairwise statistics for two rows in an EMBOSS alignment.
pub struct AlignEmbossPairCounts {
  columns : Int
  aligned : Int
  identities : Int
  mismatches : Int
  insertions : Int
  deletions : Int
  double_gap_columns : Int
  insertion_opens : Int
  deletion_opens : Int
  positives : Int
} derive(Eq, Debug)
pub extend AlignEmbossPairCounts with Eq::{equal, not_equal}
pub extend AlignEmbossPairCounts with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Construct file-level EMBOSS metadata.
pub fn AlignEmbossMetadata::create(
  program? : String = "",
  rundate? : String = "",
  report_file? : String = "",
  align_format? : String = "srspair",
  command_line? : String = "",
) -> AlignEmbossMetadata raise AlignEmbossError {
  align_emboss_validate_header_value(program, "Program")
  align_emboss_validate_header_value(rundate, "Rundate")
  align_emboss_validate_header_value(report_file, "Report_file")
  align_emboss_validate_header_value(command_line, "Command line")
  if align_format != "srspair" &&
    align_format != "pair" &&
    align_format != "simple" {
    raise AlignEmbossError(
      "EMBOSS Align_format must be srspair, pair, or simple",
    )
  }
  AlignEmbossMetadata::{
    program,
    rundate,
    report_file,
    align_format,
    command_line,
  }
}

///|
/// Construct per-alignment EMBOSS annotations.
pub fn AlignEmbossAnnotations::create(
  length : Int,
  matrix? : String = "",
  gap_penalty? : Double? = None,
  extend_penalty? : Double? = None,
  identity? : Int? = None,
  similarity? : Int? = None,
  gaps? : Int? = None,
  score? : Double? = None,
  longest_identity? : String = "",
  longest_similarity? : String = "",
  shortest_identity? : String = "",
  shortest_similarity? : String = "",
) -> AlignEmbossAnnotations raise AlignEmbossError {
  if length <= 0 {
    raise AlignEmbossError("EMBOSS alignment length must be positive")
  }
  align_emboss_validate_header_value(matrix, "Matrix")
  align_emboss_validate_optional_nonnegative_double(gap_penalty, "Gap_penalty")
  align_emboss_validate_optional_nonnegative_double(
    extend_penalty, "Extend_penalty",
  )
  align_emboss_validate_optional_count(identity, length, "Identity")
  align_emboss_validate_optional_count(similarity, length, "Similarity")
  align_emboss_validate_optional_count(gaps, length, "Gaps")
  match score {
    Some(value) =>
      if value.abs() > 1.0e300 {
        raise AlignEmbossError("EMBOSS Score must be finite")
      }
    None => ()
  }
  AlignEmbossAnnotations::{
    matrix,
    gap_penalty,
    extend_penalty,
    length,
    identity,
    similarity,
    gaps,
    score,
    longest_identity,
    longest_similarity,
    shortest_identity,
    shortest_similarity,
  }
}

///|
/// Construct one normalized EMBOSS sequence row.
pub fn AlignEmbossSequence::create(
  id : String,
  aligned_sequence : String,
  start : Int,
  end : Int,
) -> AlignEmbossSequence raise AlignEmbossError {
  align_emboss_validate_id(id)
  if aligned_sequence.length() == 0 {
    raise AlignEmbossError("EMBOSS aligned sequence must not be empty")
  }
  align_emboss_validate_aligned_fragment(aligned_sequence, "sequence row")
  if start < 0 || end < 0 {
    raise AlignEmbossError("EMBOSS sequence coordinates must be non-negative")
  }
  let sequence = align_emboss_remove_gaps(aligned_sequence)
  let span = if start >= end { start - end } else { end - start }
  if sequence.length() != span {
    raise AlignEmbossError(
      "EMBOSS ungapped sequence length does not match its coordinate span",
    )
  }
  let strand = if start <= end {
    AlignEmbossForward
  } else {
    AlignEmbossReverse
  }
  AlignEmbossSequence::{ id, sequence, aligned_sequence, start, end, strand, }
}

///|
/// Construct one EMBOSS alignment.
pub fn AlignEmbossAlignment::create(
  sequences : Array[AlignEmbossSequence],
  annotations : AlignEmbossAnnotations,
  consensus? : String = "",
) -> AlignEmbossAlignment raise AlignEmbossError {
  let copied : Array[AlignEmbossSequence] = []
  for sequence in sequences {
    copied.push(sequence)
  }
  let alignment = AlignEmbossAlignment::{
    sequences: copied,
    consensus,
    annotations,
  }
  align_emboss_validate_alignment(alignment)
  alignment
}

///|
/// Construct a complete EMBOSS report.
pub fn AlignEmbossDocument::create(
  metadata : AlignEmbossMetadata,
  alignments : Array[AlignEmbossAlignment],
) -> AlignEmbossDocument raise AlignEmbossError {
  if alignments.length() == 0 {
    raise AlignEmbossError("EMBOSS document must contain an alignment")
  }
  let copied : Array[AlignEmbossAlignment] = []
  for alignment in alignments {
    align_emboss_validate_alignment(alignment)
    copied.push(alignment)
  }
  AlignEmbossDocument::{ metadata, alignments: copied, }
}

///|
/// Parse an EMBOSS srspair, pair, or simple alignment report.
///
/// The parser accepts LF or CRLF input and validates file/alignment headers,
/// declared row order, fixed sequence-line fields, block widths, coordinate
/// continuity, final alignment width, consensus width, and reported counts.
pub fn align_emboss_parse(
  text : String,
) -> AlignEmbossDocument raise AlignEmbossError {
  let lines = align_emboss_lines(text)
  if lines.length() == 0 || (lines.length() == 1 && lines[0].length() == 0) {
    raise AlignEmbossError("Empty EMBOSS input")
  }
  let header_divider = "########################################"
  if lines[0].trim().to_owned() != header_divider {
    raise AlignEmbossError("EMBOSS file is missing its header divider")
  }
  let mut program = ""
  let mut rundate = ""
  let mut report_file = ""
  let mut align_format = "srspair"
  let mut command_line = ""
  let mut index = 1
  let mut header_closed = false
  while index < lines.length() {
    let line = lines[index]
    if line.trim().to_owned() == header_divider {
      header_closed = true
      index = index + 1
      break
    }
    if align_emboss_has_prefix(line, "#    ") {
      if command_line.length() == 0 {
        raise AlignEmbossError(
          "EMBOSS command continuation appears before Commandline",
        )
      }
      command_line = command_line + " " + line[1:].trim().to_owned()
      index = index + 1
      continue
    }
    if !align_emboss_has_prefix(line, "# ") {
      raise AlignEmbossError(
        "Unexpected EMBOSS header line " + (index + 1).to_string(),
      )
    }
    let content = line[2:].to_owned()
    let (key, value) = align_emboss_split_key_value(
      content,
      "header",
      index + 1,
    )
    match key {
      "Program" => program = value
      "Rundate" => rundate = value
      "Report_file" => report_file = value
      "Align_format" => align_format = value
      "Commandline" => command_line = value
      _ => ()
    }
    index = index + 1
  }
  if !header_closed {
    raise AlignEmbossError("Truncated EMBOSS file header")
  }
  let metadata = AlignEmbossMetadata::create(
    program~,
    rundate~,
    report_file~,
    align_format~,
    command_line~,
  )
  let alignments : Array[AlignEmbossAlignment] = []
  while index < lines.length() {
    let line = lines[index].trim().to_owned()
    if line.length() == 0 || align_emboss_is_separator(line) {
      index = index + 1
      continue
    }
    if line != "#=======================================" {
      raise AlignEmbossError(
        "Unexpected content before EMBOSS alignment at line " +
        (index + 1).to_string(),
      )
    }
    let (alignment, next_index) = align_emboss_parse_alignment(
      lines, index, metadata,
    )
    alignments.push(alignment)
    index = next_index
  }
  AlignEmbossDocument::create(metadata, alignments)
}

///|
/// Serialize an EMBOSS document in a canonical parseable form.
///
/// The upstream Biopython module is read-only; this writer is a MoonBit
/// extension for reproducible round trips. Sequence blocks default to the
/// EMBOSS width of 50 columns.
pub fn align_emboss_write(
  document : AlignEmbossDocument,
  line_width? : Int = 50,
) -> String raise AlignEmbossError {
  if line_width <= 0 || line_width > 50 {
    raise AlignEmbossError("EMBOSS line width must be between 1 and 50")
  }
  if document.alignments.length() == 0 {
    raise AlignEmbossError("EMBOSS document must contain an alignment")
  }
  let output = StringBuilder::new()
  output.write_string("########################################\n")
  if document.metadata.program.length() > 0 {
    output.write_string("# Program: " + document.metadata.program + "\n")
  }
  if document.metadata.rundate.length() > 0 {
    output.write_string("# Rundate: " + document.metadata.rundate + "\n")
  }
  if document.metadata.command_line.length() > 0 {
    output.write_string(
      "# Commandline: " + document.metadata.command_line + "\n",
    )
  }
  output.write_string(
    "# Align_format: " + document.metadata.align_format + "\n",
  )
  if document.metadata.report_file.length() > 0 {
    output.write_string(
      "# Report_file: " + document.metadata.report_file + "\n",
    )
  }
  output.write_string("########################################\n")
  for alignment in document.alignments {
    align_emboss_validate_alignment(alignment)
    output.write_string("#=======================================\n#\n")
    output.write_string(
      "# Aligned_sequences: " + alignment.sequences.length().to_string() + "\n",
    )
    for row = 0; row < alignment.sequences.length(); row = row + 1 {
      output.write_string(
        "# " + (row + 1).to_string() + ": " + alignment.sequences[row].id + "\n",
      )
    }
    if alignment.annotations.matrix.length() > 0 {
      output.write_string("# Matrix: " + alignment.annotations.matrix + "\n")
    }
    align_emboss_write_optional_double(
      output,
      "Gap_penalty",
      alignment.annotations.gap_penalty,
    )
    align_emboss_write_optional_double(
      output,
      "Extend_penalty",
      alignment.annotations.extend_penalty,
    )
    output.write_string("#\n")
    output.write_string(
      "# Length: " + alignment.annotations.length.to_string() + "\n",
    )
    align_emboss_write_optional_count(
      output,
      "Identity",
      alignment.annotations.identity,
      alignment.annotations.length,
    )
    align_emboss_write_optional_count(
      output,
      "Similarity",
      alignment.annotations.similarity,
      alignment.annotations.length,
    )
    align_emboss_write_optional_count(
      output,
      "Gaps",
      alignment.annotations.gaps,
      alignment.annotations.length,
    )
    match alignment.annotations.score {
      Some(value) => output.write_string("# Score: " + value.to_string() + "\n")
      None => ()
    }
    align_emboss_write_long_annotation(
      output,
      "Longest_Identity",
      alignment.annotations.longest_identity,
    )
    align_emboss_write_long_annotation(
      output,
      "Longest_Similarity",
      alignment.annotations.longest_similarity,
    )
    align_emboss_write_long_annotation(
      output,
      "Shortest_Identity",
      alignment.annotations.shortest_identity,
    )
    align_emboss_write_long_annotation(
      output,
      "Shortest_Similarity",
      alignment.annotations.shortest_similarity,
    )
    output.write_string("#\n#=======================================\n")
    let width = alignment.alignment_length()
    let mut column = 0
    while column < width {
      let block_end = if column + line_width < width {
        column + line_width
      } else {
        width
      }
      for row = 0; row < alignment.sequences.length(); row = row + 1 {
        let sequence = alignment.sequences[row]
        let fragment = sequence.aligned_sequence[column:block_end].to_owned()
        let boundary = align_emboss_boundary_at(sequence, column)
        let residues = align_emboss_count_residues(fragment)
        let (printed_start, printed_end) = align_emboss_printed_coordinates(
          sequence, boundary, residues,
        )
        output.write_string(
          align_emboss_sequence_prefix(sequence.id, printed_start) +
          fragment +
          " " +
          align_emboss_pad_left(printed_end.to_string(), 7) +
          "\n",
        )
      }
      if alignment.consensus.length() > 0 {
        output.write_string(
          " ".repeat(21) +
          alignment.consensus[column:block_end].to_owned() +
          "\n",
        )
      }
      output.write_char('\n')
      column = block_end
    }
    output.write_string(
      "#---------------------------------------\n" +
      "#---------------------------------------\n",
    )
  }
  output.to_string()
}

///|
/// Return the number of alignments in a report.
pub fn AlignEmbossDocument::num_alignments(self : AlignEmbossDocument) -> Int {
  self.alignments.length()
}

///|
/// Return an alignment by zero-based index.
pub fn AlignEmbossDocument::get(
  self : AlignEmbossDocument,
  index : Int,
) -> AlignEmbossAlignment? {
  if index < 0 || index >= self.alignments.length() {
    None
  } else {
    Some(self.alignments[index])
  }
}

///|
/// Return the number of sequence rows.
pub fn AlignEmbossAlignment::num_sequences(self : AlignEmbossAlignment) -> Int {
  self.sequences.length()
}

///|
/// Return the alignment width including gap columns.
pub fn AlignEmbossAlignment::alignment_length(
  self : AlignEmbossAlignment,
) -> Int {
  self.annotations.length
}

///|
/// Find the first row with an exact declared identifier.
pub fn AlignEmbossAlignment::find_sequence(
  self : AlignEmbossAlignment,
  id : String,
) -> Int? {
  for row = 0; row < self.sequences.length(); row = row + 1 {
    if self.sequences[row].id == id {
      return Some(row)
    }
  }
  None
}

///|
/// Return all row characters at an alignment column.
pub fn AlignEmbossAlignment::column(
  self : AlignEmbossAlignment,
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
/// Return the zero-based absolute residue position at an alignment column.
///
/// Gap columns and invalid indices return `None`.
pub fn AlignEmbossAlignment::column_to_sequence_position(
  self : AlignEmbossAlignment,
  row : Int,
  column : Int,
) -> Int? {
  if row < 0 ||
    row >= self.sequences.length() ||
    column < 0 ||
    column >= self.alignment_length() {
    return None
  }
  let sequence = self.sequences[row]
  if sequence.aligned_sequence.unsafe_get(column).to_int() == '-'.to_int() {
    return None
  }
  let boundary = align_emboss_boundary_at(sequence, column)
  if sequence.strand == AlignEmbossForward {
    Some(boundary)
  } else {
    Some(boundary - 1)
  }
}

///|
/// Map a zero-based absolute residue position to its alignment column.
pub fn AlignEmbossAlignment::sequence_position_to_column(
  self : AlignEmbossAlignment,
  row : Int,
  position : Int,
) -> Int? {
  if row < 0 || row >= self.sequences.length() || position < 0 {
    return None
  }
  for column = 0; column < self.alignment_length(); column = column + 1 {
    match self.column_to_sequence_position(row, column) {
      Some(value) => if value == position { return Some(column) }
      None => ()
    }
  }
  None
}

///|
/// Map a residue position from one row through the alignment to another row.
pub fn AlignEmbossAlignment::map_position(
  self : AlignEmbossAlignment,
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
/// Return per-column absolute coordinates for two rows.
///
/// Gaps are represented by `None`.
pub fn AlignEmbossAlignment::aligned_pairs(
  self : AlignEmbossAlignment,
  first_row : Int,
  second_row : Int,
) -> Array[(Int?, Int?)] raise AlignEmbossError {
  align_emboss_validate_row(self, first_row)
  align_emboss_validate_row(self, second_row)
  let pairs : Array[(Int?, Int?)] = []
  for column = 0; column < self.alignment_length(); column = column + 1 {
    pairs.push(
      (
        self.column_to_sequence_position(first_row, column),
        self.column_to_sequence_position(second_row, column),
      ),
    )
  }
  pairs
}

///|
/// Build a compact Biopython-style boundary coordinate path.
///
/// The outer array contains one coordinate array per row. Breakpoints are
/// emitted whenever the per-column movement vector changes.
pub fn AlignEmbossAlignment::coordinate_path(
  self : AlignEmbossAlignment,
) -> Array[Array[Int]] {
  let result : Array[Array[Int]] = []
  for sequence in self.sequences {
    result.push([sequence.start])
  }
  if self.alignment_length() == 0 {
    return result
  }
  for column = 1; column < self.alignment_length(); column = column + 1 {
    let mut changed = false
    for row = 0; row < self.sequences.length(); row = row + 1 {
      let sequence = self.sequences[row]
      let previous = align_emboss_column_step(sequence, column - 1)
      let current = align_emboss_column_step(sequence, column)
      if previous != current {
        changed = true
        break
      }
    }
    if changed {
      for row = 0; row < self.sequences.length(); row = row + 1 {
        result[row].push(align_emboss_boundary_at(self.sequences[row], column))
      }
    }
  }
  for row = 0; row < self.sequences.length(); row = row + 1 {
    result[row].push(self.sequences[row].end)
  }
  result
}

///|
/// Compute pairwise residue, gap, gap-open, and positive-column counts.
///
/// A gap in the first row is an insertion; a gap in the second row is a
/// deletion. `positives` uses `|` and `:` from the EMBOSS consensus when
/// available, otherwise it equals exact identities.
pub fn AlignEmbossAlignment::pair_counts(
  self : AlignEmbossAlignment,
  first_row : Int,
  second_row : Int,
) -> AlignEmbossPairCounts raise AlignEmbossError {
  align_emboss_validate_row(self, first_row)
  align_emboss_validate_row(self, second_row)
  let first = self.sequences[first_row].aligned_sequence
  let second = self.sequences[second_row].aligned_sequence
  let mut aligned = 0
  let mut identities = 0
  let mut mismatches = 0
  let mut insertions = 0
  let mut deletions = 0
  let mut double_gap_columns = 0
  let mut insertion_opens = 0
  let mut deletion_opens = 0
  let mut positives = 0
  let mut in_insertion = false
  let mut in_deletion = false
  for column = 0; column < self.alignment_length(); column = column + 1 {
    let first_code = first.unsafe_get(column).to_int()
    let second_code = second.unsafe_get(column).to_int()
    let first_gap = first_code == '-'.to_int()
    let second_gap = second_code == '-'.to_int()
    if first_gap && second_gap {
      double_gap_columns = double_gap_columns + 1
      in_insertion = false
      in_deletion = false
    } else if first_gap {
      insertions = insertions + 1
      if !in_insertion {
        insertion_opens = insertion_opens + 1
      }
      in_insertion = true
      in_deletion = false
    } else if second_gap {
      deletions = deletions + 1
      if !in_deletion {
        deletion_opens = deletion_opens + 1
      }
      in_deletion = true
      in_insertion = false
    } else {
      aligned = aligned + 1
      if first_code == second_code {
        identities = identities + 1
      } else {
        mismatches = mismatches + 1
      }
      if self.consensus.length() == self.alignment_length() {
        let marker = self.consensus.unsafe_get(column).to_int()
        if marker == '|'.to_int() || marker == ':'.to_int() {
          positives = positives + 1
        }
      } else if first_code == second_code {
        positives = positives + 1
      }
      in_insertion = false
      in_deletion = false
    }
  }
  AlignEmbossPairCounts::{
    columns: self.alignment_length(),
    aligned,
    identities,
    mismatches,
    insertions,
    deletions,
    double_gap_columns,
    insertion_opens,
    deletion_opens,
    positives,
  }
}

///|
/// Return exact identity divided by columns where both rows have residues.
pub fn AlignEmbossPairCounts::identity(self : AlignEmbossPairCounts) -> Double {
  if self.aligned == 0 {
    0.0
  } else {
    self.identities.to_double() / self.aligned.to_double()
  }
}

///|
/// Return a compact alignment description.
pub fn AlignEmbossAlignment::summary(self : AlignEmbossAlignment) -> String {
  "EMBOSS alignment(rows=" +
  self.sequences.length().to_string() +
  ", columns=" +
  self.alignment_length().to_string() +
  ", matrix=" +
  self.annotations.matrix +
  ")"
}

///|
/// Return a compact report description.
pub fn AlignEmbossDocument::summary(self : AlignEmbossDocument) -> String {
  "EMBOSS report(program=" +
  self.metadata.program +
  ", format=" +
  self.metadata.align_format +
  ", alignments=" +
  self.alignments.length().to_string() +
  ")"
}

///|
/// Small EMBOSS srspair document used by tests and examples.
pub fn align_emboss_example_text() -> String {
  "########################################\n" +
  "# Program: water\n" +
  "# Rundate: Wed Jan 16 17:23:19 2002\n" +
  "# Commandline: water\n" +
  "#    -asequence reference.fa\n" +
  "#    -bsequence query.fa\n" +
  "# Align_format: srspair\n" +
  "# Report_file: stdout\n" +
  "########################################\n" +
  "#=======================================\n" +
  "#\n" +
  "# Aligned_sequences: 2\n" +
  "# 1: reference_sequence\n" +
  "# 2: query_sequence\n" +
  "# Matrix: EDNAFULL\n" +
  "# Gap_penalty: 10.0\n" +
  "# Extend_penalty: 0.5\n" +
  "#\n" +
  "# Length: 18\n" +
  "# Identity: 16/18 (88.9%)\n" +
  "# Similarity: 17/18 (94.4%)\n" +
  "# Gaps: 1/18 ( 5.6%)\n" +
  "# Score: 72.5\n" +
  "#\n" +
  "#=======================================\n" +
  "reference_seq     79 ACGTTGAGT-CTGGGATG     95\n" +
  "                     ||||||||| ||||:|||\n" +
  "query_sequence     1 ACGTTGAGTACTGGAATG     18\n" +
  "#---------------------------------------\n" +
  "#---------------------------------------\n"
}

///|
fn align_emboss_parse_alignment(
  lines : Array[String],
  start_index : Int,
  metadata : AlignEmbossMetadata,
) -> (AlignEmbossAlignment, Int) raise AlignEmbossError {
  let mut index = start_index + 1
  let identifiers : Array[String] = []
  let mut number_of_sequences = -1
  let mut length = -1
  let mut matrix = ""
  let mut gap_penalty : Double? = None
  let mut extend_penalty : Double? = None
  let mut identity : Int? = None
  let mut similarity : Int? = None
  let mut gaps : Int? = None
  let mut score : Double? = None
  let mut longest_identity = ""
  let mut longest_similarity = ""
  let mut shortest_identity = ""
  let mut shortest_similarity = ""
  let mut metadata_closed = false
  while index < lines.length() {
    let line = lines[index]
    let trimmed = line.trim().to_owned()
    if trimmed == "#=======================================" {
      metadata_closed = true
      index = index + 1
      break
    }
    if trimmed == "#" || trimmed.length() == 0 {
      index = index + 1
      continue
    }
    if !align_emboss_has_prefix(line, "# ") {
      raise AlignEmbossError(
        "Unexpected EMBOSS alignment header line " + (index + 1).to_string(),
      )
    }
    let content = line[2:].to_owned()
    let (key, value) = align_emboss_split_annotation(content, index + 1)
    if key == "Aligned_sequences" {
      if number_of_sequences >= 0 {
        raise AlignEmbossError("Duplicate EMBOSS Aligned_sequences")
      }
      number_of_sequences = align_emboss_parse_positive_int(
        value,
        "Aligned_sequences",
        index + 1,
      )
      index = index + 1
      for expected = 1; expected <= number_of_sequences; expected = expected + 1 {
        if index >= lines.length() ||
          !align_emboss_has_prefix(lines[index], "# ") {
          raise AlignEmbossError("Truncated EMBOSS sequence identifier list")
        }
        let (number_text, identifier) = align_emboss_split_key_value(
          lines[index][2:].to_owned(),
          "sequence identifier",
          index + 1,
        )
        let number = align_emboss_parse_positive_int(
          number_text,
          "sequence number",
          index + 1,
        )
        if number != expected {
          raise AlignEmbossError(
            "Non-contiguous EMBOSS sequence identifier numbering",
          )
        }
        align_emboss_validate_id(identifier)
        identifiers.push(identifier)
        index = index + 1
      }
      continue
    }
    match key {
      "Matrix" => matrix = value
      "Gap_penalty" =>
        gap_penalty = Some(
          align_emboss_parse_number(value, "Gap_penalty", index + 1),
        )
      "Extend_penalty" =>
        extend_penalty = Some(
          align_emboss_parse_number(value, "Extend_penalty", index + 1),
        )
      "Length" =>
        length = align_emboss_parse_positive_int(value, "Length", index + 1)
      "Identity" =>
        identity = Some(
          align_emboss_parse_fraction_count(value, "Identity", index + 1),
        )
      "Similarity" =>
        similarity = Some(
          align_emboss_parse_fraction_count(value, "Similarity", index + 1),
        )
      "Gaps" =>
        gaps = Some(align_emboss_parse_fraction_count(value, "Gaps", index + 1))
      "Score" =>
        score = Some(align_emboss_parse_number(value, "Score", index + 1))
      "Longest_Identity" => longest_identity = value
      "Longest_Similarity" => longest_similarity = value
      "Shortest_Identity" => shortest_identity = value
      "Shortest_Similarity" => shortest_similarity = value
      _ =>
        raise AlignEmbossError(
          "Unknown EMBOSS alignment annotation '" + key + "'",
        )
    }
    index = index + 1
  }
  if !metadata_closed {
    raise AlignEmbossError("Truncated EMBOSS alignment header")
  }
  if number_of_sequences <= 0 || identifiers.length() != number_of_sequences {
    raise AlignEmbossError("Number of EMBOSS sequences is missing")
  }
  if length <= 0 {
    raise AlignEmbossError("Length of EMBOSS alignment is missing")
  }
  let annotations = AlignEmbossAnnotations::create(
    length,
    matrix~,
    gap_penalty~,
    extend_penalty~,
    identity~,
    similarity~,
    gaps~,
    score~,
    longest_identity~,
    longest_similarity~,
    shortest_identity~,
    shortest_similarity~,
  )
  let aligned_rows : Array[String] = []
  let ungapped_rows : Array[String] = []
  let starts : Array[Int] = []
  let ends : Array[Int] = []
  let directions : Array[Int] = []
  for _row = 0; _row < number_of_sequences; _row = _row + 1 {
    aligned_rows.push("")
    ungapped_rows.push("")
    starts.push(0)
    ends.push(0)
    directions.push(0)
  }
  let consensus_builder = StringBuilder::new(size_hint=length)
  let mut consensus_seen = false
  let mut row_index = 0
  let mut block_width = 0
  while index < lines.length() {
    let line = lines[index]
    let trimmed = line.trim().to_owned()
    if align_emboss_is_separator(trimmed) {
      index = index + 1
      break
    }
    if trimmed.length() == 0 && line.length() < 21 {
      if row_index == number_of_sequences {
        row_index = 0
        block_width = 0
      }
      index = index + 1
      continue
    }
    if line.length() < 21 {
      raise AlignEmbossError(
        "Malformed EMBOSS alignment body at line " + (index + 1).to_string(),
      )
    }
    let prefix = line[0:21].to_owned().trim().to_owned()
    if prefix.length() == 0 {
      if block_width <= 0 {
        raise AlignEmbossError(
          "EMBOSS consensus appears before a sequence block",
        )
      }
      let raw = line[21:].to_owned()
      let chunk = align_emboss_consensus_chunk(raw, block_width, index + 1)
      consensus_builder.write_string(chunk)
      consensus_seen = true
      index = index + 1
      continue
    }
    if row_index == number_of_sequences {
      row_index = 0
      block_width = 0
    }
    let prefix_words = align_emboss_split_whitespace(prefix)
    if prefix_words.length() != 2 {
      raise AlignEmbossError(
        "Malformed EMBOSS sequence prefix at line " + (index + 1).to_string(),
      )
    }
    let printed_id = prefix_words[0]
    if !align_emboss_has_prefix(identifiers[row_index], printed_id) {
      raise AlignEmbossError(
        "Unexpected EMBOSS sequence identifier at line " +
        (index + 1).to_string(),
      )
    }
    let printed_start = align_emboss_parse_nonnegative_int(
      prefix_words[1],
      "sequence start",
      index + 1,
    )
    let suffix_words = align_emboss_split_whitespace(line[21:].to_owned())
    if suffix_words.length() != 2 {
      raise AlignEmbossError(
        "Malformed EMBOSS sequence fields at line " + (index + 1).to_string(),
      )
    }
    let fragment = suffix_words[0]
    align_emboss_validate_aligned_fragment(
      fragment,
      "line " + (index + 1).to_string(),
    )
    let printed_end = align_emboss_parse_nonnegative_int(
      suffix_words[1],
      "sequence end",
      index + 1,
    )
    if block_width == 0 {
      block_width = fragment.length()
    } else if fragment.length() != block_width {
      raise AlignEmbossError(
        "Inconsistent EMBOSS sequence block width at line " +
        (index + 1).to_string(),
      )
    }
    let ungapped = align_emboss_remove_gaps(fragment)
    let consumed = ungapped.length()
    let previous_total = ungapped_rows[row_index].length()
    if previous_total == 0 && consumed > 0 {
      if printed_start < printed_end {
        let start_boundary = printed_start - 1
        if printed_start <= 0 || printed_end != start_boundary + consumed {
          raise AlignEmbossError(
            "Invalid forward EMBOSS coordinates at line " +
            (index + 1).to_string(),
          )
        }
        starts[row_index] = start_boundary
        ends[row_index] = printed_end
        directions[row_index] = 1
      } else if printed_start > printed_end {
        let end_boundary = printed_end - 1
        if printed_end <= 0 || end_boundary != printed_start - consumed {
          raise AlignEmbossError(
            "Invalid reverse EMBOSS coordinates at line " +
            (index + 1).to_string(),
          )
        }
        starts[row_index] = printed_start
        ends[row_index] = end_boundary
        directions[row_index] = -1
      } else {
        if consumed != 1 || printed_start <= 0 {
          raise AlignEmbossError(
            "Ambiguous EMBOSS coordinates at line " + (index + 1).to_string(),
          )
        }
        // A one-residue block prints the same coordinate in either direction.
        // Keep the one-based residue coordinate until another block resolves it.
        starts[row_index] = printed_start
        ends[row_index] = printed_end
        directions[row_index] = 0
      }
    } else if consumed == 0 {
      if directions[row_index] >= 0 {
        if printed_start != ends[row_index] || printed_end != printed_start {
          raise AlignEmbossError(
            "Invalid forward gap-only EMBOSS coordinates at line " +
            (index + 1).to_string(),
          )
        }
      } else if printed_start - 1 != ends[row_index] ||
        printed_end != printed_start {
        raise AlignEmbossError(
          "Invalid reverse gap-only EMBOSS coordinates at line " +
          (index + 1).to_string(),
        )
      }
    } else if directions[row_index] > 0 {
      let start_boundary = printed_start - 1
      if printed_start <= 0 ||
        start_boundary != ends[row_index] ||
        printed_end != start_boundary + consumed {
        raise AlignEmbossError(
          "Discontinuous forward EMBOSS coordinates at line " +
          (index + 1).to_string(),
        )
      }
      ends[row_index] = printed_end
    } else if directions[row_index] < 0 {
      let end_boundary = printed_end - 1
      if printed_end <= 0 ||
        printed_start != ends[row_index] ||
        end_boundary != printed_start - consumed {
        raise AlignEmbossError(
          "Discontinuous reverse EMBOSS coordinates at line " +
          (index + 1).to_string(),
        )
      }
      ends[row_index] = end_boundary
    } else if printed_start == ends[row_index] + 1 {
      let start_boundary = printed_start - 1
      if printed_end != start_boundary + consumed {
        raise AlignEmbossError(
          "Discontinuous forward EMBOSS coordinates at line " +
          (index + 1).to_string(),
        )
      }
      starts[row_index] = starts[row_index] - 1
      ends[row_index] = printed_end
      directions[row_index] = 1
    } else if printed_start + 1 == ends[row_index] {
      let end_boundary = printed_end - 1
      if printed_end <= 0 || end_boundary != printed_start - consumed {
        raise AlignEmbossError(
          "Discontinuous reverse EMBOSS coordinates at line " +
          (index + 1).to_string(),
        )
      }
      ends[row_index] = end_boundary
      directions[row_index] = -1
    } else {
      raise AlignEmbossError(
        "Cannot determine EMBOSS row orientation at line " +
        (index + 1).to_string(),
      )
    }
    aligned_rows[row_index] = aligned_rows[row_index] + fragment
    ungapped_rows[row_index] = ungapped_rows[row_index] + ungapped
    row_index = row_index + 1
    index = index + 1
  }
  let sequences : Array[AlignEmbossSequence] = []
  for row = 0; row < number_of_sequences; row = row + 1 {
    if directions[row] == 0 && ungapped_rows[row].length() == 1 {
      starts[row] = starts[row] - 1
      directions[row] = 1
    }
    if aligned_rows[row].length() != length {
      raise AlignEmbossError(
        "EMBOSS row '" +
        identifiers[row] +
        "' has width " +
        aligned_rows[row].length().to_string() +
        "; expected " +
        length.to_string(),
      )
    }
    sequences.push(
      AlignEmbossSequence::create(
        identifiers[row],
        aligned_rows[row],
        starts[row],
        ends[row],
      ),
    )
  }
  let consensus = if consensus_seen {
    let value = consensus_builder.to_string()
    if value.length() != length {
      raise AlignEmbossError(
        "EMBOSS consensus width does not match declared Length",
      )
    }
    value
  } else {
    ""
  }
  let alignment = AlignEmbossAlignment::create(
    sequences,
    annotations,
    consensus~,
  )
  align_emboss_validate_reported_counts(alignment)
  ignore(metadata)
  (alignment, index)
}

///|
fn align_emboss_validate_alignment(
  alignment : AlignEmbossAlignment,
) -> Unit raise AlignEmbossError {
  if alignment.sequences.length() == 0 {
    raise AlignEmbossError("EMBOSS alignment must contain a sequence")
  }
  if alignment.annotations.length <= 0 {
    raise AlignEmbossError("EMBOSS alignment length must be positive")
  }
  for sequence in alignment.sequences {
    align_emboss_validate_id(sequence.id)
    if sequence.aligned_sequence.length() != alignment.annotations.length {
      raise AlignEmbossError(
        "EMBOSS sequence width does not match alignment length",
      )
    }
    align_emboss_validate_aligned_fragment(
      sequence.aligned_sequence,
      "sequence row",
    )
    if align_emboss_remove_gaps(sequence.aligned_sequence) != sequence.sequence {
      raise AlignEmbossError("EMBOSS ungapped sequence is inconsistent")
    }
    let span = if sequence.start >= sequence.end {
      sequence.start - sequence.end
    } else {
      sequence.end - sequence.start
    }
    if span != sequence.sequence.length() {
      raise AlignEmbossError("EMBOSS sequence coordinate span is inconsistent")
    }
    if sequence.start <= sequence.end && sequence.strand != AlignEmbossForward {
      raise AlignEmbossError("EMBOSS sequence strand is inconsistent")
    }
    if sequence.start > sequence.end && sequence.strand != AlignEmbossReverse {
      raise AlignEmbossError("EMBOSS sequence strand is inconsistent")
    }
  }
  if alignment.consensus.length() != 0 &&
    alignment.consensus.length() != alignment.annotations.length {
    raise AlignEmbossError("EMBOSS consensus width must match alignment length")
  }
  for index = 0; index < alignment.consensus.length(); index = index + 1 {
    let code = alignment.consensus.unsafe_get(index).to_int()
    if code < 32 || code > 126 {
      raise AlignEmbossError("EMBOSS consensus contains a control character")
    }
  }
}

///|
fn align_emboss_validate_reported_counts(
  alignment : AlignEmbossAlignment,
) -> Unit raise AlignEmbossError {
  if alignment.sequences.length() != 2 {
    return
  }
  let counts = alignment.pair_counts(0, 1)
  match alignment.annotations.identity {
    Some(value) =>
      if value != counts.identities {
        raise AlignEmbossError(
          "Reported EMBOSS Identity does not match sequence rows",
        )
      }
    None => ()
  }
  match alignment.annotations.gaps {
    Some(value) =>
      if value !=
        counts.insertions + counts.deletions + counts.double_gap_columns {
        raise AlignEmbossError(
          "Reported EMBOSS Gaps does not match sequence rows",
        )
      }
    None => ()
  }
  match alignment.annotations.similarity {
    Some(value) =>
      if alignment.consensus.length() == alignment.alignment_length() &&
        value != counts.positives {
        raise AlignEmbossError(
          "Reported EMBOSS Similarity does not match consensus",
        )
      }
    None => ()
  }
}

///|
fn align_emboss_validate_row(
  alignment : AlignEmbossAlignment,
  row : Int,
) -> Unit raise AlignEmbossError {
  if row < 0 || row >= alignment.sequences.length() {
    raise AlignEmbossError("EMBOSS row index is out of range")
  }
}

///|
fn align_emboss_boundary_at(
  sequence : AlignEmbossSequence,
  column : Int,
) -> Int {
  let mut consumed = 0
  for index = 0; index < column; index = index + 1 {
    if sequence.aligned_sequence.unsafe_get(index).to_int() != '-'.to_int() {
      consumed = consumed + 1
    }
  }
  if sequence.strand == AlignEmbossForward {
    sequence.start + consumed
  } else {
    sequence.start - consumed
  }
}

///|
fn align_emboss_column_step(
  sequence : AlignEmbossSequence,
  column : Int,
) -> Int {
  if sequence.aligned_sequence.unsafe_get(column).to_int() == '-'.to_int() {
    0
  } else if sequence.strand == AlignEmbossForward {
    1
  } else {
    -1
  }
}

///|
fn align_emboss_printed_coordinates(
  sequence : AlignEmbossSequence,
  boundary : Int,
  residues : Int,
) -> (Int, Int) {
  if sequence.strand == AlignEmbossForward {
    if residues == 0 {
      (boundary, boundary)
    } else {
      (boundary + 1, boundary + residues)
    }
  } else if residues == 0 {
    (boundary + 1, boundary + 1)
  } else {
    (boundary, boundary - residues + 1)
  }
}

///|
fn align_emboss_sequence_prefix(id : String, start : Int) -> String {
  let displayed = if id.length() > 13 { id[0:13].to_owned() } else { id }
  align_emboss_pad_right(displayed, 13) +
  " " +
  align_emboss_pad_left(start.to_string(), 6) +
  " "
}

///|
fn align_emboss_pad_left(text : String, width : Int) -> String {
  if text.length() >= width {
    text
  } else {
    " ".repeat(width - text.length()) + text
  }
}

///|
fn align_emboss_pad_right(text : String, width : Int) -> String {
  if text.length() >= width {
    text
  } else {
    text + " ".repeat(width - text.length())
  }
}

///|
fn align_emboss_consensus_chunk(
  raw : String,
  width : Int,
  line_number : Int,
) -> String raise AlignEmbossError {
  if raw.length() <= width {
    raw + " ".repeat(width - raw.length())
  } else {
    let extra = raw[width:].to_owned()
    if extra.trim().length() > 0 {
      raise AlignEmbossError(
        "EMBOSS consensus exceeds block width at line " +
        line_number.to_string(),
      )
    }
    raw[0:width].to_owned()
  }
}

///|
fn align_emboss_split_annotation(
  content : String,
  line_number : Int,
) -> (String, String) raise AlignEmbossError {
  let colon = align_emboss_find_char(content, ':'.to_int())
  if colon >= 0 {
    return (
      content[0:colon].trim().to_owned(),
      content[colon + 1:].trim().to_owned(),
    )
  }
  let marker = " = "
  let equal = align_emboss_find_substring(content, marker)
  if equal >= 0 {
    return (
      content[0:equal].trim().to_owned(),
      content[equal + marker.length():].trim().to_owned(),
    )
  }
  raise AlignEmbossError(
    "Malformed EMBOSS annotation at line " + line_number.to_string(),
  )
}

///|
fn align_emboss_split_key_value(
  content : String,
  context : String,
  line_number : Int,
) -> (String, String) raise AlignEmbossError {
  let separator = align_emboss_find_char(content, ':'.to_int())
  if separator <= 0 {
    raise AlignEmbossError(
      "Malformed EMBOSS " + context + " at line " + line_number.to_string(),
    )
  }
  let key = content[0:separator].trim().to_owned()
  let value = content[separator + 1:].trim().to_owned()
  if key.length() == 0 || value.length() == 0 {
    raise AlignEmbossError(
      "Empty EMBOSS " + context + " field at line " + line_number.to_string(),
    )
  }
  (key, value)
}

///|
fn align_emboss_parse_fraction_count(
  text : String,
  field : String,
  line_number : Int,
) -> Int raise AlignEmbossError {
  let slash = align_emboss_find_char(text, '/'.to_int())
  let count_text = if slash < 0 {
    text
  } else {
    text[0:slash].trim().to_owned()
  }
  align_emboss_parse_nonnegative_int(count_text, field, line_number)
}

///|
fn align_emboss_parse_positive_int(
  text : String,
  field : String,
  line_number : Int,
) -> Int raise AlignEmbossError {
  let value = align_emboss_parse_nonnegative_int(text, field, line_number)
  if value <= 0 {
    raise AlignEmbossError(
      "EMBOSS " + field + " must be positive at line " + line_number.to_string(),
    )
  }
  value
}

///|
fn align_emboss_parse_nonnegative_int(
  text : String,
  field : String,
  line_number : Int,
) -> Int raise AlignEmbossError {
  if text.length() == 0 {
    raise AlignEmbossError(
      "Empty EMBOSS " + field + " at line " + line_number.to_string(),
    )
  }
  let mut value = 0
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise AlignEmbossError(
        "Invalid EMBOSS " + field + " at line " + line_number.to_string(),
      )
    }
    value = value * 10 + code - '0'.to_int()
  }
  value
}

///|
fn align_emboss_parse_number(
  text : String,
  field : String,
  line_number : Int,
) -> Double raise AlignEmbossError {
  if text.length() == 0 {
    raise AlignEmbossError(
      "Empty EMBOSS " + field + " at line " + line_number.to_string(),
    )
  }
  let mut has_digit = false
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code >= '0'.to_int() && code <= '9'.to_int() {
      has_digit = true
    } else if code != '-'.to_int() &&
      code != '+'.to_int() &&
      code != '.'.to_int() &&
      code != 'e'.to_int() &&
      code != 'E'.to_int() {
      raise AlignEmbossError(
        "Invalid EMBOSS " + field + " at line " + line_number.to_string(),
      )
    }
  }
  if !has_digit {
    raise AlignEmbossError(
      "Invalid EMBOSS " + field + " at line " + line_number.to_string(),
    )
  }
  match parse_double(text) {
    Some(value) =>
      if value.abs() > 1.0e300 {
        raise AlignEmbossError(
          "Non-finite EMBOSS " + field + " at line " + line_number.to_string(),
        )
      } else {
        value
      }
    None =>
      raise AlignEmbossError(
        "Invalid EMBOSS " + field + " at line " + line_number.to_string(),
      )
  }
}

///|
fn align_emboss_validate_id(id : String) -> Unit raise AlignEmbossError {
  if id.length() == 0 {
    raise AlignEmbossError("EMBOSS sequence identifier must not be empty")
  }
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    if code <= 32 || code == 127 {
      raise AlignEmbossError(
        "EMBOSS sequence identifier must not contain whitespace",
      )
    }
  }
}

///|
fn align_emboss_validate_header_value(
  value : String,
  field : String,
) -> Unit raise AlignEmbossError {
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    if code == '\n'.to_int() || code == '\r'.to_int() {
      raise AlignEmbossError("EMBOSS " + field + " must not contain a newline")
    }
  }
}

///|
fn align_emboss_validate_aligned_fragment(
  fragment : String,
  context : String,
) -> Unit raise AlignEmbossError {
  if fragment.length() == 0 {
    raise AlignEmbossError("Empty EMBOSS aligned fragment in " + context)
  }
  for index = 0; index < fragment.length(); index = index + 1 {
    let code = fragment.unsafe_get(index).to_int()
    if code <= 32 || code > 126 {
      raise AlignEmbossError("Invalid EMBOSS aligned character in " + context)
    }
  }
}

///|
fn align_emboss_validate_optional_count(
  value : Int?,
  length : Int,
  field : String,
) -> Unit raise AlignEmbossError {
  match value {
    Some(count) =>
      if count < 0 || count > length {
        raise AlignEmbossError(
          "EMBOSS " + field + " must be between zero and Length",
        )
      }
    None => ()
  }
}

///|
fn align_emboss_validate_optional_nonnegative_double(
  value : Double?,
  field : String,
) -> Unit raise AlignEmbossError {
  match value {
    Some(number) =>
      if number < 0.0 || number.abs() > 1.0e300 {
        raise AlignEmbossError(
          "EMBOSS " + field + " must be finite and non-negative",
        )
      }
    None => ()
  }
}

///|
fn align_emboss_remove_gaps(sequence : String) -> String {
  let result = StringBuilder::new(size_hint=sequence.length())
  for index = 0; index < sequence.length(); index = index + 1 {
    let code = sequence.unsafe_get(index).to_int()
    if code != '-'.to_int() {
      result.write_char(code.unsafe_to_char())
    }
  }
  result.to_string()
}

///|
fn align_emboss_count_residues(sequence : String) -> Int {
  let mut count = 0
  for index = 0; index < sequence.length(); index = index + 1 {
    if sequence.unsafe_get(index).to_int() != '-'.to_int() {
      count = count + 1
    }
  }
  count
}

///|
fn align_emboss_lines(text : String) -> Array[String] {
  let raw = text.split("\n")
  let lines : Array[String] = []
  for line in raw {
    let owned = line.to_owned()
    if owned.length() > 0 &&
      owned.unsafe_get(owned.length() - 1).to_int() == '\r'.to_int() {
      lines.push(owned[0:owned.length() - 1].to_owned())
    } else {
      lines.push(owned)
    }
  }
  lines
}

///|
fn align_emboss_split_whitespace(text : String) -> Array[String] {
  let result : Array[String] = []
  let mut index = 0
  while index < text.length() {
    while index < text.length() &&
          align_emboss_is_whitespace(text.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index >= text.length() {
      break
    }
    let start = index
    while index < text.length() &&
          !align_emboss_is_whitespace(text.unsafe_get(index).to_int()) {
      index = index + 1
    }
    result.push(text[start:index].to_owned())
  }
  result
}

///|
fn align_emboss_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\r'.to_int() ||
  code == '\n'.to_int()
}

///|
fn align_emboss_has_prefix(text : String, prefix : String) -> Bool {
  if prefix.length() > text.length() {
    return false
  }
  for index = 0; index < prefix.length(); index = index + 1 {
    if text.unsafe_get(index).to_int() != prefix.unsafe_get(index).to_int() {
      return false
    }
  }
  true
}

///|
fn align_emboss_is_separator(text : String) -> Bool {
  text == "#---------------------------------------"
}

///|
fn align_emboss_find_char(text : String, target : Int) -> Int {
  for index = 0; index < text.length(); index = index + 1 {
    if text.unsafe_get(index).to_int() == target {
      return index
    }
  }
  -1
}

///|
fn align_emboss_find_substring(text : String, target : String) -> Int {
  if target.length() == 0 {
    return 0
  }
  if target.length() > text.length() {
    return -1
  }
  for start = 0; start <= text.length() - target.length(); start = start + 1 {
    let mut equal = true
    for offset = 0; offset < target.length(); offset = offset + 1 {
      if text.unsafe_get(start + offset).to_int() !=
        target.unsafe_get(offset).to_int() {
        equal = false
        break
      }
    }
    if equal {
      return start
    }
  }
  -1
}

///|
fn align_emboss_write_optional_double(
  output : StringBuilder,
  key : String,
  value : Double?,
) -> Unit {
  match value {
    Some(number) =>
      output.write_string("# " + key + ": " + number.to_string() + "\n")
    None => ()
  }
}

///|
fn align_emboss_write_optional_count(
  output : StringBuilder,
  key : String,
  value : Int?,
  length : Int,
) -> Unit {
  match value {
    Some(count) =>
      output.write_string(
        "# " + key + ": " + count.to_string() + "/" + length.to_string() + "\n",
      )
    None => ()
  }
}

///|
fn align_emboss_write_long_annotation(
  output : StringBuilder,
  key : String,
  value : String,
) -> Unit {
  if value.length() > 0 {
    output.write_string("# " + key + " = " + value + "\n")
  }
}
