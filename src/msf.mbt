// Biopython-compatible Bio.Align.msf support.
//
// GCG MSF files contain one multiple sequence alignment with a free-form
// preamble, an MSF/Type/Check line, per-sequence metadata, and interleaved
// aligned rows. Dots, tildes, and hyphens are normalized to hyphen gaps.

///|
/// Error raised for malformed GCG MSF data or invalid alignment operations.
pub suberror MsfError {
  MsfError(String)
}

///|
/// Sequence type declared by the GCG MSF header.
pub(all) enum MsfSequenceType {
  MsfProtein
  MsfNucleotide
} derive(Eq, Debug)

///|
/// File-level metadata from one GCG MSF alignment.
///
/// `declared_length` is retained even when a third-party file reports a value
/// different from the actual alignment width, matching Biopython's warning
/// behavior while making the discrepancy queryable.
pub struct MsfMetadata {
  header : String
  title : String
  preamble : Array[String]
  declared_length : Int
  sequence_type : MsfSequenceType
  date_text : String
  checksum_label : String
  declared_checksum : Int
} derive(Eq, Debug)

///|
/// One row and its descriptor metadata in an MSF alignment.
///
/// `aligned_sequence` uses `-` for every gap. `sequence` is the ungapped
/// residue string. `checksum` is the value declared in the file; zero means
/// that the producer did not provide a checksum.
pub struct MsfSequence {
  id : String
  sequence : String
  aligned_sequence : String
  length : Int
  checksum : Int
  weight : Double
} derive(Eq, Debug)

///|
/// A coordinate-aware GCG MSF multiple sequence alignment.
pub struct MsfAlignment {
  metadata : MsfMetadata
  sequences : Array[MsfSequence]
} derive(Eq, Debug)

///|
/// Pairwise statistics for two rows in an MSF alignment.
pub struct MsfPairCounts {
  columns : Int
  aligned : Int
  identities : Int
  mismatches : Int
  gap_columns : Int
  double_gap_columns : Int
  gap_opens : Int
} derive(Eq, Debug)

///|
/// Construct validated file-level MSF metadata.
pub fn MsfMetadata::create(
  declared_length : Int,
  sequence_type : MsfSequenceType,
  header? : String = "",
  title? : String = "",
  preamble? : Array[String] = [],
  date_text? : String = "",
  checksum_label? : String = "Check:",
  declared_checksum? : Int = 0,
) -> MsfMetadata raise MsfError {
  if declared_length <= 0 {
    raise MsfError("MSF declared alignment length must be positive")
  }
  if declared_checksum < 0 || declared_checksum > 9999 {
    raise MsfError("MSF file checksum must be between 0 and 9999")
  }
  if checksum_label != "Check:" && checksum_label != "CompCheck:" {
    raise MsfError("MSF checksum label must be Check: or CompCheck:")
  }
  let canonical_header = if header.length() == 0 {
    match sequence_type {
      MsfProtein => "!!AA_MULTIPLE_ALIGNMENT 1.0"
      MsfNucleotide => "!!NA_MULTIPLE_ALIGNMENT 1.0"
    }
  } else {
    header.trim().to_owned()
  }
  msf_validate_header(canonical_header, sequence_type)
  msf_validate_plain_text(title, "MSF title")
  msf_validate_plain_text(date_text, "MSF date text")
  let copied_preamble : Array[String] = []
  for line in preamble {
    msf_validate_plain_text(line, "MSF preamble line")
    copied_preamble.push(line)
  }
  MsfMetadata::{
    header: canonical_header,
    title,
    preamble: copied_preamble,
    declared_length,
    sequence_type,
    date_text,
    checksum_label,
    declared_checksum,
  }
}

///|
/// Construct a normalized MSF sequence row.
pub fn MsfSequence::create(
  id : String,
  aligned_sequence : String,
  sequence_type : MsfSequenceType,
  weight? : Double = 1.0,
  checksum? : Int? = None,
) -> MsfSequence raise MsfError {
  msf_validate_id(id)
  if weight != weight || weight.abs() > 1.0e300 || weight < 0.0 {
    raise MsfError("MSF sequence weight must be a finite non-negative value")
  }
  if aligned_sequence.length() == 0 {
    raise MsfError("MSF aligned sequence must not be empty")
  }
  let normalized = msf_normalize_aligned_sequence(
    aligned_sequence, sequence_type,
  )
  let sequence = msf_remove_gaps(normalized)
  if sequence.length() == 0 {
    raise MsfError("MSF sequence row must contain at least one residue")
  }
  let declared_checksum = match checksum {
    Some(value) => {
      if value < 0 || value > 9999 {
        raise MsfError("MSF sequence checksum must be between 0 and 9999")
      }
      value
    }
    None => msf_gcg_checksum(sequence)
  }
  MsfSequence::{
    id,
    sequence,
    aligned_sequence: normalized,
    length: sequence.length(),
    checksum: declared_checksum,
    weight,
  }
}

///|
/// Construct a validated MSF alignment from normalized sequence rows.
pub fn MsfAlignment::create(
  metadata : MsfMetadata,
  sequences : Array[MsfSequence],
) -> MsfAlignment raise MsfError {
  let copied : Array[MsfSequence] = []
  for sequence in sequences {
    copied.push(sequence)
  }
  let alignment = MsfAlignment::{ metadata, sequences: copied }
  msf_validate_alignment(alignment)
  alignment
}

///|
/// Construct an MSF alignment directly from row identifiers and aligned rows.
pub fn msf_from_aligned(
  ids : Array[String],
  aligned_sequences : Array[String],
  sequence_type : MsfSequenceType,
  weights? : Array[Double] = [],
) -> MsfAlignment raise MsfError {
  if ids.length() == 0 {
    raise MsfError("MSF alignment must contain at least one sequence")
  }
  if ids.length() != aligned_sequences.length() {
    raise MsfError("MSF identifiers and aligned rows must have equal lengths")
  }
  if weights.length() != 0 && weights.length() != ids.length() {
    raise MsfError("MSF weights must be empty or match the sequence count")
  }
  let sequences : Array[MsfSequence] = []
  for index = 0; index < ids.length(); index = index + 1 {
    let weight = if weights.length() == 0 { 1.0 } else { weights[index] }
    sequences.push(
      MsfSequence::create(
        ids[index],
        aligned_sequences[index],
        sequence_type,
        weight~,
      ),
    )
  }
  let width = sequences[0].aligned_sequence.length()
  let metadata = MsfMetadata::create(width, sequence_type)
  MsfAlignment::create(metadata, sequences)
}

///|
/// Parse one GCG MSF alignment.
///
/// The accepted leading headers are `!!NA_MULTIPLE_ALIGNMENT`,
/// `!!AA_MULTIPLE_ALIGNMENT`, and `PileUp`. The parser accepts `Check:` and
/// EMBOSS's `CompCheck:` spelling, interleaved rows, numeric coordinate lines,
/// CRLF input, and all three common gap characters. Setting
/// `verify_checksums` checks every non-zero sequence and file checksum.
pub fn msf_parse(
  text : String,
  verify_checksums? : Bool = true,
) -> MsfAlignment raise MsfError {
  let lines = msf_normalize_lines(text)
  if lines.length() == 0 || lines[0].trim().length() == 0 {
    raise MsfError("Empty MSF input")
  }
  let header = lines[0].trim().to_owned()
  let header_words = msf_split_whitespace(header)
  if header_words.length() == 0 || !msf_is_known_header(header_words[0]) {
    raise MsfError(
      "Unknown GCG MSF header; expected !!NA_MULTIPLE_ALIGNMENT, " +
      "!!AA_MULTIPLE_ALIGNMENT, or PileUp",
    )
  }

  let preamble : Array[String] = []
  let mut index = 1
  let mut header_line = ""
  while index < lines.length() {
    let trimmed = lines[index].trim().to_owned()
    if msf_contains_word(trimmed, "MSF:") && msf_ends_with(trimmed, "..") {
      header_line = trimmed
      index = index + 1
      break
    }
    if trimmed.length() > 0 {
      preamble.push(lines[index])
    }
    index = index + 1
  }
  if header_line.length() == 0 {
    raise MsfError("Reached end of file without an MSF/Type/Check header line")
  }
  let parsed_header = msf_parse_alignment_header(header_line)
  msf_validate_header(header, parsed_header.2)

  let ids : Array[String] = []
  let lengths : Array[Int] = []
  let checksums : Array[Int] = []
  let weights : Array[Double] = []
  let mut found_terminator = false
  while index < lines.length() {
    let line = lines[index].trim().to_owned()
    index = index + 1
    if line.length() == 0 {
      continue
    }
    if line == "//" {
      found_terminator = true
      break
    }
    if !msf_starts_with(line, "Name:") {
      raise MsfError(
        "Unexpected MSF header line before // at line " + index.to_string(),
      )
    }
    let descriptor = msf_parse_name_line(line, index)
    if msf_find_string(ids, descriptor.0) >= 0 {
      raise MsfError(
        "Duplicated MSF sequence identifier '" + descriptor.0 + "'",
      )
    }
    ids.push(descriptor.0)
    lengths.push(descriptor.1)
    checksums.push(descriptor.2)
    weights.push(descriptor.3)
  }
  if !found_terminator {
    raise MsfError("End of file while looking for MSF // header terminator")
  }
  if ids.length() == 0 {
    raise MsfError("MSF header contains no Name descriptors")
  }
  if index >= lines.length() {
    raise MsfError("End of file after // line; expected sequence data")
  }
  if lines[index].trim().length() != 0 {
    raise MsfError("Expected a blank line after the MSF // header terminator")
  }
  index = index + 1

  let rows = Array::make(ids.length(), "")
  let residue_counts = Array::make(ids.length(), 0)
  let mut body_rows = 0
  while index < lines.length() {
    let line_number = index + 1
    let line = lines[index].trim().to_owned()
    index = index + 1
    if line.length() == 0 {
      continue
    }
    let words = msf_split_whitespace(line)
    if msf_all_integer_words(words) {
      continue
    }
    let row = msf_find_string(ids, words[0])
    if row < 0 {
      raise MsfError(
        "Unexpected MSF sequence identifier '" +
        words[0] +
        "' at line " +
        line_number.to_string(),
      )
    }
    if words.length() < 2 {
      raise MsfError(
        "MSF sequence row has no residues at line " + line_number.to_string(),
      )
    }
    let fragment_builder = StringBuilder::new()
    for word_index = 1; word_index < words.length(); word_index = word_index + 1 {
      fragment_builder.write_string(words[word_index])
    }
    let fragment = msf_normalize_aligned_sequence(
      fragment_builder.to_string(),
      parsed_header.2,
    )
    let residues = msf_count_residues(fragment)
    if residue_counts[row] + residues > lengths[row] {
      raise MsfError(
        "Received more residues than declared for MSF sequence '" +
        ids[row] +
        "'",
      )
    }
    rows[row] = rows[row] + fragment
    residue_counts[row] = residue_counts[row] + residues
    body_rows = body_rows + 1
  }
  if body_rows == 0 {
    raise MsfError("MSF body contains no sequence data")
  }

  let mut actual_width = 0
  for row = 0; row < ids.length(); row = row + 1 {
    if residue_counts[row] != lengths[row] {
      raise MsfError(
        "MSF sequence '" +
        ids[row] +
        "' contains " +
        residue_counts[row].to_string() +
        " residues but declares " +
        lengths[row].to_string(),
      )
    }
    if rows[row].length() > actual_width {
      actual_width = rows[row].length()
    }
  }
  if actual_width == 0 {
    raise MsfError("MSF alignment has zero columns")
  }

  let sequences : Array[MsfSequence] = []
  for row = 0; row < ids.length(); row = row + 1 {
    let padded = msf_pad_right(rows[row], actual_width, '-')
    sequences.push(
      MsfSequence::create(
        ids[row],
        padded,
        parsed_header.2,
        weight=weights[row],
        checksum=Some(checksums[row]),
      ),
    )
  }
  let metadata = MsfMetadata::create(
    parsed_header.1,
    parsed_header.2,
    header~,
    title=parsed_header.0,
    preamble~,
    date_text=parsed_header.3,
    checksum_label=parsed_header.4,
    declared_checksum=parsed_header.5,
  )
  let alignment = MsfAlignment::create(metadata, sequences)
  if verify_checksums {
    msf_verify_checksums(alignment)
  }
  alignment
}

///|
/// Serialize an alignment as canonical GCG MSF text.
///
/// Rows are emitted in 50-column blocks by default and grouped every ten
/// columns. Sequence and file checksums are always recalculated. Gap output
/// may be `.`, `~`, or `-`.
pub fn msf_write(
  alignment : MsfAlignment,
  block_width? : Int = 50,
  group_width? : Int = 10,
  gap_character? : String = ".",
) -> String raise MsfError {
  msf_validate_alignment(alignment)
  if block_width <= 0 || group_width <= 0 || block_width % group_width != 0 {
    raise MsfError(
      "MSF block width must be positive and divisible by group width",
    )
  }
  if gap_character != "." && gap_character != "~" && gap_character != "-" {
    raise MsfError("MSF output gap character must be '.', '~', or '-'")
  }
  let output = StringBuilder::new()
  output.write_string(alignment.metadata.header)
  output.write_char('\n')
  for line in alignment.metadata.preamble {
    output.write_string(line)
    output.write_char('\n')
  }
  output.write_char('\n')
  if alignment.metadata.title.length() > 0 {
    output.write_string(alignment.metadata.title)
    output.write_char(' ')
  }
  output.write_string("MSF: ")
  output.write_string(alignment.alignment_length().to_string())
  output.write_string(" Type: ")
  output.write_string(alignment.metadata.sequence_type.code())
  if alignment.metadata.date_text.length() > 0 {
    output.write_char(' ')
    output.write_string(alignment.metadata.date_text)
  }
  output.write_char(' ')
  output.write_string(alignment.metadata.checksum_label)
  output.write_char(' ')
  output.write_string(alignment.computed_checksum().to_string())
  output.write_string(" ..\n\n")
  for sequence in alignment.sequences {
    output.write_string(" Name: ")
    output.write_string(sequence.id)
    output.write_string(" Len: ")
    output.write_string(sequence.length.to_string())
    output.write_string(" Check: ")
    output.write_string(msf_gcg_checksum(sequence.sequence).to_string())
    output.write_string(" Weight: ")
    output.write_string(sequence.weight.to_string())
    output.write_char('\n')
  }
  output.write_string("//\n\n")

  let width = alignment.alignment_length()
  let gap_code = gap_character.unsafe_get(0).to_int()
  let mut block_start = 0
  while block_start < width {
    let block_end = msf_min(block_start + block_width, width)
    for sequence in alignment.sequences {
      output.write_string("  ")
      output.write_string(sequence.id)
      output.write_string("  ")
      let mut group_start = block_start
      while group_start < block_end {
        if group_start > block_start {
          output.write_char(' ')
        }
        let group_end = msf_min(group_start + group_width, block_end)
        for column = group_start; column < group_end; column = column + 1 {
          let code = sequence.aligned_sequence.unsafe_get(column).to_int()
          output.write_char(
            (if code == '-'.to_int() { gap_code } else { code }).unsafe_to_char(),
          )
        }
        group_start = group_end
      }
      output.write_char('\n')
    }
    output.write_char('\n')
    block_start = block_end
  }
  output.to_string()
}

///|
/// Return `P` or `N`, as used in an MSF `Type:` field.
pub fn MsfSequenceType::code(self : MsfSequenceType) -> String {
  match self {
    MsfProtein => "P"
    MsfNucleotide => "N"
  }
}

///|
/// Return the number of sequence rows.
pub fn MsfAlignment::num_sequences(self : MsfAlignment) -> Int {
  self.sequences.length()
}

///|
/// Return the alignment width including gaps.
pub fn MsfAlignment::alignment_length(self : MsfAlignment) -> Int {
  if self.sequences.length() == 0 {
    0
  } else {
    self.sequences[0].aligned_sequence.length()
  }
}

///|
/// Return true when the file-level MSF length differs from the actual width.
pub fn MsfAlignment::length_mismatch(self : MsfAlignment) -> Bool {
  self.metadata.declared_length != self.alignment_length()
}

///|
/// Locate a row by its exact identifier.
pub fn MsfAlignment::find_sequence(self : MsfAlignment, id : String) -> Int? {
  let index = msf_find_string(
    self.sequences.map(fn(sequence) { sequence.id }),
    id,
  )
  if index < 0 {
    None
  } else {
    Some(index)
  }
}

///|
/// Return one alignment column, or `None` for an invalid column.
pub fn MsfAlignment::column(self : MsfAlignment, column : Int) -> String? {
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
pub fn MsfAlignment::sequence_position_to_column(
  self : MsfAlignment,
  row : Int,
  position : Int,
) -> Int? {
  if row < 0 || row >= self.sequences.length() || position < 0 {
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
/// Map an alignment column to a zero-based ungapped row position.
///
/// Gap columns and out-of-range coordinates return `None`.
pub fn MsfAlignment::column_to_sequence_position(
  self : MsfAlignment,
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
/// Map one residue through the alignment to another row.
pub fn MsfAlignment::map_position(
  self : MsfAlignment,
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
/// Return per-column residue coordinates for two rows.
pub fn MsfAlignment::aligned_pairs(
  self : MsfAlignment,
  first_row : Int,
  second_row : Int,
) -> Array[(Int?, Int?)] raise MsfError {
  msf_validate_row(self, first_row)
  msf_validate_row(self, second_row)
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
/// Return a compact Biopython-style coordinate path for all sequence rows.
///
/// A boundary is emitted whenever the per-row residue/gap movement vector
/// changes, plus the initial and final boundaries.
pub fn MsfAlignment::coordinate_path(self : MsfAlignment) -> Array[Array[Int]] {
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
      msf_movement_changes(self, column, column + 1)
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
/// Compute pairwise identity, mismatch, and gap statistics.
pub fn MsfAlignment::pair_counts(
  self : MsfAlignment,
  first_row : Int,
  second_row : Int,
) -> MsfPairCounts raise MsfError {
  msf_validate_row(self, first_row)
  msf_validate_row(self, second_row)
  let first = self.sequences[first_row].aligned_sequence
  let second = self.sequences[second_row].aligned_sequence
  let mut aligned = 0
  let mut identities = 0
  let mut mismatches = 0
  let mut gap_columns = 0
  let mut double_gap_columns = 0
  let mut gap_opens = 0
  let mut in_gap = false
  for column = 0; column < self.alignment_length(); column = column + 1 {
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
  MsfPairCounts::{
    columns: self.alignment_length(),
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
pub fn MsfPairCounts::identity(self : MsfPairCounts) -> Double {
  if self.aligned == 0 {
    0.0
  } else {
    self.identities.to_double() / self.aligned.to_double()
  }
}

///|
/// Calculate a majority consensus. Gaps do not vote.
pub fn MsfAlignment::consensus(
  self : MsfAlignment,
  minimum_fraction? : Double = 0.0,
) -> String raise MsfError {
  if minimum_fraction != minimum_fraction ||
    minimum_fraction < 0.0 ||
    minimum_fraction > 1.0 {
    raise MsfError("MSF consensus minimum fraction must be between 0 and 1")
  }
  let output = StringBuilder::new(size_hint=self.alignment_length())
  for column = 0; column < self.alignment_length(); column = column + 1 {
    let counts = Array::make(128, 0)
    let mut residues = 0
    for sequence in self.sequences {
      let code = sequence.aligned_sequence.unsafe_get(column).to_int()
      if code != '-'.to_int() {
        if code >= 0 && code < counts.length() {
          counts[code] = counts[code] + 1
        }
        residues = residues + 1
      }
    }
    if residues == 0 {
      output.write_char('-')
      continue
    }
    let mut best_code = 'A'.to_int()
    let mut best_count = -1
    for code = 0; code < counts.length(); code = code + 1 {
      if counts[code] > best_count {
        best_code = code
        best_count = counts[code]
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
/// Return per-column residue occupancy as a fraction of all rows.
pub fn MsfAlignment::occupancy(self : MsfAlignment) -> Array[Double] {
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
/// Calculate the standard GCG checksum, cycling position weights from 1 to 57.
pub fn msf_gcg_checksum(sequence : String) -> Int {
  let upper = sequence.to_upper()
  let mut position = 0
  let mut checksum = 0
  for index = 0; index < upper.length(); index = index + 1 {
    position = position + 1
    checksum = (checksum + position * upper.unsafe_get(index).to_int()) % 10000
    if position == 57 {
      position = 0
    }
  }
  checksum
}

///|
/// Return the sum-of-row-checksums file checksum.
pub fn MsfAlignment::computed_checksum(self : MsfAlignment) -> Int {
  let mut checksum = 0
  for sequence in self.sequences {
    checksum = (checksum + msf_gcg_checksum(sequence.sequence)) % 10000
  }
  checksum
}

///|
/// Return whether one row's declared checksum is absent or correct.
pub fn MsfSequence::checksum_valid(self : MsfSequence) -> Bool {
  self.checksum == 0 || self.checksum == msf_gcg_checksum(self.sequence)
}

///|
/// Return whether all non-zero row and file checksums are correct.
pub fn MsfAlignment::checksums_valid(self : MsfAlignment) -> Bool {
  for sequence in self.sequences {
    if !sequence.checksum_valid() {
      return false
    }
  }
  self.metadata.declared_checksum == 0 ||
  self.metadata.declared_checksum == self.computed_checksum()
}

///|
/// Return a compact alignment summary.
pub fn MsfAlignment::summary(self : MsfAlignment) -> String {
  "MsfAlignment(type=" +
  self.metadata.sequence_type.code() +
  ", sequences=" +
  self.num_sequences().to_string() +
  ", columns=" +
  self.alignment_length().to_string() +
  ", declared=" +
  self.metadata.declared_length.to_string() +
  ", checksum=" +
  self.computed_checksum().to_string() +
  ")"
}

///|
/// Return a small interleaved nucleotide MSF example.
pub fn msf_example_text() -> String {
  "!!NA_MULTIPLE_ALIGNMENT 1.0\n" +
  "\n" +
  "MoonBit example MSF: 12 Type: N Check: 0 ..\n" +
  "\n" +
  " Name: reference Len: 12 Check: 5688 Weight: 1.0\n" +
  " Name: query_one Len: 10 Check: 3984 Weight: 0.75\n" +
  " Name: query_two Len: 11 Check: 4901 Weight: 1.25\n" +
  "//\n" +
  "\n" +
  " reference ACGTAC\n" +
  " query_one ACG.AC\n" +
  " query_two A.GTTC\n" +
  "\n" +
  " 7 12\n" +
  " reference GTACGT\n" +
  " query_one GTAC.T\n" +
  " query_two GTACGT\n"
}

///|
fn msf_validate_alignment(alignment : MsfAlignment) -> Unit raise MsfError {
  if alignment.sequences.length() == 0 {
    raise MsfError("MSF alignment must contain at least one sequence")
  }
  let width = alignment.sequences[0].aligned_sequence.length()
  if width == 0 {
    raise MsfError("MSF alignment must contain at least one column")
  }
  let ids : Array[String] = []
  for sequence in alignment.sequences {
    msf_validate_id(sequence.id)
    if msf_find_string(ids, sequence.id) >= 0 {
      raise MsfError("Duplicated MSF sequence identifier '" + sequence.id + "'")
    }
    ids.push(sequence.id)
    if sequence.aligned_sequence.length() != width {
      raise MsfError("MSF aligned rows must have equal widths")
    }
    let normalized = msf_normalize_aligned_sequence(
      sequence.aligned_sequence,
      alignment.metadata.sequence_type,
    )
    if normalized != sequence.aligned_sequence ||
      msf_remove_gaps(normalized) != sequence.sequence ||
      sequence.sequence.length() != sequence.length {
      raise MsfError("MSF sequence row contains inconsistent derived fields")
    }
    if sequence.weight != sequence.weight ||
      sequence.weight.abs() > 1.0e300 ||
      sequence.weight < 0.0 {
      raise MsfError("MSF sequence weight must be finite and non-negative")
    }
    if sequence.checksum < 0 || sequence.checksum > 9999 {
      raise MsfError("MSF sequence checksum must be between 0 and 9999")
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
      raise MsfError(
        "MSF alignment contains an all-gap column at " + column.to_string(),
      )
    }
  }
}

///|
fn msf_verify_checksums(alignment : MsfAlignment) -> Unit raise MsfError {
  for sequence in alignment.sequences {
    if !sequence.checksum_valid() {
      raise MsfError(
        "MSF checksum mismatch for sequence '" +
        sequence.id +
        "': declared " +
        sequence.checksum.to_string() +
        ", calculated " +
        msf_gcg_checksum(sequence.sequence).to_string(),
      )
    }
  }
  if alignment.metadata.declared_checksum != 0 &&
    alignment.metadata.declared_checksum != alignment.computed_checksum() {
    raise MsfError(
      "MSF file checksum mismatch: declared " +
      alignment.metadata.declared_checksum.to_string() +
      ", calculated " +
      alignment.computed_checksum().to_string(),
    )
  }
}

///|
fn msf_validate_row(alignment : MsfAlignment, row : Int) -> Unit raise MsfError {
  if row < 0 || row >= alignment.sequences.length() {
    raise MsfError("MSF row index is out of bounds")
  }
}

///|
fn msf_validate_header(
  header : String,
  sequence_type : MsfSequenceType,
) -> Unit raise MsfError {
  let words = msf_split_whitespace(header)
  if words.length() == 0 || !msf_is_known_header(words[0]) {
    raise MsfError("Unknown GCG MSF file header")
  }
  if words[0] == "!!AA_MULTIPLE_ALIGNMENT" && sequence_type != MsfProtein {
    raise MsfError("Amino-acid MSF header conflicts with Type: N")
  }
  if words[0] == "!!NA_MULTIPLE_ALIGNMENT" && sequence_type != MsfNucleotide {
    raise MsfError("Nucleotide MSF header conflicts with Type: P")
  }
}

///|
fn msf_validate_id(id : String) -> Unit raise MsfError {
  if id.length() == 0 {
    raise MsfError("MSF sequence identifier must not be empty")
  }
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    if msf_is_whitespace(code) || code < 33 || code > 126 {
      raise MsfError(
        "MSF sequence identifier contains whitespace or a control character",
      )
    }
  }
}

///|
fn msf_validate_plain_text(
  value : String,
  label : String,
) -> Unit raise MsfError {
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    if code == '\n'.to_int() || code == '\r'.to_int() || code < 9 {
      raise MsfError(label + " contains a line break or control character")
    }
  }
}

///|
fn msf_parse_alignment_header(
  line : String,
) -> (String, Int, MsfSequenceType, String, String, Int) raise MsfError {
  let words = msf_split_whitespace(line)
  let msf_index = msf_find_string(words, "MSF:")
  if msf_index < 0 ||
    msf_index + 3 >= words.length() ||
    words[msf_index + 2] != "Type:" ||
    words[words.length() - 1] != ".." {
    raise MsfError(
      "Malformed MSF header; expected MSF: <length> Type: <P|N> ... Check: <value> ..",
    )
  }
  let check_index = msf_find_last_checksum_label(words)
  if check_index < msf_index + 4 || check_index + 2 != words.length() - 1 {
    raise MsfError("Malformed MSF Check:/CompCheck: header fields")
  }
  let length = msf_parse_nonnegative_int(words[msf_index + 1], "MSF length")
  if length == 0 {
    raise MsfError("MSF declared alignment length must be positive")
  }
  let sequence_type = match words[msf_index + 3] {
    "P" => MsfProtein
    "N" => MsfNucleotide
    value => raise MsfError("MSF Type must be P or N, not '" + value + "'")
  }
  let checksum = msf_parse_checksum(words[check_index + 1], "MSF file checksum")
  let title = msf_join_words(words, 0, msf_index)
  let date_text = msf_join_words(words, msf_index + 4, check_index)
  (title, length, sequence_type, date_text, words[check_index], checksum)
}

///|
fn msf_parse_name_line(
  line : String,
  line_number : Int,
) -> (String, Int, Int, Double) raise MsfError {
  let words = msf_split_whitespace(line)
  let name_index = msf_find_string(words, "Name:")
  let length_index = msf_find_string(words, "Len:")
  let checksum_index = msf_find_string(words, "Check:")
  let weight_index = msf_find_string(words, "Weight:")
  if name_index != 0 ||
    length_index < 2 ||
    checksum_index < length_index + 2 ||
    weight_index < checksum_index + 2 ||
    weight_index + 1 != words.length() - 1 {
    raise MsfError(
      "Malformed GCG MSF Name line at line " + line_number.to_string(),
    )
  }
  if length_index == 3 && words[2] != "oo" {
    raise MsfError(
      "Only optional 'oo' is allowed between MSF Name and Len fields",
    )
  }
  if length_index > 3 {
    raise MsfError("Unexpected tokens before MSF Len field")
  }
  let id = words[1]
  msf_validate_id(id)
  let length = msf_parse_nonnegative_int(
    words[length_index + 1],
    "MSF sequence length",
  )
  if length == 0 {
    raise MsfError("MSF sequence length must be positive")
  }
  let checksum = msf_parse_checksum(
    words[checksum_index + 1],
    "MSF sequence checksum",
  )
  let weight = match parse_double(words[weight_index + 1]) {
    Some(value) =>
      if value != value || value.abs() > 1.0e300 || value < 0.0 {
        raise MsfError("MSF sequence weight must be finite and non-negative")
      } else {
        value
      }
    None => raise MsfError("Invalid MSF sequence weight")
  }
  (id, length, checksum, weight)
}

///|
fn msf_parse_nonnegative_int(
  text : String,
  label : String,
) -> Int raise MsfError {
  if text.length() == 0 {
    raise MsfError(label + " is empty")
  }
  let mut value = 0
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise MsfError(label + " must be a non-negative integer")
    }
    let digit = code - '0'.to_int()
    if value > 214748364 || (value == 214748364 && digit > 7) {
      raise MsfError(label + " is outside the supported integer range")
    }
    value = value * 10 + digit
  }
  value
}

///|
fn msf_parse_checksum(text : String, label : String) -> Int raise MsfError {
  let value = msf_parse_nonnegative_int(text, label)
  if value > 9999 {
    raise MsfError(label + " must be between 0 and 9999")
  }
  value
}

///|
fn msf_normalize_aligned_sequence(
  sequence : String,
  sequence_type : MsfSequenceType,
) -> String raise MsfError {
  let output = StringBuilder::new(size_hint=sequence.length())
  for index = 0; index < sequence.length(); index = index + 1 {
    let original = sequence.unsafe_get(index).to_int()
    let code = msf_upper_code(original)
    if code == '-'.to_int() || code == '.'.to_int() || code == '~'.to_int() {
      output.write_char('-')
    } else if msf_valid_residue(code, sequence_type) {
      output.write_char(code.unsafe_to_char())
    } else {
      raise MsfError(
        "Invalid " +
        (match sequence_type {
          MsfProtein => "protein"
          MsfNucleotide => "nucleotide"
        }) +
        " residue in MSF aligned sequence",
      )
    }
  }
  output.to_string()
}

///|
fn msf_valid_residue(code : Int, sequence_type : MsfSequenceType) -> Bool {
  match sequence_type {
    MsfProtein =>
      (code >= 'A'.to_int() && code <= 'Z'.to_int()) ||
      code == '*'.to_int() ||
      code == '?'.to_int()
    MsfNucleotide => msf_contains_code("ACGTURYSWKMBDHVNX", code)
  }
}

///|
fn msf_upper_code(code : Int) -> Int {
  if code >= 'a'.to_int() && code <= 'z'.to_int() {
    code - 'a'.to_int() + 'A'.to_int()
  } else {
    code
  }
}

///|
fn msf_remove_gaps(sequence : String) -> String {
  let output = StringBuilder::new(size_hint=sequence.length())
  for index = 0; index < sequence.length(); index = index + 1 {
    let code = sequence.unsafe_get(index).to_int()
    if code != '-'.to_int() {
      output.write_char(code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn msf_count_residues(sequence : String) -> Int {
  let mut count = 0
  for index = 0; index < sequence.length(); index = index + 1 {
    if sequence.unsafe_get(index).to_int() != '-'.to_int() {
      count = count + 1
    }
  }
  count
}

///|
fn msf_movement_changes(
  alignment : MsfAlignment,
  first_column : Int,
  second_column : Int,
) -> Bool {
  for row = 0; row < alignment.sequences.length(); row = row + 1 {
    let first = alignment.sequences[row].aligned_sequence
      .unsafe_get(first_column)
      .to_int() !=
      '-'.to_int()
    let second = alignment.sequences[row].aligned_sequence
      .unsafe_get(second_column)
      .to_int() !=
      '-'.to_int()
    if first != second {
      return true
    }
  }
  false
}

///|
fn msf_normalize_lines(text : String) -> Array[String] {
  let raw_lines = text.split("\n")
  let lines : Array[String] = []
  for raw in raw_lines {
    let owned = raw.to_owned()
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
fn msf_split_whitespace(text : String) -> Array[String] {
  let words : Array[String] = []
  let mut index = 0
  while index < text.length() {
    while index < text.length() &&
          msf_is_whitespace(text.unsafe_get(index).to_int()) {
      index = index + 1
    }
    let start = index
    while index < text.length() &&
          !msf_is_whitespace(text.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if start < index {
      words.push(text[start:index].to_owned())
    }
  }
  words
}

///|
fn msf_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\n'.to_int() ||
  code == '\r'.to_int()
}

///|
fn msf_is_known_header(value : String) -> Bool {
  value == "!!NA_MULTIPLE_ALIGNMENT" ||
  value == "!!AA_MULTIPLE_ALIGNMENT" ||
  value == "PileUp"
}

///|
fn msf_all_integer_words(words : Array[String]) -> Bool {
  if words.length() == 0 {
    return false
  }
  for word in words {
    if word.length() == 0 {
      return false
    }
    for index = 0; index < word.length(); index = index + 1 {
      let code = word.unsafe_get(index).to_int()
      if code < '0'.to_int() || code > '9'.to_int() {
        return false
      }
    }
  }
  true
}

///|
fn msf_contains_word(text : String, word : String) -> Bool {
  msf_find_string(msf_split_whitespace(text), word) >= 0
}

///|
fn msf_find_string(values : Array[String], target : String) -> Int {
  for index = 0; index < values.length(); index = index + 1 {
    if values[index] == target {
      return index
    }
  }
  -1
}

///|
fn msf_find_last_checksum_label(words : Array[String]) -> Int {
  let mut result = -1
  for index = 0; index < words.length(); index = index + 1 {
    if words[index] == "Check:" || words[index] == "CompCheck:" {
      result = index
    }
  }
  result
}

///|
fn msf_join_words(words : Array[String], start : Int, end : Int) -> String {
  let output = StringBuilder::new()
  for index = start; index < end; index = index + 1 {
    if index > start {
      output.write_char(' ')
    }
    output.write_string(words[index])
  }
  output.to_string()
}

///|
fn msf_pad_right(value : String, width : Int, character : Char) -> String {
  if value.length() >= width {
    return value
  }
  let output = StringBuilder::new(size_hint=width)
  output.write_string(value)
  for _index = value.length(); _index < width; _index = _index + 1 {
    output.write_char(character)
  }
  output.to_string()
}

///|
fn msf_contains_code(values : String, target : Int) -> Bool {
  for index = 0; index < values.length(); index = index + 1 {
    if values.unsafe_get(index).to_int() == target {
      return true
    }
  }
  false
}

///|
fn msf_starts_with(value : String, prefix : String) -> Bool {
  value.length() >= prefix.length() &&
  value[0:prefix.length()].to_owned() == prefix
}

///|
fn msf_ends_with(value : String, suffix : String) -> Bool {
  value.length() >= suffix.length() &&
  value[value.length() - suffix.length():].to_owned() == suffix
}

///|
fn msf_min(first : Int, second : Int) -> Int {
  if first < second {
    first
  } else {
    second
  }
}
