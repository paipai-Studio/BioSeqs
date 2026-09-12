// Biopython-compatible A2M multiple sequence alignment support.
//
// A2M uses upper-case residues and '-' in model match columns, and lower-case
// residues and '.' in insertion columns. The first row defines the state of
// every column. Internally rows are normalized to upper-case residues and '-';
// the separate state string preserves the information required for round-trip
// serialization.

///|
/// Error raised for malformed A2M data or invalid alignment operations.
pub suberror A2mError {
  A2mError(String)
}

///|
/// State assigned to an A2M alignment column.
pub(all) enum A2mColumnState {
  A2mMatch
  A2mInsertion
} derive(Eq, Debug)

///|
/// One normalized sequence row in an A2M alignment.
///
/// `aligned_sequence` contains upper-case residues and '-' gaps. `sequence`
/// contains the ungapped residues.
pub struct A2mSequence {
  id : String
  description : String
  sequence : String
  aligned_sequence : String
} derive(Eq, Debug)

///|
/// A state-aware A2M multiple sequence alignment.
///
/// `states` contains one `D` (model match/deletion) or `I` (insertion) for
/// each alignment column, following Biopython's column annotation convention.
pub struct A2mAlignment {
  sequences : Array[A2mSequence]
  states : String
} derive(Eq, Debug)

///|
/// A contiguous run of insertion-state columns.
///
/// `slot` is the number of match columns preceding the run. Therefore slot 0
/// is before the first match column and slot `match_columns()` is after the
/// final match column.
pub struct A2mInsertionRun {
  slot : Int
  start_column : Int
  end_column : Int
  width : Int
} derive(Eq, Debug)

///|
/// Pairwise statistics for two rows in an A2M alignment.
pub struct A2mPairCounts {
  columns : Int
  aligned : Int
  identities : Int
  mismatches : Int
  gap_columns : Int
  double_gap_columns : Int
  match_aligned : Int
  insertion_aligned : Int
} derive(Eq, Debug)

///|
/// Construct a normalized A2M sequence row.
pub fn A2mSequence::create(
  id : String,
  aligned_sequence : String,
  description? : String = "",
) -> A2mSequence raise A2mError {
  a2m_validate_id(id)
  a2m_validate_description(description)
  if aligned_sequence.length() == 0 {
    raise A2mError("A2M aligned sequence must not be empty")
  }
  let normalized = a2m_normalize_aligned_sequence(aligned_sequence)
  A2mSequence::{
    id,
    description,
    sequence: a2m_remove_gaps(normalized),
    aligned_sequence: normalized,
  }
}

///|
/// Construct an A2M alignment from normalized rows and a `D`/`I` state string.
pub fn A2mAlignment::create(
  sequences : Array[A2mSequence],
  states : String,
) -> A2mAlignment raise A2mError {
  let copied : Array[A2mSequence] = []
  for sequence in sequences {
    copied.push(sequence)
  }
  let alignment = A2mAlignment::{ sequences: copied, states, }
  a2m_validate_alignment(alignment)
  alignment
}

///|
/// Construct an A2M alignment directly from row metadata and normalized rows.
pub fn a2m_from_aligned(
  ids : Array[String],
  aligned_sequences : Array[String],
  states : String,
  descriptions? : Array[String] = [],
) -> A2mAlignment raise A2mError {
  if ids.length() == 0 {
    raise A2mError("A2M alignment must contain at least one sequence")
  }
  if ids.length() != aligned_sequences.length() {
    raise A2mError(
      "A2M identifiers and aligned sequences must have equal lengths",
    )
  }
  if descriptions.length() != 0 && descriptions.length() != ids.length() {
    raise A2mError("A2M descriptions must be empty or match the sequence count")
  }
  let sequences : Array[A2mSequence] = []
  for index = 0; index < ids.length(); index = index + 1 {
    let description = if descriptions.length() == 0 {
      ""
    } else {
      descriptions[index]
    }
    sequences.push(
      A2mSequence::create(ids[index], aligned_sequences[index], description~),
    )
  }
  A2mAlignment::create(sequences, states)
}

///|
/// Parse one A2M multiple sequence alignment.
///
/// The first row defines the state of each column. All following rows must use
/// the same character class in each column. Wrapped sequence lines and CRLF
/// input are accepted; blank lines are ignored.
pub fn a2m_parse(text : String) -> A2mAlignment raise A2mError {
  let lines = a2m_normalize_lines(text)
  let ids : Array[String] = []
  let descriptions : Array[String] = []
  let encoded_rows : Array[String] = []
  let mut current = -1
  for raw_line in lines {
    let line = raw_line.trim().to_owned()
    if line.length() == 0 {
      continue
    }
    if line.unsafe_get(0).to_int() == '>'.to_int() {
      let header = line[1:].trim().to_owned()
      if header.length() == 0 {
        raise A2mError("A2M header must contain a sequence identifier")
      }
      let (id, description) = a2m_parse_header(header)
      a2m_validate_id(id)
      a2m_validate_description(description)
      ids.push(id)
      descriptions.push(description)
      encoded_rows.push("")
      current = encoded_rows.length() - 1
    } else {
      if current < 0 {
        raise A2mError("A2M sequence data appears before the first header")
      }
      a2m_validate_encoded_fragment(line)
      encoded_rows[current] = encoded_rows[current] + line
    }
  }
  if encoded_rows.length() == 0 {
    raise A2mError("Empty A2M input")
  }
  let width = encoded_rows[0].length()
  if width == 0 {
    raise A2mError("A2M sequences must not be empty")
  }
  let states_builder = StringBuilder::new(size_hint=width)
  for column = 0; column < width; column = column + 1 {
    let code = encoded_rows[0].unsafe_get(column).to_int()
    if a2m_is_upper(code) || code == '-'.to_int() {
      states_builder.write_char('D')
    } else if a2m_is_lower(code) || code == '.'.to_int() {
      states_builder.write_char('I')
    } else {
      raise A2mError(
        "Invalid A2M character in first sequence at column " +
        column.to_string(),
      )
    }
  }
  let states = states_builder.to_string()
  let sequences : Array[A2mSequence] = []
  for row = 0; row < encoded_rows.length(); row = row + 1 {
    let encoded = encoded_rows[row]
    if encoded.length() != width {
      raise A2mError(
        "A2M row '" +
        ids[row] +
        "' has width " +
        encoded.length().to_string() +
        "; expected " +
        width.to_string(),
      )
    }
    let normalized = StringBuilder::new(size_hint=width)
    for column = 0; column < width; column = column + 1 {
      let code = encoded.unsafe_get(column).to_int()
      let state = states.unsafe_get(column).to_int()
      if state == 'D'.to_int() {
        if code == '-'.to_int() {
          normalized.write_char('-')
        } else if a2m_is_upper(code) {
          normalized.write_char(code.unsafe_to_char())
        } else {
          raise A2mError(
            "A2M row '" +
            ids[row] +
            "' uses an insertion character in match column " +
            column.to_string(),
          )
        }
      } else if code == '.'.to_int() {
        normalized.write_char('-')
      } else if a2m_is_lower(code) {
        normalized.write_char(a2m_upper_code(code).unsafe_to_char())
      } else {
        raise A2mError(
          "A2M row '" +
          ids[row] +
          "' uses a match character in insertion column " +
          column.to_string(),
        )
      }
    }
    sequences.push(A2mSequence::{
      id: ids[row],
      description: descriptions[row],
      sequence: a2m_remove_gaps(normalized.to_string()),
      aligned_sequence: normalized.to_string(),
    })
  }
  A2mAlignment::create(sequences, states)
}

///|
/// Serialize one A2M alignment.
///
/// A `line_width` of zero emits one sequence line per row. Positive values
/// wrap rows without changing column states.
pub fn a2m_write(
  alignment : A2mAlignment,
  line_width? : Int = 0,
) -> String raise A2mError {
  a2m_validate_alignment(alignment)
  if line_width < 0 {
    raise A2mError("A2M line width must be non-negative")
  }
  let output = StringBuilder::new()
  for row = 0; row < alignment.sequences.length(); row = row + 1 {
    let sequence = alignment.sequences[row]
    output.write_char('>')
    output.write_string(sequence.id)
    if sequence.description.length() > 0 {
      output.write_char(' ')
      output.write_string(sequence.description)
    }
    output.write_char('\n')
    let encoded = alignment.encoded_sequence(row).unwrap()
    if line_width == 0 {
      output.write_string(encoded)
      output.write_char('\n')
    } else {
      let mut start = 0
      while start < encoded.length() {
        let end = if start + line_width < encoded.length() {
          start + line_width
        } else {
          encoded.length()
        }
        output.write_string(encoded[start:end].to_owned())
        output.write_char('\n')
        start = end
      }
    }
  }
  output.to_string()
}

///|
/// Return the number of rows.
pub fn A2mAlignment::num_sequences(self : A2mAlignment) -> Int {
  self.sequences.length()
}

///|
/// Return the number of alignment columns.
pub fn A2mAlignment::alignment_length(self : A2mAlignment) -> Int {
  self.states.length()
}

///|
/// Return the number of model match/deletion columns.
pub fn A2mAlignment::match_columns(self : A2mAlignment) -> Int {
  let mut count = 0
  for index = 0; index < self.states.length(); index = index + 1 {
    if self.states.unsafe_get(index).to_int() == 'D'.to_int() {
      count = count + 1
    }
  }
  count
}

///|
/// Return the number of insertion-state columns.
pub fn A2mAlignment::insertion_columns(self : A2mAlignment) -> Int {
  self.alignment_length() - self.match_columns()
}

///|
/// Return a column state by zero-based alignment coordinate.
pub fn A2mAlignment::column_state(
  self : A2mAlignment,
  column : Int,
) -> A2mColumnState? {
  if column < 0 || column >= self.states.length() {
    None
  } else if self.states.unsafe_get(column).to_int() == 'D'.to_int() {
    Some(A2mMatch)
  } else {
    Some(A2mInsertion)
  }
}

///|
/// Return all row characters at an alignment column.
pub fn A2mAlignment::column(self : A2mAlignment, column : Int) -> String? {
  if column < 0 || column >= self.states.length() {
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
/// Find the first row with an exact sequence identifier.
pub fn A2mAlignment::find_sequence(self : A2mAlignment, id : String) -> Int? {
  for index = 0; index < self.sequences.length(); index = index + 1 {
    if self.sequences[index].id == id {
      return Some(index)
    }
  }
  None
}

///|
/// Return a row in canonical A2M character form.
pub fn A2mAlignment::encoded_sequence(
  self : A2mAlignment,
  row : Int,
) -> String? {
  if row < 0 || row >= self.sequences.length() {
    return None
  }
  let aligned = self.sequences[row].aligned_sequence
  let result = StringBuilder::new(size_hint=self.states.length())
  for column = 0; column < self.states.length(); column = column + 1 {
    let code = aligned.unsafe_get(column).to_int()
    if self.states.unsafe_get(column).to_int() == 'D'.to_int() {
      result.write_char(code.unsafe_to_char())
    } else if code == '-'.to_int() {
      result.write_char('.')
    } else {
      result.write_char(a2m_lower_code(code).unsafe_to_char())
    }
  }
  Some(result.to_string())
}

///|
/// Map an ungapped row coordinate to its alignment column.
pub fn A2mAlignment::sequence_position_to_column(
  self : A2mAlignment,
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
/// Map an alignment column to an ungapped row coordinate.
///
/// Gap columns return `None`.
pub fn A2mAlignment::column_to_sequence_position(
  self : A2mAlignment,
  row : Int,
  column : Int,
) -> Int? {
  if row < 0 ||
    row >= self.sequences.length() ||
    column < 0 ||
    column >= self.states.length() {
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
/// Map a residue coordinate from one row through the alignment to another.
///
/// Returns `None` for out-of-range coordinates and target gaps.
pub fn A2mAlignment::map_position(
  self : A2mAlignment,
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
/// Return per-column coordinates for two rows.
///
/// Residues are represented by zero-based coordinates and gaps by `None`.
pub fn A2mAlignment::aligned_pairs(
  self : A2mAlignment,
  first_row : Int,
  second_row : Int,
) -> Array[(Int?, Int?)] raise A2mError {
  a2m_validate_row(self, first_row)
  a2m_validate_row(self, second_row)
  let first = self.sequences[first_row].aligned_sequence
  let second = self.sequences[second_row].aligned_sequence
  let result : Array[(Int?, Int?)] = []
  let mut first_coordinate = 0
  let mut second_coordinate = 0
  for column = 0; column < self.states.length(); column = column + 1 {
    let first_gap = first.unsafe_get(column).to_int() == '-'.to_int()
    let second_gap = second.unsafe_get(column).to_int() == '-'.to_int()
    let first_value : Int? = if first_gap {
      None
    } else {
      Some(first_coordinate)
    }
    let second_value : Int? = if second_gap {
      None
    } else {
      Some(second_coordinate)
    }
    result.push((first_value, second_value))
    if !first_gap {
      first_coordinate = first_coordinate + 1
    }
    if !second_gap {
      second_coordinate = second_coordinate + 1
    }
  }
  result
}

///|
/// Compute pairwise identity, mismatch, and gap counts for two rows.
pub fn A2mAlignment::pair_counts(
  self : A2mAlignment,
  first_row : Int,
  second_row : Int,
) -> A2mPairCounts raise A2mError {
  a2m_validate_row(self, first_row)
  a2m_validate_row(self, second_row)
  let first = self.sequences[first_row].aligned_sequence
  let second = self.sequences[second_row].aligned_sequence
  let mut aligned = 0
  let mut identities = 0
  let mut mismatches = 0
  let mut gap_columns = 0
  let mut double_gap_columns = 0
  let mut match_aligned = 0
  let mut insertion_aligned = 0
  for column = 0; column < self.states.length(); column = column + 1 {
    let first_code = first.unsafe_get(column).to_int()
    let second_code = second.unsafe_get(column).to_int()
    let first_gap = first_code == '-'.to_int()
    let second_gap = second_code == '-'.to_int()
    if first_gap && second_gap {
      double_gap_columns = double_gap_columns + 1
    } else if first_gap || second_gap {
      gap_columns = gap_columns + 1
    } else {
      aligned = aligned + 1
      if first_code == second_code {
        identities = identities + 1
      } else {
        mismatches = mismatches + 1
      }
      if self.states.unsafe_get(column).to_int() == 'D'.to_int() {
        match_aligned = match_aligned + 1
      } else {
        insertion_aligned = insertion_aligned + 1
      }
    }
  }
  A2mPairCounts::{
    columns: self.states.length(),
    aligned,
    identities,
    mismatches,
    gap_columns,
    double_gap_columns,
    match_aligned,
    insertion_aligned,
  }
}

///|
/// Return identity fraction over columns containing residues in both rows.
pub fn A2mPairCounts::identity(self : A2mPairCounts) -> Double {
  if self.aligned == 0 {
    0.0
  } else {
    self.identities.to_double() / self.aligned.to_double()
  }
}

///|
/// Locate contiguous insertion-state runs and their model slots.
pub fn A2mAlignment::insertion_runs(
  self : A2mAlignment,
) -> Array[A2mInsertionRun] {
  let result : Array[A2mInsertionRun] = []
  let mut column = 0
  let mut slot = 0
  while column < self.states.length() {
    if self.states.unsafe_get(column).to_int() == 'D'.to_int() {
      slot = slot + 1
      column = column + 1
    } else {
      let start = column
      while column < self.states.length() &&
            self.states.unsafe_get(column).to_int() == 'I'.to_int() {
        column = column + 1
      }
      result.push(A2mInsertionRun::{
        slot,
        start_column: start,
        end_column: column,
        width: column - start,
      })
    }
  }
  result
}

///|
/// Calculate a majority consensus.
///
/// Gaps do not vote. Columns without residues produce `-`; columns below
/// `minimum_fraction` produce `X`. Setting `include_insertions` to false
/// returns a model match-state consensus.
pub fn A2mAlignment::consensus(
  self : A2mAlignment,
  include_insertions? : Bool = true,
  minimum_fraction? : Double = 0.0,
) -> String raise A2mError {
  if minimum_fraction != minimum_fraction ||
    minimum_fraction < 0.0 ||
    minimum_fraction > 1.0 {
    raise A2mError("A2M consensus minimum fraction must be between 0 and 1")
  }
  let result = StringBuilder::new()
  for column = 0; column < self.states.length(); column = column + 1 {
    if !include_insertions &&
      self.states.unsafe_get(column).to_int() == 'I'.to_int() {
      continue
    }
    let counts = Array::make(26, 0)
    let mut residues = 0
    for sequence in self.sequences {
      let code = sequence.aligned_sequence.unsafe_get(column).to_int()
      if a2m_is_upper(code) {
        let index = code - 'A'.to_int()
        counts[index] = counts[index] + 1
        residues = residues + 1
      }
    }
    if residues == 0 {
      result.write_char('-')
      continue
    }
    let mut best_index = 0
    let mut best_count = counts[0]
    for index = 1; index < counts.length(); index = index + 1 {
      if counts[index] > best_count {
        best_index = index
        best_count = counts[index]
      }
    }
    let fraction = best_count.to_double() / residues.to_double()
    if fraction < minimum_fraction {
      result.write_char('X')
    } else {
      result.write_char((best_index + 'A'.to_int()).unsafe_to_char())
    }
  }
  result.to_string()
}

///|
/// Return per-column residue occupancy as a fraction of row count.
pub fn A2mAlignment::occupancy(self : A2mAlignment) -> Array[Double] {
  let result : Array[Double] = []
  for column = 0; column < self.states.length(); column = column + 1 {
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
/// Return a new alignment containing only model match/deletion columns.
pub fn A2mAlignment::match_projection(
  self : A2mAlignment,
) -> A2mAlignment raise A2mError {
  let match_count = self.match_columns()
  if match_count == 0 {
    raise A2mError("A2M alignment has no match columns to project")
  }
  let states = StringBuilder::new(size_hint=match_count)
  let mut emitted = 0
  while emitted < match_count {
    states.write_char('D')
    emitted = emitted + 1
  }
  let sequences : Array[A2mSequence] = []
  for sequence in self.sequences {
    let row = StringBuilder::new(size_hint=match_count)
    for column = 0; column < self.states.length(); column = column + 1 {
      if self.states.unsafe_get(column).to_int() == 'D'.to_int() {
        row.write_char(
          sequence.aligned_sequence.unsafe_get(column).unsafe_to_char(),
        )
      }
    }
    sequences.push(
      A2mSequence::create(
        sequence.id,
        row.to_string(),
        description=sequence.description,
      ),
    )
  }
  A2mAlignment::create(sequences, states.to_string())
}

///|
/// Extract a non-empty half-open range of alignment columns.
pub fn A2mAlignment::slice_columns(
  self : A2mAlignment,
  start : Int,
  end : Int,
) -> A2mAlignment raise A2mError {
  if start < 0 || end <= start || end > self.states.length() {
    raise A2mError("A2M column slice is out of bounds or empty")
  }
  let sequences : Array[A2mSequence] = []
  for sequence in self.sequences {
    sequences.push(
      A2mSequence::create(
        sequence.id,
        sequence.aligned_sequence[start:end].to_owned(),
        description=sequence.description,
      ),
    )
  }
  A2mAlignment::create(sequences, self.states[start:end].to_owned())
}

///|
/// Return a compact alignment summary.
pub fn A2mAlignment::summary(self : A2mAlignment) -> String {
  "A2mAlignment(sequences=" +
  self.sequences.length().to_string() +
  ", columns=" +
  self.states.length().to_string() +
  ", match=" +
  self.match_columns().to_string() +
  ", insertion=" +
  self.insertion_columns().to_string() +
  ")"
}

///|
/// Return a small state-aware A2M example.
pub fn a2m_example_text() -> String {
  ">reference SAM model reference\n" +
  "ACDefG-HIkLM\n" +
  ">query_one insertion and deletion\n" +
  "ACD..GTHI.LM\n" +
  ">query_two divergent homolog\n" +
  "ATDqrG-HI.L-\n"
}

///|
fn a2m_validate_alignment(alignment : A2mAlignment) -> Unit raise A2mError {
  if alignment.sequences.length() == 0 {
    raise A2mError("A2M alignment must contain at least one sequence")
  }
  let width = alignment.states.length()
  if width == 0 {
    raise A2mError("A2M state string must not be empty")
  }
  for column = 0; column < width; column = column + 1 {
    let state = alignment.states.unsafe_get(column).to_int()
    if state != 'D'.to_int() && state != 'I'.to_int() {
      raise A2mError(
        "A2M state string contains a character other than D or I at column " +
        column.to_string(),
      )
    }
  }
  for row = 0; row < alignment.sequences.length(); row = row + 1 {
    let sequence = alignment.sequences[row]
    a2m_validate_id(sequence.id)
    a2m_validate_description(sequence.description)
    if sequence.aligned_sequence.length() != width {
      raise A2mError(
        "A2M row '" +
        sequence.id +
        "' has width " +
        sequence.aligned_sequence.length().to_string() +
        "; expected " +
        width.to_string(),
      )
    }
    for column = 0; column < width; column = column + 1 {
      let code = sequence.aligned_sequence.unsafe_get(column).to_int()
      if code != '-'.to_int() && !a2m_is_upper(code) {
        raise A2mError(
          "Normalized A2M row '" +
          sequence.id +
          "' contains an invalid character at column " +
          column.to_string(),
        )
      }
    }
    if a2m_remove_gaps(sequence.aligned_sequence) != sequence.sequence {
      raise A2mError(
        "A2M row '" + sequence.id + "' has inconsistent ungapped sequence",
      )
    }
  }
}

///|
fn a2m_validate_row(alignment : A2mAlignment, row : Int) -> Unit raise A2mError {
  if row < 0 || row >= alignment.sequences.length() {
    raise A2mError("A2M row index is out of bounds")
  }
}

///|
fn a2m_validate_id(id : String) -> Unit raise A2mError {
  if id.length() == 0 {
    raise A2mError("A2M sequence identifier must not be empty")
  }
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    if code == ' '.to_int() ||
      code == '\t'.to_int() ||
      code == '\n'.to_int() ||
      code == '\r'.to_int() {
      raise A2mError("A2M sequence identifier must not contain whitespace")
    }
  }
}

///|
fn a2m_validate_description(description : String) -> Unit raise A2mError {
  for index = 0; index < description.length(); index = index + 1 {
    let code = description.unsafe_get(index).to_int()
    if code == '\n'.to_int() || code == '\r'.to_int() {
      raise A2mError("A2M description must not contain line breaks")
    }
  }
}

///|
fn a2m_validate_encoded_fragment(fragment : String) -> Unit raise A2mError {
  for index = 0; index < fragment.length(); index = index + 1 {
    let code = fragment.unsafe_get(index).to_int()
    if !a2m_is_upper(code) &&
      !a2m_is_lower(code) &&
      code != '-'.to_int() &&
      code != '.'.to_int() {
      raise A2mError(
        "Invalid A2M sequence character at wrapped-line offset " +
        index.to_string(),
      )
    }
  }
}

///|
fn a2m_normalize_aligned_sequence(sequence : String) -> String raise A2mError {
  let result = StringBuilder::new(size_hint=sequence.length())
  for index = 0; index < sequence.length(); index = index + 1 {
    let code = sequence.unsafe_get(index).to_int()
    if code == '-'.to_int() {
      result.write_char('-')
    } else if a2m_is_upper(code) {
      result.write_char(code.unsafe_to_char())
    } else if a2m_is_lower(code) {
      result.write_char(a2m_upper_code(code).unsafe_to_char())
    } else {
      raise A2mError(
        "Normalized A2M sequence contains an invalid character at column " +
        index.to_string(),
      )
    }
  }
  result.to_string()
}

///|
fn a2m_remove_gaps(sequence : String) -> String {
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
fn a2m_parse_header(header : String) -> (String, String) {
  for index = 0; index < header.length(); index = index + 1 {
    let code = header.unsafe_get(index).to_int()
    if code == ' '.to_int() || code == '\t'.to_int() {
      return (header[0:index].to_owned(), header[index:].trim().to_owned())
    }
  }
  (header, "")
}

///|
fn a2m_normalize_lines(text : String) -> Array[String] {
  let result : Array[String] = []
  for view in text.split("\n") {
    let mut line = view.to_owned()
    if line.length() > 0 &&
      line.unsafe_get(line.length() - 1).to_int() == '\r'.to_int() {
      line = line[0:line.length() - 1].to_owned()
    }
    result.push(line)
  }
  result
}

///|
fn a2m_is_upper(code : Int) -> Bool {
  code >= 'A'.to_int() && code <= 'Z'.to_int()
}

///|
fn a2m_is_lower(code : Int) -> Bool {
  code >= 'a'.to_int() && code <= 'z'.to_int()
}

///|
fn a2m_upper_code(code : Int) -> Int {
  if a2m_is_lower(code) {
    code - ('a'.to_int() - 'A'.to_int())
  } else {
    code
  }
}

///|
fn a2m_lower_code(code : Int) -> Int {
  if a2m_is_upper(code) {
    code + ('a'.to_int() - 'A'.to_int())
  } else {
    code
  }
}
