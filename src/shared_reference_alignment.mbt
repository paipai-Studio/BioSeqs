// Merge pairwise or multiple alignments that share one reference sequence.
//
// This implements the semantics of Biopython 1.86
// Alignment.from_alignments_with_same_reference. Insertions are represented
// as slots between reference residues; each slot is widened to the largest
// insertion present in any input without realigning the query sequences.

///|
/// Error raised by shared-reference alignment operations.
pub suberror SharedReferenceAlignmentError {
  SharedReferenceAlignmentError(String)
}

///|
/// One query row in an input or merged shared-reference alignment.
pub struct SharedReferenceSequence {
  name : String
  description : String
  sequence : String
  aligned_sequence : String
  start : Int
  end : Int
} derive(Eq, Debug)
pub extend SharedReferenceSequence with Eq::{equal, not_equal}
pub extend SharedReferenceSequence with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// An alignment containing one reference row and one or more query rows.
pub struct SharedReferenceInput {
  reference_name : String
  reference_description : String
  reference : String
  aligned_reference : String
  reference_start : Int
  reference_end : Int
  queries : Array[SharedReferenceSequence]
} derive(Eq, Debug)
pub extend SharedReferenceInput with Eq::{equal, not_equal}
pub extend SharedReferenceInput with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Counts comparing one merged query row with the shared reference.
pub struct SharedReferenceCounts {
  identities : Int
  mismatches : Int
  insertions : Int
  deletions : Int
  aligned_pairs : Int
} derive(Eq, Debug)
pub extend SharedReferenceCounts with Eq::{equal, not_equal}
pub extend SharedReferenceCounts with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// A multiple sequence alignment formed without realigning the input queries.
pub struct SharedReferenceAlignment {
  reference_name : String
  reference_description : String
  reference : String
  reference_start : Int
  reference_end : Int
  aligned_reference : String
  queries : Array[SharedReferenceSequence]
  insertion_widths : Array[Int]
} derive(Eq, Debug)
pub extend SharedReferenceAlignment with Eq::{equal, not_equal}
pub extend SharedReferenceAlignment with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Remove gap characters from an aligned sequence.
fn shared_reference_ungap(sequence : String) -> String {
  let result = StringBuilder::new(size_hint=sequence.length())
  let mut index = 0
  while index < sequence.length() {
    let residue = sequence.unsafe_get(index).unsafe_to_char()
    if residue != '-' {
      result.write_char(residue)
    }
    index = index + 1
  }
  result.to_string()
}

///|
/// Compare biological sequence text without considering letter case.
fn shared_reference_same_sequence(left : String, right : String) -> Bool {
  left.to_upper() == right.to_upper()
}

///|
/// Reject whitespace and, for raw sequences, embedded gap characters.
fn shared_reference_validate_sequence(
  sequence : String,
  label : String,
  allow_gaps : Bool,
) -> Unit raise SharedReferenceAlignmentError {
  let mut index = 0
  while index < sequence.length() {
    let residue = sequence.unsafe_get(index).unsafe_to_char()
    if residue == ' ' || residue == '\t' || residue == '\n' || residue == '\r' {
      raise SharedReferenceAlignmentError(label + " contains whitespace")
    }
    if !allow_gaps && residue == '-' {
      raise SharedReferenceAlignmentError(label + " must not contain gaps")
    }
    index = index + 1
  }
}

///|
/// Return a byte-based substring after validating its bounds.
fn shared_reference_substring(
  sequence : String,
  start : Int,
  end : Int,
) -> String raise SharedReferenceAlignmentError {
  if start < 0 || end < start || end > sequence.length() {
    raise SharedReferenceAlignmentError(
      "sequence coordinates are out of bounds",
    )
  }
  sequence[start:end].to_owned()
}

///|
/// Construct a query row.
///
/// If `sequence` is omitted, the ungapped aligned row is used. Supplying the
/// full sequence and `start` preserves local-alignment query coordinates.
pub fn shared_reference_sequence(
  name : String,
  aligned_sequence : String,
  description? : String = "",
  sequence? : String = "",
  start? : Int = 0,
) -> SharedReferenceSequence raise SharedReferenceAlignmentError {
  if start < 0 {
    raise SharedReferenceAlignmentError("query start must be non-negative")
  }
  shared_reference_validate_sequence(
    aligned_sequence, "aligned query sequence", true,
  )
  let aligned_ungapped = shared_reference_ungap(aligned_sequence)
  let raw_sequence = if sequence == "" {
    aligned_ungapped
  } else {
    shared_reference_validate_sequence(sequence, "query sequence", false)
    sequence
  }
  let end = start + aligned_ungapped.length()
  if end > raw_sequence.length() {
    raise SharedReferenceAlignmentError(
      "aligned query coordinates exceed the query sequence",
    )
  }
  let observed = shared_reference_substring(raw_sequence, start, end)
  if !shared_reference_same_sequence(observed, aligned_ungapped) {
    raise SharedReferenceAlignmentError(
      "aligned query does not match the query sequence at its coordinates",
    )
  }
  SharedReferenceSequence::{
    name,
    description,
    sequence: raw_sequence,
    aligned_sequence,
    start,
    end,
  }
}

///|
/// Construct an input alignment with a shared reference row.
pub fn shared_reference_input(
  reference : String,
  aligned_reference : String,
  queries : Array[SharedReferenceSequence],
  reference_name? : String = "reference",
  reference_description? : String = "",
  reference_start? : Int = 0,
) -> SharedReferenceInput raise SharedReferenceAlignmentError {
  if reference.length() == 0 {
    raise SharedReferenceAlignmentError("reference sequence must not be empty")
  }
  if reference_start < 0 {
    raise SharedReferenceAlignmentError("reference start must be non-negative")
  }
  if queries.length() == 0 {
    raise SharedReferenceAlignmentError(
      "an input alignment must contain at least one query",
    )
  }
  shared_reference_validate_sequence(reference, "reference sequence", false)
  shared_reference_validate_sequence(
    aligned_reference, "aligned reference sequence", true,
  )
  let reference_segment = shared_reference_ungap(aligned_reference)
  let reference_end = reference_start + reference_segment.length()
  if reference_end > reference.length() {
    raise SharedReferenceAlignmentError(
      "aligned reference coordinates exceed the reference sequence",
    )
  }
  let expected = shared_reference_substring(
    reference, reference_start, reference_end,
  )
  if !shared_reference_same_sequence(expected, reference_segment) {
    raise SharedReferenceAlignmentError(
      "aligned reference does not match the reference sequence at its coordinates",
    )
  }
  let query_copy : Array[SharedReferenceSequence] = []
  for query in queries {
    if query.aligned_sequence.length() != aligned_reference.length() {
      raise SharedReferenceAlignmentError(
        "all rows in an input alignment must have the same length",
      )
    }
    query_copy.push(query)
  }
  SharedReferenceInput::{
    reference_name,
    reference_description,
    reference,
    aligned_reference,
    reference_start,
    reference_end,
    queries: query_copy,
  }
}

///|
/// Adapt an existing BioSeqs PairwiseAlignment to a shared-reference input.
pub fn shared_reference_input_from_pairwise(
  alignment : PairwiseAlignment,
  reference_name? : String = "reference",
  reference_description? : String = "",
  query_name? : String = "query",
  query_description? : String = "",
) -> SharedReferenceInput raise SharedReferenceAlignmentError {
  let query = shared_reference_sequence(
    query_name,
    alignment.aligned2(),
    description=query_description,
    sequence=alignment.query,
    start=alignment.query_start,
  )
  let input = shared_reference_input(
    alignment.target,
    alignment.aligned1(),
    [query],
    reference_name~,
    reference_description~,
    reference_start=alignment.target_start,
  )
  if input.reference_end != alignment.target_end {
    raise SharedReferenceAlignmentError(
      "pairwise target coordinates disagree with its aligned reference",
    )
  }
  if query.end != alignment.query_end {
    raise SharedReferenceAlignmentError(
      "pairwise query coordinates disagree with its aligned query",
    )
  }
  input
}

///|
/// Count reference-gap columns at every boundary in one input alignment.
fn shared_reference_slot_widths(
  aligned_reference : String,
  span : Int,
) -> Array[Int] raise SharedReferenceAlignmentError {
  let widths = Array::make(span + 1, 0)
  let mut reference_offset = 0
  let mut column = 0
  while column < aligned_reference.length() {
    let residue = aligned_reference.unsafe_get(column).unsafe_to_char()
    if residue == '-' {
      widths[reference_offset] = widths[reference_offset] + 1
    } else {
      reference_offset = reference_offset + 1
      if reference_offset > span {
        raise SharedReferenceAlignmentError(
          "aligned reference is longer than its declared span",
        )
      }
    }
    column = column + 1
  }
  if reference_offset != span {
    raise SharedReferenceAlignmentError(
      "aligned reference is shorter than its declared span",
    )
  }
  widths
}

///|
/// Render a reference row with synchronized insertion slots.
fn shared_reference_render_reference(
  reference : String,
  reference_start : Int,
  span : Int,
  insertion_widths : Array[Int],
) -> String {
  let result = StringBuilder::new()
  let mut boundary = 0
  while boundary <= span {
    result.write_string("-".repeat(insertion_widths[boundary]))
    if boundary < span {
      result.write_char(
        reference.unsafe_get(reference_start + boundary).unsafe_to_char(),
      )
    }
    boundary = boundary + 1
  }
  result.to_string()
}

///|
/// Project one query into the synchronized insertion slots.
fn shared_reference_project_query(
  aligned_reference : String,
  aligned_query : String,
  span : Int,
  insertion_widths : Array[Int],
) -> String raise SharedReferenceAlignmentError {
  let insertions : Array[StringBuilder] = []
  let bases = Array::make(span, '-')
  let mut boundary = 0
  while boundary <= span {
    insertions.push(StringBuilder::new())
    boundary = boundary + 1
  }
  let mut reference_offset = 0
  let mut column = 0
  while column < aligned_reference.length() {
    let reference_residue = aligned_reference
      .unsafe_get(column)
      .unsafe_to_char()
    let query_residue = aligned_query.unsafe_get(column).unsafe_to_char()
    if reference_residue == '-' {
      insertions[reference_offset].write_char(query_residue)
    } else {
      if reference_offset >= span {
        raise SharedReferenceAlignmentError(
          "query projection exceeded the reference span",
        )
      }
      bases[reference_offset] = query_residue
      reference_offset = reference_offset + 1
    }
    column = column + 1
  }
  if reference_offset != span {
    raise SharedReferenceAlignmentError(
      "query projection did not cover the reference span",
    )
  }
  let result = StringBuilder::new()
  boundary = 0
  while boundary <= span {
    let insertion = insertions[boundary].to_string()
    result.write_string(insertion)
    let padding = insertion_widths[boundary] - insertion.length()
    if padding < 0 {
      raise SharedReferenceAlignmentError(
        "query insertion is wider than its synchronized slot",
      )
    }
    result.write_string("-".repeat(padding))
    if boundary < span {
      result.write_char(bases[boundary])
    }
    boundary = boundary + 1
  }
  result.to_string()
}

///|
/// Merge alignments whose first row is the same reference sequence.
///
/// Queries are not realigned against each other. At each reference boundary,
/// the insertion slot is widened to the largest width found in any input.
pub fn alignments_with_same_reference(
  inputs : Array[SharedReferenceInput],
) -> SharedReferenceAlignment raise SharedReferenceAlignmentError {
  if inputs.length() == 0 {
    raise SharedReferenceAlignmentError("no alignments were provided")
  }
  let first = inputs[0]
  let reference_length = first.reference.length()
  let reference_start = first.reference_start
  let reference_end = first.reference_end
  let span = reference_end - reference_start
  let insertion_widths = Array::make(span + 1, 0)
  let mut query_count = 0
  for input in inputs {
    if input.reference.length() != reference_length {
      raise SharedReferenceAlignmentError(
        "all reference sequences must have the same length",
      )
    }
    if !shared_reference_same_sequence(input.reference, first.reference) {
      raise SharedReferenceAlignmentError("all reference sequences must match")
    }
    if input.reference_start != reference_start ||
      input.reference_end != reference_end {
      raise SharedReferenceAlignmentError(
        "reference coordinates are inconsistent across alignments",
      )
    }
    let widths = shared_reference_slot_widths(input.aligned_reference, span)
    let mut boundary = 0
    while boundary <= span {
      if widths[boundary] > insertion_widths[boundary] {
        insertion_widths[boundary] = widths[boundary]
      }
      boundary = boundary + 1
    }
    query_count = query_count + input.queries.length()
  }
  if query_count == 0 {
    raise SharedReferenceAlignmentError(
      "the alignments do not contain query sequences",
    )
  }
  let aligned_reference = shared_reference_render_reference(
    first.reference,
    reference_start,
    span,
    insertion_widths,
  )
  let merged_queries : Array[SharedReferenceSequence] = []
  for input in inputs {
    for query in input.queries {
      let aligned_sequence = shared_reference_project_query(
        input.aligned_reference,
        query.aligned_sequence,
        span,
        insertion_widths,
      )
      merged_queries.push(SharedReferenceSequence::{
        name: query.name,
        description: query.description,
        sequence: query.sequence,
        aligned_sequence,
        start: query.start,
        end: query.end,
      })
    }
  }
  SharedReferenceAlignment::{
    reference_name: first.reference_name,
    reference_description: first.reference_description,
    reference: first.reference,
    reference_start,
    reference_end,
    aligned_reference,
    queries: merged_queries,
    insertion_widths,
  }
}

///|
/// Return the number of rows, including the reference.
pub fn SharedReferenceAlignment::num_sequences(
  self : SharedReferenceAlignment,
) -> Int {
  self.queries.length() + 1
}

///|
/// Return the number of query rows.
pub fn SharedReferenceAlignment::num_queries(
  self : SharedReferenceAlignment,
) -> Int {
  self.queries.length()
}

///|
/// Return the number of alignment columns.
pub fn SharedReferenceAlignment::alignment_length(
  self : SharedReferenceAlignment,
) -> Int {
  self.aligned_reference.length()
}

///|
/// Return an aligned row by MSA row index; row zero is the reference.
pub fn SharedReferenceAlignment::row(
  self : SharedReferenceAlignment,
  index : Int,
) -> String? {
  if index == 0 {
    Some(self.aligned_reference)
  } else if index > 0 && index <= self.queries.length() {
    Some(self.queries[index - 1].aligned_sequence)
  } else {
    None
  }
}

///|
/// Return a row name by MSA row index.
pub fn SharedReferenceAlignment::row_name(
  self : SharedReferenceAlignment,
  index : Int,
) -> String? {
  if index == 0 {
    Some(self.reference_name)
  } else if index > 0 && index <= self.queries.length() {
    Some(self.queries[index - 1].name)
  } else {
    None
  }
}

///|
/// Find the first merged query by name.
pub fn SharedReferenceAlignment::find_query(
  self : SharedReferenceAlignment,
  name : String,
) -> SharedReferenceSequence? {
  for query in self.queries {
    if query.name == name {
      return Some(query)
    }
  }
  None
}

///|
/// Return a column containing the reference followed by all queries.
pub fn SharedReferenceAlignment::column(
  self : SharedReferenceAlignment,
  index : Int,
) -> String? {
  if index < 0 || index >= self.alignment_length() {
    return None
  }
  let result = StringBuilder::new(size_hint=self.num_sequences())
  result.write_char(self.aligned_reference.unsafe_get(index).unsafe_to_char())
  for query in self.queries {
    result.write_char(query.aligned_sequence.unsafe_get(index).unsafe_to_char())
  }
  Some(result.to_string())
}

///|
/// Map an absolute reference position to an alignment column.
pub fn SharedReferenceAlignment::reference_to_column(
  self : SharedReferenceAlignment,
  position : Int,
) -> Int? {
  if position < self.reference_start || position >= self.reference_end {
    return None
  }
  let mut reference_position = self.reference_start
  let mut column = 0
  while column < self.aligned_reference.length() {
    let residue = self.aligned_reference.unsafe_get(column).unsafe_to_char()
    if residue != '-' {
      if reference_position == position {
        return Some(column)
      }
      reference_position = reference_position + 1
    }
    column = column + 1
  }
  None
}

///|
/// Map an alignment column to an absolute reference position.
pub fn SharedReferenceAlignment::column_to_reference(
  self : SharedReferenceAlignment,
  column : Int,
) -> Int? {
  if column < 0 || column >= self.alignment_length() {
    return None
  }
  let mut reference_position = self.reference_start
  let mut index = 0
  while index <= column {
    let residue = self.aligned_reference.unsafe_get(index).unsafe_to_char()
    if residue != '-' {
      if index == column {
        return Some(reference_position)
      }
      reference_position = reference_position + 1
    } else if index == column {
      return None
    }
    index = index + 1
  }
  None
}

///|
/// Map an absolute query position to an alignment column.
pub fn SharedReferenceAlignment::query_to_column(
  self : SharedReferenceAlignment,
  query_index : Int,
  position : Int,
) -> Int? {
  if query_index < 0 || query_index >= self.queries.length() {
    return None
  }
  let query = self.queries[query_index]
  if position < query.start || position >= query.end {
    return None
  }
  let mut query_position = query.start
  let mut column = 0
  while column < query.aligned_sequence.length() {
    let residue = query.aligned_sequence.unsafe_get(column).unsafe_to_char()
    if residue != '-' {
      if query_position == position {
        return Some(column)
      }
      query_position = query_position + 1
    }
    column = column + 1
  }
  None
}

///|
/// Map an alignment column to an absolute query position.
pub fn SharedReferenceAlignment::column_to_query(
  self : SharedReferenceAlignment,
  query_index : Int,
  column : Int,
) -> Int? {
  if query_index < 0 ||
    query_index >= self.queries.length() ||
    column < 0 ||
    column >= self.alignment_length() {
    return None
  }
  let query = self.queries[query_index]
  let mut query_position = query.start
  let mut index = 0
  while index <= column {
    let residue = query.aligned_sequence.unsafe_get(index).unsafe_to_char()
    if residue != '-' {
      if index == column {
        return Some(query_position)
      }
      query_position = query_position + 1
    } else if index == column {
      return None
    }
    index = index + 1
  }
  None
}

///|
/// Map an absolute query position through the MSA to the reference.
pub fn SharedReferenceAlignment::query_to_reference(
  self : SharedReferenceAlignment,
  query_index : Int,
  position : Int,
) -> Int? {
  match self.query_to_column(query_index, position) {
    Some(column) => self.column_to_reference(column)
    None => None
  }
}

///|
/// Map an absolute reference position through the MSA to one query.
pub fn SharedReferenceAlignment::reference_to_query(
  self : SharedReferenceAlignment,
  query_index : Int,
  position : Int,
) -> Int? {
  match self.reference_to_column(position) {
    Some(column) => self.column_to_query(query_index, column)
    None => None
  }
}

///|
/// Count matches, mismatches, insertions, and deletions for one query.
pub fn SharedReferenceAlignment::counts(
  self : SharedReferenceAlignment,
  query_index : Int,
) -> SharedReferenceCounts? {
  if query_index < 0 || query_index >= self.queries.length() {
    return None
  }
  let query = self.queries[query_index]
  let mut identities = 0
  let mut mismatches = 0
  let mut insertions = 0
  let mut deletions = 0
  let mut column = 0
  while column < self.alignment_length() {
    let reference_residue = self.aligned_reference
      .unsafe_get(column)
      .unsafe_to_char()
    let query_residue = query.aligned_sequence
      .unsafe_get(column)
      .unsafe_to_char()
    if reference_residue == '-' {
      if query_residue != '-' {
        insertions = insertions + 1
      }
    } else if query_residue == '-' {
      deletions = deletions + 1
    } else if reference_residue.to_string().to_upper() ==
      query_residue.to_string().to_upper() {
      identities = identities + 1
    } else {
      mismatches = mismatches + 1
    }
    column = column + 1
  }
  Some(SharedReferenceCounts::{
    identities,
    mismatches,
    insertions,
    deletions,
    aligned_pairs: identities + mismatches,
  })
}

///|
/// Percent identity among residue pairs, excluding gaps.
pub fn SharedReferenceCounts::identity_percent(
  self : SharedReferenceCounts,
) -> Double {
  if self.aligned_pairs == 0 {
    0.0
  } else {
    self.identities.to_double() / self.aligned_pairs.to_double() * 100.0
  }
}

///|
/// Convert the merged alignment to the existing MultipleSeqAlignment type.
pub fn SharedReferenceAlignment::to_msa(
  self : SharedReferenceAlignment,
) -> MultipleSeqAlignment raise AlignIOError {
  let records : Array[SeqRecord] = [
    SeqRecord::new(
      Seq::new(self.aligned_reference),
      id=self.reference_name,
      name=self.reference_name,
      description=self.reference_description,
    ),
  ]
  for query in self.queries {
    records.push(
      SeqRecord::new(
        Seq::new(query.aligned_sequence),
        id=query.name,
        name=query.name,
        description=query.description,
      ),
    )
  }
  MultipleSeqAlignment::new(records)
}

///|
/// Serialize the merged alignment as aligned FASTA.
pub fn SharedReferenceAlignment::to_fasta(
  self : SharedReferenceAlignment,
) -> String {
  let result = StringBuilder::new()
  result.write_char('>')
  result.write_string(self.reference_name)
  if self.reference_description != "" {
    result.write_char(' ')
    result.write_string(self.reference_description)
  }
  result.write_char('\n')
  result.write_string(self.aligned_reference)
  result.write_char('\n')
  for query in self.queries {
    result.write_char('>')
    result.write_string(query.name)
    if query.description != "" {
      result.write_char(' ')
      result.write_string(query.description)
    }
    result.write_char('\n')
    result.write_string(query.aligned_sequence)
    result.write_char('\n')
  }
  result.to_string()
}

///|
/// Return a concise summary of a merged shared-reference alignment.
pub fn SharedReferenceAlignment::summary(
  self : SharedReferenceAlignment,
) -> String {
  "Shared-reference alignment: " +
  self.num_sequences().to_string() +
  " sequences, " +
  self.alignment_length().to_string() +
  " columns, reference " +
  self.reference_start.to_string() +
  ":" +
  self.reference_end.to_string()
}

///|
/// Deterministic sample combining one PWA and one three-row alignment.
pub fn shared_reference_alignment_sample() -> SharedReferenceAlignment raise SharedReferenceAlignmentError {
  let first_query = shared_reference_sequence(
    "seq1",
    "AC-T",
    description="one deletion",
  )
  let first = shared_reference_input(
    "ACGT",
    "ACGT",
    [first_query],
    reference_name="reference",
    reference_description="shared reference",
  )
  let second_query = shared_reference_sequence(
    "seq2",
    "ACGGT",
    description="one insertion",
  )
  let third_query = shared_reference_sequence(
    "seq3",
    "A---T",
    description="two deletions",
  )
  let second = shared_reference_input(
    "ACGT",
    "ACG-T",
    [second_query, third_query],
    reference_name="reference",
    reference_description="shared reference",
  )
  alignments_with_same_reference([first, second])
}
