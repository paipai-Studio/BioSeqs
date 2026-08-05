// Biopython-compatible Bio.Align.stockholm support.
//
// This module is intentionally separate from stockholm.mbt. The older module
// exposes a permissive block-oriented API, while this file models modern,
// coordinate-aware Stockholm alignments and their typed annotations.

///|
/// Error raised for malformed Stockholm data or invalid operations.
pub suberror AlignStockholmError {
  AlignStockholmError(String)
}

///|
/// A mapped GF or GS annotation.
///
/// `code` is the Stockholm tag and `feature` is Biopython's descriptive name.
pub struct AlignStockholmAnnotation {
  code : String
  feature : String
  value : String
} derive(Eq, Debug)

///|
/// A mapped GC annotation with one character per retained alignment column.
pub struct AlignStockholmColumnAnnotation {
  code : String
  feature : String
  value : String
} derive(Eq, Debug)

///|
/// A mapped GR annotation.
///
/// `value` contains one character per ungapped residue. `aligned_value`
/// contains one character per retained alignment column.
pub struct AlignStockholmLetterAnnotation {
  code : String
  feature : String
  value : String
  aligned_value : String
} derive(Eq, Debug)

///|
/// One structured Stockholm reference.
pub struct AlignStockholmReference {
  number : Int
  medline : String
  title : String
  author : String
  location : String
  comment : String
} derive(Eq, Debug)

///|
/// One structured database reference.
pub struct AlignStockholmDatabaseReference {
  reference : String
  comment : String
} derive(Eq, Debug)

///|
/// One nested-domain annotation.
pub struct AlignStockholmNestedDomain {
  accession : String
  location : String
} derive(Eq, Debug)

///|
/// One sequence row and its GS/GR annotations.
pub struct AlignStockholmSequence {
  id : String
  description : String
  sequence : String
  aligned_sequence : String
  annotations : Array[AlignStockholmAnnotation]
  dbxrefs : Array[String]
  letter_annotations : Array[AlignStockholmLetterAnnotation]
} derive(Eq, Debug)

///|
/// A coordinate-aware Stockholm 1.0 alignment.
///
/// `operations` contains one `M`, `D`, or `I` per retained column. A source
/// `-` gap marks a deletion column and a source `.` gap marks an insertion
/// column, matching Biopython's `Bio.Align.stockholm` behavior.
pub struct AlignStockholmAlignment {
  version : String
  sequences : Array[AlignStockholmSequence]
  annotations : Array[AlignStockholmAnnotation]
  column_annotations : Array[AlignStockholmColumnAnnotation]
  references : Array[AlignStockholmReference]
  database_references : Array[AlignStockholmDatabaseReference]
  nested_domains : Array[AlignStockholmNestedDomain]
  operations : String
  source_columns : Int
  removed_all_gap_columns : Int
} derive(Eq, Debug)

///|
/// Pairwise or all-pairs Stockholm alignment statistics.
pub struct AlignStockholmCounts {
  pairs : Int
  aligned : Int
  identities : Int
  mismatches : Int
  gap_columns : Int
  double_gap_columns : Int
  gap_opens : Int
} derive(Eq, Debug)

///|
priv struct AlignStockholmRawAnnotation {
  code : String
  value : String
}

///|
priv struct AlignStockholmRawSequenceAnnotation {
  sequence_id : String
  code : String
  value : String
}

///|
priv struct AlignStockholmReferenceBuilder {
  number : Int
  mut medline : String
  mut title : String
  mut author : String
  mut location : String
  comment : String
}

///|
priv struct AlignStockholmDatabaseReferenceBuilder {
  reference : String
  mut comment : String
}

///|
priv struct AlignStockholmNestedDomainBuilder {
  accession : String
  mut location : String
}

///|
/// Construct a generic GF or GS annotation.
pub fn AlignStockholmAnnotation::create(
  code : String,
  value : String,
  feature? : String? = None,
) -> AlignStockholmAnnotation raise AlignStockholmError {
  align_stockholm_validate_annotation_code(code)
  AlignStockholmAnnotation::{
    code,
    feature: match feature {
      Some(name) => name
      None => code
    },
    value,
  }
}

///|
/// Construct a generic GC annotation.
pub fn AlignStockholmColumnAnnotation::create(
  code : String,
  value : String,
  feature? : String? = None,
) -> AlignStockholmColumnAnnotation raise AlignStockholmError {
  align_stockholm_validate_annotation_code(code)
  AlignStockholmColumnAnnotation::{
    code,
    feature: match feature {
      Some(name) => name
      None => align_stockholm_gc_feature(code)
    },
    value,
  }
}

///|
/// Construct a generic GR annotation.
pub fn AlignStockholmLetterAnnotation::create(
  code : String,
  value : String,
  aligned_value : String,
  feature? : String? = None,
) -> AlignStockholmLetterAnnotation raise AlignStockholmError {
  align_stockholm_validate_annotation_code(code)
  AlignStockholmLetterAnnotation::{
    code,
    feature: match feature {
      Some(name) => name
      None => align_stockholm_gr_feature(code)
    },
    value,
    aligned_value,
  }
}

///|
/// Construct one standalone Stockholm row.
pub fn AlignStockholmSequence::create(
  id : String,
  aligned_sequence : String,
  description? : String = "",
  annotations? : Array[AlignStockholmAnnotation] = [],
  dbxrefs? : Array[String] = [],
  letter_annotations? : Array[AlignStockholmLetterAnnotation] = [],
) -> AlignStockholmSequence raise AlignStockholmError {
  align_stockholm_validate_id(id)
  if aligned_sequence.length() == 0 {
    raise AlignStockholmError("Stockholm aligned sequence must not be empty")
  }
  let normalized = align_stockholm_normalize_row(aligned_sequence)
  let copied_annotations : Array[AlignStockholmAnnotation] = []
  for annotation in annotations {
    copied_annotations.push(annotation)
  }
  let copied_dbxrefs : Array[String] = []
  for dbxref in dbxrefs {
    copied_dbxrefs.push(dbxref)
  }
  let copied_letter_annotations : Array[AlignStockholmLetterAnnotation] = []
  for annotation in letter_annotations {
    copied_letter_annotations.push(annotation)
  }
  let sequence = AlignStockholmSequence::{
    id,
    description,
    sequence: align_stockholm_remove_gaps(normalized),
    aligned_sequence: normalized,
    annotations: copied_annotations,
    dbxrefs: copied_dbxrefs,
    letter_annotations: copied_letter_annotations,
  }
  align_stockholm_validate_sequence_annotations(sequence)
  sequence
}

///|
/// Construct and validate one coordinate-aware Stockholm alignment.
pub fn AlignStockholmAlignment::create(
  sequences : Array[AlignStockholmSequence],
  annotations? : Array[AlignStockholmAnnotation] = [],
  column_annotations? : Array[AlignStockholmColumnAnnotation] = [],
  references? : Array[AlignStockholmReference] = [],
  database_references? : Array[AlignStockholmDatabaseReference] = [],
  nested_domains? : Array[AlignStockholmNestedDomain] = [],
  operations? : String = "",
  source_columns? : Int = 0,
  removed_all_gap_columns? : Int = 0,
) -> AlignStockholmAlignment raise AlignStockholmError {
  let copied_sequences : Array[AlignStockholmSequence] = []
  for sequence in sequences {
    copied_sequences.push(sequence)
  }
  let copied_annotations : Array[AlignStockholmAnnotation] = []
  for annotation in annotations {
    copied_annotations.push(annotation)
  }
  let copied_column_annotations : Array[AlignStockholmColumnAnnotation] = []
  for annotation in column_annotations {
    copied_column_annotations.push(annotation)
  }
  let copied_references : Array[AlignStockholmReference] = []
  for reference in references {
    copied_references.push(reference)
  }
  let copied_database_references : Array[AlignStockholmDatabaseReference] = []
  for reference in database_references {
    copied_database_references.push(reference)
  }
  let copied_nested_domains : Array[AlignStockholmNestedDomain] = []
  for domain in nested_domains {
    copied_nested_domains.push(domain)
  }
  let width = if copied_sequences.length() == 0 {
    0
  } else {
    copied_sequences[0].aligned_sequence.length()
  }
  let actual_operations = if operations.length() == 0 && width > 0 {
    align_stockholm_infer_operations(copied_sequences)
  } else {
    operations
  }
  let actual_source_columns = if source_columns == 0 {
    width + removed_all_gap_columns
  } else {
    source_columns
  }
  let alignment = AlignStockholmAlignment::{
    version: "1.0",
    sequences: copied_sequences,
    annotations: copied_annotations,
    column_annotations: copied_column_annotations,
    references: copied_references,
    database_references: copied_database_references,
    nested_domains: copied_nested_domains,
    operations: actual_operations,
    source_columns: actual_source_columns,
    removed_all_gap_columns,
  }
  align_stockholm_validate_alignment(alignment)
  alignment
}

///|
/// Construct an alignment from identifiers and printed rows.
///
/// Source `.` and `-` gaps retain their insertion/deletion operation class.
/// Columns containing only gaps are removed.
pub fn align_stockholm_from_aligned(
  ids : Array[String],
  aligned_sequences : Array[String],
) -> AlignStockholmAlignment raise AlignStockholmError {
  if ids.length() == 0 {
    raise AlignStockholmError(
      "Stockholm alignment must contain at least one sequence",
    )
  }
  if ids.length() != aligned_sequences.length() {
    raise AlignStockholmError(
      "Stockholm identifiers and aligned rows must have equal lengths",
    )
  }
  let width = aligned_sequences[0].length()
  if width == 0 {
    raise AlignStockholmError(
      "Stockholm alignment must contain at least one column",
    )
  }
  for row in aligned_sequences {
    if row.length() != width {
      raise AlignStockholmError("Stockholm aligned rows must have equal widths")
    }
  }
  let prepared = align_stockholm_prepare_rows(aligned_sequences)
  if prepared.0.length() == 0 || prepared.0[0].length() == 0 {
    raise AlignStockholmError(
      "Stockholm alignment contains no non-gap alignment column",
    )
  }
  let sequences : Array[AlignStockholmSequence] = []
  for index = 0; index < ids.length(); index = index + 1 {
    sequences.push(
      AlignStockholmSequence::create(ids[index], prepared.0[index]),
    )
  }
  AlignStockholmAlignment::create(
    sequences,
    operations=prepared.1,
    source_columns=width,
    removed_all_gap_columns=align_stockholm_count_true(prepared.2),
  )
}

///|
/// Parse exactly one Stockholm 1.0 alignment.
pub fn align_stockholm_parse(
  text : String,
) -> AlignStockholmAlignment raise AlignStockholmError {
  let alignments = align_stockholm_parse_all(text)
  if alignments.length() == 0 {
    raise AlignStockholmError("Stockholm input contains no alignment")
  }
  if alignments.length() != 1 {
    raise AlignStockholmError(
      "Expected one Stockholm alignment, found " +
      alignments.length().to_string(),
    )
  }
  alignments[0]
}

///|
/// Parse all concatenated Stockholm 1.0 alignments.
pub fn align_stockholm_parse_all(
  text : String,
) -> Array[AlignStockholmAlignment] raise AlignStockholmError {
  let lines = align_stockholm_normalize_lines(text)
  let alignments : Array[AlignStockholmAlignment] = []
  let mut index = 0
  while index < lines.length() {
    while index < lines.length() && lines[index].trim().length() == 0 {
      index = index + 1
    }
    if index >= lines.length() {
      break
    }
    if lines[index].trim() != "# STOCKHOLM 1.0" {
      raise AlignStockholmError(
        "Expected '# STOCKHOLM 1.0' at line " + (index + 1).to_string(),
      )
    }
    let start = index + 1
    index = start
    while index < lines.length() && lines[index].trim() != "//" {
      if lines[index].trim() == "# STOCKHOLM 1.0" {
        raise AlignStockholmError(
          "Stockholm alignment is missing its // terminator",
        )
      }
      index = index + 1
    }
    if index >= lines.length() {
      raise AlignStockholmError(
        "Stockholm alignment is missing its // terminator",
      )
    }
    alignments.push(align_stockholm_parse_record(lines, start, index))
    index = index + 1
  }
  alignments
}

///|
/// Write one canonical Stockholm 1.0 alignment.
pub fn align_stockholm_write(
  alignment : AlignStockholmAlignment,
) -> String raise AlignStockholmError {
  align_stockholm_validate_alignment(alignment)
  let rows = alignment.num_sequences()
  let columns = alignment.alignment_length()
  if rows == 0 {
    raise AlignStockholmError("Must have at least one Stockholm sequence")
  }
  if columns == 0 {
    raise AlignStockholmError("Non-empty Stockholm sequences are required")
  }
  let output = StringBuilder::new()
  output.write_string("# STOCKHOLM 1.0\n")
  for annotation in alignment.annotations {
    if !align_stockholm_is_known_gf(annotation.code) {
      raise AlignStockholmError(
        "Unknown Stockholm GF annotation " + annotation.code,
      )
    }
    if annotation.code == "CC" {
      output.write_string(
        align_stockholm_wrap_text("#=GF CC   ", annotation.value),
      )
    } else {
      output.write_string(
        "#=GF " + annotation.code + "   " + annotation.value + "\n",
      )
    }
  }
  for domain in alignment.nested_domains {
    if domain.accession.length() > 0 {
      output.write_string("#=GF NE   " + domain.accession + "\n")
    }
    if domain.location.length() > 0 {
      output.write_string("#=GF NL   " + domain.location + "\n")
    }
  }
  for reference in alignment.references {
    if reference.comment.length() > 0 {
      output.write_string(
        align_stockholm_wrap_text("#=GF RC   ", reference.comment),
      )
    }
    output.write_string("#=GF RN   [" + reference.number.to_string() + "]\n")
    if reference.medline.length() > 0 {
      output.write_string("#=GF RM   " + reference.medline + "\n")
    }
    if reference.title.length() > 0 {
      output.write_string(
        align_stockholm_wrap_text("#=GF RT   ", reference.title),
      )
    }
    if reference.author.length() > 0 {
      output.write_string("#=GF RA   " + reference.author + "\n")
    }
    if reference.location.length() > 0 {
      output.write_string("#=GF RL   " + reference.location + "\n")
    }
  }
  for reference in alignment.database_references {
    output.write_string("#=GF DR   " + reference.reference + "\n")
    if reference.comment.length() > 0 {
      output.write_string("#=GF DC   " + reference.comment + "\n")
    }
  }
  output.write_string("#=GF SQ   " + rows.to_string() + "\n")
  let mut name_width = 0
  for sequence in alignment.sequences {
    if sequence.id.length() > name_width {
      name_width = sequence.id.length()
    }
  }
  let start = align_stockholm_max(name_width, 20) + 12
  for sequence in alignment.sequences {
    let padded_name = align_stockholm_pad_right(sequence.id, name_width)
    for annotation in sequence.annotations {
      output.write_string(
        "#=GS " +
        padded_name +
        "  " +
        annotation.code +
        " " +
        annotation.value +
        "\n",
      )
    }
    if sequence.description.length() > 0 {
      output.write_string(
        "#=GS " + padded_name + "  DE " + sequence.description + "\n",
      )
    }
    for dbxref in sequence.dbxrefs {
      output.write_string("#=GS " + padded_name + "  DR " + dbxref + "\n")
    }
  }
  for row = 0; row < alignment.sequences.length(); row = row + 1 {
    let sequence = alignment.sequences[row]
    output.write_string(align_stockholm_pad_right(sequence.id, start))
    output.write_string(
      align_stockholm_render_row(
        sequence.aligned_sequence,
        alignment.operations,
      ),
    )
    output.write_char('\n')
    let padded_name = align_stockholm_pad_right(sequence.id, name_width)
    for annotation in sequence.letter_annotations {
      output.write_string(
        align_stockholm_pad_right(
          "#=GR " + padded_name + "  " + annotation.code + " ",
          start,
        ),
      )
      output.write_string(
        align_stockholm_render_letter_annotation(sequence, annotation.value),
      )
      output.write_char('\n')
    }
  }
  for annotation in alignment.column_annotations {
    output.write_string(
      align_stockholm_pad_right("#=GC " + annotation.code + " ", start),
    )
    output.write_string(annotation.value)
    output.write_char('\n')
  }
  output.write_string("//\n")
  output.to_string()
}

///|
/// Write concatenated Stockholm records.
pub fn align_stockholm_write_all(
  alignments : Array[AlignStockholmAlignment],
) -> String raise AlignStockholmError {
  let output = StringBuilder::new()
  for alignment in alignments {
    output.write_string(align_stockholm_write(alignment))
  }
  output.to_string()
}

///|
/// Return the number of rows.
pub fn AlignStockholmAlignment::num_sequences(
  self : AlignStockholmAlignment,
) -> Int {
  self.sequences.length()
}

///|
/// Return the retained alignment width.
pub fn AlignStockholmAlignment::alignment_length(
  self : AlignStockholmAlignment,
) -> Int {
  if self.sequences.length() == 0 {
    0
  } else {
    self.sequences[0].aligned_sequence.length()
  }
}

///|
/// Return the source width before all-gap columns were removed.
pub fn AlignStockholmAlignment::source_alignment_length(
  self : AlignStockholmAlignment,
) -> Int {
  self.source_columns
}

///|
/// Find a sequence by exact identifier.
pub fn AlignStockholmAlignment::find_sequence(
  self : AlignStockholmAlignment,
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
/// Return all GF values matching a Stockholm code or mapped feature.
pub fn AlignStockholmAlignment::annotation_values(
  self : AlignStockholmAlignment,
  name : String,
) -> Array[String] {
  let values : Array[String] = []
  for annotation in self.annotations {
    if annotation.code == name || annotation.feature == name {
      values.push(annotation.value)
    }
  }
  values
}

///|
/// Return the first GF value matching a Stockholm code or mapped feature.
pub fn AlignStockholmAlignment::annotation(
  self : AlignStockholmAlignment,
  name : String,
) -> String? {
  for annotation in self.annotations {
    if annotation.code == name || annotation.feature == name {
      return Some(annotation.value)
    }
  }
  None
}

///|
/// Return a GC value matching a Stockholm code or mapped feature.
pub fn AlignStockholmAlignment::column_annotation(
  self : AlignStockholmAlignment,
  name : String,
) -> String? {
  for annotation in self.column_annotations {
    if annotation.code == name || annotation.feature == name {
      return Some(annotation.value)
    }
  }
  None
}

///|
/// Return the first GS value matching a Stockholm code or mapped feature.
pub fn AlignStockholmSequence::annotation(
  self : AlignStockholmSequence,
  name : String,
) -> String? {
  for annotation in self.annotations {
    if annotation.code == name || annotation.feature == name {
      return Some(annotation.value)
    }
  }
  None
}

///|
/// Return a residue-level GR value matching a code or mapped feature.
pub fn AlignStockholmSequence::letter_annotation(
  self : AlignStockholmSequence,
  name : String,
) -> String? {
  for annotation in self.letter_annotations {
    if annotation.code == name || annotation.feature == name {
      return Some(annotation.value)
    }
  }
  None
}

///|
/// Return one printed alignment column.
pub fn AlignStockholmAlignment::column(
  self : AlignStockholmAlignment,
  column : Int,
) -> String? {
  if column < 0 || column >= self.alignment_length() {
    return None
  }
  let output = StringBuilder::new(size_hint=self.sequences.length())
  for sequence in self.sequences {
    output.write_char(
      sequence.aligned_sequence.unsafe_get(column).unsafe_to_char(),
    )
  }
  Some(output.to_string())
}

///|
/// Map a zero-based ungapped sequence position to an alignment column.
pub fn AlignStockholmAlignment::sequence_position_to_column(
  self : AlignStockholmAlignment,
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
/// Map a retained column to a zero-based ungapped sequence position.
pub fn AlignStockholmAlignment::column_to_sequence_position(
  self : AlignStockholmAlignment,
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
pub fn AlignStockholmAlignment::map_position(
  self : AlignStockholmAlignment,
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
pub fn AlignStockholmAlignment::aligned_pairs(
  self : AlignStockholmAlignment,
  first_row : Int,
  second_row : Int,
) -> Array[(Int?, Int?)] raise AlignStockholmError {
  align_stockholm_validate_row(self, first_row)
  align_stockholm_validate_row(self, second_row)
  let pairs : Array[(Int?, Int?)] = []
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
    pairs.push(
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
  pairs
}

///|
/// Return a compact Biopython-style coordinate path for all rows.
pub fn AlignStockholmAlignment::coordinate_path(
  self : AlignStockholmAlignment,
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
      align_stockholm_movement_changes(self, column, column + 1)
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
/// Slice retained alignment columns and all GC/GR annotations together.
pub fn AlignStockholmAlignment::slice_columns(
  self : AlignStockholmAlignment,
  start : Int,
  end : Int,
) -> AlignStockholmAlignment raise AlignStockholmError {
  if start < 0 || end <= start || end > self.alignment_length() {
    raise AlignStockholmError("Invalid Stockholm column slice")
  }
  let sequences : Array[AlignStockholmSequence] = []
  for sequence in self.sequences {
    let aligned = sequence.aligned_sequence[start:end].to_owned()
    let letter_annotations : Array[AlignStockholmLetterAnnotation] = []
    for annotation in sequence.letter_annotations {
      let aligned_value = annotation.aligned_value[start:end].to_owned()
      let value = align_stockholm_residue_annotation_value(
        aligned, aligned_value,
      )
      letter_annotations.push(AlignStockholmLetterAnnotation::{
        code: annotation.code,
        feature: annotation.feature,
        value,
        aligned_value,
      })
    }
    sequences.push(AlignStockholmSequence::{
      id: sequence.id,
      description: sequence.description,
      sequence: align_stockholm_remove_gaps(aligned),
      aligned_sequence: aligned,
      annotations: sequence.annotations,
      dbxrefs: sequence.dbxrefs,
      letter_annotations,
    })
  }
  let column_annotations : Array[AlignStockholmColumnAnnotation] = []
  for annotation in self.column_annotations {
    column_annotations.push(AlignStockholmColumnAnnotation::{
      code: annotation.code,
      feature: annotation.feature,
      value: annotation.value[start:end].to_owned(),
    })
  }
  AlignStockholmAlignment::create(
    sequences,
    annotations=self.annotations,
    column_annotations~,
    references=self.references,
    database_references=self.database_references,
    nested_domains=self.nested_domains,
    operations=self.operations[start:end].to_owned(),
    source_columns=end - start,
  )
}

///|
/// Compute statistics for one pair of rows.
pub fn AlignStockholmAlignment::pair_counts(
  self : AlignStockholmAlignment,
  first_row : Int,
  second_row : Int,
) -> AlignStockholmCounts raise AlignStockholmError {
  align_stockholm_validate_row(self, first_row)
  align_stockholm_validate_row(self, second_row)
  align_stockholm_count_pair(self, first_row, second_row)
}

///|
/// Aggregate statistics across every unordered row pair.
pub fn AlignStockholmAlignment::counts(
  self : AlignStockholmAlignment,
) -> AlignStockholmCounts {
  let mut result = AlignStockholmCounts::{
    pairs: 0,
    aligned: 0,
    identities: 0,
    mismatches: 0,
    gap_columns: 0,
    double_gap_columns: 0,
    gap_opens: 0,
  }
  for first = 0; first < self.sequences.length(); first = first + 1 {
    for second = first + 1
        second < self.sequences.length()
        second = second + 1 {
      let counts = align_stockholm_count_pair(self, first, second)
      result = AlignStockholmCounts::{
        pairs: result.pairs + counts.pairs,
        aligned: result.aligned + counts.aligned,
        identities: result.identities + counts.identities,
        mismatches: result.mismatches + counts.mismatches,
        gap_columns: result.gap_columns + counts.gap_columns,
        double_gap_columns: result.double_gap_columns +
        counts.double_gap_columns,
        gap_opens: result.gap_opens + counts.gap_opens,
      }
    }
  }
  result
}

///|
/// Return identities divided by aligned residue pairs.
pub fn AlignStockholmCounts::identity(self : AlignStockholmCounts) -> Double {
  if self.aligned == 0 {
    0.0
  } else {
    self.identities.to_double() / self.aligned.to_double()
  }
}

///|
/// Calculate a majority consensus; gaps do not vote.
pub fn AlignStockholmAlignment::consensus(
  self : AlignStockholmAlignment,
  minimum_fraction? : Double = 0.5,
  ambiguous? : String = "X",
) -> String raise AlignStockholmError {
  if minimum_fraction < 0.0 || minimum_fraction > 1.0 {
    raise AlignStockholmError(
      "Stockholm consensus minimum fraction must be between 0 and 1",
    )
  }
  if ambiguous.length() != 1 {
    raise AlignStockholmError(
      "Stockholm consensus ambiguous symbol must be one character",
    )
  }
  let output = StringBuilder::new(size_hint=self.alignment_length())
  for column = 0; column < self.alignment_length(); column = column + 1 {
    let symbols : Array[Int] = []
    let counts : Array[Int] = []
    let mut residues = 0
    for sequence in self.sequences {
      let code = sequence.aligned_sequence.unsafe_get(column).to_int()
      if code != '-'.to_int() {
        residues = residues + 1
        let normalized = align_stockholm_upper_code(code)
        let mut found = -1
        for index = 0; index < symbols.length(); index = index + 1 {
          if symbols[index] == normalized {
            found = index
            break
          }
        }
        if found < 0 {
          symbols.push(normalized)
          counts.push(1)
        } else {
          counts[found] = counts[found] + 1
        }
      }
    }
    if residues == 0 {
      output.write_string(ambiguous)
      continue
    }
    let mut best = 0
    for index = 1; index < counts.length(); index = index + 1 {
      if counts[index] > counts[best] ||
        (counts[index] == counts[best] && symbols[index] < symbols[best]) {
        best = index
      }
    }
    if counts[best].to_double() / residues.to_double() >= minimum_fraction {
      output.write_char(symbols[best].unsafe_to_char())
    } else {
      output.write_string(ambiguous)
    }
  }
  output.to_string()
}

///|
/// Return per-column non-gap occupancy.
pub fn AlignStockholmAlignment::occupancy(
  self : AlignStockholmAlignment,
) -> Array[Double] {
  let values : Array[Double] = []
  for column = 0; column < self.alignment_length(); column = column + 1 {
    let mut residues = 0
    for sequence in self.sequences {
      if sequence.aligned_sequence.unsafe_get(column).to_int() != '-'.to_int() {
        residues = residues + 1
      }
    }
    values.push(
      if self.sequences.length() == 0 {
        0.0
      } else {
        residues.to_double() / self.sequences.length().to_double()
      },
    )
  }
  values
}

///|
/// Return a compact alignment summary.
pub fn AlignStockholmAlignment::summary(
  self : AlignStockholmAlignment,
) -> String {
  "AlignStockholmAlignment(version=" +
  self.version +
  ", sequences=" +
  self.num_sequences().to_string() +
  ", columns=" +
  self.alignment_length().to_string() +
  ", source_columns=" +
  self.source_columns.to_string() +
  ", removed_all_gap=" +
  self.removed_all_gap_columns.to_string() +
  ")"
}

///|
/// Official Biopython HAT fixture used by the example and black-box tests.
pub fn align_stockholm_example_text() -> String {
  "# STOCKHOLM 1.0\n" +
  "#=GF ID   HAT\n" +
  "#=GF AC   PF02184.18\n" +
  "#=GF DE   HAT (Half-A-TPR) repeat\n" +
  "#=GF AU   SMART;\n" +
  "#=GF SE   Alignment kindly provided by SMART\n" +
  "#=GF GA   21.00 21.00;\n" +
  "#=GF TC   21.00 21.00;\n" +
  "#=GF NC   20.90 20.90;\n" +
  "#=GF BM   hmmbuild HMM.ann SEED.ann\n" +
  "#=GF SM   hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq\n" +
  "#=GF TP   Repeat\n" +
  "#=GF CL   CL0020\n" +
  "#=GF RN   [1]\n" +
  "#=GF RM   9478129\n" +
  "#=GF RT   The HAT helix, a repetitive motif implicated in RNA processing.\n" +
  "#=GF RA   Preker PJ, Keller W;\n" +
  "#=GF RL   Trends Biochem Sci 1998;23:15-16.\n" +
  "#=GF DR   INTERPRO; IPR003107;\n" +
  "#=GF DR   SMART; HAT;\n" +
  "#=GF DR   SO; 0001068; polypeptide_repeat;\n" +
  "#=GF CC   The HAT (Half A TPR) repeat is found in several RNA processing\n" +
  "#=GF CC   proteins [1].\n" +
  "#=GF SQ   3\n" +
  "#=GS CRN_DROME/191-222     AC P17886.2\n" +
  "#=GS CLF1_SCHPO/185-216    AC P87312.1\n" +
  "#=GS CLF1_SCHPO/185-216    DR PDB; 3JB9 R; 185-216;\n" +
  "#=GS O16376_CAEEL/201-233  AC O16376.2\n" +
  "CRN_DROME/191-222                KEIDRAREIYERFVYVH.PDVKNWIKFARFEES\n" +
  "CLF1_SCHPO/185-216               HENERARGIYERFVVVH.PEVTNWLRWARFEEE\n" +
  "#=GR CLF1_SCHPO/185-216    SS    --HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH\n" +
  "O16376_CAEEL/201-233             KEIDRARSVYQRFLHVHGINVQNWIKYAKFEER\n" +
  "#=GC SS_cons                     --HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH\n" +
  "#=GC seq_cons                    KEIDRARuIYERFVaVH.P-VpNWIKaARFEEc\n" +
  "//\n"
}

///|
fn align_stockholm_parse_record(
  lines : Array[String],
  start : Int,
  end : Int,
) -> AlignStockholmAlignment raise AlignStockholmError {
  let ids : Array[String] = []
  let rows : Array[String] = []
  let gf : Array[AlignStockholmRawAnnotation] = []
  let gs : Array[AlignStockholmRawSequenceAnnotation] = []
  let gc : Array[AlignStockholmRawAnnotation] = []
  let gr : Array[AlignStockholmRawSequenceAnnotation] = []
  let mut last_sequence_id = ""
  for line_index = start; line_index < end; line_index = line_index + 1 {
    let line = lines[line_index].trim().to_owned()
    if line.length() == 0 {
      continue
    }
    if align_stockholm_starts_with(line, "#=GF") {
      let rest = align_stockholm_markup_rest(line, "#=GF", line_index)
      let fields = align_stockholm_take_fields(rest, 1)
      if fields.0.length() != 1 || fields.1.length() == 0 {
        raise AlignStockholmError(
          "Malformed GF annotation at line " + (line_index + 1).to_string(),
        )
      }
      gf.push(AlignStockholmRawAnnotation::{
        code: fields.0[0],
        value: fields.1,
      })
      continue
    }
    if align_stockholm_starts_with(line, "#=GS") {
      let rest = align_stockholm_markup_rest(line, "#=GS", line_index)
      let fields = align_stockholm_take_fields(rest, 2)
      if fields.0.length() != 2 {
        raise AlignStockholmError(
          "Malformed GS annotation at line " + (line_index + 1).to_string(),
        )
      }
      gs.push(AlignStockholmRawSequenceAnnotation::{
        sequence_id: fields.0[0],
        code: fields.0[1],
        value: fields.1,
      })
      continue
    }
    if align_stockholm_starts_with(line, "#=GC") {
      let rest = align_stockholm_markup_rest(line, "#=GC", line_index)
      let fields = align_stockholm_take_fields(rest, 1)
      if fields.0.length() != 1 || fields.1.length() == 0 {
        raise AlignStockholmError(
          "Malformed GC annotation at line " + (line_index + 1).to_string(),
        )
      }
      align_stockholm_append_raw(gc, fields.0[0], fields.1)
      continue
    }
    if align_stockholm_starts_with(line, "#=GR") {
      let rest = align_stockholm_markup_rest(line, "#=GR", line_index)
      let fields = align_stockholm_take_fields(rest, 2)
      if fields.0.length() != 2 || fields.1.length() == 0 {
        raise AlignStockholmError(
          "Malformed GR annotation at line " + (line_index + 1).to_string(),
        )
      }
      if last_sequence_id.length() == 0 || fields.0[0] != last_sequence_id {
        raise AlignStockholmError("GR annotation must follow its sequence row")
      }
      for annotation in gr {
        if annotation.sequence_id == fields.0[0] &&
          annotation.code == fields.0[1] {
          raise AlignStockholmError(
            "Duplicate GR annotation " + fields.0[1] + " for " + fields.0[0],
          )
        }
      }
      gr.push(AlignStockholmRawSequenceAnnotation::{
        sequence_id: fields.0[0],
        code: fields.0[1],
        value: fields.1,
      })
      continue
    }
    if line.unsafe_get(0).to_int() == '#'.to_int() {
      continue
    }
    let tokens = align_stockholm_split_whitespace(line)
    if tokens.length() != 2 {
      raise AlignStockholmError(
        "Could not split Stockholm sequence line " +
        (line_index + 1).to_string() +
        " into identifier and aligned sequence",
      )
    }
    align_stockholm_validate_id(tokens[0])
    for existing in ids {
      if existing == tokens[0] {
        raise AlignStockholmError(
          "Duplicate Stockholm sequence identifier " + tokens[0],
        )
      }
    }
    if rows.length() > 0 && tokens[1].length() != rows[0].length() {
      raise AlignStockholmError(
        "Aligned sequence " +
        tokens[0] +
        " consists of " +
        tokens[1].length().to_string() +
        " letters, expected " +
        rows[0].length().to_string() +
        " letters",
      )
    }
    if tokens[1].length() == 0 {
      raise AlignStockholmError("Stockholm sequence row must not be empty")
    }
    align_stockholm_validate_printed_row(tokens[1])
    ids.push(tokens[0])
    rows.push(tokens[1])
    last_sequence_id = tokens[0]
  }
  if ids.length() == 0 {
    raise AlignStockholmError(
      "Stockholm alignment must contain at least one sequence",
    )
  }
  let prepared = align_stockholm_prepare_rows(rows)
  if prepared.0[0].length() == 0 {
    raise AlignStockholmError(
      "Stockholm alignment contains no non-gap alignment column",
    )
  }
  let annotations = align_stockholm_build_gf(gf, ids.length())
  let references = align_stockholm_build_references(gf)
  let database_references = align_stockholm_build_database_references(gf)
  let nested_domains = align_stockholm_build_nested_domains(gf)
  let column_annotations = align_stockholm_build_gc(
    gc,
    rows[0].length(),
    prepared.2,
  )
  let sequences : Array[AlignStockholmSequence] = []
  for row_index = 0; row_index < ids.length(); row_index = row_index + 1 {
    let sequence_annotations : Array[AlignStockholmAnnotation] = []
    let dbxrefs : Array[String] = []
    let mut description = ""
    for annotation in gs {
      if annotation.sequence_id != ids[row_index] {
        continue
      }
      if annotation.code == "DR" {
        dbxrefs.push(annotation.value)
        continue
      }
      if annotation.code == "DE" {
        if description.length() > 0 {
          raise AlignStockholmError(
            "Duplicate GS DE annotation for " + ids[row_index],
          )
        }
        description = annotation.value
        continue
      }
      for existing in sequence_annotations {
        if existing.code == annotation.code {
          raise AlignStockholmError(
            "Duplicate GS annotation " +
            annotation.code +
            " for " +
            ids[row_index],
          )
        }
      }
      sequence_annotations.push(AlignStockholmAnnotation::{
        code: annotation.code,
        feature: align_stockholm_gs_feature(annotation.code),
        value: annotation.value,
      })
    }
    let letter_annotations : Array[AlignStockholmLetterAnnotation] = []
    for annotation in gr {
      if annotation.sequence_id != ids[row_index] {
        continue
      }
      if annotation.value.length() != rows[0].length() {
        raise AlignStockholmError(
          "GR " +
          annotation.code +
          " length is " +
          annotation.value.length().to_string() +
          ", expected " +
          rows[0].length().to_string(),
        )
      }
      let aligned_value = align_stockholm_remove_columns(
        annotation.value,
        prepared.2,
      )
      let value = if annotation.code == "CSA" {
        align_stockholm_remove_character(annotation.value, '-'.to_int())
      } else {
        align_stockholm_remove_character(annotation.value, '.'.to_int())
      }
      let ungapped = align_stockholm_remove_gaps(prepared.0[row_index])
      if value.length() != ungapped.length() {
        raise AlignStockholmError(
          "GR " +
          annotation.code +
          " does not contain one value per residue for " +
          ids[row_index],
        )
      }
      letter_annotations.push(AlignStockholmLetterAnnotation::{
        code: annotation.code,
        feature: align_stockholm_gr_feature(annotation.code),
        value,
        aligned_value,
      })
    }
    sequences.push(AlignStockholmSequence::{
      id: ids[row_index],
      description,
      sequence: align_stockholm_remove_gaps(prepared.0[row_index]),
      aligned_sequence: prepared.0[row_index],
      annotations: sequence_annotations,
      dbxrefs,
      letter_annotations,
    })
  }
  for annotation in gs {
    if align_stockholm_find_string(ids, annotation.sequence_id) < 0 {
      raise AlignStockholmError(
        "Failed to find GS sequence " + annotation.sequence_id,
      )
    }
  }
  for annotation in gr {
    if align_stockholm_find_string(ids, annotation.sequence_id) < 0 {
      raise AlignStockholmError(
        "Failed to find GR sequence " + annotation.sequence_id,
      )
    }
  }
  AlignStockholmAlignment::create(
    sequences,
    annotations~,
    column_annotations~,
    references~,
    database_references~,
    nested_domains~,
    operations=prepared.1,
    source_columns=rows[0].length(),
    removed_all_gap_columns=align_stockholm_count_true(prepared.2),
  )
}

///|
fn align_stockholm_build_gf(
  raw : Array[AlignStockholmRawAnnotation],
  rows : Int,
) -> Array[AlignStockholmAnnotation] raise AlignStockholmError {
  let annotations : Array[AlignStockholmAnnotation] = []
  let codes = [
    "ID", "AC", "DE", "AU", "SE", "SS", "GA", "TC", "NC", "BM", "SM", "TP", "PI",
    "CL", "WK", "CB", "**", "CC",
  ]
  for code in codes {
    let values : Array[String] = []
    for annotation in raw {
      if annotation.code == code {
        values.push(annotation.value)
      }
    }
    if values.length() == 0 {
      continue
    }
    if code == "AU" {
      for value in values {
        annotations.push(AlignStockholmAnnotation::{
          code,
          feature: align_stockholm_gf_feature(code),
          value,
        })
      }
    } else if code == "WK" {
      let merged = align_stockholm_merge_wikipedia(values)
      for value in merged {
        annotations.push(AlignStockholmAnnotation::{
          code,
          feature: align_stockholm_gf_feature(code),
          value,
        })
      }
    } else if code == "SM" || code == "CC" || code == "**" {
      annotations.push(AlignStockholmAnnotation::{
        code,
        feature: align_stockholm_gf_feature(code),
        value: align_stockholm_join(values, " "),
      })
    } else {
      if values.length() != 1 {
        raise AlignStockholmError("GF " + code + " must occur at most once")
      }
      annotations.push(AlignStockholmAnnotation::{
        code,
        feature: align_stockholm_gf_feature(code),
        value: values[0],
      })
    }
  }
  let sq_values : Array[String] = []
  for annotation in raw {
    if annotation.code == "SQ" {
      sq_values.push(annotation.value)
    }
  }
  if sq_values.length() > 1 {
    raise AlignStockholmError("GF SQ must occur at most once")
  }
  if sq_values.length() == 1 {
    let declared = align_stockholm_parse_positive_int(
      sq_values[0],
      "Stockholm GF SQ",
    )
    if declared != rows {
      raise AlignStockholmError(
        "Inconsistent number of sequences in Stockholm alignment",
      )
    }
  }
  annotations
}

///|
fn align_stockholm_build_references(
  raw : Array[AlignStockholmRawAnnotation],
) -> Array[AlignStockholmReference] raise AlignStockholmError {
  let builders : Array[AlignStockholmReferenceBuilder] = []
  let pending_comments : Array[String] = []
  for annotation in raw {
    if annotation.code == "RC" {
      pending_comments.push(annotation.value)
    } else if annotation.code == "RN" {
      let value = annotation.value
      if value.length() < 3 ||
        value.unsafe_get(0).to_int() != '['.to_int() ||
        value.unsafe_get(value.length() - 1).to_int() != ']'.to_int() {
        raise AlignStockholmError("Malformed Stockholm GF RN annotation")
      }
      let number = align_stockholm_parse_positive_int(
        value[1:value.length() - 1].to_owned(),
        "Stockholm reference number",
      )
      builders.push(AlignStockholmReferenceBuilder::{
        number,
        medline: "",
        title: "",
        author: "",
        location: "",
        comment: align_stockholm_join(pending_comments, " "),
      })
      pending_comments.clear()
    } else if annotation.code == "RM" ||
      annotation.code == "RT" ||
      annotation.code == "RA" ||
      annotation.code == "RL" {
      if builders.length() == 0 {
        raise AlignStockholmError(
          "Stockholm GF " + annotation.code + " appears before GF RN",
        )
      }
      let last = builders.length() - 1
      if annotation.code == "RM" {
        if builders[last].medline.length() > 0 {
          raise AlignStockholmError("Duplicate Stockholm GF RM annotation")
        }
        builders[last].medline = annotation.value
      } else if annotation.code == "RT" {
        builders[last].title = align_stockholm_append_text(
          builders[last].title,
          annotation.value,
        )
      } else if annotation.code == "RA" {
        builders[last].author = align_stockholm_append_text(
          builders[last].author,
          annotation.value,
        )
      } else {
        builders[last].location = align_stockholm_append_text(
          builders[last].location,
          annotation.value,
        )
      }
    }
  }
  if pending_comments.length() > 0 {
    raise AlignStockholmError(
      "Stockholm GF RC annotation is not followed by GF RN",
    )
  }
  let references : Array[AlignStockholmReference] = []
  for builder in builders {
    references.push(AlignStockholmReference::{
      number: builder.number,
      medline: builder.medline,
      title: builder.title,
      author: builder.author,
      location: builder.location,
      comment: builder.comment,
    })
  }
  references
}

///|
fn align_stockholm_build_database_references(
  raw : Array[AlignStockholmRawAnnotation],
) -> Array[AlignStockholmDatabaseReference] raise AlignStockholmError {
  let builders : Array[AlignStockholmDatabaseReferenceBuilder] = []
  for annotation in raw {
    if annotation.code == "DR" {
      builders.push(AlignStockholmDatabaseReferenceBuilder::{
        reference: annotation.value,
        comment: "",
      })
    } else if annotation.code == "DC" {
      if builders.length() == 0 {
        raise AlignStockholmError(
          "Stockholm GF DC appears before a database reference",
        )
      }
      let last = builders.length() - 1
      if builders[last].comment.length() > 0 {
        raise AlignStockholmError("Duplicate Stockholm GF DC annotation")
      }
      builders[last].comment = annotation.value
    }
  }
  let references : Array[AlignStockholmDatabaseReference] = []
  for builder in builders {
    references.push(AlignStockholmDatabaseReference::{
      reference: builder.reference,
      comment: builder.comment,
    })
  }
  references
}

///|
fn align_stockholm_build_nested_domains(
  raw : Array[AlignStockholmRawAnnotation],
) -> Array[AlignStockholmNestedDomain] raise AlignStockholmError {
  let builders : Array[AlignStockholmNestedDomainBuilder] = []
  for annotation in raw {
    if annotation.code == "NE" {
      builders.push(AlignStockholmNestedDomainBuilder::{
        accession: annotation.value,
        location: "",
      })
    } else if annotation.code == "NL" {
      if builders.length() == 0 {
        raise AlignStockholmError(
          "Stockholm GF NL appears before a nested domain",
        )
      }
      let last = builders.length() - 1
      if builders[last].location.length() > 0 {
        raise AlignStockholmError("Duplicate Stockholm GF NL annotation")
      }
      builders[last].location = annotation.value
    }
  }
  let domains : Array[AlignStockholmNestedDomain] = []
  for builder in builders {
    domains.push(AlignStockholmNestedDomain::{
      accession: builder.accession,
      location: builder.location,
    })
  }
  domains
}

///|
fn align_stockholm_build_gc(
  raw : Array[AlignStockholmRawAnnotation],
  source_width : Int,
  removed : Array[Bool],
) -> Array[AlignStockholmColumnAnnotation] raise AlignStockholmError {
  let annotations : Array[AlignStockholmColumnAnnotation] = []
  for annotation in raw {
    if annotation.value.length() != source_width {
      raise AlignStockholmError(
        "GC " +
        annotation.code +
        " length is " +
        annotation.value.length().to_string() +
        ", expected " +
        source_width.to_string(),
      )
    }
    annotations.push(AlignStockholmColumnAnnotation::{
      code: annotation.code,
      feature: align_stockholm_gc_feature(annotation.code),
      value: align_stockholm_remove_columns(annotation.value, removed),
    })
  }
  annotations
}

///|
fn align_stockholm_prepare_rows(
  rows : Array[String],
) -> (Array[String], String, Array[Bool]) raise AlignStockholmError {
  if rows.length() == 0 {
    raise AlignStockholmError(
      "Stockholm alignment must contain at least one sequence",
    )
  }
  let width = rows[0].length()
  let operations = Array::make(width, 'M'.to_int())
  let all_gap = Array::make(width, true)
  for row in rows {
    if row.length() != width {
      raise AlignStockholmError("Stockholm aligned rows must have equal widths")
    }
    align_stockholm_validate_printed_row(row)
    for column = 0; column < width; column = column + 1 {
      let code = row.unsafe_get(column).to_int()
      if code == '-'.to_int() {
        if operations[column] == 'I'.to_int() {
          raise AlignStockholmError(
            "Stockholm column mixes insertion and deletion gap symbols",
          )
        }
        operations[column] = 'D'.to_int()
      } else if code == '.'.to_int() {
        if operations[column] == 'D'.to_int() {
          raise AlignStockholmError(
            "Stockholm column mixes insertion and deletion gap symbols",
          )
        }
        operations[column] = 'I'.to_int()
      } else {
        all_gap[column] = false
      }
    }
  }
  let normalized : Array[String] = []
  for row in rows {
    let output = StringBuilder::new(size_hint=width)
    for column = 0; column < width; column = column + 1 {
      if !all_gap[column] {
        let code = row.unsafe_get(column).to_int()
        output.write_char(
          (if code == '.'.to_int() { '-'.to_int() } else { code }).unsafe_to_char(),
        )
      }
    }
    normalized.push(output.to_string())
  }
  let operation_output = StringBuilder::new(size_hint=width)
  for column = 0; column < width; column = column + 1 {
    if !all_gap[column] {
      operation_output.write_char(operations[column].unsafe_to_char())
    }
  }
  (normalized, operation_output.to_string(), all_gap)
}

///|
fn align_stockholm_validate_alignment(
  alignment : AlignStockholmAlignment,
) -> Unit raise AlignStockholmError {
  if alignment.version != "1.0" {
    raise AlignStockholmError("Only Stockholm version 1.0 is supported")
  }
  if alignment.sequences.length() == 0 {
    raise AlignStockholmError(
      "Stockholm alignment must contain at least one sequence",
    )
  }
  let width = alignment.sequences[0].aligned_sequence.length()
  if width == 0 {
    raise AlignStockholmError(
      "Stockholm alignment must contain at least one column",
    )
  }
  let ids : Array[String] = []
  for sequence in alignment.sequences {
    align_stockholm_validate_id(sequence.id)
    if align_stockholm_find_string(ids, sequence.id) >= 0 {
      raise AlignStockholmError(
        "Duplicate Stockholm sequence identifier " + sequence.id,
      )
    }
    ids.push(sequence.id)
    if sequence.aligned_sequence.length() != width {
      raise AlignStockholmError("Stockholm aligned rows must have equal widths")
    }
    if sequence.sequence !=
      align_stockholm_remove_gaps(sequence.aligned_sequence) {
      raise AlignStockholmError(
        "Stockholm ungapped sequence does not match aligned row",
      )
    }
    align_stockholm_validate_sequence_annotations(sequence)
  }
  if alignment.operations.length() != width {
    raise AlignStockholmError(
      "Stockholm operations length must match alignment width",
    )
  }
  for column = 0; column < width; column = column + 1 {
    let operation = alignment.operations.unsafe_get(column).to_int()
    if operation != 'M'.to_int() &&
      operation != 'D'.to_int() &&
      operation != 'I'.to_int() {
      raise AlignStockholmError(
        "Stockholm operations may contain only M, D, and I",
      )
    }
    let mut all_gap = true
    let mut has_gap = false
    for sequence in alignment.sequences {
      if sequence.aligned_sequence.unsafe_get(column).to_int() == '-'.to_int() {
        has_gap = true
      } else {
        all_gap = false
      }
    }
    if all_gap {
      raise AlignStockholmError(
        "Stockholm alignment must not retain all-gap columns",
      )
    }
    if operation == 'M'.to_int() && has_gap {
      raise AlignStockholmError(
        "Stockholm match operation cannot contain a gap",
      )
    }
  }
  for annotation in alignment.column_annotations {
    if annotation.code.length() == 0 {
      raise AlignStockholmError("Stockholm GC code must not be empty")
    }
    if annotation.value.length() != width {
      raise AlignStockholmError(
        "Stockholm GC annotation length must match alignment width",
      )
    }
  }
  if alignment.removed_all_gap_columns < 0 ||
    alignment.source_columns != width + alignment.removed_all_gap_columns {
    raise AlignStockholmError("Invalid Stockholm source-column metadata")
  }
}

///|
fn align_stockholm_validate_sequence_annotations(
  sequence : AlignStockholmSequence,
) -> Unit raise AlignStockholmError {
  for annotation in sequence.annotations {
    if annotation.code.length() == 0 {
      raise AlignStockholmError("Stockholm GS code must not be empty")
    }
  }
  for annotation in sequence.letter_annotations {
    if annotation.code.length() == 0 {
      raise AlignStockholmError("Stockholm GR code must not be empty")
    }
    if annotation.value.length() != sequence.sequence.length() {
      raise AlignStockholmError(
        "Stockholm GR residue annotation length must match sequence length",
      )
    }
    if annotation.aligned_value.length() != sequence.aligned_sequence.length() {
      raise AlignStockholmError(
        "Stockholm aligned GR annotation length must match alignment width",
      )
    }
  }
}

///|
fn align_stockholm_count_pair(
  alignment : AlignStockholmAlignment,
  first_row : Int,
  second_row : Int,
) -> AlignStockholmCounts {
  let mut aligned = 0
  let mut identities = 0
  let mut mismatches = 0
  let mut gap_columns = 0
  let mut double_gap_columns = 0
  let mut gap_opens = 0
  let mut previous_gap_state = 0
  for column = 0; column < alignment.alignment_length(); column = column + 1 {
    let first = alignment.sequences[first_row].aligned_sequence
      .unsafe_get(column)
      .to_int()
    let second = alignment.sequences[second_row].aligned_sequence
      .unsafe_get(column)
      .to_int()
    let first_gap = first == '-'.to_int()
    let second_gap = second == '-'.to_int()
    let gap_state = if first_gap && second_gap {
      3
    } else if first_gap {
      1
    } else if second_gap {
      2
    } else {
      0
    }
    if gap_state == 0 {
      aligned = aligned + 1
      if align_stockholm_upper_code(first) == align_stockholm_upper_code(second) {
        identities = identities + 1
      } else {
        mismatches = mismatches + 1
      }
    } else {
      if gap_state == 3 {
        double_gap_columns = double_gap_columns + 1
      } else {
        gap_columns = gap_columns + 1
      }
      if gap_state != previous_gap_state {
        gap_opens = gap_opens + 1
      }
    }
    previous_gap_state = gap_state
  }
  AlignStockholmCounts::{
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
fn align_stockholm_movement_changes(
  alignment : AlignStockholmAlignment,
  first_column : Int,
  second_column : Int,
) -> Bool {
  for sequence in alignment.sequences {
    let first_gap = sequence.aligned_sequence.unsafe_get(first_column).to_int() ==
      '-'.to_int()
    let second_gap = sequence.aligned_sequence
      .unsafe_get(second_column)
      .to_int() ==
      '-'.to_int()
    if first_gap != second_gap {
      return true
    }
  }
  false
}

///|
fn align_stockholm_gf_feature(code : String) -> String {
  match code {
    "ID" => "identifier"
    "AC" => "accession"
    "DE" => "definition"
    "AU" => "author"
    "SE" => "source of seed"
    "SS" => "source of structure"
    "GA" => "gathering method"
    "TC" => "trusted cutoff"
    "NC" => "noise cutoff"
    "BM" => "build method"
    "SM" => "search method"
    "TP" => "type"
    "PI" => "previous identifier"
    "CC" => "comment"
    "CL" => "clan"
    "WK" => "wikipedia"
    "CB" => "calibration method"
    "**" => "**"
    _ => code
  }
}

///|
fn align_stockholm_gs_feature(code : String) -> String {
  match code {
    "AC" => "accession"
    "OS" => "organism"
    "OC" => "organism classification"
    "LO" => "look"
    _ => code
  }
}

///|
fn align_stockholm_gr_feature(code : String) -> String {
  match code {
    "SS" => "secondary structure"
    "PP" => "posterior probability"
    "CSA" => "Catalytic Site Atlas"
    "SA" => "surface accessibility"
    "TM" => "transmembrane"
    "LI" => "ligand binding"
    "AS" => "active site"
    "pAS" => "active site - Pfam predicted"
    "sAS" => "active site - from SwissProt"
    "IN" => "intron"
    _ => code
  }
}

///|
fn align_stockholm_gc_feature(code : String) -> String {
  match code {
    "RF" => "reference coordinate annotation"
    "seq_cons" => "consensus sequence"
    "scorecons" => "consensus score"
    "scorecons_70" => "consensus score 70"
    "scorecons_80" => "consensus score 80"
    "scorecons_90" => "consensus score 90"
    "MM" => "model mask"
    "SS_cons" => "consensus secondary structure"
    "PP_cons" => "consensus posterior probability"
    "CSA_cons" => "consensus Catalytic Site Atlas"
    "SA_cons" => "consensus surface accessibility"
    "TM_cons" => "consensus transmembrane"
    "LI_cons" => "consensus ligand binding"
    "AS_cons" => "consensus active site"
    "pAS_cons" => "consensus active site - Pfam predicted"
    "sAS_cons" => "consensus active site - from SwissProt"
    "IN_cons" => "consensus intron"
    "RNA_elements" => "RNA elements"
    "RNA_structural_element" => "RNA structural element"
    "RNA_structural_elements" => "RNA structural elements"
    "RNA_ligand_AdoCbl" => "RNA ligand AdoCbl"
    "RNA_ligand_AqCbl" => "RNA ligand AqCbl"
    "RNA_ligand_FMN" => "RNA ligand FMN"
    "RNA_ligand_Guanidinium" => "RNA ligand Guanidinium"
    "RNA_ligand_SAM" => "RNA ligand SAM"
    "RNA_ligand_THF_1" => "RNA ligand THF 1"
    "RNA_ligand_THF_2" => "RNA ligand THF 2"
    "RNA_ligand_TPP" => "RNA ligand TPP"
    "RNA_ligand_preQ1" => "RNA ligand preQ1"
    "RNA_motif_k_turn" => "RNA motif k turn"
    "Repeat_unit" => "Repeat unit"
    "2L3J_B_SS" => "2L3J B SS"
    "CORE" => "CORE"
    "PK" => "PK"
    "PK_SS" => "PK SS"
    "cons" => "cons"
    _ => code
  }
}

///|
fn align_stockholm_is_known_gf(code : String) -> Bool {
  code == "ID" ||
  code == "AC" ||
  code == "DE" ||
  code == "AU" ||
  code == "SE" ||
  code == "SS" ||
  code == "GA" ||
  code == "TC" ||
  code == "NC" ||
  code == "BM" ||
  code == "SM" ||
  code == "TP" ||
  code == "PI" ||
  code == "CC" ||
  code == "CL" ||
  code == "WK" ||
  code == "CB" ||
  code == "**"
}

///|
fn align_stockholm_infer_operations(
  sequences : Array[AlignStockholmSequence],
) -> String {
  let width = if sequences.length() == 0 {
    0
  } else {
    sequences[0].aligned_sequence.length()
  }
  let output = StringBuilder::new(size_hint=width)
  for column = 0; column < width; column = column + 1 {
    let mut has_gap = false
    for sequence in sequences {
      if sequence.aligned_sequence.unsafe_get(column).to_int() == '-'.to_int() {
        has_gap = true
        break
      }
    }
    output.write_char(if has_gap { 'D' } else { 'M' })
  }
  output.to_string()
}

///|
fn align_stockholm_render_row(row : String, operations : String) -> String {
  let output = StringBuilder::new(size_hint=row.length())
  for column = 0; column < row.length(); column = column + 1 {
    let code = row.unsafe_get(column).to_int()
    if code == '-'.to_int() &&
      operations.unsafe_get(column).to_int() == 'I'.to_int() {
      output.write_char('.')
    } else {
      output.write_char(code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn align_stockholm_render_letter_annotation(
  sequence : AlignStockholmSequence,
  value : String,
) -> String {
  let output = StringBuilder::new(size_hint=sequence.aligned_sequence.length())
  let mut residue = 0
  for column = 0
      column < sequence.aligned_sequence.length()
      column = column + 1 {
    if sequence.aligned_sequence.unsafe_get(column).to_int() == '-'.to_int() {
      output.write_char('.')
    } else {
      output.write_char(value.unsafe_get(residue).unsafe_to_char())
      residue = residue + 1
    }
  }
  output.to_string()
}

///|
fn align_stockholm_residue_annotation_value(
  aligned_sequence : String,
  aligned_value : String,
) -> String {
  let output = StringBuilder::new()
  for column = 0; column < aligned_sequence.length(); column = column + 1 {
    if aligned_sequence.unsafe_get(column).to_int() != '-'.to_int() {
      output.write_char(aligned_value.unsafe_get(column).unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn align_stockholm_remove_columns(
  value : String,
  removed : Array[Bool],
) -> String {
  let output = StringBuilder::new(size_hint=value.length())
  for index = 0; index < value.length(); index = index + 1 {
    if !removed[index] {
      output.write_char(value.unsafe_get(index).unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn align_stockholm_count_true(values : Array[Bool]) -> Int {
  let mut count = 0
  for value in values {
    if value {
      count = count + 1
    }
  }
  count
}

///|
fn align_stockholm_remove_gaps(sequence : String) -> String {
  let output = StringBuilder::new(size_hint=sequence.length())
  for index = 0; index < sequence.length(); index = index + 1 {
    let code = sequence.unsafe_get(index).to_int()
    if code != '-'.to_int() && code != '.'.to_int() {
      output.write_char(code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn align_stockholm_remove_character(value : String, target : Int) -> String {
  let output = StringBuilder::new(size_hint=value.length())
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    if code != target {
      output.write_char(code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn align_stockholm_normalize_row(
  row : String,
) -> String raise AlignStockholmError {
  align_stockholm_validate_printed_row(row)
  let output = StringBuilder::new(size_hint=row.length())
  for index = 0; index < row.length(); index = index + 1 {
    let code = row.unsafe_get(index).to_int()
    output.write_char(
      (if code == '.'.to_int() { '-'.to_int() } else { code }).unsafe_to_char(),
    )
  }
  output.to_string()
}

///|
fn align_stockholm_validate_printed_row(
  row : String,
) -> Unit raise AlignStockholmError {
  for index = 0; index < row.length(); index = index + 1 {
    let code = row.unsafe_get(index).to_int()
    if align_stockholm_is_whitespace(code) || code < 33 || code > 126 {
      raise AlignStockholmError(
        "Stockholm sequence rows may contain only printable non-space ASCII",
      )
    }
  }
}

///|
fn align_stockholm_validate_id(id : String) -> Unit raise AlignStockholmError {
  if id.length() == 0 {
    raise AlignStockholmError("Stockholm sequence identifier must not be empty")
  }
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    if align_stockholm_is_whitespace(code) || code < 33 || code > 126 {
      raise AlignStockholmError(
        "Stockholm sequence identifiers may not contain whitespace",
      )
    }
  }
}

///|
fn align_stockholm_validate_annotation_code(
  code : String,
) -> Unit raise AlignStockholmError {
  if code.length() == 0 {
    raise AlignStockholmError("Stockholm annotation code must not be empty")
  }
  for index = 0; index < code.length(); index = index + 1 {
    let value = code.unsafe_get(index).to_int()
    if align_stockholm_is_whitespace(value) || value < 33 || value > 126 {
      raise AlignStockholmError(
        "Stockholm annotation codes may not contain whitespace",
      )
    }
  }
}

///|
fn align_stockholm_validate_row(
  alignment : AlignStockholmAlignment,
  row : Int,
) -> Unit raise AlignStockholmError {
  if row < 0 || row >= alignment.sequences.length() {
    raise AlignStockholmError(
      "Stockholm row index " + row.to_string() + " is out of range",
    )
  }
}

///|
fn align_stockholm_append_raw(
  values : Array[AlignStockholmRawAnnotation],
  code : String,
  value : String,
) -> Unit {
  for index = 0; index < values.length(); index = index + 1 {
    if values[index].code == code {
      values[index] = AlignStockholmRawAnnotation::{
        code,
        value: values[index].value + value,
      }
      return
    }
  }
  values.push(AlignStockholmRawAnnotation::{ code, value })
}

///|
fn align_stockholm_markup_rest(
  line : String,
  prefix : String,
  line_index : Int,
) -> String raise AlignStockholmError {
  if line.length() == prefix.length() ||
    !align_stockholm_is_whitespace(line.unsafe_get(prefix.length()).to_int()) {
    raise AlignStockholmError(
      "Malformed Stockholm markup at line " + (line_index + 1).to_string(),
    )
  }
  line[prefix.length():].to_owned().trim().to_owned()
}

///|
fn align_stockholm_take_fields(
  value : String,
  count : Int,
) -> (Array[String], String) {
  let fields : Array[String] = []
  let mut index = 0
  while fields.length() < count {
    while index < value.length() &&
          align_stockholm_is_whitespace(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index >= value.length() {
      break
    }
    let start = index
    while index < value.length() &&
          !align_stockholm_is_whitespace(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    fields.push(value[start:index].to_owned())
  }
  while index < value.length() &&
        align_stockholm_is_whitespace(value.unsafe_get(index).to_int()) {
    index = index + 1
  }
  (fields, if index < value.length() { value[index:].to_owned() } else { "" })
}

///|
fn align_stockholm_split_whitespace(value : String) -> Array[String] {
  let fields : Array[String] = []
  let mut index = 0
  while index < value.length() {
    while index < value.length() &&
          align_stockholm_is_whitespace(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index >= value.length() {
      break
    }
    let start = index
    while index < value.length() &&
          !align_stockholm_is_whitespace(value.unsafe_get(index).to_int()) {
      index = index + 1
    }
    fields.push(value[start:index].to_owned())
  }
  fields
}

///|
fn align_stockholm_merge_wikipedia(values : Array[String]) -> Array[String] {
  let merged : Array[String] = []
  let mut index = 0
  while index < values.length() {
    let output = StringBuilder::new()
    let mut value = values[index]
    while value.length() > 0 &&
          value.unsafe_get(value.length() - 1).to_int() == '/'.to_int() &&
          index + 1 < values.length() {
      output.write_string(value[0:value.length() - 1].to_owned())
      index = index + 1
      value = values[index]
    }
    output.write_string(value)
    merged.push(output.to_string())
    index = index + 1
  }
  merged
}

///|
fn align_stockholm_wrap_text(prefix : String, text : String) -> String {
  if text.length() == 0 {
    return prefix + "\n"
  }
  let words = align_stockholm_split_whitespace(text)
  if words.length() == 0 {
    return prefix + "\n"
  }
  let output = StringBuilder::new()
  let mut line = prefix
  for word in words {
    let separator = if line.length() == prefix.length() { "" } else { " " }
    if line.length() > prefix.length() && line.length() + 1 + word.length() > 79 {
      output.write_string(line)
      output.write_char('\n')
      line = prefix + word
    } else {
      line = line + separator + word
    }
  }
  output.write_string(line)
  output.write_char('\n')
  output.to_string()
}

///|
fn align_stockholm_parse_positive_int(
  text : String,
  label : String,
) -> Int raise AlignStockholmError {
  if text.length() == 0 {
    raise AlignStockholmError(label + " is empty")
  }
  let mut value = 0
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise AlignStockholmError(label + " must be a positive integer")
    }
    let digit = code - '0'.to_int()
    if value > 214748364 || (value == 214748364 && digit > 7) {
      raise AlignStockholmError(
        label + " is outside the supported integer range",
      )
    }
    value = value * 10 + digit
  }
  if value <= 0 {
    raise AlignStockholmError(label + " must be positive")
  }
  value
}

///|
fn align_stockholm_normalize_lines(text : String) -> Array[String] {
  let lines : Array[String] = []
  for raw in text.split("\n") {
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
fn align_stockholm_find_string(values : Array[String], target : String) -> Int {
  for index = 0; index < values.length(); index = index + 1 {
    if values[index] == target {
      return index
    }
  }
  -1
}

///|
fn align_stockholm_join(values : Array[String], separator : String) -> String {
  let output = StringBuilder::new()
  for index = 0; index < values.length(); index = index + 1 {
    if index > 0 {
      output.write_string(separator)
    }
    output.write_string(values[index])
  }
  output.to_string()
}

///|
fn align_stockholm_append_text(existing : String, value : String) -> String {
  if existing.length() == 0 {
    value
  } else {
    existing + " " + value
  }
}

///|
fn align_stockholm_pad_right(value : String, width : Int) -> String {
  if value.length() >= width {
    return value
  }
  let output = StringBuilder::new(size_hint=width)
  output.write_string(value)
  for _index = value.length(); _index < width; _index = _index + 1 {
    output.write_char(' ')
  }
  output.to_string()
}

///|
fn align_stockholm_starts_with(value : String, prefix : String) -> Bool {
  if prefix.length() > value.length() {
    return false
  }
  for index = 0; index < prefix.length(); index = index + 1 {
    if value.unsafe_get(index).to_int() != prefix.unsafe_get(index).to_int() {
      return false
    }
  }
  true
}

///|
fn align_stockholm_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\n'.to_int() ||
  code == '\r'.to_int()
}

///|
fn align_stockholm_upper_code(code : Int) -> Int {
  if code >= 'a'.to_int() && code <= 'z'.to_int() {
    code - 'a'.to_int() + 'A'.to_int()
  } else {
    code
  }
}

///|
fn align_stockholm_max(first : Int, second : Int) -> Int {
  if first > second {
    first
  } else {
    second
  }
}
