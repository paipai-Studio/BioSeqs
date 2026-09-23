// Biopython-compatible Bio.Align.mauve support.
//
// This module is intentionally separate from mauve.mbt. The older module
// provides MAF-like rearrangement summaries, while this file models genuine
// Mauve/ProgressiveMauve XMFA documents, aligned rows, and coordinate paths.

///|
/// Error raised for malformed XMFA data or invalid coordinate operations.
pub suberror AlignMauveError {
  AlignMauveError(String)
}

///|
/// One non-sequence file-level XMFA metadata entry.
pub struct AlignMauveMetadataEntry {
  key : String
  value : String
} derive(Eq, Debug)
pub extend AlignMauveMetadataEntry with Eq::{equal, not_equal}
pub extend AlignMauveMetadataEntry with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One sequence declared by `#SequenceN*` header fields.
///
/// `number` and `entry` retain XMFA's one-based values. `entry == None`
/// denotes separate-file input without a `SequenceNEntry` line.
pub struct AlignMauveSource {
  number : Int
  file : String
  entry : Int?
  format : String
  identifier : String
} derive(Eq, Debug)
pub extend AlignMauveSource with Eq::{equal, not_equal}
pub extend AlignMauveSource with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One aligned row in a locally collinear XMFA block.
///
/// `start` and `end` are zero-based forward-axis half-open coordinates.
/// Reverse rows have a decreasing `coordinates` path from `end` to `start`.
/// `source_sequence` is always stored in the forward source orientation.
pub struct AlignMauveRow {
  source_number : Int
  identifier : String
  start : Int
  end : Int
  strand : String
  description : String
  aligned_sequence : String
  source_sequence : String
  coordinates : Array[Int]
} derive(Eq, Debug)
pub extend AlignMauveRow with Eq::{equal, not_equal}
pub extend AlignMauveRow with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One Mauve locally collinear block.
pub struct AlignMauveBlock {
  rows : Array[AlignMauveRow]
  width : Int
} derive(Eq, Debug)
pub extend AlignMauveBlock with Eq::{equal, not_equal}
pub extend AlignMauveBlock with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// A complete XMFA document.
pub struct AlignMauveDocument {
  format_version : String
  sources : Array[AlignMauveSource]
  metadata : Array[AlignMauveMetadataEntry]
  blocks : Array[AlignMauveBlock]
} derive(Eq, Debug)
pub extend AlignMauveDocument with Eq::{equal, not_equal}
pub extend AlignMauveDocument with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Pairwise statistics summed over all unordered row pairs in one or more
/// blocks. Columns where both rows contain gaps are reported separately and
/// do not contribute to gap counts.
pub struct AlignMauveCounts {
  pair_count : Int
  columns : Int
  aligned : Int
  identities : Int
  mismatches : Int
  gaps : Int
  gap_opens : Int
  double_gap_columns : Int
} derive(Eq, Debug)
pub extend AlignMauveCounts with Eq::{equal, not_equal}
pub extend AlignMauveCounts with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One interval-index entry for an aligned source segment.
pub struct AlignMauveIndexEntry {
  identifier : String
  start : Int
  end : Int
  block_index : Int
  row_index : Int
  strand : String
} derive(Eq, Debug)
pub extend AlignMauveIndexEntry with Eq::{equal, not_equal}
pub extend AlignMauveIndexEntry with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// In-memory interval index over all non-empty XMFA rows.
pub struct AlignMauveIndex {
  identifiers : Array[String]
  entries : Array[AlignMauveIndexEntry]
} derive(Eq, Debug)
pub extend AlignMauveIndex with Eq::{equal, not_equal}
pub extend AlignMauveIndex with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One exact residue-to-residue coordinate projection.
pub struct AlignMauvePositionMapping {
  block_index : Int
  column : Int
  source_position : Int
  target_position : Int
} derive(Eq, Debug)
pub extend AlignMauvePositionMapping with Eq::{equal, not_equal}
pub extend AlignMauvePositionMapping with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// One contiguous aligned part of a projected source interval.
///
/// Source and target intervals are zero-based forward-axis half-open ranges.
/// `column_start` and `column_end` identify the shared alignment columns.
pub struct AlignMauveRangeMapping {
  block_index : Int
  column_start : Int
  column_end : Int
  source_start : Int
  source_end : Int
  target_start : Int
  target_end : Int
  source_strand : String
  target_strand : String
} derive(Eq, Debug)
pub extend AlignMauveRangeMapping with Eq::{equal, not_equal}
pub extend AlignMauveRangeMapping with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Aggregate document summary.
pub struct AlignMauveSummary {
  sources : Int
  blocks : Int
  rows : Int
  alignment_columns : Int
  aligned_pairs : Int
  identities : Int
  mismatches : Int
  gaps : Int
} derive(Eq, Debug)
pub extend AlignMauveSummary with Eq::{equal, not_equal}
pub extend AlignMauveSummary with @moonbitlang/core/debug.Debug::{to_repr}

///|
priv struct AlignMauvePendingRow {
  source_number : Int
  start : Int
  end : Int
  strand : String
  description : String
  mut aligned_sequence : String
}

///|
fn align_mauve_fail(message : String) -> Unit raise AlignMauveError {
  raise AlignMauveError(message)
}

///|
fn align_mauve_copy_sources(
  values : Array[AlignMauveSource],
) -> Array[AlignMauveSource] {
  let copied : Array[AlignMauveSource] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn align_mauve_clone_row(value : AlignMauveRow) -> AlignMauveRow {
  let coordinates : Array[Int] = []
  for coordinate in value.coordinates {
    coordinates.push(coordinate)
  }
  AlignMauveRow::{
    source_number: value.source_number,
    identifier: value.identifier,
    start: value.start,
    end: value.end,
    strand: value.strand,
    description: value.description,
    aligned_sequence: value.aligned_sequence,
    source_sequence: value.source_sequence,
    coordinates,
  }
}

///|
fn align_mauve_copy_rows(values : Array[AlignMauveRow]) -> Array[AlignMauveRow] {
  let copied : Array[AlignMauveRow] = []
  for value in values {
    copied.push(align_mauve_clone_row(value))
  }
  copied
}

///|
fn align_mauve_clone_block(value : AlignMauveBlock) -> AlignMauveBlock {
  AlignMauveBlock::{
    rows: align_mauve_copy_rows(value.rows),
    width: value.width,
  }
}

///|
fn align_mauve_copy_blocks(
  values : Array[AlignMauveBlock],
) -> Array[AlignMauveBlock] {
  let copied : Array[AlignMauveBlock] = []
  for value in values {
    copied.push(align_mauve_clone_block(value))
  }
  copied
}

///|
fn align_mauve_copy_metadata(
  values : Array[AlignMauveMetadataEntry],
) -> Array[AlignMauveMetadataEntry] {
  let copied : Array[AlignMauveMetadataEntry] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn align_mauve_copy_ints(values : Array[Int]) -> Array[Int] {
  let copied : Array[Int] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn align_mauve_copy_index_entries(
  values : Array[AlignMauveIndexEntry],
) -> Array[AlignMauveIndexEntry] {
  let copied : Array[AlignMauveIndexEntry] = []
  for value in values {
    copied.push(value)
  }
  copied
}

///|
fn align_mauve_starts_with(value : String, prefix : String) -> Bool {
  if prefix.length() > value.length() {
    return false
  }
  for index in 0..<prefix.length() {
    if value.unsafe_get(index).to_int() != prefix.unsafe_get(index).to_int() {
      return false
    }
  }
  true
}

///|
fn align_mauve_ends_with(value : String, suffix : String) -> Bool {
  if suffix.length() > value.length() {
    return false
  }
  let offset = value.length() - suffix.length()
  for index in 0..<suffix.length() {
    if value.unsafe_get(offset + index).to_int() !=
      suffix.unsafe_get(index).to_int() {
      return false
    }
  }
  true
}

///|
fn align_mauve_is_space(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\r'.to_int() ||
  code == '\n'.to_int()
}

///|
fn align_mauve_trim_line(value : String) -> String {
  let mut start = 0
  let mut end = value.length()
  while start < end && align_mauve_is_space(value.unsafe_get(start).to_int()) {
    start = start + 1
  }
  while end > start && align_mauve_is_space(value.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  value[start:end].to_owned()
}

///|
fn align_mauve_validate_plain_text(
  value : String,
  label : String,
  allow_empty : Bool,
) -> Unit raise AlignMauveError {
  if !allow_empty && value.length() == 0 {
    align_mauve_fail(label + " must not be empty")
  }
  for index in 0..<value.length() {
    let code = value.unsafe_get(index).to_int()
    if code == '\n'.to_int() || code == '\r'.to_int() || code == 0 {
      align_mauve_fail(label + " must not contain line breaks or NUL")
    }
  }
}

///|
fn align_mauve_parse_positive_int(
  text : String,
  label : String,
) -> Int raise AlignMauveError {
  if text.length() == 0 {
    align_mauve_fail(label + " must not be empty")
  }
  let mut value = 0
  for index in 0..<text.length() {
    let code = text.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      align_mauve_fail(label + " must be a decimal integer")
    }
    let digit = code - '0'.to_int()
    if value > 214748364 || (value == 214748364 && digit > 7) {
      align_mauve_fail(label + " exceeds the supported integer range")
    }
    value = value * 10 + digit
  }
  value
}

///|
fn align_mauve_index_of_char(value : String, wanted : Int, start : Int) -> Int {
  let mut index = start
  while index < value.length() {
    if value.unsafe_get(index).to_int() == wanted {
      return index
    }
    index = index + 1
  }
  -1
}

///|
fn align_mauve_header_parts(
  line : String,
  line_number : Int,
) -> (String, String) raise AlignMauveError {
  let body = align_mauve_trim_line(line[1:line.length()].to_owned())
  let mut split = 0
  while split < body.length() &&
        !align_mauve_is_space(body.unsafe_get(split).to_int()) {
    split = split + 1
  }
  if split == 0 || split == body.length() {
    align_mauve_fail(
      "XMFA header line " +
      line_number.to_string() +
      " must contain a key and value",
    )
  }
  let key = body[0:split].to_owned()
  let value = align_mauve_trim_line(body[split:body.length()].to_owned())
  if value.length() == 0 {
    align_mauve_fail(
      "XMFA header line " + line_number.to_string() + " has an empty value",
    )
  }
  (key, value)
}

///|
fn align_mauve_find_metadata(
  values : Array[AlignMauveMetadataEntry],
  key : String,
) -> String? {
  for value in values {
    if value.key == key {
      return Some(value.value)
    }
  }
  None
}

///|
fn align_mauve_find_source_property(
  values : Array[AlignMauveMetadataEntry],
  number : Int,
  suffix : String,
) -> String? {
  align_mauve_find_metadata(values, "Sequence" + number.to_string() + suffix)
}

///|
fn align_mauve_parse_sequence_key(
  key : String,
) -> (Int, String)? raise AlignMauveError {
  if !align_mauve_starts_with(key, "Sequence") {
    return None
  }
  let suffix = if align_mauve_ends_with(key, "File") {
    "File"
  } else if align_mauve_ends_with(key, "Entry") {
    "Entry"
  } else if align_mauve_ends_with(key, "Format") {
    "Format"
  } else {
    align_mauve_fail("Unexpected XMFA sequence header keyword '" + key + "'")
    ""
  }
  let start = "Sequence".length()
  let end = key.length() - suffix.length()
  if end <= start {
    align_mauve_fail("Malformed XMFA sequence header keyword '" + key + "'")
  }
  let number = align_mauve_parse_positive_int(
    key[start:end].to_owned(),
    "XMFA sequence number",
  )
  if number <= 0 {
    align_mauve_fail("XMFA sequence numbers are one-based and must be positive")
  }
  Some((number, suffix))
}

///|
fn align_mauve_find_source(
  sources : Array[AlignMauveSource],
  number : Int,
) -> AlignMauveSource? {
  for source in sources {
    if source.number == number {
      return Some(source)
    }
  }
  None
}

///|
fn align_mauve_find_source_by_identifier(
  sources : Array[AlignMauveSource],
  identifier : String,
) -> AlignMauveSource? {
  for source in sources {
    if source.identifier == identifier {
      return Some(source)
    }
  }
  None
}

///|
fn align_mauve_validate_strand(strand : String) -> Unit raise AlignMauveError {
  if strand != "+" && strand != "-" {
    align_mauve_fail("XMFA strand must be '+' or '-'")
  }
}

///|
fn align_mauve_complement(base : Char) -> Char raise AlignMauveError {
  match base {
    'A' => 'T'
    'C' => 'G'
    'G' => 'C'
    'T' => 'A'
    'U' => 'A'
    'R' => 'Y'
    'Y' => 'R'
    'S' => 'S'
    'W' => 'W'
    'K' => 'M'
    'M' => 'K'
    'B' => 'V'
    'D' => 'H'
    'H' => 'D'
    'V' => 'B'
    'N' => 'N'
    'X' => 'X'
    'a' => 't'
    'c' => 'g'
    'g' => 'c'
    't' => 'a'
    'u' => 'a'
    'r' => 'y'
    'y' => 'r'
    's' => 's'
    'w' => 'w'
    'k' => 'm'
    'm' => 'k'
    'b' => 'v'
    'd' => 'h'
    'h' => 'd'
    'v' => 'b'
    'n' => 'n'
    'x' => 'x'
    _ => {
      align_mauve_fail(
        "XMFA sequence contains non-IUPAC residue '" + base.to_string() + "'",
      )
      'N'
    }
  }
}

///|
fn align_mauve_validate_aligned_sequence(
  sequence : String,
) -> Unit raise AlignMauveError {
  if sequence.length() == 0 {
    align_mauve_fail("XMFA aligned sequence must not be empty")
  }
  for index in 0..<sequence.length() {
    let base = sequence.unsafe_get(index).unsafe_to_char()
    if base != '-' {
      ignore(align_mauve_complement(base))
    }
  }
}

///|
fn align_mauve_remove_gaps(sequence : String) -> String {
  let output = StringBuilder(size_hint=sequence.length())
  for index in 0..<sequence.length() {
    let base = sequence.unsafe_get(index).unsafe_to_char()
    if base != '-' {
      output.write_char(base)
    }
  }
  output.to_string()
}

///|
fn align_mauve_reverse_complement(
  sequence : String,
) -> String raise AlignMauveError {
  let output = StringBuilder(size_hint=sequence.length())
  let mut index = sequence.length() - 1
  while index >= 0 {
    output.write_char(
      align_mauve_complement(sequence.unsafe_get(index).unsafe_to_char()),
    )
    index = index - 1
  }
  output.to_string()
}

///|
fn align_mauve_make_coordinates(
  aligned_sequence : String,
  start : Int,
  end : Int,
  strand : String,
) -> Array[Int] raise AlignMauveError {
  let coordinates : Array[Int] = []
  let mut current = if strand == "+" { start } else { end }
  coordinates.push(current)
  for column in 0..<aligned_sequence.length() {
    if aligned_sequence.unsafe_get(column).to_int() != '-'.to_int() {
      current = current + (if strand == "+" { 1 } else { -1 })
    }
    coordinates.push(current)
  }
  let expected = if strand == "+" { end } else { start }
  if current != expected {
    align_mauve_fail(
      "XMFA row coordinate span does not match its ungapped sequence length",
    )
  }
  coordinates
}

///|
fn align_mauve_parse_description(
  line : String,
  line_number : Int,
) -> AlignMauvePendingRow raise AlignMauveError {
  let body = align_mauve_trim_line(line[1:line.length()].to_owned())
  let mut first_end = 0
  while first_end < body.length() &&
        !align_mauve_is_space(body.unsafe_get(first_end).to_int()) {
    first_end = first_end + 1
  }
  if first_end == 0 || first_end == body.length() {
    align_mauve_fail(
      "Malformed XMFA row description at line " + line_number.to_string(),
    )
  }
  let locus = body[0:first_end].to_owned()
  let mut strand_start = first_end
  while strand_start < body.length() &&
        align_mauve_is_space(body.unsafe_get(strand_start).to_int()) {
    strand_start = strand_start + 1
  }
  let mut strand_end = strand_start
  while strand_end < body.length() &&
        !align_mauve_is_space(body.unsafe_get(strand_end).to_int()) {
    strand_end = strand_end + 1
  }
  if strand_start == strand_end || strand_end == body.length() {
    align_mauve_fail(
      "XMFA row description at line " +
      line_number.to_string() +
      " must include strand and description",
    )
  }
  let strand = body[strand_start:strand_end].to_owned()
  align_mauve_validate_strand(strand)
  let description = align_mauve_trim_line(
    body[strand_end:body.length()].to_owned(),
  )
  align_mauve_validate_plain_text(description, "XMFA row description", false)
  let colon = align_mauve_index_of_char(locus, ':'.to_int(), 0)
  if colon <= 0 ||
    colon + 1 >= locus.length() ||
    align_mauve_index_of_char(locus, ':'.to_int(), colon + 1) >= 0 {
    align_mauve_fail(
      "Malformed XMFA sequence locus at line " + line_number.to_string(),
    )
  }
  let dash = align_mauve_index_of_char(locus, '-'.to_int(), colon + 1)
  if dash <= colon + 1 ||
    dash + 1 >= locus.length() ||
    align_mauve_index_of_char(locus, '-'.to_int(), dash + 1) >= 0 {
    align_mauve_fail(
      "Malformed XMFA coordinate range at line " + line_number.to_string(),
    )
  }
  let source_number = align_mauve_parse_positive_int(
    locus[0:colon].to_owned(),
    "XMFA row sequence number",
  )
  if source_number <= 0 {
    align_mauve_fail("XMFA row sequence number must be positive")
  }
  let file_start = align_mauve_parse_positive_int(
    locus[colon + 1:dash].to_owned(),
    "XMFA row start",
  )
  let file_end = align_mauve_parse_positive_int(
    locus[dash + 1:locus.length()].to_owned(),
    "XMFA row end",
  )
  let (start, end) = if file_start == 0 {
    if file_end != 0 {
      align_mauve_fail("XMFA start coordinate zero is only valid for 0-0 rows")
    }
    (0, 0)
  } else {
    if file_end < file_start {
      align_mauve_fail("XMFA row end must not precede its start")
    }
    (file_start - 1, file_end)
  }
  AlignMauvePendingRow::{
    source_number,
    start,
    end,
    strand,
    description,
    aligned_sequence: "",
  }
}

///|
fn align_mauve_finalize_block(
  pending : Array[AlignMauvePendingRow],
  sources : Array[AlignMauveSource],
) -> AlignMauveBlock raise AlignMauveError {
  if pending.length() == 0 {
    align_mauve_fail("XMFA alignment block must contain at least one row")
  }
  let rows : Array[AlignMauveRow] = []
  let mut width = -1
  let mut has_residue = false
  for item in pending {
    for row in rows {
      if row.source_number == item.source_number {
        align_mauve_fail(
          "XMFA alignment block contains a duplicated sequence number",
        )
      }
    }
    let source = match align_mauve_find_source(sources, item.source_number) {
      Some(value) => value
      None => {
        align_mauve_fail(
          "XMFA row refers to undeclared sequence " +
          item.source_number.to_string(),
        )
        sources[0]
      }
    }
    let row = AlignMauveRow::create(
      source,
      item.start,
      item.end,
      item.strand,
      item.description,
      item.aligned_sequence,
    )
    if width < 0 {
      width = row.aligned_sequence.length()
    } else if row.aligned_sequence.length() != width {
      align_mauve_fail("XMFA rows in one block must have equal widths")
    }
    if row.end > row.start {
      has_residue = true
    }
    rows.push(row)
  }
  if width <= 0 {
    align_mauve_fail("XMFA alignment block width must be positive")
  }
  if !has_residue {
    align_mauve_fail("XMFA alignment block must contain at least one residue")
  }
  AlignMauveBlock::{ rows, width, }
}

///|
/// Construct one validated metadata entry.
pub fn AlignMauveMetadataEntry::create(
  key : String,
  value : String,
) -> AlignMauveMetadataEntry raise AlignMauveError {
  align_mauve_validate_plain_text(key, "XMFA metadata key", false)
  align_mauve_validate_plain_text(value, "XMFA metadata value", false)
  for index in 0..<key.length() {
    if align_mauve_is_space(key.unsafe_get(index).to_int()) {
      align_mauve_fail("XMFA metadata key must not contain whitespace")
    }
  }
  AlignMauveMetadataEntry::{ key, value, }
}

///|
/// Construct one validated XMFA source declaration.
pub fn AlignMauveSource::create(
  number : Int,
  file : String,
  format? : String = "FastA",
  entry? : Int = 0,
  identifier? : String = "",
) -> AlignMauveSource raise AlignMauveError {
  if number <= 0 {
    align_mauve_fail("XMFA source number must be positive")
  }
  align_mauve_validate_plain_text(file, "XMFA source file", false)
  align_mauve_validate_plain_text(format, "XMFA source format", false)
  if entry < 0 {
    align_mauve_fail("XMFA source entry must be non-negative")
  }
  let resolved_identifier = if identifier.length() > 0 {
    identifier
  } else if entry > 0 {
    (entry - 1).to_string()
  } else {
    file
  }
  align_mauve_validate_plain_text(
    resolved_identifier, "XMFA source identifier", false,
  )
  AlignMauveSource::{
    number,
    file,
    entry: if entry == 0 {
      None
    } else {
      Some(entry)
    },
    format,
    identifier: resolved_identifier,
  }
}

///|
/// Construct one validated aligned XMFA row.
pub fn AlignMauveRow::create(
  source : AlignMauveSource,
  start : Int,
  end : Int,
  strand : String,
  description : String,
  aligned_sequence : String,
) -> AlignMauveRow raise AlignMauveError {
  if start < 0 || end < start {
    align_mauve_fail("XMFA row interval must satisfy 0 <= start <= end")
  }
  align_mauve_validate_strand(strand)
  align_mauve_validate_plain_text(description, "XMFA row description", false)
  align_mauve_validate_aligned_sequence(aligned_sequence)
  let displayed_sequence = align_mauve_remove_gaps(aligned_sequence)
  if displayed_sequence.length() != end - start {
    align_mauve_fail(
      "XMFA row coordinate span must equal its ungapped sequence length",
    )
  }
  if start == end && displayed_sequence.length() != 0 {
    align_mauve_fail("XMFA 0-0 row must contain gaps only")
  }
  let source_sequence = if strand == "+" {
    displayed_sequence
  } else {
    align_mauve_reverse_complement(displayed_sequence)
  }
  let coordinates = align_mauve_make_coordinates(
    aligned_sequence, start, end, strand,
  )
  AlignMauveRow::{
    source_number: source.number,
    identifier: source.identifier,
    start,
    end,
    strand,
    description,
    aligned_sequence,
    source_sequence,
    coordinates,
  }
}

///|
/// Construct one validated XMFA block.
pub fn AlignMauveBlock::create(
  rows : Array[AlignMauveRow],
) -> AlignMauveBlock raise AlignMauveError {
  if rows.length() == 0 {
    align_mauve_fail("XMFA block must contain at least one row")
  }
  let copied = align_mauve_copy_rows(rows)
  let width = copied[0].aligned_sequence.length()
  let mut has_residue = false
  for index in 0..<copied.length() {
    let row = copied[index]
    if row.aligned_sequence.length() != width {
      align_mauve_fail("XMFA rows in one block must have equal widths")
    }
    if row.end > row.start {
      has_residue = true
    }
    for other in 0..<index {
      if copied[other].source_number == row.source_number {
        align_mauve_fail(
          "XMFA block contains a duplicated source sequence number",
        )
      }
    }
  }
  if width <= 0 || !has_residue {
    align_mauve_fail(
      "XMFA block must have positive width and at least one residue",
    )
  }
  AlignMauveBlock::{ rows: copied, width, }
}

///|
/// Parse one complete Mauve/ProgressiveMauve XMFA document.
///
/// File coordinates are normalized from one-based inclusive intervals to
/// zero-based half-open intervals. Sequence rows may be wrapped across lines.
/// Every alignment block must be terminated by an `=` line.
pub fn align_mauve_parse(
  text : String,
) -> AlignMauveDocument raise AlignMauveError {
  if text.length() == 0 {
    align_mauve_fail("Empty XMFA input")
  }
  let raw_lines = text.split("\n").to_array()
  let lines : Array[String] = []
  for raw in raw_lines {
    let owned = raw.to_owned()
    let line = if owned.length() > 0 &&
      owned.unsafe_get(owned.length() - 1).to_int() == '\r'.to_int() {
      owned[0:owned.length() - 1].to_owned()
    } else {
      owned
    }
    lines.push(line)
  }
  let all_metadata : Array[AlignMauveMetadataEntry] = []
  let mut index = 0
  while index < lines.length() {
    let trimmed = align_mauve_trim_line(lines[index])
    if trimmed.length() == 0 {
      index = index + 1
      continue
    }
    if !align_mauve_starts_with(trimmed, "#") {
      break
    }
    let (key, value) = align_mauve_header_parts(trimmed, index + 1)
    if align_mauve_find_metadata(all_metadata, key) is Some(_) {
      align_mauve_fail("Duplicated XMFA header keyword '" + key + "'")
    }
    all_metadata.push(AlignMauveMetadataEntry::create(key, value))
    index = index + 1
  }
  if all_metadata.length() == 0 {
    align_mauve_fail("XMFA input must begin with metadata header lines")
  }
  let format_version = match
    align_mauve_find_metadata(all_metadata, "FormatVersion") {
    Some(value) => value
    None => {
      align_mauve_fail("XMFA header is missing FormatVersion")
      ""
    }
  }
  let mut maximum_source = 0
  for entry in all_metadata {
    match align_mauve_parse_sequence_key(entry.key) {
      Some((number, _)) =>
        if number > maximum_source {
          maximum_source = number
        }
      None => ()
    }
  }
  if maximum_source == 0 {
    align_mauve_fail("XMFA header declares no source sequences")
  }
  let source_files : Array[String] = []
  let source_formats : Array[String] = []
  let source_entries : Array[Int] = []
  for number = 1; number <= maximum_source; number = number + 1 {
    let file = match
      align_mauve_find_source_property(all_metadata, number, "File") {
      Some(value) => value
      None => {
        align_mauve_fail(
          "XMFA header is missing Sequence" + number.to_string() + "File",
        )
        ""
      }
    }
    let format = match
      align_mauve_find_source_property(all_metadata, number, "Format") {
      Some(value) => value
      None => {
        align_mauve_fail(
          "XMFA header is missing Sequence" + number.to_string() + "Format",
        )
        ""
      }
    }
    let source_entry = match
      align_mauve_find_source_property(all_metadata, number, "Entry") {
      Some(value) => {
        let parsed = align_mauve_parse_positive_int(value, "XMFA source entry")
        if parsed <= 0 {
          align_mauve_fail("XMFA source entry must be positive")
        }
        parsed
      }
      None => 0
    }
    source_files.push(file)
    source_formats.push(format)
    source_entries.push(source_entry)
  }
  let mut combined = true
  for source_index in 1..<source_files.length() {
    if source_files[source_index] != source_files[0] {
      combined = false
    }
  }
  let sources : Array[AlignMauveSource] = []
  for source_index in 0..<source_files.length() {
    let identifier = if combined {
      if source_entries[source_index] == 0 {
        align_mauve_fail(
          "Combined-file XMFA headers require SequenceNEntry for every source",
        )
      }
      (source_entries[source_index] - 1).to_string()
    } else {
      source_files[source_index]
    }
    for source in sources {
      if source.identifier == identifier {
        align_mauve_fail(
          "XMFA source identifiers must be unique after header normalization",
        )
      }
    }
    sources.push(
      AlignMauveSource::create(
        source_index + 1,
        source_files[source_index],
        format=source_formats[source_index],
        entry=source_entries[source_index],
        identifier~,
      ),
    )
  }
  let metadata : Array[AlignMauveMetadataEntry] = []
  for entry in all_metadata {
    if align_mauve_parse_sequence_key(entry.key) is None {
      metadata.push(entry)
    }
  }
  let blocks : Array[AlignMauveBlock] = []
  let pending : Array[AlignMauvePendingRow] = []
  while index < lines.length() {
    let trimmed = align_mauve_trim_line(lines[index])
    let line_number = index + 1
    index = index + 1
    if trimmed.length() == 0 {
      continue
    }
    if align_mauve_starts_with(trimmed, "#") {
      align_mauve_fail("XMFA metadata is only allowed before alignment blocks")
    }
    if trimmed == "=" {
      blocks.push(align_mauve_finalize_block(pending, sources))
      pending.clear()
      continue
    }
    if align_mauve_starts_with(trimmed, ">") {
      pending.push(align_mauve_parse_description(trimmed, line_number))
      continue
    }
    if pending.length() == 0 {
      align_mauve_fail(
        "XMFA sequence data at line " +
        line_number.to_string() +
        " appears before a row description",
      )
    }
    for char_index in 0..<trimmed.length() {
      if align_mauve_is_space(trimmed.unsafe_get(char_index).to_int()) {
        align_mauve_fail(
          "XMFA sequence lines must not contain internal whitespace",
        )
      }
    }
    let last = pending.length() - 1
    pending[last].aligned_sequence = pending[last].aligned_sequence + trimmed
  }
  if pending.length() > 0 {
    align_mauve_fail("XMFA final alignment block is missing its '=' terminator")
  }
  AlignMauveDocument::{ format_version, sources, metadata, blocks, }
}

///|
/// Serialize an XMFA document in canonical Mauve form.
///
/// `line_width == 0` writes one physical line per aligned row, matching
/// Biopython's writer. Positive widths wrap sequence rows deterministically.
pub fn align_mauve_write(
  document : AlignMauveDocument,
  line_width? : Int = 0,
) -> String raise AlignMauveError {
  if line_width < 0 {
    align_mauve_fail("XMFA line width must be non-negative")
  }
  if document.sources.length() == 0 {
    align_mauve_fail("XMFA document must declare at least one source")
  }
  let output = StringBuilder()
  output.write_string("#FormatVersion ")
  output.write_string(document.format_version)
  output.write_char('\n')
  for source in document.sources {
    output.write_string("#Sequence")
    output.write_string(source.number.to_string())
    output.write_string("File\t")
    output.write_string(source.file)
    output.write_char('\n')
    match source.entry {
      Some(entry) => {
        output.write_string("#Sequence")
        output.write_string(source.number.to_string())
        output.write_string("Entry\t")
        output.write_string(entry.to_string())
        output.write_char('\n')
      }
      None => ()
    }
    output.write_string("#Sequence")
    output.write_string(source.number.to_string())
    output.write_string("Format\t")
    output.write_string(source.format)
    output.write_char('\n')
  }
  for entry in document.metadata {
    if entry.key != "FormatVersion" {
      output.write_char('#')
      output.write_string(entry.key)
      output.write_char('\t')
      output.write_string(entry.value)
      output.write_char('\n')
    }
  }
  for block in document.blocks {
    for row in block.rows {
      output.write_string("> ")
      output.write_string(row.source_number.to_string())
      output.write_char(':')
      if row.start == 0 && row.end == 0 {
        output.write_string("0-0")
      } else {
        output.write_string((row.start + 1).to_string())
        output.write_char('-')
        output.write_string(row.end.to_string())
      }
      output.write_char(' ')
      output.write_string(row.strand)
      output.write_char(' ')
      output.write_string(row.description)
      output.write_char('\n')
      if line_width == 0 {
        output.write_string(row.aligned_sequence)
        output.write_char('\n')
      } else {
        let mut offset = 0
        while offset < row.aligned_sequence.length() {
          let end = if offset + line_width < row.aligned_sequence.length() {
            offset + line_width
          } else {
            row.aligned_sequence.length()
          }
          output.write_string(row.aligned_sequence[offset:end].to_owned())
          output.write_char('\n')
          offset = end
        }
      }
    }
    output.write_string("=\n")
  }
  output.to_string()
}

///|
pub fn AlignMauveMetadataEntry::key(self : AlignMauveMetadataEntry) -> String {
  self.key
}

///|
pub fn AlignMauveMetadataEntry::value(self : AlignMauveMetadataEntry) -> String {
  self.value
}

///|
pub fn AlignMauveSource::number(self : AlignMauveSource) -> Int {
  self.number
}

///|
pub fn AlignMauveSource::file(self : AlignMauveSource) -> String {
  self.file
}

///|
pub fn AlignMauveSource::entry(self : AlignMauveSource) -> Int? {
  self.entry
}

///|
pub fn AlignMauveSource::format(self : AlignMauveSource) -> String {
  self.format
}

///|
pub fn AlignMauveSource::identifier(self : AlignMauveSource) -> String {
  self.identifier
}

///|
pub fn AlignMauveRow::source_number(self : AlignMauveRow) -> Int {
  self.source_number
}

///|
pub fn AlignMauveRow::identifier(self : AlignMauveRow) -> String {
  self.identifier
}

///|
pub fn AlignMauveRow::start(self : AlignMauveRow) -> Int {
  self.start
}

///|
pub fn AlignMauveRow::end(self : AlignMauveRow) -> Int {
  self.end
}

///|
pub fn AlignMauveRow::strand(self : AlignMauveRow) -> String {
  self.strand
}

///|
pub fn AlignMauveRow::description(self : AlignMauveRow) -> String {
  self.description
}

///|
pub fn AlignMauveRow::aligned_sequence(self : AlignMauveRow) -> String {
  self.aligned_sequence
}

///|
pub fn AlignMauveRow::source_sequence(self : AlignMauveRow) -> String {
  self.source_sequence
}

///|
pub fn AlignMauveRow::coordinates(self : AlignMauveRow) -> Array[Int] {
  align_mauve_copy_ints(self.coordinates)
}

///|
/// Return the source boundary at one alignment-column boundary.
pub fn AlignMauveRow::coordinate_at_boundary(
  self : AlignMauveRow,
  column : Int,
) -> Int raise AlignMauveError {
  if column < 0 || column >= self.coordinates.length() {
    align_mauve_fail("XMFA alignment-column boundary is out of range")
  }
  self.coordinates[column]
}

///|
/// Return the zero-based forward source coordinate represented by one residue
/// column, or `None` when this row has a gap in that column.
pub fn AlignMauveRow::residue_coordinate(
  self : AlignMauveRow,
  column : Int,
) -> Int? raise AlignMauveError {
  if column < 0 || column >= self.aligned_sequence.length() {
    align_mauve_fail("XMFA alignment column is out of range")
  }
  if self.aligned_sequence.unsafe_get(column).to_int() == '-'.to_int() {
    return None
  }
  if self.strand == "+" {
    Some(self.coordinates[column])
  } else {
    Some(self.coordinates[column] - 1)
  }
}

///|
/// Locate one forward-axis source residue in this aligned row.
pub fn AlignMauveRow::column_for_coordinate(
  self : AlignMauveRow,
  coordinate : Int,
) -> Int? raise AlignMauveError {
  if coordinate < 0 {
    align_mauve_fail("XMFA source coordinate must be non-negative")
  }
  for column in 0..<self.aligned_sequence.length() {
    match self.residue_coordinate(column) {
      Some(value) => if value == coordinate { return Some(column) }
      None => ()
    }
  }
  None
}

///|
pub fn AlignMauveBlock::rows(self : AlignMauveBlock) -> Array[AlignMauveRow] {
  align_mauve_copy_rows(self.rows)
}

///|
pub fn AlignMauveBlock::row_count(self : AlignMauveBlock) -> Int {
  self.rows.length()
}

///|
pub fn AlignMauveBlock::width(self : AlignMauveBlock) -> Int {
  self.width
}

///|
pub fn AlignMauveBlock::row_by_identifier(
  self : AlignMauveBlock,
  identifier : String,
) -> AlignMauveRow? {
  for row in self.rows {
    if row.identifier == identifier {
      return Some(row)
    }
  }
  None
}

///|
pub fn AlignMauveBlock::row_by_source_number(
  self : AlignMauveBlock,
  source_number : Int,
) -> AlignMauveRow? {
  for row in self.rows {
    if row.source_number == source_number {
      return Some(row)
    }
  }
  None
}

///|
/// Return Biopython-style compact coordinate rows.
///
/// A boundary is retained whenever any row changes between residue movement
/// and gap movement, plus the first and final boundaries.
pub fn AlignMauveBlock::compact_coordinates(
  self : AlignMauveBlock,
) -> Array[Array[Int]] {
  let boundaries : Array[Int] = [0]
  for column = 1; column < self.width; column = column + 1 {
    let mut changed = false
    for row in self.rows {
      let previous = row.coordinates[column] - row.coordinates[column - 1]
      let next = row.coordinates[column + 1] - row.coordinates[column]
      if previous != next {
        changed = true
      }
    }
    if changed {
      boundaries.push(column)
    }
  }
  boundaries.push(self.width)
  let compact : Array[Array[Int]] = []
  for row in self.rows {
    let coordinates : Array[Int] = []
    for boundary in boundaries {
      coordinates.push(row.coordinates[boundary])
    }
    compact.push(coordinates)
  }
  compact
}

///|
/// Project one source residue to a target residue through a shared block
/// column. A gap in either row yields `None`.
pub fn AlignMauveBlock::map_position(
  self : AlignMauveBlock,
  source_identifier : String,
  target_identifier : String,
  source_position : Int,
) -> Int? raise AlignMauveError {
  if source_position < 0 {
    align_mauve_fail("XMFA source position must be non-negative")
  }
  let source = match self.row_by_identifier(source_identifier) {
    Some(row) => row
    None => return None
  }
  let target = match self.row_by_identifier(target_identifier) {
    Some(row) => row
    None => return None
  }
  let column = match source.column_for_coordinate(source_position) {
    Some(value) => value
    None => return None
  }
  target.residue_coordinate(column)
}

///|
/// Calculate pairwise alignment counts over all unordered row pairs.
pub fn AlignMauveBlock::counts(self : AlignMauveBlock) -> AlignMauveCounts {
  let mut pair_count = 0
  let mut columns = 0
  let mut aligned = 0
  let mut identities = 0
  let mut mismatches = 0
  let mut gaps = 0
  let mut gap_opens = 0
  let mut double_gap_columns = 0
  for left = 0; left < self.rows.length(); left = left + 1 {
    for right = left + 1; right < self.rows.length(); right = right + 1 {
      pair_count = pair_count + 1
      columns = columns + self.width
      let mut active_gap = 0
      for column in 0..<self.width {
        let a = self.rows[left].aligned_sequence
          .unsafe_get(column)
          .unsafe_to_char()
        let b = self.rows[right].aligned_sequence
          .unsafe_get(column)
          .unsafe_to_char()
        if a == '-' && b == '-' {
          double_gap_columns = double_gap_columns + 1
        } else if a == '-' || b == '-' {
          gaps = gaps + 1
          let direction = if a == '-' { 1 } else { 2 }
          if active_gap != direction {
            gap_opens = gap_opens + 1
          }
          active_gap = direction
        } else {
          active_gap = 0
          aligned = aligned + 1
          if a.to_string().to_upper() == b.to_string().to_upper() {
            identities = identities + 1
          } else {
            mismatches = mismatches + 1
          }
        }
      }
    }
  }
  AlignMauveCounts::{
    pair_count,
    columns,
    aligned,
    identities,
    mismatches,
    gaps,
    gap_opens,
    double_gap_columns,
  }
}

///|
pub fn AlignMauveDocument::format_version(self : AlignMauveDocument) -> String {
  self.format_version
}

///|
pub fn AlignMauveDocument::sources(
  self : AlignMauveDocument,
) -> Array[AlignMauveSource] {
  align_mauve_copy_sources(self.sources)
}

///|
pub fn AlignMauveDocument::blocks(
  self : AlignMauveDocument,
) -> Array[AlignMauveBlock] {
  align_mauve_copy_blocks(self.blocks)
}

///|
pub fn AlignMauveDocument::metadata_entries(
  self : AlignMauveDocument,
) -> Array[AlignMauveMetadataEntry] {
  align_mauve_copy_metadata(self.metadata)
}

///|
pub fn AlignMauveDocument::metadata_value(
  self : AlignMauveDocument,
  key : String,
) -> String? {
  align_mauve_find_metadata(self.metadata, key)
}

///|
pub fn AlignMauveDocument::source_count(self : AlignMauveDocument) -> Int {
  self.sources.length()
}

///|
pub fn AlignMauveDocument::block_count(self : AlignMauveDocument) -> Int {
  self.blocks.length()
}

///|
pub fn AlignMauveDocument::identifiers(
  self : AlignMauveDocument,
) -> Array[String] {
  let identifiers : Array[String] = []
  for source in self.sources {
    identifiers.push(source.identifier)
  }
  identifiers
}

///|
pub fn AlignMauveDocument::source(
  self : AlignMauveDocument,
  identifier : String,
) -> AlignMauveSource? {
  align_mauve_find_source_by_identifier(self.sources, identifier)
}

///|
/// Sum pairwise counts over all XMFA blocks.
pub fn AlignMauveDocument::counts(
  self : AlignMauveDocument,
) -> AlignMauveCounts {
  let mut pair_count = 0
  let mut columns = 0
  let mut aligned = 0
  let mut identities = 0
  let mut mismatches = 0
  let mut gaps = 0
  let mut gap_opens = 0
  let mut double_gap_columns = 0
  for block in self.blocks {
    let counts = block.counts()
    pair_count = pair_count + counts.pair_count
    columns = columns + counts.columns
    aligned = aligned + counts.aligned
    identities = identities + counts.identities
    mismatches = mismatches + counts.mismatches
    gaps = gaps + counts.gaps
    gap_opens = gap_opens + counts.gap_opens
    double_gap_columns = double_gap_columns + counts.double_gap_columns
  }
  AlignMauveCounts::{
    pair_count,
    columns,
    aligned,
    identities,
    mismatches,
    gaps,
    gap_opens,
    double_gap_columns,
  }
}

///|
pub fn AlignMauveCounts::pair_count(self : AlignMauveCounts) -> Int {
  self.pair_count
}

///|
pub fn AlignMauveCounts::columns(self : AlignMauveCounts) -> Int {
  self.columns
}

///|
pub fn AlignMauveCounts::aligned(self : AlignMauveCounts) -> Int {
  self.aligned
}

///|
pub fn AlignMauveCounts::identities(self : AlignMauveCounts) -> Int {
  self.identities
}

///|
pub fn AlignMauveCounts::mismatches(self : AlignMauveCounts) -> Int {
  self.mismatches
}

///|
pub fn AlignMauveCounts::gaps(self : AlignMauveCounts) -> Int {
  self.gaps
}

///|
pub fn AlignMauveCounts::gap_opens(self : AlignMauveCounts) -> Int {
  self.gap_opens
}

///|
pub fn AlignMauveCounts::double_gap_columns(self : AlignMauveCounts) -> Int {
  self.double_gap_columns
}

///|
/// Build an in-memory interval index over all non-empty rows.
pub fn align_mauve_index(document : AlignMauveDocument) -> AlignMauveIndex {
  let entries : Array[AlignMauveIndexEntry] = []
  for block_index in 0..<document.blocks.length() {
    let block = document.blocks[block_index]
    for row_index in 0..<block.rows.length() {
      let row = block.rows[row_index]
      if row.end > row.start {
        entries.push(AlignMauveIndexEntry::{
          identifier: row.identifier,
          start: row.start,
          end: row.end,
          block_index,
          row_index,
          strand: row.strand,
        })
      }
    }
  }
  AlignMauveIndex::{ identifiers: document.identifiers(), entries, }
}

///|
pub fn AlignMauveIndex::entries(
  self : AlignMauveIndex,
) -> Array[AlignMauveIndexEntry] {
  align_mauve_copy_index_entries(self.entries)
}

///|
/// Find rows overlapping one zero-based half-open source interval.
pub fn AlignMauveIndex::query(
  self : AlignMauveIndex,
  identifier : String,
  start : Int,
  end : Int,
) -> Array[AlignMauveIndexEntry] raise AlignMauveError {
  if start < 0 || end <= start {
    align_mauve_fail("XMFA index query must satisfy 0 <= start < end")
  }
  if !self.identifiers.contains(identifier) {
    align_mauve_fail("Unknown XMFA source identifier '" + identifier + "'")
  }
  let matches : Array[AlignMauveIndexEntry] = []
  for entry in self.entries {
    if entry.identifier == identifier && entry.start < end && entry.end > start {
      matches.push(entry)
    }
  }
  matches
}

///|
pub fn AlignMauveIndexEntry::identifier(self : AlignMauveIndexEntry) -> String {
  self.identifier
}

///|
pub fn AlignMauveIndexEntry::start(self : AlignMauveIndexEntry) -> Int {
  self.start
}

///|
pub fn AlignMauveIndexEntry::end(self : AlignMauveIndexEntry) -> Int {
  self.end
}

///|
pub fn AlignMauveIndexEntry::block_index(self : AlignMauveIndexEntry) -> Int {
  self.block_index
}

///|
pub fn AlignMauveIndexEntry::row_index(self : AlignMauveIndexEntry) -> Int {
  self.row_index
}

///|
pub fn AlignMauveIndexEntry::strand(self : AlignMauveIndexEntry) -> String {
  self.strand
}

///|
/// Project one source residue through every block containing both sequences.
pub fn align_mauve_map_position(
  document : AlignMauveDocument,
  source_identifier : String,
  target_identifier : String,
  source_position : Int,
) -> Array[AlignMauvePositionMapping] raise AlignMauveError {
  if align_mauve_find_source_by_identifier(document.sources, source_identifier)
    is None {
    align_mauve_fail(
      "Unknown XMFA source identifier '" + source_identifier + "'",
    )
  }
  if align_mauve_find_source_by_identifier(document.sources, target_identifier)
    is None {
    align_mauve_fail(
      "Unknown XMFA target identifier '" + target_identifier + "'",
    )
  }
  if source_position < 0 {
    align_mauve_fail("XMFA source position must be non-negative")
  }
  let mappings : Array[AlignMauvePositionMapping] = []
  for block_index in 0..<document.blocks.length() {
    let block = document.blocks[block_index]
    let source = match block.row_by_identifier(source_identifier) {
      Some(row) => row
      None => continue
    }
    let target = match block.row_by_identifier(target_identifier) {
      Some(row) => row
      None => continue
    }
    let column = match source.column_for_coordinate(source_position) {
      Some(value) => value
      None => continue
    }
    match target.residue_coordinate(column) {
      Some(target_position) =>
        mappings.push(AlignMauvePositionMapping::{
          block_index,
          column,
          source_position,
          target_position,
        })
      None => ()
    }
  }
  mappings
}

///|
pub fn AlignMauvePositionMapping::block_index(
  self : AlignMauvePositionMapping,
) -> Int {
  self.block_index
}

///|
pub fn AlignMauvePositionMapping::column(
  self : AlignMauvePositionMapping,
) -> Int {
  self.column
}

///|
pub fn AlignMauvePositionMapping::source_position(
  self : AlignMauvePositionMapping,
) -> Int {
  self.source_position
}

///|
pub fn AlignMauvePositionMapping::target_position(
  self : AlignMauvePositionMapping,
) -> Int {
  self.target_position
}

///|
fn align_mauve_append_range_mapping(
  mappings : Array[AlignMauveRangeMapping],
  block_index : Int,
  column_start : Int,
  column_end : Int,
  first_source : Int,
  last_source : Int,
  first_target : Int,
  last_target : Int,
  source_strand : String,
  target_strand : String,
) -> Unit {
  let source_start = if first_source < last_source {
    first_source
  } else {
    last_source
  }
  let source_end = if first_source > last_source {
    first_source + 1
  } else {
    last_source + 1
  }
  let target_start = if first_target < last_target {
    first_target
  } else {
    last_target
  }
  let target_end = if first_target > last_target {
    first_target + 1
  } else {
    last_target + 1
  }
  mappings.push(AlignMauveRangeMapping::{
    block_index,
    column_start,
    column_end,
    source_start,
    source_end,
    target_start,
    target_end,
    source_strand,
    target_strand,
  })
}

///|
/// Project a source interval through all shared, residue-aligned block runs.
///
/// Gaps split the result into separate mappings. Reverse rows retain their
/// strand labels while returned coordinate intervals remain forward-axis.
pub fn align_mauve_map_range(
  document : AlignMauveDocument,
  source_identifier : String,
  target_identifier : String,
  start : Int,
  end : Int,
) -> Array[AlignMauveRangeMapping] raise AlignMauveError {
  if start < 0 || end <= start {
    align_mauve_fail("XMFA source interval must satisfy 0 <= start < end")
  }
  if align_mauve_find_source_by_identifier(document.sources, source_identifier)
    is None {
    align_mauve_fail(
      "Unknown XMFA source identifier '" + source_identifier + "'",
    )
  }
  if align_mauve_find_source_by_identifier(document.sources, target_identifier)
    is None {
    align_mauve_fail(
      "Unknown XMFA target identifier '" + target_identifier + "'",
    )
  }
  let mappings : Array[AlignMauveRangeMapping] = []
  for block_index in 0..<document.blocks.length() {
    let block = document.blocks[block_index]
    let source = match block.row_by_identifier(source_identifier) {
      Some(row) => row
      None => continue
    }
    let target = match block.row_by_identifier(target_identifier) {
      Some(row) => row
      None => continue
    }
    let source_step = if source.strand == "+" { 1 } else { -1 }
    let target_step = if target.strand == "+" { 1 } else { -1 }
    let mut active = false
    let mut column_start = 0
    let mut first_source = 0
    let mut last_source = 0
    let mut first_target = 0
    let mut last_target = 0
    let mut previous_column = -2
    for column in 0..<block.width {
      let source_position = source.residue_coordinate(column)
      let target_position = target.residue_coordinate(column)
      let usable = match (source_position, target_position) {
        (Some(source_value), Some(_)) =>
          source_value >= start && source_value < end
        _ => false
      }
      if usable {
        let source_value = source_position.unwrap()
        let target_value = target_position.unwrap()
        let contiguous = active &&
          column == previous_column + 1 &&
          source_value == last_source + source_step &&
          target_value == last_target + target_step
        if !contiguous {
          if active {
            align_mauve_append_range_mapping(
              mappings,
              block_index,
              column_start,
              previous_column + 1,
              first_source,
              last_source,
              first_target,
              last_target,
              source.strand,
              target.strand,
            )
          }
          active = true
          column_start = column
          first_source = source_value
          first_target = target_value
        }
        last_source = source_value
        last_target = target_value
        previous_column = column
      } else if active {
        align_mauve_append_range_mapping(
          mappings,
          block_index,
          column_start,
          previous_column + 1,
          first_source,
          last_source,
          first_target,
          last_target,
          source.strand,
          target.strand,
        )
        active = false
      }
    }
    if active {
      align_mauve_append_range_mapping(
        mappings,
        block_index,
        column_start,
        previous_column + 1,
        first_source,
        last_source,
        first_target,
        last_target,
        source.strand,
        target.strand,
      )
    }
  }
  mappings
}

///|
pub fn AlignMauveRangeMapping::block_index(
  self : AlignMauveRangeMapping,
) -> Int {
  self.block_index
}

///|
pub fn AlignMauveRangeMapping::column_start(
  self : AlignMauveRangeMapping,
) -> Int {
  self.column_start
}

///|
pub fn AlignMauveRangeMapping::column_end(self : AlignMauveRangeMapping) -> Int {
  self.column_end
}

///|
pub fn AlignMauveRangeMapping::source_start(
  self : AlignMauveRangeMapping,
) -> Int {
  self.source_start
}

///|
pub fn AlignMauveRangeMapping::source_end(self : AlignMauveRangeMapping) -> Int {
  self.source_end
}

///|
pub fn AlignMauveRangeMapping::target_start(
  self : AlignMauveRangeMapping,
) -> Int {
  self.target_start
}

///|
pub fn AlignMauveRangeMapping::target_end(self : AlignMauveRangeMapping) -> Int {
  self.target_end
}

///|
pub fn AlignMauveRangeMapping::source_strand(
  self : AlignMauveRangeMapping,
) -> String {
  self.source_strand
}

///|
pub fn AlignMauveRangeMapping::target_strand(
  self : AlignMauveRangeMapping,
) -> String {
  self.target_strand
}

///|
/// Reconstruct the known forward source sequence from all XMFA blocks.
///
/// Unknown positions up to the largest observed endpoint use `fill`, which
/// must be exactly one non-gap IUPAC character. Conflicting overlapping rows
/// are rejected.
pub fn align_mauve_reconstruct_sequence(
  document : AlignMauveDocument,
  identifier : String,
  fill? : String = "N",
) -> String raise AlignMauveError {
  if align_mauve_find_source_by_identifier(document.sources, identifier) is None {
    align_mauve_fail("Unknown XMFA source identifier '" + identifier + "'")
  }
  if fill.length() != 1 || fill == "-" {
    align_mauve_fail(
      "XMFA reconstruction fill must be one non-gap IUPAC character",
    )
  }
  let fill_char = fill.unsafe_get(0).unsafe_to_char()
  ignore(align_mauve_complement(fill_char))
  let mut length = 0
  for block in document.blocks {
    for row in block.rows {
      if row.identifier == identifier && row.end > length {
        length = row.end
      }
    }
  }
  let sequence : Array[Char] = []
  for _ in 0..<length {
    sequence.push(fill_char)
  }
  let assigned = Array::make(length, false)
  for block in document.blocks {
    for row in block.rows {
      if row.identifier != identifier || row.end == row.start {
        continue
      }
      for offset in 0..<row.source_sequence.length() {
        let position = row.start + offset
        let base = row.source_sequence.unsafe_get(offset).unsafe_to_char()
        if assigned[position] && sequence[position] != base {
          align_mauve_fail(
            "Conflicting XMFA residues overlap at source position " +
            position.to_string(),
          )
        }
        sequence[position] = base
        assigned[position] = true
      }
    }
  }
  let output = StringBuilder(size_hint=length)
  for base in sequence {
    output.write_char(base)
  }
  output.to_string()
}

///|
/// Return aggregate document statistics.
pub fn align_mauve_summary(document : AlignMauveDocument) -> AlignMauveSummary {
  let counts = document.counts()
  let mut rows = 0
  let mut alignment_columns = 0
  for block in document.blocks {
    rows = rows + block.rows.length()
    alignment_columns = alignment_columns + block.width
  }
  AlignMauveSummary::{
    sources: document.sources.length(),
    blocks: document.blocks.length(),
    rows,
    alignment_columns,
    aligned_pairs: counts.aligned,
    identities: counts.identities,
    mismatches: counts.mismatches,
    gaps: counts.gaps,
  }
}

///|
pub fn AlignMauveSummary::sources(self : AlignMauveSummary) -> Int {
  self.sources
}

///|
pub fn AlignMauveSummary::blocks(self : AlignMauveSummary) -> Int {
  self.blocks
}

///|
pub fn AlignMauveSummary::rows(self : AlignMauveSummary) -> Int {
  self.rows
}

///|
pub fn AlignMauveSummary::alignment_columns(self : AlignMauveSummary) -> Int {
  self.alignment_columns
}

///|
pub fn AlignMauveSummary::aligned_pairs(self : AlignMauveSummary) -> Int {
  self.aligned_pairs
}

///|
pub fn AlignMauveSummary::identities(self : AlignMauveSummary) -> Int {
  self.identities
}

///|
pub fn AlignMauveSummary::mismatches(self : AlignMauveSummary) -> Int {
  self.mismatches
}

///|
pub fn AlignMauveSummary::gaps(self : AlignMauveSummary) -> Int {
  self.gaps
}

///|
/// Small deterministic XMFA fixture adapted from Biopython's official Mauve
/// combined-file test data.
pub fn align_mauve_example_text() -> String {
  "#FormatVersion Mauve1\n" +
  "#Sequence1File\tcombined.fa\n" +
  "#Sequence1Entry\t1\n" +
  "#Sequence1Format\tFastA\n" +
  "#Sequence2File\tcombined.fa\n" +
  "#Sequence2Entry\t2\n" +
  "#Sequence2Format\tFastA\n" +
  "#Sequence3File\tcombined.fa\n" +
  "#Sequence3Entry\t3\n" +
  "#Sequence3Format\tFastA\n" +
  "#BackboneFile\tcombined.xmfa.bbcols\n" +
  "> 1:2-49 - combined.fa\n" +
  "AAGCCCTCCTAGCACACACCCGGAGTGG-CCGGGCCGTACTTTCCTTTT\n" +
  "> 2:0-0 + combined.fa\n" +
  "-------------------------------------------------\n" +
  "> 3:2-48 + combined.fa\n" +
  "AAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGCTTTCCTTTT\n" +
  "=\n" +
  "> 1:1-1 + combined.fa\n" +
  "G\n" +
  "=\n" +
  "> 1:50-50 + combined.fa\n" +
  "A\n" +
  "=\n" +
  "> 2:1-41 + combined.fa\n" +
  "GAAGAGGAAAAGTAGATCCCTGGCGTCCGGAGCTGGGACGT\n" +
  "=\n" +
  "> 3:1-1 + combined.fa\n" +
  "C\n" +
  "=\n" +
  "> 3:49-49 + combined.fa\n" +
  "C\n" +
  "=\n"
}
