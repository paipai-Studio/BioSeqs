// Biopython-compatible Bio.Align.nexus support.
//
// This module is intentionally separate from nexus.mbt. The older module
// exposes Bio.Nexus-style block access, while this file models one modern,
// coordinate-aware alignment matrix.

///|
/// Error raised for malformed NEXUS alignment data or invalid operations.
pub suberror AlignNexusError {
  AlignNexusError(String)
}

///|
/// Sequence datatype declared by a NEXUS FORMAT command.
pub(all) enum AlignNexusDataType {
  AlignNexusDna
  AlignNexusRna
  AlignNexusProtein
  AlignNexusStandard
} derive(Eq, Debug)

///|
/// File-level metadata retained from a DATA or CHARACTERS block.
///
/// `declared_characters` is the NCHAR value in the source file.
/// `removed_all_gap_columns` records columns omitted from the alignment,
/// matching Biopython's `Bio.Align.nexus` behavior.
pub struct AlignNexusMetadata {
  block_name : String
  declared_taxa : Int
  declared_characters : Int
  data_type : AlignNexusDataType
  missing_character : String
  gap_character : String
  match_character : String?
  interleaved : Bool
  respect_case : Bool
  symbols : String
  removed_all_gap_columns : Int
} derive(Eq, Debug)

///|
/// One NEXUS matrix row.
///
/// `aligned_sequence` always uses `-` for gaps. `sequence` is the same row
/// with gaps removed; missing-data symbols remain coordinate-bearing.
pub struct AlignNexusSequence {
  id : String
  sequence : String
  aligned_sequence : String
} derive(Eq, Debug)

///|
/// A coordinate-aware NEXUS multiple sequence alignment.
pub struct AlignNexusAlignment {
  metadata : AlignNexusMetadata
  sequences : Array[AlignNexusSequence]
} derive(Eq, Debug)

///|
/// Pairwise or all-pairs alignment statistics.
pub struct AlignNexusCounts {
  pairs : Int
  aligned : Int
  identities : Int
  mismatches : Int
  gap_columns : Int
  double_gap_columns : Int
  gap_opens : Int
} derive(Eq, Debug)

///|
/// Internal FORMAT state used while parsing commands.
priv struct AlignNexusFormat {
  data_type : AlignNexusDataType
  missing_character : String
  gap_character : String
  match_character : String?
  interleaved : Bool
  respect_case : Bool
  symbols : String
  labels : Bool
}

///|
/// Construct validated NEXUS alignment metadata.
pub fn AlignNexusMetadata::create(
  declared_taxa : Int,
  declared_characters : Int,
  data_type : AlignNexusDataType,
  block_name? : String = "data",
  missing_character? : String = "?",
  gap_character? : String = "-",
  match_character? : String? = None,
  interleaved? : Bool = false,
  respect_case? : Bool = false,
  symbols? : String = "",
  removed_all_gap_columns? : Int = 0,
) -> AlignNexusMetadata raise AlignNexusError {
  if declared_taxa <= 0 {
    raise AlignNexusError("NEXUS NTAX must be positive")
  }
  if declared_characters <= 0 {
    raise AlignNexusError("NEXUS NCHAR must be positive")
  }
  let normalized_block = block_name.to_lower().trim().to_owned()
  if normalized_block != "data" && normalized_block != "characters" {
    raise AlignNexusError("NEXUS alignment block must be DATA or CHARACTERS")
  }
  align_nexus_validate_format_char(
    missing_character, "NEXUS missing-data character",
  )
  align_nexus_validate_format_char(gap_character, "NEXUS gap character")
  if missing_character == gap_character {
    raise AlignNexusError(
      "NEXUS missing-data and gap characters must be different",
    )
  }
  match match_character {
    Some(value) => {
      align_nexus_validate_format_char(value, "NEXUS match character")
      if value == missing_character || value == gap_character {
        raise AlignNexusError(
          "NEXUS match character must differ from missing and gap characters",
        )
      }
    }
    None => ()
  }
  align_nexus_validate_symbols(data_type, symbols, respect_case)
  if removed_all_gap_columns < 0 ||
    removed_all_gap_columns >= declared_characters {
    raise AlignNexusError("NEXUS removed all-gap column count is outside NCHAR")
  }
  AlignNexusMetadata::{
    block_name: normalized_block,
    declared_taxa,
    declared_characters,
    data_type,
    missing_character,
    gap_character,
    match_character,
    interleaved,
    respect_case,
    symbols,
    removed_all_gap_columns,
  }
}

///|
/// Construct and normalize one standalone NEXUS row.
pub fn AlignNexusSequence::create(
  id : String,
  aligned_sequence : String,
  data_type : AlignNexusDataType,
  missing_character? : String = "?",
  gap_character? : String = "-",
  symbols? : String = "",
  respect_case? : Bool = false,
) -> AlignNexusSequence raise AlignNexusError {
  align_nexus_validate_id(id)
  if aligned_sequence.length() == 0 {
    raise AlignNexusError("NEXUS aligned sequence must not be empty")
  }
  align_nexus_validate_format_char(
    missing_character, "NEXUS missing-data character",
  )
  align_nexus_validate_format_char(gap_character, "NEXUS gap character")
  if missing_character == gap_character {
    raise AlignNexusError(
      "NEXUS missing-data and gap characters must be different",
    )
  }
  align_nexus_validate_symbols(data_type, symbols, respect_case)
  let normalized = align_nexus_normalize_row(
    aligned_sequence, data_type, missing_character, gap_character, symbols, respect_case,
  )
  AlignNexusSequence::{
    id,
    sequence: align_nexus_remove_gaps(normalized),
    aligned_sequence: normalized,
  }
}

///|
/// Construct a validated coordinate-aware NEXUS alignment.
pub fn AlignNexusAlignment::create(
  metadata : AlignNexusMetadata,
  sequences : Array[AlignNexusSequence],
) -> AlignNexusAlignment raise AlignNexusError {
  let copied : Array[AlignNexusSequence] = []
  for sequence in sequences {
    copied.push(sequence)
  }
  let alignment = AlignNexusAlignment::{ metadata, sequences: copied, }
  align_nexus_validate_alignment(alignment)
  alignment
}

///|
/// Construct an alignment directly from identifiers and printed rows.
///
/// Source gap symbols are normalized to `-`, and columns containing only gaps
/// are removed exactly as they are when parsing a NEXUS file.
pub fn align_nexus_from_aligned(
  ids : Array[String],
  aligned_sequences : Array[String],
  data_type : AlignNexusDataType,
  missing_character? : String = "?",
  gap_character? : String = "-",
  symbols? : String = "",
  respect_case? : Bool = false,
) -> AlignNexusAlignment raise AlignNexusError {
  if ids.length() == 0 {
    raise AlignNexusError("NEXUS alignment must contain at least one sequence")
  }
  if ids.length() != aligned_sequences.length() {
    raise AlignNexusError(
      "NEXUS identifiers and aligned rows must have equal lengths",
    )
  }
  let width = aligned_sequences[0].length()
  if width == 0 {
    raise AlignNexusError("NEXUS alignment must contain at least one column")
  }
  for row in aligned_sequences {
    if row.length() != width {
      raise AlignNexusError("NEXUS aligned rows must have equal widths")
    }
  }
  let normalized : Array[String] = []
  for row in aligned_sequences {
    normalized.push(
      align_nexus_normalize_row(
        row, data_type, missing_character, gap_character, symbols, respect_case,
      ),
    )
  }
  let prepared = align_nexus_remove_all_gap_columns(normalized)
  if prepared.0.length() == 0 || prepared.0[0].length() == 0 {
    raise AlignNexusError(
      "NEXUS alignment must contain a non-gap alignment column",
    )
  }
  let metadata = AlignNexusMetadata::create(
    ids.length(),
    width,
    data_type,
    missing_character~,
    gap_character~,
    respect_case~,
    symbols~,
    removed_all_gap_columns=prepared.1,
  )
  let sequences : Array[AlignNexusSequence] = []
  for index = 0; index < ids.length(); index = index + 1 {
    align_nexus_validate_id(ids[index])
    sequences.push(AlignNexusSequence::{
      id: ids[index],
      sequence: align_nexus_remove_gaps(prepared.0[index]),
      aligned_sequence: prepared.0[index],
    })
  }
  AlignNexusAlignment::create(metadata, sequences)
}

///|
/// Parse one Biopython-compatible NEXUS alignment.
///
/// The parser recognizes DATA and CHARACTERS matrices, TAXA/TAXLABELS,
/// sequential and interleaved layouts, nested comments, quoted identifiers,
/// doubled quote escaping, custom missing/gap/match symbols, and the DNA, RNA,
/// protein, nucleotide, and standard datatypes.
pub fn align_nexus_parse(
  text : String,
) -> AlignNexusAlignment raise AlignNexusError {
  let lines = align_nexus_normalize_lines(text)
  if lines.length() == 0 ||
    (lines.length() == 1 && lines[0].trim().length() == 0) {
    raise AlignNexusError("Empty NEXUS file")
  }
  if lines[0].trim() != "#NEXUS" {
    raise AlignNexusError("File does not start with NEXUS header")
  }
  let body_builder = StringBuilder::new()
  for index = 1; index < lines.length(); index = index + 1 {
    if index > 1 {
      body_builder.write_char('\n')
    }
    body_builder.write_string(lines[index])
  }
  let commands = align_nexus_scan_commands(body_builder.to_string())
  let mut active_block = ""
  let mut tax_ntax = -1
  let mut data_ntax = -1
  let mut data_nchar = -1
  let taxlabels : Array[String] = []
  let mut block_name = ""
  let mut matrix_body : String? = None
  let mut format = AlignNexusFormat::{
    data_type: AlignNexusDna,
    missing_character: "?",
    gap_character: "-",
    match_character: None,
    interleaved: false,
    respect_case: false,
    symbols: "",
    labels: true,
  }
  for command in commands {
    let split = align_nexus_split_keyword(command)
    let keyword = split.0
    let options = split.1
    if keyword == "begin" {
      if active_block.length() != 0 {
        raise AlignNexusError("NEXUS blocks must not be nested")
      }
      let words = align_nexus_tokenize(options)
      if words.length() != 1 {
        raise AlignNexusError("Malformed NEXUS BEGIN command")
      }
      active_block = words[0].to_lower()
      continue
    }
    if keyword == "end" || keyword == "endblock" {
      if active_block.length() == 0 {
        raise AlignNexusError("Unmatched NEXUS END command")
      }
      active_block = ""
      continue
    }
    if active_block.length() == 0 {
      raise AlignNexusError(
        "NEXUS command '" + keyword + "' appears outside a block",
      )
    }
    if active_block == "taxa" {
      if keyword == "dimensions" {
        let dimensions = align_nexus_parse_dimensions(options)
        if dimensions.0 >= 0 {
          tax_ntax = dimensions.0
        }
        if dimensions.1 >= 0 {
          raise AlignNexusError("NCHAR is not valid in a TAXA block")
        }
      } else if keyword == "taxlabels" {
        let parsed_labels = align_nexus_parse_taxlabels(options)
        taxlabels.clear()
        for label in parsed_labels {
          taxlabels.push(label)
        }
      } else if keyword != "title" && keyword != "link" {
        raise AlignNexusError(
          "Unsupported command '" + keyword + "' in NEXUS TAXA block",
        )
      }
      continue
    }
    if active_block == "data" || active_block == "characters" {
      if keyword == "dimensions" {
        let dimensions = align_nexus_parse_dimensions(options)
        if dimensions.0 >= 0 {
          data_ntax = dimensions.0
        }
        if dimensions.1 >= 0 {
          data_nchar = dimensions.1
        }
      } else if keyword == "format" {
        format = align_nexus_parse_format(options, format)
      } else if keyword == "taxlabels" {
        let parsed_labels = align_nexus_parse_taxlabels(options)
        taxlabels.clear()
        for label in parsed_labels {
          taxlabels.push(label)
        }
      } else if keyword == "matrix" {
        match matrix_body {
          Some(_) =>
            raise AlignNexusError(
              "A NEXUS file may contain only one alignment matrix",
            )
          None => {
            matrix_body = Some(options)
            block_name = active_block
          }
        }
      } else if keyword != "title" &&
        keyword != "link" &&
        keyword != "options" &&
        keyword != "charlabels" &&
        keyword != "charstatelabels" &&
        keyword != "statelabels" &&
        keyword != "eliminate" {
        raise AlignNexusError(
          "Unsupported command '" +
          keyword +
          "' in NEXUS " +
          active_block.to_upper() +
          " block",
        )
      }
      continue
    }
    // Alignment parsing deliberately ignores unrelated SETS, TREES, CODONS,
    // and vendor-specific blocks after their command boundaries are scanned.
  }
  if active_block.length() != 0 {
    raise AlignNexusError("NEXUS block is missing its END command")
  }
  let body = match matrix_body {
    Some(value) => value
    None =>
      raise AlignNexusError("NEXUS file does not contain an alignment matrix")
  }
  let ntax = if data_ntax >= 0 { data_ntax } else { tax_ntax }
  if ntax <= 0 {
    raise AlignNexusError("NEXUS NTAX must be specified before MATRIX")
  }
  if data_nchar <= 0 {
    raise AlignNexusError("NEXUS NCHAR must be specified before MATRIX")
  }
  if tax_ntax >= 0 && data_ntax >= 0 && tax_ntax != data_ntax {
    raise AlignNexusError("NEXUS TAXA and alignment NTAX values disagree")
  }
  if taxlabels.length() != 0 && taxlabels.length() != ntax {
    raise AlignNexusError("NEXUS TAXLABELS count does not match NTAX")
  }
  align_nexus_validate_format(format)
  let parsed = align_nexus_parse_matrix(
    body, ntax, data_nchar, format, taxlabels,
  )
  let resolved = align_nexus_resolve_and_normalize_rows(parsed.1, format)
  let prepared = align_nexus_remove_all_gap_columns(resolved)
  if prepared.0.length() == 0 || prepared.0[0].length() == 0 {
    raise AlignNexusError("NEXUS matrix contains no non-gap alignment column")
  }
  let metadata = AlignNexusMetadata::create(
    ntax,
    data_nchar,
    format.data_type,
    block_name~,
    missing_character=format.missing_character,
    gap_character=format.gap_character,
    match_character=format.match_character,
    interleaved=format.interleaved,
    respect_case=format.respect_case,
    symbols=format.symbols,
    removed_all_gap_columns=prepared.1,
  )
  let sequences : Array[AlignNexusSequence] = []
  for index = 0; index < ntax; index = index + 1 {
    sequences.push(AlignNexusSequence::{
      id: parsed.0[index],
      sequence: align_nexus_remove_gaps(prepared.0[index]),
      aligned_sequence: prepared.0[index],
    })
  }
  AlignNexusAlignment::create(metadata, sequences)
}

///|
/// Write one canonical NEXUS DATA alignment.
///
/// When `interleave` is omitted, alignments wider than 1000 columns are
/// interleaved, matching Biopython 1.86. Interleaved blocks default to 70
/// columns.
pub fn align_nexus_write(
  alignment : AlignNexusAlignment,
  interleave? : Bool? = None,
  block_width? : Int = 70,
) -> String raise AlignNexusError {
  align_nexus_validate_alignment(alignment)
  let width = alignment.alignment_length()
  if width == 0 {
    raise AlignNexusError("Non-empty NEXUS sequences are required")
  }
  if block_width <= 0 {
    raise AlignNexusError("NEXUS writer block width must be positive")
  }
  let use_interleave = match interleave {
    Some(value) => value
    None => width > 1000
  }
  let names : Array[String] = []
  let mut name_width = 0
  for sequence in alignment.sequences {
    let safe = align_nexus_safe_name(sequence.id)
    names.push(safe)
    if safe.length() > name_width {
      name_width = safe.length()
    }
  }
  let output = StringBuilder::new()
  output.write_string("#NEXUS\n")
  output.write_string("begin data;\n")
  output.write_string(
    "dimensions ntax=" +
    alignment.num_sequences().to_string() +
    " nchar=" +
    width.to_string() +
    ";\n",
  )
  output.write_string(
    "format datatype=" +
    alignment.metadata.data_type.code() +
    " missing=" +
    alignment.metadata.missing_character +
    " gap=" +
    alignment.metadata.gap_character,
  )
  if alignment.metadata.data_type == AlignNexusStandard {
    output.write_string(
      " symbols=\"" +
      align_nexus_escape_double_quotes(alignment.metadata.symbols) +
      "\"",
    )
  }
  if alignment.metadata.respect_case {
    output.write_string(" respectcase")
  }
  if use_interleave {
    output.write_string(" interleave")
  }
  output.write_string(";\n")
  output.write_string("matrix\n")
  if use_interleave {
    let mut start = 0
    while start < width {
      let end = align_nexus_min(start + block_width, width)
      for row = 0; row < alignment.sequences.length(); row = row + 1 {
        output.write_string(align_nexus_pad_right(names[row], name_width + 1))
        output.write_string(
          align_nexus_output_fragment(
            alignment.sequences[row].aligned_sequence,
            start,
            end,
            alignment.metadata.gap_character,
          ),
        )
        output.write_char('\n')
      }
      output.write_char('\n')
      start = end
    }
  } else {
    for row = 0; row < alignment.sequences.length(); row = row + 1 {
      output.write_string(align_nexus_pad_right(names[row], name_width + 1))
      output.write_string(
        align_nexus_output_fragment(
          alignment.sequences[row].aligned_sequence,
          0,
          width,
          alignment.metadata.gap_character,
        ),
      )
      output.write_char('\n')
    }
  }
  output.write_string(";\n")
  output.write_string("end;\n")
  output.to_string()
}

///|
/// Return the canonical lowercase FORMAT datatype name.
pub fn AlignNexusDataType::code(self : AlignNexusDataType) -> String {
  match self {
    AlignNexusDna => "dna"
    AlignNexusRna => "rna"
    AlignNexusProtein => "protein"
    AlignNexusStandard => "standard"
  }
}

///|
/// Return the Biopython-style molecule type annotation.
pub fn AlignNexusDataType::molecule_type(self : AlignNexusDataType) -> String {
  match self {
    AlignNexusDna => "DNA"
    AlignNexusRna => "RNA"
    AlignNexusProtein => "protein"
    AlignNexusStandard => ""
  }
}

///|
/// Return the number of sequence rows.
pub fn AlignNexusAlignment::num_sequences(self : AlignNexusAlignment) -> Int {
  self.sequences.length()
}

///|
/// Return the alignment width after removal of all-gap source columns.
pub fn AlignNexusAlignment::alignment_length(self : AlignNexusAlignment) -> Int {
  if self.sequences.length() == 0 {
    0
  } else {
    self.sequences[0].aligned_sequence.length()
  }
}

///|
/// Return the source width before all-gap columns were removed.
pub fn AlignNexusAlignment::source_alignment_length(
  self : AlignNexusAlignment,
) -> Int {
  self.metadata.declared_characters
}

///|
/// Locate the first row with an exact identifier.
pub fn AlignNexusAlignment::find_sequence(
  self : AlignNexusAlignment,
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
pub fn AlignNexusAlignment::column(
  self : AlignNexusAlignment,
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
/// Map a zero-based ungapped row position to an alignment column.
pub fn AlignNexusAlignment::sequence_position_to_column(
  self : AlignNexusAlignment,
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
/// Map a printed column to a zero-based ungapped row position.
pub fn AlignNexusAlignment::column_to_sequence_position(
  self : AlignNexusAlignment,
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
pub fn AlignNexusAlignment::map_position(
  self : AlignNexusAlignment,
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
pub fn AlignNexusAlignment::aligned_pairs(
  self : AlignNexusAlignment,
  first_row : Int,
  second_row : Int,
) -> Array[(Int?, Int?)] raise AlignNexusError {
  align_nexus_validate_row(self, first_row)
  align_nexus_validate_row(self, second_row)
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
pub fn AlignNexusAlignment::coordinate_path(
  self : AlignNexusAlignment,
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
      align_nexus_movement_changes(self, column, column + 1)
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
pub fn AlignNexusAlignment::pair_counts(
  self : AlignNexusAlignment,
  first_row : Int,
  second_row : Int,
) -> AlignNexusCounts raise AlignNexusError {
  align_nexus_validate_row(self, first_row)
  align_nexus_validate_row(self, second_row)
  align_nexus_count_pair(self, first_row, second_row)
}

///|
/// Aggregate statistics across every unordered pair of rows.
pub fn AlignNexusAlignment::counts(
  self : AlignNexusAlignment,
) -> AlignNexusCounts {
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
      let counts = align_nexus_count_pair(self, first, second)
      pairs = pairs + 1
      aligned = aligned + counts.aligned
      identities = identities + counts.identities
      mismatches = mismatches + counts.mismatches
      gap_columns = gap_columns + counts.gap_columns
      double_gap_columns = double_gap_columns + counts.double_gap_columns
      gap_opens = gap_opens + counts.gap_opens
    }
  }
  AlignNexusCounts::{
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
pub fn AlignNexusCounts::identity(self : AlignNexusCounts) -> Double {
  if self.aligned == 0 {
    0.0
  } else {
    self.identities.to_double() / self.aligned.to_double()
  }
}

///|
/// Calculate a majority consensus. Gaps do not vote.
pub fn AlignNexusAlignment::consensus(
  self : AlignNexusAlignment,
  minimum_fraction? : Double = 0.0,
) -> String raise AlignNexusError {
  if minimum_fraction != minimum_fraction ||
    minimum_fraction < 0.0 ||
    minimum_fraction > 1.0 {
    raise AlignNexusError(
      "NEXUS consensus minimum fraction must be between 0 and 1",
    )
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
    let mut best_code = 0
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
/// Return per-column non-gap occupancy.
pub fn AlignNexusAlignment::occupancy(
  self : AlignNexusAlignment,
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
/// Return a compact alignment summary.
pub fn AlignNexusAlignment::summary(self : AlignNexusAlignment) -> String {
  "AlignNexusAlignment(type=" +
  self.metadata.data_type.code() +
  ", sequences=" +
  self.num_sequences().to_string() +
  ", columns=" +
  self.alignment_length().to_string() +
  ", source_columns=" +
  self.source_alignment_length().to_string() +
  ", removed_all_gap=" +
  self.metadata.removed_all_gap_columns.to_string() +
  ")"
}

///|
/// Return a compact interleaved NEXUS example with quoting and MATCHCHAR.
pub fn align_nexus_example_text() -> String {
  "#NEXUS\n" +
  "[MoonBit NEXUS example with a [nested] comment]\n" +
  "begin data;\n" +
  "dimensions ntax=3 nchar=12;\n" +
  "format datatype=dna missing=? gap=- matchchar=. interleave=yes;\n" +
  "matrix\n" +
  "reference       ACGT--\n" +
  "'query one'     A.GT??\n" +
  "'isn''t_query'  .C-T--\n" +
  "\n" +
  "reference       ACGTAC\n" +
  "query_one       ..G-AC\n" +
  "'isn''t_query'  A-GTAC\n" +
  ";\n" +
  "end;\n"
}

///|
fn align_nexus_validate_alignment(
  alignment : AlignNexusAlignment,
) -> Unit raise AlignNexusError {
  if alignment.sequences.length() == 0 {
    raise AlignNexusError("NEXUS alignment must contain at least one sequence")
  }
  if alignment.sequences.length() != alignment.metadata.declared_taxa {
    raise AlignNexusError("NEXUS sequence count does not match declared NTAX")
  }
  let width = alignment.sequences[0].aligned_sequence.length()
  if width + alignment.metadata.removed_all_gap_columns !=
    alignment.metadata.declared_characters {
    raise AlignNexusError(
      "NEXUS alignment width and removed columns do not match NCHAR",
    )
  }
  align_nexus_validate_format_char(
    alignment.metadata.missing_character,
    "NEXUS missing-data character",
  )
  align_nexus_validate_format_char(
    alignment.metadata.gap_character,
    "NEXUS gap character",
  )
  align_nexus_validate_symbols(
    alignment.metadata.data_type,
    alignment.metadata.symbols,
    alignment.metadata.respect_case,
  )
  for sequence in alignment.sequences {
    align_nexus_validate_id(sequence.id)
    if sequence.aligned_sequence.length() != width {
      raise AlignNexusError("NEXUS aligned rows must have equal widths")
    }
    for index = 0; index < sequence.aligned_sequence.length(); index = index + 1 {
      let code = sequence.aligned_sequence.unsafe_get(index).to_int()
      if code != '-'.to_int() &&
        code != alignment.metadata.missing_character.unsafe_get(0).to_int() &&
        !align_nexus_valid_residue(
          code,
          alignment.metadata.data_type,
          alignment.metadata.symbols,
          alignment.metadata.respect_case,
        ) {
        raise AlignNexusError(
          "NEXUS sequence contains a residue inconsistent with FORMAT",
        )
      }
    }
    if align_nexus_remove_gaps(sequence.aligned_sequence) != sequence.sequence {
      raise AlignNexusError(
        "NEXUS sequence row has inconsistent derived fields",
      )
    }
  }
  for column = 0; column < width; column = column + 1 {
    let mut has_residue = false
    for sequence in alignment.sequences {
      if sequence.aligned_sequence.unsafe_get(column).to_int() != '-'.to_int() {
        has_residue = true
      }
    }
    if !has_residue {
      raise AlignNexusError(
        "NEXUS alignment contains an unremoved all-gap column",
      )
    }
  }
}

///|
fn align_nexus_validate_row(
  alignment : AlignNexusAlignment,
  row : Int,
) -> Unit raise AlignNexusError {
  if row < 0 || row >= alignment.sequences.length() {
    raise AlignNexusError("NEXUS row index is out of bounds")
  }
}

///|
fn align_nexus_count_pair(
  alignment : AlignNexusAlignment,
  first_row : Int,
  second_row : Int,
) -> AlignNexusCounts {
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
  AlignNexusCounts::{
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
fn align_nexus_scan_commands(
  body : String,
) -> Array[String] raise AlignNexusError {
  let commands : Array[String] = []
  let mut current = StringBuilder::new()
  let mut quote = 0
  let mut comment_depth = 0
  let mut index = 0
  while index < body.length() {
    let code = body.unsafe_get(index).to_int()
    if comment_depth > 0 {
      if code == '['.to_int() {
        comment_depth = comment_depth + 1
      } else if code == ']'.to_int() {
        comment_depth = comment_depth - 1
      } else if code == '\n'.to_int() {
        current.write_char('\n')
      }
      index = index + 1
      continue
    }
    if quote != 0 {
      current.write_char(code.unsafe_to_char())
      if code == quote {
        if index + 1 < body.length() &&
          body.unsafe_get(index + 1).to_int() == quote {
          current.write_char(quote.unsafe_to_char())
          index = index + 2
          continue
        }
        quote = 0
      }
      index = index + 1
      continue
    }
    if code == '['.to_int() {
      comment_depth = 1
    } else if code == ']'.to_int() {
      raise AlignNexusError("Unmatched closing NEXUS comment bracket")
    } else if code == '\''.to_int() || code == '"'.to_int() {
      quote = code
      current.write_char(code.unsafe_to_char())
    } else if code == ';'.to_int() {
      let command = current.to_string().trim().to_owned()
      if command.length() > 0 {
        commands.push(command)
      }
      current = StringBuilder::new()
    } else {
      current.write_char(code.unsafe_to_char())
    }
    index = index + 1
  }
  if comment_depth != 0 {
    raise AlignNexusError("Unterminated NEXUS comment")
  }
  if quote != 0 {
    raise AlignNexusError("Unterminated quoted NEXUS token")
  }
  if current.to_string().trim().length() != 0 {
    raise AlignNexusError("NEXUS command is missing a semicolon")
  }
  commands
}

///|
fn align_nexus_split_keyword(command : String) -> (String, String) {
  let mut index = 0
  while index < command.length() &&
        align_nexus_is_whitespace(command.unsafe_get(index).to_int()) {
    index = index + 1
  }
  let start = index
  while index < command.length() &&
        !align_nexus_is_whitespace(command.unsafe_get(index).to_int()) {
    index = index + 1
  }
  (
    command[start:index].to_owned().to_lower(),
    command[index:].to_owned().trim().to_owned(),
  )
}

///|
fn align_nexus_tokenize(text : String) -> Array[String] raise AlignNexusError {
  let tokens : Array[String] = []
  let mut index = 0
  while index < text.length() {
    while index < text.length() &&
          align_nexus_is_whitespace(text.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index >= text.length() {
      break
    }
    let code = text.unsafe_get(index).to_int()
    if code == '='.to_int() || code == ','.to_int() {
      tokens.push(code.unsafe_to_char().to_string())
      index = index + 1
      continue
    }
    if code == '\''.to_int() || code == '"'.to_int() {
      let quote = code
      let output = StringBuilder::new()
      index = index + 1
      let mut closed = false
      while index < text.length() {
        let current = text.unsafe_get(index).to_int()
        if current == quote {
          if index + 1 < text.length() &&
            text.unsafe_get(index + 1).to_int() == quote {
            output.write_char(quote.unsafe_to_char())
            index = index + 2
          } else {
            index = index + 1
            closed = true
            break
          }
        } else {
          output.write_char(current.unsafe_to_char())
          index = index + 1
        }
      }
      if !closed {
        raise AlignNexusError("Unterminated quoted NEXUS token")
      }
      tokens.push(output.to_string())
      continue
    }
    let start = index
    while index < text.length() {
      let current = text.unsafe_get(index).to_int()
      if align_nexus_is_whitespace(current) ||
        current == '='.to_int() ||
        current == ','.to_int() {
        break
      }
      index = index + 1
    }
    tokens.push(text[start:index].to_owned())
  }
  tokens
}

///|
fn align_nexus_parse_dimensions(
  text : String,
) -> (Int, Int) raise AlignNexusError {
  let tokens = align_nexus_tokenize(text)
  let mut ntax = -1
  let mut nchar = -1
  let mut index = 0
  while index < tokens.length() {
    let key = tokens[index].to_lower()
    if key == "newtaxa" {
      index = index + 1
      continue
    }
    if key != "ntax" && key != "nchar" {
      raise AlignNexusError(
        "Unknown NEXUS DIMENSIONS option '" + tokens[index] + "'",
      )
    }
    if index + 2 >= tokens.length() || tokens[index + 1] != "=" {
      raise AlignNexusError("Malformed NEXUS DIMENSIONS assignment")
    }
    let value = align_nexus_parse_positive_int(
      tokens[index + 2],
      "NEXUS " + key.to_upper(),
    )
    if key == "ntax" {
      if ntax >= 0 {
        raise AlignNexusError("Duplicated NEXUS NTAX assignment")
      }
      ntax = value
    } else {
      if nchar >= 0 {
        raise AlignNexusError("Duplicated NEXUS NCHAR assignment")
      }
      nchar = value
    }
    index = index + 3
  }
  (ntax, nchar)
}

///|
fn align_nexus_parse_format(
  text : String,
  initial : AlignNexusFormat,
) -> AlignNexusFormat raise AlignNexusError {
  let tokens = align_nexus_tokenize(text)
  let mut data_type = initial.data_type
  let mut missing = initial.missing_character
  let mut gap = initial.gap_character
  let mut match_character = initial.match_character
  let mut interleaved = initial.interleaved
  let mut respect_case = initial.respect_case
  let mut symbols = initial.symbols
  let mut labels = initial.labels
  let mut index = 0
  while index < tokens.length() {
    let key = tokens[index].to_lower()
    if key == "interleave" {
      if index + 1 < tokens.length() && tokens[index + 1] == "=" {
        if index + 2 >= tokens.length() {
          raise AlignNexusError("Missing NEXUS INTERLEAVE value")
        }
        interleaved = align_nexus_parse_yes_no(
          tokens[index + 2],
          "NEXUS INTERLEAVE",
        )
        index = index + 3
      } else {
        interleaved = true
        index = index + 1
      }
      continue
    }
    if key == "nointerleave" {
      interleaved = false
      index = index + 1
      continue
    }
    if key == "respectcase" {
      respect_case = true
      index = index + 1
      continue
    }
    if key == "labels" {
      if index + 1 < tokens.length() && tokens[index + 1] == "=" {
        if index + 2 >= tokens.length() {
          raise AlignNexusError("Missing NEXUS LABELS value")
        }
        let value = tokens[index + 2].to_lower()
        if value == "no" || value == "none" {
          labels = false
        } else if value == "yes" || value == "left" {
          labels = true
        } else {
          raise AlignNexusError("Unsupported NEXUS LABELS value")
        }
        index = index + 3
      } else {
        labels = true
        index = index + 1
      }
      continue
    }
    if key == "nolabels" {
      labels = false
      index = index + 1
      continue
    }
    if key == "transpose" ||
      key == "tokens" ||
      key == "notokens" ||
      key == "equate" {
      raise AlignNexusError(
        "NEXUS FORMAT option '" + key + "' is not supported for alignments",
      )
    }
    if key != "datatype" &&
      key != "missing" &&
      key != "gap" &&
      key != "matchchar" &&
      key != "symbols" {
      raise AlignNexusError(
        "Unknown NEXUS FORMAT option '" + tokens[index] + "'",
      )
    }
    if index + 2 >= tokens.length() || tokens[index + 1] != "=" {
      raise AlignNexusError(
        "Malformed NEXUS FORMAT assignment for " + key.to_upper(),
      )
    }
    let value = tokens[index + 2]
    if key == "datatype" {
      data_type = match value.to_lower() {
        "dna" | "nucleotide" => AlignNexusDna
        "rna" => AlignNexusRna
        "protein" => AlignNexusProtein
        "standard" => AlignNexusStandard
        other =>
          raise AlignNexusError("Unsupported NEXUS datatype '" + other + "'")
      }
    } else if key == "missing" {
      missing = value
    } else if key == "gap" {
      gap = value
    } else if key == "matchchar" {
      match_character = Some(value)
    } else {
      symbols = value
    }
    index = index + 3
  }
  let result = AlignNexusFormat::{
    data_type,
    missing_character: missing,
    gap_character: gap,
    match_character,
    interleaved,
    respect_case,
    symbols,
    labels,
  }
  align_nexus_validate_format(result)
  result
}

///|
fn align_nexus_validate_format(
  format : AlignNexusFormat,
) -> Unit raise AlignNexusError {
  align_nexus_validate_format_char(
    format.missing_character,
    "NEXUS missing-data character",
  )
  align_nexus_validate_format_char(format.gap_character, "NEXUS gap character")
  if format.missing_character == format.gap_character {
    raise AlignNexusError(
      "NEXUS missing-data and gap characters must be different",
    )
  }
  match format.match_character {
    Some(value) => {
      align_nexus_validate_format_char(value, "NEXUS match character")
      if value == format.missing_character || value == format.gap_character {
        raise AlignNexusError(
          "NEXUS match character must differ from missing and gap characters",
        )
      }
    }
    None => ()
  }
  align_nexus_validate_symbols(
    format.data_type,
    format.symbols,
    format.respect_case,
  )
}

///|
fn align_nexus_parse_taxlabels(
  text : String,
) -> Array[String] raise AlignNexusError {
  let tokens = align_nexus_tokenize(text)
  let labels : Array[String] = []
  for token in tokens {
    if token != "," {
      align_nexus_validate_id(token)
      labels.push(token)
    }
  }
  if labels.length() == 0 {
    raise AlignNexusError("NEXUS TAXLABELS must not be empty")
  }
  labels
}

///|
fn align_nexus_parse_matrix(
  body : String,
  ntax : Int,
  nchar : Int,
  format : AlignNexusFormat,
  taxlabels : Array[String],
) -> (Array[String], Array[String]) raise AlignNexusError {
  let raw_lines = align_nexus_normalize_lines(body)
  let lines : Array[String] = []
  for raw in raw_lines {
    let trimmed = raw.trim().to_owned()
    if trimmed.length() > 0 {
      lines.push(trimmed)
    }
  }
  if lines.length() == 0 {
    raise AlignNexusError("NEXUS MATRIX must not be empty")
  }
  if format.interleaved {
    align_nexus_parse_interleaved_matrix(
      lines,
      ntax,
      nchar,
      format.labels,
      taxlabels,
    )
  } else {
    align_nexus_parse_sequential_matrix(
      lines,
      ntax,
      nchar,
      format.labels,
      taxlabels,
    )
  }
}

///|
fn align_nexus_parse_interleaved_matrix(
  lines : Array[String],
  ntax : Int,
  nchar : Int,
  labels : Bool,
  taxlabels : Array[String],
) -> (Array[String], Array[String]) raise AlignNexusError {
  if !labels && taxlabels.length() != ntax {
    raise AlignNexusError("NEXUS LABELS=NO matrices require NTAX TAXLABELS")
  }
  let entry_ids : Array[String] = []
  let entry_chunks : Array[String] = []
  let mut line_index = 0
  while line_index < lines.length() {
    if labels {
      let parsed = align_nexus_parse_labeled_line(lines[line_index])
      let mut chunk = align_nexus_remove_whitespace(parsed.1)
      line_index = line_index + 1
      if chunk.length() == 0 {
        if line_index >= lines.length() {
          raise AlignNexusError("Missing sequence after NEXUS matrix label")
        }
        chunk = align_nexus_remove_whitespace(lines[line_index])
        line_index = line_index + 1
      }
      align_nexus_validate_id(parsed.0)
      entry_ids.push(parsed.0)
      entry_chunks.push(chunk)
    } else {
      entry_ids.push(taxlabels[entry_ids.length() % ntax])
      entry_chunks.push(align_nexus_remove_whitespace(lines[line_index]))
      line_index = line_index + 1
    }
  }
  if entry_ids.length() == 0 || entry_ids.length() % ntax != 0 {
    raise AlignNexusError(
      "Interleaved NEXUS MATRIX does not contain complete NTAX blocks",
    )
  }
  let ids : Array[String] = []
  let rows = Array::make(ntax, "")
  let block_count = entry_ids.length() / ntax
  for block = 0; block < block_count; block = block + 1 {
    let first_chunk_length = entry_chunks[block * ntax].length()
    if first_chunk_length == 0 {
      raise AlignNexusError("NEXUS interleaved matrix chunk must not be empty")
    }
    for row = 0; row < ntax; row = row + 1 {
      let entry = block * ntax + row
      if entry_chunks[entry].length() != first_chunk_length {
        raise AlignNexusError(
          "Rows in one NEXUS interleaved block must have equal widths",
        )
      }
      if block == 0 {
        ids.push(entry_ids[entry])
      } else if !align_nexus_labels_equivalent(ids[row], entry_ids[entry]) {
        raise AlignNexusError(
          "NEXUS interleaved taxon order or label changed between blocks",
        )
      }
      rows[row] = rows[row] + entry_chunks[entry]
    }
  }
  for row = 0; row < ntax; row = row + 1 {
    if rows[row].length() != nchar {
      raise AlignNexusError(
        "NEXUS NCHAR does not match interleaved row length for '" +
        ids[row] +
        "'",
      )
    }
  }
  (ids, rows)
}

///|
fn align_nexus_parse_sequential_matrix(
  lines : Array[String],
  ntax : Int,
  nchar : Int,
  labels : Bool,
  taxlabels : Array[String],
) -> (Array[String], Array[String]) raise AlignNexusError {
  if !labels && taxlabels.length() != ntax {
    raise AlignNexusError("NEXUS LABELS=NO matrices require NTAX TAXLABELS")
  }
  let ids : Array[String] = []
  let rows : Array[String] = []
  let mut line_index = 0
  for row = 0; row < ntax; row = row + 1 {
    if line_index >= lines.length() {
      raise AlignNexusError("Not enough taxa in NEXUS MATRIX")
    }
    let mut sequence = ""
    if labels {
      let parsed = align_nexus_parse_labeled_line(lines[line_index])
      align_nexus_validate_id(parsed.0)
      ids.push(parsed.0)
      sequence = align_nexus_remove_whitespace(parsed.1)
    } else {
      ids.push(taxlabels[row])
      sequence = align_nexus_remove_whitespace(lines[line_index])
    }
    line_index = line_index + 1
    while sequence.length() < nchar {
      if line_index >= lines.length() {
        raise AlignNexusError("NEXUS MATRIX row is shorter than declared NCHAR")
      }
      sequence = sequence + align_nexus_remove_whitespace(lines[line_index])
      line_index = line_index + 1
    }
    if sequence.length() != nchar {
      raise AlignNexusError("NEXUS MATRIX row is longer than declared NCHAR")
    }
    rows.push(sequence)
  }
  if line_index != lines.length() {
    raise AlignNexusError("Too many taxa or sequence lines in NEXUS MATRIX")
  }
  (ids, rows)
}

///|
fn align_nexus_parse_labeled_line(
  line : String,
) -> (String, String) raise AlignNexusError {
  let mut index = 0
  while index < line.length() &&
        align_nexus_is_whitespace(line.unsafe_get(index).to_int()) {
    index = index + 1
  }
  if index >= line.length() {
    raise AlignNexusError("Empty NEXUS MATRIX row")
  }
  let code = line.unsafe_get(index).to_int()
  if code == '\''.to_int() || code == '"'.to_int() {
    let quote = code
    let label = StringBuilder::new()
    index = index + 1
    let mut closed = false
    while index < line.length() {
      let current = line.unsafe_get(index).to_int()
      if current == quote {
        if index + 1 < line.length() &&
          line.unsafe_get(index + 1).to_int() == quote {
          label.write_char(quote.unsafe_to_char())
          index = index + 2
        } else {
          index = index + 1
          closed = true
          break
        }
      } else {
        label.write_char(current.unsafe_to_char())
        index = index + 1
      }
    }
    if !closed {
      raise AlignNexusError("Unterminated quoted NEXUS taxon name")
    }
    if index < line.length() &&
      !align_nexus_is_whitespace(line.unsafe_get(index).to_int()) {
      raise AlignNexusError(
        "Quoted NEXUS taxon name must be followed by whitespace",
      )
    }
    (label.to_string(), line[index:].to_owned().trim().to_owned())
  } else {
    let start = index
    while index < line.length() &&
          !align_nexus_is_whitespace(line.unsafe_get(index).to_int()) {
      index = index + 1
    }
    (line[start:index].to_owned(), line[index:].to_owned().trim().to_owned())
  }
}

///|
fn align_nexus_resolve_and_normalize_rows(
  rows : Array[String],
  format : AlignNexusFormat,
) -> Array[String] raise AlignNexusError {
  let resolved : Array[String] = []
  for row_index = 0; row_index < rows.length(); row_index = row_index + 1 {
    let output = StringBuilder::new(size_hint=rows[row_index].length())
    for column = 0; column < rows[row_index].length(); column = column + 1 {
      let mut code = rows[row_index].unsafe_get(column).to_int()
      match format.match_character {
        Some(character) =>
          if code == character.unsafe_get(0).to_int() {
            if row_index == 0 {
              raise AlignNexusError(
                "NEXUS match character is not allowed in the first row",
              )
            }
            code = resolved[0].unsafe_get(column).to_int()
          }
        None => ()
      }
      output.write_char(code.unsafe_to_char())
    }
    resolved.push(
      align_nexus_normalize_row(
        output.to_string(),
        format.data_type,
        format.missing_character,
        format.gap_character,
        format.symbols,
        format.respect_case,
      ),
    )
  }
  resolved
}

///|
fn align_nexus_normalize_row(
  row : String,
  data_type : AlignNexusDataType,
  missing_character : String,
  gap_character : String,
  symbols : String,
  respect_case : Bool,
) -> String raise AlignNexusError {
  let output = StringBuilder::new(size_hint=row.length())
  let missing_code = missing_character.unsafe_get(0).to_int()
  let gap_code = gap_character.unsafe_get(0).to_int()
  for index = 0; index < row.length(); index = index + 1 {
    let code = row.unsafe_get(index).to_int()
    if code == gap_code {
      output.write_char('-')
    } else if code == missing_code ||
      align_nexus_valid_residue(code, data_type, symbols, respect_case) {
      output.write_char(code.unsafe_to_char())
    } else {
      raise AlignNexusError(
        "Illegal " +
        data_type.code() +
        " character '" +
        code.unsafe_to_char().to_string() +
        "' in NEXUS MATRIX",
      )
    }
  }
  output.to_string()
}

///|
fn align_nexus_valid_residue(
  code : Int,
  data_type : AlignNexusDataType,
  symbols : String,
  respect_case : Bool,
) -> Bool {
  match data_type {
    AlignNexusDna =>
      align_nexus_contains_code("ACGTRYSWKMBDHVNacgtryswkmbdhvn", code)
    AlignNexusRna =>
      align_nexus_contains_code("ACGURYSWKMBDHVNacguryswkmbdhvn", code)
    AlignNexusProtein =>
      align_nexus_contains_code(
        "ACDEFGHIKLMNPQRSTVWYBZX*acdefghiklmnpqrstvywbzx", code,
      )
    AlignNexusStandard =>
      if align_nexus_contains_code(symbols, code) {
        true
      } else if !respect_case {
        align_nexus_contains_code(symbols, align_nexus_swap_case(code))
      } else {
        false
      }
  }
}

///|
fn align_nexus_validate_symbols(
  data_type : AlignNexusDataType,
  symbols : String,
  respect_case : Bool,
) -> Unit raise AlignNexusError {
  if data_type == AlignNexusStandard && symbols.length() == 0 {
    raise AlignNexusError(
      "NEXUS standard datatype requires a non-empty SYMBOLS value",
    )
  }
  let seen : Array[Int] = []
  for index = 0; index < symbols.length(); index = index + 1 {
    let code = symbols.unsafe_get(index).to_int()
    if code < 33 ||
      code > 126 ||
      align_nexus_is_whitespace(code) ||
      code == '\''.to_int() ||
      code == '"'.to_int() ||
      code == '['.to_int() ||
      code == ']'.to_int() ||
      code == '('.to_int() ||
      code == ')'.to_int() ||
      code == ','.to_int() ||
      code == ';'.to_int() ||
      code == '='.to_int() {
      raise AlignNexusError(
        "NEXUS SYMBOLS must contain distinct printable state characters",
      )
    }
    let comparable = if respect_case {
      code
    } else {
      align_nexus_upper_code(code)
    }
    for previous in seen {
      if previous == comparable {
        raise AlignNexusError("NEXUS SYMBOLS contains a duplicate state")
      }
    }
    seen.push(comparable)
  }
}

///|
fn align_nexus_validate_format_char(
  value : String,
  label : String,
) -> Unit raise AlignNexusError {
  if value.length() != 1 {
    raise AlignNexusError(label + " must be one ASCII character")
  }
  let code = value.unsafe_get(0).to_int()
  if code < 33 ||
    code > 126 ||
    code == '['.to_int() ||
    code == ']'.to_int() ||
    code == '\''.to_int() ||
    code == '"'.to_int() ||
    code == ';'.to_int() {
    raise AlignNexusError(
      label + " must be one printable non-reserved character",
    )
  }
}

///|
fn align_nexus_validate_id(id : String) -> Unit raise AlignNexusError {
  if id.length() == 0 {
    raise AlignNexusError("NEXUS taxon identifier must not be empty")
  }
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    if code < 32 ||
      code == 127 ||
      code == '\n'.to_int() ||
      code == '\r'.to_int() {
      raise AlignNexusError(
        "NEXUS taxon identifier contains a control character",
      )
    }
  }
}

///|
fn align_nexus_remove_all_gap_columns(
  rows : Array[String],
) -> (Array[String], Int) {
  if rows.length() == 0 {
    return ([], 0)
  }
  let width = rows[0].length()
  let keep = Array::make(width, false)
  let mut removed = 0
  for column = 0; column < width; column = column + 1 {
    for row in rows {
      if row.unsafe_get(column).to_int() != '-'.to_int() {
        keep[column] = true
      }
    }
    if !keep[column] {
      removed = removed + 1
    }
  }
  let result : Array[String] = []
  for row in rows {
    let output = StringBuilder::new(size_hint=width - removed)
    for column = 0; column < width; column = column + 1 {
      if keep[column] {
        output.write_char(row.unsafe_get(column).unsafe_to_char())
      }
    }
    result.push(output.to_string())
  }
  (result, removed)
}

///|
fn align_nexus_movement_changes(
  alignment : AlignNexusAlignment,
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
fn align_nexus_labels_equivalent(first : String, second : String) -> Bool {
  if first == second {
    return true
  }
  let left = StringBuilder::new(size_hint=first.length())
  let right = StringBuilder::new(size_hint=second.length())
  for index = 0; index < first.length(); index = index + 1 {
    let code = first.unsafe_get(index).to_int()
    left.write_char(
      (if code == ' '.to_int() { '_'.to_int() } else { code }).unsafe_to_char(),
    )
  }
  for index = 0; index < second.length(); index = index + 1 {
    let code = second.unsafe_get(index).to_int()
    right.write_char(
      (if code == ' '.to_int() { '_'.to_int() } else { code }).unsafe_to_char(),
    )
  }
  left.to_string() == right.to_string()
}

///|
fn align_nexus_safe_name(id : String) -> String {
  let mut quote = false
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    let safe = (code >= 'A'.to_int() && code <= 'Z'.to_int()) ||
      (code >= 'a'.to_int() && code <= 'z'.to_int()) ||
      (code >= '0'.to_int() && code <= '9'.to_int()) ||
      code == '_'.to_int() ||
      code == '-'.to_int() ||
      code == '.'.to_int()
    if !safe {
      quote = true
    }
  }
  if !quote {
    return id
  }
  let output = StringBuilder::new(size_hint=id.length() + 2)
  output.write_char('\'')
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    output.write_char(code.unsafe_to_char())
    if code == '\''.to_int() {
      output.write_char('\'')
    }
  }
  output.write_char('\'')
  output.to_string()
}

///|
fn align_nexus_output_fragment(
  row : String,
  start : Int,
  end : Int,
  gap_character : String,
) -> String {
  let output = StringBuilder::new(size_hint=end - start)
  let gap_code = gap_character.unsafe_get(0).to_int()
  for index = start; index < end; index = index + 1 {
    let code = row.unsafe_get(index).to_int()
    output.write_char(
      (if code == '-'.to_int() { gap_code } else { code }).unsafe_to_char(),
    )
  }
  output.to_string()
}

///|
fn align_nexus_escape_double_quotes(value : String) -> String {
  let output = StringBuilder::new(size_hint=value.length())
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    output.write_char(code.unsafe_to_char())
    if code == '"'.to_int() {
      output.write_char('"')
    }
  }
  output.to_string()
}

///|
fn align_nexus_remove_gaps(sequence : String) -> String {
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
fn align_nexus_remove_whitespace(value : String) -> String {
  let output = StringBuilder::new(size_hint=value.length())
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    if !align_nexus_is_whitespace(code) {
      output.write_char(code.unsafe_to_char())
    }
  }
  output.to_string()
}

///|
fn align_nexus_normalize_lines(text : String) -> Array[String] {
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
fn align_nexus_parse_positive_int(
  text : String,
  label : String,
) -> Int raise AlignNexusError {
  if text.length() == 0 {
    raise AlignNexusError(label + " is empty")
  }
  let mut value = 0
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise AlignNexusError(label + " must be a positive integer")
    }
    let digit = code - '0'.to_int()
    if value > 214748364 || (value == 214748364 && digit > 7) {
      raise AlignNexusError(label + " is outside the supported integer range")
    }
    value = value * 10 + digit
  }
  if value <= 0 {
    raise AlignNexusError(label + " must be positive")
  }
  value
}

///|
fn align_nexus_parse_yes_no(
  value : String,
  label : String,
) -> Bool raise AlignNexusError {
  match value.to_lower() {
    "yes" | "true" => true
    "no" | "false" => false
    _ => raise AlignNexusError(label + " must be YES or NO")
  }
}

///|
fn align_nexus_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\n'.to_int() ||
  code == '\r'.to_int()
}

///|
fn align_nexus_contains_code(values : String, target : Int) -> Bool {
  for index = 0; index < values.length(); index = index + 1 {
    if values.unsafe_get(index).to_int() == target {
      return true
    }
  }
  false
}

///|
fn align_nexus_upper_code(code : Int) -> Int {
  if code >= 'a'.to_int() && code <= 'z'.to_int() {
    code - 'a'.to_int() + 'A'.to_int()
  } else {
    code
  }
}

///|
fn align_nexus_swap_case(code : Int) -> Int {
  if code >= 'a'.to_int() && code <= 'z'.to_int() {
    code - 'a'.to_int() + 'A'.to_int()
  } else if code >= 'A'.to_int() && code <= 'Z'.to_int() {
    code - 'A'.to_int() + 'a'.to_int()
  } else {
    code
  }
}

///|
fn align_nexus_pad_right(value : String, width : Int) -> String {
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
fn align_nexus_min(first : Int, second : Int) -> Int {
  if first < second {
    first
  } else {
    second
  }
}
