// Alignment-aware Exonerate cigar/vulgar support.
//
// This module follows Biopython 1.86 Bio.Align.exonerate semantics. Exonerate
// coordinates are zero-based boundaries. Operations keep query and target
// steps separately so DNA/protein 3:1 alignments do not lose information.

///|
/// Error raised for malformed alignment-aware Exonerate data.
pub suberror AlignExonerateError {
  AlignExonerateError(String)
}

///|
/// Orientation or molecule type reported for an Exonerate sequence.
pub(all) enum AlignExonerateStrand {
  AlignExonerateForward
  AlignExonerateReverse
  AlignExonerateProtein
} derive(Eq, Debug)

///|
/// File-level Exonerate metadata.
pub struct AlignExonerateMetadata {
  program : String
  command_line : String
  hostname : String
} derive(Eq, Debug)

///|
/// One normalized Exonerate path operation.
///
/// Codes use alignment-aware semantics:
/// `M` match, `5`/`3` splice sites, `N` intron, `C` codon, `D` deletion,
/// `I` insertion, `U` non-equivalenced region, `S` split codon, and `F`
/// frame shift. Steps are non-negative molecular coordinates.
pub struct AlignExonerateOperation {
  code : String
  query_step : Int
  target_step : Int
} derive(Eq, Debug)

///|
/// One pairwise alignment from an Exonerate report.
pub struct AlignExonerateAlignment {
  query_id : String
  query_start : Int
  query_end : Int
  query_strand : AlignExonerateStrand
  target_id : String
  target_start : Int
  target_end : Int
  target_strand : AlignExonerateStrand
  score : Double
  operations : Array[AlignExonerateOperation]
} derive(Eq, Debug)

///|
/// A complete Exonerate report.
pub struct AlignExonerateDocument {
  metadata : AlignExonerateMetadata
  alignments : Array[AlignExonerateAlignment]
} derive(Eq, Debug)

///|
/// Operation and coordinate statistics for one Exonerate alignment.
pub struct AlignExonerateCounts {
  operations : Int
  matching_operations : Int
  aligned_query_units : Int
  aligned_target_units : Int
  insertions : Int
  deletions : Int
  insertion_opens : Int
  deletion_opens : Int
  introns : Int
  intron_units : Int
  non_equivalenced_units : Int
  split_codons : Int
  frameshifts : Int
} derive(Eq, Debug)

///|
/// Construct Exonerate file metadata.
pub fn AlignExonerateMetadata::create(
  command_line? : String = "",
  hostname? : String = "",
) -> AlignExonerateMetadata raise AlignExonerateError {
  align_exonerate_validate_header_value(command_line, "Command line")
  align_exonerate_validate_header_value(hostname, "Hostname")
  AlignExonerateMetadata::{ program: "exonerate", command_line, hostname, }
}

///|
/// Construct and validate one normalized Exonerate operation.
pub fn AlignExonerateOperation::create(
  code : String,
  query_step : Int,
  target_step : Int,
) -> AlignExonerateOperation raise AlignExonerateError {
  align_exonerate_validate_operation(code, query_step, target_step)
  AlignExonerateOperation::{ code, query_step, target_step, }
}

///|
/// Construct one alignment and validate operation totals against its bounds.
pub fn AlignExonerateAlignment::create(
  query_id : String,
  query_start : Int,
  query_end : Int,
  query_strand : AlignExonerateStrand,
  target_id : String,
  target_start : Int,
  target_end : Int,
  target_strand : AlignExonerateStrand,
  score : Double,
  operations : Array[AlignExonerateOperation],
) -> AlignExonerateAlignment raise AlignExonerateError {
  align_exonerate_validate_id(query_id, "query")
  align_exonerate_validate_id(target_id, "target")
  align_exonerate_validate_axis(query_start, query_end, query_strand, "query")
  align_exonerate_validate_axis(
    target_start, target_end, target_strand, "target",
  )
  if score.abs() > 1.0e300 {
    raise AlignExonerateError("Exonerate score must be finite")
  }
  if operations.length() == 0 {
    raise AlignExonerateError("Exonerate alignment must contain an operation")
  }
  let copied : Array[AlignExonerateOperation] = []
  let mut query_span = 0
  let mut target_span = 0
  for operation in operations {
    align_exonerate_validate_operation(
      operation.code,
      operation.query_step,
      operation.target_step,
    )
    query_span = query_span + operation.query_step
    target_span = target_span + operation.target_step
    copied.push(operation)
  }
  if query_span != align_exonerate_abs(query_end - query_start) {
    raise AlignExonerateError(
      "Exonerate query operation span does not match its coordinates",
    )
  }
  if target_span != align_exonerate_abs(target_end - target_start) {
    raise AlignExonerateError(
      "Exonerate target operation span does not match its coordinates",
    )
  }
  AlignExonerateAlignment::{
    query_id,
    query_start,
    query_end,
    query_strand,
    target_id,
    target_start,
    target_end,
    target_strand,
    score,
    operations: copied,
  }
}

///|
/// Construct a complete report. Exonerate reports with no hits are valid.
pub fn AlignExonerateDocument::create(
  metadata : AlignExonerateMetadata,
  alignments : Array[AlignExonerateAlignment],
) -> AlignExonerateDocument raise AlignExonerateError {
  if metadata.program != "exonerate" {
    raise AlignExonerateError("Exonerate metadata program must be exonerate")
  }
  let copied : Array[AlignExonerateAlignment] = []
  for alignment in alignments {
    let validated = AlignExonerateAlignment::create(
      alignment.query_id,
      alignment.query_start,
      alignment.query_end,
      alignment.query_strand,
      alignment.target_id,
      alignment.target_start,
      alignment.target_end,
      alignment.target_strand,
      alignment.score,
      alignment.operations,
    )
    copied.push(validated)
  }
  AlignExonerateDocument::{ metadata, alignments: copied, }
}

///|
pub fn AlignExonerateDocument::num_alignments(
  self : AlignExonerateDocument,
) -> Int {
  self.alignments.length()
}

///|
pub fn AlignExonerateDocument::get(
  self : AlignExonerateDocument,
  index : Int,
) -> AlignExonerateAlignment? {
  if index < 0 || index >= self.alignments.length() {
    None
  } else {
    Some(self.alignments[index])
  }
}

///|
pub fn AlignExonerateDocument::find_by_query(
  self : AlignExonerateDocument,
  query_id : String,
) -> Array[AlignExonerateAlignment] {
  let result : Array[AlignExonerateAlignment] = []
  for alignment in self.alignments {
    if alignment.query_id == query_id {
      result.push(alignment)
    }
  }
  result
}

///|
pub fn AlignExonerateDocument::best_alignment(
  self : AlignExonerateDocument,
) -> AlignExonerateAlignment? {
  if self.alignments.length() == 0 {
    return None
  }
  let mut best = self.alignments[0]
  for index = 1; index < self.alignments.length(); index = index + 1 {
    if self.alignments[index].score > best.score {
      best = self.alignments[index]
    }
  }
  Some(best)
}

///|
pub fn AlignExonerateDocument::summary(self : AlignExonerateDocument) -> String {
  "Exonerate report(host=" +
  self.metadata.hostname +
  ", alignments=" +
  self.alignments.length().to_string() +
  ")"
}

///|
pub fn AlignExonerateAlignment::query_length(
  self : AlignExonerateAlignment,
) -> Int {
  if self.query_strand == AlignExonerateReverse {
    self.query_start
  } else {
    self.query_end
  }
}

///|
pub fn AlignExonerateAlignment::target_length(
  self : AlignExonerateAlignment,
) -> Int {
  if self.target_strand == AlignExonerateReverse {
    self.target_start
  } else {
    self.target_end
  }
}

///|
pub fn AlignExonerateAlignment::num_operations(
  self : AlignExonerateAlignment,
) -> Int {
  self.operations.length()
}

///|
pub fn AlignExonerateAlignment::operation_codes(
  self : AlignExonerateAlignment,
) -> String {
  let output = StringBuilder::new(size_hint=self.operations.length())
  for operation in self.operations {
    output.write_string(operation.code)
  }
  output.to_string()
}

///|
/// Return `[target_coordinates, query_coordinates]` boundary paths.
pub fn AlignExonerateAlignment::coordinate_path(
  self : AlignExonerateAlignment,
) -> Array[Array[Int]] {
  let target = [self.target_start]
  let query = [self.query_start]
  let target_sign = align_exonerate_strand_sign(self.target_strand)
  let query_sign = align_exonerate_strand_sign(self.query_strand)
  let mut target_position = self.target_start
  let mut query_position = self.query_start
  for operation in self.operations {
    target_position = target_position + target_sign * operation.target_step
    query_position = query_position + query_sign * operation.query_step
    target.push(target_position)
    query.push(query_position)
  }
  [target, query]
}

///|
/// Map an absolute query residue to the corresponding target residue.
///
/// A target gap returns `None`. For translated alignments, a protein residue
/// maps to the first nucleotide of its codon and each nucleotide maps to its
/// containing protein residue.
pub fn AlignExonerateAlignment::query_position_to_target(
  self : AlignExonerateAlignment,
  position : Int,
) -> Int? {
  align_exonerate_map_position(self, position, true)
}

///|
/// Map an absolute target residue to the corresponding query residue.
pub fn AlignExonerateAlignment::target_position_to_query(
  self : AlignExonerateAlignment,
  position : Int,
) -> Int? {
  align_exonerate_map_position(self, position, false)
}

///|
/// Return target/query residue pairs for all query units in aligned segments.
pub fn AlignExonerateAlignment::aligned_pairs(
  self : AlignExonerateAlignment,
) -> Array[(Int, Int)] {
  let pairs : Array[(Int, Int)] = []
  let query_sign = align_exonerate_strand_sign(self.query_strand)
  let target_sign = align_exonerate_strand_sign(self.target_strand)
  let mut query_boundary = self.query_start
  let mut target_boundary = self.target_start
  for operation in self.operations {
    if operation.query_step > 0 && operation.target_step > 0 {
      for offset = 0; offset < operation.query_step; offset = offset + 1 {
        let query_position = if query_sign > 0 {
          query_boundary + offset
        } else {
          query_boundary - 1 - offset
        }
        let target_offset = offset *
          operation.target_step /
          operation.query_step
        let target_position = if target_sign > 0 {
          target_boundary + target_offset
        } else {
          target_boundary - 1 - target_offset
        }
        pairs.push((target_position, query_position))
      }
    }
    query_boundary = query_boundary + query_sign * operation.query_step
    target_boundary = target_boundary + target_sign * operation.target_step
  }
  pairs
}

///|
/// Calculate operation-aware alignment statistics.
pub fn AlignExonerateAlignment::counts(
  self : AlignExonerateAlignment,
) -> AlignExonerateCounts {
  let mut matching_operations = 0
  let mut aligned_query_units = 0
  let mut aligned_target_units = 0
  let mut insertions = 0
  let mut deletions = 0
  let mut insertion_opens = 0
  let mut deletion_opens = 0
  let mut introns = 0
  let mut intron_units = 0
  let mut non_equivalenced_units = 0
  let mut split_codons = 0
  let mut frameshifts = 0
  let mut previous_insertion = false
  let mut previous_deletion = false
  for operation in self.operations {
    if operation.query_step > 0 && operation.target_step > 0 {
      aligned_query_units = aligned_query_units + operation.query_step
      aligned_target_units = aligned_target_units + operation.target_step
      if operation.code == "M" || operation.code == "C" {
        matching_operations = matching_operations + 1
      }
    }
    let insertion = operation.target_step == 0 && operation.query_step > 0
    let deletion = operation.query_step == 0 && operation.target_step > 0
    if insertion {
      insertions = insertions + operation.query_step
      if !previous_insertion {
        insertion_opens = insertion_opens + 1
      }
    }
    if deletion {
      deletions = deletions + operation.target_step
      if !previous_deletion {
        deletion_opens = deletion_opens + 1
      }
    }
    if operation.code == "N" {
      introns = introns + 1
      intron_units = intron_units + operation.query_step + operation.target_step
    } else if operation.code == "U" {
      non_equivalenced_units = non_equivalenced_units +
        operation.query_step +
        operation.target_step
    } else if operation.code == "S" {
      split_codons = split_codons + 1
    } else if operation.code == "F" {
      frameshifts = frameshifts + 1
    }
    previous_insertion = insertion
    previous_deletion = deletion
  }
  AlignExonerateCounts::{
    operations: self.operations.length(),
    matching_operations,
    aligned_query_units,
    aligned_target_units,
    insertions,
    deletions,
    insertion_opens,
    deletion_opens,
    introns,
    intron_units,
    non_equivalenced_units,
    split_codons,
    frameshifts,
  }
}

///|
pub fn AlignExonerateAlignment::summary(
  self : AlignExonerateAlignment,
) -> String {
  "Exonerate alignment(query=" +
  self.query_id +
  ", target=" +
  self.target_id +
  ", score=" +
  self.score.to_string() +
  ", operations=" +
  self.operations.length().to_string() +
  ")"
}

///|
/// Parse an alignment-aware Exonerate cigar/vulgar report.
pub fn align_exonerate_parse(
  text : String,
) -> AlignExonerateDocument raise AlignExonerateError {
  let lines = align_exonerate_lines(text)
  if lines.length() < 3 {
    raise AlignExonerateError("Truncated Exonerate report header")
  }
  let command_line = align_exonerate_parse_header_line(
    lines[0],
    "Command line: ",
    "Command line",
    1,
  )
  let hostname = align_exonerate_parse_header_line(
    lines[1],
    "Hostname: ",
    "Hostname",
    2,
  )
  let metadata = AlignExonerateMetadata::create(command_line~, hostname~)
  let alignments : Array[AlignExonerateAlignment] = []
  let mut footer_seen = false
  for index = 2; index < lines.length(); index = index + 1 {
    let trimmed = lines[index].trim().to_owned()
    if footer_seen {
      if trimmed.length() > 0 {
        raise AlignExonerateError(
          "Additional data after Exonerate completion marker at line " +
          (index + 1).to_string(),
        )
      }
      continue
    }
    if trimmed.length() == 0 {
      continue
    }
    if trimmed == "-- completed exonerate analysis" {
      footer_seen = true
    } else if trimmed.has_prefix("vulgar: ") {
      alignments.push(
        align_exonerate_parse_vulgar_line(trimmed[8:].to_owned(), index + 1),
      )
    } else if trimmed.has_prefix("cigar: ") {
      alignments.push(
        align_exonerate_parse_cigar_line(trimmed[7:].to_owned(), index + 1),
      )
    } else {
      raise AlignExonerateError(
        "Unexpected Exonerate body line at line " + (index + 1).to_string(),
      )
    }
  }
  if !footer_seen {
    raise AlignExonerateError(
      "Failed to find completed Exonerate analysis marker",
    )
  }
  AlignExonerateDocument::create(metadata, alignments)
}

///|
/// Serialize a complete report as canonical `vulgar` or `cigar` output.
pub fn align_exonerate_write(
  document : AlignExonerateDocument,
  format? : String = "vulgar",
) -> String raise AlignExonerateError {
  if format != "vulgar" && format != "cigar" {
    raise AlignExonerateError("Exonerate output format must be vulgar or cigar")
  }
  let validated = AlignExonerateDocument::create(
    document.metadata,
    document.alignments,
  )
  let output = StringBuilder::new()
  output.write_string(
    "Command line: [" + validated.metadata.command_line + "]\n",
  )
  output.write_string("Hostname: [" + validated.metadata.hostname + "]\n")
  for alignment in validated.alignments {
    if format == "vulgar" {
      output.write_string(align_exonerate_format_vulgar(alignment))
    } else {
      output.write_string(align_exonerate_format_cigar(alignment))
    }
  }
  output.write_string("-- completed exonerate analysis\n")
  output.to_string()
}

///|
fn align_exonerate_parse_vulgar_line(
  content : String,
  line_number : Int,
) -> AlignExonerateAlignment raise AlignExonerateError {
  let words = align_exonerate_split_whitespace(content)
  if words.length() < 12 || (words.length() - 9) % 3 != 0 {
    raise AlignExonerateError(
      "Malformed Exonerate vulgar field count at line " +
      line_number.to_string(),
    )
  }
  let query_id = words[0]
  let query_start = align_exonerate_parse_int(
    words[1],
    "query start",
    line_number,
  )
  let query_end = align_exonerate_parse_int(words[2], "query end", line_number)
  let query_strand = align_exonerate_parse_strand(
    words[3],
    "query",
    line_number,
  )
  let target_id = words[4]
  let target_start = align_exonerate_parse_int(
    words[5],
    "target start",
    line_number,
  )
  let target_end = align_exonerate_parse_int(
    words[6],
    "target end",
    line_number,
  )
  let target_strand = align_exonerate_parse_strand(
    words[7],
    "target",
    line_number,
  )
  let score = align_exonerate_parse_score(words[8], line_number)
  let operations : Array[AlignExonerateOperation] = []
  let mut index = 9
  while index < words.length() {
    let raw_code = words[index]
    let query_step = align_exonerate_parse_int(
      words[index + 1],
      "query step",
      line_number,
    )
    let target_step = align_exonerate_parse_int(
      words[index + 2],
      "target step",
      line_number,
    )
    if raw_code == "G" {
      if query_step == 0 && target_step > 0 {
        operations.push(AlignExonerateOperation::create("D", 0, target_step))
      } else if target_step == 0 && query_step > 0 {
        operations.push(AlignExonerateOperation::create("I", query_step, 0))
      } else {
        raise AlignExonerateError(
          "Exonerate vulgar gap must advance exactly one sequence at line " +
          line_number.to_string(),
        )
      }
    } else if raw_code == "I" {
      operations.push(
        AlignExonerateOperation::create("N", query_step, target_step),
      )
    } else if raw_code == "N" {
      if target_step > 0 {
        operations.push(AlignExonerateOperation::create("U", 0, target_step))
      }
      if query_step > 0 {
        operations.push(AlignExonerateOperation::create("U", query_step, 0))
      }
      if query_step == 0 && target_step == 0 {
        raise AlignExonerateError(
          "Empty Exonerate non-equivalenced operation at line " +
          line_number.to_string(),
        )
      }
    } else {
      operations.push(
        AlignExonerateOperation::create(raw_code, query_step, target_step),
      )
    }
    index = index + 3
  }
  AlignExonerateAlignment::create(
    query_id, query_start, query_end, query_strand, target_id, target_start, target_end,
    target_strand, score, operations,
  )
}

///|
fn align_exonerate_parse_cigar_line(
  content : String,
  line_number : Int,
) -> AlignExonerateAlignment raise AlignExonerateError {
  let words = align_exonerate_split_whitespace(content)
  if words.length() < 11 || (words.length() - 9) % 2 != 0 {
    raise AlignExonerateError(
      "Malformed Exonerate cigar field count at line " + line_number.to_string(),
    )
  }
  let query_id = words[0]
  let query_start = align_exonerate_parse_int(
    words[1],
    "query start",
    line_number,
  )
  let query_end = align_exonerate_parse_int(words[2], "query end", line_number)
  let query_strand = align_exonerate_parse_strand(
    words[3],
    "query",
    line_number,
  )
  let target_id = words[4]
  let target_start = align_exonerate_parse_int(
    words[5],
    "target start",
    line_number,
  )
  let target_end = align_exonerate_parse_int(
    words[6],
    "target end",
    line_number,
  )
  let target_strand = align_exonerate_parse_strand(
    words[7],
    "target",
    line_number,
  )
  let score = align_exonerate_parse_score(words[8], line_number)
  let operations : Array[AlignExonerateOperation] = []
  let mut raw_query = 0
  let mut raw_target = 0
  let mut normalized_query = 0
  let mut normalized_target = 0
  let mut index = 9
  while index < words.length() {
    let code = words[index]
    let step = align_exonerate_parse_positive_int(
      words[index + 1],
      "cigar step",
      line_number,
    )
    if code == "M" {
      raw_query = raw_query + step
      raw_target = raw_target + step
    } else if code == "I" {
      if query_strand == AlignExonerateProtein &&
        target_strand != AlignExonerateProtein {
        raw_query = raw_query + step * 3
      } else {
        raw_query = raw_query + step
      }
    } else if code == "D" {
      if target_strand == AlignExonerateProtein &&
        query_strand != AlignExonerateProtein {
        raw_target = raw_target + step * 3
      } else {
        raw_target = raw_target + step
      }
    } else {
      raise AlignExonerateError(
        "Unknown Exonerate cigar operation " +
        code +
        " at line " +
        line_number.to_string(),
      )
    }
    let next_query = align_exonerate_normalize_cigar_offset(
      raw_query, query_strand, target_strand,
    )
    let next_target = align_exonerate_normalize_cigar_offset(
      raw_target, target_strand, query_strand,
    )
    let query_step = next_query - normalized_query
    let target_step = next_target - normalized_target
    operations.push(
      AlignExonerateOperation::create(code, query_step, target_step),
    )
    normalized_query = next_query
    normalized_target = next_target
    index = index + 2
  }
  AlignExonerateAlignment::create(
    query_id, query_start, query_end, query_strand, target_id, target_start, target_end,
    target_strand, score, operations,
  )
}

///|
fn align_exonerate_format_header(
  alignment : AlignExonerateAlignment,
) -> Array[String] {
  [
    alignment.query_id,
    alignment.query_start.to_string(),
    alignment.query_end.to_string(),
    align_exonerate_strand_symbol(alignment.query_strand),
    alignment.target_id,
    alignment.target_start.to_string(),
    alignment.target_end.to_string(),
    align_exonerate_strand_symbol(alignment.target_strand),
    alignment.score.to_string(),
  ]
}

///|
fn align_exonerate_format_vulgar(
  alignment : AlignExonerateAlignment,
) -> String raise AlignExonerateError {
  let words = ["vulgar:"]
  for word in align_exonerate_format_header(alignment) {
    words.push(word)
  }
  let mut index = 0
  while index < alignment.operations.length() {
    let operation = alignment.operations[index]
    if operation.code == "U" {
      let mut query_step = 0
      let mut target_step = 0
      while index < alignment.operations.length() &&
            alignment.operations[index].code == "U" {
        query_step = query_step + alignment.operations[index].query_step
        target_step = target_step + alignment.operations[index].target_step
        index = index + 1
      }
      align_exonerate_append_vulgar(words, "N", query_step, target_step)
      continue
    }
    let code = if operation.code == "N" {
      "I"
    } else if operation.code == "D" || operation.code == "I" {
      "G"
    } else {
      operation.code
    }
    align_exonerate_append_vulgar(
      words,
      code,
      operation.query_step,
      operation.target_step,
    )
    index = index + 1
  }
  words.join(" ") + "\n"
}

///|
fn align_exonerate_append_vulgar(
  words : Array[String],
  code : String,
  query_step : Int,
  target_step : Int,
) -> Unit {
  words.push(code)
  words.push(query_step.to_string())
  words.push(target_step.to_string())
}

///|
fn align_exonerate_format_cigar(
  alignment : AlignExonerateAlignment,
) -> String raise AlignExonerateError {
  let words = ["cigar:"]
  for word in align_exonerate_format_header(alignment) {
    words.push(word)
  }
  for operation in alignment.operations {
    if operation.code == "M" {
      let step = if alignment.query_strand == AlignExonerateProtein &&
        alignment.target_strand != AlignExonerateProtein {
        operation.target_step
      } else if alignment.target_strand == AlignExonerateProtein &&
        alignment.query_strand != AlignExonerateProtein {
        operation.query_step
      } else {
        if operation.query_step != operation.target_step {
          raise AlignExonerateError(
            "Cannot encode unequal non-translated match steps as cigar",
          )
        }
        operation.query_step
      }
      align_exonerate_append_cigar(words, "M", step)
    } else if operation.code == "5" || operation.code == "3" {
      align_exonerate_append_cigar_movement(words, operation)
    } else if operation.code == "N" {
      align_exonerate_append_cigar_movement(words, operation)
    } else if operation.code == "C" {
      if operation.query_step != operation.target_step {
        raise AlignExonerateError("Cannot encode unequal codon steps as cigar")
      }
      align_exonerate_append_cigar(words, "M", operation.query_step)
    } else if operation.code == "D" {
      align_exonerate_append_cigar(words, "D", operation.target_step)
    } else if operation.code == "I" {
      align_exonerate_append_cigar(words, "I", operation.query_step)
    } else if operation.code == "U" ||
      operation.code == "S" ||
      operation.code == "F" {
      if operation.target_step > 0 {
        align_exonerate_append_cigar(words, "D", operation.target_step)
      }
      if operation.query_step > 0 {
        align_exonerate_append_cigar(words, "I", operation.query_step)
      }
    } else {
      raise AlignExonerateError(
        "Cannot encode Exonerate operation " + operation.code + " as cigar",
      )
    }
  }
  words.join(" ") + "\n"
}

///|
fn align_exonerate_append_cigar_movement(
  words : Array[String],
  operation : AlignExonerateOperation,
) -> Unit raise AlignExonerateError {
  if operation.query_step == 0 {
    align_exonerate_append_cigar(words, "D", operation.target_step)
  } else if operation.target_step == 0 {
    align_exonerate_append_cigar(words, "I", operation.query_step)
  } else if operation.query_step == operation.target_step {
    align_exonerate_append_cigar(words, "M", operation.query_step)
  } else {
    raise AlignExonerateError(
      "Cannot encode two-axis splice/intron movement as cigar",
    )
  }
}

///|
fn align_exonerate_append_cigar(
  words : Array[String],
  code : String,
  step : Int,
) -> Unit raise AlignExonerateError {
  if step <= 0 {
    raise AlignExonerateError("Exonerate cigar step must be positive")
  }
  words.push(code)
  words.push(step.to_string())
}

///|
fn align_exonerate_map_position(
  alignment : AlignExonerateAlignment,
  position : Int,
  from_query : Bool,
) -> Int? {
  let query_sign = align_exonerate_strand_sign(alignment.query_strand)
  let target_sign = align_exonerate_strand_sign(alignment.target_strand)
  let mut query_boundary = alignment.query_start
  let mut target_boundary = alignment.target_start
  for operation in alignment.operations {
    let source_boundary = if from_query {
      query_boundary
    } else {
      target_boundary
    }
    let source_sign = if from_query { query_sign } else { target_sign }
    let source_step = if from_query {
      operation.query_step
    } else {
      operation.target_step
    }
    let destination_boundary = if from_query {
      target_boundary
    } else {
      query_boundary
    }
    let destination_sign = if from_query { target_sign } else { query_sign }
    let destination_step = if from_query {
      operation.target_step
    } else {
      operation.query_step
    }
    let offset = if source_sign > 0 {
      position - source_boundary
    } else {
      source_boundary - 1 - position
    }
    if offset >= 0 && offset < source_step {
      if destination_step == 0 {
        return None
      }
      let destination_offset = offset * destination_step / source_step
      if destination_sign > 0 {
        return Some(destination_boundary + destination_offset)
      }
      return Some(destination_boundary - 1 - destination_offset)
    }
    query_boundary = query_boundary + query_sign * operation.query_step
    target_boundary = target_boundary + target_sign * operation.target_step
  }
  None
}

///|
fn align_exonerate_validate_operation(
  code : String,
  query_step : Int,
  target_step : Int,
) -> Unit raise AlignExonerateError {
  if query_step < 0 || target_step < 0 {
    raise AlignExonerateError("Exonerate operation steps must be non-negative")
  }
  if query_step == 0 && target_step == 0 {
    raise AlignExonerateError("Exonerate operation must advance a sequence")
  }
  if code == "M" {
    if query_step == 0 || target_step == 0 {
      raise AlignExonerateError("Exonerate match must advance both sequences")
    }
  } else if code == "5" || code == "3" {
    if query_step != 2 && target_step != 2 {
      raise AlignExonerateError(
        "Exonerate splice-site operation must contain a two-unit step",
      )
    }
  } else if code == "N" {
    if (query_step == 0) == (target_step == 0) {
      raise AlignExonerateError(
        "Exonerate intron must advance exactly one sequence",
      )
    }
  } else if code == "C" {
    if query_step == 0 ||
      target_step == 0 ||
      query_step % 3 != 0 ||
      target_step % 3 != 0 {
      raise AlignExonerateError(
        "Exonerate codon operation must use positive multiples of three",
      )
    }
  } else if code == "D" {
    if query_step != 0 || target_step <= 0 {
      raise AlignExonerateError(
        "Exonerate deletion must advance only the target",
      )
    }
  } else if code == "I" {
    if target_step != 0 || query_step <= 0 {
      raise AlignExonerateError(
        "Exonerate insertion must advance only the query",
      )
    }
  } else if code == "U" {
    if (query_step == 0) == (target_step == 0) {
      raise AlignExonerateError(
        "Normalized Exonerate U must advance exactly one sequence",
      )
    }
  } else if code == "S" {
    ()
  } else if code == "F" {
    if (query_step == 0) == (target_step == 0) {
      raise AlignExonerateError(
        "Exonerate frame shift must advance exactly one sequence",
      )
    }
  } else {
    raise AlignExonerateError("Unknown Exonerate operation " + code)
  }
}

///|
fn align_exonerate_validate_axis(
  start : Int,
  end : Int,
  strand : AlignExonerateStrand,
  field : String,
) -> Unit raise AlignExonerateError {
  if start < 0 || end < 0 {
    raise AlignExonerateError(
      "Exonerate " + field + " coordinates must be non-negative",
    )
  }
  if strand == AlignExonerateReverse {
    if end >= start {
      raise AlignExonerateError(
        "Reverse Exonerate " + field + " coordinates must decrease",
      )
    }
  } else if end <= start {
    raise AlignExonerateError(
      "Forward/protein Exonerate " + field + " coordinates must increase",
    )
  }
}

///|
fn align_exonerate_validate_id(
  id : String,
  field : String,
) -> Unit raise AlignExonerateError {
  if id.length() == 0 {
    raise AlignExonerateError(
      "Exonerate " + field + " identifier must not be empty",
    )
  }
  for index = 0; index < id.length(); index = index + 1 {
    let code = id.unsafe_get(index).to_int()
    if align_exonerate_is_whitespace(code) || code == 127 {
      raise AlignExonerateError(
        "Exonerate " + field + " identifier must not contain whitespace",
      )
    }
  }
}

///|
fn align_exonerate_validate_header_value(
  value : String,
  field : String,
) -> Unit raise AlignExonerateError {
  for index = 0; index < value.length(); index = index + 1 {
    let code = value.unsafe_get(index).to_int()
    if code == '\n'.to_int() || code == '\r'.to_int() {
      raise AlignExonerateError(
        "Exonerate " + field + " must not contain a newline",
      )
    }
  }
}

///|
fn align_exonerate_parse_header_line(
  line : String,
  prefix : String,
  field : String,
  line_number : Int,
) -> String raise AlignExonerateError {
  if !line.has_prefix(prefix) {
    raise AlignExonerateError(
      "Missing Exonerate " + field + " at line " + line_number.to_string(),
    )
  }
  let value = line[prefix.length():].trim().to_owned()
  if value.length() < 2 ||
    value.unsafe_get(0).to_int() != '['.to_int() ||
    value.unsafe_get(value.length() - 1).to_int() != ']'.to_int() {
    raise AlignExonerateError(
      "Malformed bracketed Exonerate " +
      field +
      " at line " +
      line_number.to_string(),
    )
  }
  value[1:value.length() - 1].to_owned()
}

///|
fn align_exonerate_parse_strand(
  text : String,
  field : String,
  line_number : Int,
) -> AlignExonerateStrand raise AlignExonerateError {
  if text == "+" {
    AlignExonerateForward
  } else if text == "-" {
    AlignExonerateReverse
  } else if text == "." {
    AlignExonerateProtein
  } else {
    raise AlignExonerateError(
      "Invalid Exonerate " +
      field +
      " strand at line " +
      line_number.to_string(),
    )
  }
}

///|
fn align_exonerate_strand_symbol(strand : AlignExonerateStrand) -> String {
  if strand == AlignExonerateForward {
    "+"
  } else if strand == AlignExonerateReverse {
    "-"
  } else {
    "."
  }
}

///|
fn align_exonerate_strand_sign(strand : AlignExonerateStrand) -> Int {
  if strand == AlignExonerateReverse {
    -1
  } else {
    1
  }
}

///|
fn align_exonerate_parse_positive_int(
  text : String,
  field : String,
  line_number : Int,
) -> Int raise AlignExonerateError {
  let value = align_exonerate_parse_int(text, field, line_number)
  if value <= 0 {
    raise AlignExonerateError(
      "Exonerate " +
      field +
      " must be positive at line " +
      line_number.to_string(),
    )
  }
  value
}

///|
fn align_exonerate_parse_int(
  text : String,
  field : String,
  line_number : Int,
) -> Int raise AlignExonerateError {
  if text.length() == 0 {
    raise AlignExonerateError(
      "Empty Exonerate " + field + " at line " + line_number.to_string(),
    )
  }
  let mut value = 0
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code < '0'.to_int() || code > '9'.to_int() {
      raise AlignExonerateError(
        "Invalid Exonerate " + field + " at line " + line_number.to_string(),
      )
    }
    let digit = code - '0'.to_int()
    if value > 214748364 || (value == 214748364 && digit > 7) {
      raise AlignExonerateError(
        "Exonerate " +
        field +
        " is too large at line " +
        line_number.to_string(),
      )
    }
    value = value * 10 + digit
  }
  value
}

///|
fn align_exonerate_parse_score(
  text : String,
  line_number : Int,
) -> Double raise AlignExonerateError {
  let mut has_digit = false
  for index = 0; index < text.length(); index = index + 1 {
    let code = text.unsafe_get(index).to_int()
    if code >= '0'.to_int() && code <= '9'.to_int() {
      has_digit = true
    } else if code != '+'.to_int() &&
      code != '-'.to_int() &&
      code != '.'.to_int() &&
      code != 'e'.to_int() &&
      code != 'E'.to_int() {
      raise AlignExonerateError(
        "Invalid Exonerate score at line " + line_number.to_string(),
      )
    }
  }
  if !has_digit {
    raise AlignExonerateError(
      "Invalid Exonerate score at line " + line_number.to_string(),
    )
  }
  match parse_double(text) {
    Some(value) =>
      if value.abs() > 1.0e300 {
        raise AlignExonerateError(
          "Non-finite Exonerate score at line " + line_number.to_string(),
        )
      } else {
        value
      }
    None =>
      raise AlignExonerateError(
        "Invalid Exonerate score at line " + line_number.to_string(),
      )
  }
}

///|
fn align_exonerate_normalize_cigar_offset(
  raw : Int,
  strand : AlignExonerateStrand,
  other : AlignExonerateStrand,
) -> Int {
  if strand == AlignExonerateProtein && other != AlignExonerateProtein {
    (raw + 2) / 3
  } else {
    raw
  }
}

///|
fn align_exonerate_lines(text : String) -> Array[String] {
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
fn align_exonerate_split_whitespace(text : String) -> Array[String] {
  let result : Array[String] = []
  let mut index = 0
  while index < text.length() {
    while index < text.length() &&
          align_exonerate_is_whitespace(text.unsafe_get(index).to_int()) {
      index = index + 1
    }
    if index >= text.length() {
      break
    }
    let start = index
    while index < text.length() &&
          !align_exonerate_is_whitespace(text.unsafe_get(index).to_int()) {
      index = index + 1
    }
    result.push(text[start:index].to_owned())
  }
  result
}

///|
fn align_exonerate_is_whitespace(code : Int) -> Bool {
  code == ' '.to_int() ||
  code == '\t'.to_int() ||
  code == '\n'.to_int() ||
  code == '\r'.to_int()
}

///|
fn align_exonerate_abs(value : Int) -> Int {
  if value < 0 {
    -value
  } else {
    value
  }
}

///|
/// Small report covering splicing and translated reverse-strand coordinates.
pub fn align_exonerate_example_text() -> String {
  "Command line: [exonerate -m est2genome query.fa target.fa --showvulgar yes]\n" +
  "Hostname: [moonbit]\n" +
  "vulgar: transcript 0 18 + chromosome 100 218 + 250 M 6 6 5 0 2 I 0 96 3 0 2 M 12 12\n" +
  "vulgar: protein 0 4 . genome 500 488 - 80 M 4 12\n" +
  "-- completed exonerate analysis\n"
}
