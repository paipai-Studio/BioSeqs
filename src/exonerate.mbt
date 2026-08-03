///|
/// Bio.SearchIO.ExonerateIO module, inspired by Biopython's
/// Bio.SearchIO.ExonerateIO. Parses Exonerate output in vulgar and cigar
/// formats. Exonerate is a tool for sequence alignment, particularly for
/// spliced alignments.
///
/// Vulgar format line:
///   `vulgar: qname qstart qend qstrand tname tstart tend tstrand score vblocks`
/// where vblocks is a series of 4-field tuples: `code qlen tlen score`.
///
/// Cigar format line:
///   `cigar: qname qstart qend qstrand tname tstart tend tstrand score cblocks`
/// where cblocks is a series of 2-field pairs: `code length`.

// ============================================================================
// Data Structures
// ============================================================================

///|
/// Represents a single alignment block (vulgar/cigar operation) between
/// query and target sequences.
pub struct AlignmentBlock {
  /// Block code: "M", "I", "5", "3", "S", "G", "U", "V", "C", "N", etc.
  code : String
  /// Length consumed on the query sequence.
  query_length : Int
  /// Length consumed on the target sequence.
  target_length : Int
  /// Block score (for introns, this holds the intron length).
  block_score : Double
}

///|
/// Construct an AlignmentBlock.
pub fn AlignmentBlock::new(
  code : String,
  query_length : Int,
  target_length : Int,
  block_score : Double,
) -> AlignmentBlock {
  AlignmentBlock::{ code, query_length, target_length, block_score }
}

///|
/// Get the block code.
pub fn AlignmentBlock::code(self : AlignmentBlock) -> String {
  self.code
}

///|
/// Get the query length consumed by this block.
pub fn AlignmentBlock::query_length(self : AlignmentBlock) -> Int {
  self.query_length
}

///|
/// Get the target length consumed by this block.
pub fn AlignmentBlock::target_length(self : AlignmentBlock) -> Int {
  self.target_length
}

///|
/// Get the block score.
pub fn AlignmentBlock::block_score(self : AlignmentBlock) -> Double {
  self.block_score
}

///|
/// Represents a single exonerate alignment record (one vulgar/cigar line).
pub struct ExonerateRecord {
  query_name : String
  query_start : Int
  query_end : Int
  query_strand : String
  target_name : String
  target_start : Int
  target_end : Int
  target_strand : String
  score : Double
  alignment_blocks : Array[AlignmentBlock]
}

///|
/// Construct an ExonerateRecord.
pub fn ExonerateRecord::new(
  query_name : String,
  query_start : Int,
  query_end : Int,
  query_strand : String,
  target_name : String,
  target_start : Int,
  target_end : Int,
  target_strand : String,
  score : Double,
  alignment_blocks : Array[AlignmentBlock],
) -> ExonerateRecord {
  ExonerateRecord::{
    query_name,
    query_start,
    query_end,
    query_strand,
    target_name,
    target_start,
    target_end,
    target_strand,
    score,
    alignment_blocks,
  }
}

///|
/// Get the query sequence name.
pub fn ExonerateRecord::query_name(self : ExonerateRecord) -> String {
  self.query_name
}

///|
/// Get the query start coordinate.
pub fn ExonerateRecord::query_start(self : ExonerateRecord) -> Int {
  self.query_start
}

///|
/// Get the query end coordinate.
pub fn ExonerateRecord::query_end(self : ExonerateRecord) -> Int {
  self.query_end
}

///|
/// Get the query strand ("+" or "-").
pub fn ExonerateRecord::query_strand(self : ExonerateRecord) -> String {
  self.query_strand
}

///|
/// Get the target sequence name.
pub fn ExonerateRecord::target_name(self : ExonerateRecord) -> String {
  self.target_name
}

///|
/// Get the target start coordinate.
pub fn ExonerateRecord::target_start(self : ExonerateRecord) -> Int {
  self.target_start
}

///|
/// Get the target end coordinate.
pub fn ExonerateRecord::target_end(self : ExonerateRecord) -> Int {
  self.target_end
}

///|
/// Get the target strand ("+" or "-").
pub fn ExonerateRecord::target_strand(self : ExonerateRecord) -> String {
  self.target_strand
}

///|
/// Get the overall alignment score.
pub fn ExonerateRecord::score(self : ExonerateRecord) -> Double {
  self.score
}

///|
/// Get the list of alignment blocks.
pub fn ExonerateRecord::alignment_blocks(
  self : ExonerateRecord
) -> Array[AlignmentBlock] {
  self.alignment_blocks
}

// ============================================================================
// Internal parsing helpers
// ============================================================================

///|
/// Split a string by whitespace (space and tab), returning non-empty tokens.
fn exo_split_whitespace(s : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let len = s.length()
  let mut i = 0
  while i < len {
    while i < len && (s.unsafe_get(i) == ' ' || s.unsafe_get(i) == '\t') {
      i = i + 1
    }
    if i >= len {
      break
    }
    let start = i
    while i < len && s.unsafe_get(i) != ' ' && s.unsafe_get(i) != '\t' {
      i = i + 1
    }
    result.push(s[start:i].to_owned())
  }
  result
}

///|
/// Split content into lines by newline (handles \r\n and \n).
fn exo_split_lines(content : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let len = content.length()
  let mut start = 0
  let mut i = 0
  while i < len {
    let c = content.unsafe_get(i)
    if c == '\n' {
      let line_end = if i > 0 && content.unsafe_get(i - 1) == '\r' {
        i - 1
      } else {
        i
      }
      result.push(content[start:line_end].to_owned())
      start = i + 1
    }
    i = i + 1
  }
  if start < len {
    let line_end = if len > 0 && content.unsafe_get(len - 1) == '\r' {
      len - 1
    } else {
      len
    }
    result.push(content[start:line_end].to_owned())
  } else if start == len && len > 0 {
    let line_end = if content.unsafe_get(len - 1) == '\r' {
      len - 1
    } else {
      len
    }
    if line_end > 0 {
      result.push(content[start:line_end].to_owned())
    }
  }
  result
}

///|
/// Parse a 32-bit signed integer from a String (no error checking).
fn exo_parse_int(s : String) -> Int {
  let n = s.length()
  if n == 0 {
    return 0
  }
  let mut result = 0
  let zero = '0'.to_int()
  let mut i = 0
  let mut is_negative = false
  if s.unsafe_get(0).to_int() == 45 {
    is_negative = true
    i = 1
  }
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c >= zero && c <= zero + 9 {
      result = result * 10 + (c - zero)
    }
    i = i + 1
  }
  if is_negative {
    -result
  } else {
    result
  }
}

///|
/// Parse a Double from a String (supports scientific notation, no error).
fn exo_parse_double(s : String) -> Double {
  let n = s.length()
  if n == 0 {
    return 0.0
  }
  let mut result = 0.0
  let mut decimal = 0.0
  let mut decimal_places = 0
  let mut is_negative = false
  let mut has_exponent = false
  let mut exponent_val = 0
  let mut exp_negative = false
  let zero = '0'.to_int()
  let mut i = 0
  if s.unsafe_get(0).to_int() == 45 {
    is_negative = true
    i = 1
  }
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c == 46 {
      i = i + 1
      while i < n {
        let dc = s.unsafe_get(i).to_int()
        if dc >= zero && dc <= zero + 9 {
          decimal = decimal * 10.0 + (dc - zero).to_double()
          decimal_places = decimal_places + 1
        }
        i = i + 1
      }
      break
    }
    if c == 101 || c == 69 {
      has_exponent = true
      i = i + 1
      if i < n && s.unsafe_get(i).to_int() == 45 {
        exp_negative = true
        i = i + 1
      } else if i < n && s.unsafe_get(i).to_int() == 43 {
        i = i + 1
      }
      while i < n {
        let ec = s.unsafe_get(i).to_int()
        if ec >= zero && ec <= zero + 9 {
          exponent_val = exponent_val * 10 + (ec - zero)
        }
        i = i + 1
      }
      break
    }
    if c >= zero && c <= zero + 9 {
      result = result * 10.0 + (c - zero).to_double()
    }
    i = i + 1
  }
  while decimal_places > 0 {
    decimal = decimal / 10.0
    decimal_places = decimal_places - 1
  }
  let mut final_result = result + decimal
  if has_exponent {
    if exp_negative {
      let mut j = 0
      while j < exponent_val {
        final_result = final_result / 10.0
        j = j + 1
      }
    } else {
      let mut j = 0
      while j < exponent_val {
        final_result = final_result * 10.0
        j = j + 1
      }
    }
  }
  if is_negative {
    -final_result
  } else {
    final_result
  }
}

///|
/// Check whether a string token is a numeric value (integer or float,
/// possibly with a leading sign). Used to distinguish codes from numbers
/// when parsing variable-arity vulgar blocks.
fn exo_is_numeric_token(s : String) -> Bool {
  let n = s.length()
  if n == 0 {
    return false
  }
  let mut i = 0
  if s.unsafe_get(0).to_int() == 45 || s.unsafe_get(0).to_int() == 43 {
    i = 1
    if n == 1 {
      return false
    }
  }
  while i < n {
    let c = s.unsafe_get(i).to_int()
    let zero = '0'.to_int()
    if !(c >= zero && c <= zero + 9) && c != 46 && c != 101 && c != 69 {
      return false
    }
    i = i + 1
  }
  true
}

// ============================================================================
// Vulgar parsing
// ============================================================================

///|
/// Parse the vulgar block tokens into AlignmentBlock entries.
/// Each block is normally 4 fields: `code qlen tlen score`. To be tolerant
/// of inputs where splice-site codes (5/3) omit the score, we treat a
/// following code as a boundary and default the score to 0.
fn exo_parse_vulgar_blocks(tokens : Array[String]) -> Array[AlignmentBlock] {
  let blocks : Array[AlignmentBlock] = Array::new()
  let n = tokens.length()
  let mut i = 0
  while i < n {
    let code = tokens[i]
    if i + 2 < n {
      let qlen = exo_parse_int(tokens[i + 1])
      let tlen = exo_parse_int(tokens[i + 2])
      let (score, advance) = if i + 3 < n &&
        exo_is_numeric_token(tokens[i + 3]) {
        (exo_parse_double(tokens[i + 3]), 4)
      } else {
        (0.0, 3)
      }
      i = i + advance
      blocks.push(AlignmentBlock::new(code, qlen, tlen, score))
    } else {
      // Not enough tokens for a full block; stop.
      i = i + 1
    }
  }
  blocks
}

///|
/// Parse a single vulgar line into an ExonerateRecord.
fn exo_parse_vulgar_line(line : String) -> ExonerateRecord? {
  let tokens = exo_split_whitespace(line)
  // Need at least: "vulgar:" qname qstart qend qstrand tname tstart tend
  // tstrand score = 10 tokens, plus at least one block token.
  if tokens.length() < 10 {
    return None
  }
  if tokens[0] != "vulgar:" {
    return None
  }
  let query_name = tokens[1]
  let query_start = exo_parse_int(tokens[2])
  let query_end = exo_parse_int(tokens[3])
  let query_strand = tokens[4]
  let target_name = tokens[5]
  let target_start = exo_parse_int(tokens[6])
  let target_end = exo_parse_int(tokens[7])
  let target_strand = tokens[8]
  let score = exo_parse_double(tokens[9])
  let block_tokens : Array[String] = Array::new()
  for j in 10..<tokens.length() {
    block_tokens.push(tokens[j])
  }
  let blocks = exo_parse_vulgar_blocks(block_tokens)
  Some(
    ExonerateRecord::new(
      query_name,
      query_start,
      query_end,
      query_strand,
      target_name,
      target_start,
      target_end,
      target_strand,
      score,
      blocks,
    ),
  )
}

///|
/// Parse Exonerate vulgar format output. Returns one ExonerateRecord per
/// "vulgar:" line. Blank lines and non-vulgar lines are ignored.
pub fn parse_exonerate_vulgar(content : String) -> Array[ExonerateRecord] {
  let records : Array[ExonerateRecord] = Array::new()
  let lines = exo_split_lines(content)
  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed.has_prefix("vulgar:") {
      match exo_parse_vulgar_line(trimmed.to_owned()) {
        Some(rec) => records.push(rec)
        None => ()
      }
    }
  }
  records
}

// ============================================================================
// Cigar parsing
// ============================================================================

///|
/// Determine query and target lengths for a CIGAR operation given the
/// operation code and the single length value. For most operations both
/// sequences consume the length; introns (I) consume neither; gaps (U/V/D)
/// consume only one side.
fn exo_cigar_lengths(
  code : String,
  length : Int,
) -> (Int, Int, Double) {
  match code {
    "I" => (0, 0, length.to_double())
    "D" => (length, 0, 0.0)
    "U" => (0, length, 0.0)
    "V" => (length, 0, 0.0)
    "5" | "3" => (0, 0, 0.0)
    "G" | "N" => (0, 0, 0.0)
    _ => (length, length, 0.0)
  }
}

///|
/// Parse the cigar block tokens (code/length pairs) into AlignmentBlock
/// entries.
fn exo_parse_cigar_blocks(tokens : Array[String]) -> Array[AlignmentBlock] {
  let blocks : Array[AlignmentBlock] = Array::new()
  let n = tokens.length()
  let mut i = 0
  while i + 1 < n {
    let code = tokens[i]
    let length = exo_parse_int(tokens[i + 1])
    let (qlen, tlen, score) = exo_cigar_lengths(code, length)
    blocks.push(AlignmentBlock::new(code, qlen, tlen, score))
    i = i + 2
  }
  blocks
}

///|
/// Parse a single cigar line into an ExonerateRecord.
fn exo_parse_cigar_line(line : String) -> ExonerateRecord? {
  let tokens = exo_split_whitespace(line)
  if tokens.length() < 10 {
    return None
  }
  if tokens[0] != "cigar:" {
    return None
  }
  let query_name = tokens[1]
  let query_start = exo_parse_int(tokens[2])
  let query_end = exo_parse_int(tokens[3])
  let query_strand = tokens[4]
  let target_name = tokens[5]
  let target_start = exo_parse_int(tokens[6])
  let target_end = exo_parse_int(tokens[7])
  let target_strand = tokens[8]
  let score = exo_parse_double(tokens[9])
  let block_tokens : Array[String] = Array::new()
  for j in 10..<tokens.length() {
    block_tokens.push(tokens[j])
  }
  let blocks = exo_parse_cigar_blocks(block_tokens)
  Some(
    ExonerateRecord::new(
      query_name,
      query_start,
      query_end,
      query_strand,
      target_name,
      target_start,
      target_end,
      target_strand,
      score,
      blocks,
    ),
  )
}

///|
/// Parse Exonerate cigar format output. Returns one ExonerateRecord per
/// "cigar:" line. Blank lines and non-cigar lines are ignored.
pub fn parse_exonerate_cigar(content : String) -> Array[ExonerateRecord] {
  let records : Array[ExonerateRecord] = Array::new()
  let lines = exo_split_lines(content)
  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed.has_prefix("cigar:") {
      match exo_parse_cigar_line(trimmed.to_owned()) {
        Some(rec) => records.push(rec)
        None => ()
      }
    }
  }
  records
}

// ============================================================================
// Auto-detecting parser
// ============================================================================

///|
/// Parse Exonerate output, auto-detecting vulgar vs cigar format per line.
/// Mixed inputs (some vulgar lines, some cigar lines) are supported; each
/// line is dispatched by its prefix.
pub fn parse_exonerate(content : String) -> Array[ExonerateRecord] {
  let records : Array[ExonerateRecord] = Array::new()
  let lines = exo_split_lines(content)
  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    if trimmed.has_prefix("vulgar:") {
      match exo_parse_vulgar_line(trimmed.to_owned()) {
        Some(rec) => records.push(rec)
        None => ()
      }
    } else if trimmed.has_prefix("cigar:") {
      match exo_parse_cigar_line(trimmed.to_owned()) {
        Some(rec) => records.push(rec)
        None => ()
      }
    }
  }
  records
}

// ============================================================================
// Analysis functions
// ============================================================================

///|
/// Filter records by minimum score. Returns records whose score is greater
/// than or equal to the threshold.
pub fn exonerate_filter_by_score(
  records : Array[ExonerateRecord],
  threshold : Double,
) -> Array[ExonerateRecord] {
  let result : Array[ExonerateRecord] = Array::new()
  for rec in records {
    if rec.score >= threshold {
      result.push(rec)
    }
  }
  result
}

///|
/// Total aligned (match) length for a record: sum of query lengths of all
/// "M" blocks. Matches Biopython's notion of alignment span from vulgar M
/// blocks.
pub fn exonerate_total_alignment_length(record : ExonerateRecord) -> Int {
  let mut total = 0
  for block in record.alignment_blocks {
    if block.code == "M" {
      total = total + block.query_length
    }
  }
  total
}

///|
/// Total intron length for a record: sum of the intron lengths (block_score
/// of "I" blocks) across all intron blocks.
pub fn exonerate_total_intron_length(record : ExonerateRecord) -> Int {
  let mut total = 0
  for block in record.alignment_blocks {
    if block.code == "I" {
      total = total + block.block_score.to_int()
    }
  }
  total
}

///|
/// Count the number of intron ("I") blocks in a record.
pub fn exonerate_count_introns(record : ExonerateRecord) -> Int {
  let mut count = 0
  for block in record.alignment_blocks {
    if block.code == "I" {
      count = count + 1
    }
  }
  count
}

///|
/// Count the number of match ("M") blocks in a record.
pub fn exonerate_count_matches(record : ExonerateRecord) -> Int {
  let mut count = 0
  for block in record.alignment_blocks {
    if block.code == "M" {
      count = count + 1
    }
  }
  count
}

// ============================================================================
// Reconstruction functions
// ============================================================================

///|
/// Reconstruct the vulgar block string from a record's alignment blocks.
/// Each block is rendered as `code qlen tlen score` (score as an integer).
pub fn exonerate_to_vulgar_string(record : ExonerateRecord) -> String {
  let buf = StringBuilder::new()
  let n = record.alignment_blocks.length()
  for i in 0..<n {
    let block = record.alignment_blocks[i]
    if i > 0 {
      buf.write_char(' ')
    }
    buf.write_string(block.code)
    buf.write_char(' ')
    buf.write_string(block.query_length.to_string())
    buf.write_char(' ')
    buf.write_string(block.target_length.to_string())
    buf.write_char(' ')
    buf.write_string(block.block_score.to_int().to_string())
  }
  buf.to_string()
}

///|
/// Reconstruct the cigar block string from a record's alignment blocks.
/// Each block is rendered as `code length`, where length is the intron
/// length for "I" blocks, the query length for "V"/"D" blocks, the target
/// length for "U" blocks, and the query length otherwise.
pub fn exonerate_to_cigar_string(record : ExonerateRecord) -> String {
  let buf = StringBuilder::new()
  let n = record.alignment_blocks.length()
  for i in 0..<n {
    let block = record.alignment_blocks[i]
    if i > 0 {
      buf.write_char(' ')
    }
    let length = match block.code {
      "I" => block.block_score.to_int()
      "U" => block.target_length
      _ => block.query_length
    }
    buf.write_string(block.code)
    buf.write_char(' ')
    buf.write_string(length.to_string())
  }
  buf.to_string()
}

// ============================================================================
// Summary
// ============================================================================

///|
/// Generate a human-readable summary of a list of exonerate records.
pub fn exonerate_summary(records : Array[ExonerateRecord]) -> String {
  let buf = StringBuilder::new()
  buf.write_string("Exonerate records: ")
  buf.write_string(records.length().to_string())
  buf.write_string("\n")
  for i in 0..<records.length() {
    let rec = records[i]
    let aln_len = exonerate_total_alignment_length(rec)
    let intron_len = exonerate_total_intron_length(rec)
    let n_introns = exonerate_count_introns(rec)
    let n_matches = exonerate_count_matches(rec)
    buf.write_string("  [")
    buf.write_string((i + 1).to_string())
    buf.write_string("] ")
    buf.write_string(rec.query_name)
    buf.write_string(" (")
    buf.write_string(rec.query_start.to_string())
    buf.write_string("-")
    buf.write_string(rec.query_end.to_string())
    buf.write_string(rec.query_strand)
    buf.write_string(") -> ")
    buf.write_string(rec.target_name)
    buf.write_string(" (")
    buf.write_string(rec.target_start.to_string())
    buf.write_string("-")
    buf.write_string(rec.target_end.to_string())
    buf.write_string(rec.target_strand)
    buf.write_string(") score=")
    buf.write_string(rec.score.to_int().to_string())
    buf.write_string(" blocks=")
    buf.write_string(rec.alignment_blocks.length().to_string())
    buf.write_string(" M=")
    buf.write_string(n_matches.to_string())
    buf.write_string(" aln_len=")
    buf.write_string(aln_len.to_string())
    buf.write_string(" introns=")
    buf.write_string(n_introns.to_string())
    buf.write_string(" intron_len=")
    buf.write_string(intron_len.to_string())
    buf.write_string("\n")
  }
  buf.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Sample exonerate vulgar output for testing and demos.
pub fn exonerate_sample_data() -> String {
  "vulgar: query1 0 100 + target1 0 95 + 1850 M 80 80 80 I 0 0 2000 M 20 15 20\n" +
  "vulgar: query2 10 210 + target2 30 240 + 3200 M 100 100 100 5 0 0 0 I 0 0 1500 3 0 0 0 M 100 100 100\n" +
  "vulgar: query3 0 50 - target3 500 550 + 850 M 50 50 50\n"
}

///|
/// Sample exonerate cigar output for testing and demos.
pub fn exonerate_sample_cigar_data() -> String {
  "cigar: query1 0 100 + target1 0 95 + 1850 M 80 I 2000 M 20\n" +
  "cigar: query2 10 210 + target2 30 240 + 3200 M 100 I 1500 M 100\n" +
  "cigar: query3 0 50 - target3 500 550 + 850 M 50\n"
}
