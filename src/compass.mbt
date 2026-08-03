///|
/// COMPASS - Comparison of Multiple Protein Alignments with Assessment of
/// Statistical Significance.
///
/// References:
///   - Biopython Bio.Compass module
///   - Sadreyev RI, Baker D, Grishin NV (2003) Profile-profile comparisons by
///     COMPASS predict intricate homologies between protein families.
///     Protein Sci 12(10):2262-72.
///
/// COMPASS is a tool for profile-profile alignment. This module parses
/// COMPASS output files, which contain alignment records with statistical
/// significance scores (E-values, Smith-Waterman scores, percentage identity).

// ============================================================================
// Data structures
// ============================================================================

///|
/// A single COMPASS alignment record.
pub struct CompassRecord {
  query_name : String
  template_name : String
  query_n_seqs : Int
  template_n_seqs : Int
  query_n_cols : Int
  template_n_cols : Int
  sw_score : Double
  e_value : Double
  percentage_identity : Double
  query_start : Int
  query_end : Int
  template_start : Int
  template_end : Int
  aligned_query : String
  aligned_template : String
  consensus_line : String
}

///|
/// Create a CompassRecord.
pub fn CompassRecord::new(
  query_name : String,
  template_name : String,
  query_n_seqs : Int,
  template_n_seqs : Int,
  query_n_cols : Int,
  template_n_cols : Int,
  sw_score : Double,
  e_value : Double,
  percentage_identity : Double,
  query_start : Int,
  query_end : Int,
  template_start : Int,
  template_end : Int,
  aligned_query : String,
  aligned_template : String,
  consensus_line : String,
) -> CompassRecord {
  {
    query_name,
    template_name,
    query_n_seqs,
    template_n_seqs,
    query_n_cols,
    template_n_cols,
    sw_score,
    e_value,
    percentage_identity,
    query_start,
    query_end,
    template_start,
    template_end,
    aligned_query,
    aligned_template,
    consensus_line,
  }
}

///|
/// Get the query alignment name.
pub fn CompassRecord::query_name(self : CompassRecord) -> String {
  self.query_name
}

///|
/// Get the template alignment name.
pub fn CompassRecord::template_name(self : CompassRecord) -> String {
  self.template_name
}

///|
/// Get the number of sequences in the query alignment.
pub fn CompassRecord::query_n_seqs(self : CompassRecord) -> Int {
  self.query_n_seqs
}

///|
/// Get the number of sequences in the template alignment.
pub fn CompassRecord::template_n_seqs(self : CompassRecord) -> Int {
  self.template_n_seqs
}

///|
/// Get the number of columns in the query alignment.
pub fn CompassRecord::query_n_cols(self : CompassRecord) -> Int {
  self.query_n_cols
}

///|
/// Get the number of columns in the template alignment.
pub fn CompassRecord::template_n_cols(self : CompassRecord) -> Int {
  self.template_n_cols
}

///|
/// Get the Smith-Waterman score.
pub fn CompassRecord::sw_score(self : CompassRecord) -> Double {
  self.sw_score
}

///|
/// Get the E-value.
pub fn CompassRecord::e_value(self : CompassRecord) -> Double {
  self.e_value
}

///|
/// Get the percentage identity.
pub fn CompassRecord::percentage_identity(self : CompassRecord) -> Double {
  self.percentage_identity
}

///|
/// Get the query start position.
pub fn CompassRecord::query_start(self : CompassRecord) -> Int {
  self.query_start
}

///|
/// Get the query end position.
pub fn CompassRecord::query_end(self : CompassRecord) -> Int {
  self.query_end
}

///|
/// Get the template start position.
pub fn CompassRecord::template_start(self : CompassRecord) -> Int {
  self.template_start
}

///|
/// Get the template end position.
pub fn CompassRecord::template_end(self : CompassRecord) -> Int {
  self.template_end
}

///|
/// Get the aligned query sequence.
pub fn CompassRecord::aligned_query(self : CompassRecord) -> String {
  self.aligned_query
}

///|
/// Get the aligned template sequence.
pub fn CompassRecord::aligned_template(self : CompassRecord) -> String {
  self.aligned_template
}

///|
/// Get the consensus line.
pub fn CompassRecord::consensus_line(self : CompassRecord) -> String {
  self.consensus_line
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Extract the text after a prefix and trim surrounding whitespace.
fn compass_extract_after(line : String, prefix : String) -> String {
  let rest = line[prefix.length():line.length()].to_owned()
  trim(rest)
}

///|
/// Parse an alignment line of the form "Query: 1 MALKSLVRLFG 120".
/// Returns (start, sequence, end). The first token is the label
/// ("Query:" or "Template:"), the second token is the start position,
/// the last token is the end position, and all middle tokens form the
/// aligned sequence.
fn compass_parse_alignment_line(line : String) -> (Int, String, Int) {
  let tokens = split_by_whitespace(line)
  if tokens.length() >= 4 {
    let start = parse_int(tokens[1])
    let end_idx = tokens.length() - 1
    let end = parse_int(tokens[end_idx])
    let mut seq = ""
    for i in 2..<end_idx {
      seq = seq + tokens[i]
    }
    (start, seq, end)
  } else {
    (0, "", 0)
  }
}

///|
/// Parse a complete COMPASS output file into an array of records.
///
/// A new record begins at each "COMPASS version" line. The parser
/// extracts all field values and the alignment block (query sequence,
/// consensus line, and template sequence).
pub fn parse_compass(content : String) -> Array[CompassRecord] {
  let lines = split_by_char(content, '\n'.to_int())
  let records : Array[CompassRecord] = []

  let mut query_name = ""
  let mut template_name = ""
  let mut query_n_seqs = 0
  let mut template_n_seqs = 0
  let mut query_n_cols = 0
  let mut template_n_cols = 0
  let mut sw_score = 0.0
  let mut e_value = 0.0
  let mut percentage_identity = 0.0
  let mut query_start = 0
  let mut query_end = 0
  let mut template_start = 0
  let mut template_end = 0
  let mut aligned_query = ""
  let mut aligned_template = ""
  let mut consensus_line = ""

  let mut has_record = false
  // State: 0 = normal field parsing
  //        1 = just saw Query: alignment line, expecting consensus/template
  //        2 = just saw consensus, expecting Template: alignment line
  let mut state = 0

  for line in lines {
    let trimmed = trim(line.to_string())

    // State 2: expecting Template: alignment line after consensus
    if state == 2 {
      if trimmed.length() > 0 && starts_with(trimmed, 0, "Template:") {
        let (ts, tseq, te) = compass_parse_alignment_line(trimmed)
        template_start = ts
        template_end = te
        aligned_template = tseq
        state = 0
        continue
      }
      // Not a Template: line, reset to normal parsing
      state = 0
    }

    // State 1: expecting consensus or Template after Query: alignment line
    if state == 1 {
      if trimmed.length() > 0 && starts_with(trimmed, 0, "Template:") {
        let (ts, tseq, te) = compass_parse_alignment_line(trimmed)
        template_start = ts
        template_end = te
        aligned_template = tseq
        state = 0
        continue
      }
      // This is the consensus line
      consensus_line = trimmed
      state = 2
      continue
    }

    // State 0: normal field parsing
    if trimmed.length() == 0 {
      continue
    }

    if starts_with(trimmed, 0, "COMPASS version") {
      // Save previous record if it exists
      if has_record {
        records.push(
          CompassRecord::new(
            query_name,
            template_name,
            query_n_seqs,
            template_n_seqs,
            query_n_cols,
            template_n_cols,
            sw_score,
            e_value,
            percentage_identity,
            query_start,
            query_end,
            template_start,
            template_end,
            aligned_query,
            aligned_template,
            consensus_line,
          ),
        )
        // Reset all fields for the new record
        query_name = ""
        template_name = ""
        query_n_seqs = 0
        template_n_seqs = 0
        query_n_cols = 0
        template_n_cols = 0
        sw_score = 0.0
        e_value = 0.0
        percentage_identity = 0.0
        query_start = 0
        query_end = 0
        template_start = 0
        template_end = 0
        aligned_query = ""
        aligned_template = ""
        consensus_line = ""
      }
      has_record = true
    } else if starts_with(trimmed, 0, "Query alignment:") {
      query_name = compass_extract_after(trimmed, "Query alignment:")
    } else if starts_with(trimmed, 0, "Template alignment:") {
      template_name = compass_extract_after(trimmed, "Template alignment:")
    } else if starts_with(
      trimmed, 0, "Number of sequences in query alignment:",
    ) {
      query_n_seqs = parse_int(
        compass_extract_after(trimmed, "Number of sequences in query alignment:"),
      )
    } else if starts_with(
      trimmed, 0, "Number of sequences in template alignment:",
    ) {
      template_n_seqs = parse_int(
        compass_extract_after(
          trimmed, "Number of sequences in template alignment:",
        ),
      )
    } else if starts_with(
      trimmed, 0, "Number of columns in query alignment:",
    ) {
      query_n_cols = parse_int(
        compass_extract_after(trimmed, "Number of columns in query alignment:"),
      )
    } else if starts_with(
      trimmed, 0, "Number of columns in template alignment:",
    ) {
      template_n_cols = parse_int(
        compass_extract_after(
          trimmed, "Number of columns in template alignment:",
        ),
      )
    } else if starts_with(trimmed, 0, "Smith-Waterman score:") {
      sw_score = parse_double(
        compass_extract_after(trimmed, "Smith-Waterman score:"),
      ).unwrap_or(0.0)
    } else if starts_with(trimmed, 0, "E-value:") {
      e_value = parse_double(
        compass_extract_after(trimmed, "E-value:"),
      ).unwrap_or(0.0)
    } else if starts_with(trimmed, 0, "Percentage identity:") {
      percentage_identity = parse_double(
        compass_extract_after(trimmed, "Percentage identity:"),
      ).unwrap_or(0.0)
    } else if starts_with(trimmed, 0, "Query:") {
      // Query alignment line: "Query: 1 MALKSLVRLFG 120"
      let (qs, qseq, qe) = compass_parse_alignment_line(trimmed)
      query_start = qs
      query_end = qe
      aligned_query = qseq
      state = 1
    } else if starts_with(trimmed, 0, "Template:") {
      // Standalone Template alignment line (without preceding Query:)
      let (ts, tseq, te) = compass_parse_alignment_line(trimmed)
      template_start = ts
      template_end = te
      aligned_template = tseq
    }
  }

  // Save the last record
  if has_record {
    records.push(
      CompassRecord::new(
        query_name,
        template_name,
        query_n_seqs,
        template_n_seqs,
        query_n_cols,
        template_n_cols,
        sw_score,
        e_value,
        percentage_identity,
        query_start,
        query_end,
        template_start,
        template_end,
        aligned_query,
        aligned_template,
        consensus_line,
      ),
    )
  }

  records
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Extract the COMPASS version string from the output.
/// Returns the version (e.g., "2.4.2") or an empty string if not found.
pub fn compass_version(content : String) -> String {
  let lines = split_by_char(content, '\n'.to_int())
  for line in lines {
    let trimmed = trim(line.to_string())
    if starts_with(trimmed, 0, "COMPASS version") {
      let tokens = split_by_whitespace(trimmed)
      if tokens.length() >= 3 {
        return tokens[2]
      }
    }
  }
  ""
}

///|
/// Filter records by E-value threshold (records with E-value <= threshold).
pub fn compass_filter_by_evalue(
  records : Array[CompassRecord],
  threshold : Double,
) -> Array[CompassRecord] {
  let result : Array[CompassRecord] = []
  for r in records {
    if r.e_value <= threshold {
      result.push(r)
    }
  }
  result
}

///|
/// Filter records by percentage identity threshold
/// (records with identity >= threshold).
pub fn compass_filter_by_identity(
  records : Array[CompassRecord],
  threshold : Double,
) -> Array[CompassRecord] {
  let result : Array[CompassRecord] = []
  for r in records {
    if r.percentage_identity >= threshold {
      result.push(r)
    }
  }
  result
}

///|
/// Get the alignment length of a record (length of the aligned query sequence).
pub fn compass_alignment_length(record : CompassRecord) -> Int {
  record.aligned_query.length()
}

///|
/// Generate a human-readable summary of all COMPASS records.
pub fn compass_summary(records : Array[CompassRecord]) -> String {
  let mut s = "COMPASS Summary\n"
  s = s + "===============\n"
  s = s + "Total records: " + records.length().to_string() + "\n\n"
  for r in records {
    s = s + "Query: " + r.query_name + "\n"
    s = s + "  Template: " + r.template_name + "\n"
    s = s + "  SW score: " + r.sw_score.to_string() + "\n"
    s = s + "  E-value: " + r.e_value.to_string() + "\n"
    s = s + "  Identity: " + r.percentage_identity.to_string() + "%\n"
    s = s + "  Alignment length: " +
      compass_alignment_length(r).to_string() +
      "\n\n"
  }
  s
}

///|
/// Format a CompassRecord as a readable string.
pub fn CompassRecord::to_string(self : CompassRecord) -> String {
  let mut s = "CompassRecord(query=" + self.query_name
  s = s + ", template=" + self.template_name
  s = s + ", sw_score=" + self.sw_score.to_string()
  s = s + ", e_value=" + self.e_value.to_string()
  s = s + ", identity=" + self.percentage_identity.to_string()
  s = s + ", query_range=" +
    self.query_start.to_string() +
    "-" +
    self.query_end.to_string()
  s = s + ", template_range=" +
    self.template_start.to_string() +
    "-" +
    self.template_end.to_string()
  s = s + ")"
  s
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate sample COMPASS output text with two records for testing.
pub fn compass_sample_data() -> String {
  let mut s = ""
  // Record 1: strong match with low E-value
  s = s + "COMPASS version 2.4.2\n"
  s = s + "Query alignment: query1.msa\n"
  s = s + "Template alignment: template1.msa\n"
  s = s + " Number of sequences in query alignment: 5\n"
  s = s + " Number of sequences in template alignment: 4\n"
  s = s + " Number of columns in query alignment: 120\n"
  s = s + " Number of columns in template alignment: 110\n"
  s = s + " Smith-Waterman score: 185.5\n"
  s = s + " E-value: 3.5e-10\n"
  s = s + " Percentage identity: 28.5\n"
  s = s + " Query:      1 MALKSLVRLFG 120\n"
  s = s + "                  :  :: :::\n"
  s = s + " Template:   1 MGVKSAVKT 110\n"
  s = s + "\n"
  // Record 2: weaker match with higher E-value
  s = s + "COMPASS version 2.4.2\n"
  s = s + "Query alignment: query2.msa\n"
  s = s + "Template alignment: template2.msa\n"
  s = s + " Number of sequences in query alignment: 3\n"
  s = s + " Number of sequences in template alignment: 6\n"
  s = s + " Number of columns in query alignment: 80\n"
  s = s + " Number of columns in template alignment: 90\n"
  s = s + " Smith-Waterman score: 95.2\n"
  s = s + " E-value: 1.2e-05\n"
  s = s + " Percentage identity: 15.3\n"
  s = s + " Query:      1 ACDEFGHIKL 80\n"
  s = s + "                  : :  : ::: \n"
  s = s + " Template:   1 ACDEFGHIKL 90\n"
  s
}
