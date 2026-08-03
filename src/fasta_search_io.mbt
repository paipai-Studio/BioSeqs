///|
/// Bio.SearchIO.FastaIO - FASTA search tool tabular output parsing.
/// Mirrors Biopython's Bio.SearchIO.FastaIO parser for the `fasta36`
/// `-m 8` (compact tabular) and `-m 9` (tabular with comment header)
/// output formats.
///
/// Column layout (tab-separated, 12 columns):
///   1. query id
///   2. subject id
///   3. percent identity (0-100)
///   4. alignment length
///   5. mismatches
///   6. gap openings
///   7. query start
///   8. query end
///   9. subject start
///  10. subject end
///  11. E-value
///  12. bit score
///
/// For `-m 9` output, the file also contains comment lines beginning with
/// `#` that summarise the program invocation and per-query statistics.
/// These comment lines are parsed and exposed via FastaMeta records.

// ============================================================================
// Internal helpers (mirror the ones in hmmer_io.mbt but kept local to avoid
// cross-file name clashes)
// ============================================================================

///|
fn fasta_sio_starts_with(s : String, prefix : String) -> Bool {
  let n = s.length()
  let p = prefix.length()
  if p > n {
    return false
  }
  let mut i = 0
  while i < p {
    if s.unsafe_get(i).to_int() != prefix.unsafe_get(i).to_int() {
      return false
    }
    i = i + 1
  }
  true
}

///|
fn fasta_sio_parse_double(s : String) -> Double {
  match parse_double(s) {
    Some(v) => v
    None => 0.0
  }
}

///|
fn fasta_sio_parse_int(s : String) -> Int {
  parse_int(s)
}

///|
/// Trim leading and trailing ASCII whitespace from a string.
fn fasta_sio_trim(s : String) -> String {
  let n = s.length()
  let mut start = 0
  let mut end = n
  while start < end {
    let c = s.unsafe_get(start).to_int()
    if c == ' '.to_int() || c == '\t'.to_int() || c == '\r'.to_int() {
      start = start + 1
    } else {
      break
    }
  }
  while end > start {
    let c = s.unsafe_get(end - 1).to_int()
    if c == ' '.to_int() || c == '\t'.to_int() || c == '\r'.to_int() {
      end = end - 1
    } else {
      break
    }
  }
  s[start:end].to_string()
}

///|
/// Split a line by runs of ASCII whitespace, dropping empty tokens.
fn fasta_sio_split_whitespace(s : String) -> Array[String] {
  let n = s.length()
  let tokens : Array[String] = []
  let mut start = 0
  let mut i = 0
  let mut in_token = false
  while i < n {
    let c = s.unsafe_get(i).to_int()
    let is_ws = c == ' '.to_int() || c == '\t'.to_int() || c == '\r'.to_int()
    if is_ws {
      if in_token {
        tokens.push(s[start:i].to_string())
        in_token = false
      }
    } else if !in_token {
      start = i
      in_token = true
    }
    i = i + 1
  }
  if in_token {
    tokens.push(s[start:n].to_string())
  }
  tokens
}

// ============================================================================
// Data structures
// ============================================================================

///|
/// Metadata extracted from FASTA `-m 9` comment (`#`) header lines.
pub struct FastaMeta {
  /// Program name (e.g. "FASTA", "SSEARCH", "GGSEARCH").
  mut program : String
  /// Program version string (e.g. "36.3.8h").
  mut version : String
  /// Database name (may be empty).
  mut database : String
  /// Raw comment lines, in order.
  comments : Array[String]
}

///|
pub fn FastaMeta::new() -> FastaMeta {
  FastaMeta::{ program: "", version: "", database: "", comments: [] }
}

///|
/// Result of parsing FASTA `-m 8`/`-m 9` output: both the SearchIO
/// QueryResult objects and the optional metadata block.
pub struct FastaSearchResult {
  /// Parsed query results (one per query id).
  queries : Array[QueryResult]
  /// Metadata from `#` header lines (empty for `-m 8`).
  meta : FastaMeta
}

///|
pub fn FastaSearchResult::new(
  queries : Array[QueryResult],
  meta : FastaMeta,
) -> FastaSearchResult {
  FastaSearchResult::{ queries, meta }
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Parse FASTA `-m 8` (compact tabular) output.
///
/// Each non-empty, non-comment line is split on tabs (falling back to
/// runs of whitespace if a line contains no tabs) into the standard
/// 12 columns. Rows sharing the same query id are aggregated into a
/// single QueryResult with one Hit per (subject id) and one HSP per row.
pub fn parse_fasta_m8(content : String) -> Array[QueryResult] {
  let lines = split_by_char(content, '\n'.to_int())
  let results : Array[QueryResult] = []
  let query_order : Array[String] = []
  let by_query : Map[String, Array[FastaRow]] = Map([], capacity=0)

  let mut i = 0
  while i < lines.length() {
    let ln = lines[i]
    if ln.length() == 0 || fasta_sio_starts_with(ln, "#") {
      i = i + 1
      continue
    }
    let row = fasta_sio_parse_row(ln)
    if row.query_id.length() == 0 {
      i = i + 1
      continue
    }
    let arr = by_query.get(row.query_id)
    match arr {
      Some(a) => a.push(row)
      None => {
        query_order.push(row.query_id)
        let na : Array[FastaRow] = [row]
        by_query.set(row.query_id, na)
        ()
      }
    }
    i = i + 1
  }

  let mut qi = 0
  while qi < query_order.length() {
    let qid = query_order[qi]
    let rows = by_query.get(qid)
    let rows_arr = match rows {
      Some(a) => a
      None => {
        qi = qi + 1
        continue
      }
    }
    let hits : Array[Hit] = []
    let target_order : Array[String] = []
    let by_target : Map[String, Array[FastaRow]] = Map([], capacity=0)
    let mut ri = 0
    while ri < rows_arr.length() {
      let r = rows_arr[ri]
      let arr = by_target.get(r.subject_id)
      match arr {
        Some(a) => a.push(r)
        None => {
          target_order.push(r.subject_id)
          let na : Array[FastaRow] = [r]
          by_target.set(r.subject_id, na)
          ()
        }
      }
      ri = ri + 1
    }
    let mut ti = 0
    while ti < target_order.length() {
      let sid = target_order[ti]
      let trows = by_target.get(sid)
      let trows_arr = match trows {
        Some(a) => a
        None => {
          ti = ti + 1
          continue
        }
      }
      let hsps : Array[HSP] = []
      let mut tbs = 0.0
      let mut mev = 1.0
      let mut hi = 0
      while hi < trows_arr.length() {
        let r = trows_arr[hi]
        let qstrand = if r.query_end >= r.query_start {
          Strand::Plus
        } else {
          Strand::Minus
        }
        let sstrand = if r.subject_end >= r.subject_start {
          Strand::Plus
        } else {
          Strand::Minus
        }
        let aln_len = if r.alignment_length > 0 {
          r.alignment_length
        } else {
          0
        }
        let frag = HSPFragment::new(
          r.query_start,
          r.query_end,
          r.subject_start,
          r.subject_end,
          qstrand,
          sstrand,
          aln_len,
          r.gap_openings,
          r.bit_score,
          r.evalue,
          "",
          "",
          "",
        )
        let frags : Array[HSPFragment] = [frag]
        hsps.push(HSP::new(sid, "", frags, r.bit_score, r.evalue, r.bit_score))
        tbs = tbs + r.bit_score
        if r.evalue < mev {
          mev = r.evalue
        }
        ignore(r.percent_identity)
        ignore(r.mismatches)
        hi = hi + 1
      }
      hits.push(Hit::new(sid, "", "", hsps, tbs, mev))
      ti = ti + 1
    }
    results.push(
      QueryResult::new(qid, "", 0, hits, "fasta", "36", ""),
    )
    qi = qi + 1
  }

  results
}

///|
/// Parse FASTA `-m 9` (tabular with `#` comment header) output.
///
/// In addition to the tabular rows (parsed exactly like `-m 8`), the
/// `#`-prefixed comment lines are collected into a FastaMeta record.
/// Common comment patterns recognised:
///   * `# FASTA <version>`        - program and version
///   * `# Query: <id>`            - query id (informational)
///   * `# Database: <name>`       - database name
pub fn parse_fasta_m9(content : String) -> FastaSearchResult {
  let lines = split_by_char(content, '\n'.to_int())
  let meta = FastaMeta::new()
  let data_lines : Array[String] = []

  let mut i = 0
  while i < lines.length() {
    let ln = lines[i]
    if ln.length() == 0 {
      i = i + 1
      continue
    }
    if fasta_sio_starts_with(ln, "#") {
      let comment = fasta_sio_trim(ln)
      meta.comments.push(comment)
      // "# FASTA 36.3.8h May 10, 2023"
      if fasta_sio_starts_with(comment, "# FASTA ") {
        let rest = comment[2:comment.length()].to_string()
        let parts = fasta_sio_split_whitespace(rest)
        if parts.length() > 0 {
          meta.program = parts[0]
        }
        if parts.length() > 1 {
          meta.version = parts[1]
        }
      } else if fasta_sio_starts_with(comment, "# Database:") {
        let rest = fasta_sio_substr_after(comment, "# Database:")
        let trimmed = fasta_sio_trim(rest)
        if trimmed.length() > 0 {
          meta.database = trimmed
        }
      } else if fasta_sio_starts_with(comment, "# Query:") {
        let rest = fasta_sio_substr_after(comment, "# Query:")
        let trimmed = fasta_sio_trim(rest)
        if trimmed.length() > 0 && meta.program.length() == 0 {
          meta.program = "FASTA"
        }
      }
      i = i + 1
      continue
    }
    data_lines.push(ln)
    i = i + 1
  }

  // Parse the data rows as `-m 8` style tabular output.
  let mut data_text = ""
  let mut di = 0
  while di < data_lines.length() {
    if di > 0 {
      data_text = data_text + "\n"
    }
    data_text = data_text + data_lines[di]
    di = di + 1
  }
  let queries = parse_fasta_m8(data_text)
  FastaSearchResult::new(queries, meta)
}

///|
/// Return the substring after the first occurrence of `marker`, or "" if not found.
fn fasta_sio_substr_after(s : String, marker : String) -> String {
  let n = s.length()
  let m = marker.length()
  if m == 0 || m > n {
    return ""
  }
  let mut i = 0
  while i <= n - m {
    let mut j = 0
    let mut matched = true
    while j < m {
      if s.unsafe_get(i + j).to_int() != marker.unsafe_get(j).to_int() {
        matched = false
        break
      }
      j = j + 1
    }
    if matched {
      return s[i + m:n].to_string()
    }
    i = i + 1
  }
  ""
}

///|
/// Parsed FASTA tabular row.
struct FastaRow {
  query_id : String
  subject_id : String
  percent_identity : Double
  alignment_length : Int
  mismatches : Int
  gap_openings : Int
  query_start : Int
  query_end : Int
  subject_start : Int
  subject_end : Int
  evalue : Double
  bit_score : Double
}

///|
/// Parse a single FASTA `-m 8`/`-m 9` tabular line into a FastaRow.
/// Lines are split on tabs; if a line has no tabs, falls back to
/// whitespace splitting.
fn fasta_sio_parse_row(ln : String) -> FastaRow {
  // Prefer tab-splitting (the standard FASTA -m 8 delimiter).
  let has_tab = {
    let mut found = false
    let mut k = 0
    while k < ln.length() {
      if ln.unsafe_get(k).to_int() == '\t'.to_int() {
        found = true
        break
      }
      k = k + 1
    }
    found
  }
  let parts = if has_tab {
    split_by_char(ln, '\t'.to_int())
  } else {
    fasta_sio_split_whitespace(ln)
  }
  if parts.length() < 12 {
    // Pad short rows with empty strings so field access is safe.
    let padded : Array[String] = []
    let mut p = 0
    while p < 12 {
      if p < parts.length() {
        padded.push(parts[p])
      } else {
        padded.push("")
      }
      p = p + 1
    }
    return FastaRow::{
      query_id: padded[0],
      subject_id: padded[1],
      percent_identity: fasta_sio_parse_double(padded[2]),
      alignment_length: fasta_sio_parse_int(padded[3]),
      mismatches: fasta_sio_parse_int(padded[4]),
      gap_openings: fasta_sio_parse_int(padded[5]),
      query_start: fasta_sio_parse_int(padded[6]),
      query_end: fasta_sio_parse_int(padded[7]),
      subject_start: fasta_sio_parse_int(padded[8]),
      subject_end: fasta_sio_parse_int(padded[9]),
      evalue: fasta_sio_parse_double(padded[10]),
      bit_score: fasta_sio_parse_double(padded[11]),
    }
  }
  FastaRow::{
    query_id: parts[0],
    subject_id: parts[1],
    percent_identity: fasta_sio_parse_double(parts[2]),
    alignment_length: fasta_sio_parse_int(parts[3]),
    mismatches: fasta_sio_parse_int(parts[4]),
    gap_openings: fasta_sio_parse_int(parts[5]),
    query_start: fasta_sio_parse_int(parts[6]),
    query_end: fasta_sio_parse_int(parts[7]),
    subject_start: fasta_sio_parse_int(parts[8]),
    subject_end: fasta_sio_parse_int(parts[9]),
    evalue: fasta_sio_parse_double(parts[10]),
    bit_score: fasta_sio_parse_double(parts[11]),
  }
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate a small FASTA `-m 8` (compact tabular) sample for tests and demos.
pub fn fasta_m8_sample() -> String {
  "query1\tsubject1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.3\n" +
    "query1\tsubject2\t80.0\t80\t16\t2\t1\t80\t10\t90\t1e-30\t150.1\n" +
    "query2\tsubject1\t70.0\t60\t18\t3\t1\t60\t5\t65\t1e-20\t100.5\n"
}

///|
/// Generate a small FASTA `-m 9` (tabular with comment header) sample
/// for tests and demos.
pub fn fasta_m9_sample() -> String {
  "# FASTA 36.3.8h May 10, 2023\n" +
    "# Query: query1 - example query 1\n" +
    "# Database: uniprot_sprot.fasta\n" +
    "# Fields: query id, subject id, % identity, alignment length, mismatches, gap openings, q.start, q.end, s.start, s.end, e-value, bit score\n" +
    "query1\tsubject1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.3\n" +
    "query1\tsubject2\t80.0\t80\t16\t2\t1\t80\t10\t90\t1e-30\t150.1\n" +
    "# Query: query2 - example query 2\n" +
    "query2\tsubject1\t70.0\t60\t18\t3\t1\t60\t5\t65\t1e-20\t100.5\n"
}
