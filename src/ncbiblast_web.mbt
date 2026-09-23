///|
/// Bio.Blast.NCBIWWW - NCBI BLAST web service client.
///
/// Port of Biopython's Bio.Blast.NCBIWWW module, which submits sequence
/// similarity searches to the public NCBI BLAST URL API
/// (https://blast.ncbi.nlm.nih.gov/Blast.cgi): a PUT call registers the
/// search and returns an RID (request id) plus an estimated completion
/// time (RTOE); GET calls with CMD=Get then poll the search status
/// (WAITING / READY / FAILED) and finally download the result text.
///
/// As with the other REST clients in this library, the transport is
/// injected (`NcbiBlastClient::new(fetch~...)`), so no real networking or
/// sleeping happens inside the library; the same Biopython-style
/// three-requests-per-second rate limiting is tracked via an injectable
/// clock.

///|
/// Error raised for invalid qblast arguments, protocol failures, or
/// transport errors.
pub suberror NcbiBlastError {
  NcbiBlastError(String)
}

///|
fn[T] ncbiblast_fail(msg : String) -> T raise NcbiBlastError {
  raise NcbiBlastError(msg)
}

///|
/// Base URL of the NCBI BLAST service.
pub let ncbi_blast_url : String = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

///|
/// BLAST search programs accepted by the NCBI URL API.
pub let ncbi_blast_programs : Array[String] = [
  "blastn", "blastp", "blastx", "tblastn", "tblastx",
]

///|
/// Status of a submitted NCBI BLAST search.
pub enum NcbiBlastStatus {
  Waiting
  Ready
  Failed
  Unknown
} derive(Debug)
pub extend NcbiBlastStatus with @moonbitlang/core/debug.Debug::{to_repr}

// ============================================================================
// URL construction
// ============================================================================

///|
/// Percent-encode a value for application/x-www-form-urlencoded query
/// strings: alphanumerics and a small set of punctuation pass through,
/// spaces become '+' and every other byte becomes %XX.
pub fn qblast_quote(value : String) -> String {
  let result = StringBuilder::new()
  let mut i = 0
  while i < value.length() {
    let code = value.unsafe_get(i).to_int()
    if (code >= 'a'.to_int() && code <= 'z'.to_int()) ||
      (code >= 'A'.to_int() && code <= 'Z'.to_int()) ||
      (code >= '0'.to_int() && code <= '9'.to_int()) ||
      code == '-'.to_int() ||
      code == '_'.to_int() ||
      code == '.'.to_int() ||
      code == '~'.to_int() ||
      code == '*'.to_int() {
      result.write_string(value[i:i + 1].to_owned())
    } else if code == ' '.to_int() {
      result.write_string("+")
    } else {
      let hi = (code >> 4) & 0x0F
      let lo = code & 0x0F
      result.write_string(
        "%" + ncbi_blast_hex_digit(hi) + ncbi_blast_hex_digit(lo),
      )
    }
    i = i + 1
  }
  result.to_string()
}

///|
fn ncbi_blast_hex_digit(d : Int) -> String {
  match d {
    0 => "0"
    1 => "1"
    2 => "2"
    3 => "3"
    4 => "4"
    5 => "5"
    6 => "6"
    7 => "7"
    8 => "8"
    9 => "9"
    10 => "A"
    11 => "B"
    12 => "C"
    13 => "D"
    14 => "E"
    _ => "F"
  }
}

///|
/// Encode ordered key/value pairs as an urlencoded query string.
fn qblast_urlencode(params : Array[(String, String)]) -> String {
  let parts : Array[String] = Array::new()
  let mut i = 0
  while i < params.length() {
    let (key, value) = params[i]
    parts.push(qblast_quote(key) + "=" + qblast_quote(value))
    i = i + 1
  }
  let result = StringBuilder::new()
  i = 0
  while i < parts.length() {
    if i > 0 {
      result.write_string("&")
    }
    result.write_string(parts[i])
    i = i + 1
  }
  result.to_string()
}

///|
/// Validate the BLAST program name (mirrors Biopython's check).
pub fn qblast_validate_program(program : String) -> Unit raise NcbiBlastError {
  let mut ok = false
  let mut i = 0
  while i < ncbi_blast_programs.length() {
    if ncbi_blast_programs[i] == program {
      ok = true
    }
    i = i + 1
  }
  if !ok {
    ncbiblast_fail(
      "Program specified is " +
      program +
      " but should be one of " +
      ncbi_blast_programs.to_string(),
    )
  }
}

///|
/// Build the PUT URL that submits a BLAST search.
///
/// Optional parameters mirror Biopython qblast() arguments:
/// expect, hitlist_size, entrez_query, format_type (XML/Text/HTML),
/// descriptions, alignments, matrix_name, gap_costs, word_size,
/// threshold, nucl_reward, nucl_penalty, megablast, ncbi_gi.
pub fn qblast_put_url(
  program : String,
  database : String,
  query : String,
  expect? : Double = 10.0,
  hitlist_size? : Int = 500,
  entrez_query? : String = "(none)",
  format_type? : String = "XML",
  descriptions? : Int = 500,
  alignments? : Int = 500,
  matrix_name? : String = "",
  gap_costs? : String = "",
  word_size? : Int = 0,
  threshold? : Int = 0,
  nucl_reward? : Int = 0,
  nucl_penalty? : Int = 0,
  megablast? : Bool = false,
  ncbi_gi? : Bool = false,
) -> String raise NcbiBlastError {
  qblast_validate_program(program)
  if database.length() == 0 {
    ncbiblast_fail("No database specified for BLAST search")
  }
  if query.length() == 0 {
    ncbiblast_fail("No sequence supplied for BLAST search")
  }

  let params : Array[(String, String)] = Array::new()
  params.push(("CMD", "Put"))
  if megablast {
    // Biopython keeps PROGRAM as blastn for megablast, which just adds SERVICE.
    params.push(("PROGRAM", "blastn"))
    params.push(("SERVICE", "megablast"))
  } else {
    params.push(("PROGRAM", program))
  }
  params.push(("DATABASE", database))
  params.push(("QUERY", query))
  params.push(("DESCRIPTIONS", descriptions.to_string()))
  params.push(("ALIGNMENTS", alignments.to_string()))
  params.push(("HITLIST_SIZE", hitlist_size.to_string()))
  params.push(("EXPECT", expect.to_string()))
  if entrez_query.length() > 0 && entrez_query != "(none)" {
    params.push(("ENTREZ_QUERY", entrez_query))
  }
  if matrix_name.length() > 0 {
    params.push(("MATRIX_NAME", matrix_name))
  }
  if gap_costs.length() > 0 {
    params.push(("GAPCOSTS", gap_costs))
  }
  if word_size > 0 {
    params.push(("WORD_SIZE", word_size.to_string()))
  }
  if threshold > 0 {
    params.push(("THRESHOLD", threshold.to_string()))
  }
  if nucl_reward != 0 {
    params.push(("NUCL_REWARD", nucl_reward.to_string()))
  }
  if nucl_penalty != 0 {
    params.push(("NUCL_PENALTY", nucl_penalty.to_string()))
  }
  if ncbi_gi {
    params.push(("NCBI_GI", "Yes"))
  }
  params.push(("FORMAT_TYPE", format_type))

  ncbi_blast_url + "?" + qblast_urlencode(params)
}

///|
/// Build the GET URL used to poll an RID's status (SearchInfo object).
pub fn qblast_search_info_url(rid : String) -> String {
  let params : Array[(String, String)] = [
    ("CMD", "Get"),
    ("RID", rid),
    ("FORMAT_OBJECT", "SearchInfo"),
  ]
  ncbi_blast_url + "?" + qblast_urlencode(params)
}

///|
/// Build the GET URL used to download the finished result.
pub fn qblast_result_url(rid : String, format_type? : String = "XML") -> String {
  let params : Array[(String, String)] = Array::new()
  params.push(("CMD", "Get"))
  params.push(("RID", rid))
  params.push(("FORMAT_TYPE", format_type))
  ncbi_blast_url + "?" + qblast_urlencode(params)
}

// ============================================================================
// Response parsing
// ============================================================================

///|
/// Read the token following "label =" or "label=" in a response, stopping
/// at whitespace (mirrors Biopython's RID/RTOE/Status regular expressions).
fn ncbi_blast_scan_field(text : String, label : String) -> String? {
  let spaced = label + " = "
  let packed = label + "="
  let mut start = -1
  match text.find(spaced) {
    Some(idx) => start = idx + spaced.length()
    None =>
      match text.find(packed) {
        Some(idx) => start = idx + packed.length()
        None => ()
      }
  }
  if start < 0 {
    return None
  }
  let mut end = start
  while end < text.length() {
    let ch = text.unsafe_get(end).to_int()
    let is_end = ch == ' '.to_int() ||
      ch == '\t'.to_int() ||
      ch == '\n'.to_int() ||
      ch == '\r'.to_int() ||
      ch == ';'.to_int()
    if is_end {
      break
    }
    end = end + 1
  }
  Some(text[start:end].to_owned())
}

///|
/// Parse a non-negative base-10 integer; returns 0 when there are no digits.
fn ncbi_blast_atoi(value : String) -> Int {
  let zero = '0'.to_int()
  let mut result = 0
  let mut found = false
  let mut i = 0
  while i < value.length() {
    let c = value.unsafe_get(i).to_int()
    if c >= zero && c <= zero + 9 {
      result = result * 10 + (c - zero)
      found = true
    } else if found {
      break
    }
    i = i + 1
  }
  result
}

///|
/// Parse a PUT response: returns (RID, RTOE seconds). Raises if absent.
pub fn qblast_parse_put(
  response : String,
) -> (String, Int) raise NcbiBlastError {
  let rid = match ncbi_blast_scan_field(response, "RID") {
    Some(value) => value
    None =>
      ncbiblast_fail("No RID found in BLAST response (submission rejected?)")
  }
  let rtoe = match ncbi_blast_scan_field(response, "RTOE") {
    Some(value) => ncbi_blast_atoi(value)
    None => 0
  }
  (rid, rtoe)
}

///|
/// Parse a SearchInfo response into a status.
pub fn qblast_parse_status(response : String) -> NcbiBlastStatus {
  match ncbi_blast_scan_field(response, "Status") {
    Some(status) => {
      let s = status
      if s == "WAITING" {
        NcbiBlastStatus::Waiting
      } else if s == "READY" {
        NcbiBlastStatus::Ready
      } else if s == "FAILED" {
        NcbiBlastStatus::Failed
      } else {
        NcbiBlastStatus::Unknown
      }
    }
    None => NcbiBlastStatus::Unknown
  }
}

// ============================================================================
// Client
// ============================================================================

///|
/// NCBI BLAST client with an injectable transport.
pub struct NcbiBlastClient {
  priv fetch : (String) -> String raise NcbiBlastError
  priv now : () -> Double
  priv mut previous : Double
}

///|
/// Create a client.
///
/// - fetch: callback performing the HTTP transfer: url -> response text.
/// - now: optional wall clock in seconds; defaults to a virtual clock.
pub fn NcbiBlastClient::new(
  fetch~ : (String) -> String raise NcbiBlastError,
  now? : () -> Double,
) -> NcbiBlastClient {
  { fetch, now: now.unwrap_or(fn() { 0.0 }), previous: 0.0, }
}

///|
/// Enforce NCBI's "three requests per second" guideline and issue a GET.
fn NcbiBlastClient::call(
  self : NcbiBlastClient,
  url : String,
) -> String raise NcbiBlastError {
  let delay = 0.333333333
  let current = (self.now)()
  let wait = self.previous + delay - current
  if wait > 0 {
    self.previous = current + wait
  } else {
    self.previous = current
  }
  (self.fetch)(url)
}

///|
/// Submit a search and return (RID, RTOE).
pub fn NcbiBlastClient::submit(
  self : NcbiBlastClient,
  program : String,
  database : String,
  query : String,
  expect? : Double = 10.0,
  hitlist_size? : Int = 500,
  entrez_query? : String = "(none)",
  format_type? : String = "XML",
  matrix_name? : String = "",
  gap_costs? : String = "",
  word_size? : Int = 0,
  threshold? : Int = 0,
  nucl_reward? : Int = 0,
  nucl_penalty? : Int = 0,
  megablast? : Bool = false,
  ncbi_gi? : Bool = false,
) -> (String, Int) raise NcbiBlastError {
  let url = qblast_put_url(
    program,
    database,
    query,
    expect~,
    hitlist_size~,
    entrez_query~,
    format_type~,
    matrix_name~,
    gap_costs~,
    word_size~,
    threshold~,
    nucl_reward~,
    nucl_penalty~,
    megablast~,
    ncbi_gi~,
  )
  let response = self.call(url)
  qblast_parse_put(response)
}

///|
/// Poll an RID's status.
pub fn NcbiBlastClient::check_status(
  self : NcbiBlastClient,
  rid : String,
) -> NcbiBlastStatus raise NcbiBlastError {
  let url = qblast_search_info_url(rid)
  let response = self.call(url)
  qblast_parse_status(response)
}

///|
/// Download the result text for a finished RID.
pub fn NcbiBlastClient::fetch_result(
  self : NcbiBlastClient,
  rid : String,
  format_type? : String = "XML",
) -> String raise NcbiBlastError {
  let url = qblast_result_url(rid, format_type~)
  self.call(url)
}

///|
/// High-level workflow: submit, poll until READY (or FAILED), and return
/// the result text. The caller controls timing through its injected
/// fetch/clock; at most `max_polls` status checks are performed.
pub fn NcbiBlastClient::qblast(
  self : NcbiBlastClient,
  program : String,
  database : String,
  query : String,
  max_polls? : Int = 120,
  expect? : Double = 10.0,
  hitlist_size? : Int = 500,
  entrez_query? : String = "(none)",
  format_type? : String = "XML",
  matrix_name? : String = "",
  gap_costs? : String = "",
  word_size? : Int = 0,
  megablast? : Bool = false,
  ncbi_gi? : Bool = false,
) -> String raise NcbiBlastError {
  let (rid, rtoe) = self.submit(
    program,
    database,
    query,
    expect~,
    hitlist_size~,
    entrez_query~,
    format_type~,
    matrix_name~,
    gap_costs~,
    word_size~,
    megablast~,
    ncbi_gi~,
  )
  let limit = max_polls
  let mut polls = 0
  let mut outcome = NcbiBlastStatus::Waiting
  while polls < limit {
    let status = self.check_status(rid)
    outcome = status
    if status is NcbiBlastStatus::Ready {
      break
    }
    if status is NcbiBlastStatus::Failed {
      ncbiblast_fail("BLAST search " + rid + " failed at NCBI")
    }
    if status is NcbiBlastStatus::Unknown {
      ncbiblast_fail(
        "BLAST search " +
        rid +
        " returned an unknown status (RTOE was " +
        rtoe.to_string() +
        "s)",
      )
    }
    polls = polls + 1
  }
  if outcome is NcbiBlastStatus::Waiting {
    ncbiblast_fail(
      "BLAST search " +
      rid +
      " did not finish after " +
      limit.to_string() +
      " polls",
    )
  }
  self.fetch_result(rid, format_type~)
}
