///|
/// Bio.TogoWS - access to the TogoWS integrated web services of DBCLS, Japan.
///
/// Port of Biopython's Bio.TogoWS module. The TogoWS REST service provides
/// simple access to a range of databases, acting as a proxy to shield users
/// from all the different provider APIs. See: http://togows.dbcls.jp/
///
/// In the absence of clear usage guidelines, this port enforces the same
/// "up to three queries per second" limit as the original module: every
/// request goes through a rate limiter that records the earliest time the
/// next call is allowed (one third of a second apart). Because a pure
/// MoonBit library cannot sleep, the actual clock and the HTTP transfer are
/// injected by the host via the `now` and `fetch` callbacks of
/// `TogoWSClient::new`.

///|
/// Base URL of the TogoWS REST service.
pub const TOGOWS_BASE_URL : String = "http://togows.dbcls.jp"

///|
/// Delay enforced between two subsequent TogoWS requests (one third of a
/// second, i.e. up to three queries per second).
pub const TOGOWS_DELAY : Double = 0.333333333

///|
/// Error raised by the TogoWS client for invalid arguments, unsupported
/// databases / formats / conversions, or malformed server responses.
pub suberror TogoError {
  TogoError(String)
}

///|
fn[T] tgw_fail(msg : String) -> T raise TogoError {
  raise TogoError(msg)
}

// ---------------------------------------------------------------------------
// Character / string helpers
// ---------------------------------------------------------------------------

///|
fn tgw_is_space(c : Char) -> Bool {
  let cp = c.to_int()
  cp == 32 || (cp >= 9 && cp <= 13)
}

///|
/// Convert a string to an array of characters (code points).
fn tgw_chars(s : String) -> Array[Char] {
  let out : Array[Char] = Array(capacity=s.length())
  for ch in s {
    out.push(ch)
  }
  out
}

///|
/// Python str.strip(): remove leading and trailing whitespace.
fn tgw_strip(s : String) -> String {
  let chars = tgw_chars(s)
  let mut start = 0
  let mut end = chars.length()
  while start < end && tgw_is_space(chars[start]) {
    start = start + 1
  }
  while end > start && tgw_is_space(chars[end - 1]) {
    end = end - 1
  }
  let buf = StringBuilder::new()
  for i = start; i < end; i = i + 1 {
    buf.write_char(chars[i])
  }
  buf.to_string()
}

///|
/// Python str.split() with no arguments: split on runs of whitespace and
/// drop empty pieces.
fn tgw_split_whitespace(s : String) -> Array[String] {
  let out : Array[String] = []
  let mut current = StringBuilder::new()
  let mut has = false
  for ch in s {
    if tgw_is_space(ch) {
      if has {
        out.push(current.to_string())
        current = StringBuilder::new()
        has = false
      }
    } else {
      current.write_char(ch)
      has = true
    }
  }
  if has {
    out.push(current.to_string())
  }
  out
}

///|
/// Join identifiers with commas, like Python's ",".join(ids).
fn tgw_join_comma(ids : Array[String]) -> String {
  let buf = StringBuilder::new()
  for i, id in ids {
    if i > 0 {
      buf.write_string(",")
    }
    buf.write_string(id)
  }
  buf.to_string()
}

///|
/// Lexicographic comparison of two strings by code point, like Python's
/// str comparison used by sorted().
fn tgw_cmp_str(a : String, b : String) -> Int {
  let ac = tgw_chars(a)
  let bc = tgw_chars(b)
  let n = if ac.length() < bc.length() { ac.length() } else { bc.length() }
  for i = 0; i < n; i = i + 1 {
    let x = ac[i].to_int()
    let y = bc[i].to_int()
    if x != y {
      return if x < y { -1 } else { 1 }
    }
  }
  if ac.length() < bc.length() {
    -1
  } else if ac.length() > bc.length() {
    1
  } else {
    0
  }
}

///|
/// ", ".join(sorted(words)) for error messages.
fn tgw_sorted_join(words : Array[String]) -> String {
  let sorted = words.copy()
  sorted.sort_by(fn(a : String, b : String) -> Int { tgw_cmp_str(a, b) })
  let buf = StringBuilder::new()
  for i, w in sorted {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(w)
  }
  buf.to_string()
}

///|
/// Array equality for string lists.
fn tgw_eq(a : Array[String], b : Array[String]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if a[i] != b[i] {
      return false
    }
  }
  true
}

///|
/// Strict integer parsing (like Python int() on a stripped string):
/// optional sign followed by decimal digits.
fn tgw_parse_int(s : String) -> Int? {
  let chars = tgw_chars(s)
  let mut i = 0
  let mut sign = 1
  if chars.length() > 0 {
    let c = chars[0].to_int()
    if c == 43 {
      // '+'
      i = 1
    } else if c == 45 {
      // '-'
      sign = -1
      i = 1
    }
  }
  if i >= chars.length() {
    return None
  }
  let mut value = 0
  while i < chars.length() {
    let cp = chars[i].to_int()
    if cp < 48 || cp > 57 {
      return None
    }
    value = value * 10 + (cp - 48)
    i = i + 1
  }
  Some(sign * value)
}

///|
/// Two-digit uppercase hexadecimal representation of a byte.
fn tgw_hex2(byte : Int) -> String {
  let digits = tgw_chars("0123456789ABCDEF")
  let buf = StringBuilder::new()
  buf.write_char(digits[byte >> 4])
  buf.write_char(digits[byte & 0xF])
  buf.to_string()
}

///|
/// Percent-encode a single code point: ASCII bytes directly, larger code
/// points first encoded as UTF-8. Matches urllib.parse.quote.
fn tgw_quote_char(cp : Int, buf : StringBuilder) -> Unit {
  if cp < 0x80 {
    buf.write_string("%")
    buf.write_string(tgw_hex2(cp))
  } else if cp < 0x800 {
    tgw_quote_byte(0xC0 | (cp >> 6), buf)
    tgw_quote_byte(0x80 | (cp & 0x3F), buf)
  } else if cp < 0x10000 {
    tgw_quote_byte(0xE0 | (cp >> 12), buf)
    tgw_quote_byte(0x80 | ((cp >> 6) & 0x3F), buf)
    tgw_quote_byte(0x80 | (cp & 0x3F), buf)
  } else {
    tgw_quote_byte(0xF0 | (cp >> 18), buf)
    tgw_quote_byte(0x80 | ((cp >> 12) & 0x3F), buf)
    tgw_quote_byte(0x80 | ((cp >> 6) & 0x3F), buf)
    tgw_quote_byte(0x80 | (cp & 0x3F), buf)
  }
}

///|
fn tgw_quote_byte(byte : Int, buf : StringBuilder) -> Unit {
  buf.write_string("%")
  buf.write_string(tgw_hex2(byte))
}

///|
/// Percent-encode a string exactly like urllib.parse.quote: ASCII letters,
/// digits and '_.-~' are never encoded, characters listed in `safe` (a
/// default of "/" matches the original module's usage) are kept verbatim,
/// and everything else is encoded with UTF-8 percent escapes.
pub fn togows_quote(text~ : String, safe? : String) -> String {
  let safe_chars = safe.unwrap_or("/")
  let buf = StringBuilder::new()
  for ch in text {
    let cp = ch.to_int()
    let always_safe = (cp >= 65 && cp <= 90) ||
      (cp >= 97 && cp <= 122) ||
      (cp >= 48 && cp <= 57) ||
      cp == 95 ||
      cp == 46 ||
      cp == 45 ||
      cp == 126
    if always_safe {
      buf.write_char(ch)
    } else {
      // check the safe list (single-character entries)
      let mut is_safe = false
      for sc in safe_chars {
        if sc.to_int() == cp {
          is_safe = true
          break
        }
      }
      if is_safe {
        buf.write_char(ch)
      } else {
        tgw_quote_char(cp, buf)
      }
    }
  }
  buf.to_string()
}

// ---------------------------------------------------------------------------
// Pure URL builders
// ---------------------------------------------------------------------------

///|
/// Build the TogoWS 'entry' URL for a single identifier.
///
/// - db: database name, e.g. "pubmed", "nuccore", "uniprot", "compound".
/// - id: identifier (or a comma separated list of identifiers).
/// - format: optional return format, e.g. "fasta", "json", "xml", "ttl".
/// - field: optional field extraction, e.g. "title" for pubmed.
pub fn togows_entry_url(
  db~ : String,
  id~ : String,
  format? : String,
  field? : String,
) -> String {
  let mut url = TOGOWS_BASE_URL + "/entry/" + db + "/" + togows_quote(text=id)
  match field {
    Some(f) => url = url + "/" + f
    None => ()
  }
  match format {
    Some(f) => url = url + "." + f
    None => ()
  }
  url
}

///|
/// Build the TogoWS 'entry' URL for multiple identifiers (joined with
/// commas before quoting, like the original module).
pub fn togows_entry_multi_url(
  db~ : String,
  ids~ : Array[String],
  format? : String,
  field? : String,
) -> String {
  togows_entry_url(db~, id=tgw_join_comma(ids), format?, field?)
}

///|
/// Build the TogoWS search count URL.
pub fn togows_search_count_url(db~ : String, query~ : String) -> String {
  TOGOWS_BASE_URL + "/search/" + db + "/" + togows_quote(text=query) + "/count"
}

///|
/// Build the TogoWS search URL.
///
/// The offset (1 based) and limit must be provided together or not at all,
/// and both must be at least one - exactly like the original module.
pub fn togows_search_url(
  db~ : String,
  query~ : String,
  offset? : Int,
  limit? : Int,
  format? : String,
) -> String raise TogoError {
  let has_offset = offset is Some(_)
  let has_limit = limit is Some(_)
  if has_offset != has_limit {
    tgw_fail("Expect BOTH offset AND limit to be provided (or neither)")
  }
  let mut url = TOGOWS_BASE_URL +
    "/search/" +
    db +
    "/" +
    togows_quote(text=query)
  match (offset, limit) {
    (Some(o), Some(l)) => {
      if o <= 0 {
        tgw_fail("Offset should be at least one, not \{o}")
      }
      if l <= 0 {
        tgw_fail("Count should be at least one, not \{l}")
      }
      url = url + "/\{o},\{l}"
    }
    _ => ()
  }
  match format {
    Some(f) => url = url + "." + f
    None => ()
  }
  url
}

///|
/// Build the TogoWS conversion URL.
pub fn togows_convert_url(in_format~ : String, out_format~ : String) -> String {
  TOGOWS_BASE_URL + "/convert/" + in_format + "." + out_format
}

// ---------------------------------------------------------------------------
// Client with injectable fetch and clock
// ---------------------------------------------------------------------------

///|
/// A TogoWS REST client.
///
/// The HTTP transfer is provided by the host through `fetch`: it receives
/// the request URL and an optional POST body, and returns the response
/// body as text (raising `TogoError` on transport failures). The clock is
/// provided through `now` (seconds); it defaults to a virtual clock that
/// never advances, in which case every call is scheduled `TOGOWS_DELAY`
/// seconds after the previous one - matching the "up to three queries per
/// second" limit of the original module.
pub struct TogoWSClient {
  priv fetch : (String, String?) -> String raise TogoError
  priv now : () -> Double
  priv mut previous : Double
  priv mut entry_dbs : Array[String]?
  priv mut search_dbs : Array[String]?
  priv entry_fields : Map[String, Array[String]]
  priv entry_formats : Map[String, Array[String]]
  priv mut convert_formats : Array[(String, String)]?
}

///|
/// Create a TogoWS client.
///
/// - fetch: callback performing the HTTP transfer: (url, post_body) -> text.
/// - now: optional wall clock in seconds; defaults to a virtual clock.
pub fn TogoWSClient::new(
  fetch~ : (String, String?) -> String raise TogoError,
  now? : () -> Double,
) -> TogoWSClient {
  {
    fetch,
    now: now.unwrap_or(fn() { 0.0 }),
    previous: 0.0,
    entry_dbs: None,
    search_dbs: None,
    entry_fields: Map([], capacity=8),
    entry_formats: Map([], capacity=8),
    convert_formats: None,
  }
}

///|
/// The scheduled time of the most recent (or currently pending) request.
///
/// Between two calls at most TOGOWS_DELAY seconds may pass; the rate
/// limiter records the earliest permitted time of the next request here.
pub fn TogoWSClient::last_call_time(self : TogoWSClient) -> Double {
  self.previous
}

///|
/// Perform a request through the rate limiter.
///
/// This mirrors the private _open() of the original module: a request may
/// start once TOGOWS_DELAY seconds have elapsed since the previous one.
fn TogoWSClient::open_url(
  self : TogoWSClient,
  url : String,
  post? : String,
) -> String raise TogoError {
  let current = (self.now)()
  let wait = self.previous + TOGOWS_DELAY - current
  if wait > 0.0 {
    // The original sleeps for `wait` here and records current + wait.
    self.previous = current + wait
  } else {
    self.previous = current
  }
  // like the original: an empty POST body is treated as no body
  let body : String? = match post {
    Some(p) => if p.length() > 0 { Some(p) } else { None }
    None => None
  }
  (self.fetch)(url, body)
}

///|
/// Fetch a plain text list of values from a TogoWS URL (whitespace split).
fn TogoWSClient::get_fields(
  self : TogoWSClient,
  url : String,
) -> Array[String] raise TogoError {
  tgw_split_whitespace(tgw_strip(self.open_url(url)))
}

///|
/// The list of databases supported by the TogoWS entry service.
pub fn TogoWSClient::get_entry_dbs(
  self : TogoWSClient,
) -> Array[String] raise TogoError {
  match self.entry_dbs {
    Some(dbs) => dbs
    None => {
      let dbs = self.get_fields(TOGOWS_BASE_URL + "/entry")
      self.entry_dbs = Some(dbs)
      dbs
    }
  }
}

///|
/// The list of databases supported by the TogoWS search service.
pub fn TogoWSClient::get_search_dbs(
  self : TogoWSClient,
) -> Array[String] raise TogoError {
  match self.search_dbs {
    Some(dbs) => dbs
    None => {
      let dbs = self.get_fields(TOGOWS_BASE_URL + "/search")
      self.search_dbs = Some(dbs)
      dbs
    }
  }
}

///|
/// The fields available for field extraction within an entry database.
pub fn TogoWSClient::get_entry_fields(
  self : TogoWSClient,
  db~ : String,
) -> Array[String] raise TogoError {
  match self.entry_fields.get(db) {
    Some(fields) => fields
    None => {
      let fields = self.get_fields(TOGOWS_BASE_URL + "/entry/" + db + "?fields")
      self.entry_fields.set(db, fields)
      fields
    }
  }
}

///|
/// The return formats supported by an entry database.
pub fn TogoWSClient::get_entry_formats(
  self : TogoWSClient,
  db~ : String,
) -> Array[String] raise TogoError {
  match self.entry_formats.get(db) {
    Some(formats) => formats
    None => {
      let formats = self.get_fields(
        TOGOWS_BASE_URL + "/entry/" + db + "?formats",
      )
      self.entry_formats.set(db, formats)
      formats
    }
  }
}

///|
/// The supported format conversion pairs, e.g. ("genbank", "fasta").
pub fn TogoWSClient::get_convert_formats(
  self : TogoWSClient,
) -> Array[(String, String)] raise TogoError {
  match self.convert_formats {
    Some(pairs) => pairs
    None => {
      let tokens = self.get_fields(TOGOWS_BASE_URL + "/convert/")
      let pairs : Array[(String, String)] = []
      for token in tokens {
        let parts = tgw_split_dot(token)
        // tokens without exactly one dot can never match an (in, out) pair
        if parts.length() == 2 {
          pairs.push((parts[0], parts[1]))
        }
      }
      self.convert_formats = Some(pairs)
      pairs
    }
  }
}

///|
/// Split "a.b" on the single dot, like Python's pair.split(".").
fn tgw_split_dot(token : String) -> Array[String] {
  let out : Array[String] = []
  let mut current = StringBuilder::new()
  for ch in token {
    if ch.to_int() == 46 {
      out.push(current.to_string())
      current = StringBuilder::new()
    } else {
      current.write_char(ch)
    }
  }
  out.push(current.to_string())
  out
}

///|
/// Call TogoWS 'entry' to fetch a record.
///
/// - db: database, e.g. "pubmed", "nuccore", "uniprot", "compound", ...
/// - id: identifier (or comma separated identifiers).
/// - format: optional return format, e.g. "xml", "json", "gff", "fasta",
///   "ttl" (RDF Turtle).
/// - field: optional field extraction, e.g. "title" for pubmed.
///
/// The database name, field and format are validated against the lists
/// published by the TogoWS service (fetched once and cached). Like the
/// original module, the pubmed field alias "ti" is silently mapped to
/// "title" for backwards compatibility.
pub fn TogoWSClient::entry(
  self : TogoWSClient,
  db~ : String,
  id~ : String,
  format? : String,
  field? : String,
) -> String raise TogoError {
  self.entry_impl(db~, ids=[id], format?, field?)
}

///|
/// Call TogoWS 'entry' with multiple identifiers (joined with commas).
pub fn TogoWSClient::entry_multi(
  self : TogoWSClient,
  db~ : String,
  ids~ : Array[String],
  format? : String,
  field? : String,
) -> String raise TogoError {
  self.entry_impl(db~, ids~, format?, field?)
}

///|
fn TogoWSClient::entry_impl(
  self : TogoWSClient,
  db~ : String,
  ids~ : Array[String],
  format? : String,
  field? : String,
) -> String raise TogoError {
  let dbs = self.get_entry_dbs()
  if !dbs.contains(db) {
    tgw_fail("TogoWS entry fetch does not officially support database '\{db}'.")
  }
  let mut field = field
  match field {
    Some(f) =>
      if f.length() > 0 {
        let fields = self.get_entry_fields(db~)
        let mut fixed = f
        if db == "pubmed" && f == "ti" && fields.contains("title") {
          // Backwards compatibility fix for TogoWS change Nov/Dec 2013
          fixed = "title"
        }
        if !fields.contains(fixed) {
          tgw_fail(
            "TogoWS entry fetch does not explicitly support field '\{fixed}' for database '\{db}'. Only: \{tgw_sorted_join(fields)}",
          )
        }
        field = Some(fixed)
      } else {
        field = None
      }
    None => ()
  }
  match format {
    Some(fm) =>
      if fm.length() > 0 {
        let formats = self.get_entry_formats(db~)
        if !formats.contains(fm) {
          tgw_fail(
            "TogoWS entry fetch does not explicitly support format '\{fm}' for database '\{db}'. Only: \{tgw_sorted_join(formats)}",
          )
        }
      }
    None => ()
  }
  let id = tgw_join_comma(ids)
  let mut url = TOGOWS_BASE_URL + "/entry/" + db + "/" + togows_quote(text=id)
  match field {
    Some(f) => url = url + "/" + f
    None => ()
  }
  match format {
    Some(fm) => if fm.length() > 0 { url = url + "." + fm }
    None => ()
  }
  self.open_url(url)
}

///|
/// Call TogoWS search count to see how many matches a search gives.
pub fn TogoWSClient::search_count(
  self : TogoWSClient,
  db~ : String,
  query~ : String,
) -> Int raise TogoError {
  let url = togows_search_count_url(db~, query~)
  let data = self.open_url(url)
  if data.length() == 0 {
    tgw_fail("TogoWS returned no data from URL \{url}")
  }
  let stripped = tgw_strip(data)
  match tgw_parse_int(stripped) {
    Some(n) => n
    None => tgw_fail("Expected an integer from URL \{url}, got: '\{data}'")
  }
}

///|
/// Call TogoWS search.
///
/// - db: database, e.g. "pubmed", "ncbi-genbank", "uniprot", "taxonomy".
/// - query: search term.
/// - offset / limit: 1-based paging; both must be given together and be at
///   least one. TogoWS applies a default count limit of 100 results.
/// - format: optional return format, e.g. "json", "ttl"; by default plain
///   text with one result per line is returned.
pub fn TogoWSClient::search(
  self : TogoWSClient,
  db~ : String,
  query~ : String,
  offset? : Int,
  limit? : Int,
  format? : String,
) -> String raise TogoError {
  let url = togows_search_url(db~, query~, offset?, limit?, format?)
  self.open_url(url)
}

///|
/// Call TogoWS search iterating over the results in batches.
///
/// This first calls search_count and then pulls the identifiers in batches
/// (at most 100 per request, enforced by TogoWS). The result ordering
/// mirrors the original generator: an optional overall limit, error on a
/// batch of unexpected size, on a repeated batch, or on an identifier that
/// already appeared in the previous batch.
pub fn TogoWSClient::search_iter(
  self : TogoWSClient,
  db~ : String,
  query~ : String,
  limit? : Int,
  batch? : Int,
) -> Array[String] raise TogoError {
  let batch_size = batch.unwrap_or(100)
  let count = self.search_count(db~, query~)
  let out : Array[String] = []
  if count == 0 {
    return out
  }
  let mut remain = count
  match limit {
    Some(l) => if l < remain { remain = l }
    None => ()
  }
  let mut offset = 1 // TogoWS does not use zero based counting
  let mut prev_ids : Array[String] = []
  while remain > 0 {
    let size = if batch_size < remain { batch_size } else { remain }
    let text = self.search(db~, query~, offset~, limit=size)
    let ids = tgw_split_whitespace(tgw_strip(text))
    if ids.length() != size {
      tgw_fail("Got \{ids.length()}, expected \{size}")
    }
    if tgw_eq(ids, prev_ids) {
      tgw_fail("Same search results for previous offset")
    }
    for identifier in ids {
      if prev_ids.contains(identifier) {
        tgw_fail("Result \{identifier} was in previous batch")
      }
      out.push(identifier)
    }
    offset = offset + size
    remain = remain - size
    prev_ids = ids
  }
  out
}

///|
/// Call TogoWS for file format conversion.
///
/// - data: input record(s) as a string (sent as the POST body).
/// - in_format: input file format, e.g. "genbank".
/// - out_format: requested output format, e.g. "fasta".
///
/// The pair is validated against the conversions published by the TogoWS
/// service (fetched once and cached).
pub fn TogoWSClient::convert(
  self : TogoWSClient,
  data~ : String,
  in_format~ : String,
  out_format~ : String,
) -> String raise TogoError {
  let pairs = self.get_convert_formats()
  let mut supported = false
  for pair in pairs {
    if pair.0 == in_format && pair.1 == out_format {
      supported = true
      break
    }
  }
  if !supported {
    let buf = StringBuilder::new()
    for i, pair in pairs {
      if i > 0 {
        buf.write_string("\n")
      }
      buf.write_string("\{pair.0} -> \{pair.1}")
    }
    tgw_fail("Unsupported conversion. Choose from:\n\{buf.to_string()}")
  }
  let url = togows_convert_url(in_format~, out_format~)
  self.open_url(url, post=data)
}
