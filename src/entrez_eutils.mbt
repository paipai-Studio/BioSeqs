///|
/// Bio.Entrez E-utilities - the complete NCBI Entrez Programming Utilities
/// request layer that complements the simulated access in entrez.mbt.
///
/// Port of Biopython's Bio.Entrez request semantics (urllib-based, 1.88):
///   - epost:    POST a file of identifiers for future use; returns
///               QueryKey/WebEnv for use with efetch/esummary.
///   - esummary: retrieve document summaries for primary IDs.
///   - espell:   retrieve spelling suggestions.
///   - elink:    check for linked records / Related Articles (list ids are
///               transmitted as repeated "id" parameters, like the original
///               module's join_ids=False behaviour).
///   - ecitmatch: retrieve PMIDs for input citation strings; the citation
///               bdata uses "|" separators (not URL-encoded) and "\r"
///               between citations, exactly like the original.
///
/// Parameters are appended in the same order as the original: user
/// parameters first, then tool (default "biopython"), email and api_key if
/// configured. URL encoding follows urllib.parse.urlencode (quote_plus:
/// space -> '+').
///
/// The HTTP transfer and the clock are injected, so all behaviour can be
/// verified offline. NCBI rate limiting (3 queries/second, or 10 with an
/// api_key) is accounted on the injected clock.

///|
/// Error raised by the E-utilities client for invalid arguments, missing
/// response fields or transport failures.
pub suberror EntrezUtilError {
  EntrezUtilError(String)
}

///|
fn[T] entrez_util_fail(msg : String) -> T raise EntrezUtilError {
  raise EntrezUtilError(msg)
}

///|
/// E-utilities client with an injectable transport and clock.
pub struct EntrezClient {
  priv fetch : (String, Option[String]) -> String raise EntrezUtilError
  priv now : () -> Double
  priv mut previous : Double
  priv tool : String
  priv email : String
  priv api_key : String
}

///|
/// Create an E-utilities client.
///
/// - fetch: callback performing the HTTP transfer: (url, post_body) -> text.
/// - now: optional wall clock in seconds; defaults to a virtual clock.
/// - tool: tool identifier, default "biopython".
/// - email: optional e-mail address declared to NCBI.
/// - api_key: optional NCBI API key (raises the rate limit to 10/s).
pub fn EntrezClient::new(
  fetch~ : (String, Option[String]) -> String raise EntrezUtilError,
  now? : () -> Double,
  tool? : String,
  email? : String,
  api_key? : String,
) -> EntrezClient {
  {
    fetch,
    now: now.unwrap_or(fn() { 0.0 }),
    previous: 0.0,
    tool: tool.unwrap_or("biopython"),
    email: email.unwrap_or(""),
    api_key: api_key.unwrap_or(""),
  }
}

///|
/// Base URL of the NCBI E-utilities.
pub let entrez_eutils_base_url : String =
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

///|
/// URL-encode a string following urllib.parse.urlencode's quote_plus
/// (space becomes '+', unreserved characters kept, everything else
/// percent-encoded over UTF-8 bytes).
pub fn entrez_quote_plus(text : String) -> String {
  let buf = StringBuilder::new()
  for ch in text {
    let cp = ch.to_int()
    let unreserved = (cp >= 65 && cp <= 90) ||
      (cp >= 97 && cp <= 122) ||
      (cp >= 48 && cp <= 57) ||
      cp == 95 ||
      cp == 46 ||
      cp == 45 ||
      cp == 126
    if cp == 32 {
      buf.write_char('+')
    } else if unreserved {
      buf.write_char(ch)
    } else {
      tgw_quote_char(cp, buf)
    }
  }
  buf.to_string()
}

///|
/// Enforce NCBI rate limiting on the injected clock: at most three queries
/// per second without an API key (0.37 s pause like the original, which is
/// more reliable than the theoretical 0.333 s), 0.1 s with one.
fn EntrezClient::rate_wait(self : EntrezClient) -> Unit {
  let delay = if self.api_key.length() > 0 { 0.1 } else { 0.37 }
  let current = (self.now)()
  let wait = self.previous + delay - current
  if wait > 0 {
    self.previous = current + wait
  } else {
    self.previous = current
  }
}

///|
/// Append the standard identifying parameters (tool/email/api_key) after
/// the user parameters, mirroring _construct_params' setdefault order.
fn EntrezClient::standard_params(
  self : EntrezClient,
  params : Array[(String, String)],
) -> Array[(String, String)] {
  let out = params.copy()
  let mut has_tool = false
  let mut has_email = false
  let mut has_key = false
  for kv in out {
    let (k, _) = kv
    if k == "tool" {
      has_tool = true
    }
    if k == "email" {
      has_email = true
    }
    if k == "api_key" {
      has_key = true
    }
  }
  if !has_tool {
    out.push(("tool", self.tool))
  }
  if !has_email && self.email.length() > 0 {
    out.push(("email", self.email))
  }
  if !has_key && self.api_key.length() > 0 {
    out.push(("api_key", self.api_key))
  }
  out
}

///|
/// Build the urlencoded parameter string for a request.
fn entrez_build_query(
  params : Array[(String, String)],
  ecitmatch : Bool,
) -> String {
  let buf = StringBuilder::new()
  let mut first = true
  for kv in params {
    let (k, v) = kv
    if !first {
      buf.write_char('&')
    }
    first = false
    buf.write_string(entrez_quote_plus(k))
    buf.write_char('=')
    buf.write_string(entrez_quote_plus(v))
  }
  let mut qs = buf.to_string()
  if ecitmatch {
    // Don't URL-encode pipe ("|") characters, this is expected by the
    // ecitmatch tool.
    qs = entrez_replace(qs, "%7C", "|")
  }
  qs
}

///|
/// Replace all occurrences of a substring (used for the ecitmatch pipes).
fn entrez_replace(text : String, from : String, to : String) -> String {
  if from.length() == 0 {
    return text
  }
  let buf = StringBuilder::new()
  let mut i = 0
  while i < text.length() {
    if i + from.length() <= text.length() &&
      text[i:(i + from.length())].to_string() == from {
      buf.write_string(to)
      i = i + from.length()
    } else {
      // copy the character verbatim; char_at would render characters
      // outside its whitelist (e.g. '=', '%') as decimal code strings
      buf.write_string(text[i:(i + 1)].to_owned())
      i = i + 1
    }
  }
  buf.to_string()
}

///|
/// Issue a request through the client's transport.
///
/// - post: force HTTP POST; when omitted, POST is used if the encoded
///   parameters would be over 1000 characters long, or if there are 200 or
///   more comma-separated ids (NCBI's recommendation).
fn EntrezClient::request(
  self : EntrezClient,
  cgi : String,
  params : Array[(String, String)],
  post? : Bool,
  ecitmatch? : Bool,
) -> String raise EntrezUtilError {
  self.rate_wait()
  let all = self.standard_params(params)
  let qs = entrez_build_query(all, ecitmatch.unwrap_or(false))
  let mut do_post = post.unwrap_or(false)
  if post is None {
    if qs.length() > 1000 {
      do_post = true
    } else {
      for kv in all {
        let (k, v) = kv
        if k == "id" {
          // NCBI prefers an HTTP POST instead of an HTTP GET if there are
          // more than about 200 IDs
          let mut count = 1
          for ch in v {
            if ch == ',' {
              count = count + 1
            }
          }
          if count >= 200 {
            do_post = true
          }
        }
      }
    }
  }
  if do_post {
    (self.fetch)(cgi, Some(qs))
  } else {
    (self.fetch)(cgi + "?" + qs, None)
  }
}

// ---------------------------------------------------------------------------
// XML response parsing helpers
// ---------------------------------------------------------------------------

///|
/// Extract the text of the first occurrence of a tag, tolerating attributes
/// (e.g. `<QueryKey key="1">1</QueryKey>`).
fn entrez_extract_text(xml : String, tag : String) -> String? {
  let open = "<" + tag
  let start_pos = xml.find(open)
  if start_pos is None {
    return None
  }
  let after_open = xml[start_pos.unwrap():].to_owned()
  let gt = after_open.find(">")
  if gt is None {
    return None
  }
  let after_gt = after_open[gt.unwrap() + 1:].to_owned()
  let end_pos = after_gt.find("</" + tag + ">")
  if end_pos is None {
    return None
  }
  Some(after_gt[0:end_pos.unwrap()].to_owned())
}

///|
/// Extract the text of every occurrence of a tag.
fn entrez_extract_all(xml : String, tag : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut remaining = xml
  while remaining.length() > 0 {
    let start_pos = remaining.find("<" + tag + ">")
    if start_pos is None {
      break
    }
    let start = start_pos.unwrap() + tag.length() + 2
    let after_start = remaining[start:].to_owned()
    let end_pos = after_start.find("</" + tag + ">")
    if end_pos is None {
      break
    }
    result.push(after_start[0:end_pos.unwrap()].to_owned())
    remaining = after_start[end_pos.unwrap() + tag.length() + 3:].to_owned()
  }
  result
}

///|
/// Extract an attribute value from a tag's opening text, e.g.
/// entrez_attr("<Item Name=\"PdbDescr\" Type=\"String\">", "Name") -> "PdbDescr".
fn entrez_attr(opening : String, name : String) -> String? {
  let needle = name + "=\""
  let start_pos = opening.find(needle)
  if start_pos is None {
    return None
  }
  let start = start_pos.unwrap() + needle.length()
  let rest = opening[start:].to_owned()
  let end_pos = rest.find("\"")
  if end_pos is None {
    return None
  }
  Some(rest[0:end_pos.unwrap()].to_owned())
}

///|
/// Cut out every occurrence of a `<tag ...>...</tag>` block.
fn entrez_blocks(xml : String, tag : String) -> Array[String] {
  let result : Array[String] = Array::new()
  let mut remaining = xml
  while remaining.length() > 0 {
    let start_pos = remaining.find("<" + tag)
    if start_pos is None {
      break
    }
    let after_open = remaining[start_pos.unwrap():].to_owned()
    let gt = after_open.find(">")
    if gt is None {
      break
    }
    let body = after_open[gt.unwrap() + 1:].to_owned()
    let end_pos = body.find("</" + tag + ">")
    if end_pos is None {
      break
    }
    result.push(body[0:end_pos.unwrap()].to_owned())
    remaining = body[end_pos.unwrap():].to_owned()
  }
  result
}

// ---------------------------------------------------------------------------
// epost
// ---------------------------------------------------------------------------

///|
/// Result of an EPost request: the WebEnv session and query key to use in
/// subsequent efetch/esummary/esearch calls.
pub struct EPostResult {
  query_key : String
  web_env : String
}

///|
/// Create an EPostResult.
pub fn EPostResult::new(query_key~ : String, web_env~ : String) -> EPostResult {
  { query_key, web_env }
}

///|
/// Get the QueryKey.
pub fn EPostResult::query_key(self : EPostResult) -> String {
  self.query_key
}

///|
/// Get the WebEnv.
pub fn EPostResult::web_env(self : EPostResult) -> String {
  self.web_env
}

///|
/// EPost: post a file of identifiers for future use in the user's
/// environment. Always uses HTTP POST. The ids are transmitted as a single
/// comma-separated "id" parameter.
pub fn epost(
  client : EntrezClient,
  db : String,
  ids : Array[String],
) -> EPostResult raise EntrezUtilError {
  let params : Array[(String, String)] = [("db", db), ("id", join_comma(ids))]
  let xml = client.request(entrez_eutils_base_url + "epost.fcgi", params, post=true)
  let query_key = match entrez_extract_text(xml, "QueryKey") {
    Some(v) => v
    None => entrez_util_fail("epost response is missing QueryKey")
  }
  let web_env = match entrez_extract_text(xml, "WebEnv") {
    Some(v) => v
    None => entrez_util_fail("epost response is missing WebEnv")
  }
  EPostResult::new(query_key~, web_env~)
}

///|
/// Join ids with commas (the E-utilities id list convention).
fn join_comma(ids : Array[String]) -> String {
  let buf = StringBuilder::new()
  let mut first = true
  for id in ids {
    if !first {
      buf.write_char(',')
    }
    first = false
    buf.write_string(id)
  }
  buf.to_string()
}

// ---------------------------------------------------------------------------
// esummary
// ---------------------------------------------------------------------------

///|
/// One `<Item Name=... Type=...>value</Item>` of a document summary.
pub struct EntrezItem {
  name : String
  type_ : String
  value : String
}

///|
/// Create an EntrezItem.
pub fn EntrezItem::new(name~ : String, type_~ : String, value~ : String) -> EntrezItem {
  { name, type_, value }
}

///|
/// Get the item name.
pub fn EntrezItem::name(self : EntrezItem) -> String {
  self.name
}

///|
/// Get the item type.
pub fn EntrezItem::type_(self : EntrezItem) -> String {
  self.type_
}

///|
/// Get the item value.
pub fn EntrezItem::value(self : EntrezItem) -> String {
  self.value
}

///|
/// A document summary (DocSum) from ESummary.
pub struct EntrezDocSum {
  id : String
  items : Array[EntrezItem]
}

///|
/// Create an EntrezDocSum.
pub fn EntrezDocSum::new(id~ : String, items~ : Array[EntrezItem]) -> EntrezDocSum {
  { id, items }
}

///|
/// Get the summary's primary id.
pub fn EntrezDocSum::id(self : EntrezDocSum) -> String {
  self.id
}

///|
/// Get the summary's items (in document order).
pub fn EntrezDocSum::items(self : EntrezDocSum) -> Array[EntrezItem] {
  self.items
}

///|
/// Look up the value of an item by name.
pub fn EntrezDocSum::get(self : EntrezDocSum, name : String) -> String? {
  for item in self.items {
    if item.name == name {
      return Some(item.value)
    }
  }
  None
}

///|
/// ESummary: retrieve document summaries for a list of primary IDs.
pub fn esummary(
  client : EntrezClient,
  db : String,
  ids : Array[String],
) -> Array[EntrezDocSum] raise EntrezUtilError {
  let params : Array[(String, String)] = [("db", db), ("id", join_comma(ids))]
  let xml = client.request(entrez_eutils_base_url + "esummary.fcgi", params)
  parse_esummary(xml)
}

///|
/// Scan a DocSum block for `<Item Name=... Type=...>value</Item>` entries,
/// returning (name, type, value) triples in document order.
fn entrez_scan_items(block : String) -> Array[(String, String, String)] {
  let result : Array[(String, String, String)] = Array::new()
  let mut rest = block
  let mut go = true
  while go {
    go = false
    let pos = rest.find("<Item ")
    if pos is Some(p) {
      let after = rest[p:].to_owned()
      let gt = after.find(">")
      let end_pos = after.find("</Item>")
      if gt is Some(g) && end_pos is Some(e) && e > g {
        let opening = after[0:(g + 1)].to_owned()
        let value = after[(g + 1):e].to_owned()
        let name = entrez_attr(opening, "Name").unwrap_or("")
        let type_ = entrez_attr(opening, "Type").unwrap_or("")
        result.push((name, type_, value))
        rest = after[(e + 7):].to_owned()
        go = true
      }
    }
  }
  result
}

///|
/// Parse an ESummary XML response into document summaries.
pub fn parse_esummary(xml : String) -> Array[EntrezDocSum] raise EntrezUtilError {
  let result : Array[EntrezDocSum] = Array::new()
  for block in entrez_blocks(xml, "DocSum") {
    let id = match entrez_extract_text(block, "Id") {
      Some(v) => v
      None => entrez_util_fail("DocSum is missing Id")
    }
    let items : Array[EntrezItem] = Array::new()
    for triple in entrez_scan_items(block) {
      let (name, type_, value) = triple
      items.push(EntrezItem::new(name~, type_~, value~))
    }
    result.push(EntrezDocSum::new(id~, items~))
  }
  result
}

// ---------------------------------------------------------------------------
// espell
// ---------------------------------------------------------------------------

///|
/// Result of an ESpell request.
pub struct ESpellResult {
  database : String
  query : String
  corrected_query : String
  replaced : Array[String]
}

///|
/// Create an ESpellResult.
pub fn ESpellResult::new(
  database~ : String,
  query~ : String,
  corrected_query~ : String,
  replaced~ : Array[String],
) -> ESpellResult {
  { database, query, corrected_query, replaced }
}

///|
/// Get the corrected query (empty when no correction is available).
pub fn ESpellResult::corrected_query(self : ESpellResult) -> String {
  self.corrected_query
}

///|
/// Get the original query.
pub fn ESpellResult::query(self : ESpellResult) -> String {
  self.query
}

///|
/// Get the individual replaced terms from the SpelledQuery block.
pub fn ESpellResult::replaced(self : ESpellResult) -> Array[String] {
  self.replaced
}

///|
/// ESpell: retrieve spelling suggestions for a query term.
pub fn espell(
  client : EntrezClient,
  db? : String,
  term~ : String,
) -> ESpellResult raise EntrezUtilError {
  let params : Array[(String, String)] = []
  if db is Some(d) {
    params.push(("db", d))
  }
  params.push(("term", term))
  let xml = client.request(entrez_eutils_base_url + "espell.fcgi", params)
  parse_espell(xml)
}

///|
/// Parse an ESpell XML response.
pub fn parse_espell(xml : String) -> ESpellResult raise EntrezUtilError {
  let query = match entrez_extract_text(xml, "Query") {
    Some(v) => v
    None => entrez_util_fail("espell response is missing Query")
  }
  let corrected = entrez_extract_text(xml, "CorrectedQuery").unwrap_or("")
  let replaced = entrez_extract_all(xml, "Replaced")
  let database = entrez_extract_text(xml, "Database").unwrap_or("")
  ESpellResult::new(
    database~,
    query~,
    corrected_query=corrected,
    replaced~,
  )
}

// ---------------------------------------------------------------------------
// elink
// ---------------------------------------------------------------------------

///|
/// One LinkSetDb block of an ELink response: destination database, link
/// name and the linked ids.
pub struct EntrezLinkSetDb {
  db_to : String
  link_name : String
  ids : Array[String]
}

///|
/// Create an EntrezLinkSetDb.
pub fn EntrezLinkSetDb::new(
  db_to~ : String,
  link_name~ : String,
  ids~ : Array[String],
) -> EntrezLinkSetDb {
  { db_to, link_name, ids }
}

///|
/// Get the destination database.
pub fn EntrezLinkSetDb::db_to(self : EntrezLinkSetDb) -> String {
  self.db_to
}

///|
/// Get the link name, e.g. "pubmed_pubmed".
pub fn EntrezLinkSetDb::link_name(self : EntrezLinkSetDb) -> String {
  self.link_name
}

///|
/// Get the linked ids.
pub fn EntrezLinkSetDb::ids(self : EntrezLinkSetDb) -> Array[String] {
  self.ids
}

///|
/// One LinkSet of an ELink response.
pub struct EntrezLinkSet {
  db_from : String
  link_set_dbs : Array[EntrezLinkSetDb]
}

///|
/// Create an EntrezLinkSet.
pub fn EntrezLinkSet::new(
  db_from~ : String,
  link_set_dbs~ : Array[EntrezLinkSetDb],
) -> EntrezLinkSet {
  { db_from, link_set_dbs }
}

///|
/// Get the source database.
pub fn EntrezLinkSet::db_from(self : EntrezLinkSet) -> String {
  self.db_from
}

///|
/// Get the LinkSetDb blocks of this link set.
pub fn EntrezLinkSet::link_set_dbs(self : EntrezLinkSet) -> Array[EntrezLinkSetDb] {
  self.link_set_dbs
}

///|
/// ELink: check for linked records. Multiple source ids are transmitted as
/// repeated "id" parameters providing a one-to-one mapping from source UIDs
/// to destination UIDs (like the original module's join_ids=False).
pub fn elink(
  client : EntrezClient,
  dbfrom? : String,
  db? : String,
  ids~ : Array[String],
  cmd? : String,
  linkname? : String,
) -> Array[EntrezLinkSet] raise EntrezUtilError {
  let params : Array[(String, String)] = []
  if dbfrom is Some(f) {
    params.push(("dbfrom", f))
  }
  if db is Some(d) {
    params.push(("db", d))
  }
  for id in ids {
    params.push(("id", id))
  }
  if cmd is Some(c) {
    params.push(("cmd", c))
  }
  if linkname is Some(l) {
    params.push(("linkname", l))
  }
  let xml = client.request(entrez_eutils_base_url + "elink.fcgi", params)
  parse_elink(xml)
}

///|
/// Parse an ELink XML response.
pub fn parse_elink(xml : String) -> Array[EntrezLinkSet] raise EntrezUtilError {
  let result : Array[EntrezLinkSet] = Array::new()
  for block in entrez_blocks(xml, "LinkSet") {
    let db_from = entrez_extract_text(block, "DbFrom").unwrap_or("")
    let link_set_dbs : Array[EntrezLinkSetDb] = Array::new()
    for lsd in entrez_blocks(block, "LinkSetDb") {
      let db_to = entrez_extract_text(lsd, "DbTo").unwrap_or("")
      let link_name = entrez_extract_text(lsd, "LinkName").unwrap_or("")
      let ids : Array[String] = Array::new()
      for link in entrez_blocks(lsd, "Link") {
        match entrez_extract_text(link, "Id") {
          Some(v) => ids.push(v)
          None => ()
        }
      }
      link_set_dbs.push(EntrezLinkSetDb::new(db_to~, link_name~, ids~))
    }
    result.push(EntrezLinkSet::new(db_from~, link_set_dbs~))
  }
  result
}

// ---------------------------------------------------------------------------
// ecitmatch
// ---------------------------------------------------------------------------

///|
/// A citation for ECitMatch: the six fields are joined with "|" in the
/// given order (journal_title, year, volume, first_page, author_name, key).
pub struct EntrezCitation {
  journal_title : String
  year : String
  volume : String
  first_page : String
  author_name : String
  key : String
}

///|
/// Create an EntrezCitation.
pub fn EntrezCitation::new(
  journal_title~ : String,
  year~ : String,
  volume~ : String,
  first_page~ : String,
  author_name~ : String,
  key~ : String,
) -> EntrezCitation {
  { journal_title, year, volume, first_page, author_name, key }
}

///|
/// Format the citation for the bdata parameter.
pub fn EntrezCitation::bdata_field(self : EntrezCitation) -> String {
  self.journal_title +
  "|" +
  self.year +
  "|" +
  self.volume +
  "|" +
  self.first_page +
  "|" +
  self.author_name +
  "|" +
  self.key
}

///|
/// ECitMatch: retrieve PubMed IDs (PMIDs) for input citation strings.
/// Returns (key, pmid) pairs in response order.
pub fn ecitmatch(
  client : EntrezClient,
  db : String,
  citations~ : Array[EntrezCitation],
) -> Array[(String, String)] raise EntrezUtilError {
  // XML is the only supported value, and it actually returns TXT.
  let params : Array[(String, String)] = [("retmode", "xml"), ("db", db)]
  let buf = StringBuilder::new()
  let mut first = true
  for cit in citations {
    if !first {
      // citations are joined by a carriage return
      buf.write_char('\r')
    }
    first = false
    buf.write_string(cit.bdata_field())
  }
  params.push(("bdata", buf.to_string()))
  let text = client.request(
    entrez_eutils_base_url + "ecitmatch.cgi",
    params,
    ecitmatch=true,
  )
  parse_ecitmatch(text)
}

///|
/// Parse an ECitMatch response. Supports both the XML form
/// (`<ByteId key="..." id="..."/>`) and the plain-text pipe form
/// (`journal|year|volume|page|author|key|pmid`).
pub fn parse_ecitmatch(text : String) -> Array[(String, String)] raise EntrezUtilError {
  let result : Array[(String, String)] = Array::new()
  match text.find("<") {
    Some(_) => {
      // scan every <ByteId key="..." id="..."/> element in the document
      let mut remaining = text
      while remaining.length() > 0 {
        let pos = remaining.find("<ByteId ")
        if pos is None {
          break
        }
        let rest = remaining[pos.unwrap():].to_owned()
        let gt = rest.find(">")
        if gt is None {
          break
        }
        let opening = rest[0:(gt.unwrap() + 1)].to_owned()
        let key = entrez_attr(opening, "key").unwrap_or("")
        let id = entrez_attr(opening, "id").unwrap_or("")
        result.push((key, id))
        remaining = rest[(gt.unwrap() + 1):].to_owned()
      }
    }
    None => {
      for line in split_by_char(text, '\n'.to_int()) {
        let trimmed = line.trim().to_owned()
        if trimmed.length() == 0 {
          continue
        }
        let fields = split_by_char(trimmed, '|'.to_int())
        if fields.length() < 2 {
          entrez_util_fail("cannot parse ecitmatch line: " + trimmed)
        }
        result.push((fields[fields.length() - 2], fields[fields.length() - 1]))
      }
    }
  }
  result
}
