///|
/// Bio.Entrez E-utilities - the complete NCBI Entrez Programming Utilities
/// request layer that complements the simulated access in entrez.mbt.
///
/// Port of Biopython's Bio.Entrez request semantics (urllib-based, 1.88):
///   - einfo:    list databases, or a database's statistics/field/link names.
///   - esearch:  search a database for primary IDs, term translation and
///               optional QueryKey/WebEnv history.
///   - efetch:   retrieve records by IDs (or history), raw text returned.
///   - epost:    POST a file of identifiers for future use; returns
///               QueryKey/WebEnv for use with efetch/esummary.
///   - esummary: retrieve document summaries for primary IDs.
///   - espell:   retrieve spelling suggestions.
///   - elink:    check for linked records / Related Articles (list ids are
///               transmitted as repeated "id" parameters, like the original
///               module's join_ids=False behaviour).
///   - egquery:  global cross-database counts for one search term.
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
  priv fetch : (String, String?) -> String raise EntrezUtilError
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
  fetch~ : (String, String?) -> String raise EntrezUtilError,
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
pub let entrez_eutils_base_url : String = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

///|
/// URL-encode a string following urllib.parse.urlencode's quote_plus
/// (space becomes '+', unreserved characters kept, everything else
/// percent-encoded over UTF-8 bytes).
pub fn entrez_quote_plus(text : String) -> String {
  let buf = StringBuilder()
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
  let buf = StringBuilder()
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
  let buf = StringBuilder()
  let mut i = 0
  while i < text.length() {
    if i + from.length() <= text.length() &&
      text[i:i + from.length()].to_owned() == from {
      buf.write_string(to)
      i = i + from.length()
    } else {
      // copy the character verbatim; char_at would render characters
      // outside its whitelist (e.g. '=', '%') as decimal code strings
      buf.write_string(text[i:i + 1].to_owned())
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
  let result : Array[String] = Array(capacity=0)
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
  let result : Array[String] = Array(capacity=0)
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
  { query_key, web_env, }
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
  let xml = client.request(
    entrez_eutils_base_url + "epost.fcgi",
    params,
    post=true,
  )
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
  let buf = StringBuilder()
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
pub fn EntrezItem::new(
  name~ : String,
  type_~ : String,
  value~ : String,
) -> EntrezItem {
  { name, type_, value, }
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
pub fn EntrezDocSum::new(
  id~ : String,
  items~ : Array[EntrezItem],
) -> EntrezDocSum {
  { id, items, }
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
  let result : Array[(String, String, String)] = Array(capacity=0)
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
        let opening = after[0:g + 1].to_owned()
        let value = after[g + 1:e].to_owned()
        let name = entrez_attr(opening, "Name").unwrap_or("")
        let type_ = entrez_attr(opening, "Type").unwrap_or("")
        result.push((name, type_, value))
        rest = after[e + 7:].to_owned()
        go = true
      }
    }
  }
  result
}

///|
/// Parse an ESummary XML response into document summaries.
pub fn parse_esummary(
  xml : String,
) -> Array[EntrezDocSum] raise EntrezUtilError {
  let result : Array[EntrezDocSum] = Array(capacity=0)
  for block in entrez_blocks(xml, "DocSum") {
    let id = match entrez_extract_text(block, "Id") {
      Some(v) => v
      None => entrez_util_fail("DocSum is missing Id")
    }
    let items : Array[EntrezItem] = Array(capacity=0)
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
  { database, query, corrected_query, replaced, }
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
  ESpellResult::new(database~, query~, corrected_query=corrected, replaced~)
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
  { db_to, link_name, ids, }
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
  { db_from, link_set_dbs, }
}

///|
/// Get the source database.
pub fn EntrezLinkSet::db_from(self : EntrezLinkSet) -> String {
  self.db_from
}

///|
/// Get the LinkSetDb blocks of this link set.
pub fn EntrezLinkSet::link_set_dbs(
  self : EntrezLinkSet,
) -> Array[EntrezLinkSetDb] {
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
  let result : Array[EntrezLinkSet] = Array(capacity=0)
  for block in entrez_blocks(xml, "LinkSet") {
    let db_from = entrez_extract_text(block, "DbFrom").unwrap_or("")
    let link_set_dbs : Array[EntrezLinkSetDb] = Array(capacity=0)
    for lsd in entrez_blocks(block, "LinkSetDb") {
      let db_to = entrez_extract_text(lsd, "DbTo").unwrap_or("")
      let link_name = entrez_extract_text(lsd, "LinkName").unwrap_or("")
      let ids : Array[String] = Array(capacity=0)
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
  { journal_title, year, volume, first_page, author_name, key, }
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
  let buf = StringBuilder()
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
pub fn parse_ecitmatch(
  text : String,
) -> Array[(String, String)] raise EntrezUtilError {
  let result : Array[(String, String)] = Array(capacity=0)
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
        let opening = rest[0:gt.unwrap() + 1].to_owned()
        let key = entrez_attr(opening, "key").unwrap_or("")
        let id = entrez_attr(opening, "id").unwrap_or("")
        result.push((key, id))
        remaining = rest[gt.unwrap() + 1:].to_owned()
      }
    }
    None =>
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
  result
}

// ---------------------------------------------------------------------------
// esearch
// ---------------------------------------------------------------------------

///|
/// Result of an ESearch: the hit count, the returned window and the primary
/// IDs. With `usehistory` the QueryKey/WebEnv pair is also present.
pub struct ESearchResult {
  count : String
  retmax : String
  retstart : String
  ids : Array[String]
  translation : String?
  query_key : String?
  web_env : String?
}

///|
/// Create an ESearchResult.
pub fn ESearchResult::new(
  count~ : String,
  retmax~ : String,
  retstart~ : String,
  ids~ : Array[String],
  translation? : String? = None,
  query_key? : String? = None,
  web_env? : String? = None,
) -> ESearchResult {
  { count, retmax, retstart, ids, translation, query_key, web_env, }
}

///|
/// Get the number of hits (the Count field, as a string).
pub fn ESearchResult::count(self : ESearchResult) -> String {
  self.count
}

///|
/// Get the returned window size (RetMax).
pub fn ESearchResult::retmax(self : ESearchResult) -> String {
  self.retmax
}

///|
/// Get the window start (RetStart).
pub fn ESearchResult::retstart(self : ESearchResult) -> String {
  self.retstart
}

///|
/// Get the list of primary IDs (IdList).
pub fn ESearchResult::ids(self : ESearchResult) -> Array[String] {
  self.ids
}

///|
/// Get the term translation, if the server produced one.
pub fn ESearchResult::translation(self : ESearchResult) -> String? {
  self.translation
}

///|
/// Get the QueryKey for use with efetch/esummary (usehistory only).
pub fn ESearchResult::query_key(self : ESearchResult) -> String? {
  self.query_key
}

///|
/// Get the WebEnv for use with efetch/esummary (usehistory only).
pub fn ESearchResult::web_env(self : ESearchResult) -> String? {
  self.web_env
}

///|
/// Parse an esearch.fcgi XML response.
pub fn parse_esearch(xml : String) -> ESearchResult raise EntrezUtilError {
  let count = match entrez_extract_text(xml, "Count") {
    Some(v) => v
    None => entrez_util_fail("esearch response is missing Count")
  }
  let retmax = entrez_extract_text(xml, "RetMax").unwrap_or("")
  let retstart = entrez_extract_text(xml, "RetStart").unwrap_or("")
  let ids = entrez_extract_all(xml, "Id")
  // the real esearch DTD (20060628) spells the term translation tag
  // <QueryTranslation>; <TranslationSet>/<TranslationStack> are siblings
  let translation = entrez_extract_text(xml, "QueryTranslation")
  let query_key = entrez_extract_text(xml, "QueryKey")
  let web_env = entrez_extract_text(xml, "WebEnv")
  ESearchResult::new(
    count~,
    retmax~,
    retstart~,
    ids~,
    translation~,
    query_key~,
    web_env~,
  )
}

///|
/// ESearch: search a database and retrieve the primary IDs (for use in
/// efetch, elink and esummary), term translations, and optionally retain
/// the results in the user's environment (usehistory="y").
pub fn EntrezClient::esearch(
  self : EntrezClient,
  db : String,
  term : String,
  retstart? : String,
  retmax? : String,
  idtype? : String,
  usehistory? : String,
) -> ESearchResult raise EntrezUtilError {
  let params : Array[(String, String)] = [("db", db), ("term", term)]
  match retstart {
    Some(v) => params.push(("retstart", v))
    None => ()
  }
  match retmax {
    Some(v) => params.push(("retmax", v))
    None => ()
  }
  match idtype {
    Some(v) => params.push(("idtype", v))
    None => ()
  }
  match usehistory {
    Some(v) => params.push(("usehistory", v))
    None => ()
  }
  parse_esearch(self.request(entrez_eutils_base_url + "esearch.fcgi", params))
}

// ---------------------------------------------------------------------------
// efetch
// ---------------------------------------------------------------------------

///|
/// EFetch: retrieve records in the requested format from a list of primary
/// IDs (or from the user's environment with a QueryKey/WebEnv pair). The
/// raw response text is returned unchanged; pass rettype/retmode to select
/// the format (e.g. rettype="gb", retmode="text" for GenBank records).
/// As recommended by NCBI, the request automatically uses HTTP POST when
/// there are many identifiers.
pub fn EntrezClient::efetch(
  self : EntrezClient,
  db : String,
  ids : Array[String],
  rettype? : String,
  retmode? : String,
  query_key? : String,
  web_env? : String,
) -> String raise EntrezUtilError {
  let params : Array[(String, String)] = [("db", db), ("id", join_comma(ids))]
  match rettype {
    Some(v) => params.push(("rettype", v))
    None => ()
  }
  match retmode {
    Some(v) => params.push(("retmode", v))
    None => ()
  }
  match query_key {
    Some(v) => params.push(("query_key", v))
    None => ()
  }
  match web_env {
    Some(v) => params.push(("WebEnv", v))
    None => ()
  }
  self.request(entrez_eutils_base_url + "efetch.fcgi", params)
}

// ---------------------------------------------------------------------------
// einfo
// ---------------------------------------------------------------------------

///|
/// Result of an EInfo request. Without a database argument the result
/// carries the list of database names; with a database it carries the
/// database statistics together with the indexed field and link names.
pub struct EInfoResult {
  db_list : Array[String]
  dbname : String
  description : String
  menu : String
  count : String
  last_update : String
  field_names : Array[String]
  link_names : Array[String]
}

///|
/// Create an EInfoResult.
pub fn EInfoResult::new(
  db_list~ : Array[String],
  dbname~ : String,
  description~ : String,
  menu~ : String,
  count~ : String,
  last_update~ : String,
  field_names~ : Array[String],
  link_names~ : Array[String],
) -> EInfoResult {
  {
    db_list,
    dbname,
    description,
    menu,
    count,
    last_update,
    field_names,
    link_names,
  }
}

///|
/// Get the names of all Entrez databases (filled when no db was requested).
pub fn EInfoResult::db_list(self : EInfoResult) -> Array[String] {
  self.db_list
}

///|
/// Get the database name (db-specific requests only).
pub fn EInfoResult::dbname(self : EInfoResult) -> String {
  self.dbname
}

///|
/// Get the database description (db-specific requests only).
pub fn EInfoResult::description(self : EInfoResult) -> String {
  self.description
}

///|
/// Get the database menu name (db-specific requests only).
pub fn EInfoResult::menu(self : EInfoResult) -> String {
  self.menu
}

///|
/// Get the number of records in the database (db-specific requests only).
pub fn EInfoResult::count(self : EInfoResult) -> String {
  self.count
}

///|
/// Get the date of the last update (db-specific requests only).
pub fn EInfoResult::last_update(self : EInfoResult) -> String {
  self.last_update
}

///|
/// Get the names of the indexed fields (db-specific requests only).
pub fn EInfoResult::field_names(self : EInfoResult) -> Array[String] {
  self.field_names
}

///|
/// Get the names of the available links (db-specific requests only).
pub fn EInfoResult::link_names(self : EInfoResult) -> Array[String] {
  self.link_names
}

///|
/// Parse an einfo.fcgi XML response (either form: with or without a
/// database argument).
pub fn parse_einfo(xml : String) -> EInfoResult {
  if entrez_extract_text(xml, "DbList") is Some(_) {
    EInfoResult::new(
      db_list=entrez_extract_all(xml, "DbName"),
      dbname="",
      description="",
      menu="",
      count="",
      last_update="",
      field_names=[],
      link_names=[],
    )
  } else {
    let fields : Array[String] = Array(capacity=0)
    for block in entrez_blocks(xml, "Field") {
      match entrez_extract_text(block, "Name") {
        Some(v) => fields.push(v)
        None => ()
      }
    }
    let links : Array[String] = Array(capacity=0)
    for block in entrez_blocks(xml, "Link") {
      match entrez_extract_text(block, "Name") {
        Some(v) => links.push(v)
        None => ()
      }
    }
    EInfoResult::new(
      db_list=[],
      dbname=entrez_extract_text(xml, "DbName").unwrap_or(""),
      description=entrez_extract_text(xml, "Description").unwrap_or(""),
      // the current einfo DTD (20190110, used by Biopython 1.88 fixtures)
      // spells this tag <MenuName>
      menu=entrez_extract_text(xml, "MenuName").unwrap_or(""),
      count=entrez_extract_text(xml, "Count").unwrap_or(""),
      last_update=entrez_extract_text(xml, "LastUpdate").unwrap_or(""),
      field_names=fields,
      link_names=links,
    )
  }
}

///|
/// EInfo: provide field names, index term counts, last update and
/// available links for each Entrez database. Without `db` the result
/// lists the available databases.
pub fn EntrezClient::einfo(
  self : EntrezClient,
  db? : String,
) -> EInfoResult raise EntrezUtilError {
  let params : Array[(String, String)] = []
  match db {
    Some(v) => params.push(("db", v))
    None => ()
  }
  parse_einfo(self.request(entrez_eutils_base_url + "einfo.fcgi", params))
}

// ---------------------------------------------------------------------------
// egquery
// ---------------------------------------------------------------------------

///|
/// One row of an EGQuery result: the database name and the number of
/// matching records.
pub struct EGQueryRow {
  dbname : String
  count : Int
}

///|
/// Create an EGQueryRow.
pub fn EGQueryRow::new(dbname~ : String, count~ : Int) -> EGQueryRow {
  { dbname, count, }
}

///|
/// Get the database name.
pub fn EGQueryRow::dbname(self : EGQueryRow) -> String {
  self.dbname
}

///|
/// Get the number of matching records in the database.
pub fn EGQueryRow::count(self : EGQueryRow) -> Int {
  self.count
}

///|
/// Result of an EGQuery global search: the term and one row per database.
pub struct EGQueryResult {
  term : String
  rows : Array[EGQueryRow]
}

///|
/// Create an EGQueryResult.
pub fn EGQueryResult::new(
  term~ : String,
  rows~ : Array[EGQueryRow],
) -> EGQueryResult {
  { term, rows, }
}

///|
/// Get the searched term.
pub fn EGQueryResult::term(self : EGQueryResult) -> String {
  self.term
}

///|
/// Get the per-database rows in document order.
pub fn EGQueryResult::rows(self : EGQueryResult) -> Array[EGQueryRow] {
  self.rows
}

///|
/// Look up the count for one database (None when the database is absent).
pub fn EGQueryResult::get(self : EGQueryResult, dbname : String) -> Int? {
  for row in self.rows {
    if row.dbname == dbname {
      return Some(row.count)
    }
  }
  None
}

///|
/// Parse an egquery.fcgi XML response.
pub fn parse_egquery(xml : String) -> EGQueryResult {
  let term = entrez_extract_text(xml, "Term").unwrap_or("")
  let rows : Array[EGQueryRow] = Array(capacity=0)
  for block in entrez_blocks(xml, "Result") {
    let dbname = entrez_extract_text(block, "DbName").unwrap_or("")
    let count_text = entrez_extract_text(block, "Count").unwrap_or("0")
    let count = parse_int_or(count_text, 0)
    rows.push(EGQueryRow::new(dbname~, count~))
  }
  EGQueryResult::new(term~, rows~)
}

///|
fn parse_int_or(text : String, default : Int) -> Int {
  let trimmed = text.trim().to_owned()
  if trimmed.length() == 0 {
    return default
  }
  let mut value = 0
  let mut sign = 1
  let mut idx = 0
  if trimmed.unsafe_get(0).to_int() == 45 {
    // '-'
    sign = -1
    idx = 1
  }
  let mut any = false
  while idx < trimmed.length() {
    let c = trimmed.unsafe_get(idx).to_int()
    if c >= 48 && c <= 57 {
      value = value * 10 + (c - 48)
      any = true
      idx += 1
    } else {
      return default
    }
  }
  if !any {
    default
  } else {
    sign * value
  }
}

///|
/// EGQuery: provide Entrez database counts for a single global search term.
pub fn EntrezClient::egquery(
  self : EntrezClient,
  term : String,
) -> EGQueryResult raise EntrezUtilError {
  let params : Array[(String, String)] = [("term", term)]
  parse_egquery(self.request(entrez_eutils_base_url + "egquery.fcgi", params))
}
