///|
/// ScanProsite - Scan a protein sequence against the PROSITE motif database
/// through the ExPASy web service.
///
/// Port of Biopython's Bio.ExPASy.ScanProsite module. `scanprosite_url`
/// mirrors ScanProsite.scan()'s URL construction (query parameters encoded
/// with urllib.parse.urlencode quote_plus semantics), `scanprosite_scan`
/// performs the transfer through an injected `fetch` callback, and
/// `scanprosite_read` parses the `<scanprosite_response>` XML document into
/// a `ScanPrositeRecord`. The HTTP transfer is injected so every behaviour
/// can be exercised offline with a canned response.

///|
/// Base URL of the ScanProsite web service.
pub const SCANPROSITE_URL : String = "https://prosite.expasy.org"

///|
/// Error raised by the ScanProsite client for unsupported output formats,
/// transport failures, or malformed XML documents.
pub suberror ScanPrositeError {
  ScanPrositeError(String)
}

///|
fn[T] sp_fail(msg : String) -> T raise ScanPrositeError {
  raise ScanPrositeError(msg)
}

// ---------------------------------------------------------------------------
// Record model
// ---------------------------------------------------------------------------

///|
/// A single motif match reported by ScanProsite.
///
/// The typed fields correspond to the elements Biopython reads from each
/// ``<match>`` block. Any additional element (unknown tag) is preserved in
/// ``extras`` as a string, mirroring how Biopython keeps a plain dict entry.
pub struct ScanPrositeMatch {
  mut sequence_ac : String
  mut sequence_id : String
  mut sequence_db : String
  mut signature_ac : String
  mut start : Int
  mut stop : Int
  mut level : String
  mut level_tag : String
  extras : Map[String, String]
}

///|
pub fn ScanPrositeMatch::new() -> ScanPrositeMatch {
  ScanPrositeMatch::{
    sequence_ac: "",
    sequence_id: "",
    sequence_db: "",
    signature_ac: "",
    start: 0,
    stop: 0,
    level: "",
    level_tag: "",
    extras: Map([], capacity=4),
  }
}

///|
/// The parsed result of a ScanProsite search.
///
/// ``matches`` is the list of motif hits. ``n_match`` and ``n_seq`` are the
/// totals reported by the server. ``capped`` and ``warning`` mirror the two
/// fields on Biopython's ``Record``; Biopython leaves them unset (None) and
/// this port does the same.
pub struct ScanPrositeRecord {
  mut n_match : Int
  mut n_seq : Int
  capped : Int?
  warning : Int?
  matches : Array[ScanPrositeMatch]
}

///|
pub fn ScanPrositeRecord::new() -> ScanPrositeRecord {
  ScanPrositeRecord::{
    n_match: 0,
    n_seq: 0,
    capped: None,
    warning: None,
    matches: [],
  }
}

// ---------------------------------------------------------------------------
// Whitespace / character helpers
// ---------------------------------------------------------------------------

///|
fn sp_is_ws(cp : Int) -> Bool {
  cp == 32 || cp == 9 || cp == 10 || cp == 13 || cp == 11 || cp == 12
}

///|
/// Python str.strip(): strip whitespace from both ends.
fn sp_ws_strip(s : String) -> String {
  let mut start = 0
  while start < s.length() && sp_is_ws(s.unsafe_get(start).to_int()) {
    start = start + 1
  }
  let mut end = s.length()
  while end > start && sp_is_ws(s.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  s[start:end].to_string()
}

///|
/// Split on runs of whitespace, dropping empty pieces (like str.split()).
fn sp_split_ws(s : String) -> Array[String] {
  let out : Array[String] = []
  let mut cur = StringBuilder()
  let mut has = false
  for i = 0; i < s.length(); i = i + 1 {
    if sp_is_ws(s.unsafe_get(i).to_int()) {
      if has {
        out.push(cur.to_string())
        cur = StringBuilder()
        has = false
      }
    } else {
      cur.write_char(s.unsafe_get(i).unsafe_to_char())
      has = true
    }
  }
  if has {
    out.push(cur.to_string())
  }
  out
}

///|
/// Index of the first occurrence of code point `cp` at or after `from`, or
/// -1 when absent.
fn sp_index_of_char_from(s : String, from : Int, cp : Int) -> Int {
  let mut i = from
  while i < s.length() {
    if s.unsafe_get(i).to_int() == cp {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Index of the first occurrence of code point `cp`, or -1 when absent.
fn sp_index_of_char(s : String, cp : Int) -> Int {
  sp_index_of_char_from(s, 0, cp)
}

///|
/// Strip a single pair of surrounding single or double quotes, when both are
/// present.
fn sp_strip_quotes(s : String) -> String {
  let mut start = 0
  let mut end = s.length()
  if end > 0 {
    let c0 = s.unsafe_get(0).to_int()
    if c0 == 34 || c0 == 39 {
      start = 1
    }
  }
  if end > start {
    let cl = s.unsafe_get(end - 1).to_int()
    if cl == 34 || cl == 39 {
      end = end - 1
    }
  }
  s[start:end].to_string()
}

// ---------------------------------------------------------------------------
// URL construction
// ---------------------------------------------------------------------------

///|
/// Build a ScanProsite query URL.
///
/// - seq: query sequence, or UniProtKB (Swiss-Prot / TrEMBL) accession.
/// - output: result format; only "xml" is supported (matching Biopython).
/// - mirror: service mirror (default SCANPROSITE_URL).
/// - params: additional query parameters, passed as (key, value) pairs in
///   the order they should appear (mirrors Biopython's ``**keywords``).
///
/// Raises ScanPrositeError for an unsupported output format.
pub fn scanprosite_url(
  seq~ : String,
  output? : String,
  mirror? : String,
  params? : Array[(String, String)],
) -> String raise ScanPrositeError {
  let out = output.unwrap_or("xml")
  if out != "xml" {
    sp_fail("Only output='xml' is supported currently.")
  }
  let m = mirror.unwrap_or(SCANPROSITE_URL)
  let extra = params.unwrap_or([])
  let buf = StringBuilder()
  buf.write_string("seq=")
  buf.write_string(entrez_quote_plus(seq))
  buf.write_string("&output=")
  buf.write_string(entrez_quote_plus(out))
  for (k, v) in extra {
    buf.write_char('&')
    buf.write_string(entrez_quote_plus(k))
    buf.write_char('=')
    buf.write_string(entrez_quote_plus(v))
  }
  m + "/cgi-bin/prosite/scanprosite/PSScan.cgi?" + buf.to_string()
}

///|
/// Execute a ScanProsite search and return the raw XML response text.
///
/// This mirrors Biopython's ScanProsite.scan(), which returns a handle to
/// the search result. The transfer is provided by the injected `fetch`
/// callback (url, optional POST body) -> response text.
pub fn scanprosite_scan(
  fetch~ : (String, String?) -> String raise ScanPrositeError,
  seq~ : String,
  output? : String,
  mirror? : String,
  params? : Array[(String, String)],
) -> String raise ScanPrositeError {
  let url = scanprosite_url(seq~, output?, mirror?, params?)
  (fetch)(url, None)
}

// ---------------------------------------------------------------------------
// Minimal XML tokenizer / tree, scoped to the ScanProsite response shape
// ---------------------------------------------------------------------------

///|
enum SpXmlToken {
  Open(String, Map[String, String])
  Close(String)
  Text(String)
}

///|
/// Tokenize a well-formed XML document into open / close / text tokens.
/// Whitespace-only text is dropped, which is sufficient for the flat
/// ScanProsite response format.
fn sp_tokenize(xml : String) -> Array[SpXmlToken] raise ScanPrositeError {
  let tokens : Array[SpXmlToken] = []
  let mut i = 0
  let mut text = StringBuilder()
  let mut text_has_nonspace = false
  let n = xml.length()
  while i < n {
    let cp = xml.unsafe_get(i).to_int()
    if cp == 60 {
      // '<' : flush pending text, then read a tag
      if text_has_nonspace {
        tokens.push(Text(text.to_string()))
      }
      text = StringBuilder()
      text_has_nonspace = false
      let j = sp_index_of_char_from(xml, i, 62)
      if j < 0 {
        sp_fail("Unterminated XML tag")
      }
      if i + 1 < n && xml.unsafe_get(i + 1).to_int() == 47 {
        // closing tag </name>
        let name = sp_ws_strip(xml[i + 2:j].to_string())
        tokens.push(Close(name))
        i = j + 1
      } else {
        let (name, attrs) = sp_parse_open_tag(xml[i + 1:j].to_string())
        tokens.push(Open(name, attrs))
        i = j + 1
      }
    } else {
      text.write_char(xml.unsafe_get(i).unsafe_to_char())
      if !sp_is_ws(cp) {
        text_has_nonspace = true
      }
      i = i + 1
    }
  }
  if text_has_nonspace {
    tokens.push(Text(text.to_string()))
  }
  tokens
}

///|
/// Parse the inner content of an open tag: ``name attr="v" ...``.
fn sp_parse_open_tag(inner : String) -> (String, Map[String, String]) {
  let parts = sp_split_ws(inner)
  let name = parts[0]
  let attrs : Map[String, String] = Map([], capacity=4)
  let mut k = 1
  while k < parts.length() {
    let p = parts[k]
    let eq = sp_index_of_char(p, 61)
    if eq >= 0 {
      let key = p[0:eq].to_string()
      let value = sp_strip_quotes(p[eq + 1:p.length()].to_string())
      attrs.set(key, value)
    }
    k = k + 1
  }
  (name, attrs)
}

///|
/// A small DOM node: element name, attributes, child elements, and the text
/// accumulated directly under it (not including children).
struct SpNode {
  name : String
  attrs : Map[String, String]
  children : Array[SpNode]
  mut text : String
}

///|
/// Parse one element (and its subtree) starting at token `pos`.
fn sp_parse_node(
  tokens : Array[SpXmlToken],
  pos : Int,
) -> (SpNode, Int) raise ScanPrositeError {
  match tokens[pos] {
    Open(name, attrs) => {
      let node = SpNode::{ name, attrs, children: [], text: "" }
      let mut i = pos + 1
      let tb = StringBuilder()
      while i < tokens.length() {
        match tokens[i] {
          Close(_) => {
            node.text = tb.to_string()
            return (node, i + 1)
          }
          Text(content) => {
            tb.write_string(content)
            i = i + 1
          }
          Open(_, _) => {
            let (child, next) = sp_parse_node(tokens, i)
            node.children.push(child)
            i = next
          }
        }
      }
      sp_fail("Unexpected end of XML document")
    }
    _ => sp_fail("Expected an XML element")
  }
}

///|
/// Fill a ScanPrositeMatch from a parsed ``<match>`` element.
fn sp_build_match(node : SpNode) -> ScanPrositeMatch {
  let m = ScanPrositeMatch::new()
  for field in node.children {
    let text = sp_ws_strip(field.text)
    let name = field.name
    if name == "sequence_ac" {
      m.sequence_ac = text
    } else if name == "sequence_id" {
      m.sequence_id = text
    } else if name == "sequence_db" {
      m.sequence_db = text
    } else if name == "signature_ac" {
      m.signature_ac = text
    } else if name == "start" {
      m.start = parse_int(text)
    } else if name == "stop" {
      m.stop = parse_int(text)
    } else if name == "level" {
      m.level = text
    } else if name == "level_tag" {
      m.level_tag = text
    } else {
      m.extras.set(field.name, text)
    }
  }
  m
}

// ---------------------------------------------------------------------------
// XML result parsing
// ---------------------------------------------------------------------------

///|
/// Parse a ScanProsite response (XML) into a ScanPrositeRecord.
///
/// Like Biopython's ScanProsite.read(), the document must begin directly
/// with ``<scanprosite_response>``; otherwise the raw text is reported in a
/// ScanPrositeError (Biopython raises ValueError with the same payload).
pub fn scanprosite_read(xml : String) -> ScanPrositeRecord raise ScanPrositeError {
  if !starts_with(xml, 0, "<scanprosite_response>") {
    sp_fail(xml)
  }
  let tokens = sp_tokenize(xml)
  if tokens.length() == 0 {
    sp_fail("Empty XML document")
  }
  let (root, _) = sp_parse_node(tokens, 0)
  let record = ScanPrositeRecord::new()
  for child in root.children {
    if child.name == "matchset" {
      record.n_match = parse_int(child.attrs.get("n_match").unwrap_or("0"))
      record.n_seq = parse_int(child.attrs.get("n_seq").unwrap_or("0"))
      for match_node in child.children {
        if match_node.name == "match" {
          record.matches.push(sp_build_match(match_node))
        }
      }
    }
  }
  record
}

///|
/// A sample ScanProsite XML response for tests and demos.
pub fn scanprosite_sample_xml() -> String {
  "<scanprosite_response>\n" +
  "<matchset n_match=\"3\" n_seq=\"1\">\n" +
  "<match>\n" +
  "<sequence_ac>P00750</sequence_ac>\n" +
  "<sequence_id>TPA_HUMAN</sequence_id>\n" +
  "<sequence_db>swissprot</sequence_db>\n" +
  "<start>114</start>\n" +
  "<stop>135</stop>\n" +
  "<signature_ac>PS00001</signature_ac>\n" +
  "<level>0</level>\n" +
  "<level_tag>(0)</level_tag>\n" +
  "</match>\n" +
  "<match>\n" +
  "<sequence_ac>P00750</sequence_ac>\n" +
  "<sequence_id>TPA_HUMAN</sequence_id>\n" +
  "<sequence_db>swissprot</sequence_db>\n" +
  "<start>147</start>\n" +
  "<stop>172</stop>\n" +
  "<signature_ac>PS00134</signature_ac>\n" +
  "<level>0</level>\n" +
  "<level_tag>(0)</level_tag>\n" +
  "</match>\n" +
  "<match>\n" +
  "<sequence_ac>P00750</sequence_ac>\n" +
  "<sequence_id>TPA_HUMAN</sequence_id>\n" +
  "<sequence_db>swissprot</sequence_db>\n" +
  "<start>123</start>\n" +
  "<stop>456</stop>\n" +
  "<signature_ac>PS00002</signature_ac>\n" +
  "<level>0</level>\n" +
  "<level_tag>(0)</level_tag>\n" +
  "</match>\n" +
  "</matchset>\n" +
  "</scanprosite_response>"
}