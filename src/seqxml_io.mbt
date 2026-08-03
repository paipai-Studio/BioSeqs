// SeqXML format parser and writer, mirroring BioPython's
// Bio.SeqIO.SeqXmlIO.
//
// SeqXML is a lightweight XML format for sequence data exchange.
// Structure:
//   <?xml version="1.0" encoding="UTF-8"?>
//   <seqXML>
//     <species name="Homo sapiens" ncbiTaxID="9606"/>
//     <sourceDB name="UniProt"/>
//     <entry id="seq1" desc="Test sequence 1">
//       <dna>ATGCATGC</dna>
//     </entry>
//     <entry id="seq2" desc="Test protein">
//       <protein>MKLVGV</protein>
//     </entry>
//   </seqXML>
//
// Each <entry> has a required id and optional desc attribute, and contains
// exactly one sequence element (<dna>, <rna>, or <protein>) plus any number
// of <property name="..." value="..."/> elements.

// ============================================================================
// Data Structures
// ============================================================================

///|
/// The molecular type of a SeqXML sequence element.
pub enum SeqXmlType {
  DnaType
  RnaType
  ProteinType
  UnknownType
} derive(Eq)

///|
/// A single SeqXML entry: an identifier, an optional description, a
/// sequence of a known molecular type, and a list of name/value properties.
pub struct SeqXmlEntry {
  id : String
  description : String
  seq_type : SeqXmlType
  sequence : String
  properties : Array[(String, String)]
}

///|
/// A complete SeqXML document: a list of entries together with optional
/// species and source database metadata.
pub struct SeqXmlDocument {
  entries : Array[SeqXmlEntry]
  mut species_name : String
  mut species_taxid : String
  mut source_db : String
}

// ============================================================================
// SeqXmlType methods
// ============================================================================

///|
/// Map a sequence element tag name to a SeqXmlType.
/// "dna" -> DnaType, "rna" -> RnaType, "protein" -> ProteinType.
pub fn SeqXmlType::from_tag(tag : String) -> SeqXmlType {
  match tag {
    "dna" => DnaType
    "rna" => RnaType
    "protein" => ProteinType
    _ => UnknownType
  }
}

///|
/// Return the XML tag name for this sequence type.
pub fn SeqXmlType::to_tag(self : SeqXmlType) -> String {
  match self {
    DnaType => "dna"
    RnaType => "rna"
    ProteinType => "protein"
    UnknownType => "unknown"
  }
}

///|
/// Return a human-readable description of this sequence type.
pub fn SeqXmlType::description(self : SeqXmlType) -> String {
  match self {
    DnaType => "DNA"
    RnaType => "RNA"
    ProteinType => "Protein"
    UnknownType => "Unknown"
  }
}

// ============================================================================
// SeqXmlEntry methods
// ============================================================================

///|
/// Construct a SeqXmlEntry with the given id, sequence type, sequence, and
/// optional description.
pub fn SeqXmlEntry::new(
  id : String,
  seq_type : SeqXmlType,
  sequence : String,
  description? : String = "",
) -> SeqXmlEntry {
  SeqXmlEntry::{ id, description, seq_type, sequence, properties: [] }
}

///|
/// Return the entry identifier.
pub fn SeqXmlEntry::id(self : SeqXmlEntry) -> String {
  self.id
}

///|
/// Return the entry description.
pub fn SeqXmlEntry::description(self : SeqXmlEntry) -> String {
  self.description
}

///|
/// Return the molecular type of the entry.
pub fn SeqXmlEntry::seq_type(self : SeqXmlEntry) -> SeqXmlType {
  self.seq_type
}

///|
/// Return the sequence string of the entry.
pub fn SeqXmlEntry::sequence(self : SeqXmlEntry) -> String {
  self.sequence
}

///|
/// Return the list of (name, value) properties attached to the entry.
pub fn SeqXmlEntry::properties(self : SeqXmlEntry) -> Array[(String, String)] {
  self.properties
}

///|
/// Append a (name, value) property to this entry.
pub fn SeqXmlEntry::add_property(
  self : SeqXmlEntry,
  name : String,
  value : String,
) -> Unit {
  self.properties.push((name, value))
}

///|
/// Return a concise string representation of this entry.
pub fn SeqXmlEntry::to_string(self : SeqXmlEntry) -> String {
  let type_desc = SeqXmlType::description(self.seq_type)
  "SeqXmlEntry(id=\{self.id}, type=\{type_desc}, length=\{self.sequence.length()}, description=\{self.description}, properties=\{self.properties.length()})"
}

// ============================================================================
// SeqXmlDocument methods
// ============================================================================

///|
/// Construct an empty SeqXmlDocument.
pub fn SeqXmlDocument::new() -> SeqXmlDocument {
  SeqXmlDocument::{
    entries: [],
    species_name: "",
    species_taxid: "",
    source_db: "",
  }
}

///|
/// Return the list of entries in this document.
pub fn SeqXmlDocument::entries(self : SeqXmlDocument) -> Array[SeqXmlEntry] {
  self.entries
}

///|
/// Return the species name recorded in the document header.
pub fn SeqXmlDocument::species_name(self : SeqXmlDocument) -> String {
  self.species_name
}

///|
/// Return the NCBI taxonomy ID recorded in the document header.
pub fn SeqXmlDocument::species_taxid(self : SeqXmlDocument) -> String {
  self.species_taxid
}

///|
/// Return the source database name recorded in the document header.
pub fn SeqXmlDocument::source_db(self : SeqXmlDocument) -> String {
  self.source_db
}

///|
/// Return the number of entries in this document.
pub fn SeqXmlDocument::n_entries(self : SeqXmlDocument) -> Int {
  self.entries.length()
}

///|
/// Set the species metadata for this document.
pub fn SeqXmlDocument::set_species(
  self : SeqXmlDocument,
  name : String,
  taxid : String,
) -> Unit {
  self.species_name = name
  self.species_taxid = taxid
}

///|
/// Set the source database metadata for this document.
pub fn SeqXmlDocument::set_source_db(
  self : SeqXmlDocument,
  name : String,
) -> Unit {
  self.source_db = name
}

///|
/// Append an entry to this document.
pub fn SeqXmlDocument::add_entry(
  self : SeqXmlDocument,
  entry : SeqXmlEntry,
) -> Unit {
  self.entries.push(entry)
}

// ============================================================================
// Internal XML helpers
// ============================================================================

///|
/// Find the first occurrence of `sub` in `s` at or after position `from`.
/// Returns the index of the match, or None if not found.
fn seqxml_index_of_from(s : String, sub : String, from : Int) -> Int? {
  let n = s.length()
  let m = sub.length()
  if m == 0 || from >= n {
    return None
  }
  let start = if from < 0 { 0 } else { from }
  if start + m > n {
    return None
  }
  let mut i = start
  while i <= n - m {
    let mut j = 0
    let mut matched = true
    while j < m {
      if s.unsafe_get(i + j) != sub.unsafe_get(j) {
        matched = false
        break
      }
      j = j + 1
    }
    if matched {
      return Some(i)
    }
    i = i + 1
  }
  None
}

///|
/// Check whether a character (as UInt16) is an XML tag-name delimiter:
/// whitespace, ">", or "/".
fn seqxml_is_tag_delim(c : UInt16) -> Bool {
  let sp = ' '.to_int().to_uint16()
  let tab = '\t'.to_int().to_uint16()
  let nl = '\n'.to_int().to_uint16()
  let cr = '\r'.to_int().to_uint16()
  let gt = '>'.to_int().to_uint16()
  let sl = '/'.to_int().to_uint16()
  c == sp || c == tab || c == nl || c == cr || c == gt || c == sl
}

///|
/// Find the next opening tag named `tag_name` at or after position `from`.
/// Returns the position of the "<" or None. A tag boundary (whitespace, ">",
/// or "/") must follow the tag name to avoid prefix matches.
fn seqxml_find_tag(
  content : String,
  tag_name : String,
  from : Int,
) -> Int? {
  let pattern = "<" + tag_name
  let plen = pattern.length()
  let n = content.length()
  let mut search = from
  while search < n {
    match seqxml_index_of_from(content, pattern, search) {
      None => return None
      Some(pos) => {
        let after = pos + plen
        if after >= n {
          return Some(pos)
        }
        if seqxml_is_tag_delim(content.unsafe_get(after)) {
          return Some(pos)
        }
        search = pos + 1
      }
    }
  }
  None
}

///|
/// Given the position of a "<" that starts a tag, find the matching ">"
/// that closes the opening tag, respecting quoted attribute values.
/// Returns the position of the ">".
fn seqxml_find_tag_end(s : String, from : Int) -> Int? {
  let n = s.length()
  let dq = '"'.to_int().to_uint16()
  let sq = '\''.to_int().to_uint16()
  let gt = '>'.to_int().to_uint16()
  let mut i = from + 1
  let mut in_dq = false
  let mut in_sq = false
  while i < n {
    let c = s.unsafe_get(i)
    if in_dq {
      if c == dq {
        in_dq = false
      }
    } else if in_sq {
      if c == sq {
        in_sq = false
      }
    } else if c == dq {
      in_dq = true
    } else if c == sq {
      in_sq = true
    } else if c == gt {
      return Some(i)
    }
    i = i + 1
  }
  None
}

///|
/// Extract the value of an XML attribute from a tag body string.
/// Handles both double- and single-quoted values, and enforces a word
/// boundary before the attribute name. Returns the raw (still-escaped) value.
fn seqxml_extract_attr_value(tag : String, attr : String) -> String? {
  let n = tag.length()
  // Try double-quoted form: attr="..."
  let dq_pattern = attr + "=\""
  let dq_len = dq_pattern.length()
  let mut search = 0
  while search < n {
    match seqxml_index_of_from(tag, dq_pattern, search) {
      None => break
      Some(pos) => {
        let boundary_ok = if pos == 0 {
          true
        } else {
          let prev = tag.unsafe_get(pos - 1).to_int()
          prev == ' '.to_int() ||
          prev == '\t'.to_int() ||
          prev == '\n'.to_int() ||
          prev == '\r'.to_int()
        }
        if boundary_ok {
          let val_start = pos + dq_len
          match seqxml_index_of_from(tag, "\"", val_start) {
            Some(end) => return Some(tag[val_start:end].to_owned())
            None => return None
          }
        }
        search = pos + 1
      }
    }
  }
  // Try single-quoted form: attr='...'
  let sq_pattern = attr + "='"
  let sq_len = sq_pattern.length()
  let mut search2 = 0
  while search2 < n {
    match seqxml_index_of_from(tag, sq_pattern, search2) {
      None => break
      Some(pos) => {
        let boundary_ok = if pos == 0 {
          true
        } else {
          let prev = tag.unsafe_get(pos - 1).to_int()
          prev == ' '.to_int() ||
          prev == '\t'.to_int() ||
          prev == '\n'.to_int() ||
          prev == '\r'.to_int()
        }
        if boundary_ok {
          let val_start = pos + sq_len
          match seqxml_index_of_from(tag, "'", val_start) {
            Some(end) => return Some(tag[val_start:end].to_owned())
            None => return None
          }
        }
        search2 = pos + 1
      }
    }
  }
  None
}

///|
/// Extract an attribute value and XML-unescape it, returning `default`
/// when the attribute is absent.
fn seqxml_attr_or(
  tag : String,
  attr : String,
  default : String,
) -> String {
  match seqxml_extract_attr_value(tag, attr) {
    Some(v) => seqxml_unescape(v)
    None => default
  }
}

///|
/// Escape a string for use inside an XML attribute value.
fn seqxml_escape_attr(s : String) -> String {
  s.replace(old="&", new="&amp;")
   .replace(old="<", new="&lt;")
   .replace(old=">", new="&gt;")
   .replace(old="\"", new="&quot;")
   .replace(old="'", new="&apos;")
}

///|
/// Escape a string for use as XML text content.
fn seqxml_escape_text(s : String) -> String {
  s.replace(old="&", new="&amp;")
   .replace(old="<", new="&lt;")
   .replace(old=">", new="&gt;")
}

///|
/// Decode XML entities (&amp; &lt; &gt; &quot; &apos;) back to characters.
/// &amp; is decoded last to avoid double-decoding.
fn seqxml_unescape(s : String) -> String {
  s.replace(old="&lt;", new="<")
   .replace(old="&gt;", new=">")
   .replace(old="&quot;", new="\"")
   .replace(old="&apos;", new="'")
   .replace(old="&amp;", new="&")
}

///|
/// Trim whitespace (spaces, tabs, newlines, carriage returns) from both
/// ends of a string.
fn seqxml_trim_ws(s : String) -> String {
  let n = s.length()
  let sp = ' '.to_int()
  let tab = '\t'.to_int()
  let nl = '\n'.to_int()
  let cr = '\r'.to_int()
  let mut start = 0
  while start < n {
    let c = s.unsafe_get(start).to_int()
    if c == sp || c == tab || c == nl || c == cr {
      start = start + 1
    } else {
      break
    }
  }
  let mut end = n
  while end > start {
    let c = s.unsafe_get(end - 1).to_int()
    if c == sp || c == tab || c == nl || c == cr {
      end = end - 1
    } else {
      break
    }
  }
  s[start:end].to_owned()
}

// ============================================================================
// Internal parsing helpers
// ============================================================================

///|
/// Parse the <species> and <sourceDB> header elements and set the
/// corresponding fields on the document.
fn seqxml_parse_header(content : String, doc : SeqXmlDocument) -> Unit {
  match seqxml_find_tag(content, "species", 0) {
    Some(start) =>
      match seqxml_find_tag_end(content, start) {
        Some(gt_pos) => {
          let body = content[start + 1 : gt_pos].to_owned()
          let clean = if body.has_suffix("/") {
            body[0 : body.length() - 1].to_owned()
          } else {
            body
          }
          let name = seqxml_attr_or(clean, "name", "")
          let taxid = seqxml_attr_or(clean, "ncbiTaxID", "")
          doc.set_species(name, taxid)
        }
        None => ()
      }
    None => ()
  }
  match seqxml_find_tag(content, "sourceDB", 0) {
    Some(start) =>
      match seqxml_find_tag_end(content, start) {
        Some(gt_pos) => {
          let body = content[start + 1 : gt_pos].to_owned()
          let clean = if body.has_suffix("/") {
            body[0 : body.length() - 1].to_owned()
          } else {
            body
          }
          let name = seqxml_attr_or(clean, "name", "")
          doc.set_source_db(name)
        }
        None => ()
      }
    None => ()
  }
}

///|
/// Extract the trimmed text content of the first element named `tag` found
/// within `inner`. Returns None if no such element exists.
fn seqxml_extract_element_text(inner : String, tag : String) -> String? {
  match seqxml_find_tag(inner, tag, 0) {
    Some(open_start) =>
      match seqxml_find_tag_end(inner, open_start) {
        Some(gt_pos) => {
          let content_start = gt_pos + 1
          let close_tag = "</" + tag + ">"
          match seqxml_index_of_from(inner, close_tag, content_start) {
            Some(close_pos) => {
              let text = inner[content_start : close_pos].to_owned()
              Some(seqxml_unescape(seqxml_trim_ws(text)))
            }
            None => None
          }
        }
        None => None
      }
    None => None
  }
}

///|
/// Parse the sequence element from an entry's inner content, returning
/// the molecular type and the (unescaped) sequence string.
fn seqxml_parse_sequence(inner : String) -> (SeqXmlType, String) {
  match seqxml_extract_element_text(inner, "dna") {
    Some(text) => return (DnaType, text)
    None => ()
  }
  match seqxml_extract_element_text(inner, "rna") {
    Some(text) => return (RnaType, text)
    None => ()
  }
  match seqxml_extract_element_text(inner, "protein") {
    Some(text) => return (ProteinType, text)
    None => ()
  }
  (UnknownType, "")
}

///|
/// Parse all <property name="..." value="..."/> elements within an entry's
/// inner content.
fn seqxml_parse_properties(inner : String) -> Array[(String, String)] {
  let result : Array[(String, String)] = []
  let n = inner.length()
  let mut pos = 0
  while pos < n {
    match seqxml_find_tag(inner, "property", pos) {
      None => break
      Some(prop_start) =>
        match seqxml_find_tag_end(inner, prop_start) {
          Some(gt_pos) => {
            let body = inner[prop_start + 1 : gt_pos].to_owned()
            let clean = if body.has_suffix("/") {
              body[0 : body.length() - 1].to_owned()
            } else {
              body
            }
            let name = seqxml_attr_or(clean, "name", "")
            let value = seqxml_attr_or(clean, "value", "")
            result.push((name, value))
            pos = gt_pos + 1
          }
          None => break
        }
    }
  }
  result
}

// ============================================================================
// Public parsing functions
// ============================================================================

///|
/// Parse a full SeqXML document from a string and return a SeqXmlDocument.
///
/// Handles the XML declaration, root element, species/sourceDB header
/// elements, and all <entry> elements with their sequence and property
/// children.
pub fn seqxml_parse(content : String) -> SeqXmlDocument {
  let doc = SeqXmlDocument::new()
  seqxml_parse_header(content, doc)
  let close_tag = "</entry>"
  let close_len = close_tag.length()
  let n = content.length()
  let mut pos = 0
  while pos < n {
    match seqxml_find_tag(content, "entry", pos) {
      None => break
      Some(entry_start) =>
        match seqxml_find_tag_end(content, entry_start) {
          None => break
          Some(gt_pos) => {
            let body = content[entry_start + 1 : gt_pos].to_owned()
            let self_closing = body.has_suffix("/")
            let clean_body = if self_closing {
              body[0 : body.length() - 1].to_owned()
            } else {
              body
            }
            let id = seqxml_attr_or(clean_body, "id", "")
            let desc = seqxml_attr_or(clean_body, "desc", "")
            let inner_start = gt_pos + 1
            let close_pos_opt : Int? = if self_closing {
              None
            } else {
              seqxml_index_of_from(content, close_tag, inner_start)
            }
            let (seq_type, sequence, properties) = match close_pos_opt {
              Some(close_pos) => {
                let inner = content[inner_start : close_pos].to_owned()
                let (st, seq) = seqxml_parse_sequence(inner)
                let props = seqxml_parse_properties(inner)
                (st, seq, props)
              }
              None => (UnknownType, "", [])
            }
            let entry = SeqXmlEntry::new(
              id,
              seq_type,
              sequence,
              description=desc,
            )
            for p in properties {
              entry.add_property(p.0, p.1)
            }
            doc.add_entry(entry)
            pos = match close_pos_opt {
              Some(close_pos) => close_pos + close_len
              None => if self_closing { gt_pos + 1 } else { n }
            }
          }
        }
    }
  }
  doc
}

///|
/// Extract the id attribute from a single <entry ...> tag line.
/// Returns the unescaped id, or None if the line is not an entry tag or
/// has no id attribute.
pub fn seqxml_parse_entry_id(line : String) -> String? {
  match seqxml_find_tag(line, "entry", 0) {
    Some(start) =>
      match seqxml_find_tag_end(line, start) {
        Some(gt_pos) => {
          let body = line[start + 1 : gt_pos].to_owned()
          let clean = if body.has_suffix("/") {
            body[0 : body.length() - 1].to_owned()
          } else {
            body
          }
          match seqxml_extract_attr_value(clean, "id") {
            Some(v) => Some(seqxml_unescape(v))
            None => None
          }
        }
        None => None
      }
    None => None
  }
}

// ============================================================================
// Public writing functions
// ============================================================================

///|
/// Convert a single SeqXmlEntry to its XML representation, indented as it
/// would appear inside a <seqXML> document (2-space base indent).
pub fn seqxml_entry_to_xml(entry : SeqXmlEntry) -> String {
  let buf = StringBuilder::new()
  // Opening <entry> tag with id (and optional desc).
  buf.write_string("  <entry id=\"")
  buf.write_string(seqxml_escape_attr(entry.id))
  buf.write_string("\"")
  if entry.description.length() > 0 {
    buf.write_string(" desc=\"")
    buf.write_string(seqxml_escape_attr(entry.description))
    buf.write_string("\"")
  }
  buf.write_string(">\n")
  // Property elements (4-space indent).
  for p in entry.properties {
    buf.write_string("    <property name=\"")
    buf.write_string(seqxml_escape_attr(p.0))
    buf.write_string("\" value=\"")
    buf.write_string(seqxml_escape_attr(p.1))
    buf.write_string("\"/>\n")
  }
  // Sequence element (4-space indent).
  let tag = SeqXmlType::to_tag(entry.seq_type)
  buf.write_string("    <")
  buf.write_string(tag)
  buf.write_string(">")
  buf.write_string(seqxml_escape_text(entry.sequence))
  buf.write_string("</")
  buf.write_string(tag)
  buf.write_string(">\n")
  // Closing </entry> tag.
  buf.write_string("  </entry>\n")
  buf.to_string()
}

///|
/// Convert a full SeqXmlDocument to its SeqXML string representation,
/// including the XML declaration, header metadata, and all entries.
pub fn seqxml_to_xml(doc : SeqXmlDocument) -> String {
  let buf = StringBuilder::new()
  buf.write_string("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
  buf.write_string("<seqXML>\n")
  // Species header (if either name or taxid is set).
  if doc.species_name.length() > 0 || doc.species_taxid.length() > 0 {
    buf.write_string("  <species")
    if doc.species_name.length() > 0 {
      buf.write_string(" name=\"")
      buf.write_string(seqxml_escape_attr(doc.species_name))
      buf.write_string("\"")
    }
    if doc.species_taxid.length() > 0 {
      buf.write_string(" ncbiTaxID=\"")
      buf.write_string(seqxml_escape_attr(doc.species_taxid))
      buf.write_string("\"")
    }
    buf.write_string("/>\n")
  }
  // Source database header (if set).
  if doc.source_db.length() > 0 {
    buf.write_string("  <sourceDB name=\"")
    buf.write_string(seqxml_escape_attr(doc.source_db))
    buf.write_string("\"/>\n")
  }
  // Entries.
  for entry in doc.entries {
    buf.write_string(seqxml_entry_to_xml(entry))
  }
  buf.write_string("</seqXML>\n")
  buf.to_string()
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Auto-detect the molecular type of a sequence string.
/// If every character is in {A,C,G,T,N} (case-insensitive) it is DNA;
/// otherwise, if every character is in {A,C,G,U,N} it is RNA;
/// otherwise it is treated as protein.
fn seqxml_detect_type(seq_str : String) -> SeqXmlType {
  if seq_str.length() == 0 {
    return UnknownType
  }
  let upper = seq_str.to_upper()
  let n = upper.length()
  let a = 'A'.to_int().to_uint16()
  let c = 'C'.to_int().to_uint16()
  let g = 'G'.to_int().to_uint16()
  let t = 'T'.to_int().to_uint16()
  let u = 'U'.to_int().to_uint16()
  let nn = 'N'.to_int().to_uint16()
  let mut all_dna = true
  let mut all_rna = true
  for i in 0..<n {
    let ch = upper.unsafe_get(i)
    if !(ch == a || ch == c || ch == g || ch == t || ch == nn) {
      all_dna = false
    }
    if !(ch == a || ch == c || ch == g || ch == u || ch == nn) {
      all_rna = false
    }
  }
  if all_dna {
    DnaType
  } else if all_rna {
    RnaType
  } else {
    ProteinType
  }
}

///|
/// Convert a SeqXmlDocument into an array of SeqRecord objects.
/// Each entry's id becomes the SeqRecord id and name; the description is
/// carried over, and the sequence is built via Seq::new.
pub fn seqxml_to_seq_records(doc : SeqXmlDocument) -> Array[SeqRecord] {
  let records : Array[SeqRecord] = []
  for entry in doc.entries {
    let seq = Seq::new(entry.sequence)
    let rec = SeqRecord::new(
      seq,
      id=entry.id,
      name=entry.id,
      description=entry.description,
    )
    records.push(rec)
  }
  records
}

///|
/// Build a SeqXmlDocument from an array of SeqRecord objects.
/// When `seq_type` is UnknownType (the default), the type of each record is
/// auto-detected from its sequence content.
pub fn seqxml_from_seq_records(
  records : Array[SeqRecord],
  seq_type? : SeqXmlType = UnknownType,
) -> SeqXmlDocument {
  let doc = SeqXmlDocument::new()
  for rec in records {
    let seq_str = rec.seq.to_string()
    let st = if seq_type == UnknownType {
      seqxml_detect_type(seq_str)
    } else {
      seq_type
    }
    let entry = SeqXmlEntry::new(
      rec.id,
      st,
      seq_str,
      description=rec.description,
    )
    doc.add_entry(entry)
  }
  doc
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Return a sample SeqXmlDocument with three entries: one DNA, one RNA, and
/// one protein sequence. The protein entry also carries a sample property.
pub fn seqxml_sample_document() -> SeqXmlDocument {
  let doc = SeqXmlDocument::new()
  doc.set_species("Homo sapiens", "9606")
  doc.set_source_db("BioSeqs")
  let e1 = SeqXmlEntry::new(
    "sample_dna",
    DnaType,
    "ATGCGTACGTTGC",
    description="Sample DNA sequence",
  )
  let e2 = SeqXmlEntry::new(
    "sample_rna",
    RnaType,
    "AUGCGUACGUUGC",
    description="Sample RNA sequence",
  )
  let e3 = SeqXmlEntry::new(
    "sample_protein",
    ProteinType,
    "MKLVGV",
    description="Sample protein sequence",
  )
  e3.add_property("feature", "domain")
  doc.add_entry(e1)
  doc.add_entry(e2)
  doc.add_entry(e3)
  doc
}
