///|
/// Prodoc - parse records from the prosite.doc documentation file of PROSITE.
///
/// Port of Biopython's Bio.ExPASy.Prodoc module. A Prodoc record carries an
/// accession, a list of PROSITE cross-references, free-format text, and a
/// list of literature / electronic references. Both the single-record reader
/// and the multi-record iterator are provided, mirroring Biopython's
/// ``read`` and ``parse``.

///|
/// Error raised when a Prodoc document is malformed (bad accession line,
/// unexpected reference line, or more than one record where one is expected).
pub suberror ProdocError {
  ProdocError(String)
}

///|
fn[T] pd_fail(msg : String) -> T raise ProdocError {
  raise ProdocError(msg)
}

///|
/// A citation inside a Prodoc record.
///
/// ``number`` is the reference number (or "E1", "E2", ... for electronic
/// references). ``authors`` holds the author list (empty for electronic
/// references) and ``citation`` holds the bibliographic text or URL.
pub struct ProdocReference {
  mut number : String
  mut authors : String
  mut citation : String
}

///|
pub fn ProdocReference::new() -> ProdocReference {
  ProdocReference::{ number: "", authors: "", citation: "" }
}

///|
/// A single Prodoc (prosite.doc) record.
///
/// ``prosite_refs`` is a list of (accession, name) pairs referencing the
/// matching PROSITE entries.
pub struct ProdocRecord {
  mut accession : String
  prosite_refs : Array[(String, String)]
  mut text : String
  references : Array[ProdocReference]
}

///|
pub fn ProdocRecord::new() -> ProdocRecord {
  ProdocRecord::{
    accession: "",
    prosite_refs: [],
    text: "",
    references: [],
  }
}

// ---------------------------------------------------------------------------
// Text helpers
// ---------------------------------------------------------------------------

///|
fn pd_is_ws(cp : Int) -> Bool {
  cp == 32 || cp == 9 || cp == 10 || cp == 13 || cp == 11 || cp == 12
}

///|
/// Python str.strip(): strip whitespace from both ends.
fn pd_strip(s : String) -> String {
  let mut start = 0
  while start < s.length() && pd_is_ws(s.unsafe_get(start).to_int()) {
    start = start + 1
  }
  let mut end = s.length()
  while end > start && pd_is_ws(s.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  s[start:end].to_owned()
}

///|
/// Python str.rstrip(): strip trailing whitespace.
fn pd_rstrip(s : String) -> String {
  let mut end = s.length()
  while end > 0 && pd_is_ws(s.unsafe_get(end - 1).to_int()) {
    end = end - 1
  }
  s[0:end].to_owned()
}

///|
fn pd_ends_with(s : String, suffix : String) -> Bool {
  if suffix.length() > s.length() {
    return false
  }
  starts_with(s, s.length() - suffix.length(), suffix)
}

///|
/// Split a document into lines, mirroring Python's file iteration: each line
/// loses its trailing newline, a \r\n pair is treated as a single newline,
/// and a final trailing newline does not produce an extra empty line.
fn pd_lines(content : String) -> Array[String] {
  let raw = split_by_char(content, '\n'.to_int())
  let lines : Array[String] = []
  let n = raw.length()
  for i = 0; i < n; i = i + 1 {
    let mut line = raw[i]
    if line.length() > 0 && line.unsafe_get(line.length() - 1).to_int() == 13 {
      line = line[0:line.length() - 1].to_owned()
    }
    if i == n - 1 && line.length() == 0 {
      // trailing empty piece produced by a final newline: drop it
    } else {
      lines.push(line)
    }
  }
  lines
}

///|
/// Split on the two-character sequence "; " (semicolon + space), like
/// Python's line.split("; ").
fn pd_split_semicolon_space(s : String) -> Array[String] {
  let out : Array[String] = []
  let mut cur = StringBuilder()
  let mut i = 0
  while i < s.length() {
    if i + 1 < s.length() &&
      s.unsafe_get(i).to_int() == 59 &&
      s.unsafe_get(i + 1).to_int() == 32 {
      out.push(cur.to_string())
      cur = StringBuilder()
      i = i + 2
    } else {
      cur.write_char(s.unsafe_get(i).unsafe_to_char())
      i = i + 1
    }
  }
  out.push(cur.to_string())
  out
}

// ---------------------------------------------------------------------------
// Reference line handling
// ---------------------------------------------------------------------------

///|
/// Parse a "[NN] authors..." or "[EN] url..." reference start line.
fn pd_read_reference_start(record : ProdocRecord, line : String) -> Unit {
  let reference = ProdocReference::new()
  reference.number = pd_strip(line[1:3].to_owned())
  if starts_with(line, 1, "E") {
    reference.citation = pd_strip(line[4:line.length()].to_owned())
  } else {
    reference.authors = pd_strip(line[4:line.length()].to_owned())
  }
  record.references.push(reference)
}

///|
/// Handle a continuation line of the current reference. The line is known to
/// be non-blank. Indented lines (five leading spaces) append to the authors
/// when the author list is still continuing (ends with a comma), otherwise
/// to the citation.
fn pd_read_reference_line(
  record : ProdocRecord,
  line : String,
) -> Unit raise ProdocError {
  if starts_with(line, 0, "     ") {
    let last = record.references.length() - 1
    if record.references[last].authors.length() > 0 &&
      pd_ends_with(record.references[last].authors, ",") {
      record.references[last].authors = record.references[last].authors +
        pd_rstrip(line[4:line.length()].to_owned())
    } else {
      record.references[last].citation = record.references[last].citation +
        line[5:line.length()].to_owned()
    }
  } else {
    pd_fail("I don't understand the reference line\n\{line}")
  }
}

///|
/// Read the next record starting at line index `i` (which may fall on blank
/// lines or "//" separators). Returns None when only blank / separator lines
/// remain; otherwise the record and the index just past the record.
fn pd_next_record(
  lines : Array[String],
  start : Int,
) -> (ProdocRecord?, Int) raise ProdocError {
  let n = lines.length()
  let mut i = start
  // Skip blank lines and "//" separators between records.
  while i < n {
    let line = pd_rstrip(lines[i])
    if line.length() > 0 && !starts_with(line, 0, "//") {
      break
    }
    i = i + 1
  }
  if i >= n {
    return (None, i)
  }

  let record = ProdocRecord::new()

  // Accession line: "{PDOCnnnnn}".
  let acc_line = pd_rstrip(lines[i])
  if !starts_with(acc_line, 0, "{PDOC") {
    pd_fail("Line does not start with '{PDOC':\n\{acc_line}")
  }
  if !pd_ends_with(acc_line, "}") {
    pd_fail("I don't understand accession line\n\{acc_line}")
  }
  record.accession = acc_line[1:acc_line.length() - 1].to_owned()
  i = i + 1

  // PROSITE cross-reference lines: "{PSnnnnn; NAME}".
  while i < n {
    if starts_with(lines[i], 0, "{PS") {
      let rl = pd_rstrip(lines[i])
      if !pd_ends_with(rl, "}") {
        pd_fail("I don't understand the Prosite reference on line\n\{rl}")
      }
      let inner = rl[1:rl.length() - 1].to_owned()
      let parts = pd_split_semicolon_space(inner)
      if parts.length() != 2 {
        pd_fail("I don't understand the Prosite reference on line\n\{rl}")
      }
      record.prosite_refs.push((parts[0], parts[1]))
      i = i + 1
    } else {
      break
    }
  }
  if i >= n {
    pd_fail("Unexpected end of stream.")
  }

  // The text body must begin with "{BEGIN}".
  if !starts_with(lines[i], 0, "{BEGIN") {
    pd_fail("Line does not start with '{BEGIN':\n\{pd_rstrip(lines[i])}")
  }
  i = i + 1

  // Read text / references / copyright until "{END}".
  // mode: 1 = text, 2 = reference, 3 = copyright, 0 = none (skip).
  let mut mode = 1
  let text_buf = StringBuilder()
  while i < n {
    let line = lines[i]
    if starts_with(line, 0, "{END}") {
      for j = 0; j < record.references.length(); j = j + 1 {
        record.references[j].citation = pd_rstrip(record.references[j].citation)
        record.references[j].authors = pd_rstrip(record.references[j].authors)
      }
      record.text = text_buf.to_string()
      return (Some(record), i + 1)
    } else if starts_with(line, 0, "[") &&
      starts_with(line, 3, "]") &&
      starts_with(line, 4, " ") {
      pd_read_reference_start(record, line)
      mode = 2
      i = i + 1
    } else if starts_with(line, 0, "+----") {
      mode = 3
      i = i + 1
    } else if mode == 1 {
      text_buf.write_string(line)
      text_buf.write_char('\n')
      i = i + 1
    } else if mode == 2 {
      if pd_strip(line).length() == 0 {
        mode = 0
      } else {
        pd_read_reference_line(record, line)
      }
      i = i + 1
    } else if mode == 3 {
      if starts_with(line, 0, "+----") {
        mode = 0
      }
      i = i + 1
    } else {
      // mode == 0: skip the line (matches Biopython's None read_line)
      i = i + 1
    }
  }
  pd_fail("Unexpected end of stream.")
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

///|
/// Parse a document containing one or more Prodoc records.
///
/// This mirrors Biopython's Prodoc.parse() generator: blank lines and "//"
/// separators are skipped, and each record is returned in order.
pub fn prodoc_parse(content : String) -> Array[ProdocRecord] raise ProdocError {
  let lines = pd_lines(content)
  let out : Array[ProdocRecord] = []
  let mut i = 0
  while i < lines.length() {
    let (record, next) = pd_next_record(lines, i)
    i = next
    match record {
      Some(r) => out.push(r)
      None => break
    }
  }
  out
}

///|
/// Read a document containing exactly one Prodoc record.
///
/// Mirrors Biopython's Prodoc.read(): returns None for an empty document,
/// and raises ProdocError when anything follows the first record.
pub fn prodoc_read(content : String) -> ProdocRecord? raise ProdocError {
  let lines = pd_lines(content)
  let (record, next) = pd_next_record(lines, 0)
  match record {
    None => None
    Some(_) =>
      if next < lines.length() {
        pd_fail("More than one Prodoc record found")
      } else {
        record
      }
  }
}

///|
/// A sample two-record prosite.doc document for tests and demos.
pub fn prodoc_sample_content() -> String {
  "{PDOC00001}\n" +
  "{PS00001; ASN_GLYCOSYLATION}\n" +
  "{PS00007; PROTEIN_KINASE_ATP}\n" +
  "{BEGIN}\n" +
  "N-glycosylation site documentation.\n" +
  "Second paragraph line.\n" +
  "[ 1] Hunter, T., Plowman, G.D.\n" +
  "     In \"Methods in Enzymology\", Academic Press, New York, (1991).\n" +
  "[E1] https://www.ebi.ac.uk/interpro/entry/PDOC00001/\n" +
  "+----------------------+\n" +
  "| Copyright notice     |\n" +
  "+----------------------+\n" +
  "{END}\n" +
  "//\n" +
  "{PDOC00002}\n" +
  "{PS00002; GLYCOSAMINOGLYCAN}\n" +
  "{BEGIN}\n" +
  "Glycosaminoglycan attachment site.\n" +
  "[ 1] Bourdon, M.A.,\n" +
  "     Krusius, T.,\n" +
  "     Ruoslahti, E.\n" +
  "     Identification and synthesis of a recognition signal.\n" +
  "\n" +
  "{END}\n" +
  "//"
}