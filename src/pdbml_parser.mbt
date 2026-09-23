///|
/// PDBML (PDB XML) parser.
///
/// Port of Biopython's Bio.PDB.PDBMLParser. Parses PDBML XML into header
/// metadata plus a Model -> Chain -> Residue -> Atom hierarchy, matching the
/// structure grouping of the mmCIF parser (models are 0-indexed in document
/// order; the per-file model serial number is preserved).
///
/// NOTE: Biopython 1.88's `_parse_resolution_from` uses `if element:` to test
/// for presence, which is always False for leaf elements (xml.etree's truth
/// value reflects child count, not presence), so it always returns None. This
/// port checks presence correctly and returns the resolution when found.

///|
pub suberror PdbmlError {
  PdbmlError(String)
}

///|
fn[T] pdbml_fail(msg : String) -> T raise PdbmlError {
  raise PdbmlError(msg)
}

// ---------------------------------------------------------------------------
// Structure types
// ---------------------------------------------------------------------------

///|
pub struct PdbmlAtom {
  name : String
  fullname : String
  coord : Array[Double]
  b_factor : Double
  occupancy : Double
  altloc : String
  serial_number : Int
  element : String
}

///|
pub struct PdbmlResidue {
  component_id : String
  hetero_field : String
  sequence_id : Int
  insertion_code : String
  atoms : Array[PdbmlAtom]
}

///|
pub struct PdbmlChain {
  chain_id : String
  residues : Array[PdbmlResidue]
}

///|
pub struct PdbmlModel {
  serial_number : Int
  chains : Array[PdbmlChain]
}

///|
pub struct PdbmlHeader {
  name : String
  head : String
  idcode : String
  deposition_date : String
  structure_method : String
  resolution : Double?
}

///|
pub struct PdbmlStructure {
  idcode : String
  header : PdbmlHeader
  models : Array[PdbmlModel]
}

// ---------------------------------------------------------------------------
// Minimal namespace-aware XML tokenizer
// ---------------------------------------------------------------------------

///|
struct PdbmlXmlElem {
  name : String
  attributes : Map[String, String]
  text : String
  children : Array[PdbmlXmlElem]
}

///|
fn pdbml_is_ws(cp : Int) -> Bool {
  cp == 32 || cp == 9 || cp == 10 || cp == 13
}

///|
fn pdbml_skip_ws(s : String, i : Int) -> Int {
  let mut j = i
  while j < s.length() && pdbml_is_ws(s.unsafe_get(j).to_int()) {
    j = j + 1
  }
  j
}

///|
fn pdbml_is_name_char(cp : Int) -> Bool {
  (cp >= 65 && cp <= 90) || (cp >= 97 && cp <= 122) || (cp >= 48 && cp <= 57) ||
    cp == '_'.to_int() || cp == '.'.to_int() || cp == '-'.to_int() || cp == ':'.to_int()
}

///|
fn pdbml_read_name(s : String, i : Int) -> (String, Int) {
  let mut j = i
  while j < s.length() && pdbml_is_name_char(s.unsafe_get(j).to_int()) {
    j = j + 1
  }
  (s[i:j].to_string(), j)
}

///|
/// Local name after the namespace prefix (last `:`-separated segment).
fn pdbml_local_name(name : String) -> String {
  let parts = split_by_char(name, ':'.to_int())
  parts[parts.length() - 1]
}

///|
fn pdbml_is_ns_attr(name : String) -> Bool {
  name == "xmlns" || starts_with(name, 0, "xmlns:") || starts_with(name, 0, "xsi:")
}

///|
/// Parse an opening tag `<name attrs...>` (or `/>`). Returns local name,
/// attributes, self-closing flag, and index just past `>`.
fn pdbml_parse_open_tag(
  s : String,
  i : Int,
) -> (String, Map[String, String], Bool, Int) raise PdbmlError {
  let (tname, tk) = pdbml_read_name(s, i + 1)
  let local = pdbml_local_name(tname)
  let attrs : Map[String, String] = Map([], capacity=8)
  let mut j = tk
  while j < s.length() {
    j = pdbml_skip_ws(s, j)
    let cp = s.unsafe_get(j).to_int()
    if cp == '>'.to_int() {
      return (local, attrs, false, j + 1)
    }
    if cp == '/'.to_int() {
      if j + 1 < s.length() && s.unsafe_get(j + 1).to_int() == '>'.to_int() {
        return (local, attrs, true, j + 2)
      }
      pdbml_fail("Malformed self-closing tag")
    }
    let (aname, ak) = pdbml_read_name(s, j)
    j = pdbml_skip_ws(s, ak)
    if s.unsafe_get(j).to_int() != '='.to_int() {
      pdbml_fail("Expected '=' in attribute")
    }
    j = pdbml_skip_ws(s, j + 1)
    let q = s.unsafe_get(j).to_int()
    if q != '"'.to_int() && q != '\''.to_int() {
      pdbml_fail("Expected quoted attribute value")
    }
    j = j + 1
    let val = StringBuilder()
    while j < s.length() && s.unsafe_get(j).to_int() != q {
      val.write_char(s.unsafe_get(j).unsafe_to_char())
      j = j + 1
    }
    j = j + 1
    if !pdbml_is_ns_attr(aname) {
      attrs.set(pdbml_local_name(aname), val.to_string())
    }
  }
  pdbml_fail("Unterminated open tag")
}

///|
fn pdbml_parse_comment(s : String, i : Int) -> Int {
  // s[i..i+4] == "<!--"
  let mut k = i + 4
  while k + 2 < s.length() {
    if s.unsafe_get(k).to_int() == '-'.to_int() &&
    s.unsafe_get(k + 1).to_int() == '-'.to_int() &&
    s.unsafe_get(k + 2).to_int() == '>'.to_int() {
      return k + 3
    }
    k = k + 1
  }
  k
}

///|
/// Parse an element starting at `<`. Returns the element and index past its end.
fn pdbml_parse_element(
  s : String,
  i : Int,
) -> (PdbmlXmlElem, Int) raise PdbmlError {
  let (name, attrs, self_closing, jj) = pdbml_parse_open_tag(s, i)
  let mut j = jj
  if self_closing {
    return (
      PdbmlXmlElem::{ name, attributes: attrs, text: "", children: [] },
      j,
    )
  }
  let text = StringBuilder()
  let children : Array[PdbmlXmlElem] = []
  while j < s.length() {
    let cp = s.unsafe_get(j).to_int()
    if cp == '<'.to_int() {
      if j + 1 < s.length() && s.unsafe_get(j + 1).to_int() == '/'.to_int() {
        let (cname, ck) = pdbml_read_name(s, j + 2)
        let mut k2 = ck
        while k2 < s.length() && s.unsafe_get(k2).to_int() != '>'.to_int() {
          k2 = k2 + 1
        }
        if pdbml_local_name(cname) != name {
          pdbml_fail("Mismatched closing tag")
        }
        return (
          PdbmlXmlElem::{ name, attributes: attrs, text: trim(text.to_string()), children },
          k2 + 1,
        )
      }
      if j + 3 < s.length() &&
      s.unsafe_get(j + 1).to_int() == '!'.to_int() &&
      s.unsafe_get(j + 2).to_int() == '-'.to_int() &&
      s.unsafe_get(j + 3).to_int() == '-'.to_int() {
        let k2 = pdbml_parse_comment(s, j)
        j = k2
        continue
      }
      let (child, nj) = pdbml_parse_element(s, j)
      children.push(child)
      j = nj
      continue
    }
    let mut k2 = j
    while k2 < s.length() && s.unsafe_get(k2).to_int() != '<'.to_int() {
      k2 = k2 + 1
    }
    text.write_string(s[j:k2].to_string())
    j = k2
  }
  pdbml_fail("Unterminated element")
}

///|
fn pdbml_parse_document(xml : String) -> PdbmlXmlElem raise PdbmlError {
  let mut i = 0
  while i < xml.length() {
    i = pdbml_skip_ws(xml, i)
    if i >= xml.length() {
      break
    }
    if xml.unsafe_get(i).to_int() == '<'.to_int() {
      if i + 1 < xml.length() && xml.unsafe_get(i + 1).to_int() == '?'.to_int() {
        let mut k = i + 2
        while k + 1 < xml.length() &&
        !(xml.unsafe_get(k).to_int() == '?'.to_int() &&
        xml.unsafe_get(k + 1).to_int() == '>'.to_int()) {
          k = k + 1
        }
        i = k + 2
        continue
      }
      if i + 3 < xml.length() &&
      xml.unsafe_get(i + 1).to_int() == '!'.to_int() &&
      xml.unsafe_get(i + 2).to_int() == '-'.to_int() &&
      xml.unsafe_get(i + 3).to_int() == '-'.to_int() {
        let k = pdbml_parse_comment(xml, i)
        i = k
        continue
      }
      let (root, _) = pdbml_parse_element(xml, i)
      return root
    }
    i = i + 1
  }
  pdbml_fail("No root element")
}

// ---------------------------------------------------------------------------
// Query helpers
// ---------------------------------------------------------------------------

///|
fn pdbml_find_child(node : PdbmlXmlElem, name : String) -> PdbmlXmlElem? {
  for c in node.children {
    if c.name == name {
      return Some(c)
    }
  }
  None
}

///|
fn pdbml_find_path(node : PdbmlXmlElem, names : Array[String]) -> PdbmlXmlElem? {
  let mut cur = node
  for n in names {
    match pdbml_find_child(cur, n) {
      Some(c) => cur = c
      None => return None
    }
  }
  Some(cur)
}

///|
fn pdbml_child_text(node : PdbmlXmlElem, name : String) -> String {
  match pdbml_find_child(node, name) {
    Some(c) => c.text
    None => ""
  }
}

///|
fn pdbml_path_text(root : PdbmlXmlElem, names : Array[String]) -> String {
  match pdbml_find_path(root, names) {
    Some(e) => e.text
    None => ""
  }
}

///|
fn pdbml_attr_at_path(root : PdbmlXmlElem, names : Array[String], attr : String) -> String {
  match pdbml_find_path(root, names) {
    Some(e) => e.attributes.get(attr).unwrap_or("")
    None => ""
  }
}

// ---------------------------------------------------------------------------
// Field parsing (faithful to Bio.PDB.PDBMLParser)
// ---------------------------------------------------------------------------

///|
fn pdbml_parse_resolution(root : PdbmlXmlElem) -> Double? {
  let candidates : Array[Array[String]] = [
    ["refineCategory", "refine", "ls_d_res_high"],
    ["refine_histCategory", "refine_hist", "d_res_high"],
    ["em_3d_reconstructionCategory", "em_3d_reconstruction", "resolution"],
  ]
  for cand in candidates {
    match pdbml_find_path(root, cand) {
      Some(e) => {
        let t = trim(e.text)
        if t.length() > 0 {
          match parse_double(t) {
            Some(v) => return Some(v)
            None => ()
          }
        }
      }
      None => ()
    }
  }
  None
}

///|
fn pdbml_parse_header(root : PdbmlXmlElem) -> PdbmlHeader {
  PdbmlHeader::{
    name: pdbml_path_text(root, ["structCategory", "struct", "title"]),
    head: pdbml_path_text(root, ["struct_keywordsCategory", "struct_keywords", "text"]),
    idcode: pdbml_attr_at_path(root, ["entryCategory", "entry"], "id"),
    deposition_date: pdbml_path_text(
      root,
      ["pdbx_database_statusCategory", "pdbx_database_status", "recvd_initial_deposition_date"],
    ),
    structure_method: pdbml_attr_at_path(root, ["exptlCategory", "exptl"], "method"),
    resolution: pdbml_parse_resolution(root),
  }
}

///|
fn pdbml_parse_residue_id(elem : PdbmlXmlElem) -> (String, Int, String) {
  let group = pdbml_child_text(elem, "group_PDB")
  let component_id = pdbml_child_text(elem, "label_comp_id")
  let hetero_field = if group == "HETATM" {
    if component_id == "HOH" || component_id == "WAT" {
      "W"
    } else {
      "H"
    }
  } else {
    " "
  }
  let sequence_id = parse_int(pdbml_child_text(elem, "auth_seq_id"))
  let ins = pdbml_child_text(elem, "pdbx_PDB_ins_code")
  let insertion_code = if ins.length() == 0 { " " } else { ins }
  (hetero_field, sequence_id, insertion_code)
}

///|
fn pdbml_parse_atom(elem : PdbmlXmlElem) -> PdbmlAtom {
  let name = pdbml_child_text(elem, "label_atom_id")
  let x = parse_double(pdbml_child_text(elem, "Cartn_x")).unwrap_or(0.0)
  let y = parse_double(pdbml_child_text(elem, "Cartn_y")).unwrap_or(0.0)
  let z = parse_double(pdbml_child_text(elem, "Cartn_z")).unwrap_or(0.0)
  let altloc_raw = pdbml_child_text(elem, "label_alt_id")
  let altloc = if altloc_raw.length() == 0 { " " } else { altloc_raw }
  PdbmlAtom::{
    name,
    fullname: name,
    coord: [x, y, z],
    b_factor: parse_double(pdbml_child_text(elem, "B_iso_or_equiv")).unwrap_or(0.0),
    occupancy: parse_double(pdbml_child_text(elem, "occupancy")).unwrap_or(0.0),
    altloc,
    serial_number: parse_int(elem.attributes.get("id").unwrap_or("0")),
    element: pdbml_child_text(elem, "type_symbol"),
  }
}

// ---------------------------------------------------------------------------
// Grouping + public entry point
// ---------------------------------------------------------------------------

///|
struct PdbmlAtomEntry {
  model_number : Int
  chain_id : String
  component_id : String
  hetero_field : String
  sequence_id : Int
  insertion_code : String
  atom : PdbmlAtom
}

///|
fn pdbml_same_residue(a : PdbmlAtomEntry, b : PdbmlAtomEntry) -> Bool {
  a.hetero_field == b.hetero_field && a.sequence_id == b.sequence_id &&
    a.insertion_code == b.insertion_code && a.component_id == b.component_id
}

///|
fn pdbml_group(entries : Array[PdbmlAtomEntry]) -> Array[PdbmlModel] {
  let models : Array[PdbmlModel] = Array(capacity=0)
  let n = entries.length()
  let mut i = 0
  while i < n {
    let mn = entries[i].model_number
    let chains : Array[PdbmlChain] = Array(capacity=0)
    while i < n && entries[i].model_number == mn {
      let cid = entries[i].chain_id
      let residues : Array[PdbmlResidue] = Array(capacity=0)
      while i < n && entries[i].model_number == mn && entries[i].chain_id == cid {
        let key = entries[i]
        let atoms : Array[PdbmlAtom] = Array(capacity=0)
        while i < n && entries[i].model_number == mn && entries[i].chain_id == cid &&
        pdbml_same_residue(entries[i], key) {
          atoms.push(entries[i].atom)
          i = i + 1
        }
        residues.push(PdbmlResidue::{
          component_id: key.component_id,
          hetero_field: key.hetero_field,
          sequence_id: key.sequence_id,
          insertion_code: key.insertion_code,
          atoms,
        })
      }
      chains.push(PdbmlChain::{ chain_id: cid, residues })
    }
    models.push(PdbmlModel::{ serial_number: mn, chains })
  }
  models
}

///|
/// Parse a PDBML document into a structure (equivalent to
/// `PDBMLParser().get_structure(source)`).
pub fn pdbml_get_structure(xml : String) -> PdbmlStructure raise PdbmlError {
  let root = pdbml_parse_document(xml)
  let header = pdbml_parse_header(root)
  let entries : Array[PdbmlAtomEntry] = Array(capacity=0)
  match pdbml_find_child(root, "atom_siteCategory") {
    Some(category) => {
      for elem in category.children {
        if elem.name == "atom_site" {
          let (hetero, seq, ins) = pdbml_parse_residue_id(elem)
          entries.push(PdbmlAtomEntry::{
            model_number: parse_int(pdbml_child_text(elem, "pdbx_PDB_model_num")),
            chain_id: pdbml_child_text(elem, "auth_asym_id"),
            component_id: pdbml_child_text(elem, "label_comp_id"),
            hetero_field: hetero,
            sequence_id: seq,
            insertion_code: ins,
            atom: pdbml_parse_atom(elem),
          })
        }
      }
    }
    None => ()
  }
  PdbmlStructure::{
    idcode: header.idcode,
    header,
    models: pdbml_group(entries),
  }
}