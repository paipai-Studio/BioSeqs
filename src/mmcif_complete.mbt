///|
/// mmCIF (PDBx) I/O — complete bidirectional support for AtomArray.
///
/// This module provides:
///   1. Full-featured mmCIF **writer** that emits the standard
///      categories needed to round-trip an AtomArray faithfully:
///        _entry / _entity / _pdbx_poly_seq_scheme / _atom_site / _struct
///      with proper mmCIF quoting (double-quoted values with
///      embedded quotes doubled).
///
///   2. Direct AtomArray ↔ mmCIF convenience entry points so users
///      never have to touch the intermediate Structure layer.

// ============================================================================
// mmCIF writer helpers
// ============================================================================

fn mmcif_escape_field(value : String) -> String {
  if value.length() == 0 {
    return ".".to_string()
  }
  let first = value.unsafe_get(0).unsafe_to_char()
  let starts_digit = first >= '0' && first <= '9'
  let starts_special = first == '.' ||
    first == '-' ||
    first == '+' ||
    first == '?'
  let mut needs_quote = starts_digit || starts_special
  if !needs_quote {
    let mut j = 0
    while j < value.length() {
      let c = value.unsafe_get(j).unsafe_to_char()
      if c == '"' || c == '#' || c == ';' || c == '_' || c == ' ' || c == '\t' {
        needs_quote = true
        break
      }
      j = j + 1
    }
  }
  if needs_quote {
    let mut escaped = ""
    let mut j = 0
    while j < value.length() {
      let c = value.unsafe_get(j).unsafe_to_char()
      if c == '"' {
        escaped = escaped + "\"\""
      } else {
        escaped = escaped + c.to_string()
      }
      j = j + 1
    }
    return "\"" + escaped + "\""
  }
  value
}

///|
fn mmcif_write_loop(
  sb : StringBuilder,
  category : String,
  headers : Array[String],
  rows : Array[Array[String]],
) -> Unit {
  sb.write_string("loop_\n")
  for h in headers {
    sb.write_string("_")
    sb.write_string(category)
    sb.write_string(".")
    sb.write_string(h)
    sb.write_string("\n")
  }
  for row in rows {
    let mut first = true
    for v in row {
      if !first {
        sb.write_string(" ")
      }
      sb.write_string(mmcif_escape_field(v))
      first = false
    }
    sb.write_string("\n")
  }
  ()
}

///|
fn mmcif_write_single(
  sb : StringBuilder,
  category : String,
  attr : String,
  value : String,
) -> Unit {
  sb.write_string("_")
  sb.write_string(category)
  sb.write_string(".")
  sb.write_string(attr)
  sb.write_string("   ")
  sb.write_string(mmcif_escape_field(value))
  sb.write_string("\n")
  ()
}

///|
/// Format a Double to 6 decimals — mmCIF convention.
fn mmcif_format_double(d : Double) -> String {
  let s = d.to_string()
  // Ensure we always have a decimal point with 6 digits
  let dot_pos = {
    let mut j = 0
    let mut found = -1
    while j < s.length() {
      if s.unsafe_get(j).unsafe_to_char() == '.' {
        found = j
        break
      }
      j = j + 1
    }
    found
  }
  if dot_pos >= 0 {
    let int_part = substring(s, 0, dot_pos)
    let mut frac = substring(s, dot_pos + 1, s.length() - dot_pos - 1)
    while frac.length() < 6 {
      frac = frac + "0"
    }
    if frac.length() > 6 {
      frac = substring(frac, 0, 6)
    }
    return int_part + "." + frac
  } else {
    return s + ".000000"
  }
}

///|
/// Given a chain char, return its 1-based entity id string.
fn entity_id_for(chain_order : Array[Char], chain : Char) -> String {
  let mut j = 0
  while j < chain_order.length() {
    if chain_order[j] == chain {
      return (j + 1).to_string()
    }
    j = j + 1
  }
  "1".to_string()
}

// ============================================================================
// Top-level writer — AtomArray → complete mmCIF text
// ============================================================================

///|
/// Full-featured mmCIF writer for an AtomArray.
///
/// Emits the standard categories:
///   _entry.id
///   _entity  (one row per unique chain — all marked polymer)
///   _pdbx_poly_seq_scheme  (residue sequence per chain)
///   _atom_site  (complete atom records)
///   _struct.title
pub fn AtomArray::to_mmcif_full(self : AtomArray) -> String {
  let sb = StringBuilder()
  let n = self.len()
  let chain_order : Array[Char] = Array(capacity=8)

  // --- Discover unique chains in order ---
  if n > 0 {
    let mut prev_chain = self.chain_ids[0]
    chain_order.push(prev_chain)
    let mut i = 1
    while i < n {
      let ch = self.chain_ids[i]
      if ch != prev_chain {
        chain_order.push(ch)
        prev_chain = ch
      }
      i = i + 1
    }
  }

  // --- data block header ---
  sb.write_string("data_BIOSEQ\n#\n")

  // --- _entry ---
  mmcif_write_single(sb, "entry", "id", "BIOSEQ".to_string())
  sb.write_string("#\n")

  // --- _struct ---
  mmcif_write_single(sb, "struct", "title", "BioSeqs generated structure")
  sb.write_string("#\n")

  // --- _entity (one per chain) ---
  let entity_headers : Array[String] = [
    "id", "type", "src_method", "pdbx_description",
  ]
  let entity_rows : Array[Array[String]] = Array(capacity=chain_order.length())
  let mut ei = 0
  while ei < chain_order.length() {
    let eid : String = (ei + 1).to_string()
    let desc : String = "Chain ".to_string() + chain_order[ei].to_string()
    entity_rows.push([eid, "polymer".to_string(), "man".to_string(), desc])
    ei = ei + 1
  }
  mmcif_write_loop(sb, "entity", entity_headers, entity_rows)
  sb.write_string("#\n")

  // --- _pdbx_poly_seq_scheme (residue sequence per chain) ---
  let scheme_headers : Array[String] = [
    "asym_id", "entity_id", "seq_id", "mon_id", "pdbx_PDB_ins_code",
  ]
  let scheme_rows : Array[Array[String]] = Array(capacity=n)
  let res_starts = self.get_residue_starts()
  let n_res = res_starts.length()
  let mut ri = 0
  while ri < n_res {
    let rs = res_starts[ri]
    let chain = self.chain_ids[rs]
    let rn = self.res_names[rs]
    let rid = self.res_ids[rs]
    let entity_id = entity_id_for(chain_order, chain)
    scheme_rows.push([
      chain.to_string(),
      entity_id,
      rid.to_string(),
      rn,
      ".".to_string(), // ins_code
    ])
    ri = ri + 1
  }
  mmcif_write_loop(sb, "pdbx_poly_seq_scheme", scheme_headers, scheme_rows)
  sb.write_string("#\n")

  // --- _atom_site ---
  let atom_headers : Array[String] = [
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
    "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x",
    "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "pdbx_PDB_model_num",
  ]
  let atom_rows : Array[Array[String]] = Array(capacity=n)
  let mut serial : Int32 = 1
  let mut mi = 0
  while mi < n {
    let entity_id = entity_id_for(chain_order, self.chain_ids[mi])
    let alt = if self.altlocs[mi] == ' ' {
      ".".to_string()
    } else {
      self.altlocs[mi].to_string()
    }
    let group : String = if is_atom(self.names[mi]) { "ATOM" } else { "HETATM" }
    atom_rows.push([
      group,
      serial.to_string(),
      self.elements[mi],
      self.names[mi],
      alt,
      self.res_names[mi],
      self.chain_ids[mi].to_string(),
      entity_id,
      self.res_ids[mi].to_string(),
      ".".to_string(), // ins_code
      mmcif_format_double(self.coords[mi].x),
      mmcif_format_double(self.coords[mi].y),
      mmcif_format_double(self.coords[mi].z),
      mmcif_format_double(self.occupancies[mi]),
      mmcif_format_double(self.b_factors[mi]),
      "1".to_string(),
    ])
    serial = serial + 1
    mi = mi + 1
  }
  mmcif_write_loop(sb, "atom_site", atom_headers, atom_rows)
  sb.write_string("#\n")

  sb.write_string("#\n")
  sb.to_string()
}

///|
fn is_atom(name : String) -> Bool {
  // Accept standard atom names as ATOM; rest → HETATM
  name == "N" ||
  name == "CA" ||
  name == "C" ||
  name == "O" ||
  name == "CB" ||
  name == "CG" ||
  name == "CG1" ||
  name == "CG2" ||
  name == "CD" ||
  name == "CD1" ||
  name == "CD2" ||
  name == "CE" ||
  name == "CE1" ||
  name == "CE2" ||
  name == "CE3" ||
  name == "CZ" ||
  name == "CZ2" ||
  name == "CZ3" ||
  name == "CH2" ||
  name == "OG" ||
  name == "OG1" ||
  name == "OD1" ||
  name == "OD2" ||
  name == "OE1" ||
  name == "OE2" ||
  name == "ND1" ||
  name == "ND2" ||
  name == "NE" ||
  name == "NE1" ||
  name == "NH1" ||
  name == "NH2" ||
  name == "NZ" ||
  name == "SG" ||
  name == "SD" ||
  name == "FE" ||
  name == "OH" ||
  name == "OXT"
}

// ============================================================================
// AtomArray ↔ mmCIF convenience entry points
// ============================================================================

///|
/// Convert this AtomArray to mmCIF (PDBx/mmCIF) text.
/// Writes the full standard set of categories.
pub fn AtomArray::to_mmcif(self : AtomArray) -> String {
  self.to_mmcif_full()
}

///|
/// Parse mmCIF (PDBx/mmCIF) text into an AtomArray.
/// Picks the first data block if the input contains multiple.
pub fn AtomArray::from_mmcif(text : String) -> AtomArray {
  let block = parse_mmcif(text)
  let atoms = block.to_atoms()
  AtomArray::from_atoms(atoms)
}
