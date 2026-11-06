///|
/// Additional structure I/O format utilities.
///
/// Biotite supports multiple file formats. This module provides
/// converters and utilities for CIF, mmCIF, and PDBx formats.

///|
/// Convert PDB-style record to simplified mmCIF category string.
///
/// This is a lightweight converter for common structure records.
pub fn pdb_to_mmcif_header(pdb_id : String, title : String) -> String {
  let mut result = "data_" + pdb_id + "\n"
  result = result + "#\n"
  result = result + "_entry.id " + pdb_id + "\n"
  result = result + "_struct.title '" + title + "'\n"
  result = result + "#\n"
  result
}

///|
/// Generate mmCIF atom_site header.
pub fn mmcif_atom_site_header() -> String {
  let mut result = "loop_\n"
  result = result + "_atom_site.group_PDB\n"
  result = result + "_atom_site.id\n"
  result = result + "_atom_site.type_symbol\n"
  result = result + "_atom_site.label_atom_id\n"
  result = result + "_atom_site.label_alt_id\n"
  result = result + "_atom_site.label_comp_id\n"
  result = result + "_atom_site.label_asym_id\n"
  result = result + "_atom_site.label_seq_id\n"
  result = result + "_atom_site.Cartn_x\n"
  result = result + "_atom_site.Cartn_y\n"
  result = result + "_atom_site.Cartn_z\n"
  result = result + "_atom_site.occupancy\n"
  result = result + "_atom_site.B_iso_or_equiv\n"
  result = result + "_atom_site.pdbx_formal_charge\n"
  result = result + "_atom_site.auth_seq_id\n"
  result = result + "_atom_site.auth_comp_id\n"
  result = result + "_atom_site.auth_asym_id\n"
  result = result + "_atom_site.auth_atom_id\n"
  result
}

///|
/// Format a coordinate for CIF output (6 decimal places).
fn fmt_coord(val : Double) -> String {
  // Simple formatting without printf
  let int_part = val.floor().to_int()
  let frac = val - int_part.to_double()
  let frac_abs = if frac < 0.0 { -frac } else { frac }
  let frac_scaled = (frac_abs * 1000000.0 + 0.5).floor().to_int()
  let mut frac_str = frac_scaled.to_string()
  while frac_str.length() < 6 {
    frac_str = "0" + frac_str
  }
  let sign = if val < 0.0 && int_part == 0 { "-" } else { "" }
  sign + int_part.to_string() + "." + frac_str
}

///|
/// Format occupancy/B-factor for CIF (2 decimal places).
fn fmt_2dec(val : Double) -> String {
  let int_part = val.floor().to_int()
  let frac = val - int_part.to_double()
  let frac_abs = if frac < 0.0 { -frac } else { frac }
  let frac_scaled = (frac_abs * 100.0 + 0.5).floor().to_int()
  let mut frac_str = frac_scaled.to_string()
  while frac_str.length() < 2 {
    frac_str = "0" + frac_str
  }
  int_part.to_string() + "." + frac_str
}

///|
/// Generate an mmCIF atom_site record line.
pub fn mmcif_atom_site_record(
  group : String,
  id : Int32,
  element : String,
  atom_name : String,
  alt_id : String,
  comp_id : String,
  asym_id : String,
  seq_id : Int32,
  x : Double,
  y : Double,
  z : Double,
  occupancy : Double,
  b_factor : Double,
  auth_seq_id : Int32,
) -> String {
  let alt = if alt_id.length() == 0 { "." } else { alt_id }
  group +
  " " +
  id.to_string() +
  " " +
  element +
  " " +
  atom_name +
  " " +
  alt +
  " " +
  comp_id +
  " " +
  asym_id +
  " " +
  seq_id.to_string() +
  " " +
  fmt_coord(x) +
  " " +
  fmt_coord(y) +
  " " +
  fmt_coord(z) +
  " " +
  fmt_2dec(occupancy) +
  " " +
  fmt_2dec(b_factor) +
  " ? " +
  auth_seq_id.to_string() +
  " " +
  comp_id +
  " " +
  asym_id +
  " " +
  atom_name +
  "\n"
}

///|
/// Generate a complete mmCIF file from atom records.
///
/// Each atom is a tuple of (element, atom_name, res_name, chain_id,
/// res_seq, x, y, z, occupancy, b_factor).
pub fn atoms_to_mmcif(
  pdb_id : String,
  title : String,
  atoms : Array[
    (
      String,
      String,
      String,
      String,
      Int32,
      Double,
      Double,
      Double,
      Double,
      Double,
    ),
  ],
) -> String {
  let mut result = pdb_to_mmcif_header(pdb_id, title)
  result = result + mmcif_atom_site_header()
  for i in 0..<atoms.length() {
    let (elem, atom_name, res_name, chain, res_seq, x, y, z, occ, bfac) = atoms[i]
    result = result +
      mmcif_atom_site_record(
        "ATOM",
        i + 1,
        elem,
        atom_name,
        "",
        res_name,
        chain,
        res_seq,
        x,
        y,
        z,
        occ,
        bfac,
        res_seq,
      )
  }
  result = result + "#\n"
  result
}

///|
/// Parse a simple CIF key-value pair.
///
/// Returns the value for a given key from a CIF text block.
pub fn cif_get_value(cif_text : String, key : String) -> String {
  let prefix = "_" + key + " "
  let n = cif_text.length()
  let prefix_len = prefix.length()
  for i in 0..<=(n - prefix_len) {
    let mut is_match = true
    for j in 0..<prefix_len {
      if cif_text.unsafe_get(i + j) != prefix.unsafe_get(j) {
        is_match = false
        break
      }
    }
    if is_match {
      // Extract value after prefix
      let mut start = i + prefix_len
      while start < n && cif_text.unsafe_get(start) == 32 {
        start = start + 1
      }
      let mut end = start
      while end < n &&
            cif_text.unsafe_get(end) != 10 &&
            cif_text.unsafe_get(end) != 13 {
        end = end + 1
      }
      let mut value = ""
      for k in start..<end {
        value = value + cif_text.unsafe_get(k).unsafe_to_char().to_string()
      }
      return value
    }
  }
  ""
}
