///|
/// mmCIF (Macromolecular Crystallographic Information File) format writer.
/// Inspired by Biopython Bio.PDB.MMCIFIO.
///
/// This module writes protein Structure objects to mmCIF format strings.
/// It complements the existing mmCIF parser in mmcif.mbt.

///|
/// Format a coordinate value with 3 decimal places.
pub fn mmcif_format_coord(x : Double) -> String {
  mmcifio_format_double(x, 3)
}

///|
/// Format a B-factor or occupancy value with 2 decimal places.
pub fn mmcif_format_bfactor(b : Double) -> String {
  mmcifio_format_double(b, 2)
}

///|
/// Format a double with the specified precision (no width padding).
/// Truncates (does not round) the fractional part to the given precision.
fn mmcifio_format_double(value : Double, precision : Int) -> String {
  let s = value.to_string()
  let parts = s.split(".").to_array()
  let int_part = if parts.length() > 0 { parts[0].to_owned() } else { "0" }
  let raw_frac = if parts.length() > 1 { parts[1].to_owned() } else { "" }
  let mut padded_frac = raw_frac
  if padded_frac.length() < precision {
    padded_frac = padded_frac + "0".repeat(precision - padded_frac.length())
  } else if padded_frac.length() > precision {
    padded_frac = padded_frac[0:precision].to_owned()
  }
  int_part + "." + padded_frac
}

///|
/// Escape special characters in mmCIF values.
///
/// - Empty values become `?` (unknown).
/// - Values containing whitespace, quotes, or brackets are wrapped in
///   single quotes; embedded single quotes are doubled.
/// - Other values are returned as-is.
pub fn mmcif_escape_value(value : String) -> String {
  if value.length() == 0 {
    return "?"
  }
  if mmcifio_value_needs_quote(value) {
    let escaped = mmcifio_replace_all(value, "'", "''")
    "'" + escaped + "'"
  } else {
    value
  }
}

///|
/// Check whether an mmCIF value needs to be quoted.
fn mmcifio_value_needs_quote(value : String) -> Bool {
  let n = value.length()
  let mut i = 0
  while i < n {
    let c = value.unsafe_get(i).unsafe_to_char()
    if c == ' ' ||
      c == '\t' ||
      c == '\n' ||
      c == '\r' ||
      c == '\'' ||
      c == '"' ||
      c == '[' ||
      c == ']' ||
      c == '{' ||
      c == '}' {
      return true
    }
    i = i + 1
  }
  false
}

///|
/// Replace all occurrences of `from` with `to` in a string.
fn mmcifio_replace_all(s : String, from : String, to : String) -> String {
  let from_len = from.length()
  if from_len == 0 {
    return s
  }
  let buf = StringBuilder::new()
  let n = s.length()
  let mut i = 0
  while i < n {
    if i + from_len <= n && s[i:i + from_len].to_owned() == from {
      buf.write_string(to)
      i = i + from_len
    } else {
      buf.write_char(s.unsafe_get(i).unsafe_to_char())
      i = i + 1
    }
  }
  buf.to_string()
}

///|
/// Format a single atom line for the `_atom_site` loop using the
/// `ATOM` group_PDB value.
pub fn mmcif_format_atom(atom : Atom, atom_id : Int) -> String {
  mmcif_format_atom_with_group(atom, atom_id, "ATOM")
}

///|
/// Format a single atom line with an explicit group_PDB value
/// (`ATOM` or `HETATM`).
pub fn mmcif_format_atom_with_group(
  atom : Atom,
  atom_id : Int,
  group_pdb : String,
) -> String {
  let type_symbol = if atom.element.length() > 0 { atom.element } else { "?" }
  let label_atom_id = mmcif_escape_value(atom.name)
  let label_alt_id = if atom.altloc == ' ' { "." } else { "\{atom.altloc}" }
  let label_comp_id = mmcif_escape_value(atom.resname)
  let label_asym_id = mmcif_escape_value(atom.chainid.to_string())
  let label_entity_id = "1"
  let label_seq_id = atom.resseq.to_string()
  let pdbx_PDB_ins_code = if atom.icode == ' ' { "?" } else { "\{atom.icode}" }
  let cartn_x = mmcif_format_coord(atom.coord.x)
  let cartn_y = mmcif_format_coord(atom.coord.y)
  let cartn_z = mmcif_format_coord(atom.coord.z)
  let occupancy = mmcif_format_bfactor(atom.occupancy)
  let b_iso_or_equiv = mmcif_format_bfactor(atom.bfactor)
  let pdbx_formal_charge = if atom.charge == 0 {
    "?"
  } else if atom.charge > 0 {
    "\{atom.charge}+"
  } else {
    "\{-atom.charge}-"
  }
  let auth_seq_id = atom.resseq.to_string()
  let auth_comp_id = mmcif_escape_value(atom.resname)
  let auth_asym_id = mmcif_escape_value(atom.chainid.to_string())
  let auth_atom_id = mmcif_escape_value(atom.name)
  group_pdb +
  "   " +
  atom_id.to_string() +
  "   " +
  type_symbol +
  "   " +
  label_atom_id +
  "   " +
  label_alt_id +
  "   " +
  label_comp_id +
  "   " +
  label_asym_id +
  "   " +
  label_entity_id +
  "   " +
  label_seq_id +
  "   " +
  pdbx_PDB_ins_code +
  "   " +
  cartn_x +
  "   " +
  cartn_y +
  "   " +
  cartn_z +
  "   " +
  occupancy +
  "   " +
  b_iso_or_equiv +
  "   " +
  pdbx_formal_charge +
  "   " +
  auth_seq_id +
  "   " +
  auth_comp_id +
  "   " +
  auth_asym_id +
  "   " +
  auth_atom_id
}

///|
/// Write the mmCIF header categories (`entry` and `entity`) for a Structure.
pub fn write_mmcif_header(structure : Structure) -> String {
  let buf = StringBuilder::new()
  let sid = if structure.id.length() > 0 { structure.id } else { "XXXX" }
  buf.write_string("_entry.id   ")
  buf.write_string(sid)
  buf.write_char('\n')
  buf.write_string("#\n")
  buf.write_string("_entity.id   1\n")
  buf.write_string("_entity.type   polymer\n")
  buf.write_string("_entity.src_method   man\n")
  buf.write_string("#\n")
  buf.to_string()
}

///|
/// Write the `_atom_site` loop category for a Structure.
pub fn write_mmcif_atom_site(structure : Structure) -> String {
  let buf = StringBuilder::new()
  buf.write_string("loop_\n")
  buf.write_string("_atom_site.group_PDB\n")
  buf.write_string("_atom_site.id\n")
  buf.write_string("_atom_site.type_symbol\n")
  buf.write_string("_atom_site.label_atom_id\n")
  buf.write_string("_atom_site.label_alt_id\n")
  buf.write_string("_atom_site.label_comp_id\n")
  buf.write_string("_atom_site.label_asym_id\n")
  buf.write_string("_atom_site.label_entity_id\n")
  buf.write_string("_atom_site.label_seq_id\n")
  buf.write_string("_atom_site.pdbx_PDB_ins_code\n")
  buf.write_string("_atom_site.Cartn_x\n")
  buf.write_string("_atom_site.Cartn_y\n")
  buf.write_string("_atom_site.Cartn_z\n")
  buf.write_string("_atom_site.occupancy\n")
  buf.write_string("_atom_site.B_iso_or_equiv\n")
  buf.write_string("_atom_site.pdbx_formal_charge\n")
  buf.write_string("_atom_site.auth_seq_id\n")
  buf.write_string("_atom_site.auth_comp_id\n")
  buf.write_string("_atom_site.auth_asym_id\n")
  buf.write_string("_atom_site.auth_atom_id\n")
  let mut atom_id = 1
  for model in structure.models {
    for chain in model.chains {
      for residue in chain.residues {
        let group = if residue.is_het() { "HETATM" } else { "ATOM" }
        for atom in residue.atoms {
          buf.write_string(mmcif_format_atom_with_group(atom, atom_id, group))
          buf.write_char('\n')
          atom_id = atom_id + 1
        }
      }
    }
  }
  buf.write_string("#\n")
  buf.to_string()
}

///|
/// Write a Structure to an mmCIF format string.
pub fn write_mmcif(structure : Structure) -> String {
  let buf = StringBuilder::new()
  let sid = if structure.id.length() > 0 { structure.id } else { "XXXX" }
  buf.write_string("data_")
  buf.write_string(sid)
  buf.write_char('\n')
  buf.write_string("#\n")
  buf.write_string(write_mmcif_header(structure))
  buf.write_string(write_mmcif_atom_site(structure))
  buf.to_string()
}

///|
/// Create a sample Structure for testing and demos.
///
/// The structure contains a single model with one chain of two residues
/// (ALA and GLY), for a total of 8 atoms.
pub fn mmcif_sample_structure() -> Structure {
  let ala_n = Atom::new(
    name="N",
    coord=Vector3::new(1.0, 2.0, 3.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    occupancy=1.0,
    bfactor=20.0,
    element="N",
  )
  let ala_ca = Atom::new(
    name="CA",
    coord=Vector3::new(1.5, 2.5, 3.5),
    resname="ALA",
    chainid='A',
    resseq=1,
    occupancy=1.0,
    bfactor=20.0,
    element="C",
  )
  let ala_c = Atom::new(
    name="C",
    coord=Vector3::new(2.0, 2.5, 4.0),
    resname="ALA",
    chainid='A',
    resseq=1,
    occupancy=1.0,
    bfactor=20.0,
    element="C",
  )
  let ala_o = Atom::new(
    name="O",
    coord=Vector3::new(2.5, 2.5, 4.5),
    resname="ALA",
    chainid='A',
    resseq=1,
    occupancy=1.0,
    bfactor=20.0,
    element="O",
  )
  let ala = Residue::new(
    resname="ALA",
    chainid='A',
    resseq=1,
    atoms=[ala_n, ala_ca, ala_c, ala_o],
  )
  let gly_n = Atom::new(
    name="N",
    coord=Vector3::new(3.0, 2.5, 5.0),
    resname="GLY",
    chainid='A',
    resseq=2,
    occupancy=1.0,
    bfactor=20.0,
    element="N",
  )
  let gly_ca = Atom::new(
    name="CA",
    coord=Vector3::new(3.5, 2.5, 5.5),
    resname="GLY",
    chainid='A',
    resseq=2,
    occupancy=1.0,
    bfactor=20.0,
    element="C",
  )
  let gly_c = Atom::new(
    name="C",
    coord=Vector3::new(4.0, 2.5, 6.0),
    resname="GLY",
    chainid='A',
    resseq=2,
    occupancy=1.0,
    bfactor=20.0,
    element="C",
  )
  let gly_o = Atom::new(
    name="O",
    coord=Vector3::new(4.5, 2.5, 6.5),
    resname="GLY",
    chainid='A',
    resseq=2,
    occupancy=1.0,
    bfactor=20.0,
    element="O",
  )
  let gly = Residue::new(
    resname="GLY",
    chainid='A',
    resseq=2,
    atoms=[gly_n, gly_ca, gly_c, gly_o],
  )
  let chain = Chain::new(id='A', residues=[ala, gly])
  let model = Model::new(id=1, chains=[chain])
  Structure::new(id="1ABC", models=[model])
}
