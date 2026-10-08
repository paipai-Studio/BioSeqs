///|
/// PDBQT format — PDB with AutoDock atom types and partial charges.
/// Ported from biotite.structure.io.pdbqt.
///
/// PDBQT extends PDB with:
///   - Partial charges in the B-factor column (columns 66-76)
///   - Atom type in columns 78-79 (e.g. HD, OA, N, A)
///   - ROOT/BRANCH/ENDBRANCH/TORSDOF directives for flexible torsions

///|
/// AutoDock atom type mapping.
/// Maps PDB element + atom name to AutoDock atom type.
fn autodock_atom_type(element : String, atom_name : String) -> String {
  let e = to_uppercase(element)
  let a = to_uppercase(atom_name)
  // Hydrogen
  if e == "H" {
    return "HD"
  }
  // Carbon: aromatic vs aliphatic
  if e == "C" {
    // Aromatic carbons (ring atoms in PHE, TYR, TRP, HIS)
    if a == "CG" || a == "CD1" || a == "CD2" || a == "CE1" || a == "CE2" || a == "CZ" || a == "CE3" || a == "CH2" || a == "CZ2" || a == "CZ3" || a == "NE1" {
      return "A"
    }
    return "C"
  }
  // Oxygen: acceptor vs hydroxyl
  if e == "O" {
    // Hydroxyl oxygens
    if a == "OG" || a == "OG1" || a == "OH" {
      return "OA"
    }
    return "OA"
  }
  // Nitrogen: donor vs acceptor
  if e == "N" {
    return "NA"
  }
  // Sulfur
  if e == "S" {
    return "SA"
  }
  // F, Cl, Br, I
  if e == "F" {
    return "F"
  }
  if e == "CL" {
    return "CL"
  }
  if e == "BR" {
    return "BR"
  }
  if e == "I" {
    return "I"
  }
  // Metals
  if e == "MG" || e == "MN" || e == "ZN" || e == "CA" || e == "FE" {
    return "M"
  }
  // Default: use element
  e
}

///|
/// Estimate partial charge for an atom (very rough — for real use, read from input).
fn estimate_partial_charge(resname : String, atom_name : String) -> Double {
  let r = to_uppercase(resname)
  let a = to_uppercase(atom_name)
  // Backbone
  if a == "N" {
    return -0.3
  }
  if a == "CA" {
    return 0.1
  }
  if a == "C" {
    return 0.2
  }
  if a == "O" {
    return -0.4
  }
  if a == "OXT" {
    return -0.4
  }
  // Side chain charges
  if r == "ASP" && (a == "OD1" || a == "OD2") {
    return -0.5
  }
  if r == "GLU" && (a == "OE1" || a == "OE2") {
    return -0.5
  }
  if r == "LYS" && a == "NZ" {
    return 0.5
  }
  if r == "ARG" && (a == "NE" || a == "NH1" || a == "NH2") {
    return 0.3
  }
  0.0
}

///|
/// Write PDBQT string from atoms.
/// Includes AutoDock atom types and partial charges.
pub fn write_pdbqt(atoms : Array[Atom]) -> String {
  let sb = StringBuilder()
  sb.write_string("ROOT\n")
  for i = 0; i < atoms.length(); i = i + 1 {
    let a = atoms[i]
    let ad_type = autodock_atom_type(a.element, a.name)
    let charge = estimate_partial_charge(a.resname, a.name)
    sb.write_string(format_pdbqt_line(i + 1, a, ad_type, charge))
  }
  sb.write_string("ENDROOT\n")
  sb.write_string("TER\n")
  sb.write_string("END\n")
  sb.to_string()
}

///|
/// Format a single PDBQT atom line.
fn format_pdbqt_line(
  serial : Int,
  atom : Atom,
  ad_type : String,
  charge : Double,
) -> String {
  let sb = StringBuilder()
  // Record name: "ATOM  "
  sb.write_string("ATOM  ")
  // Serial: %5d
  sb.write_string(format_int_field(serial, 5))
  // Space
  sb.write_char(' ')
  // Atom name: %-4s
  sb.write_string(format_str_field(atom.name, 4))
  // Alternate location indicator
  sb.write_char(' ')
  // Residue name: %3s
  sb.write_string(format_str_field_right(atom.resname, 3))
  // Space
  sb.write_char(' ')
  // Chain ID
  sb.write_char(atom.chainid)
  // Residue sequence: %4d
  sb.write_string(format_int_field(atom.resseq, 4))
  // iCode
  sb.write_char(' ')
  // Space
  sb.write_char(' ')
  // X: %8.3f
  sb.write_string(format_double_field(atom.coord.x, 8, 3))
  // Y: %8.3f
  sb.write_string(format_double_field(atom.coord.y, 8, 3))
  // Z: %8.3f
  sb.write_string(format_double_field(atom.coord.z, 8, 3))
  // Occupancy: %6.2f
  sb.write_string(format_double_field(atom.occupancy, 6, 2))
  // B-factor (partial charge): %6.3f
  sb.write_string(format_double_field(charge, 6, 3))
  // Space + atom type
  sb.write_char(' ')
  sb.write_string(format_str_field_right(ad_type, 2))
  sb.write_char('\n')
  sb.to_string()
}

///|
/// Parse PDBQT string into atoms.
/// Extracts coordinates, atom types, and partial charges.
pub fn parse_pdbqt(content : String) -> Array[Atom] {
  let atoms : Array[Atom] = Array(capacity=0)
  let raw = content.split("\n").to_array()
  for sv in raw {
    let line = sv.to_owned()
    if line.length() < 6 {
      continue
    }
    // Check for ATOM or HETATM record
    let record = substring_safe(line, 0, 6).trim().to_owned()
    if record != "ATOM" && record != "HETATM" {
      continue
    }
    if line.length() < 54 {
      continue
    }
    // Parse fields
    let name = substring_safe(line, 12, 16).trim().to_owned()
    let resname = substring_safe(line, 17, 20).trim().to_owned()
    let chainid = if line.length() > 21 {
      line.unsafe_get(21).unsafe_to_char()
    } else {
      ' '
    }
    let resseq = parse_int(substring_safe(line, 22, 26).trim().to_owned())
    let x = parse_double_field(line, 30, 38)
    let y = parse_double_field(line, 38, 46)
    let z = parse_double_field(line, 46, 54)
    let coord = Vector3::new(x, y, z)
    // Partial charge (columns 70-76)
    let charge = if line.length() >= 76 {
      parse_double_field(line, 70, 76)
    } else {
      0.0
    }
    // Atom type (columns 77-79)
    let ad_type = if line.length() >= 79 {
      substring_safe(line, 77, 79).trim().to_owned()
    } else {
      ""
    }
    // Derive element from atom type
    let element = ad_type_to_element(ad_type)
    atoms.push(Atom::new(
      name=name,
      coord=coord,
      resname=resname,
      chainid=chainid,
      resseq=resseq,
      element=element,
      charge=charge.to_int(),
    ))
  }
  atoms
}

///|
/// Convert AutoDock atom type to element symbol.
fn ad_type_to_element(ad_type : String) -> String {
  let t = to_uppercase(ad_type)
  if t == "HD" {
    return "H"
  }
  if t == "C" || t == "A" {
    return "C"
  }
  if t == "OA" || t == "OS" || t == "O" {
    return "O"
  }
  if t == "NA" || t == "NS" || t == "N" {
    return "N"
  }
  if t == "SA" || t == "S" {
    return "S"
  }
  if t == "P" {
    return "P"
  }
  if t == "F" {
    return "F"
  }
  if t == "CL" {
    return "CL"
  }
  if t == "BR" {
    return "BR"
  }
  if t == "I" {
    return "I"
  }
  if t == "M" {
    return "MG"
  }
  // Default: first character
  if t.length() > 0 {
    return t[0:1].to_owned()
  }
  ""
}

///|
/// Format string right-justified in a fixed-width field.
fn format_str_field_right(s : String, width : Int) -> String {
  if s.length() >= width {
    return s.substring(start=0, end=width)
  }
  let sb = StringBuilder()
  for _i = 0; _i < width - s.length(); _i = _i + 1 {
    sb.write_char(' ')
  }
  sb.write_string(s)
  sb.to_string()
}
