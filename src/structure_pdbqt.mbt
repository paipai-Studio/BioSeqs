///|
/// PDBQT format reader/writer for AutoDock Vina.
/// Ported from biotite.structure.io.pdbqt.
///
/// PDBQT extends PDB with partial charges (Q) and AutoDock atom
/// types (T) in columns 71-79.

///|
/// A PDBQT atom record.
pub struct PdbqtAtom {
  serial : Int
  atom_name : String
  res_name : String
  chain : String
  res_id : Int
  x : Double
  y : Double
  z : Double
  partial_charge : Double
  atom_type : String // AutoDock type: A, C, HD, HS, N, NA, OA, OS, SA, etc.
  is_hetatm : Bool
  // For flexible docking: root/branch markers
  branch_marker : String // "", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"
  branch_data : (Int, Int) // For BRANCH: (start, end) atom serials
}

///|
/// Parse a PDBQT file string into atom records.
pub fn pdbqt_parse(content : String) -> Array[PdbqtAtom] {
  let lines = content.split("\n").to_array()
  let atoms : Array[PdbqtAtom] = Array(capacity=0)
  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    let rlen = trimmed.length()
    let cut = if rlen < 6 { rlen } else { 6 }
    let record = trimmed[0:cut]
    if record == "ATOM  " || record == "ATOM" || record == "HETATM" || record == "HETATM" {
      atoms.push(parse_pdbqt_atom_line(line))
    }
  }
  atoms
}

///|
/// Parse a single PDBQT ATOM/HETATM line.
fn parse_pdbqt_atom_line(line : StringView) -> PdbqtAtom {
  let len = line.length()
  let cut6 = if len < 6 { len } else { 6 }
  let is_hetatm = line[0:cut6].trim() == "HETATM"
  let serial = if len >= 11 {
    pdb_parse_int_fast(line[6:11].trim())
  } else {
    0
  }
  let atom_name = if len >= 16 {
    line[12:16].trim().to_string()
  } else {
    ""
  }
  let res_name = if len >= 20 {
    line[17:20].trim().to_string()
  } else {
    ""
  }
  let chain = if len >= 22 {
    let c = line[21:22]
    if c.trim().length() == 0 {
      " "
    } else {
      c.to_string()
    }
  } else {
    " "
  }
  let res_id = if len >= 26 {
    pdb_parse_int_fast(line[22:26].trim())
  } else {
    0
  }
  let x = if len >= 38 {
    pdb_parse_double_fast(line[30:38].trim())
  } else {
    0.0
  }
  let y = if len >= 46 {
    pdb_parse_double_fast(line[38:46].trim())
  } else {
    0.0
  }
  let z = if len >= 54 {
    pdb_parse_double_fast(line[46:54].trim())
  } else {
    0.0
  }
  // Partial charge and atom type: last two whitespace-separated fields
  let (partial_charge, atom_type) = parse_pdbqt_charge_type(line, len)
  PdbqtAtom::{
    serial,
    atom_name,
    res_name,
    chain,
    res_id,
    x,
    y,
    z,
    partial_charge,
    atom_type,
    is_hetatm,
    branch_marker: "",
    branch_data: (0, 0),
  }
}

///|
/// Parse partial charge and atom type from the tail of a PDBQT line.
fn parse_pdbqt_charge_type(line : StringView, len : Int) -> (Double, String) {
  // The last two whitespace-separated fields are charge and type
  let tokens : Array[StringView] = Array(capacity=0)
  let parts = line.split(" ")
  for p in parts {
    let t = p.trim()
    if t.length() > 0 {
      tokens.push(t)
    }
  }
  let n = tokens.length()
  if n >= 2 {
    let type_str = tokens[n - 1]
    let charge_str = tokens[n - 2]
    let charge = pdb_parse_double_fast(charge_str)
    if type_str.length() <= 4 {
      return (charge, type_str.to_string())
    }
  }
  (0.0, "")
}

///|
/// Write PDBQT atoms to a PDBQT format string.
pub fn pdbqt_write(atoms : Array[PdbqtAtom]) -> String {
  let sb = StringBuilder()
  for i = 0; i < atoms.length(); i = i + 1 {
    sb.write_string(format_pdbqt_atom_line(atoms[i], i + 1))
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// Format a PDBQT ATOM/HETATM line.
fn format_pdbqt_atom_line(a : PdbqtAtom, serial : Int) -> String {
  let sb = StringBuilder()
  // Record name
  if a.is_hetatm {
    sb.write_string("HETATM")
  } else {
    sb.write_string("ATOM  ")
  }
  // Serial (5 chars, right-aligned)
  sb.write_string(format_int_field(a.serial, 5))
  // Space
  sb.write_char(' ')
  // Atom name (4 chars)
  sb.write_string(format_atom_name_field(a.atom_name))
  // AltLoc
  sb.write_char(' ')
  // Res name (3 chars, right-aligned)
  sb.write_string(format_str_field(a.res_name, 3))
  // Space
  sb.write_char(' ')
  // Chain
  sb.write_string(a.chain)
  // Res seq (4 chars, right-aligned)
  sb.write_string(format_int_field(a.res_id, 4))
  // iCode
  sb.write_char(' ')
  // Spaces 27-30
  sb.write_string("   ")
  // x (8 chars, 3 decimal)
  sb.write_string(format_float_field(a.x, 8, 3))
  // y
  sb.write_string(format_float_field(a.y, 8, 3))
  // z
  sb.write_string(format_float_field(a.z, 8, 3))
  // Partial charge (7 chars, 3 decimal) + space
  sb.write_string(format_float_field(a.partial_charge, 7, 3))
  sb.write_char(' ')
  // Atom type (2 chars, right-aligned)
  sb.write_string(format_str_field(a.atom_type, 2))
  sb.to_string()
}

///|
/// Format integer right-aligned in width chars.
fn format_int_field(n : Int, width : Int) -> String {
  let s = n.to_string()
  if s.length() >= width {
    return s
  }
  let mut pad = ""
  for _i = 0; _i < width - s.length(); _i = _i + 1 {
    pad = pad + " "
  }
  pad + s
}

///|
/// Format string right-aligned in width chars.
fn format_str_field(s : String, width : Int) -> String {
  if s.length() >= width {
    return s
  }
  let mut pad = ""
  for _i = 0; _i < width - s.length(); _i = _i + 1 {
    pad = pad + " "
  }
  pad + s
}

///|
/// Format atom name: 4 chars, with special handling for elements.
fn format_atom_name_field(name : String) -> String {
  if name.length() >= 4 {
    return name
  }
  // If name starts with a digit (e.g. "1HE2"), left-align
  if name.length() > 0 {
    let first = name.substring(start=0, end=1)
    if first >= "0" && first <= "9" {
      let mut pad = name
      for _i = 0; _i < 4 - name.length(); _i = _i + 1 {
        pad = pad + " "
      }
      return pad
    }
  }
  // Otherwise, right-align in 4 chars
  let mut pad = ""
  for _i = 0; _i < 4 - name.length(); _i = _i + 1 {
    pad = pad + " "
  }
  pad + name
}

///|
/// Format float in width chars with precision decimal places.
fn format_float_field(x : Double, width : Int, precision : Int) -> String {
  // Format with precision decimal places
  let mut s = x.to_string()
  let dot_opt = s.find(".")
  let dot = match dot_opt { Some(d) => d; None => -1 }
  if dot < 0 {
    s = s + "."
    for _i = 0; _i < precision; _i = _i + 1 {
      s = s + "0"
    }
  } else {
    let int_part = s.substring(start=0, end=dot)
    let frac_part = s.substring(start=dot + 1, end=s.length())
    let frac_p = if frac_part.length() >= precision {
      frac_part.substring(start=0, end=precision)
    } else {
      let mut f = frac_part
      for _i = 0; _i < precision - frac_part.length(); _i = _i + 1 {
        f = f + "0"
      }
      f
    }
    s = int_part + "." + frac_p
  }
  // Pad to width
  if s.length() >= width {
    return s
  }
  let mut pad = ""
  for _i = 0; _i < width - s.length(); _i = _i + 1 {
    pad = pad + " "
  }
  pad + s
}

///|
/// Get the AutoDock atom type from element and context.
/// Simplified mapping — full PDBQT typing requires AD4 type rules.
pub fn pdbqt_atom_type(element : String, is_aromatic : Bool) -> String {
  let e = to_uppercase(element)
  match e {
    "H" => "HD" // polar hydrogen (default)
    "C" => if is_aromatic { "A" } else { "C" }
    "N" => "NA" // acceptor (default)
    "O" => "OA" // acceptor
    "S" => "SA"
    "F" => "F"
    "CL" => "CL"
    "BR" => "BR"
    "I" => "I"
    "P" => "P"
    "MG" => "MG"
    "ZN" => "ZN"
    "CA" => "CA"
    "MN" => "MN"
    "FE" => "FE"
    _ => e
  }
}
