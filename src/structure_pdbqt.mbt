///|
/// PDBQT format reader/writer for AutoDock Vina.
/// Ported from biotite.structure.io.pdbqt.
///
/// PDBQT extends PDB with partial charges (Q) and AutoDock atom
/// types (T) in columns 71-79.

///|
/// A PDBQT atom record.
pub struct PdbqtAtom {
  serial : Int32
  atom_name : String
  res_name : String
  chain : String
  res_id : Int32
  x : Double
  y : Double
  z : Double
  partial_charge : Double
  atom_type : String // AutoDock type: A, C, HD, HS, N, NA, OA, OS, SA, etc.
  is_hetatm : Bool
  // For flexible docking: root/branch markers
  branch_marker : String // "", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF"
  branch_data : (Int32, Int32) // For BRANCH: (start, end) atom serials
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
    if record == "ATOM  " ||
      record == "ATOM" ||
      record == "HETATM" ||
      record == "HETATM" {
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
  let serial = if len >= 11 { pdb_parse_int_fast(line[6:11].trim()) } else { 0 }
  let atom_name = if len >= 16 { line[12:16].trim().to_owned() } else { "" }
  let res_name = if len >= 20 { line[17:20].trim().to_owned() } else { "" }
  let chain = if len >= 22 {
    let c = line[21:22]
    if c.trim().length() == 0 {
      " "
    } else {
      c.to_owned()
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
fn parse_pdbqt_charge_type(line : StringView, _len : Int32) -> (Double, String) {
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
      return (charge, type_str.to_owned())
    }
  }
  (0.0, "")
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
