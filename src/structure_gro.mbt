///|
/// GROMACS GRO coordinate file format — read and write.
/// Ported from biotite.structure.io.gro.
///
/// GRO format specification:
///   Line 1:  title string (≤ 78 chars)
///   Line 2:  number of atoms
///   Lines 3..N+2: atom records (fixed-width columns)
///   Last line:  box vectors (3 doubles, nm)
///
/// Atom record columns (0-indexed):
///   0-4:   residue number (5 chars, right-justified)
///   5-9:   residue name (5 chars, left-justified)
///   10-14: atom name (5 chars, left-justified)
///   15-19: atom number (5 chars, right-justified)
///   20-27: x position (8 chars, 3 decimal places, nm)
///   28-35: y position
///   36-43: z position
///   44-51: x velocity (optional, nm/ps)
///   52-59: y velocity
///   60-67: z velocity

///|
/// A GRO atom record.
pub(all) struct GroAtom {
  resseq : Int
  resname : String
  atom_name : String
  atom_seq : Int
  x : Double       // nm
  y : Double       // nm
  z : Double       // nm
  vx : Double?     // nm/ps
  vy : Double?
  vz : Double?
}

///|
/// A parsed GRO file.
pub(all) struct GroFile {
  title : String
  atoms : Array[GroAtom]
  box : (Double, Double, Double)  // box dimensions (nm)
}

///|
/// Parse a GRO file from a string.
/// Returns the parsed structure with atoms and box vectors.
pub fn parse_gro(content : String) -> GroFile? {
  let raw = content.split("\n").to_array()
  // Convert StringView to String
  let lines : Array[String] = Array(capacity=raw.length())
  for sv in raw {
    lines.push(sv.to_owned())
  }
  if lines.length() < 3 {
    return None
  }
  let title = lines[0].trim().to_owned()
  let n_atoms = parse_int(lines[1].trim().to_owned())
  if lines.length() < n_atoms + 3 {
    return None
  }
  let atoms : Array[GroAtom] = Array(capacity=n_atoms)
  for i = 0; i < n_atoms; i = i + 1 {
    let line = lines[i + 2]
    if line.length() < 44 {
      continue
    }
    let resseq = parse_int_field(line, 0, 5)
    let resname = substring_safe(line, 5, 10).trim().to_owned()
    let atom_name = substring_safe(line, 10, 15).trim().to_owned()
    let atom_seq = parse_int_field(line, 15, 20)
    let x = parse_double_field(line, 20, 28)
    let y = parse_double_field(line, 28, 36)
    let z = parse_double_field(line, 36, 44)
    let vx : Double? = if line.length() >= 52 {
      Some(parse_double_field(line, 44, 52))
    } else {
      None
    }
    let vy : Double? = if line.length() >= 60 {
      Some(parse_double_field(line, 52, 60))
    } else {
      None
    }
    let vz : Double? = if line.length() >= 68 {
      Some(parse_double_field(line, 60, 68))
    } else {
      None
    }
    atoms.push(GroAtom::{
      resseq,
      resname,
      atom_name,
      atom_seq,
      x,
      y,
      z,
      vx,
      vy,
      vz,
    })
  }
  // Box line
  let box_line = lines[n_atoms + 2].trim().to_owned()
  let box_parts = box_line.split(" ").to_array()
  let mut bx = 0.0
  let mut by = 0.0
  let mut bz = 0.0
  let mut idx = 0
  for p in box_parts {
    let trimmed = p.trim().to_owned()
    if trimmed.length() == 0 {
      continue
    }
    match parse_double(trimmed) {
      Some(v) => {
        if idx == 0 {
          bx = v
        } else if idx == 1 {
          by = v
        } else if idx == 2 {
          bz = v
        }
        idx = idx + 1
      }
      None => ()
    }
  }
  Some(GroFile::{ title, atoms, box: (bx, by, bz) })
}

///|
/// Write a GRO file string from a GroFile.
pub fn write_gro(gro : GroFile) -> String {
  let sb = StringBuilder()
  // Title
  sb.write_string(gro.title)
  sb.write_char('\n')
  // Number of atoms
  sb.write_string(gro.atoms.length().to_string())
  sb.write_char('\n')
  // Atom records
  for a in gro.atoms {
    // Residue number: %5d
    sb.write_string(format_int_field(a.resseq, 5))
    // Residue name: %-5s
    sb.write_string(format_str_field(a.resname, 5))
    // Atom name: %-5s
    sb.write_string(format_str_field(a.atom_name, 5))
    // Atom number: %5d
    sb.write_string(format_int_field(a.atom_seq, 5))
    // Positions: %8.3f
    sb.write_string(format_double_field(a.x, 8, 3))
    sb.write_string(format_double_field(a.y, 8, 3))
    sb.write_string(format_double_field(a.z, 8, 3))
    // Velocities (optional)
    match a.vx {
      Some(v) => sb.write_string(format_double_field(v, 8, 4))
      None => ()
    }
    match a.vy {
      Some(v) => sb.write_string(format_double_field(v, 8, 4))
      None => ()
    }
    match a.vz {
      Some(v) => sb.write_string(format_double_field(v, 8, 4))
      None => ()
    }
    sb.write_char('\n')
  }
  // Box line
  sb.write_string(format_double_field(gro.box.0, 10, 5))
  sb.write_char(' ')
  sb.write_string(format_double_field(gro.box.1, 10, 5))
  sb.write_char(' ')
  sb.write_string(format_double_field(gro.box.2, 10, 5))
  sb.write_char('\n')
  sb.to_string()
}

///|
/// Convert GRO atoms to PDB Atom array (nm → Å, multiply by 10).
pub fn gro_to_pdb_atoms(gro : GroFile) -> Array[Atom] {
  let result : Array[Atom] = Array(capacity=gro.atoms.length())
  for ga in gro.atoms {
    let coord = Vector3::new(ga.x * 10.0, ga.y * 10.0, ga.z * 10.0)
    result.push(Atom::new(
      name=ga.atom_name,
      coord~,
      resname=ga.resname,
      resseq=ga.resseq,
    ))
  }
  result
}

///|
/// Convert PDB Atom array to GRO atoms (Å → nm, divide by 10).
pub fn pdb_to_gro_atoms(
  atoms : Array[Atom],
  title~ : String = "Generated by BioSeqs",
  box~ : (Double, Double, Double) = (0.0, 0.0, 0.0),
) -> GroFile {
  let gro_atoms : Array[GroAtom] = Array(capacity=atoms.length())
  for i = 0; i < atoms.length(); i = i + 1 {
    let a = atoms[i]
    gro_atoms.push(GroAtom::{
      resseq: a.resseq,
      resname: a.resname,
      atom_name: a.name,
      atom_seq: i + 1,
      x: a.coord.x / 10.0,
      y: a.coord.y / 10.0,
      z: a.coord.z / 10.0,
      vx: None,
      vy: None,
      vz: None,
    })
  }
  GroFile::{ title, atoms: gro_atoms, box }
}

// --- Helper functions ---

///|
/// Extract substring safely with bounds checking.
fn substring_safe(s : String, start : Int, end : Int) -> String {
  let len = s.length()
  let st = if start < 0 {
    0
  } else if start > len {
    len
  } else {
    start
  }
  let en = if end < 0 {
    0
  } else if end > len {
    len
  } else {
    end
  }
  if st >= en {
    return ""
  }
  s.substring(start=st, end=en)
}

///|
/// Parse integer from a fixed-width field.
fn parse_int_field(line : String, start : Int, end : Int) -> Int {
  let s = substring_safe(line, start, end).trim().to_owned()
  parse_int(s)
}

///|
/// Parse double from a fixed-width field.
fn parse_double_field(line : String, start : Int, end : Int) -> Double {
  let s = substring_safe(line, start, end).trim().to_owned()
  match parse_double(s) {
    Some(v) => v
    None => 0.0
  }
}

///|
/// Format integer in a fixed-width field (right-justified).
fn format_int_field(val : Int, width : Int) -> String {
  let s = val.to_string()
  if s.length() >= width {
    return s.substring(start=s.length() - width, end=s.length())
  }
  let sb = StringBuilder()
  for _i = 0; _i < width - s.length(); _i = _i + 1 {
    sb.write_char(' ')
  }
  sb.write_string(s)
  sb.to_string()
}

///|
/// Format string in a fixed-width field (left-justified).
fn format_str_field(s : String, width : Int) -> String {
  if s.length() >= width {
    return s.substring(start=0, end=width)
  }
  let sb = StringBuilder()
  sb.write_string(s)
  for _i = 0; _i < width - s.length(); _i = _i + 1 {
    sb.write_char(' ')
  }
  sb.to_string()
}

///|
/// Format double in a fixed-width field with decimal places (right-justified).
fn format_double_field(val : Double, width : Int, decimals : Int) -> String {
  let mut v = val
  let negative = v < 0.0
  if negative {
    v = -v
  }
  // Round to specified decimals
  let multiplier = pow10(decimals)
  v = (v * multiplier + 0.5).to_int().to_double() / multiplier
  // Build string
  let int_part = v.to_int()
  let frac_val = (v - int_part.to_double()) * multiplier + 0.5
  let frac_int = frac_val.to_int()
  let int_str = int_part.to_string()
  let frac_str = frac_int.to_string()
  // Pad fractional part
  let frac_padded = StringBuilder()
  for _i = frac_str.length(); _i < decimals; _i = _i + 1 {
    frac_padded.write_char('0')
  }
  frac_padded.write_string(frac_str)
  let full = if negative {
    "-" + int_str + "." + frac_padded.to_string()
  } else {
    int_str + "." + frac_padded.to_string()
  }
  // Right-justify to width
  if full.length() >= width {
    return full
  }
  let sb = StringBuilder()
  for _i = 0; _i < width - full.length(); _i = _i + 1 {
    sb.write_char(' ')
  }
  sb.write_string(full)
  sb.to_string()
}

///|
/// Integer power of 10.
fn pow10(n : Int) -> Double {
  let mut result = 1.0
  for _i = 0; _i < n; _i = _i + 1 {
    result = result * 10.0
  }
  result
}
