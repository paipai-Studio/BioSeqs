///|
/// PQR format (PDB-like topology plus per-atom charge and radius) reader and
/// writer. Commonly produced/required by APBS, DelPhi and CHARMM-GUI solvers.
///
/// Each record is an ATOM/HETATM line whose two trailing whitespace-separated
/// columns are the atomic charge and radius:
/// `ATOM serial name resname [chain] resid x y z charge radius`.

///|
/// One PQR atom record.
pub(all) struct PqrAtom {
  serial : Int32
  name : String
  resname : String
  chain : String
  resid : Int32
  x : Double
  y : Double
  z : Double
  charge : Double
  radius : Double
}

///|
/// A parsed PQR structure: an ordered list of atom records.
pub struct PqrStructure {
  atoms : Array[PqrAtom]
}

///|
/// Error raised for malformed PQR input.
pub suberror PqrError {
  PqrError(String)
}

///|
/// Default constructor for PqrStructure.
pub fn PqrStructure::new(atoms : Array[PqrAtom]) -> PqrStructure {
  PqrStructure::{ atoms, }
}

///|
/// Number of atoms in the structure.
pub fn PqrStructure::natoms(self : PqrStructure) -> Int32 {
  self.atoms.length()
}

///|
/// Sum of all partial charges (net charge of the structure).
pub fn PqrStructure::total_charge(self : PqrStructure) -> Double {
  let mut total = 0.0
  for a in self.atoms {
    total = total + a.charge
  }
  total
}

///|
/// Axis-aligned bounding box as `(min, max)` arrays of `[x, y, z]`. Both
/// arrays are empty when the structure has no atoms.
pub fn PqrStructure::bounds(
  self : PqrStructure,
) -> (Array[Double], Array[Double]) {
  let mins : Array[Double] = Array(capacity=3)
  let maxs : Array[Double] = Array(capacity=3)
  if self.atoms.length() == 0 {
    return (mins, maxs)
  }
  let first = self.atoms[0]
  mins.push(first.x)
  mins.push(first.y)
  mins.push(first.z)
  maxs.push(first.x)
  maxs.push(first.y)
  maxs.push(first.z)
  for k = 1; k < self.atoms.length(); k = k + 1 {
    let a = self.atoms[k]
    if a.x < mins[0] {
      mins[0] = a.x
    }
    if a.y < mins[1] {
      mins[1] = a.y
    }
    if a.z < mins[2] {
      mins[2] = a.z
    }
    if a.x > maxs[0] {
      maxs[0] = a.x
    }
    if a.y > maxs[1] {
      maxs[1] = a.y
    }
    if a.z > maxs[2] {
      maxs[2] = a.z
    }
  }
  (mins, maxs)
}

///|
/// Atom coordinates as an `[atom][xyz]` array.
pub fn pqr_coords(s : PqrStructure) -> Array[Array[Double]] {
  let coords : Array[Array[Double]] = Array(capacity=s.atoms.length())
  for a in s.atoms {
    coords.push([a.x, a.y, a.z])
  }
  coords
}

///|
/// Split a PQR line into whitespace-separated tokens.
fn pqr_tokens(line : StringView) -> Array[String] {
  let tokens : Array[String] = Array(capacity=0)
  let parts = line.split(" ")
  for p in parts {
    let t = p.trim()
    if t.length() > 0 {
      tokens.push(t.to_owned())
    }
  }
  tokens
}

///|
/// Parse a double token, defaulting to 0.0 when missing or malformed.
fn pqr_opt_double(s : String) -> Double {
  match parse_double(s) {
    Some(v) => v
    None => 0.0
  }
}

///|
/// Parse a single ATOM/HETATM line. Whitespace splitting is tried first:
/// with 11+ fields the chain column is present, with exactly 10 fields the
/// line has no chain id. When the token count is lower, fall back to the
/// fixed-column PDB layout (first 54 columns) with charge/radius taken from
/// the trailing tokens. Raises when the line is too short for either mode.
fn pqr_parse_atom_line(line : StringView) -> PqrAtom raise PqrError {
  let tokens = pqr_tokens(line)
  let n = tokens.length()
  if n >= 11 {
    // record serial name resname chain resid x y z charge radius [extra...]
    PqrAtom::{
      serial: parse_int(tokens[1]),
      name: tokens[2],
      resname: tokens[3],
      chain: tokens[4],
      resid: parse_int(tokens[5]),
      x: pqr_opt_double(tokens[6]),
      y: pqr_opt_double(tokens[7]),
      z: pqr_opt_double(tokens[8]),
      charge: pqr_opt_double(tokens[n - 2]),
      radius: pqr_opt_double(tokens[n - 1]),
    }
  } else if n == 10 {
    // record serial name resname resid x y z charge radius (no chain)
    PqrAtom::{
      serial: parse_int(tokens[1]),
      name: tokens[2],
      resname: tokens[3],
      chain: "",
      resid: parse_int(tokens[4]),
      x: pqr_opt_double(tokens[5]),
      y: pqr_opt_double(tokens[6]),
      z: pqr_opt_double(tokens[7]),
      charge: pqr_opt_double(tokens[8]),
      radius: pqr_opt_double(tokens[9]),
    }
  } else if line.length() >= 54 && n >= 3 {
    // Fixed-column PDB fallback: coords in columns 31-54, charge/radius in
    // the last two whitespace-separated fields.
    PqrAtom::{
      serial: pdb_parse_int_fast(line[6:11].trim()),
      name: line[12:16].trim().to_owned(),
      resname: line[17:20].trim().to_owned(),
      chain: line[21:22].trim().to_owned(),
      resid: pdb_parse_int_fast(line[22:26].trim()),
      x: pdb_parse_double_fast(line[30:38].trim()),
      y: pdb_parse_double_fast(line[38:46].trim()),
      z: pdb_parse_double_fast(line[46:54].trim()),
      charge: pqr_opt_double(tokens[n - 2]),
      radius: pqr_opt_double(tokens[n - 1]),
    }
  } else {
    raise PqrError(
      "PQR atom line has too few columns (\{n} fields): \{line.to_owned()}",
    )
  }
}

///|
/// Parse PQR text into a structure. Only ATOM/HETATM records are consumed;
/// REMARK, TER, END, CRYST1 and any other lines are skipped.
pub fn pqr_parse(text : String) -> PqrStructure raise PqrError {
  let atoms : Array[PqrAtom] = Array(capacity=0)
  let lines = text.split("\n").to_array()
  for line in lines {
    let trimmed = line.trim()
    if trimmed.length() == 0 {
      continue
    }
    let rlen = trimmed.length()
    let cut = if rlen < 6 { rlen } else { 6 }
    let record = trimmed[0:cut].trim().to_owned()
    if record == "ATOM" || record == "HETATM" {
      atoms.push(pqr_parse_atom_line(line))
    }
  }
  PqrStructure::{ atoms, }
}

///|
/// Serialize a structure to PQR text with PDB-style aligned columns:
/// 6 + 7 + 1 + 4 + 1 + 3 + 1 + 1 + 4 + 4 + 8*3 + 10*2 chars per line.
/// Atoms without a chain id keep the no-chain field order when re-read.
pub fn pqr_write(s : PqrStructure) -> String {
  let buf = StringBuilder()
  for a in s.atoms {
    buf.write_string("ATOM  ")
    buf.write_string(pad_left(a.serial.to_string(), 7, ' '))
    buf.write_char(' ')
    buf.write_string(pad_right(a.name, 4, ' '))
    buf.write_char(' ')
    buf.write_string(pad_right(a.resname, 3, ' '))
    buf.write_char(' ')
    buf.write_string(pad_right(a.chain, 1, ' '))
    buf.write_string(pad_left(a.resid.to_string(), 4, ' '))
    buf.write_string("    ")
    buf.write_string(format_double(a.x, 8, 3))
    buf.write_string(format_double(a.y, 8, 3))
    buf.write_string(format_double(a.z, 8, 3))
    buf.write_string(format_double(a.charge, 10, 4))
    buf.write_string(format_double(a.radius, 10, 4))
    buf.write_char('\n')
  }
  buf.to_string()
}
