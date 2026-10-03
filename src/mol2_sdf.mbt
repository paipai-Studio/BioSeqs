///|
/// Mol2 (Sybyl) format I/O + SDF writer + structure conversions.
///
/// Mol2 parsing/writing is new. SDF parsing already exists in chemminer.mbt
/// (chem_parse_sdf); this module adds sdf_to_string and conversions to
/// the PDB Chain/Residue/Atom model.

// ===== Mol2 =====

///|
/// A Mol2 atom.
pub(all) struct Mol2Atom {
  serial : Int
  name : String
  x : Double
  y : Double
  z : Double
  sybyl_type : String
  subst_id : Int
  subst_name : String
  charge : Double
}

///|
/// A Mol2 bond.
pub(all) struct Mol2Bond {
  serial : Int
  origin : Int
  target : Int
  bond_type : String
}

///|
/// A Mol2 molecule.
pub struct Mol2Molecule {
  name : String
  atoms : Array[Mol2Atom]
  bonds : Array[Mol2Bond]
  mol_type : String
  charge_type : String
}

///|
fn mol2_pd(s : String) -> Double {
  try {
    @string.parse_double(s)
  } catch {
    _ => 0.0
  }
}

///|
fn mol2_pi(s : String) -> Int {
  try {
    @string.parse_int(s)
  } catch {
    _ => 0
  }
}

///|
fn mol2_trim(s : StringView) -> String {
  let n = s.length()
  let mut start = 0
  while start < n {
    let ch = s[start:start + 1]
    if ch == " " || ch == "\t" || ch == "\r" || ch == "\n" {
      start = start + 1
    } else {
      break
    }
  }
  let mut end = n
  while end > start {
    let ch = s[end - 1:end]
    if ch == " " || ch == "\t" || ch == "\r" || ch == "\n" {
      end = end - 1
    } else {
      break
    }
  }
  if end <= start {
    ""
  } else {
    s[start:end].to_string()
  }
}

///|
fn split_fields_sv(line : StringView) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  let n = line.length()
  let mut i = 0
  while i < n {
    while i < n {
      let ch = line[i:i + 1]
      if ch == " " || ch == "\t" {
        i = i + 1
      } else {
        break
      }
    }
    if i >= n {
      break
    }
    let start = i
    while i < n {
      let ch = line[i:i + 1]
      if ch == " " || ch == "\t" {
        break
      }
      i = i + 1
    }
    result.push(line[start:i].to_string())
  }
  result
}

///|
/// Parse a Mol2 file content into an array of molecules.
pub fn parse_mol2(content : String) -> Array[Mol2Molecule] {
  let lines = content.split("\n").to_array()
  let molecules : Array[Mol2Molecule] = Array(capacity=0)
  let mut i = 0
  while i < lines.length() {
    let line = mol2_trim(lines[i])
    if line == "@<TRIPOS>MOLECULE" {
      i = i + 1
      let mut name = ""
      let mut mol_type = ""
      let mut charge_type = ""
      if i < lines.length() {
        name = mol2_trim(lines[i])
        i = i + 1
      }
      if i < lines.length() {
        i = i + 1 // counts line
      }
      if i < lines.length() {
        mol_type = mol2_trim(lines[i])
        i = i + 1
      }
      if i < lines.length() {
        charge_type = mol2_trim(lines[i])
        i = i + 1
      }
      // comment (skip if present and not a section)
      if i < lines.length() {
        let t = mol2_trim(lines[i])
        if t.length() > 0 && t[0] != '@' {
          i = i + 1
        }
      }
      let atoms : Array[Mol2Atom] = Array(capacity=0)
      let bonds : Array[Mol2Bond] = Array(capacity=0)
      while i < lines.length() {
        let section = mol2_trim(lines[i])
        if section == "@<TRIPOS>ATOM" {
          i = i + 1
          while i < lines.length() {
            let t = mol2_trim(lines[i])
            if t.length() == 0 || t[0] == '@' {
              break
            }
            let parts = split_fields_sv(lines[i])
            if parts.length() >= 6 {
              atoms.push(Mol2Atom::{
                serial: mol2_pi(parts[0]),
                name: parts[1],
                x: mol2_pd(parts[2]),
                y: mol2_pd(parts[3]),
                z: mol2_pd(parts[4]),
                sybyl_type: if parts.length() > 5 { parts[5] } else { "" },
                subst_id: if parts.length() > 6 { mol2_pi(parts[6]) } else { 0 },
                subst_name: if parts.length() > 7 { parts[7] } else { "" },
                charge: if parts.length() > 8 { mol2_pd(parts[8]) } else { 0.0 },
              })
            }
            i = i + 1
          }
        } else if section == "@<TRIPOS>BOND" {
          i = i + 1
          while i < lines.length() {
            let t = mol2_trim(lines[i])
            if t.length() == 0 || t[0] == '@' {
              break
            }
            let parts = split_fields_sv(lines[i])
            if parts.length() >= 4 {
              bonds.push(Mol2Bond::{
                serial: mol2_pi(parts[0]),
                origin: mol2_pi(parts[1]),
                target: mol2_pi(parts[2]),
                bond_type: parts[3],
              })
            }
            i = i + 1
          }
        } else if section == "@<TRIPOS>MOLECULE" {
          break
        } else if section.length() > 0 && section[0] == '@' {
          i = i + 1
          while i < lines.length() {
            let t = mol2_trim(lines[i])
            if t.length() > 0 && t[0] == '@' {
              break
            }
            i = i + 1
          }
        } else {
          i = i + 1
        }
      }
      molecules.push(Mol2Molecule::{ name, atoms, bonds, mol_type, charge_type })
    } else {
      i = i + 1
    }
  }
  molecules
}

///|
/// Convert a Mol2Molecule to a Mol2 file string.
pub fn mol2_to_string(mol : Mol2Molecule) -> String {
  let sb = StringBuilder::new()
  sb.write_string("@<TRIPOS>MOLECULE\n")
  sb.write_string(mol.name)
  sb.write_char('\n')
  sb.write_string(mol.atoms.length().to_string())
  sb.write_char(' ')
  sb.write_string(mol.bonds.length().to_string())
  sb.write_char('\n')
  sb.write_string(mol.mol_type)
  sb.write_char('\n')
  sb.write_string(mol.charge_type)
  sb.write_char('\n')
  sb.write_string("@<TRIPOS>ATOM\n")
  for atom in mol.atoms {
    sb.write_string(atom.serial.to_string())
    sb.write_char(' ')
    sb.write_string(atom.name)
    sb.write_char(' ')
    sb.write_string(atom.x.to_string())
    sb.write_char(' ')
    sb.write_string(atom.y.to_string())
    sb.write_char(' ')
    sb.write_string(atom.z.to_string())
    sb.write_char(' ')
    sb.write_string(atom.sybyl_type)
    sb.write_char(' ')
    sb.write_string(atom.subst_id.to_string())
    sb.write_char(' ')
    sb.write_string(atom.subst_name)
    sb.write_char(' ')
    sb.write_string(atom.charge.to_string())
    sb.write_char('\n')
  }
  if mol.bonds.length() > 0 {
    sb.write_string("@<TRIPOS>BOND\n")
    for bond in mol.bonds {
      sb.write_string(bond.serial.to_string())
      sb.write_char(' ')
      sb.write_string(bond.origin.to_string())
      sb.write_char(' ')
      sb.write_string(bond.target.to_string())
      sb.write_char(' ')
      sb.write_string(bond.bond_type)
      sb.write_char('\n')
    }
  }
  sb.to_string()
}

///|
/// Convert a Mol2Molecule to a Chain (grouping by subst_id).
pub fn mol2_to_chain(mol : Mol2Molecule) -> Chain {
  let subst_map : Map[Int, Array[Mol2Atom]] = Map([])
  for atom in mol.atoms {
    let key = atom.subst_id
    match subst_map.get(key) {
      Some(arr) => arr.push(atom)
      None => subst_map[key] = [atom]
    }
  }
  let subst_ids : Array[Int] = Array(capacity=0)
  for k, _ in subst_map {
    subst_ids.push(k)
  }
  subst_ids.sort()
  let residues : Array[Residue] = Array(capacity=0)
  for sid in subst_ids {
    match subst_map.get(sid) {
      Some(atoms) => {
        let resname = if atoms[0].subst_name.length() > 0 {
          atoms[0].subst_name
        } else {
          "UNK"
        }
        let pdb_atoms : Array[Atom] = Array(capacity=0)
        for a in atoms {
          let elem_parts = a.sybyl_type.split(".").to_array()
          let elem = if elem_parts.length() > 0 { elem_parts[0].to_string() } else { "C" }
          pdb_atoms.push(Atom::new(
            name=a.name,
            coord=Vector3::new(a.x, a.y, a.z),
            resname~,
            chainid='A',
            resseq=sid,
            element=elem,
          ))
        }
        residues.push(Residue::new(resname~, chainid='A', resseq=sid, atoms=pdb_atoms))
      }
      None => ()
    }
  }
  Chain::new(id='A', residues~)
}

// ===== SDF writer (using existing SdfMolecule from chemminer) =====

///|
fn format_int_3d(n : Int) -> String {
  if n >= 0 && n < 10 {
    "  " + n.to_string()
  } else if n >= 10 && n < 100 {
    " " + n.to_string()
  } else if n >= 0 {
    n.to_string()
  } else if n > -10 {
    " " + n.to_string()
  } else {
    n.to_string()
  }
}

///|
fn format_double_10d(d : Double) -> String {
  let s = d.to_string()
  if s.length() >= 10 {
    s.substring(start=0, end=10)
  } else {
    let mut result = ""
    let padding = 10 - s.length()
    for j = 0; j < padding; j = j + 1 {
      result = result + " "
    }
    result + s
  }
}

///|
fn mol2_pad_left(s : String, width : Int) -> String {
  if s.length() >= width {
    s
  } else {
    let mut result = ""
    for j = 0; j < width - s.length(); j = j + 1 {
      result = result + " "
    }
    result + s
  }
}

///|
/// Convert an SdfMolecule to an SDF file string.
pub fn sdf_to_string(mol : SdfMolecule) -> String {
  let sb = StringBuilder::new()
  sb.write_string(mol.name)
  sb.write_char('\n')
  sb.write_string("  BioSeqs\n")
  sb.write_char('\n')
  let na = mol.atoms.length()
  let nb = mol.bonds.length()
  sb.write_string(format_int_3d(na))
  sb.write_string(format_int_3d(nb))
  sb.write_string("  0  0  0  0  0  0  0  0  0  0999 V2000\n")
  for atom in mol.atoms {
    sb.write_string(format_double_10d(atom.x))
    sb.write_string(format_double_10d(atom.y))
    sb.write_string(format_double_10d(atom.z))
    sb.write_char(' ')
    sb.write_string(mol2_pad_left(atom.symbol, 3))
    sb.write_string("  0")
    let chg_code = match atom.charge {
      0 => 0
      3 => 1
      2 => 2
      1 => 3
      -1 => 5
      -2 => 6
      _ => 0
    }
    sb.write_string(format_int_3d(chg_code))
    sb.write_string("  0  0  0  0  0  0  0  0  0  0  0  0\n")
  }
  for bond in mol.bonds {
    sb.write_string(format_int_3d(bond.atom1 + 1))
    sb.write_string(format_int_3d(bond.atom2 + 1))
    sb.write_string(format_int_3d(bond.bond_type))
    sb.write_string("  0  0  0\n")
  }
  sb.write_string("M  END\n")
  for key, value in mol.properties {
    sb.write_string("> <")
    sb.write_string(key)
    sb.write_string(">\n")
    sb.write_string(value)
    sb.write_char('\n')
    sb.write_char('\n')
  }
  sb.write_string("$$$$\n")
  sb.to_string()
}

///|
/// Convert an SdfMolecule to a Chain (each atom → its own residue).
pub fn sdf_to_chain(mol : SdfMolecule) -> Chain {
  let residues : Array[Residue] = Array(capacity=0)
  for i = 0; i < mol.atoms.length(); i = i + 1 {
    let a = mol.atoms[i]
    let resname = a.symbol
    let atom = Atom::new(
      name=a.symbol,
      coord=Vector3::new(a.x, a.y, a.z),
      resname~,
      chainid='A',
      resseq=i + 1,
      element=a.symbol,
    )
    residues.push(Residue::new(resname~, chainid='A', resseq=i + 1, atoms=[atom]))
  }
  Chain::new(id='A', residues~)
}

///|
/// Convert a Chain to a Mol2Molecule (each residue → a substructure).
pub fn chain_to_mol2(chain : Chain, name : String) -> Mol2Molecule {
  let atoms : Array[Mol2Atom] = Array(capacity=0)
  let bonds : Array[Mol2Bond] = Array(capacity=0)
  let mut serial = 1
  let residues = chain.get_residues()
  for res in residues {
    let resname = res.resname
    let resseq = res.resseq
    for atom in res.get_atoms() {
      let sybyl = atom.element
      atoms.push(Mol2Atom::{
        serial,
        name: atom.name,
        x: atom.coord.x,
        y: atom.coord.y,
        z: atom.coord.z,
        sybyl_type: sybyl,
        subst_id: resseq,
        subst_name: resname,
        charge: 0.0,
      })
      serial = serial + 1
    }
  }
  Mol2Molecule::{ name, atoms, bonds, mol_type: "SMALL", charge_type: "NO_CHARGES" }
}

///|
/// Convert a Chain to an SdfMolecule.
pub fn chain_to_sdf(chain : Chain, name : String) -> SdfMolecule {
  let atoms : Array[SdfAtom] = Array(capacity=0)
  let residues = chain.get_residues()
  for res in residues {
    for atom in res.get_atoms() {
      atoms.push(SdfAtom::new(atom.element, atom.coord.x, atom.coord.y, atom.coord.z, 0))
    }
  }
  SdfMolecule::new(name, atoms, [])
}

// ===== Mol2Molecule accessors =====

///|
pub fn Mol2Molecule::get_name(self : Mol2Molecule) -> String {
  self.name
}

///|
pub fn Mol2Molecule::get_atoms(self : Mol2Molecule) -> Array[Mol2Atom] {
  self.atoms
}

///|
pub fn Mol2Molecule::get_bonds(self : Mol2Molecule) -> Array[Mol2Bond] {
  self.bonds
}

///|
pub fn Mol2Molecule::num_atoms(self : Mol2Molecule) -> Int {
  self.atoms.length()
}

///|
pub fn Mol2Molecule::num_bonds(self : Mol2Molecule) -> Int {
  self.bonds.length()
}

///|
pub fn Mol2Molecule::get_mol_type(self : Mol2Molecule) -> String {
  self.mol_type
}

///|
pub fn Mol2Molecule::get_charge_type(self : Mol2Molecule) -> String {
  self.charge_type
}

// ===== Mol2Atom accessors =====

///|
pub fn Mol2Atom::get_serial(self : Mol2Atom) -> Int {
  self.serial
}

///|
pub fn Mol2Atom::get_name(self : Mol2Atom) -> String {
  self.name
}

///|
pub fn Mol2Atom::get_x(self : Mol2Atom) -> Double {
  self.x
}

///|
pub fn Mol2Atom::get_y(self : Mol2Atom) -> Double {
  self.y
}

///|
pub fn Mol2Atom::get_z(self : Mol2Atom) -> Double {
  self.z
}

///|
pub fn Mol2Atom::get_sybyl_type(self : Mol2Atom) -> String {
  self.sybyl_type
}

///|
pub fn Mol2Atom::get_subst_id(self : Mol2Atom) -> Int {
  self.subst_id
}

///|
pub fn Mol2Atom::get_subst_name(self : Mol2Atom) -> String {
  self.subst_name
}

///|
pub fn Mol2Atom::get_charge(self : Mol2Atom) -> Double {
  self.charge
}

// ===== Mol2Bond accessors =====

///|
pub fn Mol2Bond::get_serial(self : Mol2Bond) -> Int {
  self.serial
}

///|
pub fn Mol2Bond::get_origin(self : Mol2Bond) -> Int {
  self.origin
}

///|
pub fn Mol2Bond::get_target(self : Mol2Bond) -> Int {
  self.target
}

///|
pub fn Mol2Bond::get_bond_type(self : Mol2Bond) -> String {
  self.bond_type
}
