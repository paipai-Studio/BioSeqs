// ===========================================================================
// component_separation.mbt — Ligand / Water / Ion separation
//
// Automatically identify and extract ligands, water molecules, and ions
// from PDB structures. Also classifies every residue into one of:
// Protein, Nucleic, Water, Ion, Ligand, Other.
// ===========================================================================

///|
/// Classification of a single residue.
pub enum ResidueClass {
  Protein
  Nucleic
  Water
  Ion
  Ligand
  Other
} derive(Eq, Show)

///|
/// Result of separating a structure into component groups.
pub(all) struct SeparationResult {
  protein : Array[Residue]
  nucleic : Array[Residue]
  ligands : Array[Residue]
  waters : Array[Residue]
  ions : Array[Residue]
  other : Array[Residue]
}

///|
/// Summary counts for a SeparationResult.
pub(all) struct SeparationSummary {
  num_protein : Int
  num_nucleic : Int
  num_ligands : Int
  num_waters : Int
  num_ions : Int
  num_other : Int
  num_ligand_atoms : Int
  num_water_atoms : Int
  num_ion_atoms : Int
}

// ---------------------------------------------------------------------------
// Residue-name classification helpers
// ---------------------------------------------------------------------------

///|
/// Check whether a residue name corresponds to a water molecule.
pub fn is_water_resname(resname : String) -> Bool {
  match resname {
    "HOH" | "WAT" | "TIP" | "TIP3" | "TIP4" | "TIP5" |
    "DOD" | "SOL" | "H2O" | "D3O" | "WAT0" | "MOH" => true
    _ => false
  }
}

///|
/// Known metal cation residue names (PDB 3-letter codes).
fn is_metal_ion(resname : String) -> Bool {
  match resname {
    "LI" | "NA" | "K" | "RB" | "CS" | "BE" | "MG" | "CA" | "SR" | "BA" |
    "ZN" | "FE" | "CU" | "MN" | "NI" | "CO" | "CD" | "HG" | "PB" |
    "AL" | "AG" | "AU" | "CR" | "MO" | "W" | "V" | "TI" | "PD" | "PT" |
    "RU" | "RH" | "IR" | "OS" | "RE" | "NB" | "TA" | "ZR" | "GA" | "IN" |
    "LA" | "CE" | "PR" | "ND" | "SM" | "EU" | "GD" | "TB" |
    "DY" | "HO" | "ER" | "TM" | "YB" | "LU" | "TH" | "U" |
    "AS" | "SE" | "TE" | "B" => true
    _ => false
  }
}

///|
/// Known halide / simple anion residue names.
fn is_halide_ion(resname : String) -> Bool {
  match resname {
    "CL" | "BR" | "I" | "F" => true
    _ => false
  }
}

///|
/// Check whether a residue name corresponds to a known ion.
/// Covers metal cations and halide anions.
pub fn is_ion_resname(resname : String) -> Bool {
  if is_metal_ion(resname) || is_halide_ion(resname) {
    return true
  }
  // Some PDB files use charge notation: NA+, K+, MG2, CA2, FE2, FE3, etc.
  let n = resname.length()
  if n >= 2 && n <= 4 {
    let last = resname.unsafe_get(n - 1).to_int()
    let plus = '+'.to_int()
    let minus = '-'.to_int()
    if last == plus || last == minus {
      // Strip trailing + or -
      let base = resname.substring(start=0, end=n - 1)
      if base != resname {
        return is_metal_ion(base) || is_halide_ion(base)
      }
    }
    // Also try stripping a trailing digit: MG2, CA2, FE3 -> MG, CA, FE
    if n >= 3 {
      let last_char = resname.unsafe_get(n - 1).to_int()
      let zero = '0'.to_int()
      let nine = '9'.to_int()
      if last_char >= zero && last_char <= nine {
        let base2 = resname.substring(start=0, end=n - 1)
        return is_metal_ion(base2) || is_halide_ion(base2)
      }
    }
  }
  false
}

///|
/// Known metal element symbols (for single-atom fallback detection).
fn is_metal_element(element : String) -> Bool {
  let e = to_uppercase(element)
  is_metal_ion(e)
}

///|
/// Known halogen element symbols.
fn is_halogen_element(element : String) -> Bool {
  let e = to_uppercase(element)
  is_halide_ion(e)
}

///|
/// Standard amino acid residue names (20 canonical + common modifications).
fn is_standard_aa(resname : String) -> Bool {
  match resname {
    "ALA" | "ARG" | "ASN" | "ASP" | "CYS" | "GLN" | "GLU" |
    "GLY" | "HIS" | "ILE" | "LEU" | "LYS" | "MET" | "PHE" |
    "PRO" | "SER" | "THR" | "TRP" | "TYR" | "VAL" |
    "MSE" | "SEC" | "PYL" | "ASX" | "GLX" | "XLE" | "XAA" => true
    _ => false
  }
}

///|
/// Standard nucleic acid residue names.
fn is_nucleic_resname(resname : String) -> Bool {
  match resname {
    "DA" | "DC" | "DG" | "DT" | "DU" |
    "A" | "C" | "G" | "T" | "U" |
    "ADE" | "CYT" | "GUA" | "THY" | "URA" |
    "PSU" | "H2U" | "5MU" | "OMU" | "1MA" | "2MA" |
    "MIA" | "OMG" | "YYG" | "7MG" => true
    _ => false
  }
}

// ---------------------------------------------------------------------------
// Residue-level classification
// ---------------------------------------------------------------------------

///|
/// Classify a single residue into Protein, Nucleic, Water, Ion, Ligand, or Other.
///
/// Decision order:
///   1. Water resname -> Water
///   2. Ion resname (or single-atom hetero with metal/halogen element) -> Ion
///   3. Standard amino acid -> Protein
///   4. Nucleic acid resname -> Nucleic
///   5. Hetero residue (not water/ion/AA/nucleic) -> Ligand
///   6. Everything else -> Other
pub fn classify_residue(residue : Residue) -> ResidueClass {
  let rn = residue.resname
  let n_atoms = residue.atoms.length()

  // 1. Water
  if is_water_resname(rn) {
    return Water
  }

  // 2. Ion — known ion resname
  if is_ion_resname(rn) {
    return Ion
  }

  // 2b. Ion — single-atom hetero with metal/halogen element
  if n_atoms == 1 && residue.is_het() {
    let atom = residue.atoms[0]
    if atom.element.length() > 0 {
      if is_metal_element(atom.element) || is_halogen_element(atom.element) {
        return Ion
      }
    }
    // Fallback: resname itself is a 1-2 letter element-like code
    let upper = to_uppercase(rn)
    if rn.length() <= 2 && is_metal_ion(upper) {
      return Ion
    }
  }

  // 3. Standard amino acid
  if is_standard_aa(rn) {
    return Protein
  }

  // 4. Nucleic acid
  if is_nucleic_resname(rn) {
    return Nucleic
  }

  // 5. Hetero residue -> Ligand
  if residue.is_het() {
    return Ligand
  }

  // 6. Non-het, non-standard -> could be modified residue or unknown
  // If it has a CA atom and multiple atoms, treat as protein-like
  if n_atoms >= 2 {
    let has_ca = match residue.find_atom("CA") {
      Some(_) => true
      None => false
    }
    if has_ca {
      return Protein
    }
  }

  Other
}

// ---------------------------------------------------------------------------
// Chain-level extraction
// ---------------------------------------------------------------------------

///|
/// Extract all water residues from a chain.
pub fn extract_waters(chain : Chain) -> Array[Residue] {
  let result = Array(capacity=0)
  for res in chain.residues {
    if classify_residue(res) is Water {
      result.push(res)
    }
  }
  result
}

///|
/// Extract all ion residues from a chain.
pub fn extract_ions(chain : Chain) -> Array[Residue] {
  let result = Array(capacity=0)
  for res in chain.residues {
    if classify_residue(res) is Ion {
      result.push(res)
    }
  }
  result
}

///|
/// Extract all ligand residues from a chain.
pub fn extract_ligands(chain : Chain) -> Array[Residue] {
  let result = Array(capacity=0)
  for res in chain.residues {
    if classify_residue(res) is Ligand {
      result.push(res)
    }
  }
  result
}

///|
/// Extract all protein residues from a chain.
pub fn extract_protein(chain : Chain) -> Array[Residue] {
  let result = Array(capacity=0)
  for res in chain.residues {
    if classify_residue(res) is Protein {
      result.push(res)
    }
  }
  result
}

///|
/// Extract all nucleic acid residues from a chain.
pub fn extract_nucleic(chain : Chain) -> Array[Residue] {
  let result = Array(capacity=0)
  for res in chain.residues {
    if classify_residue(res) is Nucleic {
      result.push(res)
    }
  }
  result
}

// ---------------------------------------------------------------------------
// Chain-building helpers
// ---------------------------------------------------------------------------

///|
/// Build a water-only chain.
pub fn water_chain(chain : Chain) -> Chain {
  Chain::new(id=chain.id, residues=extract_waters(chain))
}

///|
/// Build an ion-only chain.
pub fn ion_chain(chain : Chain) -> Chain {
  Chain::new(id=chain.id, residues=extract_ions(chain))
}

///|
/// Build a ligand-only chain.
pub fn ligand_chain(chain : Chain) -> Chain {
  Chain::new(id=chain.id, residues=extract_ligands(chain))
}

///|
/// Build a protein-only chain.
pub fn protein_chain(chain : Chain) -> Chain {
  Chain::new(id=chain.id, residues=extract_protein(chain))
}

///|
/// Build a nucleic-acid-only chain.
pub fn nucleic_chain(chain : Chain) -> Chain {
  Chain::new(id=chain.id, residues=extract_nucleic(chain))
}

// ---------------------------------------------------------------------------
// Full separation
// ---------------------------------------------------------------------------

///|
/// Separate all residues in a chain into component groups.
pub fn separate_chain(chain : Chain) -> SeparationResult {
  let protein = Array(capacity=0)
  let nucleic = Array(capacity=0)
  let ligands = Array(capacity=0)
  let waters = Array(capacity=0)
  let ions = Array(capacity=0)
  let other = Array(capacity=0)
  for res in chain.residues {
    match classify_residue(res) {
      Protein => protein.push(res)
      Nucleic => nucleic.push(res)
      Water => waters.push(res)
      Ion => ions.push(res)
      Ligand => ligands.push(res)
      Other => other.push(res)
    }
  }
  SeparationResult::{ protein, nucleic, ligands, waters, ions, other }
}

///|
/// Separate all residues across all chains in a model.
pub fn separate_model(model : Model) -> SeparationResult {
  let protein = Array(capacity=0)
  let nucleic = Array(capacity=0)
  let ligands = Array(capacity=0)
  let waters = Array(capacity=0)
  let ions = Array(capacity=0)
  let other = Array(capacity=0)
  for chain in model.chains {
    let r = separate_chain(chain)
    for res in r.protein {
      protein.push(res)
    }
    for res in r.nucleic {
      nucleic.push(res)
    }
    for res in r.ligands {
      ligands.push(res)
    }
    for res in r.waters {
      waters.push(res)
    }
    for res in r.ions {
      ions.push(res)
    }
    for res in r.other {
      other.push(res)
    }
  }
  SeparationResult::{ protein, nucleic, ligands, waters, ions, other }
}

// ---------------------------------------------------------------------------
// SeparationResult accessors
// ---------------------------------------------------------------------------

///|
/// Get protein residues.
pub fn SeparationResult::get_protein(self : SeparationResult) -> Array[Residue] {
  self.protein
}

///|
/// Get nucleic acid residues.
pub fn SeparationResult::get_nucleic(self : SeparationResult) -> Array[Residue] {
  self.nucleic
}

///|
/// Get ligand residues.
pub fn SeparationResult::get_ligands(self : SeparationResult) -> Array[Residue] {
  self.ligands
}

///|
/// Get water residues.
pub fn SeparationResult::get_waters(self : SeparationResult) -> Array[Residue] {
  self.waters
}

///|
/// Get ion residues.
pub fn SeparationResult::get_ions(self : SeparationResult) -> Array[Residue] {
  self.ions
}

///|
/// Get other residues.
pub fn SeparationResult::get_other(self : SeparationResult) -> Array[Residue] {
  self.other
}

///|
/// Total number of residues.
pub fn SeparationResult::total_residues(self : SeparationResult) -> Int {
  self.protein.length() + self.nucleic.length() + self.ligands.length() +
  self.waters.length() + self.ions.length() + self.other.length()
}

///|
/// Generate a summary of the separation.
pub fn SeparationResult::summary(self : SeparationResult) -> SeparationSummary {
  let mut num_ligand_atoms = 0
  for res in self.ligands {
    num_ligand_atoms = num_ligand_atoms + res.atoms.length()
  }
  let mut num_water_atoms = 0
  for res in self.waters {
    num_water_atoms = num_water_atoms + res.atoms.length()
  }
  let mut num_ion_atoms = 0
  for res in self.ions {
    num_ion_atoms = num_ion_atoms + res.atoms.length()
  }
  SeparationSummary::{
    num_protein: self.protein.length(),
    num_nucleic: self.nucleic.length(),
    num_ligands: self.ligands.length(),
    num_waters: self.waters.length(),
    num_ions: self.ions.length(),
    num_other: self.other.length(),
    num_ligand_atoms,
    num_water_atoms,
    num_ion_atoms,
  }
}

///|
/// Pretty-print a summary as a human-readable string.
pub fn SeparationSummary::to_string(self : SeparationSummary) -> String {
  let sb = StringBuilder::new()
  sb.write_string("Component Separation Summary:\n")
  sb.write_string("  Protein residues: ")
  sb.write_string(self.num_protein.to_string())
  sb.write_char('\n')
  sb.write_string("  Nucleic residues: ")
  sb.write_string(self.num_nucleic.to_string())
  sb.write_char('\n')
  sb.write_string("  Ligand residues:  ")
  sb.write_string(self.num_ligands.to_string())
  sb.write_string(" (")
  sb.write_string(self.num_ligand_atoms.to_string())
  sb.write_string(" atoms)\n")
  sb.write_string("  Water residues:   ")
  sb.write_string(self.num_waters.to_string())
  sb.write_string(" (")
  sb.write_string(self.num_water_atoms.to_string())
  sb.write_string(" atoms)\n")
  sb.write_string("  Ion residues:     ")
  sb.write_string(self.num_ions.to_string())
  sb.write_string(" (")
  sb.write_string(self.num_ion_atoms.to_string())
  sb.write_string(" atoms)\n")
  sb.write_string("  Other residues:   ")
  sb.write_string(self.num_other.to_string())
  sb.write_char('\n')
  sb.to_string()
}

// ---------------------------------------------------------------------------
// Ligand grouping — group contiguous hetero residues into distinct ligands
// ---------------------------------------------------------------------------

///|
/// A single ligand molecule: one or more residues that form a connected unit.
pub(all) struct LigandMolecule {
  resname : String
  chainid : Char
  resseqs : Array[Int]
  residues : Array[Residue]
  num_atoms : Int
}

///|
/// Group ligand residues into distinct ligand molecules.
/// Residues with the same resname and consecutive resseq are grouped together.
pub fn group_ligands(ligand_residues : Array[Residue]) -> Array[LigandMolecule] {
  let result = Array(capacity=0)
  let n = ligand_residues.length()
  if n == 0 {
    return result
  }
  let used = Array::make(n, false)
  for i in 0..<n {
    if used[i] {
      continue
    }
    used[i] = true
    let group = Array(capacity=1)
    group.push(ligand_residues[i])
    let rn = ligand_residues[i].resname
    let cid = ligand_residues[i].chainid
    // Greedily find connected residues (same resname, same chain, consecutive resseq)
    let mut changed = true
    while changed {
      changed = false
      for j in 0..<n {
        if used[j] {
          continue
        }
        if ligand_residues[j].resname != rn || ligand_residues[j].chainid != cid {
          continue
        }
        // Check if consecutive resseq with any in current group
        for k in 0..<group.length() {
          let diff = ligand_residues[j].resseq - group[k].resseq
          if diff == 1 || diff == -1 {
            group.push(ligand_residues[j])
            used[j] = true
            changed = true
            break
          }
        }
      }
    }
    // Build LigandMolecule
    let resseqs = Array(capacity=group.length())
    let mut num_atoms = 0
    for res in group {
      resseqs.push(res.resseq)
      num_atoms = num_atoms + res.atoms.length()
    }
    result.push(LigandMolecule::{
      resname: rn,
      chainid: cid,
      resseqs,
      residues: group,
      num_atoms,
    })
  }
  result
}

///|
/// Get the center of mass of a ligand molecule.
pub fn LigandMolecule::center_of_mass(self : LigandMolecule) -> Vector3 {
  let mut x = 0.0
  let mut y = 0.0
  let mut z = 0.0
  let mut count = 0
  for res in self.residues {
    for atom in res.atoms {
      x = x + atom.coord.x
      y = y + atom.coord.y
      z = z + atom.coord.z
      count = count + 1
    }
  }
  if count == 0 {
    return Vector3::new(0.0, 0.0, 0.0)
  }
  let d = count.to_double()
  Vector3::new(x / d, y / d, z / d)
}

///|
/// Get all atoms in a ligand molecule.
pub fn LigandMolecule::get_atoms(self : LigandMolecule) -> Array[Atom] {
  let result = Array(capacity=0)
  for res in self.residues {
    for atom in res.atoms {
      result.push(atom)
    }
  }
  result
}

///|
/// Count the number of atoms in a ligand molecule.
pub fn LigandMolecule::num_atoms_(self : LigandMolecule) -> Int {
  self.num_atoms
}

// ---------------------------------------------------------------------------
// Ligand-protein contact detection
// ---------------------------------------------------------------------------

///|
/// Find ligand residues that are within a distance cutoff of any protein atom.
/// Returns the indices of contacting ligands.
pub fn contacting_ligands(
  ligands : Array[LigandMolecule],
  protein_atoms : Array[Atom],
  cutoff : Double,
) -> Array[Int] {
  let cutoff_sq = cutoff * cutoff
  let result = Array(capacity=0)
  for i in 0..<ligands.length() {
    let lig_atoms = ligands[i].get_atoms()
    let mut found = false
    for j in 0..<lig_atoms.length() {
      if found {
        break
      }
      for k in 0..<protein_atoms.length() {
        let dx = lig_atoms[j].coord.x - protein_atoms[k].coord.x
        let dy = lig_atoms[j].coord.y - protein_atoms[k].coord.y
        let dz = lig_atoms[j].coord.z - protein_atoms[k].coord.z
        if dx * dx + dy * dy + dz * dz <= cutoff_sq {
          found = true
          break
        }
      }
    }
    if found {
      result.push(i)
    }
  }
  result
}

///|
/// Check if an ion is coordinated by protein atoms (within cutoff).
/// Returns the list of coordinating atom indices.
pub fn ion_coordination(
  ion : Residue,
  protein_atoms : Array[Atom],
  cutoff : Double,
) -> Array[Int] {
  let cutoff_sq = cutoff * cutoff
  let result = Array(capacity=0)
  if ion.atoms.length() == 0 {
    return result
  }
  let ion_atom = ion.atoms[0]
  for k in 0..<protein_atoms.length() {
    let dx = ion_atom.coord.x - protein_atoms[k].coord.x
    let dy = ion_atom.coord.y - protein_atoms[k].coord.y
    let dz = ion_atom.coord.z - protein_atoms[k].coord.z
    if dx * dx + dy * dy + dz * dz <= cutoff_sq {
      result.push(k)
    }
  }
  result
}
