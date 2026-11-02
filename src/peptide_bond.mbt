///|
/// Peptide bond detection and cis/trans classification.
///
/// Biotite provides `peptide_bond()` and related functions to identify
/// peptide bonds and classify them as cis or trans based on the omega
/// dihedral angle.

///|
/// Peptide bond geometry classification.
pub enum PeptideBondType {
  Trans // omega ~ 180° (most common, ~99.7%)
  Cis // omega ~ 0° (rare, ~0.3%)
  CisPro // cis with proline (more common, ~5% of prolines)
}

///|
/// Get bond type as string.
pub fn PeptideBondType::to_string(self : PeptideBondType) -> String {
  match self {
    Trans => "trans"
    Cis => "cis"
    CisPro => "cis-pro"
  }
}

///|
/// A peptide bond between two residues.
pub struct PeptideBond {
  residue1 : Int32 // index of first residue
  residue2 : Int32 // index of second residue
  omega : Double // dihedral angle in degrees
  bond_type : PeptideBondType
  length : Double // C-N bond length in Angstroms
}

///|
/// Get the omega angle.
pub fn PeptideBond::get_omega(self : PeptideBond) -> Double {
  self.omega
}

///|
/// Get the bond type.
pub fn PeptideBond::get_bond_type(self : PeptideBond) -> PeptideBondType {
  self.bond_type
}

///|
/// Get the C-N bond length.
pub fn PeptideBond::get_length(self : PeptideBond) -> Double {
  self.length
}

///|
/// Classify a peptide bond as cis, trans, or cis-proline.
///
/// The omega angle is the dihedral angle C(i-1)-N(i)-CA(i)-C(i).
/// - |omega| > 90°: trans
/// - |omega| <= 90°: cis (or cis-pro if residue i is proline)
pub fn classify_peptide_bond(
  omega : Double,
  is_proline : Bool,
) -> PeptideBondType {
  let abs_omega = omega.abs()
  // Normalize to [-180, 180]
  let mut a = abs_omega
  while a > 180.0 {
    a = a - 360.0
  }
  a = a.abs()
  if a > 90.0 {
    Trans
  } else if is_proline {
    CisPro
  } else {
    Cis
  }
}

///|
/// Calculate the omega dihedral angle for a peptide bond.
///
/// Omega = dihedral(C_prev, N_curr, CA_curr, C_curr)
/// Requires the C atom of the previous residue and N, CA, C of current.
pub fn calculate_omega(
  c_prev : Vector3,
  n_curr : Vector3,
  ca_curr : Vector3,
  c_curr : Vector3,
) -> Double {
  dihedral_angle(c_prev, n_curr, ca_curr, c_curr)
}

///|
/// Detect all peptide bonds in a chain.
///
/// For each pair of consecutive residues, checks if a peptide bond
/// exists (C-N distance < 1.4 Å) and classifies it.
pub fn detect_peptide_bonds(residues : Array[Residue]) -> Array[PeptideBond] {
  let bonds : Array[PeptideBond] = []
  if residues.length() < 2 {
    return bonds
  }
  for i in 1..<residues.length() {
    let prev = residues[i - 1]
    let curr = residues[i]
    // Find C atom in previous residue
    let mut c_prev : Vector3? = None
    for a in prev.get_atoms() {
      if a.name == "C" {
        c_prev = Some(a.get_coord())
        break
      }
    }
    // Find N, CA, C atoms in current residue
    let mut n_curr : Vector3? = None
    let mut ca_curr : Vector3? = None
    let mut c_curr : Vector3? = None
    for a in curr.get_atoms() {
      match a.name {
        "N" => n_curr = Some(a.get_coord())
        "CA" => ca_curr = Some(a.get_coord())
        "C" => c_curr = Some(a.get_coord())
        _ => ()
      }
    }
    match (c_prev, n_curr, ca_curr, c_curr) {
      (Some(cp), Some(n), Some(ca), Some(c)) => {
        let bond_length = cp.distance(n)
        // Peptide bond if C-N < 1.4 Å
        if bond_length < 1.4 {
          let omega = calculate_omega(cp, n, ca, c)
          let is_pro = curr.resname == "PRO"
          let bt = classify_peptide_bond(omega, is_pro)
          bonds.push(PeptideBond::{
            residue1: i - 1,
            residue2: i,
            omega,
            bond_type: bt,
            length: bond_length,
          })
        }
      }
      _ => ()
    }
  }
  bonds
}

///|
/// Count cis peptide bonds.
pub fn count_cis_bonds(bonds : Array[PeptideBond]) -> Int32 {
  let mut count = 0
  for b in bonds {
    match b.bond_type {
      Cis => count = count + 1
      CisPro => count = count + 1
      _ => ()
    }
  }
  count
}

///|
/// Find all cis-proline bonds.
pub fn find_cis_prolines(bonds : Array[PeptideBond]) -> Array[Int32] {
  let result : Array[Int32] = []
  for b in bonds {
    match b.bond_type {
      CisPro => result.push(b.residue2)
      _ => ()
    }
  }
  result
}

///|
/// Check if a peptide bond is unusual (cis non-proline).
pub fn is_unusual_bond(bond : PeptideBond) -> Bool {
  match bond.bond_type {
    Cis => true // cis non-proline is very rare
    _ => false
  }
}
