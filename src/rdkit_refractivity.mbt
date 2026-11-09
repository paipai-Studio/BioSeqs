/// Molar refractivity (MR) descriptor.
///
/// RDKit equivalent: `rdkit.Chem.Crippen.MolMR`
///
/// Molar refractivity is calculated from the Lorentz-Lorenz equation
/// using atomic contributions (Crippen fragment method).
/// MR = sum of atomic refractivity contributions

///|
/// Atomic refractivity contributions (simplified Wildman-Crippen table)
fn atomic_mr(symbol : String, aromatic : Bool) -> Double {
  match symbol {
    "H" => 1.1
    "C" => if aromatic { 2.78 } else { 3.42 }
    "N" => if aromatic { 2.86 } else { 2.81 }
    "O" => if aromatic { 1.64 } else { 1.54 }
    "F" => 1.05
    "P" => 7.84
    "S" => if aromatic { 6.36 } else { 6.72 }
    "Cl" => 6.03
    "Br" => 8.88
    "I" => 13.94
    "B" => 3.86
    "Si" => 5.62
    "Se" => 8.46
    "As" => 9.34
    _ => 2.0
  }
}

///|
/// Compute molar refractivity (MR) for a molecule.
/// Uses atomic contribution method (Crippen-style).
pub fn rdkit_molar_refractivity(mol : Molecule) -> Double {
  let mut mr = 0.0
  for atom in mol.atoms {
    mr = mr + atomic_mr(atom.symbol, atom.aromatic)
  }
  // Bond corrections: each bond adds a small contribution
  for bond in mol.bonds {
    match bond.order {
      2 => mr = mr + 1.5 // double bond correction
      3 => mr = mr + 3.0 // triple bond correction
      _ => ()
    }
  }
  mr
}

///|
/// Compute polarizability estimate from molar refractivity.
/// alpha = MR / (4 * pi * N_A / 3) ≈ MR / 4.18 (in Angstrom^3)
pub fn rdkit_polarizability(mol : Molecule) -> Double {
  rdkit_molar_refractivity(mol) / 4.18
}

///|
/// Compute molar volume estimate from MR and LogP.
/// MV = MR / (n^2 - 1) where n is refractive index ≈ 1.5
/// Simplified: MV ≈ MR * 0.8
pub fn rdkit_molar_volume(mol : Molecule) -> Double {
  rdkit_molar_refractivity(mol) * 0.8
}

///|
/// Compute parachor (surface tension related) from MR.
/// Parachor ≈ MR * 40 (rough estimate)
pub fn rdkit_parachor(mol : Molecule) -> Double {
  rdkit_molar_refractivity(mol) * 40.0
}
