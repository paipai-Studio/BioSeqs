// ============================================================================
// Molecular descriptors (RDKit/Open Babel compatible)
// ============================================================================
// Implements: LogP (Wildman-Crippen), TPSA (Ertl), HBD/HBA, rotatable bonds,
// formal charge, fraction sp3, aromatic/saturated ring count, valence electrons.
// All operate on the Molecule type from rdkit.mbt.

///|
/// Atomic contribution to LogP (Wildman-Crippen method).
/// Simplified table: (symbol, degree, aromatic) -> contribution.
fn logp_atom_contribution(
  symbol : String,
  degree : Int,
  aromatic : Bool,
  in_ring : Bool,
) -> Double {
  if aromatic {
    match symbol {
      "C" => 0.13
      "N" => -0.53
      "O" => -0.48
      "S" => 0.00
      _ => 0.00
    }
  } else {
    match symbol {
      "C" =>
        if degree == 4 { 0.20 } else if degree == 3 { 0.14 } else { 0.10 }
      "N" =>
        if degree == 3 { -1.00 } else if degree == 2 { -0.50 } else { -1.50 }
      "O" =>
        if degree == 2 { -0.45 } else { -0.60 }
      "S" => -0.20
      "P" => 0.25
      "F" => -0.25
      "Cl" => 0.06
      "Br" => 0.20
      "I" => 0.50
      "H" => 0.00
      _ => 0.00
    }
  }
}

///|
/// Calculate Wildman-Crippen LogP (octanol-water partition coefficient).
/// Sum of atom contributions + correction for ring atoms.
pub fn mol_logp(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut logp = 0.0
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    let deg = neighbors[k].length()
    let in_ring = mol_is_in_ring(mol, k)
    let contrib = logp_atom_contribution(a.symbol, deg, a.aromatic, in_ring)
    logp = logp + contrib
    // Ring correction
    if in_ring && not(a.aromatic) {
      logp = logp - 0.05
    }
  }
  // Bond corrections
  for b in mol.bonds {
    if b.order == 2 {
      logp = logp - 0.42
    } else if b.order == 3 {
      logp = logp - 0.70
    }
  }
  logp
}

///|
/// Atomic contribution to TPSA (topological polar surface area, Ertl method).
/// Based on atom type (N/O with specific neighbor patterns).
fn tpsa_atom_contribution(
  symbol : String,
  degree : Int,
  num_h_bonds : Int,
  in_ring : Bool,
) -> Double {
  // N types
  if symbol == "N" {
    if in_ring && degree == 2 {
      return 3.24 // aromatic N in ring
    }
    if degree == 1 {
      return 23.79 // primary amine
    }
    if degree == 2 {
      return 12.36 // secondary amine
    }
    if degree == 3 {
      return 3.24 // tertiary amine
    }
    return 0.0
  }
  // O types
  if symbol == "O" {
    if degree == 1 {
      return 20.23 // hydroxyl/ether
    }
    if degree == 2 {
      return 9.23 // ring O
    }
    return 0.0
  }
  // S and P contribute 0 in basic TPSA
  0.0
}

///|
/// Calculate topological polar surface area (TPSA) in Angstrom^2.
pub fn mol_tpsa(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut tpsa = 0.0
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    let deg = neighbors[k].length()
    let in_ring = mol_is_in_ring(mol, k)
    tpsa = tpsa + tpsa_atom_contribution(a.symbol, deg, 0, in_ring)
  }
  tpsa
}

///|
/// Count H-bond donors (N-H or O-H bonds).
/// An HBD is any N or O atom with at least one implicit H.
pub fn mol_hbd_count(mol : Molecule) -> Int {
  let n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut count = 0
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    if a.symbol == "N" || a.symbol == "O" {
      let deg = neighbors[k].length()
      // Implicit H count = typical_valence - degree - charge
      let typical = if a.symbol == "N" { 3 } else { 2 }
      let implicit_h = typical - deg + a.charge
      if implicit_h > 0 {
        count = count + implicit_h
      }
    }
  }
  count
}

///|
/// Count H-bond acceptors (N or O atoms).
pub fn mol_hba_count(mol : Molecule) -> Int {
  let mut count = 0
  for a in mol.atoms {
    if a.symbol == "N" || a.symbol == "O" {
      count = count + 1
    }
  }
  count
}

///|
/// Count rotatable bonds.
/// A rotatable bond is a single bond, not in a ring, not to a terminal H,
/// and both endpoints have at least 2 heavy neighbors.
pub fn mol_rotatable_bonds(mol : Molecule) -> Int {
  let n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut count = 0
  for b in mol.bonds {
    if b.order != 1 {
      continue
    }
    // Skip if either atom is H
    let a1 = mol.atoms[b.from]
    let a2 = mol.atoms[b.to]
    if a1.symbol == "H" || a2.symbol == "H" {
      continue
    }
    // Skip if in ring
    let rings = mol_find_rings(mol)
    let mut in_ring = false
    for ring in rings {
      let mut has_from = false
      let mut has_to = false
      for a in ring {
        if a == b.from {
          has_from = true
        }
        if a == b.to {
          has_to = true
        }
      }
      if has_from && has_to {
        in_ring = true
      }
    }
    if in_ring {
      continue
    }
    // Both endpoints need >= 2 neighbors
    if neighbors[b.from].length() >= 2 && neighbors[b.to].length() >= 2 {
      count = count + 1
    }
  }
  count
}

///|
/// Calculate formal charge (sum of all atom charges).
pub fn mol_formal_charge(mol : Molecule) -> Int {
  let mut charge = 0
  for a in mol.atoms {
    charge = charge + a.charge
  }
  charge
}

///|
/// Calculate fraction of sp3 carbons (C with 4 single bonds / total C).
pub fn mol_fraction_sp3(mol : Molecule) -> Double {
  let n = mol.num_atoms()
  let neighbors = mol_adjacency(mol)
  let mut total_c = 0
  let mut sp3_c = 0
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    if a.symbol == "C" {
      total_c = total_c + 1
      // Check if all bonds are single
      let mut all_single = true
      for (_, order) in neighbors[k] {
        if order != 1 {
          all_single = false
        }
      }
      if all_single && not(a.aromatic) {
        sp3_c = sp3_c + 1
      }
    }
  }
  if total_c == 0 {
    return 0.0
  }
  sp3_c.to_double() / total_c.to_double()
}

///|
/// Count valence electrons of the molecule.
pub fn mol_valence_electrons(mol : Molecule) -> Int {
  let mut total = 0
  for a in mol.atoms {
    total = total + atom_valence_electrons(a.symbol)
  }
  total
}

///|
/// Valence electrons for an element.
fn atom_valence_electrons(symbol : String) -> Int {
  match symbol {
    "H" => 1
    "C" => 4
    "N" => 5
    "O" => 6
    "F" => 7
    "P" => 5
    "S" => 6
    "Cl" => 7
    "Br" => 7
    "I" => 7
    "B" => 3
    "Si" => 4
    _ => 4
  }
}

///|
/// Count heavy atoms (non-H).
pub fn mol_heavy_atom_count(mol : Molecule) -> Int {
  let mut count = 0
  for a in mol.atoms {
    if a.symbol != "H" {
      count = count + 1
    }
  }
  count
}

///|
/// Molecular weight from atom symbols.
pub fn mol_molecular_weight(mol : Molecule) -> Double {
  let mut mw = 0.0
  for a in mol.atoms {
    mw = mw + atom_atomic_weight(a.symbol)
  }
  mw
}

///|
/// Atomic weight for common elements.
fn atom_atomic_weight(symbol : String) -> Double {
  match symbol {
    "H" => 1.008
    "C" => 12.011
    "N" => 14.007
    "O" => 15.999
    "F" => 18.998
    "P" => 30.974
    "S" => 32.065
    "Cl" => 35.453
    "Br" => 79.904
    "I" => 126.904
    "B" => 10.811
    "Si" => 28.086
    _ => 12.0
  }
}

///|
/// Count total atoms.
pub fn mol_num_atoms(mol : Molecule) -> Int {
  mol.num_atoms()
}

///|
/// Count total bonds.
pub fn mol_num_bonds(mol : Molecule) -> Int {
  mol.num_bonds()
}

///|
/// Full descriptor report as a formatted string.
pub fn mol_descriptor_report(mol : Molecule) -> String {
  let sb = StringBuilder()
  sb.write_string("Molecular Descriptors\n")
  sb.write_string("=====================\n")
  sb.write_string("Atoms: " + mol_num_atoms(mol).to_string() + "\n")
  sb.write_string("Heavy atoms: " + mol_heavy_atom_count(mol).to_string() + "\n")
  sb.write_string("Bonds: " + mol_num_bonds(mol).to_string() + "\n")
  sb.write_string("Molecular weight: " + mol_molecular_weight(mol).to_string() + "\n")
  sb.write_string("LogP: " + mol_logp(mol).to_string() + "\n")
  sb.write_string("TPSA: " + mol_tpsa(mol).to_string() + "\n")
  sb.write_string("HBD: " + mol_hbd_count(mol).to_string() + "\n")
  sb.write_string("HBA: " + mol_hba_count(mol).to_string() + "\n")
  sb.write_string("Rotatable bonds: " + mol_rotatable_bonds(mol).to_string() + "\n")
  sb.write_string("Formal charge: " + mol_formal_charge(mol).to_string() + "\n")
  sb.write_string("Fraction sp3: " + mol_fraction_sp3(mol).to_string() + "\n")
  sb.write_string("Rings: " + mol_ring_count(mol).to_string() + "\n")
  sb.write_string("Aromatic rings: " + mol_aromatic_ring_count(mol).to_string() + "\n")
  sb.write_string("Saturated rings: " + mol_saturated_ring_count(mol).to_string() + "\n")
  sb.write_string("Valence electrons: " + mol_valence_electrons(mol).to_string() + "\n")
  sb.to_string()
}
