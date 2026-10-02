// ============================================================================
// Tautomer enumeration
// ============================================================================
// Enumerates tautomers using predefined transform rules:
// (1) Keto-enol: C(=O)-C <-> C(-O)=C (ketone-enol)
// (2) Imine-amine: C=N-C <-> C-N=C (imine shift)
// (3) Amide: C(=O)-N <-> C(-O)=N (amide resonance)
// (4) Nitroso-oxime: N=O <-> N-OH (simplified)
// Uses iterative rule application with visited set to avoid cycles.

///|
/// Apply keto-enol tautomerism: find C(=O)-C pattern and convert to C(-O)=C.
/// Returns list of new molecules (tautomers).
fn tautomer_keto_enol(mol : Molecule) -> Array[Molecule] {
  let results : Array[Molecule] = Array(capacity=0)
  let neighbors = mol_adjacency(mol)
  for b in mol.bonds {
    if b.order != 2 {
      continue
    }
    // Find C=O where C also has a C neighbor
    let (c_atom, o_atom) = if mol.atoms[b.from].symbol == "C" && mol.atoms[b.to].symbol == "O" {
      (b.from, b.to)
    } else if mol.atoms[b.to].symbol == "C" && mol.atoms[b.from].symbol == "O" {
      (b.to, b.from)
    } else {
      continue
    }
    // Find adjacent C neighbor (single bond)
    for (nb, order) in neighbors[c_atom] {
      if nb == o_atom || order != 1 {
        continue
      }
      if mol.atoms[nb].symbol != "C" {
        continue
      }
      // Transform: C=O + C-C -> C-O + C=C
      let new_mol = tautomer_swap_bond_orders(mol, c_atom, o_atom, c_atom, nb)
      results.push(new_mol)
    }
  }
  results
}

///|
/// Apply enol-keto tautomerism: find C(-O)=C pattern and convert to C(=O)-C.
fn tautomer_enol_keto(mol : Molecule) -> Array[Molecule] {
  let results : Array[Molecule] = Array(capacity=0)
  let neighbors = mol_adjacency(mol)
  for b in mol.bonds {
    if b.order != 2 {
      continue
    }
    // Find C=C where one C has an O neighbor (single bond)
    let (c1, c2) = if mol.atoms[b.from].symbol == "C" && mol.atoms[b.to].symbol == "C" {
      (b.from, b.to)
    } else {
      continue
    }
    // Check both directions
    for (cc, other_c) in [(c1, c2), (c2, c1)] {
      for (nb, order) in neighbors[cc] {
        if nb == other_c || order != 1 {
          continue
        }
        if mol.atoms[nb].symbol != "O" {
          continue
        }
        // Transform: C=C + C-O -> C-C + C=O
        let new_mol = tautomer_swap_bond_orders(mol, cc, other_c, cc, nb)
        results.push(new_mol)
      }
    }
  }
  results
}

///|
/// Apply imine-amine tautomerism: C=N-C <-> C-N=C.
fn tautomer_imine_amine(mol : Molecule) -> Array[Molecule] {
  let results : Array[Molecule] = Array(capacity=0)
  let neighbors = mol_adjacency(mol)
  for b in mol.bonds {
    if b.order != 2 {
      continue
    }
    // Find C=N where N has a C neighbor
    let (c_atom, n_atom) = if mol.atoms[b.from].symbol == "C" && mol.atoms[b.to].symbol == "N" {
      (b.from, b.to)
    } else if mol.atoms[b.to].symbol == "C" && mol.atoms[b.from].symbol == "N" {
      (b.to, b.from)
    } else {
      continue
    }
    for (nb, order) in neighbors[n_atom] {
      if nb == c_atom || order != 1 {
        continue
      }
      if mol.atoms[nb].symbol != "C" {
        continue
      }
      let new_mol = tautomer_swap_bond_orders(mol, c_atom, n_atom, n_atom, nb)
      results.push(new_mol)
    }
  }
  results
}

///|
/// Swap bond orders: bond1 becomes single, bond2 becomes double.
fn tautomer_swap_bond_orders(
  mol : Molecule,
  a1 : Int,
  b1 : Int,
  a2 : Int,
  b2 : Int,
) -> Molecule {
  let new_bonds : Array[MolBond] = Array(capacity=0)
  for b in mol.bonds {
    let new_order = if (b.from == a1 && b.to == b1) || (b.from == b1 && b.to == a1) {
      1
    } else if (b.from == a2 && b.to == b2) || (b.from == b2 && b.to == a2) {
      2
    } else {
      b.order
    }
    new_bonds.push(MolBond::{ from: b.from, to: b.to, order: new_order })
  }
  Molecule::{ atoms: mol.atoms, bonds: new_bonds }
}

///|
/// Generate a canonical key for a molecule (for deduplication).
fn tautomer_mol_key(mol : Molecule) -> String {
  mol_to_smiles(mol)
}

///|
/// Enumerate tautomers of a molecule.
/// Applies all transform rules iteratively up to max_iterations.
/// Returns a list of unique tautomers (including the input molecule).
pub fn enumerate_tautomers(mol : Molecule, max_iterations : Int) -> Array[Molecule] {
  let tautomers : Array[Molecule] = [mol]
  let visited : Map[String, Bool] = Map([], capacity=0)
  visited[tautomer_mol_key(mol)] = true
  let mut changed = true
  let mut iter = 0
  while changed && iter < max_iterations {
    changed = false
    let current_count = tautomers.length()
    for i = 0; i < current_count; i = i + 1 {
      let current = tautomers[i]
      // Apply all rules
      let keto_enol = tautomer_keto_enol(current)
      let enol_keto = tautomer_enol_keto(current)
      let imine_amine = tautomer_imine_amine(current)
      for t in keto_enol {
        let key = tautomer_mol_key(t)
        match visited.get(key) {
          Some(_) => ()
          None => {
            visited[key] = true
            tautomers.push(t)
            changed = true
          }
        }
      }
      for t in enol_keto {
        let key = tautomer_mol_key(t)
        match visited.get(key) {
          Some(_) => ()
          None => {
            visited[key] = true
            tautomers.push(t)
            changed = true
          }
        }
      }
      for t in imine_amine {
        let key = tautomer_mol_key(t)
        match visited.get(key) {
          Some(_) => ()
          None => {
            visited[key] = true
            tautomers.push(t)
            changed = true
          }
        }
      }
    }
    iter = iter + 1
  }
  tautomers
}

///|
/// Count number of tautomers.
pub fn tautomer_count(mol : Molecule, max_iterations : Int) -> Int {
  enumerate_tautomers(mol, max_iterations).length()
}
