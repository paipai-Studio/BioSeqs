// ============================================================================
// Murcko scaffold extraction
// ============================================================================
// The Murcko scaffold is the ring system + linkers connecting rings.
// Algorithm: (1) iteratively prune terminal (degree-1) atoms to remove side chains,
// (2) the remaining atoms form the scaffold (rings + linkers between rings).
// Reference: Bemis & Murcko, J. Med. Chem. 1996, 39, 2887.

///|
/// Result of Murcko scaffold extraction.
pub struct MurckoResult {
  scaffold_atoms : Array[Int]    // indices of scaffold atoms in original molecule
  side_chain_atoms : Array[Int]  // indices of side chain atoms
  num_rings : Int                // number of ring systems in scaffold
}

///|
/// Extract Murcko scaffold from a molecule.
/// Returns the set of scaffold atom indices and side chain atom indices.
pub fn murcko_scaffold(mol : Molecule) -> MurckoResult {
  let n = mol.atoms.length()
  // Step 1: Iteratively prune terminal atoms (degree 1)
  let is_scaffold = Array::make(n, true)
  let neighbors = mol_adjacency(mol)
  let mut changed = true
  while changed {
    changed = false
    for k = 0; k < n; k = k + 1 {
      if not(is_scaffold[k]) {
        continue
      }
      // Count scaffold neighbors
      let mut deg = 0
      for (nb, _) in neighbors[k] {
        if is_scaffold[nb] {
          deg = deg + 1
        }
      }
      if deg <= 1 {
        is_scaffold[k] = false
        changed = true
      }
    }
  }
  // Collect results
  let scaffold_atoms : Array[Int] = Array(capacity=0)
  let side_chain_atoms : Array[Int] = Array(capacity=0)
  for k = 0; k < n; k = k + 1 {
    if is_scaffold[k] {
      scaffold_atoms.push(k)
    } else {
      side_chain_atoms.push(k)
    }
  }
  // Count ring systems in scaffold
  let rings = mol_find_rings(mol)
  let mut num_rings = 0
  for ring in rings {
    let mut all_scaffold = true
    for a in ring {
      if not(is_scaffold[a]) {
        all_scaffold = false
      }
    }
    if all_scaffold {
      num_rings = num_rings + 1
    }
  }
  MurckoResult::{ scaffold_atoms, side_chain_atoms, num_rings }
}

///|
/// Build a new Molecule containing only the scaffold atoms.
pub fn murcko_scaffold_mol(mol : Molecule) -> Molecule {
  let result = murcko_scaffold(mol)
  if result.scaffold_atoms.length() == 0 {
    return Molecule::{ atoms: [], bonds: [] }
  }
  // Map old indices to new indices
  let n = mol.atoms.length()
  let old_to_new = Array::make(n, -1)
  let new_atoms : Array[MolAtom] = Array(capacity=0)
  for i = 0; i < result.scaffold_atoms.length(); i = i + 1 {
    let old_idx = result.scaffold_atoms[i]
    old_to_new[old_idx] = i
    new_atoms.push(mol.atoms[old_idx])
  }
  // Copy bonds where both endpoints are in scaffold
  let new_bonds : Array[MolBond] = Array(capacity=0)
  for b in mol.bonds {
    let nf = old_to_new[b.from]
    let nt = old_to_new[b.to]
    if nf >= 0 && nt >= 0 {
      new_bonds.push(MolBond::{ from: nf, to: nt, order: b.order })
    }
  }
  Molecule::{ atoms: new_atoms, bonds: new_bonds }
}

///|
/// Count side chain atoms.
pub fn murcko_side_chain_count(mol : Molecule) -> Int {
  let result = murcko_scaffold(mol)
  result.side_chain_atoms.length()
}

///|
/// Check if molecule is its own scaffold (no side chains).
pub fn murcko_is_pure_scaffold(mol : Molecule) -> Bool {
  let result = murcko_scaffold(mol)
  result.side_chain_atoms.length() == 0
}

///|
/// Get Murcko scaffold as SMILES.
pub fn murcko_scaffold_smiles(mol : Molecule) -> String {
  let scaffold = murcko_scaffold_mol(mol)
  mol_to_smiles(scaffold)
}
