/// Canonical SMILES generation.
///
/// RDKit equivalent: `rdkit.Chem.MolToSmiles(mol, canonical=True)`
///
/// Canonical SMILES uses the Morgan algorithm for invariant computation
/// to produce a unique string representation regardless of atom numbering.

///|
/// Compute Morgan invariants (extended connectivity values) for all atoms.
/// Each round, each atom's invariant is updated based on its neighbors' invariants.
fn canonical_morgan_invariants(
  mol : Molecule,
  rounds~ : Int32 = 10,
) -> Array[Int32] {
  let n = mol.num_atoms()
  if n == 0 {
    return []
  }
  // Initial invariants: hash of (symbol, charge, aromatic, degree)
  let invariants : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    let atom = mol.atoms[i]
    let deg = mol_degree(mol, i)
    // Simple hash: combine properties into a single integer
    let mut h = 0
    for c in atom.symbol {
      h = h * 31 + c.to_int()
    }
    h = h * 17 + atom.charge
    h = h * 7 + (if atom.aromatic { 1 } else { 0 })
    h = h * 13 + deg
    invariants[i] = h
  }
  // Build adjacency list
  let adj : Array[Array[Int32]] = Array::make(n, [])
  for bond in mol.bonds {
    adj[bond.from].push(bond.to)
    adj[bond.to].push(bond.from)
  }
  // Iterative refinement
  for r = 0; r < rounds; r = r + 1 {
    let new_inv = invariants.copy()
    for i = 0; i < n; i = i + 1 {
      let mut h = invariants[i]
      for j in adj[i] {
        h = h * 31 + invariants[j]
      }
      new_inv[i] = h
    }
    // Check if invariants changed
    let mut changed = false
    for i = 0; i < n; i = i + 1 {
      if new_inv[i] != invariants[i] {
        changed = true
        break
      }
    }
    for i = 0; i < n; i = i + 1 {
      invariants[i] = new_inv[i]
    }
    if !changed {
      break
    }
  }
  invariants
}

///|
/// Get canonical atom ordering: sort atoms by their Morgan invariants.
/// Ties are broken by (degree, symbol, charge).
fn canonical_atom_order(mol : Molecule) -> Array[Int32] {
  let n = mol.num_atoms()
  if n == 0 {
    return []
  }
  let invariants = canonical_morgan_invariants(mol)
  let order : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    order[i] = i
  }
  // Sort by invariant, then by degree, then by symbol
  order.sort_by(fn(a, b) {
    if invariants[a] != invariants[b] {
      invariants[a] - invariants[b]
    } else {
      let da = mol_degree(mol, a)
      let db = mol_degree(mol, b)
      if da != db {
        da - db
      } else {
        let sa = mol.atoms[a].symbol
        let sb = mol.atoms[b].symbol
        if sa < sb {
          -1
        } else if sa > sb {
          1
        } else {
          0
        }
      }
    }
  })
  order
}

///|
/// Generate canonical SMILES for a molecule.
/// Reorders atoms by Morgan invariants, then generates SMILES.
pub fn rdkit_canonical_smiles(mol : Molecule) -> String {
  let n = mol.num_atoms()
  if n == 0 {
    return ""
  }
  let order = canonical_atom_order(mol)
  // Build inverse mapping: old_index -> new_index
  let rank : Array[Int32] = Array::make(n, 0)
  for new_idx = 0; new_idx < n; new_idx = new_idx + 1 {
    // avoid closure capture issue
    rank[order[new_idx]] = new_idx
  }
  // Build canonical molecule with reordered atoms
  let new_atoms : Array[MolAtom] = Array::make(n, {
    symbol: "",
    charge: 0,
    aromatic: false,
  })
  for i = 0; i < n; i = i + 1 {
    new_atoms[i] = mol.atoms[order[i]]
  }
  let new_bonds : Array[MolBond] = Array(capacity=mol.bonds.length())
  for bond in mol.bonds {
    new_bonds.push({
      from: rank[bond.from],
      to: rank[bond.to],
      order: bond.order,
    })
  }
  let canonical_mol = { atoms: new_atoms, bonds: new_bonds, }
  mol_to_smiles(canonical_mol)
}

///|
/// Check if two molecules are the same (canonical SMILES match).
pub fn rdkit_molecules_identical(mol1 : Molecule, mol2 : Molecule) -> Bool {
  rdkit_canonical_smiles(mol1) == rdkit_canonical_smiles(mol2)
}

///|
/// Get the Morgan invariant ranking for each atom.
/// Useful for understanding canonical ordering.
pub fn rdkit_atom_invariants(mol : Molecule) -> Array[Int32] {
  canonical_morgan_invariants(mol)
}

///|
/// Count distinct atom environments (number of unique invariants).
pub fn rdkit_unique_environments(mol : Molecule) -> Int32 {
  let inv = canonical_morgan_invariants(mol)
  let unique : Array[Int32] = Array(capacity=0)
  for v in inv {
    let mut found = false
    for u in unique {
      if u == v {
        found = true
        break
      }
    }
    if !found {
      unique.push(v)
    }
  }
  unique.length()
}
