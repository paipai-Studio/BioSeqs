/// R-group decomposition.
///
/// RDKit equivalent: `rdkit.Chem.rdRGroupDecomposition`
///
/// Given a scaffold and a set of molecules, decompose each molecule
/// into its R-group substituents relative to the scaffold.
/// This is core to SAR (Structure-Activity Relationship) analysis.

///|
/// An R-group substituent
pub struct RdkitRGroup {
  label : String // "R1", "R2", etc.
  smiles : String // SMILES of the substituent
  atom_indices : Array[Int32] // atoms in the original molecule
}

///|
pub fn RdkitRGroup::label(self : RdkitRGroup) -> String {
  self.label
}

///|
pub fn RdkitRGroup::smiles(self : RdkitRGroup) -> String {
  self.smiles
}

///|
pub fn RdkitRGroup::atom_indices(self : RdkitRGroup) -> Array[Int32] {
  self.atom_indices
}

///|
/// Result of R-group decomposition for one molecule
pub struct RdkitRGroupDecomp {
  scaffold_match : Bool
  r_groups : Array[RdkitRGroup]
}

///|
pub fn RdkitRGroupDecomp::scaffold_match(self : RdkitRGroupDecomp) -> Bool {
  self.scaffold_match
}

///|
pub fn RdkitRGroupDecomp::r_groups(
  self : RdkitRGroupDecomp,
) -> Array[RdkitRGroup] {
  self.r_groups
}

///|
/// Decompose a molecule into R-groups relative to a scaffold.
///
/// scaffold: the core structure (Molecule)
/// mol: the molecule to decompose
///
/// Strategy: find the scaffold as a substructure, then identify
/// the attachment points (scaffold atoms with external bonds)
/// and extract the substituents.
pub fn rdkit_rgroup_decompose(
  scaffold : Molecule,
  mol : Molecule,
) -> RdkitRGroupDecomp {
  // Use MCS to find the scaffold mapping
  let mcs = find_mcs(scaffold, mol)
  if mcs.num_atoms == 0 {
    return { scaffold_match: false, r_groups: [], }
  }
  // Build set of matched atoms in mol
  let matched_atoms : Array[Int32] = Array(capacity=0)
  for pair in mcs.mapping {
    matched_atoms.push(pair.1) // atom index in mol
  }
  // Find attachment points: matched atoms with bonds to unmatched atoms
  let r_groups : Array[RdkitRGroup] = Array(capacity=0)
  let mut r_label = 1
  for atom_idx in matched_atoms {
    // Check all bonds of this atom
    for bond in mol.bonds {
      let neighbor = if bond.from == atom_idx {
        bond.to
      } else if bond.to == atom_idx {
        bond.from
      } else {
        -1
      }
      if neighbor >= 0 {
        // Check if neighbor is unmatched
        let mut is_matched = false
        for ma in matched_atoms {
          if ma == neighbor {
            is_matched = true
            break
          }
        }
        if !is_matched {
          // Found an R-group attachment
          // Collect all atoms reachable from neighbor without going through matched atoms
          let group_atoms = rdkit_collect_substituent(
            mol, neighbor, matched_atoms,
          )
          // Build SMILES for the substituent (simplified)
          let group_smiles = rdkit_substituent_smiles(mol, group_atoms)
          r_groups.push({
            label: "R" + r_label.to_string(),
            smiles: group_smiles,
            atom_indices: group_atoms,
          })
          r_label = r_label + 1
        }
      }
    }
  }
  { scaffold_match: true, r_groups, }
}

///|
/// Collect all atoms reachable from a starting atom without crossing
/// through matched (scaffold) atoms.
fn rdkit_collect_substituent(
  mol : Molecule,
  start : Int32,
  matched : Array[Int32],
) -> Array[Int32] {
  let visited : Array[Bool] = Array::make(mol.num_atoms(), false)
  let result : Array[Int32] = Array(capacity=0)
  let queue : Array[Int32] = [start]
  visited[start] = true
  while queue.length() > 0 {
    let current = queue.remove(0)
    result.push(current)
    // Find neighbors
    for bond in mol.bonds {
      let neighbor = if bond.from == current {
        bond.to
      } else if bond.to == current {
        bond.from
      } else {
        -1
      }
      if neighbor >= 0 && !visited[neighbor] {
        // Check if neighbor is a matched atom
        let mut is_matched = false
        for m in matched {
          if m == neighbor {
            is_matched = true
            break
          }
        }
        if !is_matched {
          visited[neighbor] = true
          queue.push(neighbor)
        }
      }
    }
  }
  result
}

///|
/// Build a simplified SMILES for a set of atoms.
fn rdkit_substituent_smiles(mol : Molecule, atoms : Array[Int32]) -> String {
  let sb = StringBuilder()
  for i = 0; i < atoms.length(); i = i + 1 {
    if i > 0 {
      sb.write_char('-')
    }
    sb.write_string(mol.atoms[atoms[i]].symbol)
  }
  sb.to_string()
}

///|
/// Batch R-group decomposition for multiple molecules.
pub fn rdkit_rgroup_decompose_batch(
  scaffold : Molecule,
  mols : Array[Molecule],
) -> Array[RdkitRGroupDecomp] {
  let results : Array[RdkitRGroupDecomp] = Array(capacity=0)
  for mol in mols {
    results.push(rdkit_rgroup_decompose(scaffold, mol))
  }
  results
}

///|
/// Get all unique R-group SMILES across a set of decompositions.
pub fn rdkit_rgroup_unique_substituents(
  decomp : Array[RdkitRGroupDecomp],
) -> Array[String] {
  let unique : Array[String] = Array(capacity=0)
  for d in decomp {
    for rg in d.r_groups {
      let mut found = false
      for u in unique {
        if u == rg.smiles {
          found = true
          break
        }
      }
      if !found {
        unique.push(rg.smiles)
      }
    }
  }
  unique
}
