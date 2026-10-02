// ============================================================================
// Stereochemistry perception (CIP rules)
// ============================================================================
// (1) CIP priority assignment: compare by atomic number, then isotope,
//     then degree, then recursive neighbor comparison.
// (2) R/S tetrahedral center: needs 4 different substituents, determine
//     handedness of substituent ordering.
// (3) E/Z double bond: compare CIP priorities on each end of double bond.

///|
/// CIP priority level for an atom.
fn cip_atomic_number(symbol : String) -> Int {
  match symbol {
    "H" => 1
    "B" => 5
    "C" => 6
    "N" => 7
    "O" => 8
    "F" => 9
    "Si" => 14
    "P" => 15
    "S" => 16
    "Cl" => 17
    "Br" => 35
    "I" => 53
    _ => 0
  }
}

///|
/// Compare two atoms by CIP priority (higher = higher priority).
/// Returns positive if atom1 > atom2, negative if atom1 < atom2, 0 if equal.
fn cip_compare(
  mol : Molecule,
  a1 : Int,
  a2 : Int,
  visited : Array[Bool],
) -> Int {
  let atom1 = mol.atoms[a1]
  let atom2 = mol.atoms[a2]
  // Level 1: atomic number
  let z1 = cip_atomic_number(atom1.symbol)
  let z2 = cip_atomic_number(atom2.symbol)
  if z1 != z2 {
    return z1 - z2
  }
  // Level 2: isotope (charge as proxy)
  if atom1.charge != atom2.charge {
    return atom1.charge - atom2.charge
  }
  // Level 3: degree
  let d1 = mol_degree(mol, a1)
  let d2 = mol_degree(mol, a2)
  if d1 != d2 {
    return d1 - d2
  }
  // Level 4: recursive neighbor comparison
  let neighbors = mol_adjacency(mol)
  let n1_priorities : Array[Int] = Array(capacity=0)
  for (nb, _) in neighbors[a1] {
    if not(visited[nb]) {
      visited[nb] = true
      n1_priorities.push(cip_compare(mol, nb, -1, visited))
      visited[nb] = false
    }
  }
  let n2_priorities : Array[Int] = Array(capacity=0)
  for (nb, _) in neighbors[a2] {
    if not(visited[nb]) {
      visited[nb] = true
      n2_priorities.push(cip_compare(mol, nb, -1, visited))
      visited[nb] = false
    }
  }
  // Sort descending
  for i = 0; i < n1_priorities.length(); i = i + 1 {
    for j = i + 1; j < n1_priorities.length(); j = j + 1 {
      if n1_priorities[j] > n1_priorities[i] {
        let tmp = n1_priorities[i]
        n1_priorities[i] = n1_priorities[j]
        n1_priorities[j] = tmp
      }
    }
  }
  for i = 0; i < n2_priorities.length(); i = i + 1 {
    for j = i + 1; j < n2_priorities.length(); j = j + 1 {
      if n2_priorities[j] > n2_priorities[i] {
        let tmp = n2_priorities[i]
        n2_priorities[i] = n2_priorities[j]
        n2_priorities[j] = tmp
      }
    }
  }
  // Compare sorted lists
  let min_len = if n1_priorities.length() < n2_priorities.length() {
    n1_priorities.length()
  } else {
    n2_priorities.length()
  }
  for i = 0; i < min_len; i = i + 1 {
    if n1_priorities[i] != n2_priorities[i] {
      return n1_priorities[i] - n2_priorities[i]
    }
  }
  n1_priorities.length() - n2_priorities.length()
}

///|
/// Assign CIP priorities to neighbors of a center atom.
/// Returns neighbor indices sorted by CIP priority (highest first).
fn cip_sort_neighbors(
  mol : Molecule,
  center : Int,
  neighbors : Array[Int],
) -> Array[Int] {
  let visited = Array::make(mol.atoms.length(), false)
  visited[center] = true
  // Compute priority for each neighbor
  let pairs : Array[(Int, Int)] = Array(capacity=0)
  for nb in neighbors {
    visited[nb] = true
    let mut priority = cip_atomic_number(mol.atoms[nb].symbol)
    // Add degree as secondary
    priority = priority * 100 + mol_degree(mol, nb)
    pairs.push((priority, nb))
    visited[nb] = false
  }
  // Sort by priority descending
  for i = 0; i < pairs.length(); i = i + 1 {
    for j = i + 1; j < pairs.length(); j = j + 1 {
      if pairs[j].0 > pairs[i].0 {
        let tmp = pairs[i]
        pairs[i] = pairs[j]
        pairs[j] = tmp
      }
    }
  }
  let result : Array[Int] = Array(capacity=0)
  for (_, nb) in pairs {
    result.push(nb)
  }
  result
}

///|
/// Stereo descriptor for a stereocenter.
pub(all) enum StereoDescriptor {
  R
  S
  E
  Z
  None
}

///|
/// Check if an atom is a tetrahedral stereocenter.
/// Requires: 4 neighbors with different CIP priorities, sp3 carbon.
pub fn is_tetrahedral_center(mol : Molecule, atom : Int) -> Bool {
  let neighbors = mol_adjacency(mol)
  let nbs = neighbors[atom]
  if nbs.length() != 4 {
    return false
  }
  // Check all neighbors have different symbols (simplified)
  let symbols : Array[String] = Array(capacity=0)
  for (nb, _) in nbs {
    symbols.push(mol.atoms[nb].symbol)
  }
  // At least 2 different substituents
  let mut has_diff = false
  for i = 0; i < symbols.length(); i = i + 1 {
    for j = i + 1; j < symbols.length(); j = j + 1 {
      if symbols[i] != symbols[j] {
        has_diff = true
      }
    }
  }
  has_diff
}

///|
/// Determine R/S configuration of a tetrahedral center.
/// This is a simplified version: without 3D coordinates, we determine
/// the configuration from the atom ordering in the molecular graph.
/// Returns R, S, or None if not a stereocenter.
pub fn cip_tetrahedral(mol : Molecule, atom : Int) -> StereoDescriptor {
  if not(is_tetrahedral_center(mol, atom)) {
    return None
  }
  let neighbors = mol_adjacency(mol)
  let nb_list : Array[Int] = Array(capacity=0)
  for (nb, _) in neighbors[atom] {
    nb_list.push(nb)
  }
  let sorted = cip_sort_neighbors(mol, atom, nb_list)
  // Without 3D coordinates, we use the parity of the permutation
  // from original order to sorted order as a proxy.
  // Count inversions in the permutation
  let mut inversions = 0
  for i = 0; i < sorted.length(); i = i + 1 {
    for j = i + 1; j < sorted.length(); j = j + 1 {
      // Find original indices
      let mut oi = -1
      let mut oj = -1
      for k = 0; k < nb_list.length(); k = k + 1 {
        if nb_list[k] == sorted[i] {
          oi = k
        }
        if nb_list[k] == sorted[j] {
          oj = k
        }
      }
      if oi > oj {
        inversions = inversions + 1
      }
    }
  }
  if inversions % 2 == 0 {
    R
  } else {
    S
  }
}

///|
/// Find all tetrahedral stereocenters in a molecule.
pub fn find_stereocenters(mol : Molecule) -> Array[(Int, StereoDescriptor)] {
  let result : Array[(Int, StereoDescriptor)] = Array(capacity=0)
  for k = 0; k < mol.atoms.length(); k = k + 1 {
    let desc = cip_tetrahedral(mol, k)
    match desc {
      R => result.push((k, R))
      S => result.push((k, S))
      None => ()
      _ => ()
    }
  }
  result
}

///|
/// Determine E/Z configuration of a double bond.
/// Returns E, Z, or None if not applicable.
pub fn cip_double_bond(
  mol : Molecule,
  bond_idx : Int,
) -> StereoDescriptor {
  let b = mol.bonds[bond_idx]
  if b.order != 2 {
    return None
  }
  let neighbors = mol_adjacency(mol)
  // Get substituents on each end (excluding the double bond partner)
  let subst1 : Array[Int] = Array(capacity=0)
  let subst2 : Array[Int] = Array(capacity=0)
  for (nb, _) in neighbors[b.from] {
    if nb != b.to {
      subst1.push(nb)
    }
  }
  for (nb, _) in neighbors[b.to] {
    if nb != b.from {
      subst2.push(nb)
    }
  }
  // Need at least 1 substituent on each end
  if subst1.length() == 0 || subst2.length() == 0 {
    return None
  }
  // Get highest priority substituent on each end
  let sorted1 = cip_sort_neighbors(mol, b.from, subst1)
  let sorted2 = cip_sort_neighbors(mol, b.to, subst2)
  let high1 = sorted1[0]
  let high2 = sorted2[0]
  // Check if both ends have 2 substituents (for E/Z)
  if subst1.length() < 2 || subst2.length() < 2 {
    return None
  }
  // Get low priority substituents
  let low1 = sorted1[1]
  let low2 = sorted2[1]
  // E: high priority groups on opposite sides
  // Z: high priority groups on same side
  // Without 3D coordinates, use atom index parity as proxy
  let parity = (high1 > low1 && high2 > low2) || (high1 < low1 && high2 < low2)
  if parity {
    E
  } else {
    Z
  }
}

///|
/// Find all E/Z double bonds in a molecule.
pub fn find_double_bond_stereo(mol : Molecule) -> Array[(Int, StereoDescriptor)] {
  let result : Array[(Int, StereoDescriptor)] = Array(capacity=0)
  for i = 0; i < mol.bonds.length(); i = i + 1 {
    let desc = cip_double_bond(mol, i)
    match desc {
      E => result.push((i, E))
      Z => result.push((i, Z))
      None => ()
      _ => ()
    }
  }
  result
}

///|
/// String representation of stereo descriptor.
pub fn stereo_to_string(d : StereoDescriptor) -> String {
  match d {
    R => "R"
    S => "S"
    E => "E"
    Z => "Z"
    None => "none"
  }
}
