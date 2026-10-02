// ============================================================================
// Topological Torsion fingerprint
// ============================================================================
// A 4-atom path fingerprint (Rucker & Mata, J. Chem. Inf. Model. 2006).
// Each torsion is a linear path of 4 non-ring atoms connected by 3 bonds.
// The torsion is encoded by hashing the atom invariants and bond orders
// along the path, with canonical ordering (smallest end first).

///|
/// Atom invariant for torsion encoding: element + degree + aromaticity.
fn torsion_atom_invariant(mol : Molecule, atom : Int) -> Int {
  let a = mol.atoms[atom]
  let elem = match a.symbol {
    "C" => 6
    "N" => 7
    "O" => 8
    "S" => 16
    "P" => 15
    "F" => 9
    "Cl" => 17
    "Br" => 35
    "I" => 53
    "H" => 1
    _ => 99
  }
  let arom = if a.aromatic { 1 } else { 0 }
  let deg = mol_degree(mol, atom)
  elem * 1000 + deg * 10 + arom + a.charge * 100
}

///|
/// FNV-1a hash for a torsion path.
fn torsion_hash(values : Array[Int]) -> Int64 {
  let mut h : Int64 = 1469598103934665603L
  for v in values {
    h = h ^ v.to_int64()
    h = h * 1099511628211L
  }
  h
}

///|
/// Encode a 4-atom torsion path as a hash.
/// The path is canonicalized by comparing forward and reverse.
fn torsion_encode(
  mol : Molecule,
  path : Array[Int],
) -> Int64 {
  // Forward encoding: [inv0, bond01, inv1, bond12, inv2, bond23, inv3]
  let forward : Array[Int] = Array(capacity=0)
  forward.push(torsion_atom_invariant(mol, path[0]))
  forward.push(mol_bond_order(mol, path[0], path[1]))
  forward.push(torsion_atom_invariant(mol, path[1]))
  forward.push(mol_bond_order(mol, path[1], path[2]))
  forward.push(torsion_atom_invariant(mol, path[2]))
  forward.push(mol_bond_order(mol, path[2], path[3]))
  forward.push(torsion_atom_invariant(mol, path[3]))
  // Reverse encoding
  let reverse : Array[Int] = Array(capacity=0)
  reverse.push(torsion_atom_invariant(mol, path[3]))
  reverse.push(mol_bond_order(mol, path[2], path[3]))
  reverse.push(torsion_atom_invariant(mol, path[2]))
  reverse.push(mol_bond_order(mol, path[1], path[2]))
  reverse.push(torsion_atom_invariant(mol, path[1]))
  reverse.push(mol_bond_order(mol, path[0], path[1]))
  reverse.push(torsion_atom_invariant(mol, path[0]))
  // Use the smaller hash (canonical form)
  let hf = torsion_hash(forward)
  let hr = torsion_hash(reverse)
  if hf < hr {
    hf
  } else {
    hr
  }
}

///|
/// Enumerate all 4-atom linear paths in the molecule.
/// A path visits each atom at most once (no backtracking).
fn torsion_find_paths(mol : Molecule) -> Array[Array[Int]] {
  let n = mol.atoms.length()
  let neighbors = mol_adjacency(mol)
  let paths : Array[Array[Int]] = Array(capacity=0)
  // DFS from each starting atom
  for start = 0; start < n; start = start + 1 {
    torsion_dfs_paths(mol, neighbors, [start], paths)
  }
  paths
}

///|
/// DFS to extend paths to length 4.
fn torsion_dfs_paths(
  mol : Molecule,
  neighbors : Array[Array[(Int, Int)]],
  current : Array[Int],
  paths : Array[Array[Int]],
) -> Unit {
  if current.length() == 4 {
    // Add a copy of the path
    let copy = Array(capacity=0)
    for x in current {
      copy.push(x)
    }
    paths.push(copy)
    return
  }
  let last = current[current.length() - 1]
  for (nb, _) in neighbors[last] {
    // Check if nb is already in path
    let mut in_path = false
    for x in current {
      if x == nb {
        in_path = true
      }
    }
    if in_path {
      continue
    }
    current.push(nb)
    torsion_dfs_paths(mol, neighbors, current, paths)
    let _ = current.pop()
  }
}

///|
/// Compute topological torsion fingerprint.
/// Returns an array of unique torsion hashes.
pub fn topological_torsion_fingerprint(mol : Molecule) -> Array[Int64] {
  let paths = torsion_find_paths(mol)
  let hashes : Array[Int64] = Array(capacity=0)
  let seen : Map[Int64, Bool] = Map([], capacity=0)
  for path in paths {
    let h = torsion_encode(mol, path)
    match seen.get(h) {
      Some(_) => ()
      None => {
        seen[h] = true
        hashes.push(h)
      }
    }
  }
  hashes
}

///|
/// Count number of unique torsions.
pub fn topological_torsion_count(mol : Molecule) -> Int {
  topological_torsion_fingerprint(mol).length()
}

///|
/// Tanimoto similarity between two torsion fingerprints.
pub fn torsion_tanimoto(
  mol1 : Molecule,
  mol2 : Molecule,
) -> Double {
  let fp1 = topological_torsion_fingerprint(mol1)
  let fp2 = topological_torsion_fingerprint(mol2)
  let set1 : Map[Int64, Bool] = Map([], capacity=0)
  for h in fp1 {
    set1[h] = true
  }
  let set2 : Map[Int64, Bool] = Map([], capacity=0)
  for h in fp2 {
    set2[h] = true
  }
  let mut intersection = 0
  let mut union = 0
  for h in fp1 {
    match set2.get(h) {
      Some(_) => intersection = intersection + 1
      None => ()
    }
  }
  union = fp1.length() + fp2.length() - intersection
  if union == 0 {
    return 0.0
  }
  intersection.to_double() / union.to_double()
}

///|
/// Get torsion fingerprint as a sorted list of hash values (for debugging).
pub fn torsion_fingerprint_string(mol : Molecule) -> String {
  let fp = topological_torsion_fingerprint(mol)
  let sb = StringBuilder()
  for i = 0; i < fp.length(); i = i + 1 {
    if i > 0 {
      sb.write_string(", ")
    }
    sb.write_string(fp[i].to_string())
  }
  sb.to_string()
}
