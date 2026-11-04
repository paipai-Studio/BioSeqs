/// RDKit topological/path fingerprint (RDKFingerprint).
///
/// RDKit equivalent: `rdkit.Chem.RDKFingerprint`
///
/// Unlike Morgan (circular) fingerprints, the RDKit path fingerprint
/// enumerates all linear paths up to a maximum length and hashes each
/// path into a bit vector. This captures different structural information
/// and is complementary to Morgan/ECFP.

///|
/// Generate RDKit path fingerprint as a bit vector.
///
/// Enumerates all paths up to max_path_len (default 7) in the molecular
/// graph, hashes each path, and sets bits in a fingerprint of n_bits (default 2048).
pub fn rdkit_path_fingerprint(
  mol : Molecule,
  max_path_len~ : Int32 = 7,
  n_bits~ : Int32 = 2048,
) -> Array[Bool] {
  let n = mol.num_atoms()
  let fp : Array[Bool] = Array::make(n_bits, false)
  if n == 0 {
    return fp
  }
  // Build adjacency list with bond orders
  let adj : Array[Array[(Int32, Int32)]] = Array::make(n, [])
  for bi = 0; bi < mol.bonds.length(); bi = bi + 1 {
    let bond = mol.bonds[bi]
    adj[bond.from].push((bond.to, bond.order))
    adj[bond.to].push((bond.from, bond.order))
  }
  // Enumerate all paths starting from each atom
  for start = 0; start < n; start = start + 1 {
    rdkit_enum_paths(mol, adj, start, max_path_len, fp, n_bits)
  }
  fp
}

///|
/// Enumerate all paths from a starting atom (DFS).
fn rdkit_enum_paths(
  mol : Molecule,
  adj : Array[Array[(Int32, Int32)]],
  start : Int32,
  max_len : Int32,
  fp : Array[Bool],
  n_bits : Int32,
) -> Unit {
  // Stack: (current_atom, visited_set, path_hash, path_length)
  // Use iterative DFS to avoid deep recursion
  let visited : Array[Bool] = Array::make(mol.num_atoms(), false)
  visited[start] = true
  // Initial hash from atom symbol
  let init_hash = rdkit_atom_hash(mol, start)
  rdkit_dfs_paths(mol, adj, start, visited, init_hash, 1, max_len, fp, n_bits)
  visited[start] = false
}

///|
fn rdkit_dfs_paths(
  mol : Molecule,
  adj : Array[Array[(Int32, Int32)]],
  current : Int32,
  visited : Array[Bool],
  path_hash : Int32,
  path_len : Int32,
  max_len : Int32,
  fp : Array[Bool],
  n_bits : Int32,
) -> Unit {
  // Set bit for current path (if length >= 1)
  if path_len >= 1 {
    let bit = (path_hash % n_bits + n_bits) % n_bits
    fp[bit] = true
  }
  // Stop if max length reached
  if path_len >= max_len {
    return
  }
  // Continue DFS
  for neighbor_info in adj[current] {
    let next = neighbor_info.0
    let bond_order = neighbor_info.1
    if !visited[next] {
      visited[next] = true
      // Extend hash: combine with bond order and next atom
      let new_hash = path_hash * 31 +
        bond_order * 17 +
        rdkit_atom_hash(mol, next)
      rdkit_dfs_paths(
        mol,
        adj,
        next,
        visited,
        new_hash,
        path_len + 1,
        max_len,
        fp,
        n_bits,
      )
      visited[next] = false
    }
  }
}

///|
/// Simple hash for an atom: based on symbol, charge, aromaticity.
fn rdkit_atom_hash(mol : Molecule, idx : Int32) -> Int32 {
  let atom = mol.atoms[idx]
  let mut h = 0
  for c in atom.symbol {
    h = h * 31 + c.to_int()
  }
  h = h * 17 + atom.charge
  h = h * 7 + (if atom.aromatic { 1 } else { 0 })
  h
}

///|
/// Generate path fingerprint as integer array (bit indices that are set).
pub fn rdkit_path_fingerprint_bits(
  mol : Molecule,
  max_path_len~ : Int32 = 7,
  n_bits~ : Int32 = 2048,
) -> Array[Int32] {
  let fp = rdkit_path_fingerprint(mol, max_path_len~, n_bits~)
  let bits : Array[Int32] = Array(capacity=0)
  for i = 0; i < fp.length(); i = i + 1 {
    if fp[i] {
      bits.push(i)
    }
  }
  bits
}

///|
/// Tanimoto similarity between two path fingerprints.
pub fn rdkit_path_tanimoto(mol1 : Molecule, mol2 : Molecule) -> Double {
  let fp1 = rdkit_path_fingerprint(mol1)
  let fp2 = rdkit_path_fingerprint(mol2)
  let len = if fp1.length() < fp2.length() {
    fp1.length()
  } else {
    fp2.length()
  }
  let mut intersection = 0
  let mut union = 0
  for i = 0; i < len; i = i + 1 {
    if fp1[i] || fp2[i] {
      union = union + 1
    }
    if fp1[i] && fp2[i] {
      intersection = intersection + 1
    }
  }
  if union == 0 {
    return 0.0
  }
  intersection.to_double() / union.to_double()
}

///|
/// Dice similarity between two path fingerprints.
/// Dice = 2|A ∩ B| / (|A| + |B|)
pub fn rdkit_path_dice(mol1 : Molecule, mol2 : Molecule) -> Double {
  let fp1 = rdkit_path_fingerprint(mol1)
  let fp2 = rdkit_path_fingerprint(mol2)
  let len = if fp1.length() < fp2.length() {
    fp1.length()
  } else {
    fp2.length()
  }
  let mut intersection = 0
  let mut count1 = 0
  let mut count2 = 0
  for i = 0; i < len; i = i + 1 {
    if fp1[i] {
      count1 = count1 + 1
    }
    if fp2[i] {
      count2 = count2 + 1
    }
    if fp1[i] && fp2[i] {
      intersection = intersection + 1
    }
  }
  let total = count1 + count2
  if total == 0 {
    return 0.0
  }
  2.0 * intersection.to_double() / total.to_double()
}

///|
/// Count number of bits set in path fingerprint.
pub fn rdkit_path_bit_count(mol : Molecule) -> Int32 {
  let fp = rdkit_path_fingerprint(mol)
  let mut count = 0
  for i = 0; i < fp.length(); i = i + 1 {
    if fp[i] {
      count = count + 1
    }
  }
  count
}
