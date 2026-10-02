// ============================================================================
// Maximum Common Substructure (MCS)
// ============================================================================
// Finds the largest common substructure between two molecules using
// a product graph + maximum clique approach (Bron-Kerbosch algorithm).
// Reference: Clique-based MCS: Raymond et al., J. Chem. Inf. Comput. Sci. 2002.

///|
/// MCS result: matched atom pairs and bond count.
pub struct McsResult {
  mapping : Array[(Int, Int)]  // (atom_in_mol1, atom_in_mol2)
  num_atoms : Int
  num_bonds : Int
}

///|
/// Check if two atoms are compatible for MCS matching.
fn mcs_atom_compatible(
  mol1 : Molecule,
  a1 : Int,
  mol2 : Molecule,
  a2 : Int,
) -> Bool {
  mol1.atoms[a1].symbol == mol2.atoms[a2].symbol
}

///|
/// Check if two atom pairs are bond-compatible in the product graph.
fn mcs_bond_compatible(
  mol1 : Molecule,
  a1 : Int,
  b1 : Int,
  mol2 : Molecule,
  a2 : Int,
  b2 : Int,
) -> Bool {
  let o1 = mol_bond_order(mol1, a1, b1)
  let o2 = mol_bond_order(mol2, a2, b2)
  // -1 means no bond; both must have same bond order (including no bond = no edge)
  if o1 == -1 || o2 == -1 {
    return false
  }
  o1 == o2
}

///|
/// Build the product graph: nodes are compatible atom pairs,
/// edges connect pairs that are bond-compatible.
fn mcs_build_product_graph(
  mol1 : Molecule,
  mol2 : Molecule,
) -> (Array[(Int, Int)], Array[Array[Int]]) {
  let n1 = mol1.atoms.length()
  let n2 = mol2.atoms.length()
  // Generate candidate node pairs
  let nodes : Array[(Int, Int)] = Array(capacity=0)
  for i = 0; i < n1; i = i + 1 {
    for j = 0; j < n2; j = j + 1 {
      if mcs_atom_compatible(mol1, i, mol2, j) {
        nodes.push((i, j))
      }
    }
  }
  // Build adjacency: edge between node pairs if bond-compatible
  let nn = nodes.length()
  let adj : Array[Array[Int]] = Array::makei(nn, { _ => Array(capacity=0) })
  for i = 0; i < nn; i = i + 1 {
    for j = i + 1; j < nn; j = j + 1 {
      let (a1, a2) = nodes[i]
      let (b1, b2) = nodes[j]
      // Must not reuse same atom from either molecule
      if a1 == b1 || a2 == b2 {
        continue
      }
      if mcs_bond_compatible(mol1, a1, b1, mol2, a2, b2) {
        adj[i].push(j)
        adj[j].push(i)
      }
    }
  }
  (nodes, adj)
}

///|
/// Bron-Kerbosch maximum clique algorithm (with pivoting).
/// Returns the largest clique found.
fn mcs_bron_kerbosch(
  adj : Array[Array[Int]],
  r : Array[Int],
  p : Array[Int],
  best : Ref[Array[Int]],
) -> Unit {
  if p.length() == 0 {
    if r.length() > best.val.length() {
      best.val = Array(capacity=0)
      for x in r {
        best.val.push(x)
      }
    }
    return
  }
  // Pruning: if |R| + |P| <= |best|, cannot improve
  if r.length() + p.length() <= best.val.length() {
    return
  }
  // Pick pivot (first element of P for simplicity)
  let pivot = p[0]
  // Iterate over P \ N(pivot)
  for i = 0; i < p.length(); i = i + 1 {
    let v = p[i]
    // Skip if v is neighbor of pivot
    let mut is_neighbor = false
    for nb in adj[pivot] {
      if nb == v {
        is_neighbor = true
      }
    }
    if is_neighbor {
      continue
    }
    // R ∪ {v}
    let new_r = Array(capacity=0)
    for x in r {
      new_r.push(x)
    }
    new_r.push(v)
    // P ∩ N(v)
    let new_p = Array(capacity=0)
    for x in p {
      for nb in adj[v] {
        if nb == x {
          new_p.push(x)
        }
      }
    }
    mcs_bron_kerbosch(adj, new_r, new_p, best)
  }
}

///|
/// Find Maximum Common Substructure between two molecules.
pub fn find_mcs(mol1 : Molecule, mol2 : Molecule) -> McsResult {
  let (nodes, adj) = mcs_build_product_graph(mol1, mol2)
  if nodes.length() == 0 {
    return McsResult::{ mapping: [], num_atoms: 0, num_bonds: 0 }
  }
  // Initial P = all nodes
  let initial_p : Array[Int] = Array(capacity=0)
  for i = 0; i < nodes.length(); i = i + 1 {
    initial_p.push(i)
  }
  let best = Ref::new(Array(capacity=0))
  mcs_bron_kerbosch(adj, [], initial_p, best)
  // Build mapping from best clique
  let mapping : Array[(Int, Int)] = Array(capacity=0)
  for idx in best.val {
    mapping.push(nodes[idx])
  }
  // Count matched bonds
  let mut num_bonds = 0
  for i = 0; i < mapping.length(); i = i + 1 {
    for j = i + 1; j < mapping.length(); j = j + 1 {
      let (a1, _) = mapping[i]
      let (b1, _) = mapping[j]
      if mol_bond_order(mol1, a1, b1) > 0 {
        num_bonds = num_bonds + 1
      }
    }
  }
  McsResult::{ mapping, num_atoms: mapping.length(), num_bonds }
}

///|
/// MCS atom count (convenience).
pub fn mcs_atom_count(mol1 : Molecule, mol2 : Molecule) -> Int {
  find_mcs(mol1, mol2).num_atoms
}

///|
/// MCS bond count (convenience).
pub fn mcs_bond_count(mol1 : Molecule, mol2 : Molecule) -> Int {
  find_mcs(mol1, mol2).num_bonds
}

///|
/// Tanimoto-like similarity based on MCS size.
/// sim = |MCS| / max(|mol1|, |mol2|)
pub fn mcs_similarity(mol1 : Molecule, mol2 : Molecule) -> Double {
  let mcs = find_mcs(mol1, mol2)
  let max_atoms = if mol1.atoms.length() > mol2.atoms.length() {
    mol1.atoms.length()
  } else {
    mol2.atoms.length()
  }
  if max_atoms == 0 {
    return 0.0
  }
  mcs.num_atoms.to_double() / max_atoms.to_double()
}
