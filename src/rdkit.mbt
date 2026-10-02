///|
/// RDKit — cheminformatics core (SMILES parsing + Morgan/ECFP fingerprints).
///
/// A self-contained molecular graph + a SMILES parser for the common core
/// grammar (elements, aromatic lowercase atoms, bonds, branches, ring
/// closures, charges) plus the Morgan (circular / ECFP) fingerprint.

///|
/// A heavy atom in a molecule graph.
pub struct MolAtom {
  symbol : String
  charge : Int
  aromatic : Bool
}

///|
/// A bond between two atom indices.
pub struct MolBond {
  from : Int
  to : Int
  order : Int // 1 = single, 2 = double, 3 = triple
}

///|
/// A molecular graph.
pub struct Molecule {
  atoms : Array[MolAtom]
  bonds : Array[MolBond]
}

///|
/// Number of atoms.
pub fn Molecule::num_atoms(self : Molecule) -> Int {
  self.atoms.length()
}

///|
/// Number of bonds.
pub fn Molecule::num_bonds(self : Molecule) -> Int {
  self.bonds.length()
}

///|
/// Atom symbols in order.
pub fn Molecule::atom_symbols(self : Molecule) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for a in self.atoms {
    out.push(a.symbol)
  }
  out
}

///|
/// Bonds as `(from, to, order)` triples.
pub fn Molecule::bond_triples(self : Molecule) -> Array[(Int, Int, Int)] {
  let out : Array[(Int, Int, Int)] = Array(capacity=0)
  for b in self.bonds {
    out.push((b.from, b.to, b.order))
  }
  out
}

///|
fn mol_is_upper(c : UInt16) -> Bool {
  c >= 65 && c <= 90
}

///|
fn mol_is_lower(c : UInt16) -> Bool {
  c >= 97 && c <= 122
}

///|
fn mol_is_digit(c : UInt16) -> Bool {
  c >= 48 && c <= 57
}

///|
fn mol_digit_value(c : UInt16) -> Int {
  (c - 48).to_int()
}

///|
/// Parse an element symbol starting at `i`; returns (symbol, aromatic, length).
/// Lowercase letters are aromatic atoms (c/n/o/s/p); uppercase starts a normal
/// element (1 or 2 letters), and a bracket `[...]` is handled inline by taking
/// the first element symbol inside.
fn mol_parse_element(s : String, i : Int) -> (String, Bool, Int) {
  let n = s.length()
  let c = s.unsafe_get(i)
  if mol_is_lower(c) {
    // 芳香原子（小写）
    let sym = mol_upper_symbol(c)
    (sym, true, 1)
  } else if mol_is_upper(c) {
    // 两字母元素（Cl/Br/Si...）
    if i + 1 < n && mol_is_lower(s.unsafe_get(i + 1)) {
      let sb = StringBuilder()
      sb.write_char(c.unsafe_to_char())
      sb.write_char(s.unsafe_get(i + 1).unsafe_to_char())
      (sb.to_string(), false, 2)
    } else {
      let sb = StringBuilder()
      sb.write_char(c.unsafe_to_char())
      (sb.to_string(), false, 1)
    }
  } else if c.to_int() == 91 {
    // '[' 括号内：提取第一个大写元素（2 字母）
    let mut k = i + 1
    while k < n && s.unsafe_get(k).to_int() != 93 {
      if mol_is_upper(s.unsafe_get(k)) {
        // 尝试两字母
        if k + 1 < n && mol_is_lower(s.unsafe_get(k + 1)) {
          let sb = StringBuilder()
          sb.write_char(s.unsafe_get(k).unsafe_to_char())
          sb.write_char(s.unsafe_get(k + 1).unsafe_to_char())
          return (sb.to_string(), false, k + 2 - i)
        }
        let sb = StringBuilder()
        sb.write_char(s.unsafe_get(k).unsafe_to_char())
        return (sb.to_string(), false, k + 1 - i)
      }
      k = k + 1
    }
    ("C", false, 1)
  } else {
    ("C", false, 1)
  }
}

///|
fn mol_upper_symbol(c : UInt16) -> String {
  let upper = if c >= 97 && c <= 122 { c - 32 } else { c }
  let sb = StringBuilder()
  sb.write_char(upper.unsafe_to_char())
  sb.to_string()
}

///|
/// Parse a SMILES string into a molecular graph (core grammar).
pub fn parse_smiles(smiles : String) -> Molecule {
  let atoms : Array[MolAtom] = Array(capacity=0)
  let bonds : Array[MolBond] = Array(capacity=0)
  let ring : Map[Int, Int] = Map([], capacity=0)
  let branch_stack : Array[Int] = Array(capacity=0)
  let mut prev = -1
  let mut pending_order = 1

  let n = smiles.length()
  let mut i = 0
  while i < n {
    let c = smiles.unsafe_get(i)
    let ci = c.to_int()
    if mol_is_upper(c) || mol_is_lower(c) || ci == 91 {
      let (symbol, aromatic, len) = mol_parse_element(smiles, i)
      let atom_idx = atoms.length()
      atoms.push(MolAtom::{ symbol, charge: 0, aromatic })
      if prev >= 0 {
        bonds.push(MolBond::{ from: prev, to: atom_idx, order: pending_order })
      }
      prev = atom_idx
      pending_order = 1
      i = i + len
    } else if ci == 40 {
      // (
      branch_stack.push(prev)
      i = i + 1
    } else if ci == 41 {
      // )
      if branch_stack.length() > 0 {
        prev = branch_stack[branch_stack.length() - 1]
        let _ = branch_stack.pop()
      }
      i = i + 1
    } else if ci == 61 {
      // =
      pending_order = 2
      i = i + 1
    } else if ci == 35 {
      // #
      pending_order = 3
      i = i + 1
    } else if ci == 45 {
      // -
      pending_order = 1
      i = i + 1
    } else if mol_is_digit(c) {
      let rn = mol_digit_value(c)
      match ring.get(rn) {
        Some(open_atom) => {
          if prev >= 0 && prev != open_atom {
            bonds.push(MolBond::{ from: prev, to: open_atom, order: pending_order })
          }
          ring.remove(rn)
        }
        None => {
          if prev >= 0 {
            ring.set(rn, prev)
          }
        }
      }
      pending_order = 1
      i = i + 1
    } else {
      i = i + 1
    }
  }
  Molecule::{ atoms, bonds }
}


///|
/// FNV-1a 32-bit hash over a string's code units, reduced to a non-negative
/// 31-bit Int (Int64 intermediate avoids overflow).
fn mol_string_hash(s : String) -> Int {
  let mut h : Int64 = 2166136261L
  let n = s.length()
  let mut i = 0
  while i < n {
    h = h ^ s.unsafe_get(i).to_int().to_int64()
    h = h * 16777619L
    i = i + 1
  }
  h.to_int() & 0x7fffffff
}

///|
fn mol_insertion_sort_strings(arr : Array[String]) -> Array[String] {
  let out = arr.copy()
  let mut i = 1
  while i < out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && mol_str_gt(out[j], key) {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
    i = i + 1
  }
  out
}

///|
fn mol_str_gt(a : String, b : String) -> Bool {
  // 按字节序比较（短者为 MT 前缀？这里用最短优先的简单 lexicographic）
  let na = a.length()
  let nb = b.length()
  let m = if na < nb { na } else { nb }
  let mut i = 0
  while i < m {
    let ca = a.unsafe_get(i)
    let cb = b.unsafe_get(i)
    if ca > cb {
      return true
    } else if ca < cb {
      return false
    }
    i = i + 1
  }
  na > nb
}

///|
/// Morgan (circular / ECFP) fingerprint.
///
/// Returns the multiset of environment hashes: one hash per atom at every
/// radius from 0 up to `radius`. Identical molecular graphs yield identical
/// multisets; different (sub)structures diverge as the radius grows.
pub fn morgan_fingerprint(mol : Molecule, radius : Int) -> Array[Int] {
  let n = mol.atoms.length()
  if n == 0 {
    return []
  }

  // 邻接表（邻居 + 键序）
  let neighbors : Array[Array[(Int, Int)]] = Array(capacity=0)
  for _k = 0; _k < n; _k = _k + 1 {
    neighbors.push([])
  }
  for b in mol.bonds {
    neighbors[b.from].push((b.to, b.order))
    neighbors[b.to].push((b.from, b.order))
  }

  // 初始 invariant 字符串
  let mut atom_inv : Array[String] = Array(capacity=0)
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    let deg = neighbors[k].length()
    let arom = if a.aromatic { "A" } else { "S" }
    atom_inv.push(
      a.symbol + "|" + deg.to_string() + "|" + a.charge.to_string() + "|" + arom,
    )
  }

  let fingerprint : Array[Int] = Array(capacity=0)
  for k = 0; k < n; k = k + 1 {
    fingerprint.push(mol_string_hash(atom_inv[k]))
  }

  let mut r = 0
  while r < radius {
    let new_inv : Array[String] = Array(capacity=0)
    for k = 0; k < n; k = k + 1 {
      let parts : Array[String] = Array(capacity=0)
      for (nb, order) in neighbors[k] {
        parts.push(order.to_string() + ":" + atom_inv[nb])
      }
      let sorted = mol_insertion_sort_strings(parts)
      let sb = StringBuilder()
      sb.write_string(atom_inv[k])
      sb.write_char('(')
      for p in sorted {
        sb.write_string(p)
        sb.write_char(',')
      }
      sb.write_char(')')
      let sig = sb.to_string()
      new_inv.push(sig)
      fingerprint.push(mol_string_hash(sig))
    }
    atom_inv = new_inv
    r = r + 1
  }

  fingerprint
}

// ============================================================================
// Adjacency list, distance matrix, ring perception
// ============================================================================

///|
/// Build an adjacency list: neighbors[i] = array of (neighbor, bond_order).
pub fn mol_adjacency(mol : Molecule) -> Array[Array[(Int, Int)]] {
  let n = mol.atoms.length()
  let neighbors : Array[Array[(Int, Int)]] = Array(capacity=0)
  for _k = 0; _k < n; _k = _k + 1 {
    neighbors.push([])
  }
  for b in mol.bonds {
    neighbors[b.from].push((b.to, b.order))
    neighbors[b.to].push((b.from, b.order))
  }
  neighbors
}

///|
/// BFS shortest path from a single source. Returns dist array (-1 = unreachable).
fn mol_bfs_single(
  neighbors : Array[Array[(Int, Int)]],
  src : Int,
  n : Int,
) -> Array[Int] {
  let dist = Array::make(n, -1)
  dist[src] = 0
  let queue : Array[Int] = [src]
  let mut head = 0
  while head < queue.length() {
    let u = queue[head]
    head = head + 1
    for (v, _) in neighbors[u] {
      if dist[v] == -1 {
        dist[v] = dist[u] + 1
        queue.push(v)
      }
    }
  }
  dist
}

///|
/// Topological distance matrix (all-pairs shortest path via BFS).
/// Returns n x n matrix; dist[i][j] = number of bonds between atoms i and j.
/// -1 means unreachable (disconnected graph).
pub fn mol_distance_matrix(mol : Molecule) -> Array[Array[Int]] {
  let n = mol.atoms.length()
  let neighbors = mol_adjacency(mol)
  let matrix : Array[Array[Int]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    matrix.push(mol_bfs_single(neighbors, i, n))
  }
  matrix
}

///|
/// Get the bond order between two atoms (0 if not bonded).
pub fn mol_bond_order(mol : Molecule, a : Int, b : Int) -> Int {
  for bond in mol.bonds {
    if (bond.from == a && bond.to == b) || (bond.from == b && bond.to == a) {
      return bond.order
    }
  }
  0
}

///|
/// Degree (number of neighbors) of an atom.
pub fn mol_degree(mol : Molecule, atom : Int) -> Int {
  let mut count = 0
  for b in mol.bonds {
    if b.from == atom || b.to == atom {
      count = count + 1
    }
  }
  count
}

///|
/// Find all unique rings via cycle basis (SSSR-lite).
/// For each non-tree edge, find the fundamental cycle.
/// Returns array of rings, each ring is an array of atom indices.
pub fn mol_find_rings(mol : Molecule) -> Array[Array[Int]] {
  let n = mol.atoms.length()
  let neighbors = mol_adjacency(mol)
  let parent = Array::make(n, -1)
  let visited = Array::make(n, false)
  let tree_edges : Array[(Int, Int)] = Array(capacity=0)
  let non_tree_edges : Array[(Int, Int)] = Array(capacity=0)
  let queue : Array[Int] = Array(capacity=0)
  for start = 0; start < n; start = start + 1 {
    if visited[start] {
      continue
    }
    visited[start] = true
    queue.push(start)
    let mut head = 0
    while head < queue.length() {
      let u = queue[head]
      head = head + 1
      for (v, _) in neighbors[u] {
        if not(visited[v]) {
          visited[v] = true
          parent[v] = u
          tree_edges.push((u, v))
          queue.push(v)
        } else if v > u {
          let mut is_tree = false
          for (te1, te2) in tree_edges {
            if (te1 == u && te2 == v) || (te1 == v && te2 == u) {
              is_tree = true
              break
            }
          }
          if not(is_tree) {
            non_tree_edges.push((u, v))
          }
        }
      }
    }
  }
  let rings : Array[Array[Int]] = Array(capacity=0)
  for (u, v) in non_tree_edges {
    let path = mol_tree_path(parent, u, v)
    if path.length() > 2 {
      rings.push(path)
    }
  }
  rings
}

///|
/// Find path between two nodes in a tree via parent pointers.
fn mol_tree_path(
  parent : Array[Int],
  u : Int,
  v : Int,
) -> Array[Int] {
  let path_u : Array[Int] = Array(capacity=0)
  let mut cur = u
  while cur != -1 {
    path_u.push(cur)
    cur = parent[cur]
  }
  let path_v : Array[Int] = Array(capacity=0)
  cur = v
  while cur != -1 {
    path_v.push(cur)
    cur = parent[cur]
  }
  let set_u : Map[Int, Bool] = Map([], capacity=0)
  for a in path_u {
    set_u[a] = true
  }
  let mut lca = -1
  let mut lca_v_idx = -1
  for i = 0; i < path_v.length(); i = i + 1 {
    match set_u.get(path_v[i]) {
      Some(_) => {
        lca = path_v[i]
        lca_v_idx = i
        break
      }
      None => ()
    }
  }
  if lca == -1 {
    return [u, v]
  }
  let result : Array[Int] = Array(capacity=0)
  let mut idx_u = 0
  while path_u[idx_u] != lca {
    result.push(path_u[idx_u])
    idx_u = idx_u + 1
  }
  result.push(lca)
  let mut j = lca_v_idx - 1
  while j >= 0 {
    result.push(path_v[j])
    j = j - 1
  }
  result
}

///|
/// Check if an atom is in a ring.
pub fn mol_is_in_ring(mol : Molecule, atom : Int) -> Bool {
  let rings = mol_find_rings(mol)
  for ring in rings {
    for a in ring {
      if a == atom {
        return true
      }
    }
  }
  false
}

///|
/// Count the number of rings (ring count = number of fundamental cycles).
pub fn mol_ring_count(mol : Molecule) -> Int {
  mol_find_rings(mol).length()
}

///|
/// Check if a ring is aromatic (all atoms aromatic or Huckel 4n+2 rule).
fn mol_is_aromatic_ring(
  mol : Molecule,
  ring : Array[Int],
) -> Bool {
  let mut all_aromatic = true
  for a in ring {
    if not(mol.atoms[a].aromatic) {
      all_aromatic = false
    }
  }
  if all_aromatic {
    return true
  }
  let ring_size = ring.length()
  if ring_size < 3 || ring_size > 7 {
    return false
  }
  let mut pi_electrons = 0
  for i = 0; i < ring_size; i = i + 1 {
    let atom = ring[i]
    let a = mol.atoms[atom]
    let deg = mol_degree(mol, atom)
    if a.aromatic {
      pi_electrons += 1
    } else if a.symbol == "N" && deg == 3 {
      pi_electrons += 2
    } else if a.symbol == "O" && deg == 2 {
      pi_electrons += 2
    } else if a.symbol == "S" && deg == 2 {
      pi_electrons += 2
    } else if a.symbol == "N" && deg == 2 {
      pi_electrons += 1
    } else if a.symbol == "C" {
      pi_electrons += 1
    } else {
      pi_electrons += 1
    }
  }
  let mut n = 0
  while 4 * n + 2 <= pi_electrons {
    if 4 * n + 2 == pi_electrons {
      return true
    }
    n = n + 1
  }
  false
}

///|
/// Count aromatic rings.
pub fn mol_aromatic_ring_count(mol : Molecule) -> Int {
  let rings = mol_find_rings(mol)
  let mut count = 0
  for ring in rings {
    if mol_is_aromatic_ring(mol, ring) {
      count = count + 1
    }
  }
  count
}

///|
/// Count saturated rings (all single bonds).
pub fn mol_saturated_ring_count(mol : Molecule) -> Int {
  let rings = mol_find_rings(mol)
  let mut count = 0
  for ring in rings {
    let mut all_single = true
    let rs = ring.length()
    for i = 0; i < rs; i = i + 1 {
      let j = (i + 1) % rs
      let bo = mol_bond_order(mol, ring[i], ring[j])
      if bo != 1 && bo != 0 {
        all_single = false
      }
    }
    if all_single {
      count = count + 1
    }
  }
  count
}

// ============================================================================
// Kekulization / Aromaticity perception
// ============================================================================

///|
/// Perceive aromaticity: mark atoms as aromatic if in aromatic rings (Huckel 4n+2).
pub fn mol_perceive_aromaticity(mol : Molecule) -> Molecule {
  let rings = mol_find_rings(mol)
  let new_atoms : Array[MolAtom] = Array(capacity=0)
  for a in mol.atoms {
    new_atoms.push(MolAtom::{ symbol: a.symbol, charge: a.charge, aromatic: a.aromatic })
  }
  for ring in rings {
    if mol_is_aromatic_ring(mol, ring) {
      for a_idx in ring {
        new_atoms[a_idx] = MolAtom::{
          symbol: new_atoms[a_idx].symbol,
          charge: new_atoms[a_idx].charge,
          aromatic: true,
        }
      }
    }
  }
  Molecule::{ atoms: new_atoms, bonds: mol.bonds }
}

///|
/// Kekulize: assign alternating single/double bonds to aromatic rings.
pub fn mol_kekulize(mol : Molecule) -> Molecule {
  let rings = mol_find_rings(mol)
  let new_bonds : Array[MolBond] = Array(capacity=0)
  for b in mol.bonds {
    new_bonds.push(MolBond::{ from: b.from, to: b.to, order: b.order })
  }
  for ring in rings {
    if not(mol_is_aromatic_ring(mol, ring)) {
      continue
    }
    let rs = ring.length()
    for i = 0; i < rs; i = i + 1 {
      let j = (i + 1) % rs
      let a1 = ring[i]
      let a2 = ring[j]
      for bi = 0; bi < new_bonds.length(); bi = bi + 1 {
        let b = new_bonds[bi]
        if (b.from == a1 && b.to == a2) || (b.from == a2 && b.to == a1) {
          if i % 2 == 0 {
            new_bonds[bi] = MolBond::{ from: b.from, to: b.to, order: 2 }
          } else {
            new_bonds[bi] = MolBond::{ from: b.from, to: b.to, order: 1 }
          }
          break
        }
      }
    }
  }
  let new_atoms : Array[MolAtom] = Array(capacity=0)
  for a in mol.atoms {
    new_atoms.push(MolAtom::{ symbol: a.symbol, charge: a.charge, aromatic: false })
  }
  Molecule::{ atoms: new_atoms, bonds: new_bonds }
}

// ============================================================================
// Canonical SMILES generation
// ============================================================================

///|
/// Compute Morgan canonical ranks for atoms.
fn mol_canonical_ranks(mol : Molecule) -> Array[Int] {
  let n = mol.atoms.length()
  let neighbors = mol_adjacency(mol)
  let mut invariants : Array[Int] = Array(capacity=0)
  for k = 0; k < n; k = k + 1 {
    let a = mol.atoms[k]
    let deg = neighbors[k].length()
    let arom_bit = if a.aromatic { 1 } else { 0 }
    let elem_val = match a.symbol {
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
      "B" => 5
      _ => 99
    }
    invariants.push(elem_val * 1000 + deg * 10 + a.charge * 100 + arom_bit)
  }
  let mut changed = true
  let mut iter = 0
  while changed && iter < n + 10 {
    changed = false
    let new_inv : Array[Int] = Array(capacity=0)
    for k = 0; k < n; k = k + 1 {
      let neighbor_invs : Array[Int] = Array(capacity=0)
      for (nb, order) in neighbors[k] {
        neighbor_invs.push(invariants[nb] * 10 + order)
      }
      for i = 0; i < neighbor_invs.length(); i = i + 1 {
        for j = i + 1; j < neighbor_invs.length(); j = j + 1 {
          if neighbor_invs[j] < neighbor_invs[i] {
            let tmp = neighbor_invs[i]
            neighbor_invs[i] = neighbor_invs[j]
            neighbor_invs[j] = tmp
          }
        }
      }
      let mut combined = invariants[k] * 31
      for ni in neighbor_invs {
        combined = combined * 17 + ni
      }
      new_inv.push(combined)
    }
    for k = 0; k < n; k = k + 1 {
      if new_inv[k] != invariants[k] {
        changed = true
      }
    }
    invariants = new_inv
    iter = iter + 1
  }
  let indices : Array[Int] = Array(capacity=0)
  for k = 0; k < n; k = k + 1 {
    indices.push(k)
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if invariants[indices[j]] > invariants[indices[i]] {
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
    }
  }
  let ranks = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    ranks[indices[i]] = i
  }
  ranks
}

///|
/// SMILES symbol for an atom (lowercase if aromatic).
fn mol_smiles_symbol(a : MolAtom) -> String {
  if a.aromatic {
    let mut s = a.symbol
    if s.length() == 1 {
      let c = s[0].unsafe_to_char().to_int()
      if c >= 65 && c <= 90 {
        let lower = (c + 32).unsafe_to_char()
        s = lower.to_string()
      }
    }
    s
  } else {
    a.symbol
  }
}

///|
/// Bond symbol for SMILES output.
fn mol_bond_symbol(order : Int, aromatic : Bool) -> String {
  if aromatic {
    ""
  } else {
    match order {
      1 => ""
      2 => "="
      3 => "#"
      _ => ""
    }
  }
}
///|
/// Recursive DFS helper for SMILES generation.
fn mol_smiles_dfs(
  mol : Molecule,
  atom : Int,
  parent : Int,
  visited : Array[Bool],
  neighbors : Array[Array[(Int, Int)]],
  ranks : Array[Int],
  ring_num : Map[(Int, Int), Int],
  ring_counter : Ref[Int],
  result : StringBuilder,
) -> Unit {
  visited[atom] = true
  if parent >= 0 {
    let bo = mol_bond_order(mol, parent, atom)
    let arom = mol.atoms[parent].aromatic && mol.atoms[atom].aromatic
    result.write_string(mol_bond_symbol(bo, arom))
  }
  result.write_string(mol_smiles_symbol(mol.atoms[atom]))
  // Collect ring closures and unvisited neighbors
  let ring_neighbors : Array[Int] = Array(capacity=0)
  let unvisited : Array[Int] = Array(capacity=0)
  for (nb, _) in neighbors[atom] {
    if visited[nb] && nb != parent {
      ring_neighbors.push(nb)
    } else if not(visited[nb]) {
      unvisited.push(nb)
    }
  }
  // Handle ring closures
  for rn in ring_neighbors {
    let key = if atom < rn { (atom, rn) } else { (rn, atom) }
    match ring_num.get(key) {
      Some(num) => result.write_string(num.to_string())
      None => {
        ring_num[key] = ring_counter.val
        result.write_string(ring_counter.val.to_string())
        ring_counter.val = ring_counter.val + 1
      }
    }
  }
  // Sort unvisited by rank (ascending = best first)
  for i = 0; i < unvisited.length(); i = i + 1 {
    for j = i + 1; j < unvisited.length(); j = j + 1 {
      if ranks[unvisited[j]] < ranks[unvisited[i]] {
        let tmp = unvisited[i]
        unvisited[i] = unvisited[j]
        unvisited[j] = tmp
      }
    }
  }
  // Visit first child inline, rest as branches
  for i = 0; i < unvisited.length(); i = i + 1 {
    if i > 0 {
      result.write_string("(")
    }
    mol_smiles_dfs(mol, unvisited[i], atom, visited, neighbors, ranks, ring_num, ring_counter, result)
    if i > 0 {
      result.write_string(")")
    }
  }
}

///|
/// Generate canonical SMILES from a Molecule.
/// Uses Morgan invariant-based atom ordering + DFS traversal.
pub fn mol_to_smiles(mol : Molecule) -> String {
  let n = mol.atoms.length()
  if n == 0 {
    return ""
  }
  if n == 1 {
    return mol_smiles_symbol(mol.atoms[0])
  }
  let neighbors = mol_adjacency(mol)
  let ranks = mol_canonical_ranks(mol)
  let mut start = 0
  for k = 1; k < n; k = k + 1 {
    if ranks[k] < ranks[start] {
      start = k
    }
  }
  let visited = Array::make(n, false)
  let ring_num : Map[(Int, Int), Int] = Map([], capacity=0)
  let ring_counter = Ref::new(1)
  let result : StringBuilder = StringBuilder()
  mol_smiles_dfs(mol, start, -1, visited, neighbors, ranks, ring_num, ring_counter, result)
  result.to_string()
}
