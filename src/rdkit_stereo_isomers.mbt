/// rdkit_stereo_isomers.mbt — stereoisomer enumeration (RDKit gap fill, A1)
///
/// Built on top of mol_stereo.mbt CIP perception: potential tetrahedral
/// stereocenters (degree 3 with an implicit H, or degree 4) are detected
/// by CIP-priority distinctness of the substituents, and all R/S
/// combinations are enumerated up to `max_isomers`. Stereo-SMILES output
/// uses the [C@H]/[C@@H] convention where the atom's DFS neighbour order
/// stands in for a true 3D parity (same simplified convention as
/// mol_stereo.cip_tetrahedral).

///|
/// A single enumerated stereoisomer: center assignments, a readable
/// label ("R,S"), and an annotated SMILES.
pub struct RdkitStereoIsomer {
  assignments : Array[(Int32, Bool)] // (atom, true = R)
  label : String
  smiles : String
}

///|
fn stereo_nb_priority(mol : Molecule, nb : Int32) -> Int32 {
  cip_atomic_number(mol.atoms[nb].symbol) * 100 + mol_degree(mol, nb)
}

///|
fn stereo_all_distinct(values : Array[Int32]) -> Bool {
  for i = 0; i < values.length(); i = i + 1 {
    for j = i + 1; j < values.length(); j = j + 1 {
      if values[i] == values[j] {
        return false
      }
    }
  }
  true
}

///|
/// Find potential tetrahedral stereocenters: sp3-like atoms of degree 3
/// (implicit H) or degree 4 whose CIP-priority substituent levels are all
/// distinct (substituent distinctness is approximated by atomic number
/// and degree, matching mol_stereo's simplified CIP comparison).
pub fn rdkit_find_potential_stereocenters(mol : Molecule) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=0)
  let neighbors = mol_adjacency(mol)
  for atom = 0; atom < mol.atoms.length(); atom = atom + 1 {
    let deg = neighbors[atom].length()
    if deg != 3 && deg != 4 {
      continue
    }
    let nbs : Array[Int32] = Array(capacity=0)
    for (nb, _) in neighbors[atom] {
      nbs.push(nb)
    }
    let priorities : Array[Int32] = Array(capacity=0)
    for nb in nbs {
      priorities.push(stereo_nb_priority(mol, nb))
    }
    if stereo_all_distinct(priorities) {
      out.push(atom)
    }
  }
  out
}

///|
/// Count stereoisomers as 2^n over potential centers (1 when the
/// molecule has no potential center, matching the RDKit convention).
pub fn rdkit_count_stereoisomers(mol : Molecule) -> Int32 {
  let n = rdkit_find_potential_stereocenters(mol).length()
  1 << n
}

///|
/// Enumerate R/S combinations over the potential centers.
/// Centers with fewer than 2^max potential centers are covered fully;
/// otherwise enumeration stops at `max_isomers`.
pub fn rdkit_enumerate_stereoisomers(
  mol : Molecule,
  max_isomers~ : Int32 = 64,
) -> Array[RdkitStereoIsomer] {
  let centers = rdkit_find_potential_stereocenters(mol)
  let out : Array[RdkitStereoIsomer] = Array(capacity=0)
  let n = centers.length()
  if n == 0 {
    return out
  }
  let n_combos = 1 << n
  let limit = if n_combos > max_isomers { max_isomers } else { n_combos }
  for combo = 0; combo < limit; combo = combo + 1 {
    let assignments : Array[(Int32, Bool)] = Array(capacity=0)
    let label_parts : Array[String] = Array(capacity=0)
    for k = 0; k < n; k = k + 1 {
      let is_r = ((combo >> k) & 1) == 0
      assignments.push((centers[k], is_r))
      label_parts.push(if is_r { "R" } else { "S" })
    }
    let label = stereo_join_labels(label_parts)
    let smiles = rdkit_stereo_smiles(mol, assignments)
    out.push(RdkitStereoIsomer::{ assignments, label, smiles, })
  }
  out
}

///|
fn stereo_join_labels(parts : Array[String]) -> String {
  let sb = StringBuilder()
  for i = 0; i < parts.length(); i = i + 1 {
    if i > 0 {
      let _ = sb.write(",")
    }
    let _ = sb.write(parts[i])
  }
  sb.to_string()
}

///|
/// Generate a SMILES with tetrahedral stereo annotations for the given
/// (atom, is_R) assignments: [C@H] for R and [C@@H] for S on degree-3
/// centers (implicit H), [C@]/[C@@] on degree-4 centers. The @/@@ choice
/// is relative to the DFS neighbour order (graph parity convention).
pub fn rdkit_stereo_smiles(
  mol : Molecule,
  assignments : Array[(Int32, Bool)],
) -> String {
  let n = mol.atoms.length()
  if n == 0 {
    return ""
  }
  let is_r : Map[Int32, Bool] = Map([], capacity=0)
  for (atom, r) in assignments {
    is_r[atom] = r
  }
  let neighbors = mol_adjacency(mol)
  let ranks = mol_canonical_ranks(mol)
  let mut start = 0
  for k = 1; k < n; k = k + 1 {
    if ranks[k] < ranks[start] {
      start = k
    }
  }
  let ring_map = mol_ring_digit_map(mol, neighbors, start)
  let visited : Array[Bool] = Array::make(n, false)
  let result = StringBuilder()
  stereo_smiles_dfs(
    mol, start, -1, visited, neighbors, ranks, ring_map, is_r, result,
  )
  result.to_string()
}

///|
fn stereo_smiles_dfs(
  mol : Molecule,
  atom : Int32,
  parent : Int32,
  visited : Array[Bool],
  neighbors : Array[Array[(Int32, Int32)]],
  ranks : Array[Int32],
  ring_map : Map[Int32, Array[Int32]],
  is_r : Map[Int32, Bool],
  result : StringBuilder,
) -> Unit {
  if visited[atom] {
    return
  }
  visited[atom] = true
  if parent >= 0 {
    let bo = mol_bond_order(mol, parent, atom)
    let arom = mol.atoms[parent].aromatic && mol.atoms[atom].aromatic
    result.write_string(mol_bond_symbol(bo, arom))
  }
  // Atom symbol, with stereo annotation when assigned
  match is_r.get(atom) {
    Some(r) => {
      let deg = neighbors[atom].length()
      let sym = mol.atoms[atom].symbol
      let charge = mol.atoms[atom].charge
      let marker = if r { "@" } else { "@@" }
      result.write_string("[")
      result.write_string(sym)
      result.write_string(marker)
      if deg == 3 {
        result.write_string("H")
      }
      if charge > 0 {
        result.write_string("+")
      } else if charge < 0 {
        result.write_string("-")
      }
      result.write_string("]")
    }
    None => result.write_string(mol_smiles_symbol(mol.atoms[atom]))
  }
  // Ring closure digits touching this atom
  match ring_map.get(atom) {
    Some(nums) =>
      for rn in nums {
        result.write_string(rn.to_string())
      }
    None => ()
  }
  // Unvisited neighbours (same layout as mol_smiles_dfs)
  let unvisited : Array[Int32] = Array(capacity=0)
  for (nb, _) in neighbors[atom] {
    if !visited[nb] {
      unvisited.push(nb)
    }
  }
  for i = 0; i < unvisited.length(); i = i + 1 {
    for j = i + 1; j < unvisited.length(); j = j + 1 {
      if ranks[unvisited[j]] < ranks[unvisited[i]] {
        let tmp = unvisited[i]
        unvisited[i] = unvisited[j]
        unvisited[j] = tmp
      }
    }
  }
  for i = 0; i < unvisited.length(); i = i + 1 {
    if i > 0 {
      result.write_string("(")
    }
    stereo_smiles_dfs(
      mol,
      unvisited[i],
      atom,
      visited,
      neighbors,
      ranks,
      ring_map,
      is_r,
      result,
    )
    if i > 0 {
      result.write_string(")")
    }
  }
}
