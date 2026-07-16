///|
/// NeighborSearch - KD-tree based neighbor search for protein structures.
/// Inspired by Biopython Bio.PDB.NeighborSearch.
/// Implements KD-tree for efficient spatial neighbor queries.

///|
/// A node in the KD-tree.
struct KDNode {
  point : Vector3
  atom_index : Int
  left : KDNode?
  right : KDNode?
  axis : Int
}

///|
/// NeighborSearch class for finding atoms within a certain distance.
pub struct NeighborSearch {
  atoms : Array[Atom]
  root : KDNode?
  atom_count : Int
}

///|
/// Create a new NeighborSearch from a list of atoms.
pub fn NeighborSearch::new(atoms : Array[Atom]) -> NeighborSearch {
  let n = atoms.length()
  let indices : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    indices.push(i)
  }
  let root = build_kd_tree(atoms, indices, 0)
  NeighborSearch::{
    atoms,
    root,
    atom_count: n,
  }
}

///|
/// Build a KD-tree from atom coordinates.
fn build_kd_tree(atoms : Array[Atom], indices : Array[Int], depth : Int) -> KDNode? {
  let n = indices.length()
  if n == 0 {
    return None
  }

  let axis = depth % 3

  let sorted_indices = sort_indices_by_axis(atoms, indices, axis)

  let median = n / 2

  let left_indices : Array[Int] = Array::new()
  let right_indices : Array[Int] = Array::new()

  for i = 0; i < median; i = i + 1 {
    left_indices.push(sorted_indices[i])
  }
  for i = median + 1; i < n; i = i + 1 {
    right_indices.push(sorted_indices[i])
  }

  let mid_idx = sorted_indices[median]
  let atom = atoms[mid_idx]

  let left = build_kd_tree(atoms, left_indices, depth + 1)
  let right = build_kd_tree(atoms, right_indices, depth + 1)

  Some(KDNode::{
    point: atom.coord,
    atom_index: mid_idx,
    left,
    right,
    axis,
  })
}

///|
/// Sort indices by the value on a given axis.
fn sort_indices_by_axis(
  atoms : Array[Atom], indices : Array[Int], axis : Int,
) -> Array[Int] {
  let n = indices.length()
  let result : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    result.push(indices[i])
  }

  let mut i = 1
  while i < n {
    let mut j = i
    while j > 0 {
      let a = get_coord(atoms[result[j - 1]].coord, axis)
      let b = get_coord(atoms[result[j]].coord, axis)
      if a > b {
        let tmp = result[j]
        result[j] = result[j - 1]
        result[j - 1] = tmp
      } else {
        break
      }
      j = j - 1
    }
    i = i + 1
  }

  result
}

///|
/// Get coordinate value for a given axis (0=x, 1=y, 2=z).
fn get_coord(v : Vector3, axis : Int) -> Double {
  if axis == 0 {
    v.x
  } else if axis == 1 {
    v.y
  } else {
    v.z
  }
}

///|
/// Search for atoms within a given radius from a point.
pub fn NeighborSearch::search(
  self : NeighborSearch, center : Vector3, radius : Double,
) -> Array[Atom] {
  let result_indices : Array[Int] = Array::new()
  match self.root {
    Some(root) => {
      radius_search(root, center, radius, result_indices)
    }
    None => ()
  }
  let result : Array[Atom] = Array::new()
  for idx in result_indices {
    result.push(self.atoms[idx])
  }
  result
}

///|
/// Recursive radius search in KD-tree.
fn radius_search(
  node : KDNode,
  center : Vector3,
  radius : Double,
  result : Array[Int],
) -> Unit {
  let dist = node.point.distance(center)
  if dist <= radius {
    result.push(node.atom_index)
  }

  let axis = node.axis
  let center_coord = get_coord(center, axis)
  let node_coord = get_coord(node.point, axis)
  let diff = center_coord - node_coord

  if diff <= radius {
    match node.right {
      Some(right) => {
        radius_search(right, center, radius, result)
      }
      None => ()
    }
  }

  if diff >= -radius {
    match node.left {
      Some(left) => {
        radius_search(left, center, radius, result)
      }
      None => ()
    }
  }
}

///|
/// Search for atoms within a given radius from a query atom.
pub fn NeighborSearch::search_atom(
  self : NeighborSearch, query : Atom, radius : Double,
) -> Array[Atom] {
  self.search(query.coord, radius)
}

///|
/// Find all pairs of atoms within a given distance.
pub fn NeighborSearch::search_all(self : NeighborSearch, radius : Double) -> Array[(Atom, Atom)] {
  let result : Array[(Atom, Atom)] = Array::new()
  let n = self.atoms.length()
  for i = 0; i < n; i = i + 1 {
    let atom = self.atoms[i]
    let neighbors = self.search(atom.coord, radius)
    for neighbor in neighbors {
      if atom.name != neighbor.name || atom.resname != neighbor.resname ||
         atom.resseq != neighbor.resseq || atom.chainid != neighbor.chainid {
        result.push((atom, neighbor))
      }
    }
  }
  result
}

///|
/// Get the number of atoms in the search structure.
pub fn NeighborSearch::size(self : NeighborSearch) -> Int {
  self.atom_count
}

///|
/// Find the nearest neighbor to a point.
pub fn NeighborSearch::nearest(
  self : NeighborSearch, center : Vector3,
) -> Atom? {
  if self.atoms.length() == 0 {
    return None
  }
  let mut best_idx = 0
  let mut best_dist = self.atoms[0].coord.distance(center)
  let mut i = 1
  while i < self.atoms.length() {
    let d = self.atoms[i].coord.distance(center)
    if d < best_dist {
      best_dist = d
      best_idx = i
    }
    i = i + 1
  }
  Some(self.atoms[best_idx])
}
