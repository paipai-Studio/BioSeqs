///|
/// biotite_structure_pairwise.mbt — residue/atom matching and graph paths.

///|
/// Residue signature: chain + residue ID + residue name.
fn biotite_residue_signature(arr : AtomArray, i : Int32) -> String {
  arr.get_chain_id(i).to_string() +
  ":" +
  arr.res_ids[i].to_string() +
  ":" +
  arr.get_res_name(i)
}

///|
/// Atom signature: residue signature + atom name.
fn biotite_atom_signature(arr : AtomArray, i : Int32) -> String {
  biotite_residue_signature(arr, i) + ":" + arr.get_name(i)
}

///|
/// Match residues shared by two structures, returning (index in a, index in b).
pub fn biotite_find_shared_residues(
  a : AtomArray,
  b : AtomArray,
) -> Array[(Int32, Int32)] {
  let index : Map[String, Int32] = Map([])
  for j in 0..<b.len() {
    let sig = biotite_residue_signature(b, j)
    if !index.contains(sig) {
      index[sig] = j
    }
  }
  let result : Array[(Int32, Int32)] = []
  let seen : Map[String, Bool] = Map([])
  for i in 0..<a.len() {
    let sig = biotite_residue_signature(a, i)
    match index.get(sig) {
      Some(j) =>
        if !seen.contains(sig) {
          seen[sig] = true
          result.push((i, j))
        }
      None => ()
    }
  }
  result
}

///|
/// Match atoms shared by two structures, returning (index in a, index in b).
pub fn biotite_find_shared_atoms(
  a : AtomArray,
  b : AtomArray,
) -> Array[(Int32, Int32)] {
  let index : Map[String, Int32] = Map([])
  for j in 0..<b.len() {
    let sig = biotite_atom_signature(b, j)
    if !index.contains(sig) {
      index[sig] = j
    }
  }
  let result : Array[(Int32, Int32)] = []
  for i in 0..<a.len() {
    let sig = biotite_atom_signature(a, i)
    match index.get(sig) {
      Some(j) => result.push((i, j))
      None => ()
    }
  }
  result
}

///|
/// Breadth-first shortest path through a bond list.
/// Returns the vertex sequence including both endpoints, or [] if unreachable.
pub fn biotite_find_connected_path(
  bonds : Array[Biotite3Bond],
  n_atoms : Int32,
  start : Int32,
  end : Int32,
) -> Array[Int32] {
  if start < 0 || start >= n_atoms || end < 0 || end >= n_atoms {
    return []
  }
  let visited : Array[Bool] = Array(capacity=n_atoms)
  let previous : Array[Int32] = Array(capacity=n_atoms)
  for _i in 0..<n_atoms {
    visited.push(false)
    previous.push(-1)
  }
  let queue : Array[Int32] = [start]
  visited[start] = true
  let mut head = 0
  while head < queue.length() {
    let current = queue[head]
    head = head + 1
    if current == end {
      break
    }
    for bond in bonds {
      let neighbor = if bond.atom1 == current {
        bond.atom2
      } else if bond.atom2 == current {
        bond.atom1
      } else {
        -1
      }
      if neighbor >= 0 && neighbor < n_atoms && !visited[neighbor] {
        visited[neighbor] = true
        previous[neighbor] = current
        queue.push(neighbor)
      }
    }
  }
  if !visited[end] {
    return []
  }
  let path : Array[Int32] = []
  let mut node = end
  while node != -1 {
    path.push(node)
    if node == start {
      break
    }
    node = previous[node]
  }
  let reversed : Array[Int32] = Array(capacity=path.length())
  for i in 0..<path.length() {
    reversed.push(path[path.length() - 1 - i])
  }
  reversed
}
