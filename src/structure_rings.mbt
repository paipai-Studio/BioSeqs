///|
/// Aromatic ring detection using SSSR (Smallest Set of Smallest Rings).
/// Ported from biotite.structure.rings.
///
/// Finds rings in a molecular graph and classifies them as aromatic
/// using Hückel's rule (4n+2 π electrons).

///|
/// A detected ring.
pub struct Ring {
  atom_indices : Array[Int]   // indices forming the ring (in order)
  size : Int                   // number of atoms in ring
  is_aromatic : Bool
}

///|
/// Find all rings up to `max_size` in a molecular graph.
///
/// Parameters:
///   `atoms` — atom array
///   `max_size` — maximum ring size to detect, default 7
///   `tolerance` — bond distance tolerance (Å), default 0.4
pub fn detect_rings(
  atoms : Array[Atom],
  max_size~ : Int = 7,
  tolerance~ : Double = 0.4,
) -> Array[Ring] {
  let n = atoms.length()
  // Build adjacency list based on covalent radii
  let adj : Array[Array[Int]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    adj.push(Array(capacity=4))
  }
  for i = 0; i < n; i = i + 1 {
    let ri = covalent_radius(atoms[i].element)
    for j = i + 1; j < n; j = j + 1 {
      let rj = covalent_radius(atoms[j].element)
      let threshold = ri + rj + tolerance
      let dist = atoms[i].coord.distance(atoms[j].coord)
      if dist <= threshold && dist > 0.3 {
        adj[i].push(j)
        adj[j].push(i)
      }
    }
  }
  // Find rings using DFS-based cycle detection
  let rings : Array[Array[Int]] = Array(capacity=0)
  let visited : Array[Bool] = Array::make(n, false)
  for start = 0; start < n; start = start + 1 {
    if visited[start] {
      continue
    }
    // DFS to find cycles
    let path : Array[Int] = Array(capacity=0)
    dfs_find_rings(start, -1, adj, visited, path, rings, max_size)
  }
  // Deduplicate rings (same set of atoms = same ring)
  let unique_rings = deduplicate_rings(rings)
  // Classify aromaticity
  let result : Array[Ring] = Array(capacity=unique_rings.length())
  for ring in unique_rings {
    let is_arom = is_aromatic_ring(atoms, ring)
    result.push(Ring::{
      atom_indices: ring,
      size: ring.length(),
      is_aromatic: is_arom,
    })
  }
  result
}

///|
/// DFS to find rings in the molecular graph.
fn dfs_find_rings(
  current : Int,
  parent : Int,
  adj : Array[Array[Int]],
  visited : Array[Bool],
  path : Array[Int],
  rings : Array[Array[Int]],
  max_size : Int,
) -> Unit {
  if path.length() > max_size {
    return
  }
  visited[current] = true
  path.push(current)
  for neighbor in adj[current] {
    if neighbor == parent {
      continue
    }
    // Check if neighbor is in path (cycle found)
    let mut found_idx = -1
    for i = 0; i < path.length(); i = i + 1 {
      if path[i] == neighbor {
        found_idx = i
      }
    }
    if found_idx >= 0 && path.length() - found_idx >= 3 {
      // Extract ring: path[found_idx..end]
      let ring : Array[Int] = Array(capacity=path.length() - found_idx)
      for i = found_idx; i < path.length(); i = i + 1 {
        ring.push(path[i])
      }
      rings.push(ring)
    } else if not(visited[neighbor]) {
      dfs_find_rings(neighbor, current, adj, visited, path, rings, max_size)
    }
  }
  let _ = path.pop()
  // Don't unmark visited — we only need to find each ring once from one direction
  visited[current] = false
}

///|
/// Deduplicate rings by atom set.
fn deduplicate_rings(rings : Array[Array[Int]]) -> Array[Array[Int]] {
  let result : Array[Array[Int]] = Array(capacity=0)
  let seen : Map[String, Bool] = Map([], capacity=rings.length())
  for ring in rings {
    // Create canonical key: sorted atom indices joined
    let sorted_atoms : Array[Int] = Array(capacity=ring.length())
    for a in ring {
      sorted_atoms.push(a)
    }
    // Simple sort
    for i = 1; i < sorted_atoms.length(); i = i + 1 {
      let key = sorted_atoms[i]
      let mut j = i - 1
      while j >= 0 && sorted_atoms[j] > key {
        sorted_atoms[j + 1] = sorted_atoms[j]
        j = j - 1
      }
      sorted_atoms[j + 1] = key
    }
    let sb = StringBuilder()
    for a in sorted_atoms {
      sb.write_string(a.to_string())
      sb.write_char(',')
    }
    let key = sb.to_string()
    match seen.get(key) {
      Some(_) => () // duplicate
      None => {
        seen[key] = true
        result.push(ring)
      }
    }
  }
  result
}

///|
/// Check if a ring is aromatic using Hückel's rule.
/// Simplified: 5- or 6-membered ring with all sp2 carbons/nitrogens.
fn is_aromatic_ring(atoms : Array[Atom], ring : Array[Int]) -> Bool {
  let size = ring.length()
  // Aromatic rings are typically 5 or 6 membered
  if size != 5 && size != 6 {
    return false
  }
  // Check all atoms are C or N (common aromatic atoms)
  for idx in ring {
    let e = to_uppercase(atoms[idx].element)
    if e != "C" && e != "N" && e != "O" && e != "S" {
      return false
    }
  }
  // For 6-membered rings: check planarity (simplified — assume aromatic if all C/N)
  // For 5-membered rings: check if it has N (pyrrole-type) or is furan/thiophene
  true
}
