///|
/// Molecular graph operations for AtomArray.
///
/// Ring detection, distance matrices, contact maps, connected
/// components — all built on top of the existing BondList.

// ============================================================
// P3-1: Ring detection on AtomArray
// ============================================================

///|
/// Detect all rings up to `max_size` in the AtomArray.
///
/// Re-implements a simple SSSR-style ring detection for AtomArray directly.
pub fn AtomArray::detect_rings(
  self : AtomArray,
  max_size~ : Int32 = 7,
) -> Array[Ring] {
  let n = self.len()
  if n == 0 {
    return []
  }
  // Build adjacency via bond inference
  let bl = bondlist_from_coordinates(self.coords, self.elements)
  let adj : Array[Array[Int32]] = []
  let mut idx = 0
  while idx < n {
    adj.push(bl.get_neighbors(idx))
    idx = idx + 1
  }
  // DFS ring detection
  let rings : Array[Array[Int32]] = Array(capacity=0)
  let visited = Array::make(n, false)
  let mut start = 0
  while start < n {
    if !visited[start] {
      let path : Array[Int32] = Array(capacity=max_size)
      sg_dfs_find_rings(start, -1, adj, visited, path, rings, max_size)
    }
    start = start + 1
  }
  let unique = sg_deduplicate_rings(rings)
  let result : Array[Ring] = Array(capacity=unique.length())
  for r in unique {
    let is_arom = sg_is_aromatic(self, r)
    result.push(Ring::{
      atom_indices: r,
      size: r.length(),
      is_aromatic: is_arom,
    })
  }
  result
}

///|
/// DFS ring finder helper.
fn sg_dfs_find_rings(
  current : Int32,
  parent : Int32,
  adj : Array[Array[Int32]],
  visited : Array[Bool],
  path : Array[Int32],
  rings : Array[Array[Int32]],
  max_size : Int32,
) -> Unit {
  if path.length() >= max_size {
    return
  }
  visited[current] = true
  path.push(current)
  for neighbor in adj[current] {
    if neighbor == parent {
      continue
    }
    let mut found_idx = -1
    let mut k = 0
    while k < path.length() {
      if path[k] == neighbor {
        found_idx = k
      }
      k = k + 1
    }
    if found_idx >= 0 && path.length() - found_idx >= 3 {
      let ring : Array[Int32] = Array(capacity=path.length() - found_idx)
      let mut a = found_idx
      while a < path.length() {
        ring.push(path[a])
        a = a + 1
      }
      rings.push(ring)
    } else if !visited[neighbor] {
      sg_dfs_find_rings(neighbor, current, adj, visited, path, rings, max_size)
    }
  }
  ignore(path.pop())
  visited[current] = false
}

///|
/// Deduplicate rings by sorted atom set key.
fn sg_deduplicate_rings(rings : Array[Array[Int32]]) -> Array[Array[Int32]] {
  let result : Array[Array[Int32]] = Array(capacity=0)
  let seen : Map[String, Bool] = Map([], capacity=rings.length())
  for ring in rings {
    // Sort (insertion)
    let sorted : Array[Int32] = Array(capacity=ring.length())
    for a in ring {
      sorted.push(a)
    }
    let mut i = 1
    while i < sorted.length() {
      let key = sorted[i]
      let mut j = i - 1
      while j >= 0 && sorted[j] > key {
        sorted[j + 1] = sorted[j]
        j = j - 1
      }
      sorted[j + 1] = key
      i = i + 1
    }
    let sb = StringBuilder()
    for a in sorted {
      sb.write_string(a.to_string())
      sb.write_char(',')
    }
    let k = sb.to_string()
    match seen.get(k) {
      Some(_) => ()
      None => {
        seen[k] = true
        result.push(ring)
      }
    }
  }
  result
}

///|
/// Simplified aromaticity check for AtomArray.
fn sg_is_aromatic(arr : AtomArray, ring : Array[Int32]) -> Bool {
  let size = ring.length()
  if size != 5 && size != 6 {
    return false
  }
  for idx in ring {
    let e = arr.elements[idx]
    if e != "C" && e != "N" && e != "O" && e != "S" {
      return false
    }
  }
  true
}

///|
/// Check if an atom is part of any ring.
pub fn AtomArray::is_atom_in_ring(self : AtomArray, atom_idx : Int32) -> Bool {
  let rings = self.detect_rings()
  for ring in rings {
    for idx in ring.atom_indices {
      if idx == atom_idx {
        return true
      }
    }
  }
  false
}

///|
/// Find rotatable single bonds excluding bonds that are part of a ring.
///
/// Biotite's canonical definition:
///   1. BondType::Single
///   2. Both atoms have >= 2 neighbors (non-terminal)
///   3. At least one of the two atoms has >= 3 neighbors (prevents
///      -CH2-CH2- in long alkanes from every bond being flagged)
///   4. The bond is not inside a ring
pub fn AtomArray::find_rotatable_bonds_strict(self : AtomArray) -> Array[Bond] {
  let bl = bondlist_from_coordinates(self.coords, self.elements)
  let rings = self.detect_rings()
  // Build a set of bonds that are inside any ring (unordered pair)
  let in_ring : Array[(Int32, Int32)] = []
  for ring in rings {
    let ra = ring.atom_indices
    let n = ra.length()
    let mut i = 0
    while i < n {
      let a = ra[i]
      let b = ra[(i + 1) % n]
      // Normalize order
      let pair = if a < b { (a, b) } else { (b, a) }
      in_ring.push(pair)
      i = i + 1
    }
  }
  let result : Array[Bond] = []
  for b in bl.get_bonds() {
    if !(b.bond_type is BondType::Single) {
      continue
    }
    let n1 = bl.get_neighbors(b.atom1).length()
    let n2 = bl.get_neighbors(b.atom2).length()
    if n1 < 2 || n2 < 2 {
      continue
    }
    // Require at least one atom with >= 3 neighbors
    if n1 < 3 && n2 < 3 {
      continue
    }
    // Ring check
    let pair = if b.atom1 < b.atom2 {
      (b.atom1, b.atom2)
    } else {
      (b.atom2, b.atom1)
    }
    let mut skip = false
    for rp in in_ring {
      if rp.0 == pair.0 && rp.1 == pair.1 {
        skip = true
        break
      }
    }
    if skip {
      continue
    }
    result.push(b)
  }
  result
}

// ============================================================
// P3-2: Distance matrix + contact map
// ============================================================

///|
/// Compute the all-pairs Euclidean distance matrix (n×n).
/// Returns Array[Array[Double]], where d[i][j] = distance(coords[i], coords[j]).
/// Symmetric, diagonal = 0.
pub fn AtomArray::distance_matrix(self : AtomArray) -> Array[Array[Double]] {
  let n = self.len()
  let result : Array[Array[Double]] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let row = Array::make(n, 0.0)
    let mut j = i + 1
    while j < n {
      let d = self.coords[i].distance(self.coords[j])
      row[j] = d
      let prev_row = result[j]
      // We can't easily write result[j][i] here since we only have capacity Array
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  // Now fill lower triangle
  let mut a = 0
  while a < n {
    let mut b = 0
    while b < a {
      let d = result[b][a]
      result[a][b] = d
      b = b + 1
    }
    a = a + 1
  }
  result
}

///|
/// Compute the Cα-only distance matrix.
/// Same as distance_matrix but filtered to CA atoms first.
pub fn AtomArray::ca_distance_matrix(self : AtomArray) -> Array[Array[Double]] {
  let cas = self.filter_ca()
  cas.distance_matrix()
}

///|
/// Compute a contact map: binary matrix where contact_map[i][j] = true
/// if distance between atom i and atom j is below `cutoff` Å.
pub fn AtomArray::contact_map(
  self : AtomArray,
  cutoff~ : Double = 5.0,
) -> Array[Array[Bool]] {
  let n = self.len()
  let result : Array[Array[Bool]] = Array(capacity=n)
  let mut i = 0
  while i < n {
    let row = Array::make(n, false)
    let mut j = i + 1
    while j < n {
      let d = self.coords[i].distance(self.coords[j])
      if d < cutoff {
        row[j] = true
      }
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  // Mirror lower triangle
  let mut a = 0
  while a < n {
    let mut b = 0
    while b < a {
      let v = result[b][a]
      result[a][b] = v
      b = b + 1
    }
    a = a + 1
  }
  result
}

///|
/// Compute the residue-level contact map.
/// Two residues are considered in contact if any pair of heavy atoms
/// (one from each residue) is within `cutoff` Å.
pub fn AtomArray::residue_contact_map(
  self : AtomArray,
  cutoff~ : Double = 5.0,
) -> Array[Array[Bool]] {
  let starts = self.get_residue_starts()
  let n_res = starts.length()
  let heavy = self.filter_heavy()
  let heavy_map = self.get_heavy_atom_indices()
  // We need per-residue heavy atom ranges
  let res_heavy_ranges : Array[(Int32, Int32)] = [] // (start_in_heavy, end_in_heavy)
  let mut cur = 0
  for s in starts {
    let mut end = cur
    while end < heavy_map.length() &&
          heavy_map[end] >= s &&
          (
            end + 1 >= heavy_map.length() ||
            heavy_map[end + 1] >=
            s +
            (if cur + 1 < starts.length() { starts[cur + 1] - s } else { 0 })
          ) {
      end = end + 1
    }
    let actual_end = if cur + 1 < starts.length() {
      // All heavy atoms in [s, starts[cur+1])
      let next_s = starts[cur + 1]
      let mut k = cur
      while k < heavy_map.length() && heavy_map[k] < next_s {
        k = k + 1
      }
      k
    } else {
      heavy_map.length()
    }
    res_heavy_ranges.push((cur, actual_end))
    cur = actual_end
  }
  // Now compare all pairs of residues
  let result : Array[Array[Bool]] = Array(capacity=n_res)
  let mut i = 0
  while i < n_res {
    let row = Array::make(n_res, false)
    let (si, ei) = res_heavy_ranges[i]
    let mut j = i + 1
    while j < n_res {
      let (sj, ej) = res_heavy_ranges[j]
      let mut has_contact = false
      let mut a = si
      while a < ei && !has_contact {
        let idx_a = heavy_map[a]
        let mut b = sj
        while b < ej && !has_contact {
          let idx_b = heavy_map[b]
          if self.coords[idx_a].distance(self.coords[idx_b]) < cutoff {
            has_contact = true
          }
          b = b + 1
        }
        a = a + 1
      }
      if has_contact {
        row[j] = true
      }
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  // Mirror
  let mut a = 0
  while a < n_res {
    let mut b = 0
    while b < a {
      let v = result[b][a]
      result[a][b] = v
      b = b + 1
    }
    a = a + 1
  }
  result
}

///|
/// Internal helper: return indices of heavy (non-H) atoms.
pub fn AtomArray::get_heavy_atom_indices(self : AtomArray) -> Array[Int32] {
  let result : Array[Int32] = []
  for i in 0..<self.len() {
    let e = self.elements[i]
    if e != "H" && e != "D" {
      result.push(i)
    }
  }
  result
}

// ============================================================
// P3-3: RMS family (AtomArray native)
// ============================================================

///|
/// Compute raw RMSD (no alignment) between two AtomArrays.
/// Assumes the same atom ordering (same length).
pub fn AtomArray::rmsd(self : AtomArray, other : AtomArray) -> Double {
  let n = self.len()
  if n == 0 || n != other.len() {
    return 0.0
  }
  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    let dx = self.coords[i].x - other.coords[i].x
    let dy = self.coords[i].y - other.coords[i].y
    let dz = self.coords[i].z - other.coords[i].z
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
    i = i + 1
  }
  (sum_sq / n.to_double()).sqrt()
}

///|
/// Compute per-atom RMSF (root mean square fluctuation) across multiple
/// models / frames. `frames` should be a stacked AtomArray (same atom
/// ordering as self).
///
/// Returns an array of per-atom RMSF values.
pub fn AtomArray::rmsf(
  self : AtomArray,
  frames : AtomArrayStack,
) -> Array[Double] {
  let n_frames = frames.n_models()
  let n_atoms = self.len()
  if n_frames == 0 || n_atoms == 0 {
    return []
  }
  // Use self as reference. Center per-frame first.
  let result : Array[Double] = Array(capacity=n_atoms)
  let mut i = 0
  while i < n_atoms {
    let mut sum_sq = 0.0
    let mut f = 0
    while f < n_frames {
      let model = frames.get_model(f)
      let dx = model.coords[i].x - self.coords[i].x
      let dy = model.coords[i].y - self.coords[i].y
      let dz = model.coords[i].z - self.coords[i].z
      sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
      f = f + 1
    }
    result.push((sum_sq / n_frames.to_double()).sqrt())
    i = i + 1
  }
  result
}

// ============================================================
// P3-4: Connected component / polymer filter
// ============================================================

///|
/// Get connected component ID for each atom (0-indexed).
/// Uses inferred bond connectivity.
pub fn AtomArray::connected_components(self : AtomArray) -> Array[Int32] {
  let bl = bondlist_from_coordinates(self.coords, self.elements)
  let n = self.len()
  let visited = Array::make(n, false)
  let labels = Array::make(n, -1)
  let mut label = 0
  let mut i = 0
  while i < n {
    if !visited[i] {
      // BFS
      let queue : Array[Int32] = [i]
      visited[i] = true
      labels[i] = label
      while queue.length() > 0 {
        let cur = queue[0]
        ignore(queue.remove(0))
        let neigh = bl.get_neighbors(cur)
        for nbr in neigh {
          if nbr >= 0 && nbr < n && !visited[nbr] {
            visited[nbr] = true
            labels[nbr] = label
            queue.push(nbr)
          }
        }
      }
      label = label + 1
    }
    i = i + 1
  }
  labels
}

///|
/// Filter to only atoms in the largest connected component (polymer).
pub fn AtomArray::filter_polymer(self : AtomArray) -> AtomArray {
  let labels = self.connected_components()
  // Find max label frequency
  let counts : Array[Int32] = []
  let mut i = 0
  while i < labels.length() {
    let l = labels[i]
    while counts.length() <= l {
      counts.push(0)
    }
    counts[l] = counts[l] + 1
    i = i + 1
  }
  let mut max_label = 0
  let mut max_count = 0
  let mut j = 0
  while j < counts.length() {
    if counts[j] > max_count {
      max_count = counts[j]
      max_label = j
    }
    j = j + 1
  }
  let mask : Array[Bool] = []
  let mut k = 0
  while k < labels.length() {
    mask.push(labels[k] == max_label)
    k = k + 1
  }
  self.filter(mask)
}
