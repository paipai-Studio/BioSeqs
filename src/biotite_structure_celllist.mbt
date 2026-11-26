/// Biotite `structure.celllist` — high-level CellList API.
///
/// Wraps the repository's grid-based spatial neighbor search with the
/// Biotite-facing operations: adjacent_atoms, nearest_neighbors,
/// conflict and create_by_boxes.

///|
/// A Biotite-style cell list built from atom positions.
pub struct BiotiteCellList {
  cell_size : Double
  inner : CellList
  positions : Array[Vector3]
}

///|
/// Build a Biotite cell list from positions.
/// `cell_size` should be at least the cutoff radius used later.
pub fn biotite_celllist_new(
  positions : Array[Vector3],
  cell_size~ : Double = 2.0,
) -> BiotiteCellList {
  let flat : Array[Array[Double]] = Array(capacity=positions.length())
  for p in positions {
    flat.push([p.x, p.y, p.z])
  }
  BiotiteCellList::{
    cell_size,
    inner: celllist_new(flat, cell_size),
    positions,
  }
}

///|
/// All atoms within `radius` of atom `index`, excluding itself.
/// Indices are returned in ascending order (Biotite order is unspecified,
/// ascending is used for determinism).
pub fn BiotiteCellList::adjacent_atoms(
  self : BiotiteCellList,
  index : Int32,
  radius : Double,
) -> Array[Int32] {
  let raw = celllist_neighbors(self.inner, index, radius)
  let out = raw.copy()
  out.sort_by(fn(a, b) { Int32::compare(a, b) })
  out
}

///|
/// The `number` nearest neighbors of atom `index`.
///
/// `distance` caps the search radius; when absent the whole cell list is
/// scanned. Returns indices sorted by ascending distance together with
/// the actual distances. Equidistant atoms are all included (all=True).
pub fn BiotiteCellList::nearest_neighbors(
  self : BiotiteCellList,
  index : Int32,
  number~ : Int32 = 1,
  distance~ : Double? = None,
) -> (Array[Int32], Array[Double]) {
  let n = self.positions.length()
  let dists : Array[(Double, Int32)] = Array(capacity=n)
  if index < 0 || index >= n {
    return ([], [])
  }
  let center = self.positions[index]
  for j = 0; j < n; j = j + 1 {
    if j == index {
      continue
    }
    let d = center.distance(self.positions[j])
    match distance {
      Some(r) => if d <= r { dists.push((d, j)) }
      None => dists.push((d, j))
    }
  }
  dists.sort_by(fn(a, b) { Double::compare(a.0, b.0) })
  let idx : Array[Int32] = []
  let ds : Array[Double] = []
  let mut last_d = -1.0
  for k = 0; k < dists.length() && idx.length() < number; k = k + 1 {
    let (d, j) = dists[k]
    if idx.length() >= number && d > last_d {
      break
    }
    idx.push(j)
    ds.push(d)
    last_d = d
  }
  (idx, ds)
}

///|
/// All pairs of atoms closer than `fixed_distance`.
/// Returns (i, j) pairs with i < j.
pub fn BiotiteCellList::conflict(
  self : BiotiteCellList,
  fixed_distance : Double,
) -> Array[(Int32, Int32)] {
  let pairs = celllist_all_pairs(self.inner, fixed_distance)
  let out : Array[(Int32, Int32)] = []
  let r2 = fixed_distance * fixed_distance
  for (i, j) in pairs {
    let p1 = self.positions[i]
    let p2 = self.positions[j]
    let dx = p1.x - p2.x
    let dy = p1.y - p2.y
    let dz = p1.z - p2.z
    if dx * dx + dy * dy + dz * dz < r2 {
      out.push((i, j))
    }
  }
  out
}

///|
/// Collect the atom indices that lie inside each axis-aligned box.
/// `boxes` holds (min_corner, max_corner) pairs.
pub fn BiotiteCellList::create_by_boxes(
  self : BiotiteCellList,
  boxes : Array[(Vector3, Vector3)],
) -> Array[Array[Int32]] {
  let out : Array[Array[Int32]] = Array(capacity=boxes.length())
  for (lo, hi) in boxes {
    let inside : Array[Int32] = []
    let min_x = if lo.x < hi.x { lo.x } else { hi.x }
    let max_x = if lo.x < hi.x { hi.x } else { lo.x }
    let min_y = if lo.y < hi.y { lo.y } else { hi.y }
    let max_y = if lo.y < hi.y { hi.y } else { lo.y }
    let min_z = if lo.z < hi.z { lo.z } else { hi.z }
    let max_z = if lo.z < hi.z { hi.z } else { lo.z }
    for i = 0; i < self.positions.length(); i = i + 1 {
      let p = self.positions[i]
      if p.x >= min_x &&
        p.x <= max_x &&
        p.y >= min_y &&
        p.y <= max_y &&
        p.z >= min_z &&
        p.z <= max_z {
        inside.push(i)
      }
    }
    out.push(inside)
  }
  out
}

///|
pub fn _biotite_structure_celllist_dummy() -> Unit {
  let cl = biotite_celllist_new([], cell_size=2.0)
  let _ = cl.adjacent_atoms(0, 1.0)
  let _ = cl.nearest_neighbors(0, number=1, distance=None)
  let _ = cl.conflict(1.0)
  let _ = cl.create_by_boxes([])
  ignore(cl.cell_size)
}
