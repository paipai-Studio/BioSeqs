/// Biotite `structure.celllist` — official `CellList` API.
///
/// Mirrors the official Rust implementation
/// (`src/rust/structure/celllist.rs` in biotite-dev/biotite, backing
/// `biotite.structure.celllist.CellList`; the pre-Rust Cython version
/// exposes the same three query methods):
/// - atoms are bucketed into a 3D grid of cubic cells with edge
///   `cell_size`; the grid origin is the per-axis minimum coordinate,
///   the cell of a coordinate is `floor((c - min) / cell_size)` and
///   the grid dimensions are the maximum coordinate's cell + 1,
/// - `get_atoms(coord, radius)` returns all atoms with
///   `distance <= radius` (boundary inclusive), scanning every cell
///   within `ceil(radius / cell_size)` cells; the result is
///   independent of `cell_size` (official parametrized test),
/// - `get_atoms_in_cells(coord, cell_radius=1)` returns the atoms in
///   the cell of the coordinate and its `cell_radius` neighbors with
///   no distance filter,
/// - results are collected in the official cell traversal order
///   (axis-major loops ascending per dimension),
/// - `create_adjacency_matrix(threshold)` has
///   `m[i,j] = dist(i,j) <= threshold` semantics (the official test
///   compares against a brute-force distance matrix); negative
///   thresholds yield an all-False matrix (official: ValueError),
/// - construction with a non-positive `cell_size` or without any
///   coordinate raises in the official implementation ("Cell size
///   must be greater than 0" / "Coordinates must not be empty");
///   this port has no exceptions, so a non-positive cell size falls
///   back to 1.0 and empty coordinates yield an empty grid whose
///   queries are always empty,
/// - `selection` stores only the masked atoms while the returned
///   indices still refer to the original array (official semantics).

///|
/// Cubic-cell grid over atom positions with the official query
/// methods.
pub struct BiotiteCellList {
  /// Edge length of one cell (a valid positive value).
  cell_size : Double
  min_corner : Vector3
  dims : (Int32, Int32, Int32)
  cells : Array[Array[Int32]]
  positions : Array[Vector3]
  selected : Array[Bool]?
}

///|
/// Cell index of a coordinate along one axis: truncation of
/// `(c - min) / cell_size` (equals floor for coordinates at or above
/// the grid minimum; negative results are handled as out of bounds).
fn bcl_cell_index(
  coord : Double,
  min_coord : Double,
  cell_size : Double,
) -> Int32 {
  ((coord - min_coord) / cell_size).to_int()
}

///|
/// Ceil of a ratio as an integer (official `ceil(radius / cell_size)`).
fn bcl_ceil_ratio(ratio : Double) -> Int32 {
  let fl = ratio.to_int()
  if ratio > fl.to_double() {
    fl + 1
  } else {
    fl
  }
}

///|
/// Build a Biotite cell list from atom positions.
pub fn biotite_celllist_new(
  positions : Array[Vector3],
  cell_size~ : Double = 2.0,
  selection~ : Array[Bool]? = None,
) -> BiotiteCellList {
  let cs = if cell_size > 0.0 { cell_size } else { 1.0 }
  let n = positions.length()
  // Per-axis coordinate range (official `calculate_coord_range`).
  let mut min_x = 0.0
  let mut min_y = 0.0
  let mut min_z = 0.0
  let mut max_x = 0.0
  let mut max_y = 0.0
  let mut max_z = 0.0
  let mut first = true
  for p in positions {
    if first {
      min_x = p.x
      min_y = p.y
      min_z = p.z
      max_x = p.x
      max_y = p.y
      max_z = p.z
      first = false
    } else {
      if p.x < min_x {
        min_x = p.x
      }
      if p.y < min_y {
        min_y = p.y
      }
      if p.z < min_z {
        min_z = p.z
      }
      if p.x > max_x {
        max_x = p.x
      }
      if p.y > max_y {
        max_y = p.y
      }
      if p.z > max_z {
        max_z = p.z
      }
    }
  }
  // Grid dimensions: maximum coordinate's cell + 1 per axis.
  let di = bcl_cell_index(max_x, min_x, cs) + 1
  let dj = bcl_cell_index(max_y, min_y, cs) + 1
  let dk = bcl_cell_index(max_z, min_z, cs) + 1
  let cells : Array[Array[Int32]] = Array(capacity=di * dj * dk)
  let n_cells = di * dj * dk
  for _i = 0; _i < n_cells; _i = _i + 1 {
    cells.push([])
  }
  for i = 0; i < n; i = i + 1 {
    match selection {
      Some(sel) => if i < sel.length() && !sel[i] { continue }
      None => ()
    }
    let p = positions[i]
    let ci = bcl_cell_index(p.x, min_x, cs)
    let cj = bcl_cell_index(p.y, min_y, cs)
    let ck = bcl_cell_index(p.z, min_z, cs)
    if ci >= 0 && ci < di && cj >= 0 && cj < dj && ck >= 0 && ck < dk {
      cells[(ci * dj + cj) * dk + ck].push(i)
    }
  }
  BiotiteCellList::{
    cell_size: cs,
    min_corner: Vector3::new(min_x, min_y, min_z),
    dims: (di, dj, dk),
    cells,
    positions,
    selected: selection,
  }
}

///|
/// Cell position of a query coordinate.
fn BiotiteCellList::cell_position(
  self : BiotiteCellList,
  c : Vector3,
) -> (Int32, Int32, Int32) {
  (
    bcl_cell_index(c.x, self.min_corner.x, self.cell_size),
    bcl_cell_index(c.y, self.min_corner.y, self.cell_size),
    bcl_cell_index(c.z, self.min_corner.z, self.cell_size),
  )
}

///|
/// Whether atom `i` is stored (selection mask; all atoms when no
/// selection was given).
fn BiotiteCellList::is_selected(self : BiotiteCellList, i : Int32) -> Bool {
  match self.selected {
    Some(sel) => i < sel.length() && sel[i]
    None => true
  }
}

///|
/// Collect the atoms of all cells within `cell_radius` cells of the
/// coordinate (official cell traversal order: axis-major ascending).
fn BiotiteCellList::cells_around(
  self : BiotiteCellList,
  coord : Vector3,
  cell_radius : Int32,
) -> Array[Int32] {
  let out : Array[Int32] = []
  let (di, dj, dk) = self.dims
  if di <= 0 || dj <= 0 || dk <= 0 {
    return out
  }
  let (ci, cj, ck) = self.cell_position(coord)
  for a = ci - cell_radius; a <= ci + cell_radius; a = a + 1 {
    if a < 0 || a >= di {
      continue
    }
    for b = cj - cell_radius; b <= cj + cell_radius; b = b + 1 {
      if b < 0 || b >= dj {
        continue
      }
      for c = ck - cell_radius; c <= ck + cell_radius; c = c + 1 {
        if c < 0 || c >= dk {
          continue
        }
        for idx in self.cells[(a * dj + b) * dk + c] {
          out.push(idx)
        }
      }
    }
  }
  out
}

///|
/// All atoms within `radius` (Euclidean, boundary inclusive) of the
/// coordinate (official single-coordinate `get_atoms` with the
/// default MAPPING result).
pub fn BiotiteCellList::get_atoms(
  self : BiotiteCellList,
  coord : Vector3,
  radius : Double,
) -> Array[Int32] {
  let out : Array[Int32] = []
  if self.positions.length() == 0 {
    return out
  }
  let sq_r = radius * radius
  for idx in self.cells_around(coord, bcl_ceil_ratio(radius / self.cell_size)) {
    let p = self.positions[idx]
    let dx = p.x - coord.x
    let dy = p.y - coord.y
    let dz = p.z - coord.z
    if dx * dx + dy * dy + dz * dz <= sq_r {
      out.push(idx)
    }
  }
  out
}

///|
/// The atoms in the cell of the coordinate and its `cell_radius`
/// neighbors, with no distance filter (official `get_atoms_in_cells`,
/// default `cell_radius=1`).
pub fn BiotiteCellList::get_atoms_in_cells(
  self : BiotiteCellList,
  coord : Vector3,
  cell_radius~ : Int32 = 1,
) -> Array[Int32] {
  self.cells_around(coord, cell_radius)
}

///|
/// Neighbor pairs of multiple query coordinates as
/// `(query_idx, atom_idx)` tuples (official PAIRS result). Either a
/// single `radius` for all coordinates or `radii` per coordinate is
/// given, mirroring the official scalar/ndarray radius input.
pub fn BiotiteCellList::get_atoms_pairs(
  self : BiotiteCellList,
  coords : Array[Vector3],
  radius~ : Double? = None,
  radii~ : Array[Double]? = None,
) -> Array[(Int32, Int32)] {
  let out : Array[(Int32, Int32)] = []
  for q = 0; q < coords.length(); q = q + 1 {
    let r = match radii {
      Some(rs) => if q < rs.length() { rs[q] } else { 0.0 }
      None =>
        match radius {
          Some(rv) => rv
          None => 0.0
        }
    }
    for idx in self.get_atoms(coords[q], r) {
      out.push((q, idx))
    }
  }
  out
}

///|
/// Boolean mask of shape (queries, atoms) where
/// `mask[q][j] = dist(coords[q], atom j) <= radius` (official MASK
/// result).
pub fn BiotiteCellList::get_atoms_mask(
  self : BiotiteCellList,
  coords : Array[Vector3],
  radius : Double,
) -> Array[Array[Bool]] {
  let n = self.positions.length()
  let out : Array[Array[Bool]] = Array(capacity=coords.length())
  for q = 0; q < coords.length(); q = q + 1 {
    let row : Array[Bool] = Array(capacity=n)
    for _j = 0; _j < n; _j = _j + 1 {
      row.push(false)
    }
    for idx in self.get_atoms(coords[q], radius) {
      row[idx] = true
    }
    out.push(row)
  }
  out
}

///|
/// Adjacency matrix with `m[i,j] = dist(i,j) <= threshold_distance`
/// (official `create_adjacency_matrix`; the official test compares it
/// against a brute-force distance matrix, so the brute-force
/// computation here is semantically identical). A negative threshold
/// yields an all-False matrix (official: ValueError). Atoms not
/// masked by a constructor `selection` have all-False rows and
/// columns.
pub fn BiotiteCellList::create_adjacency_matrix(
  self : BiotiteCellList,
  threshold_distance : Double,
) -> Array[Array[Bool]] {
  let n = self.positions.length()
  let matrix : Array[Array[Bool]] = Array(capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    let row : Array[Bool] = Array(capacity=n)
    for _j = 0; _j < n; _j = _j + 1 {
      row.push(false)
    }
    matrix.push(row)
  }
  if threshold_distance < 0.0 {
    return matrix
  }
  let sq_r = threshold_distance * threshold_distance
  for i = 0; i < n; i = i + 1 {
    if !self.is_selected(i) {
      continue
    }
    let pi = self.positions[i]
    for j = 0; j < n; j = j + 1 {
      if !self.is_selected(j) {
        continue
      }
      let pj = self.positions[j]
      let dx = pi.x - pj.x
      let dy = pi.y - pj.y
      let dz = pi.z - pj.z
      if dx * dx + dy * dy + dz * dz <= sq_r {
        matrix[i][j] = true
      }
    }
  }
  matrix
}

///|
pub fn _biotite_structure_celllist_dummy() -> Unit {
  let cl = biotite_celllist_new([], cell_size=2.0)
  let _ = cl.get_atoms(Vector3::new(0.0, 0.0, 0.0), 1.0)
  let _ = cl.get_atoms_in_cells(Vector3::new(0.0, 0.0, 0.0))
  let _ = cl.get_atoms_pairs([], radius=Some(1.0))
  let _ = cl.get_atoms_mask([], 1.0)
  let _ = cl.create_adjacency_matrix(1.0)
  ignore(cl.cell_size)
}
