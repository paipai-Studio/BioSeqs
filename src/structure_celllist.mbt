///|
/// CellList — grid-based neighbor search for O(N) spatial queries.
/// Ported from biotite.structure.celllist.
///
/// Divides space into uniform cells. For each atom, finds all neighbors
/// within a cutoff radius by checking only the 3×3×3 surrounding cells.

///|
/// A cell list for fast spatial neighbor queries.
pub struct CellList {
  cell_size : Double
  origin : Vector3
  n_cells_x : Int
  n_cells_y : Int
  n_cells_z : Int
  // Flat grid: cell (i,j,k) → list of atom indices
  grid : Map[Int, Array[Int]]
}

///|
/// Build a CellList from atom coordinates.
///
/// Parameters:
///   `coords` — atom positions
///   `cell_size` — cell edge length (Å), should be ≥ cutoff radius
pub fn CellList::new(
  coords : Array[Vector3],
  cell_size~ : Double,
) -> CellList {
  if coords.length() == 0 {
    return CellList::{
      cell_size,
      origin: Vector3::new(0.0, 0.0, 0.0),
      n_cells_x: 0,
      n_cells_y: 0,
      n_cells_z: 0,
      grid: Map([], capacity=0),
    }
  }
  // Compute bounding box
  let mut min_x = coords[0].x
  let mut max_x = coords[0].x
  let mut min_y = coords[0].y
  let mut max_y = coords[0].y
  let mut min_z = coords[0].z
  let mut max_z = coords[0].z
  for c in coords {
    if c.x < min_x {
      min_x = c.x
    }
    if c.x > max_x {
      max_x = c.x
    }
    if c.y < min_y {
      min_y = c.y
    }
    if c.y > max_y {
      max_y = c.y
    }
    if c.z < min_z {
      min_z = c.z
    }
    if c.z > max_z {
      max_z = c.z
    }
  }
  let origin = Vector3::new(min_x, min_y, min_z)
  let n_cells_x = ((max_x - min_x) / cell_size).to_int() + 1
  let n_cells_y = ((max_y - min_y) / cell_size).to_int() + 1
  let n_cells_z = ((max_z - min_z) / cell_size).to_int() + 1
  // Assign atoms to cells
  let grid : Map[Int, Array[Int]] = Map([], capacity=coords.length())
  for i = 0; i < coords.length(); i = i + 1 {
    let cx = ((coords[i].x - min_x) / cell_size).to_int()
    let cy = ((coords[i].y - min_y) / cell_size).to_int()
    let cz = ((coords[i].z - min_z) / cell_size).to_int()
    let cell_id = cell_index(cx, cy, cz, n_cells_y, n_cells_z)
    match grid.get(cell_id) {
      Some(arr) => arr.push(i)
      None => {
        let arr : Array[Int] = Array(capacity=1)
        arr.push(i)
        grid[cell_id] = arr
      }
    }
  }
  CellList::{
    cell_size,
    origin,
    n_cells_x,
    n_cells_y,
    n_cells_z,
    grid,
  }
}

///|
/// Compute flat cell index from 3D cell coordinates.
fn cell_index(cx : Int, cy : Int, cz : Int, ny : Int, nz : Int) -> Int {
  cx * ny * nz + cy * nz + cz
}

///|
/// Find all atoms within `radius` of the atom at index `center_idx`.
/// Returns a list of (atom_index, distance) pairs, excluding the center itself.
pub fn CellList::get_neighbors(
  self : CellList,
  coords : Array[Vector3],
  center_idx : Int,
  radius~ : Double,
) -> Array[(Int, Double)] {
  let result : Array[(Int, Double)] = Array(capacity=0)
  if center_idx < 0 || center_idx >= coords.length() {
    return result
  }
  let center = coords[center_idx]
  let cx = ((center.x - self.origin.x) / self.cell_size).to_int()
  let cy = ((center.y - self.origin.y) / self.cell_size).to_int()
  let cz = ((center.z - self.origin.z) / self.cell_size).to_int()
  // Search 3×3×3 surrounding cells
  for dx = -1; dx <= 1; dx = dx + 1 {
    for dy = -1; dy <= 1; dy = dy + 1 {
      for dz = -1; dz <= 1; dz = dz + 1 {
        let ncx = cx + dx
        let ncy = cy + dy
        let ncz = cz + dz
        if ncx < 0 || ncx >= self.n_cells_x ||
           ncy < 0 || ncy >= self.n_cells_y ||
           ncz < 0 || ncz >= self.n_cells_z {
          continue
        }
        let cell_id = cell_index(ncx, ncy, ncz, self.n_cells_y, self.n_cells_z)
        match self.grid.get(cell_id) {
          Some(atom_indices) => {
            for idx in atom_indices {
              if idx == center_idx {
                continue
              }
              let dist = center.distance(coords[idx])
              if dist <= radius {
                result.push((idx, dist))
              }
            }
          }
          None => ()
        }
      }
    }
  }
  result
}

///|
/// Find all atom pairs within `radius` of each other.
/// Returns list of (index_a, index_b, distance) with a < b.
pub fn CellList::get_all_pairs(
  self : CellList,
  coords : Array[Vector3],
  radius~ : Double,
) -> Array[(Int, Int, Double)] {
  let result : Array[(Int, Int, Double)] = Array(capacity=0)
  let n = coords.length()
  for i = 0; i < n; i = i + 1 {
    let neighbors = self.get_neighbors(coords, i, radius~)
    for pair in neighbors {
      let j = pair.0
      let dist = pair.1
      if i < j {
        result.push((i, j, dist))
      }
    }
  }
  result
}

///|
/// Count atoms within `radius` of each atom (coordination number).
pub fn CellList::coordination_numbers(
  self : CellList,
  coords : Array[Vector3],
  radius~ : Double,
) -> Array[Int] {
  let n = coords.length()
  let result : Array[Int] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let neighbors = self.get_neighbors(coords, i, radius~)
    result.push(neighbors.length())
  }
  result
}
