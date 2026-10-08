///|
/// CellList — grid-based spatial neighbor search.
/// Ported from biotite.structure.celllist.
///
/// Divides space into uniform cells. For each atom, neighbor atoms
/// are found by searching only the surrounding cells,
/// giving O(N) average time instead of O(N²) for all-pairs.

///|
/// A CellList partitions atoms into uniform grid cells for fast
/// neighbor queries.
pub struct CellList {
  cell_size : Double
  // Cell index → list of atom indices
  cells : Map[(Int, Int, Int), Array[Int]]
  // Atom index → coordinate
  coords : Array[Array[Double]]
  // Bounding box
  min_coord : Array[Double]
  max_coord : Array[Double]
}

///|
/// Build a CellList from atom coordinates.
/// `cell_size` should be at least as large as the cutoff radius
/// for neighbor queries.
pub fn celllist_new(
  coords : Array[Array[Double]],
  cell_size : Double,
) -> CellList {
  let n = coords.length()
  if n == 0 {
    return CellList::{
      cell_size,
      cells: Map([], capacity=0),
      coords,
      min_coord: [0.0, 0.0, 0.0],
      max_coord: [0.0, 0.0, 0.0],
    }
  }
  // Compute bounding box
  let mut min_x = coords[0][0]
  let mut min_y = coords[0][1]
  let mut min_z = coords[0][2]
  let mut max_x = coords[0][0]
  let mut max_y = coords[0][1]
  let mut max_z = coords[0][2]
  for i = 1; i < n; i = i + 1 {
    let x = coords[i][0]
    let y = coords[i][1]
    let z = coords[i][2]
    if x < min_x {
      min_x = x
    }
    if y < min_y {
      min_y = y
    }
    if z < min_z {
      min_z = z
    }
    if x > max_x {
      max_x = x
    }
    if y > max_y {
      max_y = y
    }
    if z > max_z {
      max_z = z
    }
  }
  // Assign atoms to cells
  let cells : Map[(Int, Int, Int), Array[Int]] = Map([], capacity=n)
  for i = 0; i < n; i = i + 1 {
    let cx = cell_index(coords[i][0], min_x, cell_size)
    let cy = cell_index(coords[i][1], min_y, cell_size)
    let cz = cell_index(coords[i][2], min_z, cell_size)
    let key = (cx, cy, cz)
    match cells.get(key) {
      Some(arr) => arr.push(i)
      None => {
        let arr : Array[Int] = [i]
        cells[key] = arr
      }
    }
  }
  CellList::{
    cell_size,
    cells,
    coords,
    min_coord: [min_x, min_y, min_z],
    max_coord: [max_x, max_y, max_z],
  }
}

///|
/// Compute cell index for a coordinate.
fn cell_index(coord : Double, min_coord : Double, cell_size : Double) -> Int {
  ((coord - min_coord) / cell_size).to_int()
}

///|
/// Find all atoms within `radius` of atom `center_idx`.
/// Returns indices of neighbor atoms (excluding self).
pub fn celllist_neighbors(
  cl : CellList,
  center_idx : Int,
  radius : Double,
) -> Array[Int] {
  if center_idx < 0 || center_idx >= cl.coords.length() {
    return Array(capacity=0)
  }
  let cx = cl.coords[center_idx][0]
  let cy = cl.coords[center_idx][1]
  let cz = cl.coords[center_idx][2]
  celllist_neighbors_point(cl, cx, cy, cz, radius, center_idx)
}

///|
/// Find all atoms within `radius` of point (px, py, pz).
/// `exclude_idx` is excluded from results (use -1 to exclude none).
pub fn celllist_neighbors_point(
  cl : CellList,
  px : Double,
  py : Double,
  pz : Double,
  radius : Double,
  exclude_idx : Int,
) -> Array[Int] {
  let result : Array[Int] = Array(capacity=16)
  // Number of cells to search in each direction: ceil(radius / cell_size)
  let ratio = radius / cl.cell_size
  let cell_r = if ratio <= 1.0 {
    1
  } else {
    ratio.to_int() + 1
  }
  let cxi = cell_index(px, cl.min_coord[0], cl.cell_size)
  let cyi = cell_index(py, cl.min_coord[1], cl.cell_size)
  let czi = cell_index(pz, cl.min_coord[2], cl.cell_size)
  let r2 = radius * radius
  for dx = -cell_r; dx <= cell_r; dx = dx + 1 {
    for dy = -cell_r; dy <= cell_r; dy = dy + 1 {
      for dz = -cell_r; dz <= cell_r; dz = dz + 1 {
        let key = (cxi + dx, cyi + dy, czi + dz)
        match cl.cells.get(key) {
          Some(indices) => {
            for j in indices {
              if j == exclude_idx {
                continue
              }
              let dx2 = cl.coords[j][0] - px
              let dy2 = cl.coords[j][1] - py
              let dz2 = cl.coords[j][2] - pz
              let d2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2
              if d2 <= r2 {
                result.push(j)
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
/// Find all pairs of atoms within `radius` of each other.
/// Returns array of (i, j) pairs with i < j.
pub fn celllist_all_pairs(
  cl : CellList,
  radius : Double,
) -> Array[(Int, Int)] {
  let result : Array[(Int, Int)] = Array(capacity=0)
  let n = cl.coords.length()
  for i = 0; i < n; i = i + 1 {
    let neighbors = celllist_neighbors(cl, i, radius)
    for j in neighbors {
      if i < j {
        result.push((i, j))
      }
    }
  }
  result
}
