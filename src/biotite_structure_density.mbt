/// Biotite `structure.density` — electron/scattering density maps.
///
/// Each atom contributes a Gaussian density with a width derived from
/// the resolution (FWHM) and a weight given by the electron count
/// (X-ray, default) or the coherent scattering length (neutron).
/// Unknown/ambiguous elements contribute nothing (Biotite: NaN).

///|
/// Electron count per element (atomic number for the common organic
/// elements); None for unknown or pseudo atoms.
pub fn biotite_density_electron_count(element : String) -> Int32? {
  match element {
    "H" | "D" => Some(1)
    "C" => Some(6)
    "N" => Some(7)
    "O" => Some(8)
    "P" => Some(15)
    "S" => Some(16)
    "F" => Some(9)
    "Cl" => Some(17)
    "Br" => Some(35)
    "I" => Some(53)
    "Na" => Some(11)
    "Mg" => Some(12)
    "K" => Some(19)
    "Ca" => Some(20)
    "Fe" => Some(26)
    "Zn" => Some(30)
    _ => None
  }
}

///|
/// Coherent neutron scattering length (fm) per element;
/// None for unknown or pseudo atoms.
pub fn biotite_density_scattering_length(element : String) -> Double? {
  match element {
    "H" => Some(-3.739)
    "D" => Some(6.671)
    "C" => Some(6.646)
    "N" => Some(9.36)
    "O" => Some(5.803)
    "P" => Some(5.13)
    "S" => Some(2.847)
    "F" => Some(5.654)
    "Cl" => Some(9.577)
    "Br" => Some(6.795)
    "I" => Some(6.191)
    "Na" => Some(3.63)
    "Mg" => Some(4.588)
    "K" => Some(3.74)
    "Ca" => Some(4.7)
    "Fe" => Some(9.45)
    "Zn" => Some(5.68)
    _ => None
  }
}

///|
/// A 3D density grid with the box origin and the cell size.
pub struct BiotiteDensity {
  // density[z][y][x]
  density : Array[Array[Array[Double]]]
  origin : Vector3
  cell_size : Double
  n_cells_x : Int32
  n_cells_y : Int32
  n_cells_z : Int32
}

///|
/// Compute a density map over the structure's bounding box.
///
/// - `resolution`: FWHM of the per-atom Gaussian in Angstrom.
/// - `scattering_length_consistency`: use neutron scattering lengths
///   instead of electron counts as weights.
/// - `box`: optional (min_corner, max_corner); defaults to the
///   bounding box of all atoms.
pub fn biotite_density_from_atom_array(
  arr : AtomArray,
  resolution~ : Double = 1.0,
  scattering_length_consistency~ : Bool = false,
  box~ : (Vector3, Vector3)? = None,
) -> BiotiteDensity {
  let n = arr.len()
  let mut min_x = 0.0
  let mut min_y = 0.0
  let mut min_z = 0.0
  let mut max_x = 0.0
  let mut max_y = 0.0
  let mut max_z = 0.0
  if n > 0 {
    let c0 = arr.coords[0]
    min_x = c0.x
    min_y = c0.y
    min_z = c0.z
    max_x = c0.x
    max_y = c0.y
    max_z = c0.z
    for i = 1; i < n; i = i + 1 {
      let c = arr.coords[i]
      if c.x < min_x {
        min_x = c.x
      }
      if c.y < min_y {
        min_y = c.y
      }
      if c.z < min_z {
        min_z = c.z
      }
      if c.x > max_x {
        max_x = c.x
      }
      if c.y > max_y {
        max_y = c.y
      }
      if c.z > max_z {
        max_z = c.z
      }
    }
  }
  let (box_min, box_max) = match box {
    Some((lo, hi)) => (lo, hi)
    None => {
      let lo = Vector3::{ x: min_x, y: min_y, z: min_z, }
      let hi = Vector3::{ x: max_x, y: max_y, z: max_z, }
      (lo, hi)
    }
  }
  let cell_size = resolution
  let pad = 2.0 * resolution
  let origin = Vector3::{
    x: box_min.x - pad,
    y: box_min.y - pad,
    z: box_min.z - pad,
  }
  let nx = ((box_max.x - box_min.x + 2.0 * pad) / cell_size).to_int() + 1
  let ny = ((box_max.y - box_min.y + 2.0 * pad) / cell_size).to_int() + 1
  let nz = ((box_max.z - box_min.z + 2.0 * pad) / cell_size).to_int() + 1
  let grid = biotite_density_grid(nz, ny, nx)
  // FWHM -> sigma: sigma = FWHM / (2 * sqrt(2 * ln 2))
  let sigma = resolution / 2.3548200450309493
  let cutoff = 3.0 * sigma
  for i = 0; i < n; i = i + 1 {
    let weight = if scattering_length_consistency {
      biotite_density_scattering_length(arr.elements[i])
    } else {
      match biotite_density_electron_count(arr.elements[i]) {
        Some(e) => Some(e.to_double())
        None => None
      }
    }
    match weight {
      Some(w) => {
        let c = arr.coords[i]
        biotite_density_splat(
          grid, origin, cell_size, nz, ny, nx, c, w, sigma, cutoff,
        )
      }
      None => ()
    }
  }
  BiotiteDensity::{
    density: grid,
    origin,
    cell_size,
    n_cells_x: nx,
    n_cells_y: ny,
    n_cells_z: nz,
  }
}

///|
/// Allocate an all-zero grid; rows are allocated independently to
/// avoid shared-array aliasing.
fn biotite_density_grid(
  nz : Int32,
  ny : Int32,
  nx : Int32,
) -> Array[Array[Array[Double]]] {
  let grid : Array[Array[Array[Double]]] = Array(capacity=nz)
  for z = 0; z < nz; z = z + 1 {
    let row : Array[Array[Double]] = Array(capacity=ny)
    for _y = 0; _y < ny; _y = _y + 1 {
      row.push(Array::make(nx, 0.0))
    }
    grid.push(row)
  }
  grid
}

///|
/// Add one atom's Gaussian contribution to the grid (3-sigma truncated).
fn biotite_density_splat(
  grid : Array[Array[Array[Double]]],
  origin : Vector3,
  cell_size : Double,
  nz : Int32,
  ny : Int32,
  nx : Int32,
  center : Vector3,
  weight : Double,
  sigma : Double,
  cutoff : Double,
) -> Unit {
  let two_sigma2 = 2.0 * sigma * sigma
  let r_lo = ((center.x - cutoff - origin.x) / cell_size).to_int()
  let r_hi = ((center.x + cutoff - origin.x) / cell_size).to_int()
  let y_lo = ((center.y - cutoff - origin.y) / cell_size).to_int()
  let y_hi = ((center.y + cutoff - origin.y) / cell_size).to_int()
  let z_lo = ((center.z - cutoff - origin.z) / cell_size).to_int()
  let z_hi = ((center.z + cutoff - origin.z) / cell_size).to_int()
  for iz = z_lo; iz <= z_hi; iz = iz + 1 {
    if iz < 0 || iz >= nz {
      continue
    }
    for iy = y_lo; iy <= y_hi; iy = iy + 1 {
      if iy < 0 || iy >= ny {
        continue
      }
      for ix = r_lo; ix <= r_hi; ix = ix + 1 {
        if ix < 0 || ix >= nx {
          continue
        }
        let px = origin.x + ix.to_double() * cell_size
        let py = origin.y + iy.to_double() * cell_size
        let pz = origin.z + iz.to_double() * cell_size
        let dx = px - center.x
        let dy = py - center.y
        let dz = pz - center.z
        let r2 = dx * dx + dy * dy + dz * dz
        if r2 <= cutoff * cutoff {
          grid[iz][iy][ix] = grid[iz][iy][ix] +
            weight * @math.exp(-r2 / two_sigma2)
        }
      }
    }
  }
}

///|
/// Grid extents: (nz, ny, nx).
pub fn BiotiteDensity::shape(self : BiotiteDensity) -> (Int32, Int32, Int32) {
  (self.n_cells_z, self.n_cells_y, self.n_cells_x)
}

///|
/// Total density over the whole grid.
pub fn BiotiteDensity::total(self : BiotiteDensity) -> Double {
  let mut sum = 0.0
  for plane in self.density {
    for row in plane {
      for v in row {
        sum = sum + v
      }
    }
  }
  sum
}

///|
/// Maximum density value in the grid.
pub fn BiotiteDensity::max(self : BiotiteDensity) -> Double {
  let mut best = -1.7976931348623157e308
  for plane in self.density {
    for row in plane {
      for v in row {
        if v > best {
          best = v
        }
      }
    }
  }
  best
}

///|
pub fn _biotite_structure_density_dummy() -> Unit {
  let d = biotite_density_from_atom_array(AtomArray::new())
  let _ = d.shape()
  let _ = d.total()
  let _ = d.max()
  let _ = biotite_density_electron_count("C")
  let _ = biotite_density_scattering_length("H")
  ignore(d.cell_size)
  ignore(d.n_cells_x)
  ignore(d.origin)
}
