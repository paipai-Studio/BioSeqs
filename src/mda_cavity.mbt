/// Cavity detection and volume estimation.
///
/// MDAnalysis equivalent: `analysis.cavity.CavityAnalysis`
///
/// Identifies internal cavities in protein structures by grid-based
/// void space detection, complementing the existing `cavifinder.mbt`.

///|
/// A detected cavity
pub struct MdaCavity {
  id : Int32
  volume : Double
  center : Array[Double]
  n_grid_points : Int32
  grid_indices : Array[(Int32, Int32, Int32)]
}

///|
pub fn MdaCavity::id(self : MdaCavity) -> Int32 {
  self.id
}

///|
pub fn MdaCavity::volume(self : MdaCavity) -> Double {
  self.volume
}

///|
pub fn MdaCavity::center(self : MdaCavity) -> Array[Double] {
  self.center
}

///|
pub fn MdaCavity::n_grid_points(self : MdaCavity) -> Int32 {
  self.n_grid_points
}

///|
/// Detect cavities using a grid-based approach.
///
/// atom_coords: [atom][xyz]
/// atom_radii: van der Waals radii [atom]
/// grid_spacing: spacing between grid points (Å)
/// probe_radius: radius of probe sphere (Å)
pub fn mda_detect_cavities(
  atom_coords : Array[Array[Double]],
  atom_radii : Array[Double],
  grid_spacing~ : Double = 1.0,
  probe_radius~ : Double = 1.4,
) -> Array[MdaCavity] {
  let n_atoms = atom_coords.length()
  if n_atoms == 0 {
    return []
  }
  // Compute bounding box
  let mut min_x = 10000000000.0
  let mut min_y = 10000000000.0
  let mut min_z = 10000000000.0
  let mut max_x = -10000000000.0
  let mut max_y = -10000000000.0
  let mut max_z = -10000000000.0
  for i = 0; i < n_atoms; i = i + 1 {
    let r = atom_radii[i] + probe_radius
    if atom_coords[i][0] - r < min_x {
      min_x = atom_coords[i][0] - r
    }
    if atom_coords[i][1] - r < min_y {
      min_y = atom_coords[i][1] - r
    }
    if atom_coords[i][2] - r < min_z {
      min_z = atom_coords[i][2] - r
    }
    if atom_coords[i][0] + r > max_x {
      max_x = atom_coords[i][0] + r
    }
    if atom_coords[i][1] + r > max_y {
      max_y = atom_coords[i][1] + r
    }
    if atom_coords[i][2] + r > max_z {
      max_z = atom_coords[i][2] + r
    }
  }
  // Grid dimensions
  let nx = ((max_x - min_x) / grid_spacing).to_int() + 1
  let ny = ((max_y - min_y) / grid_spacing).to_int() + 1
  let nz = ((max_z - min_z) / grid_spacing).to_int() + 1
  // Mark grid points occupied by atoms
  // occupied[ix][iy][iz] = true if inside any atom sphere
  let total = nx * ny * nz
  if total > 1000000 {
    // Too large, skip
    return []
  }
  let occupied : Array[Bool] = Array::make(total, false)
  for ix = 0; ix < nx; ix = ix + 1 {
    let x = min_x + ix.to_double() * grid_spacing
    for iy = 0; iy < ny; iy = iy + 1 {
      let y = min_y + iy.to_double() * grid_spacing
      for iz = 0; iz < nz; iz = iz + 1 {
        let z = min_z + iz.to_double() * grid_spacing
        for a = 0; a < n_atoms; a = a + 1 {
          let dx = x - atom_coords[a][0]
          let dy = y - atom_coords[a][1]
          let dz = z - atom_coords[a][2]
          let r = atom_radii[a] + probe_radius
          if dx * dx + dy * dy + dz * dz < r * r {
            occupied[ix * ny * nz + iy * nz + iz] = true
            break
          }
        }
      }
    }
  }
  // Flood fill from boundary to mark exterior (free) space
  let exterior : Array[Bool] = Array::make(total, false)
  let queue : Array[Int32] = Array(capacity=0)
  // Seed from all boundary faces
  for ix = 0; ix < nx; ix = ix + 1 {
    for iy = 0; iy < ny; iy = iy + 1 {
      let idx1 = ix * ny * nz + iy * nz + 0
      let idx2 = ix * ny * nz + iy * nz + (nz - 1)
      if !occupied[idx1] && !exterior[idx1] {
        exterior[idx1] = true
        queue.push(idx1)
      }
      if !occupied[idx2] && !exterior[idx2] {
        exterior[idx2] = true
        queue.push(idx2)
      }
    }
  }
  for ix = 0; ix < nx; ix = ix + 1 {
    for iz = 0; iz < nz; iz = iz + 1 {
      let idx1 = ix * ny * nz + 0 * nz + iz
      let idx2 = ix * ny * nz + (ny - 1) * nz + iz
      if !occupied[idx1] && !exterior[idx1] {
        exterior[idx1] = true
        queue.push(idx1)
      }
      if !occupied[idx2] && !exterior[idx2] {
        exterior[idx2] = true
        queue.push(idx2)
      }
    }
  }
  for iy = 0; iy < ny; iy = iy + 1 {
    for iz = 0; iz < nz; iz = iz + 1 {
      let idx1 = 0 * ny * nz + iy * nz + iz
      let idx2 = (nx - 1) * ny * nz + iy * nz + iz
      if !occupied[idx1] && !exterior[idx1] {
        exterior[idx1] = true
        queue.push(idx1)
      }
      if !occupied[idx2] && !exterior[idx2] {
        exterior[idx2] = true
        queue.push(idx2)
      }
    }
  }
  // BFS flood fill
  let mut head = 0
  while head < queue.length() {
    let idx = queue[head]
    head = head + 1
    let ix = idx / (ny * nz)
    let rem = idx % (ny * nz)
    let iy = rem / nz
    let iz = rem % nz
    // Check 6 neighbors
    let neighbors = [
      (ix - 1, iy, iz),
      (ix + 1, iy, iz),
      (ix, iy - 1, iz),
      (ix, iy + 1, iz),
      (ix, iy, iz - 1),
      (ix, iy, iz + 1),
    ]
    for n in neighbors {
      if n.0 >= 0 && n.0 < nx && n.1 >= 0 && n.1 < ny && n.2 >= 0 && n.2 < nz {
        let nidx = n.0 * ny * nz + n.1 * nz + n.2
        if !occupied[nidx] && !exterior[nidx] {
          exterior[nidx] = true
          queue.push(nidx)
        }
      }
    }
  }
  // Remaining unoccupied, non-exterior points are cavities
  // Label connected cavity components
  let visited : Array[Bool] = Array::make(total, false)
  let cavities : Array[MdaCavity] = Array(capacity=0)
  let mut cavity_id = 0
  for idx = 0; idx < total; idx = idx + 1 {
    if occupied[idx] || exterior[idx] || visited[idx] {
      continue
    }
    // BFS to find cavity component
    let component : Array[Int32] = Array(capacity=0)
    let q2 : Array[Int32] = Array(capacity=0)
    q2.push(idx)
    visited[idx] = true
    let mut qhead = 0
    while qhead < q2.length() {
      let cur = q2[qhead]
      qhead = qhead + 1
      component.push(cur)
      let cix = cur / (ny * nz)
      let crem = cur % (ny * nz)
      let ciy = crem / nz
      let ciz = crem % nz
      let neighbors = [
        (cix - 1, ciy, ciz),
        (cix + 1, ciy, ciz),
        (cix, ciy - 1, ciz),
        (cix, ciy + 1, ciz),
        (cix, ciy, ciz - 1),
        (cix, ciy, ciz + 1),
      ]
      for n in neighbors {
        if n.0 >= 0 && n.0 < nx && n.1 >= 0 && n.1 < ny && n.2 >= 0 && n.2 < nz {
          let nidx = n.0 * ny * nz + n.1 * nz + n.2
          if !occupied[nidx] && !exterior[nidx] && !visited[nidx] {
            visited[nidx] = true
            q2.push(nidx)
          }
        }
      }
    }
    // Compute cavity properties
    let mut cx = 0.0
    let mut cy = 0.0
    let mut cz = 0.0
    let grid_indices : Array[(Int32, Int32, Int32)] = Array(capacity=0)
    for c in component {
      let cix = c / (ny * nz)
      let crem = c % (ny * nz)
      let ciy = crem / nz
      let ciz = crem % nz
      cx = cx + min_x + cix.to_double() * grid_spacing
      cy = cy + min_y + ciy.to_double() * grid_spacing
      cz = cz + min_z + ciz.to_double() * grid_spacing
      grid_indices.push((cix, ciy, ciz))
    }
    let nc = component.length()
    cx = cx / nc.to_double()
    cy = cy / nc.to_double()
    cz = cz / nc.to_double()
    let vol = nc.to_double() * grid_spacing * grid_spacing * grid_spacing
    cavities.push({
      id: cavity_id,
      volume: vol,
      center: [cx, cy, cz],
      n_grid_points: nc,
      grid_indices,
    })
    cavity_id = cavity_id + 1
  }
  cavities
}

///|
/// Total cavity volume
pub fn mda_total_cavity_volume(cavities : Array[MdaCavity]) -> Double {
  let mut total = 0.0
  for c in cavities {
    total = total + c.volume
  }
  total
}

///|
/// Largest cavity
pub fn mda_largest_cavity(cavities : Array[MdaCavity]) -> MdaCavity? {
  if cavities.length() == 0 {
    return None
  }
  let mut largest = cavities[0]
  for c in cavities {
    if c.volume > largest.volume {
      largest = c
    }
  }
  Some(largest)
}
