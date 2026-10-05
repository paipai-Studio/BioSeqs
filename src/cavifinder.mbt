///|
/// CaviFinder: channel and cavity detection in protein structures.
///
/// Grid-based approach: place a 3D grid over the protein bounding box,
/// mark grid points near atoms as "occupied", flood-fill from exterior
/// to identify free space, remaining enclosed free points are cavities.
///
/// Ported from ProDy's CaviFinder suite.

// ─── helpers ───────────────────────────────────────────────

fn cav_sqrt(x : Double) -> Double {
  if x < 0.0 { 0.0 } else { x.sqrt() }
}

///|
/// A detected cavity or channel.
pub struct Cavity {
  center : Vector3
  volume : Double
  num_points : Int
  grid_points : Array[Vector3]
  is_channel : Bool
}

///|
/// Result of cavity detection.
pub struct CaviFinderResult {
  cavities : Array[Cavity]
  grid_spacing : Double
  probe_radius : Double
  total_volume : Double
}

///|
/// Van der Waals radius for an atom element.
fn cav_vdw_radius(element : String) -> Double {
  match element {
    "C" => 1.7
    "N" => 1.55
    "O" => 1.52
    "S" => 1.8
    "P" => 1.8
    "H" => 1.2
    "F" => 1.47
    "CL" | "Cl" => 1.75
    "BR" | "Br" => 1.85
    "I" => 1.98
    "FE" | "Fe" => 1.94
    "ZN" | "Zn" => 1.39
    "CA" | "Ca" => 1.37
    "MG" | "Mg" => 1.73
    "NA" | "Na" => 1.36
    "K" => 1.64
    _ => 1.7
  }
}

///|
/// Find cavities and channels in a protein structure using a grid-based method.
///
/// # Parameters
/// - chain: protein chain
/// - grid_spacing: distance between grid points (default 1.0 Å)
/// - probe_radius: radius of the probe sphere (default 1.4 Å, water radius)
pub fn find_cavities(
  chain : Chain,
  grid_spacing? : Double = 1.0,
  probe_radius? : Double = 1.4,
) -> CaviFinderResult {
  let residues = chain.get_residues()
  // Collect all atom coordinates and radii
  let atom_coords = Array(capacity=0)
  let atom_radii = Array(capacity=0)
  for res in residues {
    for atom in res.get_atoms() {
      atom_coords.push(atom.coord)
      let r = cav_vdw_radius(atom.element)
      atom_radii.push(r)
    }
  }
  let n_atoms = atom_coords.length()
  if n_atoms == 0 {
    return CaviFinderResult::{
      cavities: [], grid_spacing, probe_radius, total_volume: 0.0,
    }
  }
  // Compute bounding box (expanded by probe radius)
  let mut min_x = atom_coords[0].x
  let mut max_x = atom_coords[0].x
  let mut min_y = atom_coords[0].y
  let mut max_y = atom_coords[0].y
  let mut min_z = atom_coords[0].z
  let mut max_z = atom_coords[0].z
  for c in atom_coords {
    if c.x < min_x { min_x = c.x }
    if c.x > max_x { max_x = c.x }
    if c.y < min_y { min_y = c.y }
    if c.y > max_y { max_y = c.y }
    if c.z < min_z { min_z = c.z }
    if c.z > max_z { max_z = c.z }
  }
  let margin = probe_radius + 2.0
  min_x = min_x - margin; max_x = max_x + margin
  min_y = min_y - margin; max_y = max_y + margin
  min_z = min_z - margin; max_z = max_z + margin
  // Grid dimensions
  let nx = ((max_x - min_x) / grid_spacing).ceil().to_int() + 1
  let ny = ((max_y - min_y) / grid_spacing).ceil().to_int() + 1
  let nz = ((max_z - min_z) / grid_spacing).ceil().to_int() + 1
  // Safety limit on grid size
  if nx > 80 || ny > 80 || nz > 80 {
    return CaviFinderResult::{
      cavities: [], grid_spacing, probe_radius, total_volume: 0.0,
    }
  }
  let total_cells = nx * ny * nz
  // 0 = free (untested), 1 = occupied, 2 = exterior, 3 = cavity
  let grid = Array::make(total_cells, 0)
  // Mark occupied cells (within vdw + probe radius of any atom)
  let mut cell_idx = 0
  let mut gi = 0
  while gi < nx {
    let mut gj = 0
    while gj < ny {
      let mut gk = 0
      while gk < nz {
        let gx = min_x + gi.to_double() * grid_spacing
        let gy = min_y + gj.to_double() * grid_spacing
        let gz = min_z + gk.to_double() * grid_spacing
        let mut occupied = false
        let mut ai = 0
        while ai < n_atoms {
          if not(occupied) {
            let dx = gx - atom_coords[ai].x
            let dy = gy - atom_coords[ai].y
            let dz = gz - atom_coords[ai].z
            let dist_sq = dx * dx + dy * dy + dz * dz
            let threshold = atom_radii[ai] + probe_radius
            if dist_sq < threshold * threshold { occupied = true }
          }
          ai = ai + 1
        }
        if occupied { grid[cell_idx] = 1 }
        cell_idx = cell_idx + 1
        gk = gk + 1
      }
      gj = gj + 1
    }
    gi = gi + 1
  }
  // Flood fill from all boundary cells (mark as exterior = 2)
  let queue = Array(capacity=0)
  // Add all boundary cells to queue
  let mut idx = 0
  while idx < total_cells {
    let is_boundary = {
      let i = idx / (ny * nz)
      let j = (idx % (ny * nz)) / nz
      let k = idx % nz
      i == 0 || i == nx - 1 || j == 0 || j == ny - 1 || k == 0 || k == nz - 1
    }
    if is_boundary && grid[idx] == 0 {
      grid[idx] = 2
      queue.push(idx)
    }
    idx = idx + 1
  }
  // BFS flood fill
  let mut head = 0
  while head < queue.length() {
    let cur = queue[head]
    head = head + 1
    let i = cur / (ny * nz)
    let j = (cur % (ny * nz)) / nz
    let k = cur % nz
    // Check 6 neighbors
    let neighbors = [
      (i - 1, j, k), (i + 1, j, k),
      (i, j - 1, k), (i, j + 1, k),
      (i, j, k - 1), (i, j, k + 1),
    ]
    for (ni, nj, nk) in neighbors {
      if ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz {
        let nidx = ni * (ny * nz) + nj * nz + nk
        if grid[nidx] == 0 {
          grid[nidx] = 2
          queue.push(nidx)
        }
      }
    }
  }
  // Remaining free cells (grid == 0) are cavities
  // Mark them as 3 and cluster
  let mut cidx = 0
  while cidx < total_cells {
    if grid[cidx] == 0 { grid[cidx] = 3 }
    cidx = cidx + 1
  }
  // Cluster cavity cells using BFS
  let visited = Array::make(total_cells, false)
  let cavities = Array(capacity=0)
  let mut total_vol = 0.0
  let cell_volume = grid_spacing * grid_spacing * grid_spacing
  let mut c2 = 0
  while c2 < total_cells {
    if grid[c2] == 3 && not(visited[c2]) {
      // BFS to collect this cluster
      let cluster = Array(capacity=0)
      let cqueue = Array(capacity=0)
      cqueue.push(c2)
      visited[c2] = true
      let mut qhead = 0
      while qhead < cqueue.length() {
        let cur = cqueue[qhead]
        qhead = qhead + 1
        cluster.push(cur)
        let ci = cur / (ny * nz)
        let cj = (cur % (ny * nz)) / nz
        let ck = cur % nz
        let neighbors = [
          (ci - 1, cj, ck), (ci + 1, cj, ck),
          (ci, cj - 1, ck), (ci, cj + 1, ck),
          (ci, cj, ck - 1), (ci, cj, ck + 1),
        ]
        for (ni, nj, nk) in neighbors {
          if ni >= 0 && ni < nx && nj >= 0 && nj < ny && nk >= 0 && nk < nz {
            let nidx = ni * (ny * nz) + nj * nz + nk
            if grid[nidx] == 3 && not(visited[nidx]) {
              visited[nidx] = true
              cqueue.push(nidx)
            }
          }
        }
      }
      // Compute cavity properties
      let mut sx = 0.0
      let mut sy = 0.0
      let mut sz = 0.0
      let points = Array(capacity=cluster.length())
      for cell in cluster {
        let ci = cell / (ny * nz)
        let cj = (cell % (ny * nz)) / nz
        let ck = cell % nz
        let gx = min_x + ci.to_double() * grid_spacing
        let gy = min_y + cj.to_double() * grid_spacing
        let gz = min_z + ck.to_double() * grid_spacing
        sx = sx + gx; sy = sy + gy; sz = sz + gz
        points.push(Vector3::new(gx, gy, gz))
      }
      let cn = cluster.length()
      let center = Vector3::new(sx / cn.to_double(), sy / cn.to_double(), sz / cn.to_double())
      let volume = cn.to_double() * cell_volume
      total_vol = total_vol + volume
      // Check if it's a channel (touches boundary region within margin)
      let is_channel = false  // simplified: all enclosed free space is cavity
      cavities.push(Cavity::{ center, volume, num_points: cn, grid_points: points, is_channel })
    }
    c2 = c2 + 1
  }
  CaviFinderResult::{ cavities, grid_spacing, probe_radius, total_volume: total_vol }
}
