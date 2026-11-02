///|
/// Atom density analysis — local and global atom density in 3D structures.
///
/// Biotite provides atom density calculations for evaluating packing
/// quality and identifying under/over-packed regions.

///|
/// Calculate local atom density around a point.
///
/// Returns the number of atoms within a given radius of the center point,
/// divided by the spherical volume (4/3 * pi * r^3).
pub fn local_atom_density(
  center : Vector3,
  coords : Array[Vector3],
  radius : Double,
) -> Double {
  let mut count = 0
  for c in coords {
    if center.distance(c) <= radius {
      count = count + 1
    }
  }
  let volume = 4.0 / 3.0 * 3.141592653589793 * radius * radius * radius
  if volume > 0.0 {
    count.to_double() / volume
  } else {
    0.0
  }
}

///|
/// Calculate atom density for each residue in a structure.
///
/// For each residue center, computes the local atom density within the
/// given radius. Returns an array of (residue_index, density) pairs.
pub fn residue_atom_density(
  residue_centers : Array[Vector3],
  all_coords : Array[Vector3],
  radius~ : Double = 5.0,
) -> Array[(Int32, Double)] {
  let result : Array[(Int32, Double)] = []
  for i in 0..<residue_centers.length() {
    let d = local_atom_density(residue_centers[i], all_coords, radius)
    result.push((i, d))
  }
  result
}

///|
/// Calculate global atom density.
///
/// Total atom count divided by the bounding box volume.
pub fn global_atom_density(coords : Array[Vector3]) -> Double {
  let n = coords.length()
  if n == 0 {
    return 0.0
  }
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
  let volume = (max_x - min_x) * (max_y - min_y) * (max_z - min_z)
  if volume > 0.0 {
    n.to_double() / volume
  } else {
    0.0
  }
}

///|
/// Find under-packed residues (density below threshold).
pub fn find_underpacked_residues(
  densities : Array[(Int32, Double)],
  threshold : Double,
) -> Array[Int32] {
  let result : Array[Int32] = []
  for pair in densities {
    if pair.1 < threshold {
      result.push(pair.0)
    }
  }
  result
}

///|
/// Find over-packed residues (density above threshold).
pub fn find_overpacked_residues(
  densities : Array[(Int32, Double)],
  threshold : Double,
) -> Array[Int32] {
  let result : Array[Int32] = []
  for pair in densities {
    if pair.1 > threshold {
      result.push(pair.0)
    }
  }
  result
}

///|
/// Voronoi-like cell volume estimate for an atom.
///
/// Uses a simple distance-based approach: the volume is approximated
/// as the cube of the distance to the nearest neighbor.
pub fn estimated_voronoi_volume(
  atom_idx : Int32,
  coords : Array[Vector3],
) -> Double {
  if atom_idx >= coords.length() {
    return 0.0
  }
  let mut min_dist = 100000000000000000000.0
  for i in 0..<coords.length() {
    if i != atom_idx {
      let d = coords[atom_idx].distance(coords[i])
      if d < min_dist {
        min_dist = d
      }
    }
  }
  min_dist * min_dist * min_dist
}
