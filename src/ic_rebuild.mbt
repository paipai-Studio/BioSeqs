///|
/// ic_rebuild - Rebuild Cartesian coordinates from internal coordinates.
/// Inspired by Biopython Bio.PDB.ic_rebuild module.
/// Implements the Natural Extension Reference Frame (NeRF) algorithm:
///   given 3 reference atoms + (bond length, bond angle, dihedral),
///   compute the 4th atom's (x, y, z) Cartesian coordinates.

///|
/// 3D point.
pub(all) struct IcPoint {
  x : Double
  y : Double
  z : Double
}

///|
/// Internal coordinate triple for one atom (relative to 3 predecessors).
pub(all) struct IcInternalCoord {
  atom_name : String
  bond_length : Double
  bond_angle : Double
  dihedral : Double
}

///|
/// Result of rebuilding: atom name + Cartesian position.
pub(all) struct IcRebuiltAtom {
  name : String
  pos : IcPoint
}

///|
/// Construct a point.
pub fn IcPoint::new(x : Double, y : Double, z : Double) -> IcPoint {
  IcPoint::{ x, y, z }
}

///|
/// Euclidean distance between two points.
pub fn ic_point_distance(a : IcPoint, b : IcPoint) -> Double {
  let dx = a.x - b.x
  let dy = a.y - b.y
  let dz = a.z - b.z
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
/// Vector subtraction (a - b).
pub fn ic_point_sub(a : IcPoint, b : IcPoint) -> IcPoint {
  IcPoint::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

///|
/// Dot product of two vectors (treated as vectors from origin).
pub fn ic_point_dot(a : IcPoint, b : IcPoint) -> Double {
  a.x * b.x + a.y * b.y + a.z * b.z
}

///|
/// Cross product a × b.
pub fn ic_point_cross(a : IcPoint, b : IcPoint) -> IcPoint {
  IcPoint::new(
    a.y * b.z - a.z * b.y,
    a.z * b.x - a.x * b.z,
    a.x * b.y - a.y * b.x
  )
}

///|
/// Vector norm (magnitude).
pub fn ic_point_norm(a : IcPoint) -> Double {
  (a.x * a.x + a.y * a.y + a.z * a.z).sqrt()
}

///|
/// Normalize a vector to unit length.
pub fn ic_point_normalize(a : IcPoint) -> IcPoint {
  let n = ic_point_norm(a)
  if n < 1.0e-12 {
    return IcPoint::new(0.0, 0.0, 0.0)
  }
  IcPoint::new(a.x / n, a.y / n, a.z / n)
}

///|
/// NeRF: place atom D given 3 reference atoms A, B, C and internal coords.
///   bond_length = |CD|
///   bond_angle  = angle B-C-D (radians, molecular convention: 180°=collinear)
///   dihedral    = torsion A-B-C-D (radians)
/// Reference: Parsons et al. (2005) "Practical conversion from torsion space to Cartesian space".
pub fn ic_nerf(
  a : IcPoint,
  b : IcPoint,
  c : IcPoint,
  bond_length : Double,
  bond_angle : Double,
  dihedral : Double
) -> IcPoint {
  // Step 1: compute local frame (bc_unit, n, m)
  let bc = ic_point_sub(c, b)
  let bc_unit = ic_point_normalize(bc)
  let ab = ic_point_sub(b, a)
  let n = ic_point_normalize(ic_point_cross(ab, bc_unit))
  let m = ic_point_cross(n, bc_unit)
  // Step 2: position in local frame
  // NeRF theta = pi - bond_angle (molecular convention: 180°=collinear → theta=0)
  let r = bond_length
  let theta = 3.14159265358979 - bond_angle
  let phi = dihedral
  let d2 = IcPoint::new(
    r * @math.cos(theta),
    r * @math.cos(phi) * @math.sin(theta),
    r * @math.sin(phi) * @math.sin(theta)
  )
  // Step 3: transform to global frame: D = C + d2.x * bc_unit + d2.y * m + d2.z * n
  IcPoint::new(
    c.x + d2.x * bc_unit.x + d2.y * m.x + d2.z * n.x,
    c.y + d2.x * bc_unit.y + d2.y * m.y + d2.z * n.y,
    c.z + d2.x * bc_unit.z + d2.y * m.z + d2.z * n.z
  )
}

///|
/// Rebuild a chain of atoms from internal coordinates.
///   The first 3 atoms must be given as Cartesian seeds (positions[0..2]).
///   Each subsequent atom i is placed via NeRF using atoms i-1, i-2, i-3
///   and internal_coords[i-3].
pub fn ic_rebuild_chain(
  seed_positions : Array[IcPoint],
  internal_coords : Array[IcInternalCoord]
) -> Array[IcRebuiltAtom] {
  let result : Array[IcRebuiltAtom] = []
  let n_seed = seed_positions.length()
  // Emit seed atoms
  let mut i = 0
  while i < n_seed {
    result.push(IcRebuiltAtom::{ name: "seed" + i.to_string(), pos: seed_positions[i] })
    i = i + 1
  }
  if n_seed < 3 {
    return result
  }
  // Rebuild remaining atoms
  let mut k = 0
  while k < internal_coords.length() {
    let ic = internal_coords[k]
    let a = seed_positions[n_seed - 3]
    let b = seed_positions[n_seed - 2]
    let c = seed_positions[n_seed - 1]
    let d = ic_nerf(a, b, c, ic.bond_length, ic.bond_angle, ic.dihedral)
    result.push(IcRebuiltAtom::{ name: ic.atom_name, pos: d })
    // Shift window: use last 3 positions as new seed
    seed_positions[n_seed - 3] = seed_positions[n_seed - 2]
    seed_positions[n_seed - 2] = seed_positions[n_seed - 1]
    seed_positions[n_seed - 1] = d
    k = k + 1
  }
  result
}

///|
/// Compute bond length from two Cartesian points.
pub fn ic_compute_bond_length(a : IcPoint, b : IcPoint) -> Double {
  ic_point_distance(a, b)
}

///|
/// Compute bond angle (radians) at vertex b for atoms a-b-c.
pub fn ic_compute_bond_angle(a : IcPoint, b : IcPoint, c : IcPoint) -> Double {
  let ba = ic_point_sub(b, a)
  let bc = ic_point_sub(b, c)
  let dot = ic_point_dot(ba, bc)
  let denom = ic_point_norm(ba) * ic_point_norm(bc)
  if denom < 1.0e-12 {
    return 0.0
  }
  let cos_theta = dot / denom
  // Clamp for numerical safety
  let clamped = if cos_theta > 1.0 { 1.0 } else if cos_theta < -1.0 { -1.0 } else { cos_theta }
  @math.acos(clamped)
}

///|
/// Compute dihedral angle (radians) for atoms a-b-c-d.
pub fn ic_compute_dihedral(
  a : IcPoint,
  b : IcPoint,
  c : IcPoint,
  d : IcPoint
) -> Double {
  // Use the same frame as NeRF to ensure sign consistency.
  let bc = ic_point_sub(c, b)
  let bc_unit = ic_point_normalize(bc)
  let ab = ic_point_sub(b, a)
  let n = ic_point_normalize(ic_point_cross(ab, bc_unit))
  let m = ic_point_cross(n, bc_unit)
  let cd = ic_point_sub(d, c)
  let cd_m = ic_point_dot(cd, m)
  let cd_n = ic_point_dot(cd, n)
  @math.atan2(cd_n, cd_m)
}

///|
/// Round-trip verification: extract internal coords from Cartesian, rebuild, compare.
///   Returns max deviation across all rebuilt atoms.
pub fn ic_verify_roundtrip(positions : Array[IcPoint]) -> Double {
  let n = positions.length()
  if n < 4 {
    return 0.0
  }
  // Extract internal coords for atoms 3..n-1
  let internal_coords : Array[IcInternalCoord] = []
  let mut i = 3
  while i < n {
    let bl = ic_compute_bond_length(positions[i - 1], positions[i])
    let ba = ic_compute_bond_angle(positions[i - 2], positions[i - 1], positions[i])
    let dh = ic_compute_dihedral(positions[i - 3], positions[i - 2], positions[i - 1], positions[i])
    internal_coords.push(IcInternalCoord::{ atom_name: "atom" + i.to_string(), bond_length: bl, bond_angle: ba, dihedral: dh })
    i = i + 1
  }
  // Rebuild
  let seeds : Array[IcPoint] = [positions[0], positions[1], positions[2]]
  let rebuilt = ic_rebuild_chain(seeds, internal_coords)
  // Compare
  let mut max_dev = 0.0
  let mut j = 3
  while j < n {
    let dev = ic_point_distance(positions[j], rebuilt[j].pos)
    if dev > max_dev {
      max_dev = dev
    }
    j = j + 1
  }
  max_dev
}

///|
/// Example: build a small tetrahedral-like chain and verify round-trip.
pub fn ic_rebuild_example() -> Array[IcRebuiltAtom] {
  // Seed: 3 atoms along x-axis
  let p0 = IcPoint::new(0.0, 0.0, 0.0)
  let p1 = IcPoint::new(1.5, 0.0, 0.0)
  let p2 = IcPoint::new(2.5, 1.0, 0.0)
  // Internal coords for 3 more atoms (ideal tetrahedral-ish)
  let ics : Array[IcInternalCoord] = [
    IcInternalCoord::{ atom_name: "CA", bond_length: 1.46, bond_angle: 1.91, dihedral: -0.52 },
    IcInternalCoord::{ atom_name: "CB", bond_length: 1.53, bond_angle: 1.91, dihedral: 1.57 },
    IcInternalCoord::{ atom_name: "CG", bond_length: 1.50, bond_angle: 1.91, dihedral: -1.05 },
  ]
  ic_rebuild_chain([p0, p1, p2], ics)
}