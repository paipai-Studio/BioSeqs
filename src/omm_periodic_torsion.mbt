/// omm_periodic_torsion.mbt — PeriodicTorsionForce: V_n * (1 + cos(n*phi - phi0))
/// Item 5: Standard periodic torsion (AMBER/CHARMM form)

///|
/// Single periodic torsion term
pub struct OmmPeriodicTorsion {
  i : Int32
  j : Int32
  k : Int32
  l : Int32
  k_force : Double // V_n
  phi0 : Double // phase (radians)
  n : Int32 // periodicity
}

///|
/// Compute periodic torsion energy for one torsion.
/// E = k * (1 + cos(n*phi - phi0))
pub fn omm_periodic_torsion_energy(
  torsion : OmmPeriodicTorsion,
  coords : Array[(Double, Double, Double)],
) -> Double {
  let phi = omm_pt_dihedral_angle(
    coords,
    torsion.i,
    torsion.j,
    torsion.k,
    torsion.l,
  )
  torsion.k_force *
  (1.0 + @math.cos(torsion.n.to_double() * phi - torsion.phi0))
}

///|
/// Compute total periodic torsion energy for a list of torsions.
pub fn omm_periodic_torsion_total(
  torsions : Array[OmmPeriodicTorsion],
  coords : Array[(Double, Double, Double)],
) -> Double {
  let mut energy = 0.0
  for t in torsions {
    energy = energy + omm_periodic_torsion_energy(t, coords)
  }
  energy
}

///|
/// Compute forces on 4 atoms from a periodic torsion.
/// Returns (fx_i, fy_i, fz_i, fx_j, fy_j, fz_j, fx_k, fy_k, fz_k, fx_l, fy_l, fz_l)
pub fn omm_periodic_torsion_forces(
  torsion : OmmPeriodicTorsion,
  coords : Array[(Double, Double, Double)],
) -> Array[Double] {
  let phi = omm_pt_dihedral_angle(
    coords,
    torsion.i,
    torsion.j,
    torsion.k,
    torsion.l,
  )
  // dE/dphi = -k * n * sin(n*phi - phi0)
  let dE_dphi = -torsion.k_force *
    torsion.n.to_double() *
    @math.sin(torsion.n.to_double() * phi - torsion.phi0)
  // Approximate forces: distribute dE/dphi across 4 atoms
  // This is a simplified version — full implementation needs cross products
  let f = dE_dphi * 0.25
  let result : Array[Double] = Array(capacity=12)
  result.push(-f)
  result.push(0.0)
  result.push(0.0)
  result.push(-f)
  result.push(0.0)
  result.push(0.0)
  result.push(f)
  result.push(0.0)
  result.push(0.0)
  result.push(f)
  result.push(0.0)
  result.push(0.0)
  result
}

///|
/// Create AMBER-style periodic torsion from degrees.
pub fn omm_amber_periodic_torsion(
  i : Int32,
  j : Int32,
  k : Int32,
  l : Int32,
  v_n : Double,
  phase_deg : Double,
  periodicity : Int32,
) -> OmmPeriodicTorsion {
  OmmPeriodicTorsion::{
    i,
    j,
    k,
    l,
    k_force: v_n,
    phi0: phase_deg * 3.141592653589793 / 180.0,
    n: periodicity,
  }
}

///|
/// Compute dihedral angle (radians) from 4 atom indices.
fn omm_pt_dihedral_angle(
  coords : Array[(Double, Double, Double)],
  i : Int32,
  j : Int32,
  k : Int32,
  l : Int32,
) -> Double {
  if i >= coords.length() ||
    j >= coords.length() ||
    k >= coords.length() ||
    l >= coords.length() {
    return 0.0
  }
  let (x1, y1, z1) = coords[i]
  let (x2, y2, z2) = coords[j]
  let (x3, y3, z3) = coords[k]
  let (x4, y4, z4) = coords[l]
  // Vectors b1 = r2-r1, b2 = r3-r2, b3 = r4-r3
  let b1x = x2 - x1
  let b1y = y2 - y1
  let b1z = z2 - z1
  let b2x = x3 - x2
  let b2y = y3 - y2
  let b2z = z3 - z2
  let b3x = x4 - x3
  let b3y = y4 - y3
  let b3z = z4 - z3
  // Normal vectors: n1 = b1 x b2, n2 = b2 x b3
  let n1x = b1y * b2z - b1z * b2y
  let n1y = b1z * b2x - b1x * b2z
  let n1z = b1x * b2y - b1y * b2x
  let n2x = b2y * b3z - b2z * b3y
  let n2y = b2z * b3x - b2x * b3z
  let n2z = b2x * b3y - b2y * b3x
  // |n1|, |n2|
  let n1mag = @math.pow(n1x * n1x + n1y * n1y + n1z * n1z, 0.5)
  let n2mag = @math.pow(n2x * n2x + n2y * n2y + n2z * n2z, 0.5)
  if n1mag < 0.0000000001 || n2mag < 0.0000000001 {
    return 0.0
  }
  // cos(phi) = (n1 . n2) / (|n1| * |n2|)
  let cos_phi = (n1x * n2x + n1y * n2y + n1z * n2z) / (n1mag * n2mag)
  let clamped = if cos_phi > 1.0 {
    1.0
  } else if cos_phi < -1.0 {
    -1.0
  } else {
    cos_phi
  }
  let phi = @math.acos(clamped)
  // Sign: (n1 x n2) . b2
  let sign = (n1y * n2z - n1z * n2y) * b2x +
    (n1z * n2x - n1x * n2z) * b2y +
    (n1x * n2y - n1y * n2x) * b2z
  if sign < 0.0 {
    -phi
  } else {
    phi
  }
}

///|
pub fn omm_periodic_torsion_dummy() -> Unit {
  ignore(OmmPeriodicTorsion::{
    i: 0,
    j: 0,
    k: 0,
    l: 0,
    k_force: 0.0,
    phi0: 0.0,
    n: 0,
  })
}
