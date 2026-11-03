///|
/// Spherical harmonics multipole expansion for electrostatic analysis.
///
/// Biotite doesn't directly provide this, but it's used in structure
/// analysis for computing electrostatic properties. This module
/// implements real spherical harmonics and multipole moments.

///|
/// Factorial helper.
fn factorial(n : Int32) -> Double {
  let mut result = 1.0
  for i in 2..<=n {
    result = result * i.to_double()
  }
  result
}

///|
/// Associated Legendre polynomial P_l^m(cos_theta).
///
/// Computed using the recurrence relation.
fn legendre_p(l : Int32, m : Int32, x : Double) -> Double {
  if m > l {
    return 0.0
  }
  // Start with P_m^m
  let mut p_mm = 1.0
  if m > 0 {
    let mut somx2 = (1.0 - x * x).sqrt()
    let mut fact = 1.0
    for _ in 0..<m {
      p_mm = p_mm * -fact * somx2
      fact = fact + 2.0
      somx2 = somx2 // keep value
    }
  }
  if l == m {
    return p_mm
  }
  // P_{m+1}^m = x * (2m+1) * P_m^m
  let mut p_mp1 = x * (2 * m + 1).to_double() * p_mm
  if l == m + 1 {
    return p_mp1
  }
  // Recurrence: P_l^m = ((2l-1)*x*P_{l-1}^m - (l+m-1)*P_{l-2}^m) / (l-m)
  let mut p_ll = 0.0
  for ll in (m + 2)..<=l {
    p_ll = (
        (2 * ll - 1).to_double() * x * p_mp1 - (ll + m - 1).to_double() * p_mm
      ) /
      (ll - m).to_double()
    p_mm = p_mp1
    p_mp1 = p_ll
  }
  p_ll
}

///|
/// Real spherical harmonic Y_l^m(theta, phi).
///
/// theta: polar angle [0, pi]
/// phi: azimuthal angle [0, 2*pi]
pub fn real_spherical_harmonic(
  l : Int32,
  m : Int32,
  theta : Double,
  phi : Double,
) -> Double {
  let abs_m = if m < 0 { -m } else { m }
  let x = @math.cos(theta)
  let p = legendre_p(l, abs_m, x)
  let norm = ((2 * l + 1).to_double() /
  (4.0 * 3.141592653589793) *
  factorial(l - abs_m) /
  factorial(l + abs_m)).sqrt()
  if m > 0 {
    norm * p * @math.cos(m.to_double() * phi) * 1.4142135623730951
  } else if m < 0 {
    norm * p * @math.sin(abs_m.to_double() * phi) * 1.4142135623730951
  } else {
    norm * p
  }
}

///|
/// Multipole moment of a set of point charges.
///
/// Computes the real spherical harmonic multipole moment
/// Q_l^m = sum_i q_i * r_i^l * Y_l^m(theta_i, phi_i)
pub fn multipole_moment(
  charges : Array[(Double, Double, Double, Double)],
  // (x, y, z, charge) for each atom
  l : Int32,
  m : Int32,
) -> Double {
  let mut moment = 0.0
  for atom in charges {
    let (x, y, z, q) = atom
    let r = (x * x + y * y + z * z).sqrt()
    if r < 0.0000000001 {
      continue
    }
    let theta = @math.acos(z / r)
    let phi = @math.atan2(y, x)
    let ylm = real_spherical_harmonic(l, m, theta, phi)
    moment = moment + q * pow_double(r, l) * ylm
  }
  moment
}

///|
/// Integer power of a double.
fn pow_double(base : Double, exp : Int32) -> Double {
  if exp == 0 {
    return 1.0
  }
  let mut result = 1.0
  for _ in 0..<exp {
    result = result * base
  }
  result
}

///|
/// Compute dipole moment (l=1 multipole) of point charges.
///
/// Returns (dipole_x, dipole_y, dipole_z, magnitude).
pub fn dipole_moment(
  charges : Array[(Double, Double, Double, Double)],
) -> (Double, Double, Double, Double) {
  let mut dx = 0.0
  let mut dy = 0.0
  let mut dz = 0.0
  for atom in charges {
    let (x, y, z, q) = atom
    dx = dx + q * x
    dy = dy + q * y
    dz = dz + q * z
  }
  let mag = (dx * dx + dy * dy + dz * dz).sqrt()
  (dx, dy, dz, mag)
}

///|
/// Compute quadrupole moment tensor (l=2) of point charges.
///
/// Returns the 3x3 quadrupole tensor as a flat array.
pub fn quadrupole_moment(
  charges : Array[(Double, Double, Double, Double)],
) -> Array[Double] {
  // Q_ij = sum_k q_k * (3 * r_i * r_j - |r|^2 * delta_ij)
  let q_tensor = Array::make(9, 0.0)
  for atom in charges {
    let (x, y, z, q) = atom
    let r2 = x * x + y * y + z * z
    q_tensor[0] = q_tensor[0] + q * (3.0 * x * x - r2) // Q_xx
    q_tensor[1] = q_tensor[1] + q * (3.0 * x * y) // Q_xy
    q_tensor[2] = q_tensor[2] + q * (3.0 * x * z) // Q_xz
    q_tensor[3] = q_tensor[3] + q * (3.0 * y * x) // Q_yx
    q_tensor[4] = q_tensor[4] + q * (3.0 * y * y - r2) // Q_yy
    q_tensor[5] = q_tensor[5] + q * (3.0 * y * z) // Q_yz
    q_tensor[6] = q_tensor[6] + q * (3.0 * z * x) // Q_zx
    q_tensor[7] = q_tensor[7] + q * (3.0 * z * y) // Q_zy
    q_tensor[8] = q_tensor[8] + q * (3.0 * z * z - r2) // Q_zz
  }
  q_tensor
}

///|
/// Electrostatic potential at a point due to point charges.
///
/// V(r) = sum_i q_i / |r - r_i|
/// Uses distance-dependent Coulomb potential (no dielectric constant).
pub fn electrostatic_potential(
  charges : Array[(Double, Double, Double, Double)],
  point : (Double, Double, Double),
) -> Double {
  let (px, py, pz) = point
  let mut potential = 0.0
  for atom in charges {
    let (x, y, z, q) = atom
    let dx = px - x
    let dy = py - y
    let dz = pz - z
    let dist = (dx * dx + dy * dy + dz * dz).sqrt()
    if dist > 0.0000000001 {
      potential = potential + q / dist
    }
  }
  potential
}
