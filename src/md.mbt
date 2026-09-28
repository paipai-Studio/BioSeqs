///|
/// Molecular dynamics — a minimal engine: harmonic bond + Lennard-Jones +
/// Coulomb potential energies and analytic forces, plus velocity-Verlet
/// integration. Units are kcal/mol (energy), Å (length) and ps (time).

///|
fn md_dist(a : Array[Double], b : Array[Double]) -> Double {
  let dx = a[0] - b[0]
  let dy = a[1] - b[1]
  let dz = a[2] - b[2]
  (dx * dx + dy * dy + dz * dz).sqrt()
}

///|
fn md_new_forces(n : Int) -> Array[Array[Double]] {
  let f : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    f[i] = [0.0, 0.0, 0.0]
  }
  f
}

///|
/// Harmonic bond energy: `sum 0.5 k (r - r0)^2`. `bonds` is `(i, j, r0, k)`.
pub fn md_bond_energy(
  positions : Array[Array[Double]],
  bonds : Array[(Int, Int, Double, Double)],
) -> Double {
  let mut e = 0.0
  for (i, j, r0, k) in bonds {
    let r = md_dist(positions[i], positions[j])
    let dr = r - r0
    e = e + 0.5 * k * dr * dr
  }
  e
}

///|
/// Analytic harmonic bond forces (negative gradient). `bonds` is
/// `(i, j, r0, k)`; returns one `[fx, fy, fz]` per atom.
pub fn md_bond_forces(
  positions : Array[Array[Double]],
  bonds : Array[(Int, Int, Double, Double)],
) -> Array[Array[Double]] {
  let forces = md_new_forces(positions.length())
  for (i, j, r0, k) in bonds {
    let r = md_dist(positions[i], positions[j])
    if r == 0.0 {
      continue
    }
    let dr = r - r0
    let magnitude = k * dr // +k沿 i→j 压缩恢复，F_i 方向背离 j
    for d = 0; d < 3; d = d + 1 {
      let dir = (positions[i][d] - positions[j][d]) / r
      let f = -magnitude * dir // 作用在 i 上的力
      forces[i][d] = forces[i][d] + f
      forces[j][d] = forces[j][d] - f
    }
  }
  forces
}

///|
/// Lennard-Jones energy: `4 eps ((sigma/r)^12 - (sigma/r)^6)`. `pairs` is
/// `(i, j, epsilon, sigma)`.
pub fn md_lj_energy(
  positions : Array[Array[Double]],
  pairs : Array[(Int, Int, Double, Double)],
) -> Double {
  let mut e = 0.0
  for (i, j, epsilon, sigma) in pairs {
    let r = md_dist(positions[i], positions[j])
    if r == 0.0 {
      continue
    }
    let sr = sigma / r
    let sr6 = sr * sr * sr * sr * sr * sr
    let sr12 = sr6 * sr6
    e = e + 4.0 * epsilon * (sr12 - sr6)
  }
  e
}

///|
/// Analytic Lennard-Jones forces (negative gradient).
pub fn md_lj_forces(
  positions : Array[Array[Double]],
  pairs : Array[(Int, Int, Double, Double)],
) -> Array[Array[Double]] {
  let forces = md_new_forces(positions.length())
  for (i, j, epsilon, sigma) in pairs {
    let r = md_dist(positions[i], positions[j])
    if r == 0.0 {
      continue
    }
    let sr = sigma / r
    let sr6 = sr * sr * sr * sr * sr * sr
    let sr12 = sr6 * sr6
    // dE/dr = 24 eps (sr6 - 2 sr12)/r；力幅 = -dE/dr = 24 eps (2 sr12 - sr6)/r
    let mag = 24.0 * epsilon * (2.0 * sr12 - sr6) / r
    for d = 0; d < 3; d = d + 1 {
      let dir = (positions[i][d] - positions[j][d]) / r
      let f = mag * dir
      forces[i][d] = forces[i][d] + f
      forces[j][d] = forces[j][d] - f
    }
  }
  forces
}

///|
/// Coulomb energy: `k_e q_i q_j / r` (k_e = 332.0637 kcal·Å / (mol·e²)).
/// `pairs` is `(i, j, qi, qj)`.
pub fn md_coulomb_energy(
  positions : Array[Array[Double]],
  pairs : Array[(Int, Int, Double, Double)],
) -> Double {
  let ke = 332.0637
  let mut e = 0.0
  for (i, j, qi, qj) in pairs {
    let r = md_dist(positions[i], positions[j])
    if r == 0.0 {
      continue
    }
    e = e + ke * qi * qj / r
  }
  e
}

///|
/// Coulomb forces (negative gradient): `k_e q_i q_j / r^2` along the bond.
pub fn md_coulomb_forces(
  positions : Array[Array[Double]],
  pairs : Array[(Int, Int, Double, Double)],
) -> Array[Array[Double]] {
  let ke = 332.0637
  let forces = md_new_forces(positions.length())
  for (i, j, qi, qj) in pairs {
    let r = md_dist(positions[i], positions[j])
    if r == 0.0 {
      continue
    }
    let mag = ke * qi * qj / (r * r)
    for d = 0; d < 3; d = d + 1 {
      let dir = (positions[i][d] - positions[j][d]) / r
      let f = mag * dir
      forces[i][d] = forces[i][d] + f
      forces[j][d] = forces[j][d] - f
    }
  }
  forces
}

///|
/// Kinetic energy `0.5 sum m v^2`.
pub fn md_kinetic_energy(
  velocities : Array[Array[Double]],
  masses : Array[Double],
) -> Double {
  let mut e = 0.0
  for i = 0; i < velocities.length(); i = i + 1 {
    let v = velocities[i]
    e = e + 0.5 * masses[i] * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
  }
  e
}


///|
/// Velocity-Verlet first half: returns `(new_positions, half_velocities)`.
/// `v_half = v + 0.5 dt F/m`, `x_new = x + dt v_half`.
pub fn md_verlet_move(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  forces : Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let n = positions.length()
  let x_new = md_new_forces(n)
  let v_half = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      v_half[i][d] = velocities[i][d] + 0.5 * dt * forces[i][d] / masses[i]
      x_new[i][d] = positions[i][d] + dt * v_half[i][d]
    }
  }
  (x_new, v_half)
}

///|
/// Velocity-Verlet second half: completes the velocities using the forces at
/// the new positions. `v_new = v_half + 0.5 dt F_new/m`.
pub fn md_verlet_finish(
  v_half : Array[Array[Double]],
  new_forces : Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
) -> Array[Array[Double]] {
  let n = v_half.length()
  let v_new = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      v_new[i][d] = v_half[i][d] + 0.5 * dt * new_forces[i][d] / masses[i]
    }
  }
  v_new
}

///|
/// Run a harmonic-bond simulation and return the total energy (kinetic +
/// potential) at each step. Under Verlet integration the total energy is
/// conserved (no thermostat), so the series should stay nearly constant.
pub fn md_simulate_harmonic(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  masses : Array[Double],
  bonds : Array[(Int, Int, Double, Double)],
  dt : Double,
  steps : Int,
) -> Array[Double] {
  let energies : Array[Double] = Array(capacity=0)
  let mut pos = positions
  let mut vel = velocities
  for _s = 0; _s < steps; _s = _s + 1 {
    let forces = md_bond_forces(pos, bonds)
    let pe = md_bond_energy(pos, bonds)
    let ke = md_kinetic_energy(vel, masses)
    energies.push(pe + ke)
    let (pos_new, v_half) = md_verlet_move(pos, vel, forces, masses, dt)
    let new_forces = md_bond_forces(pos_new, bonds)
    let vel_new = md_verlet_finish(v_half, new_forces, masses, dt)
    pos = pos_new
    vel = vel_new
  }
  energies
}


// ============================================================================
// Angle bending + dihedral (proper torsion) terms
// ============================================================================

///|
fn md_vec_sub(a : Array[Double], b : Array[Double]) -> Array[Double] {
  [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

///|
fn md_dot(a : Array[Double], b : Array[Double]) -> Double {
  a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

///|
fn md_norm(a : Array[Double]) -> Double {
  (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

///|
fn md_cross(a : Array[Double], b : Array[Double]) -> Array[Double] {
  [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ]
}

///|
/// Bond angle at `b` (center) spanned by `a` and `c`, in radians.
fn md_angle_value(a : Array[Double], b : Array[Double], c : Array[Double]) -> Double {
  let v1 = md_vec_sub(a, b)
  let v2 = md_vec_sub(c, b)
  let d1 = md_norm(v1)
  let d2 = md_norm(v2)
  if d1 == 0.0 || d2 == 0.0 {
    return 0.0
  }
  let mut cosv = md_dot(v1, v2) / (d1 * d2)
  if cosv > 1.0 {
    cosv = 1.0
  } else if cosv < -1.0 {
    cosv = -1.0
  }
  @math.acos(cosv)
}

///|
/// Harmonic angle energy: `sum 0.5 k_theta (theta - theta0)^2`. `angles` is
/// `(i, j, k, theta0, k_theta)` with `j` the central atom, `theta0` in radians.
pub fn md_angle_energy(
  positions : Array[Array[Double]],
  angles : Array[(Int, Int, Int, Double, Double)],
) -> Double {
  let mut e = 0.0
  for (i, j, k, theta0, kt) in angles {
    let theta = md_angle_value(positions[i], positions[j], positions[k])
    let d = theta - theta0
    e = e + 0.5 * kt * d * d
  }
  e
}

///|
/// Analytic angle forces (negative gradient). Force on each of `i`, `j`, `k`.
pub fn md_angle_forces(
  positions : Array[Array[Double]],
  angles : Array[(Int, Int, Int, Double, Double)],
) -> Array[Array[Double]] {
  let forces = md_new_forces(positions.length())
  for (i, j, k, theta0, kt) in angles {
    let v1 = md_vec_sub(positions[i], positions[j])
    let v2 = md_vec_sub(positions[k], positions[j])
    let d1 = md_norm(v1)
    let d2 = md_norm(v2)
    if d1 == 0.0 || d2 == 0.0 {
      continue
    }
    let mut cosv = md_dot(v1, v2) / (d1 * d2)
    if cosv > 1.0 {
      cosv = 1.0
    } else if cosv < -1.0 {
      cosv = -1.0
    }
    let theta = @math.acos(cosv)
    let mut sinv = (1.0 - cosv * cosv).sqrt()
    if sinv < 1.0e-12 {
      sinv = 1.0e-12
    }
    let u = [v1[0] / d1, v1[1] / d1, v1[2] / d1]
    let w = [v2[0] / d2, v2[1] / d2, v2[2] / d2]
    let fc = kt * (theta - theta0) / sinv
    // F_i = -fc (u cos - w) / d1；F_k = -fc (w cos - u) / d2
    for d = 0; d < 3; d = d + 1 {
      let f_i = -fc * (u[d] * cosv - w[d]) / d1
      let f_k = -fc * (w[d] * cosv - u[d]) / d2
      forces[i][d] = forces[i][d] + f_i
      forces[k][d] = forces[k][d] + f_k
      forces[j][d] = forces[j][d] - f_i - f_k
    }
  }
  forces
}

///|
/// Dihedral angle for atoms `i-j-k-l` (torsion around the `j-k` axis), in
/// radians. Uses the atan2 cross-product formulation, range `(-pi, pi]`.
fn md_dihedral_value(
  a : Array[Double],
  b : Array[Double],
  c : Array[Double],
  d : Array[Double],
) -> Double {
  let bvec = md_vec_sub(c, b)   // j -> k
  let nb = md_norm(bvec)
  if nb == 0.0 {
    return 0.0
  }
  let v1 = md_vec_sub(a, b)     // j -> i
  let v2 = md_vec_sub(d, c)     // k -> l
  let n1 = md_cross(v1, bvec)   // normal of plane i-j-k
  let n2 = md_cross(bvec, v2)   // normal of plane j-k-l
  @math.atan2(nb * md_dot(v1, n2), md_dot(n1, n2))
}

///|
/// Dihedral (proper torsion) energy: `sum k_phi (1 + cos(n phi - delta))`.
/// `dihedrals` is `(i, j, k, l, k_phi, n, delta)` with `n` the periodicity and
/// `delta` the phase in radians.
pub fn md_dihedral_energy(
  positions : Array[Array[Double]],
  dihedrals : Array[(Int, Int, Int, Int, Double, Int, Double)],
) -> Double {
  let mut e = 0.0
  for (i, j, k, l, kphi, n, delta) in dihedrals {
    let phi = md_dihedral_value(positions[i], positions[j], positions[k], positions[l])
    e = e + kphi * (1.0 + @math.cos(n.to_double() * phi - delta))
  }
  e
}

///|
fn md_copy_positions(positions : Array[Array[Double]]) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=0)
  for p in positions {
    out.push(p.copy())
  }
  out
}

///|
/// Dihedral forces computed by central-difference numerical gradient of
/// `md_dihedral_energy` (the analytic torsion derivatives are error-prone; the
/// numerical one is exact to 2nd order and used for the force). This is the
/// documented exception to the analytic-force convention in this module.
pub fn md_dihedral_forces(
  positions : Array[Array[Double]],
  dihedrals : Array[(Int, Int, Int, Int, Double, Int, Double)],
) -> Array[Array[Double]] {
  let n = positions.length()
  let forces = md_new_forces(n)
  let h = 1.0e-6
  let base = positions
  let mut i = 0
  while i < n {
    let mut d = 0
    while d < 3 {
      let plus = md_copy_positions(base)
      let minus = md_copy_positions(base)
      plus[i][d] = plus[i][d] + h
      minus[i][d] = minus[i][d] - h
      let ep = md_dihedral_energy(plus, dihedrals)
      let em = md_dihedral_energy(minus, dihedrals)
      forces[i][d] = -(ep - em) / (2.0 * h)
      d = d + 1
    }
    i = i + 1
  }
  forces
}
