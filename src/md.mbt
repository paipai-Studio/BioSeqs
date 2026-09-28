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
