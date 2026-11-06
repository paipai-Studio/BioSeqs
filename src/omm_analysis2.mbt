/// omm_analysis2.mbt — Pressure/virial tensor + moment of inertia + energy decomposition
/// Items 15-17: Advanced analysis

///|
// === Pressure / Virial Tensor ===

///|
/// Compute the virial tensor (3x3) from positions and forces.
/// V_ab = sum_i r_i,a * F_i,b
pub fn omm_virial_tensor(
  coords : Array[(Double, Double, Double)],
  forces : Array[(Double, Double, Double)],
) -> Array[Array[Double]] {
  let v : Array[Array[Double]] = Array(capacity=3)
  for a = 0; a < 3; a = a + 1 {
    let row : Array[Double] = Array(capacity=3)
    for b = 0; b < 3; b = b + 1 {
      row.push(0.0)
    }
    v.push(row)
  }
  for i = 0; i < coords.length(); i = i + 1 {
    if i >= forces.length() {
      continue
    }
    let r = coords[i]
    let f = forces[i]
    v[0][0] = v[0][0] + r.0 * f.0
    v[0][1] = v[0][1] + r.0 * f.1
    v[0][2] = v[0][2] + r.0 * f.2
    v[1][0] = v[1][0] + r.1 * f.0
    v[1][1] = v[1][1] + r.1 * f.1
    v[1][2] = v[1][2] + r.1 * f.2
    v[2][0] = v[2][0] + r.2 * f.0
    v[2][1] = v[2][1] + r.2 * f.1
    v[2][2] = v[2][2] + r.2 * f.2
  }
  v
}

///|
/// Compute pressure tensor from kinetic energy and virial.
/// P_ab = (1/V) * [sum_i m_i * v_i,a * v_i,b + virial_ab]
pub fn omm_pressure_tensor(
  masses : Array[Double],
  velocities : Array[(Double, Double, Double)],
  coords : Array[(Double, Double, Double)],
  forces : Array[(Double, Double, Double)],
  volume : Double,
) -> Array[Array[Double]] {
  let virial = omm_virial_tensor(coords, forces)
  // Kinetic contribution
  let kinetic : Array[Array[Double]] = Array(capacity=3)
  for a = 0; a < 3; a = a + 1 {
    let row : Array[Double] = Array(capacity=3)
    for b = 0; b < 3; b = b + 1 {
      row.push(0.0)
    }
    kinetic.push(row)
  }
  for i = 0; i < velocities.length(); i = i + 1 {
    if i >= masses.length() {
      continue
    }
    let m = masses[i]
    let v = velocities[i]
    kinetic[0][0] = kinetic[0][0] + m * v.0 * v.0
    kinetic[0][1] = kinetic[0][1] + m * v.0 * v.1
    kinetic[0][2] = kinetic[0][2] + m * v.0 * v.2
    kinetic[1][0] = kinetic[1][0] + m * v.1 * v.0
    kinetic[1][1] = kinetic[1][1] + m * v.1 * v.1
    kinetic[1][2] = kinetic[1][2] + m * v.1 * v.2
    kinetic[2][0] = kinetic[2][0] + m * v.2 * v.0
    kinetic[2][1] = kinetic[2][1] + m * v.2 * v.1
    kinetic[2][2] = kinetic[2][2] + m * v.2 * v.2
  }
  let inv_vol = 1.0 / (volume + 0.0000000001)
  let pressure : Array[Array[Double]] = Array(capacity=3)
  for a = 0; a < 3; a = a + 1 {
    let row : Array[Double] = Array(capacity=3)
    for b = 0; b < 3; b = b + 1 {
      row.push(inv_vol * (kinetic[a][b] + virial[a][b]))
    }
    pressure.push(row)
  }
  pressure
}

///|
/// Compute scalar (isotropic) pressure from pressure tensor.
pub fn omm_scalar_pressure(pressure : Array[Array[Double]]) -> Double {
  (pressure[0][0] + pressure[1][1] + pressure[2][2]) / 3.0
}

///|
// === Moment of Inertia Tensor ===

///|
/// Compute the 3x3 moment of inertia tensor.
/// I_ab = sum_i m_i * (r_i^2 * delta_ab - r_i,a * r_i,b)
pub fn omm_moment_of_inertia(
  masses : Array[Double],
  coords : Array[(Double, Double, Double)],
) -> Array[Array[Double]] {
  let moi : Array[Array[Double]] = Array(capacity=3)
  for a = 0; a < 3; a = a + 1 {
    let row : Array[Double] = Array(capacity=3)
    for b = 0; b < 3; b = b + 1 {
      row.push(0.0)
    }
    moi.push(row)
  }
  for i = 0; i < coords.length(); i = i + 1 {
    if i >= masses.length() {
      continue
    }
    let m = masses[i]
    let (x, y, z) = coords[i]
    let r2 = x * x + y * y + z * z
    moi[0][0] = moi[0][0] + m * (r2 - x * x)
    moi[0][1] = moi[0][1] - m * x * y
    moi[0][2] = moi[0][2] - m * x * z
    moi[1][0] = moi[1][0] - m * y * x
    moi[1][1] = moi[1][1] + m * (r2 - y * y)
    moi[1][2] = moi[1][2] - m * y * z
    moi[2][0] = moi[2][0] - m * z * x
    moi[2][1] = moi[2][1] - m * z * y
    moi[2][2] = moi[2][2] + m * (r2 - z * z)
  }
  moi
}

///|
/// Compute principal moments of inertia (eigenvalues of MOI tensor).
/// Simplified: returns diagonal elements (full eigendecomposition needs Jacobi).
pub fn omm_principal_moments(
  masses : Array[Double],
  coords : Array[(Double, Double, Double)],
) -> (Double, Double, Double) {
  let moi = omm_moment_of_inertia(masses, coords)
  // Simplified: use trace and invariants to estimate eigenvalues
  let trace = moi[0][0] + moi[1][1] + moi[2][2]
  // For diagonal-dominant case, diagonal elements are good approximations
  let i1 = moi[0][0]
  let i2 = moi[1][1]
  let i3 = moi[2][2]
  ignore(trace)
  (i1, i2, i3)
}

///|
// === Energy Decomposition ===

///|
/// Energy breakdown by force group
pub struct OmmEnergyDecomposition {
  bond_energy : Double
  angle_energy : Double
  torsion_energy : Double
  nonbonded_energy : Double
  restraint_energy : Double
  total : Double
}

///|
/// Decompose total energy into force group contributions.
pub fn omm_decompose_energy(
  bond_energy : Double,
  angle_energy : Double,
  torsion_energy : Double,
  nonbonded_energy : Double,
  restraint_energy~ : Double = 0.0,
) -> OmmEnergyDecomposition {
  let total = bond_energy +
    angle_energy +
    torsion_energy +
    nonbonded_energy +
    restraint_energy
  OmmEnergyDecomposition::{
    bond_energy,
    angle_energy,
    torsion_energy,
    nonbonded_energy,
    restraint_energy,
    total,
  }
}

///|
/// Compute per-atom energy contribution (approximate: nonbonded only).
pub fn omm_per_atom_energy(
  charges : Array[Double],
  coords : Array[(Double, Double, Double)],
  epsilon : Array[Double],
  sigma : Array[Double],
) -> Array[Double] {
  let n = coords.length()
  let energies : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let mut e_i = 0.0
    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }
      let dx = coords[i].0 - coords[j].0
      let dy = coords[i].1 - coords[j].1
      let dz = coords[i].2 - coords[j].2
      let r = @math.pow(dx * dx + dy * dy + dz * dz, 0.5)
      if r < 0.0000000001 {
        continue
      }
      // LJ: 4*eps*((sig/r)^12 - (sig/r)^6)
      let sig_r = sigma[i] / r
      let sr6 = sig_r * sig_r * sig_r * sig_r * sig_r * sig_r
      let sr12 = sr6 * sr6
      e_i = e_i + 4.0 * epsilon[i] * (sr12 - sr6)
      // Coulomb
      e_i = e_i + 138.935485 * charges[i] * charges[j] / r
    }
    energies.push(e_i * 0.5) // each pair counted twice
  }
  energies
}

///|
pub fn omm_analysis2_dummy() -> Unit {
  ignore(OmmEnergyDecomposition::{
    bond_energy: 0.0,
    angle_energy: 0.0,
    torsion_energy: 0.0,
    nonbonded_energy: 0.0,
    restraint_energy: 0.0,
    total: 0.0,
  })
}
