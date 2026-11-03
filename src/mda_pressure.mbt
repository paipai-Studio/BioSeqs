/// Pressure tensor and virial calculation.
///
/// MDAnalysis equivalent: `analysis.pressure.PressureTensor`
///
/// P = (1/V) * (sum(m_i * v_i ⊗ v_i) + virial)
/// where virial = sum(r_ij ⊗ F_ij)

///|
/// Compute kinetic contribution to pressure tensor.
/// velocities: [atom][xyz], masses: [atom]
fn mda_kinetic_pressure_tensor(
  velocities : Array[Array[Double]],
  masses : Array[Double],
) -> Array[Array[Double]] {
  let tensor : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
  ]
  for i = 0; i < velocities.length(); i = i + 1 {
    let m = if i < masses.length() { masses[i] } else { 1.0 }
    let vx = velocities[i][0]
    let vy = velocities[i][1]
    let vz = velocities[i][2]
    tensor[0][0] = tensor[0][0] + m * vx * vx
    tensor[0][1] = tensor[0][1] + m * vx * vy
    tensor[0][2] = tensor[0][2] + m * vx * vz
    tensor[1][0] = tensor[1][0] + m * vy * vx
    tensor[1][1] = tensor[1][1] + m * vy * vy
    tensor[1][2] = tensor[1][2] + m * vy * vz
    tensor[2][0] = tensor[2][0] + m * vz * vx
    tensor[2][1] = tensor[2][1] + m * vz * vy
    tensor[2][2] = tensor[2][2] + m * vz * vz
  }
  tensor
}

///|
/// Compute virial contribution from pairwise forces.
/// positions: [atom][xyz], forces: [atom][xyz]
fn mda_virial_tensor(
  positions : Array[Array[Double]],
  forces : Array[Array[Double]],
) -> Array[Array[Double]] {
  let tensor : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
  ]
  for i = 0; i < positions.length(); i = i + 1 {
    let x = positions[i][0]
    let y = positions[i][1]
    let z = positions[i][2]
    let fx = forces[i][0]
    let fy = forces[i][1]
    let fz = forces[i][2]
    // virial = sum(r_i ⊗ F_i)
    tensor[0][0] = tensor[0][0] + x * fx
    tensor[0][1] = tensor[0][1] + x * fy
    tensor[0][2] = tensor[0][2] + x * fz
    tensor[1][0] = tensor[1][0] + y * fx
    tensor[1][1] = tensor[1][1] + y * fy
    tensor[1][2] = tensor[1][2] + y * fz
    tensor[2][0] = tensor[2][0] + z * fx
    tensor[2][1] = tensor[2][1] + z * fy
    tensor[2][2] = tensor[2][2] + z * fz
  }
  tensor
}

///|
/// Compute full pressure tensor.
///
/// positions: [atom][xyz]
/// velocities: [atom][xyz]
/// forces: [atom][xyz]
/// masses: [atom]
/// volume: system volume
pub fn mda_pressure_tensor(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  forces : Array[Array[Double]],
  masses : Array[Double],
  volume : Double,
) -> Array[Array[Double]] {
  if volume <= 0.0 {
    return [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  }
  let kinetic = mda_kinetic_pressure_tensor(velocities, masses)
  let virial = mda_virial_tensor(positions, forces)
  // P = (1/V) * (kinetic + virial)
  let result : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
  ]
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      result[i][j] = (kinetic[i][j] + virial[i][j]) / volume
    }
  }
  result
}

///|
/// Compute scalar (isotropic) pressure from tensor.
/// P = trace(P_tensor) / 3
pub fn mda_scalar_pressure(tensor : Array[Array[Double]]) -> Double {
  (tensor[0][0] + tensor[1][1] + tensor[2][2]) / 3.0
}

///|
/// Compute pressure anisotropy: P_parallel - P_perpendicular.
/// For membrane systems: P_parallel = (P_xx + P_yy)/2, P_perp = P_zz
pub fn mda_pressure_anisotropy(tensor : Array[Array[Double]]) -> Double {
  let p_parallel = (tensor[0][0] + tensor[1][1]) / 2.0
  let p_perp = tensor[2][2]
  p_parallel - p_perp
}

///|
/// Compute surface tension from pressure tensor.
/// gamma = L_z * (P_parallel - P_perp) / 2
pub fn mda_surface_tension(
  tensor : Array[Array[Double]],
  box_z : Double,
) -> Double {
  let p_parallel = (tensor[0][0] + tensor[1][1]) / 2.0
  let p_perp = tensor[2][2]
  box_z * (p_parallel - p_perp) / 2.0
}

///|
/// Compute von Mises stress from pressure tensor deviatoric part.
pub fn mda_von_mises_stress(tensor : Array[Array[Double]]) -> Double {
  let p = mda_scalar_pressure(tensor)
  // Deviatoric tensor
  let sxx = tensor[0][0] - p
  let syy = tensor[1][1] - p
  let szz = tensor[2][2] - p
  let sxy = tensor[0][1]
  let sxz = tensor[0][2]
  let syz = tensor[1][2]
  (0.5 *
  (
    (sxx - syy) * (sxx - syy) +
    (syy - szz) * (syy - szz) +
    (szz - sxx) * (szz - sxx)
  ) +
  3.0 * (sxy * sxy + sxz * sxz + syz * syz)).sqrt()
}
