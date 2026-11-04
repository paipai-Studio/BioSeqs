/// Reaction field electrostatics.
///
/// OpenMM equivalent: `NonbondedForce::CutoffPeriodic` with reaction field
///
/// Alternative to PME for cutoff electrostatics in periodic systems.
/// Beyond the cutoff, the dielectric continuum approximation is used.
/// Faster than PME but less accurate for highly charged systems.

///|
/// Reaction field energy between two charges.
///
/// U_RF(r) = q_i*q_j * [1/r + k_rf*r^2 - c_rf]  for r < cutoff
/// where:
///   k_rf = (1/cutoff^3) * (epsilon_rf - 1) / (2*epsilon_rf + 1)
///   c_rf = (1/cutoff) * (3*epsilon_rf) / (2*epsilon_rf + 1)
///
/// epsilon_rf: reaction field dielectric (default 78.5 for water)
/// cutoff: nonbonded cutoff distance
pub fn omm_reaction_field_energy(
  r : Double,
  qi : Double,
  qj : Double,
  cutoff : Double,
  epsilon_rf~ : Double = 78.5,
) -> Double {
  if r >= cutoff || r == 0.0 {
    return 0.0
  }
  let k_rf = 1.0 /
    (cutoff * cutoff * cutoff) *
    (epsilon_rf - 1.0) /
    (2.0 * epsilon_rf + 1.0)
  let c_rf = 1.0 / cutoff * (3.0 * epsilon_rf) / (2.0 * epsilon_rf + 1.0)
  qi * qj * (1.0 / r + k_rf * r * r - c_rf)
}

///|
/// Reaction field force magnitude (negative gradient of energy).
pub fn omm_reaction_field_force(
  r : Double,
  qi : Double,
  qj : Double,
  cutoff : Double,
  epsilon_rf~ : Double = 78.5,
) -> Double {
  if r >= cutoff || r == 0.0 {
    return 0.0
  }
  let k_rf = 1.0 /
    (cutoff * cutoff * cutoff) *
    (epsilon_rf - 1.0) /
    (2.0 * epsilon_rf + 1.0)
  // dU/dr = q_i*q_j * (-1/r^2 + 2*k_rf*r)
  // F = -dU/dr = q_i*q_j * (1/r^2 - 2*k_rf*r)
  qi * qj * (1.0 / (r * r) - 2.0 * k_rf * r)
}

///|
/// Total reaction field energy for a system.
pub fn omm_reaction_field_total(
  positions : Array[Array[Double]],
  charges : Array[Double],
  cutoff : Double,
  epsilon_rf~ : Double = 78.5,
  box~ : Array[Double] = [0.0, 0.0, 0.0],
) -> Double {
  let n = positions.length()
  let mut energy = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let dx = positions[j][0] - positions[i][0]
      let dy = positions[j][1] - positions[i][1]
      let dz = positions[j][2] - positions[i][2]
      // PBC minimum image
      let mut rx = dx
      let mut ry = dy
      let mut rz = dz
      if box[0] > 0.0 {
        rx = rx - box[0] * (rx / box[0] + 0.5).floor()
        ry = ry - box[1] * (ry / box[1] + 0.5).floor()
        rz = rz - box[2] * (rz / box[2] + 0.5).floor()
      }
      let r = (rx * rx + ry * ry + rz * rz).sqrt()
      energy = energy +
        omm_reaction_field_energy(
          r,
          charges[i],
          charges[j],
          cutoff,
          epsilon_rf~,
        )
    }
  }
  energy
}

///|
/// Reaction field with shifted force (smooth at cutoff).
pub fn omm_reaction_field_shifted_energy(
  r : Double,
  qi : Double,
  qj : Double,
  cutoff : Double,
  epsilon_rf~ : Double = 78.5,
) -> Double {
  if r >= cutoff || r == 0.0 {
    return 0.0
  }
  let k_rf = 1.0 /
    (cutoff * cutoff * cutoff) *
    (epsilon_rf - 1.0) /
    (2.0 * epsilon_rf + 1.0)
  let c_rf = 1.0 / cutoff * (3.0 * epsilon_rf) / (2.0 * epsilon_rf + 1.0)
  // Shifted: U(r) - U(cutoff) to make energy continuous at cutoff
  let u_r = qi * qj * (1.0 / r + k_rf * r * r - c_rf)
  let u_cut = qi * qj * (1.0 / cutoff + k_rf * cutoff * cutoff - c_rf)
  u_r - u_cut
}

///|
/// Compare reaction field vs PME accuracy for a pair.
/// Returns the absolute difference.
pub fn omm_rf_vs_pme_error(
  r : Double,
  qi : Double,
  qj : Double,
  cutoff : Double,
  epsilon_rf~ : Double = 78.5,
) -> Double {
  let e_rf = omm_reaction_field_energy(r, qi, qj, cutoff, epsilon_rf~)
  // PME real-space contribution (approximate: erfc(alpha*r)/r)
  let alpha = 3.141592653589793.sqrt() / cutoff
  let arg = alpha * r
  // erfc approximation
  let t = 1.0 / (1.0 + 0.3275911 * arg)
  let erfc_val = 1.0 -
    t *
    (
      0.254829592 +
      t *
      (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))
    )
  let e_pme = qi * qj * erfc_val / r
  let diff = e_rf - e_pme
  if diff < 0.0 {
    -diff
  } else {
    diff
  }
}
