/// SHAKE/RATTLE constraint algorithm.
///
/// OpenMM equivalent: `Context::applyConstraints` (SHAKE + RATTLE)
///
/// Constrains bond lengths to fixed values after an unconstrained
/// integration step. Essential for rigid water models (TIP3P/TIP4P).

///|
/// A constraint between atoms i and j with fixed distance
pub struct OmmConstraint {
  i : Int32
  j : Int32
  distance : Double // equilibrium distance
}

///|
pub fn OmmConstraint::i(self : OmmConstraint) -> Int32 {
  self.i
}

///|
pub fn OmmConstraint::j(self : OmmConstraint) -> Int32 {
  self.j
}

///|
pub fn OmmConstraint::distance(self : OmmConstraint) -> Double {
  self.distance
}

///|
/// Apply SHAKE constraints to positions.
///
/// positions: unconstrained new positions
/// old_positions: positions before the step
/// constraints: list of (atom_i, atom_j, target_distance)
/// masses: per-atom masses
/// tolerance: convergence tolerance (default 1e-6)
/// max_iter: maximum iterations (default 100)
///
/// Returns constrained positions.
pub fn omm_shake(
  positions : Array[Array[Double]],
  old_positions : Array[Array[Double]],
  constraints : Array[OmmConstraint],
  masses : Array[Double],
  tolerance~ : Double = 0.000001,
  max_iter~ : Int32 = 100,
) -> Array[Array[Double]] {
  let n = positions.length()
  // Copy positions
  let result : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    result[i] = [positions[i][0], positions[i][1], positions[i][2]]
  }
  let n_c = constraints.length()
  if n_c == 0 {
    return result
  }
  // Lagrange multipliers (one per constraint)
  let lagrange : Array[Double] = Array::make(n_c, 0.0)
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let mut max_error = 0.0
    for c = 0; c < n_c; c = c + 1 {
      let ci = constraints[c]
      let i = ci.i
      let j = ci.j
      let d0 = ci.distance
      // Current bond vector
      let dx = result[j][0] - result[i][0]
      let dy = result[j][1] - result[i][1]
      let dz = result[j][2] - result[i][2]
      let d_sq = dx * dx + dy * dy + dz * dz
      let d_curr = d_sq.sqrt()
      let error = d_sq - d0 * d0
      let abs_err = if error < 0.0 { -error } else { error }
      if abs_err > max_error {
        max_error = abs_err
      }
      // Reference bond vector (before step)
      let rx = old_positions[j][0] - old_positions[i][0]
      let ry = old_positions[j][1] - old_positions[i][1]
      let rz = old_positions[j][2] - old_positions[i][2]
      // SHAKE: solve for Lagrange multiplier
      let r_dot_d = rx * dx + ry * dy + rz * dz
      if r_dot_d == 0.0 {
        continue
      }
      let mi = masses[i]
      let mj = masses[j]
      let inv_mi = 1.0 / mi
      let inv_mj = 1.0 / mj
      let gamma = error / (2.0 * r_dot_d * (inv_mi + inv_mj))
      // Apply correction
      result[i][0] = result[i][0] + gamma * inv_mi * rx
      result[i][1] = result[i][1] + gamma * inv_mi * ry
      result[i][2] = result[i][2] + gamma * inv_mi * rz
      result[j][0] = result[j][0] - gamma * inv_mj * rx
      result[j][1] = result[j][1] - gamma * inv_mj * ry
      result[j][2] = result[j][2] - gamma * inv_mj * rz
      lagrange[c] = lagrange[c] + gamma
    }
    if max_error < tolerance {
      break
    }
  }
  result
}

///|
/// Apply RATTLE velocity constraints.
/// Ensures that constrained bond velocities are perpendicular to the bond.
///
/// velocities: unconstrained velocities
/// positions: constrained positions
/// constraints: list of constraints
/// masses: per-atom masses
pub fn omm_rattle(
  velocities : Array[Array[Double]],
  positions : Array[Array[Double]],
  constraints : Array[OmmConstraint],
  masses : Array[Double],
  tolerance~ : Double = 0.000001,
  max_iter~ : Int32 = 100,
) -> Array[Array[Double]] {
  let n = velocities.length()
  let result : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    result[i] = [velocities[i][0], velocities[i][1], velocities[i][2]]
  }
  for iter = 0; iter < max_iter; iter = iter + 1 {
    let mut max_error = 0.0
    for c = 0; c < constraints.length(); c = c + 1 {
      let ci = constraints[c]
      let i = ci.i
      let j = ci.j
      // Bond vector at constrained positions
      let rx = positions[j][0] - positions[i][0]
      let ry = positions[j][1] - positions[i][1]
      let rz = positions[j][2] - positions[i][2]
      // Relative velocity
      let dvx = result[j][0] - result[i][0]
      let dvy = result[j][1] - result[i][1]
      let dvz = result[j][2] - result[i][2]
      let v_dot_r = dvx * rx + dvy * ry + dvz * rz
      let abs_err = if v_dot_r < 0.0 { -v_dot_r } else { v_dot_r }
      if abs_err > max_error {
        max_error = abs_err
      }
      let r_sq = rx * rx + ry * ry + rz * rz
      if r_sq == 0.0 {
        continue
      }
      let mi = masses[i]
      let mj = masses[j]
      let inv_mi = 1.0 / mi
      let inv_mj = 1.0 / mj
      let kappa = v_dot_r / (r_sq * (inv_mi + inv_mj))
      result[i][0] = result[i][0] + kappa * inv_mi * rx
      result[i][1] = result[i][1] + kappa * inv_mi * ry
      result[i][2] = result[i][2] + kappa * inv_mi * rz
      result[j][0] = result[j][0] - kappa * inv_mj * rx
      result[j][1] = result[j][1] - kappa * inv_mj * ry
      result[j][2] = result[j][2] - kappa * inv_mj * rz
    }
    if max_error < tolerance {
      break
    }
  }
  result
}

///|
/// Build TIP3P water constraints (O-H1, O-H2, H1-H2)
pub fn omm_tip3p_constraints(
  o : Int32,
  h1 : Int32,
  h2 : Int32,
) -> Array[OmmConstraint] {
  [
    { i: o, j: h1, distance: 0.9572, },
    { i: o, j: h2, distance: 0.9572, },
    { i: h1, j: h2, distance: 1.5139006545, },
  ]
}

///|
/// Build TIP4P water constraints (O-H1, O-H2, H1-H2, same as TIP3P geometry)
pub fn omm_tip4p_constraints(
  o : Int32,
  h1 : Int32,
  h2 : Int32,
) -> Array[OmmConstraint] {
  omm_tip3p_constraints(o, h1, h2)
}
