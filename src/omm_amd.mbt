/// omm_amd.mbt — accelerated molecular dynamics.
///
/// OpenMM equivalent: `CustomIntegrator`-based aMD recipes.
///
/// Implements the standard aMD boost (Hamelberg et al. 2004) and Gaussian
/// aMD (Miao et al. 2015). When the potential V(x) is below the threshold
/// E, it is boosted:
///
///   aMD :  dV = (E - V)^2 / (alpha + E - V)
///   GaMD:  dV = (E - V)^2 / (alpha + (E - V)^2)
///
/// and the effective force is scaled by d(V + dV)/dV, keeping the dynamics
/// on the modified potential. For V >= E the potential is unchanged.
///
/// Units: kcal/mol (energy), A (length), ps (time).

///|
/// aMD / GaMD boost parameters: threshold E and smoothing alpha.
pub struct OmmAMDParams {
  e_threshold : Double
  alpha : Double
}

///|
pub fn OmmAMDParams::new(
  e_threshold~ : Double,
  alpha~ : Double = 1.0,
) -> OmmAMDParams {
  { e_threshold, alpha, }
}

///|
/// Standard aMD boost. Returns (boosted potential, d(V')/dV).
pub fn omm_amd_boost(params : OmmAMDParams, v : Double) -> (Double, Double) {
  if v >= params.e_threshold {
    return (v, 1.0)
  }
  let x = params.e_threshold - v
  let denom = params.alpha + x
  let d_v = x * x / denom
  // d(dV)/dV = -x(2*alpha + x)/(alpha + x)^2
  let ratio = 1.0 - x * (2.0 * params.alpha + x) / (denom * denom)
  (v + d_v, ratio)
}

///|
/// Gaussian aMD boost. Returns (boosted potential, d(V')/dV).
pub fn omm_gamd_boost(params : OmmAMDParams, v : Double) -> (Double, Double) {
  if v >= params.e_threshold {
    return (v, 1.0)
  }
  let x = params.e_threshold - v
  let denom = params.alpha + x * x
  let d_v = x * x / denom
  // d(dV)/dV = -2*alpha*x / (alpha + x^2)^2
  let ratio = 1.0 - 2.0 * params.alpha * x / (denom * denom)
  (v + d_v, ratio)
}

///|
/// Run an aMD simulation (velocity Verlet on the boosted potential).
///
/// The boost scales the forces by dV'/dV evaluated at the current
/// configuration, so the trajectory samples the modified potential while
/// plain kinetic energy is conserved on it.
///
/// Returns (positions_history, velocities_history) with n_steps + 1
/// entries (including the initial state).
pub fn omm_amd_simulate(
  initial_positions : Array[Array[Double]],
  initial_velocities : Array[Array[Double]],
  masses : Array[Double],
  n_steps~ : Int32,
  force_fn~ : (Array[Array[Double]]) -> Array[Array[Double]],
  energy_fn~ : (Array[Array[Double]]) -> Double,
  params : OmmAMDParams,
  dt~ : Double = 0.001,
  gaussian~ : Bool = false,
) -> (Array[Array[Array[Double]]], Array[Array[Array[Double]]]) {
  let n = initial_positions.length()
  let pos : Array[Array[Double]] = md_new_forces(n)
  let vel : Array[Array[Double]] = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    pos[i] = [
      initial_positions[i][0],
      initial_positions[i][1],
      initial_positions[i][2],
    ]
    vel[i] = [
      initial_velocities[i][0],
      initial_velocities[i][1],
      initial_velocities[i][2],
    ]
  }
  let pos_hist : Array[Array[Array[Double]]] = Array(capacity=0)
  let vel_hist : Array[Array[Array[Double]]] = Array(capacity=0)
  fn snapshot() -> Unit {
    let pr : Array[Array[Double]] = md_new_forces(n)
    let vr : Array[Array[Double]] = md_new_forces(n)
    for i = 0; i < n; i = i + 1 {
      pr[i] = [pos[i][0], pos[i][1], pos[i][2]]
      vr[i] = [vel[i][0], vel[i][1], vel[i][2]]
    }
    pos_hist.push(pr)
    vel_hist.push(vr)
  }

  snapshot()
  let half_dt = 0.5 * dt
  for _s = 0; _s < n_steps; _s = _s + 1 {
    // Effective (boosted) forces at the current configuration.
    let ratio0 = if gaussian {
      omm_gamd_boost(params, energy_fn(pos)).1
    } else {
      omm_amd_boost(params, energy_fn(pos)).1
    }
    let f0 = force_fn(pos)
    // Half kick + half drift.
    let v_half : Array[Array[Double]] = md_new_forces(n)
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        v_half[i][d] = vel[i][d] + half_dt * ratio0 * f0[i][d] / masses[i]
        pos[i][d] = pos[i][d] + dt * v_half[i][d]
      }
    }
    // New boosted forces.
    let ratio1 = if gaussian {
      omm_gamd_boost(params, energy_fn(pos)).1
    } else {
      omm_amd_boost(params, energy_fn(pos)).1
    }
    let f1 = force_fn(pos)
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        vel[i][d] = v_half[i][d] + half_dt * ratio1 * f1[i][d] / masses[i]
      }
    }
    snapshot()
  }
  (pos_hist, vel_hist)
}
