/// omm_variable_integrator.mbt — variable time-step integrators.
///
/// OpenMM equivalent: `VariableVerletIntegrator`, `VariableLangevinIntegrator`.
///
/// Both integrators embed a local error estimate by comparing one full step
/// with two half steps (a second-order embedding) and adapt `dt` with the
/// standard PI-free controller dt' = 0.9 * dt * (tol/err)^(1/3), clamped to
/// [min_dt, max_dt]. The more accurate half-step result is adopted whenever
/// the error exceeds tolerance. The Langevin variant additionally rescales
/// velocities toward the target temperature (Berendsen-style coupling) so
/// that the thermostat works with a changing step size.
///
/// Units: kcal/mol, A, ps, K.

///|
/// Clamp a value into [lo, hi].
fn omm_vi_clamp(v : Double, lo : Double, hi : Double) -> Double {
  if v < lo {
    lo
  } else if v > hi {
    hi
  } else {
    v
  }
}

///|
/// One adaptive velocity-Verlet step.
///
/// Computes the embedded error estimate e = max_i |x_half2 - x_full| between
/// one full step and two half steps. If e > error_tol the half-step result is
/// adopted and dt shrinks; otherwise dt grows (capped by max_dt).
///
/// Returns (new_positions, new_velocities, new_dt).
pub fn omm_variable_verlet_step(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
  masses : Array[Double],
  dt~ : Double,
  error_tol~ : Double = 0.0001,
  min_dt~ : Double = 1.0e-7,
  max_dt~ : Double = 0.01,
) -> (Array[Array[Double]], Array[Array[Double]], Double) {
  let n = positions.length()
  let f0 = force_fn(positions)
  // --- full step (one Verlet of size dt) ---
  let xf : Array[Array[Double]] = md_new_forces(n)
  let vf : Array[Array[Double]] = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      xf[i][d] = positions[i][d] +
        velocities[i][d] * dt +
        0.5 * f0[i][d] * dt * dt / masses[i]
    }
  }
  let ff = force_fn(xf)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      vf[i][d] = velocities[i][d] + 0.5 * (f0[i][d] + ff[i][d]) * dt / masses[i]
    }
  }
  // --- two half steps (dt/2 each) ---
  let h = 0.5 * dt
  let xa : Array[Array[Double]] = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      xa[i][d] = positions[i][d] +
        velocities[i][d] * h +
        0.5 * f0[i][d] * h * h / masses[i]
    }
  }
  let fa = force_fn(xa)
  let va : Array[Array[Double]] = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      va[i][d] = velocities[i][d] + 0.5 * (f0[i][d] + fa[i][d]) * h / masses[i]
    }
  }
  let xb : Array[Array[Double]] = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      xb[i][d] = xa[i][d] + va[i][d] * h + 0.5 * fa[i][d] * h * h / masses[i]
    }
  }
  let fb = force_fn(xb)
  let vb : Array[Array[Double]] = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      vb[i][d] = va[i][d] + 0.5 * (fa[i][d] + fb[i][d]) * h / masses[i]
    }
  }
  // --- embedded error estimate + controller ---
  let mut err = 0.0
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      let e = (xf[i][d] - xb[i][d]).abs()
      if e > err {
        err = e
      }
    }
  }
  let mut dt_new = dt
  let use_half = err > error_tol
  if err > 0.0 {
    let mut factor = 0.9 * @math.pow(error_tol / err, 1.0 / 3.0)
    factor = omm_vi_clamp(factor, 0.2, 5.0)
    dt_new = omm_vi_clamp(dt * factor, min_dt, max_dt)
  } else {
    dt_new = omm_vi_clamp(dt * 1.5, min_dt, max_dt)
  }
  if use_half {
    (xb, vb, dt_new)
  } else {
    (xf, vf, dt_new)
  }
}

///|
/// Run an adaptive velocity-Verlet simulation.
/// Returns (positions_history, velocities_history, dt_history); the histories
/// contain n_steps + 1 entries (including the initial state).
pub fn omm_variable_verlet_simulate(
  initial_positions : Array[Array[Double]],
  initial_velocities : Array[Array[Double]],
  masses : Array[Double],
  n_steps~ : Int32,
  force_fn~ : (Array[Array[Double]]) -> Array[Array[Double]],
  dt0~ : Double = 0.001,
  error_tol~ : Double = 0.0001,
  min_dt~ : Double = 1.0e-7,
  max_dt~ : Double = 0.01,
) -> (Array[Array[Array[Double]]], Array[Array[Array[Double]]], Array[Double]) {
  let n = initial_positions.length()
  let pos_hist : Array[Array[Array[Double]]] = Array(capacity=0)
  let vel_hist : Array[Array[Array[Double]]] = Array(capacity=0)
  let dt_hist : Array[Double] = Array(capacity=0)
  let mut pos : Array[Array[Double]] = md_new_forces(n)
  let mut vel : Array[Array[Double]] = md_new_forces(n)
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
  pos_hist.push(pos)
  vel_hist.push(vel)
  dt_hist.push(dt0)
  let mut dt = dt0
  for _s = 0; _s < n_steps; _s = _s + 1 {
    let (x_new, v_new, dt_new) = omm_variable_verlet_step(
      pos,
      vel,
      force_fn,
      masses,
      dt~,
      error_tol~,
      min_dt~,
      max_dt~,
    )
    pos = x_new
    vel = v_new
    dt = dt_new
    pos_hist.push(pos)
    vel_hist.push(vel)
    dt_hist.push(dt)
  }
  (pos_hist, vel_hist, dt_hist)
}

///|
/// One adaptive Langevin step: an adaptive velocity-Verlet step followed by
/// a Berendsen-style velocity rescale toward `temperature` (coupling time
/// `tau`). The stochastic component is omitted from the error estimate.
/// Returns (new_positions, new_velocities, new_dt).
pub fn omm_variable_langevin_step(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
  masses : Array[Double],
  dt~ : Double,
  temperature~ : Double,
  tau~ : Double = 0.1,
  error_tol~ : Double = 0.0001,
  min_dt~ : Double = 1.0e-7,
  max_dt~ : Double = 0.01,
) -> (Array[Array[Double]], Array[Array[Double]], Double) {
  let (x_new, v_new, dt_new) = omm_variable_verlet_step(
    positions,
    velocities,
    force_fn,
    masses,
    dt~,
    error_tol~,
    min_dt~,
    max_dt~,
  )
  // Berendsen thermostat: lambda = sqrt(1 + (dt/tau) * (T0/T - 1))
  let n = x_new.length()
  let ke = omm_kinetic_energy(v_new, masses)
  let n_dof = if 3 * n - 3 > 0 { 3 * n - 3 } else { 1 }
  let temp = omm_temperature_from_energy(ke, n_dof)
  if temp > 1.0e-9 && tau > 0.0 {
    let lambda_sq = 1.0 + dt_new / tau * (temperature / temp - 1.0)
    if lambda_sq > 0.0 {
      let lambda = lambda_sq.sqrt()
      let scaled : Array[Array[Double]] = md_new_forces(n)
      for i = 0; i < n; i = i + 1 {
        for d = 0; d < 3; d = d + 1 {
          scaled[i][d] = v_new[i][d] * lambda
        }
      }
      return (x_new, scaled, dt_new)
    }
  }
  (x_new, v_new, dt_new)
}

///|
/// Run an adaptive Langevin simulation (see omm_variable_langevin_step).
/// Returns (positions_history, velocities_history, dt_history).
pub fn omm_variable_langevin_simulate(
  initial_positions : Array[Array[Double]],
  initial_velocities : Array[Array[Double]],
  masses : Array[Double],
  n_steps~ : Int32,
  force_fn~ : (Array[Array[Double]]) -> Array[Array[Double]],
  dt0~ : Double = 0.001,
  temperature~ : Double = 300.0,
  tau~ : Double = 0.1,
  error_tol~ : Double = 0.0001,
  min_dt~ : Double = 1.0e-7,
  max_dt~ : Double = 0.01,
) -> (Array[Array[Array[Double]]], Array[Array[Array[Double]]], Array[Double]) {
  let n = initial_positions.length()
  let pos_hist : Array[Array[Array[Double]]] = Array(capacity=0)
  let vel_hist : Array[Array[Array[Double]]] = Array(capacity=0)
  let dt_hist : Array[Double] = Array(capacity=0)
  let mut pos : Array[Array[Double]] = md_new_forces(n)
  let mut vel : Array[Array[Double]] = md_new_forces(n)
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
  pos_hist.push(pos)
  vel_hist.push(vel)
  dt_hist.push(dt0)
  let mut dt = dt0
  for _s = 0; _s < n_steps; _s = _s + 1 {
    let (x_new, v_new, dt_new) = omm_variable_langevin_step(
      pos,
      vel,
      force_fn,
      masses,
      dt~,
      temperature~,
      tau~,
      error_tol~,
      min_dt~,
      max_dt~,
    )
    pos = x_new
    vel = v_new
    dt = dt_new
    pos_hist.push(pos)
    vel_hist.push(vel)
    dt_hist.push(dt)
  }
  (pos_hist, vel_hist, dt_hist)
}
