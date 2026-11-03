/// Multiple time step (MTS) integrator: r-RESPA.
///
/// OpenMM equivalent: `MTSIntegrator` / `RespaIntegrator`
///
/// Separates forces into fast (bonded, short-range) and slow (long-range)
/// components, integrating them with different time steps for efficiency.

///|
/// r-RESPA (reversible Reference System Propagator Algorithm)
///
/// Positions and velocities are updated using a nested Verlet scheme:
/// - Outer step (dt): slow forces
/// - Inner steps (dt/n_inner): fast forces
///
/// positions: current positions
/// velocities: current velocities
/// fast_force_fn: fast forces (bonded, short-range nonbonded)
/// slow_force_fn: slow forces (long-range electrostatics, PME)
/// masses: per-atom masses
/// dt: outer time step
/// n_inner: number of inner steps
pub fn omm_respa_step(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  fast_force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
  slow_force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
  n_inner : Int32,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let n = positions.length()
  let dt_inner = dt / n_inner.to_double()
  let half_dt = dt * 0.5
  let half_dt_inner = dt_inner * 0.5
  // Initialize
  let pos : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  let vel : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    pos[i] = [positions[i][0], positions[i][1], positions[i][2]]
    vel[i] = [velocities[i][0], velocities[i][1], velocities[i][2]]
  }
  // Half step with slow forces
  let slow_f = slow_force_fn(pos)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      vel[i][d] = vel[i][d] + half_dt * slow_f[i][d] / masses[i]
    }
  }
  // Inner loop with fast forces (n_inner steps of velocity Verlet)
  for inner = 0; inner < n_inner; inner = inner + 1 {
    let fast_f = fast_force_fn(pos)
    // Half step velocity
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        vel[i][d] = vel[i][d] + half_dt_inner * fast_f[i][d] / masses[i]
      }
    }
    // Full step position
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        pos[i][d] = pos[i][d] + dt_inner * vel[i][d]
      }
    }
    // Half step velocity with new forces
    let fast_f_new = fast_force_fn(pos)
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        vel[i][d] = vel[i][d] + half_dt_inner * fast_f_new[i][d] / masses[i]
      }
    }
  }
  // Half step with slow forces at final positions
  let slow_f_final = slow_force_fn(pos)
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      vel[i][d] = vel[i][d] + half_dt * slow_f_final[i][d] / masses[i]
    }
  }
  (pos, vel)
}

///|
/// Run r-RESPA simulation for n_steps outer steps.
pub fn omm_respa_simulate(
  initial_positions : Array[Array[Double]],
  initial_velocities : Array[Array[Double]],
  fast_force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
  slow_force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
  n_steps : Int32,
  n_inner~ : Int32 = 2,
) -> Array[Array[Array[Double]]] {
  let n = initial_positions.length()
  let mut pos : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  let mut vel : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
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
  let history : Array[Array[Array[Double]]] = Array(capacity=0)
  history.push(pos)
  for step = 0; step < n_steps; step = step + 1 {
    let (new_pos, new_vel) = omm_respa_step(
      pos, vel, fast_force_fn, slow_force_fn, masses, dt, n_inner,
    )
    pos = new_pos
    vel = new_vel
    history.push(pos)
  }
  history
}
