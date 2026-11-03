/// Leapfrog Verlet integrator.
///
/// OpenMM equivalent: `LeapfrogVerletIntegrator`
///
/// Stores velocities at half-integer time steps (t-dt/2),
/// positions at integer steps (t). Different storage scheme from
/// velocity Verlet, though mathematically equivalent.

///|
/// Leapfrog: x(t+dt) = x(t) + dt * v(t+dt/2)
///           v(t+dt/2) = v(t-dt/2) + dt * F(t)/m
pub fn omm_leapfrog_step(
  positions : Array[Array[Double]],
  velocities_half : Array[Array[Double]],
  forces : Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let n = positions.length()
  let x_new : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  let v_new : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      // v(t+dt/2) = v(t-dt/2) + dt * F/m
      v_new[i][d] = velocities_half[i][d] + dt * forces[i][d] / masses[i]
      // x(t+dt) = x(t) + dt * v(t+dt/2)
      x_new[i][d] = positions[i][d] + dt * v_new[i][d]
    }
  }
  (x_new, v_new)
}

///|
/// Convert leapfrog half-step velocities to full-step velocities
/// for output/analysis: v(t) ≈ v(t-dt/2) + 0.5*dt*F(t)/m
pub fn omm_leapfrog_to_full_velocity(
  velocities_half : Array[Array[Double]],
  forces : Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
) -> Array[Array[Double]] {
  let n = velocities_half.length()
  let v_full : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      v_full[i][d] = velocities_half[i][d] + 0.5 * dt * forces[i][d] / masses[i]
    }
  }
  v_full
}

///|
/// Run leapfrog Verlet simulation.
/// force_fn: function that computes forces from positions
/// Returns (positions_history, velocities_half_history)
pub fn omm_leapfrog_simulate(
  initial_positions : Array[Array[Double]],
  initial_velocities : Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
  n_steps : Int32,
  force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
) -> (Array[Array[Array[Double]]], Array[Array[Array[Double]]]) {
  let n = initial_positions.length()
  let pos_history : Array[Array[Array[Double]]] = Array(capacity=0)
  let vel_history : Array[Array[Array[Double]]] = Array(capacity=0)
  let mut pos = initial_positions
  // Initial half-step velocity: v(dt/2) = v(0) + 0.5*dt*F(0)/m
  let f0 = force_fn(pos)
  let mut vel_half : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      vel_half[i][d] = initial_velocities[i][d] +
        0.5 * dt * f0[i][d] / masses[i]
    }
  }
  pos_history.push(pos)
  vel_history.push(vel_half)
  for step = 0; step < n_steps; step = step + 1 {
    let (x_new, v_new) = omm_leapfrog_step(
      pos,
      vel_half,
      force_fn(pos),
      masses,
      dt,
    )
    pos = x_new
    vel_half = v_new
    pos_history.push(pos)
    vel_history.push(vel_half)
  }
  (pos_history, vel_history)
}
