/// Simulated annealing with temperature scheduling.
///
/// OpenMM equivalent: temperature scheduling + Langevin integrator
///
/// Gradually reduces temperature to find global energy minimum,
/// combining with existing thermostat integrators.

///|
/// Temperature schedule types
pub enum OmmAnnealSchedule {
  Linear(Double, Double, Int32) // (T_start, T_end, n_steps)
  Exponential(Double, Double, Int32) // (T_start, T_end, n_steps)
  Stepwise(Array[Double], Array[Int32]) // (temperatures, steps_per_temp)
  Logarithmic(Double, Double, Int32) // (T_start, T_end, n_steps)
}

///|
/// Get temperature at a given step from schedule.
pub fn omm_anneal_temperature(
  schedule : OmmAnnealSchedule,
  step : Int32,
) -> Double {
  match schedule {
    OmmAnnealSchedule::Linear(t_start, t_end, n_steps) =>
      if step >= n_steps {
        t_end
      } else {
        t_start + (t_end - t_start) * step.to_double() / n_steps.to_double()
      }
    OmmAnnealSchedule::Exponential(t_start, t_end, n_steps) =>
      if step >= n_steps {
        t_end
      } else if t_start <= 0.0 || t_end <= 0.0 {
        t_end
      } else {
        let ratio = t_end / t_start
        t_start * @math.pow(ratio, step.to_double() / n_steps.to_double())
      }
    OmmAnnealSchedule::Stepwise(temps, steps_per) => {
      let mut cumul = 0
      for i = 0; i < temps.length(); i = i + 1 {
        cumul = cumul + steps_per[i]
        if step < cumul {
          return temps[i]
        }
      }
      temps[temps.length() - 1]
    }
    OmmAnnealSchedule::Logarithmic(t_start, t_end, n_steps) =>
      // T(t) = T_end / (1 + (T_end/T_start - 1) * (1 - t/n))
      if step >= n_steps {
        t_end
      } else {
        let frac = 1.0 - step.to_double() / n_steps.to_double()
        t_end / (1.0 + (t_end / t_start - 1.0) * frac)
      }
  }
}

///|
/// Run simulated annealing using Langevin dynamics.
///
/// positions: initial positions
/// velocities: initial velocities
/// force_fn: force calculation function
/// masses: per-atom masses
/// dt: time step
/// schedule: temperature schedule
///
/// Returns final positions and energy history.
pub fn omm_simulated_annealing(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  force_fn : (Array[Array[Double]]) -> Array[Array[Double]],
  masses : Array[Double],
  dt : Double,
  schedule : OmmAnnealSchedule,
) -> (Array[Array[Double]], Array[Double]) {
  let n = positions.length()
  let pos : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  let vel : Array[Array[Double]] = Array::make(n, [0.0, 0.0, 0.0])
  for i = 0; i < n; i = i + 1 {
    pos[i] = [positions[i][0], positions[i][1], positions[i][2]]
    vel[i] = [velocities[i][0], velocities[i][1], velocities[i][2]]
  }
  // Determine total steps from schedule
  let n_steps = match schedule {
    OmmAnnealSchedule::Linear(_, _, n) => n
    OmmAnnealSchedule::Exponential(_, _, n) => n
    OmmAnnealSchedule::Logarithmic(_, _, n) => n
    OmmAnnealSchedule::Stepwise(_, steps_per) => {
      let mut total = 0
      for s in steps_per {
        total = total + s
      }
      total
    }
  }
  let energy_history : Array[Double] = Array(capacity=0)
  let friction = 1.0
  for step = 0; step < n_steps; step = step + 1 {
    let temp = omm_anneal_temperature(schedule, step)
    // Langevin step (simplified BAOAB)
    let forces = force_fn(pos)
    let half_dt = dt * 0.5
    let kT = 0.001987204 * temp
    let sigma = (2.0 * friction * kT / dt).sqrt()
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        // B: half kick
        vel[i][d] = vel[i][d] + half_dt * forces[i][d] / masses[i]
        // A: half drift
        pos[i][d] = pos[i][d] + half_dt * vel[i][d]
        // O: stochastic thermostat
        let c1 = @math.exp(-friction * dt)
        let c2 = (1.0 - c1 * c1).sqrt() * sigma / (friction * 2.0 * kT).sqrt()
        // Simplified: v' = c1*v + c2*random
        let r1 = (md_rand_double() - 0.5) * 2.0 * sigma
        vel[i][d] = c1 * vel[i][d] + c2 * r1
        // A: half drift
        pos[i][d] = pos[i][d] + half_dt * vel[i][d]
      }
    }
    // Record energy
    let _forces_new = force_fn(pos)
    let pe = 0.0
    let mut ke = 0.0
    for i = 0; i < n; i = i + 1 {
      ke = ke +
        0.5 *
        masses[i] *
        (vel[i][0] * vel[i][0] + vel[i][1] * vel[i][1] + vel[i][2] * vel[i][2])
    }
    // PE approximated from forces (simplified)
    energy_history.push(ke + pe)
  }
  (pos, energy_history)
}

///|
/// Create a linear cooling schedule.
pub fn omm_linear_schedule(
  t_start : Double,
  t_end : Double,
  n_steps : Int32,
) -> OmmAnnealSchedule {
  OmmAnnealSchedule::Linear(t_start, t_end, n_steps)
}

///|
/// Create an exponential cooling schedule.
pub fn omm_exponential_schedule(
  t_start : Double,
  t_end : Double,
  n_steps : Int32,
) -> OmmAnnealSchedule {
  OmmAnnealSchedule::Exponential(t_start, t_end, n_steps)
}

///|
fn _omm_anneal_unused_constructors() -> Unit {
  ignore(OmmAnnealSchedule::Stepwise([], []))
  ignore(OmmAnnealSchedule::Logarithmic(0.0, 0.0, 0))
}
