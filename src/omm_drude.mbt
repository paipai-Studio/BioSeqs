/// omm_drude.mbt — Drude oscillator force field and integrator.
///
/// OpenMM equivalents: `DrudeForce` + `DrudeLangevinIntegrator`.
///
/// A Drude pair couples a massive core particle to a light, charged Drude
/// particle with a harmonic spring of constant `k_drude`
/// (kcal/mol/A^2). The induced dipole follows the core motion through a
/// strongly damped Langevin dynamics on the Drude particle, integrated
/// with `n_inner` sub-steps per main step.
///
/// The static polarizability of a Drude pair is
///
///   alpha = q_drude^2 * 332.06371 / k_drude   (A^3)
///
/// with q_drude in e and 332.06371 kcal*A/(mol*e^2) = 1/(4*pi*eps0).
///
/// Units: kcal/mol, A, ps, K, amu.

///|
/// One Drude oscillator: core particle index, Drude particle index, and
/// the harmonic spring constant (kcal/mol/A^2).
pub(all) struct OmmDrudePair {
  core : Int32
  particle : Int32
  k_drude : Double
}

///|
/// Static polarizability (A^3) of a Drude pair.
pub fn omm_drude_polarizability(q_drude : Double, k_drude : Double) -> Double {
  if k_drude <= 0.0 {
    return 0.0
  }
  q_drude * q_drude * 332.06371 / k_drude
}

///|
/// Spring force constant lookup helper: returns the harmonic force on
/// (core, particle) for one pair: F_core = -k r, F_drude = +k r with
/// r = x_drude - x_core.
pub fn omm_drude_spring_forces(
  positions : Array[Array[Double]],
  pairs : Array[OmmDrudePair],
) -> Array[Array[Double]] {
  let n = positions.length()
  let out = md_new_forces(n)
  for p in pairs {
    if p.core >= n || p.particle >= n {
      continue
    }
    let c = positions[p.core]
    let d = positions[p.particle]
    let rx = d[0] - c[0]
    let ry = d[1] - c[1]
    let rz = d[2] - c[2]
    // U = 0.5 k |x_d - x_c|^2: the spring pulls each particle toward the
    // other one (F_core = +k r, F_drude = -k r with r = x_d - x_c).
    out[p.core][0] = out[p.core][0] + p.k_drude * rx
    out[p.core][1] = out[p.core][1] + p.k_drude * ry
    out[p.core][2] = out[p.core][2] + p.k_drude * rz
    out[p.particle][0] = out[p.particle][0] - p.k_drude * rx
    out[p.particle][1] = out[p.particle][1] - p.k_drude * ry
    out[p.particle][2] = out[p.particle][2] - p.k_drude * rz
  }
  out
}

///|
/// Drude Langevin step.
///
/// Each main step `dt` is split into `n_inner` sub-steps. All particles
/// use velocity Verlet on (external + spring) forces; Drude particles are
/// additionally coupled to the bath with `friction_drude` (Ornstein-
/// Uhlenbeck, deterministic given the seeded RNG), cores get the weak
/// `friction_core` coupling after the last sub-step.
///
/// Returns (new_positions, new_velocities).
pub fn omm_drude_langevin_step(
  positions : Array[Array[Double]],
  velocities : Array[Array[Double]],
  force_fn~ : (Array[Array[Double]]) -> Array[Array[Double]],
  masses~ : Array[Double],
  pairs~ : Array[OmmDrudePair],
  dt~ : Double,
  temperature~ : Double = 300.0,
  friction_core~ : Double = 1.0,
  friction_drude~ : Double = 50.0,
  n_inner~ : Int32 = 4,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let n = positions.length()
  let pos : Array[Array[Double]] = md_new_forces(n)
  let vel : Array[Array[Double]] = md_new_forces(n)
  for i = 0; i < n; i = i + 1 {
    pos[i] = [positions[i][0], positions[i][1], positions[i][2]]
    vel[i] = [velocities[i][0], velocities[i][1], velocities[i][2]]
  }
  let h = dt / n_inner.to_double().max(1.0)
  let kT = MD_KB * temperature
  let half_h = 0.5 * h
  for _sub = 0; _sub < n_inner; _sub = _sub + 1 {
    // Plain velocity Verlet on (external + spring) forces, with the
    // forces re-evaluated at the updated positions.
    let f_ext = force_fn(pos)
    let f_spring = omm_drude_spring_forces(pos, pairs)
    let v_half : Array[Array[Double]] = md_new_forces(n)
    for i = 0; i < n; i = i + 1 {
      let m = if i < masses.length() { masses[i] } else { 1.0 }
      for d = 0; d < 3; d = d + 1 {
        let f = f_ext[i][d] + f_spring[i][d]
        v_half[i][d] = vel[i][d] + half_h * f / m
        pos[i][d] = pos[i][d] + h * v_half[i][d]
      }
    }
    let f_ext2 = force_fn(pos)
    let f_spring2 = omm_drude_spring_forces(pos, pairs)
    for i = 0; i < n; i = i + 1 {
      let m = if i < masses.length() { masses[i] } else { 1.0 }
      for d = 0; d < 3; d = d + 1 {
        vel[i][d] = v_half[i][d] + half_h * (f_ext2[i][d] + f_spring2[i][d]) / m
      }
    }
  }
  // Ornstein-Uhlenbeck coupling at the end of the main step. For each
  // Drude pair the thermostat acts on the RELATIVE velocity (reduced
  // mass mu = mc*md/(mc+md)), so the random force carries no net center-
  // of-mass momentum; unpaired cores get the weak friction_core damping.
  let c_drude = @math.exp(-friction_drude * dt)
  let c_core = @math.exp(-friction_core * dt)
  let thermostatted : Array[Bool] = Array::make(n, false)
  for p in pairs {
    if p.core < n && p.particle < n {
      thermostatted[p.core] = true
      thermostatted[p.particle] = true
      let mc = if p.core < masses.length() { masses[p.core] } else { 1.0 }
      let md = if p.particle < masses.length() {
        masses[p.particle]
      } else {
        1.0
      }
      let mt = mc + md
      let mu = mc * md / mt
      let s = @math.pow((1.0 - c_drude * c_drude) * kT / mu, 0.5)
      for d = 0; d < 3; d = d + 1 {
        let vcm = (mc * vel[p.core][d] + md * vel[p.particle][d]) / mt
        let vrel = vel[p.particle][d] - vel[p.core][d]
        let vrel_new = c_drude * vrel + s * md_gaussian()
        vel[p.core][d] = vcm - md / mt * vrel_new
        vel[p.particle][d] = vcm + mc / mt * vrel_new
      }
    }
  }
  for i = 0; i < n; i = i + 1 {
    if not(thermostatted[i]) {
      for d = 0; d < 3; d = d + 1 {
        vel[i][d] = vel[i][d] * c_core
      }
    }
  }
  (pos, vel)
}

///|
/// Run a Drude Langevin simulation for n_steps main steps.
/// Returns (positions_history, velocities_history) with n_steps + 1
/// entries.
pub fn omm_drude_langevin_simulate(
  initial_positions : Array[Array[Double]],
  initial_velocities : Array[Array[Double]],
  force_fn~ : (Array[Array[Double]]) -> Array[Array[Double]],
  masses~ : Array[Double],
  pairs~ : Array[OmmDrudePair],
  n_steps~ : Int32,
  dt~ : Double = 0.0005,
  temperature~ : Double = 300.0,
  friction_core~ : Double = 1.0,
  friction_drude~ : Double = 50.0,
  n_inner~ : Int32 = 4,
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
  for _s = 0; _s < n_steps; _s = _s + 1 {
    let (p2, v2) = omm_drude_langevin_step(
      pos,
      vel,
      force_fn~,
      masses~,
      pairs~,
      dt~,
      temperature~,
      friction_core~,
      friction_drude~,
      n_inner~,
    )
    for i = 0; i < n; i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        pos[i][d] = p2[i][d]
        vel[i][d] = v2[i][d]
      }
    }
    snapshot()
  }
  (pos_hist, vel_hist)
}
