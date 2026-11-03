/// Diffusion coefficient from MSD, and velocity autocorrelation (VACF).
///
/// MDAnalysis equivalent: `MSDAnalysis.diffusion_coefficient`,
/// `velocity_autocorrelation` (Green-Kubo relation)

///|
/// Compute diffusion coefficient from MSD using Einstein relation.
/// D = MSD / (2 * d * t) where d = dimensionality (3 for 3D).
///
/// msd_values: MSD at each time lag
/// dt: time step between frames (ps)
/// dimension: 3 for 3D, 2 for 2D, 1 for 1D
pub fn mda_diffusion_coefficient(
  msd_values : Array[Double],
  dt~ : Double = 1.0,
  dimension~ : Int32 = 3,
) -> Double {
  let n = msd_values.length()
  if n < 2 {
    return 0.0
  }
  // Linear regression: MSD = 2*d*D*t + intercept
  // Use the second half of the trajectory for better statistics
  let start = n / 4
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_xy = 0.0
  let mut sum_x2 = 0.0
  let mut count = 0.0
  for i = start; i < n; i = i + 1 {
    let t = i.to_double() * dt
    let y = msd_values[i]
    sum_x = sum_x + t
    sum_y = sum_y + y
    sum_xy = sum_xy + t * y
    sum_x2 = sum_x2 + t * t
    count = count + 1.0
  }
  if count < 2.0 {
    return 0.0
  }
  let denom = count * sum_x2 - sum_x * sum_x
  if denom == 0.0 {
    return 0.0
  }
  let slope = (count * sum_xy - sum_x * sum_y) / denom
  // D = slope / (2 * d)
  slope / (2.0 * dimension.to_double())
}

///|
/// Compute diffusion coefficient directly from trajectory.
/// traj: frames[atom][xyz] coordinates
pub fn mda_diffusion_from_trajectory(
  traj : Array[Array[Array[Double]]],
  dt~ : Double = 1.0,
  dimension~ : Int32 = 3,
) -> Double {
  let msd = mdanalysis_msd(traj)
  mda_diffusion_coefficient(msd, dt~, dimension~)
}

///|
/// Velocity autocorrelation function (VACF).
/// velocities: velocity at each frame, per atom, per xyz
/// Returns C(t) = <v(0)·v(t)> / <v(0)·v(0)>
pub fn mda_velocity_autocorrelation(
  velocities : Array[Array[Array[Double]]],
) -> Array[Double] {
  let n_frames = velocities.length()
  if n_frames == 0 {
    return []
  }
  let n_atoms = velocities[0].length()
  let result : Array[Double] = Array::make(n_frames, 0.0)
  // Normalize by VACF at t=0
  let mut vacf0 = 0.0
  for a = 0; a < n_atoms; a = a + 1 {
    let vx = velocities[0][a][0]
    let vy = velocities[0][a][1]
    let vz = velocities[0][a][2]
    vacf0 = vacf0 + vx * vx + vy * vy + vz * vz
  }
  if vacf0 == 0.0 {
    return result
  }
  for lag = 0; lag < n_frames; lag = lag + 1 {
    let mut sum = 0.0
    let mut count = 0.0
    let t_max = n_frames - lag
    for t = 0; t < t_max; t = t + 1 {
      for a = 0; a < n_atoms; a = a + 1 {
        let v0x = velocities[t][a][0]
        let v0y = velocities[t][a][1]
        let v0z = velocities[t][a][2]
        let vtx = velocities[t + lag][a][0]
        let vty = velocities[t + lag][a][1]
        let vtz = velocities[t + lag][a][2]
        sum = sum + v0x * vtx + v0y * vty + v0z * vtz
      }
      count = count + 1.0
    }
    if count > 0.0 {
      result[lag] = sum / count / vacf0 * n_atoms.to_double()
    }
  }
  result
}

///|
/// Diffusion coefficient via Green-Kubo relation.
/// D = (1/3) * integral(VACF(t) dt)
///
/// vacf: normalized VACF values
/// dt: time step
pub fn mda_green_kubo_diffusion(
  vacf : Array[Double],
  dt~ : Double = 1.0,
) -> Double {
  let n = vacf.length()
  if n < 2 {
    return 0.0
  }
  // Trapezoidal integration
  let mut integral = 0.0
  for i = 0; i < n - 1; i = i + 1 {
    integral = integral + (vacf[i] + vacf[i + 1]) * 0.5 * dt
  }
  integral / 3.0
}

///|
/// Stokes-Einstein radius from diffusion coefficient.
/// r = kT / (6 * pi * eta * D)
pub fn mda_stokes_einstein_radius(
  diffusion : Double,
  temperature~ : Double = 300.0,
  viscosity~ : Double = 0.001,
) -> Double {
  if diffusion <= 0.0 {
    return 0.0
  }
  let kT = 1.380649e-23 * temperature
  kT / (6.0 * 3.141592653589793 * viscosity * diffusion)
}
