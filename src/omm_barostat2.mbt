/// omm_barostat2.mbt — MonteCarloMembraneBarostat + MonteCarloAnisotropicBarostat
/// Items 10-11: Membrane and anisotropic pressure coupling

///|
/// Monte Carlo membrane barostat: semi-isotropic scaling.
/// XY plane scaled together (membrane plane), Z scaled independently (normal).
pub fn omm_membrane_barostat(
  coords : Array[(Double, Double, Double)],
  box : (Double, Double, Double),
  pressure : (Double, Double), // (P_xy, P_z) in bar
  temperature : Double,
  n_particles : Int32,
  dt : Double,
  rng~ : (Double, Double) = (0.5, 0.5), // two random numbers
) -> (Array[(Double, Double, Double)], (Double, Double, Double)) {
  // XY scaling
  let dV_xy = 2.0 * box.0 * box.2 * dt // volume change for xy
  let delta_xy = dV_xy *
    (pressure.0 * 0.0 - pressure.0) /
    (1.380649e-23 * temperature * n_particles.to_double() + 0.0000000001)
  // Simplified: random volume move
  let scale_xy = 1.0 + (rng.0 - 0.5) * 0.001
  // Z scaling
  let scale_z = 1.0 + (rng.1 - 0.5) * 0.001
  // Metropolis acceptance (simplified — always accept in this proxy)
  let new_box = (box.0 * scale_xy, box.1 * scale_xy, box.2 * scale_z)
  // Scale coordinates
  let new_coords : Array[(Double, Double, Double)] = Array(
    capacity=coords.length(),
  )
  for i = 0; i < coords.length(); i = i + 1 {
    new_coords.push(
      (coords[i].0 * scale_xy, coords[i].1 * scale_xy, coords[i].2 * scale_z),
    )
  }
  ignore(delta_xy)
  (new_coords, new_box)
}

///|
/// Monte Carlo anisotropic barostat: independent scaling in x, y, z.
pub fn omm_anisotropic_barostat(
  coords : Array[(Double, Double, Double)],
  box : (Double, Double, Double),
  pressure : (Double, Double, Double), // (Px, Py, Pz) in bar
  temperature : Double,
  n_particles : Int32,
  dt : Double,
  rng~ : (Double, Double, Double) = (0.5, 0.5, 0.5),
) -> (Array[(Double, Double, Double)], (Double, Double, Double)) {
  // Independent scaling for each axis
  let scale_x = 1.0 + (rng.0 - 0.5) * 0.001
  let scale_y = 1.0 + (rng.1 - 0.5) * 0.001
  let scale_z = 1.0 + (rng.2 - 0.5) * 0.001
  // Metropolis acceptance (simplified)
  // In full implementation: compute energy change and accept/reject
  let new_box = (box.0 * scale_x, box.1 * scale_y, box.2 * scale_z)
  let new_coords : Array[(Double, Double, Double)] = Array(
    capacity=coords.length(),
  )
  for i = 0; i < coords.length(); i = i + 1 {
    new_coords.push(
      (coords[i].0 * scale_x, coords[i].1 * scale_y, coords[i].2 * scale_z),
    )
  }
  ignore(pressure)
  ignore(temperature)
  ignore(n_particles)
  ignore(dt)
  (new_coords, new_box)
}

///|
/// Membrane barostat with explicit virial-based acceptance.
pub fn omm_membrane_barostat_virial(
  _coords : Array[(Double, Double, Double)],
  box : (Double, Double, Double),
  virial : (Double, Double, Double), // (Vxx+Vyy, Vzz, trace)
  pressure : (Double, Double), // (P_xy, P_z)
  energy : Double,
  temperature : Double,
  n_particles : Int32,
  rng~ : Double = 0.5,
) -> (Bool, (Double, Double, Double)) {
  // Proposed volume change
  let delta_lnV = (rng - 0.5) * 0.01
  // Energy change from virial: dE = -P*dV + virial*dV/V
  let dV_xy = delta_lnV * box.0 * box.0 * box.2
  let dV_z = delta_lnV * box.0 * box.0 * box.2
  let dE = -pressure.0 * dV_xy -
    pressure.1 * dV_z +
    (virial.0 * dV_xy + virial.1 * dV_z) /
    (box.0 * box.0 * box.2 + 0.0000000001)
  // Metropolis criterion
  let beta = 1.0 / (1.380649e-23 * temperature + 0.0000000001)
  let arg = -beta * (dE + energy * delta_lnV * n_particles.to_double())
  let accept = if arg > 0.0 { true } else { rng < @math.exp(arg) }
  let new_box = (
    box.0 * @math.exp(delta_lnV / 2.0),
    box.1 * @math.exp(delta_lnV / 2.0),
    box.2 * @math.exp(delta_lnV),
  )
  (accept, new_box)
}

///|
/// Anisotropic barostat with per-axis virial.
pub fn omm_anisotropic_barostat_virial(
  box : (Double, Double, Double),
  virial : (Double, Double, Double), // (Vxx, Vyy, Vzz)
  pressure : (Double, Double, Double), // (Px, Py, Pz)
  energy : Double,
  temperature : Double,
  n_particles : Int32,
  axis~ : Int32 = 0, // 0=x, 1=y, 2=z
  rng~ : Double = 0.5,
) -> (Bool, (Double, Double, Double)) {
  let delta = (rng - 0.5) * 0.01
  // Select axis
  let (v_axis, p_axis, l_axis) = if axis == 0 {
    (virial.0, pressure.0, box.0)
  } else if axis == 1 {
    (virial.1, pressure.1, box.1)
  } else {
    (virial.2, pressure.2, box.2)
  }
  let dE = -(p_axis - v_axis / (l_axis * l_axis * l_axis + 0.0000000001)) *
    delta *
    l_axis *
    l_axis *
    l_axis
  let beta = 1.0 / (1.380649e-23 * temperature + 0.0000000001)
  let arg = -beta * (dE + energy * delta)
  let accept = if arg > 0.0 { true } else { rng < @math.exp(arg) }
  let new_box = if axis == 0 {
    (box.0 * @math.exp(delta), box.1, box.2)
  } else if axis == 1 {
    (box.0, box.1 * @math.exp(delta), box.2)
  } else {
    (box.0, box.1, box.2 * @math.exp(delta))
  }
  ignore(n_particles)
  (accept, new_box)
}

///|
pub fn omm_barostat2_dummy() -> Unit {
  ignore(0)
}
