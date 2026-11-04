/// Lipid order parameter SCD (deuterium order parameter).
///
/// MDAnalysis equivalent: `analysis.leaflet.LipidOrderParameter`
///
/// SCD = -0.5 * <3*cos²(theta) - 1>
/// where theta is the angle between the C-H bond vector
/// and the membrane normal (typically z-axis).

///|
/// Compute SCD for a single C-H bond over trajectory.
///
/// ch_vectors: C-H bond vectors per frame [frame][xyz]
/// normal: membrane normal vector [xyz] (default [0, 0, 1])
pub fn mda_lipid_scd(
  ch_vectors : Array[Array[Double]],
  normal~ : Array[Double] = [0.0, 0.0, 1.0],
) -> Double {
  let n = ch_vectors.length()
  if n == 0 {
    return 0.0
  }
  // Normalize normal vector
  let nm_mag = (normal[0] * normal[0] +
  normal[1] * normal[1] +
  normal[2] * normal[2]).sqrt()
  if nm_mag == 0.0 {
    return 0.0
  }
  let nx = normal[0] / nm_mag
  let ny = normal[1] / nm_mag
  let nz = normal[2] / nm_mag
  // Average <3*cos²(theta) - 1>
  let mut sum = 0.0
  for f = 0; f < n; f = f + 1 {
    let v = ch_vectors[f]
    let v_mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    if v_mag == 0.0 {
      continue
    }
    let cos_theta = (v[0] * nx + v[1] * ny + v[2] * nz) / v_mag
    sum = sum + 3.0 * cos_theta * cos_theta - 1.0
  }
  -0.5 * sum / n.to_double()
}

///|
/// Compute SCD profile for a lipid tail.
///
/// ch_bonds_per_carbon: for each carbon in the tail, the C-H bond vectors over trajectory
/// Returns array of SCD values, one per carbon.
pub fn mda_lipid_order_profile(
  ch_bonds_per_carbon : Array[Array[Array[Double]]],
  normal~ : Array[Double] = [0.0, 0.0, 1.0],
) -> Array[Double] {
  let n_carbons = ch_bonds_per_carbon.length()
  let result : Array[Double] = Array::make(n_carbons, 0.0)
  for c = 0; c < n_carbons; c = c + 1 {
    result[c] = mda_lipid_scd(ch_bonds_per_carbon[c], normal~)
  }
  result
}

///|
/// Compute SCD from carbon positions and hydrogen positions.
///
/// carbon_traj: carbon positions per frame [frame][xyz]
/// hydrogen_traj: hydrogen positions per frame [frame][xyz]
pub fn mda_lipid_scd_from_positions(
  carbon_traj : Array[Array[Double]],
  hydrogen_traj : Array[Array[Double]],
  normal~ : Array[Double] = [0.0, 0.0, 1.0],
) -> Double {
  let n = if carbon_traj.length() < hydrogen_traj.length() {
    carbon_traj.length()
  } else {
    hydrogen_traj.length()
  }
  if n == 0 {
    return 0.0
  }
  // Compute C-H vectors
  let ch_vectors : Array[Array[Double]] = Array(capacity=0)
  for f = 0; f < n; f = f + 1 {
    ch_vectors.push([
      hydrogen_traj[f][0] - carbon_traj[f][0],
      hydrogen_traj[f][1] - carbon_traj[f][1],
      hydrogen_traj[f][2] - carbon_traj[f][2],
    ])
  }
  mda_lipid_scd(ch_vectors, normal~)
}

///|
/// Average SCD across all carbons (overall tail order).
pub fn mda_lipid_average_order(scd_profile : Array[Double]) -> Double {
  if scd_profile.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for s in scd_profile {
    sum = sum + s
  }
  sum / scd_profile.length().to_double()
}

///|
/// Classify lipid order: ordered (|SCD| > 0.4), disordered (< 0.2), intermediate.
pub fn mda_lipid_order_state(scd : Double) -> String {
  let abs_scd = if scd < 0.0 { -scd } else { scd }
  if abs_scd > 0.4 {
    "Ordered"
  } else if abs_scd < 0.2 {
    "Disordered"
  } else {
    "Intermediate"
  }
}
