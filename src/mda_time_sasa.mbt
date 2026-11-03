/// Time-dependent SASA over trajectory.
///
/// MDAnalysis equivalent: SASA computed per-frame in trajectory analysis
///
/// Tracks SASA changes over time, detects burial/exposure transitions.

///|
/// Compute SASA for a single frame using Shrake-Rupley algorithm.
/// coords: [atom][xyz], radii: van der Waals radii [atom]
/// n_points: number of sphere points (more = more accurate)
fn mda_frame_sasa(
  coords : Array[Array[Double]],
  radii : Array[Double],
  n_points : Int32,
) -> Array[Double] {
  let n_atoms = coords.length()
  let sasa : Array[Double] = Array::make(n_atoms, 0.0)
  if n_atoms == 0 {
    return sasa
  }
  // Generate sphere points (golden spiral)
  let points : Array[Array[Double]] = Array(capacity=0)
  let golden = 3.141592653589793 * (3.0 - 5.0.sqrt())
  for i = 0; i < n_points; i = i + 1 {
    let y = 1.0 - i.to_double() / (n_points - 1).to_double() * 2.0
    let r = (1.0 - y * y).sqrt()
    let theta = golden * i.to_double()
    points.push([@math.cos(theta) * r, y, @math.sin(theta) * r])
  }
  // For each atom, count exposed test points
  for i = 0; i < n_atoms; i = i + 1 {
    let ri = radii[i]
    let mut buried = 0
    for p = 0; p < n_points; p = p + 1 {
      let tx = coords[i][0] + ri * points[p][0]
      let ty = coords[i][1] + ri * points[p][1]
      let tz = coords[i][2] + ri * points[p][2]
      let mut is_buried = false
      for j = 0; j < n_atoms; j = j + 1 {
        if i == j {
          continue
        }
        let dx = tx - coords[j][0]
        let dy = ty - coords[j][1]
        let dz = tz - coords[j][2]
        let dist_sq = dx * dx + dy * dy + dz * dz
        if dist_sq < radii[j] * radii[j] {
          is_buried = true
          break
        }
      }
      if is_buried {
        buried = buried + 1
      }
    }
    let exposed = n_points - buried
    sasa[i] = 4.0 *
      3.141592653589793 *
      ri *
      ri *
      exposed.to_double() /
      n_points.to_double()
  }
  sasa
}

///|
/// Compute SASA over trajectory.
/// traj: [frame][atom][xyz], radii: [atom]
/// Returns [frame][atom] SASA values.
pub fn mda_time_dependent_sasa(
  traj : Array[Array[Array[Double]]],
  radii : Array[Double],
  n_points~ : Int32 = 50,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array(capacity=0)
  for f = 0; f < traj.length(); f = f + 1 {
    result.push(mda_frame_sasa(traj[f], radii, n_points))
  }
  result
}

///|
/// Compute total SASA per frame.
pub fn mda_total_sasa_trajectory(
  sasa_traj : Array[Array[Double]],
) -> Array[Double] {
  let result : Array[Double] = Array::make(sasa_traj.length(), 0.0)
  for f = 0; f < sasa_traj.length(); f = f + 1 {
    let mut sum = 0.0
    for a = 0; a < sasa_traj[f].length(); a = a + 1 {
      sum = sum + sasa_traj[f][a]
    }
    result[f] = sum
  }
  result
}

///|
/// Detect burial/exposure transitions for an atom.
/// Returns list of (frame, transition_type) where type is "burial" or "exposure".
pub fn mda_burial_exposure_events(
  sasa_per_frame : Array[Double],
  threshold~ : Double = 10.0,
) -> Array[(Int32, String)] {
  let events : Array[(Int32, String)] = Array(capacity=0)
  if sasa_per_frame.length() < 2 {
    return events
  }
  let mut was_exposed = sasa_per_frame[0] > threshold
  for f = 1; f < sasa_per_frame.length(); f = f + 1 {
    let is_exposed = sasa_per_frame[f] > threshold
    if is_exposed && !was_exposed {
      events.push((f, "Exposure"))
    }
    if !is_exposed && was_exposed {
      events.push((f, "Burial"))
    }
    was_exposed = is_exposed
  }
  events
}

///|
/// Compute SASA fluctuation (std dev over trajectory) per atom.
pub fn mda_sasa_fluctuation(sasa_traj : Array[Array[Double]]) -> Array[Double] {
  let n_frames = sasa_traj.length()
  if n_frames == 0 {
    return []
  }
  let n_atoms = sasa_traj[0].length()
  let result : Array[Double] = Array::make(n_atoms, 0.0)
  for a = 0; a < n_atoms; a = a + 1 {
    let mut mean = 0.0
    for f = 0; f < n_frames; f = f + 1 {
      mean = mean + sasa_traj[f][a]
    }
    mean = mean / n_frames.to_double()
    let mut variance = 0.0
    for f = 0; f < n_frames; f = f + 1 {
      let d = sasa_traj[f][a] - mean
      variance = variance + d * d
    }
    result[a] = (variance / n_frames.to_double()).sqrt()
  }
  result
}
