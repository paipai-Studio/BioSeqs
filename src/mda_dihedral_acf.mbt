/// Dihedral angle autocorrelation function.
///
/// MDAnalysis equivalent: `analysis.dihedrals.DihedralAutocorrelation`
///
/// C(t) = <cos(phi(0) - phi(t))> or <phi(0)*phi(t)> / <phi(0)^2>

///|
/// Compute dihedral autocorrelation using cosine form.
/// dihedrals: dihedral angle (radians) at each frame
/// Returns C(t) = <cos(phi(0) - phi(t))>
pub fn mda_dihedral_autocorrelation_cos(
  dihedrals : Array[Double],
) -> Array[Double] {
  let n = dihedrals.length()
  let result : Array[Double] = Array::make(n, 0.0)
  if n == 0 {
    return result
  }
  for lag = 0; lag < n; lag = lag + 1 {
    let mut sum = 0.0
    let mut count = 0.0
    let t_max = n - lag
    for t = 0; t < t_max; t = t + 1 {
      sum = sum + @math.cos(dihedrals[t] - dihedrals[t + lag])
      count = count + 1.0
    }
    if count > 0.0 {
      result[lag] = sum / count
    }
  }
  result
}

///|
/// Compute dihedral autocorrelation using normalized product form.
/// C(t) = <delta_phi(0) * delta_phi(t)> / <delta_phi(0)^2>
/// where delta_phi = phi - <phi>
pub fn mda_dihedral_autocorrelation_product(
  dihedrals : Array[Double],
) -> Array[Double] {
  let n = dihedrals.length()
  let result : Array[Double] = Array::make(n, 0.0)
  if n == 0 {
    return result
  }
  // Compute mean
  let mut mean = 0.0
  for i = 0; i < n; i = i + 1 {
    mean = mean + dihedrals[i]
  }
  mean = mean / n.to_double()
  // Compute variance (normalization)
  let mut variance = 0.0
  for i = 0; i < n; i = i + 1 {
    let d = dihedrals[i] - mean
    variance = variance + d * d
  }
  if variance == 0.0 {
    return result
  }
  // Compute ACF
  for lag = 0; lag < n; lag = lag + 1 {
    let mut sum = 0.0
    let t_max = n - lag
    for t = 0; t < t_max; t = t + 1 {
      sum = sum + (dihedrals[t] - mean) * (dihedrals[t + lag] - mean)
    }
    result[lag] = sum / variance
  }
  result
}

///|
/// Compute autocorrelation for multiple dihedrals (e.g., all backbone phi).
/// dihedrals_per_residue: [residue][frame] -> angle
/// Returns averaged ACF over all residues.
pub fn mda_dihedral_autocorrelation_ensemble(
  dihedrals_per_residue : Array[Array[Double]],
) -> Array[Double] {
  let n_res = dihedrals_per_residue.length()
  if n_res == 0 {
    return []
  }
  let n_frames = dihedrals_per_residue[0].length()
  let result : Array[Double] = Array::make(n_frames, 0.0)
  for r = 0; r < n_res; r = r + 1 {
    let acf = mda_dihedral_autocorrelation_cos(dihedrals_per_residue[r])
    for i = 0; i < n_frames; i = i + 1 {
      result[i] = result[i] + acf[i]
    }
  }
  for i = 0; i < n_frames; i = i + 1 {
    result[i] = result[i] / n_res.to_double()
  }
  result
}

///|
/// Estimate correlation time from ACF.
/// tau = integral(C(t) dt) until C drops below 1/e
pub fn mda_correlation_time(acf : Array[Double], dt~ : Double = 1.0) -> Double {
  let n = acf.length()
  if n < 2 {
    return 0.0
  }
  let threshold = 1.0 / 2.718281828459045
  let mut tau = 0.0
  for i = 0; i < n - 1; i = i + 1 {
    if acf[i] < threshold {
      break
    }
    tau = tau + (acf[i] + acf[i + 1]) * 0.5 * dt
  }
  tau
}

///|
/// Estimate relaxation time by exponential fit.
/// Fit C(t) = exp(-t/tau) -> ln(C(t)) = -t/tau
pub fn mda_relaxation_time(acf : Array[Double], dt~ : Double = 1.0) -> Double {
  let n = acf.length()
  if n < 2 {
    return 0.0
  }
  // Linear regression on ln(C) vs t for early part
  let max_fit = if n > 10 { n / 4 } else { n / 2 }
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_xy = 0.0
  let mut sum_x2 = 0.0
  let mut count = 0.0
  for i = 1; i < max_fit; i = i + 1 {
    if acf[i] > 0.0 {
      let t = i.to_double() * dt
      let y = @math.ln(acf[i])
      sum_x = sum_x + t
      sum_y = sum_y + y
      sum_xy = sum_xy + t * y
      sum_x2 = sum_x2 + t * t
      count = count + 1.0
    }
  }
  if count < 2.0 {
    return 0.0
  }
  let denom = count * sum_x2 - sum_x * sum_x
  if denom == 0.0 {
    return 0.0
  }
  let slope = (count * sum_xy - sum_x * sum_y) / denom
  if slope >= 0.0 {
    return 0.0
  }
  -1.0 / slope
}
