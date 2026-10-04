///|
/// ProDy-style conformation sampling from normal modes.
///
/// Generate alternative conformations by perturbing atomic coordinates
/// along normal modes.  Three sampling strategies are provided:
///
/// 1. **deform_mode** — deterministic deformation along a single mode
/// 2. **traverse_mode** — sweep a mode's amplitude from -range to +range
/// 3. **sample_anm/gnm_conformations** — stochastic sampling with
///    Gaussian random coefficients scaled by thermal energy kT/λ_k
///
/// Reference: Bakan & Bahar (2009) "Computing generation of alternative
/// conformations" in ProDy documentation.

// ============================================================================
// Configuration and result types
// ============================================================================

///|
/// Configuration for stochastic conformation sampling.
pub struct SamplingConfig {
  kT : Double // Thermal energy (default 1.0, in units of kT)
  num_modes : Int // Number of non-zero modes to use (0 = all available)
  seed : Int // Random seed for reproducibility
  scale : Double // Overall scale factor applied to displacements
}

///|
/// Default sampling config: kT=1.0, all modes, seed=42, scale=1.0.
pub fn SamplingConfig::default() -> SamplingConfig {
  SamplingConfig::{ kT: 1.0, num_modes: 0, seed: 42, scale: 1.0 }
}

///|
/// Create a SamplingConfig with custom parameters.
pub fn SamplingConfig::new(
  kT? : Double = 1.0,
  num_modes? : Int = 0,
  seed? : Int = 42,
  scale? : Double = 1.0,
) -> SamplingConfig {
  SamplingConfig::{ kT, num_modes, seed, scale }
}

///|
/// A single sampled conformation.
pub struct Conformation {
  coords : Array[Vector3] // 3D coordinates of all nodes
  coefficients : Array[Double] // Mode coefficients r_k used
  rmsd : Double // RMSD from the original structure
}

///|
/// Get the coordinates of this conformation.
pub fn Conformation::get_coords(self : Conformation) -> Array[Vector3] {
  self.coords
}

///|
/// Get the mode coefficients used to generate this conformation.
pub fn Conformation::get_coefficients(self : Conformation) -> Array[Double] {
  self.coefficients
}

///|
/// Get the RMSD from the original structure.
pub fn Conformation::get_rmsd(self : Conformation) -> Double {
  self.rmsd
}

///|
/// Construct a Conformation directly.
pub fn Conformation::new(
  coords : Array[Vector3],
  coefficients : Array[Double],
  rmsd : Double,
) -> Conformation {
  Conformation::{ coords, coefficients, rmsd }
}

// ============================================================================
// RNG — Box-Muller transform for Gaussian random numbers
// ============================================================================

///|
/// Simple LCG random number generator state.
priv struct SamplerRNG {
  state : Int64
}

///|
/// Create a new RNG with the given seed.
fn sampler_rng_new(seed : Int) -> SamplerRNG {
  let s = if seed == 0 { 1 } else { seed }
  SamplerRNG::{ state: s.to_int64() & 0x7FFFFFFFFFFFFFFF }
}

///|
/// Generate next uniform random number in [0, 1).
fn sampler_rng_uniform(state : SamplerRNG) -> (SamplerRNG, Double) {
  let next_state = (state.state * 6364136223846793005L + 1442695040888963407L) & 0x7FFFFFFFFFFFFFFF
  let value = next_state.to_double() / 9223372036854775808.0
  (SamplerRNG::{ state: next_state }, value)
}

///|
/// Box-Muller transform: generate a Gaussian N(0,1) random number.
fn sampler_rng_gaussian(state : SamplerRNG) -> (SamplerRNG, Double) {
  let (s1, u1) = sampler_rng_uniform(state)
  let (s2, u2) = sampler_rng_uniform(s1)
  let u1_safe = if u1 < 1.0e-10 { 1.0e-10 } else { u1 }
  let r = (-2.0 * @math.ln(u1_safe)).sqrt()
  let theta = 2.0 * 3.14159265358979323846 * u2
  (s2, r * @math.cos(theta))
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Calculate RMSD between two coordinate sets.
pub fn calc_rmsd(coords1 : Array[Vector3], coords2 : Array[Vector3]) -> Double {
  let n = coords1.length()
  if n == 0 || coords2.length() < n {
    return 0.0
  }
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let dx = coords1[i].x - coords2[i].x
    let dy = coords1[i].y - coords2[i].y
    let dz = coords1[i].z - coords2[i].z
    sum_sq = sum_sq + dx * dx + dy * dy + dz * dz
  }
  (sum_sq / n.to_double()).sqrt()
}

///|
/// Calculate the mean coordinate over a set of conformations.
pub fn mean_conformation(
  conformations : Array[Conformation]
) -> Array[Vector3]? {
  if conformations.length() == 0 {
    return None
  }
  let n = conformations.length()
  let n_atoms = conformations[0].coords.length()
  let result : Array[Vector3] = Array::makei(n_atoms, { _ => Vector3::zero() })
  for conf in conformations {
    for i = 0; i < n_atoms; i = i + 1 {
      result[i] = Vector3::new(
        result[i].x + conf.coords[i].x,
        result[i].y + conf.coords[i].y,
        result[i].z + conf.coords[i].z,
      )
    }
  }
  let nf = n.to_double()
  for i = 0; i < n_atoms; i = i + 1 {
    result[i] = Vector3::new(
      result[i].x / nf,
      result[i].y / nf,
      result[i].z / nf,
    )
  }
  Some(result)
}

///|
/// Calculate the RMSF of each atom across conformations.
pub fn ensemble_rmsf(conformations : Array[Conformation]) -> Array[Double] {
  if conformations.length() == 0 {
    return []
  }
  let n = conformations.length()
  let n_atoms = conformations[0].coords.length()
  let mean = match mean_conformation(conformations) {
    Some(m) => m
    None => return Array::make(n_atoms, 0.0)
  }
  let msf = Array::make(n_atoms, 0.0)
  for conf in conformations {
    for i = 0; i < n_atoms; i = i + 1 {
      let dx = conf.coords[i].x - mean[i].x
      let dy = conf.coords[i].y - mean[i].y
      let dz = conf.coords[i].z - mean[i].z
      msf[i] = msf[i] + dx * dx + dy * dy + dz * dz
    }
  }
  let nf = n.to_double()
  let rmsf = Array::make(n_atoms, 0.0)
  for i = 0; i < n_atoms; i = i + 1 {
    rmsf[i] = (msf[i] / nf).sqrt()
  }
  rmsf
}

///|
/// Calculate the RMSD matrix between all pairs of conformations.
pub fn ensemble_rmsd_matrix(
  conformations : Array[Conformation]
) -> Array[Array[Double]] {
  let n = conformations.length()
  let mat : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let r = calc_rmsd(conformations[i].coords, conformations[j].coords)
      mat[i][j] = r
      mat[j][i] = r
    }
  }
  mat
}

// ============================================================================
// Deterministic deformation
// ============================================================================

///|
/// Deform coordinates along a single mode eigenvector.
///
/// For ANM, the eigenvector has 3N elements (x,y,z per atom).
/// For GNM, the eigenvector has N elements (scalar per atom, applied
/// isotropically to all 3 coordinates).
///
/// `coords_new[i] = coords[i] + amplitude * eigenvector[i]`
pub fn deform_mode(
  coords : Array[Vector3],
  eigenvector : Array[Double],
  amplitude : Double,
) -> Array[Vector3] {
  let n = coords.length()
  let result : Array[Vector3] = Array::makei(n, { _ => Vector3::zero() })
  let ev_len = eigenvector.length()
  if ev_len == 3 * n {
    // ANM: 3N eigenvector
    for i = 0; i < n; i = i + 1 {
      let i3 = 3 * i
      result[i] = Vector3::new(
        coords[i].x + amplitude * eigenvector[i3],
        coords[i].y + amplitude * eigenvector[i3 + 1],
        coords[i].z + amplitude * eigenvector[i3 + 2],
      )
    }
  } else if ev_len == n {
    // GNM: N eigenvector (isotropic)
    for i = 0; i < n; i = i + 1 {
      let disp = amplitude * eigenvector[i]
      result[i] = Vector3::new(
        coords[i].x + disp,
        coords[i].y + disp,
        coords[i].z + disp,
      )
    }
  } else {
    for i = 0; i < n; i = i + 1 {
      result[i] = coords[i]
    }
  }
  result
}

///|
/// Traverse a mode: generate conformations sweeping the amplitude
/// from -range to +range in num_steps equal increments.
pub fn traverse_mode(
  coords : Array[Vector3],
  eigenvector : Array[Double],
  num_steps : Int,
  range : Double,
) -> Array[Array[Vector3]] {
  if num_steps <= 0 {
    return []
  }
  let result : Array[Array[Vector3]] = Array(capacity=num_steps)
  for step = 0; step < num_steps; step = step + 1 {
    let t = if num_steps == 1 {
      0.0
    } else {
      let frac = step.to_double() / (num_steps - 1).to_double()
      -range + 2.0 * range * frac
    }
    result.push(deform_mode(coords, eigenvector, t))
  }
  result
}

///|
/// Traverse a mode with explicit amplitude values.
pub fn traverse_mode_amplitudes(
  coords : Array[Vector3],
  eigenvector : Array[Double],
  amplitudes : Array[Double],
) -> Array[Array[Vector3]] {
  let result : Array[Array[Vector3]] = Array(capacity=amplitudes.length())
  for amp in amplitudes {
    result.push(deform_mode(coords, eigenvector, amp))
  }
  result
}

// ============================================================================
// Stochastic sampling — ANM
// ============================================================================

///|
/// Sample conformations from ANM normal modes.
///
/// For each sample, random Gaussian coefficients r_k ~ N(0, kT/λ_k)
/// are drawn and the new coordinates are:
///   coords_new = coords + scale * sum_k (r_k * v_k)
/// Zero modes (λ ≈ 0) are skipped.
pub fn sample_anm_conformations(
  coords : Array[Vector3],
  eigenvalues : Array[Double],
  eigenvectors : Array[Array[Double]],
  num_samples : Int,
  config : SamplingConfig,
) -> Array[Conformation] {
  let n = coords.length()
  let dim = 3 * n
  let num_modes_total = eigenvalues.length()

  let zero_threshold = 1.0e-8
  let mut start_mode = 0
  for k = 0; k < num_modes_total; k = k + 1 {
    if eigenvalues[k].abs() < zero_threshold {
      start_mode = k + 1
    } else {
      break
    }
  }
  let available_modes = num_modes_total - start_mode
  let num_modes = if config.num_modes > 0 && config.num_modes < available_modes {
    config.num_modes
  } else {
    available_modes
  }

  let mut rng = sampler_rng_new(config.seed)
  let conformations : Array[Conformation] = Array(capacity=num_samples)
  for _ in 0..<num_samples {
    let disp = Array::make(dim, 0.0)
    let coeffs = Array::make(num_modes, 0.0)
    for k = 0; k < num_modes; k = k + 1 {
      let mode_idx = start_mode + k
      let lambda = eigenvalues[mode_idx]
      if lambda.abs() < zero_threshold {
        coeffs[k] = 0.0
        continue
      }
      let (new_rng, gauss) = sampler_rng_gaussian(rng)
      rng = new_rng
      let sigma = (config.kT / lambda).sqrt()
      let r_k = gauss * sigma * config.scale
      coeffs[k] = r_k
      for i = 0; i < dim; i = i + 1 {
        disp[i] = disp[i] + r_k * eigenvectors[i][mode_idx]
      }
    }
    let new_coords : Array[Vector3] = Array::makei(n, { _ => Vector3::zero() })
    for i = 0; i < n; i = i + 1 {
      let i3 = 3 * i
      new_coords[i] = Vector3::new(
        coords[i].x + disp[i3],
        coords[i].y + disp[i3 + 1],
        coords[i].z + disp[i3 + 2],
      )
    }
    let rmsd = calc_rmsd(coords, new_coords)
    conformations.push(Conformation::{ coords: new_coords, coefficients: coeffs, rmsd })
  }
  conformations
}

// ============================================================================
// Stochastic sampling — GNM
// ============================================================================

///|
/// Sample conformations from GNM normal modes.
///
/// GNM eigenvectors are N-dimensional (one scalar per node).  Each scalar
/// is applied isotropically to all three Cartesian components.
pub fn sample_gnm_conformations(
  coords : Array[Vector3],
  eigenvalues : Array[Double],
  eigenvectors : Array[Array[Double]],
  num_samples : Int,
  config : SamplingConfig,
) -> Array[Conformation] {
  let n = coords.length()
  let num_modes_total = eigenvalues.length()

  let zero_threshold = 1.0e-8
  let mut start_mode = 0
  for k = 0; k < num_modes_total; k = k + 1 {
    if eigenvalues[k].abs() < zero_threshold {
      start_mode = k + 1
    } else {
      break
    }
  }
  let available_modes = num_modes_total - start_mode
  let num_modes = if config.num_modes > 0 && config.num_modes < available_modes {
    config.num_modes
  } else {
    available_modes
  }

  let mut rng = sampler_rng_new(config.seed)
  let conformations : Array[Conformation] = Array(capacity=num_samples)
  for _ in 0..<num_samples {
    let disp = Array::make(n, 0.0)
    let coeffs = Array::make(num_modes, 0.0)
    for k = 0; k < num_modes; k = k + 1 {
      let mode_idx = start_mode + k
      let lambda = eigenvalues[mode_idx]
      if lambda.abs() < zero_threshold {
        coeffs[k] = 0.0
        continue
      }
      let (new_rng, gauss) = sampler_rng_gaussian(rng)
      rng = new_rng
      let sigma = (config.kT / lambda).sqrt()
      let r_k = gauss * sigma * config.scale
      coeffs[k] = r_k
      for i = 0; i < n; i = i + 1 {
        disp[i] = disp[i] + r_k * eigenvectors[i][mode_idx]
      }
    }
    let new_coords : Array[Vector3] = Array::makei(n, { _ => Vector3::zero() })
    for i = 0; i < n; i = i + 1 {
      new_coords[i] = Vector3::new(
        coords[i].x + disp[i],
        coords[i].y + disp[i],
        coords[i].z + disp[i],
      )
    }
    let rmsd = calc_rmsd(coords, new_coords)
    conformations.push(Conformation::{ coords: new_coords, coefficients: coeffs, rmsd })
  }
  conformations
}

// ============================================================================
// Batch sampling from ANM/GNM objects directly
// ============================================================================

///|
/// Sample conformations directly from an ANM object.
pub fn sample_from_anm(
  anm : ANM,
  coords : Array[Vector3],
  num_samples : Int,
  config : SamplingConfig,
) -> Array[Conformation] {
  let eigenvalues = anm.get_eigenvalues()
  let eigenvectors = anm.get_eigenvectors()
  sample_anm_conformations(coords, eigenvalues, eigenvectors, num_samples, config)
}

///|
/// Sample conformations directly from a GNM object.
pub fn sample_from_gnm(
  gnm : GNM,
  coords : Array[Vector3],
  num_samples : Int,
  config : SamplingConfig,
) -> Array[Conformation] {
  let eigenvalues = gnm.get_eigenvalues()
  let eigenvectors = gnm.get_eigenvectors()
  sample_gnm_conformations(coords, eigenvalues, eigenvectors, num_samples, config)
}

// ============================================================================
// Analysis helpers
// ============================================================================

///|
/// Calculate the average RMSD of a set of conformations from the original.
pub fn average_rmsd(conformations : Array[Conformation]) -> Double {
  if conformations.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for conf in conformations {
    sum = sum + conf.rmsd
  }
  sum / conformations.length().to_double()
}

///|
/// Calculate the standard deviation of RMSD values.
pub fn rmsd_stddev(conformations : Array[Conformation]) -> Double {
  let n = conformations.length()
  if n < 2 {
    return 0.0
  }
  let avg = average_rmsd(conformations)
  let mut sum_sq = 0.0
  for conf in conformations {
    let diff = conf.rmsd - avg
    sum_sq = sum_sq + diff * diff
  }
  (sum_sq / n.to_double()).sqrt()
}

///|
/// Find the conformation with the maximum RMSD.
pub fn max_rmsd_conformation(
  conformations : Array[Conformation]
) -> Conformation? {
  if conformations.length() == 0 {
    return None
  }
  let mut best = conformations[0]
  for i = 1; i < conformations.length(); i = i + 1 {
    if conformations[i].rmsd > best.rmsd {
      best = conformations[i]
    }
  }
  Some(best)
}

///|
/// Find the conformation with the minimum RMSD.
pub fn min_rmsd_conformation(
  conformations : Array[Conformation]
) -> Conformation? {
  if conformations.length() == 0 {
    return None
  }
  let mut best = conformations[0]
  for i = 1; i < conformations.length(); i = i + 1 {
    if conformations[i].rmsd < best.rmsd {
      best = conformations[i]
    }
  }
  Some(best)
}
