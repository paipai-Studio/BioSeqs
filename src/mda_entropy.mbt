/// Quasi-harmonic configurational entropy.
///
/// MDAnalysis equivalent: `analysis.pca.QuasiHarmonicEntropy`
///
/// S_qh = (k_B/2) * sum(ln(2*pi*k_B*T / (m_i * omega_i^2)))
/// where omega_i are quasi-harmonic frequencies from PCA eigenvalues.

///|
/// Compute quasi-harmonic entropy from PCA eigenvalues.
///
/// eigenvalues: PCA eigenvalues (variance along each mode)
/// temperature: system temperature (K)
/// masses: atomic masses (for mass-weighting, optional)
/// kB: Boltzmann constant in appropriate units
///
/// Returns entropy in units of kB (i.e., S/kB)
pub fn mda_quasi_harmonic_entropy(
  eigenvalues : Array[Double],
  temperature~ : Double = 300.0,
  kB~ : Double = 1.0,
) -> Double {
  let n = eigenvalues.length()
  if n == 0 {
    return 0.0
  }
  // S_qh / kB = (1/2) * sum(ln(2*pi*kB*T / lambda_i))
  // where lambda_i are the eigenvalues (proportional to 1/omega_i^2)
  let kT = kB * temperature
  let mut entropy = 0.0
  for i = 0; i < n; i = i + 1 {
    if eigenvalues[i] > 0.0 {
      entropy = entropy +
        @math.ln(2.0 * 3.141592653589793 * kT / eigenvalues[i])
    }
  }
  0.5 * entropy
}

///|
/// Compute entropy from trajectory directly (via PCA).
/// traj: frames[atom][xyz]
pub fn mda_trajectory_entropy(
  traj : Array[Array[Array[Double]]],
  temperature~ : Double = 300.0,
) -> Double {
  let n_frames = traj.length()
  if n_frames < 2 {
    return 0.0
  }
  let n_atoms = traj[0].length()
  let dim = 3 * n_atoms
  // Compute mean structure
  let mean : Array[Double] = Array::make(dim, 0.0)
  for f = 0; f < n_frames; f = f + 1 {
    for a = 0; a < n_atoms; a = a + 1 {
      mean[3 * a] = mean[3 * a] + traj[f][a][0]
      mean[3 * a + 1] = mean[3 * a + 1] + traj[f][a][1]
      mean[3 * a + 2] = mean[3 * a + 2] + traj[f][a][2]
    }
  }
  for i = 0; i < dim; i = i + 1 {
    mean[i] = mean[i] / n_frames.to_double()
  }
  // Compute covariance matrix (diagonal approximation for efficiency)
  let variances : Array[Double] = Array::make(dim, 0.0)
  for f = 0; f < n_frames; f = f + 1 {
    for a = 0; a < n_atoms; a = a + 1 {
      let dx = traj[f][a][0] - mean[3 * a]
      let dy = traj[f][a][1] - mean[3 * a + 1]
      let dz = traj[f][a][2] - mean[3 * a + 2]
      variances[3 * a] = variances[3 * a] + dx * dx
      variances[3 * a + 1] = variances[3 * a + 1] + dy * dy
      variances[3 * a + 2] = variances[3 * a + 2] + dz * dz
    }
  }
  for i = 0; i < dim; i = i + 1 {
    variances[i] = variances[i] / (n_frames - 1).to_double()
  }
  mda_quasi_harmonic_entropy(variances, temperature~)
}

///|
/// Schlitter entropy approximation (upper bound).
/// S_schlitter = kB * sum(1/2 + ln(sqrt(2*pi*e / lambda_i)))
pub fn mda_schlitter_entropy(
  eigenvalues : Array[Double],
  kB~ : Double = 1.0,
) -> Double {
  let mut entropy = 0.0
  for i = 0; i < eigenvalues.length(); i = i + 1 {
    if eigenvalues[i] > 0.0 {
      entropy = entropy +
        0.5 +
        @math.ln(
          (2.0 * 3.141592653589793 * 2.718281828459045 / eigenvalues[i]).sqrt(),
        )
    }
  }
  kB * entropy
}

///|
/// Entropy difference between two states (e.g., bound vs unbound).
/// dS = S2 - S1
pub fn mda_entropy_difference(
  eigenvalues1 : Array[Double],
  eigenvalues2 : Array[Double],
  temperature~ : Double = 300.0,
) -> Double {
  let s1 = mda_quasi_harmonic_entropy(eigenvalues1, temperature~)
  let s2 = mda_quasi_harmonic_entropy(eigenvalues2, temperature~)
  s2 - s1
}

///|
/// Free energy estimate from entropy at given temperature.
/// dG = -T * dS
pub fn mda_entropy_free_energy(
  entropy : Double,
  temperature~ : Double = 300.0,
) -> Double {
  -temperature * entropy
}
