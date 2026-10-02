///|
/// ProDy-style normal mode analysis.
///
/// Extends GNM/ANM results with mode projections, deformation
/// vectors, collectivity, and participation ratios.

///|
/// Normal mode analysis result.
pub struct ModeAnalysis {
  n : Int
  dim : Int // n for GNM, 3n for ANM
  num_modes : Int
  eigenvalues : Array[Double]
  eigenvectors : Array[Array[Double]] // dim × num_modes
}

///|
/// Build ModeAnalysis from eigenvalues and eigenvectors.
pub fn ModeAnalysis::new(
  eigenvalues : Array[Double],
  eigenvectors : Array[Array[Double]],
  n : Int,
  dim : Int,
) -> ModeAnalysis {
  let num_modes = eigenvalues.length()
  ModeAnalysis::{ n, dim, num_modes, eigenvalues, eigenvectors }
}

///|
/// Build ModeAnalysis from a GNM (skip zero mode).
pub fn ModeAnalysis::from_gnm(gnm : GNM) -> ModeAnalysis {
  let n = gnm.num_nodes()
  let modes = gnm.slow_modes(n - 1)
  let num_modes = modes.length()
  let eigenvalues = Array::make(num_modes, 0.0)
  let eigenvectors : Array[Array[Double]] = Array::makei(n, { _ => Array::make(num_modes, 0.0) })
  for k = 0; k < num_modes; k = k + 1 {
    let (_, eval, evec) = modes[k]
    eigenvalues[k] = eval
    for i = 0; i < n; i = i + 1 {
      eigenvectors[i][k] = evec[i]
    }
  }
  ModeAnalysis::{ n, dim: n, num_modes, eigenvalues, eigenvectors }
}

///|
/// Build ModeAnalysis from an ANM (skip 6 zero modes).
pub fn ModeAnalysis::from_anm(anm : ANM) -> ModeAnalysis {
  let n = anm.num_nodes()
  let dim = 3 * n
  let modes = anm.slow_modes(dim - 6)
  let num_modes = modes.length()
  let eigenvalues = Array::make(num_modes, 0.0)
  let eigenvectors : Array[Array[Double]] = Array::makei(dim, { _ => Array::make(num_modes, 0.0) })
  for k = 0; k < num_modes; k = k + 1 {
    let (_, eval, evec) = modes[k]
    eigenvalues[k] = eval
    for i = 0; i < dim; i = i + 1 {
      eigenvectors[i][k] = evec[i]
    }
  }
  ModeAnalysis::{ n, dim, num_modes, eigenvalues, eigenvectors }
}

///|
/// Project a displacement vector onto mode k.
/// Returns the projection coefficient.
pub fn ModeAnalysis::project(
  self : ModeAnalysis,
  displacement : Array[Double],
  mode_idx : Int,
) -> Double {
  if mode_idx < 0 || mode_idx >= self.num_modes {
    return 0.0
  }
  let mut s = 0.0
  for i = 0; i < self.dim && i < displacement.length(); i = i + 1 {
    s = s + self.eigenvectors[i][mode_idx] * displacement[i]
  }
  s
}

///|
/// Deformation vector for mode k: the displacement of each node
/// when the structure moves along mode k.
/// Returns Array[Double] of length dim.
pub fn ModeAnalysis::deformation(
  self : ModeAnalysis,
  mode_idx : Int,
) -> Array[Double] {
  let result = Array::make(self.dim, 0.0)
  if mode_idx < 0 || mode_idx >= self.num_modes {
    return result
  }
  for i = 0; i < self.dim; i = i + 1 {
    result[i] = self.eigenvectors[i][mode_idx]
  }
  result
}

///|
/// Collectivity of mode k (Bruschweiler, 2003).
/// Measures how many atoms participate in the mode.
/// c_k = (1/N) * exp(-sum_i (v_k[i]^2 - <v_k^2>)^2 / (2*variance^2))
pub fn ModeAnalysis::collectivity(self : ModeAnalysis, mode_idx : Int) -> Double {
  if mode_idx < 0 || mode_idx >= self.num_modes {
    return 0.0
  }
  let n = self.n
  // For ANM, sum over 3 components per atom
  let sq_amplitudes = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut s = 0.0
    if self.dim == n {
      s = self.eigenvectors[i][mode_idx] * self.eigenvectors[i][mode_idx]
    } else {
      // ANM: 3 components
      s = self.eigenvectors[3 * i][mode_idx] * self.eigenvectors[3 * i][mode_idx]
      s = s + self.eigenvectors[3 * i + 1][mode_idx] * self.eigenvectors[3 * i + 1][mode_idx]
      s = s + self.eigenvectors[3 * i + 2][mode_idx] * self.eigenvectors[3 * i + 2][mode_idx]
    }
    sq_amplitudes[i] = s
  }
  // Mean
  let mut mean = 0.0
  for i = 0; i < n; i = i + 1 {
    mean = mean + sq_amplitudes[i]
  }
  mean = mean / n.to_double()
  if mean < 1.0e-20 {
    return 0.0
  }
  // Variance of log amplitudes
  let log_amps = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    if sq_amplitudes[i] > 1.0e-20 {
      log_amps[i] = @math.ln(sq_amplitudes[i])
    } else {
      log_amps[i] = -50.0
    }
  }
  let mut log_mean = 0.0
  for i = 0; i < n; i = i + 1 {
    log_mean = log_mean + log_amps[i]
  }
  log_mean = log_mean / n.to_double()
  let mut variance = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = log_amps[i] - log_mean
    variance = variance + diff * diff
  }
  variance = variance / n.to_double()
  // Collectivity
  let c = @math.exp(-variance) / n.to_double()
  c * n.to_double() // Normalize so it's in [0, 1]
}

///|
/// Participation ratio: PR_k = (sum_i v_k[i]^2)^2 / (N * sum_i v_k[i]^4)
/// PR = 1 for uniform participation, PR ~ 1/N for localized.
pub fn ModeAnalysis::participation_ratio(self : ModeAnalysis, mode_idx : Int) -> Double {
  if mode_idx < 0 || mode_idx >= self.num_modes {
    return 0.0
  }
  let n = self.n
  let mut sum_sq = 0.0
  let mut sum_4th = 0.0
  for i = 0; i < n; i = i + 1 {
    let mut s = 0.0
    if self.dim == n {
      s = self.eigenvectors[i][mode_idx] * self.eigenvectors[i][mode_idx]
    } else {
      s = self.eigenvectors[3 * i][mode_idx] * self.eigenvectors[3 * i][mode_idx]
      s = s + self.eigenvectors[3 * i + 1][mode_idx] * self.eigenvectors[3 * i + 1][mode_idx]
      s = s + self.eigenvectors[3 * i + 2][mode_idx] * self.eigenvectors[3 * i + 2][mode_idx]
    }
    sum_sq = sum_sq + s
    sum_4th = sum_4th + s * s
  }
  if sum_4th < 1.0e-20 {
    return 0.0
  }
  (sum_sq * sum_sq) / (n.to_double() * sum_4th)
}

///|
/// Correlation between predicted and experimental B-factors.
/// Returns Pearson correlation coefficient.
pub fn ModeAnalysis::bfactors_correlation(
  self : ModeAnalysis,
  experimental_bfactors : Array[Double],
) -> Double {
  let n = self.n
  if experimental_bfactors.length() < n {
    return 0.0
  }
  // Predict B-factors from all modes
  let predicted = Array::make(n, 0.0)
  for k = 0; k < self.num_modes; k = k + 1 {
    let lambda = self.eigenvalues[k]
    if lambda.abs() < 1.0e-10 {
      continue
    }
    for i = 0; i < n; i = i + 1 {
      let mut s = 0.0
      if self.dim == n {
        s = self.eigenvectors[i][k] * self.eigenvectors[i][k]
      } else {
        s = self.eigenvectors[3 * i][k] * self.eigenvectors[3 * i][k]
        s = s + self.eigenvectors[3 * i + 1][k] * self.eigenvectors[3 * i + 1][k]
        s = s + self.eigenvectors[3 * i + 2][k] * self.eigenvectors[3 * i + 2][k]
      }
      predicted[i] = predicted[i] + s / lambda
    }
  }
  // Pearson correlation
  let mut mean_p = 0.0
  let mut mean_e = 0.0
  for i = 0; i < n; i = i + 1 {
    mean_p = mean_p + predicted[i]
    mean_e = mean_e + experimental_bfactors[i]
  }
  mean_p = mean_p / n.to_double()
  mean_e = mean_e / n.to_double()
  let mut cov = 0.0
  let mut var_p = 0.0
  let mut var_e = 0.0
  for i = 0; i < n; i = i + 1 {
    let dp = predicted[i] - mean_p
    let de = experimental_bfactors[i] - mean_e
    cov = cov + dp * de
    var_p = var_p + dp * dp
    var_e = var_e + de * de
  }
  if var_p < 1.0e-20 || var_e < 1.0e-20 {
    return 0.0
  }
  cov / (var_p * var_e).sqrt()
}

///|
/// Number of modes.
pub fn ModeAnalysis::num_modes(self : ModeAnalysis) -> Int {
  self.num_modes
}

///|
/// Get eigenvalue for mode k.
pub fn ModeAnalysis::eigenvalue(self : ModeAnalysis, mode_idx : Int) -> Double {
  if mode_idx < 0 || mode_idx >= self.num_modes {
    return 0.0
  }
  self.eigenvalues[mode_idx]
}

///|
/// Get eigenvector for mode k.
pub fn ModeAnalysis::eigenvector(self : ModeAnalysis, mode_idx : Int) -> Array[Double] {
  self.deformation(mode_idx)
}

///|
/// Mode energy: lambda_k * |projection|^2
pub fn ModeAnalysis::mode_energy(
  self : ModeAnalysis,
  displacement : Array[Double],
  mode_idx : Int,
) -> Double {
  let proj = self.project(displacement, mode_idx)
  self.eigenvalues[mode_idx] * proj * proj
}

///|
/// Total energy across all modes.
pub fn ModeAnalysis::total_energy(
  self : ModeAnalysis,
  displacement : Array[Double],
) -> Double {
  let mut total = 0.0
  for k = 0; k < self.num_modes; k = k + 1 {
    total = total + self.mode_energy(displacement, k)
  }
  total
}
