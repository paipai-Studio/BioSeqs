///|
/// Essential Dynamics Analysis (EDA).
/// PCA-based analysis of conformational ensembles or trajectories.
/// Projects motions onto the essential subspace (top principal components).
/// Based on ProDy's `EDA` class.

///|
pub struct EssentialDynamics {
  n_conformations : Int
  n_coords : Int // 3 * n_atoms
  mean_coords : Array[Double] // mean structure
  eigenvalues : Array[Double] // sorted descending
  eigenvectors : Array[Array[Double]] // principal components
  projections : Array[Array[Double]] // conformation projections onto PCs
  total_variance : Double
}

///|
/// Create EDA from an ensemble of conformations.
/// Each conformation is a flat array of 3N coordinates.
pub fn EssentialDynamics::new(
  conformations : Array[Array[Double]],
) -> EssentialDynamics? {
  let n_conf = conformations.length()
  if n_conf < 2 {
    return None
  }
  let n_coords = conformations[0].length()
  if n_coords < 3 {
    return None
  }
  // Compute mean
  let mean = Array::make(n_coords, 0.0)
  for conf in conformations {
    for i = 0; i < n_coords; i = i + 1 {
      mean[i] = mean[i] + conf[i]
    }
  }
  for i = 0; i < n_coords; i = i + 1 {
    mean[i] = mean[i] / n_conf.to_double()
  }
  // Build covariance matrix (n_coords x n_coords)
  // For efficiency, if n_conf < n_coords, use the dual approach
  // Cov = (1/(n-1)) * X^T * X where X is centered
  let n_eff = n_coords
  let use_dual = n_conf < n_coords
  let dim = if use_dual { n_conf } else { n_coords }
  // Build the smaller matrix for eigen-decomposition
  let cov = Array::makei(dim, { _ => Array::make(dim, 0.0) })
  if use_dual {
    // Build X * X^T (n_conf x n_conf)
    for i = 0; i < n_conf; i = i + 1 {
      for j = 0; j < n_conf; j = j + 1 {
        let mut s = 0.0
        for kk = 0; kk < n_coords; kk = kk + 1 {
          let xi = conformations[i][kk] - mean[kk]
          let xj = conformations[j][kk] - mean[kk]
          s = s + xi * xj
        }
        cov[i][j] = s / (n_conf - 1).to_double()
      }
    }
  } else {
    // Build X^T * X (n_coords x n_coords)
    for i = 0; i < n_coords; i = i + 1 {
      for j = 0; j < n_coords; j = j + 1 {
        let mut s = 0.0
        for kk = 0; kk < n_conf; kk = kk + 1 {
          let xi = conformations[kk][i] - mean[i]
          let xj = conformations[kk][j] - mean[j]
          s = s + xi * xj
        }
        cov[i][j] = s / (n_conf - 1).to_double()
      }
    }
  }
  ignore(n_eff)
  // Eigen-decomposition
  let (evals, evecs) = linalg_jacobi(cov, dim, 200)
  // Sort by eigenvalue (descending) using insertion sort
  let indices = Array::makei(dim, { i => i })
  for sort_i = 1; sort_i < dim; sort_i = sort_i + 1 {
    let key = indices[sort_i]
    let key_val = evals[key]
    let mut sort_j = sort_i - 1
    while sort_j >= 0 && evals[indices[sort_j]] < key_val {
      indices[sort_j + 1] = indices[sort_j]
      sort_j = sort_j - 1
    }
    indices[sort_j + 1] = key
  }
  let sorted_evals = Array::make(dim, 0.0)
  let sorted_evecs = Array::makei(dim, { _ => Array::make(dim, 0.0) })
  for idx = 0; idx < dim; idx = idx + 1 {
    let src = indices[idx]
    sorted_evals[idx] = evals[src]
    for j = 0; j < dim; j = j + 1 {
      sorted_evecs[idx][j] = evecs[src][j]
    }
  }
  // If dual approach, convert eigenvectors back to n_coords space
  let (final_evecs, final_evals) = if use_dual {
    // V_full = X^T * V_dual / sqrt(lambda)
    let fe = Array::makei(dim, { _ => Array::make(n_coords, 0.0) })
    for k = 0; k < dim; k = k + 1 {
      let lam = sorted_evals[k]
      if lam > 1.0e-10 {
        let sqrt_lam = lam.sqrt()
        for j = 0; j < n_coords; j = j + 1 {
          let mut s = 0.0
          for i = 0; i < n_conf; i = i + 1 {
            s = s + (conformations[i][j] - mean[j]) * sorted_evecs[k][i]
          }
          fe[k][j] = s / sqrt_lam
        }
      }
    }
    (fe, sorted_evals)
  } else {
    (sorted_evecs, sorted_evals)
  }
  // Project conformations onto PCs
  let n_pcs = if dim < n_coords { dim } else { n_coords }
  let projections = Array::makei(n_conf, { _ => Array::make(n_pcs, 0.0) })
  for c = 0; c < n_conf; c = c + 1 {
    for k = 0; k < n_pcs; k = k + 1 {
      let mut s = 0.0
      for j = 0; j < n_coords; j = j + 1 {
        s = s + (conformations[c][j] - mean[j]) * final_evecs[k][j]
      }
      projections[c][k] = s
    }
  }
  // Total variance
  let mut total_var = 0.0
  for i = 0; i < final_evals.length(); i = i + 1 {
    if final_evals[i] > 0.0 {
      total_var = total_var + final_evals[i]
    }
  }
  Some(EssentialDynamics::{
    n_conformations: n_conf,
    n_coords,
    mean_coords: mean,
    eigenvalues: final_evals,
    eigenvectors: final_evecs,
    projections,
    total_variance: total_var,
  })
}

///|
/// Get eigenvalues (variances along PCs).
pub fn EssentialDynamics::get_eigenvalues(self : EssentialDynamics) -> Array[Double] {
  self.eigenvalues
}

///|
/// Get eigenvectors (principal components).
pub fn EssentialDynamics::get_eigenvectors(self : EssentialDynamics) -> Array[Array[Double]] {
  self.eigenvectors
}

///|
/// Get projections of conformations onto PCs.
pub fn EssentialDynamics::get_projections(self : EssentialDynamics) -> Array[Array[Double]] {
  self.projections
}

///|
/// Get mean coordinates.
pub fn EssentialDynamics::get_mean(self : EssentialDynamics) -> Array[Double] {
  self.mean_coords
}

///|
/// Variance explained by the first k principal components.
/// Returns fraction of total variance.
pub fn EssentialDynamics::variance_explained(
  self : EssentialDynamics,
  k : Int,
) -> Double {
  if self.total_variance < 1.0e-20 {
    return 0.0
  }
  let max_k = if k < self.eigenvalues.length() { k } else {
    self.eigenvalues.length()
  }
  let mut sum = 0.0
  for i = 0; i < max_k; i = i + 1 {
    if self.eigenvalues[i] > 0.0 {
      sum = sum + self.eigenvalues[i]
    }
  }
  sum / self.total_variance
}

///|
/// Number of essential modes (those with significant variance).
pub fn EssentialDynamics::num_essential_modes(
  self : EssentialDynamics,
  threshold? : Double = 0.95,
) -> Int {
  let total = self.total_variance
  if total < 1.0e-20 {
    return 0
  }
  let mut cumsum = 0.0
  for i = 0; i < self.eigenvalues.length(); i = i + 1 {
    cumsum = cumsum + self.eigenvalues[i]
    if cumsum / total >= threshold {
      return i + 1
    }
  }
  self.eigenvalues.length()
}

///|
/// Project a new conformation onto the essential subspace.
pub fn EssentialDynamics::project(
  self : EssentialDynamics,
  coords : Array[Double],
  k : Int,
) -> Array[Double] {
  let max_k = if k < self.eigenvalues.length() { k } else {
    self.eigenvalues.length()
  }
  let result = Array::make(max_k, 0.0)
  for i = 0; i < max_k; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < self.n_coords; j = j + 1 {
      s = s + (coords[j] - self.mean_coords[j]) * self.eigenvectors[i][j]
    }
    result[i] = s
  }
  result
}

///|
/// Reconstruct a conformation from its projections.
/// coords = mean + sum_k proj_k * v_k
pub fn EssentialDynamics::reconstruct(
  self : EssentialDynamics,
  projections : Array[Double],
) -> Array[Double] {
  let result = Array::make(self.n_coords, 0.0)
  // Start with mean
  for i = 0; i < self.n_coords; i = i + 1 {
    result[i] = self.mean_coords[i]
  }
  // Add contributions
  let n_pc = if projections.length() < self.eigenvalues.length() {
    projections.length()
  } else {
    self.eigenvalues.length()
  }
  for k = 0; k < n_pc; k = k + 1 {
    for j = 0; j < self.n_coords; j = j + 1 {
      result[j] = result[j] + projections[k] * self.eigenvectors[k][j]
    }
  }
  result
}

///|
/// RMSF (Root Mean Square Fluctuation) per atom from EDA.
pub fn EssentialDynamics::rmsf(self : EssentialDynamics) -> Array[Double] {
  let n_atoms = self.n_coords / 3
  let result = Array::make(n_atoms, 0.0)
  for k = 0; k < self.eigenvalues.length(); k = k + 1 {
    let lam = self.eigenvalues[k]
    if lam < 1.0e-20 {
      continue
    }
    for i = 0; i < n_atoms; i = i + 1 {
      let i3 = 3 * i
      result[i] = result[i] + lam * (
        self.eigenvectors[k][i3] * self.eigenvectors[k][i3] +
        self.eigenvectors[k][i3 + 1] * self.eigenvectors[k][i3 + 1] +
        self.eigenvectors[k][i3 + 2] * self.eigenvectors[k][i3 + 2]
      )
    }
  }
  for i = 0; i < n_atoms; i = i + 1 {
    result[i] = result[i].sqrt()
  }
  result
}
