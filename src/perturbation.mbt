///|
/// Perturbation response analysis (ProDy PRT / PRS).
///
/// Applies a unit perturbation force to each residue and computes
/// the displacement response of all other residues using the
/// inverse of the connectivity matrix (GNM Kirchhoff or ANM Hessian).

///|
/// Perturbation response analysis result.
pub struct PerturbationResponse {
  n : Int
  // Response matrix R[i][j] = displacement of residue i when
  // a unit perturbation is applied to residue j.
  response : Array[Array[Double]]
}

///|
/// Build PRS from a GNM: R = Kirchhoff^+ (pseudo-inverse).
/// Uses eigendecomposition: R = sum_k (1/lambda_k) * v_k * v_k^T
pub fn PerturbationResponse::from_gnm(gnm : GNM) -> PerturbationResponse {
  let n = gnm.num_nodes()
  let evals = gnm.get_eigenvalues()
  let evecs = gnm.get_eigenvectors()
  let response : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  for k = 0; k < evals.length(); k = k + 1 {
    let lambda = evals[k]
    if lambda.abs() < 1.0e-10 {
      continue // Skip zero mode
    }
    let inv_lambda = 1.0 / lambda
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        response[i][j] = response[i][j] + inv_lambda * evecs[i][k] * evecs[j][k]
      }
    }
  }
  PerturbationResponse::{ n, response }
}

///|
/// Build PRS from an ANM: R = Hessian^+ (pseudo-inverse).
/// The Hessian is 3N×3N, so the response is 3N×3N.
/// We return the per-residue response by summing x/y/z components.
pub fn PerturbationResponse::from_anm(anm : ANM) -> PerturbationResponse {
  let n = anm.num_nodes()
  let dim = 3 * n
  let evals = anm.get_eigenvalues()
  let evecs = anm.get_eigenvectors()
  // Full 3N×3N pseudo-inverse
  let full_response : Array[Array[Double]] = Array::makei(dim, { _ => Array::make(dim, 0.0) })
  for k = 0; k < evals.length(); k = k + 1 {
    let lambda = evals[k]
    if lambda.abs() < 1.0e-10 {
      continue
    }
    let inv_lambda = 1.0 / lambda
    for i = 0; i < dim; i = i + 1 {
      for j = 0; j < dim; j = j + 1 {
        full_response[i][j] = full_response[i][j] + inv_lambda * evecs[i][k] * evecs[j][k]
      }
    }
  }
  // Reduce to N×N: response[i][j] = sum of 3×3 block
  let response : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let mut s = 0.0
      for a = 0; a < 3; a = a + 1 {
        for b = 0; b < 3; b = b + 1 {
          s = s + full_response[3 * i + a][3 * j + b]
        }
      }
      response[i][j] = s
    }
  }
  PerturbationResponse::{ n, response }
}

///|
/// Get the response matrix R[i][j].
pub fn PerturbationResponse::response_matrix(
  self : PerturbationResponse,
) -> Array[Array[Double]] {
  self.response
}

///|
/// Sensitivity of residue i: average response of all residues
/// when residue i is perturbed.
/// sensitivity[i] = (1/N) * sum_j R[j][i]
pub fn PerturbationResponse::sensitivity(
  self : PerturbationResponse,
  residue_idx : Int,
) -> Double {
  if residue_idx < 0 || residue_idx >= self.n {
    return 0.0
  }
  let mut s = 0.0
  for j = 0; j < self.n; j = j + 1 {
    s = s + self.response[j][residue_idx].abs()
  }
  s / self.n.to_double()
}

///|
/// Effectiveness: average displacement of residue i when any
/// residue is perturbed.
/// effectiveness[i] = (1/N) * sum_j R[i][j]
pub fn PerturbationResponse::effectiveness(
  self : PerturbationResponse,
  residue_idx : Int,
) -> Double {
  if residue_idx < 0 || residue_idx >= self.n {
    return 0.0
  }
  let mut s = 0.0
  for j = 0; j < self.n; j = j + 1 {
    s = s + self.response[residue_idx][j].abs()
  }
  s / self.n.to_double()
}

///|
/// Perturbation profile: sensitivity for all residues.
pub fn PerturbationResponse::perturbation_profile(
  self : PerturbationResponse,
) -> Array[Double] {
  let result = Array::make(self.n, 0.0)
  for i = 0; i < self.n; i = i + 1 {
    result[i] = self.sensitivity(i)
  }
  result
}

///|
/// Identify hotspot residues: those with highest sensitivity.
/// Returns indices sorted by sensitivity descending.
pub fn PerturbationResponse::hotspot_residues(
  self : PerturbationResponse,
  top_k : Int,
) -> Array[Int] {
  let sensitivities = Array::makei(self.n, { i => (i, self.sensitivity(i)) })
  // Sort by sensitivity descending (simple selection sort)
  let result : Array[Int] = []
  let used = Array::make(self.n, false)
  for iter = 0; iter < top_k && iter < self.n; iter = iter + 1 {
    let mut best_idx = -1
    let mut best_val = -1.0
    for i = 0; i < self.n; i = i + 1 {
      if not(used[i]) && sensitivities[i].1 > best_val {
        best_val = sensitivities[i].1
        best_idx = i
      }
    }
    if best_idx >= 0 {
      used[best_idx] = true
      result.push(best_idx)
    }
  }
  result
}

///|
/// Cross-correlation between residues i and j from PRS.
/// Normalized: C[i][j] = R[i][j] / sqrt(R[i][i] * R[j][j])
pub fn PerturbationResponse::cross_correlation(
  self : PerturbationResponse,
) -> Array[Array[Double]] {
  let n = self.n
  let corr : Array[Array[Double]] = Array::makei(n, { _ => Array::make(n, 0.0) })
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      let denom = (self.response[i][i] * self.response[j][j]).abs()
      if denom > 1.0e-20 {
        corr[i][j] = self.response[i][j] / denom.sqrt()
      }
    }
  }
  corr
}

///|
/// Number of residues.
pub fn PerturbationResponse::num_residues(self : PerturbationResponse) -> Int {
  self.n
}

///|
/// Dynamic flexibility index (DFI) for residue i.
/// DFI[i] = sum_j R[i][j] / sum_jk R[j][k]
pub fn PerturbationResponse::dfi(self : PerturbationResponse, residue_idx : Int) -> Double {
  if residue_idx < 0 || residue_idx >= self.n {
    return 0.0
  }
  let mut row_sum = 0.0
  for j = 0; j < self.n; j = j + 1 {
    row_sum = row_sum + self.response[residue_idx][j].abs()
  }
  let mut total = 0.0
  for i = 0; i < self.n; i = i + 1 {
    for j = 0; j < self.n; j = j + 1 {
      total = total + self.response[i][j].abs()
    }
  }
  if total < 1.0e-20 {
    return 0.0
  }
  row_sum / total
}

///|
/// DFI profile for all residues.
pub fn PerturbationResponse::dfi_profile(self : PerturbationResponse) -> Array[Double] {
  let result = Array::make(self.n, 0.0)
  for i = 0; i < self.n; i = i + 1 {
    result[i] = self.dfi(i)
  }
  result
}
