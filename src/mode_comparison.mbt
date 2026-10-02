///|
/// Mode comparison and overlap analysis.
/// Compares normal modes from different models (e.g., ANM vs GNM,
/// or same model for different proteins).
/// Based on ProDy's `calcOverlap`, `cumOverlap`, `scalingFactor`.

///|
/// Compute overlap between two mode vectors (eigenvectors).
/// Overlap = |v1 · v2| (absolute dot product).
pub fn mode_overlap(vec1 : Array[Double], vec2 : Array[Double]) -> Double {
  let n = if vec1.length() < vec2.length() {
    vec1.length()
  } else {
    vec2.length()
  }
  let mut dot = 0.0
  for i = 0; i < n; i = i + 1 {
    dot = dot + vec1[i] * vec2[i]
  }
  dot.abs()
}

///|
/// Cumulative overlap of a mode from set 1 with the first k modes of set 2.
/// Returns sqrt(sum of squared overlaps).
pub fn cumulative_overlap(
  mode : Array[Double],
  modes_set : Array[Array[Double]],
  k : Int,
) -> Double {
  let max_k = if k < modes_set.length() { k } else { modes_set.length() }
  let mut sum_sq = 0.0
  for i = 0; i < max_k; i = i + 1 {
    let ov = mode_overlap(mode, modes_set[i])
    sum_sq = sum_sq + ov * ov
  }
  sum_sq.sqrt()
}

///|
/// Scaling factor between two mode vectors.
/// Finds the scalar s that minimizes ||v1 - s*v2||^2.
/// s = (v1 · v2) / (v2 · v2)
pub fn scaling_factor(vec1 : Array[Double], vec2 : Array[Double]) -> Double {
  let n = if vec1.length() < vec2.length() {
    vec1.length()
  } else {
    vec2.length()
  }
  let mut dot12 = 0.0
  let mut dot22 = 0.0
  for i = 0; i < n; i = i + 1 {
    dot12 = dot12 + vec1[i] * vec2[i]
    dot22 = dot22 + vec2[i] * vec2[i]
  }
  if dot22 < 1.0e-20 {
    return 0.0
  }
  dot12 / dot22
}

///|
/// Covariance overlap between two sets of modes.
/// Measures how similar two subspaces are, weighted by eigenvalues.
/// CO = sum_i sum_j (1/sqrt(lambda_i * mu_j)) * (v_i · w_j)^2
/// Normalized to [0, 1].
pub fn covariance_overlap(
  evals1 : Array[Double],
  evecs1 : Array[Array[Double]],
  evals2 : Array[Double],
  evecs2 : Array[Array[Double]],
) -> Double {
  let n1 = evals1.length()
  let n2 = evals2.length()
  let n_vec = if n1 > 0 && n2 > 0 {
    if evecs1[0].length() < evecs2[0].length() {
      evecs1[0].length()
    } else {
      evecs2[0].length()
    }
  } else {
    return 0.0
  }
  let mut total = 0.0
  let mut norm1 = 0.0
  let mut norm2 = 0.0
  for i = 0; i < n1; i = i + 1 {
    let li = if evals1[i] > 1.0e-10 { evals1[i] } else { 1.0e-10 }
    norm1 = norm1 + 1.0 / li
    for j = 0; j < n2; j = j + 1 {
      let lj = if evals2[j] > 1.0e-10 { evals2[j] } else { 1.0e-10 }
      let mut dot = 0.0
      for kk = 0; kk < n_vec; kk = kk + 1 {
        dot = dot + evecs1[i][kk] * evecs2[j][kk]
      }
      total = total + dot * dot / (li * lj).sqrt()
    }
  }
  for j = 0; j < n2; j = j + 1 {
    let lj = if evals2[j] > 1.0e-10 { evals2[j] } else { 1.0e-10 }
    norm2 = norm2 + 1.0 / lj
  }
  if norm1 < 1.0e-20 || norm2 < 1.0e-20 {
    return 0.0
  }
  total / (norm1 * norm2).sqrt()
}

///|
/// Angle between two subspaces (in degrees).
/// Computed from the singular values of the projection matrix V1^T V2.
/// Returns the largest principal angle.
pub fn subspace_angle(
  evecs1 : Array[Array[Double]],
  evecs2 : Array[Array[Double]],
) -> Double {
  let n1 = evecs1.length()
  let n2 = evecs2.length()
  if n1 == 0 || n2 == 0 {
    return 90.0
  }
  let n_vec = if evecs1[0].length() < evecs2[0].length() {
    evecs1[0].length()
  } else {
    evecs2[0].length()
  }
  // Build projection matrix P = V1^T * V2 (n1 x n2)
  let proj = Array::makei(n1, { _ => Array::make(n2, 0.0) })
  for i = 0; i < n1; i = i + 1 {
    for j = 0; j < n2; j = j + 1 {
      let mut s = 0.0
      for kk = 0; kk < n_vec; kk = kk + 1 {
        s = s + evecs1[i][kk] * evecs2[j][kk]
      }
      proj[i][j] = s
    }
  }
  // Compute P^T * P (n2 x n2) and find smallest eigenvalue
  // The smallest singular value = sqrt(smallest eigenvalue of P^T P)
  // Largest principal angle = acos(smallest singular value)
  let ptp = Array::makei(n2, { _ => Array::make(n2, 0.0) })
  for i = 0; i < n2; i = i + 1 {
    for j = 0; j < n2; j = j + 1 {
      let mut s = 0.0
      for kk = 0; kk < n1; kk = kk + 1 {
        s = s + proj[kk][i] * proj[kk][j]
      }
      ptp[i][j] = s
    }
  }
  // Eigen-decomposition of P^T P (symmetric, use Jacobi)
  let (eigenvalues, _) = linalg_jacobi(ptp, n2, 100)
  // Find smallest eigenvalue
  let mut min_ev = eigenvalues[0]
  for i = 1; i < n2; i = i + 1 {
    if eigenvalues[i] < min_ev {
      min_ev = eigenvalues[i]
    }
  }
  if min_ev < 0.0 {
    min_ev = 0.0
  }
  if min_ev > 1.0 {
    min_ev = 1.0
  }
  let cos_angle = min_ev.sqrt()
  // acos
  let angle_rad = @math.acos(cos_angle)
  angle_rad * 180.0 / 3.14159265358979323846
}

///|
/// Full overlap matrix between two sets of modes.
/// Returns n1 x n2 matrix of |v_i · w_j|.
pub fn overlap_matrix(
  modes1 : Array[Array[Double]],
  modes2 : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n1 = modes1.length()
  let n2 = modes2.length()
  let result = Array::makei(n1, { _ => Array::make(n2, 0.0) })
  for i = 0; i < n1; i = i + 1 {
    for j = 0; j < n2; j = j + 1 {
      result[i][j] = mode_overlap(modes1[i], modes2[j])
    }
  }
  result
}

///|
/// Effective overlap: a single number summarizing how well
/// the first k modes of set 2 explain mode i from set 1.
pub fn effective_overlap(
  mode : Array[Double],
  modes_set : Array[Array[Double]],
  k : Int,
) -> Double {
  cumulative_overlap(mode, modes_set, k)
}
