/// Scanpy-inspired: pp.normalize_pearson_residuals — analytic Pearson residual normalization
/// Reference: Lause et al. (2021) Genome Biology

///|
/// Compute analytic Pearson residuals for normalization.
///
/// The Pearson residual for cell i, gene j is:
///   r_ij = (x_ij - mu_ij) / sqrt(mu_ij + theta * mu_ij^2)
/// where mu_ij = s_i * g_j / (sum of all) and theta is the overdispersion.
///
/// Parameters:
/// - `x`: cells × genes count matrix
/// - `theta`: overdispersion parameter (default 100)
/// - `clip`: clip residuals to ±clip (default 30; 0 = no clip)
pub fn scanpy_normalize_pearson_residuals(
  x : Array[Array[Double]],
  theta~ : Double = 100.0,
  clip~ : Double = 30.0,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  // Compute row sums (cell sizes) and col sums (gene counts)
  let row_sums = Array(capacity=n_cells)
  let col_sums = Array(capacity=n_genes)
  for _ in 0..<n_genes {
    col_sums.push(0.0)
  }
  let mut total = 0.0
  for i in 0..<n_cells {
    let mut s = 0.0
    for j in 0..<n_genes {
      s += x[i][j]
      col_sums[j] += x[i][j]
    }
    row_sums.push(s)
    total += s
  }
  if total <= 0.0 {
    return x
  }
  // Compute residuals
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=n_genes)
    for j in 0..<n_genes {
      // mu_ij = row_sum_i * col_sum_j / total
      let mu = row_sums[i] * col_sums[j] / total
      if mu > 0.0 {
        let variance = mu + theta * mu * mu
        let std = variance.sqrt()
        let residual = (x[i][j] - mu) / std
        if clip > 0.0 {
          if residual > clip {
            row.push(clip)
          } else if residual < -clip {
            row.push(-clip)
          } else {
            row.push(residual)
          }
        } else {
          row.push(residual)
        }
      } else {
        row.push(0.0)
      }
    }
    result.push(row)
  }
  result
}

///|
/// Analyze Pearson residuals: compute the optimal theta by fitting.
/// Uses a simple grid search to find theta that minimizes residual variance.
///
/// Parameters:
/// - `x`: cells × genes count matrix
/// - `min_theta`: minimum theta to search (default 0.01)
/// - `max_theta`: maximum theta (default 1000)
/// - `n_grid`: number of grid points (default 50)
pub fn scanpy_analyze_pearson_residuals(
  x : Array[Array[Double]],
  min_theta~ : Double = 0.01,
  max_theta~ : Double = 1000.0,
  n_grid~ : Int32 = 50,
) -> (Double, Array[Double]) {
  let mut best_theta = min_theta
  let mut best_var = 1.0e100
  let variances = Array(capacity=n_grid)
  for g in 0..<n_grid {
    let t = min_theta *
      @math.pow(max_theta / min_theta, g.to_double() / (n_grid - 1).to_double())
    let residuals = scanpy_normalize_pearson_residuals(x, theta=t, clip=0.0)
    // Compute mean variance of residuals
    let n_cells = residuals.length()
    if n_cells == 0 {
      variances.push(1.0e100)
      continue
    }
    let n_genes = residuals[0].length()
    let mut total_var = 0.0
    for j in 0..<n_genes {
      let mut s = 0.0
      for i in 0..<n_cells {
        s += residuals[i][j]
      }
      let mean = s / n_cells.to_double()
      let mut v = 0.0
      for i in 0..<n_cells {
        let diff = residuals[i][j] - mean
        v += diff * diff
      }
      total_var += v / n_cells.to_double()
    }
    let mean_var = if n_genes > 0 {
      total_var / n_genes.to_double()
    } else {
      1.0e100
    }
    variances.push(mean_var)
    // Prefer theta where residual variance is closest to 1
    let dev = (mean_var - 1.0) * (mean_var - 1.0)
    if dev < best_var {
      best_var = dev
      best_theta = t
    }
  }
  (best_theta, variances)
}
