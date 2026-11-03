/// Scanpy-inspired: pp.combat — batch effect correction using ComBat algorithm
/// Reference: Johnson et al. (2007) Biostatistics

///|
/// Run ComBat batch correction on an AnnData-like matrix.
///
/// Parameters:
/// - `x`: cells × genes matrix (log-transformed)
/// - `batches`: batch label per cell
/// - `parametric`: use parametric (true) or non-parametric (false) estimation
pub fn scanpy_combat(
  x : Array[Array[Double]],
  batches : Array[String],
  parametric~ : Bool = true,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  // Identify unique batches
  let batch_set : Set[String] = Set::from_iter(batches.iter())
  let batch_labels = batch_set.to_array()
  let n_batches = batch_labels.length()
  if n_batches <= 1 {
    return x // no correction needed with single batch
  }
  // Map cells to batch indices
  let batch_idx = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let mut idx = 0
    for b in 0..<n_batches {
      if batches[i] == batch_labels[b] {
        idx = b
      }
    }
    batch_idx.push(idx)
  }
  // Compute grand mean per gene
  let grand_mean = Array(capacity=n_genes)
  for j in 0..<n_genes {
    let mut s = 0.0
    for i in 0..<n_cells {
      s += x[i][j]
    }
    grand_mean.push(s / n_cells.to_double())
  }
  // Compute batch means and variances per gene
  let batch_means = Array(capacity=n_genes)
  let batch_vars = Array(capacity=n_genes)
  let batch_sizes = Array(capacity=n_batches)
  for b in 0..<n_batches {
    let mut cnt = 0
    for i in 0..<n_cells {
      if batch_idx[i] == b {
        cnt += 1
      }
    }
    batch_sizes.push(cnt)
  }
  for j in 0..<n_genes {
    let means = Array(capacity=n_batches)
    let vars = Array(capacity=n_batches)
    for b in 0..<n_batches {
      let mut s = 0.0
      for i in 0..<n_cells {
        if batch_idx[i] == b {
          s += x[i][j]
        }
      }
      let mean = if batch_sizes[b] > 0 {
        s / batch_sizes[b].to_double()
      } else {
        0.0
      }
      means.push(mean)
      let mut v = 0.0
      for i in 0..<n_cells {
        if batch_idx[i] == b {
          let diff = x[i][j] - mean
          v += diff * diff
        }
      }
      if batch_sizes[b] > 1 {
        vars.push(v / (batch_sizes[b] - 1).to_double())
      } else {
        vars.push(1.0)
      }
    }
    batch_means.push(means)
    batch_vars.push(vars)
  }
  // Compute pooled variance (prior variance estimate)
  let pooled_var = Array(capacity=n_genes)
  for j in 0..<n_genes {
    let mut s = 0.0
    let mut total_dof = 0
    for b in 0..<n_batches {
      if batch_sizes[b] > 1 {
        s += batch_vars[j][b] * (batch_sizes[b] - 1).to_double()
        total_dof += batch_sizes[b] - 1
      }
    }
    if total_dof > 0 {
      pooled_var.push(s / total_dof.to_double())
    } else {
      pooled_var.push(1.0)
    }
  }
  // Apply ComBat correction
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let b = batch_idx[i]
    let row = Array(capacity=n_genes)
    for j in 0..<n_genes {
      let delta = batch_means[j][b] - grand_mean[j]
      let sigma_sq = if parametric { batch_vars[j][b] } else { pooled_var[j] }
      let scale = if sigma_sq > 0.0 {
        (pooled_var[j] / sigma_sq).sqrt()
      } else {
        1.0
      }
      let corrected = (x[i][j] - delta) * scale
      row.push(corrected)
    }
    result.push(row)
  }
  result
}
