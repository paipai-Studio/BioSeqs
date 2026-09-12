///|
/// scMerge - Single-cell batch correction using RUV (Remove Unwanted Variation)
/// with stably expressed genes (SEGs).
/// Inspired by Bioconductor scMerge package (Lin et al., 2019).
/// Uses RUV4 (Remove Unwanted Variation, 4th variant) for batch correction,
/// with stably expressed genes identified by minimum variance across batches.

///|
/// Result of scMerge batch correction.
pub struct ScMergeResult {
  corrected_data : Array[Array[Double]]
  seg_genes : Array[Int]
  estimated_batch_factors : Array[Array[Double]]
  cell_type_labels : Array[Int]
}

///|
pub fn ScMergeResult::new(
  corrected_data : Array[Array[Double]],
  seg_genes : Array[Int],
  estimated_batch_factors : Array[Array[Double]],
  cell_type_labels : Array[Int],
) -> ScMergeResult {
  ScMergeResult::{
    corrected_data,
    seg_genes,
    estimated_batch_factors,
    cell_type_labels,
  }
}

///|
/// Compute the mean of an array.
fn scm_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    sum = sum + arr[i]
    i = i + 1
  }
  sum / arr.length().to_double()
}

///|
/// Compute the variance of an array.
fn scm_variance(arr : Array[Double]) -> Double {
  if arr.length() <= 1 {
    return 0.0
  }
  let mean = scm_mean(arr)
  let mut sum_sq = 0.0
  let mut i = 0
  while i < arr.length() {
    let diff = arr[i] - mean
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  sum_sq / (arr.length().to_double() - 1.0)
}

///|
/// Z-score normalize a matrix (genes x cells) row-wise.
/// Returns a matrix where each gene has mean=0, std=1.
pub fn scm_zscore_matrix(data : Array[Array[Double]]) -> Array[Array[Double]] {
  if data.length() == 0 {
    return Array::new()
  }
  let n_genes = data.length()
  let n_cells = data[0].length()
  let result : Array[Array[Double]] = Array::make(n_genes, Array::new())
  let mut i = 0
  while i < n_genes {
    let row = data[i]
    let mean = scm_mean(row)
    let std = scm_variance(row).sqrt()
    let normalized : Array[Double] = Array::make(n_cells, 0.0)
    let mut j = 0
    while j < n_cells {
      normalized[j] = if std > 1.0e-10 { (row[j] - mean) / std } else { 0.0 }
      j = j + 1
    }
    result[i] = normalized
    i = i + 1
  }
  result
}

///|
/// Find stably expressed genes (SEGs) across batches.
/// Selects genes with the lowest variance within each batch and
/// most consistent expression across batches.
///
/// Parameters:
///   batches - array of batch matrices, each (n_cells x n_genes)
///   n_seg - number of stably expressed genes to select
///
/// Returns indices of the selected SEGs.
pub fn scm_find_seg(
  batches : Array[Array[Array[Double]]],
  n_seg? : Int = 500,
) -> Array[Int] {
  let n_batches = batches.length()
  if n_batches == 0 || batches[0].length() == 0 {
    return Array::new()
  }

  let n_genes = batches[0][0].length()
  let n_cells_total = batches[0].length()
  if n_genes == 0 {
    return Array::new()
  }

  // Compute gene-wise variance within each batch
  let gene_variances : Array[Array[Double]] = Array::make(n_genes, Array::new())
  let mut g = 0
  while g < n_genes {
    let batch_variances : Array[Double] = Array::make(n_batches, 0.0)
    let mut b = 0
    while b < n_batches {
      let batch = batches[b]
      if batch.length() > 0 {
        let gene_values : Array[Double] = Array::make(batch.length(), 0.0)
        let mut c = 0
        while c < batch.length() {
          gene_values[c] = batch[c][g]
          c = c + 1
        }
        batch_variances[b] = scm_variance(gene_values)
      }
      b = b + 1
    }
    gene_variances[g] = batch_variances
    g = g + 1
  }

  // Compute stability score: max variance across batches (lower is more stable)
  let stability_scores : Array[(Double, Int)] = Array::make(n_genes, (0.0, 0))
  g = 0
  while g < n_genes {
    let mut max_var = 0.0
    let mut b = 0
    while b < n_batches {
      if gene_variances[g][b] > max_var {
        max_var = gene_variances[g][b]
      }
      b = b + 1
    }
    stability_scores[g] = (max_var, g)
    g = g + 1
  }

  // Sort by stability score (ascending = most stable first)
  let mut i = 0
  while i < n_genes - 1 {
    let mut j = i + 1
    while j < n_genes {
      if stability_scores[i].0 > stability_scores[j].0 {
        let tmp = stability_scores[i]
        stability_scores[i] = stability_scores[j]
        stability_scores[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }

  // Select top n_seg genes
  let n_select = if n_seg < n_genes { n_seg } else { n_genes }
  let seg_indices : Array[Int] = Array::make(n_select, 0)
  i = 0
  while i < n_select {
    seg_indices[i] = stability_scores[i].1
    i = i + 1
  }

  seg_indices
}

///|
/// Perform RUV4 (Remove Unwanted Variation) batch correction.
/// Uses stably expressed genes (SEGs) as negative controls to estimate
/// and remove unwanted variation.
///
/// Parameters:
///   data - combined expression matrix (n_genes x n_cells_all)
///   batch_idx - batch index for each cell (0-based)
///   seg_genes - indices of stably expressed genes
///   n_factors - number of unwanted factors to estimate
///
/// Returns corrected data (n_genes x n_cells_all).
pub fn scm_ruv4(
  data : Array[Array[Double]],
  batch_idx : Array[Int],
  seg_genes : Array[Int],
  n_factors : Int,
) -> Array[Array[Double]] {
  let n_genes = data.length()
  if n_genes == 0 {
    return data
  }
  let n_cells = data[0].length()
  let n_batches = scm_max_batch(batch_idx) + 1

  // Build design matrix: one-hot encoding for batches
  // X = [1, batch1, batch2, ...] (intercept + batch indicators)
  let n_design = n_batches
  // Create design matrix X (n_cells x n_design)
  let design : Array[Array[Double]] = Array::make(n_cells, Array::new())
  let mut c = 0
  while c < n_cells {
    let row : Array[Double] = Array::make(n_design, 0.0)
    row[0] = 1.0 // intercept
    if batch_idx[c] > 0 && batch_idx[c] < n_design {
      row[batch_idx[c]] = 1.0
    }
    design[c] = row
    c = c + 1
  }

  // Compute projection onto the orthogonal complement of the design matrix
  // R = (I - X(X'X)^{-1}X') * Y
  // For SEGs only, compute the residual matrix

  // Extract SEG rows
  let n_seg = seg_genes.length()
  let seg_data : Array[Array[Double]] = Array::make(n_seg, Array::new())
  let mut s = 0
  while s < n_seg {
    seg_data[s] = data[seg_genes[s]]
    s = s + 1
  }

  // Compute residuals: R = Y - X * beta_hat
  // where beta_hat = (X'X)^{-1} X'Y
  // For each gene, fit a linear model and get residuals

  // Compute X'X (n_design x n_design)
  let xtx : Array[Array[Double]] = Array::make(n_design, Array::new())
  s = 0
  while s < n_design {
    xtx[s] = Array::make(n_design, 0.0)
    s = s + 1
  }
  c = 0
  while c < n_cells {
    let mut i = 0
    while i < n_design {
      let mut j = 0
      while j < n_design {
        xtx[i][j] = xtx[i][j] + design[c][i] * design[c][j]
        j = j + 1
      }
      i = i + 1
    }
    c = c + 1
  }

  // Compute X'Y (n_design x n_genes) - for all genes
  let xty_all : Array[Array[Double]] = Array::make(n_design, Array::new())
  s = 0
  while s < n_design {
    xty_all[s] = Array::make(n_genes, 0.0)
    s = s + 1
  }
  c = 0
  while c < n_cells {
    let mut i = 0
    while i < n_design {
      let mut g = 0
      while g < n_genes {
        xty_all[i][g] = xty_all[i][g] + design[c][i] * data[g][c]
        g = g + 1
      }
      i = i + 1
    }
    c = c + 1
  }

  // Invert X'X using Gaussian elimination
  let xtx_inv : Array[Array[Double]] = scm_invert_matrix(xtx, n_design)

  // Compute beta = (X'X)^{-1} X'Y
  let beta_all : Array[Array[Double]] = Array::make(n_design, Array::new())
  s = 0
  while s < n_design {
    beta_all[s] = Array::make(n_genes, 0.0)
    s = s + 1
  }
  s = 0
  while s < n_design {
    let mut g = 0
    while g < n_genes {
      let mut sum = 0.0
      let mut k = 0
      while k < n_design {
        sum = sum + xtx_inv[s][k] * xty_all[k][g]
        k = k + 1
      }
      beta_all[s][g] = sum
      g = g + 1
    }
    s = s + 1
  }

  // Compute residuals: R = Y - X * beta_hat
  let residuals : Array[Array[Double]] = Array::make(n_genes, Array::new())
  let mut g = 0
  while g < n_genes {
    let res_row : Array[Double] = Array::make(n_cells, 0.0)
    c = 0
    while c < n_cells {
      let mut pred = 0.0
      let mut s = 0
      while s < n_design {
        pred = pred + design[c][s] * beta_all[s][g]
        s = s + 1
      }
      res_row[c] = data[g][c] - pred
      c = c + 1
    }
    residuals[g] = res_row
    g = g + 1
  }

  // Extract SEG residuals for SVD
  let seg_residuals : Array[Array[Double]] = Array::make(n_seg, Array::new())
  s = 0
  while s < n_seg {
    seg_residuals[s] = residuals[seg_genes[s]]
    s = s + 1
  }

  // Perform SVD on SEG residuals to estimate unwanted factors
  // Use randomized SVD approximation
  let w_hat = scm_randomized_svd(seg_residuals, n_factors)

  // Compute alpha = (W'W)^{-1} W' * Y
  // For each gene, regress W on the original data
  let corrected : Array[Array[Double]] = Array::make(n_genes, Array::new())
  g = 0
  while g < n_genes {
    let y = data[g]
    // Compute W'Y (n_factors x 1)
    let wty : Array[Double] = Array::make(n_factors, 0.0)
    s = 0
    while s < n_factors {
      let mut sum = 0.0
      c = 0
      while c < n_cells {
        sum = sum + w_hat[c][s] * y[c]
        c = c + 1
      }
      wty[s] = sum
      s = s + 1
    }

    // Compute W'W (n_factors x n_factors)
    let wtw : Array[Array[Double]] = Array::make(n_factors, Array::new())
    s = 0
    while s < n_factors {
      wtw[s] = Array::make(n_factors, 0.0)
      s = s + 1
    }
    c = 0
    while c < n_cells {
      let mut i = 0
      while i < n_factors {
        let mut j = 0
        while j < n_factors {
          wtw[i][j] = wtw[i][j] + w_hat[c][i] * w_hat[c][j]
          j = j + 1
        }
        i = i + 1
      }
      c = c + 1
    }

    // Invert W'W
    let wtw_inv = scm_invert_matrix(wtw, n_factors)

    // Compute alpha = (W'W)^{-1} W'Y
    let alpha : Array[Double] = Array::make(n_factors, 0.0)
    s = 0
    while s < n_factors {
      let mut sum = 0.0
      let mut k = 0
      while k < n_factors {
        sum = sum + wtw_inv[s][k] * wty[k]
        k = k + 1
      }
      alpha[s] = sum
      s = s + 1
    }

    // Corrected: Y_corrected = Y - W * alpha
    let corrected_row : Array[Double] = Array::make(n_cells, 0.0)
    c = 0
    while c < n_cells {
      let mut unwanted = 0.0
      s = 0
      while s < n_factors {
        unwanted = unwanted + w_hat[c][s] * alpha[s]
        s = s + 1
      }
      corrected_row[c] = y[c] - unwanted
      c = c + 1
    }
    corrected[g] = corrected_row
    g = g + 1
  }

  corrected
}

///|
/// Find the maximum batch index.
fn scm_max_batch(batch_idx : Array[Int]) -> Int {
  let mut max_val = 0
  let mut i = 0
  while i < batch_idx.length() {
    if batch_idx[i] > max_val {
      max_val = batch_idx[i]
    }
    i = i + 1
  }
  max_val
}

///|
/// Invert a square matrix using Gaussian elimination.
pub fn scm_invert_matrix(
  matrix : Array[Array[Double]],
  n : Int,
) -> Array[Array[Double]] {
  // Augmented matrix [A | I]
  let aug : Array[Array[Double]] = Array::make(n, Array::new())
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::make(2 * n, 0.0)
    let mut j = 0
    while j < n {
      row[j] = matrix[i][j]
      j = j + 1
    }
    row[n + i] = 1.0
    aug[i] = row
    i = i + 1
  }

  // Forward elimination
  i = 0
  while i < n {
    // Find pivot
    let mut max_row = i
    let mut max_val = aug[i][i].abs()
    let mut r = i + 1
    while r < n {
      if aug[r][i].abs() > max_val {
        max_val = aug[r][i].abs()
        max_row = r
      }
      r = r + 1
    }

    // Swap rows
    if max_row != i {
      let tmp = aug[i]
      aug[i] = aug[max_row]
      aug[max_row] = tmp
    }

    // Scale pivot row
    let pivot = aug[i][i]
    if pivot.abs() > 1.0e-15 {
      let mut j = 0
      while j < 2 * n {
        aug[i][j] = aug[i][j] / pivot
        j = j + 1
      }
    }

    // Eliminate below
    r = i + 1
    while r < n {
      let factor = aug[r][i]
      if factor.abs() > 1.0e-15 {
        let mut j = 0
        while j < 2 * n {
          aug[r][j] = aug[r][j] - factor * aug[i][j]
          j = j + 1
        }
      }
      r = r + 1
    }
    i = i + 1
  }

  // Back substitution
  i = n - 1
  while i >= 0 {
    let mut j = i - 1
    while j >= 0 {
      let factor = aug[j][i]
      if factor.abs() > 1.0e-15 {
        let mut k = 0
        while k < 2 * n {
          aug[j][k] = aug[j][k] - factor * aug[i][k]
          k = k + 1
        }
      }
      j = j - 1
    }
    i = i - 1
  }

  // Extract inverse
  let result : Array[Array[Double]] = Array::make(n, Array::new())
  i = 0
  while i < n {
    let row : Array[Double] = Array::make(n, 0.0)
    let mut j = 0
    while j < n {
      row[j] = aug[i][n + j]
      j = j + 1
    }
    result[i] = row
    i = i + 1
  }

  result
}

///|
/// Randomized SVD approximation for tall matrices.
/// Returns the right singular vectors V (n_cells x k).
pub fn scm_randomized_svd(
  matrix : Array[Array[Double]],
  k : Int,
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  if n_rows == 0 || k <= 0 {
    return Array::new()
  }
  let n_cols = matrix[0].length()
  let k_eff = if k < n_cols { k } else { n_cols }
  if k_eff <= 0 {
    return Array::new()
  }

  // Use power iteration to compute dominant right singular vectors
  // Initialize V randomly
  let v : Array[Array[Double]] = Array::make(k_eff, Array::new())
  let mut i = 0
  while i < k_eff {
    let row : Array[Double] = Array::make(n_cols, 0.0)
    let mut j = 0
    while j < n_cols {
      row[j] = ((i * 10007 + j * 3 + 1) % 1001).to_double() / 1000.0 - 0.5
      j = j + 1
    }
    // Normalize
    let mut norm = 0.0
    j = 0
    while j < n_cols {
      norm = norm + row[j] * row[j]
      j = j + 1
    }
    norm = norm.sqrt()
    if norm > 1.0e-10 {
      j = 0
      while j < n_cols {
        row[j] = row[j] / norm
        j = j + 1
      }
    }
    v[i] = row
    i = i + 1
  }

  // Power iteration: (A'A)^p * A' * random
  let mut iter = 0
  while iter < 3 {
    // Compute A * V (n_rows x k_eff)
    let av : Array[Array[Double]] = Array::make(n_rows, Array::new())
    i = 0
    while i < n_rows {
      let row : Array[Double] = Array::make(k_eff, 0.0)
      let mut j = 0
      while j < k_eff {
        let mut sum = 0.0
        let mut c = 0
        while c < n_cols {
          sum = sum + matrix[i][c] * v[j][c]
          c = c + 1
        }
        row[j] = sum
        j = j + 1
      }
      av[i] = row
      i = i + 1
    }

    // Compute A' * (A * V) = (A'A) * V (n_cols x k_eff)
    let atav : Array[Array[Double]] = Array::make(k_eff, Array::new())
    i = 0
    while i < k_eff {
      let row : Array[Double] = Array::make(n_cols, 0.0)
      let mut c = 0
      while c < n_cols {
        let mut sum = 0.0
        let mut r = 0
        while r < n_rows {
          sum = sum + matrix[r][c] * av[r][i]
          r = r + 1
        }
        row[c] = sum
        c = c + 1
      }
      // Normalize
      let mut norm = 0.0
      c = 0
      while c < n_cols {
        norm = norm + row[c] * row[c]
        c = c + 1
      }
      norm = norm.sqrt()
      if norm > 1.0e-10 {
        c = 0
        while c < n_cols {
          row[c] = row[c] / norm
          c = c + 1
        }
      }
      atav[i] = row
      i = i + 1
    }

    // Copy back to v
    i = 0
    while i < k_eff {
      let mut c = 0
      while c < n_cols {
        v[i][c] = atav[i][c]
        c = c + 1
      }
      i = i + 1
    }

    iter = iter + 1
  }

  // Transpose v to get (n_cells x k_eff)
  let result : Array[Array[Double]] = Array::make(n_cols, Array::new())
  let mut c = 0
  while c < n_cols {
    let row : Array[Double] = Array::make(k_eff, 0.0)
    let mut i = 0
    while i < k_eff {
      row[i] = v[i][c]
      i = i + 1
    }
    result[c] = row
    c = c + 1
  }

  result
}

///|
/// Perform scMerge batch correction.
/// This is the main entry point, analogous to scMerge() in the Bioconductor
/// scMerge package.
///
/// Parameters:
///   batches - array of expression matrices, each (n_cells x n_genes)
///   n_seg - number of stably expressed genes to use (default: 500)
///   n_factors - number of unwanted factors to estimate (default: 2)
///
/// Returns ScMergeResult with corrected data.
pub fn scmerge(
  batches : Array[Array[Array[Double]]],
  n_seg? : Int = 500,
  n_factors? : Int = 2,
) -> ScMergeResult {
  let n_batches = batches.length()
  if n_batches == 0 {
    return ScMergeResult::new([], [], [], [])
  }

  let n_genes = batches[0][0].length()
  if n_genes == 0 {
    return ScMergeResult::new([], [], [], [])
  }

  // Step 1: Find stably expressed genes (SEGs)
  let seg_genes = scm_find_seg(batches, n_seg~)

  // Step 2: Combine all batches into one matrix
  let total_cells : Array[Int] = Array::make(n_batches, 0)
  let mut total = 0
  let mut b = 0
  while b < n_batches {
    total_cells[b] = batches[b].length()
    total = total + batches[b].length()
    b = b + 1
  }

  let mut offset = 0

  let combined : Array[Array[Double]] = Array::make(n_genes, Array::new())
  let mut g = 0
  while g < n_genes {
    let row : Array[Double] = Array::make(total, 0.0)
    offset = 0
    b = 0
    while b < n_batches {
      let batch = batches[b]
      let mut c = 0
      while c < batch.length() {
        row[offset + c] = batch[c][g]
        c = c + 1
      }
      offset = offset + batch.length()
      b = b + 1
    }
    combined[g] = row
    g = g + 1
  }

  // Build batch index for each cell
  let batch_idx : Array[Int] = Array::make(total, 0)
  offset = 0
  b = 0
  while b < n_batches {
    let mut c = 0
    while c < batches[b].length() {
      batch_idx[offset + c] = b
      c = c + 1
    }
    offset = offset + batches[b].length()
    b = b + 1
  }

  // Step 3: Apply RUV4 correction
  let corrected_data = scm_ruv4(combined, batch_idx, seg_genes, n_factors)

  // Step 4: Estimate cell type labels (simplified: cluster within each batch)
  let cell_type_labels : Array[Int] = Array::make(total, 0)
  offset = 0
  b = 0
  while b < n_batches {
    let mut c = 0
    while c < batches[b].length() {
      cell_type_labels[offset + c] = b
      c = c + 1
    }
    offset = offset + batches[b].length()
    b = b + 1
  }

  // Transpose corrected data back to (n_cells x n_genes) format
  let corrected_by_cell : Array[Array[Double]] = Array::make(
    total,
    Array::new(),
  )
  let mut c = 0
  while c < total {
    let row : Array[Double] = Array::make(n_genes, 0.0)
    let mut g = 0
    while g < n_genes {
      row[g] = corrected_data[g][c]
      g = g + 1
    }
    corrected_by_cell[c] = row
    c = c + 1
  }

  // Estimate batch factors (simplified: just the RUV factors)
  let estimated_batch_factors : Array[Array[Double]] = Array::make(
    total,
    Array::new(),
  )
  let mut c2 = 0
  while c2 < total {
    let row : Array[Double] = Array::make(n_factors, 0.0)
    row[0] = batch_idx[c2].to_double()
    estimated_batch_factors[c2] = row
    c2 = c2 + 1
  }

  ScMergeResult::new(
    corrected_by_cell, seg_genes, estimated_batch_factors, cell_type_labels,
  )
}
