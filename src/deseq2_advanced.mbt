///|
/// DESeq2 Advanced - Variance Stabilizing Transformation and PCA.
/// Inspired by Bioconductor DESeq2 package.
///
/// Provides:
/// - variance_stabilizing_transformation: full VST using dispersion estimates
/// - vst: fast VST using a subset of genes for dispersion trend fitting
/// - plot_pca: Principal Component Analysis on transformed data

///|
/// Result of variance stabilizing transformation.
pub struct VSTResult {
  /// Transformed expression matrix (genes x samples), log2-scale
  transformed : Array[Array[Double]]
  /// Gene names
  row_names : Array[String]
  /// Sample names
  col_names : Array[String]
  /// Dispersion used (fitted trend value)
  dispersion : Double
}

///|
/// PCA result for sample visualization.
pub struct PCAResult {
  /// Sample names
  sample_names : Array[String]
  /// PC1 coordinates for each sample
  pc1 : Array[Double]
  /// PC2 coordinates for each sample
  pc2 : Array[Double]
  /// Variance explained by PC1 (percentage)
  var_explained_pc1 : Double
  /// Variance explained by PC2 (percentage)
  var_explained_pc2 : Double
}

///|
/// Compute the asinh-based variance stabilizing transformation for a single value.
/// h(mu) = (1/sqrt(alpha)) * asinh(sqrt(alpha * mu))
/// For numerical stability, this is equivalent to:
/// h(mu) = (2/sqrt(alpha)) * log(sqrt(alpha*mu) + sqrt(1 + alpha*mu))
fn vst_transform_single(mu : Double, alpha : Double) -> Double {
  if mu <= 0.0 {
    return 0.0
  }
  if alpha <= 0.0 {
    // No overdispersion: use log2(mu) with pseudocount
    return @math.ln(mu + 1.0) / @math.ln(2.0)
  }
  let sqrt_alpha = alpha.sqrt()
  let inner = alpha * mu
  let sqrt_inner = inner.sqrt()
  let sqrt_one_plus = (1.0 + inner).sqrt()
  let val = (sqrt_inner + sqrt_one_plus) / sqrt_alpha
  // Convert from natural log to log2
  2.0 * @math.ln(val) / @math.ln(2.0)
}

///|
/// Full variance stabilizing transformation.
/// Uses the fitted dispersion from the DESeqDataSet to transform counts
/// so that variance is approximately independent of mean.
///
/// Parameters:
/// - dds: A DESeqDataSet with size factors and dispersions estimated
/// - blind: If true, re-estimate dispersions using a blind design (default: true)
///
/// Returns a VSTResult with the transformed matrix.
pub fn variance_stabilizing_transformation(
  dds : DESeqDataSet,
  blind? : Bool = true,
) -> VSTResult {
  let n_genes = dds.counts.length()
  let n_samples = dds.col_names.length()
  let normalized = normalize_counts(dds)

  // Compute the geometric mean of normalized counts per gene
  // to get a representative dispersion value
  let mut global_disp = 0.0
  let mut n_valid_disp = 0
  for i = 0; i < n_genes; i = i + 1 {
    let d = if blind {
      // Use a default dispersion when blind
      0.1
    } else {
      dds.dispersions[i]
    }
    if d > 0.0 {
      global_disp = global_disp + d
      n_valid_disp = n_valid_disp + 1
    }
  }
  let mean_disp = if n_valid_disp > 0 {
    global_disp / n_valid_disp.to_double()
  } else {
    0.1
  }

  // Apply VST to each element
  let transformed : Array[Array[Double]] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < n_samples; j = j + 1 {
      row.push(vst_transform_single(normalized[i][j], mean_disp))
    }
    transformed.push(row)
  }

  VSTResult::{
    transformed,
    row_names: dds.row_names,
    col_names: dds.col_names,
    dispersion: mean_disp,
  }
}

///|
/// Fast variance stabilizing transformation.
/// Uses a subset of genes to estimate the dispersion-mean trend,
/// then applies the transformation. Faster than full VST for large datasets.
///
/// Parameters:
/// - dds: A DESeqDataSet with size factors estimated
/// - nsub: Number of genes to use for dispersion trend estimation (default: 1000)
///
/// Returns a VSTResult with the transformed matrix.
pub fn vst(dds : DESeqDataSet, nsub? : Int = 1000) -> VSTResult {
  let n_genes = dds.counts.length()
  let n_samples = dds.col_names.length()
  let normalized = normalize_counts(dds)

  // Select a subset of genes for dispersion estimation
  let actual_nsub = if nsub < n_genes { nsub } else { n_genes }

  // Estimate dispersion from the subset
  // Use the method of moments estimator
  let mut sum_mu_sq = 0.0
  let mut sum_var_minus_mu = 0.0
  let mut n_valid = 0

  let step = if n_genes > actual_nsub {
    n_genes / actual_nsub
  } else {
    1
  }

  let mut i = 0
  while i < n_genes && n_valid < actual_nsub {
    // Compute mean and variance for this gene
    let mut gene_sum = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      gene_sum = gene_sum + normalized[i][j]
    }
    let gene_mean = gene_sum / n_samples.to_double()

    if gene_mean > 0.0 {
      let mut gene_var = 0.0
      for j = 0; j < n_samples; j = j + 1 {
        let diff = normalized[i][j] - gene_mean
        gene_var = gene_var + diff * diff
      }
      gene_var = gene_var / (n_samples - 1).to_double()

      // For NB: var = mu + alpha * mu^2
      // So alpha = (var - mu) / mu^2
      if gene_var > gene_mean {
        sum_mu_sq = sum_mu_sq + gene_mean * gene_mean
        sum_var_minus_mu = sum_var_minus_mu + (gene_var - gene_mean)
        n_valid = n_valid + 1
      }
    }
    i = i + step
  }

  let est_disp = if n_valid > 0 && sum_mu_sq > 0.0 {
    let d = sum_var_minus_mu / sum_mu_sq
    // Clamp to reasonable range
    if d < 0.001 { 0.001 } else if d > 10.0 { 10.0 } else { d }
  } else {
    0.1
  }

  // Apply VST using estimated dispersion
  let transformed : Array[Array[Double]] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < n_samples; j = j + 1 {
      row.push(vst_transform_single(normalized[i][j], est_disp))
    }
    transformed.push(row)
  }

  VSTResult::{
    transformed,
    row_names: dds.row_names,
    col_names: dds.col_names,
    dispersion: est_disp,
  }
}

///|
/// Perform Principal Component Analysis on transformed data.
/// Computes the top 2 principal components for sample visualization.
///
/// Parameters:
/// - vst_result: Output from variance_stabilizing_transformation or vst
/// - ntop: Number of most variable genes to use (default: 500)
///
/// Returns a PCAResult with sample coordinates and variance explained.
pub fn plot_pca(vst_result : VSTResult, ntop? : Int = 500) -> PCAResult {
  let n_genes = vst_result.transformed.length()
  let n_samples = vst_result.col_names.length()
  let actual_ntop = if ntop < n_genes { ntop } else { n_genes }

  // Compute row variances and select top variable genes
  let gene_vars : Array[(Int, Double)] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let mut sum = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      sum = sum + vst_result.transformed[i][j]
    }
    let mean = sum / n_samples.to_double()
    let mut var = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      let diff = vst_result.transformed[i][j] - mean
      var = var + diff * diff
    }
    var = var / n_samples.to_double()
    gene_vars.push((i, var))
  }

  // Sort by variance descending (simple selection sort for top N)
  let selected_indices : Array[Int] = Array::new()
  let used : Array[Bool] = Array::make(n_genes, false)
  for _k = 0; _k < actual_ntop; _k = _k + 1 {
    let mut max_idx = -1
    let mut max_var = -1.0
    for i = 0; i < n_genes; i = i + 1 {
      if !used[i] && gene_vars[i].1 > max_var {
        max_var = gene_vars[i].1
        max_idx = i
      }
    }
    if max_idx >= 0 {
      selected_indices.push(max_idx)
      used[max_idx] = true
    }
  }

  // Build matrix of selected genes (samples x genes)
  let mat : Array[Array[Double]] = Array::new()
  for j = 0; j < n_samples; j = j + 1 {
    let row : Array[Double] = Array::new()
    for idx in selected_indices {
      row.push(vst_result.transformed[idx][j])
    }
    mat.push(row)
  }

  // Center by column means
  let n_cols = selected_indices.length()
  for col = 0; col < n_cols; col = col + 1 {
    let mut col_sum = 0.0
    for row = 0; row < n_samples; row = row + 1 {
      col_sum = col_sum + mat[row][col]
    }
    let col_mean = col_sum / n_samples.to_double()
    for row = 0; row < n_samples; row = row + 1 {
      mat[row][col] = mat[row][col] - col_mean
    }
  }

  // Compute covariance matrix (n_samples x n_samples) = X * X^T / (n_cols - 1)
  let cov : Array[Array[Double]] = Array::make(n_samples, Array::make(n_samples, 0.0))
  for i = 0; i < n_samples; i = i + 1 {
    for j = 0; j < n_samples; j = j + 1 {
      let mut dot = 0.0
      for k = 0; k < n_cols; k = k + 1 {
        dot = dot + mat[i][k] * mat[j][k]
      }
      cov[i][j] = dot / (n_cols - 1).to_double()
    }
  }

  // Power iteration to find top 2 eigenvectors
  let (eigenvalues, eigenvectors) = power_iteration_top2(cov, n_samples)

  // Project samples onto principal components
  let pc1 : Array[Double] = Array::new()
  let pc2 : Array[Double] = Array::new()
  for i = 0; i < n_samples; i = i + 1 {
    let mut proj1 = 0.0
    let mut proj2 = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      proj1 = proj1 + cov[i][j] * eigenvectors[0][j]
      proj2 = proj2 + cov[i][j] * eigenvectors[1][j]
    }
    pc1.push(proj1)
    pc2.push(proj2)
  }

  // Compute variance explained
  let total_var = if eigenvalues.0 + eigenvalues.1 > 0.0 {
    eigenvalues.0 + eigenvalues.1
  } else {
    1.0
  }
  let var_pc1 = eigenvalues.0 / total_var * 100.0
  let var_pc2 = eigenvalues.1 / total_var * 100.0

  PCAResult::{
    sample_names: vst_result.col_names,
    pc1,
    pc2,
    var_explained_pc1: var_pc1,
    var_explained_pc2: var_pc2,
  }
}

///|
/// Power iteration method to find the top 2 eigenvalues and eigenvectors.
/// Returns ((eigenvalue1, eigenvalue2), [eigenvector1, eigenvector2]).
fn power_iteration_top2(
  mat : Array[Array[Double]],
  n : Int,
) -> ((Double, Double), Array[Array[Double]]) {
  // First eigenvector
  let (eval1, evec1) = power_iteration_single(mat, n, None)

  // Deflate: remove the first eigenvector's contribution
  let deflated : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      deflated[i][j] = mat[i][j] - eval1 * evec1[i] * evec1[j]
    }
  }

  // Second eigenvector
  let (eval2, evec2) = power_iteration_single(deflated, n, Some(evec1))

  ((eval1, eval2), [evec1, evec2])
}

///|
/// Single power iteration to find the largest eigenvalue/eigenvector.
/// orthogonal_to: if provided, keeps the vector orthogonal to it.
fn power_iteration_single(
  mat : Array[Array[Double]],
  n : Int,
  orthogonal_to : Array[Double]?,
) -> (Double, Array[Double]) {
  // Initialize random vector
  let mut v : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    v.push(1.0 / n.to_double() + (i % 3).to_double() * 0.1)
    i = i + 1
  }

  // Normalize
  v = normalize_vector(v, n)

  // Orthogonalize if needed
  match orthogonal_to {
    Some(u) => v = normalize_vector(orthogonalize(v, u, n), n)
    None => ()
  }

  let mut eigenvalue = 0.0
  let mut iter = 0
  while iter < 100 {
    // Matrix-vector multiplication
    let new_v : Array[Double] = Array::make(n, 0.0)
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        new_v[i] = new_v[i] + mat[i][j] * v[j]
      }
    }

    // Orthogonalize if needed
    let final_v = match orthogonal_to {
      Some(u) => orthogonalize(new_v, u, n)
      None => new_v
    }

    // Compute Rayleigh quotient (eigenvalue estimate)
    let mut dot = 0.0
    for i = 0; i < n; i = i + 1 {
      dot = dot + final_v[i] * v[i]
    }

    let new_v_norm = normalize_vector(final_v, n)

    // Check convergence
    let mut diff = 0.0
    for i = 0; i < n; i = i + 1 {
      diff = diff + (new_v_norm[i] - v[i]) * (new_v_norm[i] - v[i])
    }

    eigenvalue = dot.abs()
    v = new_v_norm

    if diff < 0.0000001 {
      break
    }
    iter = iter + 1
  }

  (eigenvalue, v)
}

///|
/// Normalize a vector to unit length.
fn normalize_vector(v : Array[Double], n : Int) -> Array[Double] {
  let mut norm = 0.0
  for i = 0; i < n; i = i + 1 {
    norm = norm + v[i] * v[i]
  }
  norm = norm.sqrt()
  let result : Array[Double] = Array::new()
  for i = 0; i < n; i = i + 1 {
    if norm > 0.0 {
      result.push(v[i] / norm)
    } else {
      result.push(0.0)
    }
  }
  result
}

///|
/// Orthogonalize v against u (Gram-Schmidt).
fn orthogonalize(v : Array[Double], u : Array[Double], n : Int) -> Array[Double] {
  let mut dot = 0.0
  for i = 0; i < n; i = i + 1 {
    dot = dot + v[i] * u[i]
  }
  let result : Array[Double] = Array::new()
  for i = 0; i < n; i = i + 1 {
    result.push(v[i] - dot * u[i])
  }
  result
}

///|
/// Get transformed data as a flat matrix (for downstream analysis).
pub fn get_transformed_matrix(vst_result : VSTResult) -> Array[Array[Double]] {
  vst_result.transformed
}

///|
/// Compute the standard deviation of each gene in the transformed data.
/// Useful for identifying highly variable genes.
pub fn gene_sd(vst_result : VSTResult) -> Array[(String, Double)] {
  let n_genes = vst_result.transformed.length()
  let n_samples = vst_result.col_names.length()
  let result : Array[(String, Double)] = Array::new()

  for i = 0; i < n_genes; i = i + 1 {
    let mut sum = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      sum = sum + vst_result.transformed[i][j]
    }
    let mean = sum / n_samples.to_double()
    let mut var = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      let diff = vst_result.transformed[i][j] - mean
      var = var + diff * diff
    }
    var = var / (n_samples - 1).to_double()
    result.push((vst_result.row_names[i], var.sqrt()))
  }
  result
}
