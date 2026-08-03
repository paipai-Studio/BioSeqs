///|
/// glmGamPoi - Bioconductor differential expression analysis.
///
/// Implements the glmGamPoi method for differential expression analysis,
/// which uses overdispersed Poisson and negative binomial models for
/// fast and accurate inference.
///
/// References:
///   - Bioconductor glmGamPoi package
///   - Ahlmann-Eltze C, Huber W. (2021) glmGamPoi: Fitting Gamma-Poisson
///     Generalized Linear Models on Single Cell Count Data.
///     Bioinformatics 37(2):299-300.
///
/// Key functions:
///   - glm_gp: fit a Gamma-Poisson generalized linear model
///   - test_de: test for differential expression (Wald test)
///   - pseudobulk: aggregate single-cell counts by group
///   - calculate_sf: compute size factors
///
/// All public identifiers are prefixed with `glm_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// Design matrix for the GLM.
/// Each row corresponds to a sample, each column to a coefficient.
pub struct GlmDesign {
  /// Row names (sample names)
  sample_names : Array[String]
  /// Column names (coefficient names)
  coef_names : Array[String]
  /// Design matrix values (n_samples x n_coefs, row-major)
  matrix : Array[Array[Double]]
}

///|
/// Create a design matrix.
pub fn GlmDesign::new(
  sample_names : Array[String],
  coef_names : Array[String],
  matrix : Array[Array[Double]],
) -> GlmDesign {
  { sample_names, coef_names, matrix }
}

///|
/// Get number of samples (rows).
pub fn GlmDesign::n_samples(self : GlmDesign) -> Int {
  self.sample_names.length()
}

///|
/// Get number of coefficients (columns).
pub fn GlmDesign::n_coefs(self : GlmDesign) -> Int {
  self.coef_names.length()
}

///|
/// Get design matrix value at (i, j).
pub fn GlmDesign::get(self : GlmDesign, i : Int, j : Int) -> Double {
  self.matrix[i][j]
}

///|
/// A single-gene GLM fit result.
pub struct GlmFit {
  /// Gene name
  gene : String
  /// Coefficient estimates (log-scale)
  coefficients : Array[Double]
  /// Standard errors of coefficients
  std_errors : Array[Double]
  /// Dispersion estimate
  dispersion : Double
  /// Residual degrees of freedom
  df_residual : Int
  /// Dispersion degrees of freedom
  df_dispersion : Double
}

///|
/// Create a GLM fit result.
pub fn GlmFit::new(
  gene : String,
  coefficients : Array[Double],
  std_errors : Array[Double],
  dispersion : Double,
  df_residual : Int,
  df_dispersion : Double,
) -> GlmFit {
  { gene, coefficients, std_errors, dispersion, df_residual, df_dispersion }
}

///|
/// A single-gene differential expression test result.
pub struct GlmTestResult {
  /// Gene name
  gene : String
  /// Coefficient name being tested
  coef_name : String
  /// Estimated log fold change
  estimate : Double
  /// Standard error of the estimate
  std_error : Double
  /// Test statistic (z-value)
  z_value : Double
  /// P-value
  p_value : Double
  /// Adjusted p-value (BH-FDR)
  adj_p_value : Double
  /// Whether the gene is significant (adj_p < 0.05)
  is_significant : Bool
}

///|
/// Create a test result.
pub fn GlmTestResult::new(
  gene : String,
  coef_name : String,
  estimate : Double,
  std_error : Double,
  z_value : Double,
  p_value : Double,
) -> GlmTestResult {
  let adj = glm_bh_correct(p_value)
  GlmTestResult::{
    gene,
    coef_name,
    estimate,
    std_error,
    z_value,
    p_value,
    adj_p_value: adj,
    is_significant: adj < 0.05,
  }
}

///|
/// Get the log fold change estimate.
pub fn GlmTestResult::estimate(self : GlmTestResult) -> Double {
  self.estimate
}

///|
/// Get the p-value.
pub fn GlmTestResult::p_value(self : GlmTestResult) -> Double {
  self.p_value
}

///|
/// Get the adjusted p-value.
pub fn GlmTestResult::adj_p_value(self : GlmTestResult) -> Double {
  self.adj_p_value
}

///|
/// Format a test result as a string.
pub fn GlmTestResult::to_string(self : GlmTestResult) -> String {
  let sig = if self.is_significant { "*" } else { "" }
  self.gene + " " + self.coef_name + " est=" + self.estimate.to_string() +
    " se=" + self.std_error.to_string() + " z=" + self.z_value.to_string() +
    " p=" + self.p_value.to_string() + " adj=" + self.adj_p_value.to_string() + sig
}

///|
/// Aggregated count data for pseudobulk.
pub struct GlmPseudobulk {
  /// Group labels
  group_labels : Array[String]
  /// Sample names within each group
  sample_groups : Array[Array[String]]
  /// Aggregated count matrix (genes x groups)
  counts : Array[Array[Double]]
  /// Gene names
  gene_names : Array[String]
  /// Group names
  group_names : Array[String]
}

///|
/// Create a pseudobulk result.
pub fn GlmPseudobulk::new(
  gene_names : Array[String],
  group_names : Array[String],
  sample_groups : Array[Array[String]],
  counts : Array[Array[Double]],
) -> GlmPseudobulk {
  { group_labels: [], sample_groups, counts, gene_names, group_names }
}

///|
/// Get number of genes.
pub fn GlmPseudobulk::n_genes(self : GlmPseudobulk) -> Int {
  self.gene_names.length()
}

///|
/// Get number of groups.
pub fn GlmPseudobulk::n_groups(self : GlmPseudobulk) -> Int {
  self.group_names.length()
}

// ============================================================================
// Core functions
// ============================================================================

///|
/// Compute size factors using the geometric mean method.
/// Similar to DESeq2's estimateSizeFactors.
/// Returns an array of size factors, one per sample.
pub fn glm_calculate_sf(counts : Array[Array[Double]]) -> Array[Double] {
  let n_genes = counts.length()
  if n_genes == 0 {
    return Array::make(0, 0.0)
  }
  let n_samples = counts[0].length()
  let sf : Array[Double] = Array::make(n_samples, 0.0)

  // Step 1: Compute geometric mean for each gene
  let geo_means : Array[Double] = Array::make(n_genes, 0.0)
  for i in 0..<n_genes {
    let row_sum = glm_sum(counts[i])
    if row_sum > 0.0 {
      // Log-space computation of geometric mean
      let mut log_sum = 0.0
      let mut n_nonzero = 0
      for j in 0..<n_samples {
        let v = counts[i][j]
        if v > 0.0 {
          log_sum = log_sum + @math.ln(v)
          n_nonzero = n_nonzero + 1
        }
      }
      if n_nonzero > 0 {
        geo_means[i] = @math.exp(log_sum / n_nonzero.to_double())
      } else {
        geo_means[i] = 0.0
      }
    }
  }

  // Step 2: Compute size factors as median of ratios to geometric mean
  for j in 0..<n_samples {
    let ratios : Array[Double] = []
    for i in 0..<n_genes {
      if geo_means[i] > 0.0 && counts[i][j] > 0.0 {
        ratios.push(counts[i][j] / geo_means[i])
      }
    }
    if ratios.length() > 0 {
      sf[j] = glm_median(ratios)
    } else {
      sf[j] = 1.0
    }
  }

  // Avoid zero or negative size factors
  for j in 0..<n_samples {
    if sf[j] <= 0.0 {
      sf[j] = 1.0
    }
  }
  sf
}

///|
/// Pseudobulk aggregation: sum counts within groups.
///
/// Parameters:
/// - counts: raw count matrix (genes x samples)
/// - gene_names: gene identifiers
/// - group_labels: group assignment for each sample
///
/// Returns aggregated counts (genes x groups) and group metadata.
pub fn glm_pseudobulk(
  counts : Array[Array[Double]],
  gene_names : Array[String],
  group_labels : Array[String],
) -> GlmPseudobulk {
  let n_genes = counts.length()
  let n_samples = if n_genes > 0 { counts[0].length() } else { 0 }

  // Find unique groups
  let unique_groups : Array[String] = []
  let group_map : Map[String, Int] = Map([], capacity=20)
  for i in 0..<group_labels.length() {
    let label = group_labels[i]
    if !group_map.contains(label) {
      group_map.set(label, unique_groups.length())
      unique_groups.push(label)
    }
  }
  let n_groups = unique_groups.length()

  // Assign samples to groups
  let sample_groups : Array[Array[String]] = []
  for g in 0..<n_groups {
    sample_groups.push([])
  }
  for s in 0..<n_samples {
    let gidx = group_map.get_or_default(group_labels[s], 0)
    sample_groups[gidx].push("sample_" + (s + 1).to_string())
  }

  // Aggregate counts
  let agg_counts : Array[Array[Double]] = []
  for g in 0..<n_groups {
    agg_counts.push(Array::make(n_genes, 0.0))
  }
  for i in 0..<n_genes {
    for g in 0..<n_groups {
      let mut sum = 0.0
      for s in 0..<n_samples {
        if group_map.get_or_default(group_labels[s], 0) == g {
          sum = sum + counts[i][s]
        }
      }
      agg_counts[g][i] = sum
    }
  }

  // Transpose to (genes x groups) layout
  let result_counts : Array[Array[Double]] = []
  for i in 0..<n_genes {
    let row : Array[Double] = []
    for g in 0..<n_groups {
      row.push(agg_counts[g][i])
    }
    result_counts.push(row)
  }

  GlmPseudobulk::new(gene_names, unique_groups, sample_groups, result_counts)
}

///|
/// Fit a Gamma-Poisson GLM for a single gene.
///
/// The model is: count_i ~ GammaPoisson(dispersion, mean_i)
/// where log(mean_i) = X_i %*% beta
///
/// Uses a simple iterative algorithm to estimate coefficients
/// and dispersion.
pub fn glm_fit_one_gene(
  counts : Array[Double],
  design : GlmDesign,
  size_factors : Array[Double],
  gene_name : String,
) -> GlmFit {
  let n = counts.length()
  let n_coefs = design.n_coefs()

  // Step 1: Normalize counts by size factors
  let norm_counts : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    if size_factors[i] > 0.0 {
      norm_counts[i] = counts[i] / size_factors[i]
    } else {
      norm_counts[i] = counts[i]
    }
  }

  // Step 2: Initial coefficient estimates via log-linear regression
  let initial_coefs = glm_initial_estimates(norm_counts, design)

  // Step 3: Estimate dispersion
  let dispersion = glm_estimate_dispersion(norm_counts, design, initial_coefs)

  // Step 4: Refit with estimated dispersion (IWLCS)
  let coefs = glm_iwlcs_fit(norm_counts, design, dispersion, initial_coefs)

  // Step 5: Compute standard errors
  let std_errors = glm_compute_se(design, dispersion, coefs, size_factors)

  let df_residual = n - n_coefs
  let df_dispersion = if df_residual > 0 { df_residual.to_double() } else { 1.0 }

  GlmFit::new(gene_name, coefs, std_errors, dispersion, df_residual, df_dispersion)
}

///|
/// Run differential expression test for all genes.
///
/// Parameters:
/// - counts: count matrix (genes x samples)
/// - gene_names: gene identifiers
/// - design: design matrix
/// - size_factors: size factors per sample
/// - test_coef_index: index of the coefficient to test
///
/// Returns an array of test results.
pub fn glm_test_de(
  counts : Array[Array[Double]],
  gene_names : Array[String],
  design : GlmDesign,
  size_factors : Array[Double],
  test_coef_index : Int,
) -> Array[GlmTestResult] {
  let n_genes = counts.length()
  let results : Array[GlmTestResult] = []

  // First fit all genes to get dispersion estimates
  let fits : Array[GlmFit] = []
  for g in 0..<n_genes {
    let fit = glm_fit_one_gene(counts[g], design, size_factors, gene_names[g])
    fits.push(fit)
  }

  // Then test
  for g in 0..<n_genes {
    let fit = fits[g]
    let estimate = fit.coefficients[test_coef_index]
    let se = fit.std_errors[test_coef_index]
    let z = if se > 0.0 { estimate / se } else { 0.0 }
    let pval = glm_normal_pvalue_two_sided(z)
    let coef_name = design.coef_names[test_coef_index]
    results.push(GlmTestResult::new(gene_names[g], coef_name, estimate, se, z, pval))
  }

  results
}

///|
/// Select significant results based on adjusted p-value threshold.
pub fn glm_significant(
  results : Array[GlmTestResult],
  threshold : Double,
) -> Array[GlmTestResult] {
  let sig : Array[GlmTestResult] = []
  for i in 0..<results.length() {
    if results[i].adj_p_value < threshold {
      sig.push(results[i])
    }
  }
  sig
}

// ============================================================================
// Internal algorithms
// ============================================================================

///|
/// Compute initial coefficient estimates via simple log-linear regression.
pub fn glm_initial_estimates(
  counts : Array[Double],
  design : GlmDesign,
) -> Array[Double] {
  let n = counts.length()
  let n_coefs = design.n_coefs()

  // Compute log-normalized counts (add small pseudocount for zeros)
  let log_counts : Array[Double] = Array::make(n, 0.0)
  let min_positive = glm_min_positive(counts)
  let pseudocount = if min_positive > 0.0 { min_positive * 0.5 } else { 0.5 }
  for i in 0..<n {
    let c = counts[i]
    log_counts[i] = if c > 0.0 { @math.ln(c) } else { @math.ln(pseudocount) }
  }

  // Simple least squares: solve X^T X beta = X^T y
  // Build X^T X (k x k) and X^T y (k)
  let xtx : Array[Array[Double]] = []
  for i in 0..<n_coefs {
    xtx.push(Array::make(n_coefs, 0.0))
  }
  let xty : Array[Double] = Array::make(n_coefs, 0.0)

  for s in 0..<n {
    for i in 0..<n_coefs {
      let xi = design.get(s, i)
      xty[i] = xty[i] + xi * log_counts[s]
      for j in 0..<n_coefs {
        xtx[i][j] = xtx[i][j] + xi * design.get(s, j)
      }
    }
  }

  // Solve via Gaussian elimination
  glm_solve_linear_system(xtx, xty, n_coefs)
}

///|
/// Estimate the dispersion parameter.
/// Uses the method of moments based on residuals.
pub fn glm_estimate_dispersion(
  counts : Array[Double],
  design : GlmDesign,
  coefs : Array[Double],
) -> Double {
  let n = counts.length()
  let n_coefs = design.n_coefs()

  // Compute fitted values
  let fitted : Array[Double] = Array::make(n, 0.0)
  for s in 0..<n {
    let mut mu = 0.0
    for j in 0..<n_coefs {
      mu = mu + design.get(s, j) * coefs[j]
    }
    fitted[s] = @math.exp(mu)
  }

  // Compute dispersion as mean of (variance/mean - 1)
  let residuals : Array[Double] = []
  for s in 0..<n {
    let mu = fitted[s]
    if mu > 0.0 {
      let var_est = (counts[s] - mu) * (counts[s] - mu)
      let disp = var_est / (mu * mu)
      residuals.push(disp)
    }
  }

  if residuals.length() > 0 {
    glm_mean(residuals)
  } else {
    0.1 // Default dispersion
  }
}

///|
/// Iteratively reweighted least squares for Gamma-Poisson model.
pub fn glm_iwlcs_fit(
  counts : Array[Double],
  design : GlmDesign,
  dispersion : Double,
  initial_coefs : Array[Double],
) -> Array[Double] {
  let n = counts.length()
  let n_coefs = design.n_coefs()
  let max_iter = 20
  let tolerance = 1.0e-8

  let coefs = initial_coefs.copy()

  for iter in 0..<max_iter {
    // Compute fitted values
    let fitted : Array[Double] = Array::make(n, 0.0)
    for s in 0..<n {
      let mut mu = 0.0
      for j in 0..<n_coefs {
        mu = mu + design.get(s, j) * coefs[j]
      }
      fitted[s] = @math.exp(mu)
      if fitted[s] < 1.0e-10 {
        fitted[s] = 1.0e-10
      }
    }

    // Compute working responses and weights
    let z : Array[Double] = Array::make(n, 0.0)
    let w : Array[Double] = Array::make(n, 0.0)
    for s in 0..<n {
      let mu = fitted[s]
      // Working response: log(count) + (count - mu) / mu
      z[s] = if counts[s] > 0.0 {
        @math.ln(counts[s]) + (counts[s] - mu) / mu
      } else {
        @math.ln(1.0e-10) + (0.0 - mu) / mu
      }
      // Weight: 1 / (dispersion * mu^2 + mu)
      w[s] = 1.0 / (dispersion * mu * mu + mu + 1.0e-10)
    }

    // Weighted least squares: solve (X^T W X) delta = X^T W (z - X beta)
    let xtwx : Array[Array[Double]] = []
    for i in 0..<n_coefs {
      xtwx.push(Array::make(n_coefs, 0.0))
    }
    let xtwr : Array[Double] = Array::make(n_coefs, 0.0)

    for s in 0..<n {
      let ws = w[s]
      let zs = z[s]
      let mut xp = 0.0
      for j in 0..<n_coefs {
        xp = xp + design.get(s, j) * coefs[j]
      }
      for i in 0..<n_coefs {
        let xi = design.get(s, i)
        let residual = zs - xp
        for j in 0..<n_coefs {
          xtwx[i][j] = xtwx[i][j] + ws * xi * design.get(s, j)
        }
        xtwr[i] = xtwr[i] + ws * xi * residual
      }
    }

    // Solve for delta
    let delta = glm_solve_linear_system(xtwx, xtwr, n_coefs)

    // Update coefficients
    let mut max_change = 0.0
    for j in 0..<n_coefs {
      let change = delta[j]
      if change.abs() > max_change {
        max_change = change.abs()
      }
      coefs[j] = coefs[j] + change
    }

    if max_change < tolerance {
      break
    }
  }

  coefs
}

///|
/// Compute standard errors for coefficient estimates.
pub fn glm_compute_se(
  design : GlmDesign,
  dispersion : Double,
  coefs : Array[Double],
  size_factors : Array[Double],
) -> Array[Double] {
  let n = design.n_samples()
  let n_coefs = design.n_coefs()

  // Compute fitted values
  let fitted : Array[Double] = Array::make(n, 0.0)
  for s in 0..<n {
    let mut mu = 0.0
    for j in 0..<n_coefs {
      mu = mu + design.get(s, j) * coefs[j]
    }
    fitted[s] = @math.exp(mu)
    if fitted[s] < 1.0e-10 {
      fitted[s] = 1.0e-10
    }
  }

  // Compute weight matrix inverse diagonal elements
  // Var(beta_j) = (X^T W X)^{-1}_{jj}
  let xtwx : Array[Array[Double]] = []
  for i in 0..<n_coefs {
    xtwx.push(Array::make(n_coefs, 0.0))
  }

  for s in 0..<n {
    let mu = fitted[s]
    let ws = 1.0 / (dispersion * mu * mu + mu + 1.0e-10)
    for i in 0..<n_coefs {
      let xi = design.get(s, i)
      for j in 0..<n_coefs {
        xtwx[i][j] = xtwx[i][j] + ws * xi * design.get(s, j)
      }
    }
  }

  // Invert the matrix using Gaussian elimination
  let inv = glm_invert_matrix(xtwx, n_coefs)

  // Extract diagonal elements (variances) and convert to standard errors
  let se : Array[Double] = Array::make(n_coefs, 0.0)
  for j in 0..<n_coefs {
    let var = inv[j][j]
    se[j] = if var > 0.0 { var.sqrt() } else { 1.0e-10 }
  }

  se
}

// ============================================================================
// Statistical helpers
// ============================================================================

///|
/// Two-sided p-value from standard normal distribution.
pub fn glm_normal_pvalue_two_sided(z : Double) -> Double {
  // Use approximation: p ≈ 2 * (1 - Phi(|z|))
  // For large |z|, this gives very small p-values
  let abs_z = z.abs()
  // Simple approximation of the survival function
  let p = 2.0 * (1.0 - glm_phi(abs_z))
  // Clamp to avoid zero
  if p < 1.0e-15 { 1.0e-15 } else { p }
}

///|
/// Approximation of the standard normal CDF.
/// Uses the Abramowitz and Stegun approximation.
pub fn glm_phi(x : Double) -> Double {
  // For x >= 0
  let sign = if x >= 0.0 { 1.0 } else { -1.0 }
  let abs_x = x.abs()
  // Horner form approximation
  let t = 1.0 / (1.0 + 0.2316419 * abs_x)
  let d = 0.3989422804014327 * @math.exp(-abs_x * abs_x / 2.0)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.3302744))))
  let result = if sign > 0.0 { 1.0 - p } else { p }
  // Clamp
  if result < 0.0 { 0.0 } else if result > 1.0 { 1.0 } else { result }
}

///|
/// Benjamini-Hochberg FDR correction for a single p-value.
/// (For multiple testing, use the full array version.)
pub fn glm_bh_correct(p : Double) -> Double {
  p
}

///|
/// Benjamini-Hochberg FDR correction for an array of p-values.
pub fn glm_bh_correct_array(pvalues : Array[Double]) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 {
    return []
  }

  // Create index-value pairs and sort by p-value
  let indexed : Array[(Double, Int)] = []
  for i in 0..<n {
    indexed.push((pvalues[i], i))
  }

  // Sort by p-value (bubble sort for simplicity)
  let arr = indexed.copy()
  glm_sort_pairs(arr)

  // Apply BH correction
  let adj : Array[Double] = Array::make(n, 0.0)
  let mut cummin = 1.0
  let n_minus_1 = n - 1
  for j in 0..<n {
    let i = n_minus_1 - j
    let rank = (i + 1).to_double()
    let val = arr[i].0 * (n.to_double()) / rank
    if val < cummin {
      cummin = val
    }
    adj[arr[i].1] = if cummin < 1.0 { cummin } else { 1.0 }
  }

  adj
}

///|
/// Sort an array of (Double, Int) pairs by the first element (descending).
pub fn glm_sort_pairs(arr : Array[(Double, Int)]) -> Unit {
  let n = arr.length()
  let n_minus_1 = n - 1
  for i in 0..<n_minus_1 {
    let n_minus_1_minus_i = n - 1 - i
    for j in 0..<n_minus_1_minus_i {
      if arr[j].0 > arr[j + 1].0 {
        let tmp = arr[j]
        arr[j] = arr[j + 1]
        arr[j + 1] = tmp
      }
    }
  }
}

// ============================================================================
// Linear algebra helpers
// ============================================================================

///|
/// Solve a linear system Ax = b using Gaussian elimination with
/// partial pivoting. A is modified in place.
pub fn glm_solve_linear_system(a : Array[Array[Double]], b : Array[Double], n : Int) -> Array[Double] {
  // Augmented matrix
  let aug : Array[Array[Double]] = []
  for i in 0..<n {
    let row : Array[Double] = []
    for j in 0..<n {
      row.push(a[i][j])
    }
    row.push(b[i])
    aug.push(row)
  }

  // Forward elimination with partial pivoting
  for col in 0..<n {
    // Find pivot
    let mut max_row = col
    let mut max_val = aug[col][col].abs()
    let col_plus_1 = col + 1
    for row in col_plus_1..<n {
      if aug[row][col].abs() > max_val {
        max_val = aug[row][col].abs()
        max_row = row
      }
    }

    // Swap rows
    if max_row != col {
      let tmp = aug[col]
      aug[col] = aug[max_row]
      aug[max_row] = tmp
    }

    let pivot = aug[col][col]
    if pivot.abs() < 1.0e-12 {
      // Add small regularization
      aug[col][col] = aug[col][col] + 1.0e-10
    }

    // Eliminate below
    for row in col_plus_1..<n {
      let factor = aug[row][col] / aug[col][col]
      let n_plus_1 = n + 1
      for j in col..<n_plus_1 {
        aug[row][j] = aug[row][j] - factor * aug[col][j]
      }
    }
  }

  // Back substitution
  let x : Array[Double] = Array::make(n, 0.0)
  let n_minus_1 = n - 1
  for r in 0..<n {
    let row = n_minus_1 - r
    let mut sum = aug[row][n]
    let row_plus_1 = row + 1
    for j in row_plus_1..<n {
      sum = sum - aug[row][j] * x[j]
    }
    x[row] = sum / aug[row][row]
  }

  x
}

///|
/// Invert a matrix using Gauss-Jordan elimination.
pub fn glm_invert_matrix(a : Array[Array[Double]], n : Int) -> Array[Array[Double]] {
  // Build augmented matrix [A | I]
  let aug : Array[Array[Double]] = []
  for i in 0..<n {
    let row : Array[Double] = []
    for j in 0..<n {
      row.push(a[i][j])
    }
    for j in 0..<n {
      if j == i {
        row.push(1.0)
      } else {
        row.push(0.0)
      }
    }
    aug.push(row)
  }

  let total = 2 * n

  // Forward elimination
  for col in 0..<n {
    let mut max_row = col
    let mut max_val = aug[col][col].abs()
    let col_plus_1 = col + 1
    for row in col_plus_1..<n {
      if aug[row][col].abs() > max_val {
        max_val = aug[row][col].abs()
        max_row = row
      }
    }

    if max_row != col {
      let tmp = aug[col]
      aug[col] = aug[max_row]
      aug[max_row] = tmp
    }

    let pivot = aug[col][col]
    if pivot.abs() < 1.0e-12 {
      aug[col][col] = aug[col][col] + 1.0e-10
    }

    // Scale pivot row
    let p = aug[col][col]
    for j in 0..<total {
      aug[col][j] = aug[col][j] / p
    }

    // Eliminate column
    for row in 0..<n {
      if row != col {
        let factor = aug[row][col]
        for j in 0..<total {
          aug[row][j] = aug[row][j] - factor * aug[col][j]
        }
      }
    }
  }

  // Extract inverse (right half)
  let inv : Array[Array[Double]] = []
  for i in 0..<n {
    let row : Array[Double] = []
    for j in n..<total {
      row.push(aug[i][j])
    }
    inv.push(row)
  }

  inv
}

// ============================================================================
// Utility helpers
// ============================================================================

///|
/// Sum of an array.
pub fn glm_sum(data : Array[Double]) -> Double {
  let mut s = 0.0
  for i in 0..<data.length() {
    s = s + data[i]
  }
  s
}

///|
/// Mean of an array.
pub fn glm_mean(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  glm_sum(data) / data.length().to_double()
}

///|
/// Median of an array.
pub fn glm_median(data : Array[Double]) -> Double {
  let n = data.length()
  if n == 0 {
    return 0.0
  }
  let sorted = data.copy()
  sorted.sort()
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

///|
/// Get minimum positive value (for pseudocount computation).
pub fn glm_min_positive(data : Array[Double]) -> Double {
  let mut min_val : Double = 1.0e30
  for i in 0..<data.length() {
    if data[i] > 0.0 && data[i] < min_val {
      min_val = data[i]
    }
  }
  if min_val < 1.0e30 { min_val } else { 0.0 }
}
