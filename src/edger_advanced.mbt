///|
/// edgeR Advanced - Quasi-Likelihood F-test and Gene Set Tests.
/// Inspired by Bioconductor edgeR package.
///
/// Provides:
/// - glm_qlf_fit: Fit quasi-likelihood GLM and estimate QL dispersions
/// - glm_qlf_test: Quasi-likelihood F-test for differential expression
/// - camera: Competitive gene set test
/// - roast: Self-contained gene set test

///|
/// QL dispersion fit result containing fitted GLM and QL dispersions.
pub struct QLFData {
  /// Gene names
  genes : Array[String]
  /// Coefficients from GLM fit (genes x coefficients)
  coefficients : Array[Array[Double]]
  /// Fitted values from GLM (genes x samples)
  fitted_values : Array[Array[Double]]
  /// QL dispersion for each gene (shrunk towards trend)
  ql_dispersion : Array[Double]
  /// Common dispersion trend value
  trend_dispersion : Double
  /// Sample group labels
  group : Array[String]
  /// Library sizes (with norm factors applied)
  lib_size : Array[Double]
  /// Original counts
  counts : Array[Array[Int]]
}

///|
/// QL F-test result for differential expression.
pub struct QLFResult {
  /// Gene names
  genes : Array[String]
  /// Log-fold change
  log_fc : Array[Double]
  /// F-statistic
  f_stat : Array[Double]
  /// Degrees of freedom (numerator, denominator)
  df_num : Int
  df_denom : Int
  /// P-value from F-test
  p_value : Array[Double]
  /// FDR adjusted p-value
  fdr : Array[Double]
  /// QL dispersion used
  ql_dispersion : Array[Double]
}

///|
/// Fit a quasi-likelihood negative binomial GLM.
/// Estimates QL dispersions with empirical Bayes shrinkage towards
/// a common trend.
///
/// Parameters:
/// - dge: A DGEList with counts, groups, and library sizes
///
/// Returns a QLFData with fitted values and QL dispersions.
pub fn glm_qlf_fit(dge : DGEList) -> QLFData {
  let n_genes = dge.counts.length()
  let n_samples = dge.group.length()

  // Compute effective library sizes
  let eff_lib : Array[Double] = Array::new()
  for j = 0; j < n_samples; j = j + 1 {
    eff_lib.push(dge.lib_size[j] * dge.norm_factors[j])
  }

  // Design matrix: intercept + group indicator (2-group design)
  // For simplicity, we use a 2-column design: [1, group_indicator]
  let group_labels : Array[String] = Array::new()
  let mut has_group0 = false
  let mut has_group1 = false
  for g in dge.group {
    if !group_labels.contains(g) {
      group_labels.push(g)
    }
  }
  let n_groups = group_labels.length()

  // Build design matrix
  let design : Array[Array[Double]] = Array::new()
  for j = 0; j < n_samples; j = j + 1 {
    let row : Array[Double] = [1.0]
    if n_groups > 1 {
      if dge.group[j] == group_labels[1] {
        row.push(1.0)
      } else {
        row.push(0.0)
      }
    }
    design.push(row)
  }

  let n_coef = if n_groups > 1 { 2 } else { 1 }

  // Fit NB GLM for each gene using IRLS (Iteratively Reweighted Least Squares)
  let coefficients : Array[Array[Double]] = Array::new()
  let fitted_values : Array[Array[Double]] = Array::new()
  let raw_ql_disp : Array[Double] = Array::new()

  for i = 0; i < n_genes; i = i + 1 {
    let (coefs, fitted, deviance) = fit_nb_glm(
      dge.counts[i],
      design,
      eff_lib,
      n_samples,
      n_coef,
    )
    coefficients.push(coefs)
    fitted_values.push(fitted)

    // Compute raw QL dispersion from Pearson residual deviance
    // QL dispersion = deviance / residual_df
    let residual_df = n_samples - n_coef
    let ql_disp = if residual_df > 0 {
      let d = deviance / residual_df.to_double()
      if d < 0.0 { 0.001 } else if d > 100.0 { 100.0 } else { d }
    } else {
      1.0
    }
    raw_ql_disp.push(ql_disp)
  }

  // Estimate trend dispersion (median of raw QL dispersions)
  let sorted_disps = sort_array(raw_ql_disp.copy())
  let trend = if n_genes > 0 {
    if n_genes % 2 == 0 {
      (sorted_disps[n_genes / 2 - 1] + sorted_disps[n_genes / 2]) / 2.0
    } else {
      sorted_disps[n_genes / 2]
    }
  } else {
    1.0
  }

  // Empirical Bayes shrinkage of QL dispersions towards trend
  // Using a simple shrinkage: ql_shrunk = (raw_df * raw + prior_df * trend) / (raw_df + prior_df)
  let prior_df = 2.0 // prior degrees of freedom
  let residual_df = (n_samples - n_coef).to_double()
  let shrunk_ql_disp : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let shrunk = (residual_df * raw_ql_disp[i] + prior_df * trend) /
      (residual_df + prior_df)
    shrunk_ql_disp.push(if shrunk < 0.001 { 0.001 } else { shrunk })
  }

  QLFData::{
    genes: dge.genes,
    coefficients,
    fitted_values,
    ql_dispersion: shrunk_ql_disp,
    trend_dispersion: trend,
    group: dge.group,
    lib_size: eff_lib,
    counts: dge.counts,
  }
}

///|
/// Fit a negative binomial GLM using IRLS for a single gene.
/// Returns (coefficients, fitted_values, deviance).
fn fit_nb_glm(
  counts : Array[Int],
  design : Array[Array[Double]],
  eff_lib : Array[Double],
  n_samples : Int,
  n_coef : Int,
) -> (Array[Double], Array[Double], Double) {
  // Initialize coefficients
  let coefs : Array[Double] = Array::make(n_coef, 0.0)

  // Initialize with log of mean count
  let mut total = 0.0
  for j = 0; j < n_samples; j = j + 1 {
    total = total + counts[j].to_double()
  }
  let mean_count = total / n_samples.to_double()
  let init_log_mu = if mean_count > 0.0 {
    log_nat(mean_count)
  } else {
    0.0
  }
  for k = 0; k < n_coef; k = k + 1 {
    coefs[k] = if k == 0 { init_log_mu } else { 0.0 }
  }

  // Estimate initial dispersion (NB dispersion)
  let nb_disp = 0.1 // common dispersion

  // IRLS iterations
  let mut iter = 0
  let fitted : Array[Double] = Array::make(n_samples, mean_count)
  let mut deviance = 0.0

  while iter < 50 {
    // Compute fitted values: mu = exp(X * beta) * lib_size
    for j = 0; j < n_samples; j = j + 1 {
      let mut eta = 0.0
      for k = 0; k < n_coef; k = k + 1 {
        eta = eta + design[j][k] * coefs[k]
      }
      fitted[j] = exp(eta) * eff_lib[j] / 1000000.0 // offset by library size
      if fitted[j] < 0.001 {
        fitted[j] = 0.001
      }
    }

    // IRLS weights and working response
    let working : Array[Double] = Array::make(n_samples, 0.0)
    let weights : Array[Double] = Array::make(n_samples, 0.0)

    for j = 0; j < n_samples; j = j + 1 {
      let mu = fitted[j]
      let y = counts[j].to_double()
      // NB variance: mu + disp * mu^2
      let var = mu + nb_disp * mu * mu
      let dmu_deta = mu // since link is log, dmu/deta = mu
      let deta_dmu = 1.0 / mu

      let w = dmu_deta * dmu_deta / var
      let z = eta_from_coefs(coefs, design[j], n_coef) + (y - mu) * deta_dmu

      working[j] = z
      weights[j] = w
    }

    // Weighted least squares: solve (X^T W X) beta = X^T W z
    let new_coefs = solve_wls(design, weights, working, n_samples, n_coef)

    // Check convergence
    let mut max_change = 0.0
    for k = 0; k < n_coef; k = k + 1 {
      let diff = (new_coefs[k] - coefs[k]).abs()
      if diff > max_change {
        max_change = diff
      }
      coefs[k] = new_coefs[k]
    }

    if max_change < 0.0001 {
      break
    }
    iter = iter + 1
  }

  // Compute deviance (using NB deviance)
  deviance = 0.0
  for j = 0; j < n_samples; j = j + 1 {
    let y = counts[j].to_double()
    let mu = fitted[j]
    if y > 0.0 && mu > 0.0 {
      let term = y * log_nat(y / mu) - (y + 1.0 / nb_disp) * log_nat((1.0 + nb_disp * y) / (1.0 + nb_disp * mu))
      deviance = deviance + 2.0 * term
    } else if y == 0.0 && mu > 0.0 {
      deviance = deviance + 2.0 * (1.0 / nb_disp) * log_nat(1.0 + nb_disp * mu)
    }
  }

  (coefs, fitted, deviance)
}

///|
/// Compute linear predictor eta from coefficients and design row.
fn eta_from_coefs(
  coefs : Array[Double],
  design_row : Array[Double],
  n_coef : Int,
) -> Double {
  let mut eta = 0.0
  for k = 0; k < n_coef; k = k + 1 {
    eta = eta + design_row[k] * coefs[k]
  }
  eta
}

///|
/// Solve weighted least squares: (X^T W X) beta = X^T W z
/// Using normal equations with 2x2 or 1x1 matrix inversion.
fn solve_wls(
  design : Array[Array[Double]],
  weights : Array[Double],
  z : Array[Double],
  n_samples : Int,
  n_coef : Int,
) -> Array[Double] {
  if n_coef == 1 {
    // Simple weighted mean
    let mut sum_wz = 0.0
    let mut sum_w = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      sum_wz = sum_wz + weights[j] * design[j][0] * z[j]
      sum_w = sum_w + weights[j] * design[j][0] * design[j][0]
    }
    return [if sum_w > 0.0 { sum_wz / sum_w } else { 0.0 }]
  }

  // 2x2 system: XtWX
  let mut xtwx_00 = 0.0
  let mut xtwx_01 = 0.0
  let mut xtwx_11 = 0.0
  let mut xtwz_0 = 0.0
  let mut xtwz_1 = 0.0

  for j = 0; j < n_samples; j = j + 1 {
    let x0 = design[j][0]
    let x1 = design[j][1]
    let w = weights[j]
    xtwx_00 = xtwx_00 + w * x0 * x0
    xtwx_01 = xtwx_01 + w * x0 * x1
    xtwx_11 = xtwx_11 + w * x1 * x1
    xtwz_0 = xtwz_0 + w * x0 * z[j]
    xtwz_1 = xtwz_1 + w * x1 * z[j]
  }

  let det = xtwx_00 * xtwx_11 - xtwx_01 * xtwx_01
  if det.abs() < 0.0000001 {
    return [0.0, 0.0]
  }

  let beta0 = (xtwx_11 * xtwz_0 - xtwx_01 * xtwz_1) / det
  let beta1 = (xtwx_00 * xtwz_1 - xtwx_01 * xtwz_0) / det
  [beta0, beta1]
}

///|
/// Quasi-likelihood F-test for differential expression.
/// Tests whether a coefficient is significantly different from zero.
///
/// Parameters:
/// - qlf_data: Output from glm_qlf_fit
/// - coef: Index of coefficient to test (default: 1, the group coefficient)
///
/// Returns a QLFResult with F-statistics and p-values.
pub fn glm_qlf_test(qlf_data : QLFData, coef? : Int = 1) -> QLFResult {
  let n_genes = qlf_data.genes.length()
  let n_samples = qlf_data.group.length()
  let n_coef = if qlf_data.coefficients.length() > 0 {
    qlf_data.coefficients[0].length()
  } else {
    1
  }
  let residual_df = n_samples - n_coef

  let log_fc : Array[Double] = Array::new()
  let f_stat : Array[Double] = Array::new()
  let p_value : Array[Double] = Array::new()
  let ql_disp : Array[Double] = Array::new()

  for i = 0; i < n_genes; i = i + 1 {
    // Log-fold change from the coefficient
    let beta = if coef < qlf_data.coefficients[i].length() {
      qlf_data.coefficients[i][coef]
    } else {
      0.0
    }
    log_fc.push(beta)

    // Compute LR statistic: compare full vs reduced model
    // For a simple 2-group comparison, the F-statistic is:
    // F = beta^2 / (ql_disp * var(beta))
    // where var(beta) comes from the Fisher information

    let ql = qlf_data.ql_dispersion[i]
    ql_disp.push(ql)

    // Compute variance of coefficient from fitted values
    // For log-link NB: var(beta) = (X^T W X)^{-1} where W = diag(mu / (mu + alpha*mu^2))
    let mut info_coef = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      let mu = qlf_data.fitted_values[i][j]
      if mu > 0.0 {
        let w = mu / (1.0 + ql * mu)
        // Design matrix: [1, x_j]
        let x_coef = if coef == 0 { 1.0 } else { if qlf_data.group[j] == qlf_data.group[0] { 0.0 } else { 1.0 } }
        info_coef = info_coef + w * x_coef * x_coef
      }
    }

    let var_beta = if info_coef > 0.0 { 1.0 / info_coef } else { 1.0 }
    let f = if ql * var_beta > 0.0 {
      beta * beta / (ql * var_beta)
    } else {
      0.0
    }
    f_stat.push(f)

    // P-value from F-distribution with (1, residual_df) degrees of freedom
    let df1 = 1
    let df2 = if residual_df > 0 { residual_df } else { 1 }
    let p = f_distribution_pvalue(f, df1, df2)
    p_value.push(p)
  }

  // Compute FDR using Benjamini-Hochberg
  let fdr = bh_adjust(p_value)

  QLFResult::{
    genes: qlf_data.genes,
    log_fc,
    f_stat,
    df_num: 1,
    df_denom: residual_df,
    p_value,
    fdr,
    ql_dispersion: ql_disp,
  }
}

///|
/// Compute p-value from F-distribution using the survival function.
/// F(df1, df2): P(F > f) = I_{df2/(df2 + df1*f)}(df2/2, df1/2)
/// where I is the regularized incomplete beta function.
fn f_distribution_pvalue(f : Double, df1 : Int, df2 : Int) -> Double {
  if f <= 0.0 {
    return 1.0
  }
  if df1 <= 0 || df2 <= 0 {
    return 1.0
  }

  // Use the relationship with the beta distribution:
  // P(F > f) = I_x(df2/2, df1/2) where x = df2 / (df2 + df1 * f)
  let x = df2.to_double() / (df2.to_double() + df1.to_double() * f)
  let a = df2.to_double() / 2.0
  let b = df1.to_double() / 2.0

  let p = regularized_incomplete_beta(x, a, b)
  if p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p }
}

///|
/// Regularized incomplete beta function I_x(a, b).
/// Uses continued fraction expansion (Lentz's method).
fn regularized_incomplete_beta(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }

  // Log of prefactor: B(a,b) = Gamma(a)Gamma(b)/Gamma(a+b)
  let lbeta = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
  let prefactor = exp(a * log_nat(x) + b * log_nat(1.0 - x) - lbeta)

  // Use symmetry: I_x(a,b) = 1 - I_{1-x}(b,a)
  if x < (a + 1.0) / (a + b + 2.0) {
    prefactor * edger_betacf(x, a, b) / a
  } else {
    // Use complement
    let prefactor2 = exp(a * log_nat(x) + b * log_nat(1.0 - x) - lbeta)
    1.0 - prefactor2 * edger_betacf(1.0 - x, b, a) / b
  }
}

///|
/// Continued fraction for incomplete beta function (Lentz's method).
fn edger_betacf(x : Double, a : Double, b : Double) -> Double {
  let max_iter = 200
  let eps = 0.0000000001
  let qab = a + b
  let qap = a + 1.0
  let qam = a - 1.0

  let mut c = 1.0
  let mut d = 1.0 - qab * x / qap
  if d.abs() < eps {
    d = eps
  }
  d = 1.0 / d
  let mut h = d

  let mut m = 1
  while m <= max_iter {
    let m2 = 2 * m
    let aa = m.to_double() * (m.to_double() - qam) * x / ((qam + m2.to_double()) * (a + m2.to_double()))
    d = 1.0 + aa * d
    if d.abs() < eps {
      d = eps
    }
    c = 1.0 + aa / c
    if c.abs() < eps {
      c = eps
    }
    d = 1.0 / d
    h = h * d * c

    let aa2 = -(a + m.to_double()) * (qab + m.to_double()) * x / ((a + m2.to_double()) * (qap + m2.to_double()))
    d = 1.0 + aa2 * d
    if d.abs() < eps {
      d = eps
    }
    c = 1.0 + aa2 / c
    if c.abs() < eps {
      c = eps
    }
    d = 1.0 / d
    let del = d * c
    h = h * del

    if (del - 1.0).abs() < eps {
      break
    }
    m = m + 1
  }

  h
}

///|
/// Sin function via Taylor series (Double has no sin method).
fn edger_sin(x : Double) -> Double {
  // Taylor series for sin(x), x reduced to [-pi, pi]
  let pi = 3.141592653589793
  let two_pi = 6.283185307179586
  let mut x_reduced = x
  while x_reduced > pi {
    x_reduced = x_reduced - two_pi
  }
  while x_reduced < -pi {
    x_reduced = x_reduced + two_pi
  }
  let x2 = x_reduced * x_reduced
  let mut term = x_reduced
  let mut sum = x_reduced
  let mut i = 1
  while i <= 15 {
    term = -term * x2 / ((2 * i).to_double() * (2 * i + 1).to_double())
    sum = sum + term
    i = i + 1
  }
  sum
}

///|
/// Log-gamma function using Lanczos approximation.
fn log_gamma(x : Double) -> Double {
  if x < 0.5 {
    // Use reflection formula: Gamma(x)Gamma(1-x) = pi/sin(pi*x)
    return log_nat(3.141592653589793 / edger_sin(3.141592653589793 * x)) - log_gamma(1.0 - x)
  }

  let g = 7.0
  // Lanczos coefficients
  let c = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
  ]

  let y = x - 1.0
  let tmp = y + g + 0.5
  let mut sum = c[0]
  let mut i = 1
  while i < 9 {
    sum = sum + c[i] / (y + i.to_double())
    i = i + 1
  }

  0.5 * log_nat(2.0 * 3.141592653589793) + (y + 0.5) * log_nat(tmp) - tmp + log_nat(sum)
}

///|
/// Benjamini-Hochberg FDR adjustment.
fn bh_adjust(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return []
  }

  // Get sorted indices
  let indices = sort_indices(p_values)
  let padj : Array[Double] = Array::make(n, 1.0)

  let mut min_padj = 1.0
  let mut k = n - 1
  while k >= 0 {
    let idx = indices[k]
    let rank = (k + 1).to_double()
    let adj = p_values[idx] * n.to_double() / rank
    if adj < min_padj {
      min_padj = adj
    }
    padj[idx] = if min_padj > 1.0 { 1.0 } else { min_padj }
    k = k - 1
  }

  padj
}

///|
/// Sort array and return sorted copy.
fn sort_array(arr : Array[Double]) -> Array[Double] {
  let result = arr.copy()
  let n = result.length()
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - 1 - i {
      if result[j] > result[j + 1] {
        let tmp = result[j]
        result[j] = result[j + 1]
        result[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

///|
/// Return indices that would sort the array in ascending order.
fn sort_indices(arr : Array[Double]) -> Array[Int] {
  let n = arr.length()
  let indices : Array[Int] = Array::new()
  let used : Array[Bool] = Array::make(n, false)

  let mut i = 0
  while i < n {
    let mut min_idx = -1
    let mut min_val = 1.0e300
    let mut j = 0
    while j < n {
      if !used[j] && arr[j] < min_val {
        min_val = arr[j]
        min_idx = j
      }
      j = j + 1
    }
    if min_idx >= 0 {
      indices.push(min_idx)
      used[min_idx] = true
    }
    i = i + 1
  }

  indices
}

///|
/// Camera (Competitive Gene Set Test).
/// Tests whether a gene set is differentially expressed relative to
/// genes not in the set.
///
/// Parameters:
/// - qlf_result: Output from glm_qlf_test
/// - gene_set: Array of gene names in the test set
/// - n_perm: Number of permutations (default: 1000)
///
/// Returns (mean_stat, p_value, direction).
pub fn camera(
  qlf_result : QLFResult,
  gene_set : Array[String],
  _n_perm? : Int = 1000,
) -> (Double, Double, String) {
  let n_genes = qlf_result.genes.length()

  // Identify genes in the set
  let set_indices : Array[Int] = Array::new()
  let set_stats : Array[Double] = Array::new()
  let all_stats : Array[Double] = Array::new()

  for i = 0; i < n_genes; i = i + 1 {
    all_stats.push(qlf_result.f_stat[i])
    if gene_set.contains(qlf_result.genes[i]) {
      set_indices.push(i)
      set_stats.push(qlf_result.f_stat[i])
    }
  }

  if set_stats.length() == 0 {
    return (0.0, 1.0, "none")
  }

  // Compute mean statistic for the gene set
  let set_mean = array_mean(set_stats)
  let set_var = array_variance(set_stats)
  let n_set = set_stats.length()

  // Compute mean and variance for genes NOT in the set
  let unset_stats : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    if !gene_set.contains(qlf_result.genes[i]) {
      unset_stats.push(qlf_result.f_stat[i])
    }
  }
  let unset_mean = array_mean(unset_stats)
  let unset_var = if unset_stats.length() > 1 {
    array_variance(unset_stats)
  } else {
    1.0
  }
  let n_unset = unset_stats.length()

  // Camera statistic: two-sample t-test comparing set vs unset
  let pooled_var = ((n_set - 1).to_double() * set_var +
    (n_unset - 1).to_double() * unset_var) /
    (n_set + n_unset - 2).to_double()

  let se = (pooled_var * (1.0 / n_set.to_double() + 1.0 / n_unset.to_double())).sqrt()
  let t_stat = if se > 0.0 {
    (set_mean - unset_mean) / se
  } else {
    0.0
  }

  // Two-sided p-value from t-distribution (approximate using normal)
  let p_val = 2.0 * (1.0 - edger_normal_cdf(t_stat.abs()))

  let direction = if t_stat > 0.0 { "up" } else { "down" }
  (t_stat, p_val, direction)
}

///|
/// Roast (Self-Contained Gene Set Test).
/// Tests whether any genes in the set are differentially expressed.
///
/// Parameters:
/// - qlf_result: Output from glm_qlf_test
/// - gene_set: Array of gene names in the test set
///
/// Returns (mean_stat, p_value, direction).
pub fn roast(
  qlf_result : QLFResult,
  gene_set : Array[String],
) -> (Double, Double, String) {
  let n_genes = qlf_result.genes.length()

  // Collect statistics for genes in the set
  let set_stats : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    if gene_set.contains(qlf_result.genes[i]) {
      set_stats.push(qlf_result.f_stat[i])
    }
  }

  if set_stats.length() == 0 {
    return (0.0, 1.0, "none")
  }

  // Mean F-statistic
  let mean_stat = array_mean(set_stats)

  // Use chi-square approximation: sum of F(1, df) ~ chi2(n) approximately
  // For a one-sided test, use the mean statistic
  let df = if qlf_result.df_denom > 0 { qlf_result.df_denom } else { 1 }
  let df_dbl = df.to_double()
  let denom = if df_dbl - 2.0 > 1.0 { df_dbl - 2.0 } else { 1.0 }
  let expected_mean = df_dbl / denom

  // Z-score
  let n_set = set_stats.length()
  let z = if n_set > 0 {
    (mean_stat - expected_mean) / (expected_mean * 2.0 / n_set.to_double()).sqrt()
  } else {
    0.0
  }

  // One-sided p-value
  let p_val = 1.0 - edger_normal_cdf(z)

  let direction = if z > 0.0 { "up" } else { "down" }
  (mean_stat, p_val, direction)
}

///|
/// Normal CDF (cumulative distribution function).
fn edger_normal_cdf(x : Double) -> Double {
  0.5 * (1.0 + erf(x / 1.4142135623730951))
}

///|
/// Error function (Abramowitz and Stegun approximation).
fn erf(x : Double) -> Double {
  let t = 1.0 / (1.0 + 0.3275911 * x.abs())
  let poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
  let result = 1.0 - poly * exp(-x * x)
  if x >= 0.0 { result } else { -result }
}

///|
/// Compute mean of an array.
fn array_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for v in arr {
    sum = sum + v
  }
  sum / arr.length().to_double()
}

///|
/// Compute variance of an array.
fn array_variance(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let mean = array_mean(arr)
  let mut sum_sq = 0.0
  for v in arr {
    let diff = v - mean
    sum_sq = sum_sq + diff * diff
  }
  sum_sq / (n - 1).to_double()
}

///|
/// Get significant genes from QL F-test results.
pub fn qlf_significant_genes(
  result : QLFResult,
  alpha? : Double = 0.05,
  lfc_threshold? : Double = 0.0,
) -> Array[String] {
  let genes : Array[String] = Array::new()
  for i = 0; i < result.genes.length(); i = i + 1 {
    if result.fdr[i] < alpha && result.log_fc[i].abs() > lfc_threshold {
      genes.push(result.genes[i])
    }
  }
  genes
}

///|
/// Get top genes from QL F-test results sorted by FDR.
pub fn qlf_top_genes(
  result : QLFResult,
  n? : Int = 10,
) -> Array[(String, Double, Double, Double)] {
  let n_genes = result.genes.length()
  let actual_n = if n < n_genes { n } else { n_genes }

  let indices = sort_indices(result.fdr)
  let top : Array[(String, Double, Double, Double)] = Array::new()
  for i = 0; i < actual_n; i = i + 1 {
    let idx = indices[i]
    top.push((result.genes[idx], result.log_fc[idx], result.f_stat[idx], result.fdr[idx]))
  }
  top
}
