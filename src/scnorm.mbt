///|
/// SCnorm module, inspired by Bioconductor SCnorm package.
/// Provides single-cell RNA-seq normalization using quantile regression
/// to correct for depth-dependent bias in single-cell data.

/// Normalization method
pub enum SCnormMethod {
  /// Quantile regression normalization
  Quantile
  /// LOESS-based normalization
  Loess
  /// Spline-based normalization
  Spline
}

/// Gene filtering method
pub enum SCnormFilterMethod {
  /// Filter by expression threshold
  Threshold
  /// Filter by percentage of cells expressing
  Percent
}

/// Result of quantile regression
pub struct SCnormQuantFit {
  /// Slope coefficient
  slope : Double
  /// Intercept
  intercept : Double
  /// Fitted values
  fitted : Array[Double]
  /// Residuals
  residuals : Array[Double]
  /// Quantile level
  tau : Double
}

/// Normalization result for a gene
pub struct SCnormGeneNormResult {
  /// Gene ID
  gene_id : String
  /// Quantile fit
  quant_fit : SCnormQuantFit
  /// Normalized counts
  norm_counts : Array[Double]
  /// Whether gene passed filter
  passed_filter : Bool
  /// Filter reason if not passed
  filter_reason : String
}

/// Complete SCnorm result
pub struct SCnormResult {
  /// Gene normalization results
  gene_results : Array[SCnormGeneNormResult]
  /// Normalized count matrix (genes x cells)
  norm_matrix : Array[Array[Double]]
  /// Gene IDs
  gene_ids : Array[String]
  /// Cell IDs
  cell_ids : Array[String]
  /// Library sizes (original)
  lib_sizes : Array[Double]
  /// Quantile level used
  quantile_level : Double
  /// Method used
  norm_method : SCnormMethod
}

/// Parameters for SCnorm normalization
pub struct SCnormParams {
  /// Quantile level for regression (0-1)
  quantile : Double
  /// Whether to filter genes
  filter_genes : Bool
  /// Minimum expression threshold
  min_expr : Double
  /// Minimum percent of cells expressing
  min_percent : Double
  /// Whether to use robust regression
  robust : Bool
}

/// Create default parameters
pub fn SCnormParams::new() -> SCnormParams {
  SCnormParams::{
    quantile: 0.5,
    filter_genes: true,
    min_expr: 0.1,
    min_percent: 0.1,
    robust: true
  }
}

// ===== Quantile Regression =====

///|
/// Fit a simple quantile regression model.
/// Uses the minimization approach for quantile loss.
pub fn sc_norm_quantile_regression(
  x : Array[Double],
  y : Array[Double],
  tau : Double
) -> SCnormQuantFit {
  let n = x.length()
  if n == 0 {
    return SCnormQuantFit::{
      slope: 0.0,
      intercept: 0.0,
      fitted: Array::new(),
      residuals: Array::new(),
      tau
    }
  }
  // Compute depth (x-axis) and log-counts (y-axis)
  let (slope, intercept) = sc_norm_simple_quantile_fit(x, y, tau)
  // Compute fitted values
  let fitted : Array[Double] = Array::new()
  let residuals : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let f = slope * x[i] + intercept
    fitted.push(f)
    residuals.push(y[i] - f)
    i = i + 1
  }
  SCnormQuantFit::{
    slope,
    intercept,
    fitted,
    residuals,
    tau
  }
}

///|
/// Simple quantile regression via grid search.
/// Finds slope and intercept that minimize quantile loss.
fn sc_norm_simple_quantile_fit(x : Array[Double], y : Array[Double], tau : Double) -> (Double, Double) {
  let n = x.length()
  if n < 2 { return (0.0, if n > 0 { y[0] } else { 0.0 }) }
  // Compute median-based estimates
  let x_median = sc_norm_median(x)
  let y_median = sc_norm_median(y)
  // Initial slope estimate via quantile-based ratio
  let qy_low = sc_norm_quantile(y, tau * 0.5)
  let qy_high = sc_norm_quantile(y, 1.0 - tau * 0.5)
  let qx_low = sc_norm_quantile(x, tau * 0.5)
  let qx_high = sc_norm_quantile(x, 1.0 - tau * 0.5)
  let slope_init = if qx_high > qx_low {
    (qy_high - qy_low) / (qx_high - qx_low)
  } else {
    // Fallback: use median ratio
    if x_median > 0.0 { y_median / x_median } else { 0.0 }
  }
  let intercept_init = y_median - slope_init * x_median
  // Refine via local search
  let mut best_slope = slope_init
  let mut best_intercept = intercept_init
  let mut best_loss = sc_norm_quantile_loss(x, y, best_slope, best_intercept, tau)
  // Grid search refinement
  let scale = if best_slope.abs() > 1.0 { best_slope.abs() } else { 1.0 }
  let search_range = 0.5 * scale
  let steps = 20
  let step_size = 2.0 * search_range / steps.to_double()
  let mut s = 0
  while s <= steps {
    let test_slope = best_slope - search_range + step_size * s.to_double()
    // For each slope, find best intercept
    let (best_int_for_slope, loss_for_slope) = sc_norm_find_best_intercept(x, y, test_slope, tau)
    if loss_for_slope < best_loss {
      best_loss = loss_for_slope
      best_slope = test_slope
      best_intercept = best_int_for_slope
    }
    s = s + 1
  }
  (best_slope, best_intercept)
}

///|
/// Find best intercept for a given slope.
fn sc_norm_find_best_intercept(x : Array[Double], y : Array[Double], slope : Double, tau : Double) -> (Double, Double) {
  let n = x.length()
  // Compute residuals at multiple intercepts
  // Use the median residual as starting point
  let residuals : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    residuals.push(y[i] - slope * x[i])
    i = i + 1
  }
  // Sort residuals to find the tau-th quantile
  let sorted = sc_norm_sort_array(residuals)
  let idx = (tau * (n.to_double() - 1.0)).to_int()
  let intercept = sorted[if idx < 0 { 0 } else { if idx >= n { n - 1 } else { idx } }]
  let loss = sc_norm_quantile_loss(x, y, slope, intercept, tau)
  (intercept, loss)
}

///|
/// Compute quantile loss (pinball loss).
fn sc_norm_quantile_loss(x : Array[Double], y : Array[Double], slope : Double, intercept : Double, tau : Double) -> Double {
  let n = x.length()
  let mut loss = 0.0
  let mut i = 0
  while i < n {
    let pred = slope * x[i] + intercept
    let err = y[i] - pred
    let pinball = if err >= 0.0 { tau * err } else { (tau - 1.0) * err }
    loss = loss + pinball
    i = i + 1
  }
  loss
}

/// Compute median of an array
fn sc_norm_median(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 { return 0.0 }
  let sorted = sc_norm_sort_array(arr)
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

/// Compute quantile of an array
fn sc_norm_quantile(arr : Array[Double], q : Double) -> Double {
  let n = arr.length()
  if n == 0 { return 0.0 }
  let sorted = sc_norm_sort_array(arr)
  let idx = (q * (n.to_double() - 1.0)).to_int()
  let idx_clamped = if idx < 0 { 0 } else { if idx >= n { n - 1 } else { idx } }
  sorted[idx_clamped]
}

/// Sort array (simple insertion sort for small arrays)
fn sc_norm_sort_array(arr : Array[Double]) -> Array[Double] {
  let result = arr.copy()
  let n = result.length()
  let mut i = 1
  while i < n {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 && result[j] > key {
      result[j + 1] = result[j]
      j = j - 1
    }
    result[j + 1] = key
    i = i + 1
  }
  result
}

// ===== Normalization =====

///|
/// Compute library sizes (total counts per cell).
/// Input should be genes x cells matrix.
pub fn sc_norm_library_sizes(counts : Array[Array[Double]]) -> Array[Double] {
  let n_genes = counts.length()
  let n_cells = if n_genes > 0 { counts[0].length() } else { 0 }
  let sizes : Array[Double] = Array::new()
  let mut c = 0
  while c < n_cells {
    let mut sum = 0.0
    let mut g = 0
    while g < n_genes {
      sum = sum + counts[g][c]
      g = g + 1
    }
    sizes.push(sum)
    c = c + 1
  }
  sizes
}

///|
/// Normalize a single gene across cells using quantile regression.
pub fn sc_norm_normalize_gene(
  gene_id : String,
  counts : Array[Double],
  lib_sizes : Array[Double],
  params : SCnormParams
) -> SCnormGeneNormResult {
  let n = counts.length()
  // Compute log-transformed counts
  let log_counts : Array[Double] = Array::new()
  let log_lib_sizes : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let c = if counts[i] > 0.0 { counts[i] } else { 0.5 }
    let l = if lib_sizes[i] > 0.0 { lib_sizes[i] } else { 1.0 }
    log_counts.push(@math.ln(c))
    log_lib_sizes.push(@math.ln(l))
    i = i + 1
  }
  // Fit quantile regression
  let quant_fit = sc_norm_quantile_regression(log_lib_sizes, log_counts, params.quantile)
  // Compute normalized counts
  let norm_counts : Array[Double] = Array::new()
  let mut j = 0
  while j < n {
    let predicted = quant_fit.slope * log_lib_sizes[j] + quant_fit.intercept
    let residual = log_counts[j] - predicted
    // Exponentiate back
    let norm = @math.exp(residual)
    norm_counts.push(norm)
    j = j + 1
  }
  // Filter check
  let (passed, reason) = sc_norm_check_filter(counts, norm_counts, params)
  SCnormGeneNormResult::{
    gene_id,
    quant_fit,
    norm_counts,
    passed_filter: passed,
    filter_reason: reason
  }
}

///|
/// Check if a gene passes filter criteria.
fn sc_norm_check_filter(raw_counts : Array[Double], norm_counts : Array[Double], params : SCnormParams) -> (Bool, String) {
  if !params.filter_genes { return (true, "") }
  let n = raw_counts.length()
  let mut n_expressed = 0
  let mut i = 0
  while i < n {
    if raw_counts[i] > params.min_expr { n_expressed = n_expressed + 1 }
    i = i + 1
  }
  let percent = n_expressed.to_double() / n.to_double()
  if percent < params.min_percent {
    return (false, "Not expressed in enough cells (" + n_expressed.to_string() + "/" + n.to_string() + ")")
  }
  (true, "")
}

///|
/// Run the complete SCnorm normalization pipeline.
pub fn sc_norm_run(
  count_matrix : Array[Array[Double]],
  gene_ids : Array[String],
  cell_ids : Array[String],
  params : SCnormParams
) -> SCnormResult {
  let n_genes = count_matrix.length()
  let n_cells = if n_genes > 0 { count_matrix[0].length() } else { 0 }
  // Compute library sizes
  let lib_sizes = sc_norm_library_sizes(count_matrix)
  // Normalize each gene
  let gene_results : Array[SCnormGeneNormResult] = Array::new()
  let norm_matrix : Array[Array[Double]] = Array::new()
  let result_gene_ids : Array[String] = Array::new()
  let mut g = 0
  while g < n_genes {
    let gene_counts = count_matrix[g]
    let gene_id = if g < gene_ids.length() { gene_ids[g] } else { "gene_" + g.to_string() }
    let result = sc_norm_normalize_gene(gene_id, gene_counts, lib_sizes, params)
    gene_results.push(result)
    if result.passed_filter {
      norm_matrix.push(result.norm_counts)
      result_gene_ids.push(gene_id)
    }
    g = g + 1
  }
  SCnormResult::{
    gene_results,
    norm_matrix,
    gene_ids: result_gene_ids,
    cell_ids,
    lib_sizes,
    quantile_level: params.quantile,
    norm_method: if params.robust { SCnormMethod::Quantile } else { SCnormMethod::Loess }
  }
}

// ===== LOESS-based smoothing =====

///|
/// Simple LOESS (locally weighted smoothing) for 1D data.
pub fn sc_norm_loess_smooth(x : Array[Double], y : Array[Double], span : Double) -> Array[Double] {
  let n = x.length()
  if n < 3 { return y.copy() }
  let smoothed : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    // Find neighbors within span
    let neighbors : Array[(Double, Double)] = Array::new()
    let xi = x[i]
    let half_span = span * 0.5
    let mut j = 0
    while j < n {
      let dist = (x[j] - xi).abs()
      if dist <= half_span {
        neighbors.push((x[j], y[j]))
      }
      j = j + 1
    }
    // If not enough neighbors, use all points
    let n_neighbors = neighbors.length()
    let eff_neighbors = if n_neighbors >= 2 { n_neighbors } else { n }
    // Weighted least squares fit
    let w : Array[Double] = Array::new()
    let xs : Array[Double] = Array::new()
    let ys : Array[Double] = Array::new()
    let mut k = 0
    while k < eff_neighbors {
      let idx = if k < n_neighbors { k } else { k }
      let (xn, yn) = if k < n_neighbors { neighbors[k] } else { (x[k], y[k]) }
      let dist = (xn - xi).abs()
      let weight = if half_span > 0.0 {
        let u = dist / half_span
        if u <= 1.0 { (1.0 - u * u) * (1.0 - u * u) } else { 0.0 }
      } else { 1.0 }
      w.push(weight)
      xs.push(xn)
      ys.push(yn)
      k = k + 1
    }
    // Simple weighted average as approximation
    let mut wsum = 0.0
    let mut ysum = 0.0
    let mut m = 0
    while m < eff_neighbors {
      wsum = wsum + w[m]
      ysum = ysum + w[m] * ys[m]
      m = m + 1
    }
    let val = if wsum > 0.0 { ysum / wsum } else { y[i] }
    smoothed.push(val)
    i = i + 1
  }
  smoothed
}

// ===== Spline-based =====

///|
/// Simple spline smoothing for 1D data.
pub fn sc_norm_spline_smooth(x : Array[Double], y : Array[Double], n_knots : Int) -> Array[Double] {
  let n = x.length()
  if n < 4 { return y.copy() }
  let sorted_x = sc_norm_sort_indices(x)
  // Create knot positions
  let x_min = sc_norm_min(x)
  let x_max = sc_norm_max(x)
  let knot_span = if n_knots > 2 { (x_max - x_min) / (n_knots - 1).to_double() } else { 1.0 }
  let knots : Array[Double] = Array::new()
  let mut k = 0
  while k < n_knots {
    knots.push(x_min + knot_span * k.to_double())
    k = k + 1
  }
  // Fit using basis function approach (simplified)
  let smoothed : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let xi = x[i]
    // Find nearby knots and interpolate
    let val = sc_norm_spline_interpolate(xi, knots, x, y)
    smoothed.push(val)
    i = i + 1
  }
  smoothed
}

/// Spline interpolation at a point
fn sc_norm_spline_interpolate(x_val : Double, knots : Array[Double], x : Array[Double], y : Array[Double]) -> Double {
  let n_knots = knots.length()
  let n = x.length()
  // Find the two nearest data points
  let mut best_idx = 0
  let mut best_dist = 1.0e30
  let mut i = 0
  while i < n {
    let d = (x[i] - x_val).abs()
    if d < best_dist {
      best_dist = d
      best_idx = i
    }
    i = i + 1
  }
  // Simple linear interpolation between nearest neighbors
  if best_idx == 0 { return y[0] }
  if best_idx == n - 1 { return y[n - 1] }
  let x0 = x[best_idx - 1]
  let x1 = x[best_idx]
  let y0 = y[best_idx - 1]
  let y1 = y[best_idx]
  let t = if x1 > x0 { (x_val - x0) / (x1 - x0) } else { 0.5 }
  y0 + t * (y1 - y0)
}

fn sc_norm_sort_indices(arr : Array[Double]) -> Array[Double] {
  sc_norm_sort_array(arr)
}

fn sc_norm_min(arr : Array[Double]) -> Double {
  let mut m = arr[0]
  let mut i = 1
  while i < arr.length() {
    if arr[i] < m { m = arr[i] }
    i = i + 1
  }
  m
}

fn sc_norm_max(arr : Array[Double]) -> Double {
  let mut m = arr[0]
  let mut i = 1
  while i < arr.length() {
    if arr[i] > m { m = arr[i] }
    i = i + 1
  }
  m
}

// ===== Accessor methods =====

/// Get number of genes after filtering
pub fn SCnormResult::sc_n_n_genes(self : SCnormResult) -> Int {
  self.gene_ids.length()
}

/// Get number of cells
pub fn SCnormResult::sc_n_cells(self : SCnormResult) -> Int {
  self.cell_ids.length()
}

/// Get a normalized count for a gene and cell
pub fn SCnormResult::sc_norm_get_count(self : SCnormResult, gene_idx : Int, cell_idx : Int) -> Double {
  self.norm_matrix[gene_idx][cell_idx]
}

/// Get all normalized counts for a gene
pub fn SCnormResult::sc_norm_get_gene(self : SCnormResult, gene_idx : Int) -> Array[Double] {
  self.norm_matrix[gene_idx].copy()
}

/// Get all normalized counts for a cell
pub fn SCnormResult::sc_norm_get_cell(self : SCnormResult, cell_idx : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < self.norm_matrix.length() {
    result.push(self.norm_matrix[i][cell_idx])
    i = i + 1
  }
  result
}

/// Get the quantile fit for a gene
pub fn SCnormResult::sc_norm_get_fit(self : SCnormResult, gene_id : String) -> SCnormQuantFit? {
  let mut i = 0
  while i < self.gene_results.length() {
    if self.gene_results[i].gene_id == gene_id {
      return Some(self.gene_results[i].quant_fit)
    }
    i = i + 1
  }
  None
}
