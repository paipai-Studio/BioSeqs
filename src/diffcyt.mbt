///|
/// diffcyt - Differential discovery in high-dimensional cytometry via
/// high-resolution clustering. Inspired by Bioconductor diffcyt package.
///
/// Provides FlowSOM-style clustering, differential abundance (DA) testing
/// via negative binomial GLM, and differential state (DS) testing via
/// moderated t-test, with BH-FDR multiple testing correction.
///
/// References:
///   - Weber LM & Robinson MD (2016) diffcyt: Differential discovery in
///     high-dimensional cytometry via high-resolution clustering. Nat Methods 17:792.
///   - Van Gassen S et al. (2015) FlowSOM: Using self-organizing maps for
///     visualization and interpretation of cytometry data. Cytometry A 87:636.

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// A single cytometry cell with marker expression values.
pub struct CytometryCell {
  cell_id : String
  sample_id : String
  condition : String
  marker_values : Array[Double]
  mut cluster_id : Int
}

///|
/// Create a new CytometryCell.
pub fn CytometryCell::new(
  cell_id : String,
  sample_id : String,
  condition : String,
  marker_values : Array[Double]
) -> CytometryCell {
  CytometryCell::{
    cell_id,
    sample_id,
    condition,
    marker_values,
    cluster_id: -1,
  }
}

///|
/// Accessor for cell id.
pub fn CytometryCell::cell_id(self : CytometryCell) -> String {
  self.cell_id
}

///|
/// Accessor for sample id.
pub fn CytometryCell::sample_id(self : CytometryCell) -> String {
  self.sample_id
}

///|
/// Accessor for condition.
pub fn CytometryCell::condition(self : CytometryCell) -> String {
  self.condition
}

///|
/// Accessor for marker values.
pub fn CytometryCell::marker_values(self : CytometryCell) -> Array[Double] {
  self.marker_values
}

///|
/// Accessor for cluster id.
pub fn CytometryCell::cluster_id(self : CytometryCell) -> Int {
  self.cluster_id
}

///|
/// Set the cluster id.
pub fn CytometryCell::set_cluster_id(self : CytometryCell, id : Int) -> Unit {
  self.cluster_id = id
}

///|
/// Number of markers.
pub fn CytometryCell::n_markers(self : CytometryCell) -> Int {
  self.marker_values.length()
}

///|
/// Result of differential abundance (DA) testing for one cluster.
pub struct DiffcytDAResult {
  cluster_id : Int
  log_fc : Double
  p_value : Double
  mut fdr : Double
  mut significant : Bool
}

///|
/// Accessor for cluster id.
pub fn DiffcytDAResult::cluster_id(self : DiffcytDAResult) -> Int {
  self.cluster_id
}

///|
/// Accessor for log fold change.
pub fn DiffcytDAResult::log_fc(self : DiffcytDAResult) -> Double {
  self.log_fc
}

///|
/// Accessor for p-value.
pub fn DiffcytDAResult::p_value(self : DiffcytDAResult) -> Double {
  self.p_value
}

///|
/// Accessor for FDR.
pub fn DiffcytDAResult::fdr(self : DiffcytDAResult) -> Double {
  self.fdr
}

///|
/// Accessor for significance flag.
pub fn DiffcytDAResult::significant(self : DiffcytDAResult) -> Bool {
  self.significant
}

///|
/// Result of differential state (DS) testing for one cluster-marker pair.
pub struct DiffcytDSResult {
  cluster_id : Int
  marker_idx : Int
  marker_name : String
  log_fc : Double
  p_value : Double
  mut fdr : Double
  mut significant : Bool
}

///|
/// Accessor for cluster id.
pub fn DiffcytDSResult::cluster_id(self : DiffcytDSResult) -> Int {
  self.cluster_id
}

///|
/// Accessor for marker index.
pub fn DiffcytDSResult::marker_idx(self : DiffcytDSResult) -> Int {
  self.marker_idx
}

///|
/// Accessor for marker name.
pub fn DiffcytDSResult::marker_name(self : DiffcytDSResult) -> String {
  self.marker_name
}

///|
/// Accessor for log fold change.
pub fn DiffcytDSResult::log_fc(self : DiffcytDSResult) -> Double {
  self.log_fc
}

///|
/// Accessor for p-value.
pub fn DiffcytDSResult::p_value(self : DiffcytDSResult) -> Double {
  self.p_value
}

///|
/// Accessor for FDR.
pub fn DiffcytDSResult::fdr(self : DiffcytDSResult) -> Double {
  self.fdr
}

///|
/// Accessor for significance flag.
pub fn DiffcytDSResult::significant(self : DiffcytDSResult) -> Bool {
  self.significant
}

///|
/// Summary of a diffcyt analysis run.
pub struct DiffcytSummary {
  n_cells : Int
  n_clusters : Int
  n_markers : Int
  n_samples : Int
  n_conditions : Int
  da_results : Array[DiffcytDAResult]
  ds_results : Array[DiffcytDSResult]
  fdr_threshold : Double
}

///|
/// Number of significant DA results.
pub fn DiffcytSummary::n_sig_da(self : DiffcytSummary) -> Int {
  let mut n = 0
  for r in self.da_results {
    if r.significant {
      n = n + 1
    }
  }
  n
}

///|
/// Number of significant DS results.
pub fn DiffcytSummary::n_sig_ds(self : DiffcytSummary) -> Int {
  let mut n = 0
  for r in self.ds_results {
    if r.significant {
      n = n + 1
    }
  }
  n
}

// ===========================================================================
// Internal helpers
// ===========================================================================

///|
/// Euclidean distance between two vectors.
fn dc_distance(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  let mut s = 0.0
  for i in 0..<n {
    let d = a[i] - b[i]
    s = s + d * d
  }
  dc_sqrt(s)
}

///|
/// Square root via Newton's method.
fn dc_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = x / 2.0
  if guess < 1.0 {
    guess = 1.0
  }
  for _ in 0..<50 {
    let next = (guess + x / guess) / 2.0
    let diff = next - guess
    if diff < 0.0 {
      if -diff < 0.0000000001 {
        return next
      }
    } else if diff < 0.0000000001 {
      return next
    }
    guess = next
  }
  guess
}

///|
/// Safe natural log.
fn dc_log(x : Double) -> Double {
  if x <= 0.0 {
    return -700.0
  }
  @math.ln(x)
}

///|
/// Logistic sigmoid.
fn dc_sigmoid(x : Double) -> Double {
  if x >= 40.0 {
    return 1.0
  }
  if x <= -40.0 {
    return 0.0000000001
  }
  1.0 / (1.0 + @math.exp(-x))
}

///|
/// Normal CDF approximation (Abramowitz & Stegun).
fn dc_normal_cdf(x : Double) -> Double {
  if x < -8.0 {
    return 0.0
  }
  if x > 8.0 {
    return 1.0
  }
  let t = 1.0 / (1.0 + 0.2316419 * dc_abs(x))
  let d = 0.3989423 * @math.exp(-x * x / 2.0)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
  if x >= 0.0 {
    1.0 - p
  } else {
    p
  }
}

///|
/// Absolute value.
fn dc_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

///|
/// LCG pseudo-random number generator for reproducible clustering.
fn dc_lcg_next(state : Int) -> Int {
  let s = (state * 1103515245 + 12345) & 0x7FFFFFFF
  s
}

///|
/// Random double in [0, 1) from LCG state.
fn dc_lcg_rand(state : Int) -> (Int, Double) {
  let s = dc_lcg_next(state)
  (s, (s % 1000000).to_double() / 1000000.0)
}

// ===========================================================================
// FlowSOM-style clustering
// ===========================================================================

///|
/// Cluster cytometry cells using a simplified FlowSOM approach.
/// Uses online SOM-like updates: each cell is assigned to the nearest
/// codebook vector, which is then moved toward the cell.
///
/// `cells` is the array of cytometry cells. `n_clusters` is the desired
/// number of clusters. `n_iterations` controls training epochs.
/// Returns the array of cluster codebook vectors (centroids).
pub fn diffcyt_cluster_cells(
  cells : Array[CytometryCell],
  n_clusters : Int,
  n_iterations : Int
) -> Array[Array[Double]] {
  let n_cells = cells.length()
  if n_cells == 0 || n_clusters <= 0 {
    return []
  }
  let n_markers = cells[0].marker_values.length()
  // Initialize codebook vectors by sampling random cells.
  let mut rng_state = 42
  let codebooks : Array[Array[Double]] = Array::new()
  for _ in 0..<n_clusters {
    let (s, r) = dc_lcg_rand(rng_state)
    rng_state = s
    let idx = (r * n_cells.to_double()).to_int()
    let idx = if idx >= n_cells { n_cells - 1 } else { idx }
    let cb = Array::make(n_markers, 0.0)
    for j in 0..<n_markers {
      cb[j] = cells[idx].marker_values[j]
    }
    codebooks.push(cb)
  }
  // Training: online SOM updates with decreasing learning rate.
  for iter in 0..<n_iterations {
    let lr = 0.5 * (1.0 - iter.to_double() / n_iterations.to_double())
    let lr = if lr < 0.01 { 0.01 } else { lr }
    for ci in 0..<n_cells {
      // Find best matching unit (BMU).
      let cell_vals = cells[ci].marker_values
      let mut best_idx = 0
      let mut best_dist = dc_distance(cell_vals, codebooks[0])
      for k in 1..<n_clusters {
        let d = dc_distance(cell_vals, codebooks[k])
        if d < best_dist {
          best_dist = d
          best_idx = k
        }
      }
      // Update BMU toward the cell.
      for j in 0..<n_markers {
        codebooks[best_idx][j] = codebooks[best_idx][j] + lr * (cell_vals[j] - codebooks[best_idx][j])
      }
    }
  }
  // Assign final cluster labels.
  for ci in 0..<n_cells {
    let cell_vals = cells[ci].marker_values
    let mut best_idx = 0
    let mut best_dist = dc_distance(cell_vals, codebooks[0])
    for k in 1..<n_clusters {
      let d = dc_distance(cell_vals, codebooks[k])
      if d < best_dist {
        best_dist = d
        best_idx = k
      }
    }
    cells[ci].set_cluster_id(best_idx)
  }
  codebooks
}

///|
/// Assign cells to existing codebook vectors (for new data).
pub fn diffcyt_assign_clusters(
  cells : Array[CytometryCell],
  codebooks : Array[Array[Double]]
) -> Unit {
  let n_clusters = codebooks.length()
  for ci in 0..<cells.length() {
    let cell_vals = cells[ci].marker_values
    let mut best_idx = 0
    let mut best_dist = dc_distance(cell_vals, codebooks[0])
    for k in 1..<n_clusters {
      let d = dc_distance(cell_vals, codebooks[k])
      if d < best_dist {
        best_dist = d
        best_idx = k
      }
    }
    cells[ci].set_cluster_id(best_idx)
  }
  ignore(n_clusters)
}

// ===========================================================================
// Cluster counts and medians
// ===========================================================================

///|
/// Compute per-sample-per-cluster cell counts.
/// Returns a map: (sample_id, cluster_id) -> count.
pub fn diffcyt_calc_counts(
  cells : Array[CytometryCell],
  n_clusters : Int,
  sample_ids : Array[String]
) -> Array[Array[Int]] {
  let n_samples = sample_ids.length()
  // counts[sample_idx][cluster_idx]
  let counts : Array[Array[Int]] = Array::new()
  for _ in 0..<n_samples {
    counts.push(Array::make(n_clusters, 0))
  }
  // Build sample index lookup.
  let sample_idx_map : Map[String, Int] = Map([], capacity=n_samples)
  for i in 0..<n_samples {
    sample_idx_map.set(sample_ids[i], i)
  }
  for cell in cells {
    let s_idx = sample_idx_map.get(cell.sample_id)
    if s_idx is Some(si) {
      let c = cell.cluster_id
      if c >= 0 && c < n_clusters {
        counts[si][c] = counts[si][c] + 1
      }
    }
  }
  counts
}

///|
/// Compute per-cluster-per-marker median expression values,
/// aggregated by sample. Returns medians[sample_idx][cluster_idx * n_markers + marker_idx].
pub fn diffcyt_calc_medians_by_cluster_marker(
  cells : Array[CytometryCell],
  n_clusters : Int,
  n_markers : Int,
  sample_ids : Array[String]
) -> Array[Array[Double]] {
  let n_samples = sample_ids.length()
  // Collect expression values per (sample, cluster, marker).
  // Use a nested structure: values[sample][cluster * n_markers + marker] = Array[Double]
  let total_entries = n_clusters * n_markers
  let all_values : Array[Array[Array[Double]]] = Array::new()
  for _ in 0..<n_samples {
    let per_sample : Array[Array[Double]] = Array::new()
    for _ in 0..<total_entries {
      per_sample.push(Array::new())
    }
    all_values.push(per_sample)
  }
  let sample_idx_map : Map[String, Int] = Map([], capacity=n_samples)
  for i in 0..<n_samples {
    sample_idx_map.set(sample_ids[i], i)
  }
  for cell in cells {
    let s_opt = sample_idx_map.get(cell.sample_id)
    if s_opt is Some(si) {
      let c = cell.cluster_id
      if c >= 0 && c < n_clusters {
        for m in 0..<n_markers {
          all_values[si][c * n_markers + m].push(cell.marker_values[m])
        }
      }
    }
  }
  // Compute medians.
  let medians : Array[Array[Double]] = Array::new()
  for si in 0..<n_samples {
    let row = Array::make(total_entries, 0.0)
    for j in 0..<total_entries {
      row[j] = dc_median(all_values[si][j])
    }
    medians.push(row)
  }
  medians
}

///|
/// Compute median of an array of doubles.
fn dc_median(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  // Copy and sort.
  let sorted = Array::new()
  for v in arr {
    sorted.push(v)
  }
  // Simple insertion sort (arrays are small per cluster-marker).
  for i in 1..<n {
    let key = sorted[i]
    let mut j = i - 1
    while j >= 0 && sorted[j] > key {
      sorted[j + 1] = sorted[j]
      j = j - 1
    }
    sorted[j + 1] = key
  }
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

// ===========================================================================
// Design matrix and contrast
// ===========================================================================

///|
/// Build a design matrix for differential testing.
/// Condition 0 is the reference (intercept only), other conditions get
/// indicator columns. Returns design[sample_idx][col_idx].
pub fn diffcyt_create_design_matrix(
  sample_ids : Array[String],
  conditions : Array[String]
) -> Array[Array[Double]] {
  let n = sample_ids.length()
  // Collect unique conditions preserving order.
  let unique_conditions : Array[String] = Array::new()
  for c in conditions {
    if !unique_conditions.contains(c) {
      unique_conditions.push(c)
    }
  }
  let n_cols = unique_conditions.length() // intercept + (n_conditions - 1) dummies
  let design : Array[Array[Double]] = Array::new()
  for i in 0..<n {
    let row = Array::make(n_cols, 0.0)
    row[0] = 1.0 // intercept
    // Find condition index.
    let cond = conditions[i]
    for ci in 0..<unique_conditions.length() {
      if unique_conditions[ci] == cond {
        row[ci] = 1.0
      }
    }
    design.push(row)
  }
  design
}

///|
/// Get unique conditions from cells.
pub fn diffcyt_unique_conditions(cells : Array[CytometryCell]) -> Array[String] {
  let result : Array[String] = Array::new()
  for cell in cells {
    if !result.contains(cell.condition) {
      result.push(cell.condition)
    }
  }
  result
}

///|
/// Get unique sample IDs from cells.
pub fn diffcyt_unique_samples(cells : Array[CytometryCell]) -> Array[String] {
  let result : Array[String] = Array::new()
  for cell in cells {
    if !result.contains(cell.sample_id) {
      result.push(cell.sample_id)
    }
  }
  result
}

// ===========================================================================
// Differential Abundance (DA) testing
// ===========================================================================

///|
/// Test for differential abundance using a simplified negative binomial
/// GLM approach. For each cluster, fits a log-linear model of counts
/// against condition, and tests for difference using a Wald test.
///
/// `counts` is [sample_idx][cluster_idx] cell counts.
/// `design` is the design matrix [sample_idx][col_idx].
/// `contrast_col` is the column index of the condition to test.
/// Returns DA results for each cluster.
pub fn diffcyt_testDA(
  counts : Array[Array[Int]],
  design : Array[Array[Double]],
  contrast_col : Int,
  n_clusters : Int
) -> Array[DiffcytDAResult] {
  let n_samples = counts.length()
  let results : Array[DiffcytDAResult] = Array::new()
  // Compute library sizes (total counts per sample).
  let lib_sizes = Array::make(n_samples, 0.0)
  for si in 0..<n_samples {
    let mut s = 0.0
    for ci in 0..<n_clusters {
      s = s + counts[si][ci].to_double()
    }
    lib_sizes[si] = s
  }
  let mean_lib = dc_mean(lib_sizes)
  for cluster in 0..<n_clusters {
    // Extract counts for this cluster.
    let y = Array::make(n_samples, 0.0)
    for si in 0..<n_samples {
      y[si] = counts[si][cluster].to_double()
    }
    // Offset: log(library size / mean library size).
    let offsets = Array::make(n_samples, 0.0)
    for si in 0..<n_samples {
      offsets[si] = dc_log(lib_sizes[si] / mean_lib)
    }
    // Fit NB GLM via IRLS (simplified: Poisson-like with dispersion).
    let n_params = design[0].length()
    let beta = Array::make(n_params, 0.0)
    // Estimate dispersion from mean.
    let mean_count = dc_mean(y)
    let dispersion = if mean_count > 0.0 {
      0.1 + 1.0 / mean_count // simplified dispersion estimate
    } else {
      1.0
    }
    // IRLS iterations.
    for _ in 0..<50 {
      let mu = Array::make(n_samples, 0.0)
      for si in 0..<n_samples {
        let mut eta = offsets[si]
        for j in 0..<n_params {
          eta = eta + design[si][j] * beta[j]
        }
        mu[si] = @math.exp(eta)
        if mu[si] < 0.0000001 {
          mu[si] = 0.0000001
        }
      }
      // Weights and working response.
      let w = Array::make(n_samples, 0.0)
      let z = Array::make(n_samples, 0.0)
      for si in 0..<n_samples {
        let var = mu[si] + dispersion * mu[si] * mu[si]
        w[si] = mu[si] / var
        let eta = dc_log(mu[si]) + (y[si] - mu[si]) / mu[si]
        z[si] = eta - offsets[si]
      }
      // Weighted least squares: solve (X^T W X) beta = X^T W z.
      let xtwx = Array::new()
      for _ in 0..<n_params {
        xtwx.push(Array::make(n_params, 0.0))
      }
      let xtwz = Array::make(n_params, 0.0)
      for si in 0..<n_samples {
        for j in 0..<n_params {
          xtwz[j] = xtwz[j] + design[si][j] * w[si] * z[si]
          for k in 0..<n_params {
            xtwx[j][k] = xtwx[j][k] + design[si][j] * w[si] * design[si][k]
          }
        }
      }
      match dc_solve(xtwx, xtwz) {
        Some(new_beta) => {
          let mut converged = true
          for j in 0..<n_params {
            let diff = new_beta[j] - beta[j]
            if dc_abs(diff) > 0.0001 {
              converged = false
            }
            beta[j] = new_beta[j]
          }
          if converged {
            break
          }
        }
        None => break
      }
    }
    // Compute standard error for the contrast coefficient.
    let mu = Array::make(n_samples, 0.0)
    for si in 0..<n_samples {
      let mut eta = offsets[si]
      for j in 0..<n_params {
        eta = eta + design[si][j] * beta[j]
      }
      mu[si] = @math.exp(eta)
      if mu[si] < 0.0000001 {
        mu[si] = 0.0000001
      }
    }
    let w = Array::make(n_samples, 0.0)
    for si in 0..<n_samples {
      let var = mu[si] + dispersion * mu[si] * mu[si]
      w[si] = mu[si] / var
    }
    let xtwx = Array::new()
    for _ in 0..<n_params {
      xtwx.push(Array::make(n_params, 0.0))
    }
    for si in 0..<n_samples {
      for j in 0..<n_params {
        for k in 0..<n_params {
          xtwx[j][k] = xtwx[j][k] + design[si][j] * w[si] * design[si][k]
        }
      }
    }
    // Invert (X^T W X) to get covariance; extract diagonal for SE.
    let cov_inv = dc_matrix_inverse(xtwx)
    let se = if cov_inv is Some(ci) {
      dc_sqrt(dc_abs(ci[contrast_col][contrast_col]))
    } else {
      1.0
    }
    let se = if se < 0.0000001 { 0.0000001 } else { se }
    // Wald test.
    let wald_stat = beta[contrast_col] / se
    let p_value = 2.0 * (1.0 - dc_normal_cdf(dc_abs(wald_stat)))
    let p_value = if p_value < 0.0 { 0.0 } else if p_value > 1.0 { 1.0 } else { p_value }
    results.push(DiffcytDAResult::{
      cluster_id: cluster,
      log_fc: beta[contrast_col],
      p_value,
      fdr: p_value, // will be corrected later
      significant: false,
    })
  }
  // BH-FDR correction.
  dc_bh_fdr_da(results)
}

// ===========================================================================
// Differential State (DS) testing
// ===========================================================================

///|
/// Test for differential state (marker expression differences) using a
/// simplified moderated t-test approach. For each cluster-marker pair,
/// compares expression between conditions using a Welch-like t-test
/// with empirical Bayes moderation of variance.
///
/// `cells` is the array of clustered cells.
/// `medians` is the per-sample per-cluster-marker medians.
/// `sample_ids` and `conditions` describe sample metadata.
/// `n_clusters` and `n_markers` define dimensions.
/// `marker_names` provides marker names for results.
pub fn diffcyt_testDS(
  cells : Array[CytometryCell],
  sample_ids : Array[String],
  conditions : Array[String],
  n_clusters : Int,
  n_markers : Int,
  marker_names : Array[String]
) -> Array[DiffcytDSResult] {
  let n_samples = sample_ids.length()
  // Collect per-sample-per-cluster-per-marker median expression.
  let medians = diffcyt_calc_medians_by_cluster_marker(
    cells, n_clusters, n_markers, sample_ids
  )
  // Collect unique conditions.
  let unique_conditions : Array[String] = Array::new()
  for c in conditions {
    if !unique_conditions.contains(c) {
      unique_conditions.push(c)
    }
  }
  let n_conditions = unique_conditions.length()
  let results : Array[DiffcytDSResult] = Array::new()
  // Gather all sample variances for empirical Bayes prior.
  let all_variances : Array[Double] = Array::new()
  let all_df : Array[Double] = Array::new()
  // Pre-compute per-cluster-per-marker group means and variances.
  for cluster in 0..<n_clusters {
    for marker in 0..<n_markers {
      // Collect median values per condition.
      let group_values : Array[Array[Double]] = Array::new()
      for _ in 0..<n_conditions {
        group_values.push(Array::new())
      }
      for si in 0..<n_samples {
        let val = medians[si][cluster * n_markers + marker]
        for ci in 0..<n_conditions {
          if unique_conditions[ci] == conditions[si] {
            group_values[ci].push(val)
          }
        }
      }
      // Compute group means and pooled variance.
      let group_means = Array::make(n_conditions, 0.0)
      let group_vars = Array::make(n_conditions, 0.0)
      let group_ns = Array::make(n_conditions, 0)
      for ci in 0..<n_conditions {
        let gv = group_values[ci]
        group_ns[ci] = gv.length()
        group_means[ci] = dc_mean(gv)
        group_vars[ci] = dc_variance(gv)
      }
      // Pooled variance (Welch-like).
      let mut pooled_var = 0.0
      let mut total_df = 0.0
      for ci in 0..<n_conditions {
        if group_ns[ci] > 1 {
          pooled_var = pooled_var + group_vars[ci] * (group_ns[ci] - 1).to_double()
          total_df = total_df + (group_ns[ci] - 1).to_double()
        }
      }
      let pooled_var = if total_df > 0.0 {
        pooled_var / total_df
      } else {
        1.0
      }
      if pooled_var > 0.0 {
        all_variances.push(pooled_var)
        all_df.push(total_df)
      }
    }
  }
  // Empirical Bayes prior: estimate prior mean and variance of log-variances.
  let prior_mean_var = if all_variances.length() > 0 {
    dc_mean(all_variances)
  } else {
    1.0
  }
  let prior_mean_var = if prior_mean_var < 0.0000001 { 0.0000001 } else { prior_mean_var }
  // Moderate variances: shrink toward prior.
  let prior_d0 = 10.0 // prior degrees of freedom
  for cluster in 0..<n_clusters {
    for marker in 0..<n_markers {
      let group_values : Array[Array[Double]] = Array::new()
      for _ in 0..<n_conditions {
        group_values.push(Array::new())
      }
      for si in 0..<n_samples {
        let val = medians[si][cluster * n_markers + marker]
        for ci in 0..<n_conditions {
          if unique_conditions[ci] == conditions[si] {
            group_values[ci].push(val)
          }
        }
      }
      let group_means = Array::make(n_conditions, 0.0)
      let group_vars = Array::make(n_conditions, 0.0)
      let group_ns = Array::make(n_conditions, 0)
      for ci in 0..<n_conditions {
        let gv = group_values[ci]
        group_ns[ci] = gv.length()
        group_means[ci] = dc_mean(gv)
        group_vars[ci] = dc_variance(gv)
      }
      // Pooled variance with moderation.
      let mut sse = 0.0
      let mut d0 = 0.0
      for ci in 0..<n_conditions {
        if group_ns[ci] > 1 {
          sse = sse + group_vars[ci] * (group_ns[ci] - 1).to_double()
          d0 = d0 + (group_ns[ci] - 1).to_double()
        }
      }
      let s2 = if d0 > 0.0 { sse / d0 } else { prior_mean_var }
      // Moderated variance: (d0 * prior + d * s2) / (d0 + d)
      let mod_var = (prior_d0 * prior_mean_var + d0 * s2) / (prior_d0 + d0)
      let mod_var = if mod_var < 0.0000001 { 0.0000001 } else { mod_var }
      // Compare first vs second condition (if 2 conditions).
      if n_conditions >= 2 && group_ns[0] > 0 && group_ns[1] > 0 {
        let mean_diff = group_means[1] - group_means[0]
        let se = dc_sqrt(mod_var * (1.0 / group_ns[0].to_double() + 1.0 / group_ns[1].to_double()))
        let se = if se < 0.0000001 { 0.0000001 } else { se }
        let t_stat = mean_diff / se
        let df = prior_d0 + d0
        // Convert t to p-value using normal approximation for large df.
        let p_value = 2.0 * (1.0 - dc_normal_cdf(dc_abs(t_stat)))
        let p_value = if p_value < 0.0 { 0.0 } else if p_value > 1.0 { 1.0 } else { p_value }
        let mname = if marker < marker_names.length() {
          marker_names[marker]
        } else {
          "marker_" + marker.to_string()
        }
        results.push(DiffcytDSResult::{
          cluster_id: cluster,
          marker_idx: marker,
          marker_name: mname,
          log_fc: mean_diff,
          p_value,
          fdr: p_value,
          significant: false,
        })
      }
    }
  }
  // BH-FDR correction.
  dc_bh_fdr_ds(results)
}

// ===========================================================================
// Multiple testing correction (BH-FDR)
// ===========================================================================

///|
/// Benjamini-Hochberg FDR correction for DA results.
fn dc_bh_fdr_da(results : Array[DiffcytDAResult]) -> Array[DiffcytDAResult] {
  let n = results.length()
  if n == 0 {
    return results
  }
  // Sort indices by p-value ascending.
  let indices = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  // Bubble sort by p-value.
  for i in 0..<n {
    let mut j = i + 1
    while j < n {
      if results[indices[j]].p_value < results[indices[i]].p_value {
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
      j = j + 1
    }
  }
  // Apply BH correction.
  for rank in 0..<n {
    let idx = indices[rank]
    let rank_f = (rank + 1).to_double()
    let corrected = results[idx].p_value * n.to_double() / rank_f
    let corrected = if corrected > 1.0 { 1.0 } else { corrected }
    results[idx].fdr = corrected
  }
  // Ensure monotonicity (from largest to smallest rank).
  let mut prev_fdr = 1.0
  let mut r = n - 1
  while r >= 0 {
    let idx = indices[r]
    if results[idx].fdr > prev_fdr {
      results[idx].fdr = prev_fdr
    }
    prev_fdr = results[idx].fdr
    r = r - 1
  }
  results
}

///|
/// Benjamini-Hochberg FDR correction for DS results.
fn dc_bh_fdr_ds(results : Array[DiffcytDSResult]) -> Array[DiffcytDSResult] {
  let n = results.length()
  if n == 0 {
    return results
  }
  let indices = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  for i in 0..<n {
    let mut j = i + 1
    while j < n {
      if results[indices[j]].p_value < results[indices[i]].p_value {
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
      j = j + 1
    }
  }
  for rank in 0..<n {
    let idx = indices[rank]
    let rank_f = (rank + 1).to_double()
    let corrected = results[idx].p_value * n.to_double() / rank_f
    let corrected = if corrected > 1.0 { 1.0 } else { corrected }
    results[idx].fdr = corrected
  }
  let mut prev_fdr = 1.0
  let mut r = n - 1
  while r >= 0 {
    let idx = indices[r]
    if results[idx].fdr > prev_fdr {
      results[idx].fdr = prev_fdr
    }
    prev_fdr = results[idx].fdr
    r = r - 1
  }
  results
}

///|
/// Mark significant results based on FDR threshold.
pub fn diffcyt_mark_significant(
  da_results : Array[DiffcytDAResult],
  ds_results : Array[DiffcytDSResult],
  fdr_threshold : Double
) -> (Array[DiffcytDAResult], Array[DiffcytDSResult]) {
  for r in da_results {
    r.significant = r.fdr < fdr_threshold
  }
  for r in ds_results {
    r.significant = r.fdr < fdr_threshold
  }
  (da_results, ds_results)
}

// ===========================================================================
// Linear algebra helpers
// ===========================================================================

///|
/// Solve linear system A*x = b via Gaussian elimination with partial pivoting.
fn dc_solve(matrix : Array[Array[Double]], b : Array[Double]) -> Array[Double]? {
  let n = matrix.length()
  if n == 0 {
    return Some([])
  }
  let m : Array[Array[Double]] = Array::new()
  for i in 0..<n {
    let row = Array::make(n + 1, 0.0)
    for j in 0..<n {
      row[j] = matrix[i][j]
    }
    row[n] = b[i]
    m.push(row)
  }
  for col in 0..<n {
    let mut pivot = col
    let mut max_val = dc_abs(m[col][col])
    let mut r = col + 1
    while r < n {
      let val = dc_abs(m[r][col])
      if val > max_val {
        max_val = val
        pivot = r
      }
      r = r + 1
    }
    if max_val < 0.0000000001 {
      return None
    }
    if pivot != col {
      let tmp = m[col]
      m[col] = m[pivot]
      m[pivot] = tmp
    }
    let mut r2 = col + 1
    while r2 < n {
      let factor = m[r2][col] / m[col][col]
      for c in col..<(n + 1) {
        m[r2][c] = m[r2][c] - factor * m[col][c]
      }
      r2 = r2 + 1
    }
  }
  let x = Array::make(n, 0.0)
  let mut i = n - 1
  while i >= 0 {
    let mut s = m[i][n]
    let mut j = i + 1
    while j < n {
      s = s - m[i][j] * x[j]
      j = j + 1
    }
    x[i] = s / m[i][i]
    i = i - 1
  }
  Some(x)
}

///|
/// Invert an n×n matrix via Gauss-Jordan elimination.
fn dc_matrix_inverse(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let n = matrix.length()
  if n == 0 {
    return Some([])
  }
  // Augmented [matrix | I]
  let aug : Array[Array[Double]] = Array::new()
  for i in 0..<n {
    let row = Array::make(2 * n, 0.0)
    for j in 0..<n {
      row[j] = matrix[i][j]
    }
    row[n + i] = 1.0
    aug.push(row)
  }
  for col in 0..<n {
    let mut pivot = col
    let mut max_val = dc_abs(aug[col][col])
    let mut r = col + 1
    while r < n {
      let val = dc_abs(aug[r][col])
      if val > max_val {
        max_val = val
        pivot = r
      }
      r = r + 1
    }
    if max_val < 0.0000000001 {
      return None
    }
    if pivot != col {
      let tmp = aug[col]
      aug[col] = aug[pivot]
      aug[pivot] = tmp
    }
    let piv = aug[col][col]
    for j in 0..<(2 * n) {
      aug[col][j] = aug[col][j] / piv
    }
    for r2 in 0..<n {
      if r2 != col {
        let factor = aug[r2][col]
        for j in 0..<(2 * n) {
          aug[r2][j] = aug[r2][j] - factor * aug[col][j]
        }
      }
    }
  }
  let inv : Array[Array[Double]] = Array::new()
  for i in 0..<n {
    let row = Array::make(n, 0.0)
    for j in 0..<n {
      row[j] = aug[i][n + j]
    }
    inv.push(row)
  }
  Some(inv)
}

///|
/// Mean of an array.
fn dc_mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    s = s + v
  }
  s / n.to_double()
}

///|
/// Variance of an array (sample variance with Bessel's correction).
fn dc_variance(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let m = dc_mean(arr)
  let mut s = 0.0
  for v in arr {
    let d = v - m
    s = s + d * d
  }
  s / (n - 1).to_double()
}

// ===========================================================================
// Top results extraction
// ===========================================================================

///|
/// Extract top DA clusters sorted by p-value.
pub fn diffcyt_top_clusters_da(
  da_results : Array[DiffcytDAResult],
  n_top : Int
) -> Array[DiffcytDAResult] {
  // Copy and sort by p-value.
  let sorted : Array[DiffcytDAResult] = Array::new()
  for r in da_results {
    sorted.push(r)
  }
  let n = sorted.length()
  for i in 0..<n {
    let mut j = i + 1
    while j < n {
      if sorted[j].p_value < sorted[i].p_value {
        let tmp = sorted[i]
        sorted[i] = sorted[j]
        sorted[j] = tmp
      }
      j = j + 1
    }
  }
  let limit = if n_top < n { n_top } else { n }
  let result : Array[DiffcytDAResult] = Array::new()
  for i in 0..<limit {
    result.push(sorted[i])
  }
  result
}

///|
/// Extract top DS cluster-marker pairs sorted by p-value.
pub fn diffcyt_top_table_ds(
  ds_results : Array[DiffcytDSResult],
  n_top : Int
) -> Array[DiffcytDSResult] {
  let sorted : Array[DiffcytDSResult] = Array::new()
  for r in ds_results {
    sorted.push(r)
  }
  let n = sorted.length()
  for i in 0..<n {
    let mut j = i + 1
    while j < n {
      if sorted[j].p_value < sorted[i].p_value {
        let tmp = sorted[i]
        sorted[i] = sorted[j]
        sorted[j] = tmp
      }
      j = j + 1
    }
  }
  let limit = if n_top < n { n_top } else { n }
  let result : Array[DiffcytDSResult] = Array::new()
  for i in 0..<limit {
    result.push(sorted[i])
  }
  result
}

// ===========================================================================
// Full pipeline
// ===========================================================================

///|
/// Run the complete diffcyt pipeline: cluster, compute counts/medians,
/// test DA and DS, and mark significant results.
///
/// `cells` is the array of cytometry cells. `n_clusters` is the number of
/// clusters for FlowSOM. `fdr_threshold` is the significance threshold.
/// Returns a DiffcytSummary with all results.
pub fn diffcyt_run(
  cells : Array[CytometryCell],
  n_clusters : Int,
  fdr_threshold : Double
) -> DiffcytSummary {
  let sample_ids = diffcyt_unique_samples(cells)
  let conditions_arr : Array[String] = Array::new()
  let sample_condition_map : Map[String, String] = Map([], capacity=sample_ids.length())
  for cell in cells {
    sample_condition_map.set(cell.sample_id, cell.condition)
  }
  for sid in sample_ids {
    conditions_arr.push(sample_condition_map.get(sid).unwrap_or(""))
  }
  let unique_conds = diffcyt_unique_conditions(cells)
  let n_markers = if cells.length() > 0 { cells[0].n_markers() } else { 0 }
  // 1. Cluster cells.
  let _ = diffcyt_cluster_cells(cells, n_clusters, 10)
  // 2. Compute counts.
  let counts = diffcyt_calc_counts(cells, n_clusters, sample_ids)
  // 3. Build design matrix.
  let design = diffcyt_create_design_matrix(sample_ids, conditions_arr)
  // 4. Determine contrast column (last column = stimulated condition).
  let contrast_col = design[0].length() - 1
  // 5. DA testing.
  let da_results = diffcyt_testDA(counts, design, contrast_col, n_clusters)
  // 6. DS testing.
  let marker_names = Array::new()
  for i in 0..<n_markers {
    marker_names.push("marker_" + i.to_string())
  }
  let ds_results = diffcyt_testDS(
    cells, sample_ids, conditions_arr, n_clusters, n_markers, marker_names
  )
  // 7. Mark significant.
  let (da_results, ds_results) = diffcyt_mark_significant(
    da_results, ds_results, fdr_threshold
  )
  DiffcytSummary::{
    n_cells: cells.length(),
    n_clusters,
    n_markers,
    n_samples: sample_ids.length(),
    n_conditions: unique_conds.length(),
    da_results,
    ds_results,
    fdr_threshold,
  }
}

// ===========================================================================
// Sample data
// ===========================================================================

///|
/// Create sample cytometry data for testing.
/// 4 samples (2 ctrl, 2 stim), 3 markers, 20 cells per sample.
/// Stimulated samples have higher expression of marker_2 in a subset of cells.
pub fn diffcyt_sample_data() -> Array[CytometryCell] {
  let cells : Array[CytometryCell] = Array::new()
  let samples = ["S1", "S2", "S3", "S4"]
  let conditions = ["ctrl", "ctrl", "stim", "stim"]
  let n_markers = 3
  let mut rng_state = 12345
  let mut cell_counter = 0
  for si in 0..<samples.length() {
    for _ in 0..<20 {
      let marker_vals = Array::make(n_markers, 0.0)
      for m in 0..<n_markers {
        let (s, r) = dc_lcg_rand(rng_state)
        rng_state = s
        let base = r * 5.0
        let val = if conditions[si] == "stim" && m == 1 {
          base + 3.0 // stimulated condition upregulates marker 2
        } else {
          base
        }
        marker_vals[m] = val
      }
      let cell = CytometryCell::new(
        "cell_" + cell_counter.to_string(),
        samples[si],
        conditions[si],
        marker_vals
      )
      cells.push(cell)
      cell_counter = cell_counter + 1
    }
  }
  cells
}
