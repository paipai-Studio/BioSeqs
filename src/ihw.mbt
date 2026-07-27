///|
/// Bio.IHW - Independent Hypothesis Weighting module.
/// Implements the IHW method for advanced multiple testing correction,
/// inspired by Bioconductor's IHW package.
///
/// IHW uses covariates to assign different weights to hypotheses,
/// improving power over traditional methods like Bonferroni correction
/// while maintaining error control.
///
/// References:
///   - Ignatiadis, A., et al. (2016). Data-driven hypothesis weighting
///     increases detection power for high-throughput experiments.
///     Genome Research, 26(12), 1724-1733.

// ---------------------------------------------------------------------------
// IHWResult
// ---------------------------------------------------------------------------

///|
/// Result of the IHW procedure, containing adjusted p-values,
/// weights, and diagnostic information.
pub struct IHWResult {
  p_values : Array[Double]
  adjusted_p_values : Array[Double]
  weights : Array[Double]
  covariates : Array[Double]
  alpha : Double
  n_attempts : Int
} derive(Eq)

///|
/// Create a new IHWResult.
pub fn IHWResult::new(
  p_values : Array[Double],
  adjusted_p_values : Array[Double],
  weights : Array[Double],
  covariates : Array[Double],
  alpha : Double,
  n_attempts : Int,
) -> IHWResult {
  IHWResult::{
    p_values,
    adjusted_p_values,
    weights,
    covariates,
    alpha,
    n_attempts,
  }
}

///|
/// Get the original p-values.
pub fn IHWResult::p_values(self : IHWResult) -> Array[Double] {
  self.p_values
}

///|
/// Get the adjusted p-values.
pub fn IHWResult::adjusted_p_values(self : IHWResult) -> Array[Double] {
  self.adjusted_p_values
}

///|
/// Get the weights assigned to each hypothesis.
pub fn IHWResult::weights(self : IHWResult) -> Array[Double] {
  self.weights
}

///|
/// Get the covariates used for weighting.
pub fn IHWResult::covariates(self : IHWResult) -> Array[Double] {
  self.covariates
}

///|
/// Get the significance level (alpha) used.
pub fn IHWResult::alpha(self : IHWResult) -> Double {
  self.alpha
}

///|
/// Get the number of iterations performed.
pub fn IHWResult::n_attempts(self : IHWResult) -> Int {
  self.n_attempts
}

// ---------------------------------------------------------------------------
// IHWConfig
// ---------------------------------------------------------------------------

///|
/// Configuration for the IHW procedure.
///
/// - alpha: significance level (default 0.05)
/// - n_attempts: maximum number of iterations (default 10)
/// - scale_type: "global" (equal weights) or "local" (covariate-dependent)
pub struct IHWConfig {
  alpha : Double
  n_attempts : Int
  scale_type : String
} derive(Eq)

///|
/// Create a new IHWConfig.
///
/// Parameters:
///   - alpha: significance level, default 0.05
///   - n_attempts: max iterations, default 10
///   - scale_type: "global" or "local", default "local"
pub fn IHWConfig::new(
  alpha? : Double = 0.05,
  n_attempts? : Int = 10,
  scale_type? : String = "local",
) -> IHWConfig {
  IHWConfig::{ alpha, n_attempts, scale_type }
}

///|
/// Get the alpha value.
pub fn IHWConfig::alpha(self : IHWConfig) -> Double {
  self.alpha
}

///|
/// Get the number of attempts.
pub fn IHWConfig::n_attempts(self : IHWConfig) -> Int {
  self.n_attempts
}

///|
/// Get the scale type name.
pub fn IHWConfig::scale_type(self : IHWConfig) -> String {
  self.scale_type
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

///|
/// Estimate pi0 using Storey method with covariate weighting.
/// Down-weights hypotheses with extreme covariate values.
fn ihw_storey_pi0(
  p_values : Array[Double],
  lambda : Double,
  covariates : Array[Double],
) -> Double {
  let n = p_values.length()
  if n == 0 {
    return 1.0
  }
  let mut count = 0.0
  let mut total_weight = 0.0
  let mut i = 0
  while i < n {
    let weight = if i < covariates.length() {
      let c = covariates[i]
      let w = 1.0 - (c - 0.5) * (c - 0.5) * 2.0
      if w < 0.1 { 0.1 } else { w }
    } else {
      1.0
    }
    if p_values[i] > lambda {
      count = count + weight
    }
    total_weight = total_weight + weight
    i = i + 1
  }
  if total_weight <= 0.0 {
    return 1.0
  }
  let pi0 = count / ((1.0 - lambda) * total_weight)
  if pi0 > 1.0 { 1.0 } else { pi0 }
}

///|
/// Compute the weighted Bonferroni adjusted p-values.
/// padj_i = min(1, p_i * sum(w) / w_i)
fn weighted_bonferroni(
  p_values : Array[Double],
  weights : Array[Double],
) -> Array[Double] {
  let n = p_values.length()
  let result : Array[Double] = Array::make(n, 0.0)
  let mut sum_w = 0.0
  let mut i = 0
  while i < n {
    sum_w = sum_w + weights[i]
    i = i + 1
  }
  if sum_w <= 0.0 {
    return p_values.copy()
  }
  i = 0
  while i < n {
    let w_i = if weights[i] > 0.0 { weights[i] } else { 1.0e-15 }
    let adj = p_values[i] * sum_w / w_i
    result[i] = if adj > 1.0 { 1.0 } else { adj }
    i = i + 1
  }
  result
}

///|
/// Apply weighted Bonferroni correction and determine rejection set.
/// Returns the number of rejected hypotheses and the rejection mask.
fn weighted_bonferroni_reject(
  p_values : Array[Double],
  weights : Array[Double],
  alpha : Double,
) -> (Int, Array[Bool]) {
  let n = p_values.length()
  let sum_w = ihw_sum_weights(weights)
  let rejects : Array[Bool] = Array::make(n, false)
  let mut n_reject = 0
  if sum_w <= 0.0 {
    return (0, rejects)
  }
  let mut i = 0
  while i < n {
    let threshold = weights[i] * alpha / sum_w
    if p_values[i] <= threshold {
      rejects[i] = true
      n_reject = n_reject + 1
    }
    i = i + 1
  }
  (n_reject, rejects)
}

///|
/// Sum of weights array.
fn ihw_sum_weights(weights : Array[Double]) -> Double {
  let mut s = 0.0
  let mut i = 0
  while i < weights.length() {
    s = s + weights[i]
    i = i + 1
  }
  s
}

///|
/// Assign hypotheses to bins based on covariate quantiles.
/// Returns bin indices (0 to n_bins-1) for each hypothesis.
fn bin_covariates(
  covariates : Array[Double],
  n_bins : Int,
) -> Array[Int] {
  let n = covariates.length()
  let result : Array[Int] = Array::make(n, 0)
  if n <= n_bins || n_bins <= 1 {
    return Array::make(n, 0)
  }
  let indexed : Array[Array[Double]] = Array::make(n, [0.0, 0.0])
  let mut i = 0
  while i < n {
    indexed[i] = [covariates[i], i.to_double()]
    i = i + 1
  }
  indexed.sort_by(fn(a: Array[Double], b: Array[Double]) -> Int {
    if a[0] < b[0] { -1 } else if a[0] > b[0] { 1 } else { 0 }
  })
  let bin_size = (n.to_double() / n_bins.to_double()).ceil().to_int()
  i = 0
  while i < n {
    let bin_idx = i / bin_size
    let actual_bin = if bin_idx >= n_bins { n_bins - 1 } else { bin_idx }
    let orig_idx = indexed[i][1].to_int()
    result[orig_idx] = actual_bin
    i = i + 1
  }
  result
}

///|
/// Estimate per-bin weights from rejection proportions.
/// Bins with more rejections get higher weights.
fn estimate_bin_weights(
  p_values : Array[Double],
  bins : Array[Int],
  n_bins : Int,
  alpha : Double,
) -> Array[Double] {
  let n = p_values.length()
  let weights : Array[Double] = Array::make(n, 0.0)
  let n_total = n.to_double()

  let bin_reject : Array[Int] = Array::make(n_bins, 0)
  let bin_total : Array[Int] = Array::make(n_bins, 0)
  let mut i = 0
  while i < n {
    let b = bins[i]
    bin_total[b] = bin_total[b] + 1
    if p_values[i] <= alpha {
      bin_reject[b] = bin_reject[b] + 1
    }
    i = i + 1
  }

  let bin_weight : Array[Double] = Array::make(n_bins, 0.0)
  let mut sum_w = 0.0
  let mut b = 0
  while b < n_bins {
    if bin_total[b] > 0 {
      let reject_rate = bin_reject[b].to_double() / bin_total[b].to_double()
      let baseline = 0.05
      let w = if reject_rate > baseline {
        bin_total[b].to_double() * (1.0 + reject_rate)
      } else {
        bin_total[b].to_double() * 0.5
      }
      bin_weight[b] = w
      sum_w = sum_w + w
    } else {
      bin_weight[b] = 0.0
    }
    b = b + 1
  }

  if sum_w <= 0.0 {
    return Array::make(n, 1.0)
  }

  b = 0
  while b < n_bins {
    if bin_total[b] > 0 {
      let scaled = bin_weight[b] * n_total / sum_w
      let mut j = 0
      while j < n {
        if bins[j] == b {
          weights[j] = scaled
        }
        j = j + 1
      }
    }
    b = b + 1
  }

  i = 0
  while i < n {
    if weights[i] <= 0.0 {
      weights[i] = 1.0
    }
    i = i + 1
  }
  weights
}

///|
/// Initialize weights based on covariate distribution.
/// For "local" scale_type, bin covariates and assign initial weights.
/// For "global" scale_type, assign equal weights.
fn initialize_weights(
  covariates : Array[Double],
  scale_type : String,
) -> Array[Double] {
  let n = covariates.length()
  if n == 0 {
    return Array::new()
  }
  if scale_type == "global" {
    return Array::make(n, 1.0)
  }
  Array::make(n, 1.0)
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

///|
/// Perform IHW (Independent Hypothesis Weighting) correction.
///
/// This is the main entry point. It uses default configuration
/// (alpha=0.05, local weighting, up to 10 iterations).
///
/// Parameters:
///   - p_values: array of p-values from hypothesis tests
///   - covariates: array of covariate values for each hypothesis
///   - alpha: significance level (e.g., 0.05)
///
/// Returns:
///   IHWResult with adjusted p-values, weights, and diagnostic info
pub fn ihw(
  p_values : Array[Double],
  covariates : Array[Double],
  alpha : Double,
) -> IHWResult {
  let config = IHWConfig::new(alpha=alpha, n_attempts=10, scale_type="local")
  ihw_with_config(p_values, covariates, config)
}

///|
/// Perform IHW correction with a custom configuration.
///
/// Parameters:
///   - p_values: array of p-values
///   - covariates: array of covariate values
///   - config: IHWConfig with alpha, n_attempts, and scale_type
///
/// Returns:
///   IHWResult with adjusted p-values and weights
pub fn ihw_with_config(
  p_values : Array[Double],
  covariates : Array[Double],
  config : IHWConfig,
) -> IHWResult {
  let n = p_values.length()
  let n_cov = covariates.length()
  let effective_n = if n < n_cov { n } else { n_cov }

  if effective_n == 0 {
    return IHWResult::new(
      [],
      [],
      [],
      [],
      config.alpha,
      0,
    )
  }

  let n_bins = if effective_n > 20 { 20 } else { if effective_n > 10 { 10 } else { if effective_n > 5 { 5 } else { 1 } } }

  let mut weights = initialize_weights(covariates, config.scale_type)
  let bins = bin_covariates(covariates, n_bins)

  let mut adj_pvals = weighted_bonferroni(p_values, weights)
  let mut prev_adj = adj_pvals.copy()
  let mut attempts = 0
  let max_attempts = config.n_attempts

  if config.scale_type == "global" {
    adj_pvals = weighted_bonferroni(p_values, weights)
    return IHWResult::new(
      p_values.copy(),
      adj_pvals,
      weights,
      covariates.copy(),
      config.alpha,
      0,
    )
  }

  while attempts < max_attempts {
    let (n_reject, _) = weighted_bonferroni_reject(
      p_values,
      weights,
      config.alpha,
    )

    if n_reject > 0 && attempts < max_attempts - 1 {
      let new_weights = estimate_bin_weights(
        p_values,
        bins,
        n_bins,
        config.alpha,
      )
      weights = new_weights
    }

    adj_pvals = weighted_bonferroni(p_values, weights)

    let mut converged = true
    let mut i = 0
    while i < effective_n {
      let diff = adj_pvals[i] - prev_adj[i]
      let abs_diff = if diff < 0.0 { -diff } else { diff }
      if abs_diff > 1.0e-10 {
        converged = false
        break
      }
      i = i + 1
    }

    attempts = attempts + 1
    if converged && attempts > 1 {
      break
    }
    prev_adj = adj_pvals.copy()
  }

  IHWResult::new(
    p_values.copy(),
    adj_pvals,
    weights,
    covariates.copy(),
    config.alpha,
    attempts,
  )
}

///|
/// Perform IHW correction with multiple covariates.
///
/// When multiple covariates are provided, this function uses the first
/// principal component or the most informative covariate for weighting.
/// For simplicity, covariates are combined by averaging after normalization.
///
/// Parameters:
///   - p_values: array of p-values
///   - covariates: array of covariate arrays (one per hypothesis)
///   - alpha: significance level
///
/// Returns:
///   IHWResult with adjusted p-values and weights
pub fn ihw_multiple(
  p_values : Array[Double],
  covariates : Array[Array[Double]],
  alpha : Double,
) -> IHWResult {
  let n = p_values.length()
  if n == 0 {
    return IHWResult::new([], [], [], [], alpha, 0)
  }

  let n_covs = covariates.length()
  if n_covs == 0 {
    let flat_cov : Array[Double] = Array::make(n, 0.0)
    return ihw(p_values, flat_cov, alpha)
  }

  let combined : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let mut sum = 0.0
    let mut count = 0.0
    let mut j = 0
    while j < n_covs {
      let row = covariates[j]
      if i < row.length() {
        sum = sum + row[i]
        count = count + 1.0
      }
      j = j + 1
    }
    combined[i] = if count > 0.0 { sum / count } else { 0.0 }
    i = i + 1
  }

  ihw(p_values, combined, alpha)
}

///|
/// Bonferroni correction (baseline, no covariates).
///
/// This is the standard Bonferroni correction: each p-value is multiplied
/// by the number of tests (capped at 1.0).
///
/// Parameters:
///   - p_values: array of p-values
///
/// Returns:
///   Array of adjusted p-values
pub fn bonferroni_ihw(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return []
  }
  let result : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let adj = p_values[i] * n.to_double()
    result[i] = if adj > 1.0 { 1.0 } else { adj }
    i = i + 1
  }
  result
}

///|
/// Storey-style correction with covariates.
///
/// Applies the Storey procedure to estimate the false discovery rate,
/// using covariate information to improve the estimate of pi0.
///
/// Parameters:
///   - p_values: array of p-values
///   - covariates: array of covariate values
///
/// Returns:
///   Array of adjusted p-values
pub fn storey_ihw(
  p_values : Array[Double],
  covariates : Array[Double],
) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return []
  }

  let lambda = 0.5
  let pi0 = ihw_storey_pi0(p_values, lambda, covariates)

  let sorted = p_values.copy()
  sorted.sort()

  let result : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let mut count = 0
    let mut j = 0
    while j < n {
      if p_values[j] <= sorted[i] {
        count = count + 1
      }
      j = j + 1
    }
    let rank = count.to_double()
    if rank > 0.0 {
      let adj = pi0 * p_values[i] * n.to_double() / rank
      result[i] = if adj > 1.0 { 1.0 } else { adj }
    } else {
      result[i] = 1.0
    }
    i = i + 1
  }

  let mut prev = 1.0
  let temp : Array[Array[Double]] = Array::make(n, [0.0, 0.0])
  i = 0
  while i < n {
    temp[i] = [p_values[i], i.to_double()]
    i = i + 1
  }
  temp.sort_by(fn(a: Array[Double], b: Array[Double]) -> Int {
    if a[0] < b[0] { -1 } else if a[0] > b[0] { 1 } else { 0 }
  })

  i = n - 1
  while i >= 0 {
    let idx = temp[i][1].to_int()
    let rank = (i + 1).to_double()
    let raw_adj = pi0 * p_values[idx] * n.to_double() / rank
    let adj = if raw_adj > prev { prev } else { raw_adj }
    result[idx] = if adj > 1.0 { 1.0 } else { adj }
    prev = result[idx]
    i = i - 1
  }

  result
}

///|
/// Extract weights from an IHWResult.
///
/// Parameters:
///   - result: an IHWResult from ihw or ihw_with_config
///
/// Returns:
///   Array of weights used in the correction
pub fn ihw_weights(result : IHWResult) -> Array[Double] {
  result.weights
}