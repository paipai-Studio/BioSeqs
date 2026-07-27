///|
/// Bio.QValue module - Storey's q-value FDR correction method.
/// Implements the q-value calculation as described in Storey (2002).
/// Based on Bioconductor's qvalue package.
///
/// The q-value provides a measure of significance in terms of the
/// false discovery rate (FDR). For a given p-value threshold, the
/// q-value is the minimum FDR at which that test would be declared
/// significant.

// ---------------------------------------------------------------------------
// QValueResult - Result container for q-value computations
// ---------------------------------------------------------------------------

///|
/// Container for q-value computation results.
///
/// Fields:
///   - p_values:   The input p-values
///   - q_values:   The computed q-values (same order as p_values)
///   - pi0:        Estimated proportion of true null hypotheses
///   - lambda:     The tuning parameter λ used for pi0 estimation
///   - gamma:      Bootstrap estimate of pi0 variability
///   - significance: Boolean array indicating significance at threshold
pub struct QValueResult {
  p_values : Array[Double]
  q_values : Array[Double]
  pi0 : Double
  lambda : Double
  gamma : Double
  significance : Array[Bool]
}

///|
/// Create a new QValueResult.
pub fn QValueResult::new(
  p_values : Array[Double],
  q_values : Array[Double],
  pi0 : Double,
  lambda : Double,
  gamma : Double,
  significance : Array[Bool],
) -> QValueResult {
  QValueResult::{
    p_values,
    q_values,
    pi0,
    lambda,
    gamma,
    significance,
  }
}

///|
/// Get the input p-values.
pub fn QValueResult::p_values(self : QValueResult) -> Array[Double] {
  self.p_values
}

///|
/// Get the computed q-values.
pub fn QValueResult::q_values(self : QValueResult) -> Array[Double] {
  self.q_values
}

///|
/// Get the estimated pi0 (proportion of true null hypotheses).
pub fn QValueResult::pi0(self : QValueResult) -> Double {
  self.pi0
}

///|
/// Get the tuning parameter lambda used.
pub fn QValueResult::lambda(self : QValueResult) -> Double {
  self.lambda
}

///|
/// Get the bootstrap gamma value.
pub fn QValueResult::gamma(self : QValueResult) -> Double {
  self.gamma
}

///|
/// Get the significance boolean array.
pub fn QValueResult::significance(self : QValueResult) -> Array[Bool] {
  self.significance
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

///|
/// Validate and clamp p-values to [0, 1].
fn _validate_pvalues(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  let result : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let p = p_values[i]
    if p < 0.0 {
      result[i] = 0.0
    } else if p > 1.0 {
      result[i] = 1.0
    } else {
      result[i] = p
    }
    i = i + 1
  }
  result
}

///|
/// Simple deterministic pseudo-random number generator (LCG).
/// Returns the new state and a random Double in [0, 1).
fn _lcg_next(state : Int) -> (Int, Double) {
  let new_state = (state * 1103515245 + 12345) & 0x7fffffff
  let rand_val = new_state.to_double() / 2147483647.0
  (new_state, rand_val)
}

///|
/// Generate a bootstrap resample of the given array.
/// Uses deterministic LCG seeded from the data.
fn _bootstrap_sample(p_values : Array[Double], seed : Int) -> Array[Double] {
  let n = p_values.length()
  let result : Array[Double] = Array::make(n, 0.0)
  let mut state = seed
  let mut i = 0
  while i < n {
    let (new_state, rand_val) = _lcg_next(state)
    state = new_state
    let idx = (rand_val * n.to_double()).to_int()
    result[i] = p_values[idx]
    i = i + 1
  }
  result
}

///|
/// Compute pi0 estimate for a given lambda.
/// π_0(λ) = #{p_i > λ} / (n * (1 - λ))
fn _compute_pi0_for_lambda(p_values : Array[Double], lambda : Double) -> Double {
  let n = p_values.length()
  if n == 0 {
    return 0.0
  }
  let mut count_above = 0
  let mut i = 0
  while i < n {
    if p_values[i] > lambda {
      count_above = count_above + 1
    }
    i = i + 1
  }
  let denom = n.to_double() * (1.0 - lambda)
  if denom <= 0.0 {
    return 1.0
  }
  let pi0 = count_above.to_double() / denom
  if pi0 > 1.0 {
    1.0
  } else if pi0 < 0.0 {
    0.0
  } else {
    pi0
  }
}

///|
/// Choose optimal lambda via bootstrap MSE minimization.
/// Grid: 0.05, 0.10, ..., 0.95
fn _choose_lambda(p_values : Array[Double]) -> Double {
  let lambdas = [
    0.05, 0.10, 0.15, 0.20, 0.25,
    0.30, 0.35, 0.40, 0.45, 0.50,
    0.55, 0.60, 0.65, 0.70, 0.75,
    0.80, 0.85, 0.90, 0.95,
  ]
  let n = p_values.length()
  if n == 0 {
    return 0.5
  }

  let pi0_values : Array[Double] = Array::make(lambdas.length(), 0.0)
  let mut i = 0
  while i < lambdas.length() {
    pi0_values[i] = _compute_pi0_for_lambda(p_values, lambdas[i])
    i = i + 1
  }

  let n_bootstrap = Int::min(100, n * 2)
  let mut best_idx = -1
  let mut best_mse = 0.0

  let mut li = 0
  while li < lambdas.length() {
    let lambda = lambdas[li]
    let pi0_orig = pi0_values[li]
    let mut mse_sum : Double = 0.0
    let seed_base = (n * 13 + li * 7 + 42) & 0x7fffffff

    let mut b = 0
    while b < n_bootstrap {
      let seed = (seed_base * (b + 1) + 11) & 0x7fffffff
      let bootstrap_sample = _bootstrap_sample(p_values, seed)
      let pi0_boot = _compute_pi0_for_lambda(bootstrap_sample, lambda)
      let diff = pi0_boot - pi0_orig
      mse_sum = mse_sum + diff * diff
      b = b + 1
    }

    let mse : Double = mse_sum / n_bootstrap.to_double()
    if best_idx < 0 || mse < best_mse {
      best_mse = mse
      best_idx = li
    }
    li = li + 1
  }

  if best_idx < 0 {
    lambdas[9]
  } else {
    lambdas[best_idx]
  }
}

///|
/// Compute q-values using Storey's method.
///
/// Algorithm:
///   1. Sort p-values ascending, tracking original indices
///   2. For k = n-1 down to 0:
///      q_(k) = min(q_(k+1), p_(k) * π_0 * n / (k+1))
///      where q_(n) = p_(n) * π_0 (the largest p-value)
///   3. Map q-values back to original order
fn _compute_qvalues(p_values : Array[Double], pi0 : Double) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return []
  }
  if n == 1 {
    let q = p_values[0] * pi0
    if q > 1.0 { [1.0] } else { [q] }
  } else {
    let indexed : Array[Array[Double]] = Array::new()
    let mut i = 0
    while i < n {
      indexed.push([p_values[i], i.to_double()])
      i = i + 1
    }

    indexed.sort_by(fn(a : Array[Double], b : Array[Double]) -> Int {
      if a[0] < b[0] { -1 } else if a[0] > b[0] { 1 } else { 0 }
    })

    let q_sorted : Array[Double] = Array::make(n, 1.0)
    let mut q_prev = 1.0
    let mut k = n - 1
    while k >= 0 {
      let p = indexed[k][0]
      let rank = (k + 1).to_double()
      let q_val = p * pi0 * n.to_double() / rank
      if q_val < q_prev {
        q_prev = q_val
      }
      q_sorted[k] = if q_prev > 1.0 { 1.0 } else { q_prev }
      k = k - 1
    }

    let q_values : Array[Double] = Array::make(n, 0.0)
    let mut j = 0
    while j < n {
      let orig_idx = indexed[j][1].to_int()
      q_values[orig_idx] = q_sorted[j]
      j = j + 1
    }

    q_values
  }
}

///|
/// Compute significance boolean array from q-values at a given threshold.
fn _compute_significance(q_values : Array[Double], threshold : Double) -> Array[Bool] {
  let n = q_values.length()
  let result : Array[Bool] = Array::make(n, false)
  let mut i = 0
  while i < n {
    result[i] = q_values[i] <= threshold
    i = i + 1
  }
  result
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

///|
/// Estimate pi0 (proportion of true null hypotheses) for a given lambda.
///
/// π_0(λ) = #{p_i > λ} / (n * (1 - λ))
///
/// Parameters:
///   - p_values: Array of p-values in [0, 1]
///   - lambda:   The tuning parameter λ ∈ [0, 1)
///
/// Returns:
///   - Estimated π_0, clamped to [0, 1]
pub fn qvalue_estimate_pi0(p_values : Array[Double], lambda : Double) -> Double {
  let validated = _validate_pvalues(p_values)
  _compute_pi0_for_lambda(validated, lambda)
}

///|
/// Compute q-values from p-values using Storey's method.
///
/// This is the main function that automatically selects the optimal
/// lambda via bootstrap-based MSE minimization, estimates pi0, and
/// computes q-values.
///
/// Parameters:
///   - p_values: Array of p-values in [0, 1]
///
/// Returns:
///   - QValueResult with q-values, pi0, lambda, and significance
pub fn qvalue(p_values : Array[Double]) -> QValueResult {
  let validated = _validate_pvalues(p_values)
  let n = validated.length()
  if n == 0 {
    return QValueResult::{
      p_values: [],
      q_values: [],
      pi0: 0.0,
      lambda: 0.5,
      gamma: 0.0,
      significance: [],
    }
  }

  let lambda = _choose_lambda(validated)
  let pi0 = _compute_pi0_for_lambda(validated, lambda)
  let q_values = _compute_qvalues(validated, pi0)
  let significance = _compute_significance(q_values, 0.05)

  let gamma = bootstrap_pi0(validated, 100)

  QValueResult::{
    p_values: validated,
    q_values,
    pi0,
    lambda,
    gamma,
    significance,
  }
}

///|
/// Compute q-values with a user-specified lambda.
///
/// Parameters:
///   - p_values: Array of p-values in [0, 1]
///   - lambda:   The tuning parameter λ ∈ [0, 1)
///
/// Returns:
///   - QValueResult with q-values, pi0, and significance
pub fn qvalue_with_lambda(
  p_values : Array[Double],
  lambda : Double,
) -> QValueResult {
  let validated = _validate_pvalues(p_values)
  let n = validated.length()
  if n == 0 {
    return QValueResult::{
      p_values: [],
      q_values: [],
      pi0: 0.0,
      lambda,
      gamma: 0.0,
      significance: [],
    }
  }

  let pi0 = _compute_pi0_for_lambda(validated, lambda)
  let q_values = _compute_qvalues(validated, pi0)
  let significance = _compute_significance(q_values, 0.05)

  QValueResult::{
    p_values: validated,
    q_values,
    pi0,
    lambda,
    gamma: pi0,
    significance,
  }
}

///|
/// Storey's FDR correction (Benjamini-Hochberg-like procedure).
///
/// Computes q-values using Storey's method with automatic lambda
/// selection. Returns only the array of q-values.
///
/// Parameters:
///   - p_values: Array of p-values in [0, 1]
///
/// Returns:
///   - Array of q-values in the same order as input p-values
pub fn storey_fdr(p_values : Array[Double]) -> Array[Double] {
  let result = qvalue(p_values)
  result.q_values
}

///|
/// Bootstrap estimation of pi0.
///
/// Resamples the p-values n_bootstrap times and computes the
/// median pi0 estimate across all bootstrap samples.
///
/// Parameters:
///   - p_values:     Array of p-values in [0, 1]
///   - n_bootstrap:  Number of bootstrap iterations
///
/// Returns:
///   - Bootstrap estimate of pi0 (median across samples)
pub fn bootstrap_pi0(p_values : Array[Double], n_bootstrap : Int) -> Double {
  let validated = _validate_pvalues(p_values)
  let n = validated.length()
  if n == 0 {
    return 0.0
  }
  if n_bootstrap <= 0 {
    return _compute_pi0_for_lambda(validated, 0.5)
  }

  let lambda = _choose_lambda(validated)
  let pi0_estimates : Array[Double] = Array::make(n_bootstrap, 0.0)
  let seed_base = (n * 31 + 7) & 0x7fffffff

  let mut b = 0
  while b < n_bootstrap {
    let seed = (seed_base * (b + 1) + 13) & 0x7fffffff
    let bootstrap_sample = _bootstrap_sample(validated, seed)
    pi0_estimates[b] = _compute_pi0_for_lambda(bootstrap_sample, lambda)
    b = b + 1
  }

  let sorted = pi0_estimates.copy()
  sorted.sort()
  let mid = n_bootstrap / 2
  if n_bootstrap % 2 == 0 {
    (sorted[mid - 1] + sorted[mid]) / 2.0
  } else {
    sorted[mid]
  }
}

///|
/// Get significance boolean array from a QValueResult at a given threshold.
///
/// Parameters:
///   - result:    A QValueResult from qvalue() or qvalue_with_lambda()
///   - threshold: The FDR threshold (e.g., 0.05)
///
/// Returns:
///   - Array of Bool where true indicates q-value ≤ threshold
pub fn qvalue_significance(
  result : QValueResult,
  threshold : Double,
) -> Array[Bool] {
  _compute_significance(result.q_values, threshold)
}