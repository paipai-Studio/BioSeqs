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
  QValueResult::{ p_values, q_values, pi0, lambda, gamma, significance, }
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
    0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
    0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
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
    if q > 1.0 {
      [1.0]
    } else {
      [q]
    }
  } else {
    let indexed : Array[Array[Double]] = Array(capacity=0)
    let mut i = 0
    while i < n {
      indexed.push([p_values[i], i.to_double()])
      i = i + 1
    }

    indexed.sort_by(fn(a : Array[Double], b : Array[Double]) -> Int {
      if a[0] < b[0] {
        -1
      } else if a[0] > b[0] {
        1
      } else {
        0
      }
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
fn _compute_significance(
  q_values : Array[Double],
  threshold : Double,
) -> Array[Bool] {
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

// ===========================================================================
// Bioconductor qvalue package: official pi0est / lfdr / empPvals / qvalue
//
// The functions below mirror the official qvalue package (Storey):
//   - pi0est:     pi0 over a lambda grid, with the smoother method
//                 (cubic smoothing spline, df = 3, Reinsch form) and the
//                 closed-form bootstrap method (David Robinson).
//   - lfdr:       local FDR via probit/logit transform, Gaussian KDE of the
//                 transformed p-values (Silverman/nrd0 bandwidth * adj).
//   - empPvals:   pooled / test-specific empirical p-values from observed and
//                 simulated null statistics.
//   - qvalue:     q-values with optional pfdr finite-sample correction and
//                 fdr.level significance calls.
// ===========================================================================

///|
/// Error raised by the q-value functions, mirroring Bioconductor qvalue's
/// input validation (p-values out of range, invalid lambda grid, etc.).
pub suberror QValueError {
  QValueError(String)
}

///|
/// Result of `qvalue_pi0est`: the estimated proportion of true nulls plus
/// the lambda path used for the estimate.
///
/// Fields:
///   - pi0:        Estimated proportion of true null p-values
///   - pi0_lambda: pi0 estimate at each lambda value
///   - lambdas:    The lambda grid (sorted ascending)
///   - pi0_smooth: Fitted smoother values at each lambda (empty for bootstrap)
pub struct Pi0EstResult {
  pi0 : Double
  pi0_lambda : Array[Double]
  lambdas : Array[Double]
  pi0_smooth : Array[Double]
}

///|
/// Create a new Pi0EstResult.
pub fn Pi0EstResult::new(
  pi0 : Double,
  pi0_lambda : Array[Double],
  lambdas : Array[Double],
  pi0_smooth : Array[Double],
) -> Pi0EstResult {
  Pi0EstResult::{ pi0, pi0_lambda, lambdas, pi0_smooth, }
}

///|
/// Get the estimated pi0.
pub fn Pi0EstResult::pi0(self : Pi0EstResult) -> Double {
  self.pi0
}

///|
/// Get the pi0 estimates at each lambda.
pub fn Pi0EstResult::pi0_lambda(self : Pi0EstResult) -> Array[Double] {
  self.pi0_lambda
}

///|
/// Get the lambda grid.
pub fn Pi0EstResult::lambdas(self : Pi0EstResult) -> Array[Double] {
  self.lambdas
}

///|
/// Get the smoother fitted values (empty for the bootstrap method).
pub fn Pi0EstResult::pi0_smooth(self : Pi0EstResult) -> Array[Double] {
  self.pi0_smooth
}

///|
/// Full q-value analysis result, mirroring the object returned by the
/// Bioconductor `qvalue()` function.
pub struct QValueFullResult {
  p_values : Array[Double]
  q_values : Array[Double]
  pi0 : Double
  lfdr : Array[Double]
  significant : Array[Bool]
  fdr_level : Double
  pi0_lambda : Array[Double]
  lambdas : Array[Double]
  pi0_smooth : Array[Double]
  pi0_method : String
  pfdr : Bool
}

///|
/// Create a new QValueFullResult.
pub fn QValueFullResult::new(
  p_values : Array[Double],
  q_values : Array[Double],
  pi0 : Double,
  lfdr : Array[Double],
  significant : Array[Bool],
  fdr_level : Double,
  pi0_lambda : Array[Double],
  lambdas : Array[Double],
  pi0_smooth : Array[Double],
  pi0_method : String,
  pfdr : Bool,
) -> QValueFullResult {
  QValueFullResult::{
    p_values,
    q_values,
    pi0,
    lfdr,
    significant,
    fdr_level,
    pi0_lambda,
    lambdas,
    pi0_smooth,
    pi0_method,
    pfdr,
  }
}

///|
/// Get the input p-values.
pub fn QValueFullResult::p_values(self : QValueFullResult) -> Array[Double] {
  self.p_values
}

///|
/// Get the computed q-values.
pub fn QValueFullResult::q_values(self : QValueFullResult) -> Array[Double] {
  self.q_values
}

///|
/// Get the estimated pi0.
pub fn QValueFullResult::pi0(self : QValueFullResult) -> Double {
  self.pi0
}

///|
/// Get the estimated local FDR values (empty if lfdr_out was false).
pub fn QValueFullResult::lfdr(self : QValueFullResult) -> Array[Double] {
  self.lfdr
}

///|
/// Get the significance calls at fdr.level (all false if no level given).
pub fn QValueFullResult::significant(self : QValueFullResult) -> Array[Bool] {
  self.significant
}

///|
/// Get the FDR level used for significance calls (0.0 if not set).
pub fn QValueFullResult::fdr_level(self : QValueFullResult) -> Double {
  self.fdr_level
}

///|
/// Get the pi0 estimates at each lambda (empty if pi0 supplied directly).
pub fn QValueFullResult::pi0_lambda(self : QValueFullResult) -> Array[Double] {
  self.pi0_lambda
}

///|
/// Get the lambda grid (empty if pi0 supplied directly).
pub fn QValueFullResult::lambdas(self : QValueFullResult) -> Array[Double] {
  self.lambdas
}

///|
/// Get the smoother fitted values (empty for bootstrap or supplied pi0).
pub fn QValueFullResult::pi0_smooth(self : QValueFullResult) -> Array[Double] {
  self.pi0_smooth
}

///|
/// Get the pi0 estimation method used ("smoother" or "bootstrap").
pub fn QValueFullResult::pi0_method(self : QValueFullResult) -> String {
  self.pi0_method
}

///|
/// Whether the pFDR finite-sample correction was used.
pub fn QValueFullResult::pfdr(self : QValueFullResult) -> Bool {
  self.pfdr
}

// ---------------------------------------------------------------------------
// Internal helpers: validation, lambda grid, matrix utilities
// ---------------------------------------------------------------------------

///|
/// Strict p-value validation mirroring R's stop() on out-of-range p-values.
fn _check_pvalues_strict(p : Array[Double]) -> Unit raise QValueError {
  if p.length() == 0 {
    raise QValueError("p-values must not be empty.")
  }
  let mut i = 0
  while i < p.length() {
    if p[i] < 0.0 || p[i] > 1.0 {
      raise QValueError("p-values not in valid range [0, 1].")
    }
    i = i + 1
  }
}

///|
/// Default lambda grid: seq(0.05, 0.95, 0.05) (19 values).
fn _default_lambdas() -> Array[Double] {
  let v : Array[Double] = Array(capacity=0)
  let mut i = 1
  while i <= 19 {
    v.push(i.to_double() * 0.05)
    i = i + 1
  }
  v
}

///|
/// pi0(lambda) = #{p >= lambda} / (m * (1 - lambda)) at each lambda.
/// Uses the >= convention of the official package.
fn _pi0_lambda_path(
  p : Array[Double],
  lambdas : Array[Double],
) -> Array[Double] {
  let m = p.length().to_double()
  let ll = lambdas.length()
  let out = Array::make(ll, 0.0)
  let mut j = 0
  while j < ll {
    let lam = lambdas[j]
    let mut cnt = 0
    let mut i = 0
    while i < p.length() {
      if p[i] >= lam {
        cnt = cnt + 1
      }
      i = i + 1
    }
    out[j] = cnt.to_double() / (m * (1.0 - lam))
    j = j + 1
  }
  out
}

///|
/// Allocate a dense r x c matrix of zeros (independent rows).
fn _mat_new(r : Int, c : Int) -> Array[Array[Double]] {
  let m : Array[Array[Double]] = Array(capacity=0)
  let mut i = 0
  while i < r {
    m.push(Array::make(c, 0.0))
    i = i + 1
  }
  m
}

///|
/// Deep-copy a matrix.
fn _mat_copy(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let m : Array[Array[Double]] = Array(capacity=0)
  let mut i = 0
  while i < a.length() {
    m.push(a[i].copy())
    i = i + 1
  }
  m
}

///|
/// Identity matrix of size n.
fn _mat_identity(n : Int) -> Array[Array[Double]] {
  let m = _mat_new(n, n)
  let mut i = 0
  while i < n {
    m[i][i] = 1.0
    i = i + 1
  }
  m
}

// ---------------------------------------------------------------------------
// Reinsch cubic smoothing spline with degrees-of-freedom control
//
// For knots x (sorted, general spacing), the smoothing spline fitted values
// are f_hat = (I + alpha K)^{-1} y with K = Q R^{-1} Q^T, where Q is the
// divided-difference matrix (n x (n-2)) and R the tridiagonal band matrix.
// Degrees of freedom: df(alpha) = tr(S) = 2 + sum_i 1/(1 + alpha mu_i) with
// mu_i the eigenvalues of the generalized problem Q^T Q v = mu R v.
// We solve the generalized problem via Cholesky of R + Jacobi rotations,
// then bisect on log(alpha) until df reaches the target (3 in pi0est,
// matching smooth.spline(df = 3)).
// ---------------------------------------------------------------------------

///|
/// Jacobi eigendecomposition of a symmetric matrix.
/// Returns eigenvalues sorted ascending and the eigenvector matrix whose
/// columns are the corresponding eigenvectors.
fn _jacobi_sym(
  a0 : Array[Array[Double]],
) -> (Array[Double], Array[Array[Double]]) {
  let m = a0.length()
  let a = _mat_copy(a0)
  let v = _mat_identity(m)
  let mut sweep = 0
  while sweep < 100 {
    let mut off = 0.0
    let mut p = 0
    let mut q = 1
    let mut found = false
    let mut i = 0
    while i < m {
      let mut j = i + 1
      while j < m {
        let av = a[i][j].abs()
        if av > off {
          off = av
          p = i
          q = j
          found = true
        }
        j = j + 1
      }
      i = i + 1
    }
    if !found {
      break
    }
    if off < 1.0e-14 {
      break
    }
    let apq = a[p][q]
    let tau = (a[q][q] - a[p][p]) / (2.0 * apq)
    let t = if tau >= 0.0 {
      1.0 / (tau + (1.0 + tau * tau).sqrt())
    } else {
      -1.0 / (-tau + (1.0 + tau * tau).sqrt())
    }
    let c = 1.0 / (1.0 + t * t).sqrt()
    let s = t * c
    a[p][p] = a[p][p] - t * apq
    a[q][q] = a[q][q] + t * apq
    a[p][q] = 0.0
    a[q][p] = 0.0
    let mut i = 0
    while i < m {
      if i != p && i != q {
        let aip = a[i][p]
        let aiq = a[i][q]
        a[i][p] = c * aip - s * aiq
        a[p][i] = a[i][p]
        a[i][q] = s * aip + c * aiq
        a[q][i] = a[i][q]
      }
      i = i + 1
    }
    i = 0
    while i < m {
      let vip = v[i][p]
      let viq = v[i][q]
      v[i][p] = c * vip - s * viq
      v[i][q] = s * vip + c * viq
      i = i + 1
    }
    sweep = sweep + 1
  }
  let ev0 = Array::make(m, 0.0)
  let mut i = 0
  while i < m {
    ev0[i] = if a[i][i] < 0.0 { 0.0 } else { a[i][i] }
    i = i + 1
  }
  let idx : Array[Int] = Array(capacity=0)
  i = 0
  while i < m {
    idx.push(i)
    i = i + 1
  }
  idx.sort_by(fn(a : Int, b : Int) -> Int {
    if ev0[a] < ev0[b] {
      -1
    } else if ev0[a] > ev0[b] {
      1
    } else {
      0
    }
  })
  let ev = Array::make(m, 0.0)
  let evec = _mat_new(m, m)
  let mut j = 0
  while j < m {
    ev[j] = ev0[idx[j]]
    let mut ii = 0
    while ii < m {
      evec[ii][j] = v[ii][idx[j]]
      ii = ii + 1
    }
    j = j + 1
  }
  (ev, evec)
}

///|
/// Degrees of freedom of the smoothing spline for a given alpha.
fn _smooth_df(ev : Array[Double], alpha : Double) -> Double {
  let mut d = 2.0
  let mut i = 0
  while i < ev.length() {
    d = d + 1.0 / (1.0 + alpha * ev[i])
    i = i + 1
  }
  d
}

///|
/// Fit a cubic smoothing spline (Reinsch form) with the given equivalent
/// degrees of freedom and return the fitted values at the knots.
fn _reinsch_df_fit(
  x : Array[Double],
  y : Array[Double],
  df_target : Double,
) -> Array[Double] {
  let n = x.length()
  if n < 4 || n != y.length() {
    return y.copy()
  }
  let mm = n - 2
  // spacing h_i = x[i+1] - x[i]
  let h = Array::make(n - 1, 0.0)
  let mut i = 0
  while i < n - 1 {
    h[i] = x[i + 1] - x[i]
    i = i + 1
  }
  // Q: n x (n-2) divided-difference matrix
  let qmat = _mat_new(n, mm)
  let mut j = 0
  while j < mm {
    qmat[j][j] = 1.0 / h[j]
    qmat[j + 1][j] = -(1.0 / h[j] + 1.0 / h[j + 1])
    qmat[j + 2][j] = 1.0 / h[j + 1]
    j = j + 1
  }
  // R: (n-2) x (n-2) tridiagonal band matrix
  let rmat = _mat_new(mm, mm)
  j = 0
  while j < mm {
    rmat[j][j] = (h[j] + h[j + 1]) / 3.0
    if j + 1 < mm {
      rmat[j][j + 1] = h[j + 1] / 6.0
      rmat[j + 1][j] = rmat[j][j + 1]
    }
    j = j + 1
  }
  // A = Q^T Q
  let amat = _mat_new(mm, mm)
  i = 0
  while i < n {
    let mut jj = 0
    while jj < mm {
      if qmat[i][jj] != 0.0 {
        let mut kk = 0
        while kk < mm {
          amat[jj][kk] = amat[jj][kk] + qmat[i][jj] * qmat[i][kk]
          kk = kk + 1
        }
      }
      jj = jj + 1
    }
    i = i + 1
  }
  // Cholesky of R (R is tridiagonal SPD, L is lower bidiagonal)
  let lmat = _mat_new(mm, mm)
  lmat[0][0] = rmat[0][0].sqrt()
  i = 1
  while i < mm {
    lmat[i][i - 1] = rmat[i][i - 1] / lmat[i - 1][i - 1]
    let diag = rmat[i][i] - lmat[i][i - 1] * lmat[i][i - 1]
    let safe_diag = if diag < 1.0e-300 { 1.0e-300 } else { diag }
    lmat[i][i] = safe_diag.sqrt()
    i = i + 1
  }
  // M = L^{-1} A, forward substitution per column
  let mmat = _mat_new(mm, mm)
  let mut col = 0
  while col < mm {
    let mut ii = 0
    while ii < mm {
      let mut s = amat[ii][col]
      let mut kk = 0
      while kk < ii {
        s = s - lmat[ii][kk] * mmat[kk][col]
        kk = kk + 1
      }
      mmat[ii][col] = s / lmat[ii][ii]
      ii = ii + 1
    }
    col = col + 1
  }
  // C = L^{-1} A L^{-T}; store via N = L^{-1} M^T, C = N^T (C symmetric)
  let cmat = _mat_new(mm, mm)
  col = 0
  while col < mm {
    let z = Array::make(mm, 0.0)
    let mut ii = 0
    while ii < mm {
      let mut s = mmat[col][ii]
      let mut kk = 0
      while kk < ii {
        s = s - lmat[ii][kk] * z[kk]
        kk = kk + 1
      }
      z[ii] = s / lmat[ii][ii]
      ii = ii + 1
    }
    let mut ii2 = 0
    while ii2 < mm {
      cmat[col][ii2] = z[ii2]
      ii2 = ii2 + 1
    }
    col = col + 1
  }
  let (ev, evec) = _jacobi_sym(cmat)
  // Bisect on log(alpha) so that df(alpha) = df_target
  let mut lo_alpha = 1.0e-12
  let mut hi_alpha = 1.0e12
  let mut iter = 0
  while iter < 200 {
    let mid = (lo_alpha * hi_alpha).sqrt()
    if _smooth_df(ev, mid) > df_target {
      lo_alpha = mid
    } else {
      hi_alpha = mid
    }
    iter = iter + 1
  }
  let alpha = (lo_alpha * hi_alpha).sqrt()
  // Q^T y
  let qty = Array::make(mm, 0.0)
  j = 0
  while j < mm {
    let mut s = 0.0
    i = 0
    while i < n {
      s = s + qmat[i][j] * y[i]
      i = i + 1
    }
    qty[j] = s
    j = j + 1
  }
  // z = L^{-1} Q^T y
  let z = Array::make(mm, 0.0)
  i = 0
  while i < mm {
    let mut s = qty[i]
    let mut kk = 0
    while kk < i {
      s = s - lmat[i][kk] * z[kk]
      kk = kk + 1
    }
    z[i] = s / lmat[i][i]
    i = i + 1
  }
  // z' = V diag(1/(mu + 1/alpha)) V^T z
  let zp = Array::make(mm, 0.0)
  j = 0
  while j < mm {
    let mut w = 0.0
    let mut ii = 0
    while ii < mm {
      w = w + evec[ii][j] * z[ii]
      ii = ii + 1
    }
    let wp = w / (ev[j] + 1.0 / alpha)
    let mut ii2 = 0
    while ii2 < mm {
      zp[ii2] = zp[ii2] + evec[ii2][j] * wp
      ii2 = ii2 + 1
    }
    j = j + 1
  }
  // gamma = L^{-T} z', back substitution
  let gamma = Array::make(mm, 0.0)
  i = mm - 1
  while i >= 0 {
    let mut s = zp[i]
    let mut kk = i + 1
    while kk < mm {
      s = s - lmat[kk][i] * gamma[kk]
      kk = kk + 1
    }
    gamma[i] = s / lmat[i][i]
    i = i - 1
  }
  // fitted values = y - Q gamma
  let fit = Array::make(n, 0.0)
  i = 0
  while i < n {
    let mut s = 0.0
    let mut jj = 0
    while jj < mm {
      s = s + qmat[i][jj] * gamma[jj]
      jj = jj + 1
    }
    fit[i] = y[i] - s
    i = i + 1
  }
  fit
}

// ---------------------------------------------------------------------------
// Normal distribution helpers (dnorm, qnorm via Acklam's approximation),
// type-7 quantile, and Gaussian KDE with Silverman/nrd0 bandwidth.
// ---------------------------------------------------------------------------

///|
/// Standard normal density phi(x).
fn _dnorm(x : Double) -> Double {
  @math.exp(-0.5 * x * x) * 0.3989422804014327
}

///|
/// Inverse standard normal CDF (Peter Acklam's rational approximation,
/// relative error < 1.2e-9).
fn _qnorm(p : Double) -> Double {
  let plow = 0.02425
  let phigh = 1.0 - plow
  if p < plow {
    let q = (-2.0 * @math.ln(p)).sqrt()
    let num = (
        (
          (
            (-0.00778489400243029 * q - 0.322396458041136) * q -
            2.40075827716184
          ) *
          q -
          2.54973253934373
        ) *
        q +
        4.37466414146497
      ) *
      q +
      2.93816398269878
    let den = (
        ((0.00778469570904146 * q + 0.32246712907004) * q + 2.445134137143) * q +
        3.75440866190742
      ) *
      q +
      1.0
    return num / den
  }
  if p > phigh {
    let q = (-2.0 * @math.ln(1.0 - p)).sqrt()
    let num = (
        (
          (
            (-0.00778489400243029 * q - 0.322396458041136) * q -
            2.40075827716184
          ) *
          q -
          2.54973253934373
        ) *
        q +
        4.37466414146497
      ) *
      q +
      2.93816398269878
    let den = (
        ((0.00778469570904146 * q + 0.32246712907004) * q + 2.445134137143) * q +
        3.75440866190742
      ) *
      q +
      1.0
    return -num / den
  }
  let q = p - 0.5
  let r = q * q
  let num = (
      (
        (
          ((-39.6968302866538 * r + 220.946098424521) * r - 275.928510446969) *
          r +
          138.357751867269
        ) *
        r -
        30.6647980661472
      ) *
      r +
      2.50662827745924
    ) *
    q
  let den = (
      (
        ((-54.4760987982241 * r + 161.585836858041) * r - 155.698979859887) * r +
        66.8013118877197
      ) *
      r -
      13.2806815528857
    ) *
    r +
    1.0
  num / den
}

///|
/// Type-7 quantile (R's default) of a sorted array: h = 1 + (n-1)*p with
/// linear interpolation between order statistics.
fn _quantile7_sorted(xs : Array[Double], p : Double) -> Double {
  let n = xs.length()
  if n == 0 {
    return 0.0
  }
  if n == 1 {
    return xs[0]
  }
  let h = 1.0 + (n.to_double() - 1.0) * p
  let fl = h.floor()
  let mut lo = fl.to_int()
  if lo < 1 {
    lo = 1
  }
  if lo >= n {
    return xs[n - 1]
  }
  let frac = h - fl
  xs[lo - 1] * (1.0 - frac) + xs[lo] * frac
}

///|
/// Gaussian kernel density evaluated at the data points. Mirrors R's
/// density(x, adjust=adj, kernel="gaussian") with bw.nrd0 bandwidth and a
/// 200-point evaluation grid (linear interpolation back to the data).
fn _kde_at_points(x : Array[Double], adj : Double) -> Array[Double] {
  let n = x.length()
  let out : Array[Double] = Array::make(n, 1.0)
  if n < 2 {
    return out
  }
  let nd = n.to_double()
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + x[i]
    i = i + 1
  }
  let mean = sum / nd
  let mut ss = 0.0
  i = 0
  while i < n {
    let d = x[i] - mean
    ss = ss + d * d
    i = i + 1
  }
  let sd = (ss / (nd - 1.0)).sqrt()
  let xs = x.copy()
  xs.sort()
  let iqr_s = (_quantile7_sorted(xs, 0.75) - _quantile7_sorted(xs, 0.25)) / 1.34
  let mut lo = if iqr_s < sd { iqr_s } else { sd }
  if lo == 0.0 {
    lo = sd
  }
  if lo == 0.0 {
    lo = x[0].abs()
  }
  if lo == 0.0 {
    lo = 1.0
  }
  let mut bw = 0.9 * lo * @math.pow(nd, -0.2) * adj
  if bw < 1.0e-9 || bw != bw {
    bw = 1.0
  }
  let mut xmin = x[0]
  let mut xmax = x[0]
  i = 1
  while i < n {
    if x[i] < xmin {
      xmin = x[i]
    }
    if x[i] > xmax {
      xmax = x[i]
    }
    i = i + 1
  }
  let ng = 200
  let gmin = xmin - 3.0 * bw
  let gmax = xmax + 3.0 * bw
  let step = (gmax - gmin) / (ng - 1).to_double()
  let norm = 0.3989422804014327 / (nd * bw)
  let grid : Array[Double] = Array::make(ng, 0.0)
  let mut k = 0
  while k < ng {
    let u = gmin + step * k.to_double()
    let mut s = 0.0
    let mut ii = 0
    while ii < n {
      let z = (u - x[ii]) / bw
      s = s + @math.exp(-0.5 * z * z)
      ii = ii + 1
    }
    grid[k] = norm * s
    k = k + 1
  }
  i = 0
  while i < n {
    let pos = (x[i] - gmin) / step
    let mut idx = pos.to_int()
    if idx < 0 {
      idx = 0
    }
    if idx > ng - 2 {
      idx = ng - 2
    }
    let frac = pos - idx.to_double()
    let d = grid[idx] * (1.0 - frac) + grid[idx + 1] * frac
    out[i] = if d < 1.0e-300 { 1.0e-300 } else { d }
    i = i + 1
  }
  out
}

// ---------------------------------------------------------------------------
// pi0est: proportion of true null p-values (smoother / bootstrap)
// ---------------------------------------------------------------------------

///|
/// Estimate the proportion of true null p-values pi0, mirroring the
/// Bioconductor qvalue `pi0est()` function.
///
/// pi0(lambda) = #{p >= lambda} / (m * (1 - lambda)) over the lambda grid.
/// The "smoother" method fits a cubic smoothing spline with df = 3 (Reinsch
/// form, matching R's smooth.spline(df=3)) and takes the fitted value at the
/// largest lambda; the "bootstrap" method uses the David Robinson closed-form
/// MSE minimization (no randomness).
///
/// Parameters:
///   - p_values:   Array of p-values in [0, 1]
///   - lambdas:    Tuning grid (default seq(0.05, 0.95, 0.05)); a single
///                 lambda gives the direct estimator
///   - method:     "smoother" (default) or "bootstrap"
///   - smooth_log: If true, smooth log(pi0(lambda)) then exponentiate
///
/// Returns:
///   - Pi0EstResult with pi0, pi0_lambda, lambdas, pi0_smooth
pub fn qvalue_pi0est(
  p_values : Array[Double],
  lambdas? : Array[Double] = _default_lambdas(),
  method_? : String = "smoother",
  smooth_log? : Bool = false,
) -> Pi0EstResult raise QValueError {
  _check_pvalues_strict(p_values)
  let lams = lambdas.copy()
  lams.sort()
  let ll = lams.length()
  if ll == 0 {
    raise QValueError("lambda must contain at least one value in [0, 1).")
  }
  if ll == 2 || ll == 3 {
    raise QValueError(
      "If length of lambda is greater than 1, you need at least 4 values.",
    )
  }
  if lams[0] < 0.0 || lams[ll - 1] >= 1.0 {
    raise QValueError("Lambda must be within [0, 1).")
  }
  let m = p_values.length()
  let pi0_lambda = _pi0_lambda_path(p_values, lams)
  if ll == 1 {
    let mut pi0 = pi0_lambda[0]
    if pi0 > 1.0 {
      pi0 = 1.0
    }
    if pi0 <= 0.0 {
      raise QValueError(
        "The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.",
      )
    }
    return Pi0EstResult::{ pi0, pi0_lambda, lambdas: lams, pi0_smooth: [], }
  }
  if method_ == "smoother" {
    let mut yfit : Array[Double] = []
    if smooth_log {
      let ylog = Array::make(ll, 0.0)
      let mut j = 0
      while j < ll {
        if pi0_lambda[j] <= 0.0 {
          raise QValueError(
            "smooth.log.pi0 requires pi0(lambda) > 0 at all lambdas; use smooth_log=false.",
          )
        }
        ylog[j] = @math.ln(pi0_lambda[j])
        j = j + 1
      }
      let fit_log = _reinsch_df_fit(lams, ylog, 3.0)
      yfit = Array::make(ll, 0.0)
      let mut j2 = 0
      while j2 < ll {
        yfit[j2] = @math.exp(fit_log[j2])
        j2 = j2 + 1
      }
    } else {
      yfit = _reinsch_df_fit(lams, pi0_lambda, 3.0)
    }
    let mut pi0 = yfit[ll - 1]
    if pi0 > 1.0 {
      pi0 = 1.0
    }
    if pi0 <= 0.0 {
      raise QValueError(
        "The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.",
      )
    }
    Pi0EstResult::{ pi0, pi0_lambda, lambdas: lams, pi0_smooth: yfit, }
  } else if method_ == "bootstrap" {
    let sorted_pi0 = pi0_lambda.copy()
    sorted_pi0.sort()
    let minpi0 = _quantile7_sorted(sorted_pi0, 0.1)
    let md = m.to_double()
    let mut best_mse = 0.0
    let mut best_pi0 = 0.0
    let mut have_best = false
    let mut j = 0
    while j < ll {
      let lam = lams[j]
      let mut wcount = 0
      let mut i = 0
      while i < m {
        if p_values[i] >= lam {
          wcount = wcount + 1
        }
        i = i + 1
      }
      let w = wcount.to_double()
      let one_minus = 1.0 - lam
      let diff = pi0_lambda[j] - minpi0
      let mse = w / (md * md * one_minus * one_minus) * (1.0 - w / md) +
        diff * diff
      if have_best == false || mse < best_mse {
        best_mse = mse
        best_pi0 = pi0_lambda[j]
        have_best = true
      } else if mse == best_mse && pi0_lambda[j] < best_pi0 {
        best_pi0 = pi0_lambda[j]
      }
      j = j + 1
    }
    let mut pi0 = best_pi0
    if pi0 > 1.0 {
      pi0 = 1.0
    }
    if pi0 <= 0.0 {
      raise QValueError(
        "The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.",
      )
    }
    Pi0EstResult::{ pi0, pi0_lambda, lambdas: lams, pi0_smooth: [], }
  } else {
    raise QValueError(
      "pi0 method must be one of \"smoother\" or \"bootstrap\".",
    )
  }
}

// ---------------------------------------------------------------------------
// lfdr: local false discovery rate
// ---------------------------------------------------------------------------

///|
/// Estimate local FDR values from p-values, mirroring the Bioconductor
/// qvalue `lfdr()` function.
///
/// A probit (default) or logit transformation is applied to the p-values;
/// the marginal density of the transformed values is estimated with a
/// Gaussian kernel density (Silverman/nrd0 bandwidth times `adj`), and
/// lfdr = pi0 * null_density / marginal_density. Results are optionally
/// truncated at 1 and made monotone (non-decreasing with p).
///
/// Parameters:
///   - p_values: Array of p-values in [0, 1]
///   - pi0:      Supplied pi0; estimated via qvalue_pi0est if omitted
///   - transf:   "probit" (default) or "logit"
///   - trunc:    Cap lfdr values above 1 at 1 (default true)
///   - monotone: Force non-decreasing lfdr in p (cummax, default true)
///   - adj:      Bandwidth multiplier for the density estimate (default 1.5)
///   - eps:      Tail threshold for the probit transform (default 1.0e-8)
pub fn qvalue_lfdr(
  p_values : Array[Double],
  pi0? : Double? = None,
  transf? : String = "probit",
  trunc? : Bool = true,
  monotone? : Bool = true,
  adj? : Double = 1.5,
  eps? : Double = 1.0e-8,
) -> Array[Double] raise QValueError {
  _check_pvalues_strict(p_values)
  let n = p_values.length()
  let pi0_est = match pi0 {
    Some(v) => {
      if v <= 0.0 || v > 1.0 {
        raise QValueError("pi0 is not (0,1].")
      }
      v
    }
    None => qvalue_pi0est(p_values).pi0
  }
  let x = Array::make(n, 0.0)
  let lfdr = Array::make(n, 0.0)
  if transf == "probit" {
    let mut i = 0
    while i < n {
      let mut pc = p_values[i]
      if pc < eps {
        pc = eps
      }
      if pc > 1.0 - eps {
        pc = 1.0 - eps
      }
      x[i] = _qnorm(pc)
      i = i + 1
    }
    let dens = _kde_at_points(x, adj)
    i = 0
    while i < n {
      lfdr[i] = pi0_est * _dnorm(x[i]) / dens[i]
      i = i + 1
    }
  } else if transf == "logit" {
    let mut i = 0
    while i < n {
      let pc = p_values[i]
      x[i] = @math.ln((pc + eps) / (1.0 - pc + eps))
      i = i + 1
    }
    let dens = _kde_at_points(x, adj)
    i = 0
    while i < n {
      let ex = @math.exp(x[i])
      let one_plus = 1.0 + ex
      let dx = ex / (one_plus * one_plus)
      lfdr[i] = pi0_est * dx / dens[i]
      i = i + 1
    }
  } else {
    raise QValueError("transf must be one of \"probit\" or \"logit\".")
  }
  if trunc {
    let mut i = 0
    while i < n {
      if lfdr[i] > 1.0 {
        lfdr[i] = 1.0
      }
      i = i + 1
    }
  }
  if monotone {
    let order : Array[Int] = Array(capacity=0)
    let mut i = 0
    while i < n {
      order.push(i)
      i = i + 1
    }
    order.sort_by(fn(a : Int, b : Int) -> Int {
      if p_values[a] < p_values[b] {
        -1
      } else if p_values[a] > p_values[b] {
        1
      } else {
        0
      }
    })
    let mut running = lfdr[order[0]]
    let mut k = 1
    while k < n {
      let v = lfdr[order[k]]
      if v > running {
        running = v
      }
      lfdr[order[k]] = running
      k = k + 1
    }
  }
  lfdr
}

// ---------------------------------------------------------------------------
// empPvals: empirical p-values from test statistics
// ---------------------------------------------------------------------------

///|
/// Calculate pooled empirical p-values from observed test statistics and a
/// vector of simulated null statistics, mirroring `empPvals(pool=TRUE)`.
///
/// All null statistics are pooled: the merged values are sorted decreasing
/// (observed values precede tied null values) and each observed statistic
/// gets p = (#null values greater than it) / m0, floored at 1/m0.
///
/// Parameters:
///   - stat:  Observed test statistics (larger = more extreme)
///   - stat0: Null test statistics (flattened pool)
pub fn empirical_pvalues(
  stat : Array[Double],
  stat0 : Array[Double],
) -> Array[Double] raise QValueError {
  let m = stat.length()
  let m0 = stat0.length()
  if m0 == 0 {
    raise QValueError("stat0 must contain at least one null statistic.")
  }
  let pvals : Array[Double] = Array::make(m, 0.0)
  if m == 0 {
    return pvals
  }
  // Merged entries: [value, is_observed (1.0/0.0), original index]
  let merged : Array[Array[Double]] = Array(capacity=0)
  let mut i = 0
  while i < m {
    merged.push([stat[i], 1.0, i.to_double()])
    i = i + 1
  }
  i = 0
  while i < m0 {
    merged.push([stat0[i], 0.0, -1.0])
    i = i + 1
  }
  merged.sort_by(fn(a : Array[Double], b : Array[Double]) -> Int {
    if a[0] > b[0] {
      -1
    } else if a[0] < b[0] {
      1
    } else if a[1] > b[1] {
      -1
    } else if a[1] < b[1] {
      1
    } else {
      0
    }
  })
  let floor_p = 1.0 / m0.to_double()
  let mut obs_count = 0
  let mut pos = 0
  let total = merged.length()
  while pos < total {
    let entry = merged[pos]
    if entry[1] == 1.0 {
      let mut pv = (pos.to_double() - obs_count.to_double()) / m0.to_double()
      if pv < floor_p {
        pv = floor_p
      }
      pvals[entry[2].to_int()] = pv
      obs_count = obs_count + 1
    }
    pos = pos + 1
  }
  pvals
}

///|
/// Calculate test-specific (non-pooled) empirical p-values, mirroring
/// `empPvals(pool=FALSE)`.
///
/// stat0 is a matrix with one row per observed statistic; for test i,
/// p_i = mean(stat0[i,] >= stat[i]), floored at 1 / ncol.
///
/// Parameters:
///   - stat:       Observed test statistics (length m)
///   - stat0_rows: Null statistics, m rows of equal length
pub fn empirical_pvalues_rows(
  stat : Array[Double],
  stat0_rows : Array[Array[Double]],
) -> Array[Double] raise QValueError {
  let m = stat.length()
  if stat0_rows.length() != m {
    raise QValueError("Number of rows of stat0 must equal length of stat.")
  }
  let pvals : Array[Double] = Array::make(m, 0.0)
  if m == 0 {
    return pvals
  }
  let ncol = stat0_rows[0].length()
  if ncol == 0 {
    raise QValueError(
      "Each row of stat0 must contain at least one null statistic.",
    )
  }
  let floor_p = 1.0 / ncol.to_double()
  let mut i = 0
  while i < m {
    let row = stat0_rows[i]
    if row.length() != ncol {
      raise QValueError("stat0 must be a matrix (ragged rows are not allowed).")
    }
    let mut cnt = 0
    let mut j = 0
    while j < ncol {
      if row[j] >= stat[i] {
        cnt = cnt + 1
      }
      j = j + 1
    }
    let mut pv = cnt.to_double() / ncol.to_double()
    if pv < floor_p {
      pv = floor_p
    }
    pvals[i] = pv
    i = i + 1
  }
  pvals
}

// ---------------------------------------------------------------------------
// qvalue: full q-value analysis with pfdr and fdr.level
// ---------------------------------------------------------------------------

///|
/// Compute q-values with the official Bioconductor qvalue conventions,
/// mirroring `qvalue()`.
///
/// For p-values sorted ascending, q_(j) = pi0 * pmin(1, cummin(term_j)) from
/// the largest p down, with term_j = p*m/j for the ordinary FDR and
/// term_j = p*m / (j * (1 - (1-p)^m)) for the pFDR finite-sample correction.
///
/// Parameters:
///   - p_values:  Array of p-values in [0, 1]
///   - fdr_level: FDR level for significance calls (optional, in (0, 1])
///   - pfdr:      Use the pFDR correction (default false)
///   - pi0:       Supplied pi0; estimated via qvalue_pi0est if omitted
///   - lambdas:   Lambda grid for pi0 estimation (default 0.05..0.95)
///   - method:    pi0 method "smoother" (default) or "bootstrap"
///   - smooth_log: Smooth log(pi0) in the smoother method (default false)
///   - lfdr_out:  Include local FDR values in the result (default true)
pub fn qvalue_analyze(
  p_values : Array[Double],
  fdr_level? : Double? = None,
  pfdr? : Bool = false,
  pi0? : Double? = None,
  lambdas? : Array[Double]? = None,
  method_? : String = "smoother",
  smooth_log? : Bool = false,
  lfdr_out? : Bool = true,
) -> QValueFullResult raise QValueError {
  _check_pvalues_strict(p_values)
  match fdr_level {
    Some(level) =>
      if level <= 0.0 || level > 1.0 {
        raise QValueError("'fdr.level' must be in (0, 1].")
      }
    None => ()
  }
  let n = p_values.length()
  let mut pi0_lambda : Array[Double] = []
  let mut lams : Array[Double] = []
  let mut pi0_smooth : Array[Double] = []
  let pi0_est = match pi0 {
    Some(v) => {
      if v <= 0.0 || v > 1.0 {
        raise QValueError("pi0 is not (0,1].")
      }
      v
    }
    None => {
      let est = match lambdas {
        Some(ls) => qvalue_pi0est(p_values, lambdas=ls, method_~, smooth_log~)
        None => qvalue_pi0est(p_values, method_~, smooth_log~)
      }
      pi0_lambda = est.pi0_lambda
      lams = est.lambdas
      pi0_smooth = est.pi0_smooth
      est.pi0
    }
  }
  let q_values = _compute_qvalues_official(p_values, pi0_est, pfdr)
  let mut level = 0.0
  let significant = Array::make(n, false)
  match fdr_level {
    Some(lv) => {
      level = lv
      let mut i = 0
      while i < n {
        significant[i] = q_values[i] <= lv
        i = i + 1
      }
    }
    None => ()
  }
  let lfdr_vals = if lfdr_out {
    qvalue_lfdr(p_values, pi0=Some(pi0_est))
  } else {
    []
  }
  QValueFullResult::{
    p_values: p_values.copy(),
    q_values,
    pi0: pi0_est,
    lfdr: lfdr_vals,
    significant,
    fdr_level: level,
    pi0_lambda,
    lambdas: lams,
    pi0_smooth,
    pi0_method: method_,
    pfdr,
  }
}

///|
/// Official q-value computation: sorted ascending, cummin of p*m/rank
/// (or the pFDR term) from the largest p down, clamped at 1, times pi0.
fn _compute_qvalues_official(
  p_values : Array[Double],
  pi0 : Double,
  pfdr : Bool,
) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return []
  }
  let indexed : Array[Array[Double]] = Array(capacity=0)
  let mut i = 0
  while i < n {
    indexed.push([p_values[i], i.to_double()])
    i = i + 1
  }
  indexed.sort_by(fn(a : Array[Double], b : Array[Double]) -> Int {
    if a[0] < b[0] {
      -1
    } else if a[0] > b[0] {
      1
    } else {
      0
    }
  })
  let md = n.to_double()
  let q_sorted : Array[Double] = Array::make(n, 1.0)
  let mut running = 1.0e300
  let mut k = n - 1
  while k >= 0 {
    let pv = indexed[k][0]
    let rank = (k + 1).to_double()
    let mut term = pv * md / rank
    if pfdr {
      let onemp = if pv >= 1.0 {
        0.0
      } else {
        @math.exp(md * @math.ln(1.0 - pv))
      }
      term = pv * md / (rank * (1.0 - onemp))
    }
    if term < running {
      running = term
    }
    let capped = if running > 1.0 { 1.0 } else { running }
    q_sorted[k] = pi0 * capped
    k = k - 1
  }
  let q_values : Array[Double] = Array::make(n, 0.0)
  let mut j = 0
  while j < n {
    q_values[indexed[j][1].to_int()] = q_sorted[j]
    j = j + 1
  }
  q_values
}
