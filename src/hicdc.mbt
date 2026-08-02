///|
/// hicdc - Hi-C / HiChIP direct caller for significant chromatin interactions.
/// Inspired by Bioconductor HiC-DC+ package.
///
/// Provides negative binomial background modeling of Hi-C contact frequencies,
/// significant interaction detection, TAD calling via directionality index,
/// and A/B compartment calling via PCA on contact correlation matrices.
///
/// References:
///   - Sahin M et al. (2021) HiC-DC+: systematic 3D interaction calls and
///     differential analysis for Hi-C and HiChIP. Nat Commun 12:3366.
///   - Cullen K et al. (2020) HiC-DC: scalable and fast inference of
///     significant chromatin contacts.

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// A single Hi-C contact between two genomic bins.
pub struct HiCContact {
  chr : String
  bin1 : Int // bin index (start position / bin_size)
  bin2 : Int
  mut count : Int
  distance : Int // |bin2 - bin1|
  mut gc_content : Double
  mut mappability : Double
}

///|
/// Create a new HiCContact.
pub fn HiCContact::new(
  chr : String,
  bin1 : Int,
  bin2 : Int,
  count : Int,
  bin_size : Int
) -> HiCContact {
  let distance = if bin2 > bin1 { bin2 - bin1 } else { bin1 - bin2 }
  HiCContact::{
    chr,
    bin1,
    bin2,
    count,
    distance,
    gc_content: 0.5,
    mappability: 1.0,
  }
}

///|
/// Set GC content for a contact.
pub fn HiCContact::set_gc(self : HiCContact, gc : Double) -> Unit {
  self.gc_content = gc
}

///|
/// Set mappability for a contact.
pub fn HiCContact::set_mappability(self : HiCContact, m : Double) -> Unit {
  self.mappability = m
}

///|
/// Accessor for distance.
pub fn HiCContact::distance(self : HiCContact) -> Int {
  self.distance
}

///|
/// Accessor for count.
pub fn HiCContact::count(self : HiCContact) -> Int {
  self.count
}

///|
/// Accessor for chromosome.
pub fn HiCContact::chr(self : HiCContact) -> String {
  self.chr
}

///|
/// Accessor for bin1.
pub fn HiCContact::bin1(self : HiCContact) -> Int {
  self.bin1
}

///|
/// Accessor for bin2.
pub fn HiCContact::bin2(self : HiCContact) -> Int {
  self.bin2
}

///|
/// Accessor for GC content.
pub fn HiCContact::gc_content(self : HiCContact) -> Double {
  self.gc_content
}

///|
/// Accessor for mappability.
pub fn HiCContact::mappability(self : HiCContact) -> Double {
  self.mappability
}

///|
/// Result of HiC-DC significance testing for one contact.
pub struct HiCDCResult {
  chr : String
  bin1 : Int
  bin2 : Int
  observed : Int
  expected : Double
  log_fc : Double
  p_value : Double
  mut fdr : Double
  mut significant : Bool
}

///|
/// Accessor for observed count.
pub fn HiCDCResult::observed(self : HiCDCResult) -> Int {
  self.observed
}

///|
/// Accessor for expected count.
pub fn HiCDCResult::expected(self : HiCDCResult) -> Double {
  self.expected
}

///|
/// Accessor for p-value.
pub fn HiCDCResult::p_value(self : HiCDCResult) -> Double {
  self.p_value
}

///|
/// Accessor for FDR.
pub fn HiCDCResult::fdr(self : HiCDCResult) -> Double {
  self.fdr
}

///|
/// Accessor for significance.
pub fn HiCDCResult::significant(self : HiCDCResult) -> Bool {
  self.significant
}

///|
/// Accessor for chromosome.
pub fn HiCDCResult::chr(self : HiCDCResult) -> String {
  self.chr
}

///|
/// Accessor for bin1.
pub fn HiCDCResult::bin1(self : HiCDCResult) -> Int {
  self.bin1
}

///|
/// Accessor for bin2.
pub fn HiCDCResult::bin2(self : HiCDCResult) -> Int {
  self.bin2
}

///|
/// Accessor for log fold-change.
pub fn HiCDCResult::log_fc(self : HiCDCResult) -> Double {
  self.log_fc
}

///|
/// A TAD boundary.
pub struct TADBoundary {
  chr : String
  bin : Int
  directionality_index : Double
}

///|
/// Accessor for bin.
pub fn TADBoundary::bin(self : TADBoundary) -> Int {
  self.bin
}

///|
/// Accessor for chromosome.
pub fn TADBoundary::chr(self : TADBoundary) -> String {
  self.chr
}

///|
/// Accessor for directionality index.
pub fn TADBoundary::directionality_index(self : TADBoundary) -> Double {
  self.directionality_index
}

///|
/// An A/B compartment call.
pub struct CompartmentCall {
  chr : String
  bin : Int
  pc1 : Double
  compartment : String // "A" or "B"
}

///|
/// Accessor for bin.
pub fn CompartmentCall::bin(self : CompartmentCall) -> Int {
  self.bin
}

///|
/// Accessor for chromosome.
pub fn CompartmentCall::chr(self : CompartmentCall) -> String {
  self.chr
}

///|
/// Accessor for PC1.
pub fn CompartmentCall::pc1(self : CompartmentCall) -> Double {
  self.pc1
}

///|
/// Accessor for compartment label.
pub fn CompartmentCall::compartment(self : CompartmentCall) -> String {
  self.compartment
}

// ===========================================================================
// Internal helpers
// ===========================================================================

fn hc_abs(x : Double) -> Double {
  if x < 0.0 {
    -x
  } else {
    x
  }
}

fn hc_sqrt(x : Double) -> Double {
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
    if hc_abs(diff) < 0.0000000001 {
      return next
    }
    guess = next
  }
  guess
}

fn hc_log(x : Double) -> Double {
  if x <= 0.0 {
    return -700.0
  }
  @math.ln(x)
}

fn hc_exp(x : Double) -> Double {
  if x > 700.0 {
    return 1.7976931348623157e+308
  }
  if x < -700.0 {
    return 0.0
  }
  @math.exp(x)
}

///|
/// Normal CDF approximation.
fn hc_normal_cdf(x : Double) -> Double {
  if x < -8.0 {
    return 0.0
  }
  if x > 8.0 {
    return 1.0
  }
  let t = 1.0 / (1.0 + 0.2316419 * hc_abs(x))
  let d = 0.3989423 * hc_exp(-x * x / 2.0)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
  if x >= 0.0 {
    1.0 - p
  } else {
    p
  }
}

///|
/// Lower incomplete gamma function P(a, x) via series expansion.
fn hc_lower_gamma_p(a : Double, x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x < a + 1.0 {
    // Series expansion.
    let mut term = 1.0 / a
    let mut sum = term
    let mut n = 1
    while n < 200 {
      term = term * x / (a + n.to_double())
      sum = sum + term
      if hc_abs(term) < hc_abs(sum) * 0.0000001 {
        break
      }
      n = n + 1
    }
    sum * hc_exp(-x + a * hc_log(x) - hc_log_gamma(a))
  } else {
    // Continued fraction for upper incomplete gamma Q(a, x).
    1.0 - hc_upper_gamma_q(a, x)
  }
}

///|
/// Upper incomplete gamma function Q(a, x) via continued fraction.
fn hc_upper_gamma_q(a : Double, x : Double) -> Double {
  if x <= 0.0 {
    return 1.0
  }
  let tiny = 0.000000000000000000000000000000000000001
  let mut b = x + 1.0 - a
  let mut c = 1.0 / tiny
  let mut d = 1.0 / b
  let mut h = d
  let mut i = 1
  while i < 200 {
    let an = -(i.to_double()) * (i.to_double() - a)
    b = b + 2.0
    d = an * d + b
    if hc_abs(d) < tiny {
      d = tiny
    }
    c = b + an / c
    if hc_abs(c) < tiny {
      c = tiny
    }
    d = 1.0 / d
    let del = d * c
    h = h * del
    if hc_abs(del - 1.0) < 0.0000001 {
      break
    }
    i = i + 1
  }
  h * hc_exp(-x + a * hc_log(x) - hc_log_gamma(a))
}

///|
/// Log-gamma function via Lanczos approximation.
fn hc_log_gamma(x : Double) -> Double {
  if x < 0.5 {
    return hc_log(3.141592653589793 / (hc_sin(3.141592653589793 * x) * hc_exp(hc_log_gamma(1.0 - x))))
  }
  let g = 7.0
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
  let mut a = c[0]
  let t = x + g - 0.5
  for i in 1..<9 {
    a = a + c[i] / (x + i.to_double() - 1.0)
  }
  0.5 * hc_log(2.0 * 3.141592653589793) + (x - 0.5) * hc_log(t) - t + hc_log(a)
}

///|
/// Sine function via Taylor series.
fn hc_sin(x : Double) -> Double {
  let mut v = x
  // Reduce to [-pi, pi]
  let pi = 3.141592653589793
  while v > pi {
    v = v - 2.0 * pi
  }
  while v < -pi {
    v = v + 2.0 * pi
  }
  let mut sum = v
  let mut term = v
  let v2 = v * v
  let mut i = 1
  while i < 20 {
    term = -term * v2 / ((2 * i).to_double() * (2 * i + 1).to_double())
    sum = sum + term
    if hc_abs(term) < 0.0000000000001 {
      break
    }
    i = i + 1
  }
  sum
}

// ===========================================================================
// Negative binomial background model
// ===========================================================================

///|
/// Fit a negative binomial GLM for Hi-C contact counts as a function of
/// genomic distance, GC content, and mappability.
///
/// Model: log(E[count]) = offset + beta_dist * f(distance) + beta_gc * gc + beta_map * mappability
/// where f(distance) = log(1 + distance) for distance-dependent decay.
///
/// Returns the fitted parameters: (intercept, beta_dist, beta_gc, beta_map, dispersion).
pub fn hicdc_fit_background(
  contacts : Array[HiCContact]
) -> (Double, Double, Double, Double, Double) {
  let n = contacts.length()
  if n == 0 {
    return (0.0, 0.0, 0.0, 0.0, 1.0)
  }
  // Compute total counts for offset.
  let mut total = 0.0
  for c in contacts {
    total = total + c.count.to_double()
  }
  let log_total = hc_log(total / n.to_double())
  // Design: [1, log(1+dist), gc-0.5, mappability-1]
  // IRLS for NB GLM.
  let n_params = 4
  let beta = Array::make(n_params, 0.0)
  beta[0] = log_total // start from mean
  // Build design matrix.
  let design : Array[Array[Double]] = Array::new()
  let y : Array[Double] = Array::new()
  for c in contacts {
    let row = [
      1.0,
      hc_log(1.0 + c.distance.to_double()),
      c.gc_content - 0.5,
      c.mappability - 1.0,
    ]
    design.push(row)
    y.push(c.count.to_double())
  }
  // Estimate initial dispersion.
  let mean_y = total / n.to_double()
  let mut dispersion = 0.1
  // IRLS iterations.
  for _ in 0..<50 {
    let mu = Array::make(n, 0.0)
    for i in 0..<n {
      let mut eta = 0.0
      for j in 0..<n_params {
        eta = eta + design[i][j] * beta[j]
      }
      mu[i] = hc_exp(eta)
      if mu[i] < 0.0000001 {
        mu[i] = 0.0000001
      }
    }
    // Update dispersion via method of moments.
    let mut sum_sq = 0.0
    for i in 0..<n {
      let d = y[i] - mu[i]
      sum_sq = sum_sq + d * d
    }
    let new_disp = (sum_sq / n.to_double() - mean_y) / (mean_y * mean_y)
    dispersion = if new_disp > 0.001 { new_disp } else { 0.001 }
    // Weights and working response.
    let w = Array::make(n, 0.0)
    let z = Array::make(n, 0.0)
    for i in 0..<n {
      let var = mu[i] + dispersion * mu[i] * mu[i]
      w[i] = mu[i] / var
      let eta = hc_log(mu[i]) + (y[i] - mu[i]) / mu[i]
      z[i] = eta
    }
    // Weighted least squares.
    let xtwx = Array::new()
    for _ in 0..<n_params {
      xtwx.push(Array::make(n_params, 0.0))
    }
    let xtwz = Array::make(n_params, 0.0)
    for i in 0..<n {
      for j in 0..<n_params {
        xtwz[j] = xtwz[j] + design[i][j] * w[i] * z[i]
        for k in 0..<n_params {
          xtwx[j][k] = xtwx[j][k] + design[i][j] * w[i] * design[i][k]
        }
      }
    }
    match hc_solve(xtwx, xtwz) {
      Some(new_beta) => {
        let mut converged = true
        for j in 0..<n_params {
          if hc_abs(new_beta[j] - beta[j]) > 0.001 {
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
  (beta[0], beta[1], beta[2], beta[3], dispersion)
}

///|
/// Predict expected counts for contacts given fitted background model.
pub fn hicdc_predict_expected(
  contacts : Array[HiCContact],
  intercept : Double,
  beta_dist : Double,
  beta_gc : Double,
  beta_map : Double
) -> Array[Double] {
  let expected : Array[Double] = Array::new()
  for c in contacts {
    let eta = intercept +
      beta_dist * hc_log(1.0 + c.distance.to_double()) +
      beta_gc * (c.gc_content - 0.5) +
      beta_map * (c.mappability - 1.0)
    let mu = hc_exp(eta)
    expected.push(mu)
  }
  expected
}

// ===========================================================================
// Significance testing
// ===========================================================================

///|
/// Test each contact for significance using the fitted NB background model.
/// Uses a z-score based on standardized residuals.
///
/// Returns results with p-values (uncorrected) and FDR (BH-corrected).
pub fn hicdc_test_significance(
  contacts : Array[HiCContact],
  intercept : Double,
  beta_dist : Double,
  beta_gc : Double,
  beta_map : Double,
  dispersion : Double
) -> Array[HiCDCResult] {
  let n = contacts.length()
  let results : Array[HiCDCResult] = Array::new()
  for c in contacts {
    let eta = intercept +
      beta_dist * hc_log(1.0 + c.distance.to_double()) +
      beta_gc * (c.gc_content - 0.5) +
      beta_map * (c.mappability - 1.0)
    let mu = hc_exp(eta)
    let mu = if mu < 0.0000001 { 0.0000001 } else { mu }
    let variance = mu + dispersion * mu * mu
    let std_dev = hc_sqrt(variance)
    // Standardized residual (z-score).
    let z = (c.count.to_double() - mu) / std_dev
    // One-sided p-value (test for enrichment).
    let p_value = 1.0 - hc_normal_cdf(z)
    let p_value = if p_value < 0.0 { 0.0 } else if p_value > 1.0 { 1.0 } else { p_value }
    let log_fc = hc_log(c.count.to_double() / mu)
    results.push(HiCDCResult::{
      chr: c.chr,
      bin1: c.bin1,
      bin2: c.bin2,
      observed: c.count,
      expected: mu,
      log_fc,
      p_value,
      fdr: p_value,
      significant: false,
    })
  }
  // BH-FDR correction.
  hc_bh_fdr(results)
}

///|
/// Benjamini-Hochberg FDR correction.
fn hc_bh_fdr(results : Array[HiCDCResult]) -> Array[HiCDCResult] {
  let n = results.length()
  if n == 0 {
    return results
  }
  let indices = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  // Sort by p-value ascending (bubble sort).
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
  // Enforce monotonicity.
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
/// Mark significant results at given FDR threshold.
pub fn hicdc_mark_significant(
  results : Array[HiCDCResult],
  fdr_threshold : Double
) -> Array[HiCDCResult] {
  for r in results {
    r.significant = r.fdr < fdr_threshold
  }
  results
}

// ===========================================================================
// TAD calling via directionality index
// ===========================================================================

///|
/// Compute directionality index (DI) for each bin along a chromosome.
/// DI measures whether downstream or upstream contacts dominate.
///
/// `contacts` is the array of intra-chromosomal contacts.
/// `n_bins` is the total number of bins.
/// `window` is the number of bins on each side to consider.
/// Returns array of DI values per bin.
pub fn hicdc_directionality_index(
  contacts : Array[HiCContact],
  n_bins : Int,
  window : Int
) -> Array[Double] {
  // Build per-bin downstream and upstream contact sums.
  let downstream = Array::make(n_bins, 0.0)
  let upstream = Array::make(n_bins, 0.0)
  for c in contacts {
    let b1 = c.bin1
    let b2 = c.bin2
    if b1 >= 0 && b1 < n_bins && b2 >= 0 && b2 < n_bins {
      let cnt = c.count.to_double()
      if b2 > b1 && b2 - b1 <= window {
        downstream[b1] = downstream[b1] + cnt
        upstream[b2] = upstream[b2] + cnt
      }
    }
  }
  // Compute DI: (down - up) / |down - up| * ((up + down - |down - up|) / (up + down))
  let di = Array::make(n_bins, 0.0)
  for i in 0..<n_bins {
    let d = downstream[i]
    let u = upstream[i]
    let diff = d - u
    let abs_diff = hc_abs(diff)
    let total = d + u
    if total > 0.0 {
      let sign = if diff > 0.0 { 1.0 } else if diff < 0.0 { -1.0 } else { 0.0 }
      di[i] = sign * ((total - abs_diff) / total) * (total / (total + 1.0))
    } else {
      di[i] = 0.0
    }
  }
  di
}

///|
/// Call TAD boundaries based on directionality index sign changes.
/// A boundary is detected where DI changes sign (from negative to positive).
pub fn hicdc_call_tads(
  contacts : Array[HiCContact],
  n_bins : Int,
  window : Int,
  chr : String
) -> Array[TADBoundary] {
  let di = hicdc_directionality_index(contacts, n_bins, window)
  let boundaries : Array[TADBoundary] = Array::new()
  for i in 1..<n_bins {
    // Detect sign change from negative to positive (boundary).
    if di[i - 1] < 0.0 && di[i] >= 0.0 {
      boundaries.push(TADBoundary::{ chr, bin: i, directionality_index: di[i] })
    }
  }
  boundaries
}

// ===========================================================================
// A/B compartment calling
// ===========================================================================

///|
/// Call A/B compartments from a Hi-C contact matrix.
/// Computes Pearson correlation matrix of bin-level contact profiles,
/// then applies PCA. The first principal component sign determines
/// A (positive) or B (negative) compartment.
///
/// `contacts` is the array of intra-chromosomal contacts.
/// `n_bins` is the number of bins.
/// Returns compartment calls for each bin.
pub fn hicdc_call_compartments(
  contacts : Array[HiCContact],
  n_bins : Int,
  chr : String
) -> Array[CompartmentCall] {
  // Build contact matrix.
  let matrix : Array[Array[Double]] = Array::new()
  for _ in 0..<n_bins {
    matrix.push(Array::make(n_bins, 0.0))
  }
  for c in contacts {
    let b1 = c.bin1
    let b2 = c.bin2
    if b1 >= 0 && b1 < n_bins && b2 >= 0 && b2 < n_bins {
      matrix[b1][b2] = matrix[b1][b2] + c.count.to_double()
      if b1 != b2 {
        matrix[b2][b1] = matrix[b2][b1] + c.count.to_double()
      }
    }
  }
  // Compute Pearson correlation matrix.
  let corr = Array::new()
  for i in 0..<n_bins {
    let row = Array::make(n_bins, 0.0)
    let mean_i = hc_mean(matrix[i])
    let std_i = hc_std(matrix[i], mean_i)
    let std_i = if std_i < 0.0000001 { 0.0000001 } else { std_i }
    for j in 0..<n_bins {
      let mean_j = hc_mean(matrix[j])
      let std_j = hc_std(matrix[j], mean_j)
      let std_j = if std_j < 0.0000001 { 0.0000001 } else { std_j }
      let mut cov = 0.0
      for k in 0..<n_bins {
        cov = cov + (matrix[i][k] - mean_i) * (matrix[j][k] - mean_j)
      }
      cov = cov / n_bins.to_double()
      row[j] = cov / (std_i * std_j)
    }
    corr.push(row)
  }
  // Power iteration for first principal component.
  let pc1 = hc_power_iteration(corr, n_bins, 100)
  // Assign compartments based on sign.
  let result : Array[CompartmentCall] = Array::new()
  for i in 0..<n_bins {
    let comp = if pc1[i] >= 0.0 { "A" } else { "B" }
    result.push(CompartmentCall::{ chr, bin: i, pc1: pc1[i], compartment: comp })
  }
  result
}

///|
/// Power iteration to find the first eigenvector of a matrix.
fn hc_power_iteration(
  matrix : Array[Array[Double]],
  n : Int,
  max_iter : Int
) -> Array[Double] {
  let v = Array::make(n, 1.0)
  // Normalize initial vector.
  let mut norm = 0.0
  for i in 0..<n {
    norm = norm + v[i] * v[i]
  }
  norm = hc_sqrt(norm)
  if norm > 0.0 {
    for i in 0..<n {
      v[i] = v[i] / norm
    }
  }
  for _ in 0..<max_iter {
    let new_v = Array::make(n, 0.0)
    for i in 0..<n {
      let mut s = 0.0
      for j in 0..<n {
        s = s + matrix[i][j] * v[j]
      }
      new_v[i] = s
    }
    // Normalize.
    let mut norm = 0.0
    for i in 0..<n {
      norm = norm + new_v[i] * new_v[i]
    }
    norm = hc_sqrt(norm)
    if norm > 0.0 {
      for i in 0..<n {
        new_v[i] = new_v[i] / norm
      }
    }
    // Check convergence.
    let mut converged = true
    for i in 0..<n {
      if hc_abs(new_v[i] - v[i]) > 0.0001 {
        converged = false
      }
      v[i] = new_v[i]
    }
    if converged {
      break
    }
  }
  v
}

fn hc_mean(arr : Array[Double]) -> Double {
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

fn hc_std(arr : Array[Double], mean : Double) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    let d = v - mean
    s = s + d * d
  }
  hc_sqrt(s / n.to_double())
}

// ===========================================================================
// Linear algebra helper
// ===========================================================================

fn hc_solve(matrix : Array[Array[Double]], b : Array[Double]) -> Array[Double]? {
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
    let mut max_val = hc_abs(m[col][col])
    let mut r = col + 1
    while r < n {
      let val = hc_abs(m[r][col])
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

// ===========================================================================
// Differential analysis
// ===========================================================================

///|
/// Differential interaction analysis between two conditions.
/// Compares observed counts to expected counts from a shared background
/// model, and tests for condition-specific enrichment.
///
/// `contacts1` and `contacts2` are contacts from two conditions.
/// Returns results for contacts present in both conditions.
pub fn hicdc_differential(
  contacts1 : Array[HiCContact],
  contacts2 : Array[HiCContact],
  fdr_threshold : Double
) -> Array[HiCDCResult] {
  // Build lookup for condition 2 contacts.
  let map2 : Map[String, HiCContact] = Map([], capacity=contacts2.length())
  for c in contacts2 {
    let key = c.chr + "_" + c.bin1.to_string() + "_" + c.bin2.to_string()
    map2.set(key, c)
  }
  // Fit background on pooled data.
  let pooled : Array[HiCContact] = Array::new()
  for c in contacts1 {
    pooled.push(c)
  }
  for c in contacts2 {
    pooled.push(c)
  }
  let (intercept, beta_dist, beta_gc, beta_map, dispersion) = hicdc_fit_background(pooled)
  // Test each contact from condition 1 for differential enrichment.
  let results : Array[HiCDCResult] = Array::new()
  for c1 in contacts1 {
    let key = c1.chr + "_" + c1.bin1.to_string() + "_" + c1.bin2.to_string()
    let c2_opt = map2.get(key)
    if c2_opt is Some(c2) {
      let eta = intercept +
        beta_dist * hc_log(1.0 + c1.distance.to_double()) +
        beta_gc * (c1.gc_content - 0.5) +
        beta_map * (c1.mappability - 1.0)
      let mu = hc_exp(eta)
      let mu = if mu < 0.0000001 { 0.0000001 } else { mu }
      let variance = mu + dispersion * mu * mu
      let std_dev = hc_sqrt(variance)
      // Difference in counts.
      let diff = (c1.count - c2.count).to_double()
      let z = diff / (std_dev * hc_sqrt(2.0))
      let p_value = 2.0 * (1.0 - hc_normal_cdf(hc_abs(z)))
      let p_value = if p_value < 0.0 { 0.0 } else if p_value > 1.0 { 1.0 } else { p_value }
      let c2_count = c2.count.to_double()
      let denom = if c2_count > 0.0000001 { c2_count } else { 0.0000001 }
      let log_fc = hc_log(c1.count.to_double() / denom)
      results.push(HiCDCResult::{
        chr: c1.chr,
        bin1: c1.bin1,
        bin2: c1.bin2,
        observed: c1.count,
        expected: mu,
        log_fc,
        p_value,
        fdr: p_value,
        significant: false,
      })
    }
  }
  let _ = hc_bh_fdr(results)
  let _ = hicdc_mark_significant(results, fdr_threshold)
  results
}

// ===========================================================================
// Full pipeline
// ===========================================================================

///|
/// Run the complete HiC-DC pipeline: fit background, test significance,
/// and mark significant interactions.
pub fn hicdc_run(
  contacts : Array[HiCContact],
  fdr_threshold : Double
) -> Array[HiCDCResult] {
  let (intercept, beta_dist, beta_gc, beta_map, dispersion) = hicdc_fit_background(contacts)
  let results = hicdc_test_significance(
    contacts, intercept, beta_dist, beta_gc, beta_map, dispersion
  )
  hicdc_mark_significant(results, fdr_threshold)
}

// ===========================================================================
// Sample data
// ===========================================================================

///|
/// Create sample Hi-C contact data for testing.
/// Simulates contacts on a single chromosome with 20 bins.
/// Includes distance-dependent decay and some enriched (loop) contacts.
pub fn hicdc_sample_data() -> Array[HiCContact] {
  let contacts : Array[HiCContact] = Array::new()
  let chr = "chr1"
  let bin_size = 1 // simplified bin index
  let n_bins = 20
  // Generate intra-chromosomal contacts with distance decay.
  let mut rng_state = 7777
  for b1 in 0..<n_bins {
    for b2 in b1..<n_bins {
      let distance = b2 - b1
      // Distance-dependent expected count (decays with distance).
      let base_count = 50.0 / (1.0 + distance.to_double())
      let (s, r) = hc_lcg_rand(rng_state)
      rng_state = s
      let noise = r * 5.0
      let count = (base_count + noise).to_int()
      let count = if count < 0 { 0 } else { count }
      let contact = HiCContact::new(chr, b1, b2, count, bin_size)
      contacts.push(contact)
    }
  }
  // Add some enriched contacts (loops).
  let loops = [(2, 8), (5, 14), (10, 18)]
  for pair in loops {
    let b1 = pair.0
    let b2 = pair.1
    // Find and update the contact count.
    for c in contacts {
      if c.bin1 == b1 && c.bin2 == b2 {
        c.count = c.count + 30
      }
    }
  }
  contacts
}

///|
/// LCG random for sample data.
fn hc_lcg_rand(state : Int) -> (Int, Double) {
  let s = (state * 1103515245 + 12345) & 0x7FFFFFFF
  (s, (s % 1000000).to_double() / 1000000.0)
}
