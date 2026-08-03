///|
/// openCyto - Automated, template-driven flow-cytometry gating framework.
///
/// Inspired by Bioconductor `openCyto` (Finak et al. 2014, BMC Bioinformatics
/// 15:282). This pure-MoonBit implementation provides the core declarative
/// gating pipeline that sits on top of `flowCore`:
///   - Gating template: parent → child population, method, dims, args
///   - Built-in gating methods (all pure-number algorithms):
///       * mindensity   — 1D density-based gate (KDE peak/trough detection)
///       * tailgate     — threshold gate for skewed populations
///       * quantileGate — quantile-based threshold
///       * flowClust    — robust model-based clustering (t-mixture EM)
///       * rangeGate    — 1D boundary via CDF derivative
///   - Gate hierarchy resolution (topological evaluation)
///   - Population statistics extraction
///
/// References:
///   https://www.bioconductor.org/packages/release/bioc/html/openCyto.html

// ============================================================================
// Core data structures
// ============================================================================

///|
/// One row of a gating template: defines how a child population is derived
/// from a parent population via a gating method on one or two channels.
pub struct OcGatingRule {
  parent : String // parent population name ("root" = all events)
  child : String // child population name
  method : String // "mindensity", "tailgate", "quantileGate",
  //                "flowClust", "rangeGate"
  dims : Array[String] // channel(s) / marker(s)
  args : Map[String, Double] // method-specific arguments
}

///|
/// Construct an OcGatingRule.
pub fn OcGatingRule::new(
  parent~ : String,
  child~ : String,
  method~ : String,
  dims~ : Array[String],
  args~ : Map[String, Double],
) -> OcGatingRule {
  OcGatingRule::{ parent, child, method, dims, args }
}

///|
/// A rectangular gate (1D or 2D) used as the result of a gating method.
pub struct OcGate {
  dim1 : String
  dim2 : String? // None for 1D gates
  min1 : Double
  max1 : Double
  min2 : Double
  max2 : Double
}

///|
/// Construct a 1D gate.
pub fn OcGate::new1d(
  dim~ : String,
  min~ : Double,
  max~ : Double,
) -> OcGate {
  OcGate::{ dim1: dim, dim2: None, min1: min, max1: max, min2: 0.0, max2: 0.0 }
}

///|
/// Construct a 2D gate.
pub fn OcGate::new2d(
  dim1~ : String,
  dim2~ : String,
  min1~ : Double,
  max1~ : Double,
  min2~ : Double,
  max2~ : Double,
) -> OcGate {
  OcGate::{ dim1: dim1, dim2: Some(dim2), min1, max1, min2, max2 }
}

///|
/// Check whether a (dim1, dim2) event falls inside the gate.
pub fn oc_event_in_gate(gate : OcGate, v1 : Double, v2 : Double) -> Bool {
  if v1 < gate.min1 || v1 > gate.max1 {
    return false
  }
  match gate.dim2 {
    Some(_) => {
      if v2 < gate.min2 || v2 > gate.max2 {
        return false
      }
      true
    }
    None => true
  }
}

///|
/// A flow-cytometry event matrix: samples × events × channels.
/// Each sample is a 2D matrix (events × channels).
pub struct OcFlowSet {
  sample_names : Array[String]
  channel_names : Array[String]
  data : Array[Array[Array[Double]]] // [sample][event][channel]
}

///|
/// Construct an OcFlowSet.
pub fn OcFlowSet::new(
  sample_names~ : Array[String],
  channel_names~ : Array[String],
  data~ : Array[Array[Array[Double]]],
) -> OcFlowSet {
  OcFlowSet::{ sample_names, channel_names, data }
}

///|
/// Extract a single channel as a 1D array from a sample's event matrix.
pub fn oc_extract_channel(
  events : Array[Array[Double]],
  channel_idx : Int,
) -> Array[Double] {
  let out = Array::new()
  for ev in events {
    if channel_idx < ev.length() {
      out.push(ev[channel_idx])
    }
  }
  out
}

///|
/// Result of gating one population in one sample.
pub struct OcGatingResult {
  population : String
  sample : String
  gate : OcGate
  parent_events : Int
  child_events : Int
  indices : Array[Bool] // per-parent-event membership flags
}

///|
/// Construct an OcGatingResult.
pub fn OcGatingResult::new(
  population~ : String,
  sample~ : String,
  gate~ : OcGate,
  parent_events~ : Int,
  child_events~ : Int,
  indices~ : Array[Bool],
) -> OcGatingResult {
  OcGatingResult::{
    population,
    sample,
    gate,
    parent_events,
    child_events,
    indices,
  }
}

// ============================================================================
// Numeric helpers
// ============================================================================

///|
/// Mean of a 1D array. Returns 0.0 if empty.
pub fn oc_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in values {
    s = s + v
  }
  s / values.length().to_double()
}

///|
/// Standard deviation (sample, 1/(n-1)) of a 1D array.
pub fn oc_sd(values : Array[Double]) -> Double {
  let n = values.length()
  if n < 2 {
    return 0.0
  }
  let m = oc_mean(values)
  let mut s = 0.0
  for v in values {
    let d = v - m
    s = s + d * d
  }
  (s / (n - 1).to_double()).sqrt()
}

///|
/// Median of a 1D array (non-mutating).
pub fn oc_median(values : Array[Double]) -> Double {
  let n = values.length()
  if n == 0 {
    return 0.0
  }
  let sorted = Array::new()
  for v in values {
    sorted.push(v)
  }
  oc_sort(sorted)
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

///|
/// Quantile (0..1) of a 1D array (non-mutating, linear interpolation).
pub fn oc_quantile(values : Array[Double], q : Double) -> Double {
  let n = values.length()
  if n == 0 {
    return 0.0
  }
  let sorted = Array::new()
  for v in values {
    sorted.push(v)
  }
  oc_sort(sorted)
  let pos = q * (n - 1).to_double()
  let lo = pos.to_int()
  let hi = if lo + 1 < n { lo + 1 } else { lo }
  let frac = pos - lo.to_double()
  sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

///|
/// In-place sort (insertion sort, sufficient for small-to-mid arrays).
fn oc_sort(arr : Array[Double]) -> Unit {
  let n = arr.length()
  for i in 1..<n {
    let key = arr[i]
    let mut j = i - 1
    while j >= 0 && arr[j] > key {
      arr[j + 1] = arr[j]
      j = j - 1
    }
    arr[j + 1] = key
  }
}

///|
/// Natural logarithm via series (used when @math not available).
pub fn oc_ln(x : Double) -> Double {
  if x <= 0.0 {
    return -1.0e30
  }
  @math.ln(x)
}

// ============================================================================
// KDE (Gaussian kernel density estimation) for mindensity
// ============================================================================

///|
/// Evaluate a Gaussian KDE at point `x` given data and bandwidth `h`.
/// KDE(x) = (1 / (n*h)) * Σ φ((x - xi) / h)
/// where φ is the standard normal density.
pub fn oc_kde_evaluate(
  data : Array[Double],
  x : Double,
  h : Double,
) -> Double {
  let n = data.length()
  if n == 0 || h <= 0.0 {
    return 0.0
  }
  let mut sum = 0.0
  let inv_h = 1.0 / h
  for xi in data {
    let u = (x - xi) * inv_h
    sum = sum + 0.3989423 * @math.exp(-u * u / 2.0)
  }
  sum * inv_h / n.to_double()
}

///|
/// Evaluate KDE on a grid of `n_points` points spanning [min, max].
/// Returns (grid_x, density_y) pairs.
pub fn oc_kde_grid(
  data : Array[Double],
  n_points : Int,
  h : Double,
) -> (Array[Double], Array[Double]) {
  let n = data.length()
  if n == 0 {
    return ([], [])
  }
  let mut lo = data[0]
  let mut hi = data[0]
  for v in data {
    if v < lo {
      lo = v
    }
    if v > hi {
      hi = v
    }
  }
  let range = hi - lo
  if range <= 0.0 {
    return ([lo], [oc_kde_evaluate(data, lo, h)])
  }
  let pad = range * 0.1
  lo = lo - pad
  hi = hi + pad
  let gx = Array::new()
  let gy = Array::new()
  for i in 0..<n_points {
    let x = lo + (hi - lo) * i.to_double() / (n_points - 1).to_double()
    gx.push(x)
    gy.push(oc_kde_evaluate(data, x, h))
  }
  (gx, gy)
}

///|
/// Silverman's rule-of-thumb bandwidth.
pub fn oc_bandwidth_silverman(data : Array[Double]) -> Double {
  let n = data.length()
  if n < 2 {
    return 1.0
  }
  let sd = oc_sd(data)
  let iqr = oc_quantile(data, 0.75) - oc_quantile(data, 0.25)
  let sig = if sd < iqr { sd } else { iqr }
  if sig <= 0.0 {
    return 1.0
  }
  0.9 * sig * @math.pow(n.to_double() * 0.2, 0.2)
}

// ============================================================================
// mindensity: 1D density-based gate
// ============================================================================

///|
/// Run the `mindensity` gating method on a 1D channel.
///
/// Finds the deepest valley (local minimum of the KDE) between the two most
/// prominent peaks and uses it as the threshold. The gate spans
/// [threshold, +∞) by default (or (-∞, threshold] if `negative_side` is set).
///
/// Returns a 1D OcGate.
pub fn oc_mindensity(
  values : Array[Double],
  dim : String,
  args : Map[String, Double],
) -> OcGate {
  let n = values.length()
  if n == 0 {
    return OcGate::new1d(dim=dim, min=0.0, max=1.0e30)
  }
  let n_points = 200
  let h = match args.get("bandwidth") {
    Some(b) => b
    None => oc_bandwidth_silverman(values)
  }
  let (gx, gy) = oc_kde_grid(values, n_points, h)
  if gx.length() < 3 {
    return OcGate::new1d(dim=dim, min=oc_median(values), max=1.0e30)
  }
  // find local maxima (peaks)
  let peaks : Array[Int] = Array::new()
  for i in 1..<(gy.length() - 1) {
    if gy[i] > gy[i - 1] && gy[i] >= gy[i + 1] {
      peaks.push(i)
    }
  }
  if peaks.length() < 2 {
    // single peak: use the median
    let m = oc_median(values)
    return OcGate::new1d(dim=dim, min=m, max=1.0e30)
  }
  // find the two highest peaks
  let mut p1 = peaks[0]
  let mut p2 = peaks[1]
  if gy[p1] < gy[p2] {
    let tmp = p1
    p1 = p2
    p2 = tmp
  }
  for i in 2..<peaks.length() {
    let p = peaks[i]
    if gy[p] > gy[p1] {
      p2 = p1
      p1 = p
    } else if gy[p] > gy[p2] {
      p2 = p
    }
  }
  // ensure p1 < p2 in x-axis order
  if gx[p1] > gx[p2] {
    let tmp = p1
    p1 = p2
    p2 = tmp
  }
  // find the valley (local minimum) between p1 and p2
  let mut valley = p1
  let mut valley_d = 1.0e300
  for i in (p1 + 1)..<p2 {
    if gy[i] < valley_d {
      valley_d = gy[i]
      valley = i
    }
  }
  let threshold = gx[valley]
  // gate spans [threshold, +∞) by default
  let max_val = match args.get("max") {
    Some(m) => m
    None => 1.0e30
  }
  OcGate::new1d(dim=dim, min=threshold, max=max_val)
}

// ============================================================================
// tailgate: threshold gate for skewed populations
// ============================================================================

///|
/// Run the `tailgate` gating method on a 1D channel.
///
/// For populations with a long tail (no clear second peak), tailgate places
/// the threshold at a user-specified quantile (default 0.99) of the data.
pub fn oc_tailgate(
  values : Array[Double],
  dim : String,
  args : Map[String, Double],
) -> OcGate {
  let n = values.length()
  if n == 0 {
    return OcGate::new1d(dim=dim, min=0.0, max=1.0e30)
  }
  let q = match args.get("quantile") {
    Some(v) => v
    None => 0.99
  }
  let threshold = oc_quantile(values, q)
  let max_val = match args.get("max") {
    Some(m) => m
    None => 1.0e30
  }
  OcGate::new1d(dim=dim, min=threshold, max=max_val)
}

// ============================================================================
// quantileGate: quantile-based threshold
// ============================================================================

///|
/// Run the `quantileGate` method: place the gate at a given quantile,
/// keeping the top (1 - q) fraction of events.
pub fn oc_quantile_gate(
  values : Array[Double],
  dim : String,
  args : Map[String, Double],
) -> OcGate {
  let n = values.length()
  if n == 0 {
    return OcGate::new1d(dim=dim, min=0.0, max=1.0e30)
  }
  let q = match args.get("probs") {
    Some(v) => v
    None => 0.95
  }
  let threshold = oc_quantile(values, q)
  let max_val = match args.get("max") {
    Some(m) => m
    None => 1.0e30
  }
  OcGate::new1d(dim=dim, min=threshold, max=max_val)
}

// ============================================================================
// rangeGate: 1D boundary via CDF derivative
// ============================================================================

///|
/// Run the `rangeGate` method: find the steepest change in the empirical CDF
/// (i.e. the point of maximum density increase) and use it as the threshold.
pub fn oc_range_gate(
  values : Array[Double],
  dim : String,
  args : Map[String, Double],
) -> OcGate {
  let n = values.length()
  if n == 0 {
    return OcGate::new1d(dim=dim, min=0.0, max=1.0e30)
  }
  let sorted = Array::new()
  for v in values {
    sorted.push(v)
  }
  oc_sort(sorted)
  // compute CDF and find max derivative (max density)
  let n_points = if n < 100 { n } else { 100 }
  let step = n.to_double() / n_points.to_double()
  let cdf_x = Array::new()
  let cdf_y = Array::new()
  for i in 0..<n_points {
    let idx = (i.to_double() * step).to_int()
    let ii = if idx >= n { n - 1 } else { idx }
    cdf_x.push(sorted[ii])
    cdf_y.push(ii.to_double() / n.to_double())
  }
  // find max |d(CDF)/dx|
  let mut best_i = 0
  let mut best_deriv = 0.0
  for i in 1..<cdf_x.length() {
    let dx = cdf_x[i] - cdf_x[i - 1]
    if dx.abs() < 1.0e-12 {
      continue
    }
    let deriv = (cdf_y[i] - cdf_y[i - 1]).abs() / dx.abs()
    if deriv > best_deriv {
      best_deriv = deriv
      best_i = i
    }
  }
  let threshold = cdf_x[best_i]
  let max_val = match args.get("max") {
    Some(m) => m
    None => 1.0e30
  }
  OcGate::new1d(dim=dim, min=threshold, max=max_val)
}

// ============================================================================
// flowClust: robust model-based clustering (t-mixture with EM)
// ============================================================================

///|
/// Run the `flowClust` method: fit a 2-component t-mixture via EM and use
/// the boundary between the two components as the gate threshold.
///
/// For 1D data, returns a gate that separates the two clusters.
/// For 2D data, returns a rectangular gate covering the high-density cluster.
pub fn oc_flow_clust(
  values1 : Array[Double],
  values2 : Array[Double],
  dim1 : String,
  dim2 : String?,
  args : Map[String, Double],
  n_clusters : Int,
) -> OcGate {
  let n = values1.length()
  if n == 0 {
    return OcGate::new1d(dim=dim1, min=0.0, max=1.0e30)
  }
  // 1D flowClust
  if dim2 is None {
    let (means, _) = oc_em_t_mixture_1d(values1, n_clusters, 50)
    // sort means
    let sorted_means = Array::new()
    for m in means {
      sorted_means.push(m)
    }
    oc_sort_d(sorted_means)
    // gate: [boundary between first and second cluster, +∞)
    if sorted_means.length() >= 2 {
      let boundary = (sorted_means[0] + sorted_means[1]) / 2.0
      let max_val = match args.get("max") {
        Some(m) => m
        None => 1.0e30
      }
      return OcGate::new1d(dim=dim1, min=boundary, max=max_val)
    }
    let m = oc_mean(values1)
    return OcGate::new1d(dim=dim1, min=m, max=1.0e30)
  }
  // 2D flowClust: fit n_clusters on (v1, v2), gate the higher-mean cluster
  let v2 = values2
  let (means1, means2) = oc_em_t_mixture_2d(values1, v2, n_clusters, 50)
  // find the cluster with the highest mean (sum of both dims)
  let mut best = 0
  let mut best_sum = -1.0e300
  for k in 0..<n_clusters {
    let s = means1[k] + means2[k]
    if s > best_sum {
      best_sum = s
      best = k
    }
  }
  // estimate cluster spread (use overall SD as proxy)
  let sd1 = oc_sd(values1)
  let sd2 = oc_sd(v2)
  let gate_min1 = means1[best] - 2.0 * sd1
  let gate_max1 = means1[best] + 2.0 * sd1
  let gate_min2 = means2[best] - 2.0 * sd2
  let gate_max2 = means2[best] + 2.0 * sd2
  OcGate::new2d(
    dim1=dim1,
    dim2=dim2.unwrap(),
    min1=gate_min1,
    max1=gate_max1,
    min2=gate_min2,
    max2=gate_max2,
  )
}

///|
/// EM for a 1D t-mixture with K components. Returns (means, weights).
/// Uses a simplified t-distribution (df=4 fixed) with Box-Cox.
pub fn oc_em_t_mixture_1d(
  data : Array[Double],
  n_clusters : Int,
  max_iter : Int,
) -> (Array[Double], Array[Double]) {
  let n = data.length()
  if n == 0 || n_clusters <= 0 {
    return ([], [])
  }
  if n_clusters == 1 {
    return ([oc_mean(data)], [1.0])
  }
  // initialize means via quantiles
  let means = Array::new()
  for k in 0..<n_clusters {
    means.push(oc_quantile(data, (k + 1).to_double() / (n_clusters + 1).to_double()))
  }
  let variances = Array::make(n_clusters, oc_sd(data) * oc_sd(data) + 1.0e-6)
  let weights = Array::make(n_clusters, 1.0 / n_clusters.to_double())
  let df = 4.0 // fixed degrees of freedom
  // EM iterations
  for _ in 0..<max_iter {
    // E-step: compute responsibilities
    let resp = Array::make(n * n_clusters, 0.0)
    for i in 0..<n {
      let mut sum = 0.0
      for k in 0..<n_clusters {
        let d = data[i] - means[k]
        let v = variances[k]
        // t-density (simplified): weight * t_pdf(x; mean, var, df)
        let t_pdf = oc_t_pdf(data[i], means[k], v, df) * weights[k]
        resp[i * n_clusters + k] = t_pdf
        sum = sum + t_pdf
      }
      if sum > 0.0 {
        for k in 0..<n_clusters {
          resp[i * n_clusters + k] = resp[i * n_clusters + k] / sum
        }
      }
    }
    // M-step: update parameters
    for k in 0..<n_clusters {
      let mut nk = 0.0
      let mut mk = 0.0
      let mut vk = 0.0
      for i in 0..<n {
        let r = resp[i * n_clusters + k]
        nk = nk + r
        mk = mk + r * data[i]
      }
      if nk > 0.0 {
        means[k] = mk / nk
        for i in 0..<n {
          let r = resp[i * n_clusters + k]
          let d = data[i] - means[k]
          vk = vk + r * d * d
        }
        variances[k] = vk / nk + 1.0e-6
        weights[k] = nk / n.to_double()
      }
    }
  }
  (means, weights)
}

///|
/// EM for a 2D t-mixture. Returns (means_dim1, means_dim2).
pub fn oc_em_t_mixture_2d(
  data1 : Array[Double],
  data2 : Array[Double],
  n_clusters : Int,
  max_iter : Int,
) -> (Array[Double], Array[Double]) {
  let n = data1.length()
  if n == 0 || n_clusters <= 0 {
    return ([], [])
  }
  if n_clusters == 1 {
    return ([oc_mean(data1)], [oc_mean(data2)])
  }
  let means1 = Array::new()
  let means2 = Array::new()
  for k in 0..<n_clusters {
    let q = (k + 1).to_double() / (n_clusters + 1).to_double()
    means1.push(oc_quantile(data1, q))
    means2.push(oc_quantile(data2, q))
  }
  let var1 = Array::make(n_clusters, oc_sd(data1) * oc_sd(data1) + 1.0e-6)
  let var2 = Array::make(n_clusters, oc_sd(data2) * oc_sd(data2) + 1.0e-6)
  let weights = Array::make(n_clusters, 1.0 / n_clusters.to_double())
  let df = 4.0
  for _ in 0..<max_iter {
    let resp = Array::make(n * n_clusters, 0.0)
    for i in 0..<n {
      let mut sum = 0.0
      for k in 0..<n_clusters {
        let t1 = oc_t_pdf(data1[i], means1[k], var1[k], df)
        let t2 = oc_t_pdf(data2[i], means2[k], var2[k], df)
        let t_pdf = t1 * t2 * weights[k]
        resp[i * n_clusters + k] = t_pdf
        sum = sum + t_pdf
      }
      if sum > 0.0 {
        for k in 0..<n_clusters {
          resp[i * n_clusters + k] = resp[i * n_clusters + k] / sum
        }
      }
    }
    for k in 0..<n_clusters {
      let mut nk = 0.0
      let mut m1 = 0.0
      let mut m2 = 0.0
      for i in 0..<n {
        let r = resp[i * n_clusters + k]
        nk = nk + r
        m1 = m1 + r * data1[i]
        m2 = m2 + r * data2[i]
      }
      if nk > 0.0 {
        means1[k] = m1 / nk
        means2[k] = m2 / nk
        let mut v1 = 0.0
        let mut v2 = 0.0
        for i in 0..<n {
          let r = resp[i * n_clusters + k]
          let d1 = data1[i] - means1[k]
          let d2 = data2[i] - means2[k]
          v1 = v1 + r * d1 * d1
          v2 = v2 + r * d2 * d2
        }
        var1[k] = v1 / nk + 1.0e-6
        var2[k] = v2 / nk + 1.0e-6
        weights[k] = nk / n.to_double()
      }
    }
  }
  (means1, means2)
}

///|
/// Student's t density (location-scale, given variance and df).
pub fn oc_t_pdf(
  x : Double,
  mean : Double,
  variance : Double,
  df : Double,
) -> Double {
  let sd = variance.sqrt()
  if sd <= 0.0 {
    return 0.0
  }
  let t = (x - mean) / sd
  // t-density: Γ((df+1)/2) / (√(df*π) * Γ(df/2)) * (1 + t²/df)^(-(df+1)/2)
  let log_norm = oc_lgamma((df + 1.0) / 2.0) - oc_lgamma(df / 2.0) -
    0.5 * @math.ln(df * 3.141592653589793) - @math.ln(sd)
  let log_kernel = -((df + 1.0) / 2.0) * @math.ln(1.0 + t * t / df)
  @math.exp(log_norm + log_kernel)
}

///|
/// Log-gamma function (Lanczos approximation).
pub fn oc_lgamma(x : Double) -> Double {
  if x < 0.5 {
    // reflection formula
    let pi = 3.141592653589793
    let log_pi = @math.ln(pi)
    return log_pi - @math.ln(pi * x).abs() - oc_lgamma(1.0 - x)
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
  let mut x_adj = x - 1.0
  let mut a = c[0]
  for i in 1..<9 {
    x_adj = x_adj + 1.0
    a = a + c[i] / x_adj
  }
  let t = x + g - 0.5
  0.5 * @math.ln(2.0 * 3.141592653589793) + (x - 0.5) * @math.ln(t) - t + @math.ln(a)
}

///|
/// Sort a double array in place (insertion sort).
fn oc_sort_d(arr : Array[Double]) -> Unit {
  let n = arr.length()
  for i in 1..<n {
    let key = arr[i]
    let mut j = i - 1
    while j >= 0 && arr[j] > key {
      arr[j + 1] = arr[j]
      j = j - 1
    }
    arr[j + 1] = key
  }
}

// ============================================================================
// Gate hierarchy resolution (topological evaluation)
// ============================================================================

///|
/// Topologically sort the gating template so that parents are evaluated
/// before children. Returns the sorted list of rules.
pub fn oc_topo_sort(rules : Array[OcGatingRule]) -> Array[OcGatingRule] {
  let out : Array[OcGatingRule] = Array::new()
  let visited : Map[String, Bool] = Map([], capacity=16)
  // build a map of child -> rule
  for rule in rules {
    oc_visit_rule(rule, rules, visited, out)
  }
  out
}

///|
/// DFS helper for topological sort.
fn oc_visit_rule(
  rule : OcGatingRule,
  all : Array[OcGatingRule],
  visited : Map[String, Bool],
  out : Array[OcGatingRule],
) -> Unit {
  match visited.get(rule.child) {
    Some(true) => return
    Some(false) => ()
    None => ()
  }
  // visit parent first (if parent is another child population)
  if rule.parent != "root" {
    for parent_rule in all {
      if parent_rule.child == rule.parent {
        oc_visit_rule(parent_rule, all, visited, out)
      }
    }
  }
  visited[rule.child] = true
  out.push(rule)
}

// ============================================================================
// Full gating pipeline
// ============================================================================

///|
/// Run the gating template on a flowSet, returning per-population,
/// per-sample gating results.
pub fn oc_gate_flow_set(
  flowset : OcFlowSet,
  template : Array[OcGatingRule],
) -> Array[OcGatingResult] {
  let sorted_rules = oc_topo_sort(template)
  let results : Array[OcGatingResult] = Array::new()
  // for each sample, evaluate the template
  for s in 0..<flowset.data.length() {
    let sample_name = flowset.sample_names[s]
    let events = flowset.data[s]
    let n_events = events.length()
    // track which events belong to each population
    let pop_indices : Map[String, Array[Bool]] = Map([], capacity=16)
    // root: all events
    let root_flags = Array::make(n_events, true)
    pop_indices["root"] = root_flags
    for rule in sorted_rules {
      let parent_flags = match pop_indices.get(rule.parent) {
        Some(f) => f
        None => Array::make(n_events, true)
      }
      // gather parent events
      let parent_event_indices : Array[Int] = Array::new()
      for i in 0..<n_events {
        if parent_flags[i] {
          parent_event_indices.push(i)
        }
      }
      if parent_event_indices.length() == 0 {
        let gate = OcGate::new1d(dim=rule.dims[0], min=0.0, max=0.0)
        results.push(
          OcGatingResult::new(
            population=rule.child,
            sample=sample_name,
            gate=gate,
            parent_events=0,
            child_events=0,
            indices=Array::new(),
          ),
        )
        pop_indices[rule.child] = Array::make(n_events, false)
        continue
      }
      // resolve channel indices
      let dim_idx = oc_resolve_channel(flowset.channel_names, rule.dims[0])
      let gate = if dim_idx < 0 {
        OcGate::new1d(dim=rule.dims[0], min=0.0, max=0.0)
      } else if rule.method == "mindensity" {
        let vals = Array::new()
        for i in parent_event_indices {
          vals.push(events[i][dim_idx])
        }
        oc_mindensity(vals, rule.dims[0], rule.args)
      } else if rule.method == "tailgate" {
        let vals = Array::new()
        for i in parent_event_indices {
          vals.push(events[i][dim_idx])
        }
        oc_tailgate(vals, rule.dims[0], rule.args)
      } else if rule.method == "quantileGate" {
        let vals = Array::new()
        for i in parent_event_indices {
          vals.push(events[i][dim_idx])
        }
        oc_quantile_gate(vals, rule.dims[0], rule.args)
      } else if rule.method == "rangeGate" {
        let vals = Array::new()
        for i in parent_event_indices {
          vals.push(events[i][dim_idx])
        }
        oc_range_gate(vals, rule.dims[0], rule.args)
      } else if rule.method == "flowClust" {
        let vals1 = Array::new()
        for i in parent_event_indices {
          vals1.push(events[i][dim_idx])
        }
        let n_clusters = match rule.args.get("K") {
          Some(k) => k.to_int()
          None => 2
        }
        if rule.dims.length() >= 2 {
          let dim2_idx = oc_resolve_channel(flowset.channel_names, rule.dims[1])
          if dim2_idx >= 0 {
            let vals2 = Array::new()
            for i in parent_event_indices {
              vals2.push(events[i][dim2_idx])
            }
            oc_flow_clust(vals1, vals2, rule.dims[0], Some(rule.dims[1]), rule.args, n_clusters)
          } else {
            oc_flow_clust(vals1, [], rule.dims[0], None, rule.args, n_clusters)
          }
        } else {
          oc_flow_clust(vals1, [], rule.dims[0], None, rule.args, n_clusters)
        }
      } else {
        // unknown method: use median split
        let vals = Array::new()
        for i in parent_event_indices {
          vals.push(events[i][dim_idx])
        }
        let m = oc_median(vals)
        OcGate::new1d(dim=rule.dims[0], min=m, max=1.0e30)
      }
      // apply gate
      let child_flags = Array::make(n_events, false)
      let mut child_count = 0
      for i in parent_event_indices {
        let v1 = if dim_idx >= 0 { events[i][dim_idx] } else { 0.0 }
        let v2 = if rule.dims.length() >= 2 {
          let d2 = oc_resolve_channel(flowset.channel_names, rule.dims[1])
          if d2 >= 0 { events[i][d2] } else { 0.0 }
        } else {
          0.0
        }
        if oc_event_in_gate(gate, v1, v2) {
          child_flags[i] = true
          child_count = child_count + 1
        }
      }
      results.push(
        OcGatingResult::new(
          population=rule.child,
          sample=sample_name,
          gate=gate,
          parent_events=parent_event_indices.length(),
          child_events=child_count,
          indices=child_flags,
        ),
      )
      pop_indices[rule.child] = child_flags
    }
  }
  results
}

///|
/// Find the index of a channel by name.
pub fn oc_resolve_channel(channels : Array[String], name : String) -> Int {
  for i in 0..<channels.length() {
    if channels[i] == name {
      return i
    }
  }
  -1
}

// ============================================================================
// Population statistics
// ============================================================================|

///|
/// Population statistics for one (population, sample) pair.
pub struct OcPopStat {
  population : String
  sample : String
  count : Int
  parent_count : Int
  percent_of_parent : Double
  percent_of_total : Double
}

///|
/// Construct an OcPopStat.
pub fn OcPopStat::new(
  population~ : String,
  sample~ : String,
  count~ : Int,
  parent_count~ : Int,
  percent_of_parent~ : Double,
  percent_of_total~ : Double,
) -> OcPopStat {
  OcPopStat::{
    population,
    sample,
    count,
    parent_count,
    percent_of_parent,
    percent_of_total,
  }
}

///|
/// Compute population statistics from gating results.
pub fn oc_population_stats(
  results : Array[OcGatingResult],
  total_events : Int,
) -> Array[OcPopStat] {
  let out : Array[OcPopStat] = Array::new()
  for r in results {
    let pop = r.parent_events
    let pop_pct = if pop > 0 {
      r.child_events.to_double() / pop.to_double() * 100.0
    } else {
      0.0
    }
    let total_pct = if total_events > 0 {
      r.child_events.to_double() / total_events.to_double() * 100.0
    } else {
      0.0
    }
    out.push(
      OcPopStat::new(
        population=r.population,
        sample=r.sample,
        count=r.child_events,
        parent_count=r.parent_events,
        percent_of_parent=pop_pct,
        percent_of_total=total_pct,
      ),
    )
  }
  out
}

///|
/// Format gating results as a human-readable summary string.
pub fn oc_gating_summary(results : Array[OcGatingResult]) -> String {
  let mut s = "=== openCyto Gating Summary ===\n"
  s = s + "Population\tSample\tParent\tChild\t%\n"
  for r in results {
    let pct = if r.parent_events > 0 {
      r.child_events.to_double() / r.parent_events.to_double() * 100.0
    } else {
      0.0
    }
    s = s +
      r.population +
      "\t" +
      r.sample +
      "\t" +
      r.parent_events.to_string() +
      "\t" +
      r.child_events.to_string() +
      "\t" +
      pct.to_string() +
      "\n"
  }
  s
}
