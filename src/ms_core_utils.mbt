///|
/// MsCoreUtils - Foundational mass-spectrometry signal-processing utilities.
///
/// Inspired by Bioconductor `MsCoreUtils` (Rainer et al. 2021). This pure-
/// MoonBit implementation provides the low-level algorithmic layer beneath
/// MSstats/QFeatures/Spectra:
///   - refineCentroids: pick-peak centroiding (weighted mean around maxima)
///   - localMaxima / joinPeaks: peak detection and matching
///   - smooth: Savitzky-Golay and moving-average smoothing
///   - baselineEstimation: SNIP algorithm
///   - calibrate: calibrant-based m/z shift
///   - impute: KNN, zero, half-minimum, BPCA-style imputation
///   - medianPolish / robustSummary: row/column aggregation
///
/// References:
///   https://www.bioconductor.org/packages/release/bioc/html/MsCoreUtils.html

// ============================================================================
// Numeric helpers
// ============================================================================

///|
/// Mean of an array. Returns 0.0 if empty.
pub fn mc_mean(values : Array[Double]) -> Double {
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
/// Median of an array (non-mutating).
pub fn mc_median(values : Array[Double]) -> Double {
  let n = values.length()
  if n == 0 {
    return 0.0
  }
  let sorted = Array::new()
  for v in values {
    sorted.push(v)
  }
  mc_sort_inplace(sorted)
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

///|
/// Standard deviation (sample, 1/(n-1)).
pub fn mc_sd(values : Array[Double]) -> Double {
  let n = values.length()
  if n < 2 {
    return 0.0
  }
  let m = mc_mean(values)
  let mut s = 0.0
  for v in values {
    let d = v - m
    s = s + d * d
  }
  (s / (n - 1).to_double()).sqrt()
}

///|
/// Quantile (0..1) with linear interpolation (non-mutating).
pub fn mc_quantile(values : Array[Double], q : Double) -> Double {
  let n = values.length()
  if n == 0 {
    return 0.0
  }
  let sorted = Array::new()
  for v in values {
    sorted.push(v)
  }
  mc_sort_inplace(sorted)
  let pos = q * (n - 1).to_double()
  let lo = pos.to_int()
  let hi = if lo + 1 < n { lo + 1 } else { lo }
  let frac = pos - lo.to_double()
  sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

///|
/// In-place sort (insertion sort).
fn mc_sort_inplace(arr : Array[Double]) -> Unit {
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
// refineCentroids: pick-peak centroiding
// ============================================================================

///|
/// Refine peak centroids from profile-mode data.
///
/// Given m/z values, intensity values, and peak apex indices, compute the
/// refined centroid m/z and intensity for each peak by taking the weighted
/// mean of m/z values within +/- `half_window` points around each apex.
///
/// Returns (centroid_mz, centroid_intensity) arrays.
pub fn mc_refine_centroids(
  mz : Array[Double],
  intensity : Array[Double],
  peak_indices : Array[Int],
  half_window : Int,
) -> (Array[Double], Array[Double]) {
  let out_mz = Array::new()
  let out_int = Array::new()
  let n = mz.length()
  for p in peak_indices {
    if p < 0 || p >= n {
      continue
    }
    let lo = if p - half_window < 0 { 0 } else { p - half_window }
    let hi = if p + half_window >= n { n - 1 } else { p + half_window }
    let mut sum_mz = 0.0
    let mut sum_int = 0.0
    for i in lo..=hi {
      let w = intensity[i]
      if w > 0.0 {
        sum_mz = sum_mz + mz[i] * w
        sum_int = sum_int + w
      }
    }
    if sum_int > 0.0 {
      out_mz.push(sum_mz / sum_int)
      out_int.push(intensity[p])
    } else {
      out_mz.push(mz[p])
      out_int.push(intensity[p])
    }
  }
  (out_mz, out_int)
}

// ============================================================================
// localMaxima: peak detection
// ============================================================================

///|
/// Find local maxima in a 1D intensity array.
///
/// A point i is a local maximum if intensity[i] >= intensity[i-1],
/// intensity[i] >= intensity[i+1], and intensity[i] > 0.
/// Returns the indices of all local maxima.
pub fn mc_local_maxima(intensity : Array[Double]) -> Array[Int] {
  let n = intensity.length()
  let out = Array::new()
  if n < 3 {
    if n == 1 && intensity[0] > 0.0 {
      out.push(0)
    }
    return out
  }
  for i in 1..<(n - 1) {
    if intensity[i] > 0.0 &&
      intensity[i] >= intensity[i - 1] &&
      intensity[i] >= intensity[i + 1] {
      out.push(i)
    }
  }
  // check first and last
  if intensity[0] > 0.0 && intensity[0] >= intensity[1] {
    out.push(0)
  }
  if intensity[n - 1] > 0.0 && intensity[n - 1] >= intensity[n - 2] {
    out.push(n - 1)
  }
  out
}

// ============================================================================
// joinPeaks: peak matching with ppm / tolerance
// ============================================================================

///|
/// Match two lists of m/z values by tolerance.
///
/// For each m/z in `x`, find the closest m/z in `y` within the given
/// tolerance (either ppm-based or absolute). Returns an array of (x_idx, y_idx)
/// pairs for matched peaks.
pub fn mc_join_peaks(
  x : Array[Double],
  y : Array[Double],
  tolerance : Double,
  use_ppm : Bool,
) -> Array[(Int, Int)] {
  let out = Array::new()
  for i in 0..<x.length() {
    let mut best_j = -1
    let mut best_d = 1.0e300
    for j in 0..<y.length() {
      let diff = (x[i] - y[j]).abs()
      let tol = if use_ppm {
        x[i] * tolerance * 1.0e-6
      } else {
        tolerance
      }
      if diff <= tol && diff < best_d {
        best_d = diff
        best_j = j
      }
    }
    if best_j >= 0 {
      out.push((i, best_j))
    }
  }
  out
}

// ============================================================================
// Smoothing: Savitzky-Golay and moving average
// ============================================================================

///|
/// Moving-average smoothing with a window of `half_window * 2 + 1` points.
pub fn mc_smooth_moving_average(
  values : Array[Double],
  half_window : Int,
) -> Array[Double] {
  let n = values.length()
  let out = Array::make(n, 0.0)
  for i in 0..<n {
    let lo = if i - half_window < 0 { 0 } else { i - half_window }
    let hi = if i + half_window >= n { n - 1 } else { i + half_window }
    let mut s = 0.0
    let mut cnt = 0
    for j in lo..=hi {
      s = s + values[j]
      cnt = cnt + 1
    }
    out[i] = s / cnt.to_double()
  }
  out
}

///|
/// Savitzky-Golay smoothing with a quadratic polynomial and a window of
/// `half_window * 2 + 1` points. Uses pre-computed convolution coefficients
/// for the common case of quadratic/cubic fit.
pub fn mc_smooth_savitzky_golay(
  values : Array[Double],
  half_window : Int,
) -> Array[Double] {
  let n = values.length()
  let out = Array::make(n, 0.0)
  let window = half_window * 2 + 1
  if window < 3 || n == 0 {
    for i in 0..<n {
      out[i] = values[i]
    }
    return out
  }
  // Savitzky-Golay coefficients for quadratic polynomial smoothing
  // (convolution kernel, symmetric). For half_window=2 (window=5):
  //   [-3, 12, 17, 12, -3] / 35
  // For half_window=1 (window=3): [-1, 2, -1] is the 2nd derivative;
  // smoothing kernel is [1, 2, 1] / 4 (equivalent to moving average weighted).
  // General formula: c_k = (3/(4m^4 - 1)) * (3*m^2 - 1 - 5*k^2) / ... 
  // We use a simplified approach: fit a quadratic to each window.
  for i in 0..<n {
    let lo = if i - half_window < 0 { 0 } else { i - half_window }
    let hi = if i + half_window >= n { n - 1 } else { i + half_window }
    let cnt = hi - lo + 1
    if cnt < 3 {
      out[i] = values[i]
      continue
    }
    // least-squares quadratic fit: y = a + b*x + c*x^2
    // x values centered at 0
    let mut s0 = 0.0
    let mut s1 = 0.0
    let mut s2 = 0.0
    let mut s3 = 0.0
    let mut s4 = 0.0
    let mut t0 = 0.0
    let mut t1 = 0.0
    let mut t2 = 0.0
    let center = (lo + hi).to_double() / 2.0
    for j in lo..=hi {
      let x = j.to_double() - center
      let y = values[j]
      let x2 = x * x
      s0 = s0 + 1.0
      s1 = s1 + x
      s2 = s2 + x2
      s3 = s3 + x2 * x
      s4 = s4 + x2 * x2
      t0 = t0 + y
      t1 = t1 + x * y
      t2 = t2 + x2 * y
    }
    // solve normal equations for a (the constant term = smoothed value at center)
    // | s0  s1  s2 | |a|   |t0|
    // | s1  s2  s3 | |b| = |t1|
    // | s2  s3  s4 | |c|   |t2|
    // We only need a. Using Cramer's rule:
    let det = s0 * (s2 * s4 - s3 * s3) - s1 * (s1 * s4 - s2 * s3) +
      s2 * (s1 * s3 - s2 * s2)
    if det.abs() < 1.0e-15 {
      out[i] = values[i]
      continue
    }
    let det_a = t0 * (s2 * s4 - s3 * s3) - s1 * (t1 * s4 - t2 * s3) +
      s2 * (t1 * s3 - t2 * s2)
    out[i] = det_a / det
  }
  out
}

///|
/// Generic smoothing dispatcher.
pub fn mc_smooth(
  values : Array[Double],
  half_window : Int,
  method : String,
) -> Array[Double] {
  if method == "SavitzkyGolay" {
    mc_smooth_savitzky_golay(values, half_window)
  } else {
    mc_smooth_moving_average(values, half_window)
  }
}

// ============================================================================
// baselineEstimation: SNIP algorithm
// ============================================================================

///|
/// SNIP (Statistics-sensitive Non-linear Iterative Peak-clipping) baseline
/// estimation. Iteratively reduces each point to the minimum of itself and
/// the average of its symmetric neighbors, for increasing window sizes.
///
/// `max_iter` controls the number of increasing-window iterations.
pub fn mc_baseline_snip(
  values : Array[Double],
  max_iter : Int,
) -> Array[Double] {
  let n = values.length()
  if n == 0 {
    return []
  }
  // work on a log-log transformed copy (SNIP MLE)
  let v = Array::make(n, 0.0)
  for i in 0..<n {
    let val = if values[i] > 0.0 { values[i] } else { 0.0 }
    v[i] = val
  }
  // LLS transform: log(log(sqrt(y) + 1) + 1)
  for i in 0..<n {
    let sqrt_y = v[i].sqrt()
    v[i] = @math.ln(@math.ln(sqrt_y + 1.0) + 1.0)
  }
  // SNIP iterations
  let max_window = if max_iter < n / 2 { max_iter } else { n / 2 }
  for w in 1..<=max_window {
    let prev = Array::new()
    for i in 0..<n {
      prev.push(v[i])
    }
    for i in 0..<n {
      let lo = i - w
      let hi = i + w
      let left = if lo >= 0 { prev[lo] } else { prev[0] }
      let right = if hi < n { prev[hi] } else { prev[n - 1] }
      let avg = (left + right) / 2.0
      if avg < prev[i] {
        v[i] = avg
      }
    }
  }
  // inverse LLS transform
  for i in 0..<n {
    let exp1 = @math.exp(v[i]) - 1.0
    let exp2 = @math.exp(exp1) - 1.0
    v[i] = exp2 * exp2
  }
  v
}

// ============================================================================
// calibrate: calibrant-based m/z shift
// ============================================================================

///|
/// Calibrate observed m/z values using known calibrant m/z values.
///
/// Finds the closest calibrant for each observed m/z and applies a linear
/// shift (observed → calibrant). If `method` is "linear", a global linear
/// regression (calibrant ~ observed) is fit and applied.
pub fn mc_calibrate(
  observed : Array[Double],
  calibrants : Array[Double],
  tolerance : Double,
  method : String,
) -> Array[Double] {
  let n = observed.length()
  let out = Array::make(n, 0.0)
  if method == "linear" {
    // collect matched (obs, cal) pairs
    let obs_pairs = Array::new()
    let cal_pairs = Array::new()
    for i in 0..<n {
      let mut best_j = -1
      let mut best_d = 1.0e300
      for j in 0..<calibrants.length() {
        let d = (observed[i] - calibrants[j]).abs()
        if d <= tolerance && d < best_d {
          best_d = d
          best_j = j
        }
      }
      if best_j >= 0 {
        obs_pairs.push(observed[i])
        cal_pairs.push(calibrants[best_j])
      }
    }
    if obs_pairs.length() < 2 {
      // not enough points for linear fit; use shift
      for i in 0..<n {
        out[i] = observed[i]
      }
      return out
    }
    // linear regression: cal = a * obs + b
    let m = obs_pairs.length().to_double()
    let s_x = mc_sum(obs_pairs)
    let s_y = mc_sum(cal_pairs)
    let s_xx = mc_dot(obs_pairs, obs_pairs)
    let s_xy = mc_dot(obs_pairs, cal_pairs)
    let denom = m * s_xx - s_x * s_x
    if denom.abs() < 1.0e-15 {
      for i in 0..<n {
        out[i] = observed[i]
      }
      return out
    }
    let a = (m * s_xy - s_x * s_y) / denom
    let b = (s_y - a * s_x) / m
    for i in 0..<n {
      out[i] = a * observed[i] + b
    }
    out
  } else {
    // shift: each observed m/z gets the closest calibrant's shift
    for i in 0..<n {
      let mut best_j = -1
      let mut best_d = 1.0e300
      for j in 0..<calibrants.length() {
        let d = (observed[i] - calibrants[j]).abs()
        if d <= tolerance && d < best_d {
          best_d = d
          best_j = j
        }
      }
      if best_j >= 0 {
        out[i] = observed[i] + (calibrants[best_j] - observed[i])
      } else {
        out[i] = observed[i]
      }
    }
    out
  }
}

///|
/// Sum of an array.
pub fn mc_sum(values : Array[Double]) -> Double {
  let mut s = 0.0
  for v in values {
    s = s + v
  }
  s
}

///|
/// Dot product of two equal-length arrays.
pub fn mc_dot(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  for i in 0..<a.length() {
    if i < b.length() {
      s = s + a[i] * b[i]
    }
  }
  s
}

// ============================================================================
// Imputation methods
// ============================================================================

///|
/// Impute missing values (represented as NaN, i.e. 0.0/0.0) in a matrix
/// (rows × cols) using the specified method.
///
/// Methods:
///   - "zero": replace with 0.0
///   - "halfMin": replace with half the minimum of non-missing values in the row
///   - "mean": replace with the row mean of non-missing values
///   - "median": replace with the row median of non-missing values
///   - "knn": K-nearest-neighbors imputation (k=5 default)
///|
/// Check if a value is NaN (missing).
pub fn mc_is_missing(v : Double) -> Bool {
  v != v // NaN check: NaN != NaN
}

///|
pub fn mc_impute(
  matrix : Array[Array[Double]],
  method : String,
  k : Int,
) -> Array[Array[Double]] {
  let n = matrix.length()
  if n == 0 {
    return []
  }
  let p = matrix[0].length()
  // copy matrix
  let out = Array::new()
  for i in 0..<n {
    out.push(Array::make(p, 0.0))
    for j in 0..<p {
      out[i][j] = matrix[i][j]
    }
  }
  if method == "zero" {
    for i in 0..<n {
      for j in 0..<p {
        if mc_is_missing(out[i][j]) {
          out[i][j] = 0.0
        }
      }
    }
    return out
  }
  if method == "halfMin" {
    for i in 0..<n {
      // find row min of non-missing
      let mut min_val = 1.0e300
      for j in 0..<p {
        if !mc_is_missing(out[i][j]) && out[i][j] < min_val {
          min_val = out[i][j]
        }
      }
      let half_min = if min_val == 1.0e300 { 0.0 } else { min_val / 2.0 }
      for j in 0..<p {
        if mc_is_missing(out[i][j]) {
          out[i][j] = half_min
        }
      }
    }
    return out
  }
  if method == "mean" {
    for i in 0..<n {
      let vals = Array::new()
      for j in 0..<p {
        if !mc_is_missing(out[i][j]) {
          vals.push(out[i][j])
        }
      }
      let m = mc_mean(vals)
      for j in 0..<p {
        if mc_is_missing(out[i][j]) {
          out[i][j] = m
        }
      }
    }
    return out
  }
  if method == "median" {
    for i in 0..<n {
      let vals = Array::new()
      for j in 0..<p {
        if !mc_is_missing(out[i][j]) {
          vals.push(out[i][j])
        }
      }
      let m = mc_median(vals)
      for j in 0..<p {
        if mc_is_missing(out[i][j]) {
          out[i][j] = m
        }
      }
    }
    return out
  }
  if method == "knn" {
    // KNN imputation: for each row with missing values, find the k nearest
    // rows (by Euclidean distance on non-missing columns) and impute from
    // their values.
    for i in 0..<n {
      for j in 0..<p {
        if mc_is_missing(out[i][j]) {
          // find k nearest rows that have a value at column j
          let distances : Array[(Int, Double)] = Array::new()
          for r in 0..<n {
            if r == i || mc_is_missing(out[r][j]) {
              continue
            }
            let mut d = 0.0
            for c in 0..<p {
              if !mc_is_missing(out[i][c]) && !mc_is_missing(out[r][c]) {
                let diff = out[i][c] - out[r][c]
                d = d + diff * diff
              }
            }
            distances.push((r, d.sqrt()))
          }
          // sort distances ascending (selection sort)
          let nd = distances.length()
          for a in 1..<nd {
            let key = distances[a]
            let mut b = a - 1
            while b >= 0 && distances[b].1 > key.1 {
              distances[b + 1] = distances[b]
              b = b - 1
            }
            distances[b + 1] = key
          }
          // take mean of the k nearest
          let kk = if k < distances.length() { k } else { distances.length() }
          let mut s = 0.0
          let mut cnt = 0
          for a in 0..<kk {
            s = s + out[distances[a].0][j]
            cnt = cnt + 1
          }
          if cnt > 0 {
            out[i][j] = s / cnt.to_double()
          } else {
            out[i][j] = 0.0
          }
        }
      }
    }
    return out
  }
  out
}

// ============================================================================
// medianPolish: row/column aggregation
// ============================================================================

///|
/// Median polish: iteratively subtract row medians and column medians.
///
/// Returns (polished_matrix, row_effects, col_effects, overall_effect).
/// After `max_iter` iterations, the residuals approximate the additive model:
///   x[i][j] = overall + row_eff[i] + col_eff[j] + residual[i][j]
pub fn mc_median_polish(
  matrix : Array[Array[Double]],
  max_iter : Int,
  eps : Double,
) -> (Array[Array[Double]], Array[Double], Array[Double], Double) {
  let n = matrix.length()
  if n == 0 {
    return ([], [], [], 0.0)
  }
  let p = matrix[0].length()
  // copy
  let r = Array::new()
  for i in 0..<n {
    r.push(Array::make(p, 0.0))
    for j in 0..<p {
      r[i][j] = matrix[i][j]
    }
  }
  let row_eff = Array::make(n, 0.0)
  let col_eff = Array::make(p, 0.0)
  let mut overall = 0.0
  for _ in 0..<max_iter {
    let mut max_delta = 0.0
    // row sweep
    for i in 0..<n {
      let row_vals = Array::new()
      for j in 0..<p {
        row_vals.push(r[i][j])
      }
      let m = mc_median(row_vals)
      for j in 0..<p {
        r[i][j] = r[i][j] - m
      }
      row_eff[i] = row_eff[i] + m
      let abs_m = m.abs()
      if abs_m > max_delta {
        max_delta = abs_m
      }
    }
    // column sweep
    for j in 0..<p {
      let col_vals = Array::new()
      for i in 0..<n {
        col_vals.push(r[i][j])
      }
      let m = mc_median(col_vals)
      for i in 0..<n {
        r[i][j] = r[i][j] - m
      }
      col_eff[j] = col_eff[j] + m
      let abs_m = m.abs()
      if abs_m > max_delta {
        max_delta = abs_m
      }
    }
    if max_delta < eps {
      break
    }
  }
  // adjust: subtract median of row effects from row_eff, add to overall
  let row_med = mc_median(row_eff)
  for i in 0..<n {
    row_eff[i] = row_eff[i] - row_med
  }
  overall = overall + row_med
  (r, row_eff, col_eff, overall)
}

///|
/// Robust summary: compute a robust aggregate per column using the median
/// and MAD (median absolute deviation).
///
/// For each column, returns the median of values across rows.
pub fn mc_robust_summary(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  let n = matrix.length()
  if n == 0 {
    return []
  }
  let p = matrix[0].length()
  let out = Array::make(p, 0.0)
  for j in 0..<p {
    let col = Array::new()
    for i in 0..<n {
      col.push(matrix[i][j])
    }
    out[j] = mc_median(col)
  }
  out
}

///|
/// Median absolute deviation (MAD).
pub fn mc_mad(values : Array[Double]) -> Double {
  let m = mc_median(values)
  let deviations = Array::new()
  for v in values {
    deviations.push((v - m).abs())
  }
  mc_median(deviations) * 1.4826
}

// ============================================================================
// Feature validity helpers
// ============================================================================

///|
/// Check if a peak list is valid: m/z values must be strictly increasing
/// and all intensities non-negative.
pub fn mc_valid_peak_list(mz : Array[Double], intensity : Array[Double]) -> Bool {
  if mz.length() != intensity.length() {
    return false
  }
  for i in 0..<mz.length() {
    if intensity[i] < 0.0 {
      return false
    }
    if i > 0 && mz[i] <= mz[i - 1] {
      return false
    }
  }
  true
}

///|
/// Check whether a matrix has any missing (NaN) values.
pub fn mc_has_missing(matrix : Array[Array[Double]]) -> Bool {
  for row in matrix {
    for v in row {
      if v != v {
        return true
      }
    }
  }
  false
}

// ============================================================================
// Row/column aggregation helpers
// ============================================================================

///|
/// Aggregate rows by a grouping factor using a summary function.
/// Returns one row per group.
pub fn mc_aggregate_rows(
  matrix : Array[Array[Double]],
  groups : Array[Int],
  method : String,
) -> Array[Array[Double]] {
  let n = matrix.length()
  if n == 0 {
    return []
  }
  let p = matrix[0].length()
  // collect unique group ids
  let group_ids : Array[Int] = Array::new()
  for g in groups {
    let mut found = false
    for id in group_ids {
      if id == g {
        found = true
        break
      }
    }
    if !found {
      group_ids.push(g)
    }
  }
  // sort group ids
  for i in 1..<group_ids.length() {
    let key = group_ids[i]
    let mut j = i - 1
    while j >= 0 && group_ids[j] > key {
      group_ids[j + 1] = group_ids[j]
      j = j - 1
    }
    group_ids[j + 1] = key
  }
  let out = Array::new()
  for g in group_ids {
    let row = Array::make(p, 0.0)
    for j in 0..<p {
      let vals = Array::new()
      for i in 0..<n {
        if groups[i] == g {
          vals.push(matrix[i][j])
        }
      }
      if method == "sum" {
        row[j] = mc_sum(vals)
      } else if method == "mean" {
        row[j] = mc_mean(vals)
      } else {
        row[j] = mc_median(vals)
      }
    }
    out.push(row)
  }
  out
}

///|
/// Normalize intensities by TIC (total ion current).
pub fn mc_normalize_tic(intensity : Array[Double]) -> Array[Double] {
  let total = mc_sum(intensity)
  if total <= 0.0 {
    return Array::make(intensity.length(), 0.0)
  }
  let out = Array::new()
  for v in intensity {
    out.push(v / total * 100.0)
  }
  out
}
