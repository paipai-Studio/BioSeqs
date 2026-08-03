///|
/// Bio.Statistics.lowess module - LOWESS smoothing algorithm.
/// Inspired by Biopython's Bio.Statistics.lowess module.
/// Implements Cleveland's (1979) locally weighted regression for
/// robust scatterplot smoothing.

///|
/// Tricube weight function: (1 - u^3)^3 for u in [0, 1], 0 otherwise.
/// Used as the kernel weight in LOWESS.
pub fn lowess_tricube_weight(u : Double) -> Double {
  if u < 0.0 {
    return 0.0
  }
  if u >= 1.0 {
    return 0.0
  }
  let tmp = 1.0 - u * u * u
  tmp * tmp * tmp
}

///|
/// Bisquare weight function: (1 - u^2)^2 for u in [0, 1], 0 otherwise.
/// Used for robustness weights in iterative LOWESS.
pub fn lowess_bisquare_weight(u : Double) -> Double {
  if u < 0.0 {
    return 0.0
  }
  if u >= 1.0 {
    return 0.0
  }
  let tmp = 1.0 - u * u
  tmp * tmp
}

///|
/// Calculate the median of an array (used for robustness scale estimate).
pub fn lowess_median(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let sorted = arr.copy()
  sorted.sort()
  let n = sorted.length()
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

///|
/// Weighted linear regression. Returns (intercept, slope).
/// Fits y = intercept + slope * x minimizing weighted sum of squared residuals.
pub fn lowess_weighted_linear_regression(
  x : Array[Double],
  y : Array[Double],
  w : Array[Double],
) -> (Double, Double) {
  let n = x.length()
  if n == 0 {
    return (0.0, 0.0)
  }
  if n == 1 {
    return (y[0], 0.0)
  }
  let mut sum_w = 0.0
  let mut sum_wx = 0.0
  let mut sum_wy = 0.0
  let mut sum_wxx = 0.0
  let mut sum_wxy = 0.0
  let mut i = 0
  while i < n {
    let wi = w[i]
    let xi = x[i]
    let yi = y[i]
    sum_w = sum_w + wi
    sum_wx = sum_wx + wi * xi
    sum_wy = sum_wy + wi * yi
    sum_wxx = sum_wxx + wi * xi * xi
    sum_wxy = sum_wxy + wi * xi * yi
    i = i + 1
  }
  if sum_w <= 0.0 {
    return (0.0, 0.0)
  }
  let mean_x = sum_wx / sum_w
  let mean_y = sum_wy / sum_w
  let var_x = sum_wxx / sum_w - mean_x * mean_x
  let cov_xy = sum_wxy / sum_w - mean_x * mean_y
  let slope = if var_x > 1.0e-30 { cov_xy / var_x } else { 0.0 }
  let intercept = mean_y - slope * mean_x
  (intercept, slope)
}

///|
/// Find the window [left, right] containing the r nearest neighbors of x[i].
/// Assumes x is sorted in ascending order.
fn lowess_find_window(
  x : Array[Double],
  n : Int,
  i : Int,
  r : Int,
) -> (Int, Int) {
  let mut left = i
  let mut right = i
  let mut count = 1
  while count < r {
    if left > 0 && right < n - 1 {
      if x[i] - x[left - 1] <= x[right + 1] - x[i] {
        left = left - 1
      } else {
        right = right + 1
      }
    } else if left > 0 {
      left = left - 1
    } else if right < n - 1 {
      right = right + 1
    } else {
      break
    }
    count = count + 1
  }
  (left, right)
}

///|
/// Compute the LOWESS smoothed value at a single point i.
/// Uses tricube kernel weights multiplied by robustness weights.
fn lowess_compute_point(
  x : Array[Double],
  y : Array[Double],
  n : Int,
  i : Int,
  r : Int,
  rw : Array[Double],
) -> Double {
  let (left, right) = lowess_find_window(x, n, i, r)
  // Bandwidth: max distance from x[i] to window boundaries
  let dist_left = x[i] - x[left]
  let dist_right = x[right] - x[i]
  let h_i = if dist_left > dist_right { dist_left } else { dist_right }
  // Weighted linear regression with tricube * robustness weights
  let mut sum_w = 0.0
  let mut sum_wx = 0.0
  let mut sum_wy = 0.0
  let mut sum_wxx = 0.0
  let mut sum_wxy = 0.0
  let mut j = left
  while j <= right {
    let dist = (x[j] - x[i]).abs()
    let u = if h_i > 0.0 { dist / h_i } else { 0.0 }
    let w = lowess_tricube_weight(u) * rw[j]
    let xj = x[j]
    let yj = y[j]
    sum_w = sum_w + w
    sum_wx = sum_wx + w * xj
    sum_wy = sum_wy + w * yj
    sum_wxx = sum_wxx + w * xj * xj
    sum_wxy = sum_wxy + w * xj * yj
    j = j + 1
  }
  if sum_w <= 0.0 {
    return y[i]
  }
  let mean_x = sum_wx / sum_w
  let mean_y = sum_wy / sum_w
  let var_x = sum_wxx / sum_w - mean_x * mean_x
  let cov_xy = sum_wxy / sum_w - mean_x * mean_y
  let slope = if var_x > 1.0e-30 { cov_xy / var_x } else { 0.0 }
  let intercept = mean_y - slope * mean_x
  intercept + slope * x[i]
}

///|
/// LOWESS smoothing algorithm (Cleveland 1979).
///
/// Parameters:
/// - x: independent variable array (must be sorted in ascending order)
/// - y: dependent variable array
/// - f: smoothing window fraction (0-1), controls the proportion of data used
/// - n_steps: number of robustness iterations (0 = no robustness iteration)
/// - delta: performance optimization parameter; when delta > 0, skip nearby points
///
/// Returns the smoothed y values.
pub fn lowess(
  x : Array[Double],
  y : Array[Double],
  f~ : Double = 0.6667,
  n_steps~ : Int = 0,
  delta~ : Double = 0.0,
) -> Array[Double] {
  let n = x.length()
  if n == 0 {
    return Array::new()
  }
  if n == 1 {
    return [y[0]]
  }
  // Window size: number of nearest neighbors to use
  let mut r = (f * n.to_double()).to_int()
  if r < 2 {
    r = 2
  }
  if r > n {
    r = n
  }
  let ys : Array[Double] = Array::make(n, 0.0)
  // Robustness weights (initialized to 1.0, updated each iteration)
  let rw : Array[Double] = Array::make(n, 1.0)
  let mut step = 0
  while step <= n_steps {
    // Smoothing pass with optional delta optimization
    if delta > 0.0 {
      // Delta optimization: compute at anchor points spaced ~delta apart,
      // linearly interpolate the points in between.
      ys[0] = lowess_compute_point(x, y, n, 0, r, rw)
      let mut prev_anchor = 0
      while prev_anchor < n - 1 {
        // Find next anchor: first point where x[j] - x[prev] >= delta
        let mut next_anchor = prev_anchor + 1
        while next_anchor < n - 1 && x[next_anchor] - x[prev_anchor] < delta {
          next_anchor = next_anchor + 1
        }
        ys[next_anchor] = lowess_compute_point(x, y, n, next_anchor, r, rw)
        // Linearly interpolate between prev_anchor and next_anchor
        let mut k = prev_anchor + 1
        while k < next_anchor {
          let denom = x[next_anchor] - x[prev_anchor]
          let t = if denom > 0.0 {
            (x[k] - x[prev_anchor]) / denom
          } else {
            0.0
          }
          ys[k] = ys[prev_anchor] + t * (ys[next_anchor] - ys[prev_anchor])
          k = k + 1
        }
        prev_anchor = next_anchor
      }
    } else {
      // No delta optimization: compute every point directly
      let mut i = 0
      while i < n {
        ys[i] = lowess_compute_point(x, y, n, i, r, rw)
        i = i + 1
      }
    }
    // Compute robustness weights for the next iteration
    if step < n_steps {
      let residuals : Array[Double] = Array::new()
      let mut k = 0
      while k < n {
        residuals.push((y[k] - ys[k]).abs())
        k = k + 1
      }
      let s = lowess_median(residuals)
      if s > 0.0 {
        let mut k = 0
        while k < n {
          let u = (y[k] - ys[k]).abs() / (6.0 * s)
          rw[k] = lowess_bisquare_weight(u)
          k = k + 1
        }
      }
    }
    step = step + 1
  }
  ys
}

///|
/// Robust LOWESS with default 2 robustness iterations.
/// Less sensitive to outliers than standard LOWESS.
pub fn lowess_robust(
  x : Array[Double],
  y : Array[Double],
  f~ : Double = 0.6667,
  n_steps~ : Int = 2,
  delta~ : Double = 0.0,
) -> Array[Double] {
  lowess(x, y, f=f, n_steps=n_steps, delta=delta)
}

///|
/// Predict y for a new x value based on existing LOWESS smoothed results.
/// Inside the data range: linear interpolation of smoothed_y.
/// Outside the data range: local LOWESS using the original y data and
/// window fraction f for more robust extrapolation.
pub fn lowess_predict(
  x_new : Double,
  x : Array[Double],
  y : Array[Double],
  smoothed_y : Array[Double],
  f~ : Double = 0.6667,
) -> Double {
  let n = x.length()
  if n == 0 {
    return 0.0
  }
  if n == 1 {
    return smoothed_y[0]
  }
  // Inside the data range: interpolate the smoothed values
  if x_new > x[0] && x_new < x[n - 1] {
    let mut i = 0
    while i < n - 1 && x[i + 1] < x_new {
      i = i + 1
    }
    let denom = x[i + 1] - x[i]
    let t = if denom > 0.0 { (x_new - x[i]) / denom } else { 0.0 }
    return smoothed_y[i] + t * (smoothed_y[i + 1] - smoothed_y[i])
  }
  // At or beyond boundaries: clamp to nearest smoothed value
  if x_new == x[0] {
    return smoothed_y[0]
  }
  if x_new == x[n - 1] {
    return smoothed_y[n - 1]
  }
  // Outside the data range: local LOWESS extrapolation using y and f
  let mut r = (f * n.to_double()).to_int()
  if r < 2 {
    r = 2
  }
  if r > n {
    r = n
  }
  // Find the r nearest neighbors of x_new (x is sorted ascending)
  let mut center = 0
  while center < n && x[center] < x_new {
    center = center + 1
  }
  if center >= n {
    center = n - 1
  }
  let mut left = center
  let mut right = center
  let mut count = 1
  while count < r {
    if left > 0 && right < n - 1 {
      if x_new - x[left - 1] <= x[right + 1] - x_new {
        left = left - 1
      } else {
        right = right + 1
      }
    } else if left > 0 {
      left = left - 1
    } else if right < n - 1 {
      right = right + 1
    } else {
      break
    }
    count = count + 1
  }
  let dist_left = (x_new - x[left]).abs()
  let dist_right = (x[right] - x_new).abs()
  let h = if dist_left > dist_right { dist_left } else { dist_right }
  // Weighted linear regression with tricube weights centered at x_new
  let mut sum_w = 0.0
  let mut sum_wx = 0.0
  let mut sum_wy = 0.0
  let mut sum_wxx = 0.0
  let mut sum_wxy = 0.0
  let mut j = left
  while j <= right {
    let dist = (x[j] - x_new).abs()
    let u = if h > 0.0 { dist / h } else { 0.0 }
    let w = lowess_tricube_weight(u)
    let xj = x[j]
    let yj = y[j]
    sum_w = sum_w + w
    sum_wx = sum_wx + w * xj
    sum_wy = sum_wy + w * yj
    sum_wxx = sum_wxx + w * xj * xj
    sum_wxy = sum_wxy + w * xj * yj
    j = j + 1
  }
  if sum_w <= 0.0 {
    if x_new < x[0] {
      return smoothed_y[0]
    }
    return smoothed_y[n - 1]
  }
  let mean_x = sum_wx / sum_w
  let mean_y = sum_wy / sum_w
  let var_x = sum_wxx / sum_w - mean_x * mean_x
  let cov_xy = sum_wxy / sum_w - mean_x * mean_y
  let slope = if var_x > 1.0e-30 { cov_xy / var_x } else { 0.0 }
  let intercept = mean_y - slope * mean_x
  intercept + slope * x_new
}

///|
/// Smooth a time series / signal using LOWESS.
/// Convenience function where x = 0, 1, 2, ... (indices as x values).
pub fn lowess_smooth_series(
  values : Array[Double],
  window_frac~ : Double = 0.3,
) -> Array[Double] {
  let n = values.length()
  if n == 0 {
    return Array::new()
  }
  let x : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    x[i] = i.to_double()
    i = i + 1
  }
  lowess(x, values, f=window_frac, n_steps=0, delta=0.0)
}
