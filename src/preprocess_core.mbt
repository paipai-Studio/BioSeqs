///|
/// preprocessCore - Bioconductor preprocessing core module.
/// Provides fundamental normalization and preprocessing algorithms
/// for microarray and RNA-seq data.
///
/// Key algorithms:
/// - Quantile normalization (normalize_quantiles)
/// - Quantile normalization with target distribution (normalize_quantiles_use_target)
/// - Quantile target determination (normalize_quantiles_determine_target)
/// - Invariant set normalization (normalize_invariantset)
/// - Linear interpolation for data alignment
/// - Subset quantile normalization

// ============================================================================
// Data structures
// ============================================================================

///|
/// Configuration for quantile normalization.
pub struct QuantileNormConfig {
  copy : Bool
  keep_names : Bool
  method : String
}

///|
/// Create default quantile normalization config.
pub fn QuantileNormConfig::new() -> QuantileNormConfig {
  QuantileNormConfig::{
    copy: true,
    keep_names: true,
    method: "quantile",
  }
}

///|
/// Result from invariant set normalization.
pub struct InvariantSetResult {
  normalized_data : Array[Array[Double]]
  invariant_set_indices : Array[Int]
  reference_index : Int
  fit_coefficients : Array[Double] // [intercept, slope]
}

///|
/// Create invariant set result.
pub fn InvariantSetResult::new(
  normalized_data : Array[Array[Double]],
  invariant_set_indices : Array[Int],
  reference_index : Int,
  fit_coefficients : Array[Double],
) -> InvariantSetResult {
  InvariantSetResult::{
    normalized_data,
    invariant_set_indices,
    reference_index,
    fit_coefficients,
  }
}

// ============================================================================
// Helper: Basic math utilities
// ============================================================================

///|
/// Absolute value for double.
fn pc_abs(x : Double) -> Double {
  if x < 0.0 { -x } else { x }
}

///|
/// Sort array of doubles in ascending order (bubble sort for simplicity).
fn pc_sort_doubles(arr : Array[Double]) -> Array[Double] {
  let n = arr.length()
  let result = arr.copy()
  for i in 0..<n {
    let start = i + 1
    for j in start..<n {
      if result[i] > result[j] {
        let tmp = result[i]
        result[i] = result[j]
        result[j] = tmp
      }
    }
  }
  result
}

///|
/// Get sorted indices (argsort) - returns indices that would sort array.
fn pc_argsort(arr : Array[Double]) -> Array[Int] {
  let n = arr.length()
  let indices : Array[Int] = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  // Sort indices by value
  for i in 0..<n {
    let start = i + 1
    for j in start..<n {
      if arr[indices[i]] > arr[indices[j]] {
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
    }
  }
  indices
}

///|
/// Linear interpolation: given sorted xs, sorted ys, and query x, return interpolated y.
pub fn pc_interp_linear(
  xs : Array[Double],
  ys : Array[Double],
  x : Double,
) -> Double {
  let n = xs.length()
  if n == 0 { return 0.0 }
  if n == 1 { return ys[0] }
  // Handle extrapolation below minimum
  if x <= xs[0] {
    return ys[0]
  }
  // Handle extrapolation above maximum
  if x >= xs[n - 1] {
    return ys[n - 1]
  }
  // Binary search for the right interval
  let mut lo = 0
  let mut hi = n - 1
  while lo < hi - 1 {
    let mid = (lo + hi) / 2
    if xs[mid] <= x {
      lo = mid
    } else {
      hi = mid
    }
  }
  // Linear interpolation between xs[lo] and xs[hi]
  let t = (x - xs[lo]) / (xs[hi] - xs[lo])
  ys[lo] + t * (ys[hi] - ys[lo])
}

///|
/// Linear model fit via least squares (simple y = a + bx).
/// Returns [intercept, slope].
fn pc_linear_fit(
  xs : Array[Double],
  ys : Array[Double],
) -> Array[Double] {
  let n = xs.length()
  if n < 2 { return [0.0, 0.0] }
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_xx = 0.0
  let mut sum_xy = 0.0
  for i in 0..<n {
    let x = xs[i]
    let y = ys[i]
    sum_x = sum_x + x
    sum_y = sum_y + y
    sum_xx = sum_xx + x * x
    sum_xy = sum_xy + x * y
  }
  let denom = n.to_double() * sum_xx - sum_x * sum_x
  if pc_abs(denom) < 0.0000000001 {
    return [sum_y / n.to_double(), 0.0]
  }
  let slope = (n.to_double() * sum_xy - sum_x * sum_y) / denom
  let intercept = (sum_y - slope * sum_x) / n.to_double()
  [intercept, slope]
}

// ============================================================================
// Quantile Normalization (core algorithm)
// ============================================================================

///|
/// Core quantile normalization algorithm.
/// For each row (feature) across columns (samples):
/// 1. Sort each column independently
/// 2. Take the mean of each row in the sorted matrix
/// 3. Rearrange each column to match the original order, but using the mean values
///
/// This makes the distributions of all samples identical.
pub fn normalize_quantiles(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  if n_rows == 0 { return matrix }
  let n_cols = matrix[0].length()
  if n_cols == 0 { return matrix }

  // Transpose: work with columns as arrays
  // column_arrays[col] = [row0, row1, ..., rowN]
  let column_arrays : Array[Array[Double]] = []
  for _col in 0..<n_cols {
    column_arrays.push([])
  }
  for col in 0..<n_cols {
    let col_data : Array[Double] = Array::make(n_rows, 0.0)
    for row in 0..<n_rows {
      col_data[row] = matrix[row][col]
    }
    column_arrays[col] = col_data
  }

  // Step 1: Sort each column, keep track of original indices
  let sorted_cols : Array[Array[Double]] = []
  for _col in 0..<n_cols {
    sorted_cols.push([])
  }
  let order_indices : Array[Array[Int]] = []
  for _col in 0..<n_cols {
    order_indices.push([])
  }
  for col in 0..<n_cols {
    let idx = pc_argsort(column_arrays[col])
    let sorted : Array[Double] = Array::make(n_rows, 0.0)
    for r in 0..<n_rows {
      sorted[r] = column_arrays[col][idx[r]]
    }
    sorted_cols[col] = sorted
    order_indices[col] = idx
  }

  // Step 2: Compute mean per row across sorted columns (the target distribution)
  let means : Array[Double] = Array::make(n_rows, 0.0)
  for r in 0..<n_rows {
    let mut s = 0.0
    for col in 0..<n_cols {
      s = s + sorted_cols[col][r]
    }
    means[r] = s / n_cols.to_double()
  }

  // Step 3: Put mean values back into original positions
  let result : Array[Array[Double]] = []
  for _row in 0..<n_rows {
    result.push(Array::make(n_cols, 0.0))
  }
  for col in 0..<n_cols {
    for r in 0..<n_rows {
      let original_pos = order_indices[col][r]
      result[original_pos][col] = means[r]
    }
  }

  result
}

///|
/// Quantile normalization using a specified target distribution.
/// Each column is normalized to match the target distribution.
///
/// target should be a sorted array of length = n_rows.
pub fn normalize_quantiles_use_target(
  matrix : Array[Array[Double]],
  target : Array[Double],
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  if n_rows == 0 { return matrix }
  let n_cols = matrix[0].length()
  if n_cols == 0 { return matrix }
  if target.length() != n_rows { return matrix }

  let sorted_target = pc_sort_doubles(target)

  let result : Array[Array[Double]] = []
  for _row in 0..<n_rows {
    result.push(Array::make(n_cols, 0.0))
  }

  for col in 0..<n_cols {
    // Extract column
    let col_data : Array[Double] = Array::make(n_rows, 0.0)
    for row in 0..<n_rows {
      col_data[row] = matrix[row][col]
    }
    let order_idx = pc_argsort(col_data)
    // Place sorted target values back to original positions
    for r in 0..<n_rows {
      let original_pos = order_idx[r]
      result[original_pos][col] = sorted_target[r]
    }
  }

  result
}

///|
/// Determine target quantile distribution from a matrix.
/// Takes the mean of sorted values across columns.
pub fn normalize_quantiles_determine_target(
  matrix : Array[Array[Double]],
) -> Array[Double] {
  let n_rows = matrix.length()
  if n_rows == 0 { return [] }
  let n_cols = matrix[0].length()
  if n_cols == 0 { return [] }

  // Sort each column
  let sorted_cols : Array[Array[Double]] = []
  for _col in 0..<n_cols {
    sorted_cols.push([])
  }
  for col in 0..<n_cols {
    let col_data : Array[Double] = Array::make(n_rows, 0.0)
    for row in 0..<n_rows {
      col_data[row] = matrix[row][col]
    }
    sorted_cols[col] = pc_sort_doubles(col_data)
  }

  // Mean per row
  let means : Array[Double] = Array::make(n_rows, 0.0)
  for r in 0..<n_rows {
    let mut s = 0.0
    for col in 0..<n_cols {
      s = s + sorted_cols[col][r]
    }
    means[r] = s / n_cols.to_double()
  }
  means
}

///|
/// Subset quantile normalization (normalize just a subset and apply to all).
/// Useful for normalizing by group or using a subset of features.
pub fn normalize_quantiles_subset(
  matrix : Array[Array[Double]],
  subset_indices : Array[Int],
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  if n_rows == 0 { return matrix }
  let n_cols = matrix[0].length()
  if n_cols == 0 { return matrix }
  if subset_indices.length() == 0 { return normalize_quantiles(matrix) }

  // Build subset matrix
  let sub_n = subset_indices.length()
  let sub_matrix : Array[Array[Double]] = []
  for _s_idx in 0..<sub_n {
    sub_matrix.push([])
  }
  for s_idx in 0..<sub_n {
    let orig_row = subset_indices[s_idx]
    if orig_row < 0 || orig_row >= n_rows { continue }
    sub_matrix[s_idx] = matrix[orig_row].copy()
  }

  // Determine target from subset
  let target_distribution = normalize_quantiles_determine_target(sub_matrix)

  // Apply target to all rows of original matrix
  // Resize target if needed by interpolation
  if target_distribution.length() == n_rows {
    normalize_quantiles_use_target(matrix, target_distribution)
  } else {
    // Interpolate target to correct length
    let interp_target : Array[Double] = Array::make(n_rows, 0.0)
    let from_len = target_distribution.length()
    let xs : Array[Double] = []
    for i in 0..<from_len {
      xs.push(i.to_double())
    }
    for r in 0..<n_rows {
      let pct = r.to_double() / (n_rows - 1).to_double()
      let from_pos = pct * (from_len - 1).to_double()
      interp_target[r] = pc_interp_linear(
        xs,
        pc_sort_doubles(target_distribution),
        from_pos,
      )
    }
    normalize_quantiles_use_target(matrix, interp_target)
  }
}

// ============================================================================
// Invariant Set Normalization
// ============================================================================

///|
/// Find invariant set of features between reference and target array.
/// An invariant set is a set of features whose rank order is preserved
/// (roughly same ranks in both arrays). Used for normalization between arrays.
///
/// Returns indices of invariant features.
pub fn find_invariant_set(
  reference : Array[Double],
  target : Array[Double],
  threshold~ : Double = 0.005,
  max_iterations~ : Int = 10,
) -> Array[Int] {
  let n = reference.length()
  if n != target.length() || n == 0 { return [] }

  // Rank each array
  let ref_order = pc_argsort(reference)
  let tgt_order = pc_argsort(target)

  let ref_ranks : Array[Double] = Array::make(n, 0.0)
  let tgt_ranks : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    ref_ranks[ref_order[i]] = i.to_double()
    tgt_ranks[tgt_order[i]] = i.to_double()
  }

  // Start with all indices as invariant set
  let mut invariant : Array[Int] = []
  for i in 0..<n {
    invariant.push(i)
  }

  let mut iter = 0
  while iter < max_iterations {
    // Fit line on current invariant set
    let xs : Array[Double] = invariant.map(fn(i) { ref_ranks[i] })
    let ys : Array[Double] = invariant.map(fn(i) { tgt_ranks[i] })
    let coef = pc_linear_fit(xs, ys)
    let intercept = coef[0]
    let slope = coef[1]

    // Compute residuals and filter
    let new_invariant : Array[Int] = []
    let mut max_residual = 0.0
    for idx in invariant {
      let predicted = intercept + slope * ref_ranks[idx]
      let residual = pc_abs(tgt_ranks[idx] - predicted)
      if max_residual < residual { max_residual = residual }
    }
    let cutoff = threshold * n.to_double()
    if max_residual < cutoff || invariant.length() < 3 {
      break
    }
    let rel_cutoff = if cutoff > 0.0001 { cutoff } else { max_residual * 0.5 }
    for idx in invariant {
      let predicted = intercept + slope * ref_ranks[idx]
      let residual = pc_abs(tgt_ranks[idx] - predicted)
      if residual <= rel_cutoff {
        new_invariant.push(idx)
      }
    }
    if new_invariant.length() == invariant.length() {
      break
    }
    invariant = new_invariant
    if invariant.length() < 3 { break }
    iter = iter + 1
  }

  invariant
}

///|
/// Normalize all columns to a reference column using invariant set.
/// Reference column is specified by reference_index.
pub fn normalize_invariantset(
  matrix : Array[Array[Double]],
  reference_index~ : Int = 0,
  threshold~ : Double = 0.005,
) -> InvariantSetResult {
  let n_rows = matrix.length()
  if n_rows == 0 {
    return InvariantSetResult::new([], [], -1, [0.0, 0.0])
  }
  let n_cols = matrix[0].length()
  if reference_index < 0 || reference_index >= n_cols {
    return InvariantSetResult::new(matrix.copy(), [], 0, [0.0, 1.0])
  }

  // Extract reference column
  let ref_col : Array[Double] = Array::make(n_rows, 0.0)
  for row in 0..<n_rows {
    ref_col[row] = matrix[row][reference_index]
  }

  let result_matrix : Array[Array[Double]] = []
  for row in 0..<n_rows {
    let r : Array[Double] = Array::make(n_cols, 0.0)
    r[reference_index] = ref_col[row]
    result_matrix.push(r)
  }

  // Find invariant set using first non-reference column for global set
  let first_other = if reference_index == 0 { 1 } else { 0 }
  let first_col : Array[Double] = Array::make(n_rows, 0.0)
  for row in 0..<n_rows {
    first_col[row] = matrix[row][first_other]
  }
  let invariant_set = find_invariant_set(ref_col, first_col, threshold=threshold)

  // Get invariant set values for fitting
  let ref_inv_vals : Array[Double] = invariant_set.map(fn(i) { ref_col[i] })

  // Reference slope/intercept for reference column
  let ref_fit = pc_linear_fit(ref_inv_vals, ref_inv_vals)

  // Normalize each column
  for col in 0..<n_cols {
    if col == reference_index { continue }
    // Extract column
    let col_data : Array[Double] = Array::make(n_rows, 0.0)
    for row in 0..<n_rows {
      col_data[row] = matrix[row][col]
    }
    // Get invariant set values from this column
    let col_inv_vals : Array[Double] = invariant_set.map(fn(i) { col_data[i] })
    // Fit: col ~ ref (on invariant set)
    let fit = pc_linear_fit(col_inv_vals, ref_inv_vals)
    let intercept = fit[0]
    let slope = fit[1]
    // Apply: normalized = intercept + slope * col
    for row in 0..<n_rows {
      result_matrix[row][col] = intercept + slope * col_data[row]
    }
  }

  InvariantSetResult::new(
    result_matrix,
    invariant_set,
    reference_index,
    ref_fit,
  )
}

// ============================================================================
// Additional preprocessing utilities
// ============================================================================

///|
/// Log2 transformation with an offset to handle zeros.
/// log2(x + offset)
pub fn log2_transform(
  matrix : Array[Array[Double]],
  offset~ : Double = 1.0,
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  if n_rows == 0 { return matrix }
  let n_cols = matrix[0].length()
  let result : Array[Array[Double]] = []
  for _row in 0..<n_rows {
    result.push(Array::make(n_cols, 0.0))
  }
  for row in 0..<n_rows {
    for col in 0..<n_cols {
      let v = matrix[row][col] + offset
      if v <= 0.0 {
        result[row][col] = 0.0
      } else {
        result[row][col] = @math.log2(v)
      }
    }
  }
  result
}

///|
/// Background correction: subtract background estimate, floor at minimum.
/// Uses a simple percentile-based background estimate from each column.
pub fn background_correct_percentile(
  matrix : Array[Array[Double]],
  percentile~ : Double = 0.02,
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  if n_rows == 0 { return matrix }
  let n_cols = matrix[0].length()

  let result : Array[Array[Double]] = []
  for _row in 0..<n_rows {
    result.push(Array::make(n_cols, 0.0))
  }

  for col in 0..<n_cols {
    let col_data : Array[Double] = Array::make(n_rows, 0.0)
    for row in 0..<n_rows {
      col_data[row] = matrix[row][col]
    }
    let sorted = pc_sort_doubles(col_data)
    let idx = (percentile * (n_rows - 1).to_double()).to_int()
    let bg = sorted[if idx < 0 { 0 } else { if idx >= n_rows { n_rows - 1 } else { idx } }]
    for row in 0..<n_rows {
      let corrected = matrix[row][col] - bg
      result[row][col] = if corrected < 0.0 { 0.0 } else { corrected }
    }
  }
  result
}

///|
/// Median center each column (subtract column median).
pub fn median_center_columns(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  if n_rows == 0 { return matrix }
  let n_cols = matrix[0].length()

  let result : Array[Array[Double]] = []
  for _row in 0..<n_rows {
    let r : Array[Double] = Array::make(n_cols, 0.0)
    result.push(r)
  }

  for col in 0..<n_cols {
    let col_data : Array[Double] = Array::make(n_rows, 0.0)
    for row in 0..<n_rows {
      col_data[row] = matrix[row][col]
    }
    let sorted = pc_sort_doubles(col_data)
    let median = if n_rows % 2 == 0 {
      (sorted[n_rows / 2 - 1] + sorted[n_rows / 2]) / 2.0
    } else {
      sorted[n_rows / 2]
    }
    for row in 0..<n_rows {
      result[row][col] = matrix[row][col] - median
    }
  }
  result
}

///|
/// Compute column summary statistics.
/// Returns [mean, median, sd, min, max] for each column.
pub fn column_summary(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_rows = matrix.length()
  if n_rows == 0 { return [] }
  let n_cols = matrix[0].length()
  let result : Array[Array[Double]] = []
  for _col in 0..<n_cols {
    result.push([])
  }

  for col in 0..<n_cols {
    let col_data : Array[Double] = Array::make(n_rows, 0.0)
    let mut sum = 0.0
    for row in 0..<n_rows {
      let v = matrix[row][col]
      col_data[row] = v
      sum = sum + v
    }
    let mean = sum / n_rows.to_double()
    let sorted = pc_sort_doubles(col_data)
    let min_v = sorted[0]
    let max_v = sorted[n_rows - 1]
    let median = if n_rows % 2 == 0 {
      (sorted[n_rows / 2 - 1] + sorted[n_rows / 2]) / 2.0
    } else {
      sorted[n_rows / 2]
    }
    let mut sum_sq = 0.0
    for row in 0..<n_rows {
      let diff = matrix[row][col] - mean
      sum_sq = sum_sq + diff * diff
    }
    let sd = if n_rows > 1 {
      (sum_sq / (n_rows - 1).to_double()).sqrt()
    } else {
      0.0
    }
    result[col] = [mean, median, sd, min_v, max_v]
  }
  result
}
