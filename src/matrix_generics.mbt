///|
/// MatrixGenerics - Bioconductor row/col summary statistics.
///
/// Provides efficient row and column summary statistics for matrix-like
/// objects, mirroring the Bioconductor MatrixGenerics package.
///
/// Core functions:
///   - rowMeans / colMeans: compute row or column means
///   - rowSums / colSums: compute row or column sums
///   - rowVars / colVars: compute row or column variances
///   - rowSds / colSds: compute row or column standard deviations
///   - rowMedians / colMedians: compute row or column medians
///   - rowMins / colMins: compute row or column minima
///   - rowMaxs / colMaxs: compute row or column maxima
///   - rowRanges / colRanges: compute row or column ranges
///   - rowMad / colMad: median absolute deviation
///   - rowCounts: count non-NA values per row
///   - rowAnys / colAnys: any value in row/col meets condition
///   - rowAlls / colAlls: all values in row/col meet condition
///
/// All public identifiers are prefixed with `mg_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// A simple 2D matrix backed by a flat array (row-major order).
pub struct MgMatrix {
  n_rows : Int
  n_cols : Int
  data : Array[Double]
}

///|
/// Create a matrix from a 2D array (rows x cols).
pub fn MgMatrix::new(data : Array[Array[Double]]) -> MgMatrix {
  let n_rows = data.length()
  if n_rows == 0 {
    return { n_rows: 0, n_cols: 0, data: Array::make(0, 0.0) }
  }
  let n_cols = data[0].length()
  let flat : Array[Double] = []
  for i in 0..<n_rows {
    for j in 0..<n_cols {
      flat.push(data[i][j])
    }
  }
  { n_rows, n_cols, data: flat }
}

///|
/// Create a matrix filled with a value.
pub fn MgMatrix::zeros(n_rows : Int, n_cols : Int) -> MgMatrix {
  let data = Array::make(n_rows * n_cols, 0.0)
  { n_rows, n_cols, data }
}

///|
/// Get number of rows.
pub fn MgMatrix::dim_rows(self : MgMatrix) -> Int {
  self.n_rows
}

///|
/// Get number of columns.
pub fn MgMatrix::dim_cols(self : MgMatrix) -> Int {
  self.n_cols
}

///|
/// Get element at (i, j).
pub fn MgMatrix::get(self : MgMatrix, i : Int, j : Int) -> Double {
  self.data[i * self.n_cols + j]
}

///|
/// Set element at (i, j).
pub fn MgMatrix::set(self : MgMatrix, i : Int, j : Int, val : Double) -> Unit {
  self.data[i * self.n_cols + j] = val
}

///|
/// Get a row as an array.
pub fn MgMatrix::get_row(self : MgMatrix, i : Int) -> Array[Double] {
  let row : Array[Double] = []
  for j in 0..<self.n_cols {
    row.push(self.data[i * self.n_cols + j])
  }
  row
}

///|
/// Get a column as an array.
pub fn MgMatrix::get_col(self : MgMatrix, j : Int) -> Array[Double] {
  let col : Array[Double] = []
  for i in 0..<self.n_rows {
    col.push(self.data[i * self.n_cols + j])
  }
  col
}

// ============================================================================
// Row-level statistics
// ============================================================================

///|
/// Compute the mean of each row.
pub fn mg_row_means(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_rows
  let result : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let row = mat.get_row(i)
    result[i] = mg_mean(row)
  }
  result
}

///|
/// Compute the sum of each row.
pub fn mg_row_sums(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_rows
  let result : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let row = mat.get_row(i)
    let mut s = 0.0
    for j in 0..<row.length() {
      s = s + row[j]
    }
    result[i] = s
  }
  result
}

///|
/// Compute the variance of each row.
pub fn mg_row_vars(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_rows
  let result : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let row = mat.get_row(i)
    result[i] = mg_variance(row)
  }
  result
}

///|
/// Compute the standard deviation of each row.
pub fn mg_row_sds(mat : MgMatrix) -> Array[Double] {
  let vars = mg_row_vars(mat)
  let n = vars.length()
  let result : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    result[i] = vars[i].sqrt()
  }
  result
}

///|
/// Compute the median of each row.
pub fn mg_row_medians(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_rows
  let result : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let row = mat.get_row(i)
    result[i] = mg_median(row)
  }
  result
}

///|
/// Compute the minimum of each row.
pub fn mg_row_mins(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_rows
  let result : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let row = mat.get_row(i)
    result[i] = mg_min(row)
  }
  result
}

///|
/// Compute the maximum of each row.
pub fn mg_row_maxs(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_rows
  let result : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let row = mat.get_row(i)
    result[i] = mg_max(row)
  }
  result
}

///|
/// Compute the range (max - min) of each row.
pub fn mg_row_ranges(mat : MgMatrix) -> Array[Double] {
  let mins = mg_row_mins(mat)
  let maxs = mg_row_maxs(mat)
  let n = mins.length()
  let result : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    result[i] = maxs[i] - mins[i]
  }
  result
}

///|
/// Count non-zero values in each row.
pub fn mg_row_counts(mat : MgMatrix) -> Array[Int] {
  let n = mat.n_rows
  let result : Array[Int] = Array::make(n, 0)
  for i in 0..<n {
    let row = mat.get_row(i)
    let mut count = 0
    for j in 0..<row.length() {
      if row[j] != 0.0 {
        count = count + 1
      }
    }
    result[i] = count
  }
  result
}

///|
/// Compute the median absolute deviation for each row.
pub fn mg_row_mads(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_rows
  let result : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let row = mat.get_row(i)
    result[i] = mg_mad(row)
  }
  result
}

///|
/// Check if any value in each row exceeds a threshold.
pub fn mg_row_anys(mat : MgMatrix, threshold : Double) -> Array[Bool] {
  let n = mat.n_rows
  let result : Array[Bool] = Array::make(n, false)
  for i in 0..<n {
    let row = mat.get_row(i)
    let mut found = false
    for j in 0..<row.length() {
      if row[j] > threshold {
        found = true
        break
      }
    }
    result[i] = found
  }
  result
}

///|
/// Check if all values in each row exceed a threshold.
pub fn mg_row_alls(mat : MgMatrix, threshold : Double) -> Array[Bool] {
  let n = mat.n_rows
  let result : Array[Bool] = Array::make(n, false)
  for i in 0..<n {
    let row = mat.get_row(i)
    let mut all_above = true
    for j in 0..<row.length() {
      if row[j] <= threshold {
        all_above = false
        break
      }
    }
    result[i] = all_above
  }
  result
}

// ============================================================================
// Column-level statistics
// ============================================================================

///|
/// Compute the mean of each column.
pub fn mg_col_means(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_cols
  let result : Array[Double] = Array::make(n, 0.0)
  for j in 0..<n {
    let col = mat.get_col(j)
    result[j] = mg_mean(col)
  }
  result
}

///|
/// Compute the sum of each column.
pub fn mg_col_sums(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_cols
  let result : Array[Double] = Array::make(n, 0.0)
  for j in 0..<n {
    let col = mat.get_col(j)
    let mut s = 0.0
    for i in 0..<col.length() {
      s = s + col[i]
    }
    result[j] = s
  }
  result
}

///|
/// Compute the variance of each column.
pub fn mg_col_vars(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_cols
  let result : Array[Double] = Array::make(n, 0.0)
  for j in 0..<n {
    let col = mat.get_col(j)
    result[j] = mg_variance(col)
  }
  result
}

///|
/// Compute the standard deviation of each column.
pub fn mg_col_sds(mat : MgMatrix) -> Array[Double] {
  let vars = mg_col_vars(mat)
  let n = vars.length()
  let result : Array[Double] = Array::make(n, 0.0)
  for j in 0..<n {
    result[j] = vars[j].sqrt()
  }
  result
}

///|
/// Compute the median of each column.
pub fn mg_col_medians(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_cols
  let result : Array[Double] = Array::make(n, 0.0)
  for j in 0..<n {
    let col = mat.get_col(j)
    result[j] = mg_median(col)
  }
  result
}

///|
/// Compute the minimum of each column.
pub fn mg_col_mins(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_cols
  let result : Array[Double] = Array::make(n, 0.0)
  for j in 0..<n {
    let col = mat.get_col(j)
    result[j] = mg_min(col)
  }
  result
}

///|
/// Compute the maximum of each column.
pub fn mg_col_maxs(mat : MgMatrix) -> Array[Double] {
  let n = mat.n_cols
  let result : Array[Double] = Array::make(n, 0.0)
  for j in 0..<n {
    let col = mat.get_col(j)
    result[j] = mg_max(col)
  }
  result
}

///|
/// Compute the range of each column.
pub fn mg_col_ranges(mat : MgMatrix) -> Array[Double] {
  let mins = mg_col_mins(mat)
  let maxs = mg_col_maxs(mat)
  let n = mins.length()
  let result : Array[Double] = Array::make(n, 0.0)
  for j in 0..<n {
    result[j] = maxs[j] - mins[j]
  }
  result
}

///|
/// Check if any value in each column exceeds a threshold.
pub fn mg_col_anys(mat : MgMatrix, threshold : Double) -> Array[Bool] {
  let n = mat.n_cols
  let result : Array[Bool] = Array::make(n, false)
  for j in 0..<n {
    let col = mat.get_col(j)
    let mut found = false
    for i in 0..<col.length() {
      if col[i] > threshold {
        found = true
        break
      }
    }
    result[j] = found
  }
  result
}

///|
/// Check if all values in each column exceed a threshold.
pub fn mg_col_alls(mat : MgMatrix, threshold : Double) -> Array[Bool] {
  let n = mat.n_cols
  let result : Array[Bool] = Array::make(n, false)
  for j in 0..<n {
    let col = mat.get_col(j)
    let mut all_above = true
    for i in 0..<col.length() {
      if col[i] <= threshold {
        all_above = false
        break
      }
    }
    result[j] = all_above
  }
  result
}

// ============================================================================
// Block-level processing
// ============================================================================

///|
/// Process rows in blocks of a given size.
pub fn mg_block_apply_rows(
  mat : MgMatrix,
  block_size : Int,
  block_fn : (MgMatrix) -> Array[Double],
) -> Array[Double] {
  let n_rows = mat.n_rows
  let n_cols = mat.n_cols
  let result : Array[Double] = Array::make(n_rows, 0.0)
  let mut start = 0
  while start < n_rows {
    let end = if start + block_size < n_rows { start + block_size } else { n_rows }
    let block_data : Array[Array[Double]] = []
    let mut i = start
    while i < end {
      let row_data : Array[Double] = []
      let mut j = 0
      while j < n_cols {
        row_data.push(mat.get(i, j))
        j = j + 1
      }
      block_data.push(row_data)
      i = i + 1
    }
    let block = MgMatrix::new(block_data)
    let block_result = block_fn(block)
    let mut k = 0
    while k < block_result.length() {
      result[start + k] = block_result[k]
      k = k + 1
    }
    start = end
  }
  result
}

///|
/// Process columns in blocks of a given size.
pub fn mg_block_apply_cols(
  mat : MgMatrix,
  block_size : Int,
  block_fn : (MgMatrix) -> Array[Double],
) -> Array[Double] {
  let n_rows = mat.n_rows
  let n_cols = mat.n_cols
  let result : Array[Double] = Array::make(n_cols, 0.0)
  let mut start = 0
  while start < n_cols {
    let end = if start + block_size < n_cols { start + block_size } else { n_cols }
    let block_data : Array[Array[Double]] = []
    let mut i = 0
    while i < n_rows {
      let row_data : Array[Double] = []
      let mut j = start
      while j < end {
        row_data.push(mat.get(i, j))
        j = j + 1
      }
      block_data.push(row_data)
      i = i + 1
    }
    let block = MgMatrix::new(block_data)
    let block_result = block_fn(block)
    let mut k = 0
    while k < block_result.length() {
      result[start + k] = block_result[k]
      k = k + 1
    }
    start = end
  }
  result
}

// ============================================================================
// Helper functions
// ============================================================================

///|
/// Compute mean of an array.
pub fn mg_mean(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < data.length() {
    sum = sum + data[i]
    i = i + 1
  }
  sum / data.length().to_double()
}

///|
/// Compute variance of an array (sample variance).
pub fn mg_variance(data : Array[Double]) -> Double {
  let n = data.length()
  if n <= 1 {
    return 0.0
  }
  let mean = mg_mean(data)
  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    let diff = data[i] - mean
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  sum_sq / (n.to_double() - 1.0)
}

///|
/// Compute median of an array.
pub fn mg_median(data : Array[Double]) -> Double {
  let n = data.length()
  if n == 0 {
    return 0.0
  }
  let sorted = data.copy()
  sorted.sort()
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

///|
/// Compute minimum of an array.
pub fn mg_min(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let mut min_val = data[0]
  let mut i = 1
  while i < data.length() {
    if data[i] < min_val {
      min_val = data[i]
    }
    i = i + 1
  }
  min_val
}

///|
/// Compute maximum of an array.
pub fn mg_max(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let mut max_val = data[0]
  let mut i = 1
  while i < data.length() {
    if data[i] > max_val {
      max_val = data[i]
    }
    i = i + 1
  }
  max_val
}

///|
/// Compute median absolute deviation (MAD).
pub fn mg_mad(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let med = mg_median(data)
  let abs_devs : Array[Double] = []
  let mut i = 0
  while i < data.length() {
    abs_devs.push((data[i] - med).abs())
    i = i + 1
  }
  1.4826 * mg_median(abs_devs)
}
