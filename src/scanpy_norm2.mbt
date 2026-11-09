/// Scanpy second-round: Advanced normalization and scaling
/// Inspired by scanpy.pp: normalize_total (exclude_highly_expressed, target_sum),
/// scale (zero_center, max_value), log1p (base)

///|
/// Normalize total counts with exclude_highly_expressed option.
/// When exclude_highly_expressed=true, the target sum is computed from
/// genes that are NOT in the top 99th percentile.
pub fn scanpy_normalize_total_exclude(
  x : Array[Array[Double]],
  exclude_highly_expressed~ : Bool = false,
  target_sum~ : Double = 0.0,
  percentile~ : Double = 0.99,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=n_genes)
    // Compute total or exclude high genes
    let mut total = 0.0
    if exclude_highly_expressed {
      // Find threshold for top percentile genes
      let sorted_vals = Array(capacity=n_genes)
      for g in 0..<n_genes {
        sorted_vals.push(x[i][g])
      }
      sorted_vals.sort_by(fn(a, b) {
        if a < b {
          -1
        } else if a > b {
          1
        } else {
          0
        }
      })
      let thresh_idx = (percentile * n_genes.to_double()).floor().to_int()
      let thresh = if thresh_idx < n_genes && thresh_idx >= 0 {
        sorted_vals[thresh_idx]
      } else {
        0.0
      }
      for g in 0..<n_genes {
        if x[i][g] < thresh {
          total = total + x[i][g]
        }
      }
    } else {
      for g in 0..<n_genes {
        total = total + x[i][g]
      }
    }
    // Determine target
    let target = if target_sum > 0.0 {
      target_sum
    } else if total > 0.0 {
      total // use median-like estimate (here just use this cell's total)
    } else {
      1.0
    }
    // Normalize
    let scale = if total > 0.0 { target / total } else { 1.0 }
    for g in 0..<n_genes {
      row.push(x[i][g] * scale)
    }
    result.push(row)
  }
  result
}

///|
/// Estimate target_sum from median of total counts across cells.
pub fn scanpy_estimate_target_sum(x : Array[Array[Double]]) -> Double {
  let totals = Array(capacity=x.length())
  for i in 0..<x.length() {
    let mut total = 0.0
    for g in 0..<x[i].length() {
      total = total + x[i][g]
    }
    totals.push(total)
  }
  totals.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let n = totals.length()
  if n == 0 {
    return 0.0
  }
  if n % 2 == 1 {
    totals[n / 2]
  } else {
    (totals[n / 2 - 1] + totals[n / 2]) / 2.0
  }
}

///|
/// Scale data with zero_center and max_value options.
/// If zero_center=true (default), subtract mean and divide by std.
/// If zero_center=false, divide by std without centering.
/// If max_value > 0, clip values to [-max_value, max_value].
pub fn scanpy_scale_ext(
  x : Array[Array[Double]],
  zero_center~ : Bool = true,
  max_value~ : Double = 10.0,
  axis~ : Int32 = 0, // 0=genes (columns), 1=cells (rows)
) -> Array[Array[Double]] {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=n_genes)
    for g in 0..<n_genes {
      row.push(x[i][g])
    }
    result.push(row)
  }
  if axis == 0 {
    // Per-gene (column) scaling
    for g in 0..<n_genes {
      let mut mean = 0.0
      for i in 0..<n_cells {
        mean = mean + result[i][g]
      }
      mean = mean / n_cells.to_double()
      let mut variance = 0.0
      for i in 0..<n_cells {
        let diff = if zero_center { result[i][g] - mean } else { result[i][g] }
        variance = variance + diff * diff
      }
      let std = @math.pow(variance / n_cells.to_double(), 0.5)
      if std > 0.0 {
        for i in 0..<n_cells {
          if zero_center {
            result[i][g] = (result[i][g] - mean) / std
          } else {
            result[i][g] = result[i][g] / std
          }
          // Clip
          if max_value > 0.0 {
            if result[i][g] > max_value {
              result[i][g] = max_value
            } else if result[i][g] < -max_value {
              result[i][g] = -max_value
            }
          }
        }
      }
    }
  } else {
    // Per-cell (row) scaling
    for i in 0..<n_cells {
      let mut mean = 0.0
      for g in 0..<n_genes {
        mean = mean + result[i][g]
      }
      mean = mean / n_genes.to_double()
      let mut variance = 0.0
      for g in 0..<n_genes {
        let diff = if zero_center { result[i][g] - mean } else { result[i][g] }
        variance = variance + diff * diff
      }
      let std = @math.pow(variance / n_genes.to_double(), 0.5)
      if std > 0.0 {
        for g in 0..<n_genes {
          if zero_center {
            result[i][g] = (result[i][g] - mean) / std
          } else {
            result[i][g] = result[i][g] / std
          }
          if max_value > 0.0 {
            if result[i][g] > max_value {
              result[i][g] = max_value
            } else if result[i][g] < -max_value {
              result[i][g] = -max_value
            }
          }
        }
      }
    }
  }
  result
}

///|
/// Log1p with configurable base.
/// log_b(1 + x) = ln(1 + x) / ln(b)
pub fn scanpy_log1p_base(
  x : Array[Array[Double]],
  base~ : Double = 2.718281828459045, // natural log by default
) -> Array[Array[Double]] {
  let log_base = @math.ln(base)
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let row = Array(capacity=x[i].length())
    for g in 0..<x[i].length() {
      row.push(@math.ln(1.0 + x[i][g]) / log_base)
    }
    result.push(row)
  }
  result
}
