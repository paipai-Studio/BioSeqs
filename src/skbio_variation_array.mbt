/// Variation array for compositional data.
///
/// scikit-bio equivalent: `composition.variation_array()`
///
/// The variation array stores the log-ratio of every pair of components.
/// VA[i][j] = ln(x_i / x_j)
/// It is a compact representation of compositional data and the basis
/// for many CoDA visualization techniques.

///|
/// Compute the variation array for a composition.
/// Returns an n×n matrix where entry [i][j] = ln(x_i / x_j).
pub fn skbio_variation_array(x : Array[Double]) -> Array[Array[Double]] {
  let n = x.length()
  let result : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    result[i] = Array::make(n, 0.0)
  }
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      if i != j && x[i] > 0.0 && x[j] > 0.0 {
        result[i][j] = @math.ln(x[i]) - @math.ln(x[j])
      }
    }
  }
  result
}

///|
/// Compute the variation array for multiple compositions (averaged).
/// Returns an n×n matrix of mean log-ratios across all samples.
pub fn skbio_variation_array_mean(
  compositions : Array[Array[Double]],
) -> Array[Array[Double]] {
  if compositions.length() == 0 {
    return []
  }
  let n = compositions[0].length()
  let sum_va : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    sum_va[i] = Array::make(n, 0.0)
  }
  for comp in compositions {
    let va = skbio_variation_array(comp)
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        sum_va[i][j] = sum_va[i][j] + va[i][j]
      }
    }
  }
  let count = compositions.length().to_double()
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      sum_va[i][j] = sum_va[i][j] / count
    }
  }
  sum_va
}

///|
/// Compute the variation array variance across multiple compositions.
/// Returns an n×n matrix of variances of log-ratios.
pub fn skbio_variation_array_variance(
  compositions : Array[Array[Double]],
) -> Array[Array[Double]] {
  if compositions.length() == 0 {
    return []
  }
  let n = compositions[0].length()
  let mean_va = skbio_variation_array_mean(compositions)
  let sum_sq : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    sum_sq[i] = Array::make(n, 0.0)
  }
  for comp in compositions {
    let va = skbio_variation_array(comp)
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        let d = va[i][j] - mean_va[i][j]
        sum_sq[i][j] = sum_sq[i][j] + d * d
      }
    }
  }
  let count = compositions.length().to_double()
  if count > 1.0 {
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < n; j = j + 1 {
        sum_sq[i][j] = sum_sq[i][j] / (count - 1.0)
      }
    }
  }
  sum_sq
}

///|
/// Extract the upper triangle of a variation array as a flat array.
/// Useful for plotting and clustering.
pub fn skbio_variation_array_upper(va : Array[Array[Double]]) -> Array[Double] {
  let n = va.length()
  let result : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      result.push(va[i][j])
    }
  }
  result
}

///|
/// Find the pair of components with maximum variation (max |log-ratio|).
/// Returns (index_i, index_j, value).
pub fn skbio_variation_array_max(x : Array[Double]) -> (Int32, Int32, Double) {
  let va = skbio_variation_array(x)
  let n = va.length()
  let mut best_i : Int32 = 0
  let mut best_j : Int32 = 0
  let mut best_val : Double = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if va[i][j].abs() > best_val {
        best_val = va[i][j].abs()
        best_i = i
        best_j = j
      }
    }
  }
  (best_i, best_j, best_val)
}
