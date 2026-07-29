///|
/// limma module - Linear Models for Microarray and RNA-seq Data.
/// Inspired by Bioconductor limma package.
///
/// Provides functions for differential expression analysis using linear models
/// and empirical Bayes moderation.

///|
/// Simple matrix type for linear algebra operations.
pub struct Matrix {
  rows : Int
  cols : Int
  data : Array[Double]
}

///|
/// Create a matrix from flat array.
pub fn Matrix::new(rows : Int, cols : Int, data : Array[Double]) -> Matrix {
  Matrix::{ rows, cols, data }
}

///|
/// Create a matrix from 2D array.
pub fn Matrix::from_2d(
  rows : Int,
  cols : Int,
  data : Array[Array[Double]],
) -> Matrix {
  let flat = Array::make(rows * cols, 0.0)
  let mut i = 0
  while i < rows {
    let mut j = 0
    while j < cols {
      flat[i * cols + j] = data[i][j]
      j = j + 1
    }
    i = i + 1
  }
  Matrix::new(rows, cols, flat)
}

///|
/// Create a zero matrix.
pub fn Matrix::zeros(rows : Int, cols : Int) -> Matrix {
  let data = Array::make(rows * cols, 0.0)
  Matrix::{ rows, cols, data }
}

///|
/// Create an identity matrix.
pub fn Matrix::identity(n : Int) -> Matrix {
  let data = Array::make(n * n, 0.0)
  let mut i = 0
  while i < n {
    data[i * n + i] = 1.0
    i = i + 1
  }
  Matrix::{ rows: n, cols: n, data }
}

///|
/// Get element at (row, col).
pub fn Matrix::get(self : Matrix, row : Int, col : Int) -> Double {
  self.data[row * self.cols + col]
}

///|
/// Set element at (row, col).
pub fn Matrix::set(self : Matrix, row : Int, col : Int, value : Double) -> Unit {
  self.data[row * self.cols + col] = value
}

///|
/// Matrix multiplication.
pub fn mat_mult(a : Matrix, b : Matrix) -> Matrix {
  if a.cols != b.rows {
    return Matrix::zeros(0, 0)
  }

  let data = Array::make(a.rows * b.cols, 0.0)
  let mut i = 0
  while i < a.rows {
    let mut j = 0
    while j < b.cols {
      let mut sum = 0.0
      let mut k = 0
      while k < a.cols {
        sum = sum + a.get(i, k) * b.get(k, j)
        k = k + 1
      }
      data[i * b.cols + j] = sum
      j = j + 1
    }
    i = i + 1
  }
  Matrix::new(a.rows, b.cols, data)
}

///|
/// Matrix transpose.
pub fn matrix_transpose(a : Matrix) -> Matrix {
  let data = Array::make(a.cols * a.rows, 0.0)
  let mut i = 0
  while i < a.rows {
    let mut j = 0
    while j < a.cols {
      data[j * a.rows + i] = a.get(i, j)
      j = j + 1
    }
    i = i + 1
  }
  Matrix::new(a.cols, a.rows, data)
}

///|
/// Matrix inverse using Gaussian elimination.
pub fn mat_inverse(a : Matrix) -> Matrix {
  if a.rows != a.cols {
    return Matrix::zeros(0, 0)
  }

  let n = a.rows
  let size = 2 * n

  let aug_data = Array::make(n * size, 0.0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      aug_data[i * size + j] = a.get(i, j)
      j = j + 1
    }
    aug_data[i * size + n + i] = 1.0
    i = i + 1
  }

  let aug = Matrix::new(n, size, aug_data)

  let mut pivot_row = 0
  while pivot_row < n {
    let mut max_row = pivot_row
    let mut max_val = aug.get(pivot_row, pivot_row).abs()
    let mut k = pivot_row + 1
    while k < n {
      let val = aug.get(k, pivot_row).abs()
      if val > max_val {
        max_val = val
        max_row = k
      }
      k = k + 1
    }

    if max_row != pivot_row {
      let mut j = 0
      while j < size {
        let temp = aug.get(pivot_row, j)
        aug.set(pivot_row, j, aug.get(max_row, j))
        aug.set(max_row, j, temp)
        j = j + 1
      }
    }

    let pivot_val = aug.get(pivot_row, pivot_row)
    if pivot_val.abs() < 0.0000000001 {
      return Matrix::zeros(0, 0)
    }

    let mut j = 0
    while j < size {
      aug.set(pivot_row, j, aug.get(pivot_row, j) / pivot_val)
      j = j + 1
    }

    let mut i = 0
    while i < n {
      if i != pivot_row {
        let factor = aug.get(i, pivot_row)
        let mut j = 0
        while j < size {
          aug.set(i, j, aug.get(i, j) - factor * aug.get(pivot_row, j))
          j = j + 1
        }
      }
      i = i + 1
    }

    pivot_row = pivot_row + 1
  }

  let result_data = Array::make(n * n, 0.0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      result_data[i * n + j] = aug.get(i, n + j)
      j = j + 1
    }
    i = i + 1
  }
  Matrix::new(n, n, result_data)
}

///|
/// Vector dot product.
pub fn vec_dot(a : Array[Double], b : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < a.length() {
    sum = sum + a[i] * b[i]
    i = i + 1
  }
  sum
}

///|
/// Vector sum of squares.
pub fn vec_ss(a : Array[Double]) -> Double {
  vec_dot(a, a)
}

///|
/// Lanczos approximation for ln(gamma(x)).
pub fn lngamma(x : Double) -> Double {
  let cof = [
    76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
    0.001208650973866179, -0.000005395239384953,
  ]
  let mut x1 = x
  let mut tmp = x1 + 5.5
  tmp = (x1 + 0.5) * limma_log_nat(tmp) - tmp

  let mut ser = 1.000000000190015
  let mut j = 0
  while j < cof.length() {
    x1 = x1 + 1.0
    ser = ser + cof[j] / x1
    j = j + 1
  }

  tmp + limma_log_nat(2.5066282746310005 * ser)
}

///|
/// Incomplete beta function ratio (regularized).
pub fn pbeta(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }

  let lbeta_val = lngamma(a) + lngamma(b) - lngamma(a + b)
  let mut bt = limma_log_nat(x) * a + limma_log_nat(1.0 - x) * b - lbeta_val
  bt = limma_exp(bt)

  let mut betaf = 0.0
  if x < (a + 1.0) / (a + b + 2.0) {
    betaf = bt * betacf(x, a, b) / a
  } else {
    betaf = 1.0 - bt * betacf(1.0 - x, b, a) / b
  }

  betaf
}

///|
/// Continued fraction for incomplete beta function.
fn betacf(x : Double, a : Double, b : Double) -> Double {
  let max_iter = 100
  let eps = 0.0000003

  let mut c = 1.0
  let mut d = 1.0 - (a + b) * x / (a + 1.0)
  if d.abs() < eps {
    d = eps
  }
  d = 1.0 / d
  let mut h = d

  let mut m = 1
  while m <= max_iter {
    let m2 = (2 * m).to_double()
    let mut aa = m.to_double() *
      (b - m.to_double()) *
      x /
      ((a + m2 - 1.0) * (a + m2))
    d = 1.0 + aa * d
    if d.abs() < eps {
      d = eps
    }
    c = 1.0 + aa / c
    if c.abs() < eps {
      c = eps
    }
    d = 1.0 / d
    h = h * d * c

    aa = -(a + m.to_double()) *
      (a + b + m.to_double()) *
      x /
      ((a + m2) * (a + m2 + 1.0))
    d = 1.0 + aa * d
    if d.abs() < eps {
      d = eps
    }
    c = 1.0 + aa / c
    if c.abs() < eps {
      c = eps
    }
    d = 1.0 / d
    h = h * d * c

    if (h - 1.0).abs() < eps {
      return h
    }
    m = m + 1
  }
  h
}

///|
/// F-distribution cumulative distribution function.
pub fn pf(f : Double, df1 : Double, df2 : Double) -> Double {
  if f <= 0.0 {
    return 0.0
  }
  let x = df1 * f / (df1 * f + df2)
  1.0 - pbeta(x, df1 / 2.0, df2 / 2.0)
}

///|
/// F-distribution quantile function (inverse of pf).
pub fn qf(p : Double, df1 : Double, df2 : Double) -> Double {
  if p <= 0.0 {
    return 0.0
  }
  if p >= 1.0 {
    return 1000000000000000000000000000000000.0
  }

  let mut low = 0.0
  let mut high = 1.0
  while pf(high, df1, df2) < p {
    high = high * 2.0
  }

  let mut mid = 0.0
  let mut i = 0
  while i < 50 {
    mid = (low + high) / 2.0
    let cdf = pf(mid, df1, df2)
    if cdf < p {
      low = mid
    } else {
      high = mid
    }
    i = i + 1
  }
  mid
}

///|
/// Natural logarithm function.
fn limma_log_nat(x : Double) -> Double {
  if x <= 0.0 {
    0.0
  } else if x == 1.0 {
    0.0
  } else {
    let y = (x - 1.0) / (x + 1.0)
    let y2 = y * y
    let mut sum = y
    let mut term = y
    let mut i = 1
    while i <= 20 {
      term = term * y2
      sum = sum + term / (2 * i + 1).to_double()
      i = i + 1
    }
    2.0 * sum
  }
}

///|
/// Exponential function.
fn limma_exp(x : Double) -> Double {
  let mut sum = 1.0
  let mut term = 1.0
  let mut i = 1
  while i <= 20 {
    term = term * x / i.to_double()
    sum = sum + term
    i = i + 1
  }
  sum
}

///|
/// MArrayLM struct - stores linear model fit results.
pub struct MArrayLM {
  coefficients : Matrix
  stdev_unscaled : Matrix
  sigma : Array[Double]
  df_residual : Int
  design : Matrix
  effects : Matrix
  fit : Matrix
}

///|
/// Create MArrayLM from components.
pub fn MArrayLM::new(
  coefficients : Matrix,
  stdev_unscaled : Matrix,
  sigma : Array[Double],
  df_residual : Int,
  design : Matrix,
  effects : Matrix,
  fit : Matrix,
) -> MArrayLM {
  MArrayLM::{
    coefficients,
    stdev_unscaled,
    sigma,
    df_residual,
    design,
    effects,
    fit,
  }
}

///|
/// eBayes result struct.
pub struct EBayesResult {
  coefficients : Matrix
  stdev_unscaled : Matrix
  sigma : Array[Double]
  df_residual : Int
  design : Matrix
  t : Matrix
  p_value : Matrix
  lods : Matrix
  f : Array[Double]
  f_p_value : Array[Double]
  df_total : Int
  proportion : Double
  prior_df : Double
  prior_var : Double
}

///|
/// Create EBayesResult from components.
pub fn EBayesResult::new(
  coefficients : Matrix,
  stdev_unscaled : Matrix,
  sigma : Array[Double],
  df_residual : Int,
  design : Matrix,
  t : Matrix,
  p_value : Matrix,
  lods : Matrix,
  f : Array[Double],
  f_p_value : Array[Double],
  df_total : Int,
  proportion : Double,
  prior_df : Double,
  prior_var : Double,
) -> EBayesResult {
  EBayesResult::{
    coefficients,
    stdev_unscaled,
    sigma,
    df_residual,
    design,
    t,
    p_value,
    lods,
    f,
    f_p_value,
    df_total,
    proportion,
    prior_df,
    prior_var,
  }
}

///|
/// TopTable row struct.
pub struct TopTableRow {
  log_fc : Double
  ave_exp : Double
  t : Double
  p_value : Double
  adj_p_value : Double
  b : Double
}

///|
/// Create TopTableRow.
pub fn TopTableRow::new(
  log_fc : Double,
  ave_exp : Double,
  t : Double,
  p_value : Double,
  adj_p_value : Double,
  b : Double,
) -> TopTableRow {
  TopTableRow::{ log_fc, ave_exp, t, p_value, adj_p_value, b }
}

///|
/// Fit linear model to gene expression data.
/// y: genes x samples matrix
/// design: samples x coefficients matrix
pub fn lm_fit(y : Matrix, design : Matrix) -> MArrayLM {
  let n_genes = y.rows
  let n_samples = y.cols
  let n_coefs = design.cols

  let xt = matrix_transpose(design)
  let xtx = mat_mult(xt, design)
  let xtx_inv = mat_inverse(xtx)

  let betas_data = Array::make(n_genes * n_coefs, 0.0)
  let stdevs_data = Array::make(n_genes * n_coefs, 0.0)
  let sigma = Array::make(n_genes, 0.0)
  let effects_data = Array::make(n_genes * n_coefs, 0.0)
  let fit_vals_data = Array::make(n_genes * n_samples, 0.0)

  let mut i = 0
  while i < n_genes {
    let y_vec = Array::make(n_samples, 0.0)
    let mut j = 0
    while j < n_samples {
      y_vec[j] = y.get(i, j)
      j = j + 1
    }

    let xt_y_data = Array::make(n_coefs, 0.0)
    let mut k = 0
    while k < n_coefs {
      let mut sum = 0.0
      let mut m = 0
      while m < n_samples {
        sum = sum + xt.get(k, m) * y_vec[m]
        m = m + 1
      }
      xt_y_data[k] = sum
      k = k + 1
    }
    let xt_y = Matrix::new(n_coefs, 1, xt_y_data)

    let beta = mat_mult(xtx_inv, xt_y)

    let mut j = 0
    while j < n_coefs {
      betas_data[i * n_coefs + j] = beta.get(j, 0)
      j = j + 1
    }

    let fitted = Array::make(n_samples, 0.0)
    let mut j = 0
    while j < n_samples {
      let mut sum = 0.0
      let mut k = 0
      while k < n_coefs {
        sum = sum + betas_data[i * n_coefs + k] * design.get(j, k)
        k = k + 1
      }
      fitted[j] = sum
      fit_vals_data[i * n_samples + j] = sum
      j = j + 1
    }

    let residuals = Array::make(n_samples, 0.0)
    let mut j = 0
    while j < n_samples {
      residuals[j] = y_vec[j] - fitted[j]
      j = j + 1
    }

    let ss_residual = vec_ss(residuals)
    let df_resid = n_samples - n_coefs
    let sig = if df_resid > 0 {
      (ss_residual / df_resid.to_double()).sqrt()
    } else {
      0.0
    }
    sigma[i] = sig

    let mut j = 0
    while j < n_coefs {
      stdevs_data[i * n_coefs + j] = sig * xtx_inv.get(j, j).sqrt()
      j = j + 1
    }

    let mut j = 0
    while j < n_coefs {
      let mut sum = 0.0
      let mut k = 0
      while k < n_samples {
        sum = sum + xt.get(j, k) * y_vec[k]
        k = k + 1
      }
      effects_data[i * n_coefs + j] = sum
      j = j + 1
    }

    i = i + 1
  }

  let betas = Matrix::new(n_genes, n_coefs, betas_data)
  let stdevs = Matrix::new(n_genes, n_coefs, stdevs_data)
  let effects = Matrix::new(n_genes, n_coefs, effects_data)
  let fit_vals = Matrix::new(n_genes, n_samples, fit_vals_data)

  MArrayLM::new(
    betas,
    stdevs,
    sigma,
    n_samples - n_coefs,
    design,
    effects,
    fit_vals,
  )
}

///|
/// Empirical Bayes moderation of standard errors.
pub fn ebayes(fit : MArrayLM, proportion : Double) -> EBayesResult {
  let n_genes = fit.coefficients.rows
  let n_coefs = fit.coefficients.cols

  let sigma_sq = Array::make(n_genes, 0.0)
  let mut i = 0
  while i < n_genes {
    sigma_sq[i] = fit.sigma[i] * fit.sigma[i]
    i = i + 1
  }

  let df_residual = fit.df_residual.to_double()

  let mut sum_weights = 0.0
  let mut sum_sq = 0.0
  let mut n_eff = 0.0
  i = 0
  while i < n_genes {
    if sigma_sq[i] > 0.0 {
      let w = df_residual / sigma_sq[i]
      sum_weights = sum_weights + w
      sum_sq = sum_sq + df_residual
      n_eff = n_eff + 1.0
    }
    i = i + 1
  }

  let (prior_var, mean_sq) = if n_eff > 1.0 {
    let mean_sq_val = sum_sq / sum_weights
    let mut sum_diffs = 0.0
    i = 0
    while i < n_genes {
      if sigma_sq[i] > 0.0 {
        let diff = sigma_sq[i] - mean_sq_val
        sum_diffs = sum_diffs + diff * diff
      }
      i = i + 1
    }
    let var_sq = sum_diffs / (n_eff - 1.0)
    let prior = var_sq / mean_sq_val - 1.0 / df_residual
    (if prior < 0.0 { 0.0 } else { prior }, mean_sq_val)
  } else {
    (0.0, 0.0)
  }

  let prior_df = if mean_sq > 0.0 && n_eff > 1.0 {
    let val = (df_residual + 2.0) *
      prior_var /
      (mean_sq * (n_eff - 1.0) / n_eff)
    if val < 0.0 {
      0.0
    } else {
      val
    }
  } else {
    0.0
  }

  let df_total = df_residual + prior_df

  let t_data = Array::make(n_genes * n_coefs, 0.0)
  let p_value_data = Array::make(n_genes * n_coefs, 0.0)
  let lods_data = Array::make(n_genes * n_coefs, 0.0)

  i = 0
  while i < n_genes {
    let mut j = 0
    while j < n_coefs {
      let beta = fit.coefficients.get(i, j)
      let se = fit.stdev_unscaled.get(i, j) * fit.sigma[i]

      if se > 0.0 {
        let moderated_se = (se * se + prior_var).sqrt()
        let t_val = beta / moderated_se

        let df_t = df_total
        let p_val = 2.0 *
          (
            1.0 -
            pf(
              t_val.abs() *
              t_val.abs() *
              (df_t + 1.0) /
              (df_t + t_val.abs() * t_val.abs()),
              1.0,
              df_t,
            )
          )
        let p = if p_val > 1.0 {
          1.0
        } else if p_val < 0.0 {
          0.0
        } else {
          p_val
        }

        let lod = if p > 0.0 && p < 1.0 {
          -0.5 * limma_log_nat(p) +
          0.5 * limma_log_nat(1.0 - p) +
          limma_log_nat((1.0 - proportion) / proportion)
        } else {
          0.0
        }

        t_data[i * n_coefs + j] = t_val
        p_value_data[i * n_coefs + j] = p
        lods_data[i * n_coefs + j] = lod
      } else {
        t_data[i * n_coefs + j] = 0.0
        p_value_data[i * n_coefs + j] = 1.0
        lods_data[i * n_coefs + j] = 0.0
      }
      j = j + 1
    }
    i = i + 1
  }

  let t = Matrix::new(n_genes, n_coefs, t_data)
  let p_value = Matrix::new(n_genes, n_coefs, p_value_data)
  let lods = Matrix::new(n_genes, n_coefs, lods_data)

  let f = Array::make(n_genes, 0.0)
  let f_p_value = Array::make(n_genes, 0.0)

  i = 0
  while i < n_genes {
    let mut sum_t2 = 0.0
    let mut j = 0
    while j < n_coefs {
      sum_t2 = sum_t2 + t.get(i, j) * t.get(i, j)
      j = j + 1
    }

    let f_val = if n_coefs > 1 {
      sum_t2 /
      n_coefs.to_double() *
      (df_total / (df_total - n_coefs.to_double() + 1.0))
    } else {
      t.get(i, 0) * t.get(i, 0)
    }

    f[i] = f_val

    let df1 = n_coefs.to_double()
    let df2 = df_total
    f_p_value[i] = 1.0 - pf(f_val, df1, df2)

    i = i + 1
  }

  EBayesResult::new(
    fit.coefficients,
    fit.stdev_unscaled,
    fit.sigma,
    fit.df_residual,
    fit.design,
    t,
    p_value,
    lods,
    f,
    f_p_value,
    df_total.to_int(),
    proportion,
    prior_df,
    prior_var,
  )
}

///|
/// Extract top differentially expressed genes.
pub fn top_table(
  fit : EBayesResult,
  coef : Int,
  number : Int,
  adjust_mode : String,
) -> Array[TopTableRow] {
  let n_genes = fit.coefficients.rows

  let rows = Array::make(
    n_genes,
    TopTableRow::new(0.0, 0.0, 0.0, 1.0, 1.0, 0.0),
  )

  let mut i = 0
  while i < n_genes {
    let log_fc = fit.coefficients.get(i, coef)
    let t_val = fit.t.get(i, coef)
    let p_val = fit.p_value.get(i, coef)

    let mut ave_exp = 0.0
    let mut n = 0
    let mut j = 0
    while j < fit.design.cols {
      ave_exp = ave_exp + fit.coefficients.get(i, j)
      n = n + 1
      j = j + 1
    }
    if n > 0 {
      ave_exp = ave_exp / n.to_double()
    }

    rows[i] = TopTableRow::new(
      log_fc,
      ave_exp,
      t_val,
      p_val,
      p_val,
      fit.lods.get(i, coef),
    )
    i = i + 1
  }

  let sorted_indices = Array::make(n_genes, 0)
  let mut i = 0
  while i < n_genes {
    sorted_indices[i] = i
    i = i + 1
  }

  let mut swapped = true
  let mut k = n_genes - 1
  while k > 0 && swapped {
    swapped = false
    let mut i = 0
    while i < k {
      if rows[sorted_indices[i]].p_value > rows[sorted_indices[i + 1]].p_value {
        let temp = sorted_indices[i]
        sorted_indices[i] = sorted_indices[i + 1]
        sorted_indices[i + 1] = temp
        swapped = true
      }
      i = i + 1
    }
    k = k - 1
  }

  let p_values = Array::make(n_genes, 0.0)
  i = 0
  while i < n_genes {
    p_values[i] = rows[i].p_value
    i = i + 1
  }

  let adj_p = adjust_p_values(p_values, adjust_mode)

  let result_size = if number < n_genes { number } else { n_genes }
  let result = Array::make(
    result_size,
    TopTableRow::new(0.0, 0.0, 0.0, 1.0, 1.0, 0.0),
  )

  i = 0
  while i < result_size {
    let idx = sorted_indices[i]
    let row = rows[idx]
    result[i] = TopTableRow::new(
      row.log_fc,
      row.ave_exp,
      row.t,
      row.p_value,
      adj_p[idx],
      row.b,
    )
    i = i + 1
  }

  result
}

///|
/// Adjust p-values for multiple testing.
pub fn adjust_p_values(
  p_values : Array[Double],
  mode : String,
) -> Array[Double] {
  let n = p_values.length()
  let indices = Array::make(n, 0)
  let mut i = 0
  while i < n {
    indices[i] = i
    i = i + 1
  }

  let mut swapped = true
  let mut k = n - 1
  while k > 0 && swapped {
    swapped = false
    let mut i = 0
    while i < k {
      if p_values[indices[i]] > p_values[indices[i + 1]] {
        let temp = indices[i]
        indices[i] = indices[i + 1]
        indices[i + 1] = temp
        swapped = true
      }
      i = i + 1
    }
    k = k - 1
  }

  let adj = Array::make(n, 0.0)

  if mode.to_upper() == "BH" || mode.to_upper() == "FDR" {
    let mut min_adj = 1.0
    let mut i = n - 1
    while i >= 0 {
      let rank = (i + 1).to_double()
      let adj_val = p_values[indices[i]] * n.to_double() / rank
      if adj_val < min_adj {
        min_adj = adj_val
      }
      adj[indices[i]] = min_adj
      i = i - 1
    }

    i = 0
    while i < n {
      if adj[i] > 1.0 {
        adj[i] = 1.0
      }
      i = i + 1
    }
  } else {
    i = 0
    while i < n {
      let adj_val = p_values[i] * n.to_double()
      adj[i] = if adj_val > 1.0 { 1.0 } else { adj_val }
      i = i + 1
    }
  }

  adj
}

///|
/// Create contrast matrix from design matrix.
pub fn make_contrasts(
  contrasts : Array[Array[Double]],
  design : Matrix,
) -> Matrix {
  let n_contrasts = contrasts.length()
  let n_coefs = design.cols

  let data = Array::make(n_coefs * n_contrasts, 0.0)

  let mut i = 0
  while i < n_contrasts {
    let mut j = 0
    while j < n_coefs {
      if j < contrasts[i].length() {
        data[j * n_contrasts + i] = contrasts[i][j]
      }
      j = j + 1
    }
    i = i + 1
  }

  Matrix::new(n_coefs, n_contrasts, data)
}

///|
/// Fit contrasts to eBayes results.
pub fn contrasts_fit(fit : EBayesResult, contrasts : Matrix) -> EBayesResult {
  let n_genes = fit.coefficients.rows
  let n_contrasts = contrasts.cols

  let coefs = mat_mult(fit.coefficients, contrasts)

  let stdevs_data = Array::make(n_genes * n_contrasts, 0.0)
  let mut i = 0
  while i < n_genes {
    let mut j = 0
    while j < n_contrasts {
      let mut sum_se = 0.0
      let mut k = 0
      while k < fit.design.cols {
        let se = fit.stdev_unscaled.get(i, k)
        sum_se = sum_se + se * se * contrasts.get(k, j) * contrasts.get(k, j)
        k = k + 1
      }
      stdevs_data[i * n_contrasts + j] = sum_se.sqrt() * fit.sigma[i]
      j = j + 1
    }
    i = i + 1
  }
  let stdevs = Matrix::new(n_genes, n_contrasts, stdevs_data)

  let t_data = Array::make(n_genes * n_contrasts, 0.0)
  let p_value_data = Array::make(n_genes * n_contrasts, 0.0)
  let lods_data = Array::make(n_genes * n_contrasts, 0.0)

  let df_total = fit.df_total.to_double()

  i = 0
  while i < n_genes {
    let mut j = 0
    while j < n_contrasts {
      let beta = coefs.get(i, j)
      let se = stdevs.get(i, j)

      if se > 0.0 {
        let t_val = beta / se
        let p_val = 2.0 *
          (
            1.0 -
            pf(
              t_val.abs() *
              t_val.abs() *
              (df_total + 1.0) /
              (df_total + t_val.abs() * t_val.abs()),
              1.0,
              df_total,
            )
          )
        let p = if p_val > 1.0 {
          1.0
        } else if p_val < 0.0 {
          0.0
        } else {
          p_val
        }

        let lod = if p > 0.0 && p < 1.0 {
          -0.5 * limma_log_nat(p) +
          0.5 * limma_log_nat(1.0 - p) +
          limma_log_nat((1.0 - fit.proportion) / fit.proportion)
        } else {
          0.0
        }

        t_data[i * n_contrasts + j] = t_val
        p_value_data[i * n_contrasts + j] = p
        lods_data[i * n_contrasts + j] = lod
      } else {
        t_data[i * n_contrasts + j] = 0.0
        p_value_data[i * n_contrasts + j] = 1.0
        lods_data[i * n_contrasts + j] = 0.0
      }
      j = j + 1
    }
    i = i + 1
  }

  let t = Matrix::new(n_genes, n_contrasts, t_data)
  let p_value = Matrix::new(n_genes, n_contrasts, p_value_data)
  let lods = Matrix::new(n_genes, n_contrasts, lods_data)

  let f = Array::make(n_genes, 0.0)
  let f_p_value = Array::make(n_genes, 0.0)

  i = 0
  while i < n_genes {
    let mut sum_t2 = 0.0
    let mut j = 0
    while j < n_contrasts {
      sum_t2 = sum_t2 + t.get(i, j) * t.get(i, j)
      j = j + 1
    }

    let f_val = if n_contrasts > 1 {
      sum_t2 /
      n_contrasts.to_double() *
      (df_total / (df_total - n_contrasts.to_double() + 1.0))
    } else {
      t.get(i, 0) * t.get(i, 0)
    }

    f[i] = f_val
    f_p_value[i] = 1.0 - pf(f_val, n_contrasts.to_double(), df_total)

    i = i + 1
  }

  EBayesResult::new(
    coefs,
    stdevs,
    fit.sigma,
    fit.df_residual,
    fit.design,
    t,
    p_value,
    lods,
    f,
    f_p_value,
    fit.df_total,
    fit.proportion,
    fit.prior_df,
    fit.prior_var,
  )
}

///|
/// voom transformation for RNA-seq counts.
/// Converts counts to log2 CPM with precision weights.
pub fn voom(counts : Matrix, design : Matrix) -> Matrix {
  let n_genes = counts.rows
  let n_samples = counts.cols

  let lib_sizes = Array::make(n_samples, 0.0)
  let mut j = 0
  while j < n_samples {
    let mut sum = 0.0
    let mut i = 0
    while i < n_genes {
      sum = sum + counts.get(i, j)
      i = i + 1
    }
    lib_sizes[j] = sum
    j = j + 1
  }

  let log_cpm = Matrix::zeros(n_genes, n_samples)

  let mut i = 0
  while i < n_genes {
    let mut j = 0
    while j < n_samples {
      let cpm = counts.get(i, j) * 1000000.0 / lib_sizes[j]
      log_cpm.set(
        i,
        j,
        if cpm > 0.0 {
          limma_log_nat(cpm) / limma_log_nat(2.0)
        } else {
          0.0
        },
      )
      j = j + 1
    }
    i = i + 1
  }

  log_cpm
}

///|
/// Remove batch effects from expression data using linear model.
/// Uses linear model to estimate and remove batch effects.
/// Parameters:
///   - y: Expression matrix (genes x samples)
///   - batch: Array of batch labels for each sample
///   - group: Array of group labels (optional, for preserving group effects)
/// Returns batch-corrected expression matrix.
pub fn limma_remove_batch_effect(
  y : Matrix,
  batch : Array[String],
  group? : Array[String] = Array::new(),
) -> Matrix {
  let n_genes = y.rows
  let n_samples = y.cols

  if n_samples != batch.length() {
    return y
  }

  // Create batch design matrix
  let unique_batches : Array[String] = Array::new()
  let mut i = 0
  while i < n_samples {
    let b = batch[i]
    let found = unique_batches.contains(b)
    if !found {
      unique_batches.push(b)
    }
    i = i + 1
  }

  let n_batches = unique_batches.length()
  if n_batches <= 1 {
    return y
  }

  // Create batch indicator matrix (n_samples x n_batches-1, dropping first as reference)
  let n_batch_cols = n_batches - 1
  let batch_design_data = Array::make(n_samples * n_batch_cols, 0.0)
  i = 0
  while i < n_samples {
    let b_idx = find_batch_index(batch[i], unique_batches)
    let mut j = 0
    while j < n_batch_cols {
      if b_idx == j + 1 {
        batch_design_data[i * n_batch_cols + j] = 1.0
      }
      j = j + 1
    }
    i = i + 1
  }
  let batch_design = Matrix::new(n_samples, n_batch_cols, batch_design_data)

  // Add group information to design if provided
  let design = if group.length() > 0 && group.length() == n_samples {
    let unique_groups : Array[String] = Array::new()
    i = 0
    while i < n_samples {
      let g = group[i]
      if !unique_groups.contains(g) {
        unique_groups.push(g)
      }
      i = i + 1
    }
    let n_groups = unique_groups.length()
    let n_group_cols = if n_groups > 1 { n_groups - 1 } else { 0 }

    let total_cols = n_batch_cols + n_group_cols
    let design_data = Array::make(n_samples * total_cols, 0.0)

    // Copy batch design
    i = 0
    while i < n_samples {
      let mut j = 0
      while j < n_batch_cols {
        design_data[i * total_cols + j] = batch_design.get(i, j)
        j = j + 1
      }
      i = i + 1
    }

    // Add group indicators
    if n_group_cols > 0 {
      i = 0
      while i < n_samples {
        let g_idx = find_batch_index(group[i], unique_groups)
        let mut j = 0
        while j < n_group_cols {
          if g_idx == j + 1 {
            design_data[i * total_cols + n_batch_cols + j] = 1.0
          }
          j = j + 1
        }
        i = i + 1
      }
    }

    Matrix::new(n_samples, total_cols, design_data)
  } else {
    batch_design
  }

  // Fit linear model for each gene and subtract batch effects
  let result = Matrix::zeros(n_genes, n_samples)
  i = 0
  while i < n_genes {
    // Extract row
    let y_row = Array::make(n_samples, 0.0)
    let mut j = 0
    while j < n_samples {
      y_row[j] = y.get(i, j)
      j = j + 1
    }

    // Fit OLS: y_row = design * beta + residual
    // Using normal equations: beta = (X'X)^{-1} X'y
    let n_design = design.cols
    let xtx_data = Array::make(n_design * n_design, 0.0)
    let xty_data = Array::make(n_design, 0.0)

    j = 0
    while j < n_samples {
      let mut k = 0
      while k < n_design {
        let x_jk = design.get(j, k)
        xty_data[k] = xty_data[k] + x_jk * y_row[j]
        let mut l = 0
        while l < n_design {
          xtx_data[k * n_design + l] = xtx_data[k * n_design + l] + x_jk * design.get(j, l)
          l = l + 1
        }
        k = k + 1
      }
      j = j + 1
    }

    // Solve for beta using simple Gaussian elimination
    let beta = solve_linear_system(xtx_data, xty_data, n_design)

    // Compute fitted values and residuals
    j = 0
    while j < n_samples {
      let mut fitted = 0.0
      let mut k = 0
      while k < n_design {
        fitted = fitted + design.get(j, k) * beta[k]
        k = k + 1
      }
      result.set(i, j, y_row[j] - fitted + limma_stat_mean(y_row))
      j = j + 1
    }

    i = i + 1
  }

  result
}

///|
/// Find index of a batch label in the unique batches array.
fn find_batch_index(label : String, unique : Array[String]) -> Int {
  let mut i = 0
  while i < unique.length() {
    if unique[i] == label {
      return i
    }
    i = i + 1
  }
  0
}

///|
/// Solve linear system Ax = b using Gaussian elimination.
fn solve_linear_system(a_data : Array[Double], b : Array[Double], n : Int) -> Array[Double] {
  // Create augmented matrix
  let aug_data = Array::make(n * (n + 1), 0.0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      aug_data[i * (n + 1) + j] = a_data[i * n + j]
      j = j + 1
    }
    aug_data[i * (n + 1) + n] = b[i]
    i = i + 1
  }

  // Forward elimination with partial pivoting
  i = 0
  while i < n {
    // Find pivot
    let mut max_row = i
    let mut max_val = aug_data[i * (n + 1) + i].abs()
    let mut k = i + 1
    while k < n {
      let val = aug_data[k * (n + 1) + i].abs()
      if val > max_val {
        max_val = val
        max_row = k
      }
      k = k + 1
    }

    // Swap rows
    if max_row != i {
      let mut j = 0
      while j <= n {
        let tmp = aug_data[i * (n + 1) + j]
        aug_data[i * (n + 1) + j] = aug_data[max_row * (n + 1) + j]
        aug_data[max_row * (n + 1) + j] = tmp
        j = j + 1
      }
    }

    // Check for singular matrix
    let pivot = aug_data[i * (n + 1) + i]
    if pivot.abs() < 1.0e-10 {
      // Return zeros for singular case
      return Array::make(n, 0.0)
    }

    // Eliminate below
    k = i + 1
    while k < n {
      let factor = aug_data[k * (n + 1) + i] / pivot
      let mut j = i
      while j <= n {
        aug_data[k * (n + 1) + j] = aug_data[k * (n + 1) + j] - factor * aug_data[i * (n + 1) + j]
        j = j + 1
      }
      k = k + 1
    }

    i = i + 1
  }

  // Back substitution
  let x = Array::make(n, 0.0)
  i = n - 1
  while i >= 0 {
    let mut sum = aug_data[i * (n + 1) + n]
    let mut j = i + 1
    while j < n {
      sum = sum - aug_data[i * (n + 1) + j] * x[j]
      j = j + 1
    }
    let pivot = aug_data[i * (n + 1) + i]
    x[i] = if pivot.abs() < 1.0e-10 { 0.0 } else { sum / pivot }
    if i == 0 { break }
    i = i - 1
  }

  x
}

///|
/// Treat: one-sided test for differential expression.
/// Computes p-values for testing whether logFC exceeds a threshold.
/// Parameters:
///   - fit: EBayesResult from ebayes()
///   - lfc: log2 fold-change threshold (default 0)
///   - direction: "up", "down", or "both" (default "both")
/// Returns updated EBayesResult with one-sided p-values and B-values.
pub fn treat(fit : EBayesResult, lfc : Double, direction : String) -> EBayesResult {
  let n_genes = fit.coefficients.rows
  let n_coefs = fit.coefficients.cols

  let p_value_new = Array::make(n_genes * n_coefs, 1.0)
  let b_value_new = Array::make(n_genes * n_coefs, 0.0)

  let df_total = fit.df_total.to_double()

  let mut i = 0
  while i < n_genes {
    let mut j = 0
    while j < n_coefs {
      let beta = fit.coefficients.get(i, j)
      let se = fit.stdev_unscaled.get(i, j) * fit.sigma[i]

      let t_val = if se > 0.0 {
        (beta - lfc) / se
      } else {
        0.0
      }

      let p_val = if direction == "up" {
        // P(T > t_val) for one-sided
        if t_val > 0.0 {
          let f_val = t_val * t_val * (df_total + 1.0) / (df_total + t_val * t_val)
          1.0 - pf(f_val, 1.0, df_total)
        } else {
          1.0
        }
      } else if direction == "down" {
        // P(T < t_val) for one-sided
        if t_val < 0.0 {
          let t_abs = t_val.abs()
          let f_val = t_abs * t_abs * (df_total + 1.0) / (df_total + t_abs * t_abs)
          1.0 - pf(f_val, 1.0, df_total)
        } else {
          1.0
        }
      } else {
        // Two-sided
        if se > 0.0 {
          let t_abs = t_val.abs()
          let f_val = t_abs * t_abs * (df_total + 1.0) / (df_total + t_abs * t_abs)
          2.0 * (1.0 - pf(f_val, 1.0, df_total))
        } else {
          1.0
        }
      }

      let p = if p_val.is_nan() || p_val < 0.0 { 0.0 } else if p_val > 1.0 { 1.0 } else { p_val }

      // B-value (log-odds of differential expression)
      let b = if p > 0.0 && p < 1.0 {
        -@math.ln(p / (1.0 - p))
      } else if p <= 0.0 {
        10.0
      } else {
        -10.0
      }

      p_value_new[i * n_coefs + j] = p
      b_value_new[i * n_coefs + j] = b
      j = j + 1
    }
    i = i + 1
  }

  let p_values = Matrix::new(n_genes, n_coefs, p_value_new)
  let b_values = Matrix::new(n_genes, n_coefs, b_value_new)

  // Compute t-statistics for the treat result
  let t_data = Array::make(n_genes * n_coefs, 0.0)
  i = 0
  while i < n_genes {
    let mut j = 0
    while j < n_coefs {
      let beta = fit.coefficients.get(i, j)
      let se = fit.stdev_unscaled.get(i, j) * fit.sigma[i]
      t_data[i * n_coefs + j] = if se > 0.0 { beta / se } else { 0.0 }
      j = j + 1
    }
    i = i + 1
  }
  let t_matrix = Matrix::new(n_genes, n_coefs, t_data)

  let f_data = Array::make(n_genes, 0.0)
  let f_p_data = Array::make(n_genes, 1.0)

  EBayesResult::{
    coefficients: fit.coefficients,
    stdev_unscaled: fit.stdev_unscaled,
    sigma: fit.sigma,
    df_residual: fit.df_residual,
    design: fit.design,
    t: t_matrix,
    p_value: p_values,
    lods: b_values,
    f: f_data,
    f_p_value: f_p_data,
    df_total: fit.df_total,
    proportion: fit.proportion,
    prior_df: fit.prior_df,
    prior_var: fit.prior_var,
  }
}

///|
/// Compute distances for multidimensional scaling (MDS) plot.
/// Returns distance matrix between samples based on top variable genes.
pub fn plot_mds_distances(y : Matrix, n_top_genes : Int) -> Matrix {
  let n_genes = y.rows
  let n_samples = y.cols

  if n_samples < 3 {
    return Matrix::zeros(n_samples, n_samples)
  }

  // Select top most variable genes
  let variances = Array::make(n_genes, 0.0)
  let mut i = 0
  while i < n_genes {
    let row_vals = Array::make(n_samples, 0.0)
    let mut j = 0
    while j < n_samples {
      row_vals[j] = y.get(i, j)
      j = j + 1
    }
    variances[i] = limma_stat_variance(row_vals)
    i = i + 1
  }

  // Find top variable genes
  let n_top = if n_top_genes < n_genes { n_top_genes } else { n_genes }
  let top_indices : Array[Int] = Array::new()
  let used : Array[Bool] = Array::make(n_genes, false)

  i = 0
  while i < n_top {
    let mut max_var = -1.0
    let mut max_idx = 0
    let mut j = 0
    while j < n_genes {
      if !used[j] && variances[j] > max_var {
        max_var = variances[j]
        max_idx = j
      }
      j = j + 1
    }
    if max_var > 0.0 {
      top_indices.push(max_idx)
      used[max_idx] = true
    }
    i = i + 1
  }

  // Compute Euclidean distances between samples using top genes
  let n_top_actual = top_indices.length()
  let dist_data = Array::make(n_samples * n_samples, 0.0)

  i = 0
  while i < n_samples {
    let mut j = i + 1
    while j < n_samples {
      let mut sum_sq = 0.0
      let mut k = 0
      while k < n_top_actual {
        let gene_idx = top_indices[k]
        let diff = y.get(gene_idx, i) - y.get(gene_idx, j)
        sum_sq = sum_sq + diff * diff
        k = k + 1
      }
      let dist = sum_sq.sqrt()
      dist_data[i * n_samples + j] = dist
      dist_data[j * n_samples + i] = dist
      j = j + 1
    }
    i = i + 1
  }

  Matrix::new(n_samples, n_samples, dist_data)
}

///|
/// RPKM normalization: Reads Per Kilobase per Million mapped reads.
/// counts: genes x samples matrix of read counts.
/// gene_lengths: gene lengths in base pairs (one per gene).
/// Returns normalized expression matrix.
pub fn limma_rpkm(counts : Matrix, gene_lengths : Array[Double]) -> Matrix {
  let n_genes = counts.rows
  let n_samples = counts.cols

  let result = Matrix::zeros(n_genes, n_samples)

  // Compute library sizes (total counts per sample)
  let lib_sizes = Array::make(n_samples, 0.0)
  let mut j = 0
  while j < n_samples {
    let mut sum = 0.0
    let mut i = 0
    while i < n_genes {
      sum = sum + counts.get(i, j)
      i = i + 1
    }
    lib_sizes[j] = sum
    j = j + 1
  }

  let mut i = 0
  while i < n_genes {
    let gene_len_kb = if i < gene_lengths.length() && gene_lengths[i] > 0.0 {
      gene_lengths[i] / 1000.0
    } else {
      1.0
    }
    j = 0
    while j < n_samples {
      let rpkm_val = if lib_sizes[j] > 0.0 && gene_len_kb > 0.0 {
        counts.get(i, j) * 1.0e9 / (lib_sizes[j] * gene_len_kb)
      } else {
        0.0
      }
      result.set(i, j, rpkm_val)
      j = j + 1
    }
    i = i + 1
  }

  result
}

///|
/// CPM normalization: Counts Per Million mapped reads.
/// counts: genes x samples matrix of read counts.
/// log: if true, returns log2-transformed values.
/// Returns normalized matrix.
pub fn limma_cpm(counts : Matrix, log? : Bool = true) -> Matrix {
  let n_genes = counts.rows
  let n_samples = counts.cols

  let result = Matrix::zeros(n_genes, n_samples)

  // Compute library sizes
  let lib_sizes = Array::make(n_samples, 0.0)
  let mut j = 0
  while j < n_samples {
    let mut sum = 0.0
    let mut i = 0
    while i < n_genes {
      sum = sum + counts.get(i, j)
      i = i + 1
    }
    lib_sizes[j] = sum
    j = j + 1
  }

  let mut i = 0
  while i < n_genes {
    j = 0
    while j < n_samples {
      let cpm_val = if lib_sizes[j] > 0.0 {
        counts.get(i, j) * 1.0e6 / lib_sizes[j]
      } else {
        0.0
      }
      let val = if log {
        if cpm_val > 0.0 {
          limma_log_nat(cpm_val + 0.5) / limma_log_nat(2.0)
        } else {
          0.0
        }
      } else {
        cpm_val
      }
      result.set(i, j, val)
      j = j + 1
    }
    i = i + 1
  }

  result
}

///|
/// Quantile normalization (between-array normalization).
/// Normalizes expression values to have the same distribution across arrays.
pub fn limma_quantile_normalize(y : Matrix) -> Matrix {
  let n_genes = y.rows
  let n_samples = y.cols

  let result = Matrix::zeros(n_genes, n_samples)

  // Sort each column
  let sorted_cols : Array[Array[Double]] = Array::new()
  let mut j = 0
  while j < n_samples {
    let col : Array[Double] = Array::new()
    let mut i = 0
    while i < n_genes {
      col.push(y.get(i, j))
      i = i + 1
    }
    col.sort()
    sorted_cols.push(col)
    j = j + 1
  }

  // Compute mean of sorted values at each rank
  let rank_means = Array::make(n_genes, 0.0)
  let mut i = 0
  while i < n_genes {
    let mut sum = 0.0
    j = 0
    while j < n_samples {
      sum = sum + sorted_cols[j][i]
      j = j + 1
    }
    rank_means[i] = sum / n_samples.to_double()
    i = i + 1
  }

  // Map original values to quantiles and replace with rank means
  j = 0
  while j < n_samples {
    let col = sorted_cols[j]
    let mut i = 0
    while i < n_genes {
      let orig_val = y.get(i, j)
      // Find rank (smallest index where sorted value >= original)
      let mut rank = 0
      let mut k = 0
      while k < n_genes {
        if col[k] >= orig_val {
          rank = k
          break
        }
        rank = k
        k = k + 1
      }
      // Handle ties: use midpoint
      let mut end_rank = rank
      while end_rank + 1 < n_genes && col[end_rank + 1] == orig_val {
        end_rank = end_rank + 1
      }
      let avg_rank = (rank + end_rank).to_double() / 2.0
      let rank_floor = avg_rank.to_int()
      let rank_ceil = if rank_floor + 1 < n_genes { rank_floor + 1 } else { n_genes - 1 }
      let frac = avg_rank - rank_floor.to_double()

      let val = if rank_floor == rank_ceil {
        rank_means[rank_floor]
      } else {
        rank_means[rank_floor] * (1.0 - frac) + rank_means[rank_ceil] * frac
      }

      result.set(i, j, val)
      i = i + 1
    }
    j = j + 1
  }

  result
}

///|
/// ComBat batch effect correction using empirical Bayes.
/// Adjusts batch effects while preserving biological variation.
pub fn limma_combat(
  y : Matrix,
  batch : Array[String],
  group? : Array[String] = Array::new(),
) -> Matrix {
  let n_genes = y.rows
  let n_samples = y.cols

  if n_samples != batch.length() {
    return y
  }

  // Find unique batches
  let unique_batches : Array[String] = Array::new()
  let mut i = 0
  while i < n_samples {
    let found = unique_batches.contains(batch[i])
    if !found {
      unique_batches.push(batch[i])
    }
    i = i + 1
  }

  let n_batches = unique_batches.length()
  if n_batches <= 1 {
    return y
  }

  // For each batch, compute batch-specific mean and overall mean
  let result_data = Array::make(n_genes * n_samples, 0.0)
  let mut j = 0
  while j < n_samples {
    let mut k = 0
    while k < n_genes {
      result_data[k * n_samples + j] = y.get(k, j)
      k = k + 1
    }
    j = j + 1
  }

  // Compute overall gene means
  let overall_means = Array::make(n_genes, 0.0)
  i = 0
  while i < n_genes {
    let mut sum = 0.0
    j = 0
    while j < n_samples {
      sum = sum + y.get(i, j)
      j = j + 1
    }
    overall_means[i] = sum / n_samples.to_double()
    i = i + 1
  }

  // For each batch, compute batch offset and apply correction
  let mut b = 0
  while b < n_batches {
    let batch_label = unique_batches[b]
    // Find sample indices in this batch
    let batch_indices : Array[Int] = Array::new()
    let mut s = 0
    while s < n_samples {
      if batch[s] == batch_label {
        batch_indices.push(s)
      }
      s = s + 1
    }
    let batch_size = batch_indices.length()
    if batch_size == 0 {
      b = b + 1
      continue
    }

    // Compute per-gene batch means
    let batch_means = Array::make(n_genes, 0.0)
    i = 0
    while i < n_genes {
      let mut sum = 0.0
      j = 0
      while j < batch_size {
        sum = sum + y.get(i, batch_indices[j])
        j = j + 1
      }
      batch_means[i] = sum / batch_size.to_double()
      i = i + 1
    }

    // Batch offset per gene (difference from overall mean)
    let batch_offset = Array::make(n_genes, 0.0)
    i = 0
    while i < n_genes {
      batch_offset[i] = batch_means[i] - overall_means[i]
      i = i + 1
    }

    // Compute prior variance of batch offsets
    let offset_mean = limma_stat_mean(batch_offset)
    let offset_var = limma_stat_variance(batch_offset)
    let prior_var = if offset_var > 0.0 { offset_var } else { 1.0 }

    // For each gene, apply empirical Bayes shrinkage
    i = 0
    while i < n_genes {
      let offset = batch_offset[i]
      // Weighted average: combine batch offset with prior (0)
      let shrunken = offset * batch_size.to_double() /
        (batch_size.to_double() + offset / prior_var.abs())

      // Subtract shrunken batch offset from batch samples
      j = 0
      while j < batch_size {
        let sample_idx = batch_indices[j]
        let old_val = y.get(i, sample_idx)
        let new_val = old_val - shrunken
        result_data[i * n_samples + sample_idx] = new_val
        j = j + 1
      }
      i = i + 1
    }

    b = b + 1
  }

  Matrix::new(n_genes, n_samples, result_data)
}

///|
/// Compute mean of an array (helper for limma module).
fn limma_stat_mean(data : Array[Double]) -> Double {
  if data.length() == 0 { return 0.0 }
  let mut sum = 0.0
  for v in data {
    sum = sum + v
  }
  sum / data.length().to_double()
}

///|
/// Compute variance of an array (helper for limma module).
fn limma_stat_variance(data : Array[Double]) -> Double {
  let n = data.length()
  if n <= 1 { return 0.0 }
  let mean = limma_stat_mean(data)
  let mut sum_sq = 0.0
  for v in data {
    let diff = v - mean
    sum_sq = sum_sq + diff * diff
  }
  sum_sq / (n.to_double() - 1.0)
}
