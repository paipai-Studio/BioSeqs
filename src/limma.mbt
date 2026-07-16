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
pub fn Matrix::from_2d(rows : Int, cols : Int, data : Array[Array[Double]]) -> Matrix {
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
pub fn mat_transpose(a : Matrix) -> Matrix {
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
  let cof = [76.18009172947146, -86.50532032941677, 24.01409824083091,
             -1.231739572450155, 0.001208650973866179, -0.000005395239384953]
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
  if x <= 0.0 { return 0.0 }
  if x >= 1.0 { return 1.0 }
  
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
  if d.abs() < eps { d = eps }
  d = 1.0 / d
  let mut h = d
  
  let mut m = 1
  while m <= max_iter {
    let m2 = (2 * m).to_double()
    let mut aa = m.to_double() * (b - m.to_double()) * x / ((a + m2 - 1.0) * (a + m2))
    d = 1.0 + aa * d
    if d.abs() < eps { d = eps }
    c = 1.0 + aa / c
    if c.abs() < eps { c = eps }
    d = 1.0 / d
    h = h * d * c
    
    aa = -(a + m.to_double()) * (a + b + m.to_double()) * x / ((a + m2) * (a + m2 + 1.0))
    d = 1.0 + aa * d
    if d.abs() < eps { d = eps }
    c = 1.0 + aa / c
    if c.abs() < eps { c = eps }
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
  if f <= 0.0 { return 0.0 }
  let x = df1 * f / (df1 * f + df2)
  1.0 - pbeta(x, df1 / 2.0, df2 / 2.0)
}

///|
/// F-distribution quantile function (inverse of pf).
pub fn qf(p : Double, df1 : Double, df2 : Double) -> Double {
  if p <= 0.0 { return 0.0 }
  if p >= 1.0 { return 1000000000000000000000000000000000.0 }
  
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
  TopTableRow::{
    log_fc,
    ave_exp,
    t,
    p_value,
    adj_p_value,
    b,
  }
}

///|
/// Fit linear model to gene expression data.
/// y: genes x samples matrix
/// design: samples x coefficients matrix
pub fn lm_fit(y : Matrix, design : Matrix) -> MArrayLM {
  let n_genes = y.rows
  let n_samples = y.cols
  let n_coefs = design.cols
  
  let xt = mat_transpose(design)
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
    let sig = if df_resid > 0 { (ss_residual / df_resid.to_double()).sqrt() } else { 0.0 }
    sigma[i] = sig
    
    let mut j = 0
    while j < n_coefs {
      stdevs_data[i * n_coefs + j] = sig * (xtx_inv.get(j, j)).sqrt()
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
  
  MArrayLM::new(betas, stdevs, sigma, n_samples - n_coefs, design, effects, fit_vals)
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
    let val = (df_residual + 2.0) * prior_var / (mean_sq * (n_eff - 1.0) / n_eff)
    if val < 0.0 { 0.0 } else { val }
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
        let p_val = 2.0 * (1.0 - pf(t_val.abs() * t_val.abs() * (df_t + 1.0) / (df_t + t_val.abs() * t_val.abs()), 1.0, df_t))
        let p = if p_val > 1.0 { 1.0 } else if p_val < 0.0 { 0.0 } else { p_val }
        
        let lod = if p > 0.0 && p < 1.0 {
          -0.5 * limma_log_nat(p) + 0.5 * limma_log_nat(1.0 - p) + limma_log_nat((1.0 - proportion) / proportion)
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
      (sum_t2 / n_coefs.to_double()) * (df_total / (df_total - n_coefs.to_double() + 1.0))
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
  
  let rows = Array::make(n_genes, TopTableRow::new(0.0, 0.0, 0.0, 1.0, 1.0, 0.0))
  
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
    
    rows[i] = TopTableRow::new(log_fc, ave_exp, t_val, p_val, p_val, fit.lods.get(i, coef))
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
  let result = Array::make(result_size, TopTableRow::new(0.0, 0.0, 0.0, 1.0, 1.0, 0.0))
  
  i = 0
  while i < result_size {
    let idx = sorted_indices[i]
    let row = rows[idx]
    result[i] = TopTableRow::new(row.log_fc, row.ave_exp, row.t, row.p_value, adj_p[idx], row.b)
    i = i + 1
  }
  
  result
}

///|
/// Adjust p-values for multiple testing.
pub fn adjust_p_values(p_values : Array[Double], mode : String) -> Array[Double] {
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
pub fn make_contrasts(contrasts : Array[Array[Double]], design : Matrix) -> Matrix {
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
      stdevs_data[i * n_contrasts + j] = (sum_se).sqrt() * fit.sigma[i]
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
        let p_val = 2.0 * (1.0 - pf(t_val.abs() * t_val.abs() * (df_total + 1.0) / (df_total + t_val.abs() * t_val.abs()), 1.0, df_total))
        let p = if p_val > 1.0 { 1.0 } else if p_val < 0.0 { 0.0 } else { p_val }
        
        let lod = if p > 0.0 && p < 1.0 {
          -0.5 * limma_log_nat(p) + 0.5 * limma_log_nat(1.0 - p) + limma_log_nat((1.0 - fit.proportion) / fit.proportion)
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
      (sum_t2 / n_contrasts.to_double()) * (df_total / (df_total - n_contrasts.to_double() + 1.0))
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
      log_cpm.set(i, j, if cpm > 0.0 { limma_log_nat(cpm) / limma_log_nat(2.0) } else { 0.0 })
      j = j + 1
    }
    i = i + 1
  }
  
  log_cpm
}
