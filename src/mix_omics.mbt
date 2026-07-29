///|
/// mixOmics - Multi-omics integration methods.
/// Provides PLS, sparse PLS (sPLS), and DIABLO algorithms
/// for integrating multiple omics data modalities.

// ============================================================================
// Data Structures
// ============================================================================

///|
/// Result of Partial Least Squares (PLS) regression.
pub struct PLSResult {
  x_scores : Array[Array[Double]]
  y_scores : Array[Array[Double]]
  x_loadings : Array[Array[Double]]
  y_loadings : Array[Array[Double]]
  x_weights : Array[Array[Double]]
  y_weights : Array[Array[Double]]
  explained_variance : Array[(Double, Double)]
  ncomp : Int
  n : Int
  p : Int
  q : Int
}

///|
/// Result of Sparse PLS (sPLS) with variable selection.
pub struct SPLSResult {
  x_scores : Array[Array[Double]]
  y_scores : Array[Array[Double]]
  x_loadings : Array[Array[Double]]
  y_loadings : Array[Array[Double]]
  x_weights : Array[Array[Double]]
  y_weights : Array[Array[Double]]
  x_selected : Array[Array[Bool]]
  y_selected : Array[Array[Bool]]
  explained_variance : Array[(Double, Double)]
  keep_x : Array[Int]
  keep_y : Array[Int]
  ncomp : Int
  n : Int
  p : Int
  q : Int
}

///|
/// Result of DIABLO (multi-block PLS-DA) integration.
pub struct DIABLOResult {
  block_scores : Array[Array[Array[Double]]]
  block_loadings : Array[Array[Array[Double]]]
  block_weights : Array[Array[Array[Double]]]
  block_selected : Array[Array[Array[Bool]]]
  ncomp : Int
  n_blocks : Int
  block_names : Array[String]
  block_sizes : Array[Int]
}

///|
/// Options for PLS computation.
pub struct PLSOptions {
  ncomp : Int
  mode : String
  scale : Bool
  max_iter : Int
  tol : Double
} derive(Eq, Debug)

///|
/// Options for sparse PLS computation.
pub struct SPLSOptions {
  ncomp : Int
  keep_x : Array[Int]
  keep_y : Array[Int]
  scale : Bool
  max_iter : Int
  tol : Double
} derive(Eq, Debug)

///|
/// Options for DIABLO multi-block integration.
pub struct DIABLOOptions {
  ncomp : Int
  keep_variables : Array[Array[Int]]
  design : Array[Array[Double]]
  scale : Bool
  max_iter : Int
  tol : Double
} derive(Eq, Debug)

// ============================================================================
// Helper Functions
// ============================================================================

///|
/// Center and scale matrix to zero mean and unit variance.
/// Returns scaled data, means and standard deviations.
pub fn mix_scale_matrix(
  data : Array[Array[Double]],
  n : Int,
  p : Int,
) -> (Array[Array[Double]], Array[Double], Array[Double]) {
  let means = Array::make(p, 0.0)
  let stds = Array::make(p, 0.0)

  for j = 0; j < p; j = j + 1 {
    let mut sum = 0.0
    for i = 0; i < n; i = i + 1 {
      sum = sum + data[i][j]
    }
    means[j] = sum / n.to_double()

    let mut var_sum = 0.0
    for i = 0; i < n; i = i + 1 {
      let diff = data[i][j] - means[j]
      var_sum = var_sum + diff * diff
    }
    let variance = if n > 1 { var_sum / (n.to_double() - 1.0) } else { 0.0 }
    stds[j] = if variance > 0.0 { variance.sqrt() } else { 1.0 }
  }

  let scaled : Array[Array[Double]] = []
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      row[j] = if stds[j] > 0.0 {
        (data[i][j] - means[j]) / stds[j]
      } else {
        data[i][j] - means[j]
      }
    }
    scaled.push(row)
  }

  (scaled, means, stds)
}

///|
/// Matrix-vector multiply: data (n x p) * x (p) -> result (n).
pub fn mix_mat_vec(
  data : Array[Array[Double]],
  n : Int,
  p : Int,
  x : Array[Double],
) -> Array[Double] {
  let result = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut sum = 0.0
    for j = 0; j < p; j = j + 1 {
      sum = sum + data[i][j] * x[j]
    }
    result[i] = sum
  }
  result
}

///|
/// Transpose matrix-vector multiply: data^T (p x n) * x (n) -> result (p).
pub fn mix_mat_t_vec(
  data : Array[Array[Double]],
  n : Int,
  p : Int,
  x : Array[Double],
) -> Array[Double] {
  let result = Array::make(p, 0.0)
  for j = 0; j < p; j = j + 1 {
    let mut sum = 0.0
    for i = 0; i < n; i = i + 1 {
      sum = sum + data[i][j] * x[i]
    }
    result[j] = sum
  }
  result
}

///|
/// Compute Euclidean norm of a vector.
pub fn mix_vector_norm(v : Array[Double]) -> Double {
  let mut sum = 0.0
  for i = 0; i < v.length(); i = i + 1 {
    sum = sum + v[i] * v[i]
  }
  sum.sqrt()
}

///|
/// Normalize vector to unit length. Returns zero vector if norm is zero.
pub fn mix_normalize(v : Array[Double]) -> Array[Double] {
  let norm = mix_vector_norm(v)
  let result = Array::make(v.length(), 0.0)
  if norm > 0.0 {
    for i = 0; i < v.length(); i = i + 1 {
      result[i] = v[i] / norm
    }
  }
  result
}

///|
/// Deflate matrix: data = data - scores * loadings^T.
pub fn mix_matrix_deflate(
  data : Array[Array[Double]],
  scores : Array[Double],
  loadings : Array[Double],
  n : Int,
  p : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      row[j] = data[i][j] - scores[i] * loadings[j]
    }
    result.push(row)
  }
  result
}

///|
/// Compute covariance matrix: X^T Y / n (p x q).
pub fn mix_covariance_matrix(
  x : Array[Array[Double]],
  y : Array[Array[Double]],
  n : Int,
  p : Int,
  q : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for j = 0; j < p; j = j + 1 {
    let row = Array::make(q, 0.0)
    for k = 0; k < q; k = k + 1 {
      let mut sum = 0.0
      for i = 0; i < n; i = i + 1 {
        sum = sum + x[i][j] * y[i][k]
      }
      row[k] = sum / n.to_double()
    }
    result.push(row)
  }
  result
}

///|
/// Soft thresholding operator: S(x, lambda) = sign(x) * max(0, |x| - lambda).
pub fn soft_threshold(x : Double, lambda : Double) -> Double {
  let abs_x = if x < 0.0 { -x } else { x }
  let shrunk = abs_x - lambda
  if shrunk <= 0.0 {
    0.0
  } else {
    if x < 0.0 { -shrunk } else { shrunk }
  }
}

///|
/// Select top k variables by absolute loading value.
/// Returns array of selected indices (length k or fewer).
pub fn select_top_k(values : Array[Double], k : Int) -> Array[Int] {
  let n = values.length()
  if k <= 0 || n == 0 {
    return []
  }
  let effective_k = if k < n { k } else { n }

  let indexed : Array[(Int, Double)] = []
  for i = 0; i < n; i = i + 1 {
    let abs_val = if values[i] < 0.0 { -values[i] } else { values[i] }
    indexed.push((i, abs_val))
  }

  indexed.sort_by(fn(a : (Int, Double), b : (Int, Double)) -> Int {
    if a.1 > b.1 { -1 } else if a.1 < b.1 { 1 } else { 0 }
  })

  let result = Array::make(effective_k, 0)
  for i = 0; i < effective_k; i = i + 1 {
    result[i] = indexed[i].0
  }
  result
}

///|
/// Compute total variance of a matrix (sum of column variances).
fn mix_total_variance(data : Array[Array[Double]], n : Int, p : Int) -> Double {
  let mut total = 0.0
  for j = 0; j < p; j = j + 1 {
    let mut mean = 0.0
    for i = 0; i < n; i = i + 1 {
      mean = mean + data[i][j]
    }
    mean = mean / n.to_double()

    let mut var_sum = 0.0
    for i = 0; i < n; i = i + 1 {
      let diff = data[i][j] - mean
      var_sum = var_sum + diff * diff
    }
    total = total + var_sum / (n.to_double() - 1.0)
  }
  total
}

///|
/// Dot product of two vectors.
fn nipals_dot(a : Array[Double], b : Array[Double]) -> Double {
  let len = if a.length() < b.length() { a.length() } else { b.length() }
  let mut sum = 0.0
  for i = 0; i < len; i = i + 1 {
    sum = sum + a[i] * b[i]
  }
  sum
}

///|
/// Euclidean distance between two vectors.
fn nipals_vec_norm(a : Array[Double], b : Array[Double]) -> Double {
  let len = if a.length() < b.length() { a.length() } else { b.length() }
  let mut sum = 0.0
  for i = 0; i < len; i = i + 1 {
    let d = a[i] - b[i]
    sum = sum + d * d
  }
  sum.sqrt()
}

///|
/// Create a zero-initialized 2D array (n x p).
fn make_zero_matrix(n : Int, p : Int) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for i = 0; i < n; i = i + 1 {
    result.push(Array::make(p, 0.0))
  }
  result
}

///|
/// Create a zero-initialized 3D array (dim1 x dim2 x dim3).
fn make_zero_tensor(
  dim1 : Int,
  dim2 : Int,
  dim3 : Int,
) -> Array[Array[Array[Double]]] {
  let result : Array[Array[Array[Double]]] = []
  for i = 0; i < dim1; i = i + 1 {
    let inner : Array[Array[Double]] = []
    for j = 0; j < dim2; j = j + 1 {
      inner.push(Array::make(dim3, 0.0))
    }
    result.push(inner)
  }
  result
}

// ============================================================================
// PLS Implementation
// ============================================================================

///|
/// Create default PLS options.
pub fn pls_new() -> PLSOptions {
  PLSOptions::{
    ncomp: 2,
    mode: "regression",
    scale: true,
    max_iter: 500,
    tol: 0.000001,
  }
}

///|
/// Create PLS options with all parameters.
pub fn pls_options_full(
  ncomp : Int,
  mode : String,
  scale : Bool,
  max_iter : Int,
  tol : Double,
) -> PLSOptions {
  PLSOptions::{ ncomp, mode, scale, max_iter, tol }
}

///|
/// Run PLS on X and Y data matrices.
/// X is (n x p), Y is (n x q).
pub fn run_pls(
  x_data : Array[Array[Double]],
  y_data : Array[Array[Double]],
  options : PLSOptions,
) -> PLSResult {
  let n = x_data.length()
  let p = if n > 0 { x_data[0].length() } else { 0 }
  let q = if y_data.length() > 0 {
    y_data[0].length()
  } else {
    0
  }

  let ncomp = if options.ncomp > n {
    n
  } else if options.ncomp > p {
    p
  } else if options.ncomp > q {
    q
  } else {
    options.ncomp
  }
  let effective_ncomp = if ncomp < 1 { 1 } else { ncomp }

  let (x_scaled, _, _) = if options.scale {
    mix_scale_matrix(x_data, n, p)
  } else {
    (x_data, Array::make(p, 0.0), Array::make(p, 1.0))
  }

  let (y_scaled, _, _) = if options.scale {
    mix_scale_matrix(y_data, n, q)
  } else {
    (y_data, Array::make(q, 0.0), Array::make(q, 1.0))
  }

  let total_x_var = mix_total_variance(x_scaled, n, p)
  let total_y_var = mix_total_variance(y_scaled, n, q)

  let mut x_current = x_scaled
  let mut y_current = y_scaled

  let x_scores = make_zero_matrix(n, effective_ncomp)
  let y_scores = make_zero_matrix(n, effective_ncomp)
  let x_loadings = make_zero_matrix(p, effective_ncomp)
  let y_loadings = make_zero_matrix(q, effective_ncomp)
  let x_weights = make_zero_matrix(p, effective_ncomp)
  let y_weights = make_zero_matrix(q, effective_ncomp)
  let explained_variance = Array::make(effective_ncomp, (0.0, 0.0))

  for comp = 0; comp < effective_ncomp; comp = comp + 1 {
    let (t, u, w_x, w_y, c_x, c_y) = pls_component(
      x_current,
      y_current,
      n,
      p,
      q,
      options.max_iter,
      options.tol,
    )

    for i = 0; i < n; i = i + 1 {
      x_scores[i][comp] = t[i]
      y_scores[i][comp] = u[i]
    }
    for j = 0; j < p; j = j + 1 {
      x_loadings[j][comp] = c_x[j]
      x_weights[j][comp] = w_x[j]
    }
    for k = 0; k < q; k = k + 1 {
      y_loadings[k][comp] = c_y[k]
      y_weights[k][comp] = w_y[k]
    }

    let mut x_var_explained = 0.0
    for j = 0; j < p; j = j + 1 {
      let loading_var = c_x[j] * c_x[j]
      x_var_explained = x_var_explained + loading_var
    }
    x_var_explained = if total_x_var > 0.0 {
      x_var_explained / total_x_var
    } else {
      0.0
    }

    let mut y_var_explained = 0.0
    for k = 0; k < q; k = k + 1 {
      let loading_var = c_y[k] * c_y[k]
      y_var_explained = y_var_explained + loading_var
    }
    y_var_explained = if total_y_var > 0.0 {
      y_var_explained / total_y_var
    } else {
      0.0
    }

    explained_variance[comp] = (x_var_explained, y_var_explained)

    x_current = mix_matrix_deflate(x_current, t, c_x, n, p)
    y_current = mix_matrix_deflate(y_current, t, c_y, n, q)
  }

  PLSResult::{
    x_scores,
    y_scores,
    x_loadings,
    y_loadings,
    x_weights,
    y_weights,
    explained_variance,
    ncomp: effective_ncomp,
    n,
    p,
    q,
  }
}

///|
/// Compute single PLS component using NIPALS algorithm.
/// Returns (t, u, w_x, w_y, c_x, c_y).
pub fn pls_component(
  x : Array[Array[Double]],
  y : Array[Array[Double]],
  n : Int,
  p : Int,
  q : Int,
  max_iter : Int,
  tol : Double,
) -> (
  Array[Double],
  Array[Double],
  Array[Double],
  Array[Double],
  Array[Double],
  Array[Double],
) {
  let (t, u, w_x, w_y) = nipals_pls(x, y, n, p, q, max_iter, tol)

  let t_norm_sq = nipals_dot(t, t)

  let c_x = Array::make(p, 0.0)
  let c_y = Array::make(q, 0.0)

  let xt = mix_mat_t_vec(x, n, p, t)
  if t_norm_sq > 0.0 {
    for j = 0; j < p; j = j + 1 {
      c_x[j] = xt[j] / t_norm_sq
    }
  }

  let yt = mix_mat_t_vec(y, n, q, t)
  if t_norm_sq > 0.0 {
    for k = 0; k < q; k = k + 1 {
      c_y[k] = yt[k] / t_norm_sq
    }
  } else {
    let yt_norm = nipals_dot(yt, yt)
    if yt_norm > 0.0 {
      for k = 0; k < q; k = k + 1 {
        c_y[k] = yt[k] / yt_norm
      }
    }
  }

  (t, u, w_x, w_y, c_x, c_y)
}

///|
/// NIPALS algorithm for PLS: computes scores and weights.
pub fn nipals_pls(
  x : Array[Array[Double]],
  y : Array[Array[Double]],
  n : Int,
  p : Int,
  q : Int,
  max_iter : Int,
  tol : Double,
) -> (Array[Double], Array[Double], Array[Double], Array[Double]) {
  let mut u = Array::make(n, 0.0)
  let mut w_x = Array::make(p, 0.0)
  let mut w_y = Array::make(q, 0.0)
  let mut t = Array::make(n, 0.0)

  for i = 0; i < n; i = i + 1 {
    u[i] = if q > 0 { y[i][0] } else { 0.0 }
  }

  let mut converged = false
  let mut iter = 0

  while !converged && iter < max_iter {
    let u_norm_sq = nipals_dot(u, u)
    if u_norm_sq < 0.0000000000000001 {
      break
    }

    let wx = mix_mat_t_vec(x, n, p, u)
    for j = 0; j < p; j = j + 1 {
      w_x[j] = wx[j] / u_norm_sq
    }
    w_x = mix_normalize(w_x)

    let t_new = mix_mat_vec(x, n, p, w_x)

    let t_norm_sq = nipals_dot(t_new, t_new)
    if t_norm_sq < 0.0000000000000001 {
      t = t_new
      break
    }

    let wy = mix_mat_t_vec(y, n, q, t_new)
    for k = 0; k < q; k = k + 1 {
      w_y[k] = wy[k] / t_norm_sq
    }
    w_y = mix_normalize(w_y)

    let u_new = mix_mat_vec(y, n, q, w_y)

    let diff = nipals_vec_norm(t_new, t)
    t = t_new
    u = u_new

    if diff < tol {
      converged = true
    }
    iter = iter + 1
  }

  (t, u, w_x, w_y)
}

// ============================================================================
// Sparse PLS Implementation
// ============================================================================

///|
/// Create default sparse PLS options.
pub fn spls_new() -> SPLSOptions {
  SPLSOptions::{
    ncomp: 2,
    keep_x: [],
    keep_y: [],
    scale: true,
    max_iter: 500,
    tol: 0.000001,
  }
}

///|
/// Create sparse PLS options with all parameters.
pub fn spls_options_full(
  ncomp : Int,
  keep_x : Array[Int],
  keep_y : Array[Int],
  scale : Bool,
  max_iter : Int,
  tol : Double,
) -> SPLSOptions {
  SPLSOptions::{ ncomp, keep_x, keep_y, scale, max_iter, tol }
}

///|
/// Run sparse PLS on X and Y data matrices.
pub fn run_spls(
  x_data : Array[Array[Double]],
  y_data : Array[Array[Double]],
  options : SPLSOptions,
) -> SPLSResult {
  let n = x_data.length()
  let p = if n > 0 { x_data[0].length() } else { 0 }
  let q = if y_data.length() > 0 {
    y_data[0].length()
  } else {
    0
  }

  let ncomp = if options.ncomp > n {
    n
  } else if options.ncomp > p {
    p
  } else if options.ncomp > q {
    q
  } else {
    options.ncomp
  }
  let effective_ncomp = if ncomp < 1 { 1 } else { ncomp }

  let (x_scaled, _, _) = if options.scale {
    mix_scale_matrix(x_data, n, p)
  } else {
    (x_data, Array::make(p, 0.0), Array::make(p, 1.0))
  }

  let (y_scaled, _, _) = if options.scale {
    mix_scale_matrix(y_data, n, q)
  } else {
    (y_data, Array::make(q, 0.0), Array::make(q, 1.0))
  }

  let total_x_var = mix_total_variance(x_scaled, n, p)
  let total_y_var = mix_total_variance(y_scaled, n, q)

  let mut x_current = x_scaled
  let mut y_current = y_scaled

  let x_scores = make_zero_matrix(n, effective_ncomp)
  let y_scores = make_zero_matrix(n, effective_ncomp)
  let x_loadings = make_zero_matrix(p, effective_ncomp)
  let y_loadings = make_zero_matrix(q, effective_ncomp)
  let x_weights = make_zero_matrix(p, effective_ncomp)
  let y_weights = make_zero_matrix(q, effective_ncomp)
  let x_selected = Array::make(p, Array::make(effective_ncomp, false))
  let y_selected = Array::make(q, Array::make(effective_ncomp, false))
  let explained_variance = Array::make(effective_ncomp, (0.0, 0.0))

  for comp = 0; comp < effective_ncomp; comp = comp + 1 {
    let kx = if comp < options.keep_x.length() {
      options.keep_x[comp]
    } else {
      p
    }
    let ky = if comp < options.keep_y.length() {
      options.keep_y[comp]
    } else {
      q
    }

    let (t, u, w_x, w_y, c_x, c_y, sel_x, sel_y) = spls_component(
      x_current,
      y_current,
      n,
      p,
      q,
      kx,
      ky,
      options.max_iter,
      options.tol,
    )

    for i = 0; i < n; i = i + 1 {
      x_scores[i][comp] = t[i]
      y_scores[i][comp] = u[i]
    }
    for j = 0; j < p; j = j + 1 {
      x_loadings[j][comp] = c_x[j]
      x_weights[j][comp] = w_x[j]
      x_selected[j][comp] = sel_x[j]
    }
    for k = 0; k < q; k = k + 1 {
      y_loadings[k][comp] = c_y[k]
      y_weights[k][comp] = w_y[k]
      y_selected[k][comp] = sel_y[k]
    }

    let mut x_var_explained = 0.0
    for j = 0; j < p; j = j + 1 {
      let loading_var = c_x[j] * c_x[j]
      x_var_explained = x_var_explained + loading_var
    }
    x_var_explained = if total_x_var > 0.0 {
      x_var_explained / total_x_var
    } else {
      0.0
    }

    let mut y_var_explained = 0.0
    for k = 0; k < q; k = k + 1 {
      let loading_var = c_y[k] * c_y[k]
      y_var_explained = y_var_explained + loading_var
    }
    y_var_explained = if total_y_var > 0.0 {
      y_var_explained / total_y_var
    } else {
      0.0
    }

    explained_variance[comp] = (x_var_explained, y_var_explained)

    x_current = mix_matrix_deflate(x_current, t, c_x, n, p)
    y_current = mix_matrix_deflate(y_current, t, c_y, n, q)
  }

  SPLSResult::{
    x_scores,
    y_scores,
    x_loadings,
    y_loadings,
    x_weights,
    y_weights,
    x_selected,
    y_selected,
    explained_variance,
    keep_x: options.keep_x,
    keep_y: options.keep_y,
    ncomp: effective_ncomp,
    n,
    p,
    q,
  }
}

///|
/// Compute single sparse PLS component with variable selection.
/// Returns scores, weights, loadings and selection vectors.
pub fn spls_component(
  x : Array[Array[Double]],
  y : Array[Array[Double]],
  n : Int,
  p : Int,
  q : Int,
  keep_x : Int,
  keep_y : Int,
  max_iter : Int,
  tol : Double,
) -> (
  Array[Double],
  Array[Double],
  Array[Double],
  Array[Double],
  Array[Double],
  Array[Double],
  Array[Bool],
  Array[Bool],
) {
  let (t_init, u_init, w_x_init, w_y_init) = nipals_pls(
    x,
    y,
    n,
    p,
    q,
    max_iter,
    tol,
  )

  let x_loadings = mix_mat_t_vec(x, n, p, t_init)
  let y_loadings = mix_mat_t_vec(y, n, q, t_init)

  let x_selected_idx = select_top_k(x_loadings, keep_x)
  let y_selected_idx = select_top_k(y_loadings, keep_y)

  let sel_x = Array::make(p, false)
  for i = 0; i < x_selected_idx.length(); i = i + 1 {
    sel_x[x_selected_idx[i]] = true
  }

  let sel_y = Array::make(q, false)
  for i = 0; i < y_selected_idx.length(); i = i + 1 {
    sel_y[y_selected_idx[i]] = true
  }

  let mut w_x = w_x_init
  let mut w_y = w_y_init
  let mut final_t = t_init
  let mut final_u = u_init

  let mut iter = 0
  let mut converged = false

  while !converged && iter < max_iter {
    let u_norm_sq = nipals_dot(final_u, final_u)
    if u_norm_sq < 0.0000000000000001 {
      break
    }

    let wx = mix_mat_t_vec(x, n, p, final_u)
    for j = 0; j < p; j = j + 1 {
      wx[j] = wx[j] / u_norm_sq
      if !sel_x[j] {
        wx[j] = soft_threshold(wx[j], 0.001)
      }
    }
    w_x = mix_normalize(wx)

    let t_new = mix_mat_vec(x, n, p, w_x)

    let t_norm_sq = nipals_dot(t_new, t_new)
    if t_norm_sq < 0.0000000000000001 {
      final_t = t_new
      break
    }

    let wy = mix_mat_t_vec(y, n, q, t_new)
    for k = 0; k < q; k = k + 1 {
      wy[k] = wy[k] / t_norm_sq
      if !sel_y[k] {
        wy[k] = soft_threshold(wy[k], 0.001)
      }
    }
    w_y = mix_normalize(wy)

    let u_new = mix_mat_vec(y, n, q, w_y)

    let diff = nipals_vec_norm(t_new, final_t)
    final_t = t_new
    final_u = u_new

    if diff < tol {
      converged = true
    }
    iter = iter + 1
  }

  let c_x = mix_mat_t_vec(x, n, p, final_t)
  let t_norm_sq = nipals_dot(final_t, final_t)
  if t_norm_sq > 0.0 {
    for j = 0; j < p; j = j + 1 {
      c_x[j] = c_x[j] / t_norm_sq
    }
  }

  let c_y = mix_mat_t_vec(y, n, q, final_t)
  if t_norm_sq > 0.0 {
    for k = 0; k < q; k = k + 1 {
      c_y[k] = c_y[k] / t_norm_sq
    }
  }

  (final_t, final_u, w_x, w_y, c_x, c_y, sel_x, sel_y)
}

// ============================================================================
// DIABLO Implementation
// ============================================================================

///|
/// Create default DIABLO options.
pub fn diablo_new() -> DIABLOOptions {
  DIABLOOptions::{
    ncomp: 2,
    keep_variables: [],
    design: [],
    scale: true,
    max_iter: 500,
    tol: 0.000001,
  }
}

///|
/// Create DIABLO options with all parameters.
pub fn diablo_options_full(
  ncomp : Int,
  keep_variables : Array[Array[Int]],
  design : Array[Array[Double]],
  scale : Bool,
  max_iter : Int,
  tol : Double,
) -> DIABLOOptions {
  DIABLOOptions::{ ncomp, keep_variables, design, scale, max_iter, tol }
}

///|
/// Run DIABLO multi-block integration.
/// blocks is array of (n x p_i) matrices.
pub fn run_diablo(
  blocks : Array[Array[Array[Double]]],
  block_names : Array[String],
  options : DIABLOOptions,
) -> DIABLOResult {
  let n_blocks = blocks.length()
  let n = if n_blocks > 0 { blocks[0].length() } else { 0 }
  let block_sizes = Array::make(n_blocks, 0)
  for b = 0; b < n_blocks; b = b + 1 {
    block_sizes[b] = if blocks[b].length() > 0 {
      blocks[b][0].length()
    } else {
      0
    }
  }

  let min_size = if block_sizes.length() > 0 { block_sizes[0] } else { 0 }
  let ncomp = if options.ncomp > n {
    n
  } else if options.ncomp > min_size {
    min_size
  } else {
    options.ncomp
  }
  let effective_ncomp = if ncomp < 1 { 1 } else { ncomp }

  let scaled_blocks : Array[Array[Array[Double]]] = []
  for b = 0; b < n_blocks; b = b + 1 {
    let pb = block_sizes[b]
    if options.scale {
      let (scaled, _, _) = mix_scale_matrix(blocks[b], n, pb)
      scaled_blocks.push(scaled)
    } else {
      scaled_blocks.push(blocks[b])
    }
  }

  let current_blocks = scaled_blocks

  let block_scores = make_zero_tensor(n_blocks, n, effective_ncomp)
  let block_loadings : Array[Array[Array[Double]]] = []
  let block_weights : Array[Array[Array[Double]]] = []
  let block_selected : Array[Array[Array[Bool]]] = []

  for b = 0; b < n_blocks; b = b + 1 {
    let pb = block_sizes[b]
    block_loadings.push(make_zero_matrix(pb, effective_ncomp))
    block_weights.push(make_zero_matrix(pb, effective_ncomp))
    let sel_block = Array::make(pb, Array::make(effective_ncomp, false))
    block_selected.push(sel_block)
  }

  for comp = 0; comp < effective_ncomp; comp = comp + 1 {
    let keep_for_comp = Array::make(n_blocks, 0)
    for b = 0; b < n_blocks; b = b + 1 {
      if comp < options.keep_variables.length() &&
        b < options.keep_variables[comp].length() {
        keep_for_comp[b] = options.keep_variables[comp][b]
      } else {
        keep_for_comp[b] = block_sizes[b]
      }
    }

    let (scores, weights, loadings, selected) = diablo_component(
      current_blocks,
      n_blocks,
      n,
      block_sizes,
      keep_for_comp,
      options.design,
      options.max_iter,
      options.tol,
    )

    for b = 0; b < n_blocks; b = b + 1 {
      let pb = block_sizes[b]
      for i = 0; i < n; i = i + 1 {
        block_scores[b][i][comp] = scores[b][i]
      }
      for j = 0; j < pb; j = j + 1 {
        block_loadings[b][j][comp] = loadings[b][j]
        block_weights[b][j][comp] = weights[b][j]
        block_selected[b][j][comp] = selected[b][j]
      }
    }

    for b = 0; b < n_blocks; b = b + 1 {
      let pb = block_sizes[b]
      current_blocks[b] = mix_matrix_deflate(
        current_blocks[b],
        scores[b],
        loadings[b],
        n,
        pb,
      )
    }
  }

  DIABLOResult::{
    block_scores,
    block_loadings,
    block_weights,
    block_selected,
    ncomp: effective_ncomp,
    n_blocks,
    block_names,
    block_sizes,
  }
}

///|
/// Compute single DIABLO component across all blocks.
pub fn diablo_component(
  blocks : Array[Array[Array[Double]]],
  n_blocks : Int,
  n : Int,
  block_sizes : Array[Int],
  keep_variables : Array[Int],
  design : Array[Array[Double]],
  max_iter : Int,
  tol : Double,
) -> (
  Array[Array[Double]],
  Array[Array[Double]],
  Array[Array[Double]],
  Array[Array[Bool]],
) {
  let block_weights : Array[Array[Double]] = []
  let block_scores : Array[Array[Double]] = []
  let block_loadings : Array[Array[Double]] = []
  let block_selected : Array[Array[Bool]] = []

  for b = 0; b < n_blocks; b = b + 1 {
    let pb = block_sizes[b]

    let w = Array::make(pb, 0.0)
    for j = 0; j < pb; j = j + 1 {
      w[j] = 1.0 / (pb.to_double()).sqrt()
    }
    block_weights.push(w)

    block_scores.push(Array::make(n, 0.0))

    block_loadings.push(Array::make(pb, 0.0))

    let sel = Array::make(pb, false)
    block_selected.push(sel)
  }

  for b = 0; b < n_blocks; b = b + 1 {
    block_scores[b] = mix_mat_vec(blocks[b], n, block_sizes[b], block_weights[b])
  }

  let mut global_score = Array::make(n, 0.0)
  for b = 0; b < n_blocks; b = b + 1 {
    for i = 0; i < n; i = i + 1 {
      global_score[i] = global_score[i] + block_scores[b][i]
    }
  }
  let gs_norm = mix_vector_norm(global_score)
  if gs_norm > 0.0 {
    for i = 0; i < n; i = i + 1 {
      global_score[i] = global_score[i] / gs_norm
    }
  }

  let mut converged = false
  let mut iter = 0

  while !converged && iter < max_iter {
    let new_global = Array::make(n, 0.0)

    for b = 0; b < n_blocks; b = b + 1 {
      let pb = block_sizes[b]

      let w = mix_mat_t_vec(blocks[b], n, pb, block_scores[b])

      let keep_b = if b < keep_variables.length() {
        keep_variables[b]
      } else {
        pb
      }

      if keep_b > 0 && keep_b < pb {
        let sel_idx = select_top_k(w, keep_b)
        for j = 0; j < pb; j = j + 1 {
          block_selected[b][j] = false
        }
        for i = 0; i < sel_idx.length(); i = i + 1 {
          block_selected[b][sel_idx[i]] = true
        }

        for j = 0; j < pb; j = j + 1 {
          if !block_selected[b][j] {
            w[j] = 0.0
          }
        }
      } else {
        for j = 0; j < pb; j = j + 1 {
          block_selected[b][j] = true
        }
      }

      let w_norm = mix_vector_norm(w)
      if w_norm > 0.0 {
        for j = 0; j < pb; j = j + 1 {
          w[j] = w[j] / w_norm
        }
      }

      block_weights[b] = w

      let new_score = mix_mat_vec(blocks[b], n, pb, w)
      block_scores[b] = new_score

      let design_weight = if b < design.length() {
        let sum_row = design[b]
        let mut dw = 0.0
        for c = 0; c < n_blocks; c = c + 1 {
          dw = dw + sum_row[c]
        }
        if dw > 0.0 { dw } else { 1.0 }
      } else {
        1.0
      }

      for i = 0; i < n; i = i + 1 {
        new_global[i] = new_global[i] + design_weight * new_score[i]
      }
    }

    let global_norm = mix_vector_norm(new_global)
    if global_norm > 0.0 {
      for i = 0; i < n; i = i + 1 {
        new_global[i] = new_global[i] / global_norm
      }
    }

    let diff = nipals_vec_norm(new_global, global_score)
    global_score = new_global

    if diff < tol {
      converged = true
    }
    iter = iter + 1
  }

  for b = 0; b < n_blocks; b = b + 1 {
    let pb = block_sizes[b]
    let loadings = mix_mat_t_vec(blocks[b], n, pb, block_scores[b])
    let scores_norm_sq = nipals_dot(block_scores[b], block_scores[b])
    if scores_norm_sq > 0.0 {
      for j = 0; j < pb; j = j + 1 {
        loadings[j] = loadings[j] / scores_norm_sq
      }
    }
    block_loadings[b] = loadings
  }

  (block_scores, block_weights, block_loadings, block_selected)
}

///|
/// Update block scores as weighted combination of individual block scores.
pub fn diablo_update_scores(
  blocks : Array[Array[Array[Double]]],
  weights : Array[Array[Double]],
  _scores : Array[Array[Double]],
  n_blocks : Int,
  n : Int,
) -> Array[Array[Double]] {
  let block_sizes = Array::make(n_blocks, 0)
  for b = 0; b < n_blocks; b = b + 1 {
    block_sizes[b] = if blocks[b].length() > 0 {
      blocks[b][0].length()
    } else {
      0
    }
  }

  let result : Array[Array[Double]] = []
  for b = 0; b < n_blocks; b = b + 1 {
    let pb = block_sizes[b]
    result.push(mix_mat_vec(blocks[b], n, pb, weights[b]))
  }
  result
}