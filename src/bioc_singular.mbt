///|
/// BiocSingular - SVD for dimensionality reduction.
/// Inspired by Bioconductor BiocSingular package.
/// Provides Singular Value Decomposition (SVD) for genomic data analysis.
/// Supports exact SVD, IRLBA, and randomized SVD methods.

///|
/// Results from Singular Value Decomposition.
pub struct SVResult {
  d : Array[Double]
  u : Array[Array[Double]]
  v : Array[Array[Double]]
  nrow : Int
  ncol : Int
  rank : Int
} derive(Eq)

///|
/// SVD computation method.
pub enum SVDMethod {
  ExactSVD
  IRLBA
  Randomized
} derive(Eq, Debug)

///|
/// Create ExactSVD method variant.
pub fn exact_svd_method() -> SVDMethod { ExactSVD }

///|
/// Create IRLBA method variant.
pub fn irlba_method() -> SVDMethod { IRLBA }

///|
/// Create Randomized method variant.
pub fn randomized_method() -> SVDMethod { Randomized }

///|
/// Options for SVD computation.
pub struct SVDOptions {
  method : SVDMethod
  rank : Int
  tol : Double
  max_iter : Int
  n_oversamples : Int
} derive(Eq, Debug)

// ============================================================================
// Random number generation utilities
// ============================================================================

///|
/// Simple Linear Congruential Generator state.
priv struct RNGState {
  state : Int
}

///|
/// Create a new RNG state with seed.
fn rng_new(seed : Int) -> RNGState {
  RNGState::{ state: seed }
}

///|
/// Generate next pseudo-random number in [0, 1).
fn rng_next(state : RNGState) -> (RNGState, Double) {
  let next_state = (state.state * 1103515245 + 12345) & 0x7FFFFFFF
  let value = next_state.to_double() / 2147483648.0
  (RNGState::{ state: next_state }, value)
}

///|
/// Generate standard normal random number using Box-Muller transform.
fn rng_normal(state : RNGState) -> (RNGState, Double) {
  let (s1, u1) = rng_next(state)
  let (s2, u2) = rng_next(s1)
  let pi = 3.14159265358979323846
  let z = (-2.0 * @math.ln(u1)).sqrt() * @math.cos(2.0 * pi * u2)
  (s2, z)
}

// ============================================================================
// Matrix operations
// ============================================================================

///|
/// Matrix-vector multiply: result = matrix * x.
/// Matrix is stored in row-major order: data[i * ncol + j].
pub fn svd_mat_vec(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
  x : Array[Double],
) -> Array[Double] {
  let result = Array::make(nrow, 0.0)
  for i = 0; i < nrow; i = i + 1 {
    let mut sum = 0.0
    for j = 0; j < ncol; j = j + 1 {
      sum = sum + matrix[i * ncol + j] * x[j]
    }
    result[i] = sum
  }
  result
}

///|
/// Transpose matrix-vector multiply: result = A^T * x.
/// Matrix is stored in row-major order: data[i * ncol + j].
pub fn svd_mat_t_vec(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
  x : Array[Double],
) -> Array[Double] {
  let result = Array::make(ncol, 0.0)
  for j = 0; j < ncol; j = j + 1 {
    let mut sum = 0.0
    for i = 0; i < nrow; i = i + 1 {
      sum = sum + matrix[i * ncol + j] * x[i]
    }
    result[j] = sum
  }
  result
}

///|
/// Matrix-matrix multiply: result = matrix * other.
/// matrix is row-major (nrow x ncol), other is (ncol x other_cols).
pub fn svd_mat_mat(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
  other : Array[Array[Double]],
  other_cols : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for i = 0; i < nrow; i = i + 1 {
    let row = Array::make(other_cols, 0.0)
    for j = 0; j < other_cols; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < ncol; k = k + 1 {
        sum = sum + matrix[i * ncol + k] * other[k][j]
      }
      row[j] = sum
    }
    result.push(row)
  }
  result
}

///|
/// Transpose flat matrix to 2D array.
/// Matrix is stored in row-major order: data[i * ncol + j].
pub fn svd_transpose(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for i = 0; i < ncol; i = i + 1 {
    let row = Array::make(nrow, 0.0)
    for j = 0; j < nrow; j = j + 1 {
      row[j] = matrix[j * ncol + i]
    }
    result.push(row)
  }
  result
}

///|
/// Full matrix multiply: C = A * B.
/// A is (nrow_a x ncol_a), B is (ncol_a x ncol_b).
/// Returns C as (nrow_a x ncol_b) 2D array.
pub fn svd_multiply(
  a : Array[Array[Double]],
  nrow_a : Int,
  ncol_a : Int,
  b : Array[Array[Double]],
  ncol_b : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for i = 0; i < nrow_a; i = i + 1 {
    let row = Array::make(ncol_b, 0.0)
    for j = 0; j < ncol_b; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < ncol_a; k = k + 1 {
        sum = sum + a[i][k] * b[k][j]
      }
      row[j] = sum
    }
    result.push(row)
  }
  result
}

///|
/// Compute A^T * A for a row-major matrix.
/// Returns an ncol x ncol symmetric matrix as 2D array.
fn svd_ata(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for i = 0; i < ncol; i = i + 1 {
    let row = Array::make(ncol, 0.0)
    for j = i; j < ncol; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < nrow; k = k + 1 {
        sum = sum + matrix[k * ncol + i] * matrix[k * ncol + j]
      }
      row[j] = sum
      row[i] = sum
    }
    result.push(row)
  }
  result
}

///|
/// Compute A * A^T for a row-major matrix.
/// Returns an nrow x nrow symmetric matrix as 2D array.
fn svd_aat(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = Array::new()
  for i = 0; i < nrow; i = i + 1 {
    let row = Array::make(nrow, 0.0)
    for j = i; j < nrow; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < ncol; k = k + 1 {
        sum = sum + matrix[i * ncol + k] * matrix[j * ncol + k]
      }
      row[j] = sum
      row[i] = sum
    }
    result.push(row)
  }
  result
}

// ============================================================================
// Jacobi eigenvalue algorithm for symmetric matrices
// ============================================================================

///|
/// Jacobi eigenvalue algorithm for symmetric n×n matrix.
/// Returns (eigenvalues, eigenvectors) where eigenvectors are columns of V.
/// Convergence tolerance: off-diagonal elements < tol.
pub fn jacobi_rotation(
  matrix : Array[Array[Double]],
  n : Int,
  max_iter : Int,
) -> (Array[Double], Array[Array[Double]]) {
  let a : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = matrix[i][j]
    }
    a.push(row)
  }

  let v : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    row[i] = 1.0
    v.push(row)
  }

  let tol = 1.0e-12
  let mut iter = 0

  while iter < max_iter {
    let mut max_off = 0.0
    let mut p = 0
    let mut q = 1

    for i = 0; i < n; i = i + 1 {
      for j = i + 1; j < n; j = j + 1 {
        if a[i][j].abs() > max_off {
          max_off = a[i][j].abs()
          p = i
          q = j
        }
      }
    }

    if max_off < tol {
      break
    }

    let app = a[p][p]
    let aqq = a[q][q]
    let apq = a[p][q]

    let theta = 0.5 * @math.atan2(2.0 * apq, app - aqq)
    let c = @math.cos(theta)
    let s = @math.sin(theta)

    let new_app = app * c * c + aqq * s * s + 2.0 * apq * c * s
    let new_aqq = app * s * s + aqq * c * c - 2.0 * apq * c * s
    let new_apq = (app - aqq) * c * s + apq * (c * c - s * s)

    a[p][p] = new_app
    a[q][q] = new_aqq
    a[p][q] = new_apq
    a[q][p] = new_apq

    for k = 0; k < n; k = k + 1 {
      if k != p && k != q {
        let aik = a[k][p]
        let ajk = a[k][q]
        a[k][p] = aik * c + ajk * s
        a[p][k] = a[k][p]
        a[k][q] = -aik * s + ajk * c
        a[q][k] = a[k][q]
      }
    }

    for k = 0; k < n; k = k + 1 {
      let vik = v[k][p]
      let vjk = v[k][q]
      v[k][p] = vik * c + vjk * s
      v[k][q] = -vik * s + vjk * c
    }

    iter = iter + 1
  }

  let eigenvalues = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    eigenvalues[i] = a[i][i]
  }

  (eigenvalues, v)
}

///|
/// Compute eigenvalues and eigenvectors of symmetric matrix.
/// Uses Jacobi eigenvalue algorithm.
pub fn eigh_symmetric(
  matrix : Array[Array[Double]],
  n : Int,
) -> (Array[Double], Array[Array[Double]]) {
  jacobi_rotation(matrix, n, 100)
}

///|
/// Sort eigenvalues and eigenvectors by descending eigenvalue.
fn sort_eigen(
  eigenvalues : Array[Double],
  eigenvectors : Array[Array[Double]],
  n : Int,
) -> (Array[Double], Array[Array[Double]]) {
  let indices = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    indices[i] = i
  }

  for i = 0; i < n - 1; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if eigenvalues[indices[j]] > eigenvalues[indices[i]] {
        let tmp = indices[i]
        indices[i] = indices[j]
        indices[j] = tmp
      }
    }
  }

  let sorted_eigenvalues = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    sorted_eigenvalues[i] = eigenvalues[indices[i]]
  }

  let sorted_eigenvectors : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = eigenvectors[i][indices[j]]
    }
    sorted_eigenvectors.push(row)
  }

  (sorted_eigenvalues, sorted_eigenvectors)
}

// ============================================================================
// Randomized SVD helpers
// ============================================================================

///|
/// Generate a random Gaussian matrix (rows x cols).
/// Uses Box-Muller transform with LCG.
pub fn svd_random_matrix(rows : Int, cols : Int) -> Array[Array[Double]] {
  let mut rng = rng_new(42)
  let result : Array[Array[Double]] = Array::new()
  for i = 0; i < rows; i = i + 1 {
    let row = Array::make(cols, 0.0)
    for j = 0; j < cols; j = j + 1 {
      let (new_rng, val) = rng_normal(rng)
      rng = new_rng
      row[j] = val
    }
    result.push(row)
  }
  result
}

///|
/// QR decomposition via modified Gram-Schmidt.
/// A is (n x n) matrix. Returns (Q, R) where Q is orthogonal and R is upper triangular.
pub fn svd_qr_decompose(
  a : Array[Array[Double]],
  n : Int,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let r : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    r.push(Array::make(n, 0.0))
  }

  let q : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let col = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      col[j] = a[j][i]
    }
    q.push(col)
  }

  for j = 0; j < n; j = j + 1 {
    for i = 0; i < j; i = i + 1 {
      let mut dot = 0.0
      for k = 0; k < n; k = k + 1 {
        dot = dot + q[i][k] * q[j][k]
      }
      r[i][j] = dot
      for k = 0; k < n; k = k + 1 {
        q[j][k] = q[j][k] - dot * q[i][k]
      }
    }

    let mut norm_sq_acc = 0.0
    for k = 0; k < n; k = k + 1 {
      norm_sq_acc = norm_sq_acc + q[j][k] * q[j][k]
    }
    let norm_val = norm_sq_acc.sqrt()
    r[j][j] = norm_val

    if norm_val > 1.0e-15 {
      for k = 0; k < n; k = k + 1 {
        q[j][k] = q[j][k] / norm_val
      }
    }
  }

  let q_t : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = q[j][i]
    }
    q_t.push(row)
  }

  (q_t, r)
}

///|
/// QR decomposition for n×k matrix (n >= k).
/// Returns Q (n×k, orthonormal columns) and R (k×k, upper triangular).
fn svd_qr_decompose_rect(
  a : Array[Array[Double]],
  n : Int,
  k : Int,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let r : Array[Array[Double]] = Array::new()
  for i = 0; i < k; i = i + 1 {
    r.push(Array::make(k, 0.0))
  }

  let q : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    q.push(Array::make(k, 0.0))
  }

  for col_idx = 0; col_idx < k; col_idx = col_idx + 1 {
    for i = 0; i < n; i = i + 1 {
      q[i][col_idx] = a[i][col_idx]
    }

    for j = 0; j < col_idx; j = j + 1 {
      let mut dot = 0.0
      for i = 0; i < n; i = i + 1 {
        dot = dot + q[i][j] * q[i][col_idx]
      }
      r[j][col_idx] = dot
      for i = 0; i < n; i = i + 1 {
        q[i][col_idx] = q[i][col_idx] - dot * q[i][j]
      }
    }

    let mut norm_sq = 0.0
    for i = 0; i < n; i = i + 1 {
      norm_sq = norm_sq + q[i][col_idx] * q[i][col_idx]
    }
    let norm_val = norm_sq.sqrt()
    r[col_idx][col_idx] = norm_val

    if norm_val > 1.0e-15 {
      for i = 0; i < n; i = i + 1 {
        q[i][col_idx] = q[i][col_idx] / norm_val
      }
    }
  }

  (q, r)
}

// ============================================================================
// Options constructors
// ============================================================================

///|
/// Create default SVD options.
/// Method: ExactSVD, rank: 0 (auto), tol: 1e-7, max_iter: 1000, n_oversamples: 10.
pub fn svd_options_new() -> SVDOptions {
  SVDOptions::{
    method: ExactSVD,
    rank: 0,
    tol: 1.0e-7,
    max_iter: 1000,
    n_oversamples: 10,
  }
}

///|
/// Create SVD options with specified method.
pub fn svd_options_with_method(method : SVDMethod) -> SVDOptions {
  let opts = svd_options_new()
  SVDOptions::{
    method,
    rank: opts.rank,
    tol: opts.tol,
    max_iter: opts.max_iter,
    n_oversamples: opts.n_oversamples,
  }
}

///|
/// Create SVD options with all parameters.
pub fn svd_options_full(
  rank : Int,
  tol : Double,
  max_iter : Int,
  n_oversamples : Int,
  method : SVDMethod,
) -> SVDOptions {
  SVDOptions::{
    method,
    rank,
    tol,
    max_iter,
    n_oversamples,
  }
}

// ============================================================================
// Exact SVD via eigendecomposition
// ============================================================================

///|
/// Compute exact SVD using eigendecomposition of A^T A (or A A^T).
/// Handles all edge cases: empty matrix, nrow=1, ncol=1, etc.
/// Returns SVResult with singular values sorted descending.
pub fn run_exact_svd(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
  rank : Int,
) -> SVResult {
  let mn = if nrow < ncol { nrow } else { ncol }
  let effective_rank = if rank <= 0 || rank > mn { mn } else { rank }

  if nrow == 0 || ncol == 0 {
    return SVResult::{
      d: Array::new(),
      u: Array::new(),
      v: Array::new(),
      nrow,
      ncol,
      rank: 0,
    }
  }

  if ncol <= nrow {
    let ata = svd_ata(matrix, nrow, ncol)
    let (eigenvalues, eigenvectors) = eigh_symmetric(ata, ncol)
    let (sorted_evals, sorted_evecs) = sort_eigen(eigenvalues, eigenvectors, ncol)

    let d = Array::make(effective_rank, 0.0)
    let v : Array[Array[Double]] = Array::new()
    let u : Array[Array[Double]] = Array::new()

    for i = 0; i < effective_rank; i = i + 1 {
      let sigma = if sorted_evals[i] > 0.0 {
        sorted_evals[i].sqrt()
      } else {
        0.0
      }
      d[i] = sigma

      let v_col = Array::make(ncol, 0.0)
      for j = 0; j < ncol; j = j + 1 {
        v_col[j] = sorted_evecs[j][i]
      }
      v.push(v_col)

      let u_col = Array::make(nrow, 0.0)
      let av_result = svd_mat_vec(matrix, nrow, ncol, v_col)
      if sigma > 1.0e-15 {
        for j = 0; j < nrow; j = j + 1 {
          u_col[j] = av_result[j] / sigma
        }
      } else {
        u_col[i % nrow] = 1.0
      }
      u.push(u_col)
    }

    SVResult::{ d, u, v, nrow, ncol, rank: effective_rank }
  } else {
    let aat = svd_aat(matrix, nrow, ncol)
    let (eigenvalues, eigenvectors) = eigh_symmetric(aat, nrow)
    let (sorted_evals, sorted_evecs) = sort_eigen(eigenvalues, eigenvectors, nrow)

    let d = Array::make(effective_rank, 0.0)
    let u : Array[Array[Double]] = Array::new()
    let v : Array[Array[Double]] = Array::new()

    for i = 0; i < effective_rank; i = i + 1 {
      let sigma = if sorted_evals[i] > 0.0 {
        sorted_evals[i].sqrt()
      } else {
        0.0
      }
      d[i] = sigma

      let u_col = Array::make(nrow, 0.0)
      for j = 0; j < nrow; j = j + 1 {
        u_col[j] = sorted_evecs[j][i]
      }
      u.push(u_col)

      let v_col = Array::make(ncol, 0.0)
      let atu_result = svd_mat_t_vec(matrix, nrow, ncol, u_col)
      if sigma > 1.0e-15 {
        for j = 0; j < ncol; j = j + 1 {
          v_col[j] = atu_result[j] / sigma
        }
      } else {
        v_col[i % ncol] = 1.0
      }
      v.push(v_col)
    }

    SVResult::{ d, u, v, nrow, ncol, rank: effective_rank }
  }
}

// ============================================================================
// IRLBA (Implicitly Restarted Lanczos Bidiagonalization Algorithm)
// ============================================================================

///|
/// Simple SVD of small bidiagonal matrix (n x n).
/// Returns singular values and left/right singular vectors.
fn svd_bidiagonal_small(
  diag : Array[Double],
  offdiag : Array[Double],
  n : Int,
) -> (Array[Double], Array[Array[Double]], Array[Array[Double]]) {
  let ncol = n
  let nrow = n

  let a : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    a.push(row)
  }
  for i = 0; i < n; i = i + 1 {
    a[i][i] = diag[i]
  }
  for i = 0; i < n - 1; i = i + 1 {
    a[i + 1][i] = offdiag[i]
  }

  let ata = svd_multiply(a, nrow, ncol, a, ncol)
  let (eigenvalues, eigenvectors) = eigh_symmetric(ata, n)
  let (sorted_evals, sorted_evecs) = sort_eigen(eigenvalues, eigenvectors, n)

  let d = Array::make(n, 0.0)
  let v : Array[Array[Double]] = Array::new()
  let u : Array[Array[Double]] = Array::new()

  for i = 0; i < n; i = i + 1 {
    let sigma = if sorted_evals[i] > 0.0 {
      sorted_evals[i].sqrt()
    } else {
      0.0
    }
    d[i] = sigma

    let v_col = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      v_col[j] = sorted_evecs[j][i]
    }
    v.push(v_col)

    let u_col = Array::make(n, 0.0)
    let av_mat = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < n; k = k + 1 {
        sum = sum + a[j][k] * v_col[k]
      }
      av_mat[j] = sum
    }
    if sigma > 1.0e-15 {
      for j = 0; j < n; j = j + 1 {
        u_col[j] = av_mat[j] / sigma
      }
    } else {
      u_col[i % n] = 1.0
    }
    u.push(u_col)
  }

  (d, u, v)
}

///|
/// IRLBA: Implicitly Restarted Lanczos Bidiagonalization Algorithm.
/// Computes the largest singular values of a matrix.
pub fn run_irlba(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
  rank : Int,
  options : SVDOptions,
) -> SVResult {
  let mn = if nrow < ncol { nrow } else { ncol }
  let effective_rank = if rank <= 0 || rank > mn { mn } else { rank }

  if nrow == 0 || ncol == 0 {
    return SVResult::{
      d: Array::new(),
      u: Array::new(),
      v: Array::new(),
      nrow,
      ncol,
      rank: 0,
    }
  }

  if effective_rank >= mn {
    return run_exact_svd(matrix, nrow, ncol, effective_rank)
  }

  let max_iter = options.max_iter
  let tol = options.tol
  let l = if effective_rank + 10 > mn - 1 {
    mn - 1
  } else {
    effective_rank + 10
  }

  let mut rng = rng_new(123)
  let u_vecs : Array[Array[Double]] = Array::new()
  let v_vecs : Array[Array[Double]] = Array::new()

  let x = Array::make(ncol, 0.0)
  for j = 0; j < ncol; j = j + 1 {
    let (new_rng, val) = rng_normal(rng)
    rng = new_rng
    x[j] = val
  }

  let mut norm_x = 0.0
  for j = 0; j < ncol; j = j + 1 {
    norm_x = norm_x + x[j] * x[j]
  }
  norm_x = norm_x.sqrt()
  if norm_x > 1.0e-15 {
    for j = 0; j < ncol; j = j + 1 {
      x[j] = x[j] / norm_x
    }
  }

  for i = 0; i < l; i = i + 1 {
    let u_col = Array::make(nrow, 0.0)
    u_vecs.push(u_col)
    let v_col = Array::make(ncol, 0.0)
    v_vecs.push(v_col)
  }

  let alphas = Array::make(l, 0.0)
  let betas = Array::make(l, 0.0)

  let mut v_prev = x.copy()
  let mut beta_prev = 0.0
  let prev_d = Array::make(effective_rank, 0.0)
  let mut converged = false

  for iter = 0; iter < max_iter; iter = iter + 1 {
    for j = 0; j < l; j = j + 1 {
      let u_new = svd_mat_vec(matrix, nrow, ncol, v_prev)
      for k = 0; k < nrow; k = k + 1 {
        u_new[k] = u_new[k] - beta_prev * u_vecs[if j > 0 { j - 1 } else { 0 }][k]
      }

      let mut alpha = 0.0
      for k = 0; k < nrow; k = k + 1 {
        alpha = alpha + u_new[k] * u_new[k]
      }
      alpha = alpha.sqrt()

      if alpha < 1.0e-15 {
        let (new_rng, val) = rng_normal(rng)
        rng = new_rng
        for k = 0; k < nrow; k = k + 1 {
          u_new[k] = val
        }
        let mut nrm = 0.0
        for k = 0; k < nrow; k = k + 1 {
          nrm = nrm + u_new[k] * u_new[k]
        }
        nrm = nrm.sqrt()
        if nrm > 1.0e-15 {
          for k = 0; k < nrow; k = k + 1 {
            u_new[k] = u_new[k] / nrm
          }
        }
        alpha = 0.0
      } else {
        for k = 0; k < nrow; k = k + 1 {
          u_new[k] = u_new[k] / alpha
        }
      }

      for k = 0; k < nrow; k = k + 1 {
        u_vecs[j][k] = u_new[k]
      }
      alphas[j] = alpha

      let v_new = svd_mat_t_vec(matrix, nrow, ncol, u_vecs[j])
      for k = 0; k < ncol; k = k + 1 {
        v_new[k] = v_new[k] - alpha * v_prev[k]
      }

      let mut beta = 0.0
      for k = 0; k < ncol; k = k + 1 {
        beta = beta + v_new[k] * v_new[k]
      }
      beta = beta.sqrt()

      if beta < 1.0e-15 {
        let (new_rng, val) = rng_normal(rng)
        rng = new_rng
        for k = 0; k < ncol; k = k + 1 {
          v_new[k] = val
        }
        let mut nrm = 0.0
        for k = 0; k < ncol; k = k + 1 {
          nrm = nrm + v_new[k] * v_new[k]
        }
        nrm = nrm.sqrt()
        if nrm > 1.0e-15 {
          for k = 0; k < ncol; k = k + 1 {
            v_new[k] = v_new[k] / nrm
          }
        }
        beta = 0.0
      } else {
        for k = 0; k < ncol; k = k + 1 {
          v_new[k] = v_new[k] / beta
        }
      }

      for k = 0; k < ncol; k = k + 1 {
        v_vecs[j][k] = v_new[k]
      }
      betas[j] = beta

      v_prev = v_new.copy()
      beta_prev = beta
    }

    let (s, u_hat, v_hat) = svd_bidiagonal_small(alphas, betas, l)

    let d = Array::make(effective_rank, 0.0)
    let u_result : Array[Array[Double]] = Array::new()
    let v_result : Array[Array[Double]] = Array::new()

    for i = 0; i < effective_rank; i = i + 1 {
      d[i] = s[i]

      let u_col = Array::make(nrow, 0.0)
      for k = 0; k < nrow; k = k + 1 {
        for j = 0; j < l; j = j + 1 {
          u_col[k] = u_col[k] + u_vecs[j][k] * u_hat[j][i]
        }
      }
      u_result.push(u_col)

      let v_col = Array::make(ncol, 0.0)
      for k = 0; k < ncol; k = k + 1 {
        for j = 0; j < l; j = j + 1 {
          v_col[k] = v_col[k] + v_vecs[j][k] * v_hat[j][i]
        }
      }
      v_result.push(v_col)
    }

    if iter > 0 {
      let mut d_diff = 0.0
      for i = 0; i < effective_rank; i = i + 1 {
        let diff = (d[i] - prev_d[i]).abs()
        d_diff = d_diff + diff
      }
      if d_diff < tol {
        converged = true
      }
    }

    if converged && iter > 0 {
      return SVResult::{
        d,
        u: u_result,
        v: v_result,
        nrow,
        ncol,
        rank: effective_rank,
      }
    }

    for i = 0; i < effective_rank; i = i + 1 {
      prev_d[i] = d[i]
    }

    v_prev = v_result[effective_rank - 1].copy()
    let mut nrm = 0.0
    for k = 0; k < ncol; k = k + 1 {
      nrm = nrm + v_prev[k] * v_prev[k]
    }
    nrm = nrm.sqrt()
    if nrm > 1.0e-15 {
      for k = 0; k < ncol; k = k + 1 {
        v_prev[k] = v_prev[k] / nrm
      }
    }
    beta_prev = 0.0
  }

  SVResult::{
    d: Array::new(),
    u: Array::new(),
    v: Array::new(),
    nrow,
    ncol,
    rank: 0,
  }
}

// ============================================================================
// Randomized SVD
// ============================================================================

///|
/// Randomized SVD for large matrices.
/// Uses random projection to reduce dimensions before computing SVD.
pub fn run_randomized_svd(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
  rank : Int,
  options : SVDOptions,
) -> SVResult {
  let mn = if nrow < ncol { nrow } else { ncol }
  let effective_rank = if rank <= 0 || rank > mn { mn } else { rank }

  if nrow == 0 || ncol == 0 {
    return SVResult::{
      d: Array::new(),
      u: Array::new(),
      v: Array::new(),
      nrow,
      ncol,
      rank: 0,
    }
  }

  if effective_rank >= mn {
    return run_exact_svd(matrix, nrow, ncol, effective_rank)
  }

  let n_oversamples = options.n_oversamples
  let sketch_size = effective_rank + n_oversamples
  let l = if sketch_size > mn { mn } else { sketch_size }

  let omega = svd_random_matrix(ncol, l)
  let y : Array[Array[Double]] = Array::new()
  for i = 0; i < nrow; i = i + 1 {
    let row = Array::make(l, 0.0)
    for j = 0; j < l; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < ncol; k = k + 1 {
        sum = sum + matrix[i * ncol + k] * omega[k][j]
      }
      row[j] = sum
    }
    y.push(row)
  }

  let (q, _r) = svd_qr_decompose_rect(y, nrow, l)

  let b : Array[Array[Double]] = Array::new()
  for i = 0; i < l; i = i + 1 {
    let row = Array::make(ncol, 0.0)
    for j = 0; j < ncol; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < nrow; k = k + 1 {
        sum = sum + q[k][i] * matrix[k * ncol + j]
      }
      row[j] = sum
    }
    b.push(row)
  }

  let b_flat = Array::make(l * ncol, 0.0)
  for i = 0; i < l; i = i + 1 {
    for j = 0; j < ncol; j = j + 1 {
      b_flat[i * ncol + j] = b[i][j]
    }
  }

  let b_result = run_exact_svd(b_flat, l, ncol, effective_rank)

  if b_result.d.length() == 0 {
    return SVResult::{
      d: Array::new(),
      u: Array::new(),
      v: Array::new(),
      nrow,
      ncol,
      rank: 0,
    }
  }

  let d = b_result.d
  let v = b_result.v

  let u : Array[Array[Double]] = Array::new()
  for i = 0; i < effective_rank; i = i + 1 {
    let u_col = Array::make(nrow, 0.0)
    for j = 0; j < nrow; j = j + 1 {
      let mut sum = 0.0
      for k = 0; k < l; k = k + 1 {
        sum = sum + q[j][k] * b_result.u[k][i]
      }
      u_col[j] = sum
    }
    u.push(u_col)
  }

  SVResult::{ d, u, v, nrow, ncol, rank: effective_rank }
}

// ============================================================================
// Main SVD entry point
// ============================================================================

///|
/// Run SVD with specified options.
/// Dispatches to appropriate SVD method based on options.method.
pub fn run_svd(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
  options : SVDOptions,
) -> SVResult {
  match options.method {
    ExactSVD => {
      run_exact_svd(matrix, nrow, ncol, options.rank)
    }
    IRLBA => {
      run_irlba(matrix, nrow, ncol, options.rank, options)
    }
    Randomized => {
      run_randomized_svd(matrix, nrow, ncol, options.rank, options)
    }
  }
}

///|
/// Create an SVResult with specified components.
pub fn sv_result_new(
  d : Array[Double],
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  nrow : Int,
  ncol : Int,
  rank : Int,
) -> SVResult {
  SVResult::{ d, u, v, nrow, ncol, rank }
}

///|
/// Get singular values from SVResult.
pub fn sv_result_d(result : SVResult) -> Array[Double] {
  result.d
}

///|
/// Get left singular vectors from SVResult.
pub fn sv_result_u(result : SVResult) -> Array[Array[Double]] {
  result.u
}

///|
/// Get right singular vectors from SVResult.
pub fn sv_result_v(result : SVResult) -> Array[Array[Double]] {
  result.v
}

///|
/// Get number of rows from SVResult.
pub fn sv_result_nrow(result : SVResult) -> Int {
  result.nrow
}

///|
/// Get number of columns from SVResult.
pub fn sv_result_ncol(result : SVResult) -> Int {
  result.ncol
}

///|
/// Get rank from SVResult.
pub fn sv_result_rank(result : SVResult) -> Int {
  result.rank
}

///|
/// Create example matrix for testing SVD.
pub fn create_test_matrix() -> (Array[Double], Int, Int) {
  let nrow = 5
  let ncol = 4
  let data = [
    1.0, 2.0, 3.0, 4.0,
    2.0, 3.0, 4.0, 5.0,
    3.0, 4.0, 5.0, 6.0,
    4.0, 5.0, 6.0, 7.0,
    5.0, 6.0, 7.0, 8.0,
  ]
  (data, nrow, ncol)
}

///|
/// Compute truncated SVD returning only top k components.
pub fn run_svd_truncated(
  matrix : Array[Double],
  nrow : Int,
  ncol : Int,
  k : Int,
) -> SVResult {
  let opts = svd_options_full(k, 1.0e-7, 1000, 10, ExactSVD)
  run_svd(matrix, nrow, ncol, opts)
}