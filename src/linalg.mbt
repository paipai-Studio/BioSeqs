///|
/// Shared dense linear-algebra primitives.
///
/// The canonical implementations of the operations that had previously been
/// re-implemented across many modules: symmetric eigen-decomposition (Jacobi),
/// matrix inverse (Gauss-Jordan with partial pivoting), matrix multiply and
/// transpose. Other modules that need these should delegate here rather than
/// keep a local copy — a duplicated Jacobi with a sign error was the root cause
/// of a `2.67e14` eigenvalue bug in the past.

///|
/// Jacobi eigenvalue decomposition of a symmetric n×n matrix.
/// Returns (eigenvalues, eigenvector-matrix) where the columns of the returned
/// matrix are the eigenvectors. Eigenvalues are in whatever order the sweep
/// leaves them (not sorted).
pub fn linalg_jacobi(
  matrix : Array[Array[Double]],
  n : Int,
  max_iter : Int,
) -> (Array[Double], Array[Array[Double]]) {
  let a : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = matrix[i][j]
    }
    a.push(row)
  }

  let v : Array[Array[Double]] = Array(capacity=0)
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

    a[p][p] = new_app
    a[q][q] = new_aqq
    a[p][q] = 0.0
    a[q][p] = 0.0

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
/// Sort eigenvalues descending and reorder eigenvector columns accordingly.
pub fn linalg_sort_eigen(
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

  let sorted_eigenvectors : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = eigenvectors[i][indices[j]]
    }
    sorted_eigenvectors.push(row)
  }

  (sorted_eigenvalues, sorted_eigenvectors)
}

///|
/// Matrix inverse via Gauss-Jordan elimination with partial pivoting.
/// Assumes the matrix is invertible.
pub fn linalg_inverse(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = a.length()
  let m : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(2 * n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = a[i][j]
    }
    row[n + i] = 1.0
    m.push(row)
  }
  for col = 0; col < n; col = col + 1 {
    let mut pivot = col
    let mut best = m[col][col].abs()
    for r = col + 1; r < n; r = r + 1 {
      if m[r][col].abs() > best {
        best = m[r][col].abs()
        pivot = r
      }
    }
    let tmp = m[col]
    m[col] = m[pivot]
    m[pivot] = tmp
    let pv = m[col][col]
    for j = 0; j < 2 * n; j = j + 1 {
      m[col][j] = m[col][j] / pv
    }
    for r = 0; r < n; r = r + 1 {
      if r != col {
        let factor = m[r][col]
        for j = 0; j < 2 * n; j = j + 1 {
          m[r][j] = m[r][j] - factor * m[col][j]
        }
      }
    }
  }
  let inv : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = m[i][n + j]
    }
    inv.push(row)
  }
  inv
}

///|
/// Dense matrix multiply (a: m×k, b: k×n → m×n).
pub fn linalg_mat_mul(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let m = a.length()
  let k = if m > 0 { a[0].length() } else { 0 }
  let n = if b.length() > 0 { b[0].length() } else { 0 }
  let c : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < m; i = i + 1 {
    let row = Array::make(n, 0.0)
    for kk = 0; kk < k; kk = kk + 1 {
      let aik = a[i][kk]
      for j = 0; j < n; j = j + 1 {
        row[j] = row[j] + aik * b[kk][j]
      }
    }
    c.push(row)
  }
  c
}

///|
/// Matrix transpose.
pub fn linalg_transpose(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let m = a.length()
  let n = if m > 0 { a[0].length() } else { 0 }
  let t : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < n; j = j + 1 {
    let row = Array::make(m, 0.0)
    for i = 0; i < m; i = i + 1 {
      row[i] = a[i][j]
    }
    t.push(row)
  }
  t
}