///|
/// Lanczos iteration for the top eigenvalues/eigenvectors of a symmetric
/// linear operator `A` (given as a matrix-vector closure). Runs `O(k·n²)` for a
/// dense operator / `O(k·nnz)` for a sparse one, instead of the `O(n³)` of the
/// full Jacobi sweep.

///|
pub fn linalg_lanczos_topk(
  mat_vec : (Array[Double]) -> Array[Double],
  n : Int,
  n_comps : Int,
) -> (Array[Double], Array[Array[Double]]) {
  if n == 0 || n_comps <= 0 {
    return ([], [])
  }
  // a few extra iterations so the top `n_comps` Ritz values converge
  let n_iter = if n_comps + 15 < n { n_comps + 15 } else { n }
  let want = if n_comps < n { n_comps } else { n }

  // 确定性起始向量: 第一个基向量 e0。非退化(与常数/Fiedler 都不正交),
  // 能张成一般矩阵的完整 Krylov 子空间。
  let mut q : Array[Double] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    q.push(if i == 0 { 1.0 } else { 0.0 })
  }

  let v : Array[Array[Double]] = Array(capacity=n_iter)
  let alpha : Array[Double] = Array(capacity=n_iter)
  let beta : Array[Double] = Array(capacity=n_iter)

  let mut prev_beta = 0.0
  let mut m = 0
  for j = 0; j < n_iter; j = j + 1 {
    v.push(q.copy())
    let aq = mat_vec(q)
    // r = A q - beta_{j-1} v_{j-1} - alpha_j q
    let r : Array[Double] = Array(capacity=n)
    let mut a_j = 0.0
    for i = 0; i < n; i = i + 1 {
      let val = aq[i] - (if j > 0 { prev_beta * v[j - 1][i] } else { 0.0 })
      a_j = a_j + q[i] * val
      r.push(val)
    }
    for i = 0; i < n; i = i + 1 {
      r[i] = r[i] - a_j * q[i]
    }
    // full re-orthogonalization against all previous Lanczos vectors
    for t = 0; t < j; t = t + 1 {
      let vv = v[t]
      let mut dot = 0.0
      for i = 0; i < n; i = i + 1 {
        dot = dot + vv[i] * r[i]
      }
      for i = 0; i < n; i = i + 1 {
        r[i] = r[i] - dot * vv[i]
      }
    }
    let mut nr = 0.0
    for i = 0; i < n; i = i + 1 {
      nr = nr + r[i] * r[i]
    }
    let bj = nr.sqrt()
    alpha.push(a_j)
    m = j + 1
    if bj < 1.0e-12 {
      beta.push(0.0)
      break
    }
    beta.push(bj)
    let nq : Array[Double] = Array(capacity=n)
    for i = 0; i < n; i = i + 1 {
      nq.push(r[i] / bj)
    }
    q = nq
    prev_beta = bj
  }

  // tridiagonal T (m × m)
  let tmat : Array[Array[Double]] = Array(capacity=m)
  for i = 0; i < m; i = i + 1 {
    let row = Array::make(m, 0.0)
    row[i] = alpha[i]
    tmat.push(row)
  }
  for i = 0; i + 1 < m; i = i + 1 {
    tmat[i][i + 1] = beta[i]
    tmat[i + 1][i] = beta[i]
  }

  // Ritz values/vectors from the small T
  let (te, tv) = linalg_jacobi(tmat, m, 100 * m * m)
  let (se, sv) = linalg_sort_eigen(te, tv, m)

  // 若 Krylov 子空间提前耗尽(m < want), 剩余特征值补 0(秩亏), 特征向量补零
  let eigen_values : Array[Double] = Array(capacity=want)
  for t = 0; t < want; t = t + 1 {
    eigen_values.push(if t < m { se[t] } else { 0.0 })
  }
  // eigen_vectors[i][t] = Σ_j sv[j][t] · v[j][i]
  let eigen_vectors : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=want)
    for t = 0; t < want; t = t + 1 {
      if t < m {
        let mut s = 0.0
        for j = 0; j < m; j = j + 1 {
          s = s + sv[j][t] * v[j][i]
        }
        row.push(s)
      } else {
        row.push(0.0)
      }
    }
    eigen_vectors.push(row)
  }
  (eigen_values, eigen_vectors)
}