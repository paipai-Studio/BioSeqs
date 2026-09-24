///|
/// Tests for the Lanczos top-k eigensolver (`linalg_lanczos_topk`): its top-k
/// eigenvalues must match a full Jacobi sweep on a well-separated spectrum, and
/// the Ritz vectors must be (approximate) eigenvectors.

///|
test "linalg_lanczos_topk_matches_jacobi" {
  let n = 40
  // 对角降序 + 小耦合的三对角对称矩阵, 特征值良好分离
  let a : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    row[i] = (n - i).to_double() + 1.0
    a.push(row)
  }
  for i = 0; i + 1 < n; i = i + 1 {
    a[i][i + 1] = 0.5
    a[i + 1][i] = 0.5
  }
  let mat_vec = fn(v : Array[Double]) -> Array[Double] {
    let out = Array::make(n, 0.0)
    for i = 0; i < n; i = i + 1 {
      let mut s = 0.0
      for j = 0; j < n; j = j + 1 {
        s = s + a[i][j] * v[j]
      }
      out[i] = s
    }
    out
  }

  // Jacobi 全谱(参照)
  let (je, jv) = @src.linalg_jacobi(a, n, 100 * n * n)
  let (se, _) = @src.linalg_sort_eigen(je, jv, n)

  // Lanczos top-3（n_iter = 3+15 = 18 < 40, 真正走近似）
  let (le, lv) = @src.linalg_lanczos_topk(mat_vec, n, 3)
  assert_eq(le.length(), 3)
  assert_eq(lv.length(), n)

  // 特征值一致
  for t = 0; t < 3; t = t + 1 {
    assert_true(assert_close(le[t], se[t], 1.0e-2))
  }

  // Ritz 向量残差 ||A·x - λ·x|| 小(符号无关的正确性检查)
  for t = 0; t < 3; t = t + 1 {
    let x : Array[Double] = Array::make(n, 0.0)
    for i = 0; i < n; i = i + 1 {
      x[i] = lv[i][t]
    }
    let ax = mat_vec(x)
    let mut norm = 0.0
    for i = 0; i < n; i = i + 1 {
      let diff = ax[i] - le[t] * x[i]
      norm = norm + diff * diff
    }
    assert_true(norm.sqrt() < 1.0e-2)
  }
}