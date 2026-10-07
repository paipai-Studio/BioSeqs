///|
/// Tests for scanpy `pp.pca`. Property-based checks (orthonormal loadings,
/// centred reconstruction, variance-ratio ordering) that hold for any input.

///|
test "scanpy_pca" {
  let x = [[1.0, 2.0, 3.0], [2.0, 5.0, 1.0], [0.0, 1.0, 4.0], [3.0, 2.0, 0.0]]
  let pca = @src.scanpy_pca(x, 3)
  assert_eq(pca.coords.length(), 4)
  assert_eq(pca.coords[0].length(), 3)
  assert_eq(pca.loadings.length(), 3)
  assert_eq(pca.loadings[0].length(), 3)

  // variance ratio sums to 1 and is non-increasing
  let r = pca.variance_ratio
  let mut s = 0.0
  for v in r {
    s = s + v
  }
  assert_true(assert_close(s, 1.0, 1.0e-6))
  assert_true(r[0] >= r[1] - 1.0e-9)
  assert_true(r[1] >= r[2] - 1.0e-9)

  // loadings are orthonormal rows
  for a = 0; a < 3; a = a + 1 {
    for b = 0; b < 3; b = b + 1 {
      let mut dot = 0.0
      for j = 0; j < 3; j = j + 1 {
        dot = dot + pca.loadings[a][j] * pca.loadings[b][j]
      }
      assert_true(assert_close(dot, if a == b { 1.0 } else { 0.0 }, 1.0e-6))
    }
  }

  // reconstruction: coords @ loadings == centred X
  let means = [1.5, 2.5, 2.0]
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      let mut re = 0.0
      for c = 0; c < 3; c = c + 1 {
        re = re + pca.coords[i][c] * pca.loadings[c][j]
      }
      assert_true(assert_close(re, x[i][j] - means[j], 1.0e-6))
    }
  }
}
///|
/// Empty matrix (0 rows) → all output arrays empty.
test "scanpy_pca_empty_matrix" {
  let x : Array[Array[Double]] = []
  let pca = @src.scanpy_pca(x, 3)
  assert_eq(pca.coords.length(), 0)
  assert_eq(pca.loadings.length(), 0)
  assert_eq(pca.variance_ratio.length(), 0)
}

///|
/// Single sample (1 row): after column centring every entry is 0, so the
/// covariance is the zero matrix and all variance ratios are 0.
test "scanpy_pca_single_sample" {
  let x = [[1.0, 5.0, 3.0]]
  let pca = @src.scanpy_pca(x, 2)
  // coords has 1 row, but n_comps clamped to min(2, 3) = 2
  assert_eq(pca.coords.length(), 1)
  assert_eq(pca.coords[0].length(), 2)
  // with zero covariance, all coords are 0 and variance_ratio is 0
  for c = 0; c < 2; c = c + 1 {
    assert_true(assert_close(pca.coords[0][c], 0.0, 1.0e-10))
    assert_true(assert_close(pca.variance_ratio[c], 0.0, 1.0e-10))
  }
}

///|
/// n_comps > n_vars: the function should silently clamp to n_vars
/// rather than indexing out of bounds.
test "scanpy_pca_n_comps_exceeds_n_vars" {
  let x = [[1.0, 2.0], [3.0, 1.0], [0.0, 5.0]]
  let pca = @src.scanpy_pca(x, 10) // n_comps=10, n_vars=2
  // clamped to n_vars = 2
  assert_eq(pca.loadings.length(), 2)
  assert_eq(pca.loadings[0].length(), 2)
  assert_eq(pca.coords.length(), 3)
  assert_eq(pca.coords[0].length(), 2)
  assert_eq(pca.variance_ratio.length(), 2)
  // variance ratio still sums to 1 (two non-degenerate components)
  let mut s = 0.0
  for v in pca.variance_ratio {
    s = s + v
  }
  assert_true(assert_close(s, 1.0, 1.0e-6))
}
