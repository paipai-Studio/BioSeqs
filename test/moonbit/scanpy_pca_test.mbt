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