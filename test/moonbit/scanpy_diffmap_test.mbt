///|
/// Tests for scanpy `tl.diffmap`.

///|
test "scanpy_diffmap_two_clusters" {
  // two cliques joined by a weak bridge -> Fiedler vector separates clusters
  let conn = [
    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.1, 1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
  ]
  let dm = @src.scanpy_diffmap(conn, 6)
  assert_eq(dm.eigen_basis.length(), 6)
  assert_eq(dm.eigen_values.length(), 6)
  // top eigenvalue ~1 (constant)
  assert_true(assert_close(dm.eigen_values[0], 1.0, 1.0e-6))
  // Fiedler vector (index 1) separates the two clusters by sign
  let v = dm.eigen_basis
  assert_true(v[0][1] * v[4][1] < 0.0)
  assert_true(assert_close(v[0][1], v[1][1], 1.0e-6))
  assert_true(v[3][1] * v[4][1] > 0.0)
}