///|
/// Tests for scanpy `pp.neighbors` (kNN + UMAP connectivities).

///|
test "scanpy_neighbors" {
  // two clear clusters: 3 points near origin, 3 near (5,5)
  let x = [
    [0.0, 0.0],
    [0.0, 0.1],
    [0.0, 0.2],
    [5.0, 5.0],
    [5.0, 5.1],
    [5.0, 5.2],
  ]
  let conn = @src.scanpy_neighbors(x, 3)
  assert_eq(conn.length(), 6)
  assert_eq(conn[0].length(), 6)

  // symmetric + diagonal 1 + values in [0,1]
  for i = 0; i < 6; i = i + 1 {
    assert_true(assert_close(conn[i][i], 1.0, 1.0e-9))
    for j = 0; j < 6; j = j + 1 {
      assert_true(assert_close(conn[i][j], conn[j][i], 1.0e-9))
      assert_true(conn[i][j] >= -1.0e-12)
      assert_true(conn[i][j] <= 1.0 + 1.0e-12)
    }
  }

  // within-cluster connectivities are positive, cross-cluster are ~0
  assert_true(conn[0][1] > 0.0)
  assert_true(conn[3][4] > 0.0)
  assert_true(conn[0][3] < 1.0e-9)
  assert_true(conn[1][5] < 1.0e-9)
}