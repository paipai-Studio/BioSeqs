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
  let g = @src.scanpy_neighbors(x, 3)
  let conn = g.to_dense()
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

///|
test "scanpy_neighbors_with_distances" {
  // Two clusters
  let x = [
    [0.0, 0.0],
    [0.0, 0.1],
    [5.0, 5.0],
    [5.0, 5.1],
  ]
  let result = @src.scanpy_neighbors_with_distances(x, 2)
  // Distances graph
  let dist = result.distances.to_dense()
  assert_eq(dist.length(), 4)
  // Diagonal should be 0 (no self-distance in kNN)
  for i = 0; i < 4; i = i + 1 {
    assert_true(dist[i][i] < 1.0e-9)
  }
  // Within-cluster distances should be small and positive
  assert_true(dist[0][1] > 0.0)
  assert_true(dist[0][1] < 1.0)
  assert_true(dist[2][3] > 0.0)
  assert_true(dist[2][3] < 1.0)
  // Cross-cluster distances should be 0 (not kNN neighbors)
  assert_true(dist[0][2] < 1.0e-9)
  // Distances should be symmetric
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 4; j = j + 1 {
      assert_true(assert_close(dist[i][j], dist[j][i], 1.0e-9))
    }
  }
  // Connectivities graph
  let conn = result.connectivities.to_dense()
  assert_eq(conn.length(), 4)
  // Diagonal should be 1
  for i = 0; i < 4; i = i + 1 {
    assert_true(assert_close(conn[i][i], 1.0, 1.0e-9))
  }
  // Within-cluster connectivities positive
  assert_true(conn[0][1] > 0.0)
  assert_true(conn[2][3] > 0.0)
}

///|
test "scanpy_neighbors_gauss_symmetric" {
  let x = [
    [0.0, 0.0],
    [1.0, 0.0],
    [2.0, 0.0],
    [5.0, 0.0],
  ]
  let g = @src.scanpy_neighbors_gauss(x, 3)
  let conn = g.to_dense()
  assert_eq(conn.length(), 4)
  // Diagonal = 1
  for i = 0; i < 4; i = i + 1 {
    assert_true(assert_close(conn[i][i], 1.0, 1.0e-9))
  }
  // Symmetric
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 4; j = j + 1 {
      assert_true(assert_close(conn[i][j], conn[j][i], 1.0e-9))
    }
  }
  // Values in [0, 1]
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 4; j = j + 1 {
      assert_true(conn[i][j] >= -1.0e-12)
      assert_true(conn[i][j] <= 1.0 + 1.0e-12)
    }
  }
  // Nearby points have higher connectivity than distant points
  assert_true(conn[0][1] >= conn[0][3])
}

///|
test "scanpy_neighbors_gauss_vs_umap" {
  // Both methods should produce symmetric matrices with 1 on diagonal,
  // but the actual connectivity values differ.
  let x = [
    [0.0, 0.0],
    [0.5, 0.0],
    [1.0, 0.0],
  ]
  let g_umap = @src.scanpy_neighbors(x, 2)
  let g_gauss = @src.scanpy_neighbors_gauss(x, 2)
  let umap = g_umap.to_dense()
  let gauss = g_gauss.to_dense()
  // Both symmetric
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      assert_true(assert_close(umap[i][j], umap[j][i], 1.0e-9))
      assert_true(assert_close(gauss[i][j], gauss[j][i], 1.0e-9))
    }
  }
  // Both have 1 on diagonal
  for i = 0; i < 3; i = i + 1 {
    assert_true(assert_close(umap[i][i], 1.0, 1.0e-9))
    assert_true(assert_close(gauss[i][i], 1.0, 1.0e-9))
  }
  // Both have positive within-cluster connectivity
  assert_true(umap[0][1] > 0.0)
  assert_true(gauss[0][1] > 0.0)
}

///|
test "scanpy_neighbors_with_distances_empty" {
  let result = @src.scanpy_neighbors_with_distances([], 3)
  assert_eq(result.distances.num_vertices(), 0)
  assert_eq(result.connectivities.num_vertices(), 0)
}

///|
test "scanpy_neighbors_gauss_empty" {
  let g = @src.scanpy_neighbors_gauss([], 3)
  assert_eq(g.num_vertices(), 0)
}
