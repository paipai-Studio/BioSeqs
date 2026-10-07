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
  let dm = @src.scanpy_diffmap(@src.sparse_graph_from_dense(conn), 6)
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
///|
/// Empty graph (0 vertices) → empty eigen_basis and eigen_values.
test "scanpy_diffmap_empty_graph" {
  let g = @src.SparseGraph::new(0, [], [])
  let dm = @src.scanpy_diffmap(g, 3)
  assert_eq(dm.eigen_basis.length(), 0)
  assert_eq(dm.eigen_values.length(), 0)
}

///|
/// Single vertex with self-weight 1.0: the transition matrix is [1.0],
/// so the sole eigenvalue is 1.0 and the eigenvector is [1.0].
test "scanpy_diffmap_single_vertex" {
  let g = @src.SparseGraph::new(1, [1.0], [[]])
  let dm = @src.scanpy_diffmap(g, 1)
  assert_eq(dm.eigen_basis.length(), 1)
  assert_eq(dm.eigen_values.length(), 1)
  // top eigenvalue of a 1×1 transition matrix with self-weight 1 is 1.0
  assert_true(assert_close(dm.eigen_values[0], 1.0, 1.0e-6))
}

///|
/// Two isolated vertices (no edges, only diagonal self-weights):
/// the transition matrix is the identity. The dominant eigenvalue is 1.0.
/// Note: Lanczos can only resolve one copy of a degenerate eigenvalue,
/// so the second value may be 0 — we only assert the dominant one.
test "scanpy_diffmap_two_isolated_vertices" {
  let g = @src.SparseGraph::new(2, [1.0, 1.0], [[], []])
  let dm = @src.scanpy_diffmap(g, 2)
  assert_eq(dm.eigen_basis.length(), 2)
  assert_eq(dm.eigen_values.length(), 2)
  // dominant eigenvalue of the identity is 1.0
  assert_true(assert_close(dm.eigen_values[0], 1.0, 1.0e-6))
}
