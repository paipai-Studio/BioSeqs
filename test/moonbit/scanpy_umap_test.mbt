///|
/// Tests for scanpy `tl.umap` (2D embedding of a connectivity graph).

///|
test "scanpy_umap_separates_clusters" {
  // two cliques with a weak bridge -> UMAP should separate the two clusters
  let conn = [
    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.1, 1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
  ]
  let emb = @src.scanpy_umap(@src.sparse_graph_from_dense(conn), 200)
  assert_eq(emb.length(), 6)
  assert_eq(emb[0].length(), 2)
  // cluster centroid distance should be much larger than within-cluster spread
  let c1 = [(emb[0][0] + emb[1][0] + emb[2][0]) / 3.0, (emb[0][1] + emb[1][1] + emb[2][1]) / 3.0]
  let c2 = [(emb[3][0] + emb[4][0] + emb[5][0]) / 3.0, (emb[3][1] + emb[4][1] + emb[5][1]) / 3.0]
  let d_c = (c1[0] - c2[0]) * (c1[0] - c2[0]) + (c1[1] - c2[1]) * (c1[1] - c2[1])
  let d_1 = (emb[0][0] - emb[1][0]) * (emb[0][0] - emb[1][0]) + (emb[0][1] - emb[1][1]) * (emb[0][1] - emb[1][1])
  let d_2 = (emb[3][0] - emb[4][0]) * (emb[3][0] - emb[4][0]) + (emb[3][1] - emb[4][1]) * (emb[3][1] - emb[4][1])
  let within = if d_1 > d_2 { d_1 } else { d_2 }
  assert_true(d_c > within + 0.5)
}
///|
/// Empty graph (0 vertices) → empty embedding.
test "scanpy_umap_empty_graph" {
  let g = @src.SparseGraph::new(0, [], [])
  let emb = @src.scanpy_umap(g, 100)
  assert_eq(emb.length(), 0)
}

///|
/// Single vertex (n < 3): returns the trivial layout [[0.0, 0.0]].
test "scanpy_umap_single_vertex" {
  let g = @src.SparseGraph::new(1, [1.0], [[]])
  let emb = @src.scanpy_umap(g, 100)
  assert_eq(emb.length(), 1)
  assert_eq(emb[0].length(), 2)
  assert_true(assert_close(emb[0][0], 0.0, 1.0e-10))
  assert_true(assert_close(emb[0][1], 0.0, 1.0e-10))
}

///|
/// Two vertices (n < 3): returns the trivial layout [[0,0], [1,0]].
/// The function short-circuits for n < 3 without running SGD.
test "scanpy_umap_two_vertices" {
  let g = @src.SparseGraph::new(2, [1.0, 1.0], [[(1, 0.5)], [(0, 0.5)]])
  let emb = @src.scanpy_umap(g, 100)
  assert_eq(emb.length(), 2)
  assert_eq(emb[0].length(), 2)
  // trivial layout: vertex 0 at (0,0), vertex 1 at (1,0)
  assert_true(assert_close(emb[0][0], 0.0, 1.0e-10))
  assert_true(assert_close(emb[0][1], 0.0, 1.0e-10))
  assert_true(assert_close(emb[1][0], 1.0, 1.0e-10))
  assert_true(assert_close(emb[1][1], 0.0, 1.0e-10))
}
