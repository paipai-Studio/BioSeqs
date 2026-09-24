///|
/// Tests for scanpy `tl.leiden`, `tl.dpt` and `tl.paga`.

///|
test "scanpy_leiden_connected" {
  // two disconnected cliques -> two connected communities
  let a = [
    [0.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 0.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 0.0],
  ]
  let labels = @src.scanpy_leiden(@src.sparse_graph_from_dense(a))
  assert_eq(labels.length(), 6)
  assert_eq((labels[0], labels[1], labels[2]), (labels[0], labels[0], labels[0]))
  assert_eq((labels[3], labels[4], labels[5]), (labels[3], labels[3], labels[3]))
  assert_true(labels[0] != labels[3])
}

///|
test "scanpy_dpt" {
  // a path 0-1-2-3-4 (chain); root = 0, pseudotime should increase away from root
  let conn = [
    [1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 0.0, 0.0],
    [0.0, 1.0, 1.0, 1.0, 0.0],
    [0.0, 0.0, 1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0],
  ]
  let dm = @src.scanpy_diffmap(@src.sparse_graph_from_dense(conn), 5)
  let dpt = @src.scanpy_dpt(dm, 0)
  assert_eq(dpt.length(), 5)
  assert_true(dpt[0] < 1.0e-9) // root has (near) zero pseudotime
  assert_true(dpt[1] < dpt[2])
  assert_true(dpt[2] < dpt[3])
  assert_true(dpt[3] < dpt[4])
}

///|
test "scanpy_paga" {
  // 6 cells in 2 tight clusters -> high within-cluster, low cross-cluster connectivity
  let conn = [
    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
  ]
  let labels = [0, 0, 0, 1, 1, 1]
  let paga = @src.scanpy_paga(@src.sparse_graph_from_dense(conn), labels)
  assert_eq(paga.length(), 2)
  assert_eq(paga[0].length(), 2)
  // diagonal stays 0, off-diagonal is the (low) cross-cluster connectivity
  assert_true(assert_close(paga[0][0], 0.0, 1.0e-12))
  assert_true(paga[0][1] < 0.1) // near-zero connectivity across disconnected clusters
}