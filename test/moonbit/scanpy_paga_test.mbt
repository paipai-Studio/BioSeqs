///|
/// Tests for scanpy PAGA (partition-based graph abstraction).
///
/// PAGA computes connectivity between clusters based on the kNN graph.
/// Given a SparseGraph and cluster labels, it returns a cluster-level
/// connectivity matrix normalized by a random-graph null model.

///|
/// Build a simple SparseGraph: path graph 0-1-2-3-4 with unit weights.
fn paga_path_graph() -> @src.SparseGraph {
  // 0 -- 1 -- 2 -- 3 -- 4
  let n = 5
  let diag = [1.0, 1.0, 1.0, 1.0, 1.0]
  let adj : Array[Array[(Int, Double)]] = [
    [(1, 1.0)],
    [(0, 1.0), (2, 1.0)],
    [(1, 1.0), (3, 1.0)],
    [(2, 1.0), (4, 1.0)],
    [(3, 1.0)],
  ]
  @src.SparseGraph::new(n, diag, adj)
}

///|
test "paga_two_clusters_connected" {
  // Graph: 0-1-2-3-4, labels: [0,0,0,1,1]
  // Cluster 0 = {0,1,2}, cluster 1 = {3,4}
  // Edge 2-3 connects the two clusters.
  let g = paga_path_graph()
  let labels = [0, 0, 0, 1, 1]
  let conn = @src.scanpy_paga(g, labels)
  assert_eq(conn.length(), 2)
  // Off-diagonal should be > 0 (clusters are connected)
  assert_true(conn[0][1] > 0.0)
  assert_true(conn[1][0] > 0.0)
  // Symmetric
  assert_eq(conn[0][1], conn[1][0])
  // Diagonal is 0 (PAGA skips i==j)
  assert_eq(conn[0][0], 0.0)
  assert_eq(conn[1][1], 0.0)
}

///|
test "paga_three_clusters_chain" {
  // Graph: 0-1-2-3-4-5, labels: [0,0,1,1,2,2]
  // Cluster 0 = {0,1}, cluster 1 = {2,3}, cluster 2 = {4,5}
  // Edges 1-2 and 3-4 connect the clusters in a chain.
  let n = 6
  let diag = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  let adj : Array[Array[(Int, Double)]] = [
    [(1, 1.0)],
    [(0, 1.0), (2, 1.0)],
    [(1, 1.0), (3, 1.0)],
    [(2, 1.0), (4, 1.0)],
    [(3, 1.0), (5, 1.0)],
    [(4, 1.0)],
  ]
  let g = @src.SparseGraph::new(n, diag, adj)
  let labels = [0, 0, 1, 1, 2, 2]
  let conn = @src.scanpy_paga(g, labels)
  assert_eq(conn.length(), 3)
  // 0-1 connected, 1-2 connected, 0-2 not directly connected
  assert_true(conn[0][1] > 0.0)
  assert_true(conn[1][2] > 0.0)
  // Symmetric
  assert_eq(conn[0][1], conn[1][0])
  assert_eq(conn[1][2], conn[2][1])
  assert_eq(conn[0][2], conn[2][0])
}

///|
test "paga_disconnected_clusters" {
  // Graph: 0-1, 2-3 (two disconnected components)
  let n = 4
  let diag = [1.0, 1.0, 1.0, 1.0]
  let adj : Array[Array[(Int, Double)]] = [
    [(1, 1.0)],
    [(0, 1.0)],
    [(3, 1.0)],
    [(2, 1.0)],
  ]
  let g = @src.SparseGraph::new(n, diag, adj)
  let labels = [0, 0, 1, 1]
  let conn = @src.scanpy_paga(g, labels)
  assert_eq(conn.length(), 2)
  // No inter-cluster edges → connectivity should be 0
  assert_eq(conn[0][1], 0.0)
  assert_eq(conn[1][0], 0.0)
}

///|
test "paga_single_cluster" {
  // All nodes in one cluster → only 1×1 matrix with 0
  let g = paga_path_graph()
  let labels = [0, 0, 0, 0, 0]
  let conn = @src.scanpy_paga(g, labels)
  assert_eq(conn.length(), 1)
  assert_eq(conn[0][0], 0.0)
}

///|
test "paga_empty_graph" {
  let g = @src.SparseGraph::new(0, [], [])
  let conn = @src.scanpy_paga(g, [])
  assert_eq(conn.length(), 0)
}

///|
test "paga_symmetry" {
  // Random-ish graph with 3 clusters, verify conn matrix is symmetric
  let n = 6
  let diag = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  let adj : Array[Array[(Int, Double)]] = [
    [(1, 0.5), (2, 0.3), (3, 0.8)],
    [(0, 0.5), (2, 0.7), (4, 0.2)],
    [(0, 0.3), (1, 0.7), (5, 0.6)],
    [(0, 0.8), (4, 0.4), (5, 0.9)],
    [(1, 0.2), (3, 0.4), (5, 0.3)],
    [(2, 0.6), (3, 0.9), (4, 0.3)],
  ]
  let g = @src.SparseGraph::new(n, diag, adj)
  let labels = [0, 0, 1, 1, 2, 2]
  let conn = @src.scanpy_paga(g, labels)
  assert_eq(conn.length(), 3)
  // Verify symmetry: conn[i][j] == conn[j][i] for all i,j
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      let diff = conn[i][j] - conn[j][i]
      assert_true(diff < 0.0001)
      assert_true(diff > -0.0001)
    }
  }
}

///|
test "paga_clamped_to_unit" {
  // With very high weights, connectivity should be clamped to [0, 1]
  let n = 4
  let diag = [1.0, 1.0, 1.0, 1.0]
  let adj : Array[Array[(Int, Double)]] = [
    [(1, 100.0), (2, 100.0)],
    [(0, 100.0), (2, 100.0)],
    [(0, 100.0), (1, 100.0), (3, 100.0)],
    [(2, 100.0)],
  ]
  let g = @src.SparseGraph::new(n, diag, adj)
  let labels = [0, 0, 0, 1]
  let conn = @src.scanpy_paga(g, labels)
  assert_eq(conn.length(), 2)
  // Connectivity should be at most 1.0
  assert_true(conn[0][1] <= 1.0)
  assert_true(conn[1][0] <= 1.0)
  assert_true(conn[0][1] >= 0.0)
}

///|
test "paga_each_node_own_cluster" {
  // Each node is its own cluster → PAGA reduces to the original graph
  let g = paga_path_graph()
  let labels = [0, 1, 2, 3, 4]
  let conn = @src.scanpy_paga(g, labels)
  assert_eq(conn.length(), 5)
  // Adjacent nodes should have connectivity > 0
  assert_true(conn[0][1] > 0.0)
  assert_true(conn[1][2] > 0.0)
  assert_true(conn[2][3] > 0.0)
  assert_true(conn[3][4] > 0.0)
  // Non-adjacent nodes should have connectivity 0
  assert_eq(conn[0][2], 0.0)
  assert_eq(conn[0][3], 0.0)
  assert_eq(conn[0][4], 0.0)
}
