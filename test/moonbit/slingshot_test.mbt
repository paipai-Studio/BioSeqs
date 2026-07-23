///|
test "slingshot_node_create" {
  let coords = [1.0, 2.0, 3.0]
  let node = @src.SlingshotNode::new("c1", "cluster1", coords, 10)
  assert_eq(node.id, "c1")
  assert_eq(node.cluster, "cluster1")
  assert_eq(node.size, 10)
  assert_eq(node.sling_coord(0), 1.0)
  assert_eq(node.sling_coord(1), 2.0)
}

///|
test "slingshot_node_terminal" {
  let coords = [1.0, 2.0]
  let node = @src.SlingshotNode::new("c1", "cluster1", coords, 5)
  let terminal = node.sling_set_terminal(true)
  assert_true(terminal.terminal)
  let start_node = node.sling_set_start(true)
  assert_true(start_node.start)
}

///|
test "slingshot_edge_create" {
  let edge = @src.SlingshotEdge::new("c1", "c2", 1.5)
  assert_eq(edge.from_id, "c1")
  assert_eq(edge.to_id, "c2")
  assert_eq(edge.weight, 1.5)
}

///|
test "slingshot_distance" {
  let node1 = @src.SlingshotNode::new("c1", "cluster1", [0.0, 0.0], 5)
  let node2 = @src.SlingshotNode::new("c2", "cluster2", [3.0, 4.0], 5)
  let dist = node1.sling_distance(node2, @src.sling_euclidean_metric())
  assert_true((dist - 5.0).abs() < 1.0e-6)
}

///|
test "slingshot_distance_matrix" {
  let nodes = [
    @src.SlingshotNode::new("c1", "cluster1", [0.0, 0.0], 5),
    @src.SlingshotNode::new("c2", "cluster2", [1.0, 0.0], 5),
    @src.SlingshotNode::new("c3", "cluster3", [0.0, 1.0], 5)
  ]
  let matrix = @src.sling_distance_matrix(nodes, @src.sling_euclidean_metric())
  assert_eq(matrix.length(), 3)
  assert_eq(matrix[0].length(), 3)
  assert_true(matrix[0][0].abs() < 1.0e-10)
  assert_true(matrix[1][1].abs() < 1.0e-10)
}

///|
test "slingshot_mst" {
  let nodes = [
    @src.SlingshotNode::new("c1", "cluster1", [0.0, 0.0], 5),
    @src.SlingshotNode::new("c2", "cluster2", [1.0, 0.0], 5),
    @src.SlingshotNode::new("c3", "cluster3", [0.5, 1.0], 5)
  ]
  let edges = @src.sling_build_mst(nodes, @src.sling_euclidean_metric())
  assert_eq(edges.length(), 2)
}

///|
test "slingshot_terminals" {
  let nodes = [
    @src.SlingshotNode::new("c1", "cluster1", [0.0, 0.0], 5),
    @src.SlingshotNode::new("c2", "cluster2", [1.0, 0.0], 5),
    @src.SlingshotNode::new("c3", "cluster3", [0.5, 1.0], 5)
  ]
  let edges = @src.sling_build_mst(nodes, @src.sling_euclidean_metric())
  let terminals = @src.sling_find_terminals(edges, nodes)
  assert_true(terminals.length() >= 2)
}

///|
test "slingshot_principal_curve" {
  let points : Array[Array[Double]] = [
    [0.0, 0.0],
    [1.0, 0.5],
    [2.0, 1.0],
    [3.0, 1.5],
    [4.0, 2.0]
  ]
  let curve = @src.sling_fit_principal_curve(points, 5)
  assert_true(curve.length() >= 2)
}

///|
test "slingshot_pseudotime" {
  let cell_coords : Array[Array[Double]] = [
    [0.0, 0.0],
    [1.0, 0.5],
    [2.0, 1.0],
    [3.0, 1.5]
  ]
  let curve : Array[Array[Double]] = [
    [0.0, 0.0],
    [1.0, 0.5],
    [2.0, 1.0],
    [3.0, 1.5]
  ]
  let pseudotime = @src.sling_compute_pseudotime(cell_coords, curve)
  assert_eq(pseudotime.length(), 4)
}

///|
test "slingshot_full_run" {
  let cell_coords : Array[Array[Double]] = [
    [0.0, 0.0],
    [0.5, 0.2],
    [1.0, 0.5],
    [2.0, 1.0],
    [2.5, 1.2],
    [3.0, 1.5]
  ]
  let cluster_labels = ["c1", "c1", "c2", "c2", "c3", "c3"]
  let params = @src.SlingshotParams::new()
  let result = @src.slingshot_run(cell_coords, cluster_labels, params)
  assert_true(result.cell_pseudotime.length() > 0)
  assert_true(result.curves.length() > 0)
}
