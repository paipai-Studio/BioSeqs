///|
test "graph_create_example" {
  let graph = @src.create_example_graph()
  assert_eq(@src.graph_node_count(graph), 8)
  assert_eq(@src.graph_edge_count(graph), 11)
}

///|
test "graph_add_node_and_edge" {
  let graph = @src.Graph::new("Test Graph", @src.GraphType::new_undirected())
  @src.graph_add_node(graph, @src.Node::new("A", "Node A"))
  @src.graph_add_node(graph, @src.Node::new("B", "Node B"))
  @src.graph_add_edge_by_id(graph, "A", "B", 1.0)
  assert_eq(@src.graph_node_count(graph), 2)
  assert_eq(@src.graph_edge_count(graph), 1)
}

///|
test "graph_get_neighbors" {
  let graph = @src.create_example_graph()
  let neighbors = @src.graph_get_neighbors(graph, "TP53")
  assert_eq(neighbors.length(), 5)
}

///|
test "graph_get_degree" {
  let graph = @src.create_example_graph()
  let degree = @src.graph_get_degree(graph, "TP53")
  assert_eq(degree, 5)
}

///|
test "graph_has_edge" {
  let graph = @src.create_example_graph()
  assert_true(@src.graph_has_edge(graph, "TP53", "BRCA1"))
  assert_true(@src.graph_has_edge(graph, "BRCA1", "TP53"))
}

///|
test "graph_shortest_path" {
  let graph = @src.create_example_graph()
  let path = @src.graph_shortest_path(graph, "TP53", "CHEK2")
  assert_true(path.is_some())
  if path.is_some() {
    assert_true(path.unwrap().length() >= 2)
  }
}

///|
test "graph_connected_components" {
  let graph = @src.create_example_graph()
  let components = @src.graph_connected_components(graph)
  assert_eq(components.length(), 1)
}

///|
test "graph_summarize" {
  let graph = @src.create_example_graph()
  let summary = @src.graph_summarize(graph)
  assert_true(@src.contains(summary, "Nodes"))
}

///|
test "graph_to_dot" {
  let graph = @src.create_example_graph()
  let dot = @src.graph_to_dot(graph)
  assert_true(@src.contains(dot, "graph"))
}

///|
test "graph_directed_pathway" {
  let graph = @src.create_pathway_graph()
  assert_eq(@src.graph_node_count(graph), 9)
  let path = @src.graph_shortest_path(graph, "G1", "M")
  assert_true(path.is_some())
}
