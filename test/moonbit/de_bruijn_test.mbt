///|
/// Tests for De Bruijn Graph.

///|
test "de_bruijn_create" {
  let graph = @src.DeBruijnGraph::new(3)
  assert_eq(graph.k(), 3)
  assert_eq(graph.node_count(), 0)
  assert_eq(graph.edge_count(), 0)
}

///|
test "de_bruijn_add_kmer" {
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_kmer("ABC")
  assert_eq(graph.node_count(), 2)
  assert_eq(graph.edge_count(), 1)
  assert_true(graph.contains_node("AB"))
  assert_true(graph.contains_node("BC"))
  assert_eq(graph.in_degree("AB"), 0)
  assert_eq(graph.out_degree("AB"), 1)
  assert_eq(graph.in_degree("BC"), 1)
  assert_eq(graph.out_degree("BC"), 0)
}

///|
test "de_bruijn_add_sequence" {
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_sequence("ABCDE")
  assert_eq(graph.node_count(), 4)
  assert_eq(graph.edge_count(), 3)
}

///|
test "de_bruijn_start_end_nodes" {
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_sequence("ABCDE")
  let starts = graph.start_nodes()
  let ends = graph.end_nodes()
  assert_eq(starts.length(), 1)
  assert_eq(starts[0], "AB")
  assert_eq(ends.length(), 1)
  assert_eq(ends[0], "DE")
}

///|
test "de_bruijn_find_eulerian_path_simple" {
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_sequence("ABCDE")
  let path = graph.find_eulerian_path()
  match path {
    Some(seq) => assert_eq(seq, "ABCDE")
    None => assert_true(false)
  }
}

///|
test "de_bruijn_find_eulerian_path_circular" {
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_sequence("ABCA")
  let path = graph.find_eulerian_path()
  assert_true(path is Some(_))
}

///|
test "de_bruijn_assemble_single" {
  let result = @src.assemble_single("ABCDEFG", 4)
  match result {
    Some(seq) => assert_eq(seq, "ABCDEFG")
    None => assert_true(false)
  }
}

///|
test "de_bruijn_assemble_sequences" {
  let reads = ["ABCDE", "CDEFG", "EFGHI"]
  let result = @src.assemble_sequences(reads, 3)
  match result {
    Some(seq) => assert_eq(seq, "ABCDEFGHI")
    None => assert_true(false)
  }
}

///|
test "de_bruijn_branch_nodes" {
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_sequence("ABCDE")
  graph.add_sequence("ABXDE")
  let branches = graph.branch_nodes()
  assert_eq(branches.length(), 1)
  assert_eq(branches[0], "AB")
}

///|
test "de_bruijn_outgoing" {
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_sequence("ABCDE")
  let outgoing = graph.outgoing("BC")
  assert_eq(outgoing.length(), 1)
  assert_eq(outgoing[0], "CD")
}

///|
test "de_bruijn_all_nodes" {
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_sequence("ABC")
  let nodes = graph.all_nodes()
  assert_eq(nodes.length(), 2)
}

///|
test "de_bruijn_to_string" {
  let graph = @src.DeBruijnGraph::new(3)
  graph.add_sequence("ABC")
  let str = graph.to_string()
  assert_true(str.length() > 0)
}
