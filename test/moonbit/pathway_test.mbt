/// Test suite for Bio.Pathway (Reaction / System / Network) and
/// Bio.Pathway.Rep (PathwayGraph / PathwayMultiGraph + BFS / DFS).
///
/// All expected values were frozen from Biopython 1.88 + Python 3.12.

// ---- helpers ----------------------------------------------------------

fn _pathway_double_eq(a : Double, b : Double) -> Bool {
  (a - b).abs() < 0.000000000001
}

fn _pathway_sorted(arr : Array[String]) -> Array[String] {
  let copy = arr.copy()
  copy.sort_by(fn(a : String, b : String) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  copy
}

fn _pathway_sorted_strings_equal(a : Array[String], b : Array[String]) -> Bool {
  let sa = _pathway_sorted(a)
  let sb = _pathway_sorted(b)
  if sa.length() != sb.length() { return false }
  let mut i = 0
  while i < sa.length() {
    if sa[i] != sb[i] { return false }
    i = i + 1
  }
  true
}

// =======================================================================
// 1. Reaction  — coefficients are signed, negative = substrate,
//                positive = product.
// =======================================================================

test "reaction_coef_basics" {
  let r1 = @src.Reaction::new("r1", Map([("S", -1.0), ("A", 1.0)]))
  assert_true(_pathway_double_eq(r1.reactants["S"], -1.0))
  assert_true(_pathway_double_eq(r1.reactants["A"], 1.0))
  assert_false(r1.reversible)
}

test "reaction_reverse_flips_signs" {
  let r1 = @src.Reaction::new("r1", Map([("S", -1.0), ("A", 1.0)]))
  let rr = r1.reverse()
  assert_true(_pathway_double_eq(rr.reactants["S"], 1.0))
  assert_true(_pathway_double_eq(rr.reactants["A"], -1.0))
}

test "reaction_reversible_flag" {
  let r = @src.Reaction::new(
    "r5",
    Map([("X", -2.0), ("Y", -1.0), ("Z", 3.0), ("W", 1.0)]),
    reversible=true,
  )
  assert_true(r.reversible)
  assert_true(_pathway_double_eq(r.reactants["X"], -2.0))
  assert_true(_pathway_double_eq(r.reactants["Z"], 3.0))
}

test "reaction_zeros_dropped" {
  let r = @src.Reaction::new(
    "rz",
    Map([("A", -1.0), ("B", 0.0), ("C", 2.0)]),
  )
  assert_true(_pathway_sorted_strings_equal(r.species(), ["A", "C"]))
}

test "reaction_species" {
  let r = @src.Reaction::new("rx", Map([("S", -1.0), ("A", 1.0)]))
  assert_true(_pathway_sorted_strings_equal(r.species(), ["A", "S"]))
}

test "reaction_to_string_irreversible" {
  let r = @src.Reaction::new("t", Map([("a", -2.0), ("b", 1.0)]))
  let s = r.to_string()
  assert_true(s.contains("-->"))
  assert_true(s.contains("a"))
  assert_true(s.contains("b"))
}

test "reaction_to_string_reversible" {
  let r = @src.Reaction::new("t", Map([("x", -1.0), ("y", 1.0)]), reversible=true)
  let s = r.to_string()
  assert_true(s.contains("<->"))
}

// =======================================================================
// 2. System + stoichiometry matrix  (species × reactions)
// =======================================================================

test "stoichiometry_simple_3reactions" {
  let r1 = @src.Reaction::new("r1", Map([("S", -1.0), ("A", 1.0)]))
  let r2 = @src.Reaction::new("r2", Map([("A", -1.0), ("B", 1.0)]))
  let r3 = @src.Reaction::new("r3", Map([("B", -1.0), ("P", 1.0)]))
  let sys = @src.System::new("linear", [r1, r2, r3])
  let (species, rxns, stoch) = sys.stochiometry()

  assert_true(_pathway_sorted_strings_equal(species, ["A", "B", "P", "S"]))
  assert_eq(rxns.length(), 3)

  let mut ri = 0
  while ri < species.length() {
    if species[ri] == "S" { break }
    ri = ri + 1
  }
  assert_true(_pathway_double_eq(stoch[ri][0], -1.0))
  assert_true(_pathway_double_eq(stoch[ri][1], 0.0))
  assert_true(_pathway_double_eq(stoch[ri][2], 0.0))

  let mut ri2 = 0
  while ri2 < species.length() {
    if species[ri2] == "A" { break }
    ri2 = ri2 + 1
  }
  assert_true(_pathway_double_eq(stoch[ri2][0], 1.0))
  assert_true(_pathway_double_eq(stoch[ri2][1], -1.0))
  assert_true(_pathway_double_eq(stoch[ri2][2], 0.0))
}

test "stoichiometry_single_reaction" {
  let r5 = @src.Reaction::new(
    "r5",
    Map([("X", -2.0), ("Y", -1.0), ("Z", 3.0), ("W", 1.0)]),
    reversible=true,
  )
  let sys = @src.System::new("single", [r5])
  let (species, rxns, stoch) = sys.stochiometry()
  assert_true(_pathway_sorted_strings_equal(species, ["W", "X", "Y", "Z"]))
  assert_eq(rxns.length(), 1)

  let mut i = 0
  while i < species.length() {
    let sp = species[i]
    match sp {
      "X" => assert_true(_pathway_double_eq(stoch[i][0], -2.0))
      "Y" => assert_true(_pathway_double_eq(stoch[i][0], -1.0))
      "Z" => assert_true(_pathway_double_eq(stoch[i][0], 3.0))
      "W" => assert_true(_pathway_double_eq(stoch[i][0], 1.0))
      _ => assert_true(false)
    }
    i = i + 1
  }
}

test "stoichiometry_glycolysis_core" {
  let sys = @src.create_example_pathway()
  let (species, rxns, _stoch) = sys.stochiometry()
  assert_eq(species.length(), 6)
  assert_eq(rxns.length(), 3)
  assert_true(_pathway_sorted_strings_equal(
    species,
    ["adp", "atp", "f16bp", "f6p", "g6p", "glucose"]
  ))
}

test "stoichiometry_empty_system" {
  let sys = @src.System::new("empty", [])
  let (species, rxns, stoch) = sys.stochiometry()
  assert_eq(species.length(), 0)
  assert_eq(rxns.length(), 0)
  assert_eq(stoch.length(), 0)
}

// =======================================================================
// 3. Network  (source / sink / interactions)
// =======================================================================

test "network_linear_source_sink" {
  let net = @src.create_example_network()
  assert_eq(net.source("S").length(), 0)
  assert_true(_pathway_sorted_strings_equal(net.sink("S"), ["A"]))

  assert_true(_pathway_sorted_strings_equal(net.source("A"), ["S"]))
  assert_true(_pathway_sorted_strings_equal(net.sink("A"), ["B"]))

  assert_true(_pathway_sorted_strings_equal(net.source("B"), ["A"]))
  assert_true(_pathway_sorted_strings_equal(net.sink("B"), ["P"]))

  assert_true(_pathway_sorted_strings_equal(net.source("P"), ["B"]))
  assert_eq(net.sink("P").length(), 0)
}

test "network_interactions" {
  let net = @src.create_example_network()
  let src_a = net.source_interactions("A")
  assert_eq(src_a.length(), 1)
  assert_eq(src_a[0].0, "S")
  assert_eq(src_a[0].1, "r1")

  let snk_a = net.sink_interactions("A")
  assert_eq(snk_a.length(), 1)
  assert_eq(snk_a[0].0, "B")
  assert_eq(snk_a[0].1, "r2")
}

test "network_multi_source" {
  let net = @src.Network::new()
  net.add_interaction("X", "Y", "e1")
  net.add_interaction("Z", "Y", "e2")
  net.add_interaction("Y", "W", "e3")
  assert_true(_pathway_sorted_strings_equal(net.source("Y"), ["X", "Z"]))
  assert_true(_pathway_sorted_strings_equal(net.sink("Y"), ["W"]))
}

test "network_unknown_species_empty" {
  let net = @src.create_example_network()
  assert_eq(net.source("unknown").length(), 0)
  assert_eq(net.sink("unknown").length(), 0)
}

// =======================================================================
// 4. PathwayGraph  (simple directed graph with labelled edges)
// =======================================================================

test "graph_basic_nodes_edges" {
  let g = @src.PathwayGraph::new()
  g.add_node("A"); g.add_node("B"); g.add_node("C")
  g.add_edge("A", "B", "e1")
  g.add_edge("B", "C", "e2")

  assert_true(_pathway_sorted_strings_equal(g.nodes(), ["A", "B", "C"]))
  assert_true(_pathway_sorted_strings_equal(g.labels(), ["e1", "e2"]))

  let (f1, t1) = g.edges("e1")
  assert_eq(f1, "A")
  assert_eq(t1, "B")
}

test "graph_children_parents_labels" {
  let g = @src.PathwayGraph::new()
  for n in ["n1", "n2", "n3", "n4", "n5"] {
    g.add_node(n)
  }
  g.add_edge("n1", "n2", "e1")
  g.add_edge("n2", "n3", "e2")
  g.add_edge("n3", "n4", "e3")
  g.add_edge("n4", "n5", "e4")
  g.add_edge("n1", "n3", "e5")
  g.add_edge("n2", "n4", "e6")

  assert_true(_pathway_sorted_strings_equal(g.children("n1"), ["e1", "e5"]))
  assert_true(_pathway_sorted_strings_equal(g.parents("n5"), ["e4"]))
}

test "graph_child_edges_parent_edges" {
  let g = @src.PathwayGraph::new()
  g.add_node("A"); g.add_node("B"); g.add_node("C")
  g.add_edge("A", "B", "e1")
  g.add_edge("A", "C", "e2")

  assert_eq(g.child_edges("A").length(), 2)

  g.add_node("D")
  g.add_edge("D", "B", "e3")
  assert_eq(g.parent_edges("B").length(), 2)
}

test "graph_duplicate_edge_label_ignored" {
  let g = @src.PathwayGraph::new()
  g.add_node("A"); g.add_node("B"); g.add_node("C")
  g.add_edge("A", "B", "shared")
  g.add_edge("A", "C", "shared")
  assert_eq(g.labels().length(), 1)
  let (f, t) = g.edges("shared")
  assert_eq(f, "A")
  assert_eq(t, "B")
}

// =======================================================================
// 5. PathwayMultiGraph  (parallel edges allowed)
// =======================================================================

test "multigraph_parallel_edges" {
  let mg = @src.PathwayMultiGraph::new()
  mg.add_node("X"); mg.add_node("Y"); mg.add_node("Z")
  mg.add_edge("X", "Y", "p1")
  mg.add_edge("X", "Y", "p2")
  mg.add_edge("Y", "Z", "yz")

  assert_true(_pathway_sorted_strings_equal(mg.labels(), ["p1", "p2", "yz"]))
  assert_eq(mg.edges("p1").length(), 1)
  assert_eq(mg.edges("p1")[0].0, "X")
  assert_eq(mg.edges("p1")[0].1, "Y")
  assert_true(_pathway_sorted_strings_equal(mg.children("X"), ["p1", "p2"]))
}

test "multigraph_empty" {
  let mg = @src.PathwayMultiGraph::new()
  assert_eq(mg.nodes().length(), 0)
  assert_eq(mg.labels().length(), 0)
}

// =======================================================================
// 6. BFS / DFS
// =======================================================================

test "bf_search_multigraph_linear" {
  let mg = @src.PathwayMultiGraph::new()
  for n in ["A", "B", "C", "D", "E"] { mg.add_node(n) }
  mg.add_edge("A", "B", "e1")
  mg.add_edge("B", "C", "e2")
  mg.add_edge("C", "D", "e3")
  mg.add_edge("D", "E", "e4")
  mg.add_edge("A", "C", "e5")

  let order = @src.bf_search(mg, root="A")
  assert_eq(order.length(), 5)
  assert_eq(order[0], "A")
  assert_true(_pathway_sorted_strings_equal(order, ["A", "B", "C", "D", "E"]))
}

test "bf_search_disconnected" {
  let mg = @src.PathwayMultiGraph::new()
  for n in ["X", "Y", "Z", "W"] { mg.add_node(n) }
  mg.add_edge("X", "Y", "xy")
  mg.add_edge("Y", "Z", "yz")

  assert_true(_pathway_sorted_strings_equal(@src.bf_search(mg, root="X"), ["X", "Y", "Z"]))
  assert_true(_pathway_sorted_strings_equal(@src.bf_search(mg, root="W"), ["W"]))
}

test "bf_search_single_node" {
  let mg = @src.PathwayMultiGraph::new()
  mg.add_node("solo")
  assert_true(_pathway_sorted_strings_equal(@src.bf_search(mg, root="solo"), ["solo"]))
}

test "bf_search_empty_graph" {
  let mg = @src.PathwayMultiGraph::new()
  assert_eq(@src.bf_search(mg).length(), 0)
}

test "df_search_multigraph_covers_all" {
  let mg = @src.PathwayMultiGraph::new()
  for n in ["A", "B", "C", "D", "E"] { mg.add_node(n) }
  mg.add_edge("A", "B", "e1")
  mg.add_edge("B", "C", "e2")
  mg.add_edge("C", "D", "e3")
  mg.add_edge("D", "E", "e4")
  mg.add_edge("A", "C", "e5")

  let order = @src.df_search(mg, root="A")
  assert_eq(order.length(), 5)
  assert_eq(order[0], "A")
  assert_true(_pathway_sorted_strings_equal(order, ["A", "B", "C", "D", "E"]))
}

test "bf_search_pathway_graph" {
  let g = @src.PathwayGraph::new()
  for n in ["n1", "n2", "n3", "n4", "n5"] { g.add_node(n) }
  g.add_edge("n1", "n2", "e1")
  g.add_edge("n2", "n3", "e2")
  g.add_edge("n3", "n4", "e3")
  g.add_edge("n4", "n5", "e4")
  g.add_edge("n1", "n3", "e5")
  g.add_edge("n2", "n4", "e6")

  let order = @src.bf_search_graph(g, root="n1")
  assert_eq(order.length(), 5)
  assert_eq(order[0], "n1")
  assert_true(_pathway_sorted_strings_equal(order, ["n1", "n2", "n3", "n4", "n5"]))
}

test "df_search_graph_single" {
  let g = @src.PathwayGraph::new()
  g.add_node("solo")
  assert_eq(@src.df_search_graph(g, root="solo").length(), 1)
}
