///|
/// Tests for GOSemSim GO semantic similarity module.

test "gosemsim_example_graph" {
  let g = @src.gosemsim_example_graph()
  assert_true(g.root_id == "GO:0008150")
  match g.get_term("GO:0008150") {
    Some(term) => assert_true(term.go_id == "GO:0008150")
    None => assert_true(false)
  }
}

test "gosemsim_term_self_similarity" {
  let g = @src.gosemsim_example_graph()
  // A term compared to itself should have maximal similarity.
  let resnik = @src.gosemsim_term_sim(g, "GO:0006810", "GO:0006810", @src.resnik_measure())
  assert_true(resnik > 0.0)
  let lin = @src.gosemsim_term_sim(g, "GO:0006810", "GO:0006810", @src.lin_measure())
  assert_true(lin > 0.99 && lin <= 1.0 + 1.0e-6)
}

test "gosemsim_resnik_similarity" {
  let g = @src.gosemsim_example_graph()
  // Parent-child should have high similarity
  let sim = @src.gosemsim_term_sim(g, "GO:0006810", "GO:0009987", @src.resnik_measure())
  assert_true(sim > 0.0)
}

test "gosemsim_lin_similarity" {
  let g = @src.gosemsim_example_graph()
  let sim = @src.gosemsim_term_sim(g, "GO:0006810", "GO:0007165", @src.lin_measure())
  assert_true(sim >= 0.0 && sim <= 1.0)
}

test "gosemsim_rel_similarity" {
  let g = @src.gosemsim_example_graph()
  let sim = @src.gosemsim_term_sim(g, "GO:0006810", "GO:0007165", @src.rel_measure())
  assert_true(sim >= 0.0 && sim <= 1.0)
}

test "gosemsim_jiang_similarity" {
  let g = @src.gosemsim_example_graph()
  let sim = @src.gosemsim_term_sim(g, "GO:0006810", "GO:0007165", @src.jiang_measure())
  assert_true(sim >= 0.0)
}

test "gosemsim_wang_similarity" {
  let g = @src.gosemsim_example_graph()
  let sim = @src.gosemsim_term_sim(g, "GO:0006810", "GO:0007165", @src.wang_measure())
  assert_true(sim >= 0.0 && sim <= 1.0)
}

test "gosemsim_missing_term" {
  let g = @src.gosemsim_example_graph()
  let sim = @src.gosemsim_term_sim(g, "GO:9999999", "GO:0006810", @src.resnik_measure())
  assert_true(sim == 0.0)
}

test "gosemsim_gen_sim_max" {
  let g = @src.gosemsim_example_graph()
  let go1 = ["GO:0006810", "GO:0007165"]
  let go2 = ["GO:0006810", "GO:0009987"]
  let sim = @src.gosemsim_gen_sim(g, go1, go2, @src.resnik_measure(), combine = "max")
  assert_true(sim > 0.0)
}

test "gosemsim_gen_sim_avg" {
  let g = @src.gosemsim_example_graph()
  let go1 = ["GO:0006810", "GO:0007165"]
  let go2 = ["GO:0006810", "GO:0009987"]
  let sim = @src.gosemsim_gen_sim(g, go1, go2, @src.lin_measure(), combine = "avg")
  assert_true(sim >= 0.0 && sim <= 1.0)
}

test "gosemsim_gen_sim_rcmax" {
  let g = @src.gosemsim_example_graph()
  let go1 = ["GO:0006810", "GO:0007165"]
  let go2 = ["GO:0006810", "GO:0009987"]
  let sim = @src.gosemsim_gen_sim(g, go1, go2, @src.rel_measure(), combine = "rcmax")
  assert_true(sim >= 0.0 && sim <= 1.0)
}

test "gosemsim_graph_operations" {
  let g = @src.GOGraph::new()
  let g1 = g.add_root("TEST:001")
  let g2 = g1.add_term(@src.GOTermNode::new("TEST:001", ic = 0.0))
  let g3 = g2.add_term(@src.GOTermNode::new("TEST:002", ic = 0.5, parents = ["TEST:001"]))
  match g3.get_term("TEST:002") {
    Some(term) => assert_true(term.ic == 0.5)
    None => assert_true(false)
  }
}
