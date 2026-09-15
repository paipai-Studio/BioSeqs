///|
/// Tests for rrvgo module (GO term semantic reduction).

test "rrvgo_sim_matrix_diagonal" {
  let g = @src.gosemsim_example_graph()
  let terms = ["GO:0009987", "GO:0006810", "GO:0009058"]
  let sim = @src.rvgo_sim_matrix(g, terms)
  assert_eq(sim.length(), 3)
  assert_true((sim[0][0] - 1.0).abs() < 1.0e-9)
  assert_true((sim[1][1] - 1.0).abs() < 1.0e-9)
}

test "rrvgo_sim_matrix_symmetric" {
  let g = @src.gosemsim_example_graph()
  let terms = ["GO:0009987", "GO:0006810", "GO:0009058", "GO:0023046"]
  let sim = @src.rvgo_sim_matrix(g, terms)
  assert_true((sim[0][1] - sim[1][0]).abs() < 1.0e-9)
  assert_true((sim[0][2] - sim[2][0]).abs() < 1.0e-9)
}

test "rrvgo_parent_child_similar" {
  let g = @src.gosemsim_example_graph()
  // cellular_process vs transport (parent-child) should be similar
  let sim = @src.gosemsim_wang_sim(g, "GO:0009987", "GO:0006810")
  assert_true(sim > 0.3)
}

test "rrvgo_dist_matrix_range" {
  let g = @src.gosemsim_example_graph()
  let sim = @src.rvgo_sim_matrix(g, ["GO:0009987", "GO:0006810"])
  let dist = @src.rvgo_dist_matrix(sim)
  assert_true(dist[0][0] < 1.0e-9)
  assert_true(dist[0][1] >= 0.0 && dist[0][1] <= 1.0)
}

test "rrvgo_hclust_merges" {
  let g = @src.gosemsim_example_graph()
  let terms = ["GO:0009987", "GO:0006810", "GO:0009058", "GO:0023046"]
  let sim = @src.rvgo_sim_matrix(g, terms)
  let dist = @src.rvgo_dist_matrix(sim)
  let merges = @src.rvgo_hclust(dist)
  assert_eq(merges.length(), 3)
}

test "rrvgo_reduce_groups" {
  let g = @src.gosemsim_example_graph()
  let terms = ["GO:0009987", "GO:0006810", "GO:0023046", "GO:0009058"]
  let scores = [0.01, 0.02, 0.03, 0.04]
  let reduced = @src.rvgo_reduce(g, terms, scores, threshold=0.5)
  assert_true(@src.rvgo_n_clusters(reduced) >= 2)
}

test "rrvgo_representative_lowest_score" {
  let g = @src.gosemsim_example_graph()
  // two very similar terms (parent-child): rep should be the one with lower score
  let terms = ["GO:0009987", "GO:0006810"]
  let scores = [0.05, 0.01]  // GO:0006810 has lower (better) score
  let reduced = @src.rvgo_reduce(g, terms, scores, threshold=0.5)
  // both should collapse to 1 cluster; representative = lowest score
  assert_eq(@src.rvgo_n_clusters(reduced), 1)
  assert_eq(@src.rvgo_representative(reduced, 0), "GO:0006810")
}

test "rrvgo_export" {
  let g = @src.gosemsim_example_graph()
  let terms = ["GO:0009987", "GO:0006810"]
  let scores = [0.01, 0.02]
  let reduced = @src.rvgo_reduce(g, terms, scores, threshold=0.5)
  let tsv = @src.rvgo_export(reduced)
  assert_true(tsv.contains("term_id"))
  assert_true(tsv.contains("GO:0009987"))
}

test "rrvgo_mds_2d" {
  let g = @src.gosemsim_example_graph()
  let terms = ["GO:0009987", "GO:0006810", "GO:0009058"]
  let sim = @src.rvgo_sim_matrix(g, terms)
  let dist = @src.rvgo_dist_matrix(sim)
  let coords = @src.rvgo_mds_2d(dist)
  assert_eq(coords.length(), 3)
}
