///|
test "gostats_log_choose_C5_2" {
  // C(5,2) = 10, log(10) ≈ 2.302585
  let result = @src.gostats_log_choose(5, 2)
  let expected = @math.ln(10.0)
  assert_true((result - expected).abs() < 1.0e-6)
}

///|
test "gostats_log_choose_C10_3" {
  // C(10,3) = 120, log(120) ≈ 4.787492
  let result = @src.gostats_log_choose(10, 3)
  let expected = @math.ln(120.0)
  assert_true((result - expected).abs() < 1.0e-6)
}

///|
test "gostats_log_choose_C6_3" {
  // C(6,3) = 20, log(20) ≈ 2.995732
  let result = @src.gostats_log_choose(6, 3)
  let expected = @math.ln(20.0)
  assert_true((result - expected).abs() < 1.0e-6)
}

///|
test "gostats_log_choose_k0_is_0" {
  // C(n,0) = 1, log(1) = 0.0
  assert_true(@src.gostats_log_choose(5, 0) == 0.0)
  assert_true(@src.gostats_log_choose(10, 0) == 0.0)
  assert_true(@src.gostats_log_choose(100, 0) == 0.0)
}

///|
test "gostats_log_choose_kn_is_0" {
  // C(n,n) = 1, log(1) = 0.0
  assert_true(@src.gostats_log_choose(5, 5) == 0.0)
  assert_true(@src.gostats_log_choose(10, 10) == 0.0)
  assert_true(@src.gostats_log_choose(100, 100) == 0.0)
}

///|
test "gostats_log_choose_symmetry" {
  // C(n,k) == C(n,n-k)
  let n = 10
  assert_true(
    (@src.gostats_log_choose(n, 2) - @src.gostats_log_choose(n, 8)).abs() < 1.0e-10
  )
  assert_true(
    (@src.gostats_log_choose(n, 3) - @src.gostats_log_choose(n, 7)).abs() < 1.0e-10
  )
}

///|
test "gostats_bh_fdr_basic" {
  let p_values = [0.01, 0.04, 0.05, 0.20]
  let fdr = @src.gostats_bh_fdr(p_values)
  assert_eq(fdr.length(), 4)
  // All adjusted p-values should be between 0 and 1
  let mut i = 0
  while i < fdr.length() {
    assert_true(fdr[i] >= 0.0 && fdr[i] <= 1.0)
    i = i + 1
  }
}

///|
test "gostats_bh_fdr_monotone" {
  // When sorted by original p-value, adjusted p-values should be non-decreasing
  let p_values = [0.005, 0.01, 0.03, 0.05, 0.1]
  let fdr = @src.gostats_bh_fdr(p_values)
  // p_values are already sorted, so fdr should be non-decreasing
  let mut i = 0
  while i < fdr.length() - 1 {
    assert_true(fdr[i] <= fdr[i + 1] + 1.0e-10)
    i = i + 1
  }
}

///|
test "gostats_bh_fdr_unsorted" {
  // BH FDR should work correctly even with unsorted input
  let p_values = [0.1, 0.01, 0.05, 0.005]
  let fdr = @src.gostats_bh_fdr(p_values)
  assert_eq(fdr.length(), 4)
  let mut i = 0
  while i < fdr.length() {
    assert_true(fdr[i] >= 0.0 && fdr[i] <= 1.0)
    i = i + 1
  }
}

///|
test "gostats_bh_fdr_empty" {
  let fdr = @src.gostats_bh_fdr([])
  assert_eq(fdr.length(), 0)
}

///|
test "gostats_hypergeometric_test_basic" {
  // Build a simple GO graph
  let term1 = @src.GOStatsTerm::new(
    "GO:0008150",
    "biological_process",
    "BP",
    [],
    ["g1", "g2", "g3", "g4", "g5"],
  )
  let term2 = @src.GOStatsTerm::new(
    "GO:0003674",
    "molecular_function",
    "MF",
    [],
    ["g3", "g4", "g5", "g6", "g7", "g8"],
  )
  let graph = @src.GOStatsGraph::new(Map([]))
  let graph = graph.add_term(term1)
  let graph = graph.add_term(term2)

  let universe = ["g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10"]
  let selected = ["g1", "g2", "g3"]

  let hg_test = @src.GOHyperGeoTest::new(universe, selected, graph)
  let results = @src.gostats_hypergeometric_test(hg_test)

  // Should have results for terms with overlap > 0
  assert_true(results.length() > 0)

  // All p-values should be valid
  let mut i = 0
  while i < results.length() {
    assert_true(results[i].p_value >= 0.0 && results[i].p_value <= 1.0)
    assert_true(results[i].fdr >= 0.0 && results[i].fdr <= 1.0)
    i = i + 1
  }
}

///|
test "gostats_hypergeometric_test_with_ontology_filter" {
  let term1 = @src.GOStatsTerm::new(
    "GO:0008150",
    "biological_process",
    "BP",
    [],
    ["g1", "g2", "g3"],
  )
  let term2 = @src.GOStatsTerm::new(
    "GO:0003674",
    "molecular_function",
    "MF",
    [],
    ["g3", "g4", "g5"],
  )
  let graph = @src.GOStatsGraph::new(Map([]))
  let graph = graph.add_term(term1)
  let graph = graph.add_term(term2)

  let universe = ["g1", "g2", "g3", "g4", "g5", "g6"]
  let selected = ["g1", "g2", "g3"]

  let hg_test = @src.GOHyperGeoTest::new(universe, selected, graph)
  let results = @src.gostats_hypergeometric_test(hg_test, ontology="BP")

  // Should only include BP terms
  let mut i = 0
  while i < results.length() {
    assert_eq(results[i].ontology, "BP")
    i = i + 1
  }
}

///|
test "gostats_conditional_test_basic" {
  // Build a GO graph with parent-child relationships
  let child_term = @src.GOStatsTerm::new(
    "GO:0006915",
    "apoptotic process",
    "BP",
    ["GO:0008150"],
    ["g1", "g2"],
  )
  let parent_term = @src.GOStatsTerm::new(
    "GO:0008150",
    "biological_process",
    "BP",
    [],
    ["g1", "g2", "g3", "g4"],
  )
  let graph = @src.GOStatsGraph::new(Map([]))
  let graph = graph.add_term(child_term)
  let graph = graph.add_term(parent_term)

  let universe = [
    "g1", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "g10",
  ]
  let selected = ["g1", "g2", "g3"]

  let hg_test = @src.GOHyperGeoTest::new(universe, selected, graph)
  let results = @src.gostats_conditional_test(hg_test)

  // Results should have valid p-values
  let mut i = 0
  while i < results.length() {
    assert_true(results[i].p_value >= 0.0 && results[i].p_value <= 1.0)
    assert_true(results[i].fdr >= 0.0 && results[i].fdr <= 1.0)
    i = i + 1
  }
}

///|
test "gostats_filter_results_by_fdr" {
  let r1 = @src.GOHyperGeoResult::new(
    "GO:0008150", "bp", "BP", 0.001, 0.01, 3.0, 1.5, 3, 5, 100,
  )
  let r2 = @src.GOHyperGeoResult::new(
    "GO:0003674", "mf", "MF", 0.05, 0.10, 1.5, 2.0, 2, 5, 100,
  )
  let r3 = @src.GOHyperGeoResult::new(
    "GO:0005515", "cc", "CC", 0.20, 0.40, 0.8, 3.0, 1, 5, 100,
  )

  let results = [r1, r2, r3]

  // FDR <= 0.05, min_observed = 1
  let filtered = @src.gostats_filter_results(results, 0.05, 1)
  assert_eq(filtered.length(), 1)
  assert_eq(filtered[0].go_id, "GO:0008150")

  // FDR <= 0.10, min_observed = 1
  let filtered2 = @src.gostats_filter_results(results, 0.10, 1)
  assert_eq(filtered2.length(), 2)

  // FDR <= 0.50, min_observed = 2
  let filtered3 = @src.gostats_filter_results(results, 0.50, 2)
  assert_eq(filtered3.length(), 2) // r1 (obs=3) and r2 (obs=2), r3 has obs=1

  // FDR <= 1.0, min_observed = 0
  let filtered4 = @src.gostats_filter_results(results, 1.0, 0)
  assert_eq(filtered4.length(), 3)
}

///|
test "gostats_graph_add_and_get_term" {
  let term = @src.GOStatsTerm::new(
    "GO:0008150", "biological_process", "BP", [], ["g1"],
  )
  let graph = @src.GOStatsGraph::new(Map([]))
  let graph = graph.add_term(term)

  assert_true(graph.has_term("GO:0008150"))
  assert_true(!graph.has_term("GO:9999999"))

  let retrieved = graph.get_term("GO:0008150")
  match retrieved {
    Some(t) => {
      assert_eq(t.go_id, "GO:0008150")
      assert_eq(t.name, "biological_process")
      assert_eq(t.ontology, "BP")
    }
    None => assert_true(false)
  }
}

///|
test "gostats_summary" {
  let r1 = @src.GOHyperGeoResult::new(
    "GO:0008150", "bp", "BP", 0.001, 0.01, 3.0, 1.5, 3, 5, 100,
  )
  let r2 = @src.GOHyperGeoResult::new(
    "GO:0003674", "mf", "MF", 0.05, 0.10, 1.5, 2.0, 2, 5, 100,
  )
  let results = [r1, r2]
  let s = @src.gostats_summary(results)
  assert_true(s.length() > 0)
}
