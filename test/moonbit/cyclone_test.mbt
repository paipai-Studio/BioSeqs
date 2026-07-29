///|
/// Tests for Cyclone module - Cell cycle scoring using gene pairs.

///|
test "cyclone_get_gene_pairs" {
  let pairs = @src.cyclone_get_gene_pairs()
  assert_true(pairs.length() > 0)
  // Should have pairs for all 4 phases
  let phases = @src.cyclone_get_phases(pairs)
  assert_true(phases.contains("G1"))
  assert_true(phases.contains("S"))
  assert_true(phases.contains("G2"))
  assert_true(phases.contains("M"))
}

///|
test "cyclone_get_phases" {
  let pairs = @src.cyclone_get_gene_pairs()
  let phases = @src.cyclone_get_phases(pairs)
  assert_true(phases.length() == 4)
}

///|
test "cyclone_count_pairs_by_phase" {
  let pairs = @src.cyclone_get_gene_pairs()
  let counts = @src.cyclone_count_pairs_by_phase(pairs)
  assert_true(counts.contains("G1"))
  assert_true(counts.contains("S"))
  assert_true(counts.contains("G2"))
  assert_true(counts.contains("M"))
  assert_true(counts.get("G1").unwrap_or(0) > 0)
  assert_true(counts.get("S").unwrap_or(0) > 0)
}

///|
test "cyclone_score_phase" {
  let pairs = @src.cyclone_get_gene_pairs()
  let test_data = @src.cyclone_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  // Score the first cell (should be G1 phase)
  let cell_expr : Array[Double] = Array::new()
  let n_genes = counts.length()
  let mut g = 0
  while g < n_genes {
    cell_expr.push(counts[g][0])
    g = g + 1
  }
  
  let g1_score = @src.cyclone_score_phase(cell_expr, gene_names, pairs, "G1")
  assert_true(g1_score >= 0.0)
  assert_true(g1_score <= 1.0)
}

///|
test "cyclone_score_cell" {
  let pairs = @src.cyclone_get_gene_pairs()
  let test_data = @src.cyclone_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let cell_expr : Array[Double] = Array::new()
  let n_genes = counts.length()
  let mut g = 0
  while g < n_genes {
    cell_expr.push(counts[g][0])
    g = g + 1
  }
  
  let scores = @src.cyclone_score_cell(cell_expr, gene_names, pairs)
  assert_true(scores.size() > 0)
  assert_true(scores.contains("G1"))
  assert_true(scores.contains("S"))
  assert_true(scores.contains("G2"))
  assert_true(scores.contains("M"))
}

///|
test "cyclone_assign_phase" {
  let scores = Map([("G1", 0.8), ("S", 0.3), ("G2", 0.5), ("M", 0.2)])
  let phase = @src.cyclone_assign_phase(scores)
  assert_eq(phase, "G1")
}

///|
test "cyclone_assign_phase_best" {
  let scores = Map([("G1", 0.2), ("S", 0.9), ("G2", 0.4), ("M", 0.3)])
  let phase = @src.cyclone_assign_phase(scores)
  assert_eq(phase, "S")
}

///|
test "cyclone_score_single_cell" {
  let pairs = @src.cyclone_get_gene_pairs()
  let test_data = @src.cyclone_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let cell_expr : Array[Double] = Array::new()
  let n_genes = counts.length()
  let mut g = 0
  while g < n_genes {
    cell_expr.push(counts[g][0])
    g = g + 1
  }
  
  let result = @src.cyclone_score_single_cell(cell_expr, "test_cell", gene_names, pairs)
  assert_eq(result.cell_id, "test_cell")
  assert_true(result.scores.size() == 4)
  assert_true(result.phases.length() == 4)
  assert_true(result.assigned_phase.length() > 0)
  assert_true(result.g1_score >= 0.0 && result.g1_score <= 1.0)
  assert_true(result.s_score >= 0.0 && result.s_score <= 1.0)
  assert_true(result.g2_score >= 0.0 && result.g2_score <= 1.0)
  assert_true(result.m_score >= 0.0 && result.m_score <= 1.0)
}

///|
test "cyclone_score_cells" {
  let test_data = @src.cyclone_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let results = @src.cyclone_score_cells(counts, cell_ids, gene_names)
  assert_eq(results.length(), 20)
  
  // Check first cell
  assert_eq(results[0].cell_id, "cell_1")
  assert_true(results[0].assigned_phase.length() > 0)
}

///|
test "cyclone_phase_distribution" {
  let test_data = @src.cyclone_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let results = @src.cyclone_score_cells(counts, cell_ids, gene_names)
  let dist = @src.cyclone_phase_distribution(results)
  
  assert_true(dist.size() > 0)
  
  // Sum of proportions should be approximately 1.0
  let keys = dist.keys()
  let mut sum = 0.0
  for key in keys {
    sum = sum + dist.get(key).unwrap_or(0.0)
  }
  assert_true((sum - 1.0).abs() < 0.01)
}

///|
test "cyclone_create_test_data" {
  let test_data = @src.cyclone_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  assert_eq(cell_ids.length(), 20)
  assert_true(gene_names.length() > 0)
  assert_true(counts.length() > 0)
  assert_true(counts[0].length() == 20)
}

///|
test "cyclone_params_new" {
  let params = @src.CycloneParams::new()
  assert_true((params.min_score - 0.0).abs() < 0.001)
  assert_eq(params.n_top_pairs, 1)
}

///|
test "cyclone_params_with_params" {
  let params = @src.CycloneParams::with_params(0.1, 5)
  assert_true((params.min_score - 0.1).abs() < 0.001)
  assert_eq(params.n_top_pairs, 5)
}

///|
test "cyclone_score_empty_cells" {
  let gene_names : Array[String] = Array::new()
  let cell_ids : Array[String] = Array::new()
  let counts : Array[Array[Double]] = Array::new()
  
  let results = @src.cyclone_score_cells(counts, cell_ids, gene_names)
  assert_eq(results.length(), 0)
}

///|
test "cyclone_phase_distribution_empty" {
  let results : Array[@src.CycloneResult] = Array::new()
  let dist = @src.cyclone_phase_distribution(results)
  assert_eq(dist.size(), 0)
}

///|
test "cyclone_average_scores" {
  let test_data = @src.cyclone_create_test_data()
  let (counts, cell_ids, gene_names) = test_data
  
  let results = @src.cyclone_score_cells(counts, cell_ids, gene_names)
  let avg = @src.cyclone_average_scores(results)
  
  assert_true(avg.size() == 4)
  assert_true(avg.contains("G1"))
  assert_true(avg.contains("S"))
  assert_true(avg.contains("G2"))
  assert_true(avg.contains("M"))
}

///|
test "cyclone_average_scores_empty" {
  let results : Array[@src.CycloneResult] = Array::new()
  let avg = @src.cyclone_average_scores(results)
  assert_eq(avg.size(), 0)
}