///|
/// Tests for SingleR module - Cell type annotation by reference profiles.

///|
test "single_r_create_reference_data" {
  let ref_data = @src.single_r_create_reference_data()
  assert_eq(ref_data.n_profiles, 10)
  assert_eq(ref_data.n_genes, 66)
  assert_true(ref_data.cell_types.length() > 0)
  assert_true(ref_data.profiles.length() > 0)
}

///|
test "single_r_create_test_data" {
  let test_data = @src.single_r_create_test_data()
  assert_eq(test_data.n_cells, 50)
  assert_eq(test_data.n_genes, 66)
  assert_eq(test_data.cell_names.length(), 50)
  assert_true(test_data.counts.length() > 0)
  assert_true(test_data.counts[0].length() > 0)
}

///|
test "single_r_spearman_correlation" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let corr = @src.single_r_spearman_correlation(x, y)
  assert_true(corr > 0.99)
  assert_true(corr <= 1.0)
}

///|
test "single_r_spearman_correlation_negative" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [5.0, 4.0, 3.0, 2.0, 1.0]
  let corr = @src.single_r_spearman_correlation(x, y)
  assert_true(corr < -0.99)
  assert_true(corr >= -1.0)
}

///|
test "single_r_pearson_correlation" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let corr = @src.single_r_pearson_correlation(x, y)
  assert_true(corr > 0.99)
  assert_true(corr <= 1.0)
}

///|
test "single_r_pearson_correlation_negative" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [5.0, 4.0, 3.0, 2.0, 1.0]
  let corr = @src.single_r_pearson_correlation(x, y)
  assert_true(corr < -0.99)
  assert_true(corr >= -1.0)
}

///|
test "single_r_correlation_noise" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [3.0, 3.5, 3.0, 3.5, 3.0]
  let corr = @src.single_r_pearson_correlation(x, y)
  assert_true(corr.abs() < 0.5)
}

///|
test "single_r_rank_values" {
  let vals = [3.0, 1.0, 2.0]
  let ranks = @src.single_r_rank_values(vals)
  assert_eq(ranks.length(), 3)
  assert_eq(ranks[1], 1.0)
  assert_eq(ranks[2], 2.0)
  assert_eq(ranks[0], 3.0)
}

///|
test "single_r_rank_values_ties" {
  let vals = [2.0, 2.0, 1.0]
  let ranks = @src.single_r_rank_values(vals)
  assert_eq(ranks.length(), 3)
  assert_eq(ranks[2], 1.0)
  assert_true((ranks[0] - 2.5).abs() < 0.01)
  assert_true((ranks[1] - 2.5).abs() < 0.01)
}

///|
test "single_r_compute_correlations" {
  let ref_data = @src.single_r_create_reference_data()
  let cell_expr = ref_data.profiles[0].expression
  let scores = @src.single_r_compute_correlations(cell_expr, ref_data, "spearman")
  assert_eq(scores.length(), ref_data.n_profiles)
  assert_true(scores[0] > 0.0)
}

///|
test "single_r_get_top_scores" {
  let scores = [0.1, 0.5, 0.3, 0.8, 0.2]
  let top = @src.single_r_get_top_scores(scores, 3)
  assert_eq(top.length(), 3)
  assert_eq(top[0].0, 3) // Index of max score (0.8)
  assert_eq(top[0].1, 0.8)
}

///|
test "single_r_aggregate_scores_by_type" {
  let ref_data = @src.single_r_create_reference_data()
  let scores = @src.single_r_compute_correlations(ref_data.profiles[0].expression, ref_data, "spearman")
  let agg = @src.single_r_aggregate_scores_by_type(scores, ref_data, 1)
  assert_true(agg.size() > 0)
  assert_true(agg.contains("T cells"))
}

///|
test "single_r_get_best_label" {
  let scores = Map([("T cells", 0.85), ("B cells", 0.3), ("NK cells", 0.5)])
  let (label, score) = @src.single_r_get_best_label(scores)
  assert_eq(label, "T cells")
  assert_eq(score, 0.85)
}

///|
test "single_r_compute_delta_score" {
  let scores = Map([("T cells", 0.85), ("B cells", 0.60), ("NK cells", 0.3)])
  let delta = @src.single_r_compute_delta_score(scores)
  assert_true(delta > 0.0)
  assert_true((delta - 0.25).abs() < 0.01)
}

///|
test "single_r_annotate_cell" {
  let ref_data = @src.single_r_create_reference_data()
  let params = @src.SingleRParams::new()
  
  // Test with a T cell-like expression
  let cell_expr = ref_data.profiles[0].expression
  let result = @src.single_r_annotate_cell(cell_expr, "test_t_cell", ref_data, params)
  
  assert_eq(result.cell_id, "test_t_cell")
  assert_true(result.scores.length() > 0)
  assert_true(result.first_annotation_fine.length() > 0)
  assert_true(result.annotations_fine.length() > 0)
}

///|
test "single_r_annotate_cell_b_cell" {
  let ref_data = @src.single_r_create_reference_data()
  let params = @src.SingleRParams::new()
  
  // Test with a B cell-like expression
  let cell_expr = ref_data.profiles[1].expression
  let result = @src.single_r_annotate_cell(cell_expr, "test_b_cell", ref_data, params)
  
  assert_true(result.first_annotation_fine.length() > 0)
  assert_true(result.scores[0] > 0.0)
}

///|
test "single_r_annotate_all_cells" {
  let ref_data = @src.single_r_create_reference_data()
  let test_data = @src.single_r_create_test_data()
  let params = @src.SingleRParams::new()
  
  let result = @src.single_r_annotate_cells(test_data, ref_data, params)
  
  assert_eq(result.cell_ids.length(), 50)
  assert_eq(result.labels.length(), 50)
  assert_eq(result.scores.length(), 50)
  assert_true(result.labels.length() > 0)
  assert_true(result.scores.length() > 0)
}

///|
test "single_r_annotation_summary" {
  let ref_data = @src.single_r_create_reference_data()
  let test_data = @src.single_r_create_test_data()
  let params = @src.SingleRParams::new()
  
  let result = @src.single_r_annotate_cells(test_data, ref_data, params)
  let summary = @src.single_r_annotation_summary(result)
  
  assert_true(summary.size() > 0)
  assert_true(summary.contains("n_cells"))
  assert_eq(summary.get("n_cells").unwrap_or(0.0), 50.0)
  assert_true(summary.contains("avg_score"))
  assert_true(summary.get("avg_score").unwrap_or(0.0) > 0.0)
}

///|
test "single_r_params_default" {
  let params = @src.SingleRParams::new()
  assert_eq(params.method, "spearman")
  assert_eq(params.n_top_scores, 1)
  assert_true(params.fine_tune)
  assert_true((params.min_score - 0.0).abs() < 0.001)
}

///|
test "single_r_params_pearson" {
  let params = @src.SingleRParams::with_method("pearson")
  assert_eq(params.method, "pearson")
}

///|
test "single_r_pearson_annotation" {
  let ref_data = @src.single_r_create_reference_data()
  let test_data = @src.single_r_create_test_data()
  let params = @src.SingleRParams::with_method("pearson")
  
  let result = @src.single_r_annotate_cells(test_data, ref_data, params)
  
  assert_eq(result.cell_ids.length(), 50)
  assert_true(result.scores.length() > 0)
  assert_true(result.scores[0] > 0.0)
}

///|
test "single_r_reference_dataset_from_matrix" {
  let expression = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
  ]
  let cell_types = ["TypeA", "TypeB", "TypeC"]
  let gene_names = ["Gene1", "Gene2", "Gene3"]
  
  let reference = @src.ReferenceDataset::from_matrix(expression, cell_types, gene_names)
  
  assert_eq(reference.n_profiles, 3)
  assert_eq(reference.n_genes, 3)
  assert_eq(reference.cell_types.length(), 3)
  assert_eq(reference.profiles.length(), 3)
}

///|
test "single_r_fine_tune_disabled" {
  let ref_data = @src.single_r_create_reference_data()
  let params = @src.SingleRParams::with_all(false, "spearman", 1, 0.0, false)
  
  let cell_expr = ref_data.profiles[0].expression
  let result_no_finetune = @src.single_r_annotate_cell(cell_expr, "test", ref_data, params)
  
  assert_true(result_no_finetune.first_annotation_fine.length() > 0)
}

///|
test "single_r_fine_tune_enabled" {
  let ref_data = @src.single_r_create_reference_data()
  let params = @src.SingleRParams::with_all(true, "spearman", 1, 0.0, false)
  
  let cell_expr = ref_data.profiles[0].expression
  let result_finetune = @src.single_r_annotate_cell(cell_expr, "test", ref_data, params)
  
  assert_true(result_finetune.first_annotation_fine.length() > 0)
  assert_true(result_finetune.scores[0] > 0.0)
}

///|
test "single_r_correlation_short_vectors" {
  let x = [1.0, 2.0]
  let y = [3.0, 4.0]
  let corr = @src.single_r_spearman_correlation(x, y)
  assert_eq(corr, 0.0)
  let corr2 = @src.single_r_pearson_correlation(x, y)
  assert_eq(corr2, 0.0)
}

///|
test "single_r_identical_profiles" {
  let ref_data = @src.single_r_create_reference_data()
  // Annotate a cell with the same profile as reference profile 0
  let cell_expr = ref_data.profiles[0].expression
  let result = @src.single_r_annotate_cell(cell_expr, "identical", ref_data, @src.SingleRParams::new())
  
  assert_true(result.scores[0] > 0.99)
  assert_eq(result.first_annotation_fine, "T cells")
}