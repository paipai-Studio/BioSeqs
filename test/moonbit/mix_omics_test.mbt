///|
/// mixOmics Tests - Multi-omics data integration

///|
test "pls_new_defaults" {
  let opts = @src.pls_new()
  assert_eq(opts.ncomp, 2)
  assert_eq(opts.mode, "regression")
  assert_true(opts.scale)
  assert_eq(opts.max_iter, 500)
  assert_true((opts.tol - 0.000001).abs() < 1.0e-15)
}

///|
test "pls_options_full_custom" {
  let opts = @src.pls_options_full(5, "canonical", false, 200, 1.0e-8)
  assert_eq(opts.ncomp, 5)
  assert_eq(opts.mode, "canonical")
  assert_false(opts.scale)
  assert_eq(opts.max_iter, 200)
  assert_true((opts.tol - 1.0e-8).abs() < 1.0e-15)
}

///|
test "spls_new_defaults" {
  let opts = @src.spls_new()
  assert_eq(opts.ncomp, 2)
  assert_eq(opts.keep_x.length(), 0)
  assert_eq(opts.keep_y.length(), 0)
  assert_true(opts.scale)
}

///|
test "spls_options_full_custom" {
  let opts = @src.spls_options_full(3, [2, 4], [1], false, 100, 1.0e-6)
  assert_eq(opts.ncomp, 3)
  assert_eq(opts.keep_x.length(), 2)
  assert_eq(opts.keep_x[0], 2)
  assert_eq(opts.keep_x[1], 4)
  assert_eq(opts.keep_y.length(), 1)
  assert_eq(opts.keep_y[0], 1)
  assert_false(opts.scale)
}

///|
test "diablo_new_defaults" {
  let opts = @src.diablo_new()
  assert_eq(opts.ncomp, 2)
  assert_eq(opts.keep_variables.length(), 0)
  assert_eq(opts.design.length(), 0)
  assert_true(opts.scale)
}

///|
test "diablo_options_full_custom" {
  let design = [[0.0, 1.0], [1.0, 0.0]]
  let opts = @src.diablo_options_full(3, [[2, 3], [1]], design, false, 200, 1.0e-7)
  assert_eq(opts.ncomp, 3)
  assert_eq(opts.keep_variables.length(), 2)
  assert_eq(opts.keep_variables[0].length(), 2)
  assert_eq(opts.keep_variables[1].length(), 1)
  assert_eq(opts.design.length(), 2)
  assert_false(opts.scale)
}

///|
test "pls_basic_dimensions" {
  let x = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.pls_new()
  let result = @src.run_pls(x, y, opts)
  assert_true(result.ncomp >= 1)
  assert_eq(result.n, 5)
  assert_eq(result.p, 4)
  assert_eq(result.q, 1)
}

///|
test "pls_component_extraction" {
  let x = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
    [4.0, 5.0, 6.0],
    [5.0, 6.0, 7.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.pls_new()
  let result = @src.run_pls(x, y, opts)
  assert_true(result.ncomp >= 1)
  assert_true(result.x_scores.length() >= 1)
  assert_true(result.y_scores.length() >= 1)
  assert_true(result.x_loadings.length() >= 1)
  assert_true(result.y_loadings.length() >= 1)
  assert_true(result.x_weights.length() >= 1)
  assert_true(result.x_scores.length() == 5)
  assert_true(result.x_loadings.length() == 3)
  assert_true(result.x_weights.length() == 3)
}

///|
test "pls_single_component" {
  let x = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 4.0],
    [4.0, 5.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0]]
  let opts = @src.pls_options_full(1, "regression", true, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.ncomp >= 1)
  assert_true(result.x_scores.length() >= 1)
  assert_true(result.x_loadings.length() >= 1)
}

///|
test "pls_nipals_convergence" {
  let x = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
    [4.0, 5.0, 6.0],
    [5.0, 6.0, 7.0],
    [6.0, 7.0, 8.0],
    [7.0, 8.0, 9.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0], [6.0], [7.0]]
  let opts = @src.pls_options_full(3, "regression", true, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.ncomp >= 1)
  assert_true(result.x_scores.length() >= 1)
}

///|
test "pls_data_scaling" {
  let x = [
    [1.0, 100.0],
    [2.0, 200.0],
    [3.0, 300.0],
    [4.0, 400.0],
    [5.0, 500.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts_scaled = @src.pls_options_full(1, "regression", true, 500, 0.000001)
  let opts_unscaled = @src.pls_options_full(1, "regression", false, 500, 0.000001)
  let result_scaled = @src.run_pls(x, y, opts_scaled)
  let result_unscaled = @src.run_pls(x, y, opts_unscaled)
  assert_true(result_scaled.ncomp >= 1)
  assert_true(result_unscaled.ncomp >= 1)
}

///|
test "pls_explained_variance" {
  let x = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
    [4.0, 5.0, 6.0],
    [5.0, 6.0, 7.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.pls_new()
  let result = @src.run_pls(x, y, opts)
  assert_true(result.explained_variance.length() >= 1)
}

///|
test "pls_ncomp_exceeds_samples" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let y = [[1.0], [2.0]]
  let opts = @src.pls_options_full(10, "regression", true, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.ncomp >= 1)
  assert_true(result.ncomp <= 2)
}

///|
test "pls_ncomp_exceeds_variables" {
  let x = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
  ]
  let y = [[1.0], [2.0], [3.0]]
  let opts = @src.pls_options_full(10, "regression", true, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.ncomp >= 1)
  assert_true(result.ncomp <= 3)
}

///|
test "pls_empty_data" {
  let x : Array[Array[Double]] = []
  let y : Array[Array[Double]] = []
  let opts = @src.pls_new()
  let result = @src.run_pls(x, y, opts)
  assert_eq(result.n, 0)
  assert_eq(result.p, 0)
  assert_eq(result.q, 0)
}

///|
test "pls_single_sample" {
  let x = [[1.0, 2.0, 3.0]]
  let y = [[1.0]]
  let opts = @src.pls_new()
  let result = @src.run_pls(x, y, opts)
  assert_eq(result.n, 1)
  assert_eq(result.p, 3)
  assert_eq(result.q, 1)
}

///|
test "pls_multivariate_y" {
  let x = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 4.0],
    [4.0, 5.0],
    [5.0, 6.0],
  ]
  let y = [[1.0, 0.5], [2.0, 1.0], [3.0, 1.5], [4.0, 2.0], [5.0, 2.5]]
  let opts = @src.pls_new()
  let result = @src.run_pls(x, y, opts)
  assert_eq(result.q, 2)
  assert_true(result.y_loadings[0].length() >= 2)
}

///|
test "pls_canonical_mode" {
  let x = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 4.0],
    [4.0, 5.0],
    [5.0, 6.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.pls_options_full(2, "canonical", true, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.ncomp >= 1)
}

///|
test "spls_variable_selection" {
  let x = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 3.0, 4.0, 5.0, 6.0],
    [3.0, 4.0, 5.0, 6.0, 7.0],
    [4.0, 5.0, 6.0, 7.0, 8.0],
    [5.0, 6.0, 7.0, 8.0, 9.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.spls_options_full(1, [2], [], true, 500, 0.000001)
  let result = @src.run_spls(x, y, opts)
  assert_true(result.ncomp >= 1)
}

///|
test "spls_sparse_component" {
  let x = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
    [4.0, 5.0, 6.0],
    [5.0, 6.0, 7.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.spls_options_full(2, [2, 2], [1, 1], true, 500, 0.000001)
  let result = @src.run_spls(x, y, opts)
  assert_true(result.ncomp >= 1)
}

///|
test "spls_keep_x_parameter" {
  let x = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.spls_options_full(1, [3], [], true, 500, 0.000001)
  let result = @src.run_spls(x, y, opts)
  assert_true(result.ncomp >= 1)
}

///|
test "spls_keep_y_parameter" {
  let x = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
    [4.0, 5.0, 6.0],
    [5.0, 6.0, 7.0],
  ]
  let y = [[1.0, 0.5], [2.0, 1.0], [3.0, 1.5], [4.0, 2.0], [5.0, 2.5]]
  let opts = @src.spls_options_full(1, [], [1], true, 500, 0.000001)
  let result = @src.run_spls(x, y, opts)
  assert_true(result.ncomp >= 1)
}

///|
test "spls_empty_data" {
  let x : Array[Array[Double]] = []
  let y : Array[Array[Double]] = []
  let opts = @src.spls_new()
  let result = @src.run_spls(x, y, opts)
  assert_eq(result.n, 0)
  assert_eq(result.p, 0)
  assert_eq(result.q, 0)
}

///|
test "spls_default_behavior" {
  let x = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
    [4.0, 5.0, 6.0],
    [5.0, 6.0, 7.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.spls_new()
  let result = @src.run_spls(x, y, opts)
  assert_true(result.ncomp >= 1)
}

///|
test "soft_thresholding_operator" {
  let result1 = @src.soft_threshold(3.0, 1.0)
  assert_true((result1 - 2.0).abs() < 1.0e-15)
  let result2 = @src.soft_threshold(-3.0, 1.0)
  assert_true((result2 - (-2.0)).abs() < 1.0e-15)
  let result3 = @src.soft_threshold(0.5, 1.0)
  assert_true(result3.abs() < 1.0e-15)
  let result4 = @src.soft_threshold(0.0, 0.0)
  assert_true(result4.abs() < 1.0e-15)
}

///|
test "diablo_two_blocks" {
  let block1 = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
    [4.0, 5.0, 6.0],
    [5.0, 6.0, 7.0],
  ]
  let block2 = [
    [1.0, 0.5],
    [2.0, 1.0],
    [3.0, 1.5],
    [4.0, 2.0],
    [5.0, 2.5],
  ]
  let blocks = [block1, block2]
  let names = ["block1", "block2"]
  let design = [[0.0, 1.0], [1.0, 0.0]]
  let opts = @src.diablo_options_full(2, [], design, true, 500, 0.000001)
  let result = @src.run_diablo(blocks, names, opts)
  assert_true(result.ncomp >= 1)
  assert_eq(result.n_blocks, 2)
}

///|
test "diablo_variable_selection" {
  let block1 = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0],
  ]
  let block2 = [
    [1.0, 0.5, 0.3],
    [2.0, 1.0, 0.6],
    [3.0, 1.5, 0.9],
    [4.0, 2.0, 1.2],
    [5.0, 2.5, 1.5],
  ]
  let blocks = [block1, block2]
  let names = ["block1", "block2"]
  let design = [[0.0, 1.0], [1.0, 0.0]]
  let opts = @src.diablo_options_full(1, [[2], [1]], design, true, 500, 0.000001)
  let result = @src.run_diablo(blocks, names, opts)
  assert_true(result.ncomp >= 1)
}

///|
test "diablo_default_behavior" {
  let block1 = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 4.0],
    [4.0, 5.0],
  ]
  let block2 = [
    [1.0, 0.5],
    [2.0, 1.0],
    [3.0, 1.5],
    [4.0, 2.0],
  ]
  let blocks = [block1, block2]
  let names = ["block1", "block2"]
  let opts = @src.diablo_new()
  let result = @src.run_diablo(blocks, names, opts)
  assert_true(result.ncomp >= 1)
  assert_eq(result.n_blocks, 2)
}

///|
test "diablo_single_block" {
  let block1 = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
  ]
  let blocks = [block1]
  let names = ["block1"]
  let design = [[0.0]]
  let opts = @src.diablo_options_full(1, [], design, true, 500, 0.000001)
  let result = @src.run_diablo(blocks, names, opts)
  assert_eq(result.n_blocks, 1)
  assert_true(result.ncomp >= 1)
}

///|
test "diablo_empty_blocks" {
  let blocks : Array[Array[Array[Double]]] = []
  let names : Array[String] = []
  let opts = @src.diablo_new()
  let result = @src.run_diablo(blocks, names, opts)
  assert_eq(result.n_blocks, 0)
}

///|
test "diablo_three_blocks" {
  let block1 = [[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]]
  let block2 = [[1.0, 0.5], [2.0, 1.0], [3.0, 1.5]]
  let block3 = [[1.0], [2.0], [3.0]]
  let blocks = [block1, block2, block3]
  let names = ["block1", "block2", "block3"]
  let design = [[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]]
  let opts = @src.diablo_options_full(1, [], design, true, 500, 0.000001)
  let result = @src.run_diablo(blocks, names, opts)
  assert_eq(result.n_blocks, 3)
  assert_true(result.ncomp >= 1)
}

///|
test "pls_options_equality" {
  let opts1 = @src.pls_new()
  let opts2 = @src.pls_new()
  assert_eq(opts1, opts2)
}

///|
test "spls_options_equality" {
  let opts1 = @src.spls_new()
  let opts2 = @src.spls_new()
  assert_eq(opts1, opts2)
}

///|
test "diablo_options_equality" {
  let opts1 = @src.diablo_new()
  let opts2 = @src.diablo_new()
  assert_eq(opts1, opts2)
}

///|
test "pls_result_dimensions_consistency" {
  let x = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0],
    [6.0, 7.0, 8.0, 9.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0], [6.0]]
  let opts = @src.pls_options_full(3, "regression", true, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.x_scores.length() >= 1)
  assert_true(result.x_loadings.length() >= 1)
  assert_true(result.x_weights.length() >= 1)
}

///|
test "pls_with_two_components" {
  let x = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.pls_options_full(2, "regression", true, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.ncomp >= 1)
  assert_true(result.x_scores.length() >= 1)
  assert_true(result.y_scores.length() >= 1)
}

///|
test "diablo_result_structure" {
  let block1 = [[1.0, 2.0], [2.0, 3.0], [3.0, 4.0], [4.0, 5.0]]
  let block2 = [[1.0, 0.5], [2.0, 1.0], [3.0, 1.5], [4.0, 2.0]]
  let blocks = [block1, block2]
  let names = ["block1", "block2"]
  let design = [[0.0, 1.0], [1.0, 0.0]]
  let opts = @src.diablo_options_full(2, [], design, true, 500, 0.000001)
  let result = @src.run_diablo(blocks, names, opts)
  assert_eq(result.block_scores.length(), 2)
  assert_eq(result.block_loadings.length(), 2)
  assert_eq(result.block_weights.length(), 2)
}

///|
test "pls_mode_regression" {
  let x = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 4.0],
    [4.0, 5.0],
    [5.0, 6.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.pls_options_full(2, "regression", true, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.ncomp >= 1)
  assert_true(result.ncomp <= 2)
}

///|
test "pls_mode_invariant_scores" {
  let x = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 4.0],
    [4.0, 5.0],
    [5.0, 6.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.pls_options_full(1, "regression", false, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.x_scores.length() >= 1)
  assert_true(result.x_scores[0].length() > 0)
}

///|
test "diablo_block_names" {
  let block1 = [[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]]
  let block2 = [[1.0, 0.5], [2.0, 1.0], [3.0, 1.5]]
  let blocks = [block1, block2]
  let names = ["rna", "protein"]
  let design = [[0.0, 1.0], [1.0, 0.0]]
  let opts = @src.diablo_options_full(1, [], design, true, 500, 0.000001)
  let result = @src.run_diablo(blocks, names, opts)
  assert_eq(result.block_names.length(), 2)
  assert_eq(result.block_names[0], "rna")
  assert_eq(result.block_names[1], "protein")
}

///|
test "pls_result_has_explained_variance" {
  let x = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [3.0, 4.0, 5.0],
    [4.0, 5.0, 6.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0]]
  let opts = @src.pls_options_full(2, "regression", true, 500, 0.000001)
  let result = @src.run_pls(x, y, opts)
  assert_true(result.explained_variance.length() >= 1)
}

///|
test "spls_selected_variables" {
  let x = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0],
  ]
  let y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
  let opts = @src.spls_options_full(1, [2], [1], true, 500, 0.000001)
  let result = @src.run_spls(x, y, opts)
  assert_true(result.ncomp >= 1)
  assert_true(result.x_selected.length() >= 1)
  assert_true(result.y_selected.length() >= 1)
}