///|
test "scnorm_params" {
  let params = @src.SCnormParams::new()
  assert_true(params.quantile > 0.0)
  assert_true(params.quantile < 1.0)
  assert_true(params.filter_genes)
  assert_true(params.robust)
}

///|
test "scnorm_library_sizes" {
  let counts : Array[Array[Double]] = [
    [10.0, 20.0, 30.0],
    [5.0, 15.0, 25.0],
    [50.0, 100.0, 150.0]
  ]
  let sizes = @src.sc_norm_library_sizes(counts)
  assert_eq(sizes.length(), 3)
  assert_true(sizes[0] > 0.0)
  assert_true(sizes[2] > sizes[0])
}

///|
test "scnorm_quantile_regression" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let fit = @src.sc_norm_quantile_regression(x, y, 0.5)
  assert_true(fit.slope > 0.0)
  assert_eq(fit.fitted.length(), 5)
  assert_eq(fit.residuals.length(), 5)
}

///|
test "scnorm_loess_smooth" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [2.0, 3.5, 5.0, 6.5, 8.0, 9.5, 11.0, 12.5, 14.0, 15.5]
  let smoothed = @src.sc_norm_loess_smooth(x, y, 0.5)
  assert_eq(smoothed.length(), 10)
}

///|
test "scnorm_spline_smooth" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let y = [2.0, 3.0, 5.0, 7.0, 8.5, 9.5, 11.0, 12.0]
  let smoothed = @src.sc_norm_spline_smooth(x, y, 3)
  assert_true(smoothed.length() > 0)
}

///|
test "scnorm_normalize_gene" {
  let counts = [10.0, 20.0, 30.0, 40.0, 50.0]
  let lib_sizes = [100.0, 200.0, 300.0, 400.0, 500.0]
  let params = @src.SCnormParams::new()
  let result = @src.sc_norm_normalize_gene("gene1", counts, lib_sizes, params)
  assert_eq(result.gene_id, "gene1")
  assert_eq(result.norm_counts.length(), 5)
}

///|
test "scnorm_full_run" {
  let count_matrix : Array[Array[Double]] = [
    [10.0, 20.0, 30.0],
    [50.0, 60.0, 70.0],
    [100.0, 200.0, 300.0]
  ]
  let gene_ids = ["gene1", "gene2", "gene3"]
  let cell_ids = ["cell1", "cell2", "cell3"]
  let params = @src.SCnormParams::new()
  let result = @src.sc_norm_run(count_matrix, gene_ids, cell_ids, params)
  assert_true(result.sc_n_cells() > 0)
}

///|
test "scnorm_result_accessors" {
  let count_matrix : Array[Array[Double]] = [
    [10.0, 20.0, 30.0],
    [50.0, 60.0, 70.0]
  ]
  let gene_ids = ["gene1", "gene2"]
  let cell_ids = ["cell1", "cell2", "cell3"]
  let params = @src.SCnormParams::new()
  let result = @src.sc_norm_run(count_matrix, gene_ids, cell_ids, params)
  assert_true(result.sc_norm_get_count(0, 0) > 0.0)
  assert_eq(result.sc_norm_get_gene(0).length(), 3)
}
