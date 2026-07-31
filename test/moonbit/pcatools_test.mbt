///|
test "biplot_options_new" {
  let opts = @src.BiplotOptions::new()
  assert_eq(opts.pc_x, 1)
  assert_eq(opts.pc_y, 2)
  assert_eq(opts.width, 60)
  assert_eq(opts.height, 20)
  assert_eq(opts.max_loadings, 8)
}

///|
test "pcatools_run_pca_simple_2d" {
  let data = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 4.0],
    [4.0, 5.0],
  ]
  let res = @src.pcatools_run_pca(data, 2, false)
  assert_eq(res.n_samples, 4)
  assert_eq(res.n_variables, 2)
  assert_eq(res.n_components, 2)
  assert_eq(res.scores.length(), 4)
  assert_eq(res.scores[0].length(), 2)
  assert_eq(res.eigenvalues.length(), 2)
  assert_true(res.variance_explained[0] > 0.9)
  let last_idx = res.cumulative_variance.length() - 1
  assert_true(res.cumulative_variance[last_idx] > 0.99)
}

///|
test "pcatools_run_pca_ncomponents_limit" {
  let data = [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
  ]
  let res = @src.pcatools_run_pca(data, 10, false)
  // 3 samples, 4 variables => p=4, so limited to 4
  assert_eq(res.n_components, 4)
}

///|
test "pcatools_run_pca_scale_true" {
  let data = [
    [1.0, 100.0],
    [2.0, 200.0],
    [3.0, 300.0],
    [4.0, 400.0],
  ]
  let res = @src.pcatools_run_pca(data, 2, true)
  assert_eq(res.used_scaling, true)
  assert_eq(res.scale.length(), 2)
  assert_eq(res.center.length(), 2)
  assert_true(res.scale[0] > 0.0)
  assert_true(res.scale[1] > 0.0)
}

///|
test "pcatools_run_pca_empty" {
  let data : Array[Array[Double]] = []
  let res = @src.pcatools_run_pca(data, 2, false)
  assert_eq(res.n_samples, 0)
  assert_eq(res.scores.length(), 0)
}

///|
test "scree_plot_ascii_returns_string" {
  let data = [
    [1.0, 2.0, 3.0],
    [2.0, 4.0, 5.0],
    [3.0, 5.0, 7.0],
    [4.0, 6.0, 9.0],
    [5.0, 7.0, 10.0],
    [6.0, 8.0, 11.0],
  ]
  let res = @src.pcatools_run_pca(data, 3, false)
  let plot = @src.scree_plot_ascii(res, height=8, width_per_pc=3)
  assert_true(plot.length() > 0)
  assert_true(plot.contains("PCs"))
}

///|
test "pca_biplot_ascii_returns_string" {
  let data = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 5.0],
    [4.0, 4.0],
    [5.0, 6.0],
  ]
  let res = @src.pcatools_run_pca(data, 2, false)
  let opts = @src.BiplotOptions::new()
  let plot = @src.pca_biplot_ascii(res, opts=opts)
  assert_true(plot.length() > 0)
  assert_true(plot.contains("PCA biplot"))
  assert_true(plot.contains("PC1"))
  assert_true(plot.contains("PC2"))
}

///|
test "find_pca_outliers_basic" {
  let data = [
    [0.0, 0.0],
    [0.1, 0.1],
    [-0.1, 0.1],
    [0.1, -0.1],
    [-0.1, -0.1],
    [0.05, 0.02],
    [-0.05, 0.03],
    [10.0, 10.0],
  ]
  let res = @src.pcatools_run_pca(data, 2, false)
  let out = @src.find_pca_outliers(res, n_pcs=2, quantile_cutoff=0.95)
  assert_eq(out.n_tested, 8)
  assert_true(out.indices.length() >= 1)
  let mut found7 = false
  for i in out.indices {
    if i == 7 { found7 = true }
  }
  assert_true(found7)
}

///|
test "find_pca_outliers_cutoff_positive" {
  let data = [
    [1.0, 1.0], [1.1, 1.0], [0.9, 1.0],
    [1.0, 1.1], [1.0, 0.9],
    [20.0, 20.0],
  ]
  let res = @src.pcatools_run_pca(data, 2, false)
  let out = @src.find_pca_outliers(res, n_pcs=2, quantile_cutoff=0.9)
  assert_true(out.cutoff > 0.0)
  assert_eq(out.n_tested, 6)
  // the extreme point [20,20] should be flagged
  assert_true(out.indices.length() >= 1)
}

///|
test "variable_correlations_shape" {
  let data = [
    [1.0, 2.0, 3.0],
    [2.0, 4.0, 6.0],
    [3.0, 6.0, 9.0],
    [4.0, 8.0, 12.0],
  ]
  let res = @src.pcatools_run_pca(data, 2, true)
  let corr = @src.variable_correlations(data, res)
  assert_eq(corr.length(), res.n_variables)
  if corr.length() > 0 {
    assert_eq(corr[0].length(), res.n_components)
    // correlations should be in [-1, 1] approx
    let mut i = 0
    while i < corr.length() {
      let mut j = 0
      while j < corr[i].length() {
        assert_true(corr[i][j] >= -1.0001 && corr[i][j] <= 1.0001)
        j = j + 1
      }
      i = i + 1
    }
  }
}

///|
test "pcatools_summary_contains_keywords" {
  let data = [
    [1.0, 2.0],
    [2.0, 3.0],
    [3.0, 4.0],
  ]
  let res = @src.pcatools_run_pca(data, 2, false)
  let s = @src.pcatools_summary(res)
  assert_true(s.contains("PCAResult"))
  assert_true(s.contains("samples"))
  assert_true(s.contains("variables"))
  assert_true(s.contains("components"))
}
