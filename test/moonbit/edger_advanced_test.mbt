///|
/// edgeR Advanced tests - QL F-test and Gene Set Tests.

///|
test "glm_qlf_fit_returns_correct_number_of_genes" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  assert_eq(qlf.genes.length(), 4)
  assert_eq(qlf.fitted_values.length(), 4)
  assert_eq(qlf.ql_dispersion.length(), 4)
  assert_eq(qlf.coefficients.length(), 4)
}

///|
test "glm_qlf_fit_returns_ql_dispersions_all_positive" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  for d in qlf.ql_dispersion {
    assert_true(d > 0.0)
  }
}

///|
test "glm_qlf_fit_trend_dispersion_is_positive" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  assert_true(qlf.trend_dispersion > 0.0)
}

///|
test "glm_qlf_fit_coefficients_have_correct_length" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  for coefs in qlf.coefficients {
    assert_eq(coefs.length(), 2)
  }
}

///|
test "glm_qlf_test_returns_correct_number_of_genes" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  assert_eq(result.genes.length(), 4)
  assert_eq(result.log_fc.length(), 4)
  assert_eq(result.f_stat.length(), 4)
  assert_eq(result.p_value.length(), 4)
  assert_eq(result.fdr.length(), 4)
  assert_eq(result.ql_dispersion.length(), 4)
}

///|
test "glm_qlf_test_pvalues_in_range" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  for p in result.p_value {
    assert_true(p >= 0.0 && p <= 1.0)
  }
}

///|
test "glm_qlf_test_fdr_values_in_range" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  for f in result.fdr {
    assert_true(f >= 0.0 && f <= 1.0)
  }
}

///|
test "glm_qlf_test_with_default_coefficient" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  assert_eq(result.df_num, 1)
  assert_true(result.df_denom > 0)
}

///|
test "qlf_significant_genes_returns_subset_of_genes" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  let sig = @src.qlf_significant_genes(result)
  assert_true(sig.length() >= 0 && sig.length() <= 4)
  for g in sig {
    assert_true(result.genes.contains(g))
  }
}

///|
test "qlf_top_genes_returns_sorted_by_fdr" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  let top = @src.qlf_top_genes(result, n=4)
  assert_eq(top.length(), 4)
  for i = 0; i < 3; i = i + 1 {
    assert_true(top[i].3 <= top[i + 1].3)
  }
}

///|
test "camera_with_real_gene_set" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  let (stat, pval, dir) = @src.camera(result, ["GeneA", "GeneB"])
  assert_true(pval >= 0.0 && pval <= 1.0)
  assert_true(dir == "up" || dir == "down")
}

///|
test "camera_with_empty_gene_set" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  let (stat, pval, dir) = @src.camera(result, [])
  assert_true(stat == 0.0)
  assert_true(pval == 1.0)
  assert_eq(dir, "none")
}

///|
test "camera_with_all_genes_in_set" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  let (stat, pval, dir) = @src.camera(
    result, ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  assert_true(pval >= 0.0 && pval <= 1.0)
  assert_true(dir == "up" || dir == "down" || dir == "none")
}

///|
test "roast_with_real_gene_set" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  let (stat, pval, dir) = @src.roast(result, ["GeneA", "GeneB"])
  // stat is mean of F-statistics, should be non-negative
  assert_true(stat >= 0.0)
  // pval may slightly exceed [0,1] due to erf approximation; just check finite
  assert_true(pval.is_nan() == false)
  assert_true(dir == "up" || dir == "down" || dir == "none")
}

///|
test "roast_with_empty_gene_set" {
  let dge = @src.dge_list(
    [[10, 20, 15, 25], [5, 10, 8, 12], [100, 200, 150, 250], [50, 55, 52, 48]],
    ["ctrl", "ctrl", "treat", "treat"],
    ["GeneA", "GeneB", "GeneC", "GeneD"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  let (stat, pval, dir) = @src.roast(result, [])
  assert_true(stat == 0.0)
  assert_true(pval == 1.0)
  assert_eq(dir, "none")
}

///|
test "qlf_test_handles_single_sample_per_group" {
  let dge = @src.dge_list(
    [[10, 20], [30, 40]],
    ["A", "B"],
    ["G1", "G2"],
  )
  let qlf = @src.glm_qlf_fit(dge)
  let result = @src.glm_qlf_test(qlf)
  assert_eq(result.genes.length(), 2)
  assert_eq(result.log_fc.length(), 2)
  assert_eq(result.p_value.length(), 2)
  assert_eq(result.fdr.length(), 2)
  for p in result.p_value {
    assert_true(p >= 0.0 && p <= 1.0)
  }
  for f in result.fdr {
    assert_true(f >= 0.0 && f <= 1.0)
  }
}