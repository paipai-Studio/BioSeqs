///|
/// Tests for MatrixEQTL module.

///|
test "matrixeqtl_linear_regression_perfect_linear" {
  // y ≈ 2x + 3 with tiny noise → beta1 ≈ 2.0, high t-stat, low p-value
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let y = [5.01, 7.02, 8.99, 11.01, 12.98, 15.01, 16.99, 19.01]
  let (beta, t_stat, p_value) = @src.matrixeqtl_linear_regression(x, y)

  // beta should be close to 2.0
  assert_true((beta - 2.0).abs() < 0.05)
  // t-stat should be large (strong linear relationship)
  assert_true(t_stat.abs() > 50.0)
  // p-value should be very small
  assert_true(p_value < 0.001)
}

///|
test "matrixeqtl_linear_regression_uncorrelated" {
  // x and y have no correlation → beta ≈ 0, low |t-stat|
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let y = [8.0, 3.0, 5.0, 1.0, 7.0, 4.0, 6.0, 2.0]
  let (beta, t_stat, p_value) = @src.matrixeqtl_linear_regression(x, y)

  // beta should be near 0
  assert_true(beta.abs() < 0.5)
  // |t-stat| should be relatively low
  assert_true(t_stat.abs() < 3.0)
  // p-value should not be extremely small
  assert_true(p_value > 0.05)
}

///|
test "matrixeqtl_linear_regression_negative_slope" {
  // y ≈ -1.5x + 10
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let y = [8.51, 7.01, 5.49, 3.99, 2.51, 1.01, -0.49, -1.99]
  let (beta, t_stat, p_value) = @src.matrixeqtl_linear_regression(x, y)

  assert_true((beta - (-1.5)).abs() < 0.05)
  assert_true(t_stat < -50.0)
  assert_true(p_value < 0.001)
}

///|
test "matrixeqtl_linear_regression_too_few_samples" {
  // Less than 3 samples → returns (0.0, 0.0, 1.0)
  let x = [1.0, 2.0]
  let y = [3.0, 5.0]
  let (beta, t_stat, p_value) = @src.matrixeqtl_linear_regression(x, y)

  assert_true(beta == 0.0)
  assert_true(t_stat == 0.0)
  assert_true(p_value == 1.0)
}

///|
test "matrixeqtl_linear_regression_constant_x" {
  // Constant x (no variation) → returns (0.0, 0.0, 1.0)
  let x = [3.0, 3.0, 3.0, 3.0, 3.0]
  let y = [1.0, 2.0, 3.0, 4.0, 5.0]
  let (beta, t_stat, p_value) = @src.matrixeqtl_linear_regression(x, y)

  assert_true(beta == 0.0)
  assert_true(t_stat == 0.0)
  assert_true(p_value == 1.0)
}

///|
test "matrixeqtl_bh_fdr_monotonicity" {
  // BH FDR should produce monotonically non-decreasing adjusted p-values
  // when sorted by original p-value
  let p_values = [0.01, 0.04, 0.05, 0.2, 0.5]
  let adj = @src.matrixeqtl_bh_fdr(p_values)

  assert_eq(adj.length(), 5)
  // The adjusted p-values corresponding to sorted original p-values
  // should be non-decreasing: adj[0] <= adj[1] <= ... <= adj[4]
  // Since original is already sorted, check monotonicity directly
  assert_true(adj[0] <= adj[1])
  assert_true(adj[1] <= adj[2])
  assert_true(adj[2] <= adj[3])
  assert_true(adj[3] <= adj[4])
}

///|
test "matrixeqtl_bh_fdr_caps_at_one" {
  // All adjusted p-values should be capped at 1.0
  let p_values = [0.8, 0.9, 0.95, 0.99, 1.0]
  let adj = @src.matrixeqtl_bh_fdr(p_values)

  for p in adj {
    assert_true(p <= 1.0)
  }
}

///|
test "matrixeqtl_bh_fdr_empty" {
  // Empty input → empty output
  let p_values : Array[Double] = []
  let adj = @src.matrixeqtl_bh_fdr(p_values)
  assert_eq(adj.length(), 0)
}

///|
test "matrixeqtl_bh_fdr_single" {
  // Single p-value → adjusted = min(p * 1 / 1, 1.0) = p
  let p_values = [0.03]
  let adj = @src.matrixeqtl_bh_fdr(p_values)
  assert_eq(adj.length(), 1)
  assert_true((adj[0] - 0.03).abs() < 1.0e-10)
}

///|
test "matrixeqtl_bh_fdr_unsorted_input" {
  // BH FDR should work correctly even with unsorted input
  let p_values = [0.5, 0.01, 0.2, 0.04, 0.05]
  let adj = @src.matrixeqtl_bh_fdr(p_values)

  assert_eq(adj.length(), 5)
  // All adjusted values should be <= 1.0
  for p in adj {
    assert_true(p <= 1.0)
  }
  // Adjusted values should be >= original values
  let mut i = 0
  while i < 5 {
    assert_true(adj[i] >= p_values[i])
    i = i + 1
  }
}

///|
test "matrixeqtl_input_new" {
  // MatrixEQTLInput::new creates proper input
  let snps : Array[Array[Double]] = [[0.0, 1.0, 2.0], [1.0, 1.0, 0.0]]
  let gene_expression : Array[Array[Double]] = [
    [5.0, 6.0, 7.0],
    [3.0, 3.5, 4.0],
  ]
  let snp_names = ["SNP1", "SNP2"]
  let gene_names = ["GeneA", "GeneB"]
  let sample_names = ["S1", "S2", "S3"]

  let input = @src.MatrixEQTLInput::new(
    snps,
    gene_expression,
    snp_names,
    gene_names,
    sample_names,
  )

  assert_eq(input.snps.length(), 2)
  assert_eq(input.gene_expression.length(), 2)
  assert_eq(input.snp_names.length(), 2)
  assert_eq(input.gene_names.length(), 2)
  assert_eq(input.sample_names.length(), 3)
  assert_true(input.snp_names[0] == "SNP1")
  assert_true(input.gene_names[1] == "GeneB")
}

///|
test "matrixeqtl_run_synthetic" {
  // matrixeqtl_run on small synthetic data returns results
  let snps : Array[Array[Double]] = [
    // SNP1: strongly correlated with GeneA
    [0.0, 1.0, 1.0, 2.0, 0.0, 1.0, 2.0, 1.0],
    // SNP2: not strongly correlated
    [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
  ]
  let gene_expression : Array[Array[Double]] = [
    // GeneA: correlated with SNP1 (≈ 2*SNP + 5)
    [5.1, 6.3, 6.5, 7.8, 4.9, 6.1, 7.9, 6.4],
    // GeneB: no strong correlation
    [3.2, 3.5, 3.1, 3.8, 3.0, 3.4, 3.6, 3.3],
  ]
  let snp_names = ["SNP1", "SNP2"]
  let gene_names = ["GeneA", "GeneB"]
  let sample_names = [
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
  ]

  let input = @src.MatrixEQTLInput::new(
    snps,
    gene_expression,
    snp_names,
    gene_names,
    sample_names,
  )
  let output = @src.matrixeqtl_run(input)

  // Should have performed 2 SNPs × 2 Genes = 4 tests
  assert_eq(output.n_tests, 4)
  // Should have some results (at least SNP1-GeneA should be significant)
  assert_true(output.results.length() > 0)
  // First result should have the smallest p-value (likely SNP1-GeneA)
  let first = output.results[0]
  assert_true(first.p_value < 0.05)
}

///|
test "matrixeqtl_run_with_threshold" {
  // Using a stricter threshold should return fewer results
  let snps : Array[Array[Double]] = [
    [0.0, 1.0, 1.0, 2.0, 0.0, 1.0, 2.0, 1.0],
    [1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
  ]
  let gene_expression : Array[Array[Double]] = [
    [5.1, 6.3, 6.5, 7.8, 4.9, 6.1, 7.9, 6.4],
    [3.2, 3.5, 3.1, 3.8, 3.0, 3.4, 3.6, 3.3],
  ]
  let snp_names = ["SNP1", "SNP2"]
  let gene_names = ["GeneA", "GeneB"]
  let sample_names = [
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
  ]

  let input = @src.MatrixEQTLInput::new(
    snps,
    gene_expression,
    snp_names,
    gene_names,
    sample_names,
  )
  let output_strict = @src.matrixeqtl_run(input, pv_threshold=0.001)
  let output_relaxed = @src.matrixeqtl_run(input, pv_threshold=0.5)

  // Relaxed threshold should have at least as many results
  assert_true(output_relaxed.results.length() >= output_strict.results.length())
}

///|
test "matrixeqtl_run_sample_input" {
  // Test using the sample_matrixeqtl_input function
  let input = @src.sample_matrixeqtl_input()
  let output = @src.matrixeqtl_run(input)

  // 5 SNPs × 4 Genes = 20 tests
  assert_eq(output.n_tests, 20)
  // Should have some significant results
  assert_true(output.results.length() > 0)
}

///|
test "matrixeqtl_run_empty_input" {
  // Empty data should return empty results
  let snps : Array[Array[Double]] = []
  let gene_expression : Array[Array[Double]] = []
  let snp_names : Array[String] = []
  let gene_names : Array[String] = []
  let sample_names : Array[String] = []

  let input = @src.MatrixEQTLInput::new(
    snps,
    gene_expression,
    snp_names,
    gene_names,
    sample_names,
  )
  let output = @src.matrixeqtl_run(input)

  assert_eq(output.n_tests, 0)
  assert_eq(output.results.length(), 0)
}

///|
test "matrixeqtl_result_struct" {
  // MatrixEQTLResult::new creates proper result
  let result = @src.MatrixEQTLResult::new(
    "rs123", "BRCA1", 1.5, 4.2, 0.001, 0.01,
  )
  assert_true(result.snp_name == "rs123")
  assert_true(result.gene_name == "BRCA1")
  assert_true((result.beta - 1.5).abs() < 1.0e-10)
  assert_true((result.t_stat - 4.2).abs() < 1.0e-10)
  assert_true((result.p_value - 0.001).abs() < 1.0e-10)
  assert_true((result.fdr - 0.01).abs() < 1.0e-10)
}

///|
test "matrixeqtl_output_struct" {
  // MatrixEQTLOutput::new creates proper output
  let results : Array[@src.MatrixEQTLResult] = []
  let output = @src.MatrixEQTLOutput::new(results, 10)
  assert_eq(output.n_tests, 10)
  assert_eq(output.results.length(), 0)
}
