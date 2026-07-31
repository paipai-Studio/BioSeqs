///|
/// Tests for Bioconductor baySeq module - Bayesian differential expression.

test "estimate_dispersion_basic" {
  // Two groups with different means
  let a = [100.0, 120.0, 110.0, 90.0]
  let b = [200.0, 220.0, 210.0, 190.0]
  let disp = @src.bayseq_estimate_dispersion_gene(a, b)
  assert_true(disp > 0.0)
  assert_true(disp <= 100.0)
}

test "estimate_dispersion_low_counts" {
  let a = [1.0, 2.0, 1.5, 1.2]
  let b = [1.0, 1.8, 1.3, 1.1]
  let disp = @src.bayseq_estimate_dispersion_gene(a, b)
  // With very low counts, should default to reasonable value
  assert_true(disp > 0.0)
}

test "estimate_dispersion_single_sample" {
  let a = [100.0]
  let b = [200.0]
  let disp = @src.bayseq_estimate_dispersion_gene(a, b)
  assert_true(disp == 1.0)
}

test "estimate_prior_basic" {
  let (counts, groups, _) = @src.bayseq_sample_data()
  let (shape, rate) = @src.bayseq_estimate_prior(counts, groups)
  assert_true(shape > 0.0)
  assert_true(rate > 0.0)
}

test "estimate_prior_empty" {
  let counts : Array[Array[Double]] = []
  let (shape, rate) = @src.bayseq_estimate_prior(counts, [])
  assert_true(shape == 2.0)
  assert_true(rate == 0.5)
}

test "log_likelihood_ratio_de" {
  let a = [100.0, 120.0, 110.0, 90.0]
  let b = [500.0, 520.0, 510.0, 490.0]
  let disp = 10.0
  let llr = @src.bayseq_log_likelihood_ratio(a, b, disp)
  // Higher LL = more evidence for DE
  assert_true(llr > 0.0)
}

test "log_likelihood_ratio_same" {
  let a = [100.0, 120.0, 110.0, 90.0]
  let b = [100.0, 120.0, 110.0, 90.0]
  let disp = 10.0
  let llr = @src.bayseq_log_likelihood_ratio(a, b, disp)
  // Same distribution -> LL should be close to 0
  assert_true(llr >= -1.0 && llr <= 1.0)
}

test "lgamma_positive" {
  let val = @src.bayseq_lgamma(5.0)
  assert_true(val > 0.0)
}

test "lgamma_one" {
  // lgamma(1) = 0 since gamma(1) = 1, log(1) = 0
  let val = @src.bayseq_lgamma(1.0)
  assert_true(val >= -0.1 && val <= 0.1)
}

test "lgamma_large" {
  let val = @src.bayseq_lgamma(100.0)
  // lgamma(100) ≈ 363.739
  assert_true(val > 100.0)
}

test "nb_log_likelihood_zero_mean" {
  let ll = @src.bayseq_nb_log_likelihood([1.0, 2.0, 3.0], 0.0, 10.0)
  assert_true(ll == 0.0)
}

test "test_basic" {
  let (counts, groups, names) = @src.bayseq_sample_data()
  let result = @src.bayseq_test(counts, groups, names)
  assert_eq(result.gene_results.length(), 30)
  assert_true(result.n_de_genes > 0)
  assert_true(result.n_de_genes <= 30)
}

test "test_lfc_calculation" {
  let (counts, groups, names) = @src.bayseq_sample_data()
  let result = @src.bayseq_test(counts, groups, names)
  // Genes with higher expression in group 1 (first 8) should have positive LFC
  let g0 = result.gene_results[0]
  assert_true(g0.log_fold_change > 0.0)
}

test "test_posterior_prob" {
  let (counts, groups, names) = @src.bayseq_sample_data()
  let result = @src.bayseq_test(counts, groups, names)
  let g0 = result.gene_results[0]
  assert_true(g0.posterior_prob_de >= 0.0 && g0.posterior_prob_de <= 1.0)
}

test "get_de_genes" {
  let (counts, groups, names) = @src.bayseq_sample_data()
  let result = @src.bayseq_test(counts, groups, names)
  let de_genes = @src.bayseq_get_de_genes(result)
  assert_eq(de_genes.length(), result.n_de_genes)
}

test "get_top_de" {
  let (counts, groups, names) = @src.bayseq_sample_data()
  let result = @src.bayseq_test(counts, groups, names)
  let top = @src.bayseq_get_top_de(result, 5)
  assert_eq(top.length(), 5)
}

test "get_top_de_fewer_than_n" {
  let (counts, groups, names) = @src.bayseq_sample_data()
  let result = @src.bayseq_test(counts, groups, names)
  let top = @src.bayseq_get_top_de(result, 100)
  assert_eq(top.length(), 30)
}

test "sample_data_dimensions" {
  let (counts, groups, names) = @src.bayseq_sample_data()
  assert_eq(counts.length(), 30)
  assert_eq(counts[0].length(), 8)
  assert_eq(groups.length(), 8)
  assert_eq(names.length(), 30)
}

test "summary" {
  let (counts, groups, names) = @src.bayseq_sample_data()
  let result = @src.bayseq_test(counts, groups, names)
  let summary = @src.bayseq_summary(result)
  assert_true(summary.length() > 0)
}
