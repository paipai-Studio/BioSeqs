///|
test "gsva_create_example_data" {
  let data = @src.gsva_create_example_data()
  assert_eq(data.expression.length(), 100)
  assert_eq(data.gene_names.length(), 100)
  assert_eq(data.sample_names.length(), 10)
}

///|
test "gsva_create_example_gene_sets" {
  let gene_sets = @src.gsva_create_example_gene_sets()
  assert_eq(gene_sets.length(), 5)
}

///|
test "gsva_params_new" {
  let params = @src.GSVAParams::new()
  assert_eq(params.method, "ssgsea")
}

///|
test "gsva_run" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  assert_eq(scores.scores.length(), 5)
  assert_eq(scores.scores[0].length(), 10)
}

///|
test "gsva_score_ssgsea" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let ssgsea_scores = @src.gsva_score_ssgsea(
    data.expression,
    gene_sets[0],
    data.gene_names,
  )
  assert_eq(ssgsea_scores.length(), 10)
}

///|
test "gsva_score_zscore" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let zscore_scores = @src.gsva_score_zscore(
    data.expression,
    gene_sets[0],
    data.gene_names,
  )
  assert_eq(zscore_scores.length(), 10)
}

///|
test "gsva_score_plage" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let plage_scores = @src.gsva_score_plage(
    data.expression,
    gene_sets[0],
    data.gene_names,
  )
  assert_eq(plage_scores.length(), 10)
}

///|
test "gsva_score_gsva" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let gsva_scores = @src.gsva_score_gsva(
    data.expression,
    gene_sets[0],
    data.gene_names,
    "Gaussian",
  )
  assert_eq(gsva_scores.length(), 10)
}

///|
test "gsva_rank_genes" {
  let data = @src.gsva_create_example_data()
  let ranks = @src.gsva_rank_genes(data.expression, "mean")
  assert_eq(ranks.length(), 10)
  assert_eq(ranks[0].length(), 100)
}

///|
test "gsva_filter_gene_sets" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let filtered = @src.gsva_filter_gene_sets(gene_sets, data.gene_names, 5, 100)
  assert_eq(filtered.length(), 5)
}

///|
test "gsva_statistics" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let stats = @src.gsva_statistics(scores.scores)
  assert_eq(stats.mean_score.length(), 5)
}

///|
test "gsva_summarize" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let summary = @src.gsva_summarize(scores)
  assert_true(summary.length() > 0)
}

///|
test "gsva_differential_analysis" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
  let diff = @src.gsva_differential_analysis(scores.scores, groups)
  assert_eq(diff.length(), 5)
}

///|
test "gsva_permutation_test" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let groups = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]
  let perm_result = @src.gsva_permutation_test(scores.scores, groups, 10)
  assert_eq(perm_result.length(), 5)
}

///|
test "gsva_utils" {
  let arr = [1.0, 2.0, 3.0, 4.0, 5.0]
  let mean_val = @src.gsva_mean(arr)
  assert_eq(mean_val, 3.0)

  let var_val = @src.gsva_variance(arr)
  assert_true(var_val > 0.0)

  let std_val = @src.gsva_std(arr)
  assert_true(std_val > 0.0)

  let median_val = @src.gsva_median(arr)
  assert_eq(median_val, 3.0)

  let sqrt_val = @src.gsva_sqrt(16.0)
  assert_eq(sqrt_val, 4.0)
}

///|
test "gsva_enrichment_map" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let map = @src.gsva_enrichment_map(scores, gene_sets, 0.1)
  assert_true(map.length() > 0)
}

///|
test "gsva_phenotype_correlation" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let phenotypes = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let corr = @src.gsva_phenotype_correlation(scores, phenotypes)
  assert_eq(corr.length(), 5)
}

///|
test "gsva_phenotype_report" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let phenotypes = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let report = @src.gsva_phenotype_report(scores, phenotypes, "Age")
  assert_true(report.length() > 0)
}

///|
test "gsva_survival_analysis" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let survival_time = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
  let event = [1, 1, 1, 0, 0, 1, 1, 0, 1, 0]
  let surv = @src.gsva_survival_analysis(scores, survival_time, event)
  assert_eq(surv.length(), 5)
}

///|
test "gsva_survival_report" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let survival_time = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
  let event = [1, 1, 1, 0, 0, 1, 1, 0, 1, 0]
  let report = @src.gsva_survival_report(scores, survival_time, event)
  assert_true(report.length() > 0)
}

///|
test "gsva_score_distribution" {
  let data = @src.gsva_create_example_data()
  let gene_sets = @src.gsva_create_example_gene_sets()
  let params = @src.GSVAParams::new()
  let scores = @src.gsva_run(data, gene_sets, params)
  let dist = @src.gsva_score_distribution(scores)
  assert_true(dist.length() > 0)
}
