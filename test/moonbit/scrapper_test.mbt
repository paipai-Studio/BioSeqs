///|
/// Tests for Bioconductor scrapper-inspired single-cell preprocessing.

///|
fn scrapper_test_assert_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort("values differ beyond tolerance")
  }
}

///|
fn scrapper_test_sce() -> @src.SingleCellExperiment {
  @src.SCEBuilder::new()
  |> @src.SCEBuilder::add_assay("counts", [
    [3.0, 6.0, 9.0],
    [1.0, 2.0, 3.0],
    [0.0, 1.0, 0.0],
  ])
  |> @src.SCEBuilder::add_row_data("symbol", ["A", "B", "MT-C"])
  |> @src.SCEBuilder::add_col_data("batch", ["x", "x", "y"])
  |> @src.SCEBuilder::add_reduced_dim("PCA", [
    [0.1, 0.2],
    [0.3, 0.4],
    [0.5, 0.6],
  ])
  |> @src.SCEBuilder::build
}

///|
test "scrapper: sample count matrix uses feature by cell orientation" {
  let counts = @src.scrapper_sample_counts()
  assert_eq(counts.length(), 6)
  assert_eq(counts[0].length(), 6)
  assert_eq(counts[0][3], 60.0)
}

///|
test "scrapper: computes RNA QC sums detected genes and subset proportions" {
  let metrics = @src.scrapper_compute_rna_qc_metrics(
    [[10.0, 0.0, 1.0], [0.0, 0.0, 1.0], [5.0, 0.0, 0.0]],
    subsets=[@src.ScrapperNamedSubset::new("mito", [0, 2])],
  ) catch {
    _ => abort("valid RNA QC input should be accepted")
  }
  assert_eq(metrics.sums, [15.0, 0.0, 2.0])
  assert_eq(metrics.detected, [2, 0, 2])
  assert_eq(metrics.subset_names, ["mito"])
  assert_eq(metrics.subset_proportions[0], [1.0, 0.0, 0.5])
}

///|
test "scrapper: RNA QC detection limit uses strict greater than comparison" {
  let metrics = @src.scrapper_compute_rna_qc_metrics(
    [[2.0, 1.0], [1.0, 1.0], [3.0, 0.0]],
    detection_limit=1.0,
  ) catch {
    _ => abort("valid detection limit should be accepted")
  }
  assert_eq(metrics.sums, [6.0, 2.0])
  assert_eq(metrics.detected, [2, 0])
}

///|
test "scrapper: RNA QC rejects malformed count matrices" {
  let ragged = try {
    ignore(@src.scrapper_compute_rna_qc_metrics([[1.0], [2.0, 3.0]]))
    false
  } catch {
    ScrapperError(_) => true
  }
  let negative = try {
    ignore(@src.scrapper_compute_rna_qc_metrics([[1.0, -1.0]]))
    false
  } catch {
    ScrapperError(_) => true
  }
  let non_finite = try {
    ignore(@src.scrapper_compute_rna_qc_metrics([[1.0, @double.not_a_number]]))
    false
  } catch {
    ScrapperError(_) => true
  }
  assert_true(ragged)
  assert_true(negative)
  assert_true(non_finite)
}

///|
test "scrapper: RNA QC validates named feature subsets" {
  let duplicate = try {
    ignore(
      @src.scrapper_compute_rna_qc_metrics([[1.0], [2.0]], subsets=[
        @src.ScrapperNamedSubset::new("mito", [0]),
        @src.ScrapperNamedSubset::new("mito", [1]),
      ]),
    )
    false
  } catch {
    ScrapperError(_) => true
  }
  let out_of_bounds = try {
    ignore(
      @src.scrapper_compute_rna_qc_metrics([[1.0]], subsets=[
        @src.ScrapperNamedSubset::new("mito", [1]),
      ]),
    )
    false
  } catch {
    ScrapperError(_) => true
  }
  assert_true(duplicate)
  assert_true(out_of_bounds)
}

///|
test "scrapper: unblocked MAD thresholds filter low-library cells" {
  let metrics = @src.scrapper_compute_rna_qc_metrics([[1.0, 10.0, 100.0]]) catch {
    _ => abort("valid RNA QC input should be accepted")
  }
  let thresholds = @src.scrapper_suggest_rna_qc_thresholds(
    metrics,
    sum_num_mads=0.0,
    detected_num_mads=0.0,
  ) catch {
    _ => abort("valid RNA QC metrics should produce thresholds")
  }
  let keep = @src.scrapper_filter_rna_qc_metrics(thresholds, metrics) catch {
    _ => abort("thresholds should apply to their source metrics")
  }
  assert_eq(thresholds.block_names, ["all"])
  scrapper_test_assert_close(thresholds.sum_lower[0], 10.0, 1.0e-10)
  assert_eq(keep, [false, true, true])
}

///|
test "scrapper: blocked MAD thresholds are estimated independently" {
  let metrics = @src.scrapper_compute_rna_qc_metrics([[1.0, 2.0, 100.0, 200.0]]) catch {
    _ => abort("valid RNA QC input should be accepted")
  }
  let blocks = ["A", "A", "B", "B"]
  let thresholds = @src.scrapper_suggest_rna_qc_thresholds(
    metrics,
    blocks~,
    sum_num_mads=0.0,
    detected_num_mads=0.0,
  ) catch {
    _ => abort("valid blocks should produce thresholds")
  }
  let keep = @src.scrapper_filter_rna_qc_metrics(thresholds, metrics, blocks~) catch {
    _ => abort("matching blocks should be accepted")
  }
  assert_eq(thresholds.block_names, ["A", "B"])
  assert_true(thresholds.sum_lower[1] > thresholds.sum_lower[0])
  assert_eq(keep, [false, true, false, true])
}

///|
test "scrapper: subset MAD upper threshold removes high-proportion cells" {
  let metrics = @src.scrapper_compute_rna_qc_metrics(
    [[9.0, 1.0, 5.0], [1.0, 9.0, 5.0]],
    subsets=[@src.ScrapperNamedSubset::new("mito", [0])],
  ) catch {
    _ => abort("valid RNA QC input should be accepted")
  }
  let thresholds = @src.scrapper_suggest_rna_qc_thresholds(
    metrics,
    sum_num_mads=100.0,
    detected_num_mads=100.0,
    subset_proportion_num_mads=0.0,
  ) catch {
    _ => abort("valid RNA QC metrics should produce thresholds")
  }
  let keep = @src.scrapper_filter_rna_qc_metrics(thresholds, metrics) catch {
    _ => abort("thresholds should apply to their source metrics")
  }
  assert_eq(thresholds.subset_upper[0], [0.5])
  assert_eq(keep, [false, true, true])
}

///|
test "scrapper: filter rejects a different block level ordering" {
  let metrics = @src.scrapper_compute_rna_qc_metrics([[1.0, 2.0]]) catch {
    _ => abort("valid RNA QC input should be accepted")
  }
  let thresholds = @src.scrapper_suggest_rna_qc_thresholds(metrics, blocks=[
    "A", "B",
  ]) catch {
    _ => abort("valid blocks should produce thresholds")
  }
  let raised = try {
    ignore(
      @src.scrapper_filter_rna_qc_metrics(thresholds, metrics, blocks=["B", "A"]),
    )
    false
  } catch {
    ScrapperError(_) => true
  }
  assert_true(raised)
}

///|
test "scrapper: quick RNA QC returns metrics thresholds and summary" {
  let result = @src.scrapper_quick_rna_qc(
    [[1.0, 10.0, 100.0]],
    sum_num_mads=0.0,
    detected_num_mads=0.0,
  ) catch {
    _ => abort("quick RNA QC should accept valid counts")
  }
  assert_eq(result.keep, [false, true, true])
  assert_eq(result.retained_count(), 2)
  assert_eq(result.summary(), "ScrapperRNAQC(cells=3, retained=2, subsets=0)")
}

///|
test "scrapper: sanitizes invalid size factors without mutating input" {
  let input = [2.0, 0.0, -1.0, @double.not_a_number]
  let output = @src.scrapper_sanitize_size_factors(input, fallback=3.0) catch {
    _ => abort("positive fallback should be accepted")
  }
  assert_eq(output, [2.0, 3.0, 3.0, 3.0])
  assert_eq(input[0], 2.0)
  assert_eq(input[1], 0.0)
  assert_true(input[3].is_nan())
}

///|
test "scrapper: size-factor sanitization validates fallback" {
  let raised = try {
    ignore(@src.scrapper_sanitize_size_factors([1.0], fallback=0.0))
    false
  } catch {
    ScrapperError(_) => true
  }
  assert_true(raised)
}

///|
test "scrapper: centers size factors across all cells" {
  let centered = @src.scrapper_center_size_factors([1.0, 2.0, 3.0]) catch {
    _ => abort("positive size factors should be centered")
  }
  assert_eq(centered, [0.5, 1.0, 1.5])
}

///|
test "scrapper: per-block centering gives each block unit mean" {
  let centered = @src.scrapper_center_size_factors(
    [1.0, 3.0, 10.0, 20.0],
    blocks=["A", "A", "B", "B"],
    mode=@src.scrapper_center_per_block(),
  ) catch {
    _ => abort("valid blocked size factors should be centered")
  }
  assert_eq(centered[0], 0.5)
  assert_eq(centered[1], 1.5)
  scrapper_test_assert_close(centered[2], 2.0 / 3.0, 1.0e-12)
  scrapper_test_assert_close(centered[3], 4.0 / 3.0, 1.0e-12)
}

///|
test "scrapper: lowest-block centering preserves between-block scale" {
  let centered = @src.scrapper_center_size_factors(
    [1.0, 3.0, 10.0, 20.0],
    blocks=["A", "A", "B", "B"],
    mode=@src.scrapper_center_lowest(),
  ) catch {
    _ => abort("valid blocked size factors should be centered")
  }
  assert_eq(centered, [0.5, 1.5, 5.0, 10.0])
}

///|
test "scrapper: computes centered library size factors" {
  let factors = @src.scrapper_library_size_factors([[1.0, 2.0], [3.0, 6.0]]) catch {
    _ => abort("valid counts should produce library size factors")
  }
  scrapper_test_assert_close(factors[0], 2.0 / 3.0, 1.0e-12)
  scrapper_test_assert_close(factors[1], 4.0 / 3.0, 1.0e-12)
}

///|
test "scrapper: zero libraries receive finite fallback factors" {
  let factors = @src.scrapper_library_size_factors([[0.0, 0.0]]) catch {
    _ => abort("zero libraries should use sanitized factors")
  }
  assert_eq(factors, [1.0, 1.0])
}

///|
test "scrapper: performs linear count scaling" {
  let normalized = @src.scrapper_normalize_counts(
    [[2.0, 8.0], [4.0, 4.0]],
    [2.0, 4.0],
    log_transform=false,
  ) catch {
    _ => abort("valid size factors should normalize counts")
  }
  assert_eq(normalized, [[1.0, 2.0], [2.0, 1.0]])
}

///|
test "scrapper: performs configurable log normalization" {
  let normalized = @src.scrapper_normalize_counts(
    [[3.0, 8.0]],
    [1.0, 2.0],
    pseudo_count=1.0,
    log_base=2.0,
  ) catch {
    _ => abort("valid size factors should normalize counts")
  }
  scrapper_test_assert_close(normalized[0][0], 2.0, 1.0e-12)
  scrapper_test_assert_close(normalized[0][1], 2.321928094887362, 1.0e-12)
}

///|
test "scrapper: normalization validates size factors and logarithm base" {
  let wrong_length = try {
    ignore(@src.scrapper_normalize_counts([[1.0, 2.0]], [1.0]))
    false
  } catch {
    ScrapperError(_) => true
  }
  let non_positive = try {
    ignore(@src.scrapper_normalize_counts([[1.0]], [0.0]))
    false
  } catch {
    ScrapperError(_) => true
  }
  let invalid_base = try {
    ignore(@src.scrapper_normalize_counts([[1.0]], [1.0], log_base=1.0))
    false
  } catch {
    ScrapperError(_) => true
  }
  assert_true(wrong_length)
  assert_true(non_positive)
  assert_true(invalid_base)
}

///|
test "scrapper: models per-gene means and sample variances" {
  let model = @src.scrapper_model_gene_variances(
    [[1.0, 2.0, 3.0], [2.0, 2.0, 2.0], [0.0, 0.0, 0.0]],
    mean_filter=false,
    transform=false,
    span=1.0,
    min_window_count=1,
  ) catch {
    _ => abort("valid expression should produce a variance model")
  }
  assert_eq(model.means, [2.0, 2.0, 0.0])
  assert_eq(model.variances, [1.0, 0.0, 0.0])
  let mut feature = 0
  while feature < model.means.length() {
    scrapper_test_assert_close(
      model.residuals[feature],
      model.variances[feature] - model.fitted[feature],
      1.0e-12,
    )
    feature = feature + 1
  }
}

///|
test "scrapper: LOWESS reproduces a linear mean-variance trend" {
  let trend = @src.scrapper_fit_variance_trend(
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [2.0, 4.0, 6.0, 8.0, 10.0],
    mean_filter=false,
    transform=false,
    span=1.0,
    min_window_count=5,
  ) catch {
    _ => abort("valid trend inputs should be fitted")
  }
  let mut index = 0
  while index < trend.fitted.length() {
    scrapper_test_assert_close(
      trend.fitted[index],
      (index + 1).to_double() * 2.0,
      1.0e-9,
    )
    scrapper_test_assert_close(trend.residuals[index], 0.0, 1.0e-9)
    index = index + 1
  }
}

///|
test "scrapper: variance trend extrapolates below the mean filter" {
  let trend = @src.scrapper_fit_variance_trend(
    [0.05, 1.0, 2.0, 3.0],
    [0.01, 1.0, 2.0, 3.0],
    min_mean=1.0,
    transform=false,
    span=1.0,
    min_window_count=3,
  ) catch {
    _ => abort("valid filtered trend should be fitted")
  }
  assert_true(trend.fitted[0] >= 0.0)
  assert_true(trend.fitted[0] < trend.fitted[1])
}

///|
test "scrapper: variance trend validates fitting parameters" {
  let mismatched = try {
    ignore(@src.scrapper_fit_variance_trend([1.0], [1.0, 2.0]))
    false
  } catch {
    ScrapperError(_) => true
  }
  let no_genes = try {
    ignore(
      @src.scrapper_fit_variance_trend(
        [0.0],
        [0.0],
        mean_filter=true,
        min_mean=1.0,
      ),
    )
    false
  } catch {
    ScrapperError(_) => true
  }
  assert_true(mismatched)
  assert_true(no_genes)
}

///|
test "scrapper: selects top HVGs with stable tie handling" {
  let selected = @src.scrapper_choose_highly_variable_genes(
    [0.5, 2.0, 2.0, -1.0, 1.0],
    top=1,
    keep_ties=true,
  ) catch {
    _ => abort("finite statistics should be ranked")
  }
  assert_eq(selected, [1, 2])
}

///|
test "scrapper: HVG selection supports smaller statistics and no bound" {
  let selected = @src.scrapper_choose_highly_variable_genes(
    [0.5, 2.0, 2.0, -1.0, 1.0],
    top=2,
    larger=false,
    keep_ties=false,
    use_bound=false,
  ) catch {
    _ => abort("finite statistics should be ranked")
  }
  assert_eq(selected, [3, 0])
}

///|
test "scrapper: variance model convenience method uses positive residuals" {
  let model = @src.scrapper_model_gene_variances(
    [[1.0, 2.0, 8.0], [2.0, 2.0, 2.0], [1.0, 2.0, 3.0]],
    mean_filter=false,
    transform=false,
    span=1.0,
    min_window_count=1,
  ) catch {
    _ => abort("valid expression should produce a variance model")
  }
  let selected = model.highly_variable_genes(top=3) catch {
    _ => abort("finite residuals should be ranked")
  }
  let expected = @src.scrapper_choose_highly_variable_genes(
    model.residuals,
    top=3,
    use_bound=true,
    bound=0.0,
  ) catch {
    _ => abort("finite residuals should be ranked")
  }
  assert_eq(selected, expected)
}

///|
test "scrapper: aggregates sums detection medians and means by one factor" {
  let result = @src.scrapper_aggregate_across_cells(
    [[1.0, 2.0, 3.0, 4.0], [0.0, 5.0, 0.0, 7.0]],
    [@src.ScrapperFactor::new("cluster", ["A", "A", "B", "B"])],
    compute_median=true,
  ) catch {
    _ => abort("valid grouping factor should aggregate expression")
  }
  let means = result.means() catch {
    _ => abort("aggregate sums should support means")
  }
  assert_eq(result.factor_names, ["cluster"])
  assert_eq(result.group_names, ["cluster=A", "cluster=B"])
  assert_eq(result.counts, [2, 2])
  assert_eq(result.index, [0, 0, 1, 1])
  assert_eq(result.sums, [[3.0, 7.0], [5.0, 7.0]])
  assert_eq(result.detected, [[2, 2], [1, 1]])
  assert_eq(result.medians, [[1.5, 3.5], [2.5, 3.5]])
  assert_eq(means, [[1.5, 3.5], [2.5, 3.5]])
  assert_eq(result.n_groups(), 2)
  assert_eq(
    result.summary(),
    "ScrapperAggregate(groups=2, cells=4, features=2)",
  )
}

///|
test "scrapper: aggregates unique combinations of multiple factors" {
  let result = @src.scrapper_aggregate_across_cells([[1.0, 2.0, 3.0, 4.0]], [
    @src.ScrapperFactor::new("cluster", ["A", "A", "B", "B"]),
    @src.ScrapperFactor::new("batch", ["X", "Y", "X", "Y"]),
  ]) catch {
    _ => abort("valid grouping factors should aggregate expression")
  }
  assert_eq(result.factor_names, ["cluster", "batch"])
  assert_eq(result.n_groups(), 4)
  assert_eq(result.group_names, [
    "cluster=A,batch=X", "cluster=A,batch=Y", "cluster=B,batch=X", "cluster=B,batch=Y",
  ])
  assert_eq(result.index, [0, 1, 2, 3])
}

///|
test "scrapper: aggregation can omit optional summaries" {
  let result = @src.scrapper_aggregate_across_cells(
    [[1.0, 2.0]],
    [@src.ScrapperFactor::new("group", ["A", "A"])],
    compute_sum=false,
    compute_detected=false,
    compute_median=true,
  ) catch {
    _ => abort("median-only aggregation should be accepted")
  }
  assert_eq(result.sums.length(), 0)
  assert_eq(result.detected.length(), 0)
  assert_eq(result.medians, [[1.5]])
  let means_missing = try {
    ignore(result.means())
    false
  } catch {
    ScrapperError(_) => true
  }
  assert_true(means_missing)
}

///|
test "scrapper: aggregation validates grouping factors" {
  let no_factors = try {
    ignore(@src.scrapper_aggregate_across_cells([[1.0]], []))
    false
  } catch {
    ScrapperError(_) => true
  }
  let duplicate_names = try {
    ignore(
      @src.scrapper_aggregate_across_cells([[1.0]], [
        @src.ScrapperFactor::new("group", ["A"]),
        @src.ScrapperFactor::new("group", ["B"]),
      ]),
    )
    false
  } catch {
    ScrapperError(_) => true
  }
  let wrong_length = try {
    ignore(
      @src.scrapper_aggregate_across_cells([[1.0, 2.0]], [
        @src.ScrapperFactor::new("group", ["A"]),
      ]),
    )
    false
  } catch {
    ScrapperError(_) => true
  }
  assert_true(no_factors)
  assert_true(duplicate_names)
  assert_true(wrong_length)
}

///|
test "scrapper: aggregation validates a finite detection limit" {
  let raised = try {
    ignore(
      @src.scrapper_aggregate_across_cells(
        [[1.0]],
        [@src.ScrapperFactor::new("group", ["A"])],
        detection_limit=@double.not_a_number,
      ),
    )
    false
  } catch {
    ScrapperError(_) => true
  }
  assert_true(raised)
}

///|
test "scrapper: SCE normalization adds assay and size factors immutably" {
  let original = scrapper_test_sce()
  let normalized = @src.scrapper_normalize_rna_counts_sce(original, size_factors=[
    1.0, 1.0, 1.0,
  ]) catch {
    _ => abort("valid SCE counts should be normalized")
  }
  assert_eq(@src.sce_get_assay(original, "logcounts").length(), 0)
  let logcounts = @src.sce_get_assay(normalized, "logcounts")
  assert_eq(logcounts.length(), 3)
  scrapper_test_assert_close(logcounts[0][0], 2.0, 1.0e-12)
  assert_eq(@src.sce_get_col_data(normalized, "size_factor"), ["1", "1", "1"])
  assert_eq(@src.sce_get_col_data(original, "size_factor").length(), 0)
  assert_eq(@src.sce_get_row_data(normalized, "symbol"), ["A", "B", "MT-C"])
  assert_eq(@src.sce_get_reduced_dim(normalized, "PCA").length(), 3)
}

///|
test "scrapper: SCE normalization deep-copies assay matrices" {
  let original = scrapper_test_sce()
  let normalized = @src.scrapper_normalize_rna_counts_sce(original, size_factors=[
    1.0, 1.0, 1.0,
  ]) catch {
    _ => abort("valid SCE counts should be normalized")
  }
  let copied_counts = @src.sce_get_assay(normalized, "counts")
  copied_counts[0][0] = 999.0
  assert_eq(@src.sce_get_assay(original, "counts")[0][0], 3.0)
}

///|
test "scrapper: SCE quick QC adds copied cell metadata" {
  let original = scrapper_test_sce()
  let (annotated, result) = @src.scrapper_quick_rna_qc_sce(
    original,
    subsets=[@src.ScrapperNamedSubset::new("mito", [2])],
    blocks=["x", "x", "y"],
  ) catch {
    _ => abort("valid SCE counts should support quick RNA QC")
  }
  assert_eq(result.metrics.sums, [4.0, 9.0, 12.0])
  assert_eq(@src.sce_get_col_data(annotated, "scrapper_sum").length(), 3)
  assert_eq(@src.sce_get_col_data(annotated, "scrapper_detected").length(), 3)
  assert_eq(@src.sce_get_col_data(annotated, "scrapper_keep").length(), 3)
  assert_eq(
    @src.sce_get_col_data(annotated, "scrapper_subset_mito_proportion").length(),
    3,
  )
  assert_eq(@src.sce_get_col_data(original, "scrapper_sum").length(), 0)
  assert_eq(@src.sce_get_col_data(annotated, "batch"), ["x", "x", "y"])
}
