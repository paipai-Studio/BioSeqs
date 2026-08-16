///|
fn tradeseq_adv_test_counts() -> Array[Array[Double]] {
  let flat : Array[Double] = []
  let association : Array[Double] = []
  let endpoint : Array[Double] = []
  let early_branch : Array[Double] = []
  let zero : Array[Double] = []
  for cell in 0..<24 {
    flat.push((10 + cell % 3).to_double())
    if cell < 8 {
      association.push((3 + cell).to_double())
      endpoint.push((6 + cell % 2).to_double())
      early_branch.push((8 + cell % 2).to_double())
    } else if cell < 16 {
      let step = cell - 8
      association.push((12 + 3 * step).to_double())
      endpoint.push((10 + 4 * step).to_double())
      early_branch.push(
        if step < 4 {
          (26 + step).to_double()
        } else {
          (14 + step % 2).to_double()
        },
      )
    } else {
      let step = cell - 16
      association.push((12 + 3 * step).to_double())
      endpoint.push((10 + step / 3).to_double())
      early_branch.push(
        if step < 4 {
          (4 + step).to_double()
        } else {
          (14 + step % 2).to_double()
        },
      )
    }
    zero.push(0.0)
  }
  [flat, association, endpoint, early_branch, zero]
}

///|
fn tradeseq_adv_test_trajectory() -> (
  Array[Array[Double]],
  Array[Array[Double]],
) {
  let pseudotime : Array[Array[Double]] = []
  let weights : Array[Array[Double]] = []
  for cell in 0..<24 {
    if cell < 8 {
      let time = cell.to_double() * 0.05
      pseudotime.push([time, time])
      weights.push([2.0, 2.0])
    } else if cell < 16 {
      let time = 0.4 + (cell - 8).to_double() * 0.08
      pseudotime.push([time, 0.0])
      weights.push([3.0, 0.0])
    } else {
      let time = 0.4 + (cell - 16).to_double() * 0.08
      pseudotime.push([0.0, time])
      weights.push([0.0, 4.0])
    }
  }
  (pseudotime, weights)
}

///|
fn tradeseq_adv_test_config() -> @src.TradeSeqAdvancedConfig {
  @src.TradeSeqAdvancedConfig::create(
    n_knots=4,
    smoothing_penalty=1.0,
    max_iterations=80,
    tolerance=1.0e-5,
    minimum_dispersion=1.0e-4,
    maximum_dispersion=20.0,
    ridge=1.0e-6,
    fdr_threshold=0.1,
    test_points=6,
  ) catch {
    _ => abort("valid tradeSeq advanced configuration should build")
  }
}

///|
fn tradeseq_adv_test_fit() -> @src.TradeSeqAdvancedFit {
  let (pseudotime, weights) = tradeseq_adv_test_trajectory()
  @src.tradeseq_fit_advanced(
    tradeseq_adv_test_counts(),
    pseudotime,
    weights,
    gene_names=["flat", "association", "endpoint", "early", "zero"],
    lineage_names=["left", "right"],
    offsets=Array::make(24, 0.0),
    config=tradeseq_adv_test_config(),
  ) catch {
    _ => abort("valid tradeSeq advanced data should fit")
  }
}

///|
fn tradeseq_adv_test_single_lineage_fit() -> @src.TradeSeqAdvancedFit {
  let pseudotime : Array[Array[Double]] = []
  let weights : Array[Array[Double]] = []
  let counts : Array[Double] = []
  for cell in 0..<12 {
    pseudotime.push([cell.to_double()])
    weights.push([1.0])
    counts.push((cell + 1).to_double())
  }
  @src.tradeseq_fit_advanced(
    [counts],
    pseudotime,
    weights,
    gene_names=["gene"],
    offsets=Array::make(12, 0.0),
    config=@src.TradeSeqAdvancedConfig::create(
      n_knots=3,
      max_iterations=40,
      ridge=1.0e-5,
      test_points=4,
    ),
  ) catch {
    _ => abort("valid single-lineage tradeSeq data should fit")
  }
}

///|
fn tradeseq_adv_test_sce() -> @src.SingleCellExperiment {
  let names : Array[String] = []
  for cell in 0..<24 {
    names.push("cell" + (cell + 1).to_string())
  }
  let experiment = @src.SingleCellExperiment::new(
    tradeseq_adv_test_counts(),
    ["flat", "association", "endpoint", "early", "zero"],
    names,
  )
  let (pseudotime, weights) = tradeseq_adv_test_trajectory()
  experiment.reduced_dims["slingshot.pseudotime"] = pseudotime
  experiment.reduced_dims["slingshot.weights"] = weights
  experiment.row_data["symbol"] = ["F", "A", "E", "B", "Z"]
  experiment.col_data["batch"] = Array::make(24, "one")
  experiment.metadata["source"] = "test"
  experiment.alternative_experiments["spike"] = @src.SingleCellExperiment::new(
    [[1.0, 2.0]],
    ["spike1"],
    ["s1", "s2"],
  )
  experiment
}

///|
fn tradeseq_adv_test_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.0e299
}

///|
fn tradeseq_adv_test_slingshot_result() -> @src.SlingshotAdvancedResult {
  let coordinates = [
    [-0.1, 0.0],
    [0.0, 0.1],
    [0.1, -0.1],
    [0.9, 0.0],
    [1.0, 0.1],
    [1.1, -0.1],
    [1.9, 0.0],
    [2.0, 0.1],
    [2.1, -0.1],
    [2.9, 0.9],
    [3.0, 1.0],
    [3.1, 1.1],
    [2.9, -0.9],
    [3.0, -1.0],
    [3.1, -1.1],
  ]
  let labels = [
    "A", "A", "A", "B", "B", "B", "C", "C", "C", "D", "D", "D", "E", "E", "E",
  ]
  let config = @src.SlingshotAdvancedConfig::create(
    start_clusters=["A"],
    end_clusters=["D", "E"],
    distance=@src.slingshot_center_euclidean(),
    extension=@src.slingshot_extend_none(),
    max_iterations=6,
    tolerance=1.0e-4,
    smoother_span=0.4,
    curve_points=24,
  ) catch {
    _ => abort("valid Slingshot configuration should build")
  }
  @src.slingshot_advanced(coordinates, labels, config) catch {
    _ => abort("valid Slingshot data should fit")
  }
}

///|
test "tradeSeq advanced default configuration matches fitGAM workflow" {
  let config = @src.TradeSeqAdvancedConfig::default()
  assert_eq(config.n_knots, 6)
  assert_eq(config.smoothing_penalty, 1.0)
  assert_eq(config.max_iterations, 60)
  assert_eq(config.tolerance, 1.0e-6)
  assert_eq(config.minimum_dispersion, 1.0e-4)
  assert_eq(config.maximum_dispersion, 100.0)
  assert_eq(config.ridge, 1.0e-8)
  assert_eq(config.fdr_threshold, 0.05)
  assert_eq(config.test_points, 12)
}

///|
test "tradeSeq advanced configuration preserves explicit controls" {
  let config = tradeseq_adv_test_config()
  assert_eq(config.n_knots, 4)
  assert_eq(config.smoothing_penalty, 1.0)
  assert_eq(config.max_iterations, 80)
  assert_eq(config.tolerance, 1.0e-5)
  assert_eq(config.maximum_dispersion, 20.0)
  assert_eq(config.ridge, 1.0e-6)
  assert_eq(config.fdr_threshold, 0.1)
  assert_eq(config.test_points, 6)
}

///|
test "tradeSeq advanced rejects invalid knots and smoothing penalty" {
  let mut failures = 0
  ignore(@src.TradeSeqAdvancedConfig::create(n_knots=2)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.TradeSeqAdvancedConfig::create(smoothing_penalty=-1.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.TradeSeqAdvancedConfig::create(smoothing_penalty=0.0 / 0.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "tradeSeq advanced rejects invalid iteration controls" {
  let mut failures = 0
  ignore(@src.TradeSeqAdvancedConfig::create(max_iterations=0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.TradeSeqAdvancedConfig::create(tolerance=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.TradeSeqAdvancedConfig::create(tolerance=0.0 / 0.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "tradeSeq advanced rejects invalid dispersion and ridge controls" {
  let mut failures = 0
  ignore(@src.TradeSeqAdvancedConfig::create(minimum_dispersion=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.TradeSeqAdvancedConfig::create(
      minimum_dispersion=2.0,
      maximum_dispersion=1.0,
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(@src.TradeSeqAdvancedConfig::create(ridge=0.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "tradeSeq advanced rejects invalid testing controls" {
  let mut failures = 0
  ignore(@src.TradeSeqAdvancedConfig::create(fdr_threshold=0.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.TradeSeqAdvancedConfig::create(fdr_threshold=1.1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.TradeSeqAdvancedConfig::create(test_points=1)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "tradeSeq advanced validates count matrix shape and values" {
  let (pseudotime, weights) = tradeseq_adv_test_trajectory()
  let mut failures = 0
  ignore(@src.tradeseq_fit_advanced([], pseudotime, weights)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_fit_advanced([[]], pseudotime, weights)) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced([[1.0, 2.0], [1.0]], [[0.0], [1.0]], [
      [1.0],
      [1.0],
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced([[-1.0, 2.0]], [[0.0], [1.0]], [[1.0], [1.0]]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced([[1.5, 2.0]], [[0.0], [1.0]], [[1.0], [1.0]]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced([[0.0 / 0.0, 2.0]], [[0.0], [1.0]], [
      [1.0],
      [1.0],
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 6)
}

///|
test "tradeSeq advanced validates trajectory row dimensions" {
  let mut failures = 0
  ignore(@src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0]], [[1.0]])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_fit_advanced([[1.0, 2.0]], [[], []], [[], []])) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0], [1.0, 2.0]], [[1.0], [1.0]]),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "tradeSeq advanced validates lineage weight semantics" {
  let mut failures = 0
  ignore(
    @src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0], [1.0]], [[-1.0], [1.0]]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0], [1.0]], [[0.0], [1.0]]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0], [1.0]], [
      [0.0 / 0.0],
      [1.0],
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "tradeSeq advanced requires finite pseudotime for active lineages" {
  let failed = try {
    ignore(
      @src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0 / 0.0], [1.0]], [
        [1.0],
        [1.0],
      ]),
    )
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "tradeSeq advanced permits nonfinite inactive pseudotime" {
  let fit = @src.tradeseq_fit_advanced(
    [[1.0, 2.0, 3.0, 4.0]],
    [[0.0, 0.0 / 0.0], [1.0, 0.0 / 0.0], [0.0 / 0.0, 3.0], [0.0 / 0.0, 4.0]],
    [[1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]],
    offsets=[0.0, 0.0, 0.0, 0.0],
    config=@src.TradeSeqAdvancedConfig::create(n_knots=3),
  ) catch {
    _ => abort("inactive nonfinite pseudotime should be ignored")
  }
  assert_eq(fit.n_lineages, 2)
}

///|
test "tradeSeq advanced requires pseudotime range in every lineage" {
  let failed = try {
    ignore(
      @src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0], [0.0]], [[1.0], [1.0]]),
    )
    false
  } catch {
    _ => true
  }
  assert_true(failed)
}

///|
test "tradeSeq advanced validates gene and lineage name counts" {
  let (pseudotime, weights) = tradeseq_adv_test_trajectory()
  let mut failures = 0
  ignore(
    @src.tradeseq_fit_advanced(
      tradeseq_adv_test_counts(),
      pseudotime,
      weights,
      gene_names=["one"],
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced(
      tradeseq_adv_test_counts(),
      pseudotime,
      weights,
      lineage_names=["one"],
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "tradeSeq advanced validates nonempty unique names" {
  let mut failures = 0
  ignore(
    @src.tradeseq_fit_advanced(
      [[1.0, 2.0], [2.0, 3.0]],
      [[0.0], [1.0]],
      [[1.0], [1.0]],
      gene_names=["gene", "gene"],
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0], [1.0]], [[1.0], [1.0]], lineage_names=[
      " ",
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "tradeSeq advanced validates supplied offsets" {
  let mut failures = 0
  ignore(
    @src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0], [1.0]], [[1.0], [1.0]], offsets=[
      0.0,
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_fit_advanced([[1.0, 2.0]], [[0.0], [1.0]], [[1.0], [1.0]], offsets=[
      0.0,
      0.0 / 0.0,
    ]),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "tradeSeq advanced fit exposes feature by cell dimensions" {
  let fit = tradeseq_adv_test_fit()
  assert_eq(fit.n_genes, 5)
  assert_eq(fit.n_cells, 24)
  assert_eq(fit.n_lineages, 2)
  assert_eq(fit.gene_names, ["flat", "association", "endpoint", "early", "zero"])
  assert_eq(fit.lineage_names, ["left", "right"])
  assert_eq(fit.models.length(), 5)
  assert_eq(fit.models[0].fitted_values.length(), 24)
}

///|
test "tradeSeq advanced normalizes each cell lineage weight row" {
  let fit = tradeseq_adv_test_fit()
  assert_eq(fit.cell_weights[0], [0.5, 0.5])
  assert_eq(fit.cell_weights[8], [1.0, 0.0])
  assert_eq(fit.cell_weights[16], [0.0, 1.0])
  for row in fit.cell_weights {
    let mut total = 0.0
    for value in row {
      total = total + value
    }
    assert_true((total - 1.0).abs() < 1.0e-12)
  }
}

///|
test "tradeSeq advanced derives centered library size offsets" {
  let (pseudotime, weights) = tradeseq_adv_test_trajectory()
  let fit = @src.tradeseq_fit_advanced(
    tradeseq_adv_test_counts(),
    pseudotime,
    weights,
    config=tradeseq_adv_test_config(),
  ) catch {
    _ => abort("valid automatic offsets should fit")
  }
  let mut total = 0.0
  for value in fit.offsets {
    total = total + value
  }
  assert_true(total.abs() < 1.0e-10)
}

///|
test "tradeSeq advanced preserves caller matrices" {
  let (pseudotime, weights) = tradeseq_adv_test_trajectory()
  let fit = @src.tradeseq_fit_advanced(
    tradeseq_adv_test_counts(),
    pseudotime,
    weights,
    config=tradeseq_adv_test_config(),
  ) catch {
    _ => abort("valid caller matrices should fit")
  }
  fit.pseudotime[0][0] = 99.0
  fit.cell_weights[0][0] = 0.0
  assert_eq(pseudotime[0][0], 0.0)
  assert_eq(weights[0][0], 2.0)
}

///|
test "tradeSeq advanced models have finite penalized NB statistics" {
  let fit = tradeseq_adv_test_fit()
  for model in fit.models {
    assert_eq(model.coefficients.length(), 8)
    assert_eq(model.covariance.length(), 8)
    assert_eq(model.covariance[0].length(), 8)
    assert_true(model.dispersion >= fit.config.minimum_dispersion)
    assert_true(model.dispersion <= fit.config.maximum_dispersion)
    assert_true(tradeseq_adv_test_finite(model.log_likelihood))
    assert_true(tradeseq_adv_test_finite(model.aic))
    assert_eq(model.effective_df, 8.0)
    assert_true(model.iterations >= 1)
    for value in model.coefficients {
      assert_true(tradeseq_adv_test_finite(value))
    }
    for value in model.fitted_values {
      assert_true(tradeseq_adv_test_finite(value))
      assert_true(value >= 0.0)
    }
  }
}

///|
test "tradeSeq advanced handles an all-zero gene" {
  let fit = tradeseq_adv_test_fit()
  let model = fit.models[4]
  assert_eq(model.gene_id, "zero")
  assert_true(model.dispersion >= fit.config.minimum_dispersion)
  for value in model.fitted_values {
    assert_true(value < 1.0e-5)
  }
}

///|
test "tradeSeq advanced fitted trend follows increasing counts" {
  let fit = tradeseq_adv_test_fit()
  let fitted = fit.models[1].fitted_values
  assert_true(fitted[15] > fitted[0])
  assert_true(fitted[23] > fitted[0])
}

///|
test "tradeSeq advanced associationTest returns gene-level Wald results" {
  let table = @src.tradeseq_association_test_advanced(tradeseq_adv_test_fit())
  assert_eq(table.test_name, "associationTest")
  assert_eq(table.results.length(), 5)
  assert_eq(table.fdr_threshold, 0.1)
  for result in table.results {
    assert_true(result.degrees_freedom > 0)
    assert_true(result.wald_statistic >= 0.0)
    assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
    assert_true(
      result.adjusted_p_value >= 0.0 && result.adjusted_p_value <= 1.0,
    )
    assert_true(result.adjusted_p_value + 1.0e-12 >= result.p_value)
  }
}

///|
test "tradeSeq advanced associationTest ranks temporal signal above flat gene" {
  let table = @src.tradeseq_association_test_advanced(tradeseq_adv_test_fit())
  assert_true(table.results[1].wald_statistic > table.results[0].wald_statistic)
  assert_true(table.results[1].p_value < table.results[0].p_value)
}

///|
test "tradeSeq advanced startVsEndTest detects endpoint change" {
  let table = @src.tradeseq_start_vs_end_test_advanced(tradeseq_adv_test_fit())
  assert_eq(table.test_name, "startVsEndTest")
  assert_eq(table.results.length(), 5)
  assert_eq(table.results[0].degrees_freedom, 2)
  assert_true(table.results[1].wald_statistic > table.results[0].wald_statistic)
}

///|
test "tradeSeq advanced diffEndTest detects lineage endpoint separation" {
  let table = @src.tradeseq_diff_end_test_advanced(tradeseq_adv_test_fit())
  assert_eq(table.test_name, "diffEndTest")
  assert_eq(table.results[0].degrees_freedom, 1)
  assert_true(table.results[2].wald_statistic > table.results[0].wald_statistic)
}

///|
test "tradeSeq advanced patternTest compares complete lineage smooths" {
  let table = @src.tradeseq_pattern_test_advanced(
    tradeseq_adv_test_fit(),
    n_points=5,
  )
  assert_eq(table.test_name, "patternTest")
  assert_eq(table.results.length(), 5)
  assert_true(table.results[2].degrees_freedom >= 2)
  assert_true(table.results[2].wald_statistic > table.results[0].wald_statistic)
}

///|
test "tradeSeq advanced earlyDETest restricts the comparison interval" {
  let table = @src.tradeseq_early_de_test_advanced(
    tradeseq_adv_test_fit(),
    0.35,
    0.75,
    n_points=5,
  )
  assert_eq(table.test_name, "earlyDETest")
  assert_eq(table.results.length(), 5)
  assert_true(table.results[3].degrees_freedom >= 2)
  assert_true(table.results[3].wald_statistic > table.results[0].wald_statistic)
}

///|
test "tradeSeq advanced fold-change threshold shrinks Wald evidence" {
  let fit = tradeseq_adv_test_fit()
  let unthresholded = @src.tradeseq_start_vs_end_test_advanced(fit)
  let thresholded = @src.tradeseq_start_vs_end_test_advanced(fit, l2fc=2.0)
  for gene in 0..<fit.n_genes {
    assert_true(
      thresholded.results[gene].wald_statistic <=
      unthresholded.results[gene].wald_statistic + 1.0e-10,
    )
    assert_true(
      thresholded.results[gene].p_value + 1.0e-10 >=
      unthresholded.results[gene].p_value,
    )
  }
}

///|
test "tradeSeq advanced tests validate fold-change threshold" {
  let fit = tradeseq_adv_test_fit()
  let mut failures = 0
  ignore(@src.tradeseq_association_test_advanced(fit, l2fc=-1.0)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_start_vs_end_test_advanced(fit, l2fc=0.0 / 0.0)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "tradeSeq advanced lineage comparisons require multiple lineages" {
  let fit = tradeseq_adv_test_single_lineage_fit()
  let mut failures = 0
  ignore(@src.tradeseq_diff_end_test_advanced(fit)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_pattern_test_advanced(fit)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_early_de_test_advanced(fit, 0.0, 0.5)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "tradeSeq advanced validates pattern and earlyDE grids" {
  let fit = tradeseq_adv_test_fit()
  let mut failures = 0
  ignore(@src.tradeseq_pattern_test_advanced(fit, n_points=1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_early_de_test_advanced(fit, -0.1, 0.5)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_early_de_test_advanced(fit, 0.5, 0.5)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_early_de_test_advanced(fit, 0.5, 1.1)) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_early_de_test_advanced(fit, 0.0, 0.5, n_points=1)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 5)
}

///|
test "tradeSeq advanced predicts smooths on each lineage" {
  let prediction = @src.tradeseq_predict_smooth_advanced(
    tradeseq_adv_test_fit(),
    "association",
    n_points=9,
  )
  assert_eq(prediction.gene_id, "association")
  assert_eq(prediction.pseudotime.length(), 2)
  assert_eq(prediction.pseudotime[0].length(), 9)
  assert_eq(prediction.fitted.length(), 2)
  assert_eq(prediction.standard_errors[1].length(), 9)
  assert_eq(prediction.pseudotime[0][0], 0.0)
  assert_true((prediction.pseudotime[0][8] - 0.96).abs() < 1.0e-12)
}

///|
test "tradeSeq advanced smooth predictions are finite and ordered" {
  let prediction = @src.tradeseq_predict_smooth_advanced(
    tradeseq_adv_test_fit(),
    "association",
    n_points=9,
  )
  for lineage in 0..<2 {
    assert_true(prediction.fitted[lineage][8] > prediction.fitted[lineage][0])
    for point in 0..<9 {
      assert_true(tradeseq_adv_test_finite(prediction.fitted[lineage][point]))
      assert_true(
        tradeseq_adv_test_finite(prediction.standard_errors[lineage][point]),
      )
      assert_true(prediction.standard_errors[lineage][point] >= 0.0)
    }
  }
}

///|
test "tradeSeq advanced smooth prediction validates gene and grid" {
  let fit = tradeseq_adv_test_fit()
  let mut failures = 0
  ignore(@src.tradeseq_predict_smooth_advanced(fit, "missing")) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_predict_smooth_advanced(fit, "flat", n_points=1)) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 2)
}

///|
test "tradeSeq advanced evaluateK returns gene by candidate AIC" {
  let (pseudotime, weights) = tradeseq_adv_test_trajectory()
  let counts = tradeseq_adv_test_counts()
  let evaluation = @src.tradeseq_evaluate_k_advanced(
    [counts[0], counts[1]],
    pseudotime,
    weights,
    [3, 4],
    gene_names=["flat", "association"],
    offsets=Array::make(24, 0.0),
    config=tradeseq_adv_test_config(),
  ) catch {
    _ => abort("valid knot candidates should be evaluated")
  }
  assert_eq(evaluation.candidates, [3, 4])
  assert_eq(evaluation.gene_aic.length(), 2)
  assert_eq(evaluation.gene_aic[0].length(), 2)
  assert_eq(evaluation.mean_aic.length(), 2)
  assert_true(
    evaluation.selected_n_knots == 3 || evaluation.selected_n_knots == 4,
  )
  for value in evaluation.mean_aic {
    assert_true(tradeseq_adv_test_finite(value))
  }
}

///|
test "tradeSeq advanced evaluateK validates candidate set" {
  let (pseudotime, weights) = tradeseq_adv_test_trajectory()
  let counts = tradeseq_adv_test_counts()
  let mut failures = 0
  ignore(@src.tradeseq_evaluate_k_advanced(counts, pseudotime, weights, [3])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_evaluate_k_advanced(counts, pseudotime, weights, [2, 3])) catch {
    _ => failures = failures + 1
  }
  ignore(@src.tradeseq_evaluate_k_advanced(counts, pseudotime, weights, [3, 3])) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 3)
}

///|
test "tradeSeq advanced fits directly from Slingshot optional pseudotime" {
  let slingshot = tradeseq_adv_test_slingshot_result()
  let first : Array[Double] = []
  let second : Array[Double] = []
  for cell in 0..<15 {
    first.push((cell + 2).to_double())
    second.push((5 + cell % 4).to_double())
  }
  let fit = @src.tradeseq_fit_from_slingshot(
    [first, second],
    ["trend", "stable"],
    slingshot,
    offsets=Array::make(15, 0.0),
    config=@src.TradeSeqAdvancedConfig::create(
      n_knots=3,
      max_iterations=50,
      ridge=1.0e-5,
      test_points=4,
    ),
  ) catch {
    _ => abort("Slingshot output should feed tradeSeq directly")
  }
  assert_eq(fit.n_cells, 15)
  assert_eq(fit.n_lineages, slingshot.curves.length())
  assert_eq(fit.lineage_names[0], slingshot.curves[0].name)
  for row in fit.cell_weights {
    let mut total = 0.0
    for value in row {
      total = total + value
    }
    assert_true((total - 1.0).abs() < 1.0e-10)
  }
}

///|
test "tradeSeq advanced SCE writes fitted assays and row statistics" {
  let output = @src.tradeseq_advanced_sce(
    tradeseq_adv_test_sce(),
    config=tradeseq_adv_test_config(),
  ) catch {
    _ => abort("valid SCE trajectory should fit")
  }
  assert_eq(output.fit.n_genes, 5)
  assert_eq(output.experiment.assays["tradeSeq.fitted"].length(), 5)
  assert_eq(output.experiment.assays["tradeSeq.fitted"][0].length(), 24)
  assert_eq(
    output.experiment.row_data["tradeSeq.association.pvalue"].length(),
    5,
  )
  assert_eq(output.experiment.row_data["tradeSeq.startVsEnd.padj"].length(), 5)
  assert_eq(output.experiment.row_data["tradeSeq.dispersion"].length(), 5)
  assert_eq(output.experiment.metadata["tradeSeq.nKnots"], "4")
  assert_eq(
    output.experiment.metadata["tradeSeq.lineages"],
    "Lineage1,Lineage2",
  )
}

///|
test "tradeSeq advanced SCE preserves input immutability" {
  let experiment = tradeseq_adv_test_sce()
  let output = @src.tradeseq_advanced_sce(
    experiment,
    config=tradeseq_adv_test_config(),
  ) catch {
    _ => abort("valid SCE trajectory should fit")
  }
  output.experiment.assays["counts"][0][0] = 99.0
  output.experiment.reduced_dims["slingshot.pseudotime"][0][0] = 99.0
  output.experiment.row_data["symbol"][0] = "changed"
  output.experiment.col_data["batch"][0] = "changed"
  assert_eq(experiment.assays["counts"][0][0], 10.0)
  assert_eq(experiment.reduced_dims["slingshot.pseudotime"][0][0], 0.0)
  assert_eq(experiment.row_data["symbol"][0], "F")
  assert_eq(experiment.col_data["batch"][0], "one")
  assert_false(experiment.assays.contains("tradeSeq.fitted"))
}

///|
test "tradeSeq advanced SCE recursively copies alternative experiments" {
  let experiment = tradeseq_adv_test_sce()
  let output = @src.tradeseq_advanced_sce(
    experiment,
    config=tradeseq_adv_test_config(),
  ) catch {
    _ => abort("valid SCE trajectory should fit")
  }
  output.experiment.alternative_experiments["spike"].assays["counts"][0][0] = 8.0
  assert_eq(
    experiment.alternative_experiments["spike"].assays["counts"][0][0],
    1.0,
  )
}

///|
test "tradeSeq advanced SCE supports custom input and output names" {
  let experiment = tradeseq_adv_test_sce()
  experiment.assays["raw"] = experiment.assays["counts"].map(fn(row) {
    row.copy()
  })
  experiment.reduced_dims["ptime"] = experiment.reduced_dims["slingshot.pseudotime"].map(fn(
      row,
    ) {
      row.copy()
    },
  )
  experiment.reduced_dims["lineageWeight"] = experiment.reduced_dims["slingshot.weights"].map(fn(
      row,
    ) {
      row.copy()
    },
  )
  let output = @src.tradeseq_advanced_sce(
    experiment,
    assay_name="raw",
    pseudotime_name="ptime",
    weights_name="lineageWeight",
    output_prefix="ts",
    config=tradeseq_adv_test_config(),
  ) catch {
    _ => abort("custom SCE names should be supported")
  }
  assert_true(output.experiment.assays.contains("ts.fitted"))
  assert_true(output.experiment.row_data.contains("ts.association.padj"))
  assert_eq(output.experiment.metadata["ts.nKnots"], "4")
}

///|
test "tradeSeq advanced SCE validates required entries and output prefix" {
  let experiment = tradeseq_adv_test_sce()
  let mut failures = 0
  ignore(
    @src.tradeseq_advanced_sce(
      experiment,
      assay_name="missing",
      config=tradeseq_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_advanced_sce(
      experiment,
      pseudotime_name="missing",
      config=tradeseq_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_advanced_sce(
      experiment,
      weights_name="missing",
      config=tradeseq_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  ignore(
    @src.tradeseq_advanced_sce(
      experiment,
      output_prefix=" ",
      config=tradeseq_adv_test_config(),
    ),
  ) catch {
    _ => failures = failures + 1
  }
  assert_eq(failures, 4)
}

///|
test "tradeSeq advanced summary reports fit dimensions" {
  let summary = @src.tradeseq_advanced_summary(tradeseq_adv_test_fit())
  assert_true(summary.contains("tradeSeq 1.27.0 advanced fit"))
  assert_true(summary.contains("genes: 5"))
  assert_true(summary.contains("cells: 24"))
  assert_true(summary.contains("lineages: 2"))
  assert_true(summary.contains("knots: 4"))
  assert_true(summary.contains("mean dispersion:"))
}
