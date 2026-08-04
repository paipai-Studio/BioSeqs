// Tests for the Bioconductor zinbwave-inspired ZINB factor model.

///|
fn zinbwave_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "zinbwave value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn zinbwave_test_config() -> @src.ZinbWaveConfig {
  @src.ZinbWaveConfig::create(
    n_factors=2,
    max_iterations=20,
    inner_iterations=1,
    tolerance=1.0e-4,
    ridge=0.01,
  ) catch {
    _ => abort("test configuration should be valid")
  }
}

///|
fn zinbwave_test_model() -> @src.ZinbWaveModel {
  let (counts, genes, cells, batch) = @src.zinbwave_example_data()
  let config = zinbwave_test_config()
  @src.zinbwave_fit(
    counts,
    config~,
    cell_covariates=batch,
    gene_names=genes,
    cell_names=cells,
  ) catch {
    _ => abort("example zinbwave fit should succeed")
  }
}

///|
fn zinbwave_test_matrix(
  rows : Int,
  columns : Int,
  value : Double,
) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = []
  for _ in 0..<rows {
    matrix.push(Array::make(columns, value))
  }
  matrix
}

///|
fn zinbwave_test_sce() -> @src.SingleCellExperiment {
  let (counts, genes, cells, _) = @src.zinbwave_example_data()
  @src.SingleCellExperiment::new(counts, genes, cells)
}

///|
test "zinbwave: default configuration exposes model controls" {
  let config = @src.ZinbWaveConfig::default()
  assert_eq(config.n_factors, 2)
  assert_eq(config.max_iterations, 50)
  assert_eq(config.inner_iterations, 2)
  assert_eq(config.tolerance, 1.0e-5)
  assert_eq(config.dispersion_shrinkage, 5.0)
}

///|
test "zinbwave: custom configuration preserves parameters" {
  let config = @src.ZinbWaveConfig::create(
    n_factors=1,
    max_iterations=12,
    inner_iterations=3,
    tolerance=1.0e-6,
    ridge=0.2,
    dispersion_shrinkage=8.0,
    minimum_mean=1.0e-7,
    minimum_probability=1.0e-5,
  ) catch {
    _ => abort("custom configuration should be valid")
  }
  assert_eq(config.n_factors, 1)
  assert_eq(config.max_iterations, 12)
  assert_eq(config.inner_iterations, 3)
  assert_eq(config.ridge, 0.2)
  assert_eq(config.minimum_probability, 1.0e-5)
}

///|
test "zinbwave: configuration rejects invalid factor and iteration counts" {
  let factors = try {
    ignore(@src.ZinbWaveConfig::create(n_factors=-1))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let iterations = try {
    ignore(@src.ZinbWaveConfig::create(max_iterations=0))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let inner = try {
    ignore(@src.ZinbWaveConfig::create(inner_iterations=0))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(factors)
  assert_true(iterations)
  assert_true(inner)
}

///|
test "zinbwave: configuration rejects invalid numeric controls" {
  let tolerance = try {
    ignore(@src.ZinbWaveConfig::create(tolerance=0.0))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let ridge = try {
    ignore(@src.ZinbWaveConfig::create(ridge=-1.0))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let shrinkage = try {
    ignore(@src.ZinbWaveConfig::create(dispersion_shrinkage=-1.0))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(tolerance)
  assert_true(ridge)
  assert_true(shrinkage)
}

///|
test "zinbwave: configuration rejects invalid probability bounds" {
  let mean = try {
    ignore(@src.ZinbWaveConfig::create(minimum_mean=0.0))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let probability = try {
    ignore(@src.ZinbWaveConfig::create(minimum_probability=0.5))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(mean)
  assert_true(probability)
}

///|
test "zinbwave: example data uses gene by cell orientation" {
  let (counts, genes, cells, batch) = @src.zinbwave_example_data()
  assert_eq(counts.length(), 6)
  assert_eq(counts[0].length(), 8)
  assert_eq(genes.length(), 6)
  assert_eq(cells.length(), 8)
  assert_eq(batch.length(), 8)
  assert_eq(batch[0].length(), 1)
}

///|
test "zinbwave: fit dimensions match the source matrix" {
  let model = zinbwave_test_model()
  assert_eq(model.n_genes(), 6)
  assert_eq(model.n_cells(), 8)
  assert_eq(model.n_factors(), 2)
  assert_eq(model.means.length(), 6)
  assert_eq(model.means[0].length(), 8)
  assert_eq(model.factors.length(), 8)
  assert_eq(model.factors[0].length(), 2)
}

///|
test "zinbwave: custom identifiers are retained" {
  let model = zinbwave_test_model()
  assert_eq(model.gene_names[0], "A_marker_1")
  assert_eq(model.gene_names[5], "dropout_gene")
  assert_eq(model.cell_names[0], "A1")
  assert_eq(model.cell_names[7], "B4")
}

///|
test "zinbwave: omitted identifiers are generated deterministically" {
  let config = @src.ZinbWaveConfig::create(n_factors=0, max_iterations=2) catch {
    _ => abort("configuration should succeed")
  }
  let model = @src.zinbwave_fit([[4.0, 3.0], [2.0, 1.0]], config~) catch {
    _ => abort("small fit should succeed")
  }
  assert_eq(model.gene_names, ["Gene1", "Gene2"])
  assert_eq(model.cell_names, ["Cell1", "Cell2"])
}

///|
test "zinbwave: cell and gene designs include intercepts" {
  let model = zinbwave_test_model()
  assert_eq(model.cell_design.length(), 8)
  assert_eq(model.cell_design[0].length(), 2)
  assert_eq(model.cell_design[0][0], 1.0)
  assert_eq(model.gene_design.length(), 6)
  assert_eq(model.gene_design[0], [1.0])
}

///|
test "zinbwave: latent factors are centered and scaled" {
  let model = zinbwave_test_model()
  for factor in 0..<model.n_factors() {
    let mut mean = 0.0
    let mut second_moment = 0.0
    for cell in 0..<model.n_cells() {
      mean = mean + model.factors[cell][factor]
      second_moment = second_moment +
        model.factors[cell][factor] * model.factors[cell][factor]
    }
    mean = mean / model.n_cells().to_double()
    second_moment = second_moment / model.n_cells().to_double()
    zinbwave_test_close(mean, 0.0, 1.0e-8)
    zinbwave_test_close(second_moment, 1.0, 1.0e-6)
  }
}

///|
test "zinbwave: deterministic initialization gives repeatable fit" {
  let first = zinbwave_test_model()
  let second = zinbwave_test_model()
  for cell in 0..<first.n_cells() {
    for factor in 0..<first.n_factors() {
      zinbwave_test_close(
        first.factors[cell][factor],
        second.factors[cell][factor],
        1.0e-10,
      )
    }
  }
}

///|
test "zinbwave: model likelihood is finite and improves" {
  let model = zinbwave_test_model()
  assert_true(model.log_likelihood.abs() < 1.0e300)
  assert_true(model.diagnostics.initial_log_likelihood.abs() < 1.0e300)
  assert_true(
    model.diagnostics.final_log_likelihood >
    model.diagnostics.initial_log_likelihood,
  )
}

///|
test "zinbwave: diagnostics track every outer iteration" {
  let model = zinbwave_test_model()
  assert_true(model.diagnostics.iterations > 0)
  assert_eq(
    model.diagnostics.log_likelihoods.length(),
    model.diagnostics.iterations + 1,
  )
  assert_eq(
    model.diagnostics.final_log_likelihood,
    model.diagnostics.log_likelihoods[model.diagnostics.log_likelihoods.length() -
    1],
  )
}

///|
test "zinbwave: fitted means are finite and positive" {
  let model = zinbwave_test_model()
  for row in model.means {
    for value in row {
      assert_true(value > 0.0)
      assert_true(value.abs() < 1.0e300)
    }
  }
}

///|
test "zinbwave: zero probabilities are strictly bounded" {
  let model = zinbwave_test_model()
  for row in model.zero_probabilities {
    for value in row {
      assert_true(value > 0.0)
      assert_true(value < 1.0)
    }
  }
  assert_true(model.mean_zero_probability() > 0.0)
  assert_true(model.mean_zero_probability() < 1.0)
}

///|
test "zinbwave: gene dispersions are finite and positive" {
  let model = zinbwave_test_model()
  assert_eq(model.dispersions.length(), model.n_genes())
  for value in model.dispersions {
    assert_true(value >= 1.0e-4)
    assert_true(value <= 100.0)
  }
}

///|
test "zinbwave: observation weights are positive and at most one" {
  let model = zinbwave_test_model()
  for row in model.observational_weights {
    for value in row {
      assert_true(value > 0.0)
      assert_true(value <= 1.0)
    }
  }
}

///|
test "zinbwave: positive counts receive unit observation weight" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let model = zinbwave_test_model()
  for gene in 0..<counts.length() {
    for cell in 0..<counts[gene].length() {
      if counts[gene][cell] > 0.0 {
        assert_eq(model.observational_weights[gene][cell], 1.0)
      }
    }
  }
}

///|
test "zinbwave: dropout-like zero receives downweighted NB responsibility" {
  let model = zinbwave_test_model()
  assert_true(model.observational_weights[5][0] < 0.5)
  assert_true(model.zero_probabilities[5][0] > 0.5)
}

///|
test "zinbwave: deviance residuals are finite with source dimensions" {
  let model = zinbwave_test_model()
  assert_eq(model.deviance_residuals.length(), 6)
  assert_eq(model.deviance_residuals[0].length(), 8)
  for row in model.deviance_residuals {
    for value in row {
      assert_true(value.abs() < 1.0e300)
    }
  }
}

///|
test "zinbwave: normalized values remove library offsets" {
  let model = zinbwave_test_model()
  let normalized = model.normalized_values()
  assert_eq(normalized.length(), model.n_genes())
  assert_eq(normalized[0].length(), model.n_cells())
  for row in normalized {
    for value in row {
      assert_true(value > 0.0)
    }
  }
}

///|
test "zinbwave: imputation preserves observed positive counts" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let model = zinbwave_test_model()
  let imputed = model.impute_zeros(counts) catch {
    _ => abort("imputation should succeed")
  }
  for gene in 0..<counts.length() {
    for cell in 0..<counts[gene].length() {
      if counts[gene][cell] > 0.0 {
        assert_eq(imputed[gene][cell], counts[gene][cell])
      }
    }
  }
}

///|
test "zinbwave: imputation fills structural zeros with non-negative signal" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let model = zinbwave_test_model()
  let imputed = model.impute_zeros(counts) catch {
    _ => abort("imputation should succeed")
  }
  assert_true(imputed[5][0] > 0.0)
  assert_true(imputed[5][0] < model.means[5][0])
}

///|
test "zinbwave: summary reports shape and latent dimension" {
  let summary = zinbwave_test_model().summary()
  assert_true(summary.contains("6 genes x 8 cells"))
  assert_true(summary.contains("K=2"))
  assert_true(summary.contains("logLik="))
}

///|
test "zinbwave: information criteria are finite" {
  let model = zinbwave_test_model()
  assert_true(model.n_parameters() > 0)
  assert_true(model.aic().abs() < 1.0e300)
  assert_true(model.bic().abs() < 1.0e300)
  assert_true(model.bic() > model.aic())
}

///|
test "zinbwave: expected mean query handles bounds" {
  let model = zinbwave_test_model()
  match model.expected_mean(0, 0) {
    Some(value) => assert_eq(value, model.means[0][0])
    None => abort("valid fitted mean should exist")
  }
  assert_true(model.expected_mean(-1, 0) is None)
  assert_true(model.expected_mean(0, 8) is None)
}

///|
test "zinbwave: zero probability query handles bounds" {
  let model = zinbwave_test_model()
  match model.zero_probability_at(5, 7) {
    Some(value) => assert_eq(value, model.zero_probabilities[5][7])
    None => abort("valid zero probability should exist")
  }
  assert_true(model.zero_probability_at(6, 0) is None)
}

///|
test "zinbwave: cell factor query returns an independent copy" {
  let model = zinbwave_test_model()
  match model.cell_factors(0) {
    Some(values) => {
      let original = model.factors[0][0]
      values[0] = values[0] + 10.0
      assert_eq(model.factors[0][0], original)
    }
    None => abort("valid cell factors should exist")
  }
  assert_true(model.cell_factors(8) is None)
}

///|
test "zinbwave: gene summary reports observed and fitted zeros" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let model = zinbwave_test_model()
  match model.gene_summary(counts, 5) {
    Some(summary) => {
      assert_eq(summary.gene_name, "dropout_gene")
      assert_eq(summary.observed_zero_fraction, 0.5)
      assert_true(summary.fitted_zero_fraction > 0.0)
      assert_true(summary.fitted_zero_fraction < 1.0)
      assert_true(summary.inverse_dispersion > 0.0)
    }
    None => abort("valid gene summary should exist")
  }
}

///|
test "zinbwave: gene summary rejects invalid dimensions and indices" {
  let model = zinbwave_test_model()
  assert_true(model.gene_summary([[1.0, 2.0]], 0) is None)
  let (counts, _, _, _) = @src.zinbwave_example_data()
  assert_true(model.gene_summary(counts, -1) is None)
  assert_true(model.gene_summary(counts, 6) is None)
}

///|
test "zinbwave: zero-factor model remains fully operational" {
  let (counts, genes, cells, _) = @src.zinbwave_example_data()
  let config = @src.ZinbWaveConfig::create(
    n_factors=0,
    max_iterations=5,
    inner_iterations=1,
  ) catch {
    _ => abort("zero-factor configuration should succeed")
  }
  let model = @src.zinbwave_fit(
    counts,
    config~,
    gene_names=genes,
    cell_names=cells,
  ) catch {
    _ => abort("zero-factor fit should succeed")
  }
  assert_eq(model.n_factors(), 0)
  assert_eq(model.factors.length(), 8)
  assert_eq(model.factors[0].length(), 0)
  assert_true(model.log_likelihood.abs() < 1.0e300)
}

///|
test "zinbwave: explicit mean offsets are retained" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let offsets = zinbwave_test_matrix(6, 8, 0.25)
  let config = @src.ZinbWaveConfig::create(n_factors=1, max_iterations=3) catch {
    _ => abort("configuration should succeed")
  }
  let model = @src.zinbwave_fit(counts, config~, mean_offsets=offsets) catch {
    _ => abort("offset fit should succeed")
  }
  assert_eq(model.mean_offsets[0][0], 0.25)
  assert_eq(model.mean_offsets[5][7], 0.25)
}

///|
test "zinbwave: explicit zero offsets are retained" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let offsets = zinbwave_test_matrix(6, 8, -0.5)
  let config = @src.ZinbWaveConfig::create(n_factors=1, max_iterations=3) catch {
    _ => abort("configuration should succeed")
  }
  let model = @src.zinbwave_fit(counts, config~, zero_offsets=offsets) catch {
    _ => abort("zero-offset fit should succeed")
  }
  assert_eq(model.zero_offsets[0][0], -0.5)
  assert_eq(model.zero_offsets[5][7], -0.5)
}

///|
test "zinbwave: gene-level covariates are fitted with intercept" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let gene_covariates = [[0.30], [0.35], [0.40], [0.45], [0.50], [0.55]]
  let config = @src.ZinbWaveConfig::create(n_factors=1, max_iterations=4) catch {
    _ => abort("configuration should succeed")
  }
  let model = @src.zinbwave_fit(counts, config~, gene_covariates~) catch {
    _ => abort("gene covariate fit should succeed")
  }
  assert_eq(model.gene_design[0], [1.0, 0.30])
  assert_eq(model.gamma_mu.length(), 8)
  assert_eq(model.gamma_mu[0].length(), 2)
}

///|
test "zinbwave: empty count matrix is rejected" {
  let rejected = try {
    ignore(@src.zinbwave_fit([]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}

///|
test "zinbwave: one-cell count matrix is rejected" {
  let rejected = try {
    ignore(@src.zinbwave_fit([[1.0], [2.0]]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}

///|
test "zinbwave: ragged count matrix is rejected" {
  let rejected = try {
    ignore(@src.zinbwave_fit([[1.0, 2.0], [3.0]]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}

///|
test "zinbwave: negative and non-finite counts are rejected" {
  let negative = try {
    ignore(@src.zinbwave_fit([[1.0, -1.0], [2.0, 3.0]]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let non_finite = try {
    ignore(@src.zinbwave_fit([[1.0, @double.not_a_number], [2.0, 3.0]]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(negative)
  assert_true(non_finite)
}

///|
test "zinbwave: fractional counts are rejected" {
  let rejected = try {
    ignore(@src.zinbwave_fit([[1.5, 2.0], [3.0, 4.0]]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}

///|
test "zinbwave: empty cell libraries are rejected for automatic offsets" {
  let rejected = try {
    ignore(@src.zinbwave_fit([[1.0, 0.0], [2.0, 0.0]]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}

///|
test "zinbwave: excessive factor count is rejected" {
  let config = @src.ZinbWaveConfig::create(n_factors=2) catch {
    _ => abort("configuration itself should succeed")
  }
  let rejected = try {
    ignore(@src.zinbwave_fit([[2.0, 1.0], [1.0, 2.0]], config~))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}

///|
test "zinbwave: malformed cell covariates are rejected" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let wrong_rows = try {
    ignore(@src.zinbwave_fit(counts, cell_covariates=[[0.0], [1.0]]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let ragged = try {
    ignore(
      @src.zinbwave_fit(counts, cell_covariates=[
        [0.0],
        [0.0, 1.0],
        [0.0],
        [0.0],
        [1.0],
        [1.0],
        [1.0],
        [1.0],
      ]),
    )
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(wrong_rows)
  assert_true(ragged)
}

///|
test "zinbwave: malformed gene covariates are rejected" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let rejected = try {
    ignore(@src.zinbwave_fit(counts, gene_covariates=[[0.0], [1.0], [2.0]]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}

///|
test "zinbwave: malformed offsets are rejected" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let rows = try {
    ignore(
      @src.zinbwave_fit(counts, mean_offsets=zinbwave_test_matrix(5, 8, 0.0)),
    )
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let columns = try {
    ignore(
      @src.zinbwave_fit(counts, zero_offsets=zinbwave_test_matrix(6, 7, 0.0)),
    )
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rows)
  assert_true(columns)
}

///|
test "zinbwave: non-finite offsets are rejected" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let offsets = zinbwave_test_matrix(6, 8, 0.0)
  offsets[2][3] = @double.not_a_number
  let rejected = try {
    ignore(@src.zinbwave_fit(counts, mean_offsets=offsets))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}

///|
test "zinbwave: identifier length mismatch is rejected" {
  let (counts, _, _, _) = @src.zinbwave_example_data()
  let genes = try {
    ignore(@src.zinbwave_fit(counts, gene_names=["A"]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let cells = try {
    ignore(@src.zinbwave_fit(counts, cell_names=["A"]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(genes)
  assert_true(cells)
}

///|
test "zinbwave: duplicate and empty identifiers are rejected" {
  let counts = [[2.0, 1.0], [1.0, 2.0]]
  let config = @src.ZinbWaveConfig::create(n_factors=0) catch {
    _ => abort("configuration should succeed")
  }
  let duplicate = try {
    ignore(@src.zinbwave_fit(counts, config~, gene_names=["same", "same"]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  let empty = try {
    ignore(@src.zinbwave_fit(counts, config~, cell_names=["", "cell"]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(duplicate)
  assert_true(empty)
}

///|
test "zinbwave: imputation rejects matrix shape mismatch" {
  let model = zinbwave_test_model()
  let rejected = try {
    ignore(model.impute_zeros([[1.0, 2.0], [2.0, 1.0]]))
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}

///|
test "zinbwave: SCE integration adds assays and reduced dimensions" {
  let (_, _, _, batch) = @src.zinbwave_example_data()
  let config = zinbwave_test_config()
  let output = @src.zinbwave_sce(
    zinbwave_test_sce(),
    config~,
    cell_covariates=batch,
  ) catch {
    _ => abort("SCE integration should succeed")
  }
  assert_eq(
    @src.sce_get_assay(output.experiment, "zinbwave_weights").length(),
    6,
  )
  assert_eq(
    @src.sce_get_assay(output.experiment, "zinbwave_residuals").length(),
    6,
  )
  assert_eq(
    @src.sce_get_assay(output.experiment, "zinbwave_normalized").length(),
    6,
  )
  assert_eq(
    @src.sce_get_assay(output.experiment, "zinbwave_imputed").length(),
    6,
  )
  assert_eq(@src.sce_get_reduced_dim(output.experiment, "zinbwave").length(), 8)
}

///|
test "zinbwave: SCE integration preserves source object" {
  let source = zinbwave_test_sce()
  let config = zinbwave_test_config()
  ignore(@src.zinbwave_sce(source, config~)) catch {
    _ => abort("SCE integration should succeed")
  }
  assert_eq(@src.sce_get_assay(source, "zinbwave_weights").length(), 0)
  assert_eq(@src.sce_get_reduced_dim(source, "zinbwave").length(), 0)
}

///|
test "zinbwave: SCE integration supports custom reduced dimension name" {
  let config = zinbwave_test_config()
  let output = @src.zinbwave_sce(
    zinbwave_test_sce(),
    reduced_dim_name="ZINB",
    config~,
  ) catch {
    _ => abort("custom reduced dimension should succeed")
  }
  assert_eq(@src.sce_get_reduced_dim(output.experiment, "ZINB").length(), 8)
  assert_eq(@src.sce_get_reduced_dim(output.experiment, "zinbwave").length(), 0)
}

///|
test "zinbwave: SCE integration retains row and cell identifiers" {
  let config = zinbwave_test_config()
  let output = @src.zinbwave_sce(zinbwave_test_sce(), config~) catch {
    _ => abort("SCE integration should succeed")
  }
  assert_eq(output.model.gene_names[0], "A_marker_1")
  assert_eq(output.model.cell_names[7], "B4")
  assert_eq(output.experiment.row_names[5], "dropout_gene")
}

///|
test "zinbwave: SCE integration rejects missing assay" {
  let config = zinbwave_test_config()
  let rejected = try {
    ignore(
      @src.zinbwave_sce(zinbwave_test_sce(), input_assay="missing", config~),
    )
    false
  } catch {
    ZinbWaveError(_) => true
  }
  assert_true(rejected)
}
