// Tests for Bioconductor celda_CG-inspired joint cell and feature clustering.

///|
fn celda_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "celda_CG value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn celda_test_sum(values : Array[Double]) -> Double {
  let mut total = 0.0
  for value in values {
    total = total + value
  }
  total
}

///|
fn celda_test_config() -> @src.CeldaCGConfig {
  @src.CeldaCGConfig::create(
    3,
    4,
    max_iterations=15,
    stop_iterations=3,
    chains=2,
    seed=2026,
  ) catch {
    _ => abort("valid celda_CG configuration should build")
  }
}

///|
fn celda_test_result() -> @src.CeldaCGResult {
  let (counts, features, cells, samples) = @src.celda_cg_example_data()
  @src.celda_cg_fit(
    counts,
    samples,
    celda_test_config(),
    feature_names=features,
    cell_names=cells,
  ) catch {
    _ => abort("example celda_CG analysis should succeed")
  }
}

///|
fn celda_test_sce() -> @src.SingleCellExperiment {
  let (counts, features, cells, samples) = @src.celda_cg_example_data()
  let sce = @src.SingleCellExperiment::new(counts, features, cells)
  sce.col_data["donor"] = samples
  sce.metadata["source"] = "celda_test"
  sce
}

///|
test "celda_CG: default configuration preserves upstream hyperparameters" {
  let config = @src.CeldaCGConfig::create(3, 4) catch {
    _ => abort("default configuration should build")
  }
  assert_eq(config.cell_populations, 3)
  assert_eq(config.feature_modules, 4)
  assert_eq(config.alpha, 1.0)
  assert_eq(config.beta, 1.0)
  assert_eq(config.delta, 1.0)
  assert_eq(config.gamma, 1.0)
  assert_eq(config.algorithm, @src.CeldaEm)
  assert_eq(config.max_iterations, 200)
  assert_eq(config.stop_iterations, 10)
  assert_eq(config.chains, 3)
  assert_eq(config.tolerance, 1.0e-8)
  assert_eq(config.seed, 12345)
}

///|
test "celda_CG: custom configuration preserves every option" {
  let config = @src.CeldaCGConfig::create(
    2,
    3,
    alpha=0.5,
    beta=0.25,
    delta=0.75,
    gamma=2.0,
    algorithm=@src.CeldaGibbs,
    max_iterations=20,
    stop_iterations=4,
    chains=5,
    tolerance=1.0e-6,
    seed=99,
  ) catch {
    _ => abort("custom configuration should build")
  }
  assert_eq(config.alpha, 0.5)
  assert_eq(config.beta, 0.25)
  assert_eq(config.delta, 0.75)
  assert_eq(config.gamma, 2.0)
  assert_eq(config.algorithm, @src.CeldaGibbs)
  assert_eq(config.max_iterations, 20)
  assert_eq(config.stop_iterations, 4)
  assert_eq(config.chains, 5)
  assert_eq(config.tolerance, 1.0e-6)
  assert_eq(config.seed, 99)
}

///|
test "celda_CG: with_dimensions retains inference controls" {
  let original = celda_test_config()
  let updated = original.with_dimensions(2, 5, seed=77) catch {
    _ => abort("valid dimensions should build")
  }
  assert_eq(updated.cell_populations, 2)
  assert_eq(updated.feature_modules, 5)
  assert_eq(updated.max_iterations, original.max_iterations)
  assert_eq(updated.stop_iterations, original.stop_iterations)
  assert_eq(updated.chains, original.chains)
  assert_eq(updated.seed, 77)
}

///|
test "celda_CG: configuration rejects non-positive dimensions" {
  let bad_k = try {
    ignore(@src.CeldaCGConfig::create(0, 2))
    false
  } catch {
    CeldaError(_) => true
  }
  let bad_l = try {
    ignore(@src.CeldaCGConfig::create(2, 0))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(bad_k)
  assert_true(bad_l)
}

///|
test "celda_CG: configuration rejects invalid Dirichlet priors" {
  let alpha = try {
    ignore(@src.CeldaCGConfig::create(2, 2, alpha=0.0))
    false
  } catch {
    CeldaError(_) => true
  }
  let beta = try {
    ignore(@src.CeldaCGConfig::create(2, 2, beta=-1.0))
    false
  } catch {
    CeldaError(_) => true
  }
  let delta = try {
    ignore(@src.CeldaCGConfig::create(2, 2, delta=Double::nan()))
    false
  } catch {
    CeldaError(_) => true
  }
  let gamma = try {
    ignore(@src.CeldaCGConfig::create(2, 2, gamma=1.0e301))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(alpha)
  assert_true(beta)
  assert_true(delta)
  assert_true(gamma)
}

///|
test "celda_CG: configuration rejects invalid iteration limits" {
  let maximum = try {
    ignore(@src.CeldaCGConfig::create(2, 2, max_iterations=0))
    false
  } catch {
    CeldaError(_) => true
  }
  let stopping = try {
    ignore(@src.CeldaCGConfig::create(2, 2, stop_iterations=0))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(maximum)
  assert_true(stopping)
}

///|
test "celda_CG: configuration rejects invalid chains and tolerance" {
  let chains = try {
    ignore(@src.CeldaCGConfig::create(2, 2, chains=0))
    false
  } catch {
    CeldaError(_) => true
  }
  let tolerance = try {
    ignore(@src.CeldaCGConfig::create(2, 2, tolerance=0.0))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(chains)
  assert_true(tolerance)
}

///|
test "celda_CG: configuration rejects non-positive random seeds" {
  let rejected = try {
    ignore(@src.CeldaCGConfig::create(2, 2, seed=0))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: example data uses feature by cell orientation" {
  let (counts, features, cells, samples) = @src.celda_cg_example_data()
  assert_eq(counts.length(), 12)
  assert_eq(counts[0].length(), 12)
  assert_eq(features.length(), 12)
  assert_eq(cells.length(), 12)
  assert_eq(samples.length(), 12)
  assert_eq(features[0], "T_marker_1")
  assert_eq(cells[8], "M1")
}

///|
test "celda_CG: fit rejects an empty count matrix" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects a matrix without cells" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects a non-rectangular matrix" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, 2.0], [3.0]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects negative counts" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, -1.0]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects fractional counts" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, 1.5]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects NaN counts" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, Double::nan()]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects effectively infinite counts" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, 1.0e301]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects all-zero features" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, 2.0], [0.0, 0.0]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects all-zero cells" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, 0.0], [2.0, 0.0]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects more populations than cells" {
  let config = @src.CeldaCGConfig::create(3, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, 2.0]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: fit rejects more modules than features" {
  let config = @src.CeldaCGConfig::create(1, 2) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, 2.0]], [], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: sample labels must match the number of cells" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, 2.0]], ["S1"], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: sample labels cannot be blank" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_fit([[1.0, 2.0]], ["S1", " "], config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: sample order follows first appearance" {
  let config = @src.CeldaCGConfig::create(
    1,
    1,
    max_iterations=1,
    stop_iterations=1,
    chains=1,
  ) catch {
    _ => abort("configuration should build")
  }
  let result = @src.celda_cg_fit([[2.0, 3.0, 4.0]], ["B", "A", "B"], config) catch {
    _ => abort("valid sample labels should fit")
  }
  assert_eq(result.sample_names, ["B", "A"])
  assert_eq(result.sample_indices, [0, 1, 0])
}

///|
test "celda_CG: omitted identifiers are generated deterministically" {
  let config = @src.CeldaCGConfig::create(
    1,
    1,
    max_iterations=1,
    stop_iterations=1,
    chains=1,
  ) catch {
    _ => abort("configuration should build")
  }
  let result = @src.celda_cg_fit([[2.0, 3.0], [3.0, 2.0]], [], config) catch {
    _ => abort("generated names should be accepted")
  }
  assert_eq(result.feature_names, ["Feature1", "Feature2"])
  assert_eq(result.cell_names, ["Cell1", "Cell2"])
  assert_eq(result.sample_names, ["Sample1"])
}

///|
test "celda_CG: feature identifiers must match the matrix" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(
      @src.celda_cg_fit([[1.0, 2.0], [2.0, 1.0]], [], config, feature_names=[
        "G1",
      ]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: cell identifiers must match the matrix" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(
      @src.celda_cg_fit([[1.0, 2.0], [2.0, 1.0]], [], config, cell_names=["C1"]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: identifiers cannot be blank" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(
      @src.celda_cg_fit([[1.0, 2.0], [2.0, 1.0]], [], config, feature_names=[
        "G1", " ",
      ]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: identifiers must be unique after trimming" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(
      @src.celda_cg_fit([[1.0, 2.0], [2.0, 1.0]], [], config, feature_names=[
        "G1", " G1 ",
      ]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: predefined cell labels must match and use all populations" {
  let config = @src.CeldaCGConfig::create(2, 1) catch {
    _ => abort("configuration should build")
  }
  let wrong_length = try {
    ignore(
      @src.celda_cg_fit([[2.0, 3.0, 4.0]], [], config, initial_cell_clusters=[
        0, 1,
      ]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  let missing = try {
    ignore(
      @src.celda_cg_fit([[2.0, 3.0, 4.0]], [], config, initial_cell_clusters=[
        0, 0, 0,
      ]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(wrong_length)
  assert_true(missing)
}

///|
test "celda_CG: predefined cell labels reject out-of-range values" {
  let config = @src.CeldaCGConfig::create(2, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(
      @src.celda_cg_fit([[2.0, 3.0, 4.0]], [], config, initial_cell_clusters=[
        0, 1, 2,
      ]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: predefined feature labels must match and use all modules" {
  let config = @src.CeldaCGConfig::create(1, 2) catch {
    _ => abort("configuration should build")
  }
  let wrong_length = try {
    ignore(
      @src.celda_cg_fit([[2.0, 3.0], [3.0, 2.0], [1.0, 1.0]], [], config, initial_feature_modules=[
        0, 1,
      ]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  let missing = try {
    ignore(
      @src.celda_cg_fit([[2.0, 3.0], [3.0, 2.0], [1.0, 1.0]], [], config, initial_feature_modules=[
        0, 0, 0,
      ]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(wrong_length)
  assert_true(missing)
}

///|
test "celda_CG: predefined feature labels reject out-of-range values" {
  let config = @src.CeldaCGConfig::create(1, 2) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(
      @src.celda_cg_fit([[2.0, 3.0], [3.0, 2.0], [1.0, 1.0]], [], config, initial_feature_modules=[
        0, 1, -1,
      ]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: one-population one-module likelihood matches closed form" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let value = @src.celda_cg_log_likelihood(
    [[2.0, 1.0], [1.0, 2.0]],
    [],
    [0, 0],
    [0, 0],
    config,
  ) catch {
    _ => abort("valid labels should produce a likelihood")
  }
  celda_test_close(value, -6.327936783729195, 1.0e-9)
}

///|
test "celda_CG: omitted and explicit default samples have equal likelihood" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let counts = [[2.0, 1.0], [1.0, 2.0]]
  let implicit = @src.celda_cg_log_likelihood(
    counts,
    [],
    [0, 0],
    [0, 0],
    config,
  ) catch {
    _ => abort("implicit sample should work")
  }
  let explicit = @src.celda_cg_log_likelihood(
    counts,
    ["Sample1", "Sample1"],
    [0, 0],
    [0, 0],
    config,
  ) catch {
    _ => abort("explicit sample should work")
  }
  celda_test_close(implicit, explicit, 1.0e-12)
}

///|
test "celda_CG: likelihood validates assignment occupancy" {
  let config = @src.CeldaCGConfig::create(2, 2) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(
      @src.celda_cg_log_likelihood(
        [[2.0, 1.0], [1.0, 2.0]],
        [],
        [0, 0],
        [0, 1],
        config,
      ),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: result dimensions match the example matrix" {
  let result = celda_test_result()
  assert_eq(result.n_features(), 12)
  assert_eq(result.n_cells(), 12)
  assert_eq(result.n_samples(), 2)
  assert_eq(result.cell_clusters.length(), 12)
  assert_eq(result.feature_modules.length(), 12)
}

///|
test "celda_CG: model selection metrics are finite" {
  let result = celda_test_result()
  assert_true(result.log_likelihood == result.log_likelihood)
  assert_true(result.log_likelihood.abs() < 1.0e300)
  assert_true(result.perplexity > 1.0)
  assert_true(result.perplexity < 1.0e300)
  assert_true(result.aic > 0.0)
  assert_true(result.bic > 0.0)
}

///|
test "celda_CG: diagnostics record bounded optimization work" {
  let result = celda_test_result()
  assert_true(result.diagnostics.iterations >= 1)
  assert_true(result.diagnostics.iterations <= result.config.max_iterations)
  assert_eq(
    result.diagnostics.log_likelihoods.length(),
    result.diagnostics.iterations + 1,
  )
  assert_true(result.diagnostics.best_chain >= 0)
  assert_true(result.diagnostics.best_chain < result.config.chains)
}

///|
test "celda_CG: diagnostics include one finite score per chain" {
  let result = celda_test_result()
  assert_eq(result.diagnostics.chain_scores.length(), result.config.chains)
  for score in result.diagnostics.chain_scores {
    assert_true(score == score)
    assert_true(score.abs() < 1.0e300)
  }
}

///|
test "celda_CG: hard EM is reproducible for a fixed seed" {
  let (counts, features, cells, samples) = @src.celda_cg_example_data()
  let first = @src.celda_cg_fit(
    counts,
    samples,
    celda_test_config(),
    feature_names=features,
    cell_names=cells,
  ) catch {
    _ => abort("first fit should succeed")
  }
  let second = @src.celda_cg_fit(
    counts,
    samples,
    celda_test_config(),
    feature_names=features,
    cell_names=cells,
  ) catch {
    _ => abort("second fit should succeed")
  }
  assert_eq(first.cell_clusters, second.cell_clusters)
  assert_eq(first.feature_modules, second.feature_modules)
  assert_eq(first.log_likelihood, second.log_likelihood)
}

///|
test "celda_CG: seeded Gibbs inference is reproducible" {
  let (counts, features, cells, samples) = @src.celda_cg_example_data()
  let config = @src.CeldaCGConfig::create(
    3,
    4,
    algorithm=@src.CeldaGibbs,
    max_iterations=6,
    stop_iterations=3,
    chains=2,
    seed=456,
  ) catch {
    _ => abort("Gibbs configuration should build")
  }
  let first = @src.celda_cg_fit(
    counts,
    samples,
    config,
    feature_names=features,
    cell_names=cells,
  ) catch {
    _ => abort("first Gibbs fit should succeed")
  }
  let second = @src.celda_cg_fit(
    counts,
    samples,
    config,
    feature_names=features,
    cell_names=cells,
  ) catch {
    _ => abort("second Gibbs fit should succeed")
  }
  assert_eq(first.cell_clusters, second.cell_clusters)
  assert_eq(first.feature_modules, second.feature_modules)
  assert_eq(first.log_likelihood, second.log_likelihood)
}

///|
test "celda_CG: every inferred population remains occupied" {
  let result = celda_test_result()
  let sizes = result.population_sizes()
  assert_eq(sizes.length(), 3)
  for size in sizes {
    assert_true(size > 0)
  }
  assert_eq(celda_test_sum(sizes.map(fn(value) { value.to_double() })), 12.0)
}

///|
test "celda_CG: every inferred feature module remains occupied" {
  let result = celda_test_result()
  let sizes = result.module_sizes()
  assert_eq(sizes.length(), 4)
  let mut total = 0
  for size in sizes {
    assert_true(size > 0)
    total = total + size
  }
  assert_eq(total, 12)
}

///|
test "celda_CG: example cell types form three distinct populations" {
  let result = celda_test_result()
  for index in 1..<4 {
    assert_eq(result.cell_clusters[index], result.cell_clusters[0])
  }
  for index in 5..<8 {
    assert_eq(result.cell_clusters[index], result.cell_clusters[4])
  }
  for index in 9..<12 {
    assert_eq(result.cell_clusters[index], result.cell_clusters[8])
  }
  assert_true(result.cell_clusters[0] != result.cell_clusters[4])
  assert_true(result.cell_clusters[0] != result.cell_clusters[8])
  assert_true(result.cell_clusters[4] != result.cell_clusters[8])
}

///|
test "celda_CG: example marker families form distinct modules" {
  let result = celda_test_result()
  assert_eq(result.feature_modules[0], result.feature_modules[1])
  assert_eq(result.feature_modules[1], result.feature_modules[2])
  assert_eq(result.feature_modules[3], result.feature_modules[4])
  assert_eq(result.feature_modules[4], result.feature_modules[5])
  assert_eq(result.feature_modules[6], result.feature_modules[7])
  assert_eq(result.feature_modules[7], result.feature_modules[8])
  assert_eq(result.feature_modules[9], result.feature_modules[10])
  assert_eq(result.feature_modules[10], result.feature_modules[11])
}

///|
test "celda_CG: theta is population by sample and column-normalized" {
  let result = celda_test_result()
  assert_eq(result.theta.length(), 3)
  assert_eq(result.theta[0].length(), 2)
  for sample in 0..<result.n_samples() {
    let mut total = 0.0
    for population in 0..<result.config.cell_populations {
      total = total + result.theta[population][sample]
    }
    celda_test_close(total, 1.0, 1.0e-12)
  }
}

///|
test "celda_CG: phi is module by population and column-normalized" {
  let result = celda_test_result()
  assert_eq(result.phi.length(), 4)
  assert_eq(result.phi[0].length(), 3)
  for population in 0..<result.config.cell_populations {
    let mut total = 0.0
    for module_ in 0..<result.config.feature_modules {
      total = total + result.phi[module_][population]
    }
    celda_test_close(total, 1.0, 1.0e-12)
  }
}

///|
test "celda_CG: psi stores probability only in the assigned module" {
  let result = celda_test_result()
  assert_eq(result.psi.length(), 12)
  assert_eq(result.psi[0].length(), 4)
  for feature in 0..<result.n_features() {
    for module_ in 0..<result.config.feature_modules {
      if module_ == result.feature_modules[feature] {
        assert_true(result.psi[feature][module_] > 0.0)
      } else {
        assert_eq(result.psi[feature][module_], 0.0)
      }
    }
  }
}

///|
test "celda_CG: psi reserves pseudogene probability in every module" {
  let result = celda_test_result()
  for module_ in 0..<result.config.feature_modules {
    let mut total = 0.0
    for feature in 0..<result.n_features() {
      total = total + result.psi[feature][module_]
    }
    assert_true(total > 0.0)
    assert_true(total < 1.0)
  }
}

///|
test "celda_CG: eta is a normalized feature-module abundance" {
  let result = celda_test_result()
  assert_eq(result.eta.length(), 4)
  for probability in result.eta {
    assert_true(probability > 0.0)
  }
  celda_test_close(celda_test_sum(result.eta), 1.0, 1.0e-12)
}

///|
test "celda_CG: fitted counts match source dimensions" {
  let result = celda_test_result()
  assert_eq(result.fitted_counts.length(), 12)
  assert_eq(result.fitted_counts[0].length(), 12)
  for row in result.fitted_counts {
    for value in row {
      assert_true(value >= 0.0)
      assert_true(value.abs() < 1.0e300)
    }
  }
}

///|
test "celda_CG: fitted counts preserve every cell library size" {
  let result = celda_test_result()
  for cell in 0..<result.n_cells() {
    let mut total = 0.0
    for feature in 0..<result.n_features() {
      total = total + result.fitted_counts[feature][cell]
    }
    celda_test_close(total, result.cell_totals[cell], 1.0e-9)
  }
}

///|
test "celda_CG: feature totals reconstruct the source rows" {
  let (counts, _, _, _) = @src.celda_cg_example_data()
  let result = celda_test_result()
  for feature in 0..<counts.length() {
    celda_test_close(
      result.feature_totals[feature],
      celda_test_sum(counts[feature]),
      1.0e-12,
    )
  }
}

///|
test "celda_CG: cell totals reconstruct the source columns" {
  let (counts, _, _, _) = @src.celda_cg_example_data()
  let result = celda_test_result()
  for cell in 0..<counts[0].length() {
    let mut expected = 0.0
    for feature in 0..<counts.length() {
      expected = expected + counts[feature][cell]
    }
    celda_test_close(result.cell_totals[cell], expected, 1.0e-12)
  }
}

///|
test "celda_CG: name queries return assignments and None for misses" {
  let result = celda_test_result()
  assert_eq(result.cell_population("T1"), Some(result.cell_clusters[0]))
  assert_eq(
    result.feature_module("B_marker_1"),
    Some(result.feature_modules[3]),
  )
  assert_eq(result.cell_population("missing"), None)
  assert_eq(result.feature_module("missing"), None)
}

///|
test "celda_CG: population member query agrees with population sizes" {
  let result = celda_test_result()
  for population in 0..<result.config.cell_populations {
    assert_eq(
      result.cells_in_population(population).length(),
      result.population_sizes()[population],
    )
  }
}

///|
test "celda_CG: module member query agrees with module sizes" {
  let result = celda_test_result()
  for module_ in 0..<result.config.feature_modules {
    assert_eq(
      result.features_in_module(module_).length(),
      result.module_sizes()[module_],
    )
  }
}

///|
test "celda_CG: member queries reject invalid indices" {
  let result = celda_test_result()
  let population = try {
    ignore(result.cells_in_population(3))
    false
  } catch {
    CeldaError(_) => true
  }
  let module_ = try {
    ignore(result.features_in_module(-1))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(population)
  assert_true(module_)
}

///|
test "celda_CG: top features are sorted by posterior probability" {
  let result = celda_test_result()
  let module_ = result.feature_modules[0]
  let scores = result.top_features(module_, limit=3) catch {
    _ => abort("valid module should have top features")
  }
  assert_true(scores.length() >= 1)
  for index in 1..<scores.length() {
    assert_true(scores[index - 1].probability >= scores[index].probability)
  }
  for score in scores {
    assert_eq(result.feature_modules[score.feature_index], module_)
    assert_eq(result.feature_names[score.feature_index], score.feature_name)
  }
}

///|
test "celda_CG: top feature query validates module and limit" {
  let result = celda_test_result()
  let module_ = try {
    ignore(result.top_features(4))
    false
  } catch {
    CeldaError(_) => true
  }
  let limit = try {
    ignore(result.top_features(0, limit=0))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(module_)
  assert_true(limit)
}

///|
test "celda_CG: summary reports dimensions and model family" {
  let summary = celda_test_result().summary()
  assert_true(summary.contains("celda_CG"))
  assert_true(summary.contains("12 features x 12 cells"))
  assert_true(summary.contains("K=3"))
  assert_true(summary.contains("L=4"))
}

///|
test "celda_CG: prediction returns cell by population probabilities" {
  let (counts, _, _, _) = @src.celda_cg_example_data()
  let result = celda_test_result()
  let prediction = @src.celda_cg_predict_cells(result, counts) catch {
    _ => abort("training matrix should be predictable")
  }
  assert_eq(prediction.probabilities.length(), 12)
  assert_eq(prediction.probabilities[0].length(), 3)
  assert_eq(prediction.assignments.length(), 12)
}

///|
test "celda_CG: prediction probabilities are normalized" {
  let (counts, _, _, _) = @src.celda_cg_example_data()
  let prediction = @src.celda_cg_predict_cells(celda_test_result(), counts) catch {
    _ => abort("training matrix should be predictable")
  }
  for probabilities in prediction.probabilities {
    celda_test_close(celda_test_sum(probabilities), 1.0, 1.0e-12)
    for probability in probabilities {
      assert_true(probability >= 0.0)
      assert_true(probability <= 1.0)
    }
  }
}

///|
test "celda_CG: training cells predict their inferred populations" {
  let (counts, _, _, _) = @src.celda_cg_example_data()
  let result = celda_test_result()
  let prediction = @src.celda_cg_predict_cells(result, counts) catch {
    _ => abort("training matrix should be predictable")
  }
  assert_eq(prediction.assignments, result.cell_clusters)
}

///|
test "celda_CG: prediction rejects a mismatched feature space" {
  let rejected = try {
    ignore(
      @src.celda_cg_predict_cells(celda_test_result(), [[1.0, 2.0], [2.0, 1.0]]),
    )
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: grid search evaluates the Cartesian K by L grid" {
  let (counts, features, cells, samples) = @src.celda_cg_example_data()
  let base = @src.CeldaCGConfig::create(
    1,
    1,
    max_iterations=5,
    stop_iterations=2,
    chains=1,
    seed=41,
  ) catch {
    _ => abort("base grid configuration should build")
  }
  let grid = @src.celda_cg_grid_search(
    counts,
    samples,
    [2, 3],
    [3, 4],
    base,
    feature_names=features,
    cell_names=cells,
  ) catch {
    _ => abort("valid grid should fit")
  }
  assert_eq(grid.entries.length(), 4)
  assert_true(grid.best_index >= 0)
  assert_true(grid.best_index < 4)
  assert_eq(grid.criterion, @src.CeldaBic)
}

///|
test "celda_CG: grid best_result matches the selected entry" {
  let (counts, _, _, samples) = @src.celda_cg_example_data()
  let base = @src.CeldaCGConfig::create(
    1,
    1,
    max_iterations=3,
    stop_iterations=1,
    chains=1,
  ) catch {
    _ => abort("base grid configuration should build")
  }
  let grid = @src.celda_cg_grid_search(
    counts,
    samples,
    [2, 3],
    [4],
    base,
    criterion=@src.CeldaLogLikelihood,
  ) catch {
    _ => abort("valid likelihood grid should fit")
  }
  let best = grid.best_result()
  assert_eq(best.log_likelihood, grid.entries[grid.best_index].score)
}

///|
test "celda_CG: grid search rejects empty dimensions" {
  let (counts, _, _, samples) = @src.celda_cg_example_data()
  let base = celda_test_config()
  let rejected = try {
    ignore(@src.celda_cg_grid_search(counts, samples, [], [2], base))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: grid search rejects duplicate K-L combinations" {
  let (counts, _, _, samples) = @src.celda_cg_example_data()
  let base = celda_test_config()
  let rejected = try {
    ignore(@src.celda_cg_grid_search(counts, samples, [2, 2], [3], base))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: SCE integration adds fitted assay and assignments" {
  let output = @src.celda_cg_sce(
    celda_test_sce(),
    celda_test_config(),
    sample_column="donor",
  ) catch {
    _ => abort("valid SCE integration should succeed")
  }
  assert_eq(@src.sce_get_assay(output.experiment, "celda_fitted").length(), 12)
  assert_eq(
    @src.sce_get_row_data(output.experiment, "celda_feature_module").length(),
    12,
  )
  assert_eq(
    @src.sce_get_col_data(output.experiment, "celda_cell_population").length(),
    12,
  )
}

///|
test "celda_CG: SCE integration writes sample and model metadata" {
  let output = @src.celda_cg_sce(
    celda_test_sce(),
    celda_test_config(),
    sample_column="donor",
  ) catch {
    _ => abort("valid SCE integration should succeed")
  }
  assert_eq(@src.sce_get_col_data(output.experiment, "celda_sample_label"), [
    "Donor1", "Donor1", "Donor2", "Donor2", "Donor1", "Donor1", "Donor2", "Donor2",
    "Donor1", "Donor1", "Donor2", "Donor2",
  ])
  assert_eq(output.experiment.metadata["celda_model"], "celda_CG")
  assert_eq(output.experiment.metadata["celda_K"], "3")
  assert_eq(output.experiment.metadata["celda_L"], "4")
  assert_eq(output.experiment.metadata["source"], "celda_test")
}

///|
test "celda_CG: SCE integration does not mutate the source object" {
  let source = celda_test_sce()
  ignore(
    @src.celda_cg_sce(source, celda_test_config(), sample_column="donor") catch {
      _ => abort("valid SCE integration should succeed")
    },
  )
  assert_eq(@src.sce_get_assay(source, "celda_fitted").length(), 0)
  assert_eq(@src.sce_get_row_data(source, "celda_feature_module").length(), 0)
  assert_eq(@src.sce_get_col_data(source, "celda_cell_population").length(), 0)
}

///|
test "celda_CG: SCE integration supports custom assay names" {
  let output = @src.celda_cg_sce(
    celda_test_sce(),
    celda_test_config(),
    sample_column="donor",
    fitted_assay="celda_reconstructed",
  ) catch {
    _ => abort("custom fitted assay should work")
  }
  assert_eq(
    @src.sce_get_assay(output.experiment, "celda_reconstructed").length(),
    12,
  )
  assert_eq(@src.sce_get_assay(output.experiment, "celda_fitted").length(), 0)
}

///|
test "celda_CG: SCE integration validates the input assay" {
  let source = @src.SCEBuilder::new()
    |> @src.SCEBuilder::add_assay("other", [[1.0, 2.0], [2.0, 1.0]])
    |> @src.SCEBuilder::build
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_sce(source, config))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}

///|
test "celda_CG: SCE integration validates sample column metadata" {
  let config = @src.CeldaCGConfig::create(1, 1) catch {
    _ => abort("configuration should build")
  }
  let rejected = try {
    ignore(@src.celda_cg_sce(celda_test_sce(), config, sample_column="missing"))
    false
  } catch {
    CeldaError(_) => true
  }
  assert_true(rejected)
}
