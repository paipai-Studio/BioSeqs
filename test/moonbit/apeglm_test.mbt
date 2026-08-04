// Tests for the Bioconductor apeglm-inspired heavy-tailed LFC shrinkage.

///|
fn apeglm_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "apeglm value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn apeglm_test_config() -> @src.ApeglmConfig {
  @src.ApeglmConfig::create(
    coefficient=2,
    threshold=@math.ln(2.0),
  ) catch {
    _ => abort("test configuration should be valid")
  }
}

///|
fn apeglm_test_result() -> @src.ApeglmResult {
  let (counts, design, dispersions, genes, coefficients) =
    @src.apeglm_example_data()
  let config = apeglm_test_config()
  @src.apeglm_fit(
    counts,
    design,
    dispersions,
    config~,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("example apeglm fit should succeed")
  }
}

///|
fn apeglm_test_gene(
  result : @src.ApeglmResult,
  name : String,
) -> @src.ApeglmGeneResult {
  match result.find_gene(name) {
    Some(gene) => gene
    None => abort("expected apeglm gene: " + name)
  }
}

///|
fn apeglm_test_zero_matrix(
  rows : Int,
  columns : Int,
) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = []
  for _ in 0..<rows {
    matrix.push(Array::make(columns, 0.0))
  }
  matrix
}

///|
fn apeglm_test_experiment() -> @src.SummarizedExperiment {
  let (counts, _, _, _, _) = @src.apeglm_example_data()
  let assays : Map[String, Array[Array[Double]]] = Map([("counts", counts)])
  let metadata : Map[String, String] = Map([("source", "apeglm_test")])
  @src.summarized_experiment(assays, [], [], metadata)
}

///|
test "apeglm: default configuration exposes inference controls" {
  let config = @src.ApeglmConfig::default()
  assert_eq(config.coefficient, -1)
  assert_eq(config.threshold, 0.0)
  assert_eq(config.interval_level, 0.95)
  assert_eq(config.prior_scale, 1.0)
  assert_eq(config.prior_df, 1.0)
  assert_eq(config.prior_no_shrink_scale, 15.0)
  assert_true(config.adaptive_prior)
  assert_eq(config.random_starts, 5)
}

///|
test "apeglm: custom configuration preserves parameters" {
  let config = @src.ApeglmConfig::create(
    coefficient=1,
    threshold=0.5,
    interval_level=0.9,
    prior_scale=0.7,
    prior_df=3.0,
    prior_no_shrink_scale=12.0,
    adaptive_prior=false,
    multiplier=1.5,
    max_iterations=40,
    tolerance=1.0e-6,
    max_step=1.0,
    random_starts=9,
  ) catch {
    _ => abort("custom configuration should be valid")
  }
  assert_eq(config.coefficient, 1)
  assert_eq(config.threshold, 0.5)
  assert_eq(config.interval_level, 0.9)
  assert_eq(config.prior_df, 3.0)
  assert_false(config.adaptive_prior)
  assert_eq(config.max_iterations, 40)
  assert_eq(config.random_starts, 9)
}

///|
test "apeglm: configuration rejects invalid coefficient indices" {
  let raised = try {
    ignore(@src.ApeglmConfig::create(coefficient=-2))
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(raised)
}

///|
test "apeglm: configuration rejects invalid probability controls" {
  let threshold = try {
    ignore(@src.ApeglmConfig::create(threshold=-0.1))
    false
  } catch {
    ApeglmError(_) => true
  }
  let interval = try {
    ignore(@src.ApeglmConfig::create(interval_level=1.0))
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(threshold)
  assert_true(interval)
}

///|
test "apeglm: configuration rejects invalid prior controls" {
  let scale = try {
    ignore(@src.ApeglmConfig::create(prior_scale=0.0))
    false
  } catch {
    ApeglmError(_) => true
  }
  let degrees = try {
    ignore(@src.ApeglmConfig::create(prior_df=-1.0))
    false
  } catch {
    ApeglmError(_) => true
  }
  let nuisance = try {
    ignore(@src.ApeglmConfig::create(prior_no_shrink_scale=0.0))
    false
  } catch {
    ApeglmError(_) => true
  }
  let multiplier = try {
    ignore(@src.ApeglmConfig::create(multiplier=0.0))
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(scale)
  assert_true(degrees)
  assert_true(nuisance)
  assert_true(multiplier)
}

///|
test "apeglm: configuration rejects invalid optimizer controls" {
  let iterations = try {
    ignore(@src.ApeglmConfig::create(max_iterations=0))
    false
  } catch {
    ApeglmError(_) => true
  }
  let tolerance = try {
    ignore(@src.ApeglmConfig::create(tolerance=0.0))
    false
  } catch {
    ApeglmError(_) => true
  }
  let step = try {
    ignore(@src.ApeglmConfig::create(max_step=-1.0))
    false
  } catch {
    ApeglmError(_) => true
  }
  let starts = try {
    ignore(@src.ApeglmConfig::create(random_starts=0))
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(iterations)
  assert_true(tolerance)
  assert_true(step)
  assert_true(starts)
}

///|
test "apeglm: example data has gene by sample orientation" {
  let (counts, design, dispersions, genes, coefficients) =
    @src.apeglm_example_data()
  assert_eq(counts.length(), 8)
  assert_eq(counts[0].length(), 8)
  assert_eq(design.length(), 8)
  assert_eq(design[0].length(), 3)
  assert_eq(dispersions.length(), 8)
  assert_eq(genes.length(), 8)
  assert_eq(coefficients.length(), 3)
}

///|
test "apeglm: example data identifies intended signals" {
  let (_, design, _, genes, coefficients) = @src.apeglm_example_data()
  assert_eq(genes[0], "strong_up")
  assert_eq(genes[1], "strong_down")
  assert_eq(genes[3], "low_noisy")
  assert_eq(coefficients, ["intercept", "batch", "condition"])
  assert_eq(design[0], [1.0, 0.0, 0.0])
  assert_eq(design[4], [1.0, 0.0, 1.0])
}

///|
test "apeglm: fit dimensions match inputs" {
  let result = apeglm_test_result()
  assert_eq(result.genes.length(), 8)
  assert_eq(result.gene_names.length(), 8)
  assert_eq(result.coefficient_names.length(), 3)
  assert_eq(result.coefficient, 2)
  assert_eq(result.genes[0].map.length(), 3)
  assert_eq(result.genes[0].posterior_sd.length(), 3)
}

///|
test "apeglm: omitted names are generated deterministically" {
  let result = @src.apeglm_fit(
    [[10.0, 11.0, 20.0, 22.0]],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
    [0.1],
  ) catch {
    _ => abort("small apeglm fit should succeed")
  }
  assert_eq(result.gene_names, ["feature_1"])
  assert_eq(result.coefficient_names, ["coefficient_1", "coefficient_2"])
}

///|
test "apeglm: default target is the last design coefficient" {
  let result = @src.apeglm_fit(
    [[5.0, 6.0, 10.0, 12.0]],
    [
      [1.0, 0.0, 0.0],
      [1.0, 1.0, 0.0],
      [1.0, 0.0, 1.0],
      [1.0, 1.0, 1.0],
    ],
    [0.2],
  ) catch {
    _ => abort("default coefficient fit should succeed")
  }
  assert_eq(result.coefficient, 2)
}

///|
test "apeglm: nuisance coefficients are marked no-shrink" {
  let result = apeglm_test_result()
  assert_eq(result.prior_control.no_shrink, [0, 1])
  assert_eq(result.prior_control.prior_mean, 0.0)
  assert_eq(result.prior_control.prior_df, 1.0)
}

///|
test "apeglm: example genes converge" {
  let result = apeglm_test_result()
  for gene in result.genes {
    assert_true(gene.converged)
    assert_true(gene.iterations > 0)
    assert_true(gene.iterations <= result.config.max_iterations)
  }
}

///|
test "apeglm: strong effects retain their directions" {
  let result = apeglm_test_result()
  let up = apeglm_test_gene(result, "strong_up")
  let down = apeglm_test_gene(result, "strong_down")
  assert_true(up.map[result.coefficient] > 1.0)
  assert_true(down.map[result.coefficient] < -1.0)
}

///|
test "apeglm: strong effects are only mildly shrunk" {
  let result = apeglm_test_result()
  for name in ["strong_up", "strong_down"] {
    let gene = apeglm_test_gene(result, name)
    let mle = gene.mle[result.coefficient].abs()
    let map = gene.map[result.coefficient].abs()
    assert_true(map <= mle + 1.0e-8)
    assert_true(map >= 0.8 * mle)
  }
}

///|
test "apeglm: low-count noisy effect is strongly shrunk" {
  let result = apeglm_test_result()
  let gene = apeglm_test_gene(result, "low_noisy")
  let mle = gene.mle[result.coefficient].abs()
  let map = gene.map[result.coefficient].abs()
  assert_true(mle > 1.0)
  assert_true(map < 0.4 * mle)
}

///|
test "apeglm: Cauchy tail preserves large effects" {
  let result = apeglm_test_result()
  let strong = apeglm_test_gene(result, "strong_up")
  let noisy = apeglm_test_gene(result, "low_noisy")
  let strong_ratio = strong.map[result.coefficient].abs() /
    strong.mle[result.coefficient].abs()
  let noisy_ratio = noisy.map[result.coefficient].abs() /
    noisy.mle[result.coefficient].abs()
  assert_true(strong_ratio > noisy_ratio)
  assert_true(strong_ratio > 0.9)
}

///|
test "apeglm: stable gene remains near zero" {
  let result = apeglm_test_result()
  let stable = apeglm_test_gene(result, "stable")
  assert_true(stable.map[result.coefficient].abs() < 0.2)
  assert_true(stable.fsr > 0.3)
}

///|
test "apeglm: abundant moderate signal remains detectable" {
  let result = apeglm_test_result()
  let gene = apeglm_test_gene(result, "abundant_up")
  assert_true(gene.map[result.coefficient] > 0.5)
  assert_true(gene.fsr < 1.0e-3)
}

///|
test "apeglm: posterior standard deviations are finite and positive" {
  let result = apeglm_test_result()
  for gene in result.genes {
    for value in gene.posterior_sd {
      assert_true(value > 0.0)
      assert_true(value < 1.0e300)
      assert_true(value == value)
    }
  }
}

///|
test "apeglm: posterior intervals contain MAP estimates" {
  let result = apeglm_test_result()
  for gene in result.genes {
    let estimate = gene.map[result.coefficient]
    assert_true(gene.interval_lower < estimate)
    assert_true(gene.interval_upper > estimate)
  }
}

///|
test "apeglm: local posterior probabilities are bounded" {
  let result = apeglm_test_result()
  for gene in result.genes {
    assert_true(gene.fsr >= 0.0 && gene.fsr <= 0.5)
    assert_true(gene.s_value >= 0.0 && gene.s_value <= 0.5)
    assert_true(
      gene.threshold_probability >= 0.0 &&
      gene.threshold_probability <= 1.0,
    )
  }
}

///|
test "apeglm: strong signals have lower FSR than null genes" {
  let result = apeglm_test_result()
  let up = apeglm_test_gene(result, "strong_up")
  let down = apeglm_test_gene(result, "strong_down")
  let null = apeglm_test_gene(result, "null_low")
  assert_true(up.fsr < null.fsr)
  assert_true(down.fsr < null.fsr)
}

///|
test "apeglm: threshold probability distinguishes large effects" {
  let result = apeglm_test_result()
  let up = apeglm_test_gene(result, "strong_up")
  let stable = apeglm_test_gene(result, "stable")
  assert_true(up.threshold_probability < 0.01)
  assert_true(stable.threshold_probability > 0.9)
}

///|
test "apeglm: s-values use cumulative mean in FSR order" {
  let values = @src.apeglm_s_values([0.3, 0.1, 0.2])
  apeglm_test_close(values[1], 0.1, 1.0e-12)
  apeglm_test_close(values[2], 0.15, 1.0e-12)
  apeglm_test_close(values[0], 0.2, 1.0e-12)
}

///|
test "apeglm: s-values support empty and tied inputs" {
  assert_eq(@src.apeglm_s_values([]), [])
  assert_eq(@src.apeglm_s_values([0.2, 0.2]), [0.2, 0.2])
}

///|
test "apeglm: empirical Bayes prior variance is positive" {
  let variance = @src.apeglm_prior_variance(
    [1.5, -1.2, 0.1, 0.0],
    [0.2, 0.3, 0.5, 0.6],
  ) catch {
    _ => abort("prior variance estimation should succeed")
  }
  assert_true(variance >= 1.0e-6)
  assert_true(variance <= 400.0)
}

///|
test "apeglm: empirical prior variance responds to effect spread" {
  let narrow = @src.apeglm_prior_variance(
    [0.1, -0.1, 0.0],
    [0.5, 0.5, 0.5],
  ) catch {
    _ => abort("narrow prior estimate should succeed")
  }
  let broad = @src.apeglm_prior_variance(
    [3.0, -3.0, 2.5],
    [0.5, 0.5, 0.5],
  ) catch {
    _ => abort("broad prior estimate should succeed")
  }
  assert_true(broad > narrow)
}

///|
test "apeglm: prior variance respects custom bounds" {
  let variance = @src.apeglm_prior_variance(
    [100.0, -100.0],
    [0.1, 0.1],
    minimum=0.01,
    maximum=5.0,
  ) catch {
    _ => abort("bounded prior estimate should succeed")
  }
  apeglm_test_close(variance, 5.0, 1.0e-10)
}

///|
test "apeglm: prior variance rejects malformed inputs" {
  let length = try {
    ignore(@src.apeglm_prior_variance([1.0], [0.1, 0.2]))
    false
  } catch {
    ApeglmError(_) => true
  }
  let standard_error = try {
    ignore(@src.apeglm_prior_variance([1.0], [0.0]))
    false
  } catch {
    ApeglmError(_) => true
  }
  let bounds = try {
    ignore(@src.apeglm_prior_variance([1.0], [0.1], minimum=2.0, maximum=1.0))
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(length)
  assert_true(standard_error)
  assert_true(bounds)
}

///|
test "apeglm: adaptive fit reports learned prior" {
  let result = apeglm_test_result()
  assert_true(result.prior_control.adaptive)
  assert_true(result.prior_control.prior_variance > 0.0)
  assert_true(result.prior_control.prior_scale >= 0.001)
  assert_true(result.prior_control.prior_scale <= 1.0)
}

///|
test "apeglm: fixed prior preserves configured scale" {
  let (counts, design, dispersions, genes, coefficients) =
    @src.apeglm_example_data()
  let config = @src.ApeglmConfig::create(
    coefficient=2,
    prior_scale=0.4,
    adaptive_prior=false,
  ) catch {
    _ => abort("fixed prior configuration should succeed")
  }
  let result = @src.apeglm_fit(
    counts,
    design,
    dispersions,
    config~,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("fixed prior fit should succeed")
  }
  assert_false(result.prior_control.adaptive)
  apeglm_test_close(result.prior_control.prior_scale, 0.4, 1.0e-12)
  apeglm_test_close(result.prior_control.prior_variance, 0.16, 1.0e-12)
}

///|
test "apeglm: adaptive multiplier obeys scale cap" {
  let (counts, design, dispersions, _, _) = @src.apeglm_example_data()
  let config = @src.ApeglmConfig::create(
    coefficient=2,
    multiplier=100.0,
  ) catch {
    _ => abort("multiplier configuration should succeed")
  }
  let result = @src.apeglm_fit(counts, design, dispersions, config~) catch {
    _ => abort("multiplier fit should succeed")
  }
  apeglm_test_close(result.prior_control.prior_scale, 1.0, 1.0e-12)
}

///|
test "apeglm: NB log-likelihood favors group-specific means" {
  let counts = [10.0, 10.0, 40.0, 40.0]
  let design = [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]]
  let null = @src.apeglm_log_likelihood_nb(
    counts,
    design,
    [@math.ln(25.0), 0.0],
    0.1,
  ) catch {
    _ => abort("null likelihood should succeed")
  }
  let fitted = @src.apeglm_log_likelihood_nb(
    counts,
    design,
    [@math.ln(10.0), @math.ln(4.0)],
    0.1,
  ) catch {
    _ => abort("fitted likelihood should succeed")
  }
  assert_true(fitted > null)
}

///|
test "apeglm: zero likelihood weight excludes an observation" {
  let weighted = @src.apeglm_log_likelihood_nb(
    [5.0, 100.0],
    [[1.0], [1.0]],
    [@math.ln(5.0)],
    0.2,
    weights=[1.0, 0.0],
  ) catch {
    _ => abort("weighted likelihood should succeed")
  }
  let reference = @src.apeglm_log_likelihood_nb(
    [5.0],
    [[1.0]],
    [@math.ln(5.0)],
    0.2,
  ) catch {
    _ => abort("reference likelihood should succeed")
  }
  apeglm_test_close(weighted, reference, 1.0e-12)
}

///|
test "apeglm: likelihood offsets shift the linear predictor" {
  let with_offset = @src.apeglm_log_likelihood_nb(
    [2.0, 3.0],
    [[1.0], [1.0]],
    [0.0],
    0.2,
    offsets=[@math.ln(2.0), @math.ln(2.0)],
  ) catch {
    _ => abort("offset likelihood should succeed")
  }
  let with_beta = @src.apeglm_log_likelihood_nb(
    [2.0, 3.0],
    [[1.0], [1.0]],
    [@math.ln(2.0)],
    0.2,
  ) catch {
    _ => abort("shifted beta likelihood should succeed")
  }
  apeglm_test_close(with_offset, with_beta, 1.0e-12)
}

///|
test "apeglm: likelihood rejects incompatible dimensions" {
  let rows = try {
    ignore(
      @src.apeglm_log_likelihood_nb(
        [1.0, 2.0],
        [[1.0]],
        [0.0],
        0.1,
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let columns = try {
    ignore(
      @src.apeglm_log_likelihood_nb(
        [1.0],
        [[1.0, 0.0]],
        [0.0],
        0.1,
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(rows)
  assert_true(columns)
}

///|
test "apeglm: likelihood rejects invalid dispersion and auxiliaries" {
  let dispersion = try {
    ignore(@src.apeglm_log_likelihood_nb([1.0], [[1.0]], [0.0], 0.0))
    false
  } catch {
    ApeglmError(_) => true
  }
  let offsets = try {
    ignore(
      @src.apeglm_log_likelihood_nb(
        [1.0, 2.0],
        [[1.0], [1.0]],
        [0.0],
        0.1,
        offsets=[0.0],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let weights = try {
    ignore(
      @src.apeglm_log_likelihood_nb(
        [1.0, 2.0],
        [[1.0], [1.0]],
        [0.0],
        0.1,
        weights=[1.0],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(dispersion)
  assert_true(offsets)
  assert_true(weights)
}

///|
test "apeglm: explicit zero offsets match omitted offsets" {
  let (counts, design, dispersions, genes, coefficients) =
    @src.apeglm_example_data()
  let config = apeglm_test_config()
  let baseline = @src.apeglm_fit(
    counts,
    design,
    dispersions,
    config~,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("baseline fit should succeed")
  }
  let shifted = @src.apeglm_fit(
    counts,
    design,
    dispersions,
    config~,
    offsets=apeglm_test_zero_matrix(8, 8),
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("zero offset fit should succeed")
  }
  apeglm_test_close(
    baseline.genes[0].map[2],
    shifted.genes[0].map[2],
    1.0e-12,
  )
}

///|
test "apeglm: zero observation weights suppress an outlier" {
  let config = @src.ApeglmConfig::create(
    coefficient=1,
    adaptive_prior=false,
    prior_scale=1.0,
  ) catch {
    _ => abort("weight test configuration should succeed")
  }
  let result = @src.apeglm_fit(
    [[10.0, 10.0, 20.0, 1000.0]],
    [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]],
    [0.1],
    config~,
    weights=[[1.0, 1.0, 1.0, 0.0]],
  ) catch {
    _ => abort("weighted fit should succeed")
  }
  assert_true(result.genes[0].map[1] > 0.0)
  assert_true(result.genes[0].map[1] < 2.0)
}

///|
test "apeglm: custom target coefficient is honored" {
  let (counts, design, dispersions, genes, coefficients) =
    @src.apeglm_example_data()
  let config = @src.ApeglmConfig::create(coefficient=1) catch {
    _ => abort("batch target configuration should succeed")
  }
  let result = @src.apeglm_fit(
    counts,
    design,
    dispersions,
    config~,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("batch target fit should succeed")
  }
  assert_eq(result.coefficient, 1)
  assert_eq(result.prior_control.no_shrink, [0, 2])
  assert_eq(result.summary().coefficient_name, "batch")
}

///|
test "apeglm: MAP and SD matrix accessors preserve dimensions" {
  let result = apeglm_test_result()
  let maps = result.map_matrix()
  let standard_errors = result.sd_matrix()
  assert_eq(maps.length(), 8)
  assert_eq(maps[0].length(), 3)
  assert_eq(standard_errors.length(), 8)
  assert_eq(standard_errors[0].length(), 3)
  apeglm_test_close(maps[0][2], result.genes[0].map[2], 1.0e-12)
}

///|
test "apeglm: matrix accessors return independent arrays" {
  let result = apeglm_test_result()
  let maps = result.map_matrix()
  let original = result.genes[0].map[0]
  maps[0][0] = original + 100.0
  apeglm_test_close(result.genes[0].map[0], original, 1.0e-12)
}

///|
test "apeglm: index lookup handles bounds" {
  let result = apeglm_test_result()
  match result.gene(0) {
    Some(gene) => assert_eq(gene.gene_name, "strong_up")
    None => abort("first gene should exist")
  }
  assert_true(result.gene(-1) is None)
  assert_true(result.gene(8) is None)
}

///|
test "apeglm: name lookup finds and rejects genes" {
  let result = apeglm_test_result()
  match result.find_gene("stable") {
    Some(gene) => assert_eq(gene.gene_name, "stable")
    None => abort("stable gene should exist")
  }
  assert_true(result.find_gene("missing") is None)
}

///|
test "apeglm: ranking prioritizes low directional error" {
  let result = apeglm_test_result()
  let ranked = result.ranked()
  assert_eq(ranked.length(), 8)
  assert_true(ranked[0].s_value <= ranked[1].s_value)
  assert_true(
    ranked[0].gene_name == "strong_up" ||
    ranked[0].gene_name == "strong_down",
  )
}

///|
test "apeglm: selection applies s-value and effect filters" {
  let result = apeglm_test_result()
  let broad = result.select(maximum_s_value=0.1)
  let strict = result.select(
    maximum_s_value=0.1,
    minimum_absolute_effect=1.0,
  )
  assert_true(broad.length() >= strict.length())
  for gene in strict {
    assert_true(gene.s_value <= 0.1)
    assert_true(gene.map[result.coefficient].abs() >= 1.0)
  }
}

///|
test "apeglm: threshold selection uses FSOS probabilities" {
  let result = apeglm_test_result()
  let selected = result.select(
    maximum_s_value=0.01,
    use_threshold_probability=true,
  )
  assert_true(selected.length() >= 2)
  for gene in selected {
    assert_true(gene.threshold_probability <= 0.01)
  }
}

///|
test "apeglm: summary reports fit diagnostics" {
  let result = apeglm_test_result()
  let summary = result.summary(maximum_s_value=0.1)
  assert_eq(summary.gene_count, 8)
  assert_eq(summary.converged_count, 8)
  assert_eq(summary.coefficient_name, "condition")
  assert_eq(summary.selected_count, result.select(maximum_s_value=0.1).length())
  assert_true(summary.prior_scale > 0.0)
  assert_true(summary.median_absolute_effect >= 0.0)
}

///|
test "apeglm: effect accessor converts natural log to log2" {
  let result = apeglm_test_result()
  let gene = result.genes[0]
  apeglm_test_close(
    gene.effect(2, log2_scale=true),
    gene.effect(2) / @math.ln(2.0),
    1.0e-12,
  )
}

///|
test "apeglm: standard error accessor converts to log2" {
  let result = apeglm_test_result()
  let gene = result.genes[0]
  apeglm_test_close(
    gene.standard_error(2, log2_scale=true),
    gene.standard_error(2) / @math.ln(2.0),
    1.0e-12,
  )
}

///|
test "apeglm: TSV export labels scales and includes genes" {
  let result = apeglm_test_result()
  let log2_table = result.to_tsv(log2_scale=true)
  let natural_table = result.to_tsv(log2_scale=false)
  assert_true(log2_table.contains("gene\tbaseMean\tlog2MLE\tlog2MAP"))
  assert_true(natural_table.contains("gene\tbaseMean\tlnMLE\tlnMAP"))
  assert_true(log2_table.contains("strong_up"))
  assert_true(log2_table.contains("low_noisy"))
}

///|
test "apeglm: repeated fitting is deterministic" {
  let first = apeglm_test_result()
  let second = apeglm_test_result()
  for gene in 0..<first.genes.length() {
    for coefficient in 0..<first.coefficient_names.length() {
      apeglm_test_close(
        first.genes[gene].map[coefficient],
        second.genes[gene].map[coefficient],
        1.0e-12,
      )
    }
  }
}

///|
test "apeglm: deterministic start counts preserve signal modes" {
  let (counts, design, dispersions, genes, coefficients) =
    @src.apeglm_example_data()
  let one = @src.ApeglmConfig::create(coefficient=2, random_starts=1) catch {
    _ => abort("single start configuration should succeed")
  }
  let many = @src.ApeglmConfig::create(coefficient=2, random_starts=12) catch {
    _ => abort("many start configuration should succeed")
  }
  let first = @src.apeglm_fit(
    counts,
    design,
    dispersions,
    config=one,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("single start fit should succeed")
  }
  let second = @src.apeglm_fit(
    counts,
    design,
    dispersions,
    config=many,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("many start fit should succeed")
  }
  assert_true(first.genes[0].map[2] > 0.0)
  assert_true(second.genes[0].map[2] > 0.0)
  assert_true(first.genes[1].map[2] < 0.0)
  assert_true(second.genes[1].map[2] < 0.0)
}

///|
test "apeglm: DESeq2 wrapper uses counts names and offsets" {
  let dds = @src.DESeqDataSet::new(
    [
      [48, 52, 45, 54, 198, 220, 205, 230],
      [160, 148, 171, 155, 39, 44, 35, 41],
    ],
    ["up", "down"],
    ["s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8"],
    [
      [1.0, 0.0, 0.0],
      [1.0, 1.0, 0.0],
      [1.0, 0.0, 0.0],
      [1.0, 1.0, 0.0],
      [1.0, 0.0, 1.0],
      [1.0, 1.0, 1.0],
      [1.0, 0.0, 1.0],
      [1.0, 1.0, 1.0],
    ],
  )
  let config = apeglm_test_config()
  let result = @src.apeglm_from_deseq2(
    dds,
    config~,
    coefficient_names=["intercept", "batch", "condition"],
  ) catch {
    _ => abort("DESeq2 integration should succeed")
  }
  assert_eq(result.gene_names, ["up", "down"])
  assert_true(result.genes[0].map[2] > 0.0)
  assert_true(result.genes[1].map[2] < 0.0)
}

///|
test "apeglm: SummarizedExperiment output adds posterior assays" {
  let (_, design, dispersions, genes, coefficients) =
    @src.apeglm_example_data()
  let config = apeglm_test_config()
  let output = @src.apeglm_summarized_experiment(
    apeglm_test_experiment(),
    design,
    dispersions,
    config~,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("SummarizedExperiment integration should succeed")
  }
  for name in [
    "apeglm_map",
    "apeglm_sd",
    "apeglm_fsr",
    "apeglm_svalue",
    "apeglm_fsos",
  ] {
    match @src.se_assay(output.experiment, name) {
      Some(assay) => {
        assert_eq(assay.length(), 8)
        assert_eq(assay[0].length(), 1)
      }
      None => abort("missing apeglm output assay: " + name)
    }
  }
}

///|
test "apeglm: SummarizedExperiment integration is immutable" {
  let (_, design, dispersions, genes, coefficients) =
    @src.apeglm_example_data()
  let experiment = apeglm_test_experiment()
  let config = apeglm_test_config()
  let output = @src.apeglm_summarized_experiment(
    experiment,
    design,
    dispersions,
    config~,
    gene_names=genes,
    coefficient_names=coefficients,
  ) catch {
    _ => abort("SummarizedExperiment integration should succeed")
  }
  assert_true(@src.se_assay(experiment, "apeglm_map") is None)
  assert_true(experiment.metadata.get("apeglm_scale") is None)
  assert_eq(output.experiment.metadata["apeglm_scale"], "natural-log")
  assert_eq(output.experiment.metadata["apeglm_coefficient"], "condition")
  assert_eq(output.experiment.metadata["source"], "apeglm_test")
}

///|
test "apeglm: SummarizedExperiment rejects missing assay" {
  let (_, design, dispersions, _, _) = @src.apeglm_example_data()
  let raised = try {
    ignore(
      @src.apeglm_summarized_experiment(
        apeglm_test_experiment(),
        design,
        dispersions,
        assay_name="missing",
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(raised)
}

///|
test "apeglm: fit rejects empty and single-sample counts" {
  let empty = try {
    ignore(@src.apeglm_fit([], [], []))
    false
  } catch {
    ApeglmError(_) => true
  }
  let single = try {
    ignore(@src.apeglm_fit([[1.0]], [[1.0, 0.0]], [0.1]))
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(empty)
  assert_true(single)
}

///|
test "apeglm: fit rejects ragged and fractional counts" {
  let ragged = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0], [3.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1, 0.1],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let fractional = try {
    ignore(
      @src.apeglm_fit(
        [[1.5, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(ragged)
  assert_true(fractional)
}

///|
test "apeglm: fit rejects negative and non-finite counts" {
  let negative = try {
    ignore(
      @src.apeglm_fit(
        [[-1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let overflow = try {
    ignore(
      @src.apeglm_fit(
        [[1.0e301, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(negative)
  assert_true(overflow)
}

///|
test "apeglm: fit rejects malformed design matrices" {
  let rows = try {
    ignore(@src.apeglm_fit([[1.0, 2.0]], [[1.0, 0.0]], [0.1]))
    false
  } catch {
    ApeglmError(_) => true
  }
  let columns = try {
    ignore(@src.apeglm_fit([[1.0, 2.0]], [[1.0], [1.0]], [0.1]))
    false
  } catch {
    ApeglmError(_) => true
  }
  let ragged = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0]],
        [0.1],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(rows)
  assert_true(columns)
  assert_true(ragged)
}

///|
test "apeglm: fit rejects intercept and out-of-range targets" {
  let intercept_config = @src.ApeglmConfig::create(coefficient=0) catch {
    _ => abort("intercept index is syntactically valid")
  }
  let range_config = @src.ApeglmConfig::create(coefficient=2) catch {
    _ => abort("range index is syntactically valid")
  }
  let intercept = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
        config=intercept_config,
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let range = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
        config=range_config,
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(intercept)
  assert_true(range)
}

///|
test "apeglm: fit rejects malformed dispersions" {
  let length = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0], [3.0, 4.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let value = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.0],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(length)
  assert_true(value)
}

///|
test "apeglm: fit validates offsets and observation weights" {
  let offsets = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
        offsets=[[0.0]],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let weight_shape = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
        weights=[[1.0]],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let all_zero = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
        weights=[[0.0, 0.0]],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(offsets)
  assert_true(weight_shape)
  assert_true(all_zero)
}

///|
test "apeglm: fit validates identifier lengths and uniqueness" {
  let gene_length = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
        gene_names=["a", "b"],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let duplicate = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0], [2.0, 3.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1, 0.1],
        gene_names=["same", "same"],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let coefficient_length = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
        coefficient_names=["intercept"],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  let empty = try {
    ignore(
      @src.apeglm_fit(
        [[1.0, 2.0]],
        [[1.0, 0.0], [1.0, 1.0]],
        [0.1],
        gene_names=[""],
      ),
    )
    false
  } catch {
    ApeglmError(_) => true
  }
  assert_true(gene_length)
  assert_true(duplicate)
  assert_true(coefficient_length)
  assert_true(empty)
}
