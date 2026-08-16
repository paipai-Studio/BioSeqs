// Tests for the Bioconductor DirichletMultinomial-inspired module.

///|
fn dmn_test_close(
  actual : Double,
  expected : Double,
  tolerance : Double,
) -> Unit {
  if (actual - expected).abs() > tolerance {
    abort(
      "DirichletMultinomial value " +
      actual.to_string() +
      " differs from " +
      expected.to_string(),
    )
  }
}

///|
fn dmn_test_row_sum(values : Array[Double]) -> Double {
  let mut total = 0.0
  for value in values {
    total = total + value
  }
  total
}

///|
fn dmn_test_counts() -> Array[Array[Int]] {
  [
    [40, 3, 2],
    [35, 4, 1],
    [42, 2, 3],
    [38, 5, 2],
    [45, 3, 1],
    [36, 2, 4],
    [2, 40, 3],
    [4, 35, 2],
    [3, 43, 1],
    [5, 38, 2],
    [1, 44, 3],
    [3, 36, 4],
  ]
}

///|
fn dmn_test_groups() -> Array[String] {
  ["A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B"]
}

///|
fn dmn_test_config(components : Int) -> @src.DmnConfig {
  @src.DmnConfig::create(
    components~,
    max_iterations=40,
    optimizer_iterations=100,
    soft_kmeans_iterations=50,
    tolerance=1.0e-5,
    optimizer_tolerance=1.0e-5,
    seed=2026,
  ) catch {
    DirichletMultinomialError(message) =>
      abort("valid DirichletMultinomial configuration failed: " + message)
  }
}

///|
fn dmn_test_fit(components : Int) -> @src.DmnFit {
  @src.dirichlet_multinomial_fit(
    dmn_test_counts(),
    taxon_names=["taxon_a", "taxon_b", "taxon_c"],
    config=dmn_test_config(components),
  ) catch {
    DirichletMultinomialError(message) =>
      abort("DirichletMultinomial test fit failed: " + message)
  }
}

///|
test "dirichlet_multinomial: default configuration matches upstream controls" {
  let config = @src.DmnConfig::default()
  assert_eq(config.components, 1)
  assert_eq(config.max_iterations, 100)
  assert_eq(config.optimizer_iterations, 250)
  assert_eq(config.soft_kmeans_iterations, 200)
  assert_eq(config.tolerance, 1.0e-6)
  assert_eq(config.optimizer_tolerance, 1.0e-5)
  assert_eq(config.soft_beta, 50.0)
  assert_eq(config.seed, 42)
  assert_eq(config.prior_shape, 0.1)
  assert_eq(config.prior_rate, 0.1)
  assert_eq(config.min_alpha, 1.0e-8)
}

///|
test "dirichlet_multinomial: custom configuration preserves controls" {
  let config = @src.DmnConfig::create(
    components=3,
    max_iterations=9,
    optimizer_iterations=11,
    soft_kmeans_iterations=13,
    tolerance=0.01,
    optimizer_tolerance=0.02,
    soft_beta=25.0,
    seed=99,
    prior_shape=0.2,
    prior_rate=0.3,
    min_alpha=1.0e-7,
  ) catch {
    _ => abort("custom DirichletMultinomial configuration should be valid")
  }
  assert_eq(config.components, 3)
  assert_eq(config.max_iterations, 9)
  assert_eq(config.optimizer_iterations, 11)
  assert_eq(config.soft_kmeans_iterations, 13)
  assert_eq(config.soft_beta, 25.0)
  assert_eq(config.seed, 99)
  assert_eq(config.prior_shape, 0.2)
  assert_eq(config.prior_rate, 0.3)
}

///|
test "dirichlet_multinomial: with_components preserves optimization controls" {
  let original = dmn_test_config(1)
  let changed = original.with_components(2, seed=77) catch {
    _ => abort("with_components should accept positive component counts")
  }
  assert_eq(changed.components, 2)
  assert_eq(changed.seed, 77)
  assert_eq(changed.max_iterations, original.max_iterations)
  assert_eq(changed.optimizer_iterations, original.optimizer_iterations)
  assert_eq(changed.soft_beta, original.soft_beta)
}

///|
test "dirichlet_multinomial: configuration rejects nonpositive components" {
  let raised = try {
    ignore(@src.DmnConfig::create(components=0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(raised)
}

///|
test "dirichlet_multinomial: configuration rejects invalid iteration limits" {
  let em = try {
    ignore(@src.DmnConfig::create(max_iterations=0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let optimizer = try {
    ignore(@src.DmnConfig::create(optimizer_iterations=0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let kmeans = try {
    ignore(@src.DmnConfig::create(soft_kmeans_iterations=0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(em)
  assert_true(optimizer)
  assert_true(kmeans)
}

///|
test "dirichlet_multinomial: configuration rejects invalid tolerances" {
  let em = try {
    ignore(@src.DmnConfig::create(tolerance=0.0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let optimizer = try {
    ignore(@src.DmnConfig::create(optimizer_tolerance=-0.1))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(em)
  assert_true(optimizer)
}

///|
test "dirichlet_multinomial: configuration rejects invalid soft beta and seed" {
  let beta = try {
    ignore(@src.DmnConfig::create(soft_beta=0.0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let seed = try {
    ignore(@src.DmnConfig::create(seed=0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(beta)
  assert_true(seed)
}

///|
test "dirichlet_multinomial: configuration rejects invalid Gamma priors" {
  let shape = try {
    ignore(@src.DmnConfig::create(prior_shape=0.0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let rate = try {
    ignore(@src.DmnConfig::create(prior_rate=-0.1))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(shape)
  assert_true(rate)
}

///|
test "dirichlet_multinomial: configuration rejects invalid alpha floor" {
  let zero = try {
    ignore(@src.DmnConfig::create(min_alpha=0.0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let one = try {
    ignore(@src.DmnConfig::create(min_alpha=1.0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(zero)
  assert_true(one)
}

///|
test "dirichlet_multinomial: known beta-binomial probability is exact" {
  let value = @src.dirichlet_multinomial_log_pmf([2, 1], [1.0, 1.0]) catch {
    _ => abort("valid Dirichlet-multinomial PMF should succeed")
  }
  dmn_test_close(value, -1.3862943611198906, 1.0e-10)
}

///|
test "dirichlet_multinomial: PMF normalizes over two-count compositions" {
  let mut total = 0.0
  for left in 0..<=2 {
    let log_probability = @src.dirichlet_multinomial_log_pmf([left, 2 - left], [
      1.0, 1.0,
    ]) catch {
      _ => abort("valid Dirichlet-multinomial PMF should succeed")
    }
    total = total + @math.exp(log_probability)
  }
  dmn_test_close(total, 1.0, 1.0e-10)
}

///|
test "dirichlet_multinomial: PMF is symmetric under paired permutation" {
  let left = @src.dirichlet_multinomial_log_pmf([4, 1], [2.0, 3.0]) catch {
    _ => abort("valid Dirichlet-multinomial PMF should succeed")
  }
  let right = @src.dirichlet_multinomial_log_pmf([1, 4], [3.0, 2.0]) catch {
    _ => abort("valid Dirichlet-multinomial PMF should succeed")
  }
  dmn_test_close(left, right, 1.0e-12)
}

///|
test "dirichlet_multinomial: all-zero count vector has unit probability" {
  let value = @src.dirichlet_multinomial_log_pmf([0, 0, 0], [0.5, 1.5, 2.5]) catch {
    _ => abort("zero-count Dirichlet-multinomial PMF should succeed")
  }
  dmn_test_close(value, 0.0, 1.0e-12)
}

///|
test "dirichlet_multinomial: PMF rejects dimension mismatch" {
  let raised = try {
    ignore(@src.dirichlet_multinomial_log_pmf([1, 2], [1.0]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(raised)
}

///|
test "dirichlet_multinomial: PMF rejects negative counts" {
  let raised = try {
    ignore(@src.dirichlet_multinomial_log_pmf([1, -1], [1.0, 1.0]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(raised)
}

///|
test "dirichlet_multinomial: PMF rejects nonpositive alpha" {
  let zero = try {
    ignore(@src.dirichlet_multinomial_log_pmf([1, 1], [1.0, 0.0]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let huge = try {
    ignore(@src.dirichlet_multinomial_log_pmf([1, 1], [1.0, 1.0e301]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(zero)
  assert_true(huge)
}

///|
test "dirichlet_multinomial: mean follows alpha proportions" {
  let mean = @src.dirichlet_multinomial_mean(10, [2.0, 3.0]) catch {
    _ => abort("valid Dirichlet-multinomial mean should succeed")
  }
  assert_eq(mean.length(), 2)
  dmn_test_close(mean[0], 4.0, 1.0e-12)
  dmn_test_close(mean[1], 6.0, 1.0e-12)
}

///|
test "dirichlet_multinomial: zero-total mean is zero" {
  let mean = @src.dirichlet_multinomial_mean(0, [2.0, 3.0]) catch {
    _ => abort("zero-total Dirichlet-multinomial mean should succeed")
  }
  assert_eq(mean, [0.0, 0.0])
}

///|
test "dirichlet_multinomial: mean validates total and alpha" {
  let total = try {
    ignore(@src.dirichlet_multinomial_mean(-1, [1.0]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let alpha = try {
    ignore(@src.dirichlet_multinomial_mean(2, [1.0, -1.0]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(total)
  assert_true(alpha)
}

///|
test "dirichlet_multinomial: covariance includes overdispersion inflation" {
  let covariance = @src.dirichlet_multinomial_covariance(10, [2.0, 3.0]) catch {
    _ => abort("valid Dirichlet-multinomial covariance should succeed")
  }
  dmn_test_close(covariance[0][0], 6.0, 1.0e-12)
  dmn_test_close(covariance[1][1], 6.0, 1.0e-12)
  dmn_test_close(covariance[0][1], -6.0, 1.0e-12)
  dmn_test_close(covariance[1][0], -6.0, 1.0e-12)
}

///|
test "dirichlet_multinomial: zero-total covariance is zero" {
  let covariance = @src.dirichlet_multinomial_covariance(0, [1.0, 2.0, 3.0]) catch {
    _ => abort("zero-total covariance should succeed")
  }
  for row in covariance {
    assert_eq(row, [0.0, 0.0, 0.0])
  }
}

///|
test "dirichlet_multinomial: single-component fit exposes dimensions and names" {
  let fit = dmn_test_fit(1)
  assert_eq(fit.component_count(), 1)
  assert_eq(fit.sample_count(), 12)
  assert_eq(fit.taxon_count(), 3)
  assert_eq(fit.sample_names[0], "sample1")
  assert_eq(fit.sample_names[11], "sample12")
  assert_eq(fit.taxon_names, ["taxon_a", "taxon_b", "taxon_c"])
}

///|
test "dirichlet_multinomial: single-component posterior and weight are one" {
  let fit = dmn_test_fit(1)
  assert_eq(fit.mixture_weights(), [1.0])
  for row in fit.mixture() {
    dmn_test_close(row[0], 1.0, 1.0e-12)
  }
  assert_eq(fit.assignments(), Array::make(12, 0))
}

///|
test "dirichlet_multinomial: estimates and intervals remain positive" {
  let fit = dmn_test_fit(1)
  for taxon in 0..<fit.taxon_count() {
    assert_true(fit.estimates[0][taxon] > 0.0)
    assert_true(fit.lower[0][taxon] > 0.0)
    assert_true(fit.upper[0][taxon] > 0.0)
    assert_true(fit.lower[0][taxon] <= fit.upper[0][taxon])
  }
}

///|
test "dirichlet_multinomial: component proportions and fitted scale normalize" {
  let fit = dmn_test_fit(1)
  let proportions = fit.component_proportions()
  dmn_test_close(dmn_test_row_sum(proportions[0]), 1.0, 1.0e-12)
  let scaled = fit.fitted(scale=true)
  assert_eq(scaled.length(), 3)
  assert_eq(scaled[0].length(), 1)
  dmn_test_close(scaled[0][0] + scaled[1][0] + scaled[2][0], 1.0, 1.0e-12)
}

///|
test "dirichlet_multinomial: fitted matrix transposes component alpha" {
  let fit = dmn_test_fit(1)
  let fitted = fit.fitted()
  for taxon in 0..<fit.taxon_count() {
    dmn_test_close(fitted[taxon][0], fit.estimates[0][taxon], 1.0e-12)
  }
}

///|
test "dirichlet_multinomial: two-component weights are sorted descending" {
  let fit = dmn_test_fit(2)
  let weights = fit.mixture_weights()
  assert_eq(weights.length(), 2)
  assert_true(weights[0] >= weights[1])
  dmn_test_close(weights[0] + weights[1], 1.0, 1.0e-10)
}

///|
test "dirichlet_multinomial: two-component posterior rows normalize" {
  let fit = dmn_test_fit(2)
  for row in fit.mixture() {
    assert_eq(row.length(), 2)
    dmn_test_close(dmn_test_row_sum(row), 1.0, 1.0e-10)
    assert_true(row[0] >= 0.0 && row[0] <= 1.0)
    assert_true(row[1] >= 0.0 && row[1] <= 1.0)
  }
}

///|
test "dirichlet_multinomial: synthetic clusters receive distinct assignments" {
  let assignments = dmn_test_fit(2).assignments()
  let first = assignments[0]
  let second = assignments[6]
  assert_true(first != second)
  for sample in 0..<6 {
    assert_eq(assignments[sample], first)
  }
  for sample in 6..<12 {
    assert_eq(assignments[sample], second)
  }
}

///|
test "dirichlet_multinomial: synthetic components recover dominant taxa" {
  let fit = dmn_test_fit(2)
  let assignments = fit.assignments()
  let proportions = fit.component_proportions()
  let first_component = assignments[0]
  let second_component = assignments[6]
  assert_true(proportions[first_component][0] > proportions[first_component][1])
  assert_true(
    proportions[second_component][1] > proportions[second_component][0],
  )
}

///|
test "dirichlet_multinomial: prediction normalizes and classifies novel samples" {
  let fit = dmn_test_fit(2)
  let posterior = fit.predict([[50, 1, 1], [1, 50, 1]]) catch {
    _ => abort("valid DirichletMultinomial prediction should succeed")
  }
  for row in posterior {
    dmn_test_close(dmn_test_row_sum(row), 1.0, 1.0e-10)
  }
  let assignments = fit.predict_assignments([[50, 1, 1], [1, 50, 1]]) catch {
    _ => abort("valid DirichletMultinomial assignment should succeed")
  }
  assert_true(assignments[0] != assignments[1])
}

///|
test "dirichlet_multinomial: component evidence has sample by component shape" {
  let fit = dmn_test_fit(2)
  let evidence = fit.negative_log_evidence([[50, 1, 1], [1, 50, 1]]) catch {
    _ => abort("valid DirichletMultinomial evidence should succeed")
  }
  assert_eq(evidence.length(), 2)
  assert_eq(evidence[0].length(), 2)
  for row in evidence {
    for value in row {
      assert_true(value.abs() < 1.0e300)
    }
  }
}

///|
test "dirichlet_multinomial: concentrations equal alpha row sums" {
  let fit = dmn_test_fit(2)
  let concentrations = fit.concentrations()
  for component in 0..<fit.component_count() {
    dmn_test_close(
      concentrations[component],
      dmn_test_row_sum(fit.estimates[component]),
      1.0e-10,
    )
  }
}

///|
test "dirichlet_multinomial: goodness-of-fit criteria use upstream parameter count" {
  let fit = dmn_test_fit(2)
  let parameters = 7.0
  dmn_test_close(
    fit.goodness.aic,
    fit.goodness.negative_log_evidence + parameters,
    1.0e-10,
  )
  dmn_test_close(
    fit.goodness.bic,
    fit.goodness.negative_log_evidence + 0.5 * @math.ln(12.0) * parameters,
    1.0e-10,
  )
  assert_true(fit.goodness.laplace.abs() < 1.0e300)
  assert_true(fit.goodness.log_determinant.abs() < 1.0e300)
}

///|
test "dirichlet_multinomial: fit summary reports dimensions and criteria" {
  let summary = dmn_test_fit(2).summary()
  assert_true(summary.contains("k=2"))
  assert_true(summary.contains("samples=12"))
  assert_true(summary.contains("taxa=3"))
  assert_true(summary.contains("Laplace="))
  assert_true(summary.contains("BIC="))
  assert_true(summary.contains("AIC="))
}

///|
test "dirichlet_multinomial: prediction rejects taxon mismatch" {
  let fit = dmn_test_fit(1)
  let raised = try {
    ignore(fit.predict([[1, 2]]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(raised)
}

///|
test "dirichlet_multinomial: model selection fits every requested K" {
  let selection = @src.dirichlet_multinomial_select(
    dmn_test_counts(),
    2,
    taxon_names=["taxon_a", "taxon_b", "taxon_c"],
    config=dmn_test_config(1),
  ) catch {
    _ => abort("DirichletMultinomial model selection should succeed")
  }
  assert_eq(selection.fits.length(), 2)
  assert_eq(selection.scores.length(), 2)
  assert_eq(selection.fits[0].component_count(), 1)
  assert_eq(selection.fits[1].component_count(), 2)
  assert_eq(selection.criterion, @src.DmnLaplace)
}

///|
test "dirichlet_multinomial: model selection best index minimizes score" {
  let selection = @src.dirichlet_multinomial_select(
    dmn_test_counts(),
    2,
    criterion=@src.DmnBic,
    config=dmn_test_config(1),
  ) catch {
    _ => abort("DirichletMultinomial BIC selection should succeed")
  }
  assert_eq(selection.criterion, @src.DmnBic)
  assert_eq(
    selection.best_fit().component_count(),
    selection.best_component_count(),
  )
  assert_true(selection.scores[selection.best_index] <= selection.scores[0])
  assert_true(selection.scores[selection.best_index] <= selection.scores[1])
}

///|
test "dirichlet_multinomial: model selection validates component range" {
  let zero = try {
    ignore(@src.dirichlet_multinomial_select(dmn_test_counts(), 0))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let excessive = try {
    ignore(@src.dirichlet_multinomial_select(dmn_test_counts(), 13))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(zero)
  assert_true(excessive)
}

///|
test "dirichlet_multinomial: fit rejects empty and ragged matrices" {
  let empty = try {
    ignore(@src.dirichlet_multinomial_fit([]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let ragged = try {
    ignore(@src.dirichlet_multinomial_fit([[1, 2], [3]]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(empty)
  assert_true(ragged)
}

///|
test "dirichlet_multinomial: fit rejects negative and zero-library samples" {
  let negative = try {
    ignore(@src.dirichlet_multinomial_fit([[1, -1], [2, 3]]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let zero = try {
    ignore(@src.dirichlet_multinomial_fit([[1, 2], [0, 0]]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(negative)
  assert_true(zero)
}

///|
test "dirichlet_multinomial: fit rejects too many components" {
  let raised = try {
    ignore(
      @src.dirichlet_multinomial_fit(
        [[1, 2], [2, 1]],
        config=dmn_test_config(3),
      ),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(raised)
}

///|
test "dirichlet_multinomial: fit validates sample names" {
  let count = try {
    ignore(
      @src.dirichlet_multinomial_fit([[1, 2], [2, 1]], sample_names=["one"]),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let duplicate = try {
    ignore(
      @src.dirichlet_multinomial_fit([[1, 2], [2, 1]], sample_names=[
        "same", "same",
      ]),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(count)
  assert_true(duplicate)
}

///|
test "dirichlet_multinomial: fit validates taxon names" {
  let empty = try {
    ignore(
      @src.dirichlet_multinomial_fit([[1, 2], [2, 1]], taxon_names=[
        "taxon", " ",
      ]),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let duplicate = try {
    ignore(
      @src.dirichlet_multinomial_fit([[1, 2], [2, 1]], taxon_names=[
        "same", "same",
      ]),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(empty)
  assert_true(duplicate)
}

///|
test "dirichlet_multinomial: returned mixture arrays are defensive copies" {
  let fit = dmn_test_fit(1)
  let weights = fit.mixture_weights()
  weights[0] = 0.0
  let mixture = fit.mixture()
  mixture[0][0] = 0.0
  assert_eq(fit.mixture_weights(), [1.0])
  assert_eq(fit.mixture()[0], [1.0])
}

///|
test "dirichlet_multinomial: group fit preserves labels and empirical priors" {
  let classifier = @src.dirichlet_multinomial_group_fit(
    dmn_test_counts(),
    dmn_test_groups(),
    components_by_group=[1, 1],
    taxon_names=["taxon_a", "taxon_b", "taxon_c"],
    config=dmn_test_config(1),
  ) catch {
    _ => abort("DirichletMultinomial group fit should succeed")
  }
  assert_eq(classifier.group_names, ["A", "B"])
  assert_eq(classifier.models.length(), 2)
  assert_eq(classifier.priors, [0.5, 0.5])
  assert_eq(classifier.taxon_names, ["taxon_a", "taxon_b", "taxon_c"])
}

///|
test "dirichlet_multinomial: group prediction normalizes and separates classes" {
  let classifier = @src.dirichlet_multinomial_group_fit(
    dmn_test_counts(),
    dmn_test_groups(),
    components_by_group=[1, 1],
    config=dmn_test_config(1),
  ) catch {
    _ => abort("DirichletMultinomial group fit should succeed")
  }
  let probabilities = classifier.predict([[50, 1, 1], [1, 50, 1]]) catch {
    _ => abort("DirichletMultinomial group prediction should succeed")
  }
  for row in probabilities {
    dmn_test_close(dmn_test_row_sum(row), 1.0, 1.0e-10)
  }
  let predictions = classifier.predict_assignments([[50, 1, 1], [1, 50, 1]]) catch {
    _ => abort("DirichletMultinomial group assignment should succeed")
  }
  assert_eq(predictions, ["A", "B"])
}

///|
test "dirichlet_multinomial: group fitting validates labels and component counts" {
  let labels = try {
    ignore(
      @src.dirichlet_multinomial_group_fit([[1, 2], [2, 1]], ["A"], components_by_group=[
        1, 1,
      ]),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let one_group = try {
    ignore(
      @src.dirichlet_multinomial_group_fit([[1, 2], [2, 1]], ["A", "A"], components_by_group=[
        1,
      ]),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let components = try {
    ignore(
      @src.dirichlet_multinomial_group_fit([[1, 2], [2, 1]], ["A", "B"], components_by_group=[
        1,
      ]),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(labels)
  assert_true(one_group)
  assert_true(components)
}

///|
test "dirichlet_multinomial: stratified cross-validation returns held-out results" {
  let result = @src.dirichlet_multinomial_cross_validate(
    dmn_test_counts(),
    dmn_test_groups(),
    3,
    components_by_group=[1, 1],
    taxon_names=["taxon_a", "taxon_b", "taxon_c"],
    config=dmn_test_config(1),
  ) catch {
    _ => abort("DirichletMultinomial cross-validation should succeed")
  }
  assert_eq(result.probabilities.length(), 12)
  assert_eq(result.predictions.length(), 12)
  assert_eq(result.truth, dmn_test_groups())
  assert_eq(result.group_names, ["A", "B"])
  assert_eq(result.fold_ids, [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2])
  for row in result.probabilities {
    dmn_test_close(dmn_test_row_sum(row), 1.0, 1.0e-10)
  }
  assert_true(result.accuracy >= 0.9)
  assert_eq(result.correct.to_double() / 12.0, result.accuracy)
}

///|
test "dirichlet_multinomial: cross-validation validates folds and group sizes" {
  let folds = try {
    ignore(
      @src.dirichlet_multinomial_cross_validate(
        dmn_test_counts(),
        dmn_test_groups(),
        1,
      ),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let sizes = try {
    ignore(
      @src.dirichlet_multinomial_cross_validate(
        [[5, 1], [4, 1], [1, 5], [1, 4]],
        ["A", "A", "B", "B"],
        3,
      ),
    )
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(folds)
  assert_true(sizes)
}

///|
test "dirichlet_multinomial: perfect ROC has unit AUC" {
  let roc = @src.dirichlet_multinomial_roc([true, true, false, false], [
    0.9, 0.8, 0.2, 0.1,
  ]) catch {
    _ => abort("valid DirichletMultinomial ROC should succeed")
  }
  assert_eq(roc.auc, 1.0)
  assert_eq(roc.best_threshold, 0.8)
  assert_eq(roc.sensitivity[0], 0.0)
  assert_eq(roc.specificity[0], 1.0)
  assert_eq(roc.sensitivity[roc.sensitivity.length() - 1], 1.0)
  assert_eq(roc.specificity[roc.specificity.length() - 1], 0.0)
}

///|
test "dirichlet_multinomial: tied ROC scores contribute half concordance" {
  let roc = @src.dirichlet_multinomial_roc([true, false], [0.5, 0.5]) catch {
    _ => abort("tied DirichletMultinomial ROC should succeed")
  }
  assert_eq(roc.auc, 0.5)
  assert_eq(roc.thresholds.length(), 2)
}

///|
test "dirichlet_multinomial: ROC rejects invalid labels and scores" {
  let lengths = try {
    ignore(@src.dirichlet_multinomial_roc([true], [0.5, 0.2]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let classes = try {
    ignore(@src.dirichlet_multinomial_roc([true, true], [0.5, 0.2]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  let score = try {
    ignore(@src.dirichlet_multinomial_roc([true, false], [1.0e301, 0.2]))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(lengths)
  assert_true(classes)
  assert_true(score)
}

///|
test "dirichlet_multinomial: SummarizedExperiment entry transposes feature by sample assay" {
  let assays : Map[String, Array[Array[Double]]] = Map([
    (
      "counts",
      [
        [40.0, 35.0, 42.0, 2.0, 4.0, 3.0],
        [3.0, 4.0, 2.0, 40.0, 35.0, 43.0],
        [2.0, 1.0, 3.0, 3.0, 2.0, 1.0],
      ],
    ),
  ])
  let experiment = @src.summarized_experiment(assays, [], [], Map([]))
  let fit = @src.dirichlet_multinomial_fit_se(
    experiment,
    "counts",
    sample_names=["s1", "s2", "s3", "s4", "s5", "s6"],
    taxon_names=["taxon_a", "taxon_b", "taxon_c"],
    config=dmn_test_config(2),
  ) catch {
    _ => abort("DirichletMultinomial SummarizedExperiment fit should succeed")
  }
  assert_eq(fit.sample_count(), 6)
  assert_eq(fit.taxon_count(), 3)
  assert_eq(fit.sample_names[0], "s1")
  assert_eq(fit.taxon_names[2], "taxon_c")
  assert_true(fit.assignments()[0] != fit.assignments()[3])
}

///|
test "dirichlet_multinomial: SummarizedExperiment entry rejects missing assay" {
  let experiment = @src.summarized_experiment(Map([]), [], [], Map([]))
  let raised = try {
    ignore(@src.dirichlet_multinomial_fit_se(experiment, "counts"))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(raised)
}

///|
test "dirichlet_multinomial: SummarizedExperiment entry rejects noninteger assay" {
  let assays : Map[String, Array[Array[Double]]] = Map([
    ("counts", [[1.5, 2.0], [3.0, 4.0]]),
  ])
  let experiment = @src.summarized_experiment(assays, [], [], Map([]))
  let raised = try {
    ignore(@src.dirichlet_multinomial_fit_se(experiment, "counts"))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(raised)
}

///|
test "dirichlet_multinomial: SummarizedExperiment entry rejects negative assay" {
  let assays : Map[String, Array[Array[Double]]] = Map([
    ("counts", [[1.0, -2.0], [3.0, 4.0]]),
  ])
  let experiment = @src.summarized_experiment(assays, [], [], Map([]))
  let raised = try {
    ignore(@src.dirichlet_multinomial_fit_se(experiment, "counts"))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(raised)
}

///|
test "dirichlet_multinomial: SummarizedExperiment entry rejects ragged assay" {
  let assays : Map[String, Array[Array[Double]]] = Map([
    ("counts", [[1.0, 2.0], [3.0]]),
  ])
  let experiment = @src.summarized_experiment(assays, [], [], Map([]))
  let raised = try {
    ignore(@src.dirichlet_multinomial_fit_se(experiment, "counts"))
    false
  } catch {
    DirichletMultinomialError(_) => true
  }
  assert_true(raised)
}
