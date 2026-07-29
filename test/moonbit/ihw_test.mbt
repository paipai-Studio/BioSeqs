///|
/// Tests for Bio.IHW module - Independent Hypothesis Weighting.

// ---------------------------------------------------------------------------
// IHWResult tests
// ---------------------------------------------------------------------------

///|
test "ihw_result_creation" {
  let pvals = [0.01, 0.04, 0.03]
  let adj = [0.03, 0.05, 0.04]
  let weights = [1.0, 1.0, 1.0]
  let covs = [0.5, 0.3, 0.8]
  let result = @src.IHWResult::new(pvals, adj, weights, covs, 0.05, 5)
  assert_eq(result.p_values().length(), 3)
  assert_eq(result.adjusted_p_values().length(), 3)
  assert_eq(result.weights().length(), 3)
  assert_eq(result.alpha(), 0.05)
  assert_eq(result.n_attempts(), 5)
}

///|
test "ihw_result_accessors" {
  let pvals = [0.01]
  let adj = [0.05]
  let weights = [2.0]
  let covs = [0.5]
  let result = @src.IHWResult::new(pvals, adj, weights, covs, 0.01, 3)
  assert_true(result.p_values()[0] > 0.009 && result.p_values()[0] < 0.011)
  assert_true(result.adjusted_p_values()[0] > 0.049 && result.adjusted_p_values()[0] < 0.051)
  assert_true(result.weights()[0] > 1.9 && result.weights()[0] < 2.1)
  assert_true(result.covariates()[0] > 0.4 && result.covariates()[0] < 0.6)
  assert_true(result.alpha() > 0.009 && result.alpha() < 0.011)
  assert_eq(result.n_attempts(), 3)
}

// ---------------------------------------------------------------------------
// IHWConfig tests
// ---------------------------------------------------------------------------

///|
test "ihw_config_defaults" {
  let config = @src.IHWConfig::new()
  assert_true(config.alpha() > 0.049 && config.alpha() < 0.051)
  assert_eq(config.n_attempts(), 10)
  assert_eq(config.scale_type(), "local")
}

///|
test "ihw_config_custom" {
  let config = @src.IHWConfig::new(alpha=0.01, n_attempts=20, scale_type="global")
  assert_true(config.alpha() > 0.009 && config.alpha() < 0.011)
  assert_eq(config.n_attempts(), 20)
  assert_eq(config.scale_type(), "global")
}

// ---------------------------------------------------------------------------
// bonferroni_ihw tests
// ---------------------------------------------------------------------------

///|
test "bonferroni_ihw_basic" {
  let pvals = [0.01, 0.02, 0.03, 0.04]
  let result = @src.bonferroni_ihw(pvals)
  assert_eq(result.length(), 4)
  assert_true(result[0] > 0.03 && result[0] < 0.05)
  assert_true(result[1] > 0.07 && result[1] < 0.09)
  assert_true(result[2] > 0.11 && result[2] < 0.13)
  assert_true(result[3] > 0.15 && result[3] < 0.17)
}

///|
test "bonferroni_ihw_clamped_at_1" {
  let pvals = [0.3, 0.4, 0.5]
  let result = @src.bonferroni_ihw(pvals)
  assert_eq(result.length(), 3)
  assert_true(result[0] <= 1.0)
  assert_true(result[1] <= 1.0)
  assert_true(result[2] <= 1.0)
}

///|
test "bonferroni_ihw_empty" {
  let pvals : Array[Double] = Array::new()
  let result = @src.bonferroni_ihw(pvals)
  assert_eq(result.length(), 0)
}

///|
test "bonferroni_ihw_single" {
  let pvals = [0.05]
  let result = @src.bonferroni_ihw(pvals)
  assert_eq(result.length(), 1)
  assert_true(result[0] > 0.04 && result[0] < 0.06)
}

///|
test "bonferroni_ihw_all_small" {
  let pvals = [0.001, 0.002, 0.003]
  let result = @src.bonferroni_ihw(pvals)
  assert_eq(result.length(), 3)
  assert_true(result[0] > 0.002 && result[0] < 0.004)
  assert_true(result[1] > 0.005 && result[1] < 0.007)
  assert_true(result[2] > 0.008 && result[2] < 0.01)
}

// ---------------------------------------------------------------------------
// Basic IHW correction tests
// ---------------------------------------------------------------------------

///|
test "ihw_basic_correction" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  assert_eq(result.p_values().length(), 5)
  assert_eq(result.adjusted_p_values().length(), 5)
  assert_eq(result.weights().length(), 5)
  assert_true(result.alpha() > 0.04 && result.alpha() < 0.06)
  assert_true(result.n_attempts() >= 0)
}

///|
test "ihw_adjusted_pvals_non_decreasing" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  let adj = result.adjusted_p_values()
  let mut i = 1
  while i < adj.length() {
    assert_true(adj[i] >= adj[i - 1] - 1.0e-10)
    i = i + 1
  }
}

///|
test "ihw_adjusted_pvals_bounded" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  let adj = result.adjusted_p_values()
  let mut i = 0
  while i < adj.length() {
    assert_true(adj[i] >= 0.0)
    assert_true(adj[i] <= 1.0)
    i = i + 1
  }
}

///|
test "ihw_weights_sum_positive" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  let w = result.weights()
  let mut sum = 0.0
  let mut i = 0
  while i < w.length() {
    sum = sum + w[i]
    i = i + 1
  }
  assert_true(sum > 0.0)
}

// ---------------------------------------------------------------------------
// IHW with covariates tests
// ---------------------------------------------------------------------------

///|
test "ihw_with_covariates_different_weights" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  let w = result.weights()
  let w0 = w[0]
  let w4 = w[4]
  assert_true(w0 > 0.0)
  assert_true(w4 > 0.0)
}

///|
test "ihw_with_extreme_covariates" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.001, 0.002, 0.998, 0.999, 1.0]
  let result = @src.ihw(pvals, covs, 0.05)
  assert_eq(result.weights().length(), 5)
  let w = result.weights()
  assert_true(w[0] > 0.0)
  assert_true(w[4] > 0.0)
}

///|
test "ihw_with_constant_covariates" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.5, 0.5, 0.5, 0.5, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  assert_eq(result.weights().length(), 5)
}

// ---------------------------------------------------------------------------
// ihw_with_config tests
// ---------------------------------------------------------------------------

///|
test "ihw_with_config_global" {
  let pvals = [0.01, 0.04, 0.03, 0.02]
  let covs = [0.5, 0.3, 0.8, 0.2]
  let config = @src.IHWConfig::new(alpha=0.05, n_attempts=5, scale_type="global")
  let result = @src.ihw_with_config(pvals, covs, config)
  assert_eq(result.adjusted_p_values().length(), 4)
  let adj = result.adjusted_p_values()
  assert_true(adj[0] > 0.03 && adj[0] < 0.05)
}

///|
test "ihw_with_config_custom_alpha" {
  let pvals = [0.01, 0.04, 0.03, 0.02]
  let covs = [0.5, 0.3, 0.8, 0.2]
  let config = @src.IHWConfig::new(alpha=0.01, n_attempts=5, scale_type="local")
  let result = @src.ihw_with_config(pvals, covs, config)
  assert_true(result.alpha() > 0.009 && result.alpha() < 0.011)
}

///|
test "ihw_with_config_many_attempts" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let config = @src.IHWConfig::new(alpha=0.05, n_attempts=50, scale_type="local")
  let result = @src.ihw_with_config(pvals, covs, config)
  assert_true(result.n_attempts() <= 50)
}

// ---------------------------------------------------------------------------
// ihw_multiple tests
// ---------------------------------------------------------------------------

///|
test "ihw_multiple_basic" {
  let pvals = [0.001, 0.01, 0.03, 0.04]
  let cov1 = [0.1, 0.2, 0.3, 0.4]
  let cov2 = [0.5, 0.6, 0.7, 0.8]
  let covs = [cov1, cov2]
  let result = @src.ihw_multiple(pvals, covs, 0.05)
  assert_eq(result.adjusted_p_values().length(), 4)
  assert_eq(result.weights().length(), 4)
}

///|
test "ihw_multiple_single_covariate" {
  let pvals = [0.001, 0.01, 0.03]
  let cov1 = [0.1, 0.2, 0.3]
  let covs = [cov1]
  let result = @src.ihw_multiple(pvals, covs, 0.05)
  assert_eq(result.adjusted_p_values().length(), 3)
}

///|
test "ihw_multiple_no_covariates" {
  let pvals = [0.001, 0.01, 0.03]
  let covs : Array[Array[Double]] = Array::new()
  let result = @src.ihw_multiple(pvals, covs, 0.05)
  assert_eq(result.adjusted_p_values().length(), 3)
}

// ---------------------------------------------------------------------------
// storey_ihw tests
// ---------------------------------------------------------------------------

///|
test "storey_ihw_basic" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.storey_ihw(pvals, covs)
  assert_eq(result.length(), 5)
  let mut i = 0
  while i < result.length() {
    assert_true(result[i] >= 0.0)
    assert_true(result[i] <= 1.0)
    i = i + 1
  }
}

///|
test "storey_ihw_empty" {
  let pvals : Array[Double] = Array::new()
  let covs : Array[Double] = Array::new()
  let result = @src.storey_ihw(pvals, covs)
  assert_eq(result.length(), 0)
}

///|
test "storey_ihw_single" {
  let pvals = [0.05]
  let covs = [0.5]
  let result = @src.storey_ihw(pvals, covs)
  assert_eq(result.length(), 1)
  assert_true(result[0] >= 0.0)
  assert_true(result[0] <= 1.0)
}

///|
test "storey_ihw_many" {
  let pvals = [
    0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9,
  ]
  let covs = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
  let result = @src.storey_ihw(pvals, covs)
  assert_eq(result.length(), 12)
  let mut i = 0
  while i < result.length() {
    assert_true(result[i] >= 0.0)
    assert_true(result[i] <= 1.0)
    i = i + 1
  }
}

// ---------------------------------------------------------------------------
// ihw_weights tests
// ---------------------------------------------------------------------------

///|
test "ihw_weights_extraction" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  let w = @src.ihw_weights(result)
  assert_eq(w.length(), 5)
  let mut i = 0
  while i < w.length() {
    assert_true(w[i] > 0.0)
    i = i + 1
  }
}

///|
test "ihw_weights_via_accessor" {
  let pvals = [0.01, 0.04, 0.03]
  let covs = [0.5, 0.3, 0.8]
  let result = @src.ihw(pvals, covs, 0.05)
  let w1 = result.weights()
  let w2 = @src.ihw_weights(result)
  assert_eq(w1.length(), w2.length())
  let mut i = 0
  while i < w1.length() {
    let diff = w1[i] - w2[i]
    let abs_diff = if diff < 0.0 { -diff } else { diff }
    assert_true(abs_diff < 1.0e-15)
    i = i + 1
  }
}

// ---------------------------------------------------------------------------
// Edge case tests
// ---------------------------------------------------------------------------

///|
test "ihw_empty_input" {
  let pvals : Array[Double] = Array::new()
  let covs : Array[Double] = Array::new()
  let result = @src.ihw(pvals, covs, 0.05)
  assert_eq(result.p_values().length(), 0)
  assert_eq(result.adjusted_p_values().length(), 0)
  assert_eq(result.weights().length(), 0)
  assert_eq(result.n_attempts(), 0)
}

///|
test "ihw_single_input" {
  let pvals = [0.01]
  let covs = [0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  assert_eq(result.p_values().length(), 1)
  assert_eq(result.adjusted_p_values().length(), 1)
  assert_eq(result.weights().length(), 1)
}

///|
test "ihw_all_significant" {
  let pvals = [0.0001, 0.0005, 0.001, 0.002, 0.003]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  let adj = result.adjusted_p_values()
  let n_significant = count_less_than(adj, 0.05)
  assert_true(n_significant >= 1)
}

///|
test "ihw_none_significant" {
  let pvals = [0.6, 0.7, 0.8, 0.9, 0.95]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  let adj = result.adjusted_p_values()
  let n_significant = count_less_than(adj, 0.05)
  assert_eq(n_significant, 0)
}

///|
test "ihw_very_small_alpha" {
  let pvals = [0.001, 0.01, 0.03, 0.04, 0.5]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.001)
  assert_true(result.alpha() > 0.0005 && result.alpha() < 0.0015)
}

// ---------------------------------------------------------------------------
// Multiple testing scenario tests
// ---------------------------------------------------------------------------

///|
test "ihw_many_tests" {
  let n = 100
  let pvals : Array[Double] = Array::make(n, 0.0)
  let covs : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    pvals[i] = ((i + 1).to_double() / (n.to_double() + 1.0)) * 0.1
    covs[i] = (i + 1).to_double() / n.to_double()
    i = i + 1
  }
  let result = @src.ihw(pvals, covs, 0.05)
  assert_eq(result.adjusted_p_values().length(), n)
  assert_eq(result.weights().length(), n)
}

///|
test "ihw_all_same_pvalue" {
  let pvals = [0.05, 0.05, 0.05, 0.05, 0.05]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  assert_eq(result.adjusted_p_values().length(), 5)
}

///|
test "ihw_very_small_pvalues" {
  let pvals = [1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2]
  let covs = [0.1, 0.2, 0.3, 0.4, 0.5]
  let result = @src.ihw(pvals, covs, 0.05)
  assert_eq(result.adjusted_p_values().length(), 5)
  let adj = result.adjusted_p_values()
  assert_true(adj[0] <= 1.0)
  assert_true(adj[4] <= 1.0)
}

// ---------------------------------------------------------------------------
// Helper function for tests
// ---------------------------------------------------------------------------

///|
fn count_less_than(arr : Array[Double], threshold : Double) -> Int {
  let mut count = 0
  let mut i = 0
  while i < arr.length() {
    if arr[i] < threshold {
      count = count + 1
    }
    i = i + 1
  }
  count
}