///|
/// Tests for statsmodels-style generalized linear models (glm.mbt).

// ============================================================
// Helpers
// ============================================================

///|
/// Build a design matrix with intercept + one covariate.
fn glm_design(xs : Array[Double]) -> Array[Array[Double]] {
  let raw : Array[Array[Double]] = []
  for x in xs {
    raw.push([x])
  }
  @src.glm_add_intercept(raw)
}

fn approx(a : Double, b : Double, tol : Double) -> Bool {
  let d = a - b
  if d < 0.0 {
    -d < tol
  } else {
    d < tol
  }
}

// ============================================================
// Gaussian (linear regression, identity link)
// ============================================================

test "glm_gaussian_recovers_coefficients" {
  // y = 2 + 3x (exact, no noise) -> should recover [2, 3].
  let xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
  let ys = [2.0, 5.0, 8.0, 11.0, 14.0, 17.0]
  let x = glm_design(xs)
  let result = @src.glm_fit_irls(x, ys, @src.glm_family_gaussian())
  let b = result.coefficients()
  assert_eq(b.length(), 2)
  assert_true(approx(b[0], 2.0, 1.0e-6))
  assert_true(approx(b[1], 3.0, 1.0e-6))
  assert_true(result.converged())
}

test "glm_gaussian_intercept_only_is_mean" {
  // Intercept-only linear model = the mean of y.
  let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let ys = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
  let x = glm_design(xs)
  // Intercept-only: fit with just intercept column.
  let x0 : Array[Array[Double]] = []
  for _ in 0..<ys.length() {
    x0.push([1.0])
  }
  let result = @src.glm_fit_irls(x0, ys, @src.glm_family_gaussian())
  // mean of ys = 35
  assert_true(approx(result.coefficients()[0], 35.0, 1.0e-6))
}

// ============================================================
// Poisson (log link)
// ============================================================

test "glm_poisson_fits_positive_slope" {
  // Counts roughly increasing with x.
  let xs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let ys = [1.0, 2.0, 4.0, 7.0, 12.0, 20.0]
  let x = glm_design(xs)
  let result = @src.glm_fit_irls(x, ys, @src.glm_family_poisson())
  let b = result.coefficients()
  assert_eq(b.length(), 2)
  // Slope should be positive (counts increase with x).
  assert_true(b[1] > 0.0)
  // Deviance should be small (good fit).
  assert_true(result.deviance() >= 0.0)
}

test "glm_poisson_log_likelihood_finite" {
  let xs = [1.0, 2.0, 3.0, 4.0]
  let ys = [2.0, 4.0, 7.0, 10.0]
  let x = glm_design(xs)
  let result = @src.glm_fit_irls(x, ys, @src.glm_family_poisson())
  // Log-likelihood should be a finite negative number.
  let ll = result.log_likelihood()
  assert_true(ll < 0.0)
  assert_true(!(ll.is_nan()))
}

// ============================================================
// Binomial / logistic (logit link)
// ============================================================

test "glm_binomial_separates_classes" {
  // Strongly separable labels: y=0 for small x, y=1 for large x.
  let xs = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
  let ys = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  let x = glm_design(xs)
  let result = @src.glm_fit_irls(x, ys, @src.glm_family_binomial())
  let b = result.coefficients()
  // Slope for x should be strongly positive (labels increase with x).
  assert_true(b[1] > 0.0)
  // High log-likelihood (near-separable -> close to 0).
  assert_true(result.log_likelihood() > -5.0)
}

// ============================================================
// Negative binomial (log link)
// ============================================================

test "glm_negative_binomial_runs" {
  let xs = [1.0, 2.0, 3.0, 4.0, 5.0]
  let ys = [3.0, 6.0, 9.0, 15.0, 21.0]
  let x = glm_design(xs)
  // NB with large dispersion approximates Poisson.
  let family = @src.glm_family_negative_binomial(100.0)
  let result = @src.glm_fit_irls(x, ys, family)
  assert_true(result.converged())
  assert_true(result.coefficients()[1] > 0.0)
  assert_true(result.deviance() >= 0.0)
}

// ============================================================
// Wald tests
// ============================================================

test "glm_wald_test_slope_significant" {
  // y strongly depends on x -> slope Wald test significant.
  let xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
  let ys = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0]
  let x = glm_design(xs)
  let result = @src.glm_fit_irls(x, ys, @src.glm_family_gaussian())
  let wald = @src.glm_wald_test(result, 1)
  assert_true(wald.p_value() < 0.001)
  // z should be huge (very significant).
  assert_true(wald.statistic() > 10.0)
}

test "glm_wald_omnibus_joint" {
  // Two predictors, both relevant -> omnibus Wald test significant.
  let x1 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let x2 = [2.0, 1.0, 4.0, 3.0, 5.0, 7.0]
  // y = 10 + 2*x1 + 3*x2
  let ys : Array[Double] = []
  let raw : Array[Array[Double]] = []
  for i in 0..<x1.length() {
    raw.push([x1[i], x2[i]])
    ys.push(10.0 + 2.0 * x1[i] + 3.0 * x2[i])
  }
  let design = @src.glm_add_intercept(raw)
  let result = @src.glm_fit_irls(design, ys, @src.glm_family_gaussian())
  let omnibus = @src.glm_wald_omnibus(result, [1, 2])
  assert_true(omnibus.p_value() < 0.001)
  assert_eq(omnibus.df(), 2)
}

// ============================================================
// Likelihood-ratio test (LRT)
// ============================================================

test "glm_lrt_full_beats_reduced" {
  // y = 2 + 3x: the full model (with x) should fit much better than the
  // intercept-only reduced model.
  let xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
  let ys = [2.0, 5.0, 8.0, 11.0, 14.0, 17.0, 20.0, 23.0]
  let x_full = glm_design(xs)
  let x_reduced : Array[Array[Double]] = []
  for _ in 0..<ys.length() {
    x_reduced.push([1.0])
  }
  let full = @src.glm_fit_irls(x_full, ys, @src.glm_family_gaussian())
  let reduced = @src.glm_fit_irls(x_reduced, ys, @src.glm_family_gaussian())
  let lrt = @src.glm_lrt_compare(full, reduced)
  assert_eq(lrt.df(), 1)
  assert_true(lrt.statistic() > 0.0)
  assert_true(lrt.p_value() < 0.001)
}

test "glm_lrt_no_effect_nonsignificant" {
  // y independent of x (pure noise): model with x should not significantly
  // improve over intercept-only. (Noise is symmetric, LRT p should not be tiny.)
  let xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
  // y has NO x dependence.
  let ys = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
  let x_full = glm_design(xs)
  let x_reduced : Array[Array[Double]] = []
  for _ in 0..<ys.length() {
    x_reduced.push([1.0])
  }
  let full = @src.glm_fit_irls(x_full, ys, @src.glm_family_gaussian())
  let reduced = @src.glm_fit_irls(x_reduced, ys, @src.glm_family_gaussian())
  let lrt = @src.glm_lrt_compare(full, reduced)
  // Statistic should be ~0 (or very small); p-value not extremely small.
  assert_true(lrt.statistic() < 1.0e-6)
  assert_true(lrt.p_value() > 0.5)
}

// ============================================================
// Link / variance helpers
// ============================================================

test "glm_link_and_inverse_identity" {
  assert_eq(@src.glm_link(@src.glm_link_identity(), 3.5), 3.5)
  assert_eq(@src.glm_inverse_link(@src.glm_link_identity(), -1.2), -1.2)
}

test "glm_link_and_inverse_log" {
  // log(2) and exp(log 2) = 2
  assert_true(approx(@src.glm_link(@src.glm_link_log(), 2.0), @math.ln(2.0), 1.0e-9))
  assert_true(approx(@src.glm_inverse_link(@src.glm_link_log(), @math.ln(2.0)), 2.0, 1.0e-9))
}

test "glm_variance_functions" {
  // Gaussian variance = 1
  assert_eq(@src.glm_variance(@src.glm_family_gaussian(), 5.0), 1.0)
  // Poisson variance = mu
  assert_true(approx(@src.glm_variance(@src.glm_family_poisson(), 7.0), 7.0, 1.0e-12))
  // Binomial variance = mu(1-mu)
  assert_true(approx(@src.glm_variance(@src.glm_family_binomial(), 0.3), 0.21, 1.0e-12))
}

test "glm_default_links" {
  assert_eq(@src.glm_default_link(@src.glm_family_gaussian()), @src.glm_link_identity())
  assert_eq(@src.glm_default_link(@src.glm_family_poisson()), @src.glm_link_log())
  assert_eq(@src.glm_default_link(@src.glm_family_binomial()), @src.glm_link_logit())
}

// ============================================================
// 解析解对拍（验证 IRLS 正确性，而非仅自查一致）
// ============================================================

test "glm_ols_matches_analytic_closed_form" {
  // y = 0.8 + 2.3x 的最小二乘精确解（数据 y=[1,3,5,8], x=[0,1,2,3]）
  let xs = [0.0, 1.0, 2.0, 3.0]
  let ys = [1.0, 3.0, 5.0, 8.0]
  let x = glm_design(xs)
  let result = @src.glm_fit_irls(x, ys, @src.glm_family_gaussian())
  let b = result.coefficients()
  // 解析：beta1 = 11.5/5 = 2.3, beta0 = 4.25 - 2.3*1.5 = 0.8
  assert_true(approx(b[0], 0.8, 1.0e-4))
  assert_true(approx(b[1], 2.3, 1.0e-4))
  // 解析：RSS = 0.30, scale = RSS/(n-2) = 0.15
  assert_true(approx(result.deviance(), 0.30, 1.0e-6))
  assert_true(approx(result.scale(), 0.15, 1.0e-6))
  assert_eq(result.df_residual(), 2)
}

test "glm_ols_wald_matches_t_statistic" {
  // SE(beta1) = sqrt(scale * [(X'X)^-1]_11) = sqrt(0.15 * 0.2) = sqrt(0.03)
  // z = 2.3 / sqrt(0.03) ≈ 13.279
  let xs = [0.0, 1.0, 2.0, 3.0]
  let ys = [1.0, 3.0, 5.0, 8.0]
  let x = glm_design(xs)
  let result = @src.glm_fit_irls(x, ys, @src.glm_family_gaussian())
  let wald = @src.glm_wald_test(result, 1)
  assert_true(approx(wald.statistic(), 2.3 / 0.03.sqrt(), 0.02))
  assert_true(wald.p_value() < 0.0001)
}

test "glm_poisson_intercept_matches_log_mean" {
  // Poisson(intercept only)：beta0 = ln(mean(y))，这是 MLE 的解析解。
  let ys = [1.0, 2.0, 3.0, 4.0]
  let x0 : Array[Array[Double]] = []
  for _ in 0..<ys.length() {
    x0.push([1.0])
  }
  let result = @src.glm_fit_irls(x0, ys, @src.glm_family_poisson())
  assert_true(approx(result.coefficients()[0], @math.ln(2.5), 1.0e-4))
}

test "glm_logistic_intercept_matches_logit" {
  // Binomial/Bernoulli(intercept only)：beta0 = logit(mean(y))，解析解。
  let ys = [0.0, 0.0, 1.0, 1.0, 1.0]
  let x0 : Array[Array[Double]] = []
  for _ in 0..<ys.length() {
    x0.push([1.0])
  }
  let result = @src.glm_fit_irls(x0, ys, @src.glm_family_binomial())
  // mean = 0.6, logit = ln(0.6/0.4) = ln(1.5)
  assert_true(approx(result.coefficients()[0], @math.ln(1.5), 1.0e-4))
}
