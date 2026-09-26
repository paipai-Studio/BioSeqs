///|
/// Tests for linear mixed models (lmm.mbt).

fn lmm_approx(a : Double, b : Double, tol : Double) -> Bool {
  let d = a - b
  if d < 0.0 {
    -d < tol
  } else {
    d < tol
  }
}

///|
/// n x 1 intercept-only design matrix.
fn lmm_intercept(n : Int) -> Array[Array[Double]] {
  let x : Array[Array[Double]] = []
  for _ in 0..<n {
    x.push([1.0])
  }
  x
}

///|
/// n x n identity matrix.
fn lmm_identity(n : Int) -> Array[Array[Double]] {
  let m : Array[Array[Double]] = []
  for i in 0..<n {
    let row = Array::make(n, 0.0)
    row[i] = 1.0
    m.push(row)
  }
  m
}

// ============================================================
// GRM computation
// ============================================================

test "lmm_compute_grm_symmetric" {
  let g = [
    [0.0, 1.0, 2.0],
    [0.0, 1.0, 2.0],
    [2.0, 1.0, 0.0],
  ]
  let k = @src.lmm_compute_grm(g)
  assert_eq(k.length(), 3)
  // Symmetry.
  assert_true(lmm_approx(k[0][1], k[1][0], 1.0e-9))
  assert_true(lmm_approx(k[0][2], k[2][0], 1.0e-9))
  assert_true(lmm_approx(k[1][2], k[2][1], 1.0e-9))
  // Diagonal (self-relatedness) positive.
  assert_true(k[0][0] > 0.0)
  assert_true(k[1][1] > 0.0)
  assert_true(k[2][2] > 0.0)
}

test "lmm_compute_grm_relatedness_direction" {
  let g = [
    [0.0, 1.0, 2.0],
    [0.0, 1.0, 2.0], // identical to sample 0
    [2.0, 1.0, 0.0], // opposite of sample 0
  ]
  let k = @src.lmm_compute_grm(g)
  // Identical samples -> high positive relatedness.
  assert_true(k[0][1] > 0.5)
  assert_true(lmm_approx(k[0][1], k[0][0], 1.0e-9))
  // Opposite samples -> negative relatedness.
  assert_true(k[0][2] < 0.0)
}

// ============================================================
// LMM fit
// ============================================================

test "lmm_fit_intercept_only_mean" {
  // With identity kinship, the GLS fixed effect reduces to OLS: beta = mean(y).
  let n = 4
  let y = [10.0, 20.0, 30.0, 40.0]
  let x = lmm_intercept(n)
  let k = lmm_identity(n)
  let model = @src.lmm_fit(y, x, k)
  assert_eq(model.beta().length(), 1)
  // mean of y = 25
  assert_true(lmm_approx(model.beta()[0], 25.0, 1.0e-6))
}

test "lmm_fit_variance_components_finite" {
  // Structured kinship from a small GRM; components must be finite and non-negative.
  let g = [
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 2.0],
    [2.0, 2.0],
    [2.0, 1.0],
    [1.0, 0.0],
  ]
  let k = @src.lmm_compute_grm(g)
  let n = 6
  let y = [1.0, 1.2, 2.1, 3.0, 2.2, 1.1]
  let x = lmm_intercept(n)
  let model = @src.lmm_fit(y, x, k)
  assert_true(!(model.sigma2_g().is_nan()))
  assert_true(!(model.sigma2_e().is_nan()))
  assert_true(model.sigma2_g() >= 0.0)
  assert_true(model.sigma2_e() >= 0.0)
  let h = model.heritability()
  assert_true(h >= 0.0 && h <= 1.0)
  assert_true(!(model.log_likelihood().is_nan()))
}

// ============================================================
// Association testing (EMMAX)
// ============================================================

test "lmm_associate_finds_causal_signal" {
  // Strong linear effect of the causal SNP; a null SNP is uncorrelated.
  let n = 14
  let causal = [
    0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0,
  ]
  let null_snp = [
    0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0,
  ]
  let noise = [
    1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0,
  ]
  let y : Array[Double] = []
  for i in 0..<n {
    y.push(1.0 + 2.0 * causal[i] + 0.5 * noise[i])
  }

  let x = lmm_intercept(n)
  let k = lmm_identity(n)
  let model = @src.lmm_fit(y, x, k)

  let results = @src.lmm_associate(model, [causal, null_snp])
  assert_eq(results.length(), 2)

  // Causal SNP: strongly significant.
  assert_true(results[0].p_value() < 0.01)
  // Null SNP: not significant.
  assert_true(results[1].p_value() > 0.05)
  // Causal effect estimate close to the true slope (2.0).
  assert_true(lmm_approx(results[0].beta(), 2.0, 0.5))
}

test "lmm_associate_more_samples_reduces_pvalue" {
  // A larger sample size with the same effect should give a smaller p-value.
  let n = 60
  let causal : Array[Double] = []
  let null_snp : Array[Double] = []
  let y : Array[Double] = []
  for i in 0..<n {
    let c = (i % 3).to_double() // 0,1,2 cycling
    let nu = ((i / 30) % 2).to_double() // block structure
    causal.push(c)
    null_snp.push(nu)
    y.push(1.0 + 1.5 * c + 0.3 * ((i % 5) - 2).to_double())
  }
  let x = lmm_intercept(n)
  let k = lmm_identity(n)
  let model = @src.lmm_fit(y, x, k)
  let results = @src.lmm_associate(model, [causal, null_snp])
  assert_true(results[0].p_value() < 0.001)
  assert_true(results[1].p_value() > 0.05)
}

test "lmm_associate_structured_kinship" {
  // With a structured GRM, the association machinery should still run and
  // identify the causal SNP.
  let g = [
    [0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0, 1.0],
    [1.0, 1.0, 1.0, 1.0],
    [2.0, 2.0, 2.0, 2.0],
    [2.0, 2.0, 2.0, 2.0],
  ]
  let k = @src.lmm_compute_grm(g)
  let n = 6
  let causal = [0.0, 1.0, 0.0, 1.0, 2.0, 2.0]
  let null_snp = [0.0, 0.0, 1.0, 1.0, 1.0, 0.0]
  let y : Array[Double] = []
  for i in 0..<n {
    y.push(2.0 * causal[i] + 0.2 * g[i][0])
  }
  let x = lmm_intercept(n)
  let model = @src.lmm_fit(y, x, k)
  let results = @src.lmm_associate(model, [causal, null_snp])
  assert_eq(results.length(), 2)
  assert_true(!(results[0].p_value().is_nan()))
  assert_true(!(results[1].p_value().is_nan()))
  assert_true(results[0].p_value() < results[1].p_value())
}

// ============================================================
// 解析解对拍：K=I 时 LMM beta 退化为 OLS
// ============================================================

test "lmm_identity_kinship_matches_ols" {
  // 同一组数据 y=0.8+2.3x；K=I 时 GLS 固定效应应退化为 OLS 解。
  let xs = [0.0, 1.0, 2.0, 3.0]
  let ys = [1.0, 3.0, 5.0, 8.0]
  let x : Array[Array[Double]] = []
  for v in xs {
    x.push([1.0, v])
  }
  let k = lmm_identity(4)
  let model = @src.lmm_fit(ys, x, k)
  let b = model.beta()
  assert_true(lmm_approx(b[0], 0.8, 1.0e-3))
  assert_true(lmm_approx(b[1], 2.3, 1.0e-3))
}
