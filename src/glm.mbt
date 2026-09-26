///|
/// Generalized Linear Models (statsmodels-style foundation).
///
/// A general GLM substrate for differential expression, eQTL and GWAS style
/// analyses: link functions, exponential-family variance functions, iterative
/// reweighted least squares (IRLS) fitting, Wald tests and likelihood-ratio
/// tests (LRT).
///
/// Mirrors the core of Python statsmodels' `GLM` (with `GLMResults` /
/// `wald_test` / `compare_lr_test`) and R's `glm`/`anova`.
///
/// Supported families: Gaussian, Poisson, negative binomial (NB2, fixed
/// dispersion), binomial (Bernoulli). Supported links: identity, log, logit.
///
/// Mixed models (the GEMMA/EMMAX LMM) are a separate layer built on this
/// foundation and are intentionally out of scope here.

// ============================================================
// Model primitives
// ============================================================

///|
/// Link function g: eta = g(mu) where mu = mean of the response.
pub enum GlmLink {
  Identity
  Log
  Logit
} derive(Eq, Debug)

///|
/// Exponential-family distribution (with variance function).
pub enum GlmFamily {
  Gaussian
  Poisson
  /// Negative binomial NB2 with fixed dispersion phi (V(mu) = mu + mu^2/phi).
  NegativeBinomial(Double)
  /// Binomial with a single trial (Bernoulli), V(mu) = mu(1-mu).
  Binomial
} derive(Eq, Debug)

///|
/// Constructor functions for link and family (enum variants cannot be
/// constructed directly from external packages in this MoonBit version).
pub fn glm_link_identity() -> GlmLink {
  GlmLink::Identity
}

pub fn glm_link_log() -> GlmLink {
  GlmLink::Log
}

pub fn glm_link_logit() -> GlmLink {
  GlmLink::Logit
}

pub fn glm_family_gaussian() -> GlmFamily {
  GlmFamily::Gaussian
}

pub fn glm_family_poisson() -> GlmFamily {
  GlmFamily::Poisson
}

pub fn glm_family_negative_binomial(phi : Double) -> GlmFamily {
  GlmFamily::NegativeBinomial(phi)
}

pub fn glm_family_binomial() -> GlmFamily {
  GlmFamily::Binomial
}

///|
/// Clamp a mean to a small positive minimum (avoids log(0) / division by 0).
fn glm_pos(x : Double) -> Double {
  if x > 0.0 {
    x
  } else {
    1.0e-12
  }
}

///|
/// Clamp a probability to (1e-6, 1 - 1e-6).
fn glm_clamp01(x : Double) -> Double {
  if x < 1.0e-6 {
    1.0e-6
  } else if x > 0.999999 {
    0.999999
  } else {
    x
  }
}

///|
/// Link function value: eta = g(mu).
pub fn glm_link(link : GlmLink, mu : Double) -> Double {
  match link {
    Identity => mu
    Log => @math.ln(glm_pos(mu))
    Logit => {
      let m = glm_clamp01(mu)
      @math.ln(m / (1.0 - m))
    }
  }
}

///|
/// Inverse link: mu = g^{-1}(eta).
pub fn glm_inverse_link(link : GlmLink, eta : Double) -> Double {
  match link {
    Identity => eta
    Log => @math.exp(if eta > 709.0 { 709.0 } else { eta })
    Logit => 1.0 / (1.0 + @math.exp(if eta > 709.0 { 709.0 } else { -eta }))
  }
}

///|
/// mu_eta: d mu / d eta (used in the IRLS weight and working response).
pub fn glm_mu_eta(link : GlmLink, mu : Double) -> Double {
  match link {
    Identity => 1.0
    Log => mu
    Logit => {
      let m = glm_clamp01(mu)
      m * (1.0 - m)
    }
  }
}

///|
/// Variance function V(mu) of the family.
pub fn glm_variance(family : GlmFamily, mu : Double) -> Double {
  match family {
    Gaussian => 1.0
    Poisson => glm_pos(mu)
    NegativeBinomial(phi) => {
      let ph = if phi > 0.0 { phi } else { 1.0 }
      let m = glm_pos(mu)
      m + m * m / ph
    }
    Binomial => {
      let m = glm_clamp01(mu)
      m * (1.0 - m)
    }
  }
}

///|
/// Default link for a family (statsmodels/R convention).
pub fn glm_default_link(family : GlmFamily) -> GlmLink {
  match family {
    Gaussian => Identity
    Poisson => Log
    NegativeBinomial(_) => Log
    Binomial => Logit
  }
}

// ============================================================
// Small linear-algebra helpers (avoid aliasing bug)
// ============================================================

///|
/// Allocate a rows x cols matrix of Double (distinct inner arrays).
fn glm_dmatrix(rows : Int, cols : Int) -> Array[Array[Double]] {
  let m : Array[Array[Double]] = Array::make(rows, Array::make(0, 0.0))
  for i in 0..<rows {
    m[i] = Array::make(cols, 0.0)
  }
  m
}

///|
/// Xᵀ W X for design matrix X (n×p) and diagonal weights W (n).
fn glm_xtwx(
  x : Array[Array[Double]],
  w : Array[Double]
) -> Array[Array[Double]] {
  let n = x.length()
  let p = if n > 0 { x[0].length() } else { 0 }
  let xtwx = glm_dmatrix(p, p)
  for i in 0..<p {
    for j in 0..<p {
      let mut s = 0.0
      for k in 0..<n {
        s = s + x[k][i] * w[k] * x[k][j]
      }
      xtwx[i][j] = s
    }
  }
  xtwx
}

///|
/// Xᵀ W z for design matrix X (n×p), diagonal weights W (n), vector z (n).
fn glm_xtwz(
  x : Array[Array[Double]],
  w : Array[Double],
  z : Array[Double]
) -> Array[Double] {
  let n = x.length()
  let p = if n > 0 { x[0].length() } else { 0 }
  let xtwz = Array::make(p, 0.0)
  for i in 0..<p {
    let mut s = 0.0
    for k in 0..<n {
      s = s + x[k][i] * w[k] * z[k]
    }
    xtwz[i] = s
  }
  xtwz
}

///|
/// Solve the linear system A x = b via matrix inversion + matrix-vector
/// multiplication.
fn glm_solve(a : Array[Array[Double]], b : Array[Double]) -> Array[Double] {
  let inv = linalg_inverse(a)
  let p = inv.length()
  let x = Array::make(p, 0.0)
  for i in 0..<p {
    let mut s = 0.0
    for j in 0..<p {
      s = s + inv[i][j] * b[j]
    }
    x[i] = s
  }
  x
}

// ============================================================
// Deviance and log-likelihood
// ============================================================

///|
/// Per-observation deviance contribution (unscaled). Summed over i gives the
/// residual deviance.
fn glm_deviance_one(family : GlmFamily, y : Double, mu : Double) -> Double {
  match family {
    Gaussian => (y - mu) * (y - mu)
    Poisson => {
      if y <= 0.0 {
        2.0 * mu
      } else {
        2.0 * (y * @math.ln(y / glm_pos(mu)) - (y - mu))
      }
    }
    NegativeBinomial(phi) => {
      let ph = if phi > 0.0 { phi } else { 1.0 }
      let m = glm_pos(mu)
      if y <= 0.0 {
        2.0 * ph * @math.ln((m + ph) / ph)
      } else {
        2.0 * (
          y * @math.ln(y / m) -
          (y + ph) * @math.ln((y + ph) / (m + ph))
        )
      }
    }
    Binomial => {
      let m = glm_clamp01(mu)
      let ym = glm_pos(y)
      2.0 * (
        ym * @math.ln(ym / m) +
        (1.0 - ym) * @math.ln((1.0 - ym) / (1.0 - m))
      )
    }
  }
}

///|
/// Log-likelihood contribution (omits y-only additive constants that cancel
/// in likelihood-ratio comparisons; Gaussian uses the MLE variance).
fn glm_loglik_one(family : GlmFamily, y : Double, mu : Double) -> Double {
  match family {
    Gaussian => {
      // loglik under MLE sigma^2 is handled at the aggregate level; return
      // the quadratic part here (scaled in glm_total_loglik).
      -(y - mu) * (y - mu)
    }
    Poisson => {
      if y <= 0.0 {
        -mu
      } else {
        y * @math.ln(glm_pos(mu)) -
        mu -
        bsseq_ln_gamma(y + 1.0)
      }
    }
    NegativeBinomial(phi) => {
      let ph = if phi > 0.0 { phi } else { 1.0 }
      let m = glm_pos(mu)
      let yy = glm_pos(y)
      bsseq_ln_gamma(yy + ph) -
      bsseq_ln_gamma(ph) -
      bsseq_ln_gamma(yy + 1.0) +
      yy * @math.ln(m / (m + ph)) +
      ph * @math.ln(ph / (m + ph))
    }
    Binomial => {
      let m = glm_clamp01(mu)
      let ym = glm_pos(y)
      ym * @math.ln(m) + (1.0 - ym) * @math.ln(1.0 - m)
    }
  }
}

// ============================================================
// IRLS fitting
// ============================================================

///|
/// Result of a GLM fit.
pub struct GlmResult {
  coefficients : Array[Double]
  std_errors : Array[Double]
  covariance : Array[Array[Double]]
  log_likelihood : Double
  deviance : Double
  df_residual : Int
  n_obs : Int
  converged : Bool
  iterations : Int
  /// Estimated dispersion (scale) parameter.
  scale : Double
}

///|
/// Accessors.
pub fn GlmResult::coefficients(self : GlmResult) -> Array[Double] {
  self.coefficients
}

pub fn GlmResult::std_errors(self : GlmResult) -> Array[Double] {
  self.std_errors
}

pub fn GlmResult::covariance(self : GlmResult) -> Array[Array[Double]] {
  self.covariance
}

pub fn GlmResult::log_likelihood(self : GlmResult) -> Double {
  self.log_likelihood
}

pub fn GlmResult::deviance(self : GlmResult) -> Double {
  self.deviance
}

pub fn GlmResult::df_residual(self : GlmResult) -> Int {
  self.df_residual
}

pub fn GlmResult::converged(self : GlmResult) -> Bool {
  self.converged
}

pub fn GlmResult::scale(self : GlmResult) -> Double {
  self.scale
}

///|
/// Add an intercept (leading column of 1.0) to a raw design matrix.
pub fn glm_add_intercept(x : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = x.length()
  let p = if n > 0 { x[0].length() } else { 0 }
  let out : Array[Array[Double]] = Array(capacity=0)
  for i in 0..<n {
    let row = Array::make(p + 1, 1.0)
    for j in 0..<p {
      row[j + 1] = x[i][j]
    }
    out.push(row)
  }
  out
}

///|
/// Fit a generalized linear model by iteratively reweighted least squares.
///
/// `x` is the design matrix (n rows × p columns). Include an intercept column
/// explicitly (via `glm_add_intercept`) — statsmodels-style implicit intercept
/// is not added here. `y` is the response (length n). `family` is the
/// distribution with its variance function; `link` is the link function
/// (defaults to the family's canonical link).
///
/// - `max_iter`: IRLS iterations cap (default 25).
/// - `tol`: convergence threshold on the coefficient change (default 1e-8).
pub fn glm_fit_irls(
  x : Array[Array[Double]],
  y : Array[Double],
  family : GlmFamily,
  link~ : GlmLink? = None,
  max_iter~ : Int = 25,
  tol~ : Double = 1.0e-8,
) -> GlmResult {
  let link = match link {
    Some(l) => l
    None => glm_default_link(family)
  }
  let n = y.length()
  let p = if x.length() > 0 { x[0].length() } else { 0 }

  // Initialize mu (mean) from y.
  let mu : Array[Double] = Array::make(n, 0.0)
  match family {
    Gaussian => {
      for i in 0..<n {
        mu[i] = y[i]
      }
    }
    Binomial => {
      for i in 0..<n {
        let ym = if y[i] <= 0.0 { 0.1 } else { if y[i] >= 1.0 { 0.9 } else { y[i] } }
        mu[i] = (ym + 0.5) / 2.0
      }
    }
    _ => {
      for i in 0..<n {
        mu[i] = if y[i] <= 0.0 { 0.1 } else { y[i] }
      }
    }
  }

  let mut beta : Array[Double] = Array::make(p, 0.0)
  let mut converged = false
  let mut iters = 0
  let mut w_final : Array[Double] = Array::make(n, 1.0)

  for it in 1..=max_iter {
    // eta = X beta, mu = g^-1(eta).
    let eta = glm_linear_predictor(x, beta)
    for i in 0..<n {
      mu[i] = glm_inverse_link(link, eta[i])
      // clamp mu for families that require positivity
      match family {
        Poisson | NegativeBinomial(_) => {
          if mu[i] <= 0.0 {
            mu[i] = 1.0e-12
          }
        }
        Binomial => {
          if mu[i] < 0.000001 {
            mu[i] = 0.000001
          } else if mu[i] > 0.999999 {
            mu[i] = 0.999999
          }
        }
        Gaussian => ()
      }
    }

    // Build working response and weights.
    let z : Array[Double] = Array::make(n, 0.0)
    let w : Array[Double] = Array::make(n, 0.0)
    for i in 0..<n {
      let me = glm_mu_eta(link, mu[i])
      z[i] = eta[i] + (y[i] - mu[i]) / me
      w[i] = me * me / glm_variance(family, mu[i])
      if w[i] > 1.0e12 {
        w[i] = 1.0e12
      }
    }

    // Weighted least squares update.
    let xtwx = glm_xtwx(x, w)
    let xtwz = glm_xtwz(x, w, z)
    let beta_new = glm_solve(xtwx, xtwz)

    // Convergence check.
    let mut max_delta = 0.0
    for j in 0..<p {
      let d = beta_new[j] - beta[j]
      let ad = if d < 0.0 { -d } else { d }
      if ad > max_delta {
        max_delta = ad
      }
    }
    beta = beta_new
    iters = it
    w_final = w
    if max_delta < tol {
      converged = true
      break
    }
  }

  // Final fitted mean, covariance, scale, deviance, log-likelihood.
  let eta = glm_linear_predictor(x, beta)
  for i in 0..<n {
    mu[i] = glm_inverse_link(link, eta[i])
    match family {
      Poisson | NegativeBinomial(_) => if mu[i] <= 0.0 {
        mu[i] = 1.0e-12
      }
      Binomial => {
        if mu[i] < 0.000001 {
          mu[i] = 0.000001
        } else if mu[i] > 0.999999 {
          mu[i] = 0.999999
        }
      }
      Gaussian => ()
    }
  }

  let df_residual = n - p
  // Dispersion (scale) parameter.
  let mut scale = 1.0
  let deviance = glm_total_deviance(family, y, mu)
  match family {
    Gaussian => {
      if df_residual > 0 {
        scale = deviance / df_residual.to_double()
      } else {
        scale = 1.0
      }
    }
    _ => scale = 1.0
  }

  // Covariance = scale * (X^T W X)^{-1}.
  let xtwx = glm_xtwx(x, w_final)
  let xtwx_inv = linalg_inverse(xtwx)
  let covariance = glm_matrix_scale(xtwx_inv, scale)
  let std_errors : Array[Double] = Array::make(p, 0.0)
  for j in 0..<p {
    let v = covariance[j][j]
    std_errors[j] = if v > 0.0 { v.sqrt() } else { 0.0 }
  }

  let log_likelihood = glm_total_loglik(family, y, mu, deviance)

  GlmResult::{
    coefficients: beta,
    std_errors,
    covariance,
    log_likelihood,
    deviance,
    df_residual,
    n_obs: n,
    converged,
    iterations: iters,
    scale,
  }
}

///|
/// X beta.
fn glm_linear_predictor(x : Array[Array[Double]], beta : Array[Double]) -> Array[Double] {
  let n = x.length()
  let eta = Array::make(n, 0.0)
  for i in 0..<n {
    let mut s = 0.0
    let row = x[i]
    for j in 0..<row.length() {
      s = s + row[j] * beta[j]
    }
    eta[i] = s
  }
  eta
}

///|
/// Sum of per-observation deviance.
fn glm_total_deviance(
  family : GlmFamily,
  y : Array[Double],
  mu : Array[Double]
) -> Double {
  let mut d = 0.0
  for i in 0..<y.length() {
    d = d + glm_deviance_one(family, y[i], mu[i])
  }
  d
}

///|
/// Total log-likelihood (Gaussian uses the MLE variance closed form).
fn glm_total_loglik(
  family : GlmFamily,
  y : Array[Double],
  mu : Array[Double],
  deviance : Double
) -> Double {
  match family {
    Gaussian => {
      let n = y.length()
      if n <= 1 || deviance <= 0.0 {
        return 0.0
      }
      let v = deviance / n.to_double()
      -n.to_double() / 2.0 * @math.ln(2.0 * 3.141592653589793 * v) - n.to_double() / 2.0
    }
    _ => {
      let mut ll = 0.0
      for i in 0..<y.length() {
        ll = ll + glm_loglik_one(family, y[i], mu[i])
      }
      ll
    }
  }
}

///|
/// Scale a matrix by a scalar.
fn glm_matrix_scale(a : Array[Array[Double]], s : Double) -> Array[Array[Double]] {
  let n = a.length()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i in 0..<n {
    let row = Array::make(n, 0.0)
    for j in 0..<n {
      row[j] = a[i][j] * s
    }
    out.push(row)
  }
  out
}

// ============================================================
// Wald tests
// ============================================================

///|
/// A Wald test result (single coefficient or omnibus).
pub struct WaldResult {
  statistic : Double
  df : Int
  p_value : Double
}

pub fn WaldResult::statistic(self : WaldResult) -> Double {
  self.statistic
}

pub fn WaldResult::df(self : WaldResult) -> Int {
  self.df
}

pub fn WaldResult::p_value(self : WaldResult) -> Double {
  self.p_value
}

///|
/// Single-coefficient Wald test: z = beta_j / se_j, two-sided normal p-value.
/// Tests H0: beta_j = 0.
pub fn glm_wald_test(result : GlmResult, coef_idx : Int) -> WaldResult {
  let b = result.coefficients[coef_idx]
  let se = result.std_errors[coef_idx]
  let z = if se > 0.0 { b / se } else { 0.0 }
  let p = stat_pvalue_zscore(z)
  WaldResult::{ statistic: z, df: 1, p_value: p }
}

///|
/// Omnibus Wald test for a subset of coefficients: tests H0: beta[indices] = 0.
///
/// Wald statistic = beta_Cᵀ cov_C⁻¹ beta_C ~ chi-square with k = indices.length()
/// degrees of freedom.
pub fn glm_wald_omnibus(result : GlmResult, indices : Array[Int]) -> WaldResult {
  let k = indices.length()
  if k == 0 {
    return WaldResult::{ statistic: 0.0, df: 0, p_value: 1.0 }
  }
  // Build sub-vector beta_C and sub-matrix cov_C.
  let beta_c = Array::make(k, 0.0)
  let cov_c = glm_dmatrix(k, k)
  for a in 0..<k {
    beta_c[a] = result.coefficients[indices[a]]
    for b in 0..<k {
      cov_c[a][b] = result.covariance[indices[a]][indices[b]]
    }
  }
  let cov_inv = linalg_inverse(cov_c)
  // Quadratic form beta_Cᵀ cov_C⁻¹ beta_C.
  let mut stat = 0.0
  for a in 0..<k {
    let mut row_sum = 0.0
    for b in 0..<k {
      row_sum = row_sum + beta_c[b] * cov_inv[a][b]
    }
    stat = stat + beta_c[a] * row_sum
  }
  if stat < 0.0 {
    stat = 0.0
  }
  let p = 1.0 - stat_chi2_cdf(stat, k)
  WaldResult::{ statistic: stat, df: k, p_value: p }
}

// ============================================================
// Likelihood-ratio test (LRT)
// ============================================================

///|
/// Likelihood-ratio test result comparing two nested models.
pub struct LrtResult {
  statistic : Double
  df : Int
  p_value : Double
}

pub fn LrtResult::statistic(self : LrtResult) -> Double {
  self.statistic
}

pub fn LrtResult::df(self : LrtResult) -> Int {
  self.df
}

pub fn LrtResult::p_value(self : LrtResult) -> Double {
  self.p_value
}

///|
/// Likelihood-ratio test of a full model against a reduced (nested) model.
///
/// `full` and `reduced` must be fitted on the same data with the same family;
/// `reduced` must be nested in `full` (fewer parameters). The test statistic
/// is -2*(ll_reduced - ll_full) ~ chi-square with
/// df = n_params(full) - n_params(reduced).
pub fn glm_lrt_compare(full : GlmResult, reduced : GlmResult) -> LrtResult {
  let stat = 2.0 * (full.log_likelihood - reduced.log_likelihood)
  let df_full = full.coefficients.length()
  let df_reduced = reduced.coefficients.length()
  let df = df_full - df_reduced
  let s = if stat < 0.0 { 0.0 } else { stat }
  let p = if df > 0 { 1.0 - stat_chi2_cdf(s, df) } else { 1.0 }
  LrtResult::{ statistic: s, df, p_value: p }
}