///|
/// scanpy `experimental.pp.sctransform` (regularized negative-binomial
/// normalization, Hafemeister & Satija 2019).  For every gene a negative
/// binomial GLM `y ~ 1 + log_umi` (log link) is fitted by IRLS; the dispersion
/// theta is estimated by the Pearson moment method.  Gene-wise parameters are
/// then shrunk toward the pooled (median) estimates — a deterministic,
/// simplified replacement for sctransform's regularized-model step — and
/// Pearson residuals `(y - mu) / sqrt(mu + mu^2 / theta)` are computed and
/// clipped.

///|
/// Result of sctransform: Pearson residuals (clipped and unclipped) for the
/// genes that passed the `min_cells` filter (`kept_genes` indexes into the
/// input columns), the regularized per-gene parameters and the pooled
/// estimates used for shrinkage.
pub struct ScanpySctransform {
  residuals : Array[Array[Double]]
  residuals_clipped : Array[Array[Double]]
  theta : Array[Double]
  intercept : Array[Double]
  beta : Array[Double]
  pooled_intercept : Double
  pooled_beta : Double
  theta_median : Double
  kept_genes : Array[Int32]
  n_cells : Int32
}

///|
/// Median of a non-empty slice.
fn sct_median(a : Array[Double]) -> Double {
  let n = a.length()
  if n == 0 {
    return 0.0
  }
  let s = a.copy()
  s.sort()
  if n % 2 == 1 {
    s[n / 2]
  } else {
    (s[n / 2 - 1] + s[n / 2]) / 2.0
  }
}

///|
/// Fit an NB GLM `y ~ 1 + t` (log link) by IRLS with fixed theta; returns
/// (intercept, slope).  Starting values come from an OLS fit on log(y + 1).
fn sct_fit_glm(
  y : Array[Double],
  t : Array[Double],
  theta : Double,
  max_iter : Int32,
) -> (Double, Double) {
  let n = y.length()
  if n == 0 {
    return (0.0, 0.0)
  }
  // OLS init on log(y + 1)
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sxx = 0.0
  let mut sxy = 0.0
  for i = 0; i < n; i = i + 1 {
    let lx = t[i]
    let ly = @math.ln(y[i] + 1.0)
    sx = sx + lx
    sy = sy + ly
    sxx = sxx + lx * lx
    sxy = sxy + lx * ly
  }
  let denom = n.to_double() * sxx - sx * sx
  let mut b1 = if denom.abs() > 1.0e-12 {
    (n.to_double() * sxy - sx * sy) / denom
  } else {
    0.0
  }
  let mut b0 = sy / n.to_double() - b1 * sx / n.to_double()
  // IRLS
  for _iter = 0; _iter < max_iter; _iter = _iter + 1 {
    let mut sw = 0.0
    let mut swx = 0.0
    let mut swz = 0.0
    let mut swzx = 0.0
    let mut swxx = 0.0
    for i = 0; i < n; i = i + 1 {
      let eta = b0 + b1 * t[i]
      let mu = if eta > 30.0 { @math.exp(30.0) } else { @math.exp(eta) }
      let w = mu / (1.0 + mu / theta)
      if w <= 0.0 {
        continue
      }
      let z = eta + (y[i] - mu) / mu
      sw = sw + w
      swx = swx + w * t[i]
      swz = swz + w * z
      swzx = swzx + w * z * t[i]
      swxx = swxx + w * t[i] * t[i]
    }
    if sw <= 0.0 {
      break
    }
    // 2x2 solve: [sw, swx; swx, swxx] [b0; b1] = [swz; swzx]
    let det2 = sw * swxx - swx * swx
    if det2.abs() < 1.0e-12 {
      break
    }
    let new_b0 = (swz * swxx - swzx * swx) / det2
    let new_b1 = (sw * swzx - swx * swz) / det2
    if (new_b0 - b0).abs() < 1.0e-10 && (new_b1 - b1).abs() < 1.0e-10 {
      b0 = new_b0
      b1 = new_b1
      break
    }
    b0 = new_b0
    b1 = new_b1
  }
  (b0, b1)
}

///|
/// Negative-binomial dispersion theta by the Pearson moment method:
/// `sum(mu^2) / (sum((y - mu)^2) - sum(mu))`, clipped to [0.01, 1e10] (small
/// overdispersion clips toward Poisson at 1e10).
fn sct_theta_moment(
  y : Array[Double],
  t : Array[Double],
  b0 : Double,
  b1 : Double,
) -> Double {
  let n = y.length()
  let mut sum_mu2 = 0.0
  let mut sum_res2 = 0.0
  let mut sum_mu = 0.0
  for i = 0; i < n; i = i + 1 {
    let eta = b0 + b1 * t[i]
    let mu = if eta > 30.0 { @math.exp(30.0) } else { @math.exp(eta) }
    let r = y[i] - mu
    sum_mu2 = sum_mu2 + mu * mu
    sum_res2 = sum_res2 + r * r
    sum_mu = sum_mu + mu
  }
  let denom = sum_res2 - sum_mu
  if denom <= 0.0 {
    return 1.0e10
  }
  let theta = sum_mu2 / denom
  if theta < 0.01 {
    0.01
  } else if theta > 1.0e10 {
    1.0e10
  } else {
    theta
  }
}

///|
/// Regularized sctransform.  For each gene expressed in at least `min_cells`
/// cells an NB GLM against `log10(count depth)` is fitted, the dispersion is
/// estimated by the Pearson moment method, and the parameters are shrunk
/// toward the pooled estimates with weight `n_cells / (n_cells + 30)` per
/// gene (deterministic stand-in for the regularized model).  Pearson
/// residuals `(y - mu) / sqrt(mu + mu^2 / theta_reg)` are returned clipped to
/// `[-clip_range, clip_range]` and unclipped.
pub fn scanpy_sctransform(
  x : Array[Array[Double]],
  min_cells~ : Int32 = 5,
  clip_range~ : Double = 10.0,
  max_iter~ : Int32 = 5,
) -> ScanpySctransform {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  if n_cells == 0 || n_genes == 0 {
    return ScanpySctransform::{
      residuals: [],
      residuals_clipped: [],
      theta: [],
      intercept: [],
      beta: [],
      pooled_intercept: 0.0,
      pooled_beta: 0.0,
      theta_median: 0.0,
      kept_genes: [],
      n_cells: 0,
    }
  }
  // log10 count depth per cell (guard against empty cells)
  let log_umi : Array[Double] = Array(capacity=n_cells)
  for i = 0; i < n_cells; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < n_genes; j = j + 1 {
      s = s + x[i][j]
    }
    log_umi.push(@math.log10(if s > 1.0 { s } else { 1.0 }))
  }
  // gene filter
  let kept : Array[Int32] = Array(capacity=n_genes)
  for j = 0; j < n_genes; j = j + 1 {
    let mut cnt = 0
    for i = 0; i < n_cells; i = i + 1 {
      if x[i][j] > 0.0 {
        cnt = cnt + 1
      }
    }
    if cnt >= min_cells {
      kept.push(j)
    }
  }
  let nk = kept.length()
  if nk == 0 {
    return ScanpySctransform::{
      residuals: [],
      residuals_clipped: [],
      theta: [],
      intercept: [],
      beta: [],
      pooled_intercept: 0.0,
      pooled_beta: 0.0,
      theta_median: 0.0,
      kept_genes: [],
      n_cells,
    }
  }
  // per-gene GLM fits
  let intercept = Array::make(nk, 0.0)
  let beta = Array::make(nk, 0.0)
  let theta = Array::make(nk, 0.0)
  for t = 0; t < nk; t = t + 1 {
    let j = kept[t]
    let y : Array[Double] = Array(capacity=n_cells)
    for i = 0; i < n_cells; i = i + 1 {
      y.push(x[i][j])
    }
    let (b0, b1) = sct_fit_glm(y, log_umi, 100.0, max_iter)
    intercept[t] = b0
    beta[t] = b1
    theta[t] = sct_theta_moment(y, log_umi, b0, b1)
  }
  // pooled (median) estimates and shrinkage
  let pooled_intercept = sct_median(intercept)
  let pooled_beta = sct_median(beta)
  let theta_med = sct_median(theta)
  let intercept_reg = Array::make(nk, 0.0)
  let beta_reg = Array::make(nk, 0.0)
  let theta_reg = Array::make(nk, 0.0)
  let n_f = n_cells.to_double()
  let shrink = n_f / (n_f + 30.0)
  for t = 0; t < nk; t = t + 1 {
    intercept_reg[t] = shrink * intercept[t] + (1.0 - shrink) * pooled_intercept
    beta_reg[t] = shrink * beta[t] + (1.0 - shrink) * pooled_beta
    theta_reg[t] = theta_med
  }
  // Pearson residuals
  let residuals : Array[Array[Double]] = Array(capacity=n_cells)
  let residuals_clipped : Array[Array[Double]] = Array(capacity=n_cells)
  for i = 0; i < n_cells; i = i + 1 {
    let row : Array[Double] = Array(capacity=nk)
    let rowc : Array[Double] = Array(capacity=nk)
    for t = 0; t < nk; t = t + 1 {
      let j = kept[t]
      let y = x[i][j]
      let eta = intercept_reg[t] + beta_reg[t] * log_umi[i]
      let mu = if eta > 30.0 { @math.exp(30.0) } else { @math.exp(eta) }
      let nb_var = mu + mu * mu / theta_reg[t]
      let r = if nb_var > 0.0 { (y - mu) / nb_var.sqrt() } else { 0.0 }
      row.push(r)
      rowc.push(
        if r > clip_range {
          clip_range
        } else if r < -clip_range {
          -clip_range
        } else {
          r
        },
      )
    }
    residuals.push(row)
    residuals_clipped.push(rowc)
  }
  ScanpySctransform::{
    residuals,
    residuals_clipped,
    theta: theta_reg,
    intercept: intercept_reg,
    beta: beta_reg,
    pooled_intercept,
    pooled_beta,
    theta_median: theta_med,
    kept_genes: kept,
    n_cells,
  }
}
