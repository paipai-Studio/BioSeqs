///|
/// ashr - Adaptive Shrinkage of effect sizes.
/// Inspired by Bioconductor {ashr} (Stephens M., 2017).
///
/// Implements empirical-Bayes shrinkage of noisy effect-size estimates under a
/// unimodal prior expressed as a finite mixture of zero-mean normals:
///
///   beta_hat_i ~ N(beta_i, se_i^2),   beta_i ~ sum_k pi_k N(0, sd_k^2)
///
/// The mixture weights pi_k are estimated by EM, giving per-observation
/// posterior means, posterior standard deviations and the local false sign
/// rate (lfsr). This mirrors the shrinkage used downstream of DESeq2 / edgeR /
/// limma to stabilise log fold-change estimates.

///|
/// Result of an adaptive shrinkage fit.
pub struct AshrFit {
  posterior_mean : Array[Double]
  posterior_sd : Array[Double]
  lfsr : Array[Double]
  prior_grid : Array[Double]
  prior_weights : Array[Double]
  log_likelihood : Double
  iterations : Int
  n : Int
}

// ---- univariate normal helpers -------------------------------------------

///|
/// Error function (Abramowitz & Stegun 7.1.26 approximation).
pub fn ashr_erf(x : Double) -> Double {
  let p = 0.3275911
  let a1 = 0.254829592
  let a2 = -0.284496736
  let a3 = 1.421413741
  let a4 = -1.453152027
  let a5 = 1.061405429
  let sign = if x < 0.0 { -1.0 } else { 1.0 }
  let ax = x.abs()
  let t = 1.0 / (1.0 + p * ax)
  let poly = ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t
  sign * (1.0 - poly * @math.exp(-ax * ax))
}

///|
/// Standard normal cumulative distribution function.
pub fn ashr_pnorm(x : Double) -> Double {
  0.5 * (1.0 + ashr_erf(x / (2.0).sqrt()))
}

///|
/// Normal probability density function.
pub fn ashr_dnorm(x : Double, mean : Double, sd : Double) -> Double {
  let s = if sd > 1.0e-12 { sd } else { 1.0e-12 }
  let z = (x - mean) / s
  @math.exp(-0.5 * z * z) / (s * (2.0 * 3.141592653589793).sqrt())
}

///|
/// Log normal density (numerically stable form).
fn ashr_dnorm_log(x : Double, mean : Double, sd : Double) -> Double {
  let s = if sd > 1.0e-12 { sd } else { 1.0e-12 }
  let z = (x - mean) / s
  -0.5 * z * z - @math.ln(s) - 0.5 * @math.ln(2.0 * 3.141592653589793)
}

///|
/// Default unimodal prior grid (standard deviations of the normal components).
/// The first entry (0.0) is a point mass at zero.
pub fn ashr_default_grid() -> Array[Double] {
  [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
}

// ---- EM fit ---------------------------------------------------------------

///|
/// Fit the unimodal normal-mixture prior by EM and return posterior summaries.
pub fn ashr_shrink(
  estimates : Array[Double],
  ses : Array[Double],
  grid_sd~ : Array[Double] = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0],
  max_iter~ : Int = 200,
  tol~ : Double = 1.0e-8,
) -> AshrFit {
  let n = estimates.length()
  let k = grid_sd.length()
  if n == 0 || k == 0 {
    return AshrFit::{
      posterior_mean: [],
      posterior_sd: [],
      lfsr: [],
      prior_grid: grid_sd,
      prior_weights: [],
      log_likelihood: 0.0,
      iterations: 0,
      n,
    }
  }

  // Responsibilities w[i * k + j]
  let w = Array::make(n * k, 0.0)
  let pi = Array::make(k, 1.0 / k.to_double())
  let mut log_likelihood = 0.0
  let mut iterations = 0
  let mut prev_ll = -1.0e300

  let mut iter = 0
  while iter < max_iter {
    let pi_sum = Array::make(k, 0.0)
    let mut new_ll = 0.0

    // ---- E-step ---------------------------------------------------------
    let mut i = 0
    while i < n {
      let est = estimates[i]
      let se = if ses[i] > 1.0e-12 { ses[i] } else { 1.0e-12 }
      let llk = Array::make(k, 0.0)
      let mut max_ll = -1.0e300
      let mut j = 0
      while j < k {
        let tot = (se * se + grid_sd[j] * grid_sd[j]).sqrt()
        let v = ashr_dnorm_log(est, 0.0, tot) + @math.ln(pi[j])
        llk[j] = v
        if v > max_ll {
          max_ll = v
        }
        j = j + 1
      }
      let mut sum_exp = 0.0
      j = 0
      while j < k {
        let e = @math.exp(llk[j] - max_ll)
        w[i * k + j] = e
        sum_exp = sum_exp + e
        j = j + 1
      }
      new_ll = new_ll + max_ll + @math.ln(sum_exp)
      j = 0
      while j < k {
        let r = w[i * k + j] / sum_exp
        w[i * k + j] = r
        pi_sum[j] = pi_sum[j] + r
        j = j + 1
      }
      i = i + 1
    }

    // ---- M-step ---------------------------------------------------------
    let mut jm = 0
    while jm < k {
      pi[jm] = pi_sum[jm] / n.to_double()
      jm = jm + 1
    }

    log_likelihood = new_ll
    iterations = iter + 1
    if (new_ll - prev_ll).abs() < tol {
      break
    }
    prev_ll = new_ll
    iter = iter + 1
  }

  // ---- posterior summaries ---------------------------------------------
  let post_mean = Array::make(n, 0.0)
  let post_sd = Array::make(n, 0.0)
  let lfsr = Array::make(n, 0.0)

  let mut i = 0
  while i < n {
    let est = estimates[i]
    let se = if ses[i] > 1.0e-12 { ses[i] } else { 1.0e-12 }
    let mut pm = 0.0
    let mut e_var = 0.0
    let mut e_mean_sq = 0.0
    let mut p_le = 0.0
    let mut p_ge = 0.0
    let mut j = 0
    while j < k {
      let wij = w[i * k + j]
      let sd = grid_sd[j]
      let denom = se * se + sd * sd
      let m = est * sd * sd / denom
      let v = se * se * sd * sd / denom
      pm = pm + wij * m
      e_var = e_var + wij * v
      e_mean_sq = e_mean_sq + wij * m * m
      let s = v.sqrt()
      let le = if s < 1.0e-12 { 0.5 } else { ashr_pnorm((0.0 - m) / s) }
      let ge = if s < 1.0e-12 { 0.5 } else { 1.0 - ashr_pnorm((0.0 - m) / s) }
      p_le = p_le + wij * le
      p_ge = p_ge + wij * ge
      j = j + 1
    }
    post_mean[i] = pm
    let total_var = e_var + e_mean_sq - pm * pm
    post_sd[i] = if total_var > 0.0 { total_var.sqrt() } else { 0.0 }
    lfsr[i] = if p_le < p_ge { p_le } else { p_ge }
    i = i + 1
  }

  AshrFit::{
    posterior_mean: post_mean,
    posterior_sd: post_sd,
    lfsr,
    prior_grid: grid_sd,
    prior_weights: pi,
    log_likelihood,
    iterations,
    n,
  }
}

///|
/// Convenience wrapper for z-scores (unit standard errors).
pub fn ashr_shrink_z(
  z : Array[Double],
  grid_sd~ : Array[Double] = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0],
  max_iter~ : Int = 200,
  tol~ : Double = 1.0e-8,
) -> AshrFit {
  let ses = Array::make(z.length(), 1.0)
  ashr_shrink(z, ses, grid_sd~, max_iter~, tol~)
}

// ---- post-fit utilities ---------------------------------------------------

///|
/// Number of observations whose lfsr is below the given cutoff.
pub fn ashr_count_significant(fit : AshrFit, lfsr_cutoff : Double) -> Int {
  let mut count = 0
  for v in fit.lfsr {
    if v < lfsr_cutoff {
      count = count + 1
    }
  }
  count
}

///|
/// Per-observation shrinkage factor posterior_mean / estimate (in [0, 1]).
/// Estimates of exactly zero yield a factor of 1.0.
pub fn ashr_shrinkage_factor(
  fit : AshrFit,
  estimates : Array[Double],
) -> Array[Double] {
  let out : Array[Double] = []
  let n = fit.n
  let mut i = 0
  while i < n {
    let est = estimates[i]
    if est.abs() < 1.0e-12 {
      out.push(1.0)
    } else {
      out.push(fit.posterior_mean[i] / est)
    }
    i = i + 1
  }
  out
}

///|
/// Index of the prior component carrying the largest weight.
pub fn ashr_mode_component(fit : AshrFit) -> Int {
  let mut best = 0
  let mut i = 1
  while i < fit.prior_weights.length() {
    if fit.prior_weights[i] > fit.prior_weights[best] {
      best = i
    }
    i = i + 1
  }
  best
}

///|
/// Human readable summary of an adaptive shrinkage fit.
pub fn ashr_summary(fit : AshrFit) -> String {
  let buf = StringBuilder()
  buf.write_string("AshrFit: n=\{fit.n}  iterations=\{fit.iterations}\n")
  buf.write_string("  log-likelihood: \{fit.log_likelihood}\n")
  buf.write_string("  prior components (sd, weight):\n")
  let mut i = 0
  while i < fit.prior_grid.length() {
    buf.write_string("    sd=\{fit.prior_grid[i]}  pi=\{fit.prior_weights[i]}\n")
    i = i + 1
  }
  let mut n_sig = 0
  for v in fit.lfsr {
    if v < 0.05 {
      n_sig = n_sig + 1
    }
  }
  buf.write_string("  lfsr < 0.05: \{n_sig} / \{fit.n}\n")
  buf.to_string()
}

///|
/// Small illustrative dataset: strong effects, weak effects and nulls.
pub fn ashr_example_data() -> (Array[Double], Array[Double]) {
  let estimates = [2.0, -1.8, 2.5, 0.3, -0.25, 0.05, 1.5, -0.1]
  let ses = [0.2, 0.2, 0.25, 0.2, 0.22, 0.18, 0.3, 0.2]
  (estimates, ses)
}
