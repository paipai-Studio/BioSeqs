// scvi_ppc.mbt - Posterior Predictive Checks
// Inspired by scvi-tools criticism module (PosteriorPredictiveCheck)
// Model criticism by comparing observed vs generated data

///|
pub struct PPCConfig {
  n_samples : Int32
  gene_likelihood : String
}

///|
pub struct PPCResult {
  // Per-gene mean expression: observed vs expected
  observed_mean : Array[Double]
  expected_mean : Array[Double]
  // Per-gene variance: observed vs expected
  observed_var : Array[Double]
  expected_var : Array[Double]
  // Per-gene zero fraction: observed vs expected
  observed_zero_frac : Array[Double]
  expected_zero_frac : Array[Double]
  // Coefficient of variation
  observed_cv : Array[Double]
  expected_cv : Array[Double]
  // p-values per gene
  p_values : Array[Double]
}

///|
/// Compute observed statistics from count data
pub fn ppc_observed_stats(
  x : Array[Array[Double]],
) -> (Array[Double], Array[Double], Array[Double], Array[Double]) {
  let n_cells = x.length()
  let n_genes = x[0].length()
  let means = Array(capacity=n_genes)
  let vars = Array(capacity=n_genes)
  let zero_fracs = Array(capacity=n_genes)
  let cvs = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut sum = 0.0
    let mut sum_sq = 0.0
    let mut n_zero = 0
    for i in 0..<n_cells {
      let val = x[i][g]
      sum += val
      sum_sq += val * val
      if val < 0.5 {
        n_zero += 1
      }
    }
    let mean = sum / n_cells.to_double()
    let variance = sum_sq / n_cells.to_double() - mean * mean
    means.push(mean)
    vars.push(variance)
    zero_fracs.push(n_zero.to_double() / n_cells.to_double())
    if mean > 0.0 {
      cvs.push(variance.sqrt() / mean)
    } else {
      cvs.push(0.0)
    }
  }
  (means, vars, zero_fracs, cvs)
}

///|
/// Generate posterior predictive samples from Negative Binomial
pub fn ppc_generate_nb(
  means : Array[Double],
  _dispersions : Array[Double],
  n_samples : Int32,
) -> Array[Array[Double]] {
  let n_genes = means.length()
  let result = Array(capacity=n_samples)
  for s in 0..<n_samples {
    let sample = Array(capacity=n_genes)
    for g in 0..<n_genes {
      let mu = means[g]
      let theta = _dispersions[g]
      // Simple NB sampling via Poisson-Gamma mixture
      let r = theta
      let p = theta / (theta + mu)
      // Generate from gamma-poisson
      let gamma_val = ppc_gamma_sample(r, p, s * n_genes + g)
      let count = ppc_poisson_sample(gamma_val, s * n_genes + g + 1000)
      sample.push(count.to_double())
    }
    result.push(sample)
  }
  result
}

///|
fn ppc_gamma_sample(shape : Double, scale : Double, seed : Int32) -> Double {
  // Marsaglia-Tsang method for gamma sampling
  if shape < 1.0 {
    let u = scvi_rand(seed)
    return ppc_gamma_sample(shape + 1.0, scale, seed + 1) *
      @math.pow(u, 1.0 / shape)
  }
  let d = shape - 1.0 / 3.0
  let c = 1.0 / (9.0 * d).sqrt()
  let mut x = 0.0
  for _i in 0..<100 {
    let u1 = scvi_rand(seed + _i * 2)
    let u2 = scvi_rand(seed + _i * 2 + 1)
    let v = @math.pow(1.0 + c * (2.0 * u1 - 1.0), 3.0)
    if v > 0.0 {
      let log_u2 = @math.ln(u2 + 0.0000000001)
      if log_u2 < 0.5 * v * v - v + d {
        x = d * v * scale
        break
      }
    }
  }
  if x <= 0.0 {
    shape * scale
  } else {
    x
  }
}

///|
fn ppc_poisson_sample(lam : Double, seed : Int32) -> Int32 {
  // Knuth's algorithm
  let exp_neg_lam = @math.exp(-lam)
  let mut k = 0
  let mut p = 1.0
  while true {
    k = k + 1
    p = p * scvi_rand(seed + k)
    if p <= exp_neg_lam {
      break
    }
    if k > 10000 {
      break
    }
  }
  k - 1
}

///|
/// Run posterior predictive check
pub fn ppc_run(
  x : Array[Array[Double]],
  decoded_expr : Array[Array[Double]],
  _dispersions : Array[Double],
  _config : PPCConfig,
) -> PPCResult {
  let (obs_mean, obs_var, obs_zero, obs_cv) = ppc_observed_stats(x)
  let n_genes = obs_mean.length()
  // Generate posterior predictive samples
  let exp_means = Array(capacity=n_genes)
    |> a => {
      for _i in 0..<n_genes {
        a.push(0.0)
      }
      a
    }
  let exp_vars = Array(capacity=n_genes)
    |> a => {
      for _i in 0..<n_genes {
        a.push(0.0)
      }
      a
    }
  let exp_zeros = Array(capacity=n_genes)
    |> a => {
      for _i in 0..<n_genes {
        a.push(0.0)
      }
      a
    }
  let exp_cvs = Array(capacity=n_genes)
    |> a => {
      for _i in 0..<n_genes {
        a.push(0.0)
      }
      a
    }
  // Average decoded expression as expected means
  let n_cells = decoded_expr.length()
  for g in 0..<n_genes {
    let mut sum = 0.0
    let mut sum_sq = 0.0
    let mut n_zero = 0
    for i in 0..<n_cells {
      let val = decoded_expr[i][g]
      sum += val
      sum_sq += val * val
      if val < 0.5 {
        n_zero += 1
      }
    }
    exp_means[g] = sum / n_cells.to_double()
    exp_vars[g] = sum_sq / n_cells.to_double() - exp_means[g] * exp_means[g]
    exp_zeros[g] = n_zero.to_double() / n_cells.to_double()
    if exp_means[g] > 0.0 {
      exp_cvs[g] = exp_vars[g].sqrt() / exp_means[g]
    }
  }
  // Compute p-values
  let p_values = Array(capacity=n_genes)
  for g in 0..<n_genes {
    // Two-sided test: how extreme is observed mean
    let diff = (obs_mean[g] - exp_means[g]).abs()
    let se = exp_vars[g] / n_cells.to_double().sqrt() + 0.0000000001
    let z = diff / (se + 0.0000000001)
    // Approximate p-value from z-score
    let p = 2.0 * (1.0 - ppc_normal_cdf(z))
    p_values.push(p)
  }
  {
    observed_mean: obs_mean,
    expected_mean: exp_means,
    observed_var: obs_var,
    expected_var: exp_vars,
    observed_zero_frac: obs_zero,
    expected_zero_frac: exp_zeros,
    observed_cv: obs_cv,
    expected_cv: exp_cvs,
    p_values,
  }
}

///|
fn ppc_normal_cdf(x : Double) -> Double {
  // Approximate normal CDF using error function
  0.5 * (1.0 + ppc_erf(x / 2.0.sqrt()))
}

///|
fn ppc_erf(x : Double) -> Double {
  // Abramowitz and Stegun approximation
  let t = 1.0 / (1.0 + 0.3275911 * x.abs())
  let poly = t *
    (
      0.254829592 +
      t *
      (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))
    )
  let result = 1.0 - poly * @math.exp(-x * x)
  if x >= 0.0 {
    result
  } else {
    -result
  }
}

///|
/// Create criticism report summary
pub fn ppc_create_report(result : PPCResult) -> String {
  let n_genes = result.p_values.length()
  let mut _n_significant = 0
  for g in 0..<n_genes {
    if result.p_values[g] < 0.05 {
      _n_significant += 1
    }
  }
  let mut mean_pval = 0.0
  for g in 0..<n_genes {
    mean_pval += result.p_values[g]
  }
  mean_pval /= n_genes.to_double()
  "Posterior Predictive Check Report:\n  Genes tested: {n_genes}\n  Significant (p<0.05): {n_significant}\n  Mean p-value: {mean_pval}\n  Model fit: " +
  (if mean_pval > 0.3 {
    "Good"
  } else if mean_pval > 0.1 {
    "Moderate"
  } else {
    "Poor"
  })
}

///|
pub fn PPCConfig::default() -> PPCConfig {
  { n_samples: 0, gene_likelihood: "", }
}
