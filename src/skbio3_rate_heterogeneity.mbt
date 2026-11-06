// Scikit-bio third-round: rate heterogeneity (items 9-11)
// discrete_gamma_rates, site_specific_rates, invariant_sites_fraction

///|
pub fn skbio3_discrete_gamma_rates(
  alpha : Double,
  n_categories : Int,
) -> Array[Double] {
  // Discrete gamma approximation: divide Gamma(alpha, 1/beta) into n equal-probability categories
  // Return the mean rate of each category
  // Uses simple quantile approximation via inverse CDF
  if n_categories <= 0 || alpha <= 0.0 {
    return [1.0]
  }
  let rates : Array[Double] = []
  // For discrete gamma, we divide [0,1] into n equal intervals and compute
  // the conditional mean in each interval
  // Gamma CDF approximation using series expansion
  for i in 0..<n_categories {
    let lo = i.to_double() / n_categories.to_double()
    let hi = (i + 1).to_double() / n_categories.to_double()
    // Approximate the mean rate in [lo, hi] quantile range
    // Use the midpoint quantile and invert the CDF
    let mid_q = (lo + hi) / 2.0
    let rate = skbio3_gamma_inv_cdf(mid_q, alpha, 1.0 / alpha)
    rates.push(rate)
  }
  // Normalize so mean = 1
  let mut sum = 0.0
  for i in 0..<rates.length() {
    sum = sum + rates[i]
  }
  if sum > 0.0 {
    for i in 0..<rates.length() {
      rates[i] = rates[i] * rates.length().to_double() / sum
    }
  }
  rates
}

fn skbio3_gamma_inv_cdf(
  p : Double,
  alpha : Double,
  beta : Double,
) -> Double {
  // Approximate inverse CDF of Gamma(alpha, beta) using Newton's method
  // Start with a reasonable initial guess: alpha * beta (the mean)
  let mut x = alpha * beta
  if x <= 0.0 {
    x = 0.01
  }
  // Newton's method: x_{n+1} = x_n - (CDF(x) - p) / PDF(x)
  for _ in 0..<50 {
    let cdf = skbio3_gamma_cdf(x, alpha, beta)
    let pdf = skbio3_gamma_pdf(x, alpha, beta)
    if pdf < 0.0000000000000001 {
      break
    }
    let delta = (cdf - p) / pdf
    x = x - delta
    if x <= 0.0 {
      x = 0.001
    }
    if delta < 0.0 && delta > -0.0000000001 {
      break
    }
    if delta >= 0.0 && delta < 0.0000000001 {
      break
    }
  }
  x
}

fn skbio3_gamma_pdf(
  x : Double,
  alpha : Double,
  beta : Double,
) -> Double {
  // Gamma PDF: (beta^alpha / Gamma(alpha)) * x^(alpha-1) * exp(-beta*x)
  if x <= 0.0 {
    return 0.0
  }
  let log_pdf = alpha * @math.ln(beta) - skbio3_log_gamma(alpha) +
    (alpha - 1.0) * @math.ln(x) - beta * x
  @math.exp(log_pdf)
}

fn skbio3_gamma_cdf(
  x : Double,
  alpha : Double,
  beta : Double,
) -> Double {
  // Gamma CDF via series expansion (lower incomplete gamma / complete gamma)
  if x <= 0.0 {
    return 0.0
  }
  // Use series expansion: gamma(s, x) / Gamma(s) = x^s * e^(-x) * sum_{k=0}^inf x^k / (s*(s+1)*...*(s+k))
  let scaled_x = beta * x
  let mut sum = 1.0 / alpha
  let mut term = 1.0 / alpha
  for k in 1..<100 {
    term = term * scaled_x / (alpha + k.to_double())
    sum = sum + term
    if term < 0.000000000000001 * sum {
      break
    }
  }
  let result = @math.exp(alpha * @math.ln(scaled_x) - scaled_x - skbio3_log_gamma(alpha)) * sum
  if result > 1.0 {
    1.0
  } else if result < 0.0 {
    0.0
  } else {
    result
  }
}

fn skbio3_log_gamma(
  x : Double,
) -> Double {
  // Lanczos approximation for log(Gamma(x))
  let g = 7.0
  let c = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61501716214001,
    12.5073436283956,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
  ]
  if x < 0.5 {
    // Reflection formula
    @math.ln(3.141592653589793 / (@math.sin(3.141592653589793 * x))) - skbio3_log_gamma(1.0 - x)
  } else {
    let z = x - 1.0
    let mut x0 = c[0]
    for i in 1..<9 {
      x0 = x0 + c[i] / (z + i.to_double())
    }
    let t = z + g + 0.5
    0.5 * @math.ln(2.0 * 3.141592653589793) + (z + 0.5) * @math.ln(t) - t + @math.ln(x0)
  }
}

///|
pub fn skbio3_site_specific_rates(
  alignment : Array[String],
  _tree_distances : Array[Double],
) -> Array[Double] {
  // Estimate relative rate for each site
  // Simple approach: rate ∝ variability at that site across sequences
  let n_sites = if alignment.is_empty() { 0 } else { alignment[0].length() }
  if n_sites == 0 {
    return []
  }
  let n_seqs = alignment.length()
  let rates : Array[Double] = []
  for site in 0..<n_sites {
    // Count unique characters at this site
    let chars : Array[UInt16] = []
    for seq_idx in 0..<n_seqs {
      let c = alignment[seq_idx][site]
      let mut found = false
      for k in 0..<chars.length() {
        if chars[k] == c {
          found = true
          break
        }
      }
      if not(found) {
        chars.push(c)
      }
    }
    // Rate proportional to number of unique chars (more variation = higher rate)
    let rate = (chars.length() - 1).to_double()
    rates.push(rate)
  }
  // Normalize to mean 1
  let mut sum = 0.0
  for i in 0..<rates.length() {
    sum = sum + rates[i]
  }
  if sum > 0.0 {
    for i in 0..<rates.length() {
      rates[i] = rates[i] * rates.length().to_double() / sum
    }
  }
  rates
}

///|
pub fn skbio3_invariant_sites_fraction(
  alignment : Array[String],
) -> Double {
  // Estimate fraction of invariant sites
  let n_sites = if alignment.is_empty() { 0 } else { alignment[0].length() }
  if n_sites == 0 {
    return 0.0
  }
  let n_seqs = alignment.length()
  let mut invariant_count = 0
  for site in 0..<n_sites {
    let first_char = alignment[0][site]
    let mut is_invariant = true
    for seq_idx in 1..<n_seqs {
      if alignment[seq_idx][site] != first_char {
        is_invariant = false
        break
      }
    }
    if is_invariant {
      invariant_count = invariant_count + 1
    }
  }
  invariant_count.to_double() / n_sites.to_double()
}
