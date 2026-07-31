///|
/// Bioconductor baySeq - Bayesian analysis of RNA-seq data.
/// Implements empirical Bayes estimation of dispersions and
/// posterior probability calculation for differential expression.
///
/// Reference: Bioconductor baySeq package
/// Model: Negative Binomial with conjugate Gamma prior on dispersion.

///|
/// Result for a single gene from baySeq analysis.
pub struct BayesGeneResult {
  gene_index : Int
  gene_name : String
  log_fold_change : Double
  dispersion : Double
  posterior_prob_de : Double
  map_expression_a : Double
  map_expression_b : Double
  is_de : Bool
}

///|
/// Full baySeq analysis result.
pub struct BayesResult {
  gene_results : Array[BayesGeneResult]
  n_de_genes : Int
  dispersion_prior_shape : Double
  dispersion_prior_rate : Double
  alpha_level : Double
}

///|
/// Estimate the dispersion parameter for a single gene using the method of moments.
/// For Negative Binomial: Var = mean + mean^2 / dispersion
/// => dispersion = mean^2 / (var - mean)
pub fn bayseq_estimate_dispersion_gene(counts_a : Array[Double], counts_b : Array[Double]) -> Double {
  let n_a = counts_a.length()
  let n_b = counts_b.length()
  if n_a < 2 || n_b < 2 {
    return 1.0
  }

  // Means for each group
  let mut sum_a = 0.0
  let mut sum_b = 0.0
  let mut i = 0
  while i < n_a {
    sum_a = sum_a + counts_a[i]
    i = i + 1
  }
  i = 0
  while i < n_b {
    sum_b = sum_b + counts_b[i]
    i = i + 1
  }
  let mean_a = sum_a / n_a.to_double()
  let mean_b = sum_b / n_b.to_double()
  let mean_pool = (sum_a + sum_b) / (n_a + n_b).to_double()

  if mean_pool < 0.5 {
    return 1.0
  }

  // Variances
  let mut var_sum_a = 0.0
  i = 0
  while i < n_a {
    let diff = counts_a[i] - mean_a
    var_sum_a = var_sum_a + diff * diff
    i = i + 1
  }

  let mut var_sum_b = 0.0
  i = 0
  while i < n_b {
    let diff = counts_b[i] - mean_b
    var_sum_b = var_sum_b + diff * diff
    i = i + 1
  }

  let var_pool = (var_sum_a + var_sum_b) / (n_a + n_b - 2).to_double()
  let disp = if var_pool > mean_pool {
    mean_pool * mean_pool / (var_pool - mean_pool)
  } else {
    mean_pool
  }

  // Bound dispersion to reasonable range
  if disp < 0.01 { 0.01 } else if disp > 100.0 { 100.0 } else { disp }
}

///|
/// Estimate dispersion prior parameters (shape, rate) from all genes.
pub fn bayseq_estimate_prior(counts : Array[Array[Double]], groups : Array[Int]) -> (Double, Double) {
  let n_genes = counts.length()
  if n_genes == 0 {
    return (2.0, 0.5)
  }
  let n_samples = counts[0].length()

  // Validate groups
  let mut n_group_0 = 0
  let mut n_group_1 = 0
  let mut gi = 0
  while gi < groups.length() {
    if groups[gi] == 0 { n_group_0 = n_group_0 + 1 }
    else if groups[gi] == 1 { n_group_1 = n_group_1 + 1 }
    gi = gi + 1
  }

  if n_group_0 < 2 || n_group_1 < 2 {
    return (2.0, 0.5)
  }

  let dispersions : Array[Double] = []
  let mut g = 0
  while g < n_genes {
    let c_a : Array[Double] = []
    let c_b : Array[Double] = []
    let mut s = 0
    while s < n_samples {
      if groups[s] == 0 {
        c_a.push(counts[g][s])
      } else {
        c_b.push(counts[g][s])
      }
      s = s + 1
    }
    let d = bayseq_estimate_dispersion_gene(c_a, c_b)
    if d > 0.0 {
      dispersions.push(d)
    }
    g = g + 1
  }

  // Method of moments for Gamma prior: mean=shape/rate, var=shape/rate^2
  let n = dispersions.length()
  if n < 2 {
    return (2.0, 0.5)
  }
  let mut sum = 0.0
  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    sum = sum + dispersions[i]
    sum_sq = sum_sq + dispersions[i] * dispersions[i]
    i = i + 1
  }
  let mean_d = sum / n.to_double()
  let var_d = sum_sq / n.to_double() - mean_d * mean_d

  if var_d > 0.0 && mean_d > 0.0 {
    let shape = mean_d * mean_d / var_d
    let rate = mean_d / var_d
    (shape, rate)
  } else {
    (2.0, 0.5)
  }
}

///|
/// Count number of samples in a given group.
pub fn hts_count_group(groups : Array[Int], group_id : Int) -> Int {
  let mut count = 0
  let mut i = 0
  while i < groups.length() {
    if groups[i] == group_id { count = count + 1 }
    i = i + 1
  }
  count
}

///|
/// Compute the log-likelihood ratio for differential expression under NB model.
/// Higher values = more evidence for DE.
pub fn bayseq_log_likelihood_ratio(
  counts_a : Array[Double],
  counts_b : Array[Double],
  disp : Double,
) -> Double {
  let n_a = counts_a.length()
  let n_b = counts_b.length()

  let mut sum_a = 0.0
  let mut sum_b = 0.0
  let mut i = 0
  while i < n_a {
    sum_a = sum_a + counts_a[i]
    i = i + 1
  }
  i = 0
  while i < n_b {
    sum_b = sum_b + counts_b[i]
    i = i + 1
  }

  let mean_a = sum_a / n_a.to_double()
  let mean_b = sum_b / n_b.to_double()

  if mean_a < 0.001 && mean_b < 0.001 {
    return 0.0
  }

  // Pooled mean (null model: no DE)
  let mean_pool = (sum_a + sum_b) / (n_a + n_b).to_double()

  // Log-likelihood under null (common mean)
  let ll_null = bayseq_nb_log_likelihood(counts_a, mean_pool, disp) +
    bayseq_nb_log_likelihood(counts_b, mean_pool, disp)

  // Log-likelihood under alternative (different means)
  let ll_alt = bayseq_nb_log_likelihood(counts_a, mean_a, disp) +
    bayseq_nb_log_likelihood(counts_b, mean_b, disp)

  ll_alt - ll_null
}

///|
/// Negative Binomial log-likelihood for a set of counts given mean and dispersion.
pub fn bayseq_nb_log_likelihood(counts : Array[Double], mean : Double, disp : Double) -> Double {
  if mean <= 0.0 || disp <= 0.0 {
    return 0.0
  }
  let mut ll = 0.0
  let r = disp // size parameter
  let p = disp / (disp + mean) // success probability
  let mut i = 0
  while i < counts.length() {
    let x = counts[i]
    // NB PMF: C(x + r - 1, x) * p^r * (1-p)^x
    // log form: lgamma(x+r) - lgamma(x+1) - lgamma(r) + r*ln(p) + x*ln(1-p)
    if x > 0.0 {
      let lgamma_xr = bayseq_lgamma(x + r)
      let lgamma_x1 = bayseq_lgamma(x + 1.0)
      let lgamma_r = bayseq_lgamma(r)
      let ln_p = if p > 0.0 { @math.ln(p) } else { -10000000000.0 }
      let ln_1mp = if (1.0 - p) > 0.0 { @math.ln(1.0 - p) } else { -10000000000.0 }
      ll = ll + lgamma_xr - lgamma_x1 - lgamma_r + r * ln_p + x * ln_1mp
    } else {
      let lgamma_r = bayseq_lgamma(r)
      let ln_p = if p > 0.0 { @math.ln(p) } else { -10000000000.0 }
      ll = ll - lgamma_r + r * ln_p
    }
    i = i + 1
  }
  ll
}

///|
/// Approximation of log-gamma function using Lanczos with g=5.
pub fn bayseq_lgamma(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  // For integer x <= 20, use exact recursive relation
  if x == 1.0 || x == 2.0 {
    return 0.0
  }
  if x < 1.0 {
    // Reflection formula: lgamma(x) = lgamma(x+1) - ln(x)
    return bayseq_lgamma(x + 1.0) - @math.ln(x)
  }
  // For x in [1, 2], interpolate
  if x < 2.0 {
    // lgamma(x) = lgamma(x+1) - ln(x)
    return bayseq_lgamma(x + 1.0) - @math.ln(x)
  }
  // Lanczos approximation for x >= 2
  // Coefficients for g=5
  let c = [
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5
  ]
  let y = x
  let mut tmp = x + 5.5
  tmp = tmp - (x + 0.5) * @math.ln(tmp)
  let mut ser = 1.000000000190015
  let mut i = 0
  while i < 6 {
    ser = ser + c[i] / (y + i.to_double())
    i = i + 1
  }
  -tmp + @math.ln(2.5066282746310005 * ser / y)
}

///|
/// Test for differential expression using baySeq approach.
/// Returns BayesGeneResult for each gene.
pub fn bayseq_test(
  counts : Array[Array[Double]],
  groups : Array[Int],
  gene_names : Array[String],
  dispersion_shape? : Double = 2.0,
  dispersion_rate? : Double = 0.5,
  alpha? : Double = 0.05,
) -> BayesResult {
  let n_genes = counts.length()
  let n_samples = counts[0].length()
  let n_group_0 = hts_count_group(groups, 0)
  let n_group_1 = hts_count_group(groups, 1)

  let results : Array[BayesGeneResult] = []
  let mut de_count = 0
  let mut g = 0
  while g < n_genes {
    let c_a : Array[Double] = []
    let c_b : Array[Double] = []
    let mut s = 0
    while s < n_samples {
      if groups[s] == 0 {
        c_a.push(counts[g][s])
      } else {
        c_b.push(counts[g][s])
      }
      s = s + 1
    }

    // Estimate dispersion for this gene
    let disp = bayseq_estimate_dispersion_gene(c_a, c_b)

    // Compute means
    let mut sum_a = 0.0
    let mut sum_b = 0.0
    let mut i = 0
    while i < c_a.length() { sum_a = sum_a + c_a[i]; i = i + 1 }
    i = 0
    while i < c_b.length() { sum_b = sum_b + c_b[i]; i = i + 1 }

    let mean_a = sum_a / n_group_0.to_double()
    let mean_b = sum_b / n_group_1.to_double()
    let lfc = if mean_a > 0.0 && mean_b > 0.0 {
      @math.ln(mean_b / mean_a)
    } else if mean_a <= 0.0 && mean_b > 0.0 {
      10.0
    } else if mean_a > 0.0 && mean_b <= 0.0 {
      -10.0
    } else {
      0.0
    }

    // Log-likelihood ratio as evidence for DE
    let llr = bayseq_log_likelihood_ratio(c_a, c_b, disp)

    // Convert to posterior probability via sigmoid
    let prob_de = 1.0 / (1.0 + @math.exp(-llr))

    let is_de = prob_de >= alpha
    if is_de { de_count = de_count + 1 }

    let name = if g < gene_names.length() { gene_names[g] } else { "Gene" + (g + 1).to_string() }

    results.push(BayesGeneResult::{
      gene_index: g,
      gene_name: name,
      log_fold_change: lfc,
      dispersion: disp,
      posterior_prob_de: prob_de,
      map_expression_a: mean_a,
      map_expression_b: mean_b,
      is_de: is_de,
    })
    g = g + 1
  }

  BayesResult::{
    gene_results: results,
    n_de_genes: de_count,
    dispersion_prior_shape: dispersion_shape,
    dispersion_prior_rate: dispersion_rate,
    alpha_level: alpha,
  }
}

///|
/// Get significant DE genes from a BayesResult.
pub fn bayseq_get_de_genes(result : BayesResult) -> Array[BayesGeneResult] {
  let de : Array[BayesGeneResult] = []
  let mut i = 0
  while i < result.gene_results.length() {
    if result.gene_results[i].is_de {
      de.push(result.gene_results[i])
    }
    i = i + 1
  }
  de
}

///|
/// Get top DE genes sorted by absolute log fold change.
pub fn bayseq_get_top_de(result : BayesResult, n : Int) -> Array[BayesGeneResult] {
  let all = result.gene_results
  // Sort by absolute log_fold_change descending (bubble sort)
  let mut i = 0
  while i < all.length() {
    let mut j = i + 1
    while j < all.length() {
      let abs_i = if all[i].log_fold_change >= 0.0 { all[i].log_fold_change } else { -all[i].log_fold_change }
      let abs_j = if all[j].log_fold_change >= 0.0 { all[j].log_fold_change } else { -all[j].log_fold_change }
      if abs_j > abs_i {
        let tmp = all[i]
        all[i] = all[j]
        all[j] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  let top : Array[BayesGeneResult] = []
  let count = if n < all.length() { n } else { all.length() }
  i = 0
  while i < count {
    top.push(all[i])
    i = i + 1
  }
  top
}

///|
/// Generate sample RNA-seq data with clear differential expression signal.
/// Returns (counts, groups, gene_names) with some genes truly DE.
pub fn bayseq_sample_data() -> (Array[Array[Double]], Array[Int], Array[String]) {
  let n_genes = 30
  let n_samples = 8
  let counts : Array[Array[Double]] = []
  let gene_names : Array[String] = []
  let groups = [0, 0, 0, 0, 1, 1, 1, 1]

  let mut g = 0
  while g < n_genes {
    let row : Array[Double] = []
    let mut s = 0
    while s < n_samples {
      let base_a = 500.0
      let base_b = if g < 8 { 1500.0 } else if g < 16 { 200.0 } else { 500.0 }
      let base = if groups[s] == 0 { base_a } else { base_b }
      let noise = (((g * 7 + s * 17) % 100).to_double() / 100.0) * base * 0.3
      let val = base + noise
      row.push(val)
      s = s + 1
    }
    counts.push(row)
    gene_names.push("Gene" + (g + 1).to_string())
    g = g + 1
  }
  (counts, groups, gene_names)
}

///|
/// Summarize the baySeq result.
pub fn bayseq_summary(result : BayesResult) -> String {
  "baySeq Analysis Summary:\n" +
    "  Total genes: " + result.gene_results.length().to_string() + "\n" +
    "  DE genes (posterior >= " + result.alpha_level.to_string() + "): " +
    result.n_de_genes.to_string() + "\n" +
    "  Dispersion prior: Gamma(shape=" +
    result.dispersion_prior_shape.to_string() +
    ", rate=" + result.dispersion_prior_rate.to_string() + ")"
}
