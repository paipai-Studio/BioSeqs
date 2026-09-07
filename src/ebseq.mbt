///|
/// EBSeq - Empirical Bayesian approach for identifying differentially expressed genes.
/// Inspired by Bioconductor EBSeq package.
///
/// Model: Gene counts ~ NB(mean, dispersion) where mean = size_factor * q_gene
///   and variance = mean + dispersion * mean^2.
///
/// Two-condition analysis computes posterior probability of DE (PPDE)
/// using empirical Bayes with marginal likelihoods under EE and DE models.

///|
/// Dataset for EBSeq analysis containing count matrix and experimental design.
pub struct EBSeqDataSet {
  counts : Array[Array[Int]]
  gene_names : Array[String]
  sample_names : Array[String]
  conditions : Array[Int]
  size_factors : Array[Double]
  dispersions : Array[Double]
}

///|
/// Result of EBSeq differential expression analysis.
pub struct EBSeqResult {
  gene_names : Array[String]
  ppde : Array[Double]
  log2_fc : Array[Double]
  mean_cond1 : Array[Double]
  mean_cond2 : Array[Double]
  flag_de : Array[Bool]
}

///|
pub fn EBSeqDataSet::new(
  counts : Array[Array[Int]],
  gene_names : Array[String],
  sample_names : Array[String],
  conditions : Array[Int],
) -> EBSeqDataSet {
  let n_genes = counts.length()
  let n_samples = if counts.length() > 0 { counts[0].length() } else { 0 }

  EBSeqDataSet::{
    counts,
    gene_names,
    sample_names,
    conditions,
    size_factors: Array::make(n_samples, 1.0),
    dispersions: Array::make(n_genes, 0.1),
  }
}

///|
/// Log-gamma function using Stirling's approximation with correction terms.
/// lgamma(x) ≈ 0.5*ln(2π) + (x-0.5)*ln(x) - x + 1/(12x) - 1/(360x^3)
fn ebseq_lgamma(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x == 1.0 || x == 2.0 {
    return 0.0
  }
  // Use recurrence lgamma(x) = lgamma(x+1) - ln(x) for small x
  if x < 8.0 {
    return ebseq_lgamma(x + 1.0) - @math.ln(x)
  }
  // Stirling's approximation with correction terms for x >= 8
  let ln_x = @math.ln(x)
  0.5 * 0.9189385332046727 + (x - 0.5) * ln_x - x + 1.0 / (12.0 * x) -
    1.0 / (360.0 * x * x * x)
}

///|
/// Negative Binomial log-PMF.
/// log P(k|mu,r) = lgamma(k+r) - lgamma(k+1) - lgamma(r)
///                + r*log(r/(r+mu)) + k*log(mu/(r+mu))
/// where r = 1/dispersion is the size parameter.
fn ebseq_nb_log_pmf(k : Double, mu : Double, r : Double) -> Double {
  if mu <= 0.0 || r <= 0.0 {
    return -1.0e300
  }
  let r_plus_mu = r + mu
  let log_r_over_rmu = @math.ln(r / r_plus_mu)
  let log_mu_over_rmu = @math.ln(mu / r_plus_mu)
  ebseq_lgamma(k + r) - ebseq_lgamma(k + 1.0) - ebseq_lgamma(r) +
    r * log_r_over_rmu +
    k * log_mu_over_rmu
}

///|
/// Compute geometric mean of an array of positive values.
fn ebseq_geometric_mean(vals : Array[Double]) -> Double {
  let n = vals.length()
  if n == 0 {
    return 0.0
  }
  let mut sum_log = 0.0
  let mut count = 0
  for v in vals {
    if v > 0.0 {
      sum_log = sum_log + @math.ln(v)
      count = count + 1
    }
  }
  if count == 0 {
    0.0
  } else {
    @math.exp(sum_log / count.to_double())
  }
}

///|
/// Compute median of an array of doubles.
fn ebseq_median(vals : Array[Double]) -> Double {
  let n = vals.length()
  if n == 0 {
    return 0.0
  }
  // Copy and sort
  let sorted : Array[Double] = Array::new()
  for v in vals {
    sorted.push(v)
  }
  for i = 0; i < n - 1; i = i + 1 {
    for j = 0; j < n - 1 - i; j = j + 1 {
      if sorted[j] > sorted[j + 1] {
        let tmp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = tmp
      }
    }
  }
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

///|
/// Estimate size factors using median-of-ratios method (like DESeq).
/// For each sample, compute the median of ratios to geometric mean across genes.
pub fn ebseq_estimate_size_factors(data : EBSeqDataSet) -> EBSeqDataSet {
  let n_genes = data.counts.length()
  let n_samples = data.sample_names.length()

  // Compute geometric mean per gene
  let geo_means : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let gene_vals : Array[Double] = Array::new()
    for j = 0; j < n_samples; j = j + 1 {
      gene_vals.push(data.counts[i][j].to_double())
    }
    geo_means.push(ebseq_geometric_mean(gene_vals))
  }

  // Compute size factor per sample as median of ratios
  let size_factors : Array[Double] = Array::new()
  for j = 0; j < n_samples; j = j + 1 {
    let ratios : Array[Double] = Array::new()
    for i = 0; i < n_genes; i = i + 1 {
      if geo_means[i] > 0.0 && data.counts[i][j] > 0 {
        ratios.push(data.counts[i][j].to_double() / geo_means[i])
      }
    }
    let median_val = if ratios.length() > 0 { ebseq_median(ratios) } else {
      1.0
    }
    size_factors.push(median_val)
  }

  // Normalize so geometric mean of size factors is 1
  let gm = ebseq_geometric_mean(size_factors)
  let norm_sf : Array[Double] = Array::new()
  for sf in size_factors {
    if gm > 0.0 {
      norm_sf.push(sf / gm)
    } else {
      norm_sf.push(1.0)
    }
  }

  EBSeqDataSet::{
    counts: data.counts,
    gene_names: data.gene_names,
    sample_names: data.sample_names,
    conditions: data.conditions,
    size_factors: norm_sf,
    dispersions: data.dispersions,
  }
}

///|
/// Estimate dispersions using method of moments.
/// For NB: Var = mean + dispersion * mean^2
/// => dispersion = (Var - mean) / mean^2
/// Pool estimates across genes and apply shrinkage.
pub fn ebseq_estimate_dispersion(data : EBSeqDataSet) -> EBSeqDataSet {
  let n_genes = data.counts.length()
  let n_samples = data.sample_names.length()

  // Get indices for each condition
  let cond1_idx : Array[Int] = Array::new()
  let cond2_idx : Array[Int] = Array::new()
  for j = 0; j < n_samples; j = j + 1 {
    if data.conditions[j] == 0 {
      cond1_idx.push(j)
    } else {
      cond2_idx.push(j)
    }
  }

  let raw_dispersions : Array[Double] = Array::new()
  let norm_means : Array[Double] = Array::new()

  for i = 0; i < n_genes; i = i + 1 {
    // Compute normalized counts for this gene across all samples
    let norm_counts : Array[Double] = Array::new()
    for j = 0; j < n_samples; j = j + 1 {
      if data.size_factors[j] > 0.0 {
        norm_counts.push(
          data.counts[i][j].to_double() / data.size_factors[j],
        )
      } else {
        norm_counts.push(data.counts[i][j].to_double())
      }
    }

    // Compute mean and variance of normalized counts
    let mut sum_val = 0.0
    for v in norm_counts {
      sum_val = sum_val + v
    }
    let mean_val = sum_val / n_samples.to_double()

    let mut var_sum = 0.0
    for v in norm_counts {
      let diff = v - mean_val
      var_sum = var_sum + diff * diff
    }
    let var_val = var_sum / (n_samples - 1).to_double()

    norm_means.push(if mean_val > 0.0 { @math.ln(mean_val) } else { -10.0 })

    // Method of moments: dispersion = (var - mean) / mean^2
    if mean_val > 0.0 && var_val > mean_val {
      let disp_est = (var_val - mean_val) / (mean_val * mean_val)
      raw_dispersions.push(if disp_est > 0.001 { disp_est } else { 0.001 })
    } else {
      raw_dispersions.push(0.1)
    }
  }

  // Fit dispersion trend: disp ~ a + b / mean
  // Use simple parametric fit for shrinkage
  let alpha_trend = 0.1
  let beta_trend = 0.5

  let final_dispersions : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let raw_disp = raw_dispersions[i]
    let prior_disp = alpha_trend +
      beta_trend * @math.exp(-norm_means[i])
    // Shrinkage: weight raw estimate based on sample size
    let weight = 1.0 / (1.0 + 2.0 / n_samples.to_double())
    let shrunk = (1.0 - weight) * prior_disp + weight * raw_disp
    final_dispersions.push(if shrunk > 0.001 { shrunk } else { 0.001 })
  }

  EBSeqDataSet::{
    counts: data.counts,
    gene_names: data.gene_names,
    sample_names: data.sample_names,
    conditions: data.conditions,
    size_factors: data.size_factors,
    dispersions: final_dispersions,
  }
}

///|
/// Compute marginal likelihood for a gene under the EE (equal expression) model.
/// Both conditions share a common mean q.
fn ebseq_marginal_ee(
  counts_gene : Array[Int],
  cond1_idx : Array[Int],
  cond2_idx : Array[Int],
  size_factors : Array[Double],
  dispersion : Double,
) -> Double {
  let r = 1.0 / dispersion

  // Compute sum of counts/sf and sum of 1/sf for each condition
  let mut sum_c1 = 0.0
  let mut sum_inv_sf_c1 = 0.0
  for j in cond1_idx {
    sum_c1 = sum_c1 + counts_gene[j].to_double() / size_factors[j]
    sum_inv_sf_c1 = sum_inv_sf_c1 + 1.0 / size_factors[j]
  }

  let mut sum_c2 = 0.0
  let mut sum_inv_sf_c2 = 0.0
  for j in cond2_idx {
    sum_c2 = sum_c2 + counts_gene[j].to_double() / size_factors[j]
    sum_inv_sf_c2 = sum_inv_sf_c2 + 1.0 / size_factors[j]
  }

  let total_sum = sum_c1 + sum_c2
  let total_inv_sf = sum_inv_sf_c1 + sum_inv_sf_c2

  if total_sum <= 0.0 {
    return -1.0e300
  }

  // MLE of q under EE: q_hat = total_sum / total_inv_sf
  let q_hat = total_sum / total_inv_sf

  // Compute log-likelihood at MLE
  let mut ll = 0.0
  for j in cond1_idx {
    let mu = size_factors[j] * q_hat
    ll = ll + ebseq_nb_log_pmf(counts_gene[j].to_double(), mu, r)
  }
  for j in cond2_idx {
    let mu = size_factors[j] * q_hat
    ll = ll + ebseq_nb_log_pmf(counts_gene[j].to_double(), mu, r)
  }

  // Laplace approximation: integrate out q
  // log marginal ≈ loglik(q_hat) - 0.5 * log(det(I(q_hat))) + const
  // Fisher information for q: I = sum_j sf_j^2 * r / (mu_j * (r + mu_j))
  let mut fisher = 0.0
  for j in cond1_idx {
    let mu = size_factors[j] * q_hat
    if mu > 0.0 {
      fisher = fisher +
        size_factors[j] * size_factors[j] * r / (mu * (r + mu))
    }
  }
  for j in cond2_idx {
    let mu = size_factors[j] * q_hat
    if mu > 0.0 {
      fisher = fisher +
        size_factors[j] * size_factors[j] * r / (mu * (r + mu))
    }
  }

  if fisher > 0.0 {
    ll - 0.5 * @math.ln(fisher)
  } else {
    ll
  }
}

///|
/// Compute marginal likelihood for a gene under the DE (differential expression) model.
/// Each condition has its own mean.
fn ebseq_marginal_de(
  counts_gene : Array[Int],
  cond1_idx : Array[Int],
  cond2_idx : Array[Int],
  size_factors : Array[Double],
  dispersion : Double,
) -> Double {
  let r = 1.0 / dispersion

  // Condition 1
  let mut sum_c1 = 0.0
  let mut sum_inv_sf_c1 = 0.0
  for j in cond1_idx {
    sum_c1 = sum_c1 + counts_gene[j].to_double() / size_factors[j]
    sum_inv_sf_c1 = sum_inv_sf_c1 + 1.0 / size_factors[j]
  }

  // Condition 2
  let mut sum_c2 = 0.0
  let mut sum_inv_sf_c2 = 0.0
  for j in cond2_idx {
    sum_c2 = sum_c2 + counts_gene[j].to_double() / size_factors[j]
    sum_inv_sf_c2 = sum_inv_sf_c2 + 1.0 / size_factors[j]
  }

  // MLE of q1 and q2
  let q1_hat = if sum_c1 > 0.0 && sum_inv_sf_c1 > 0.0 {
    sum_c1 / sum_inv_sf_c1
  } else {
    0.001
  }
  let q2_hat = if sum_c2 > 0.0 && sum_inv_sf_c2 > 0.0 {
    sum_c2 / sum_inv_sf_c2
  } else {
    0.001
  }

  // Log-likelihood at MLE for each condition
  let mut ll_c1 = 0.0
  for j in cond1_idx {
    let mu = size_factors[j] * q1_hat
    ll_c1 = ll_c1 + ebseq_nb_log_pmf(counts_gene[j].to_double(), mu, r)
  }

  let mut ll_c2 = 0.0
  for j in cond2_idx {
    let mu = size_factors[j] * q2_hat
    ll_c2 = ll_c2 + ebseq_nb_log_pmf(counts_gene[j].to_double(), mu, r)
  }

  let ll = ll_c1 + ll_c2

  // Laplace approximation: Fisher information for (q1, q2)
  let mut fisher1 = 0.0
  for j in cond1_idx {
    let mu = size_factors[j] * q1_hat
    if mu > 0.0 {
      fisher1 = fisher1 +
        size_factors[j] * size_factors[j] * r / (mu * (r + mu))
    }
  }
  let mut fisher2 = 0.0
  for j in cond2_idx {
    let mu = size_factors[j] * q2_hat
    if mu > 0.0 {
      fisher2 = fisher2 +
        size_factors[j] * size_factors[j] * r / (mu * (r + mu))
    }
  }

  // For DE model, we integrate out 2 parameters (q1, q2)
  // log marginal ≈ loglik - 0.5*log(det(I)) + const
  // det(I) = fisher1 * fisher2 (block diagonal since q1, q2 independent)
  if fisher1 > 0.0 && fisher2 > 0.0 {
    ll - 0.5 * @math.ln(fisher1 * fisher2)
  } else {
    ll
  }
}

///|
/// Compute posterior probabilities of DE (PPDE) for all genes.
/// PPDE = P(DE|data) = p_de * L_DE / (p_ee * L_EE + p_de * L_DE)
/// where p_ee and p_de are prior probabilities and L_EE, L_DE are marginal likelihoods.
pub fn ebseq_ppde(data : EBSeqDataSet) -> Array[Double] {
  let n_genes = data.counts.length()
  let n_samples = data.sample_names.length()

  // Get condition indices
  let cond1_idx : Array[Int] = Array::new()
  let cond2_idx : Array[Int] = Array::new()
  for j = 0; j < n_samples; j = j + 1 {
    if data.conditions[j] == 0 {
      cond1_idx.push(j)
    } else {
      cond2_idx.push(j)
    }
  }

  // Compute marginal likelihoods for each gene
  let log_ee : Array[Double] = Array::new()
  let log_de : Array[Double] = Array::new()

  for i = 0; i < n_genes; i = i + 1 {
    let l_ee = ebseq_marginal_ee(
      data.counts[i],
      cond1_idx,
      cond2_idx,
      data.size_factors,
      data.dispersions[i],
    )
    let l_de = ebseq_marginal_de(
      data.counts[i],
      cond1_idx,
      cond2_idx,
      data.size_factors,
      data.dispersions[i],
    )
    log_ee.push(l_ee)
    log_de.push(l_de)
  }

  // Estimate prior probability of DE from data (empirical Bayes)
  // Use the proportion of genes where L_DE > L_EE as initial estimate
  let mut n_de_favored = 0
  for i = 0; i < n_genes; i = i + 1 {
    if log_de[i] > log_ee[i] {
      n_de_favored = n_de_favored + 1
    }
  }
  let p_de_prior = if n_genes > 0 {
    0.5 * n_de_favored.to_double() / n_genes.to_double() + 0.1
  } else {
    0.5
  }
  // Clamp to [0.1, 0.9]
  let p_de = if p_de_prior < 0.1 {
    0.1
  } else if p_de_prior > 0.9 {
    0.9
  } else {
    p_de_prior
  }
  let p_ee = 1.0 - p_de

  // Compute PPDE for each gene
  let ppde : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    // Work in log space to avoid overflow/underflow
    let log_num = @math.ln(p_de) + log_de[i]
    let log_den_ee = @math.ln(p_ee) + log_ee[i]
    let log_den_de = @math.ln(p_de) + log_de[i]

    // log_den = log(exp(log_den_ee) + exp(log_den_de))
    //         = log_den_max + log(1 + exp(log_den_min - log_den_max))
    let log_den_max = if log_den_ee > log_den_de { log_den_ee } else {
      log_den_de
    }
    let log_den_min = if log_den_ee > log_den_de { log_den_de } else {
      log_den_ee
    }
    let diff = log_den_min - log_den_max
    let log_den = log_den_max + @math.ln(1.0 + @math.exp(diff))

    let log_ppde = log_num - log_den
    let ppde_val = @math.exp(log_ppde)

    // Clamp to [0, 1]
    if ppde_val > 1.0 {
      ppde.push(1.0)
    } else if ppde_val < 0.0 {
      ppde.push(0.0)
    } else {
      ppde.push(ppde_val)
    }
  }

  ppde
}

///|
/// Run full EBSeq differential expression test.
/// Returns EBSeqResult with PPDE, log2 fold change, condition means, and DE flags.
pub fn ebseq_test(
  data : EBSeqDataSet,
  threshold? : Double = 0.95,
) -> EBSeqResult {
  let n_genes = data.counts.length()
  let n_samples = data.sample_names.length()

  // Get condition indices
  let cond1_idx : Array[Int] = Array::new()
  let cond2_idx : Array[Int] = Array::new()
  for j = 0; j < n_samples; j = j + 1 {
    if data.conditions[j] == 0 {
      cond1_idx.push(j)
    } else {
      cond2_idx.push(j)
    }
  }

  // Compute PPDE
  let ppde = ebseq_ppde(data)

  // Compute condition means and log2 fold change
  let mean_cond1 : Array[Double] = Array::new()
  let mean_cond2 : Array[Double] = Array::new()
  let log2_fc : Array[Double] = Array::new()

  for i = 0; i < n_genes; i = i + 1 {
    // Mean for condition 1 (normalized)
    let mut sum_c1 = 0.0
    for j in cond1_idx {
      if data.size_factors[j] > 0.0 {
        sum_c1 = sum_c1 +
          data.counts[i][j].to_double() / data.size_factors[j]
      }
    }
    let m1 = if cond1_idx.length() > 0 {
      sum_c1 / cond1_idx.length().to_double()
    } else {
      0.0
    }

    // Mean for condition 2 (normalized)
    let mut sum_c2 = 0.0
    for j in cond2_idx {
      if data.size_factors[j] > 0.0 {
        sum_c2 = sum_c2 +
          data.counts[i][j].to_double() / data.size_factors[j]
      }
    }
    let m2 = if cond2_idx.length() > 0 {
      sum_c2 / cond2_idx.length().to_double()
    } else {
      0.0
    }

    mean_cond1.push(m1)
    mean_cond2.push(m2)

    // Log2 fold change: log2(cond2 / cond1)
    let lfc = if m1 > 0.0 && m2 > 0.0 {
      @math.ln(m2 / m1) / @math.ln(2.0)
    } else if m1 <= 0.0 && m2 > 0.0 {
      10.0
    } else if m1 > 0.0 && m2 <= 0.0 {
      -10.0
    } else {
      0.0
    }
    log2_fc.push(lfc)
  }

  // Flag DE genes
  let flag_de : Array[Bool] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    flag_de.push(ppde[i] >= threshold)
  }

  EBSeqResult::{
    gene_names: data.gene_names,
    ppde,
    log2_fc,
    mean_cond1,
    mean_cond2,
    flag_de,
  }
}

///|
/// Count number of differentially expressed genes in a result.
pub fn ebseq_count_de(result : EBSeqResult) -> Int {
  let mut count = 0
  for flag in result.flag_de {
    if flag {
      count = count + 1
    }
  }
  count
}

///|
/// Get names of differentially expressed genes.
pub fn ebseq_get_de_genes(result : EBSeqResult) -> Array[String] {
  let genes : Array[String] = Array::new()
  for i = 0; i < result.gene_names.length(); i = i + 1 {
    if result.flag_de[i] {
      genes.push(result.gene_names[i])
    }
  }
  genes
}

///|
/// Get top DE genes sorted by PPDE descending.
pub fn ebseq_top_genes(
  result : EBSeqResult,
  n? : Int = 10,
) -> Array[(String, Double, Double)] {
  let n_genes = result.gene_names.length()
  let actual_n = if n < n_genes { n } else { n_genes }

  // Sort indices by PPDE descending (bubble sort)
  let indices : Array[Int] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    indices.push(i)
  }
  for i = 0; i < n_genes - 1; i = i + 1 {
    for j = 0; j < n_genes - 1 - i; j = j + 1 {
      if result.ppde[indices[j]] < result.ppde[indices[j + 1]] {
        let tmp = indices[j]
        indices[j] = indices[j + 1]
        indices[j + 1] = tmp
      }
    }
  }

  let top : Array[(String, Double, Double)] = Array::new()
  for i = 0; i < actual_n; i = i + 1 {
    let idx = indices[i]
    top.push((result.gene_names[idx], result.ppde[idx], result.log2_fc[idx]))
  }
  top
}

///|
/// Generate a sample EBSeq dataset with known differential expression.
/// Genes 0-4 are up-regulated in condition 2, genes 5-9 are down-regulated,
/// and genes 10-19 are equally expressed.
pub fn ebseq_sample_dataset() -> EBSeqDataSet {
  let gene_names = [
    "GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH",
    "GeneI", "GeneJ", "GeneK", "GeneL", "GeneM", "GeneN", "GeneO", "GeneP",
    "GeneQ", "GeneR", "GeneS", "GeneT",
  ]
  let sample_names = [
    "Cond1_Rep1", "Cond1_Rep2", "Cond1_Rep3", "Cond2_Rep1", "Cond2_Rep2",
    "Cond2_Rep3",
  ]
  let conditions = [0, 0, 0, 1, 1, 1]

  let counts : Array[Array[Int]] = Array::new()
  for i = 0; i < 20; i = i + 1 {
    let row : Array[Int] = Array::new()
    let base = 100.0 + i.to_double() * 50.0

    // Determine fold change for condition 2
    let fold_change = if i < 5 {
      3.0 // up-regulated
    } else if i < 10 {
      0.3 // down-regulated
    } else {
      1.0 // equal expression
    }

    // Condition 1 replicates
    for j = 0; j < 3; j = j + 1 {
      let noise = ((i * 7 + j * 13) % 20 - 10).to_double() * 0.1
      row.push((base * (1.0 + noise)).to_int())
    }

    // Condition 2 replicates
    for j = 3; j < 6; j = j + 1 {
      let noise = ((i * 11 + j * 17) % 20 - 10).to_double() * 0.1
      row.push((base * fold_change * (1.0 + noise)).to_int())
    }

    counts.push(row)
  }

  EBSeqDataSet::new(counts, gene_names, sample_names, conditions)
}
