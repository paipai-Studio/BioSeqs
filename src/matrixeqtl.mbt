///|
/// MatrixEQTL - Fast eQTL Analysis for Large-Scale Data.
/// Inspired by Bioconductor MatrixEQTL package.
///
/// Performs efficient SNP-gene association testing using linear regression
/// for all SNP-gene pairs, with Benjamini-Hochberg FDR correction.

///|
/// Input data for MatrixEQTL analysis.
/// snps: SNP genotype matrix (snps x samples), each row is a SNP
/// gene_expression: Gene expression matrix (genes x samples), each row is a gene
/// snp_names: Names of SNPs (one per row in snps matrix)
/// gene_names: Names of genes (one per row in gene_expression matrix)
/// sample_names: Names of samples (one per column)
pub struct MatrixEQTLInput {
  snps : Array[Array[Double]]
  gene_expression : Array[Array[Double]]
  snp_names : Array[String]
  gene_names : Array[String]
  sample_names : Array[String]
}

///|
/// Create a MatrixEQTLInput struct.
pub fn MatrixEQTLInput::new(
  snps : Array[Array[Double]],
  gene_expression : Array[Array[Double]],
  snp_names : Array[String],
  gene_names : Array[String],
  sample_names : Array[String],
) -> MatrixEQTLInput {
  MatrixEQTLInput::{
    snps,
    gene_expression,
    snp_names,
    gene_names,
    sample_names,
  }
}

///|
/// Result for a single SNP-gene association test.
pub struct MatrixEQTLResult {
  snp_name : String
  gene_name : String
  beta : Double
  t_stat : Double
  p_value : Double
  fdr : Double
}

///|
/// Create a MatrixEQTLResult struct.
pub fn MatrixEQTLResult::new(
  snp_name : String,
  gene_name : String,
  beta : Double,
  t_stat : Double,
  p_value : Double,
  fdr : Double,
) -> MatrixEQTLResult {
  MatrixEQTLResult::{ snp_name, gene_name, beta, t_stat, p_value, fdr, }
}

///|
/// Output of MatrixEQTL analysis containing all significant results.
pub struct MatrixEQTLOutput {
  results : Array[MatrixEQTLResult]
  n_tests : Int
}

///|
/// Create a MatrixEQTLOutput struct.
pub fn MatrixEQTLOutput::new(
  results : Array[MatrixEQTLResult],
  n_tests : Int,
) -> MatrixEQTLOutput {
  MatrixEQTLOutput::{ results, n_tests, }
}

///|
/// Compute mean of a Double array.
fn meqtl_mean(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for v in data {
    sum = sum + v
  }
  sum / data.length().to_double()
}



///|
/// Approximate cumulative distribution function of the standard normal.
fn meqtl_pnorm(x : Double) -> Double {
  let t = 1.0 / (1.0 + 0.2316419 * x.abs())
  let d = 0.3989423 * @math.exp(-x * x / 2.0)
  let prob = d *
    t *
    (
      0.3193815 +
      t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274)))
    )
  if x >= 0.0 {
    1.0 - prob
  } else {
    prob
  }
}

///|
/// Approximate the two-tailed p-value from a t-statistic with given degrees of freedom.
/// Uses the normal approximation for large df, and a refined approximation for small df.
fn meqtl_t_dist_p_value(t_stat : Double, df : Double) -> Double {
  if df <= 0.0 {
    return 1.0
  }

  // For large df, t-distribution approaches normal
  if df > 30.0 {
    let z = t_stat.abs()
    return 2.0 * (1.0 - meqtl_pnorm(z))
  }

  // For small df, use the relationship: t^2 ~ F(1, df)
  // P(|T| > |t|) = I(df/(df+t^2), df/2, 1/2) using regularized incomplete beta
  let x = df / (df + t_stat * t_stat)
  let a = df / 2.0
  let b = 0.5
  let p = meqtl_incomplete_beta(x, a, b)
  if p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Regularized incomplete beta function I_x(a, b).
/// Used to compute t-distribution p-values.
fn meqtl_incomplete_beta(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }

  let lbeta_val = meqtl_lngamma(a) + meqtl_lngamma(b) - meqtl_lngamma(a + b)
  let mut bt = @math.ln(x) * a + @math.ln(1.0 - x) * b - lbeta_val
  bt = @math.exp(bt)

  if x < (a + 1.0) / (a + b + 2.0) {
    bt * meqtl_betacf(x, a, b) / a
  } else {
    1.0 - bt * meqtl_betacf(1.0 - x, b, a) / b
  }
}

///|
/// Continued fraction representation for the incomplete beta function.
fn meqtl_betacf(x : Double, a : Double, b : Double) -> Double {
  let max_iter = 200
  let eps = 0.0000003
  let fpmin = 1.0e-30

  let mut c = 1.0
  let mut d = 1.0 - (a + b) * x / (a + 1.0)
  if d.abs() < fpmin {
    d = fpmin
  }
  d = 1.0 / d
  let mut h = d

  let mut m = 1
  while m <= max_iter {
    let m2 = (2 * m).to_double()
    let mut aa = m.to_double() *
      (b - m.to_double()) *
      x /
      ((a + m2 - 1.0) * (a + m2))
    d = 1.0 + aa * d
    if d.abs() < fpmin {
      d = fpmin
    }
    c = 1.0 + aa / c
    if c.abs() < fpmin {
      c = fpmin
    }
    d = 1.0 / d
    let del = d * c
    h = h * del

    aa = -(a + m.to_double()) *
      (a + b + m.to_double()) *
      x /
      ((a + m2) * (a + m2 + 1.0))
    d = 1.0 + aa * d
    if d.abs() < fpmin {
      d = fpmin
    }
    c = 1.0 + aa / c
    if c.abs() < fpmin {
      c = fpmin
    }
    d = 1.0 / d
    let del2 = d * c
    h = h * del2

    if (del2 - 1.0).abs() < eps {
      return h
    }
    m = m + 1
  }
  h
}

///|
/// Log of the gamma function using the Lanczos approximation.
fn meqtl_lngamma(x : Double) -> Double {
  let cof = [
    76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
    0.001208650973866179, -0.000005395239384953,
  ]
  let mut y = x
  let mut tmp = x + 5.5
  tmp = (x + 0.5) * @math.ln(tmp) - tmp

  let mut ser = 1.000000000190015
  let mut j = 0
  while j < cof.length() {
    y = y + 1.0
    ser = ser + cof[j] / y
    j = j + 1
  }

  tmp + @math.ln(2.5066282746310005 * ser / x)
}

///|
/// Perform simple linear regression: y = beta0 + beta1 * x + error.
/// Returns (beta1, t_statistic, p_value).
///
/// Computation:
///   beta1 = cov(x,y) / var(x)
///   beta0 = mean(y) - beta1 * mean(x)
///   residual_se = sqrt(sum(residuals^2) / (n-2))
///   se_beta1 = residual_se / sqrt(sum((x - mean_x)^2))
///   t = beta1 / se_beta1
///   p-value from t-distribution with n-2 degrees of freedom
pub fn matrixeqtl_linear_regression(
  x : Array[Double],
  y : Array[Double],
) -> (Double, Double, Double) {
  let n = x.length()
  if n != y.length() || n < 3 {
    return (0.0, 0.0, 1.0)
  }

  let mean_x = meqtl_mean(x)
  let mean_y = meqtl_mean(y)

  // Compute sum of squared deviations for x
  let mut ssx = 0.0
  let mut i = 0
  while i < n {
    let diff = x[i] - mean_x
    ssx = ssx + diff * diff
    i = i + 1
  }

  // If x is constant, no regression possible
  if ssx < 1.0e-15 {
    return (0.0, 0.0, 1.0)
  }

  // Compute sum of cross-deviations
  let mut sxy = 0.0
  i = 0
  while i < n {
    sxy = sxy + (x[i] - mean_x) * (y[i] - mean_y)
    i = i + 1
  }

  let beta1 = sxy / ssx
  let beta0 = mean_y - beta1 * mean_x

  // Compute residual sum of squares
  let mut rss = 0.0
  i = 0
  while i < n {
    let residual = y[i] - (beta0 + beta1 * x[i])
    rss = rss + residual * residual
    i = i + 1
  }

  let df = (n - 2).to_double()
  let residual_se = (rss / df).sqrt()
  let se_beta1 = residual_se / ssx.sqrt()

  if se_beta1 < 1.0e-15 {
    return (beta1, 0.0, 1.0)
  }

  let t_stat = beta1 / se_beta1
  let p_value = meqtl_t_dist_p_value(t_stat, df)

  (beta1, t_stat, p_value)
}

///|
/// Benjamini-Hochberg FDR correction for multiple testing.
/// Adjusts p-values to control the false discovery rate.
///
/// Algorithm:
///   1. Sort p-values in ascending order
///   2. Adjust: p_adj[i] = p[i] * m / rank[i]
///   3. Enforce monotonicity (cumulative minimum from largest rank)
///   4. Cap at 1.0
pub fn matrixeqtl_bh_fdr(p_values : Array[Double]) -> Array[Double] {
  let m = p_values.length()
  if m == 0 {
    return []
  }

  // Create index array and sort by p-value
  let sorted_indices : Array[Int] = Array(capacity=0)
  for i = 0; i < m; i = i + 1 {
    sorted_indices.push(i)
  }

  // Bubble sort indices by p-value
  for i = 0; i < m - 1; i = i + 1 {
    for j = 0; j < m - 1 - i; j = j + 1 {
      if p_values[sorted_indices[j]] > p_values[sorted_indices[j + 1]] {
        let tmp = sorted_indices[j]
        sorted_indices[j] = sorted_indices[j + 1]
        sorted_indices[j + 1] = tmp
      }
    }
  }

  // Compute adjusted p-values with monotonicity enforcement
  let padj : Array[Double] = Array::make(m, 1.0)

  let mut min_padj = 1.0
  let mut i = m - 1
  while i >= 0 {
    let rank = (i + 1).to_double()
    let adj = p_values[sorted_indices[i]] * m.to_double() / rank
    if adj < min_padj {
      min_padj = adj
    }
    padj[sorted_indices[i]] = min_padj
    if i == 0 {
      break
    }
    i = i - 1
  }

  // Cap at 1.0
  let final_padj : Array[Double] = Array(capacity=0)
  for p in padj {
    final_padj.push(if p > 1.0 { 1.0 } else { p })
  }
  final_padj
}

///|
/// Run full MatrixEQTL analysis: test all SNP-gene pairs.
///
/// For each SNP-gene pair, fits the linear model:
///   gene_expression = beta0 + beta1 * SNP + error
///
/// Only results with p-value below pv_threshold are retained.
/// FDR correction is applied to all tests.
pub fn matrixeqtl_run(
  input : MatrixEQTLInput,
  pv_threshold? : Double = 0.05,
) -> MatrixEQTLOutput {
  let n_snps = input.snps.length()
  let n_genes = input.gene_expression.length()

  if n_snps == 0 || n_genes == 0 {
    return MatrixEQTLOutput::new([], 0)
  }

  let n_samples = input.sample_names.length()
  if n_samples < 3 {
    return MatrixEQTLOutput::new([], 0)
  }

  // Collect all SNP-gene test results
  let all_snp_names : Array[String] = Array(capacity=0)
  let all_gene_names : Array[String] = Array(capacity=0)
  let all_betas : Array[Double] = Array(capacity=0)
  let all_t_stats : Array[Double] = Array(capacity=0)
  let all_p_values : Array[Double] = Array(capacity=0)

  // Pre-compute SNP statistics for efficiency
  // For each SNP, pre-compute mean, ssx (sum of squared deviations)
  let snp_means : Array[Double] = Array::make(n_snps, 0.0)
  let snp_ssx : Array[Double] = Array::make(n_snps, 0.0)

  let mut s = 0
  while s < n_snps {
    let snp_data = input.snps[s]
    let mean_s = meqtl_mean(snp_data)
    snp_means[s] = mean_s
    let mut ssx = 0.0
    let mut i = 0
    while i < n_samples {
      let diff = snp_data[i] - mean_s
      ssx = ssx + diff * diff
      i = i + 1
    }
    snp_ssx[s] = ssx
    s = s + 1
  }

  // Pre-compute gene statistics
  let gene_means : Array[Double] = Array::make(n_genes, 0.0)

  let mut g = 0
  while g < n_genes {
    gene_means[g] = meqtl_mean(input.gene_expression[g])
    g = g + 1
  }

  // Test each SNP-gene pair
  let mut n_tests = 0
  s = 0
  while s < n_snps {
    let snp_data = input.snps[s]
    let mean_x = snp_means[s]
    let ssx = snp_ssx[s]

    // Skip constant SNPs (no variation)
    if ssx < 1.0e-15 {
      // Record trivial results for constant SNPs
      g = 0
      while g < n_genes {
        n_tests = n_tests + 1
        all_snp_names.push(input.snp_names[s])
        all_gene_names.push(input.gene_names[g])
        all_betas.push(0.0)
        all_t_stats.push(0.0)
        all_p_values.push(1.0)
        g = g + 1
      }
      s = s + 1
      continue
    }

    g = 0
    while g < n_genes {
      n_tests = n_tests + 1
      let gene_data = input.gene_expression[g]
      let mean_y = gene_means[g]

      // Compute cross-deviation sum
      let mut sxy = 0.0
      let mut i = 0
      while i < n_samples {
        sxy = sxy + (snp_data[i] - mean_x) * (gene_data[i] - mean_y)
        i = i + 1
      }

      let beta1 = sxy / ssx
      let beta0 = mean_y - beta1 * mean_x

      // Compute residual sum of squares
      let mut rss = 0.0
      i = 0
      while i < n_samples {
        let residual = gene_data[i] - (beta0 + beta1 * snp_data[i])
        rss = rss + residual * residual
        i = i + 1
      }

      let df = (n_samples - 2).to_double()
      let residual_se = (rss / df).sqrt()
      let se_beta1 = residual_se / ssx.sqrt()

      let (t_stat, p_value) = if se_beta1 < 1.0e-15 {
        (0.0, 1.0)
      } else {
        let t = beta1 / se_beta1
        let p = meqtl_t_dist_p_value(t, df)
        (t, p)
      }

      all_snp_names.push(input.snp_names[s])
      all_gene_names.push(input.gene_names[g])
      all_betas.push(beta1)
      all_t_stats.push(t_stat)
      all_p_values.push(p_value)

      g = g + 1
    }

    s = s + 1
  }

  // Apply BH FDR correction to all p-values
  let fdr_values = matrixeqtl_bh_fdr(all_p_values)

  // Filter results by p-value threshold
  let results : Array[MatrixEQTLResult] = Array(capacity=0)
  let mut i = 0
  while i < all_p_values.length() {
    if all_p_values[i] <= pv_threshold {
      results.push(
        MatrixEQTLResult::new(
          all_snp_names[i],
          all_gene_names[i],
          all_betas[i],
          all_t_stats[i],
          all_p_values[i],
          fdr_values[i],
        ),
      )
    }
    i = i + 1
  }

  // Sort results by p-value (ascending)
  let n_results = results.length()
  for a = 0; a < n_results - 1; a = a + 1 {
    for b = 0; b < n_results - 1 - a; b = b + 1 {
      if results[b].p_value > results[b + 1].p_value {
        let tmp = results[b]
        results[b] = results[b + 1]
        results[b + 1] = tmp
      }
    }
  }

  MatrixEQTLOutput::new(results, n_tests)
}

///|
/// Create a sample MatrixEQTL input dataset for testing.
pub fn sample_matrixeqtl_input() -> MatrixEQTLInput {
  let snp_names = ["SNP1", "SNP2", "SNP3", "SNP4", "SNP5"]
  let gene_names = ["GeneA", "GeneB", "GeneC", "GeneD"]
  let sample_names = [
    "Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6", "Sample7",
    "Sample8",
  ]

  // SNP genotype data (0, 1, 2 for genotype dosage)
  let snps : Array[Array[Double]] = [
    [0.0, 1.0, 1.0, 2.0, 0.0, 1.0, 2.0, 1.0],
    [1.0, 1.0, 0.0, 0.0, 2.0, 1.0, 1.0, 0.0],
    [0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2.0],
    [2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0],
    [1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0],
  ]

  // Gene expression data with some associations to SNPs
  let gene_expression : Array[Array[Double]] = [
    // GeneA correlated with SNP1
    [5.1, 6.3, 6.5, 7.8, 4.9, 6.1, 7.9, 6.4],
    // GeneB correlated with SNP2
    [10.5, 10.2, 8.1, 8.3, 12.7, 10.4, 10.1, 8.0],
    // GeneC no strong correlation
    [3.2, 3.5, 3.1, 3.8, 3.0, 3.4, 3.6, 3.3],
    // GeneD correlated with SNP4
    [12.0, 10.5, 8.2, 10.3, 11.8, 10.6, 8.0, 10.4],
  ]

  MatrixEQTLInput::new(
    snps, gene_expression, snp_names, gene_names, sample_names,
  )
}
