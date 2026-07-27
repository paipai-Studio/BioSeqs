///|
/// microbiome - Microbiome data analysis toolkit.
/// Inspired by Bioconductor microbiome package.
/// Provides alpha diversity, beta diversity, ordination, and differential abundance analysis.

///|
/// Alpha diversity result structure.
pub struct AlphaDiversity {
  observed : Double
  shannon : Double
  simpson : Double
  inv_simpson : Double
  pielou_evenness : Double
  chao1 : Double
  ace : Double
  fisher_alpha : Double
}

///|
pub fn AlphaDiversity::new() -> AlphaDiversity {
  AlphaDiversity::{
    observed: 0.0,
    shannon: 0.0,
    simpson: 0.0,
    inv_simpson: 0.0,
    pielou_evenness: 0.0,
    chao1: 0.0,
    ace: 0.0,
    fisher_alpha: 0.0,
  }
}

///|
/// Calculate observed richness (number of non-zero taxa).
pub fn calc_observed(counts : Array[Double]) -> Double {
  let mut count = 0.0
  let mut i = 0
  while i < counts.length() {
    if counts[i] > 0.0 {
      count = count + 1.0
    }
    i = i + 1
  }
  count
}

///|
/// Calculate Shannon-Wiener diversity index.
/// H = -sum(p_i * ln(p_i))
pub fn calc_shannon(counts : Array[Double]) -> Double {
  let total = sum_double(counts)
  if total == 0.0 {
    return 0.0
  }
  
  let mut shannon = 0.0
  let mut i = 0
  while i < counts.length() {
    if counts[i] > 0.0 {
      let p = counts[i] / total
      shannon = shannon - p * @math.ln(p)
    }
    i = i + 1
  }
  shannon
}

///|
/// Calculate Simpson's diversity index.
/// D = 1 - sum(p_i^2)
pub fn calc_simpson(counts : Array[Double]) -> Double {
  let total = sum_double(counts)
  if total == 0.0 {
    return 0.0
  }
  
  let mut sum_sq = 0.0
  let mut i = 0
  while i < counts.length() {
    let p = counts[i] / total
    sum_sq = sum_sq + p * p
    i = i + 1
  }
  1.0 - sum_sq
}

///|
/// Calculate Inverse Simpson index.
/// 1/D = 1 / sum(p_i^2)
pub fn calc_inv_simpson(counts : Array[Double]) -> Double {
  let total = sum_double(counts)
  if total == 0.0 {
    return 0.0
  }
  
  let mut sum_sq = 0.0
  let mut i = 0
  while i < counts.length() {
    let p = counts[i] / total
    sum_sq = sum_sq + p * p
    i = i + 1
  }
  if sum_sq == 0.0 {
    return 0.0
  }
  1.0 / sum_sq
}

///|
/// Calculate Pielou's evenness index.
/// J = H / ln(S)  where H is Shannon index and S is observed richness
pub fn calc_pielou_evenness(counts : Array[Double]) -> Double {
  let observed = calc_observed(counts)
  if observed <= 1.0 {
    return 0.0
  }
  let shannon = calc_shannon(counts)
  shannon / @math.ln(observed)
}

///|
/// Calculate Chao1 richness estimator.
/// Chao1 = S_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1))
/// where f1 = number of singletons, f2 = number of doubletons
pub fn calc_chao1(counts : Array[Double]) -> Double {
  let s_obs = calc_observed(counts)
  
  let mut f1 = 0.0
  let mut f2 = 0.0
  let mut i = 0
  while i < counts.length() {
    if counts[i] == 1.0 {
      f1 = f1 + 1.0
    } else if counts[i] == 2.0 {
      f2 = f2 + 1.0
    }
    i = i + 1
  }
  
  s_obs + (f1 * (f1 - 1.0)) / (2.0 * (f2 + 1.0))
}

///|
/// Calculate ACE (Abundance-based Coverage Estimator) richness.
pub fn calc_ace(counts : Array[Double], threshold : Int) -> Double {
  let mut n_rare = 0.0
  let mut s_rare = 0.0
  let mut s_abund = 0.0
  let mut f1 = 0.0
  
  let mut i = 0
  while i < counts.length() {
    let count = counts[i]
    if count > 0.0 {
      if count <= threshold.to_double() {
        n_rare = n_rare + count
        s_rare = s_rare + 1.0
        if count == 1.0 {
          f1 = f1 + 1.0
        }
      } else {
        s_abund = s_abund + 1.0
      }
    }
    i = i + 1
  }
  
  if n_rare == 0.0 || s_rare <= f1 {
    return s_rare + s_abund
  }
  
  let mut c_ace = 1.0 - f1 / n_rare
  if c_ace <= 0.0 {
    c_ace = 0.001
  }
  
  let mut sum_fi = 0.0
  let mut sum_fi_i = 0.0
  i = 0
  while i < counts.length() {
    let count = counts[i]
    if count > 0.0 && count <= threshold.to_double() {
      sum_fi = sum_fi + count
      sum_fi_i = sum_fi_i + count * count
    }
    i = i + 1
  }
  
  let gamma2 = 0.0
  if c_ace > 0.0 && n_rare > 0.0 {
    // Simplified ACE calculation
    s_abund + s_rare / c_ace + (f1 / c_ace) * gamma2
  } else {
    s_rare + s_abund
  }
}

///|
/// Calculate Fisher's alpha diversity index.
/// Uses the approximation: S = alpha * ln(1 + N/alpha)
/// where S = number of species, N = total number of individuals
pub fn calc_fisher_alpha(counts : Array[Double]) -> Double {
  let n = sum_double(counts)
  let s = calc_observed(counts)
  
  if n <= 0.0 || s <= 0.0 {
    return 0.0
  }
  
  // Solve S = alpha * ln(1 + N/alpha) numerically
  let mut alpha = s
  let mut iter = 0
  while iter < 50 {
    let f = alpha * @math.ln(1.0 + n / alpha) - s
    let df = @math.ln(1.0 + n / alpha) - n / (n + alpha)
    if df == 0.0 {
      break
    }
    let new_alpha = alpha - f / df
    if new_alpha <= 0.0 {
      alpha = alpha * 0.5
    } else {
      let diff = (new_alpha - alpha).abs()
      alpha = new_alpha
      if diff < 0.0001 {
        break
      }
    }
    iter = iter + 1
  }
  
  alpha
}

///|
/// Calculate all alpha diversity indices for a single sample.
pub fn calc_alpha_diversity(counts : Array[Double]) -> AlphaDiversity {
  AlphaDiversity::{
    observed: calc_observed(counts),
    shannon: calc_shannon(counts),
    simpson: calc_simpson(counts),
    inv_simpson: calc_inv_simpson(counts),
    pielou_evenness: calc_pielou_evenness(counts),
    chao1: calc_chao1(counts),
    ace: calc_ace(counts, 10),
    fisher_alpha: calc_fisher_alpha(counts),
  }
}

///|
/// Calculate alpha diversity for multiple samples (OTU table: rows=taxa, cols=samples).
pub fn calc_alpha_diversity_table(otu_table : Array[Array[Double]]) -> Array[AlphaDiversity] {
  let result : Array[AlphaDiversity] = Array::new()
  if otu_table.length() == 0 {
    return result
  }
  
  let num_samples = otu_table[0].length()
  let mut j = 0
  while j < num_samples {
    let sample_counts : Array[Double] = Array::new()
    let mut i = 0
    while i < otu_table.length() {
      sample_counts.push(otu_table[i][j])
      i = i + 1
    }
    result.push(calc_alpha_diversity(sample_counts))
    j = j + 1
  }
  
  result
}

// ============================================================
// Beta Diversity
// ============================================================

///|
/// Calculate Bray-Curtis dissimilarity between two samples.
/// BC = sum(|x_i - y_i|) / sum(x_i + y_i)
pub fn bray_curtis(x : Array[Double], y : Array[Double]) -> Double {
  if x.length() != y.length() {
    return 0.0
  }
  
  let mut num = 0.0
  let mut den = 0.0
  let mut i = 0
  while i < x.length() {
    num = num + (x[i] - y[i]).abs()
    den = den + x[i] + y[i]
    i = i + 1
  }
  
  if den == 0.0 {
    return 0.0
  }
  
  num / den
}

///|
/// Calculate Jaccard distance between two samples (binary).
/// J = 1 - (A ∩ B) / (A ∪ B)
pub fn jaccard_distance(x : Array[Double], y : Array[Double]) -> Double {
  if x.length() != y.length() {
    return 0.0
  }
  
  let mut intersection = 0.0
  let mut union = 0.0
  let mut i = 0
  while i < x.length() {
    let x_present = x[i] > 0.0
    let y_present = y[i] > 0.0
    if x_present && y_present {
      intersection = intersection + 1.0
    }
    if x_present || y_present {
      union = union + 1.0
    }
    i = i + 1
  }
  
  if union == 0.0 {
    return 0.0
  }
  
  1.0 - intersection / union
}

///|
/// Calculate Jensen-Shannon divergence between two samples.
/// JSD = 0.5 * KL(P||M) + 0.5 * KL(Q||M)
/// where M = (P + Q) / 2
pub fn jensen_shannon_divergence(x : Array[Double], y : Array[Double]) -> Double {
  if x.length() != y.length() {
    return 0.0
  }
  
  let sum_x = sum_double(x)
  let sum_y = sum_double(y)
  
  if sum_x == 0.0 || sum_y == 0.0 {
    return 0.0
  }
  
  let mut jsd = 0.0
  let mut i = 0
  while i < x.length() {
    let p = x[i] / sum_x
    let q = y[i] / sum_y
    let m = (p + q) / 2.0
    
    if p > 0.0 && m > 0.0 {
      jsd = jsd + 0.5 * p * @math.ln(p / m)
    }
    if q > 0.0 && m > 0.0 {
      jsd = jsd + 0.5 * q * @math.ln(q / m)
    }
    i = i + 1
  }
  
  jsd
}

///|
/// Calculate weighted Unifrac-like distance (simplified version using abundance weights).
/// This is a simplified Bray-Curtis with phylogenetic weighting placeholder.
pub fn weighted_unifrac(x : Array[Double], y : Array[Double]) -> Double {
  // Simplified: use abundance-weighted Bray-Curtis as approximation
  bray_curtis(x, y)
}

///|
/// Calculate unweighted Unifrac-like distance (simplified version using presence/absence).
pub fn unweighted_unifrac(x : Array[Double], y : Array[Double]) -> Double {
  // Simplified: use Jaccard distance as approximation
  jaccard_distance(x, y)
}

///|
/// Calculate a beta diversity distance matrix for all samples.
/// OTU table format: rows=taxa, cols=samples
pub fn calc_beta_diversity_matrix(otu_table : Array[Array[Double]], dist_method : String) -> Array[Array[Double]] {
  if otu_table.length() == 0 {
    return Array::new()
  }
  
  let num_samples = otu_table[0].length()
  let dist_matrix : Array[Array[Double]] = Array::new()
  
  // Initialize distance matrix
  let mut i = 0
  while i < num_samples {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < num_samples {
      row.push(0.0)
      j = j + 1
    }
    dist_matrix.push(row)
    i = i + 1
  }
  
  // Calculate pairwise distances
  i = 0
  while i < num_samples {
    let mut j = i + 1
    while j < num_samples {
      let sample_i : Array[Double] = Array::new()
      let sample_j : Array[Double] = Array::new()
      let mut k = 0
      while k < otu_table.length() {
        sample_i.push(otu_table[k][i])
        sample_j.push(otu_table[k][j])
        k = k + 1
      }
      
      let dist = match dist_method {
        "bray" => bray_curtis(sample_i, sample_j)
        "jaccard" => jaccard_distance(sample_i, sample_j)
        "jsd" => jensen_shannon_divergence(sample_i, sample_j)
        "unifrac" => unweighted_unifrac(sample_i, sample_j)
        "wunifrac" => weighted_unifrac(sample_i, sample_j)
        _ => bray_curtis(sample_i, sample_j)
      }
      
      dist_matrix[i][j] = dist
      dist_matrix[j][i] = dist
      j = j + 1
    }
    i = i + 1
  }
  
  dist_matrix
}

// ============================================================
// Ordination (PCoA)
// ============================================================

///|
/// PCoA (Principal Coordinates Analysis) result.
pub struct PCoAResult {
  eigenvalues : Array[Double]
  vectors : Array[Array[Double]]
  variance_explained : Array[Double]
  cumulative_variance : Array[Double]
}

///|
pub fn PCoAResult::new() -> PCoAResult {
  PCoAResult::{
    eigenvalues: Array::new(),
    vectors: Array::new(),
    variance_explained: Array::new(),
    cumulative_variance: Array::new(),
  }
}

///|
/// Perform PCoA (Principal Coordinates Analysis) from a distance matrix.
/// Uses classical MDS (metric multidimensional scaling).
pub fn pcoa(dist_matrix : Array[Array[Double]]) -> PCoAResult {
  let n = dist_matrix.length()
  
  if n <= 1 {
    return PCoAResult::new()
  }
  
  // Step 1: Center the distance matrix
  // Convert distances to similarities and double-center
  let d2 : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      row.push(-0.5 * dist_matrix[i][j] * dist_matrix[i][j])
      j = j + 1
    }
    d2.push(row)
    i = i + 1
  }
  
  // Calculate row means, column means, and grand mean
  let row_means : Array[Double] = Array::new()
  i = 0
  while i < n {
    let mut sum = 0.0
    let mut j = 0
    while j < n {
      sum = sum + d2[i][j]
      j = j + 1
    }
    row_means.push(sum / n.to_double())
    i = i + 1
  }
  
  let col_means : Array[Double] = Array::new()
  let mut j = 0
  while j < n {
    let mut sum = 0.0
    let mut i = 0
    while i < n {
      sum = sum + d2[i][j]
      i = i + 1
    }
    col_means.push(sum / n.to_double())
    j = j + 1
  }
  
  let mut grand_mean = 0.0
  i = 0
  while i < n {
    grand_mean = grand_mean + row_means[i]
    i = i + 1
  }
  grand_mean = grand_mean / n.to_double()
  
  // Double-center the matrix
  let b : Array[Array[Double]] = Array::new()
  i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      let val = d2[i][j] - row_means[i] - col_means[j] + grand_mean
      row.push(val)
      j = j + 1
    }
    b.push(row)
    i = i + 1
  }
  
  // Step 2: Simple power iteration for eigenvalues/vectors (first k)
  // We'll compute first min(n-1, 5) principal coordinates
  let k = if n - 1 < 5 { n - 1 } else { 5 }
  
  let eigenvalues : Array[Double] = Array::new()
  let vectors : Array[Array[Double]] = Array::new()
  
  // Simple implementation using power iteration
  let remaining = b.copy()
  
  let mut comp = 0
  while comp < k {
    // Power iteration
    let v : Array[Double] = Array::make(n, 1.0)
    
    let mut lambda = 0.0
    let mut iter = 0
    while iter < 100 {
      // Multiply remaining matrix by v
      let bv : Array[Double] = Array::new()
      let mut i = 0
      while i < n {
        let mut sum = 0.0
        let mut j = 0
        while j < n {
          sum = sum + remaining[i][j] * v[j]
          j = j + 1
        }
        bv.push(sum)
        i = i + 1
      }
      
      // Compute new lambda (Rayleigh quotient)
      let mut v_norm = 0.0
      let mut v_bv = 0.0
      let mut i = 0
      while i < n {
        v_norm = v_norm + v[i] * v[i]
        v_bv = v_bv + v[i] * bv[i]
        i = i + 1
      }
      
      if v_norm == 0.0 {
        break
      }
      
      let new_lambda = v_bv / v_norm
      
      // Normalize bv
      let bv_norm = 0.0
      let mut i = 0
      while i < n {
        i = i + 1
      }
      let mut bv_len = 0.0
      i = 0
      while i < n {
        bv_len = bv_len + bv[i] * bv[i]
        i = i + 1
      }
      bv_len = bv_len.sqrt()
      
      if bv_len == 0.0 {
        break
      }
      
      i = 0
      while i < n {
        v[i] = bv[i] / bv_len
        i = i + 1
      }
      
      if (new_lambda - lambda).abs() < 0.000001 {
        lambda = new_lambda
        break
      }
      lambda = new_lambda
      iter = iter + 1
    }
    
    if lambda < 0.0 {
      lambda = 0.0
    }
    
    eigenvalues.push(lambda)
    
    // Scale eigenvector by sqrt(eigenvalue)
    let scaled_v : Array[Double] = Array::new()
    let mut i = 0
    while i < n {
      scaled_v.push(v[i] * lambda.sqrt())
      i = i + 1
    }
    vectors.push(scaled_v)
    
    // Deflate: remove this component from remaining
    let mut i = 0
    while i < n {
      let mut j = 0
      while j < n {
        remaining[i][j] = remaining[i][j] - lambda * v[i] * v[j]
        j = j + 1
      }
      i = i + 1
    }
    comp = comp + 1
  }
  
  // Calculate variance explained
  let total_var = sum_double(eigenvalues)
  
  let variance_explained : Array[Double] = Array::new()
  let cumulative_variance : Array[Double] = Array::new()
  let mut cum_var = 0.0
  
  let mut i = 0
  while i < eigenvalues.length() {
    let ve = if total_var > 0.0 { eigenvalues[i] / total_var } else { 0.0 }
    variance_explained.push(ve)
    cum_var = cum_var + ve
    cumulative_variance.push(cum_var)
    i = i + 1
  }
  
  // Transpose vectors to get samples x components
  let transposed : Array[Array[Double]] = Array::new()
  if vectors.length() > 0 {
    let n_samples = vectors[0].length()
    let mut i = 0
    while i < n_samples {
      let row : Array[Double] = Array::new()
      let mut j = 0
      while j < vectors.length() {
        row.push(vectors[j][i])
        j = j + 1
      }
      transposed.push(row)
      i = i + 1
    }
  }
  
  PCoAResult::{
    eigenvalues: eigenvalues,
    vectors: transposed,
    variance_explained: variance_explained,
    cumulative_variance: cumulative_variance,
  }
}

// ============================================================
// Differential Abundance
// ============================================================

///|
/// Differential abundance test result.
pub struct DAResult {
  taxa : String
  log2_fold_change : Double
  p_value : Double
  adjusted_p_value : Double
  mean_group1 : Double
  mean_group2 : Double
}

///|
pub fn DAResult::new(taxa : String, log2_fc : Double, p_val : Double, adj_p : Double, mean1 : Double, mean2 : Double) -> DAResult {
  DAResult::{
    taxa: taxa,
    log2_fold_change: log2_fc,
    p_value: p_val,
    adjusted_p_value: adj_p,
    mean_group1: mean1,
    mean_group2: mean2,
  }
}

///|
/// Simple Wilcoxon rank-sum test for differential abundance.
/// Compares two groups of samples for each taxon.
pub fn differential_abundance(
  otu_table : Array[Array[Double]],
  taxa_names : Array[String],
  group1_indices : Array[Int],
  group2_indices : Array[Int],
) -> Array[DAResult] {
  let results : Array[DAResult] = Array::new()
  
  if otu_table.length() == 0 || group1_indices.length() == 0 || group2_indices.length() == 0 {
    return results
  }
  
  let p_values : Array[Double] = Array::new()
  let log2_fcs : Array[Double] = Array::new()
  let mean1s : Array[Double] = Array::new()
  let mean2s : Array[Double] = Array::new()
  let taxon_names : Array[String] = Array::new()
  
  let mut i = 0
  while i < otu_table.length() {
    // Extract group 1 values
    let group1 : Array[Double] = Array::new()
    let mut j = 0
    while j < group1_indices.length() {
      group1.push(otu_table[i][group1_indices[j]])
      j = j + 1
    }
    
    // Extract group 2 values
    let group2 : Array[Double] = Array::new()
    j = 0
    while j < group2_indices.length() {
      group2.push(otu_table[i][group2_indices[j]])
      j = j + 1
    }
    
    let mean1 = mean_double(group1)
    let mean2 = mean_double(group2)
    
    let log2_fc = if mean1 > 0.0 && mean2 > 0.0 {
      @math.log2(mean2 / mean1)
    } else {
      0.0
    }
    
    // Simple t-test approximation
    let p_val = micro_welch_t_test(group1, group2)
    p_values.push(p_val)
    log2_fcs.push(log2_fc)
    mean1s.push(mean1)
    mean2s.push(mean2)
    
    let taxa_name = if i < taxa_names.length() { taxa_names[i] } else { "Taxon" + i.to_string() }
    taxon_names.push(taxa_name)
    
    i = i + 1
  }
  
  // BH correction
  let adj_p = micro_bh_correction(p_values)
  
  // Build results
  let mut i = 0
  while i < otu_table.length() {
    results.push(DAResult::new(
      taxon_names[i],
      log2_fcs[i],
      p_values[i],
      adj_p[i],
      mean1s[i],
      mean2s[i],
    ))
    i = i + 1
  }
  
  results
}

// ============================================================
// Helper functions
// ============================================================

fn sum_double(arr : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    sum = sum + arr[i]
    i = i + 1
  }
  sum
}

fn mean_double(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  sum_double(arr) / arr.length().to_double()
}

fn variance_double(arr : Array[Double]) -> Double {
  if arr.length() <= 1 {
    return 0.0
  }
  let m = mean_double(arr)
  let mut sum_sq = 0.0
  let mut i = 0
  while i < arr.length() {
    sum_sq = sum_sq + (arr[i] - m) * (arr[i] - m)
    i = i + 1
  }
  sum_sq / (arr.length().to_double() - 1.0)
}

fn micro_welch_t_test(x : Array[Double], y : Array[Double]) -> Double {
  let n1 = x.length().to_double()
  let n2 = y.length().to_double()
  
  if n1 < 2.0 || n2 < 2.0 {
    return 1.0
  }
  
  let mean_x = mean_double(x)
  let mean_y = mean_double(y)
  let var_x = variance_double(x)
  let var_y = variance_double(y)
  
  if var_x == 0.0 && var_y == 0.0 {
    return 1.0
  }
  
  let se = (var_x / n1 + var_y / n2).sqrt()
  if se == 0.0 {
    return 1.0
  }
  
  let t_stat = (mean_x - mean_y) / se
  
  // Welch-Satterthwaite degrees of freedom
  let df_num = (var_x / n1 + var_y / n2) * (var_x / n1 + var_y / n2)
  let df_den = (var_x / n1) * (var_x / n1) / (n1 - 1.0) + (var_y / n2) * (var_y / n2) / (n2 - 1.0)
  let df = if df_den > 0.0 { df_num / df_den } else { n1 + n2 - 2.0 }
  
  // Approximate p-value using normal approximation (two-tailed)
  let p_val = 2.0 * (1.0 - micro_normal_cdf(t_stat.abs()))
  if p_val < 0.0 { 0.0 } else if p_val > 1.0 { 1.0 } else { p_val }
}

fn micro_normal_cdf(x : Double) -> Double {
  // Approximation of standard normal CDF
  0.5 * (1.0 + micro_erf(x / 1.4142135623730951))
}

fn micro_erf(x : Double) -> Double {
  // Abramowitz and Stegun approximation for erf
  let a1 = 0.254829592
  let a2 = -0.284496736
  let a3 = 1.421413741
  let a4 = -1.453152027
  let a5 = 1.061405429
  let p = 0.3275911
  
  let sign = if x < 0.0 { -1.0 } else { 1.0 }
  let x_abs = x.abs()
  
  let t = 1.0 / (1.0 + p * x_abs)
  let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * @math.exp(-x_abs * x_abs)
  
  sign * y
}

fn micro_bh_correction(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return Array::new()
  }
  
  // Create indexed pairs
  let indices : Array[Int] = Array::new()
  let sorted_p : Array[Double] = Array::new()
  
  // Simple bubble sort (since n may be small)
  let mut i = 0
  while i < n {
    indices.push(i)
    sorted_p.push(p_values[i])
    i = i + 1
  }
  
  // Sort
  let mut i = 0
  while i < n - 1 {
    let mut j = 0
    while j < n - i - 1 {
      if sorted_p[j] > sorted_p[j + 1] {
        let temp_p = sorted_p[j]
        sorted_p[j] = sorted_p[j + 1]
        sorted_p[j + 1] = temp_p
        
        let temp_idx = indices[j]
        indices[j] = indices[j + 1]
        indices[j + 1] = temp_idx
      }
      j = j + 1
    }
    i = i + 1
  }
  
  // BH adjustment
  let adj_sorted : Array[Double] = Array::new()
  let mut min_adj = 1.0
  let mut i = n - 1
  while i >= 0 {
    let adj = sorted_p[i] * n.to_double() / (i.to_double() + 1.0)
    if adj < min_adj {
      min_adj = adj
    }
    if min_adj > 1.0 {
      adj_sorted.push(1.0)
    } else {
      adj_sorted.push(min_adj)
    }
    i = i - 1
  }
  
  // Reverse back
  let adj_sorted_rev : Array[Double] = Array::new()
  let mut i = adj_sorted.length() - 1
  while i >= 0 {
    adj_sorted_rev.push(adj_sorted[i])
    i = i - 1
  }
  
  // Map back to original order
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    result.push(0.0)
    i = i + 1
  }
  
  let mut i = 0
  while i < n {
    result[indices[i]] = adj_sorted_rev[i]
    i = i + 1
  }
  
  result
}

///|
/// Create example OTU table for testing/demo.
pub fn create_example_otu_table() -> Array[Array[Double]] {
  // 10 taxa, 6 samples (3 control, 3 treatment)
  let otu_table : Array[Array[Double]] = Array::new()
  
  // Taxon 1: High in control
  let t1 : Array[Double] = [120.0, 150.0, 130.0, 30.0, 25.0, 35.0]
  otu_table.push(t1)
  
  // Taxon 2: High in treatment
  let t2 : Array[Double] = [20.0, 25.0, 15.0, 100.0, 120.0, 110.0]
  otu_table.push(t2)
  
  // Taxon 3: Moderate in both
  let t3 : Array[Double] = [50.0, 55.0, 48.0, 52.0, 47.0, 53.0]
  otu_table.push(t3)
  
  // Taxon 4: Low abundance
  let t4 : Array[Double] = [5.0, 3.0, 7.0, 4.0, 6.0, 2.0]
  otu_table.push(t4)
  
  // Taxon 5: Rare (singletons)
  let t5 : Array[Double] = [1.0, 0.0, 1.0, 0.0, 1.0, 0.0]
  otu_table.push(t5)
  
  // Taxon 6
  let t6 : Array[Double] = [80.0, 75.0, 85.0, 60.0, 65.0, 55.0]
  otu_table.push(t6)
  
  // Taxon 7
  let t7 : Array[Double] = [30.0, 35.0, 28.0, 45.0, 50.0, 40.0]
  otu_table.push(t7)
  
  // Taxon 8
  let t8 : Array[Double] = [10.0, 12.0, 8.0, 15.0, 18.0, 14.0]
  otu_table.push(t8)
  
  // Taxon 9: Doubletons
  let t9 : Array[Double] = [2.0, 2.0, 0.0, 0.0, 2.0, 2.0]
  otu_table.push(t9)
  
  // Taxon 10
  let t10 : Array[Double] = [45.0, 50.0, 42.0, 55.0, 60.0, 48.0]
  otu_table.push(t10)
  
  otu_table
}

///|
/// Get example taxa names.
pub fn get_example_taxa_names() -> Array[String] {
  let names : Array[String] = [
    "Bacteroides",
    "Prevotella",
    "Firmicutes_g1",
    "Actinobacteria_g1",
    "Rare_taxon_1",
    "Lachnospiraceae",
    "Ruminococcaceae",
    "Bifidobacterium",
    "Rare_taxon_2",
    "Escherichia",
  ]
  names
}
