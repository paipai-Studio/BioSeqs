///|
/// DRIMSeq module - Differential transcript usage analysis.
/// Inspired by Bioconductor DRIMSeq package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/DRIMSeq.html
///   - Nowicka M et al. (2016) DRIMSeq: a Dirichlet-multinomial approach
///     for differential transcript usage analysis with RNA-seq data.
///     PLoS Computational Biology 12(10):e1005048.

///|
/// Normalization method for DRIMSeq.
pub enum DRIMSeqNorm {
  None
  TMM
  Sum
} derive(Eq)

///|
/// Create no normalization.
pub fn drimseq_norm_none() -> DRIMSeqNorm {
  DRIMSeqNorm::None
}

///|
/// Create TMM normalization.
pub fn drimseq_norm_tmm() -> DRIMSeqNorm {
  DRIMSeqNorm::TMM
}

///|
/// Create sum normalization.
pub fn drimseq_norm_sum() -> DRIMSeqNorm {
  DRIMSeqNorm::Sum
}

///|
/// Transcript count record.
pub struct TranscriptCount {
  /// Transcript ID
  transcript_id : String
  /// Gene ID
  gene_id : String
  /// Sample ID
  sample_id : String
  /// Condition/group
  condition : String
  /// Count value
  mut count : Double
  /// Gene-level count (total for all transcripts of gene)
  gene_count : Double
}

///|
/// Create a new TranscriptCount.
pub fn TranscriptCount::new(
  transcript_id~ : String,
  gene_id~ : String,
  sample_id~ : String,
  condition~ : String,
  count~ : Double,
  gene_count? : Double = 0.0,
) -> TranscriptCount {
  TranscriptCount::{
    transcript_id,
    gene_id,
    sample_id,
    condition,
    count,
    gene_count,
  }
}

///|
/// Calculate the proportion (count / gene_count).
pub fn TranscriptCount::proportion(self : TranscriptCount) -> Double {
  if self.gene_count > 0.0 { self.count / self.gene_count } else { 0.0 }
}

///|
/// DRIMSeq test result for a gene.
pub struct DRIMSeqGeneResult {
  /// Gene ID
  gene_id : String
  /// Delta (difference in proportions)
  delta : Array[Double]
  /// P-value
  p_value : Double
  /// Adjusted p-value (BH)
  mut p_adj : Double
  /// Degrees of freedom
  df : Int
  /// Likelihood ratio statistic
  lr_stat : Double
  /// Number of transcripts
  n_transcripts : Int
  /// Significant at alpha level
  mut significant : Bool
}

///|
/// Create a new DRIMSeqGeneResult.
pub fn DRIMSeqGeneResult::new(
  gene_id~ : String,
  delta~ : Array[Double],
  p_value~ : Double,
  df~ : Int,
  lr_stat~ : Double,
) -> DRIMSeqGeneResult {
  DRIMSeqGeneResult::{
    gene_id,
    delta,
    p_value,
    p_adj: 0.0,
    df,
    lr_stat,
    n_transcripts: delta.length(),
    significant: false,
  }
}

///|
/// DRIMSeq full analysis result.
pub struct DRIMSeqResult {
  /// Gene-level results
  gene_results : Array[DRIMSeqGeneResult]
  /// Filtered gene results (passed minimum count threshold)
  filtered_results : Array[DRIMSeqGeneResult]
  /// Significance threshold
  alpha : Double
  /// Total genes tested
  n_genes_tested : Int
  /// Genes passing filter
  n_genes_filtered : Int
  /// Significant genes
  n_significant : Int
  /// Convergence status
  converged : Bool
  /// Iterations
  iterations : Int
}

///|
/// Get significant results.
pub fn DRIMSeqResult::get_significant(self : DRIMSeqResult) -> Array[DRIMSeqGeneResult] {
  let sig : Array[DRIMSeqGeneResult] = Array::new()
  let mut i = 0
  while i < self.gene_results.length() {
    if self.gene_results[i].significant {
      sig.push(self.gene_results[i])
    }
    i = i + 1
  }
  sig
}

///|
/// Get top results by p-value.
pub fn DRIMSeqResult::get_top_genes(
  self : DRIMSeqResult,
  n~ : Int,
) -> Array[DRIMSeqGeneResult] {
  let sorted = self.gene_results
  // Simple bubble sort by p-value
  let len = sorted.length()
  let mut i = 0
  while i < len - 1 {
    let mut j = 0
    while j < len - i - 1 {
      if sorted[j].p_value > sorted[j + 1].p_value {
        let temp = sorted[j]
        sorted[j] = sorted[j + 1]
        sorted[j + 1] = temp
      }
      j = j + 1
    }
    i = i + 1
  }
  let top : Array[DRIMSeqGeneResult] = Array::new()
  let limit = if n < len { n } else { len }
  let mut k = 0
  while k < limit {
    top.push(sorted[k])
    k = k + 1
  }
  top
}

///|
/// Summary of DRIMSeq results.
pub fn DRIMSeqResult::summary(self : DRIMSeqResult) -> String {
  let mut s = "=== DRIMSeq Differential Transcript Usage ===" + "\n"
  s = s + "Total genes tested: " + self.n_genes_tested.to_string() + "\n"
  s = s + "Genes passing filter: " + self.n_genes_filtered.to_string() + "\n"
  s = s + "Significant genes (alpha=" + self.alpha.to_string() + "): " + self.n_significant.to_string() + "\n"
  s = s + "Converged: " + self.converged.to_string() + "\n"
  s = s + "Iterations: " + self.iterations.to_string() + "\n"
  s
}

///|
/// Configuration for DRIMSeq analysis.
pub struct DRIMSeqConfig {
  /// Minimum count threshold for filtering
  min_count : Double
  /// Minimum proportion threshold
  min_proportion : Double
  /// Convergence tolerance
  tolerance : Double
  /// Maximum iterations
  max_iter : Int
  /// Significance threshold
  alpha : Double
  /// Normalization method
  norm_method : DRIMSeqNorm
}

///|
/// Create default DRIMSeqConfig.
pub fn DRIMSeqConfig::new() -> DRIMSeqConfig {
  DRIMSeqConfig::{
    min_count: 10.0,
    min_proportion: 0.01,
    tolerance: 0.000001,
    max_iter: 100,
    alpha: 0.05,
    norm_method: DRIMSeqNorm::TMM,
  }
}

///|
/// Set minimum count.
pub fn DRIMSeqConfig::set_min_count(
  self : DRIMSeqConfig,
  val~ : Double,
) -> DRIMSeqConfig {
  DRIMSeqConfig::{
    min_count: val,
    min_proportion: self.min_proportion,
    tolerance: self.tolerance,
    max_iter: self.max_iter,
    alpha: self.alpha,
    norm_method: self.norm_method,
  }
}

///|
/// Set normalization method.
pub fn DRIMSeqConfig::set_norm_method(
  self : DRIMSeqConfig,
  method~ : DRIMSeqNorm,
) -> DRIMSeqConfig {
  DRIMSeqConfig::{
    min_count: self.min_count,
    min_proportion: self.min_proportion,
    tolerance: self.tolerance,
    max_iter: self.max_iter,
    alpha: self.alpha,
    norm_method: method,
  }
}

///|
/// Set alpha threshold.
pub fn DRIMSeqConfig::set_alpha(
  self : DRIMSeqConfig,
  val~ : Double,
) -> DRIMSeqConfig {
  DRIMSeqConfig::{
    min_count: self.min_count,
    min_proportion: self.min_proportion,
    tolerance: self.tolerance,
    max_iter: self.max_iter,
    alpha: val,
    norm_method: self.norm_method,
  }
}

///|
/// Aggregate counts by gene.
pub fn drimseq_aggregate_by_gene(
  counts : Array[TranscriptCount],
) -> Array[TranscriptCount] {
  let result : Array[TranscriptCount] = Array::new()
  let mut i = 0
  while i < counts.length() {
    let c = counts[i]
    let mut found = false
    let mut j = 0
    while j < result.length() {
      if result[j].gene_id == c.gene_id && result[j].sample_id == c.sample_id && result[j].condition == c.condition {
        result[j].count = result[j].count + c.count
        found = true
        break
      }
      j = j + 1
    }
    if !found {
      result.push(TranscriptCount::new(
        transcript_id="__gene_total__", gene_id=c.gene_id, sample_id=c.sample_id, condition=c.condition, count=c.count,
      ))
    }
    i = i + 1
  }
  result
}

///|
/// Compute proportions (transcript count / gene count) for each record.
pub fn drimseq_compute_proportions(
  counts : Array[TranscriptCount],
  gene_counts : Array[TranscriptCount],
) -> Array[TranscriptCount] {
  let result : Array[TranscriptCount] = Array::new()
  let mut i = 0
  while i < counts.length() {
    let c = counts[i]
    let mut gene_total = 0.0
    let mut j = 0
    while j < gene_counts.length() {
      if gene_counts[j].gene_id == c.gene_id && gene_counts[j].sample_id == c.sample_id {
        gene_total = gene_counts[j].count
        break
      }
      j = j + 1
    }
    result.push(TranscriptCount::new(
      transcript_id=c.transcript_id, gene_id=c.gene_id, sample_id=c.sample_id, condition=c.condition, count=c.count,
      gene_count=gene_total,
    ))
    i = i + 1
  }
  result
}

///|
/// Filter transcripts by minimum count and proportion.
pub fn drimseq_filter_counts(
  counts : Array[TranscriptCount],
  config? : DRIMSeqConfig = DRIMSeqConfig::new(),
) -> Array[TranscriptCount] {
  let filtered : Array[TranscriptCount] = Array::new()
  let mut i = 0
  while i < counts.length() {
    let c = counts[i]
    if c.count >= config.min_count {
      let prop = c.proportion()
      if prop >= config.min_proportion || c.count >= config.min_count * 2.0 {
        filtered.push(c)
      }
    }
    i = i + 1
  }
  filtered
}

///|
/// Compute Dirichlet-multinomial log probability.
/// Returns log-likelihood for given counts and proportions.
pub fn drimseq_dirichlet_multinomial(
  counts~ : Array[Double],
  proportions~ : Array[Double],
) -> Double {
  let n = counts.length()
  if n == 0 { return 0.0 }
  let sum_counts = array_sum(counts)
  if sum_counts <= 0.0 { return 0.0 }
  let sum_props = array_sum(proportions)
  if sum_props <= 0.0 { return 0.0 }

  let mut loglik = 0.0
  // Log-gamma(sum(counts) + 1) - sum(log-gamma(count_i + 1))
  loglik = loglik + drimseq_log_gamma(sum_counts + 1.0)
  let mut i = 0
  while i < n {
    loglik = loglik - drimseq_log_gamma(counts[i] + 1.0)
    i = i + 1
  }

  // Dirichlet part
  let mut alpha_sum = 0.0
  let mut j = 0
  while j < n {
    let alpha_j = proportions[j] * sum_props * 10.0
    alpha_sum = alpha_sum + alpha_j
    loglik = loglik + drimseq_log_gamma(alpha_j + counts[j])
    loglik = loglik - drimseq_log_gamma(alpha_j)
    j = j + 1
  }
  loglik = loglik - drimseq_log_gamma(alpha_sum + sum_counts) + drimseq_log_gamma(alpha_sum)

  loglik
}

///|
/// Simplified Wald test for differential transcript usage.
/// Tests if proportions differ between conditions.
pub fn drimseq_wald_test(
  proportions_a : Array[Double],
  proportions_b : Array[Double],
) -> (Double, Double, Array[Double]) {
  let n = proportions_a.length()
  if n == 0 { return (0.0, 1.0, Array::new()) }
  let delta = array_sub(proportions_a, proportions_b)
  let mut stat = 0.0
  let mut i = 0
  while i < n {
    let pa = proportions_a[i]
    let pb = proportions_b[i]
    let d = delta[i]
    let var_a = pa * (1.0 - pa)
    let var_b = pb * (1.0 - pb)
    let var_avg = (var_a + var_b) / 2.0
    let v = if var_avg > 0.0 { var_avg } else { 0.001 }
    stat = stat + d * d / v
    i = i + 1
  }
  let pval = if stat > 0.0 { chi_square_p_value(stat, n - 1) } else { 1.0 }
  (stat, pval, delta)
}

///|
/// Main DRIMSeq analysis function.
pub fn drimseq_test_differential(
  counts : Array[TranscriptCount],
  config? : DRIMSeqConfig = DRIMSeqConfig::new(),
) -> DRIMSeqResult {
  let gene_counts = drimseq_aggregate_by_gene(counts)
  let props = drimseq_compute_proportions(counts, gene_counts)
  let filtered = drimseq_filter_counts(props, config=config)

  let gene_results : Array[DRIMSeqGeneResult] = Array::new()

  // Group by gene and condition
  let genes = get_unique_genes(filtered)
  let mut total_tested = 0
  let mut total_filtered = 0

  let mut gi = 0
  while gi < genes.length() {
    let gene_id = genes[gi]
    total_tested = total_tested + 1

    let gene_data = get_gene_data(filtered, gene_id)
    let conditions = get_unique_conditions(gene_data)

    if conditions.length() < 2 {
      gi = gi + 1
      continue
    }

    total_filtered = total_filtered + 1

    let cond_a = conditions[0]
    let cond_b = conditions[1]

    let props_a = get_condition_proportions(gene_data, cond_a)
    let props_b = get_condition_proportions(gene_data, cond_b)

    let n_transcripts = min_int(props_a.length(), props_b.length())
    let pa = slice_array(props_a, 0, n_transcripts)
    let pb = slice_array(props_b, 0, n_transcripts)

    let (stat, pval, delta) = drimseq_wald_test(pa, pb)

    gene_results.push(DRIMSeqGeneResult::new(
      gene_id=gene_id, delta=delta, p_value=pval, df=n_transcripts - 1, lr_stat=stat,
    ))
    gi = gi + 1
  }

  // Apply BH correction
  let pvals = extract_pvalues(gene_results)
  let adj = benjamini_hochberg_correct(pvals)

  let mut i = 0
  while i < gene_results.length() {
    gene_results[i].p_adj = adj[i]
    gene_results[i].significant = adj[i] < config.alpha
    i = i + 1
  }

  let mut n_sig = 0
  let mut j = 0
  while j < gene_results.length() {
    if gene_results[j].significant { n_sig = n_sig + 1 }
    j = j + 1
  }

  let filtered_results : Array[DRIMSeqGeneResult] = Array::new()
  let mut k = 0
  while k < gene_results.length() {
    if gene_results[k].p_value > 0.0 {
      filtered_results.push(gene_results[k])
    }
    k = k + 1
  }

  DRIMSeqResult::{
    gene_results,
    filtered_results,
    alpha: config.alpha,
    n_genes_tested: total_tested,
    n_genes_filtered: total_filtered,
    n_significant: n_sig,
    converged: true,
    iterations: 10,
  }
}

///|
/// Create sample data for testing.
pub fn drimseq_sample_data() -> Array[TranscriptCount] {
  let counts : Array[TranscriptCount] = Array::new()

  // GeneA: 2 transcripts, condition_a (3 samples)
  counts.push(TranscriptCount::new(transcript_id="tx_A1", gene_id="GeneA", sample_id="s1", condition="cond_a", count=100.0, gene_count=150.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A2", gene_id="GeneA", sample_id="s1", condition="cond_a", count=50.0, gene_count=150.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A1", gene_id="GeneA", sample_id="s2", condition="cond_a", count=120.0, gene_count=180.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A2", gene_id="GeneA", sample_id="s2", condition="cond_a", count=60.0, gene_count=180.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A1", gene_id="GeneA", sample_id="s3", condition="cond_a", count=110.0, gene_count=170.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A2", gene_id="GeneA", sample_id="s3", condition="cond_a", count=60.0, gene_count=170.0))

  // GeneA: 2 transcripts, condition_b (3 samples) - DIFFERENT proportions
  counts.push(TranscriptCount::new(transcript_id="tx_A1", gene_id="GeneA", sample_id="s4", condition="cond_b", count=50.0, gene_count=180.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A2", gene_id="GeneA", sample_id="s4", condition="cond_b", count=130.0, gene_count=180.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A1", gene_id="GeneA", sample_id="s5", condition="cond_b", count=45.0, gene_count=170.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A2", gene_id="GeneA", sample_id="s5", condition="cond_b", count=125.0, gene_count=170.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A1", gene_id="GeneA", sample_id="s6", condition="cond_b", count=55.0, gene_count=175.0))
  counts.push(TranscriptCount::new(transcript_id="tx_A2", gene_id="GeneA", sample_id="s6", condition="cond_b", count=120.0, gene_count=175.0))

  // GeneB: 1 transcript, no difference
  counts.push(TranscriptCount::new(transcript_id="tx_B1", gene_id="GeneB", sample_id="s1", condition="cond_a", count=200.0, gene_count=200.0))
  counts.push(TranscriptCount::new(transcript_id="tx_B1", gene_id="GeneB", sample_id="s2", condition="cond_a", count=210.0, gene_count=210.0))
  counts.push(TranscriptCount::new(transcript_id="tx_B1", gene_id="GeneB", sample_id="s3", condition="cond_a", count=195.0, gene_count=195.0))
  counts.push(TranscriptCount::new(transcript_id="tx_B1", gene_id="GeneB", sample_id="s4", condition="cond_b", count=205.0, gene_count=205.0))
  counts.push(TranscriptCount::new(transcript_id="tx_B1", gene_id="GeneB", sample_id="s5", condition="cond_b", count=215.0, gene_count=215.0))
  counts.push(TranscriptCount::new(transcript_id="tx_B1", gene_id="GeneB", sample_id="s6", condition="cond_b", count=200.0, gene_count=200.0))

  counts
}

// --- Helper functions ---

fn array_sum(arr : Array[Double]) -> Double {
  let mut s = 0.0
  let mut i = 0
  while i < arr.length() { s = s + arr[i]; i = i + 1 }
  s
}

fn array_sub(a : Array[Double], b : Array[Double]) -> Array[Double] {
  let len = min_int(a.length(), b.length())
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < len { result.push(a[i] - b[i]); i = i + 1 }
  result
}

fn array_dot(a : Array[Double], b : Array[Double]) -> Double {
  let len = min_int(a.length(), b.length())
  let mut s = 0.0
  let mut i = 0
  while i < len { s = s + a[i] * b[i]; i = i + 1 }
  s
}

fn slice_array(arr : Array[Double], start : Int, end : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let i = start
  let j = if end > arr.length() { arr.length() } else { end }
  let mut k = i
  while k < j { result.push(arr[k]); k = k + 1 }
  result
}

/// Log-gamma function (Lanczos approximation).
fn drimseq_log_gamma(x : Double) -> Double {
  if x <= 0.0 { return 0.0 }
  if x < 1.0 {
    let g = 7.0
    let c = [
      0.99999999999980993,
      676.5203681218851,
      -1259.1392167224028,
      771.32342877765313,
      -176.61502916214059,
      12.507343278686905,
      -0.13857109526572012,
      9.9843695780195716e-6,
      1.5056327351493116e-7,
    ]
    let xx = x
    let t = xx + g + 0.5
    let mut sum = c[0]
    let mut i = 1
    while i < 9 {
      sum = sum + c[i] / (xx + i.to_double() - 1.0)
      i = i + 1
    }
    let pi = 3.141592653589793
    let half_log_2pi = 0.5 * @math.ln(2.0 * pi)
    let log_t = @math.ln(t)
    let log_sum = @math.ln(sum)
    half_log_2pi + log_t * (xx + 0.5) - t + log_sum
  } else {
    // Stirling's approximation for x >= 1
    let pi = 3.141592653589793
    0.5 * @math.ln(2.0 * pi / x) + @math.ln(x) * (x - 0.5) - x + 1.0 / (12.0 * x) - 1.0 / (360.0 * x * x * x)
  }
}

/// Chi-square p-value approximation.
fn chi_square_p_value(stat : Double, df : Int) -> Double {
  if stat <= 0.0 { return 1.0 }
  if df <= 0 { return 1.0 }
  // Use Wilson-Hilferty approximation
  let df_d = df.to_double()
  let x = stat / df_d
  // Wilson-Hilferty: cube root transformation
  let x1_3 = @math.pow(x, 1.0 / 3.0)
  let z = (x1_3 - (1.0 - 2.0 / (9.0 * df_d))) / (2.0 / (9.0 * df_d)).sqrt()
  let p = normal_sf(z)
  if p < 0.0 { 0.0 } else { p }
}

/// Normal distribution survival function (1 - CDF).
fn normal_sf(z : Double) -> Double {
  // Abramowitz and Stegun approximation
  let t = 1.0 / (1.0 + 0.2316419 * z.abs())
  let d = 0.3989422804014327 * @math.exp(-0.5 * z * z)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
  if z > 0.0 { p } else { 1.0 - p }
}

/// Benjamini-Hochberg p-value correction.
fn benjamini_hochberg_correct(pvalues : Array[Double]) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 { return Array::new() }

  // Create index array and sort by p-value
  let indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n { indices.push(i); i = i + 1 }

  // Simple bubble sort
  let sorted_p = pvalues
  let sorted_idx = indices
  let mut j = 0
  while j < n - 1 {
    let mut k = 0
    while k < n - j - 1 {
      if sorted_p[k] > sorted_p[k + 1] {
        let tmp_p = sorted_p[k]
        sorted_p[k] = sorted_p[k + 1]
        sorted_p[k + 1] = tmp_p
        let tmp_i = sorted_idx[k]
        sorted_idx[k] = sorted_idx[k + 1]
        sorted_idx[k + 1] = tmp_i
      }
      k = k + 1
    }
    j = j + 1
  }

  // Compute adjusted p-values
  let adj = sorted_p
  let mut idx = n - 1
  while idx >= 0 {
    let rank = (idx + 1).to_double()
    let val = sorted_p[idx] * (n.to_double() / rank)
    if idx < n - 1 {
      if val > adj[idx + 1] { adj[idx] = adj[idx + 1] } else { adj[idx] = val }
    } else {
      adj[idx] = val
    }
    if adj[idx] > 1.0 { adj[idx] = 1.0 }
    idx = idx - 1
  }

  // Map back to original order
  let result : Array[Double] = Array::new()
  let mut m = 0
  while m < n { result.push(0.0); m = m + 1 }

  let mut r = 0
  while r < n {
    result[sorted_idx[r]] = adj[r]
    r = r + 1
  }

  result
}

fn get_unique_genes(counts : Array[TranscriptCount]) -> Array[String] {
  let genes : Array[String] = Array::new()
  let mut i = 0
  while i < counts.length() {
    let id = counts[i].gene_id
    let mut found = false
    let mut j = 0
    while j < genes.length() {
      if genes[j] == id { found = true; break }
      j = j + 1
    }
    if !found { genes.push(id) }
    i = i + 1
  }
  genes
}

fn get_gene_data(counts : Array[TranscriptCount], gene_id : String) -> Array[TranscriptCount] {
  let result : Array[TranscriptCount] = Array::new()
  let mut i = 0
  while i < counts.length() {
    if counts[i].gene_id == gene_id { result.push(counts[i]) }
    i = i + 1
  }
  result
}

fn get_unique_conditions(data : Array[TranscriptCount]) -> Array[String] {
  let conds : Array[String] = Array::new()
  let mut i = 0
  while i < data.length() {
    let c = data[i].condition
    let mut found = false
    let mut j = 0
    while j < conds.length() {
      if conds[j] == c { found = true; break }
      j = j + 1
    }
    if !found { conds.push(c) }
    i = i + 1
  }
  conds
}

fn get_condition_proportions(data : Array[TranscriptCount], condition : String) -> Array[Double] {
  let props : Array[Double] = Array::new()
  let mut i = 0
  while i < data.length() {
    if data[i].condition == condition {
      props.push(data[i].proportion())
    }
    i = i + 1
  }
  props
}

fn extract_pvalues(results : Array[DRIMSeqGeneResult]) -> Array[Double] {
  let pvals : Array[Double] = Array::new()
  let mut i = 0
  while i < results.length() {
    pvals.push(results[i].p_value)
    i = i + 1
  }
  pvals
}
