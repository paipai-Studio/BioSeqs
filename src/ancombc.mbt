// Bioconductor ANCOMBC - Analysis of Composition of Microbiomes with Bias Correction.
//
// Inspired by the Bioconductor ANCOMBC package (Lin and Peddada 2020,
// Nat Commun 11:5144), which performs differential abundance (DA) testing
// for microbiome data while correcting for the compositional bias inherent
// in sequencing count data. This pure-MoonBit module implements:
// - OTU/ASV count table container
// - Sampling fraction estimation
// - Log-ratio bias correction (relative to a reference feature)
// - Two-group differential abundance test (Welch's t-test on bias-corrected
//   log-abundances) with ANCOM-style W statistics
// - Benjamini-Hochberg FDR correction
//
// Reference:
//   https://bioconductor.org/packages/ANCOMBC/

// ============================================================================
// Private statistical helpers
// ============================================================================

///|
/// Natural log, clamped to avoid log(0).
fn abc_safe_ln(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  @math.ln(x)
}

///|
/// Arithmetic mean of an array of doubles. Returns 0.0 when empty.
fn abc_mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    s = s + v
  }
  s / n.to_double()
}

///|
/// Sample variance (n-1 denominator). Returns 0.0 when fewer than 2 elements.
fn abc_variance(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let m = abc_mean(arr)
  let mut s = 0.0
  for v in arr {
    let d = v - m
    s = s + d * d
  }
  s / (n - 1).to_double()
}

///|
/// Welch's t-statistic for two groups. Returns (t, df) where df is the
/// Welch-Satterthwaite degrees of freedom.
fn abc_welch_t(
  group1 : Array[Double],
  group2 : Array[Double]
) -> (Double, Double) {
  let n1 = group1.length()
  let n2 = group2.length()
  if n1 == 0 || n2 == 0 {
    return (0.0, 0.0)
  }
  let m1 = abc_mean(group1)
  let m2 = abc_mean(group2)
  let v1 = abc_variance(group1)
  let v2 = abc_variance(group2)
  let se = (v1 / n1.to_double() + v2 / n2.to_double()).sqrt()
  if se == 0.0 {
    return (0.0, 0.0)
  }
  let t = (m1 - m2) / se
  // Welch-Satterthwaite df.
  let num = (v1 / n1.to_double() + v2 / n2.to_double())
  let denom1 = v1 * v1 / ((n1.to_double() - 1.0) * n1.to_double() * n1.to_double())
  let denom2 = v2 * v2 / ((n2.to_double() - 1.0) * n2.to_double() * n2.to_double())
  let df = if denom1 + denom2 > 0.0 {
    num * num / (denom1 + denom2)
  } else {
    (n1 + n2 - 2).to_double()
  }
  (t, df)
}

///|
/// Two-tailed p-value from a t-statistic using a normal approximation
/// (sufficient for large samples; microbiome studies typically have
// dozens of samples).
fn abc_t_pvalue(t : Double, df : Double) -> Double {
  if df <= 0.0 {
    return 1.0
  }
  // Use normal approximation: p = 2 * (1 - Phi(|t|))
  let z = t.abs()
  // Standard normal CDF approximation (Abramowitz & Stegun 26.2.17).
  let phi = abc_normal_cdf(z)
  let p = 2.0 * (1.0 - phi)
  if p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Standard normal CDF approximation (Abramowitz & Stegun 26.2.17).
fn abc_normal_cdf(z : Double) -> Double {
  let sign = if z < 0.0 { -1.0 } else { 1.0 }
  let x = z.abs() / 2.0.sqrt()
  let t = 1.0 / (1.0 + 0.3275911 * x)
  let exp_val = @math.exp(0.0 - x * x / 2.0)
  let y2 = 1.0 -
    (((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t -
    0.284496736) *
    t +
    0.254829592) *
    t) *
    exp_val
  0.5 * (1.0 + sign * y2)
}

///|
/// Benjamini-Hochberg FDR correction. Returns adjusted p-values (q-values)
/// in the same order as the input.
fn abc_bh_fdr(pvalues : Array[Double]) -> Array[Double] {
  let m = pvalues.length()
  if m == 0 {
    return []
  }
  // Create index-value pairs and sort by p-value ascending.
  let indexed : Array[(Int, Double)] = Array::new()
  for i in 0..<m {
    indexed.push((i, pvalues[i]))
  }
  indexed.sort_by(fn(a : (Int, Double), b : (Int, Double)) -> Int {
    if a.1 < b.1 { -1 } else if a.1 > b.1 { 1 } else { 0 }
  })
  // Compute adjusted p-values from largest to smallest.
  let adjusted : Array[Double] = Array::make(m, 0.0)
  let mut prev = 1.0
  let mut k = m - 1
  while k >= 0 {
    let i = k + 1 // 1-based rank
    let (orig_idx, p) = indexed[k]
    let raw_adj = p * m.to_double() / i.to_double()
    let adj = if raw_adj < prev { raw_adj } else { prev }
    adjusted[orig_idx] = if adj > 1.0 { 1.0 } else { adj }
    prev = adj
    k = k - 1
  }
  adjusted
}

// ============================================================================
// Data containers
// ============================================================================

///|
/// A single OTU/ASV feature: its name and per-sample raw counts.
pub struct AncombcFeature {
  name : String
  counts : Array[Int] // per-sample counts
}

///|
/// Construct an AncombcFeature.
pub fn AncombcFeature::new(name : String, counts : Array[Int]) -> AncombcFeature {
  AncombcFeature::{ name, counts }
}

///|
pub fn AncombcFeature::name(self : AncombcFeature) -> String {
  self.name
}

///|
pub fn AncombcFeature::counts(self : AncombcFeature) -> Array[Int] {
  self.counts
}

///|
/// Total count for this feature across all samples.
pub fn AncombcFeature::total(self : AncombcFeature) -> Int {
  let mut s = 0
  for c in self.counts {
    s = s + c
  }
  s
}

///|
/// A single sample: its name and group label.
pub struct AncombcSample {
  name : String
  group : String
}

///|
/// Construct an AncombcSample.
pub fn AncombcSample::new(name : String, group : String) -> AncombcSample {
  AncombcSample::{ name, group }
}

///|
pub fn AncombcSample::name(self : AncombcSample) -> String {
  self.name
}

///|
pub fn AncombcSample::group(self : AncombcSample) -> String {
  self.group
}

///|
/// The main ANCOMBC data container: features (taxa) × samples with
/// sample metadata (group labels).
pub struct AncombcData {
  sample_names : Array[String]
  sample_groups : Array[String]
  feature_names : Array[String]
  counts : Array[Array[Int]] // [feature][sample] -> count
}

///|
/// Construct an AncombcData from parallel arrays.
pub fn AncombcData::new(
  sample_names : Array[String],
  sample_groups : Array[String],
  feature_names : Array[String],
  counts : Array[Array[Int]]
) -> AncombcData {
  AncombcData::{ sample_names, sample_groups, feature_names, counts }
}

///|
pub fn AncombcData::n_samples(self : AncombcData) -> Int {
  self.sample_names.length()
}

///|
pub fn AncombcData::n_features(self : AncombcData) -> Int {
  self.feature_names.length()
}

///|
pub fn AncombcData::sample_names(self : AncombcData) -> Array[String] {
  self.sample_names
}

///|
pub fn AncombcData::sample_groups(self : AncombcData) -> Array[String] {
  self.sample_groups
}

///|
pub fn AncombcData::feature_names(self : AncombcData) -> Array[String] {
  self.feature_names
}

///|
pub fn AncombcData::counts(self : AncombcData) -> Array[Array[Int]] {
  self.counts
}

///|
/// Get the count for a specific feature and sample.
pub fn AncombcData::get_count(
  self : AncombcData,
  feature_idx : Int,
  sample_idx : Int
) -> Int {
  self.counts[feature_idx][sample_idx]
}

///|
/// Unique group labels in the data.
pub fn AncombcData::groups(self : AncombcData) -> Array[String] {
  let seen : Map[String, Bool] = Map([], capacity=4)
  let result : Array[String] = Array::new()
  for g in self.sample_groups {
    if !seen.contains(g) {
      seen[g] = true
      result.push(g)
    }
  }
  result
}

///|
/// Total library size for each sample (column sums).
pub fn AncombcData::library_sizes(self : AncombcData) -> Array[Int] {
  let n = self.n_samples()
  let sizes = Array::make(n, 0)
  for f in self.counts {
    for j in 0..<n {
      sizes[j] = sizes[j] + f[j]
    }
  }
  sizes
}

// ============================================================================
// Sampling fraction estimation
// ============================================================================

///|
/// Estimate the sampling fraction for each sample. The sampling fraction
/// s_j is the proportion of sample j's library size relative to the total
/// across all samples. Returns an array of doubles (one per sample).
pub fn ancombc_sampling_fractions(data : AncombcData) -> Array[Double] {
  let sizes = data.library_sizes()
  let mut total = 0.0
  for s in sizes {
    total = total + s.to_double()
  }
  if total <= 0.0 {
    return Array::make(data.n_samples(), 0.0)
  }
  let fracs = Array::new()
  for s in sizes {
    fracs.push(s.to_double() / total)
  }
  fracs
}

// ============================================================================
// Bias-corrected log-abundance
// ============================================================================

///|
/// Compute bias-corrected log-abundances for each feature in each sample.
///
/// The bias-corrected value for feature i in sample j is:
///   y[i,j] = log(count[i,j] + 0.5) - log(s_j * N + 0.5)
/// where s_j is the sampling fraction and N is the number of features.
/// This adjusts for differences in sequencing depth (library size).
pub fn ancombc_bias_corrected_log(
  data : AncombcData
) -> Array[Array[Double]] {
  let n_features = data.n_features()
  let n_samples = data.n_samples()
  let fracs = ancombc_sampling_fractions(data)
  let n = n_features.to_double()
  let result : Array[Array[Double]] = Array::new()
  for i in 0..<n_features {
    let row : Array[Double] = Array::new()
    for j in 0..<n_samples {
      let count = data.counts[i][j].to_double()
      let adjusted_total = fracs[j] * n
      let log_count = abc_safe_ln(count + 0.5)
      let log_total = abc_safe_ln(adjusted_total + 0.5)
      row.push(log_count - log_total)
    }
    result.push(row)
  }
  result
}

///|
/// Select a reference feature for log-ratio computation. The reference
/// is the feature with the lowest inter-sample variance in bias-corrected
/// log-abundance (i.e., the most "stable" feature). Returns the feature
/// index.
pub fn ancombc_reference_feature(data : AncombcData) -> Int {
  let log_data = ancombc_bias_corrected_log(data)
  let n_features = data.n_features()
  if n_features == 0 {
    return -1
  }
  let mut best_idx = 0
  let mut best_var = 1000000000.0 * 1000000000.0
  for i in 0..<n_features {
    let v = abc_variance(log_data[i])
    if v < best_var {
      best_var = v
      best_idx = i
    }
  }
  best_idx
}

// ============================================================================
// Differential abundance test
// ============================================================================

///|
/// Result of a differential abundance test for a single feature.
pub struct AncombcResult {
  feature : String
  log_fold_change : Double
  w_stat : Int
  p_value : Double
  q_value : Double
  significant : Bool
}

///|
/// Construct an AncombcResult.
pub fn AncombcResult::new(
  feature : String,
  log_fold_change : Double,
  w_stat : Int,
  p_value : Double,
  q_value : Double,
  significant : Bool
) -> AncombcResult {
  AncombcResult::{ feature, log_fold_change, w_stat, p_value, q_value, significant }
}

///|
pub fn AncombcResult::feature(self : AncombcResult) -> String {
  self.feature
}

///|
pub fn AncombcResult::log_fold_change(self : AncombcResult) -> Double {
  self.log_fold_change
}

///|
pub fn AncombcResult::w_stat(self : AncombcResult) -> Int {
  self.w_stat
}

///|
pub fn AncombcResult::p_value(self : AncombcResult) -> Double {
  self.p_value
}

///|
pub fn AncombcResult::q_value(self : AncombcResult) -> Double {
  self.q_value
}

///|
pub fn AncombcResult::significant(self : AncombcResult) -> Bool {
  self.significant
}

///|
/// The full ANCOMBC test result: an array of per-feature results plus
/// metadata.
pub struct AncombcTestResult {
  results : Array[AncombcResult]
  reference_feature : String
  group1 : String
  group2 : String
  n_samples_group1 : Int
  n_samples_group2 : Int
  alpha : Double
}

///|
pub fn AncombcTestResult::results(self : AncombcTestResult) -> Array[AncombcResult] {
  self.results
}

///|
pub fn AncombcTestResult::reference_feature(
  self : AncombcTestResult
) -> String {
  self.reference_feature
}

///|
pub fn AncombcTestResult::group1(self : AncombcTestResult) -> String {
  self.group1
}

///|
pub fn AncombcTestResult::group2(self : AncombcTestResult) -> String {
  self.group2
}

///|
pub fn AncombcTestResult::alpha(self : AncombcTestResult) -> Double {
  self.alpha
}

///|
/// Number of significant features.
pub fn AncombcTestResult::n_significant(self : AncombcTestResult) -> Int {
  let mut n = 0
  for r in self.results {
    if r.significant {
      n = n + 1
    }
  }
  n
}

///|
/// Run the ANCOMBC differential abundance test between two groups.
///
/// The algorithm:
/// 1. Compute bias-corrected log-abundances.
/// 2. For each feature, subtract the log-abundance of the reference feature
///    (log-ratio transformation) to remove compositional bias.
/// 3. Perform Welch's t-test on the log-ratios between the two groups.
/// 4. Compute the ANCOM W statistic: for each feature, count how many
///    other features have a significantly different log-ratio relative to
///    this feature (at the given alpha level).
/// 5. Apply BH-FDR correction to the p-values.
/// 6. Declare a feature significant if its q-value <= alpha AND its W
///    statistic exceeds the 0.7 * n_features threshold.
pub fn ancombc_test(
  data : AncombcData,
  group1 : String,
  group2 : String,
  alpha : Double
) -> AncombcTestResult {
  let n_features = data.n_features()
  let n_samples = data.n_samples()
  // Assign samples to groups.
  let g1_indices : Array[Int] = Array::new()
  let g2_indices : Array[Int] = Array::new()
  for j in 0..<n_samples {
    if data.sample_groups[j] == group1 {
      g1_indices.push(j)
    } else if data.sample_groups[j] == group2 {
      g2_indices.push(j)
    }
  }
  // Compute bias-corrected log-abundances.
  let log_data = ancombc_bias_corrected_log(data)
  // Select reference feature.
  let ref_idx = ancombc_reference_feature(data)
  let ref_name = if ref_idx >= 0 {
    data.feature_names[ref_idx]
  } else {
    ""
  }
  // Compute log-ratios relative to reference: y_i - y_ref for each sample.
  let log_ratios : Array[Array[Double]] = Array::new()
  for i in 0..<n_features {
    let row : Array[Double] = Array::new()
    for j in 0..<n_samples {
      if ref_idx >= 0 {
        row.push(log_data[i][j] - log_data[ref_idx][j])
      } else {
        row.push(log_data[i][j])
      }
    }
    log_ratios.push(row)
  }
  // Compute pairwise log-ratios between all features for W statistic.
  // For each feature i, W_i = number of features j (j != i) where the
  // log-ratio (y_i - y_j) has a significant t-test between groups.
  let pvalues : Array[Double] = Array::make(n_features, 1.0)
  let w_stats : Array[Int] = Array::make(n_features, 0)
  let log_fcs : Array[Double] = Array::make(n_features, 0.0)
  // For each feature, compute t-test on log-ratios relative to reference.
  for i in 0..<n_features {
    let g1_vals : Array[Double] = Array::new()
    let g2_vals : Array[Double] = Array::new()
    for j in g1_indices {
      g1_vals.push(log_ratios[i][j])
    }
    for j in g2_indices {
      g2_vals.push(log_ratios[i][j])
    }
    let (t, df) = abc_welch_t(g1_vals, g2_vals)
    let p = abc_t_pvalue(t, df)
    pvalues[i] = p
    let m1 = abc_mean(g1_vals)
    let m2 = abc_mean(g2_vals)
    log_fcs[i] = m1 - m2
  }
  // Compute W statistics using pairwise comparisons.
  for i in 0..<n_features {
    for j in 0..<n_features {
      if i == j {
        continue
      }
      // Log-ratio of feature i relative to feature j.
      let g1_vals : Array[Double] = Array::new()
      let g2_vals : Array[Double] = Array::new()
      for k in g1_indices {
        g1_vals.push(log_data[i][k] - log_data[j][k])
      }
      for k in g2_indices {
        g2_vals.push(log_data[i][k] - log_data[j][k])
      }
      let (t, df) = abc_welch_t(g1_vals, g2_vals)
      let p = abc_t_pvalue(t, df)
      if p < alpha {
        w_stats[i] = w_stats[i] + 1
      }
    }
  }
  // BH-FDR correction.
  let qvalues = abc_bh_fdr(pvalues)
  // W threshold: 0.7 * (n_features - 1).
  let w_threshold = ((n_features - 1).to_double() * 0.7).to_int()
  // Build results.
  let results : Array[AncombcResult] = Array::new()
  for i in 0..<n_features {
    let sig = qvalues[i] <= alpha && w_stats[i] > w_threshold
    results.push(
      AncombcResult::new(
        data.feature_names[i],
        log_fcs[i],
        w_stats[i],
        pvalues[i],
        qvalues[i],
        sig,
      )
    )
  }
  AncombcTestResult::{
    results,
    reference_feature: ref_name,
    group1,
    group2,
    n_samples_group1: g1_indices.length(),
    n_samples_group2: g2_indices.length(),
    alpha,
  }
}

// ============================================================================
// Summary
// ============================================================================

///|
/// Text summary of the ANCOMBC test result.
pub fn ancombc_result_summary(result : AncombcTestResult) -> String {
  let sb = StringBuilder::new()
  sb.write_string("=== ANCOMBC Differential Abundance Test ===\n")
  sb.write_string("Comparison: ")
  sb.write_string(result.group1)
  sb.write_string(" vs ")
  sb.write_string(result.group2)
  sb.write_char('\n')
  sb.write_string("Reference feature: ")
  sb.write_string(result.reference_feature)
  sb.write_char('\n')
  sb.write_string("Samples: ")
  sb.write_string(result.n_samples_group1.to_string())
  sb.write_string(" vs ")
  sb.write_string(result.n_samples_group2.to_string())
  sb.write_char('\n')
  sb.write_string("Alpha: ")
  sb.write_string(result.alpha.to_string())
  sb.write_char('\n')
  sb.write_string("Significant features: ")
  sb.write_string(result.n_significant().to_string())
  sb.write_char('\n')
  sb.write_string("--- Per-feature results ---\n")
  sb.write_string(
    "feature\tlogFC\tW\tp_value\tq_value\tsignificant\n",
  )
  for r in result.results {
    sb.write_string(r.feature)
    sb.write_char('\t')
    sb.write_string(r.log_fold_change.to_string())
    sb.write_char('\t')
    sb.write_string(r.w_stat.to_string())
    sb.write_char('\t')
    sb.write_string(r.p_value.to_string())
    sb.write_char('\t')
    sb.write_string(r.q_value.to_string())
    sb.write_char('\t')
    sb.write_string(r.significant.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Build sample microbiome count data: 6 features (taxa) × 8 samples
/// in 2 groups (4 "control" + 4 "treatment"). The treatment group has
/// elevated counts for features 0 and 1 (truly differential) and similar
/// counts for the rest.
pub fn ancombc_sample_data() -> AncombcData {
  let sample_names = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"]
  let sample_groups = [
    "control", "control", "control", "control", "treatment", "treatment",
    "treatment", "treatment",
  ]
  let feature_names = ["Bacteroides", "Prevotella", "Faecalibacterium", "Roseburia", "Eubacterium", "Ruminococcus"]
  // Counts: [feature][sample]
  // Features 0,1 are elevated in treatment; features 2-5 are similar.
  let counts = [
    [120, 110, 130, 115, 250, 240, 260, 245], // Bacteroides - differential
    [80, 85, 75, 90, 180, 170, 190, 175], // Prevotella - differential
    [50, 55, 48, 52, 55, 50, 58, 53], // Faecalibacterium - not differential
    [30, 35, 28, 32, 33, 30, 35, 31], // Roseburia
    [40, 42, 38, 41, 43, 40, 45, 42], // Eubacterium
    [25, 28, 22, 26, 27, 25, 29, 24], // Ruminococcus
  ]
  AncombcData::new(sample_names, sample_groups, feature_names, counts)
}
