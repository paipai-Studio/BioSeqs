///|
/// lefse - LEfSe (Linear Discriminant Analysis Effect Size) algorithm
/// for microbiome biomarker discovery.
/// Inspired by the Huttenhower lab LEfSe tool (Segata et al. 2011,
/// Genome Biol 12:R60), commonly used via Galaxy/lefse and conda lefse.
/// Implements the full LEfSe pipeline:
/// 1. Kruskal-Wallis test across classes
/// 2. Optional class-subclass Wilcoxon pairwise tests
/// 3. LDA effect size computation
/// 4. Threshold filtering (KW p-value < alpha and |LDA| > threshold)

// ============================================================================
// Private statistical helpers
// ============================================================================

///|
/// Sum of an array of doubles.
fn lefse_sum(arr : Array[Double]) -> Double {
  let mut s = 0.0
  let mut i = 0
  while i < arr.length() {
    s = s + arr[i]
    i = i + 1
  }
  s
}

///|
/// Arithmetic mean of an array of doubles. Returns 0.0 when empty.
fn lefse_mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  lefse_sum(arr) / n.to_double()
}

///|
/// Compute ranks with average tie-breaking for a combined array.
/// Returns an array of ranks in the same order as the input values.
fn lefse_rank(data : Array[Double]) -> Array[Double] {
  let n = data.length()
  let result : Array[Double] = Array::make(n, 0.0)

  // Build indexed array: [value, original_index]
  let indexed : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((data[i], i))
    i = i + 1
  }

  // Sort by value ascending
  indexed.sort_by(fn(a : (Double, Int), b : (Double, Int)) -> Int {
    if a.0 < b.0 {
      -1
    } else if a.0 > b.0 {
      1
    } else {
      0
    }
  })

  // Assign ranks with average tie-breaking
  let mut pos = 0
  while pos < n {
    let mut end = pos
    while end < n - 1 && indexed[end + 1].0 == indexed[pos].0 {
      end = end + 1
    }
    // Average rank for tied group (1-based)
    let avg_rank = (pos + end).to_double() / 2.0 + 1.0
    let mut j = pos
    while j <= end {
      result[indexed[j].1] = avg_rank
      j = j + 1
    }
    pos = end + 1
  }

  result
}

///|
/// Standard normal CDF using Abramowitz & Stegun approximation.
fn lefse_normal_cdf(x : Double) -> Double {
  let sign = if x < 0.0 { -1.0 } else { 1.0 }
  let abs_x = x.abs()
  let t = 1.0 / (1.0 + 0.2316419 * abs_x)
  let d = 0.3989423 * @math.exp(-abs_x * abs_x / 2.0)
  let p = d *
    t *
    (
      0.3193815 +
      t * (-0.3565638 + t * (1.781478 + t * (-1.8212560 + t * 1.3302744)))
    )
  if x > 0.0 {
    1.0 - p
  } else {
    p
  }
}

///|
/// Two-tailed p-value from a Z-score.
fn lefse_pvalue_z(z : Double) -> Double {
  let p = 2.0 * (1.0 - lefse_normal_cdf(z.abs()))
  if p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Chi-square survival function (upper tail p-value) using the
/// Wilson-Hilferty approximation for chi-square distribution with df degrees
/// of freedom.
fn lefse_chi2_sf(x : Double, df : Int) -> Double {
  if x <= 0.0 || df <= 0 {
    return 1.0
  }
  let df_d = df.to_double()
  let z = ((x / df_d).pow(1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df_d))) /
    (2.0 / (9.0 * df_d)).sqrt()
  lefse_pvalue_z(z)
}

///|
/// Get unique string values from an array, preserving first-occurrence order.
fn lefse_unique_strings(arr : Array[String]) -> Array[String] {
  let seen : Map[String, Bool] = Map([], capacity=arr.length())
  let result : Array[String] = Array::new()
  let mut i = 0
  while i < arr.length() {
    if !seen.contains(arr[i]) {
      seen[arr[i]] = true
      result.push(arr[i])
    }
    i = i + 1
  }
  result
}

// ============================================================================
// Kruskal-Wallis test
// ============================================================================

///|
/// Perform the Kruskal-Wallis H test for comparing k independent groups.
///
/// The test ranks all observations jointly, then computes:
///   H = (12 / (N * (N+1))) * sum_i(n_i * (R_i_bar - (N+1)/2)^2)
/// with a tie correction factor applied.
/// The p-value is obtained from the chi-square distribution with k-1 df.
///
/// Returns (h_statistic, p_value).
pub fn lefse_kruskal_wallis(groups : Array[Array[Double]]) -> (Double, Double) {
  let k = groups.length()
  if k < 2 {
    return (0.0, 1.0)
  }

  // Check each group has at least one observation
  let mut i = 0
  while i < k {
    if groups[i].length() == 0 {
      return (0.0, 1.0)
    }
    i = i + 1
  }

  // Combine all values and compute ranks
  let combined : Array[Double] = Array::new()
  let group_sizes : Array[Int] = Array::new()
  let mut total_n = 0
  i = 0
  while i < k {
    let n_i = groups[i].length()
    group_sizes.push(n_i)
    let mut j = 0
    while j < n_i {
      combined.push(groups[i][j])
      j = j + 1
    }
    total_n = total_n + n_i
    i = i + 1
  }

  if total_n < 2 {
    return (0.0, 1.0)
  }

  let ranks = lefse_rank(combined)
  let n_dbl = total_n.to_double()

  // Compute sum of ranks for each group
  let group_rank_sums : Array[Double] = Array::make(k, 0.0)
  let mut offset = 0
  i = 0
  while i < k {
    let mut r_sum = 0.0
    let mut j = 0
    while j < group_sizes[i] {
      r_sum = r_sum + ranks[offset + j]
      j = j + 1
    }
    group_rank_sums[i] = r_sum
    offset = offset + group_sizes[i]
    i = i + 1
  }

  // Compute H statistic
  // H = (12 / (N*(N+1))) * sum_i(n_i * (R_i_bar - (N+1)/2)^2)
  let expected_rank_mean = (n_dbl + 1.0) / 2.0
  let mut h = 0.0
  i = 0
  while i < k {
    let n_i = group_sizes[i].to_double()
    let r_bar = group_rank_sums[i] / n_i
    let diff = r_bar - expected_rank_mean
    h = h + n_i * diff * diff
    i = i + 1
  }
  h = h * 12.0 / (n_dbl * (n_dbl + 1.0))

  // Tie correction factor
  // C = 1 - sum(t_j^3 - t_j) / (N^3 - N)
  // where t_j is the size of each tied group
  let sorted_vals = combined.copy()
  sorted_vals.sort()

  let mut tie_sum = 0.0
  let mut pos = 0
  while pos < total_n {
    let mut end = pos
    while end < total_n - 1 && sorted_vals[end + 1] == sorted_vals[pos] {
      end = end + 1
    }
    let tie_size = (end - pos + 1).to_double()
    if tie_size > 1.0 {
      tie_sum = tie_sum + tie_size * tie_size * tie_size - tie_size
    }
    pos = end + 1
  }

  let correction = 1.0 - tie_sum / (n_dbl * n_dbl * n_dbl - n_dbl)
  if correction > 0.0 {
    h = h / correction
  }

  // p-value from chi-square with k-1 degrees of freedom
  let df = k - 1
  let p = lefse_chi2_sf(h, df)

  (h, if p.is_nan() { 0.0 } else { p })
}

// ============================================================================
// Wilcoxon rank-sum test
// ============================================================================

///|
/// Perform the Wilcoxon rank-sum (Mann-Whitney U) test for two independent
/// samples.
///
/// Ranks all observations jointly, then computes:
///   U1 = R1 - n1*(n1+1)/2
///   U2 = n1*n2 - U1
///   U = min(U1, U2)
/// For large samples (n > 20), uses normal approximation with continuity
/// correction for the p-value.
///
/// Returns (u_statistic, p_value).
pub fn lefse_wilcoxon_rank_sum(
  x : Array[Double],
  y : Array[Double],
) -> (Double, Double) {
  let n1 = x.length()
  let n2 = y.length()
  if n1 == 0 || n2 == 0 {
    return (0.0, 1.0)
  }

  // Combine and rank
  let combined : Array[Double] = Array::new()
  let mut i = 0
  while i < n1 {
    combined.push(x[i])
    i = i + 1
  }
  i = 0
  while i < n2 {
    combined.push(y[i])
    i = i + 1
  }

  let ranks = lefse_rank(combined)

  // Sum of ranks for the first group (x)
  let mut r1 = 0.0
  i = 0
  while i < n1 {
    r1 = r1 + ranks[i]
    i = i + 1
  }

  // Compute U statistics
  let n1_d = n1.to_double()
  let n2_d = n2.to_double()
  let u1 = r1 - n1_d * (n1_d + 1.0) / 2.0
  let u2 = n1_d * n2_d - u1
  let u = if u1 < u2 { u1 } else { u2 }

  // Normal approximation
  let mean_u = n1_d * n2_d / 2.0

  // Tie correction for variance
  let sorted_vals = combined.copy()
  sorted_vals.sort()
  let total_n = n1 + n2
  let n_d = total_n.to_double()

  let mut tie_sum = 0.0
  let mut pos = 0
  while pos < total_n {
    let mut end = pos
    while end < total_n - 1 && sorted_vals[end + 1] == sorted_vals[pos] {
      end = end + 1
    }
    let tie_size = (end - pos + 1).to_double()
    if tie_size > 1.0 {
      tie_sum = tie_sum + tie_size * tie_size * tie_size - tie_size
    }
    pos = end + 1
  }

  let var_u = n1_d * n2_d / 12.0 * (n_d + 1.0 - tie_sum / (n_d * (n_d - 1.0)))

  if var_u <= 0.0 {
    return (u, 1.0)
  }

  // Continuity correction
  let z = if u > mean_u {
    (u - mean_u - 0.5) / var_u.sqrt()
  } else if u < mean_u {
    (u - mean_u + 0.5) / var_u.sqrt()
  } else {
    0.0
  }

  // Two-tailed p-value
  let p = lefse_pvalue_z(z)

  (u, if p.is_nan() { 0.0 } else { p })
}

// ============================================================================
// LDA score
// ============================================================================

///|
/// Compute the LDA effect size score for a feature.
///
/// Given the class means and the index of the class with the highest mean,
/// the LDA score is:
///   lda = log10(max_mean / mean_of_other_classes_combined)
///
/// where mean_of_other_classes_combined is the arithmetic mean of all
/// other class means (excluding the best class). The sign is positive
/// if the best class has the highest mean, negative otherwise.
/// If any mean is zero or negative, a small pseudocount is added.
pub fn lefse_lda_score(
  class_means : Array[Double],
  best_class_idx : Int,
) -> Double {
  let k = class_means.length()
  if k < 2 || best_class_idx < 0 || best_class_idx >= k {
    return 0.0
  }

  // Add pseudocount to avoid log(0)
  let pseudocount = 1.0e-10
  let mut i = 0
  let safe_means : Array[Double] = Array::new()
  while i < k {
    let v = class_means[i]
    safe_means.push(if v > 0.0 { v } else { pseudocount })
    i = i + 1
  }

  let max_mean = safe_means[best_class_idx]

  // Mean of other classes
  let mut other_sum = 0.0
  let mut other_count = 0
  i = 0
  while i < k {
    if i != best_class_idx {
      other_sum = other_sum + safe_means[i]
      other_count = other_count + 1
    }
    i = i + 1
  }

  if other_count == 0 || other_sum <= 0.0 {
    return 0.0
  }

  let other_mean = other_sum / other_count.to_double()

  if other_mean <= 0.0 || max_mean <= 0.0 {
    return 0.0
  }

  // LDA score = log10(max_mean / other_mean)
  // Convert from natural log: log10(x) = ln(x) / ln(10)
  let ln10 = 2.302585092994046
  @math.ln(max_mean / other_mean) / ln10
}

// ============================================================================
// Data structures
// ============================================================================

///|
/// Input data for LEfSe analysis.
///
/// features: abundance table, rows = features (taxa), cols = samples
/// feature_names: name of each feature (taxon)
/// sample_names: name of each sample
/// classes: class label for each sample (e.g. "Healthy", "Disease")
/// subclasses: subclass label for each sample (e.g. subject ID or
///   batch; can be empty strings if no subclass structure)
pub struct LEfSeInput {
  features : Array[Array[Double]]
  feature_names : Array[String]
  sample_names : Array[String]
  classes : Array[String]
  subclasses : Array[String]
}

///|
/// Construct an LEfSeInput.
pub fn LEfSeInput::new(
  features : Array[Array[Double]],
  feature_names : Array[String],
  sample_names : Array[String],
  classes : Array[String],
  subclasses : Array[String],
) -> LEfSeInput {
  LEfSeInput::{ features, feature_names, sample_names, classes, subclasses, }
}

///|
/// Result for a single feature from LEfSe analysis.
pub struct LEfSeResult {
  feature_name : String
  class_name : String
  kw_pvalue : Double
  lda_score : Double
  is_significant : Bool
}

///|
/// Construct an LEfSeResult.
pub fn LEfSeResult::new(
  feature_name : String,
  class_name : String,
  kw_pvalue : Double,
  lda_score : Double,
  is_significant : Bool,
) -> LEfSeResult {
  LEfSeResult::{
    feature_name,
    class_name,
    kw_pvalue,
    lda_score,
    is_significant,
  }
}

///|
pub fn LEfSeResult::feature_name(self : LEfSeResult) -> String {
  self.feature_name
}

///|
pub fn LEfSeResult::class_name(self : LEfSeResult) -> String {
  self.class_name
}

///|
pub fn LEfSeResult::kw_pvalue(self : LEfSeResult) -> Double {
  self.kw_pvalue
}

///|
pub fn LEfSeResult::lda_score(self : LEfSeResult) -> Double {
  self.lda_score
}

///|
pub fn LEfSeResult::is_significant(self : LEfSeResult) -> Bool {
  self.is_significant
}

///|
/// Output of LEfSe analysis containing all per-feature results.
pub struct LEfSeOutput {
  results : Array[LEfSeResult]
  n_significant : Int
}

///|
/// Construct an LEfSeOutput.
pub fn LEfSeOutput::new(
  results : Array[LEfSeResult],
  n_significant : Int,
) -> LEfSeOutput {
  LEfSeOutput::{ results, n_significant, }
}

///|
pub fn LEfSeOutput::results(self : LEfSeOutput) -> Array[LEfSeResult] {
  self.results
}

///|
pub fn LEfSeOutput::n_significant(self : LEfSeOutput) -> Int {
  self.n_significant
}

// ============================================================================
// LEfSe pipeline
// ============================================================================

///|
/// Run the full LEfSe pipeline for microbiome biomarker discovery.
///
/// Steps:
/// 1. For each feature, perform Kruskal-Wallis test across classes.
/// 2. If subclasses are provided (non-empty and not all identical),
///    perform pairwise Wilcoxon rank-sum tests within each subclass
///    pair for features that pass the KW test. A feature must be
///    significant in ALL subclass pair comparisons to proceed.
/// 3. For features passing steps 1-2, compute the LDA effect size.
/// 4. A feature is declared significant if:
///    - KW p-value < alpha (default 0.05)
///    - |LDA score| > lda_threshold (default 2.0)
///
/// Returns an LEfSeOutput with per-feature results and the count of
/// significant biomarkers.
pub fn lefse_run(
  input : LEfSeInput,
  alpha? : Double = 0.05,
  lda_threshold? : Double = 2.0,
) -> LEfSeOutput {
  let n_features = input.features.length()
  let n_samples = input.sample_names.length()

  if n_features == 0 || n_samples == 0 {
    return LEfSeOutput::new(Array::new(), 0)
  }

  // Get unique class labels
  let unique_classes = lefse_unique_strings(input.classes)
  let n_classes = unique_classes.length()

  if n_classes < 2 {
    return LEfSeOutput::new(Array::new(), 0)
  }

  // Check if subclasses are meaningful (not all empty or identical)
  let has_subclasses = lefse_has_subclasses(input.subclasses)

  // Build class-to-sample-index mapping
  let class_indices : Array[Array[Int]] = Array::new()
  let mut ci = 0
  while ci < n_classes {
    let indices : Array[Int] = Array::new()
    let mut j = 0
    while j < n_samples {
      if input.classes[j] == unique_classes[ci] {
        indices.push(j)
      }
      j = j + 1
    }
    class_indices.push(indices)
    ci = ci + 1
  }

  // If subclasses exist, build subclass-to-sample-index mapping per class
  let subclass_info = if has_subclasses {
    lefse_build_subclass_info(input.classes, input.subclasses, unique_classes)
  } else {
    Array::new()
  }

  // Process each feature
  let results : Array[LEfSeResult] = Array::new()
  let mut n_sig = 0

  let mut fi = 0
  while fi < n_features {
    let feature_name = if fi < input.feature_names.length() {
      input.feature_names[fi]
    } else {
      "Feature" + fi.to_string()
    }

    // Step 1: Kruskal-Wallis test
    // Build groups for this feature
    let kw_groups : Array[Array[Double]] = Array::new()
    let mut ci2 = 0
    while ci2 < n_classes {
      let group_vals : Array[Double] = Array::new()
      let mut j = 0
      while j < class_indices[ci2].length() {
        group_vals.push(input.features[fi][class_indices[ci2][j]])
        j = j + 1
      }
      kw_groups.push(group_vals)
      ci2 = ci2 + 1
    }

    let (h_stat, kw_p) = lefse_kruskal_wallis(kw_groups)

    // If KW test is not significant, skip this feature
    if kw_p >= alpha {
      results.push(LEfSeResult::new(feature_name, "", kw_p, 0.0, false))
      fi = fi + 1
      continue
    }

    // Step 2: Subclass Wilcoxon tests (if applicable)
    let wilcoxon_pass = if has_subclasses {
      lefse_subclass_wilcoxon(input.features[fi], subclass_info, alpha)
    } else {
      true
    }

    if !wilcoxon_pass {
      results.push(LEfSeResult::new(feature_name, "", kw_p, 0.0, false))
      fi = fi + 1
      continue
    }

    // Step 3: Compute class means and LDA score
    let class_means : Array[Double] = Array::new()
    let mut ci3 = 0
    while ci3 < n_classes {
      let group_vals : Array[Double] = Array::new()
      let mut j = 0
      while j < class_indices[ci3].length() {
        group_vals.push(input.features[fi][class_indices[ci3][j]])
        j = j + 1
      }
      class_means.push(lefse_mean(group_vals))
      ci3 = ci3 + 1
    }

    // Find the class with the highest mean
    let mut best_idx = 0
    let mut best_mean = class_means[0]
    let mut ci4 = 1
    while ci4 < n_classes {
      if class_means[ci4] > best_mean {
        best_mean = class_means[ci4]
        best_idx = ci4
      }
      ci4 = ci4 + 1
    }

    let lda = lefse_lda_score(class_means, best_idx)

    // Step 4: Threshold filtering
    let is_significant = kw_p < alpha && lda.abs() > lda_threshold

    let best_class = unique_classes[best_idx]

    results.push(
      LEfSeResult::new(feature_name, best_class, kw_p, lda, is_significant),
    )

    if is_significant {
      n_sig = n_sig + 1
    }

    fi = fi + 1
  }

  LEfSeOutput::new(results, n_sig)
}

// ============================================================================
// Subclass helpers
// ============================================================================

///|
/// Check if subclasses are meaningful (not all empty and not all identical).
fn lefse_has_subclasses(subclasses : Array[String]) -> Bool {
  if subclasses.length() == 0 {
    return false
  }
  let first = subclasses[0]
  // Check if all are empty or all are the same
  let mut all_empty = true
  let mut all_same = true
  let mut i = 0
  while i < subclasses.length() {
    if subclasses[i].length() > 0 {
      all_empty = false
    }
    if subclasses[i] != first {
      all_same = false
    }
    i = i + 1
  }
  // Subclasses are meaningful if not all empty and not all identical
  !all_empty && !all_same
}

///|
/// Subclass pair information for Wilcoxon testing.
/// Each entry is (subclass_label, class_label, sample_indices).
/// We group samples by (subclass, class) for pairwise Wilcoxon tests
/// across classes within each subclass.
fn lefse_build_subclass_info(
  classes : Array[String],
  subclasses : Array[String],
  unique_classes : Array[String],
) -> Array[(String, String, Array[Int])] {
  let n = classes.length()
  let info : Array[(String, String, Array[Int])] = Array::new()

  // Collect all unique subclass labels
  let unique_subclasses = lefse_unique_strings(subclasses)

  // For each subclass, for each class, collect sample indices
  let mut si = 0
  while si < unique_subclasses.length() {
    let sub = unique_subclasses[si]
    let mut ci = 0
    while ci < unique_classes.length() {
      let cls = unique_classes[ci]
      let indices : Array[Int] = Array::new()
      let mut j = 0
      while j < n {
        if subclasses[j] == sub && classes[j] == cls {
          indices.push(j)
        }
        j = j + 1
      }
      if indices.length() > 0 {
        info.push((sub, cls, indices))
      }
      ci = ci + 1
    }
    si = si + 1
  }

  info
}

///|
/// Perform pairwise Wilcoxon rank-sum tests within each subclass.
///
/// For each subclass, test all pairs of classes. A feature passes if it
/// is significant (p < alpha) in ALL pairwise class comparisons within
/// EVERY subclass that has samples in both classes of the pair.
fn lefse_subclass_wilcoxon(
  feature_values : Array[Double],
  subclass_info : Array[(String, String, Array[Int])],
  alpha : Double,
) -> Bool {
  // Group by subclass
  let unique_subclasses : Array[String] = Array::new()
  let mut i = 0
  while i < subclass_info.length() {
    let sub = subclass_info[i].0
    let mut found = false
    let mut j = 0
    while j < unique_subclasses.length() {
      if unique_subclasses[j] == sub {
        found = true
      }
      j = j + 1
    }
    if !found {
      unique_subclasses.push(sub)
    }
    i = i + 1
  }

  // Get unique classes
  let unique_classes : Array[String] = Array::new()
  i = 0
  while i < subclass_info.length() {
    let cls = subclass_info[i].1
    let mut found = false
    let mut j = 0
    while j < unique_classes.length() {
      if unique_classes[j] == cls {
        found = true
      }
      j = j + 1
    }
    if !found {
      unique_classes.push(cls)
    }
    i = i + 1
  }

  // For each subclass, for each pair of classes, do Wilcoxon test
  let mut si = 0
  while si < unique_subclasses.length() {
    let sub = unique_subclasses[si]

    // For each pair of classes
    let mut ci1 = 0
    while ci1 < unique_classes.length() {
      let mut ci2 = ci1 + 1
      while ci2 < unique_classes.length() {
        // Get sample indices for this (subclass, class) combination
        let indices1 = lefse_get_subclass_class_indices(
          subclass_info,
          sub,
          unique_classes[ci1],
        )
        let indices2 = lefse_get_subclass_class_indices(
          subclass_info,
          sub,
          unique_classes[ci2],
        )

        // Only test if both groups have samples
        if indices1.length() >= 1 && indices2.length() >= 1 {
          // Extract feature values
          let vals1 : Array[Double] = Array::new()
          let mut k = 0
          while k < indices1.length() {
            vals1.push(feature_values[indices1[k]])
            k = k + 1
          }
          let vals2 : Array[Double] = Array::new()
          k = 0
          while k < indices2.length() {
            vals2.push(feature_values[indices2[k]])
            k = k + 1
          }

          // Need at least 3 samples per group for meaningful Wilcoxon
          if vals1.length() >= 3 && vals2.length() >= 3 {
            let (_, p_val) = lefse_wilcoxon_rank_sum(vals1, vals2)
            // Feature must be significant in ALL pairwise tests
            if p_val >= alpha {
              return false
            }
          }
        }

        ci2 = ci2 + 1
      }
      ci1 = ci1 + 1
    }
    si = si + 1
  }

  true
}

///|
/// Get sample indices for a given (subclass, class) pair from subclass_info.
fn lefse_get_subclass_class_indices(
  subclass_info : Array[(String, String, Array[Int])],
  subclass : String,
  class_ : String,
) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let mut i = 0
  while i < subclass_info.length() {
    if subclass_info[i].0 == subclass && subclass_info[i].1 == class_ {
      let indices = subclass_info[i].2
      let mut j = 0
      while j < indices.length() {
        result.push(indices[j])
        j = j + 1
      }
    }
    i = i + 1
  }
  result
}

// ============================================================================
// Summary
// ============================================================================

///|
/// Generate a text summary of LEfSe analysis results.
pub fn lefse_result_summary(output : LEfSeOutput) -> String {
  let sb = StringBuilder::new()
  sb.write_string("=== LEfSe Analysis Results ===\n")
  sb.write_string("Total features tested: ")
  sb.write_string(output.results.length().to_string())
  sb.write_char('\n')
  sb.write_string("Significant biomarkers: ")
  sb.write_string(output.n_significant.to_string())
  sb.write_char('\n')
  sb.write_string("--- Per-feature results ---\n")
  sb.write_string("feature\tclass\tKW_pvalue\tLDA\tsignificant\n")
  let mut i = 0
  while i < output.results.length() {
    let r = output.results[i]
    sb.write_string(r.feature_name)
    sb.write_char('\t')
    sb.write_string(r.class_name)
    sb.write_char('\t')
    sb.write_string(r.kw_pvalue.to_string())
    sb.write_char('\t')
    sb.write_string(r.lda_score.to_string())
    sb.write_char('\t')
    sb.write_string(r.is_significant.to_string())
    sb.write_char('\n')
    i = i + 1
  }
  sb.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Build sample LEfSe input data: 6 features (taxa) × 10 samples
/// in 2 classes (5 "Healthy" + 5 "Disease"), with 2 subclasses
/// (subjects A and B). Features 0 and 1 are truly differential
/// (elevated in Disease), while the rest are similar across classes.
pub fn lefse_sample_data() -> LEfSeInput {
  let feature_names = [
    "Bacteroides", "Prevotella", "Faecalibacterium", "Roseburia", "Eubacterium",
    "Ruminococcus",
  ]
  let sample_names = [
    "H_A1", "H_A2", "H_B1", "H_B2", "H_B3", "D_A1", "D_A2", "D_B1", "D_B2", "D_B3",
  ]
  let classes = [
    "Healthy", "Healthy", "Healthy", "Healthy", "Healthy", "Disease", "Disease",
    "Disease", "Disease", "Disease",
  ]
  let subclasses = [
    "SubjectA", "SubjectA", "SubjectB", "SubjectB", "SubjectB", "SubjectA", "SubjectA",
    "SubjectB", "SubjectB", "SubjectB",
  ]
  // Features: rows = taxa, cols = samples
  // Features 0,1 are elevated in Disease; features 2-5 are similar.
  let features = [
    [120.0, 110.0, 130.0, 115.0, 125.0, 250.0, 240.0, 260.0, 245.0, 255.0], // Bacteroides
    [80.0, 85.0, 75.0, 90.0, 82.0, 180.0, 170.0, 190.0, 175.0, 185.0], // Prevotella
    [50.0, 55.0, 48.0, 52.0, 51.0, 55.0, 50.0, 58.0, 53.0, 54.0], // Faecalibacterium
    [30.0, 35.0, 28.0, 32.0, 31.0, 33.0, 30.0, 35.0, 31.0, 32.0], // Roseburia
    [40.0, 42.0, 38.0, 41.0, 40.0, 43.0, 40.0, 45.0, 42.0, 41.0], // Eubacterium
    [25.0, 28.0, 22.0, 26.0, 24.0, 27.0, 25.0, 29.0, 24.0, 26.0], // Ruminococcus
  ]
  LEfSeInput::new(features, feature_names, sample_names, classes, subclasses)
}
