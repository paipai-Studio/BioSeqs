///|
/// scikit-bio Stats module - Statistical functions for bioinformatics.

// ============================================================
// Helper: Int seed to 32-byte Bytes for chacha8 RNG
// ============================================================

///|
fn seed_to_bytes(seed : Int) -> Bytes {
  Bytes::makei(32, fn(i) { ((seed >> (i % 4 * 8)) & 0xFF).to_byte() })
}

// ============================================================
// Subsample counts
// ============================================================

///|
/// Randomly subsample from a vector of counts, with or without replacement.
/// Returns a subsampled vector where the sum equals n.
pub fn subsample_counts(
  counts : Array[Double],
  n : Int,
  replace : Bool,
  seed : Int,
) -> Array[Double] {
  let total = sum_double(counts)
  let total_int = total.to_int()

  if n < 0 {
    return Array::make(counts.length(), 0.0)
  }
  if !replace && n > total_int {
    return Array::make(counts.length(), 0.0)
  }

  let rng = @random.Rand::chacha8(seed=seed_to_bytes(seed))
  let result = Array::make(counts.length(), 0.0)

  if replace {
    let total_d = total
    let mut i = 0
    while i < n {
      let r = rng.double() * total_d
      let mut cumulative = 0.0
      let mut idx = 0
      while idx < counts.length() && cumulative + counts[idx] < r {
        cumulative = cumulative + counts[idx]
        idx = idx + 1
      }
      if idx >= counts.length() {
        idx = counts.length() - 1
      }
      result[idx] = result[idx] + 1.0
      i = i + 1
    }
  } else {
    let pool : Array[Int] = Array::new()
    let mut i = 0
    while i < counts.length() {
      let cnt = counts[i].to_int()
      let mut j = 0
      while j < cnt {
        pool.push(i)
        j = j + 1
      }
      i = i + 1
    }

    // Fisher-Yates shuffle
    let m = pool.length()
    let mut i = 0
    while i < m {
      let j = i + rng.int(limit=m - i)
      let tmp = pool[i]
      pool[i] = pool[j]
      pool[j] = tmp
      i = i + 1
    }

    let take = if n < m { n } else { m }
    let mut i = 0
    while i < take {
      let idx = pool[i]
      result[idx] = result[idx] + 1.0
      i = i + 1
    }
  }

  result
}

///|
fn sum_double(arr : Array[Double]) -> Double {
  let mut s = 0.0
  for x in arr {
    s = s + x
  }
  s
}

// ============================================================
// P-value correction (multiple comparison)
// ============================================================

///|
/// Bonferroni correction.
/// p_adj = p * m, clipped to [0, 1]
pub fn p_adjust_bonferroni(pvalues : Array[Double]) -> Array[Double] {
  let m = pvalues.length().to_double()
  let result = Array::make(pvalues.length(), 0.0)
  let mut i = 0
  while i < pvalues.length() {
    let adj = pvalues[i] * m
    result[i] = if adj > 1.0 { 1.0 } else { adj }
    i = i + 1
  }
  result
}

///|
/// Holm-Bonferroni (step-down) correction.
pub fn p_adjust_holm(pvalues : Array[Double]) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 {
    return Array::new()
  }

  let indexed : Array[(Int, Double)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((i, pvalues[i]))
    i = i + 1
  }

  sort_pairs_by_p(indexed, 0, n - 1)

  let adjusted_sorted = Array::make(n, 0.0)
  let mut prev = 0.0
  let mut i = 0
  while i < n {
    let factor = (n - i).to_double()
    let mut adj = indexed[i].1 * factor
    if adj < prev {
      adj = prev
    }
    if adj > 1.0 {
      adj = 1.0
    }
    adjusted_sorted[i] = adj
    prev = adj
    i = i + 1
  }

  let result = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    result[indexed[i].0] = adjusted_sorted[i]
    i = i + 1
  }
  result
}

///|
/// Benjamini-Hochberg (FDR) correction.
pub fn p_adjust_fdr_bh(pvalues : Array[Double]) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 {
    return Array::new()
  }

  let indexed : Array[(Int, Double)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((i, pvalues[i]))
    i = i + 1
  }

  sort_pairs_by_p(indexed, 0, n - 1)

  let adjusted_sorted = Array::make(n, 0.0)
  let mut prev = 1.0
  let mut i = n - 1
  while i >= 0 {
    let factor = n.to_double() / (i + 1).to_double()
    let mut adj = indexed[i].1 * factor
    if adj > prev {
      adj = prev
    }
    if adj > 1.0 {
      adj = 1.0
    }
    adjusted_sorted[i] = adj
    prev = adj
    i = i - 1
  }

  let result = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    result[indexed[i].0] = adjusted_sorted[i]
    i = i + 1
  }
  result
}

///|
/// Quick sort for (index, pvalue) pairs by p-value ascending.
fn sort_pairs_by_p(arr : Array[(Int, Double)], lo : Int, hi : Int) -> Unit {
  if lo < hi {
    let pivot = arr[hi].1
    let mut i = lo - 1
    let mut j = lo
    while j < hi {
      if arr[j].1 <= pivot {
        i = i + 1
        let tmp = arr[i]
        arr[i] = arr[j]
        arr[j] = tmp
      }
      j = j + 1
    }
    let tmp = arr[i + 1]
    arr[i + 1] = arr[hi]
    arr[hi] = tmp
    let pi = i + 1
    sort_pairs_by_p(arr, lo, pi - 1)
    sort_pairs_by_p(arr, pi + 1, hi)
  }
}

// ============================================================
// Permutation test
// ============================================================

///|
/// Result of a permutation test.
pub struct PermutationTestResult {
  statistic : Double
  p_value : Double
  permutations : Int
}

///|
/// Perform a permutation test on two independent samples.
/// p-value is proportion of permuted stats >= original stat (one-sided, greater).
pub fn permutation_test(
  sample1 : Array[Double],
  sample2 : Array[Double],
  test_fn : (Array[Double], Array[Double]) -> Double,
  permutations : Int,
  seed : Int,
) -> PermutationTestResult {
  let orig_stat = test_fn(sample1, sample2)

  if permutations <= 0 {
    return PermutationTestResult::{
      statistic: orig_stat,
      p_value: 0.0 / 0.0,
      permutations: 0,
    }
  }

  let combined : Array[Double] = Array::new()
  for x in sample1 {
    combined.push(x)
  }
  for x in sample2 {
    combined.push(x)
  }

  let n1 = sample1.length()
  let n_total = combined.length()

  let rng = @random.Rand::chacha8(seed=seed_to_bytes(seed))
  let mut count_ge = 0

  let mut i = 0
  while i < permutations {
    let shuffled = combined.copy()
    let mut j = 0
    while j < n_total {
      let k = j + rng.int(limit=n_total - j)
      let tmp = shuffled[j]
      shuffled[j] = shuffled[k]
      shuffled[k] = tmp
      j = j + 1
    }

    let s1 = Array::make(n1, 0.0)
    let s2 = Array::make(n_total - n1, 0.0)
    let mut j = 0
    while j < n1 {
      s1[j] = shuffled[j]
      j = j + 1
    }
    let mut j = n1
    while j < n_total {
      s2[j - n1] = shuffled[j]
      j = j + 1
    }

    let perm_stat = test_fn(s1, s2)
    if perm_stat >= orig_stat {
      count_ge = count_ge + 1
    }
    i = i + 1
  }

  let p_value = count_ge.to_double() / permutations.to_double()

  PermutationTestResult::{ statistic: orig_stat, p_value, permutations }
}

///|
/// Mean of a double array.
fn mean_double(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    0.0
  } else {
    sum_double(arr) / arr.length().to_double()
  }
}

///|
/// Difference of means test statistic (mean1 - mean2).
pub fn diff_of_means(s1 : Array[Double], s2 : Array[Double]) -> Double {
  mean_double(s1) - mean_double(s2)
}

// ============================================================
// Distance matrix utilities (for ANOSIM / PERMANOVA)
// ============================================================

///|
/// A distance matrix (square, symmetric, zero diagonal).
pub struct DistanceMatrix {
  data : Array[Double]
  size : Int
  ids : Array[String]
}

///|
/// Create a DistanceMatrix from a flat array (row-major order).
pub fn DistanceMatrix::from_flat(
  data : Array[Double],
  size : Int,
  ids : Array[String],
) -> DistanceMatrix {
  DistanceMatrix::{ data, size, ids }
}

///|
/// Create a DistanceMatrix from a 2D array.
pub fn DistanceMatrix::from_2d(
  matrix : Array[Array[Double]],
  ids : Array[String],
) -> DistanceMatrix {
  let n = matrix.length()
  let data = Array::make(n * n, 0.0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      data[i * n + j] = matrix[i][j]
      j = j + 1
    }
    i = i + 1
  }
  DistanceMatrix::{ data, size: n, ids }
}

///|
/// Get distance between sample i and sample j.
pub fn DistanceMatrix::get(self : DistanceMatrix, i : Int, j : Int) -> Double {
  self.data[i * self.size + j]
}

// ============================================================
// ANOSIM (Analysis of Similarities)
// ============================================================

///|
/// Result of an ANOSIM test.
pub struct ANOSIMResult {
  r_statistic : Double
  p_value : Double
  permutations : Int
  num_groups : Int
  sample_size : Int
}

///|
/// Perform ANOSIM (Analysis of Similarities) test.
/// groups: array of group labels (integers), one per sample.
pub fn anosim(
  dm : DistanceMatrix,
  groups : Array[Int],
  permutations : Int,
  seed : Int,
) -> ANOSIMResult {
  let n = dm.size
  let num_groups = count_unique_int(groups)
  let orig_r = anosim_r_statistic(dm, groups)

  if permutations <= 0 {
    return ANOSIMResult::{
      r_statistic: orig_r,
      p_value: 0.0 / 0.0,
      permutations: 0,
      num_groups,
      sample_size: n,
    }
  }

  let rng = @random.Rand::chacha8(seed=seed_to_bytes(seed))
  let mut count_ge = 0

  let mut i = 0
  while i < permutations {
    let perm_groups = groups.copy()
    let mut j = 0
    while j < n {
      let k = j + rng.int(limit=n - j)
      let tmp = perm_groups[j]
      perm_groups[j] = perm_groups[k]
      perm_groups[k] = tmp
      j = j + 1
    }

    let perm_r = anosim_r_statistic(dm, perm_groups)
    if perm_r >= orig_r {
      count_ge = count_ge + 1
    }
    i = i + 1
  }

  let p_value = count_ge.to_double() / permutations.to_double()

  ANOSIMResult::{
    r_statistic: orig_r,
    p_value,
    permutations,
    num_groups,
    sample_size: n,
  }
}

///|
/// Compute ANOSIM R statistic.
fn anosim_r_statistic(dm : DistanceMatrix, groups : Array[Int]) -> Double {
  let n = dm.size
  let num_pairs = n * (n - 1) / 2

  let distances = Array::make(num_pairs, 0.0)
  let is_between_arr = Array::make(num_pairs, false)

  let mut idx = 0
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      distances[idx] = dm.get(i, j)
      is_between_arr[idx] = groups[i] != groups[j]
      idx = idx + 1
      j = j + 1
    }
    i = i + 1
  }

  let ranks = rank_data(distances)

  let mut sum_within = 0.0
  let mut count_within = 0.0
  let mut sum_between = 0.0
  let mut count_between = 0.0

  let mut i = 0
  while i < num_pairs {
    if is_between_arr[i] {
      sum_between = sum_between + ranks[i]
      count_between = count_between + 1.0
    } else {
      sum_within = sum_within + ranks[i]
      count_within = count_within + 1.0
    }
    i = i + 1
  }

  let r_b = if count_between > 0.0 { sum_between / count_between } else { 0.0 }
  let r_w = if count_within > 0.0 { sum_within / count_within } else { 0.0 }

  let denom = n.to_double() * (n.to_double() - 1.0) / 4.0
  if denom == 0.0 {
    0.0
  } else {
    (r_b - r_w) / denom
  }
}

///|
/// Compute ranks of data (average ranks for ties).
fn rank_data(data : Array[Double]) -> Array[Double] {
  let n = data.length()
  let result = Array::make(n, 0.0)

  let indexed : Array[(Int, Double)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((i, data[i]))
    i = i + 1
  }

  sort_pairs_by_p(indexed, 0, n - 1)

  let mut i = 0
  while i < n {
    let mut j = i
    while j + 1 < n && indexed[j + 1].1 == indexed[i].1 {
      j = j + 1
    }
    let avg_rank = (i.to_double() + j.to_double()) / 2.0 + 1.0
    let mut k = i
    while k <= j {
      result[indexed[k].0] = avg_rank
      k = k + 1
    }
    i = j + 1
  }

  result
}

///|
/// Count unique integer values.
fn count_unique_int(arr : Array[Int]) -> Int {
  let seen : @hashmap.HashMap[Int, Bool] = @hashmap.HashMap([])
  for x in arr {
    if !seen.contains(x) {
      seen.set(x, true)
    }
  }
  seen.length()
}

// ============================================================
// PERMANOVA (Permutational Multivariate Analysis of Variance)
// ============================================================

///|
/// Result of a PERMANOVA test.
pub struct PERMANOVAResult {
  f_statistic : Double
  p_value : Double
  permutations : Int
  num_groups : Int
  sample_size : Int
  df_between : Int
  df_within : Int
  ss_between : Double
  ss_within : Double
  ss_total : Double
}

///|
/// Perform PERMANOVA test using pseudo-F statistic on distance matrix.
pub fn permanova(
  dm : DistanceMatrix,
  groups : Array[Int],
  permutations : Int,
  seed : Int,
) -> PERMANOVAResult {
  let n = dm.size
  let num_groups = count_unique_int(groups)

  let (f_stat, ss_bet, ss_with, ss_tot, df_bet, df_with) = permanova_f_statistic(
    dm, groups,
  )

  if permutations <= 0 {
    return PERMANOVAResult::{
      f_statistic: f_stat,
      p_value: 0.0 / 0.0,
      permutations: 0,
      num_groups,
      sample_size: n,
      df_between: df_bet,
      df_within: df_with,
      ss_between: ss_bet,
      ss_within: ss_with,
      ss_total: ss_tot,
    }
  }

  let rng = @random.Rand::chacha8(seed=seed_to_bytes(seed))
  let mut count_ge = 0

  let mut i = 0
  while i < permutations {
    let perm_groups = groups.copy()
    let mut j = 0
    while j < n {
      let k = j + rng.int(limit=n - j)
      let tmp = perm_groups[j]
      perm_groups[j] = perm_groups[k]
      perm_groups[k] = tmp
      j = j + 1
    }

    let (perm_f, _, _, _, _, _) = permanova_f_statistic(dm, perm_groups)
    if perm_f >= f_stat {
      count_ge = count_ge + 1
    }
    i = i + 1
  }

  let p_value = count_ge.to_double() / permutations.to_double()

  PERMANOVAResult::{
    f_statistic: f_stat,
    p_value,
    permutations,
    num_groups,
    sample_size: n,
    df_between: df_bet,
    df_within: df_with,
    ss_between: ss_bet,
    ss_within: ss_with,
    ss_total: ss_tot,
  }
}

///|
/// Compute PERMANOVA pseudo-F statistic.
fn permanova_f_statistic(
  dm : DistanceMatrix,
  groups : Array[Int],
) -> (Double, Double, Double, Double, Int, Int) {
  let n = dm.size

  let mut ss_total = 0.0
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n {
      let d = dm.get(i, j)
      ss_total = ss_total + d * d
      j = j + 1
    }
    i = i + 1
  }
  ss_total = ss_total / n.to_double()

  let group_ids : Array[Int] = Array::new()
  let group_seen : @hashmap.HashMap[Int, Bool] = @hashmap.HashMap([])
  for g in groups {
    if !group_seen.contains(g) {
      group_seen.set(g, true)
      group_ids.push(g)
    }
  }
  let num_groups = group_ids.length()

  let mut ss_within = 0.0
  for g in group_ids {
    let mut ss_g = 0.0
    let mut n_g = 0
    let mut i = 0
    while i < n {
      if groups[i] == g {
        let mut j = i + 1
        while j < n {
          if groups[j] == g {
            let d = dm.get(i, j)
            ss_g = ss_g + d * d
          }
          j = j + 1
        }
        n_g = n_g + 1
      }
      i = i + 1
    }
    if n_g > 0 {
      ss_within = ss_within + ss_g / n_g.to_double()
    }
  }

  let ss_between = ss_total - ss_within

  let df_between = num_groups - 1
  let df_within = n - num_groups

  let f_stat = if df_within > 0 && ss_within > 0.0 {
    ss_between / df_between.to_double() / (ss_within / df_within.to_double())
  } else {
    0.0
  }

  (f_stat, ss_between, ss_within, ss_total, df_between, df_within)
}

// ============================================================
// Confidence bound (normal approximation)
// ============================================================

///|
/// Calculate a confidence bound assuming a normal distribution.
/// Returns (lower, upper) confidence interval bounds.
pub fn confidence_bound(
  mean : Double,
  std : Double,
  n : Int,
  confidence : Double,
) -> (Double, Double) {
  if n <= 1 || std == 0.0 {
    (mean, mean)
  } else {
    let z = normal_ppf(0.5 + confidence / 2.0)
    let se = std / n.to_double().sqrt()
    (mean - z * se, mean + z * se)
  }
}

///|
/// Approximate inverse normal CDF (probit function).
/// Uses Beasley-Springer-Moro algorithm.
fn normal_ppf(p : Double) -> Double {
  if p <= 0.0 {
    return -1000000000000000000000000000000.0
  }
  if p >= 1.0 {
    return 1000000000000000000000000000000.0
  }

  let a0 = -39.69683028665376
  let a1 = 220.9460984245205
  let a2 = -275.9285104469687
  let a3 = 138.3577518672690
  let a4 = -30.66479806614716
  let a5 = 2.506628277459239
  let b0 = -54.47609879822406
  let b1 = 161.5858368580409
  let b2 = -155.6989798598866
  let b3 = 66.80131188771972
  let b4 = -13.28068155288572
  let c0 = -0.007784894002430293
  let c1 = -0.3223964580411365
  let c2 = -2.400758277161838
  let c3 = -2.549732539343734
  let c4 = 4.374664141464968
  let c5 = 2.938163982698783
  let d0 = 0.007784695709041462
  let d1 = 0.3224671290700398
  let d2 = 2.445134137142996
  let d3 = 3.754408661907416

  let p_low = 0.02425
  let p_high = 1.0 - p_low

  if p < p_low {
    let q = (-2.0 * @math.ln(p)).sqrt()
    (((((c0 * q + c1) * q + c2) * q + c3) * q + c4) * q + c5) /
    ((((d0 * q + d1) * q + d2) * q + d3) * q + 1.0)
  } else if p <= p_high {
    let q = p - 0.5
    let r = q * q
    (((((a0 * r + a1) * r + a2) * r + a3) * r + a4) * r + a5) *
    q /
    (((((b0 * r + b1) * r + b2) * r + b3) * r + b4) * r + 1.0)
  } else {
    let q = (-2.0 * @math.ln(1.0 - p)).sqrt()
    -(((((c0 * q + c1) * q + c2) * q + c3) * q + c4) * q + c5) /
    ((((d0 * q + d1) * q + d2) * q + d3) * q + 1.0)
  }
}
