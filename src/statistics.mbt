///|
/// Bio.Statistics module - Statistical distributions and tests.
/// Inspired by Biopython's Bio.Statistics module.
/// Provides probability distributions, hypothesis tests, and correlation analysis.

///|
/// Calculate the mean of an array.
pub fn stat_mean(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < data.length() {
    sum = sum + data[i]
    i = i + 1
  }
  sum / data.length().to_double()
}

///|
/// Calculate the variance of an array (population variance).
pub fn stat_variance(data : Array[Double]) -> Double {
  if data.length() <= 1 {
    return 0.0
  }
  let mean = stat_mean(data)
  let mut sum_sq = 0.0
  let mut i = 0
  while i < data.length() {
    let diff = data[i] - mean
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  sum_sq / (data.length().to_double() - 1.0)
}

///|
/// Calculate the standard deviation.
pub fn stat_std(data : Array[Double]) -> Double {
  stat_variance(data).sqrt()
}

///|
/// Calculate the median of an array.
pub fn stat_median(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let sorted = data.copy()
  sorted.sort()
  let n = sorted.length()
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

///|
/// Calculate the mode (most frequent value) of an array.
pub fn stat_mode(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let counts : Map[Double, Int] = Map([], capacity=data.length())
  let mut i = 0
  while i < data.length() {
    let val = data[i]
    let current = counts.get_or_default(val, 0)
    counts.set(val, current + 1)
    i = i + 1
  }
  
  let mut max_count = 0
  let mut mode_val = data[0]
  for key in counts.keys() {
    let count = counts.get_or_default(key, 0)
    if count > max_count {
      max_count = count
      mode_val = key
    }
  }
  mode_val
}

///|
/// Calculate the minimum value.
pub fn stat_min(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let mut min_val = data[0]
  let mut i = 1
  while i < data.length() {
    if data[i] < min_val {
      min_val = data[i]
    }
    i = i + 1
  }
  min_val
}

///|
/// Calculate the maximum value.
pub fn stat_max(data : Array[Double]) -> Double {
  if data.length() == 0 {
    return 0.0
  }
  let mut max_val = data[0]
  let mut i = 1
  while i < data.length() {
    if data[i] > max_val {
      max_val = data[i]
    }
    i = i + 1
  }
  max_val
}

///|
/// Calculate the range (max - min).
pub fn stat_range(data : Array[Double]) -> Double {
  stat_max(data) - stat_min(data)
}

///|
/// Calculate the sum of an array.
pub fn stat_sum(data : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < data.length() {
    sum = sum + data[i]
    i = i + 1
  }
  sum
}

///|
/// Calculate the cumulative sum.
pub fn stat_cumsum(data : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut cumsum = 0.0
  let mut i = 0
  while i < data.length() {
    cumsum = cumsum + data[i]
    result.push(cumsum)
    i = i + 1
  }
  result
}

///|
/// Calculate the dot product of two arrays.
pub fn stat_dot_product(a : Array[Double], b : Array[Double]) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + a[i] * b[i]
    i = i + 1
  }
  sum
}

///|
/// Calculate the Pearson correlation coefficient.
pub fn stat_pearson_correlation(x : Array[Double], y : Array[Double]) -> Double {
  let n = if x.length() < y.length() { x.length() } else { y.length() }
  if n < 3 {
    return 0.0
  }
  
  let mean_x = stat_mean(x)
  let mean_y = stat_mean(y)
  
  let mut num = 0.0
  let mut den_x = 0.0
  let mut den_y = 0.0
  let mut i = 0
  
  while i < n {
    let dx = x[i] - mean_x
    let dy = y[i] - mean_y
    num = num + dx * dy
    den_x = den_x + dx * dx
    den_y = den_y + dy * dy
    i = i + 1
  }
  
  let denom = (den_x * den_y).sqrt()
  if denom == 0.0 {
    0.0
  } else {
    num / denom
  }
}

///|
/// Calculate the Spearman rank correlation coefficient.
pub fn stat_spearman_correlation(x : Array[Double], y : Array[Double]) -> Double {
  let ranks_x = stat_rank(x)
  let ranks_y = stat_rank(y)
  stat_pearson_correlation(ranks_x, ranks_y)
}

///|
/// Calculate ranks for an array (average ranking for ties).
pub fn stat_rank(data : Array[Double]) -> Array[Double] {
  let n = data.length()
  let result : Array[Double] = Array::make(n, 0.0)
  
  let indexed : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row = [data[i], i.to_double()]
    indexed.push(row)
    i = i + 1
  }
  
  indexed.sort_by(fn(a: Array[Double], b: Array[Double]) -> Int { if a[0] < b[0] { -1 } else if a[0] > b[0] { 1 } else { 0 } })
  
  let mut pos = 0
  while pos < n {
    let mut end = pos
    while end < n - 1 && indexed[end + 1][0] == indexed[pos][0] {
      end = end + 1
    }
    let avg_rank = ((pos + end).to_double() / 2.0) + 1.0
    let mut j = pos
    while j <= end {
      let idx = indexed[j][1].to_int()
      result[idx] = avg_rank
      j = j + 1
    }
    pos = end + 1
  }
  
  result
}

///|
/// Calculate the Z-score (standardized value).
pub fn stat_zscore(value : Double, mean : Double, std : Double) -> Double {
  if std == 0.0 {
    0.0
  } else {
    (value - mean) / std
  }
}

///|
/// Calculate the Z-scores for an array.
pub fn stat_zscores(data : Array[Double]) -> Array[Double] {
  let mean = stat_mean(data)
  let std = stat_std(data)
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < data.length() {
    result.push(stat_zscore(data[i], mean, std))
    i = i + 1
  }
  result
}

///|
/// Calculate the p-value from Z-score (two-tailed).
pub fn stat_pvalue_zscore(z : Double) -> Double {
  let abs_z = if z < 0.0 { -z } else { z }
  2.0 * (1.0 - stat_normal_cdf(abs_z))
}

///|
/// Approximate normal CDF using Abramowitz and Stegun approximation.
pub fn stat_normal_cdf(x : Double) -> Double {
  let t = 1.0 / (1.0 + 0.2316419 * (if x < 0.0 { -x } else { x }))
  let d = 0.3989423 * @math.exp(-x * x / 2.0)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.8212560 + t * 1.3302744))))
  if x > 0.0 {
    1.0 - p
  } else {
    p
  }
}

///|
/// Calculate the t-statistic for a one-sample t-test.
pub fn stat_t_statistic(sample : Array[Double], mu_0 : Double) -> Double {
  let n = sample.length()
  if n < 2 {
    return 0.0
  }
  let mean = stat_mean(sample)
  let std = stat_std(sample)
  if std == 0.0 {
    0.0
  } else {
    (mean - mu_0) / (std / (n.to_double()).sqrt())
  }
}

///|
/// Calculate the p-value from t-distribution (approximate).
pub fn stat_pvalue_t(t : Double, df : Int) -> Double {
  let abs_t = if t < 0.0 { -t } else { t }
  // Approximation using normal distribution for large df
  let z = abs_t * (1.0 - 1.0 / (4.0 * df.to_double())) / (1.0 + t * t / (df.to_double())).sqrt()
  2.0 * (1.0 - stat_normal_cdf(z))
}

///|
/// Calculate confidence interval for the mean.
pub fn stat_confidence_interval(data : Array[Double], confidence? : Double = 0.95) -> Array[Double] {
  let n = data.length()
  if n < 2 {
    return [stat_mean(data), stat_mean(data)]
  }
  let mean = stat_mean(data)
  let std = stat_std(data)
  
  let z = if confidence == 0.99 {
    2.576
  } else if confidence == 0.90 {
    1.645
  } else {
    1.96
  }
  
  let margin = z * std / (n.to_double()).sqrt()
  [mean - margin, mean + margin]
}

///|
/// Calculate the chi-square statistic.
pub fn stat_chi_square(observed : Array[Double], expected : Array[Double]) -> Double {
  let n = if observed.length() < expected.length() { observed.length() } else { expected.length() }
  let mut chi2 = 0.0
  let mut i = 0
  while i < n {
    if expected[i] > 0.0 {
      let diff = observed[i] - expected[i]
      chi2 = chi2 + diff * diff / expected[i]
    }
    i = i + 1
  }
  chi2
}

///|
/// Expected frequencies for chi-square test (uniform distribution).
pub fn stat_uniform_expected(n_categories : Int, total : Double) -> Array[Double] {
  let expected : Array[Double] = Array::make(n_categories, total / n_categories.to_double())
  expected
}

///|
/// Calculate the frequency of each unique value.
pub fn stat_frequency_counts(data : Array[Double]) -> Map[Double, Int] {
  let counts : Map[Double, Int] = Map([], capacity=data.length())
  let mut i = 0
  while i < data.length() {
    let val = data[i]
    let current = counts.get_or_default(val, 0)
    counts.set(val, current + 1)
    i = i + 1
  }
  counts
}

///|
/// Get unique values in an array.
pub fn stat_unique_values(data : Array[Double]) -> Array[Double] {
  let map : Map[Double, Bool] = Map([], capacity=data.length())
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < data.length() {
    let val = data[i]
    let is_present = map.get_or_default(val, false)
    if !is_present {
      map.set(val, true)
      result.push(val)
    }
    i = i + 1
  }
  result
}

///|
/// Calculate the natural log (ln) of a number.
fn stat_ln(x : Double) -> Double {
  if x <= 0.0 { return 0.0 }
  if x == 1.0 { return 0.0 }
  let y = (x - 1.0) / (x + 1.0)
  let y2 = y * y
  let mut sum = y
  let mut term = y
  let mut i = 1
  while i <= 20 {
    term = term * y2
    sum = sum + term / (2 * i + 1).to_double()
    i = i + 1
  }
  2.0 * sum
}

///|
/// Calculate the exponential of a number.
fn stat_exp(x : Double) -> Double {
  let mut sum = 1.0
  let mut term = 1.0
  let mut i = 1
  while i <= 20 {
    term = term * x / i.to_double()
    sum = sum + term
    i = i + 1
  }
  sum
}

///|
/// Calculate the geometric mean of an array (only positive values).
pub fn stat_geometric_mean(data : Array[Double]) -> Double {
  let mut sum_log = 0.0
  let mut n_valid = 0
  for val in data {
    if val > 0.0 {
      sum_log = sum_log + stat_ln(val)
      n_valid = n_valid + 1
    }
  }
  if n_valid > 0 {
    stat_exp(sum_log / n_valid.to_double())
  } else {
    0.0
  }
}

///|
/// Sort an array of (score, index) pairs by score descending (quicksort).
pub fn stat_sort_indices_desc(data : Array[(Double, Int)]) -> Array[(Double, Int)] {
  if data.length() <= 1 {
    return data.copy()
  }
  let pivot = data[0]
  let left : Array[(Double, Int)] = Array::new()
  let right : Array[(Double, Int)] = Array::new()
  let mut i = 1
  while i < data.length() {
    if data[i].0 >= pivot.0 {
      left.push(data[i])
    } else {
      right.push(data[i])
    }
    i = i + 1
  }
  let result = stat_sort_indices_desc(left)
  result.push(pivot)
  let right_sorted = stat_sort_indices_desc(right)
  let mut j = 0
  while j < right_sorted.length() {
    result.push(right_sorted[j])
    j = j + 1
  }
  result
}

///|
/// Sort an array of (score, index) pairs by score ascending (quicksort).
pub fn stat_sort_indices_asc(data : Array[(Double, Int)]) -> Array[(Double, Int)] {
  if data.length() <= 1 {
    return data.copy()
  }
  let pivot = data[0]
  let left : Array[(Double, Int)] = Array::new()
  let right : Array[(Double, Int)] = Array::new()
  let mut i = 1
  while i < data.length() {
    if data[i].0 <= pivot.0 {
      left.push(data[i])
    } else {
      right.push(data[i])
    }
    i = i + 1
  }
  let result = stat_sort_indices_asc(left)
  result.push(pivot)
  let right_sorted = stat_sort_indices_asc(right)
  let mut j = 0
  while j < right_sorted.length() {
    result.push(right_sorted[j])
    j = j + 1
  }
  result
}

///|
/// In-place quicksort for Double arrays (ascending).
pub fn stat_sort_double(data : Array[Double]) -> Array[Double] {
  let result = data.copy()
  result.sort()
  result
}

///|
/// Compute mean and variance in a single pass.
pub fn stat_mean_var(data : Array[Double]) -> (Double, Double) {
  let n = data.length()
  if n == 0 {
    return (0.0, 0.0)
  }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + data[i]
    i = i + 1
  }
  let mean = sum / n.to_double()
  let mut sum_sq = 0.0
  i = 0
  while i < n {
    let diff = data[i] - mean
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  let var = if n > 1 { sum_sq / (n.to_double() - 1.0) } else { 0.0 }
  (mean, var)
}

///|
/// Standardize an array to zero mean and unit variance (Z-score normalization).
pub fn stat_standardize(data : Array[Double]) -> Array[Double] {
  let (mean, var) = stat_mean_var(data)
  let std = if var > 0.0 { var.sqrt() } else { 1.0 }
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < data.length() {
    result.push((data[i] - mean) / std)
    i = i + 1
  }
  result
}

///|
/// Wilcoxon signed-rank test for paired samples.
/// Returns (statistic, p_value).
pub fn stat_wilcoxon_signed_rank(x : Array[Double], y : Array[Double]) -> (Double, Double) {
  let n = if x.length() < y.length() { x.length() } else { y.length() }
  if n < 6 {
    return (0.0, 1.0)
  }

  let diffs : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let d = x[i] - y[i]
    if d.abs() > 1.0e-15 {
      diffs.push(d)
    }
    i = i + 1
  }

  let n_valid = diffs.length()
  if n_valid < 6 {
    return (0.0, 1.0)
  }

  // Rank absolute differences
  let abs_diffs : Array[Double] = Array::new()
  for d in diffs {
    abs_diffs.push(d.abs())
  }
  let ranks = stat_rank(abs_diffs)

  // Compute W (sum of ranks of positive differences)
  let mut w_plus = 0.0
  i = 0
  while i < n_valid {
    if diffs[i] > 0.0 {
      w_plus = w_plus + ranks[i]
    }
    i = i + 1
  }

  // Mean and variance under null
  let n_dbl = n_valid.to_double()
  let mean_w = n_dbl * (n_dbl + 1.0) / 4.0
  let var_w = n_dbl * (n_dbl + 1.0) * (2.0 * n_dbl + 1.0) / 24.0

  // Z-statistic with continuity correction
  let z = if w_plus > mean_w {
    (w_plus - mean_w - 0.5) / var_w.sqrt()
  } else if w_plus < mean_w {
    (w_plus - mean_w + 0.5) / var_w.sqrt()
  } else {
    0.0
  }

  // Two-tailed p-value
  let p = 2.0 * stat_pvalue_zscore(z)
  let stat = if w_plus > (n_dbl * (n_dbl + 1.0) / 4.0) {
    w_plus
  } else {
    n_dbl * (n_dbl + 1.0) / 2.0 - w_plus
  }
  (stat, if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p })
}

///|
/// Mann-Whitney U test (Wilcoxon rank-sum test) for two independent samples.
/// Returns (statistic, p_value).
pub fn stat_mann_whitney(x : Array[Double], y : Array[Double]) -> (Double, Double) {
  let n1 = x.length()
  let n2 = y.length()
  if n1 < 3 || n2 < 3 {
    return (0.0, 1.0)
  }

  // Combine and rank
  let combined : Array[Double] = Array::new()
  for v in x {
    combined.push(v)
  }
  for v in y {
    combined.push(v)
  }
  let ranks = stat_rank(combined)

  let n = n1 + n2
  let mut r1 = 0.0
  let mut i = 0
  while i < n1 {
    r1 = r1 + ranks[i]
    i = i + 1
  }

  // U statistic
  let u1 = r1 - n1.to_double() * (n1.to_double() + 1.0) / 2.0
  let u2 = (n1.to_double() * n2.to_double()) - u1
  let u = if u1 < u2 { u1 } else { u2 }

  // Normal approximation with tie correction
  let n1_d = n1.to_double()
  let n2_d = n2.to_double()
  let mean_u = n1_d * n2_d / 2.0
  let var_u = n1_d * n2_d * (n1_d + n2_d + 1.0) / 12.0

  let z = if var_u > 0.0 {
    (u - mean_u) / var_u.sqrt()
  } else {
    0.0
  }

  let p = 2.0 * stat_pvalue_zscore(z)
  (u, if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p })
}

///|
/// Fisher's exact test for a 2x2 contingency table.
/// Table format: [[a, b], [c, d]]
/// Returns (odds_ratio, p_value).
pub fn stat_fisher_exact(table : Array[Array[Int]]) -> (Double, Double) {
  let a = if table.length() > 0 && table[0].length() > 0 { table[0][0] } else { 0 }
  let b = if table.length() > 0 && table[0].length() > 1 { table[0][1] } else { 0 }
  let c = if table.length() > 1 && table[1].length() > 0 { table[1][0] } else { 0 }
  let d = if table.length() > 1 && table[1].length() > 1 { table[1][1] } else { 0 }

  let n = (a + b + c + d).to_double()
  if n <= 0.0 {
    return (1.0, 1.0)
  }

  // Odds ratio
  let or = if c.to_double() * b.to_double() == 0.0 {
    if a.to_double() * d.to_double() == 0.0 { 1.0 } else { 1.0e15 }
  } else {
    (a.to_double() * d.to_double()) / (b.to_double() * c.to_double())
  }

  // Hypergeometric probability
  let p_val = stat_fisher_exact_pvalue(a, b, c, d)
  (or, p_val)
}

///|
/// Compute Fisher's exact test p-value using hypergeometric distribution.
fn stat_fisher_exact_pvalue(a : Int, b : Int, c : Int, d : Int) -> Double {
  let n = (a + b + c + d).to_double()
  let a_d = a.to_double()
  let b_d = b.to_double()
  let c_d = c.to_double()
  let d_d = d.to_double()

  // Probability of the observed table
  let p_obs = stat_mult_hypergeom_pmf(a, b, c, d)

  // Enumerate all tables with same margins
  let row1_total = a + b
  let row2_total = c + d
  let col1_total = a + c
  let col2_total = b + d

  let mut p_value = 0.0
  let min_a = if (row1_total - col2_total) > 0 { row1_total - col2_total } else { 0 }
  let max_a = if row1_total < col1_total { row1_total } else { col1_total }

  let mut a_curr = min_a
  while a_curr <= max_a {
    let b_curr = row1_total - a_curr
    let c_curr = col1_total - a_curr
    let d_curr = col2_total - b_curr
    if c_curr >= 0 && d_curr >= 0 {
      let p_curr = stat_mult_hypergeom_pmf(a_curr, b_curr, c_curr, d_curr)
      if p_curr <= p_obs + 1.0e-10 {
        p_value = p_value + p_curr
      }
    }
    a_curr = a_curr + 1
  }

  if p_value > 1.0 { 1.0 } else if p_value < 0.0 { 0.0 } else { p_value }
}

///|
/// Multivariate hypergeometric PMF for a 2x2 table.
fn stat_mult_hypergeom_pmf(a : Int, b : Int, c : Int, d : Int) -> Double {
  let n = (a + b + c + d).to_double()
  if n <= 0.0 { return 0.0 }

  let a_d = a.to_double()
  let b_d = b.to_double()
  let c_d = c.to_double()
  let d_d = d.to_double()

  let row1 = a_d + b_d
  let row2 = c_d + d_d
  let col1 = a_d + c_d
  let col2 = b_d + d_d

  // log probability using gamma functions
  let ln_p = stat_ln_factorial(row1) + stat_ln_factorial(row2) +
    stat_ln_factorial(col1) + stat_ln_factorial(col2) -
    stat_ln_factorial(n) -
    stat_ln_factorial(a_d) - stat_ln_factorial(b_d) -
    stat_ln_factorial(c_d) - stat_ln_factorial(d_d)

  let p = @math.exp(ln_p)
  if p.is_nan() || p < 0.0 { 0.0 } else { p }
}

///|
/// Natural log of factorial (uses Stirling's approximation for large values).
fn stat_ln_factorial(x : Double) -> Double {
  if x <= 0.0 { return 0.0 }
  if x < 20.0 {
    // Direct computation for small values
    let mut result = 0.0
    let mut i = 1.0
    while i <= x {
      result = result + @math.ln(i)
      i = i + 1.0
    }
    result
  } else {
    // Stirling's approximation
    let pi = 3.141592653589793
    let e = 2.718281828459045
    let x1 = x + 1.0
    (x1 - 0.5) * @math.ln(x1) - x1 + 0.5 * @math.ln(2.0 * pi)
  }
}

///|
/// Kolmogorov-Smirnov test for comparing two samples.
/// Returns (statistic, p_value).
pub fn stat_ks_test(sample1 : Array[Double], sample2 : Array[Double]) -> (Double, Double) {
  let n1 = sample1.length()
  let n2 = sample2.length()
  if n1 < 2 || n2 < 2 {
    return (0.0, 1.0)
  }

  // Combine all values and compute ECDF
  let combined : Array[Double] = Array::new()
  for v in sample1 { combined.push(v) }
  for v in sample2 { combined.push(v) }
  let sorted = stat_sort_double(combined)
  let n = sorted.length()

  // Compute ECDF for each sample
  let mut max_diff = 0.0
  let mut i = 0
  while i < n {
    let x = sorted[i]
    // ECDF for sample1 at x
    let mut ecdf1 = 0.0
    for v in sample1 {
      if v <= x { ecdf1 = ecdf1 + 1.0 }
    }
    ecdf1 = ecdf1 / n1.to_double()

    // ECDF for sample2 at x
    let mut ecdf2 = 0.0
    for v in sample2 {
      if v <= x { ecdf2 = ecdf2 + 1.0 }
    }
    ecdf2 = ecdf2 / n2.to_double()

    let diff = (ecdf1 - ecdf2).abs()
    if diff > max_diff {
      max_diff = diff
    }
    i = i + 1
  }

  // Approximate p-value using Kolmogorov approximation
  let n_eff = (n1.to_double() * n2.to_double()) / (n1.to_double() + n2.to_double())
  let sqrt_n_eff = n_eff.sqrt()
  let lambda = (sqrt_n_eff + 0.12 + 0.11 / sqrt_n_eff) * max_diff

  // Kolmogorov distribution approximation
  let p = if lambda <= 0.0 {
    1.0
  } else {
    2.0 * stat_kolmogorov_cdf(lambda)
  }

  (max_diff, if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p })
}

///|
/// Kolmogorov distribution CDF (approximation).
fn stat_kolmogorov_cdf(lambda : Double) -> Double {
  if lambda <= 0.0 { return 0.0 }
  if lambda >= 10.0 { return 1.0 }

  let mut sum = 0.0
  let mut k = 1
  while k <= 100 {
    let k_dbl = k.to_double()
    let term = (-1.0).pow(k_dbl - 1.0) * @math.exp(-2.0 * k_dbl * k_dbl * lambda * lambda)
    sum = sum + term
    if term.abs() < 1.0e-15 { break }
    k = k + 1
  }
  2.0 * sum
}

///|
/// Bonferroni correction for multiple testing.
pub fn stat_bonferroni(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let corrected = p_values[i] * n.to_double()
    result.push(if corrected > 1.0 { 1.0 } else if corrected < 0.0 { 0.0 } else { corrected })
    i = i + 1
  }
  result
}

///|
/// Holm-Bonferroni step-down procedure for multiple testing.
pub fn stat_holm(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 { return Array::new() }

  // Create indexed pairs and sort by p-value
  let indexed : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((p_values[i], i))
    i = i + 1
  }

  // Sort ascending
  let sorted = stat_sort_indices_asc(indexed)

  // Apply Holm-Bonferroni
  let adjusted : Array[Double] = Array::make(n, 1.0)
  let mut cummin = 1.0
  let mut j = n - 1
  while j >= 0 {
    let rank = n - j
    let orig_idx = sorted[j].1
    let p = sorted[j].0
    let corrected = p * rank.to_double()
    let val = if corrected > 1.0 { 1.0 } else { corrected }
    cummin = if cummin < val { cummin } else { val }
    adjusted[orig_idx] = cummin
    j = j - 1
  }

  adjusted
}

///|
/// Benjamini-Yekutieli correction for multiple testing (conservative version).
pub fn stat_by(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 { return Array::new() }

  // Harmonic number for BY correction
  let mut h = 0.0
  let mut k = 1
  while k <= n {
    h = h + 1.0 / k.to_double()
    k = k + 1
  }

  let indexed : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((p_values[i], i))
    i = i + 1
  }

  let sorted = stat_sort_indices_asc(indexed)
  let adjusted : Array[Double] = Array::make(n, 1.0)
  let mut cummin = 1.0

  let mut j = n - 1
  while j >= 0 {
    let rank = n - j
    let orig_idx = sorted[j].1
    let p = sorted[j].0
    let corrected = p * (n.to_double() * h) / rank.to_double()
    let val = if corrected > 1.0 { 1.0 } else { corrected }
    cummin = if cummin < val { cummin } else { val }
    adjusted[orig_idx] = cummin
    j = j - 1
  }

  adjusted
}

///|
/// Beta distribution PDF.
pub fn stat_beta_pdf(x : Double, alpha : Double, beta : Double) -> Double {
  if x <= 0.0 || x >= 1.0 {
    return 0.0
  }
  if alpha <= 0.0 || beta <= 0.0 {
    return 0.0
  }

  // log(1/B(alpha,beta)) + (alpha-1)*log(x) + (beta-1)*log(1-x)
  let ln_beta = stat_lngamma(alpha) + stat_lngamma(beta) - stat_lngamma(alpha + beta)
  let ln_x = @math.ln(x)
  let ln_1mx = @math.ln(1.0 - x)
  let ln_pdf = -ln_beta + (alpha - 1.0) * ln_x + (beta - 1.0) * ln_1mx

  let pdf = @math.exp(ln_pdf)
  if pdf.is_nan() || pdf < 0.0 { 0.0 } else { pdf }
}

///|
/// Gamma distribution PDF.
pub fn stat_gamma_pdf(x : Double, shape : Double, rate : Double) -> Double {
  if x <= 0.0 || shape <= 0.0 || rate <= 0.0 {
    return 0.0
  }

  let ln_gamma = stat_lngamma(shape)
  let ln_pdf = shape * @math.ln(rate) - ln_gamma + (shape - 1.0) * @math.ln(x) - rate * x

  let pdf = @math.exp(ln_pdf)
  if pdf.is_nan() || pdf < 0.0 { 0.0 } else { pdf }
}

///|
/// Compute ln(gamma(x)) using Lanczos approximation.
fn stat_lngamma(x : Double) -> Double {
  if x <= 0.0 {
    return 1.0e300
  }

  let cof = [
    76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
    0.001208650973866179, -0.000005395239384953,
  ]

  let mut x1 = x
  let mut tmp = x1 + 5.5
  tmp = (x1 + 0.5) * @math.ln(tmp) - tmp

  let mut ser = 1.000000000190015
  let mut j = 0
  while j < cof.length() {
    x1 = x1 + 1.0
    ser = ser + cof[j] / x1
    j = j + 1
  }

  tmp + @math.ln(2.5066282746310005 * ser)
}

///|
/// Bartlett's test for homogeneity of variances across groups.
/// Returns (statistic, p_value).
pub fn stat_bartlett_test(samples : Array[Array[Double]]) -> (Double, Double) {
  let k = samples.length()
  if k < 2 {
    return (0.0, 1.0)
  }

  let variances : Array[Double] = Array::new()
  let sizes : Array[Int] = Array::new()
  let mut total_n = 0
  let mut i = 0
  while i < k {
    let n = samples[i].length()
    if n < 2 { return (0.0, 1.0) }
    let var = stat_variance(samples[i])
    variances.push(var)
    sizes.push(n)
    total_n = total_n + n
    i = i + 1
  }

  // Pooled variance
  let mut pooled_var_num = 0.0
  let mut pooled_var_den = 0.0
  i = 0
  while i < k {
    let n = sizes[i].to_double()
    pooled_var_num = pooled_var_num + (n - 1.0) * variances[i]
    pooled_var_den = pooled_var_den + n - 1.0
    i = i + 1
  }
  let pooled_var = pooled_var_num / pooled_var_den

  // Compute B statistic
  let mut b = 0.0
  i = 0
  while i < k {
    let n = sizes[i].to_double()
    let vi = variances[i]
    if vi > 0.0 {
      b = b + (n - 1.0) * (@math.ln(pooled_var) - @math.ln(vi))
    }
    i = i + 1
  }
  b = 2.0 * b

  // Approximate p-value using chi-square distribution
  let df = k - 1
  let p = stat_chi2_sf(b, df)
  (b, if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p })
}

///|
/// Chi-square survival function (p-value from chi-square distribution).
fn stat_chi2_sf(x : Double, df : Int) -> Double {
  if x <= 0.0 { return 1.0 }
  if df <= 0 { return 1.0 }

  // Wilson-Hilferty approximation
  let df_d = df.to_double()
  let z = ((x / df_d).pow(1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df_d))) / ((2.0 / (9.0 * df_d)).sqrt())
  stat_pvalue_zscore(z)
}

///|
/// Levene's test for homogeneity of variances.
/// Returns (statistic, p_value).
pub fn stat_levene_test(samples : Array[Array[Double]]) -> (Double, Double) {
  let k = samples.length()
  if k < 2 {
    return (0.0, 1.0)
  }

  // Compute group medians
  let medians : Array[Double] = Array::new()
  let sizes : Array[Int] = Array::new()
  let all_deviations : Array[Array[Double]] = Array::new()
  let mut total_n = 0

  let mut i = 0
  while i < k {
    let n = samples[i].length()
    if n < 2 { return (0.0, 1.0) }
    let med = stat_median(samples[i])
    medians.push(med)
    sizes.push(n)
    total_n = total_n + n

    // Absolute deviations from median
    let devs : Array[Double] = Array::new()
    for v in samples[i] {
      devs.push((v - med).abs())
    }
    all_deviations.push(devs)
    i = i + 1
  }

  // Compute overall mean of deviations
  let all_devs_flat : Array[Double] = Array::new()
  for devs in all_deviations {
    for d in devs {
      all_devs_flat.push(d)
    }
  }
  let grand_mean = stat_mean(all_devs_flat)

  // Compute group means and variances of deviations
  let group_means : Array[Double] = Array::new()
  let group_vars : Array[Double] = Array::new()
  i = 0
  while i < k {
    group_means.push(stat_mean(all_deviations[i]))
    group_vars.push(stat_variance(all_deviations[i]))
    i = i + 1
  }

  // Levene statistic: ratio of between-group to within-group variation
  let mut numerator = 0.0
  i = 0
  while i < k {
    let n = sizes[i].to_double()
    numerator = numerator + n * (group_means[i] - grand_mean).pow(2.0)
    i = i + 1
  }

  let mut denominator = 0.0
  i = 0
  while i < k {
    let n = sizes[i].to_double()
    denominator = denominator + (n - 1.0) * group_vars[i]
    i = i + 1
  }

  let levene_stat = if denominator > 0.0 {
    (numerator / (k.to_double() - 1.0)) / (denominator / (total_n.to_double() - k.to_double()))
  } else {
    0.0
  }

  // Approximate p-value using F distribution
  let df1 = k - 1
  let df2 = total_n - k
  let p = stat_f_sf(levene_stat, df1, df2)
  (levene_stat, if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p })
}

///|
/// F-distribution survival function (p-value from F distribution).
fn stat_f_sf(f : Double, df1 : Int, df2 : Int) -> Double {
  if f <= 0.0 { return 1.0 }
  if df1 <= 0 || df2 <= 0 { return 1.0 }

  // Use incomplete beta function relationship
  let x = df1.to_double() * f / (df1.to_double() * f + df2.to_double())
  if x >= 1.0 { return 0.0 }
  let a = df1.to_double() / 2.0
  let b = df2.to_double() / 2.0
  let p = stat_beta_inc(x, a, b)
  if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p }
}

///|
/// Regularized incomplete beta function I_x(a, b).
fn stat_beta_inc(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 { return 0.0 }
  if x >= 1.0 { return 1.0 }

  // Use continued fraction expansion
  let lbeta = stat_lngamma(a) + stat_lngamma(b) - stat_lngamma(a + b)
  let front = @math.exp(@math.ln(x) * a + @math.ln(1.0 - x) * b - lbeta) / a

  // Lentz's continued fraction algorithm
  let mut f = 1.0
  let mut c = 1.0
  let mut d = 0.0
  let max_iter = 200
  let eps = 1.0e-10

  let mut i = 0
  while i < max_iter {
    let m = i / 2
    let m_d = m.to_double()

    let mut numerator : Double = 0.0
    if i == 0 {
      numerator = 1.0
    } else if i % 2 == 0 {
      // Even step
      numerator = m_d * (b - m_d) * x / ((a + 2.0 * m_d - 1.0) * (a + 2.0 * m_d))
    } else {
      // Odd step
      numerator = -(a + m_d) * (a + b + m_d) * x / ((a + 2.0 * m_d) * (a + 2.0 * m_d + 1.0))
    }

    d = 1.0 + numerator * d
    if d.abs() < eps { d = eps }
    d = 1.0 / d

    c = 1.0 + numerator / c
    if c.abs() < eps { c = eps }

    let delta = c * d
    f = f * delta
    if (delta - 1.0).abs() < eps { break }

    i = i + 1
  }

  front * (f - 1.0)
}

///|
/// Chi-square distribution CDF.
pub fn stat_chi2_cdf(x : Double, df : Int) -> Double {
  if x <= 0.0 { return 0.0 }
  if df <= 0 { return 0.0 }

  let df_d = df.to_double()
  // Wilson-Hilferty approximation
  let z = ((x / df_d).pow(1.0 / 3.0) - (1.0 - 2.0 / (9.0 * df_d))) / ((2.0 / (9.0 * df_d)).sqrt())
  let p = stat_normal_cdf(z)
  if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p }
}

///|
/// Chi-square distribution quantile (inverse CDF).
pub fn stat_chi2_quantile(p : Double, df : Int) -> Double {
  if p <= 0.0 { return 0.0 }
  if p >= 1.0 { return 1.0e300 }
  if df <= 0 { return 0.0 }

  // Newton's method on stat_chi2_cdf(x, df) = p
  let df_d = df.to_double()
  let mut x = df_d * (1.0 - 2.0 / (9.0 * df_d))
  let mut i = 0
  while i < 100 {
    let cdf_val = stat_chi2_cdf(x, df)
    let diff = cdf_val - p
    if diff.abs() < 1.0e-12 { break }
    // Numerical derivative
    let h = 1.0e-7
    let deriv = (stat_chi2_cdf(x + h, df) - stat_chi2_cdf(x - h, df)) / (2.0 * h)
    if deriv.abs() < 1.0e-30 { break }
    x = x - diff / deriv
    if x < 0.0 { x = 1.0e-15 }
    i = i + 1
  }
  if x < 0.0 { 0.0 } else { x }
}

///|
/// Student's t-distribution CDF.
pub fn stat_t_cdf(t : Double, df : Double) -> Double {
  if df <= 0.0 { return 0.5 }
  let x = df / (df + t * t)
  let p = 0.5 * stat_beta_inc(x, df / 2.0, 0.5)
  if t >= 0.0 { 1.0 - p } else { p }
}

///|
/// Student's t-distribution quantile (inverse CDF).
pub fn stat_t_quantile(p : Double, df : Double) -> Double {
  if p <= 0.0 { return -1.0e300 }
  if p >= 1.0 { return 1.0e300 }
  if (p - 0.5).abs() < 1.0e-15 { return 0.0 }

  // Newton's method on stat_t_cdf(t, df) = p
  let mut t = if p < 0.5 { -1.0 } else { 1.0 }
  let mut i = 0
  while i < 100 {
    let cdf_val = stat_t_cdf(t, df)
    let diff = cdf_val - p
    if diff.abs() < 1.0e-12 { break }
    // Numerical derivative
    let h = 1.0e-7
    let deriv = (stat_t_cdf(t + h, df) - stat_t_cdf(t - h, df)) / (2.0 * h)
    if deriv.abs() < 1.0e-30 { break }
    t = t - diff / deriv
    i = i + 1
  }
  t
}

///|
/// One-sample t-test.
/// Returns (statistic, p_value).
pub fn stat_t_test_one_sample(sample : Array[Double], mu_0 : Double) -> (Double, Double) {
  let n = sample.length()
  if n < 3 {
    return (0.0, 1.0)
  }
  let mean = stat_mean(sample)
  let std = stat_std(sample)
  if std <= 0.0 {
    return (0.0, 1.0)
  }
  let t = (mean - mu_0) / (std / (n.to_double()).sqrt())
  let df = (n - 1).to_double()
  let p = 2.0 * (1.0 - stat_t_cdf(t.abs(), df))
  (t, if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p })
}

///|
/// Two-sample t-test (equal variances).
/// Returns (statistic, p_value).
pub fn stat_t_test_two_sample(x : Array[Double], y : Array[Double]) -> (Double, Double) {
  let n1 = x.length()
  let n2 = y.length()
  if n1 < 3 || n2 < 3 {
    return (0.0, 1.0)
  }
  let mean1 = stat_mean(x)
  let mean2 = stat_mean(y)
  let var1 = stat_variance(x)
  let var2 = stat_variance(y)

  // Pooled variance (equal variances assumption)
  let pooled_var = ((n1 - 1).to_double() * var1 + (n2 - 1).to_double() * var2) /
    ((n1 + n2 - 2).to_double())
  if pooled_var <= 0.0 {
    return (0.0, 1.0)
  }
  let se = (pooled_var * (1.0 / n1.to_double() + 1.0 / n2.to_double())).sqrt()
  if se <= 0.0 {
    return (0.0, 1.0)
  }
  let t = (mean1 - mean2) / se
  let df = (n1 + n2 - 2).to_double()
  let p = 2.0 * (1.0 - stat_t_cdf(t.abs(), df))
  (t, if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p })
}

///|
/// Benjamini-Hochberg (BH) step-up procedure for multiple testing.
pub fn stat_bh(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 { return Array::new() }

  let indexed : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    let p = if p_values[i].is_nan() { 1.0 } else { p_values[i] }
    indexed.push((p, i))
    i = i + 1
  }

  // Sort ascending
  let sorted = stat_sort_indices_asc(indexed)

  let adjusted : Array[Double] = Array::make(n, 1.0)
  let mut cummin = 1.0

  let mut j = n - 1
  while j >= 0 {
    let rank = (n - j).to_double()
    let orig_idx = sorted[j].1
    let p = sorted[j].0
    let corrected = p * n.to_double() / rank
    let val = if corrected > 1.0 { 1.0 } else { corrected }
    cummin = if cummin < val { cummin } else { val }
    adjusted[orig_idx] = cummin
    j = j - 1
  }

  i = 0
  while i < n {
    if adjusted[i] > 1.0 { adjusted[i] = 1.0 }
    i = i + 1
  }

  adjusted
}

///|
/// Normal distribution quantile (inverse CDF) using rational approximation.
pub fn stat_normal_quantile(p : Double) -> Double {
  if p <= 0.0 { return -1.0e300 }
  if p >= 1.0 { return 1.0e300 }

  // Peter Acklam's rational approximation
  let a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
    1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
  let b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
    6.680131188771972e+01, -1.328068155288572e+01]
  let c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
    -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
  let d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
    3.754408661907416e+00]

  let plow = 0.02425
  let phigh = 1.0 - plow

  if p < plow {
    let q = (-2.0 * @math.ln(p)).sqrt()
    c[0] * (q + c[1] * (q + c[2] * (q + c[3] * (q + c[4] * q + c[5])))) /
      (1.0 + d[0] * (q + d[1] * (q + d[2] * (q + d[3] * q))))
  } else if p <= phigh {
    let q = p - 0.5
    let r = q * q
    (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
      (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
  } else {
    let q = (-2.0 * @math.ln(1.0 - p)).sqrt()
    -c[0] * (q + c[1] * (q + c[2] * (q + c[3] * (q + c[4] * q + c[5])))) /
      (1.0 + d[0] * (q + d[1] * (q + d[2] * (q + d[3] * q))))
  }
}

///|
/// Compute the mean absolute deviation (MAD) from median.
pub fn stat_mad(data : Array[Double], constant? : Double = 1.4826) -> Double {
  if data.length() == 0 { return 0.0 }
  let med = stat_median(data)
  let devs : Array[Double] = Array::new()
  let mut i = 0
  while i < data.length() {
    devs.push((data[i] - med).abs())
    i = i + 1
  }
  stat_median(devs) * constant
}

///|
/// Log-rank p-value approximation (for two groups).
/// Returns (chi2_statistic, p_value).
pub fn stat_logrank_test(
  time : Array[Double],
  event : Array[Bool],
  group : Array[Int],
) -> (Double, Double) {
  let n = time.length()
  if n < 3 {
    return (0.0, 1.0)
  }

  // Sort by time
  let indexed : Array[(Double, Bool, Int)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((time[i], event[i], group[i]))
    i = i + 1
  }

  // Simple ordering
  let mut swapped = true
  let mut k = n - 1
  while k > 0 && swapped {
    swapped = false
    let mut j = 0
    while j < k {
      if indexed[j].0 > indexed[j + 1].0 {
        let tmp = indexed[j]
        indexed[j] = indexed[j + 1]
        indexed[j + 1] = tmp
        swapped = true
      }
      j = j + 1
    }
    k = k - 1
  }

  let mut n1_total = 0
  let mut n2_total = 0
  i = 0
  while i < n {
    if indexed[i].2 == 0 { n1_total = n1_total + 1 } else { n2_total = n2_total + 1 }
    i = i + 1
  }

  let mut n1 = n1_total.to_double()
  let mut n2 = n2_total.to_double()
  let mut o1 = 0.0
  let mut o2 = 0.0
  let mut e1 = 0.0
  let mut v = 0.0

  i = 0
  while i < n {
    let t_i = indexed[i].0
    let d = if indexed[i].1 { 1.0 } else { 0.0 }
    let at_risk = n1 + n2

    if at_risk > 0.0 {
      e1 = e1 + d * n1 / at_risk
      let v_i = d * n1 * n2 / (at_risk * at_risk * (at_risk - 1.0))
      v = v + v_i
    }

    if indexed[i].2 == 0 {
      o1 = o1 + d
    } else {
      o2 = o2 + d
    }

    n1 = n1 - 1.0
    n2 = n2 - 1.0
    i = i + 1
  }

  let z = if v > 0.0 {
    (o1 - e1) / v.sqrt()
  } else {
    0.0
  }

  let p = 2.0 * (1.0 - stat_normal_cdf(z.abs()))
  let chi2 = z * z
  (chi2, if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p })
}

///|
/// Analysis of Variance (one-way ANOVA).
/// Returns (f_statistic, p_value).
pub fn stat_anova(samples : Array[Array[Double]]) -> (Double, Double) {
  let k = samples.length()
  if k < 2 {
    return (0.0, 1.0)
  }

  let sizes : Array[Int] = Array::new()
  let means : Array[Double] = Array::new()
  let variances : Array[Double] = Array::new()
  let mut total_n = 0

  let mut i = 0
  while i < k {
    let n = samples[i].length()
    if n < 2 { return (0.0, 1.0) }
    sizes.push(n)
    means.push(stat_mean(samples[i]))
    variances.push(stat_variance(samples[i]))
    total_n = total_n + n
    i = i + 1
  }

  // Overall mean
  let all_vals : Array[Double] = Array::new()
  i = 0
  while i < k {
    for v in samples[i] {
      all_vals.push(v)
    }
    i = i + 1
  }
  let grand_mean = stat_mean(all_vals)

  // Between-group sum of squares
  let mut ss_between = 0.0
  i = 0
  while i < k {
    ss_between = ss_between + sizes[i].to_double() * (means[i] - grand_mean).pow(2.0)
    i = i + 1
  }

  // Within-group sum of squares
  let mut ss_within = 0.0
  i = 0
  while i < k {
    let m = means[i]
    for v in samples[i] {
      ss_within = ss_within + (v - m).pow(2.0)
    }
    i = i + 1
  }

  let df_between = (k - 1).to_double()
  let df_within = (total_n - k).to_double()

  if df_between <= 0.0 || df_within <= 0.0 || ss_within <= 0.0 {
    return (0.0, 1.0)
  }

  let ms_between = ss_between / df_between
  let ms_within = ss_within / df_within
  let f = ms_between / ms_within

  let p = stat_f_sf(f, k - 1, total_n - k)
  (f, if p.is_nan() || p < 0.0 { 0.0 } else if p > 1.0 { 1.0 } else { p })
}
