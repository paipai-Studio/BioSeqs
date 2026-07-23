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
