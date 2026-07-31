///|
/// Hmisc - Harrell Miscellaneous Statistical Utilities for BioSeqs.
/// Inspired by Bioconductor's Hmisc package by Frank Harrell.
/// Provides correlation analysis, variable clustering, descriptive statistics,
/// and various statistical utility functions commonly used in bioinformatics.

///|
/// Correlation type: Pearson or Spearman.
pub enum HmiscCorType {
  HC_Pearson
  HC_Spearman
}

/// Convert correlation type to string.
pub fn hmisc_cor_type_to_string(ct : HmiscCorType) -> String {
  match ct {
    HC_Pearson => "Pearson"
    HC_Spearman => "Spearman"
  }
}

///|
/// Result of a correlation matrix computation.
pub struct HmiscCorrResult {
  /// Correlation matrix (row-major)
  matrix : Array[Array[Double]]
  /// P-value matrix
  p_values : Array[Array[Double]]
  /// Number of observations used
  n : Array[Array[Int]]
  /// Variable names
  names : Array[String]
  /// Correlation type
  cor_type : HmiscCorType
}

///|
/// Descriptive statistics for a single variable.
pub struct HmiscDescriptiveStats {
  /// Variable name
  name : String
  /// Number of non-missing observations
  n : Int
  /// Number of missing observations
  n_missing : Int
  /// Mean
  mean : Double
  /// Standard deviation
  sd : Double
  /// Min
  min : Double
  /// 1st quartile
  q1 : Double
  /// Median
  median : Double
  /// 3rd quartile
  q3 : Double
  /// Max
  max : Double
}

///|
/// Variable clustering result.
pub struct HmiscVarClusResult {
  /// Dendrogram-like merge steps
  merges : Array[Array[Int]]
  /// Heights at which clusters merge
  heights : Array[Double]
  /// Variable assignments to clusters
  cluster_assignments : Array[Int]
  /// Number of clusters
  n_clusters : Int
  /// Variable names
  names : Array[String]
}

///|
/// Somers' d result with confidence limits.
pub struct HmiscSomersResult {
  /// Somers' d value
  d : Double
  /// Lower confidence limit
  lower : Double
  /// Upper confidence limit
  upper : Double
  /// Number of observations
  n : Int
}

///|
/// Compute descriptive statistics for an array of values.
/// Missing values (NaN) are automatically excluded.
pub fn hmisc_describe(values : Array[Double], name? : String = "") -> HmiscDescriptiveStats {
  let valid : Array[Double] = Array::new()
  let mut n_missing = 0
  let mut i = 0
  while i < values.length() {
    if !values[i].is_nan() {
      valid.push(values[i])
    } else {
      n_missing = n_missing + 1
    }
    i = i + 1
  }
  let n = valid.length()
  if n == 0 {
    return HmiscDescriptiveStats::{
      name,
      n: 0,
      n_missing,
      mean: 0.0,
      sd: 0.0,
      min: 0.0,
      q1: 0.0,
      median: 0.0,
      q3: 0.0,
      max: 0.0,
    }
  }
  let mut sum = 0.0
  i = 0
  while i < n {
    sum = sum + valid[i]
    i = i + 1
  }
  let mean_val = sum / n.to_double()
  let mut sum_sq = 0.0
  i = 0
  while i < n {
    let diff = valid[i] - mean_val
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  let sd = if n > 1 { (sum_sq / (n.to_double() - 1.0)).sqrt() } else { 0.0 }
  let sorted = valid.copy()
  sorted.sort()
  let min_val = sorted[0]
  let max_val = sorted[n - 1]
  let q1 = hmisc_percentile(sorted, 0.25)
  let median_val = hmisc_percentile(sorted, 0.50)
  let q3 = hmisc_percentile(sorted, 0.75)
  HmiscDescriptiveStats::{
    name,
    n,
    n_missing,
    mean: mean_val,
    sd,
    min: min_val,
    q1,
    median: median_val,
    q3,
    max: max_val,
  }
}

///|
/// Compute the percentile of a sorted array using linear interpolation.
fn hmisc_percentile(sorted : Array[Double], p : Double) -> Double {
  let n = sorted.length()
  if n == 0 {
    return 0.0
  }
  if p <= 0.0 {
    return sorted[0]
  }
  if p >= 1.0 {
    return sorted[n - 1]
  }
  let idx = p * (n.to_double() - 1.0)
  let lo = idx.floor().to_int()
  let hi = idx.ceil().to_int()
  if lo == hi {
    return sorted[lo]
  }
  let frac = idx - lo.to_double()
  sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

///|
/// Compute Pearson correlation coefficient between two arrays.
pub fn hmisc_pearson(x : Array[Double], y : Array[Double]) -> Double {
  let n = x.length()
  if n != y.length() || n < 3 {
    return @double.not_a_number
  }
  let mean_x = hmisc_array_mean(x)
  let mean_y = hmisc_array_mean(y)
  let mut sum_xy = 0.0
  let mut sum_x2 = 0.0
  let mut sum_y2 = 0.0
  let mut i = 0
  while i < n {
    let dx = x[i] - mean_x
    let dy = y[i] - mean_y
    sum_xy = sum_xy + dx * dy
    sum_x2 = sum_x2 + dx * dx
    sum_y2 = sum_y2 + dy * dy
    i = i + 1
  }
  let denom = (sum_x2 * sum_y2).sqrt()
  if denom < 1.0e-15 {
    return @double.not_a_number
  }
  sum_xy / denom
}

///|
/// Compute Spearman rank correlation between two arrays.
pub fn hmisc_spearman(x : Array[Double], y : Array[Double]) -> Double {
  let n = x.length()
  if n != y.length() || n < 3 {
    return @double.not_a_number
  }
  let mut i = 0
  while i < n {
    if x[i].is_nan() || y[i].is_nan() {
      return @double.not_a_number
    }
    i = i + 1
  }
  let rx = hmisc_rank_array(x)
  let ry = hmisc_rank_array(y)
  hmisc_pearson(rx, ry)
}

///|
/// Compute ranks of an array (average ranks for ties).
fn hmisc_rank_array(data : Array[Double]) -> Array[Double] {
  let n = data.length()
  let indices = Array::make(n, 0)
  let mut i = 0
  while i < n {
    indices[i] = i
    i = i + 1
  }
  indices.sort_by(fn(a : Int, b : Int) -> Int { if data[a] < data[b] { -1 } else if data[a] > data[b] { 1 } else { 0 } })
  let ranks = Array::make(n, 0.0)
  i = 0
  while i < n {
    let mut j = i
    while j < n - 1 && data[indices[j]] == data[indices[j + 1]] {
      j = j + 1
    }
    let avg = ((i + j).to_double() / 2.0) + 1.0
    let mut k = i
    while k <= j {
      ranks[indices[k]] = avg
      k = k + 1
    }
    i = j + 1
  }
  ranks
}

///|
/// Compute a correlation matrix for multiple variables.
/// Missing values (NaN) are handled by pairwise deletion.
pub fn hmisc_rcorr(
  data : Array[Array[Double]],
  names? : Array[String] = [],
  cor_type? : HmiscCorType = HmiscCorType::HC_Pearson,
) -> HmiscCorrResult {
  let n_vars = data.length()
  let n_obs = if n_vars > 0 { data[0].length() } else { 0 }
  let var_names = if names.length() == n_vars {
    names
  } else {
    let default_names = Array::make(n_vars, "")
    let mut i = 0
    while i < n_vars {
      default_names[i] = "V\{(i + 1)}"
      i = i + 1
    }
    default_names
  }
  let matrix : Array[Array[Double]] = Array::new()
  let p_values : Array[Array[Double]] = Array::new()
  let n_matrix : Array[Array[Int]] = Array::new()
  let mut i = 0
  while i < n_vars {
    matrix.push(Array::make(n_vars, 0.0))
    p_values.push(Array::make(n_vars, 0.0))
    n_matrix.push(Array::make(n_vars, 0))
    i = i + 1
  }
  let mut i = 0
  while i < n_vars {
    let mut j = 0
    while j < n_vars {
      if i == j {
        matrix[i][j] = 1.0
        p_values[i][j] = 0.0
        n_matrix[i][j] = n_obs
      } else if i < j {
        let x : Array[Double] = Array::new()
        let y : Array[Double] = Array::new()
        let mut k = 0
        while k < n_obs {
          if !data[i][k].is_nan() && !data[j][k].is_nan() {
            x.push(data[i][k])
            y.push(data[j][k])
          }
          k = k + 1
        }
        let n_pair = x.length()
        let r = match cor_type {
          HmiscCorType::HC_Pearson => hmisc_pearson(x, y)
          HmiscCorType::HC_Spearman => hmisc_spearman(x, y)
        }
        let p = if n_pair > 2 && !r.is_nan() {
          let t_stat = r * ((n_pair.to_double() - 2.0) / (1.0 - r * r)).sqrt()
          2.0 * hmisc_t_dist_p((t_stat.abs()), n_pair - 2)
        } else {
          @double.not_a_number
        }
        matrix[i][j] = r
        matrix[j][i] = r
        p_values[i][j] = p
        p_values[j][i] = p
        n_matrix[i][j] = n_pair
        n_matrix[j][i] = n_pair
      }
      j = j + 1
    }
    i = i + 1
  }
  HmiscCorrResult::{
    matrix,
    p_values,
    n: n_matrix,
    names: var_names,
    cor_type,
  }
}

///|
/// Approximate two-tailed p-value from Student's t-distribution.
fn hmisc_t_dist_p(t : Double, df : Int) -> Double {
  if df <= 0 || t.is_nan() {
    return @double.not_a_number
  }
  let x = df.to_double() / (df.to_double() + t * t)
  hmisc_ibeta(df.to_double() / 2.0, 0.5, x)
}

///|
/// Incomplete beta function I_x(a, b).
fn hmisc_ibeta(a : Double, b : Double, x : Double) -> Double {
  if x < 0.0 || x > 1.0 {
    return 0.0
  }
  if x == 0.0 || x == 1.0 {
    return x
  }
  let result = hmisc_log_gamma(a + b) - hmisc_log_gamma(a) - hmisc_log_gamma(b)
  let front = @math.exp(result + a * @math.ln(x) + b * @math.ln(1.0 - x)) / a
  let max_iter = 200
  let eps = 1.0e-10
  let mut f = 1.0
  let mut c_val = 1.0
  let mut d = 1.0 - (a + b) * x / (a + 1.0)
  if d.abs() < 1.0e-30 {
    d = 1.0e-30
  }
  d = 1.0 / d
  f = d
  let mut m = 1
  while m <= max_iter {
    let num = m.to_double() * (b - m.to_double()) * x / ((a + 2.0 * m.to_double() - 1.0) * (a + 2.0 * m.to_double()))
    d = 1.0 + num * d
    if d.abs() < 1.0e-30 { d = 1.0e-30 }
    c_val = 1.0 + num / c_val
    if c_val.abs() < 1.0e-30 { c_val = 1.0e-30 }
    d = 1.0 / d
    f = f * d * c_val
    let num2 = -(a + m.to_double()) * (a + b + m.to_double()) * x / ((a + 2.0 * m.to_double()) * (a + 2.0 * m.to_double() + 1.0))
    d = 1.0 + num2 * d
    if d.abs() < 1.0e-30 { d = 1.0e-30 }
    c_val = 1.0 + num2 / c_val
    if c_val.abs() < 1.0e-30 { c_val = 1.0e-30 }
    d = 1.0 / d
    let del = d * c_val
    f = f * del
    if (del - 1.0).abs() < eps {
      break
    }
    m = m + 1
  }
  front * f
}

///|
/// Natural log of the gamma function using Lanczos approximation.
fn hmisc_log_gamma(x : Double) -> Double {
  if x <= 0.0 {
    return @double.not_a_number
  }
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
  let pi = 3.14159265358979323846
  if x < 0.5 {
    return @math.ln(pi / @math.sin(pi * x)) - hmisc_log_gamma(1.0 - x)
  }
  let mut sum = c[0]
  let mut i = 1
  while i < 9 {
    sum = sum + c[i] / (x + i.to_double() - 1.0)
    i = i + 1
  }
  let t = x + g + 0.5
  0.5 * @math.ln(2.0 * pi) + (t * @math.ln(t) - t) + @math.ln(sum)
}

///|
/// Compute variable clustering using hierarchical clustering.
pub fn hmisc_varclus(
  data : Array[Array[Double]],
  names? : Array[String] = [],
  min_cluster_size? : Int = 2,
) -> HmiscVarClusResult {
  let n_vars = data.length()
  let var_names = if names.length() == n_vars {
    names
  } else {
    let default = Array::make(n_vars, "")
    let mut i = 0
    while i < n_vars {
      default[i] = "V\{(i + 1)}"
      i = i + 1
    }
    default
  }
  if n_vars <= min_cluster_size {
    let assignments = Array::make(n_vars, 0)
    let mut i = 0
    while i < n_vars {
      assignments[i] = 0
      i = i + 1
    }
    return HmiscVarClusResult::{
      merges: Array::new(),
      heights: Array::new(),
      cluster_assignments: assignments,
      n_clusters: 1,
      names: var_names,
    }
  }
  let cor_result = hmisc_rcorr(data, names=var_names, cor_type=HmiscCorType::HC_Pearson)
  let assignments = Array::make(n_vars, 0)
  let mut i = 0
  while i < n_vars {
    assignments[i] = i
    i = i + 1
  }
  let merges : Array[Array[Int]] = Array::new()
  let heights : Array[Double] = Array::new()
  let mut current_n = n_vars
  let max_clusters = hmisc_max_cluster_count(n_vars, min_cluster_size)
  while current_n > max_clusters {
    let mut best_i = 0
    let mut best_j = 0
    let mut best_sim = -1.0
    let mut a = 0
    while a < current_n {
      let mut b = a + 1
      while b < current_n {
        let sim = hmisc_cluster_similarity(assignments, cor_result.matrix, a, b)
        if sim > best_sim {
          best_sim = sim
          best_i = a
          best_j = b
        }
        b = b + 1
      }
      a = a + 1
    }
    let merge_height = 1.0 - best_sim
    heights.push(merge_height)
    let mut j = 0
    while j < n_vars {
      if assignments[j] == best_j {
        assignments[j] = best_i
      } else if assignments[j] > best_j {
        assignments[j] = assignments[j] - 1
      }
      j = j + 1
    }
    let merge_info = [best_i, best_j]
    merges.push(merge_info)
    current_n = current_n - 1
  }
  let unique_clusters = hmisc_unique_values(assignments)
  let n_clusters = unique_clusters.length()
  let final_assignments = Array::make(n_vars, 0)
  let mut i = 0
  while i < n_vars {
    let val = assignments[i]
    let mut new_val = 0
    let mut k = 0
    while k < n_clusters {
      if unique_clusters[k] == val {
        new_val = k
        break
      }
      k = k + 1
    }
    final_assignments[i] = new_val
    i = i + 1
  }
  HmiscVarClusResult::{
    merges,
    heights,
    cluster_assignments: final_assignments,
    n_clusters,
    names: var_names,
  }
}

///|
/// Compute maximum number of clusters.
fn hmisc_max_cluster_count(n_vars : Int, min_size : Int) -> Int {
  if min_size <= 1 {
    return 1
  }
  let n_clusters = n_vars / min_size
  if n_clusters < 1 { 1 } else { n_clusters }
}

///|
/// Compute average correlation between two clusters.
fn hmisc_cluster_similarity(
  assignments : Array[Int],
  cor_matrix : Array[Array[Double]],
  cluster_a : Int,
  cluster_b : Int,
) -> Double {
  let members_a : Array[Int] = Array::new()
  let members_b : Array[Int] = Array::new()
  let mut i = 0
  while i < assignments.length() {
    if assignments[i] == cluster_a {
      members_a.push(i)
    } else if assignments[i] == cluster_b {
      members_b.push(i)
    }
    i = i + 1
  }
  if members_a.length() == 0 || members_b.length() == 0 {
    return -1.0
  }
  let mut sum = 0.0
  let mut count = 0.0
  let mut i = 0
  while i < members_a.length() {
    let mut j = 0
    while j < members_b.length() {
      let r = cor_matrix[members_a[i]][members_b[j]]
      if !r.is_nan() {
        sum = sum + r.abs()
        count = count + 1.0
      }
      j = j + 1
    }
    i = i + 1
  }
  if count > 0.0 { sum / count } else { 0.0 }
}

///|
/// Extract unique values from an array (sorted).
fn hmisc_unique_values(arr : Array[Int]) -> Array[Int] {
  let sorted = arr.copy()
  sorted.sort()
  let result = Array::new()
  let mut i = 0
  while i < sorted.length() {
    if i == 0 || sorted[i] != sorted[i - 1] {
      result.push(sorted[i])
    }
    i = i + 1
  }
  result
}

///|
/// Compute Somers' d statistic for comparing two groups.
pub fn hmisc_somers_d(group1 : Array[Double], group2 : Array[Double]) -> HmiscSomersResult {
  let n1 = group1.length()
  let n2 = group2.length()
  if n1 == 0 || n2 == 0 {
    return HmiscSomersResult::{
      d: @double.not_a_number,
      lower: @double.not_a_number,
      upper: @double.not_a_number,
      n: 0,
    }
  }
  let mut concordant = 0.0
  let mut discordant = 0.0
  let mut i = 0
  while i < n1 {
    let mut j = 0
    while j < n2 {
      if group1[i] > group2[j] {
        concordant = concordant + 1.0
      } else if group1[i] < group2[j] {
        discordant = discordant + 1.0
      }
      j = j + 1
    }
    i = i + 1
  }
  let total = (n1 * n2).to_double()
  let d = if total > 0.0 {
    (concordant - discordant) / total
  } else {
    0.0
  }
  let se = (4.0 * concordant * discordant / (total * total * total)).sqrt()
  HmiscSomersResult::{
    d,
    lower: d - 1.96 * se,
    upper: d + 1.96 * se,
    n: n1 + n2,
  }
}

///|
/// Impute missing values using column means.
pub fn hmisc_impute(data : Array[Array[Double]]) -> Array[Array[Double]] {
  let n_rows = data.length()
  if n_rows == 0 {
    return data
  }
  let n_cols = data[0].length()
  let col_means = Array::make(n_cols, 0.0)
  let mut col = 0
  while col < n_cols {
    let mut sum = 0.0
    let mut count = 0
    let mut row = 0
    while row < n_rows {
      if !data[row][col].is_nan() {
        sum = sum + data[row][col]
        count = count + 1
      }
      row = row + 1
    }
    col_means[col] = if count > 0 { sum / count.to_double() } else { 0.0 }
    col = col + 1
  }
  let result : Array[Array[Double]] = Array::new()
  let mut row = 0
  while row < n_rows {
    let new_row = Array::make(n_cols, 0.0)
    let mut col = 0
    while col < n_cols {
      if data[row][col].is_nan() {
        new_row[col] = col_means[col]
      } else {
        new_row[col] = data[row][col]
      }
      col = col + 1
    }
    result.push(new_row)
    row = row + 1
  }
  result
}

///|
/// Compute mean of an array (ignoring NaN).
fn hmisc_array_mean(data : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut count = 0
  let mut i = 0
  while i < data.length() {
    if !data[i].is_nan() {
      sum = sum + data[i]
      count = count + 1
    }
    i = i + 1
  }
  if count > 0 { sum / count.to_double() } else { 0.0 }
}

///|
/// Convert HmiscCorType to string.
fn hmisc_cor_type_string(t : HmiscCorType) -> String {
  match t {
    HmiscCorType::HC_Pearson => "Pearson"
    HmiscCorType::HC_Spearman => "Spearman"
  }
}

///|
/// Format correlation result as a string table.
pub fn hmisc_rcorr_summary(result : HmiscCorrResult) -> String {
  let n = result.names.length()
  let mut s = ""
  s = s + "Correlation Matrix (\{hmisc_cor_type_string(result.cor_type)})\n"
  s = s + "Variables: \{n}\n"
  s = s + "\n"
  s = s + "        "
  let mut j = 0
  while j < n {
    s = s + " \{result.names[j]}"
    j = j + 1
  }
  s = s + "\n"
  let mut i = 0
  while i < n {
    s = s + "\{result.names[i]} "
    let mut j = 0
    while j < n {
      let val = result.matrix[i][j]
      s = s + " \{hmisc_fmt_dbl(val, 3)}"
      j = j + 1
    }
    s = s + "\n"
    i = i + 1
  }
  s
}

///|
/// Format DescriptiveStats as a string.
pub fn hmisc_describe_summary(stats : HmiscDescriptiveStats) -> String {
  "Variable: \{stats.name}\n  n: \{stats.n}, missing: \{stats.n_missing}\n  Mean: \{hmisc_fmt_dbl(stats.mean, 4)}\n  SD: \{hmisc_fmt_dbl(stats.sd, 4)}\n  Min: \{hmisc_fmt_dbl(stats.min, 4)}\n  Q1: \{hmisc_fmt_dbl(stats.q1, 4)}\n  Median: \{hmisc_fmt_dbl(stats.median, 4)}\n  Q3: \{hmisc_fmt_dbl(stats.q3, 4)}\n  Max: \{hmisc_fmt_dbl(stats.max, 4)}"
}

///|
/// Format VarClusResult as a string.
pub fn hmisc_varclus_summary(result : HmiscVarClusResult) -> String {
  let mut s = "Variable Clustering Result\n"
  s = s + "  Number of clusters: \{result.n_clusters}\n"
  s = s + "\n"
  let mut i = 0
  while i < result.names.length() {
    s = s + "  \{result.names[i]} -> Cluster \{result.cluster_assignments[i]}\n"
    i = i + 1
  }
  s
}

///|
/// Format double with decimals.
fn hmisc_fmt_dbl(value : Double, decimals : Int) -> String {
  if value.is_nan() {
    return "     NA"
  }
  let rounded = (value * @math.pow(10.0, decimals.to_double())).round() / @math.pow(10.0, decimals.to_double())
  let int_part = rounded.to_int()
  let dec_part = ((rounded - int_part.to_double()) * @math.pow(10.0, decimals.to_double())).round().to_int()
  if decimals == 0 {
    return "\{int_part}"
  }
  let sign = if value < 0.0 { "-" } else { "" }
  let dec_str = hmisc_pad_zeros(dec_part, decimals)
  "\{sign}\{int_part}.\{dec_str}"
}

///|
/// Pad integer with leading zeros.
fn hmisc_pad_zeros(value : Int, length : Int) -> String {
  let abs_val = value.abs()
  let s = "\{abs_val}"
  if s.length() >= length {
    return s
  }
  let mut result = ""
  let mut i = 0
  while i < length - s.length() {
    result = result + "0"
    i = i + 1
  }
  result + s
}

/// Create a Pearson correlation type.
pub fn hmisc_cor_type_pearson() -> HmiscCorType { HmiscCorType::HC_Pearson }

/// Create a Spearman correlation type.
pub fn hmisc_cor_type_spearman() -> HmiscCorType { HmiscCorType::HC_Spearman }

///|
/// Sample data for testing.
pub fn hmisc_sample_data() -> Array[Array[Double]] {
  [
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],
    [2.1, 3.9, 6.1, 8.0, 9.9, 12.0, 14.1, 16.0, 18.1, 20.0],
    [5.0, 4.0, 6.0, 3.0, 7.0, 5.0, 8.0, 4.0, 9.0, 6.0],
  ]
}

///|
/// Sample variable names.
pub fn hmisc_sample_names() -> Array[String] {
  ["X", "Y", "Z"]
}
