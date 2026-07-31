///|
/// rstatix - Tidy Data-Friendly Statistical Tests for BioSeqs.
/// Inspired by Bioconductor's rstatix package.
/// Provides pipe-friendly statistical tests that output tidy data frames
/// compatible with the tidyverse ecosystem.

///|
/// Result of a tidy statistical test.
pub struct RstatixTTestResult {
  test_name : String
  group1 : String
  group2 : String
  statistic : Double
  df : Double
  p_value : Double
  p_adj : Double
  ci_low : Double
  ci_high : Double
  estimate : Double
  se : Double
  n1 : Int
  n2 : Int
  method_name : String
  alternative : String
}

///|
/// Result of a correlation test.
pub struct RstatixCorTestResult {
  cor_type : String
  statistic : Double
  p_value : Double
  p_adj : Double
  ci_low : Double
  ci_high : Double
  n : Int
  method_name : String
}

///|
/// Result of ANOVA test.
pub struct RstatixAnovaResult {
  source : String
  df : Int
  ss : Double
  ms : Double
  f : Double
  p_value : Double
  p_adj : Double
}

///|
/// Perform a one-sample or two-sample t-test (Welch's).
pub fn rstatix_t_test(
  x : Array[Double],
  y? : Array[Double] = [],
  mu? : Double = 0.0,
  alternative? : String = "two.sided",
  paired? : Bool = false,
  conf_level? : Double = 0.95,
) -> RstatixTTestResult {
  let n1 = x.length()
  let is_two_sample = y.length() > 0
  if !is_two_sample {
    rstatix_t_test_one_sample(x, mu, alternative, conf_level, "Sample", "", "One-sample t-test")
  } else if paired {
    let n = rstatix_min_int(n1, y.length())
    let diffs = Array::make(n, 0.0)
    let mut i = 0
    while i < n {
      diffs[i] = x[i] - y[i]
      i = i + 1
    }
    rstatix_t_test_one_sample(diffs, mu, alternative, conf_level, "Difference", "", "Paired t-test")
  } else {
    let n2 = y.length()
    let mean1 = rstatix_mean_double(x)
    let mean2 = rstatix_mean_double(y)
    let sd1 = rstatix_sd_double(x)
    let sd2 = rstatix_sd_double(y)
    let se = (sd1 * sd1 / n1.to_double() + sd2 * sd2 / n2.to_double()).sqrt()
    let t_stat = (mean1 - mean2) / se
    let df_num = (sd1 * sd1 / n1.to_double() + sd2 * sd2 / n2.to_double())
    let df_den = @math.pow(sd1 * sd1 / n1.to_double(), 2.0) / (n1.to_double() - 1.0) + @math.pow(sd2 * sd2 / n2.to_double(), 2.0) / (n2.to_double() - 1.0)
    let df = df_num * df_num / df_den
    let p_val = rstatix_compute_t_p(t_stat, df.to_int(), alternative)
    let tcrit = rstatix_t_critical(1.0 - conf_level, df.to_int(), alternative)
    let diff = mean1 - mean2
    RstatixTTestResult::{
      test_name: "t_test",
      group1: "Group1",
      group2: "Group2",
      statistic: t_stat,
      df,
      p_value: p_val,
      p_adj: p_val,
      ci_low: diff - tcrit * se,
      ci_high: diff + tcrit * se,
      estimate: diff,
      se,
      n1,
      n2,
      method_name: "Welch's two-sample t-test",
      alternative,
    }
  }
}

///|
/// One-sample t-test helper.
fn rstatix_t_test_one_sample(
  x : Array[Double],
  mu : Double,
  alternative : String,
  conf_level : Double,
  g1 : String,
  g2 : String,
  method_name : String,
) -> RstatixTTestResult {
  let n = x.length()
  let mean_x = rstatix_mean_double(x)
  let sd_x = rstatix_sd_double(x)
  let se = sd_x / (n.to_double()).sqrt()
  let t_stat = (mean_x - mu) / se
  let df = n - 1
  let p_val = rstatix_compute_t_p(t_stat, df, alternative)
  let tcrit = rstatix_t_critical(1.0 - conf_level, df, alternative)
  let ci_half = tcrit * se
  RstatixTTestResult::{
    test_name: "t_test",
    group1: g1,
    group2: g2,
    statistic: t_stat,
    df: df.to_double(),
    p_value: p_val,
    p_adj: p_val,
    ci_low: mean_x - ci_half,
    ci_high: mean_x + ci_half,
    estimate: mean_x - mu,
    se,
    n1: n,
    n2: 0,
    method_name,
    alternative,
  }
}

///|
/// Perform a Wilcoxon rank-sum test.
pub fn rstatix_wilcox_test(
  x : Array[Double],
  y? : Array[Double] = [],
  alternative? : String = "two.sided",
  mu? : Double = 0.0,
  paired? : Bool = false,
) -> RstatixTTestResult {
  let is_two_sample = y.length() > 0
  if !is_two_sample {
    let adjusted = rstatix_sub_scalar(x, mu)
    rstatix_wilcox_signed(adjusted, alternative, "One-sample Wilcoxon")
  } else if paired {
    let n = rstatix_min_int(x.length(), y.length())
    let diffs = Array::make(n, 0.0)
    let mut i = 0
    while i < n {
      diffs[i] = x[i] - y[i] - mu
      i = i + 1
    }
    rstatix_wilcox_signed(diffs, alternative, "Paired Wilcoxon")
  } else {
    rstatix_wilcox_rank_sum(x, y, alternative)
  }
}

///|
/// Wilcoxon signed-rank test.
fn rstatix_wilcox_signed(x : Array[Double], alternative : String, method_name : String) -> RstatixTTestResult {
  let non_zero = Array::new()
  let mut i = 0
  while i < x.length() {
    if x[i].abs() > 1.0e-15 {
      non_zero.push(x[i])
    }
    i = i + 1
  }
  let n = non_zero.length()
  if n < 6 {
    return RstatixTTestResult::{
      test_name: "wilcox_test",
      group1: "",
      group2: "",
      statistic: @double.not_a_number,
      df: @double.not_a_number,
      p_value: @double.not_a_number,
      p_adj: @double.not_a_number,
      ci_low: @double.not_a_number,
      ci_high: @double.not_a_number,
      estimate: @double.not_a_number,
      se: @double.not_a_number,
      n1: n,
      n2: 0,
      method_name: method_name + " (too few observations)",
      alternative,
    }
  }
  let abs_vals = Array::make(n, 0.0)
  i = 0
  while i < n {
    abs_vals[i] = non_zero[i].abs()
    i = i + 1
  }
  let ranks = rstatix_rank_array(abs_vals)
  let mut t_plus = 0.0
  let mut t_minus = 0.0
  i = 0
  while i < n {
    if non_zero[i] > 0.0 {
      t_plus = t_plus + ranks[i]
    } else {
      t_minus = t_minus + ranks[i]
    }
    i = i + 1
  }
  let w = rstatix_min_double(t_plus, t_minus)
  let mu_w = (n.to_double() * (n.to_double() + 1.0)) / 4.0
  let sigma_w = ((n.to_double() * (n.to_double() + 1.0) * (2.0 * n.to_double() + 1.0)) / 24.0).sqrt()
  let z = (w - mu_w) / sigma_w
  let p_val = rstatix_norm_p(z, alternative)
  let med_est = rstatix_median_double(non_zero)
  RstatixTTestResult::{
    test_name: "wilcox_test",
    group1: "Sample",
    group2: "",
    statistic: w,
    df: n.to_double(),
    p_value: p_val,
    p_adj: p_val,
    ci_low: med_est - 1.96 * sigma_w,
    ci_high: med_est + 1.96 * sigma_w,
    estimate: med_est,
    se: sigma_w,
    n1: n,
    n2: 0,
    method_name,
    alternative,
  }
}

///|
/// Wilcoxon rank-sum test (two independent samples).
fn rstatix_wilcox_rank_sum(
  x : Array[Double],
  y : Array[Double],
  alternative : String,
) -> RstatixTTestResult {
  let n1 = x.length()
  let n2 = y.length()
  let n = n1 + n2
  if n1 < 3 || n2 < 3 {
    return RstatixTTestResult::{
      test_name: "wilcox_test",
      group1: "Group1",
      group2: "Group2",
      statistic: @double.not_a_number,
      df: @double.not_a_number,
      p_value: @double.not_a_number,
      p_adj: @double.not_a_number,
      ci_low: @double.not_a_number,
      ci_high: @double.not_a_number,
      estimate: @double.not_a_number,
      se: @double.not_a_number,
      n1,
      n2,
      method_name: "Wilcoxon rank-sum (too few observations)",
      alternative,
    }
  }
  let combined = Array::make(n, 0.0)
  let mut i = 0
  while i < n1 {
    combined[i] = x[i]
    i = i + 1
  }
  let mut j = 0
  while j < n2 {
    combined[n1 + j] = y[j]
    j = j + 1
  }
  let ranks = rstatix_rank_array(combined)
  let mut r1 = 0.0
  i = 0
  while i < n1 {
    r1 = r1 + ranks[i]
    i = i + 1
  }
  let u1 = r1 - (n1.to_double() * (n1.to_double() + 1.0)) / 2.0
  let u2 = n1.to_double() * n2.to_double() - u1
  let u = rstatix_min_double(u1, u2)
  let mu_u = (n1.to_double() * n2.to_double()) / 2.0
  let sigma_u = ((n1.to_double() * n2.to_double() * (n.to_double() + 1.0)) / 12.0).sqrt()
  let z = (u - mu_u) / sigma_u
  let p_val = rstatix_norm_p(z, alternative)
  let median1 = rstatix_median_double(x)
  let median2 = rstatix_median_double(y)
  RstatixTTestResult::{
    test_name: "wilcox_test",
    group1: "Group1",
    group2: "Group2",
    statistic: u,
    df: (n1 + n2).to_double(),
    p_value: p_val,
    p_adj: p_val,
    ci_low: median1 - median2 - 1.96 * sigma_u,
    ci_high: median1 - median2 + 1.96 * sigma_u,
    estimate: median1 - median2,
    se: sigma_u,
    n1,
    n2,
    method_name: "Wilcoxon rank-sum test",
    alternative,
  }
}

///|
/// Compute correlation test.
pub fn rstatix_cor_test(
  x : Array[Double],
  y : Array[Double],
  cor_method? : String = "pearson",
  conf_level? : Double = 0.95,
) -> RstatixCorTestResult {
  let n = x.length()
  if n != y.length() || n < 3 {
    return RstatixCorTestResult::{
      cor_type: cor_method,
      statistic: @double.not_a_number,
      p_value: @double.not_a_number,
      p_adj: @double.not_a_number,
      ci_low: @double.not_a_number,
      ci_high: @double.not_a_number,
      n,
      method_name: cor_method,
    }
  }
  let (r, p, mname) = match cor_method {
    "spearman" => {
      let sr = rstatix_spearman_val(x, y)
      let sp = if n > 3 {
        let t_stat = sr * ((n.to_double() - 2.0) / (1.0 - sr * sr)).sqrt()
        2.0 * (1.0 - rstatix_norm_cdf(t_stat.abs()))
      } else { @double.not_a_number }
      (sr, sp, "Spearman's rank correlation")
    }
    "kendall" => {
      let kt = rstatix_kendall_tau(x, y)
      let np = n.to_double()
      let var_kt = 2.0 * (2.0 * np + 5.0) / (9.0 * np * (np - 1.0))
      let z = kt / var_kt.sqrt()
      let kp = 2.0 * (1.0 - rstatix_norm_cdf(z.abs()))
      (kt, kp, "Kendall's tau correlation")
    }
    _ => {
      let pr = rstatix_pearson_val(x, y)
      let pp = if n > 2 {
        let t_stat = pr * ((n.to_double() - 2.0) / (1.0 - pr * pr)).sqrt()
        2.0 * rstatix_t_dist_p(t_stat.abs(), n - 2)
      } else { @double.not_a_number }
      (pr, pp, "Pearson's product-moment correlation")
    }
  }
  let (ci_low, ci_high) = if r.abs() >= 1.0 - 1.0e-15 {
    (r, r)
  } else {
    let z = 0.5 * @math.ln((1.0 + r) / (1.0 - r))
    let se_z = 1.0 / ((n.to_double() - 3.0).sqrt())
    let z_crit = 1.96
    (@math.exp(z - z_crit * se_z), @math.exp(z + z_crit * se_z))
  }
  RstatixCorTestResult::{
    cor_type: cor_method,
    statistic: r,
    p_value: p,
    p_adj: p,
    ci_low,
    ci_high,
    n,
    method_name: mname,
  }
}

///|
/// One-way ANOVA.
pub fn rstatix_anova_test(
  groups : Array[Array[Double]],
  group_names? : Array[String] = [],
) -> RstatixAnovaResult {
  let k = groups.length()
  if k < 3 {
    return RstatixAnovaResult::{
      source: "Between",
      df: k - 1,
      ss: @double.not_a_number,
      ms: @double.not_a_number,
      f: @double.not_a_number,
      p_value: @double.not_a_number,
      p_adj: @double.not_a_number,
    }
  }
  let all_data = Array::new()
  let n_groups = Array::make(k, 0)
  let group_means = Array::make(k, 0.0)
  let mut total_n = 0
  let mut i = 0
  while i < k {
    let g = groups[i]
    let ng = g.length()
    n_groups[i] = ng
    group_means[i] = rstatix_mean_double(g)
    let mut j = 0
    while j < ng {
      all_data.push(g[j])
      j = j + 1
    }
    total_n = total_n + ng
    i = i + 1
  }
  let grand_mean = rstatix_mean_double(all_data)
  let mut ss_between = 0.0
  i = 0
  while i < k {
    ss_between = ss_between + (n_groups[i].to_double() * (group_means[i] - grand_mean) * (group_means[i] - grand_mean))
    i = i + 1
  }
  let mut ss_within = 0.0
  i = 0
  while i < k {
    let g = groups[i]
    let ng = g.length()
    let mut j = 0
    while j < ng {
      ss_within = ss_within + (g[j] - group_means[i]) * (g[j] - group_means[i])
      j = j + 1
    }
    i = i + 1
  }
  let df_between = k - 1
  let df_within = total_n - k
  let ms_between = ss_between / df_between.to_double()
  let ms_within = if df_within > 0 { ss_within / df_within.to_double() } else { 0.0 }
  let f_stat = if ms_within < 1.0e-15 {
    if ss_between < 1.0e-15 { 0.0 } else { @double.infinity }
  } else {
    ms_between / ms_within
  }
  let p_val = if f_stat == @double.infinity { 0.0 } else if f_stat < 1.0e-15 { 1.0 } else { rstatix_f_p(f_stat, df_between, df_within) }
  RstatixAnovaResult::{
    source: "Between",
    df: df_between,
    ss: ss_between,
    ms: ms_between,
    f: f_stat,
    p_value: p_val,
    p_adj: p_val,
  }
}

///|
/// Kruskal-Wallis rank sum test.
pub fn rstatix_kruskal_test(
  groups : Array[Array[Double]],
  _group_names? : Array[String] = [],
) -> RstatixAnovaResult {
  let k = groups.length()
  if k < 3 {
    return RstatixAnovaResult::{
      source: "Kruskal-Wallis",
      df: k - 1,
      ss: @double.not_a_number,
      ms: @double.not_a_number,
      f: @double.not_a_number,
      p_value: @double.not_a_number,
      p_adj: @double.not_a_number,
    }
  }
  let all_data = Array::new()
  let group_ids = Array::new()
  let mut n_total = 0
  let mut i = 0
  while i < k {
    let g = groups[i]
    let mut j = 0
    while j < g.length() {
      all_data.push(g[j])
      group_ids.push(i)
      j = j + 1
    }
    n_total = n_total + g.length()
    i = i + 1
  }
  let ranks = rstatix_rank_array(all_data)
  let mut sum_r_sq = 0.0
  i = 0
  while i < k {
    let n_i = groups[i].length()
    let mut sum_r = 0.0
    let mut j = 0
    while j < n_total {
      if group_ids[j] == i {
        sum_r = sum_r + ranks[j]
      }
      j = j + 1
    }
    sum_r_sq = sum_r_sq + (sum_r * sum_r) / n_i.to_double()
    i = i + 1
  }
  let h = (12.0 / (n_total.to_double() * (n_total.to_double() + 1.0))) * sum_r_sq - 3.0 * (n_total.to_double() + 1.0)
  let df = k - 1
  let p_val = rstatix_chi2_p(h, df)
  RstatixAnovaResult::{
    source: "Kruskal-Wallis",
    df,
    ss: h,
    ms: h / df.to_double(),
    f: h,
    p_value: p_val,
    p_adj: p_val,
  }
}

///|
/// Friedman rank sum test.
pub fn rstatix_friedman_test(
  groups : Array[Array[Double]],
  _group_names? : Array[String] = [],
) -> RstatixAnovaResult {
  let k = groups.length()
  if k < 3 {
    return RstatixAnovaResult::{
      source: "Friedman",
      df: k - 1,
      ss: @double.not_a_number,
      ms: @double.not_a_number,
      f: @double.not_a_number,
      p_value: @double.not_a_number,
      p_adj: @double.not_a_number,
    }
  }
  let n_blocks = groups[0].length()
  let rank_sums = Array::make(k, 0.0)
  let mut i = 0
  while i < n_blocks {
    let block_vals = Array::make(k, 0.0)
    let mut j = 0
    while j < k {
      block_vals[j] = groups[j][i]
      j = j + 1
    }
    let block_ranks = rstatix_rank_array(block_vals)
    let mut r_i = 0
    while r_i < k {
      rank_sums[r_i] = rank_sums[r_i] + block_ranks[r_i]
      r_i = r_i + 1
    }
    i = i + 1
  }
  let mut total = 0.0
  let mut j = 0
  while j < k {
    total = total + rank_sums[j] * rank_sums[j]
    j = j + 1
  }
  let q = (12.0 / (n_blocks.to_double() * k.to_double() * (k.to_double() + 1.0))) * total - 3.0 * n_blocks.to_double() * (k.to_double() + 1.0)
  let df = k - 1
  let p_val = rstatix_chi2_p(q, df)
  RstatixAnovaResult::{
    source: "Friedman",
    df,
    ss: q,
    ms: q / df.to_double(),
    f: q,
    p_value: p_val,
    p_adj: p_val,
  }
}

///|
/// Apply Benjamini-Hochberg FDR correction.
pub fn rstatix_bh_correct(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 { return Array::new() }
  let result = Array::make(n, 0.0)
  let indices = Array::make(n, 0)
  let mut i = 0
  while i < n {
    indices[i] = i
    i = i + 1
  }
  indices.sort_by(fn(a : Int, b : Int) -> Int { if p_values[a] < p_values[b] { -1 } else if p_values[a] > p_values[b] { 1 } else { 0 } })
  let mut min_so_far = 1.0
  i = n - 1
  while i >= 0 {
    let rank = (i + 1).to_double()
    let idx = indices[i]
    let adj = p_values[idx] * (n.to_double() / rank)
    min_so_far = rstatix_min_double(min_so_far, adj)
    result[idx] = rstatix_min_double(1.0, min_so_far)
    i = i - 1
  }
  result
}

///|
/// Apply Bonferroni correction.
pub fn rstatix_bonferroni(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  let result = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    result[i] = rstatix_min_double(1.0, p_values[i] * n.to_double())
    i = i + 1
  }
  result
}

///|
/// Format t-test result.
pub fn rstatix_t_test_summary(r : RstatixTTestResult) -> String {
  "T-test: \{r.test_name}\n  Groups: \{r.group1} vs \{r.group2}\n  Statistic: \{rstatix_fmt(r.statistic)}\n  df: \{rstatix_fmt(r.df)}\n  p-value: \{rstatix_fmt(r.p_value)}\n  95% CI: [\{rstatix_fmt(r.ci_low)}, \{rstatix_fmt(r.ci_high)}]\n  Estimate: \{rstatix_fmt(r.estimate)}\n  Method: \{r.method_name}\n  Alternative: \{r.alternative}"
}

///|
/// Format correlation test result.
pub fn rstatix_cor_test_summary(r : RstatixCorTestResult) -> String {
  "Correlation Test: \{r.cor_type}\n  Statistic: \{rstatix_fmt(r.statistic)}\n  p-value: \{rstatix_fmt(r.p_value)}\n  95% CI: [\{rstatix_fmt(r.ci_low)}, \{rstatix_fmt(r.ci_high)}]\n  n: \{r.n}\n  Method: \{r.method_name}"
}

///|
/// Format ANOVA result.
pub fn rstatix_anova_summary(r : RstatixAnovaResult) -> String {
  "ANOVA Result\n  Source: \{r.source}\n  df: \{r.df}\n  SS: \{rstatix_fmt(r.ss)}\n  MS: \{rstatix_fmt(r.ms)}\n  F: \{rstatix_fmt(r.f)}\n  p-value: \{rstatix_fmt(r.p_value)}"
}

///|
/// Sample data for testing.
pub fn rstatix_sample_data() -> Array[Array[Double]] {
  [
    [10.0, 12.0, 15.0, 18.0, 20.0, 22.0, 25.0, 28.0, 30.0, 32.0],
    [8.0, 9.0, 11.0, 12.0, 14.0, 15.0, 17.0, 19.0, 20.0, 22.0],
    [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0],
  ]
}

// ============ Private Helpers ============

fn rstatix_min_int(a : Int, b : Int) -> Int { if a < b { a } else { b } }
fn rstatix_min_double(a : Double, b : Double) -> Double { if a < b { a } else { b } }

fn rstatix_mean_double(data : Array[Double]) -> Double {
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

fn rstatix_sd_double(data : Array[Double]) -> Double {
  let n = data.length()
  if n <= 1 { return 0.0 }
  let mean = rstatix_mean_double(data)
  let mut sum_sq = 0.0
  let mut count = 0
  let mut i = 0
  while i < n {
    if !data[i].is_nan() {
      let diff = data[i] - mean
      sum_sq = sum_sq + diff * diff
      count = count + 1
    }
    i = i + 1
  }
  if count > 1 { (sum_sq / (count.to_double() - 1.0)).sqrt() } else { 0.0 }
}

fn rstatix_median_double(data : Array[Double]) -> Double {
  if data.length() == 0 { return 0.0 }
  let sorted = data.copy()
  sorted.sort()
  let n = sorted.length()
  if n % 2 == 0 { (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0 } else { sorted[n / 2] }
}

fn rstatix_sub_scalar(data : Array[Double], scalar : Double) -> Array[Double] {
  let result = Array::make(data.length(), 0.0)
  let mut i = 0
  while i < data.length() {
    result[i] = data[i] - scalar
    i = i + 1
  }
  result
}

fn rstatix_pearson_val(x : Array[Double], y : Array[Double]) -> Double {
  let n = x.length()
  let mean_x = rstatix_mean_double(x)
  let mean_y = rstatix_mean_double(y)
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
  if denom < 1.0e-15 { @double.not_a_number } else { sum_xy / denom }
}

fn rstatix_spearman_val(x : Array[Double], y : Array[Double]) -> Double {
  let rx = rstatix_rank_array(x)
  let ry = rstatix_rank_array(y)
  rstatix_pearson_val(rx, ry)
}

fn rstatix_rank_array(data : Array[Double]) -> Array[Double] {
  let n = data.length()
  let indices = Array::make(n, 0)
  let mut i = 0
  while i < n { indices[i] = i; i = i + 1 }
  indices.sort_by(fn(a : Int, b : Int) -> Int { if data[a] < data[b] { -1 } else if data[a] > data[b] { 1 } else { 0 } })
  let ranks = Array::make(n, 0.0)
  i = 0
  while i < n {
    let mut j = i
    while j < n - 1 && data[indices[j]] == data[indices[j + 1]] { j = j + 1 }
    let avg = ((i + j).to_double() / 2.0) + 1.0
    let mut k = i
    while k <= j { ranks[indices[k]] = avg; k = k + 1 }
    i = j + 1
  }
  ranks
}

fn rstatix_kendall_tau(x : Array[Double], y : Array[Double]) -> Double {
  let n = x.length()
  let mut conc = 0
  let mut disc = 0
  let mut i = 0
  while i < n - 1 {
    let mut j = i + 1
    while j < n {
      let dx = x[j] - x[i]
      let dy = y[j] - y[i]
      if dx * dy > 1.0e-15 { conc = conc + 1 } else if dx * dy < -1.0e-15 { disc = disc + 1 }
      j = j + 1
    }
    i = i + 1
  }
  let total = (n * (n - 1)) / 2
  if total > 0 { (conc.to_double() - disc.to_double()) / total.to_double() } else { 0.0 }
}

fn rstatix_compute_t_p(t : Double, df : Int, alt : String) -> Double {
  let p_two = 2.0 * rstatix_t_dist_p(t.abs(), df)
  match alt {
    "greater" => if t > 0.0 { p_two / 2.0 } else { 1.0 - p_two / 2.0 }
    "less" => if t < 0.0 { p_two / 2.0 } else { 1.0 - p_two / 2.0 }
    _ => p_two
  }
}

fn rstatix_t_dist_p(t : Double, df : Int) -> Double {
  if df <= 0 || t.is_nan() { return @double.not_a_number }
  let x = df.to_double() / (df.to_double() + t * t)
  rstatix_ibeta(df.to_double() / 2.0, 0.5, x)
}

fn rstatix_ibeta(a : Double, b : Double, x : Double) -> Double {
  if x < 0.0 || x > 1.0 { return 0.0 }
  if x == 0.0 || x == 1.0 { return x }
  let g = 7.0
  let c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012,
    9.9843695780195716e-6, 1.5056327351493116e-7]
  let pi = 3.14159265358979323846
  let lbeta = rstatix_lgamma(a) + rstatix_lgamma(b) - rstatix_lgamma(a + b)
  let front = @math.exp(a * @math.ln(x) + b * @math.ln(1.0 - x) - lbeta) / a
  let max_iter = 100
  let eps = 1.0e-10
  let mut f = 1.0
  let mut c_val = 1.0
  let mut d = 1.0 - (a + b) * x / (a + 1.0)
  if d.abs() < 1.0e-30 { d = 1.0e-30 }
  d = 1.0 / d
  f = d
  let mut m = 1
  while m <= max_iter {
    let n1 = m.to_double() * (b - m.to_double()) * x / ((a + 2.0 * m.to_double() - 1.0) * (a + 2.0 * m.to_double()))
    d = 1.0 + n1 * d
    if d.abs() < 1.0e-30 { d = 1.0e-30 }
    c_val = 1.0 + n1 / c_val
    if c_val.abs() < 1.0e-30 { c_val = 1.0e-30 }
    d = 1.0 / d
    f = f * d * c_val
    let n2 = -(a + m.to_double()) * (a + b + m.to_double()) * x / ((a + 2.0 * m.to_double()) * (a + 2.0 * m.to_double() + 1.0))
    d = 1.0 + n2 * d
    if d.abs() < 1.0e-30 { d = 1.0e-30 }
    c_val = 1.0 + n2 / c_val
    if c_val.abs() < 1.0e-30 { c_val = 1.0e-30 }
    d = 1.0 / d
    let del = d * c_val
    f = f * del
    if (del - 1.0).abs() < eps { break }
    m = m + 1
  }
  front * f
}

fn rstatix_lgamma(x : Double) -> Double {
  if x <= 0.0 { return @double.not_a_number }
  let g = 7.0
  let c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012,
    9.9843695780195716e-6, 1.5056327351493116e-7]
  let pi = 3.14159265358979323846
  if x < 0.5 { return @math.ln(pi / @math.sin(pi * x)) - rstatix_lgamma(1.0 - x) }
  let mut sum = c[0]
  let mut i = 1
  while i < 9 { sum = sum + c[i] / (x + i.to_double() - 1.0); i = i + 1 }
  let t = x + g + 0.5
  0.5 * @math.ln(2.0 * pi) + (t * @math.ln(t) - t) + @math.ln(sum)
}

fn rstatix_t_critical(alpha : Double, df : Int, alt : String) -> Double {
  let a = if alt == "two.sided" { alpha / 2.0 } else { alpha }
  rstatix_approx_tc(a, df)
}

fn rstatix_approx_tc(alpha : Double, df : Int) -> Double {
  let z = rstatix_norm_q(1.0 - alpha)
  let df_d = df.to_double()
  z + (z * z * z + z) / (4.0 * df_d) + (5.0 * z + 16.0 * z * z * z + 3.0 * @math.pow(z, 5.0)) / (96.0 * df_d * df_d)
}

fn rstatix_norm_cdf(x : Double) -> Double {
  0.5 * (1.0 + rstatix_erf(x / 1.4142135623730951))
}

fn rstatix_erf(x : Double) -> Double {
  let t = 1.0 / (1.0 + 0.3275911 * x.abs())
  let a1 = 0.254829592
  let a2 = -0.284496736
  let a3 = 1.421413741
  let a4 = -1.453152027
  let a5 = 1.061405429
  let res = 1.0 - (a1 * t + a2 * t * t + a3 * t * t * t + a4 * t * t * t * t + a5 * t * t * t * t * t) * @math.exp(-x * x)
  if x >= 0.0 { res } else { -res }
}

fn rstatix_norm_q(p : Double) -> Double {
  if p <= 0.0 || p >= 1.0 { return @double.not_a_number }
  if p < 0.5 {
    let eta = (-2.0 * @math.ln(p)).sqrt()
    -(eta - (2.515517 + 0.802853 * eta + 0.010328 * eta * eta) / (1.0 + 1.432788 * eta + 0.189269 * eta * eta + 0.001308 * eta * eta * eta))
  } else {
    let eta = (-2.0 * @math.ln(1.0 - p)).sqrt()
    eta - (2.515517 + 0.802853 * eta + 0.010328 * eta * eta) / (1.0 + 1.432788 * eta + 0.189269 * eta * eta + 0.001308 * eta * eta * eta)
  }
}

fn rstatix_norm_p(z : Double, alt : String) -> Double {
  let p_two = 2.0 * (1.0 - rstatix_norm_cdf(z.abs()))
  match alt {
    "greater" => if z > 0.0 { p_two / 2.0 } else { 1.0 - p_two / 2.0 }
    "less" => if z < 0.0 { p_two / 2.0 } else { 1.0 - p_two / 2.0 }
    _ => p_two
  }
}

fn rstatix_f_p(f : Double, df1 : Int, df2 : Int) -> Double {
  if f <= 0.0 { return 1.0 }
  let x = df1.to_double() / (df1.to_double() + df2.to_double() * f)
  1.0 - rstatix_ibeta(df1.to_double() / 2.0, df2.to_double() / 2.0, x)
}

fn rstatix_chi2_p(chi2 : Double, df : Int) -> Double {
  if chi2 <= 0.0 { return 1.0 }
  let x = chi2 / (chi2 + (df.to_double() + 2.0))
  1.0 - rstatix_ibeta(df.to_double() / 2.0, 0.5, x)
}

fn rstatix_fmt(v : Double) -> String {
  if v.is_nan() { return "NA" }
  if v == @double.infinity || v == @double.neg_infinity { return "Inf" }
  let r = (v * 10000.0).round() / 10000.0
  "\{r}"
}
