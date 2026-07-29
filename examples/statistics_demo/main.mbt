// Demo for Bio.Statistics module

fn main {
  println("=== Bio.Statistics Module Demo ===\n")
  
  // Descriptive statistics
  let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
  println("Data: \{data}")
  println("Mean: \{@bio.stat_mean(data)}")
  println("Variance: \{@bio.stat_variance(data)}")
  println("Std Dev: \{@bio.stat_std(data)}")
  println("Median: \{@bio.stat_median(data)}")
  println("Min: \{@bio.stat_min(data)}")
  println("Max: \{@bio.stat_max(data)}")
  println("Sum: \{@bio.stat_sum(data)}")
  
  // Correlation analysis
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y_pos = [2.0, 4.0, 6.0, 8.0, 10.0]
  let y_neg = [10.0, 8.0, 6.0, 4.0, 2.0]
  
  println("\n=== Correlation Analysis ===")
  println("Positive correlation (x vs 2x):")
  println("  Pearson: \{@bio.stat_pearson_correlation(x, y_pos)}")
  println("  Spearman: \{@bio.stat_spearman_correlation(x, y_pos)}")
  
  println("Negative correlation (x vs 20-2x):")
  println("  Pearson: \{@bio.stat_pearson_correlation(x, y_neg)}")
  println("  Spearman: \{@bio.stat_spearman_correlation(x, y_neg)}")
  
  // Hypothesis testing
  let sample = [2.0, 3.0, 4.0, 5.0, 6.0]
  println("\n=== Hypothesis Testing ===")
  println("Sample: \{sample}")
  println("One-sample t-test (mu0=3.0):")
  let t_stat = @bio.stat_t_statistic(sample, 3.0)
  println("  t-statistic: \{t_stat}")
  
  // Confidence interval
  println("95% Confidence Interval:")
  let ci = @bio.stat_confidence_interval(sample)
  println("  [\{ci[0]}, \{ci[1]}]")
  
  // Z-scores
  println("\n=== Z-Scores ===")
  let zscores = @bio.stat_zscores([1.0, 2.0, 3.0, 4.0, 5.0])
  println("Z-scores of [1,2,3,4,5] (mean=3, std=~1.58):")
  let mut i = 0
  while i < zscores.length() {
    println("  z[\{i}] = \{zscores[i]}")
    i = i + 1
  }

  // Non-parametric tests
  println("\n=== Non-parametric Tests ===")
  let x1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y1 = [2.0, 4.0, 6.0, 8.0, 10.0]
  println("Mann-Whitney U test (x vs y):")
  let mw_result = @bio.stat_mann_whitney(x1, y1)
  println("  Statistic: \{mw_result.0}, p-value: \{mw_result.1}")

  println("Wilcoxon signed-rank test (x vs y):")
  let wx_result = @bio.stat_wilcoxon_signed_rank(x1, y1)
  println("  Statistic: \{wx_result.0}, p-value: \{wx_result.1}")

  // KS test
  println("\n=== Kolmogorov-Smirnov Test ===")
  let sample1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let sample2 = [2.0, 4.0, 6.0, 8.0, 10.0]
  println("KS test comparing two samples:")
  let ks_result = @bio.stat_ks_test(sample1, sample2)
  println("  Statistic: \{ks_result.0}, p-value: \{ks_result.1}")

  // Fisher exact test
  println("\n=== Fisher Exact Test ===")
  let table = [[10, 2], [3, 5]]
  println("Contingency table: \{table}")
  let fisher_result = @bio.stat_fisher_exact(table)
  println("  Odds ratio: \{fisher_result.0}, p-value: \{fisher_result.1}")

  // Multiple test correction
  println("\n=== Multiple Test Correction ===")
  let p_values = [0.01, 0.04, 0.03, 0.005, 0.2]
  println("Original p-values: \{p_values}")
  let bonf = @bio.stat_bonferroni(p_values)
  println("Bonferroni corrected: \{bonf}")
  let holm = @bio.stat_holm(p_values)
  println("Holm corrected: \{holm}")
  let by = @bio.stat_by(p_values)
  println("Benjamini-Yekutieli corrected: \{by}")

  // Distribution functions
  println("\n=== Distribution Functions ===")
  println("Chi-square CDF(3.84, df=1): \{@bio.stat_chi2_cdf(3.84, 1)}")
  println("Chi-square quantile(0.95, df=2): \{@bio.stat_chi2_quantile(0.95, 2)}")
  println("t-CDF(0, df=10): \{@bio.stat_t_cdf(0.0, 10.0)}")
  println("t-quantile(0.975, df=10): \{@bio.stat_t_quantile(0.975, 10.0)}")
  println("Normal quantile(0.975): \{@bio.stat_normal_quantile(0.975)}")

  // Hypothesis testing (updated)
  println("\n=== Hypothesis Testing (Updated) ===")
  println("One-sample t-test (mu0=3.0):")
  let t1_result = @bio.stat_t_test_one_sample(sample, 3.0)
  println("  t-statistic: \{t1_result.0}, p-value: \{t1_result.1}")

  println("Two-sample t-test (x vs y):")
  let t2_result = @bio.stat_t_test_two_sample([2.1, 2.3, 2.2, 2.4, 2.0], [3.1, 3.3, 3.2, 3.4, 3.0])
  println("  t-statistic: \{t2_result.0}, p-value: \{t2_result.1}")

  // BH correction
  println("\n=== Benjamini-Hochberg Correction ===")
  let bh = @bio.stat_bh(p_values)
  println("BH corrected p-values: \{bh}")

  // ANOVA
  println("\n=== One-way ANOVA ===")
  let anova_samples = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [8.0, 9.0, 10.0],
  ]
  let anova_result = @bio.stat_anova(anova_samples)
  println("  F-statistic: \{anova_result.0}, p-value: \{anova_result.1}")

  // MAD
  println("\n=== Robust Statistics ===")
  let mad_data = [1.0, 2.0, 3.0, 4.0, 5.0]
  println("MAD of [1,2,3,4,5]: \{@bio.stat_mad(mad_data)}")

  // Log-rank test
  println("\n=== Survival Analysis ===")
  let lr_time = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let lr_event = [true, true, true, false, true, false]
  let lr_group = [0, 0, 0, 1, 1, 1]
  let lr_result = @bio.stat_logrank_test(lr_time, lr_event, lr_group)
  println("Log-rank test: chi2=\{lr_result.0}, p-value=\{lr_result.1}")

  println("\nDemo completed successfully!")
}
