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
  
  println("\nDemo completed successfully!")
}
