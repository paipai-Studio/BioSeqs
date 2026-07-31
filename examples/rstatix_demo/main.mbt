///|
/// rstatix Module Demo - Tidy Data-Friendly Statistical Tests
///
/// This demo showcases the rstatix module functionality including
/// t-tests, Wilcoxon tests, ANOVA, Kruskal-Wallis, correlation tests,
/// and multiple testing correction.
fn main {
  println("=== rstatix Statistical Tests Demo ===")
  println("")

  println("1. Sample Data:")
  let data = @bio.rstatix_sample_data()
  println("   Groups: \{data.length()}")
  println("   Group 1 (n=\{data[0].length()}): \{data[0]}")
  println("   Group 2 (n=\{data[1].length()}): \{data[1]}")
  println("   Group 3 (n=\{data[2].length()}): \{data[2]}")
  println("")

  println("2. Two-Sample T-Test (Group 1 vs Group 2):")
  let t_test_result = @bio.rstatix_t_test(data[0], y=data[1])
  println("   \{@bio.rstatix_t_test_summary(t_test_result)}")
  println("")

  println("3. Two-Sample Wilcoxon Test (Group 1 vs Group 2):")
  let wilcox_result = @bio.rstatix_wilcox_test(data[0], y=data[1])
  println("   \{@bio.rstatix_t_test_summary(wilcox_result)}")
  println("")

  println("4. One-Way ANOVA (all 3 groups):")
  let anova_result = @bio.rstatix_anova_test(data)
  println("   \{@bio.rstatix_anova_summary(anova_result)}")
  println("")

  println("5. Kruskal-Wallis Test (all 3 groups):")
  let kruskal_result = @bio.rstatix_kruskal_test(data)
  println("   \{@bio.rstatix_anova_summary(kruskal_result)}")
  println("")

  println("6. Pearson Correlation Test (Group 1 vs Group 2):")
  let pearson_result = @bio.rstatix_cor_test(data[0], data[1], cor_method="pearson")
  println("   \{@bio.rstatix_cor_test_summary(pearson_result)}")
  println("")

  println("7. Spearman Correlation Test (Group 1 vs Group 2):")
  let spearman_result = @bio.rstatix_cor_test(data[0], data[1], cor_method="spearman")
  println("   \{@bio.rstatix_cor_test_summary(spearman_result)}")
  println("")

  println("8. Benjamini-Hochberg (BH) Correction:")
  let p_values = [t_test_result.p_value, wilcox_result.p_value, anova_result.p_value, kruskal_result.p_value, pearson_result.p_value, spearman_result.p_value]
  let labels = ["t-test", "Wilcoxon", "ANOVA", "Kruskal-Wallis", "Pearson", "Spearman"]
  let bh_corrected = @bio.rstatix_bh_correct(p_values)
  println("   Raw p-values vs BH-adjusted p-values:")
  let mut i = 0
  while i < p_values.length() {
    println("     \{labels[i]}: raw = \{p_values[i]}, adj = \{bh_corrected[i]}")
    i = i + 1
  }
  println("")

  println("9. Bonferroni Correction (for comparison):")
  let bonf_corrected = @bio.rstatix_bonferroni(p_values)
  i = 0
  while i < p_values.length() {
    println("     \{labels[i]}: raw = \{p_values[i]}, adj = \{bonf_corrected[i]}")
    i = i + 1
  }
  println("")

  println("10. Paired T-Test (Group 1 vs Group 2):")
  let paired_t_result = @bio.rstatix_t_test(data[0], y=data[1], paired=true)
  println("   \{@bio.rstatix_t_test_summary(paired_t_result)}")
  println("")

  println("11. One-Sample T-Test (Group 1, mu = 20):")
  let one_sample_t = @bio.rstatix_t_test(data[0], mu=20.0)
  println("   \{@bio.rstatix_t_test_summary(one_sample_t)}")
  println("")

  println("12. One-Sample Wilcoxon Test (Group 1, mu = 18):")
  let one_sample_wilcox = @bio.rstatix_wilcox_test(data[0], mu=18.0)
  println("   \{@bio.rstatix_t_test_summary(one_sample_wilcox)}")
  println("")

  println("=== Demo Complete ===")
}