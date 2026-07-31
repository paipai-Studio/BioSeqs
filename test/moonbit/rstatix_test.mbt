// Tests for rstatix module - tidy data-friendly statistical tests

// ===== rstatix_t_test =====

test "rstatix_t_test one-sample basic" {
  let x = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0]
  let result = @src.rstatix_t_test(x, mu=15.0)
  assert_eq(result.test_name, "t_test")
  assert_eq(result.method_name, "One-sample t-test")
  assert_true(result.statistic > 0.0)
  assert_true(result.df > 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(result.ci_low < result.ci_high)
  assert_true(result.n1 == 10)
  assert_eq(result.alternative, "two.sided")
}

test "rstatix_t_test one-sample greater alternative" {
  let x = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0]
  let result = @src.rstatix_t_test(x, mu=15.0, alternative="greater")
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_eq(result.alternative, "greater")
}

test "rstatix_t_test one-sample less alternative" {
  let x = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0]
  let result = @src.rstatix_t_test(x, mu=25.0, alternative="less")
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_eq(result.alternative, "less")
}

test "rstatix_t_test one-sample small array" {
  let x = [5.0, 6.0]
  let result = @src.rstatix_t_test(x, mu=0.0)
  assert_true(result.df == 1.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_t_test two-sample Welch's basic" {
  let x = [10.0, 12.0, 15.0, 18.0, 20.0, 22.0, 25.0, 28.0, 30.0, 32.0]
  let y = [8.0, 9.0, 11.0, 12.0, 14.0, 15.0, 17.0, 19.0, 20.0, 22.0]
  let result = @src.rstatix_t_test(x, y=y)
  assert_eq(result.test_name, "t_test")
  assert_eq(result.method_name, "Welch's two-sample t-test")
  assert_true(result.statistic > 0.0)
  assert_true(result.df > 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(result.ci_low < result.ci_high)
  assert_true(result.n1 == 10)
  assert_true(result.n2 == 10)
  assert_true(result.se > 0.0)
}

test "rstatix_t_test two-sample equal means" {
  let x = [5.0, 5.0, 5.0, 5.0, 5.0]
  let y = [5.0, 5.0, 5.0, 5.0, 5.0]
  let result = @src.rstatix_t_test(x, y=y)
  assert_true(result.p_value.is_nan())
}

test "rstatix_t_test paired basic" {
  let x = [10.0, 12.0, 15.0, 18.0, 20.0, 22.0, 25.0, 28.0, 30.0, 32.0]
  let y = [8.0, 9.0, 11.0, 12.0, 14.0, 15.0, 17.0, 19.0, 20.0, 22.0]
  let result = @src.rstatix_t_test(x, y=y, paired=true)
  assert_eq(result.method_name, "Paired t-test")
  assert_true(result.statistic > 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(result.ci_low < result.ci_high)
}

test "rstatix_t_test paired different length" {
  let x = [10.0, 12.0, 15.0, 18.0, 20.0, 22.0, 25.0]
  let y = [8.0, 9.0, 11.0, 12.0, 14.0]
  let result = @src.rstatix_t_test(x, y=y, paired=true)
  assert_true(result.n1 == 5)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_t_test confidence interval" {
  let x = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0]
  let result99 = @src.rstatix_t_test(x, mu=15.0, conf_level=0.99)
  let result95 = @src.rstatix_t_test(x, mu=15.0, conf_level=0.95)
  assert_true(result99.ci_low < result95.ci_low)
  assert_true(result99.ci_high > result95.ci_high)
}

// ===== rstatix_wilcox_test =====

test "rstatix_wilcox_test one-sample basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let result = @src.rstatix_wilcox_test(x, mu=5.0)
  assert_eq(result.test_name, "wilcox_test")
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(!result.statistic.is_nan())
  assert_true(result.n1 == 9)
}

test "rstatix_wilcox_test one-sample too few observations" {
  let x = [1.0, 2.0, 3.0]
  let result = @src.rstatix_wilcox_test(x, mu=2.0)
  assert_true(result.p_value.is_nan())
  assert_true(result.statistic.is_nan())
  assert_true(result.n1 == 2)
}

test "rstatix_wilcox_test one-sample greater alternative" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let result = @src.rstatix_wilcox_test(x, mu=4.0, alternative="greater")
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_eq(result.alternative, "greater")
}

test "rstatix_wilcox_test two-sample rank-sum basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let y = [8.0, 9.0, 10.0, 11.0, 12.0, 13.0]
  let result = @src.rstatix_wilcox_test(x, y=y)
  assert_eq(result.test_name, "wilcox_test")
  assert_eq(result.method_name, "Wilcoxon rank-sum test")
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(!result.statistic.is_nan())
  assert_true(result.n1 == 6)
  assert_true(result.n2 == 6)
}

test "rstatix_wilcox_test two-sample too few observations" {
  let x = [1.0, 2.0]
  let y = [3.0, 4.0]
  let result = @src.rstatix_wilcox_test(x, y=y)
  assert_true(result.p_value.is_nan())
  assert_eq(result.n1, 2)
  assert_eq(result.n2, 2)
}

test "rstatix_wilcox_test two-sample one group too small" {
  let x = [1.0, 2.0, 3.0]
  let y = [4.0, 5.0]
  let result = @src.rstatix_wilcox_test(x, y=y)
  assert_true(result.p_value.is_nan())
}

test "rstatix_wilcox_test paired basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
  let result = @src.rstatix_wilcox_test(x, y=y, paired=true)
  assert_eq(result.method_name, "Paired Wilcoxon")
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(!result.statistic.is_nan())
}

test "rstatix_wilcox_test paired with mu" {
  let x = [5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
  let y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let result = @src.rstatix_wilcox_test(x, y=y, paired=true, mu=2.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

// ===== rstatix_cor_test =====

test "rstatix_cor_test pearson basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
  let result = @src.rstatix_cor_test(x, y, cor_method="pearson")
  let n_val = result.n
  let ct = result.cor_type
  let mn = result.method_name
  let stat = result.statistic
  let pv = result.p_value
  let ci_lo = result.ci_low
  let ci_hi = result.ci_high
  assert_eq(ct, "pearson")
  assert_eq(mn, "Pearson's product-moment correlation")
  assert_true(stat > 0.99 && stat <= 1.0)
  assert_true(pv >= 0.0 && pv <= 1.0)
  assert_true(ci_lo <= ci_hi)
  assert_true(n_val == 10)
}

test "rstatix_cor_test pearson negative" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [20.0, 18.0, 16.0, 14.0, 12.0, 10.0, 8.0, 6.0, 4.0, 2.0]
  let result = @src.rstatix_cor_test(x, y)
  assert_true(result.statistic < -0.99 && result.statistic >= -1.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_cor_test pearson no correlation" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
  let result = @src.rstatix_cor_test(x, y)
  assert_true(result.statistic.is_nan())
  assert_true(result.p_value.is_nan())
}

test "rstatix_cor_test spearman basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
  let result = @src.rstatix_cor_test(x, y, cor_method="spearman")
  assert_eq(result.cor_type, "spearman")
  assert_eq(result.method_name, "Spearman's rank correlation")
  assert_true(result.statistic > 0.99 && result.statistic <= 1.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_cor_test spearman nonlinear" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0]
  let result = @src.rstatix_cor_test(x, y, cor_method="spearman")
  assert_true(result.statistic > 0.9)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_cor_test kendall basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
  let result = @src.rstatix_cor_test(x, y, cor_method="kendall")
  assert_eq(result.cor_type, "kendall")
  assert_eq(result.method_name, "Kendall's tau correlation")
  assert_true(result.statistic > 0.9)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_cor_test kendall negative" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [20.0, 18.0, 16.0, 14.0, 12.0, 10.0, 8.0, 6.0, 4.0, 2.0]
  let result = @src.rstatix_cor_test(x, y, cor_method="kendall")
  assert_true(result.statistic < -0.9)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_cor_test mismatched length" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [1.0, 2.0, 3.0, 4.0]
  let result = @src.rstatix_cor_test(x, y)
  assert_true(result.statistic.is_nan())
  assert_true(result.p_value.is_nan())
}

test "rstatix_cor_test too small" {
  let x = [1.0, 2.0]
  let y = [3.0, 4.0]
  let result = @src.rstatix_cor_test(x, y)
  assert_true(result.statistic.is_nan())
  assert_true(result.p_value.is_nan())
}

// ===== rstatix_anova_test =====

test "rstatix_anova_test basic" {
  let groups = [
    [10.0, 12.0, 15.0, 18.0, 20.0],
    [8.0, 9.0, 11.0, 12.0, 14.0],
    [5.0, 6.0, 7.0, 8.0, 9.0],
  ]
  let result = @src.rstatix_anova_test(groups)
  assert_eq(result.source, "Between")
  assert_true(result.df == 2)
  assert_true(result.f > 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(result.ss > 0.0)
  assert_true(result.ms > 0.0)
}

test "rstatix_anova_test equal groups" {
  let groups = [
    [5.0, 5.0, 5.0, 5.0, 5.0],
    [5.0, 5.0, 5.0, 5.0, 5.0],
    [5.0, 5.0, 5.0, 5.0, 5.0],
  ]
  let result = @src.rstatix_anova_test(groups)
  assert_true(result.f >= 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_anova_test too few groups" {
  let groups = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]
  let result = @src.rstatix_anova_test(groups)
  assert_true(result.p_value.is_nan())
  assert_true(result.f.is_nan())
}

test "rstatix_anova_test two groups" {
  let groups = [
    [1.0, 2.0, 3.0],
  ]
  let result = @src.rstatix_anova_test(groups)
  assert_true(result.p_value.is_nan())
}

test "rstatix_anova_test four groups" {
  let groups = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [6.0, 7.0, 8.0, 9.0, 10.0],
    [11.0, 12.0, 13.0, 14.0, 15.0],
    [16.0, 17.0, 18.0, 19.0, 20.0],
  ]
  let result = @src.rstatix_anova_test(groups)
  assert_true(result.df == 3)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

// ===== rstatix_kruskal_test =====

test "rstatix_kruskal_test basic" {
  let groups = [
    [10.0, 12.0, 15.0, 18.0, 20.0],
    [8.0, 9.0, 11.0, 12.0, 14.0],
    [5.0, 6.0, 7.0, 8.0, 9.0],
  ]
  let result = @src.rstatix_kruskal_test(groups)
  assert_eq(result.source, "Kruskal-Wallis")
  assert_true(result.df == 2)
  assert_true(result.f > 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(result.ss > 0.0)
}

test "rstatix_kruskal_test equal groups" {
  let groups = [
    [5.0, 5.0, 5.0, 5.0, 5.0],
    [5.0, 5.0, 5.0, 5.0, 5.0],
    [5.0, 5.0, 5.0, 5.0, 5.0],
  ]
  let result = @src.rstatix_kruskal_test(groups)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_kruskal_test too few groups" {
  let groups = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]
  let result = @src.rstatix_kruskal_test(groups)
  assert_true(result.p_value.is_nan())
  assert_true(result.f.is_nan())
}

// ===== rstatix_friedman_test =====

test "rstatix_friedman_test basic" {
  let groups = [
    [10.0, 12.0, 15.0, 18.0, 20.0],
    [8.0, 9.0, 11.0, 12.0, 14.0],
    [5.0, 6.0, 7.0, 8.0, 9.0],
  ]
  let result = @src.rstatix_friedman_test(groups)
  assert_eq(result.source, "Friedman")
  assert_true(result.df == 2)
  assert_true(result.f > 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(result.ss > 0.0)
}

test "rstatix_friedman_test equal groups" {
  let groups = [
    [5.0, 5.0, 5.0, 5.0, 5.0],
    [5.0, 5.0, 5.0, 5.0, 5.0],
    [5.0, 5.0, 5.0, 5.0, 5.0],
  ]
  let result = @src.rstatix_friedman_test(groups)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_friedman_test too few groups" {
  let groups = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]
  let result = @src.rstatix_friedman_test(groups)
  assert_true(result.p_value.is_nan())
  assert_true(result.f.is_nan())
}

test "rstatix_friedman_test small blocks" {
  let groups = [
    [1.0, 4.0],
    [2.0, 5.0],
    [3.0, 6.0],
  ]
  let result = @src.rstatix_friedman_test(groups)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

// ===== rstatix_bh_correct =====

test "rstatix_bh_correct basic" {
  let p_values = [0.01, 0.04, 0.03, 0.005, 0.2]
  let result = @src.rstatix_bh_correct(p_values)
  assert_true(result.length() == 5)
  let mut i = 0
  while i < 5 {
    assert_true(result[i] >= p_values[i] - 1.0e-10)
    assert_true(result[i] <= 1.0)
    i = i + 1
  }
}

test "rstatix_bh_correct empty" {
  let p_values : Array[Double] = Array::new()
  let result = @src.rstatix_bh_correct(p_values)
  assert_true(result.length() == 0)
}

test "rstatix_bh_correct single" {
  let p_values = [0.05]
  let result = @src.rstatix_bh_correct(p_values)
  assert_true(result.length() == 1)
  assert_true(result[0] >= 0.05 - 1.0e-10 && result[0] <= 0.05 + 1.0e-10)
}

test "rstatix_bh_correct all significant" {
  let p_values = [0.001, 0.002, 0.003, 0.004, 0.005]
  let result = @src.rstatix_bh_correct(p_values)
  assert_true(result.length() == 5)
  let mut i = 0
  while i < 5 {
    assert_true(result[i] <= 1.0)
    i = i + 1
  }
}

test "rstatix_bh_correct all non-significant" {
  let p_values = [0.8, 0.9, 0.95, 0.98, 0.99]
  let result = @src.rstatix_bh_correct(p_values)
  assert_true(result.length() == 5)
  let mut i = 0
  while i < 5 {
    assert_true(result[i] >= p_values[i] - 1.0e-10)
    assert_true(result[i] <= 1.0)
    i = i + 1
  }
}

test "rstatix_bh_correct preserves order" {
  let p_values = [0.2, 0.005, 0.04, 0.01, 0.03]
  let result = @src.rstatix_bh_correct(p_values)
  assert_true(result.length() == 5)
  assert_true(result[0] >= p_values[0] - 1.0e-10)
  assert_true(result[1] >= p_values[1] - 1.0e-10)
}

// ===== rstatix_bonferroni =====

test "rstatix_bonferroni basic" {
  let p_values = [0.01, 0.04, 0.03, 0.005, 0.2]
  let result = @src.rstatix_bonferroni(p_values)
  assert_true(result.length() == 5)
  assert_true(result[0] >= 0.05 - 1.0e-10)
  assert_true(result[1] >= 0.2 - 1.0e-10)
  assert_true(result[4] >= 1.0 - 1.0e-10)
}

test "rstatix_bonferroni empty" {
  let p_values : Array[Double] = Array::new()
  let result = @src.rstatix_bonferroni(p_values)
  assert_true(result.length() == 0)
}

test "rstatix_bonferroni single" {
  let p_values = [0.05]
  let result = @src.rstatix_bonferroni(p_values)
  assert_true(result.length() == 1)
  assert_true(result[0] >= 0.05 - 1.0e-10 && result[0] <= 0.05 + 1.0e-10)
}

test "rstatix_bonferroni capped at 1" {
  let p_values = [0.3, 0.4, 0.5]
  let result = @src.rstatix_bonferroni(p_values)
  assert_true(result[0] >= 0.9 - 1.0e-10)
  assert_true(result[1] >= 1.0 - 1.0e-10)
  assert_true(result[2] >= 1.0 - 1.0e-10)
}

// ===== Summary functions =====

test "rstatix_t_test_summary basic" {
  let x = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0]
  let result = @src.rstatix_t_test(x, mu=15.0)
  let summary = @src.rstatix_t_test_summary(result)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("T-test"))
  assert_true(summary.contains("Groups"))
  assert_true(summary.contains("p-value"))
  assert_true(summary.contains("95% CI"))
  assert_true(summary.contains("Method"))
}

test "rstatix_t_test_summary two-sample" {
  let x = [10.0, 12.0, 15.0, 18.0, 20.0]
  let y = [8.0, 9.0, 11.0, 12.0, 14.0]
  let result = @src.rstatix_t_test(x, y=y)
  let summary = @src.rstatix_t_test_summary(result)
  assert_true(summary.contains("Group1"))
  assert_true(summary.contains("Group2"))
}

test "rstatix_cor_test_summary basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
  let result = @src.rstatix_cor_test(x, y)
  let summary = @src.rstatix_cor_test_summary(result)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("Correlation Test"))
  assert_true(summary.contains("Statistic"))
  assert_true(summary.contains("p-value"))
  assert_true(summary.contains("n"))
  assert_true(summary.contains("Method"))
}

test "rstatix_cor_test_summary spearman" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
  let result = @src.rstatix_cor_test(x, y, cor_method="spearman")
  let summary = @src.rstatix_cor_test_summary(result)
  assert_true(summary.contains("spearman"))
}

test "rstatix_anova_summary basic" {
  let groups = [
    [10.0, 12.0, 15.0, 18.0, 20.0],
    [8.0, 9.0, 11.0, 12.0, 14.0],
    [5.0, 6.0, 7.0, 8.0, 9.0],
  ]
  let result = @src.rstatix_anova_test(groups)
  let summary = @src.rstatix_anova_summary(result)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("ANOVA Result"))
  assert_true(summary.contains("Source"))
  assert_true(summary.contains("df"))
  assert_true(summary.contains("SS"))
  assert_true(summary.contains("MS"))
  assert_true(summary.contains("F"))
  assert_true(summary.contains("p-value"))
}

test "rstatix_anova_summary kruskal" {
  let groups = [
    [10.0, 12.0, 15.0, 18.0, 20.0],
    [8.0, 9.0, 11.0, 12.0, 14.0],
    [5.0, 6.0, 7.0, 8.0, 9.0],
  ]
  let result = @src.rstatix_kruskal_test(groups)
  let summary = @src.rstatix_anova_summary(result)
  assert_true(summary.contains("Kruskal-Wallis"))
}

test "rstatix_anova_summary friedman" {
  let groups = [
    [10.0, 12.0, 15.0, 18.0, 20.0],
    [8.0, 9.0, 11.0, 12.0, 14.0],
    [5.0, 6.0, 7.0, 8.0, 9.0],
  ]
  let result = @src.rstatix_friedman_test(groups)
  let summary = @src.rstatix_anova_summary(result)
  assert_true(summary.contains("Friedman"))
}

// ===== rstatix_sample_data =====

test "rstatix_sample_data basic" {
  let data = @src.rstatix_sample_data()
  assert_true(data.length() == 3)
  assert_true(data[0].length() == 10)
  assert_true(data[1].length() == 10)
  assert_true(data[2].length() == 10)
}

test "rstatix_sample_data values" {
  let data = @src.rstatix_sample_data()
  assert_true(data[0][0] == 10.0)
  assert_true(data[0][9] == 32.0)
  assert_true(data[1][0] == 8.0)
  assert_true(data[1][9] == 22.0)
  assert_true(data[2][0] == 5.0)
  assert_true(data[2][9] == 14.0)
}

test "rstatix_sample_data used in anova" {
  let data = @src.rstatix_sample_data()
  let result = @src.rstatix_anova_test(data)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(result.f > 0.0)
}

test "rstatix_sample_data used in kruskal" {
  let data = @src.rstatix_sample_data()
  let result = @src.rstatix_kruskal_test(data)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_sample_data used in friedman" {
  let data = @src.rstatix_sample_data()
  let result = @src.rstatix_friedman_test(data)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

// ===== Edge cases and integration =====

test "rstatix_t_test constant values" {
  let x = [5.0, 5.0, 5.0, 5.0, 5.0]
  let result = @src.rstatix_t_test(x, mu=5.0)
  assert_true(result.statistic.is_nan())
  assert_true(result.p_value.is_nan())
}

test "rstatix_wilcox_test all zeros" {
  let x = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  let result = @src.rstatix_wilcox_test(x)
  assert_true(result.p_value.is_nan())
}

test "rstatix_cor_test identical values" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let result = @src.rstatix_cor_test(x, x)
  assert_true(result.statistic > 0.99 && result.statistic <= 1.0)
}

test "rstatix_bh_correct preserves non-decreasing order" {
  let p_values = [0.001, 0.01, 0.03, 0.05, 0.1]
  let result = @src.rstatix_bh_correct(p_values)
  let mut i = 1
  while i < result.length() {
    assert_true(result[i] >= result[i - 1] - 1.0e-10)
    i = i + 1
  }
}

test "rstatix_t_test_two_sample_same_values" {
  let x = [10.0, 10.0, 10.0, 10.0, 10.0]
  let y = [10.0, 10.0, 10.0, 10.0, 10.0]
  let result = @src.rstatix_t_test(x, y=y)
  assert_true(result.statistic.is_nan())
}

test "rstatix_anova_test with sample data from function" {
  let data = @src.rstatix_sample_data()
  let result = @src.rstatix_anova_test(data, group_names=["GroupA", "GroupB", "GroupC"])
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_t_test paired with mu default" {
  let x = [10.0, 12.0, 15.0, 18.0, 20.0, 22.0, 25.0, 28.0, 30.0, 32.0]
  let y = [8.0, 9.0, 11.0, 12.0, 14.0, 15.0, 17.0, 19.0, 20.0, 22.0]
  let result = @src.rstatix_t_test(x, y=y, paired=true)
  assert_true(result.statistic > 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_wilcox_test two-sample with different lengths" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let y = [9.0, 10.0, 11.0, 12.0, 13.0, 14.0]
  let result = @src.rstatix_wilcox_test(x, y=y)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

test "rstatix_cor_test pearson default method" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
  let result = @src.rstatix_cor_test(x, y)
  assert_eq(result.method_name, "Pearson's product-moment correlation")
}

test "rstatix_kruskal_test with sample data" {
  let data = @src.rstatix_sample_data()
  let result = @src.rstatix_kruskal_test(data, _group_names=["A", "B", "C"])
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(result.source == "Kruskal-Wallis")
}

test "rstatix_friedman_test with sample data" {
  let data = @src.rstatix_sample_data()
  let result = @src.rstatix_friedman_test(data, _group_names=["A", "B", "C"])
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
  assert_true(result.source == "Friedman")
}