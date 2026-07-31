// Tests for Bio.Statistics module

test "stat_mean" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.stat_mean(data)
  assert_true(result > 2.9 && result < 3.1)
}

test "stat_mean empty" {
  let data : Array[Double] = Array::new()
  let result = @src.stat_mean(data)
  assert_true(result == 0.0)
}

test "stat_variance" {
  let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
  let result = @src.stat_variance(data)
  assert_true(result > 4.0 && result < 5.0)
}

test "stat_std" {
  let data = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
  let result = @src.stat_std(data)
  assert_true(result > 2.0 && result < 2.5)
}

test "stat_median odd" {
  let data = [1.0, 3.0, 2.0, 5.0, 4.0]
  let result = @src.stat_median(data)
  assert_true(result == 3.0)
}

test "stat_median even" {
  let data = [1.0, 2.0, 3.0, 4.0]
  let result = @src.stat_median(data)
  assert_true(result == 2.5)
}

test "stat_min" {
  let data = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0]
  let result = @src.stat_min(data)
  assert_true(result == 1.0)
}

test "stat_max" {
  let data = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0]
  let result = @src.stat_max(data)
  assert_true(result == 9.0)
}

test "stat_sum" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.stat_sum(data)
  assert_true(result == 15.0)
}

test "stat_cumsum" {
  let data = [1.0, 2.0, 3.0]
  let result = @src.stat_cumsum(data)
  assert_true(result.length() == 3)
  assert_true(result[0] == 1.0)
  assert_true(result[1] == 3.0)
  assert_true(result[2] == 6.0)
}

test "stat_pearson_correlation positive" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let result = @src.stat_pearson_correlation(x, y)
  assert_true(result > 0.99 && result < 1.01)
}

test "stat_pearson_correlation negative" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [10.0, 8.0, 6.0, 4.0, 2.0]
  let result = @src.stat_pearson_correlation(x, y)
  assert_true(result < -0.99 && result > -1.01)
}

test "stat_zscore" {
  let result = @src.stat_zscore(3.0, 2.0, 0.5)
  assert_true(result > 1.9 && result < 2.1)
}

test "stat_normal_cdf at 0" {
  let result = @src.stat_normal_cdf(0.0)
  assert_true(result > 0.49 && result < 0.51)
}

test "stat_t_statistic" {
  let sample = [2.0, 3.0, 4.0, 5.0, 6.0]
  let result = @src.stat_t_statistic(sample, 3.0)
  assert_true(result > 1.4 && result < 1.6)
}

test "stat_confidence_interval" {
  let data = [2.0, 3.0, 4.0, 5.0, 6.0]
  let result = @src.stat_confidence_interval(data)
  assert_true(result.length() == 2)
  assert_true(result[0] < result[1])
}

test "stat_mann_whitney_u_basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let result = @src.stat_mann_whitney(x, y)
  assert_true(result.1 >= 0.0 && result.1 <= 1.0)
}

test "stat_wilcoxon_signed_rank_basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [2.0, 4.0, 6.0, 8.0, 10.0]
  let result = @src.stat_wilcoxon_signed_rank(x, y)
  assert_true(result.1 >= 0.0 && result.1 <= 1.0)
}

test "stat_ks_test_basic" {
  let sample1 = [1.0, 2.0, 3.0, 4.0, 5.0]
  let sample2 = [2.0, 4.0, 6.0, 8.0, 10.0]
  let result = @src.stat_ks_test(sample1, sample2)
  assert_true(result.1 >= 0.0 && result.1 <= 1.0)
}

test "stat_fisher_exact_basic" {
  let table = [[10, 2], [3, 5]]
  let result = @src.stat_fisher_exact(table)
  assert_true(result.1 >= 0.0 && result.1 <= 1.0)
}

test "stat_bonferroni_correct_basic" {
  let p_values = [0.01, 0.04, 0.03, 0.005, 0.2]
  let result = @src.stat_bonferroni(p_values)
  assert_true(result.length() == 5)
  assert_true(result[0] >= 0.01)
}

test "stat_holm_correct_basic" {
  let p_values = [0.01, 0.04, 0.03, 0.005, 0.2]
  let result = @src.stat_holm(p_values)
  assert_true(result.length() == 5)
}

test "stat_by_correct_basic" {
  let p_values = [0.01, 0.04, 0.03, 0.005, 0.2]
  let result = @src.stat_by(p_values)
  assert_true(result.length() == 5)
}

test "stat_chi2_cdf_basic" {
  // For df=3, P(chi2 <= 7.81) ~ 0.95
  let result = @src.stat_chi2_cdf(7.81, 3)
  assert_true(result > 0.90 && result < 0.99)
}

test "stat_chi2_cdf_zero" {
  let result = @src.stat_chi2_cdf(0.0, 3)
  assert_true(result == 0.0)
}

test "stat_chi2_cdf_negative" {
  let result = @src.stat_chi2_cdf(-1.0, 5)
  assert_true(result == 0.0)
}

test "stat_chi2_quantile_basic" {
  // For df=2, chi2_0.95 ~ 5.99
  let result = @src.stat_chi2_quantile(0.95, 2)
  assert_true(result > 5.0 && result < 7.0)
}

test "stat_chi2_quantile_zero_one" {
  assert_true(@src.stat_chi2_quantile(0.0, 2) == 0.0)
  assert_true(@src.stat_chi2_quantile(1.0, 2) > 1.0e100)
}

test "stat_t_cdf_central" {
  let result = @src.stat_t_cdf(0.0, 10.0)
  assert_true(result > 0.49 && result < 0.51)
}

test "stat_t_cdf_positive" {
  // t with df=10 at 2.228 ~ 0.975
  let result = @src.stat_t_cdf(2.228, 10.0)
  assert_true(result > 0.94 && result < 0.99)
}

test "stat_t_quantile_central" {
  let result = @src.stat_t_quantile(0.5, 10.0)
  assert_true(result > -0.01 && result < 0.01)
}

test "stat_t_quantile_basic" {
  // t_0.975,df=10 ~ 2.228
  let result = @src.stat_t_quantile(0.975, 10.0)
  assert_true(result > 1.5 && result < 3.0)
}

test "stat_t_test_one_sample" {
  let sample = [2.1, 2.3, 2.2, 2.4, 2.0]
  let result = @src.stat_t_test_one_sample(sample, 2.0)
  assert_true(result.0 > 0.0)
  assert_true(result.1 > 0.0 && result.1 <= 1.0)
}

test "stat_t_test_two_sample" {
  let x = [2.1, 2.3, 2.2, 2.4, 2.0]
  let y = [3.1, 3.3, 3.2, 3.4, 3.0]
  let result = @src.stat_t_test_two_sample(x, y)
  assert_true(result.1 >= 0.0 && result.1 <= 1.0)
}

test "stat_bh_basic" {
  let p_values = [0.01, 0.04, 0.03, 0.005, 0.2]
  let result = @src.stat_bh(p_values)
  assert_true(result.length() == 5)
  // BH should produce values <= original p-values and <= bonferroni
  let bonf = @src.stat_bonferroni(p_values)
  let mut i = 0
  while i < 5 {
    assert_true(result[i] <= bonf[i] + 1.0e-10)
    i = i + 1
  }
}

test "stat_bh_empty" {
  let p_values : Array[Double] = Array::new()
  let result = @src.stat_bh(p_values)
  assert_true(result.length() == 0)
}

test "stat_normal_quantile_half" {
  let result = @src.stat_normal_quantile(0.5)
  assert_true(result > -0.01 && result < 0.01)
}

test "stat_normal_quantile_standard" {
  // qnorm(0.975) ~ 1.96
  let result = @src.stat_normal_quantile(0.975)
  assert_true(result > 1.9 && result < 2.1)
}

test "stat_mad_basic" {
  let data = [1.0, 2.0, 3.0, 4.0, 5.0]
  let result = @src.stat_mad(data)
  // Median is 3, deviations are [2,1,0,1,2], median=1, with constant ~1.48
  assert_true(result > 1.0 && result < 2.0)
}

test "stat_anova_basic" {
  let samples = [
    [1.0, 2.0, 3.0],
    [2.0, 3.0, 4.0],
    [8.0, 9.0, 10.0],
  ]
  let result = @src.stat_anova(samples)
  assert_true(result.0 > 0.0)
  assert_true(result.1 >= 0.0 && result.1 <= 1.0)
}

test "stat_anova_equal" {
  let samples = [
    [1.0, 1.0, 1.0],
    [1.0, 1.0, 1.0],
    [1.0, 1.0, 1.0],
  ]
  let result = @src.stat_anova(samples)
  assert_true(result.0 >= 0.0)
}

test "stat_logrank_test_basic" {
  let time = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let event = [true, true, true, false, true, false]
  let group = [0, 0, 0, 1, 1, 1]
  let result = @src.stat_logrank_test(time, event, group)
  assert_true(result.0 >= 0.0)
  assert_true(result.1 >= 0.0 && result.1 <= 1.0)
}
