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
