///|
/// Tests for the LOWESS smoothing module.

///|
test "lowess_tricube_weight_at_zero" {
  let result = @src.lowess_tricube_weight(0.0)
  assert_eq(result, 1.0)
}

///|
test "lowess_tricube_weight_at_one" {
  let result = @src.lowess_tricube_weight(1.0)
  assert_eq(result, 0.0)
}

///|
test "lowess_tricube_weight_above_one" {
  let result = @src.lowess_tricube_weight(1.5)
  assert_eq(result, 0.0)
}

///|
test "lowess_tricube_weight_at_half" {
  let result = @src.lowess_tricube_weight(0.5)
  let expected = (1.0 - 0.5 * 0.5 * 0.5)
  let expected_cubed = expected * expected * expected
  assert_true((result - expected_cubed).abs() < 0.001)
}

///|
test "lowess_tricube_weight_negative" {
  let result = @src.lowess_tricube_weight(-0.5)
  assert_eq(result, 0.0)
}

///|
test "lowess_bisquare_weight_at_zero" {
  let result = @src.lowess_bisquare_weight(0.0)
  assert_eq(result, 1.0)
}

///|
test "lowess_bisquare_weight_at_one" {
  let result = @src.lowess_bisquare_weight(1.0)
  assert_eq(result, 0.0)
}

///|
test "lowess_bisquare_weight_above_one" {
  let result = @src.lowess_bisquare_weight(2.0)
  assert_eq(result, 0.0)
}

///|
test "lowess_bisquare_weight_at_half" {
  let result = @src.lowess_bisquare_weight(0.5)
  let expected = (1.0 - 0.5 * 0.5)
  let expected_squared = expected * expected
  assert_true((result - expected_squared).abs() < 0.001)
}

///|
test "lowess_median_odd_length" {
  let result = @src.lowess_median([1.0, 3.0, 2.0])
  assert_eq(result, 2.0)
}

///|
test "lowess_median_even_length" {
  let result = @src.lowess_median([1.0, 2.0, 3.0, 4.0])
  assert_eq(result, 2.5)
}

///|
test "lowess_median_empty" {
  let result = @src.lowess_median([])
  assert_eq(result, 0.0)
}

///|
test "lowess_median_single_element" {
  let result = @src.lowess_median([42.0])
  assert_eq(result, 42.0)
}

///|
test "lowess_weighted_linear_regression_linear" {
  let x = [0.0, 1.0, 2.0, 3.0]
  let y = [1.0, 3.0, 5.0, 7.0]
  let w = [1.0, 1.0, 1.0, 1.0]
  let (intercept, slope) = @src.lowess_weighted_linear_regression(x, y, w)
  assert_true((intercept - 1.0).abs() < 0.001)
  assert_true((slope - 2.0).abs() < 0.001)
}

///|
test "lowess_weighted_linear_regression_weighted" {
  let x = [0.0, 1.0, 2.0, 3.0]
  let y = [1.0, 3.0, 5.0, 7.0]
  let w = [1.0, 1.0, 1.0, 1.0]
  let (intercept, slope) = @src.lowess_weighted_linear_regression(x, y, w)
  let (intercept2, slope2) = @src.lowess_weighted_linear_regression(x, y, [2.0, 2.0, 2.0, 2.0])
  assert_true((intercept - intercept2).abs() < 0.001)
  assert_true((slope - slope2).abs() < 0.001)
}

///|
test "lowess_weighted_linear_regression_empty" {
  let (intercept, slope) = @src.lowess_weighted_linear_regression([], [], [])
  assert_eq(intercept, 0.0)
  assert_eq(slope, 0.0)
}

///|
test "lowess_weighted_linear_regression_single_point" {
  let (intercept, slope) = @src.lowess_weighted_linear_regression([1.0], [5.0], [1.0])
  assert_eq(intercept, 5.0)
  assert_eq(slope, 0.0)
}

///|
test "lowess_linear_data" {
  let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
  let smoothed = @src.lowess(x, y)
  assert_eq(smoothed.length(), 6)
  for i = 0; i < 6; i = i + 1 {
    assert_true((smoothed[i] - y[i]).abs() < 0.01)
  }
}

///|
test "lowess_empty_arrays" {
  let smoothed = @src.lowess([], [])
  assert_eq(smoothed.length(), 0)
}

///|
test "lowess_single_point" {
  let smoothed = @src.lowess([1.0], [5.0])
  assert_eq(smoothed.length(), 1)
  assert_eq(smoothed[0], 5.0)
}

///|
test "lowess_nonlinear_data" {
  let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
  let y = [0.0, 1.0, 4.0, 9.0, 16.0, 25.0, 36.0, 49.0, 64.0, 81.0]
  let smoothed = @src.lowess(x, y)
  assert_eq(smoothed.length(), 10)
  for i = 0; i < 10; i = i + 1 {
    assert_true((smoothed[i] - y[i]).abs() < 2.0)
  }
}

///|
test "lowess_with_custom_f" {
  let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
  let smoothed = @src.lowess(x, y, f=0.5)
  assert_eq(smoothed.length(), 6)
  for i = 0; i < 6; i = i + 1 {
    assert_true((smoothed[i] - y[i]).abs() < 0.01)
  }
}

///|
test "lowess_robust_with_outlier" {
  let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
  let y = [0.0, 2.0, 4.0, 6.0, 8.0, 100.0, 12.0, 14.0, 16.0, 18.0]
  let non_robust = @src.lowess(x, y, n_steps=0)
  let robust = @src.lowess_robust(x, y)
  assert_eq(robust.length(), 10)
  let nr_err = (non_robust[5] - 10.0).abs()
  let r_err = (robust[5] - 10.0).abs()
  assert_true(r_err <= nr_err)
}

///|
test "lowess_robust_linear_data" {
  let x = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0]
  let smoothed = @src.lowess_robust(x, y)
  assert_eq(smoothed.length(), 6)
  for i = 0; i < 6; i = i + 1 {
    assert_true((smoothed[i] - y[i]).abs() < 0.1)
  }
}

///|
test "lowess_robust_empty" {
  let smoothed = @src.lowess_robust([], [])
  assert_eq(smoothed.length(), 0)
}

///|
test "lowess_smooth_series_linear" {
  let values = [1.0, 2.0, 3.0, 4.0, 5.0]
  let smoothed = @src.lowess_smooth_series(values, window_frac=0.6)
  assert_eq(smoothed.length(), 5)
  for i = 0; i < 5; i = i + 1 {
    assert_true((smoothed[i] - values[i]).abs() < 0.5)
  }
}

///|
test "lowess_smooth_series_empty" {
  let smoothed = @src.lowess_smooth_series([])
  assert_eq(smoothed.length(), 0)
}

///|
test "lowess_smooth_series_single" {
  let smoothed = @src.lowess_smooth_series([42.0])
  assert_eq(smoothed.length(), 1)
  assert_eq(smoothed[0], 42.0)
}

///|
test "lowess_predict_interpolation" {
  let x = [0.0, 1.0, 2.0, 3.0, 4.0]
  let y = [0.0, 2.0, 4.0, 6.0, 8.0]
  let smoothed_y = [0.0, 2.0, 4.0, 6.0, 8.0]
  let result = @src.lowess_predict(2.5, x, y, smoothed_y)
  assert_true((result - 5.0).abs() < 0.001)
}

///|
test "lowess_predict_at_boundary_left" {
  let x = [0.0, 1.0, 2.0, 3.0, 4.0]
  let y = [0.0, 2.0, 4.0, 6.0, 8.0]
  let smoothed_y = [0.0, 2.0, 4.0, 6.0, 8.0]
  let result = @src.lowess_predict(0.0, x, y, smoothed_y)
  assert_eq(result, 0.0)
}

///|
test "lowess_predict_at_boundary_right" {
  let x = [0.0, 1.0, 2.0, 3.0, 4.0]
  let y = [0.0, 2.0, 4.0, 6.0, 8.0]
  let smoothed_y = [0.0, 2.0, 4.0, 6.0, 8.0]
  let result = @src.lowess_predict(4.0, x, y, smoothed_y)
  assert_eq(result, 8.0)
}

///|
test "lowess_predict_extrapolation" {
  let x = [0.0, 1.0, 2.0, 3.0, 4.0]
  let y = [0.0, 2.0, 4.0, 6.0, 8.0]
  let smoothed_y = [0.0, 2.0, 4.0, 6.0, 8.0]
  let result = @src.lowess_predict(5.0, x, y, smoothed_y)
  assert_true((result - 10.0).abs() < 1.0)
}

///|
test "lowess_predict_empty" {
  let result = @src.lowess_predict(1.0, [], [], [])
  assert_eq(result, 0.0)
}

///|
test "lowess_predict_single_point" {
  let result = @src.lowess_predict(1.0, [0.0], [5.0], [5.0])
  assert_eq(result, 5.0)
}
