///|
test "impute_param_new" {
  let p = @src.KNNImputeParam::new()
  assert_eq(p.k, 10)
  assert_eq(p.by_row, true)
  assert_eq(p.eps, 0.000001)
  assert_true(p.max_value.is_none())
  assert_true(p.min_value.is_none())
}

///|
test "impute_param_with_k" {
  let p = @src.KNNImputeParam::with_k(5)
  assert_eq(p.k, 5)
}

///|
test "impute_by_row_mean_basic" {
  let data = [
    [1.0, 2.0, Double::nan(), 4.0],
    [5.0, Double::nan(), 7.0, 8.0],
    [9.0, 10.0, 11.0, Double::nan()],
  ]
  let imputed = @src.impute_by_row_mean(data)
  assert_true(imputed[0][2] > 2.33 && imputed[0][2] < 2.334)
  assert_true(imputed[1][1] > 6.666 && imputed[1][1] < 6.667)
  assert_eq(imputed[2][3], 10.0)
}

///|
test "impute_by_col_median_basic" {
  let data = [
    [1.0, 2.0, 3.0, Double::nan()],
    [Double::nan(), 4.0, 5.0, 6.0],
    [7.0, 8.0, Double::nan(), 9.0],
  ]
  let imputed = @src.impute_by_col_median(data)
  assert_eq(imputed[1][0], 4.0)
  assert_eq(imputed[0][3], 7.5)
}

///|
test "impute_locf_basic" {
  let data = [
    [Double::nan(), 1.0, Double::nan(), 2.0, Double::nan()],
  ]
  let imputed = @src.impute_locf(data, by_row=true)
  assert_true(imputed[0][0].is_nan())
  assert_eq(imputed[0][2], 1.0)
  assert_eq(imputed[0][4], 2.0)
}

///|
test "impute_nocb_basic" {
  let data = [
    [Double::nan(), 1.0, Double::nan(), 2.0, Double::nan()],
  ]
  let imputed = @src.impute_nocb(data, by_row=true)
  assert_eq(imputed[0][0], 1.0)
  assert_eq(imputed[0][2], 2.0)
  assert_true(imputed[0][4].is_nan())
}

///|
test "impute_na_by_zero" {
  let data = [
    [1.0, Double::nan()],
    [Double::nan(), 2.0],
  ]
  let imputed = @src.impute_na_by_zero(data)
  assert_eq(imputed[0][1], 0.0)
  assert_eq(imputed[1][0], 0.0)
  assert_eq(imputed[0][0], 1.0)
  assert_eq(imputed[1][1], 2.0)
}

///|
test "impute_na_stats" {
  let data = [
    [1.0, Double::nan(), 3.0],
    [Double::nan(), Double::nan(), 6.0],
  ]
  let stats = @src.impute_na_stats(data)
  assert_eq(stats.total_rows, 2)
  assert_eq(stats.total_cols, 3)
  assert_eq(stats.total_na, 3)
  assert_true((stats.na_fraction - 0.5).abs() < 0.001)
  assert_eq(stats.row_na_counts.length(), 2)
  assert_eq(stats.row_na_counts[0], 1)
  assert_eq(stats.row_na_counts[1], 2)
}

///|
test "impute_na_summary" {
  let data = [
    [1.0, Double::nan()],
    [Double::nan(), 2.0],
  ]
  let stats = @src.impute_na_stats(data)
  let s = @src.impute_na_summary(stats)
  assert_true(s.contains("total NAs"))
}

///|
test "impute_by_knn_no_na" {
  let data = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]
  let p = @src.KNNImputeParam::new()
  let imputed = @src.impute_by_knn(data, p)
  assert_eq(imputed[0][0], 1.0)
  assert_eq(imputed[1][2], 6.0)
}

///|
test "impute_by_knn_simple" {
  let data = [
    [1.0, Double::nan(), 3.0],
    [2.0, 5.0, 4.0],
    [3.0, 6.0, 7.0],
  ]
  let imputed = @src.impute_by_knn_simple(data, 2)
  assert_true(!imputed[0][1].is_nan())
}

///|
test "impute_na_mask" {
  let data = [
    [1.0, Double::nan()],
    [Double::nan(), 2.0],
  ]
  let mask = @src.make_na_mask(data)
  assert_eq(mask[0][1], true)
  assert_eq(mask[1][0], true)
  assert_eq(mask[0][0], false)
  assert_eq(mask[1][1], false)
}

///|
test "impute_by_row_median_and_min" {
  let data = [
    [1.0, Double::nan(), 3.0, 100.0],
  ]
  let med = @src.impute_by_row_median(data)
  assert_true(!med[0][1].is_nan())
  let mn = @src.impute_by_row_min(data)
  assert_eq(mn[0][1], 1.0)
}
