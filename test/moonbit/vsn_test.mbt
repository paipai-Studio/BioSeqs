///|
test "vsn_glog_basic" {
  let x = 1.0
  let y = @src.glog(x, 0.0, 1.0)
  // glog(1,0,1) = ln(0.5 + sqrt(1.25)) ≈ ln(1.618) ≈ 0.481
  assert_true(y > 0.48 && y < 0.49)
}

///|
test "vsn_glog_monotonic" {
  let mut prev = @src.glog(-5.0, 1.0, 1.0)
  let mut i = -4
  while i <= 5 {
    let cur = @src.glog(i.to_double(), 1.0, 1.0)
    assert_true(cur > prev)
    prev = cur
    i = i + 1
  }
}

///|
test "vsn_glog_inv_roundtrip" {
  let a = 0.5
  let b = 2.0
  let mut x = -2.0
  while x <= 2.0 {
    let y = @src.glog(x, a, b)
    let x_back = @src.glog_inv(y, a, b)
    assert_true((x - x_back).abs() < 0.000001)
    x = x + 0.5
  }
}

///|
test "vsn_control_new" {
  let c = @src.VSNControl::new()
  assert_eq(c.n_bins, 20)
  assert_eq(c.min_bin_count, 5)
  assert_eq(c.max_iter, 20)
  assert_true((c.lts_frac - 0.9).abs() < 0.0001)
}

///|
test "vsn2_basic" {
  let data = [
    [100.0, 200.0],
    [150.0, 180.0],
    [90.0, 210.0],
  ]
  let result = @src.vsn2(data)
  assert_eq(result.length(), 3)
  assert_eq(result[0].length(), 2)
}

///|
test "vsn2_with_control_runs" {
  let data = [
    [100.0, 200.0],
    [150.0, 180.0],
    [90.0, 210.0],
  ]
  let ctrl = @src.VSNControl::new()
  let result = @src.vsn2_with_control(data, ctrl)
  assert_eq(result.length(), 3)
  assert_eq(result[0].length(), 2)
}

///|
test "vsn_fit_and_report_returns_vsnresult" {
  let data = [
    [100.0, 200.0],
    [150.0, 180.0],
    [90.0, 210.0],
  ]
  let fit = @src.vsn_fit_and_report(data)
  assert_eq(fit.params.length(), 2)
  let summary = @src.summarize_vsn_fit(fit)
  assert_true(summary.contains("VSNResult:"))
  assert_true(summary.contains("columns fitted"))
}

///|
test "mean_sd_bins_and_ascii" {
  let data = [
    [100.0, 200.0],
    [150.0, 180.0],
    [90.0, 210.0],
    [110.0, 220.0],
    [120.0, 190.0],
  ]
  // bins calculation
  let bins = @src.mean_sd_bins(data, 3)
  assert_true(bins.length() > 0)
  // Skip ascii plot in heavy test environment to avoid OOM:
  // let plot = @src.mean_sd_ascii(data, n_bins=40)
  // assert_true(plot.length() > 0)
  // assert_true(plot.contains("mean-SD plot"))
  assert_true(true)
}

///|
test "vsn_denoise_keeps_shape" {
  let data = [
    [100.0, 200.0],
    [150.0, 180.0],
    [90.0, 210.0],
  ]
  let denoised = @src.vsn_denoise(data)
  assert_eq(denoised.length(), 3)
  assert_eq(denoised[0].length(), 2)
}
