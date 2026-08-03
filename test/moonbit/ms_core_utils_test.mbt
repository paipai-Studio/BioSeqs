///|
/// Test file for MsCoreUtils module.

// ============================================================================
// Numeric helpers: mc_mean
// ============================================================================

test "mc_mean_basic" {
  let values = [1.0, 2.0, 3.0, 4.0, 5.0]
  let m = @src.mc_mean(values)
  assert_true((m - 3.0).abs() < 0.001)
}

test "mc_mean_empty" {
  let empty : Array[Double] = []
  let m = @src.mc_mean(empty)
  assert_eq(m, 0.0)
}

test "mc_mean_two_values" {
  let m = @src.mc_mean([2.5, 3.5])
  assert_true((m - 3.0).abs() < 0.001)
}

// ============================================================================
// Numeric helpers: mc_median
// ============================================================================

test "mc_median_odd" {
  let m = @src.mc_median([1.0, 3.0, 5.0])
  assert_eq(m, 3.0)
}

test "mc_median_even" {
  let m = @src.mc_median([1.0, 2.0, 3.0, 4.0])
  assert_true((m - 2.5).abs() < 0.001)
}

test "mc_median_unsorted" {
  let m = @src.mc_median([5.0, 1.0, 3.0, 2.0, 4.0])
  assert_eq(m, 3.0)
}

test "mc_median_empty" {
  let empty : Array[Double] = []
  let m = @src.mc_median(empty)
  assert_eq(m, 0.0)
}

// ============================================================================
// Numeric helpers: mc_sd
// ============================================================================

test "mc_sd_basic" {
  // mean=3, var = (4+1+0+1+4)/4 = 2.5, sd = sqrt(2.5) ~ 1.5811
  let s = @src.mc_sd([1.0, 2.0, 3.0, 4.0, 5.0])
  assert_true((s - 1.5811).abs() < 0.001)
}

test "mc_sd_constant" {
  let s = @src.mc_sd([5.0, 5.0, 5.0])
  assert_eq(s, 0.0)
}

test "mc_sd_single" {
  // n < 2 returns 0.0
  let s = @src.mc_sd([7.0])
  assert_eq(s, 0.0)
}

// ============================================================================
// Numeric helpers: mc_quantile
// ============================================================================

test "mc_quantile_median" {
  let q = @src.mc_quantile([1.0, 2.0, 3.0, 4.0, 5.0], 0.5)
  assert_eq(q, 3.0)
}

test "mc_quantile_min_max" {
  let qmin = @src.mc_quantile([1.0, 2.0, 3.0, 4.0, 5.0], 0.0)
  let qmax = @src.mc_quantile([1.0, 2.0, 3.0, 4.0, 5.0], 1.0)
  assert_eq(qmin, 1.0)
  assert_eq(qmax, 5.0)
}

test "mc_quantile_quarter" {
  // pos = 0.25 * 4 = 1, lo=1, frac=0 -> sorted[1] = 2.0
  let q = @src.mc_quantile([1.0, 2.0, 3.0, 4.0, 5.0], 0.25)
  assert_eq(q, 2.0)
}

// ============================================================================
// Numeric helpers: mc_sum, mc_dot
// ============================================================================

test "mc_sum_basic" {
  let s = @src.mc_sum([1.0, 2.0, 3.0])
  assert_eq(s, 6.0)
}

test "mc_sum_empty" {
  let empty : Array[Double] = []
  let s = @src.mc_sum(empty)
  assert_eq(s, 0.0)
}

test "mc_dot_basic" {
  let d = @src.mc_dot([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])
  // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
  assert_eq(d, 32.0)
}

// ============================================================================
// refineCentroids: mc_refine_centroids
// ============================================================================

test "mc_refine_centroids_symmetric_peak" {
  // Symmetric peak centered at index 2 with m/z 100.2
  let mz = [100.0, 100.1, 100.2, 100.3, 100.4]
  let intensity = [10.0, 20.0, 100.0, 20.0, 10.0]
  let (out_mz, out_int) = @src.mc_refine_centroids(mz, intensity, [2], 2)
  assert_eq(out_mz.length(), 1)
  assert_eq(out_int.length(), 1)
  // weighted mean of m/z should be exactly 100.2 (symmetric)
  assert_true((out_mz[0] - 100.2).abs() < 0.001)
  // intensity is the apex intensity
  assert_eq(out_int[0], 100.0)
}

test "mc_refine_centroids_asymmetric_peak" {
  // Asymmetric: more weight on the right side -> refined m/z shifts right
  let mz = [100.0, 100.1, 100.2, 100.3, 100.4]
  let intensity = [1.0, 1.0, 10.0, 50.0, 50.0]
  let (out_mz, _out_int) = @src.mc_refine_centroids(mz, intensity, [2], 2)
  // weighted mean should be > 100.2 (skewed right)
  assert_true(out_mz[0] > 100.2)
}

test "mc_refine_centroids_out_of_range_index" {
  let mz = [100.0, 100.1, 100.2]
  let intensity = [10.0, 20.0, 30.0]
  // index 5 is out of range -> should be skipped, output empty
  let (out_mz, out_int) = @src.mc_refine_centroids(mz, intensity, [5], 1)
  assert_eq(out_mz.length(), 0)
  assert_eq(out_int.length(), 0)
}

// ============================================================================
// localMaxima: mc_local_maxima
// ============================================================================

test "mc_local_maxima_basic" {
  let idx = @src.mc_local_maxima([1.0, 3.0, 2.0])
  assert_eq(idx, [1])
}

test "mc_local_maxima_zeros" {
  let idx = @src.mc_local_maxima([0.0, 0.0, 0.0])
  assert_eq(idx.length(), 0)
}

test "mc_local_maxima_multiple_peaks" {
  let idx = @src.mc_local_maxima([1.0, 2.0, 1.0, 2.0, 1.0])
  assert_eq(idx, [1, 3])
}

// ============================================================================
// joinPeaks: mc_join_peaks
// ============================================================================

test "mc_join_peaks_absolute_tolerance" {
  let x = [100.0, 200.0, 300.0]
  let y = [100.01, 199.99, 300.5]
  let pairs = @src.mc_join_peaks(x, y, 0.05, false)
  // 100 matches 100.01 (diff 0.01 <= 0.05)
  // 200 matches 199.99 (diff 0.01 <= 0.05)
  // 300 does not match 300.5 (diff 0.5 > 0.05)
  assert_eq(pairs.length(), 2)
  assert_eq(pairs[0], (0, 0))
  assert_eq(pairs[1], (1, 1))
}

test "mc_join_peaks_ppm_tolerance" {
  let x = [1000.0]
  let y = [1000.005]
  // tol = 1000 * 10 * 1e-6 = 0.01; diff = 0.005 <= 0.01 -> match
  let pairs = @src.mc_join_peaks(x, y, 10.0, true)
  assert_eq(pairs.length(), 1)
  assert_eq(pairs[0], (0, 0))
}

test "mc_join_peaks_no_match" {
  let x = [100.0, 200.0]
  let y = [150.0, 250.0]
  let pairs = @src.mc_join_peaks(x, y, 0.01, false)
  assert_eq(pairs.length(), 0)
}

// ============================================================================
// Smoothing: mc_smooth_moving_average
// ============================================================================

test "mc_smooth_moving_average_basic" {
  let out = @src.mc_smooth_moving_average([1.0, 2.0, 3.0, 4.0, 5.0], 1)
  assert_eq(out.length(), 5)
  assert_true((out[0] - 1.5).abs() < 0.001)
  assert_true((out[1] - 2.0).abs() < 0.001)
  assert_true((out[2] - 3.0).abs() < 0.001)
  assert_true((out[3] - 4.0).abs() < 0.001)
  assert_true((out[4] - 4.5).abs() < 0.001)
}

test "mc_smooth_moving_average_reduces_noise" {
  // Noisy signal: moving average should reduce total variation
  let noisy = [1.0, 5.0, 1.0, 5.0, 1.0, 5.0, 1.0]
  let smoothed = @src.mc_smooth_moving_average(noisy, 1)
  let mut rough_in = 0.0
  let mut rough_out = 0.0
  for i in 1..<noisy.length() {
    let d_in = noisy[i] - noisy[i - 1]
    let d_out = smoothed[i] - smoothed[i - 1]
    rough_in = rough_in + d_in * d_in
    rough_out = rough_out + d_out * d_out
  }
  assert_true(rough_out < rough_in)
}

// ============================================================================
// Smoothing: mc_smooth_savitzky_golay
// ============================================================================

test "mc_smooth_savitzky_golay_linear_signal" {
  // A quadratic/SG fit on perfectly linear data reproduces the value at the
  // window center. For the interior point with a full symmetric window, the
  // smoothed value equals the input.
  let lin = [1.0, 2.0, 3.0, 4.0, 5.0]
  let out = @src.mc_smooth_savitzky_golay(lin, 2)
  // Center point i=2 has a full symmetric window (indices 0..4) -> exact.
  assert_true((out[2] - lin[2]).abs() < 0.01)
}

test "mc_smooth_savitzky_golay_reduces_noise" {
  let noisy = [1.0, 5.0, 1.0, 5.0, 1.0, 5.0, 1.0]
  let smoothed = @src.mc_smooth_savitzky_golay(noisy, 2)
  let mut rough_in = 0.0
  let mut rough_out = 0.0
  for i in 1..<noisy.length() {
    let d_in = noisy[i] - noisy[i - 1]
    let d_out = smoothed[i] - smoothed[i - 1]
    rough_in = rough_in + d_in * d_in
    rough_out = rough_out + d_out * d_out
  }
  assert_true(rough_out < rough_in)
}

// ============================================================================
// Smoothing: mc_smooth dispatcher
// ============================================================================

test "mc_smooth_dispatcher_moving_average" {
  let out = @src.mc_smooth([1.0, 2.0, 3.0, 4.0, 5.0], 1, "MovingAverage")
  assert_true((out[2] - 3.0).abs() < 0.001)
}

test "mc_smooth_dispatcher_savitzky_golay" {
  let out = @src.mc_smooth([1.0, 2.0, 3.0, 4.0, 5.0], 2, "SavitzkyGolay")
  // Linear signal -> preserved
  assert_true((out[2] - 3.0).abs() < 0.01)
}

// ============================================================================
// Baseline: mc_baseline_snip
// ============================================================================

test "mc_baseline_snip_flat_signal" {
  // Flat positive signal -> baseline should be <= signal and >= 0
  let signal = [10.0, 10.0, 10.0, 10.0, 10.0]
  let baseline = @src.mc_baseline_snip(signal, 2)
  assert_eq(baseline.length(), 5)
  for i in 0..<signal.length() {
    assert_true(baseline[i] >= 0.0)
    assert_true(baseline[i] <= 10.0 + 0.001)
  }
}

test "mc_baseline_snip_below_peaks" {
  // Signal with a peak in the middle -> baseline at peak should be below peak
  let signal = [1.0, 1.0, 100.0, 1.0, 1.0]
  let baseline = @src.mc_baseline_snip(signal, 2)
  assert_true(baseline[2] < 100.0)
}

// ============================================================================
// Calibration: mc_calibrate
// ============================================================================

test "mc_calibrate_shift" {
  let observed = [100.0, 200.0, 300.0]
  let calibrants = [100.1, 200.1, 300.1]
  let out = @src.mc_calibrate(observed, calibrants, 0.5, "shift")
  assert_eq(out.length(), 3)
  assert_true((out[0] - 100.1).abs() < 0.001)
  assert_true((out[1] - 200.1).abs() < 0.001)
  assert_true((out[2] - 300.1).abs() < 0.001)
}

test "mc_calibrate_linear" {
  // matched pairs: (100, 100.5), (200, 200.5), (300, 300.5)
  // linear fit -> cal = 1.0 * obs + 0.5
  let observed = [100.0, 200.0, 300.0]
  let calibrants = [100.5, 200.5, 300.5]
  let out = @src.mc_calibrate(observed, calibrants, 1.0, "linear")
  assert_eq(out.length(), 3)
  assert_true((out[0] - 100.5).abs() < 0.01)
  assert_true((out[1] - 200.5).abs() < 0.01)
  assert_true((out[2] - 300.5).abs() < 0.01)
}

test "mc_calibrate_no_match_returns_input" {
  // No calibrants within tolerance -> shift method leaves values unchanged
  let observed = [100.0, 200.0]
  let calibrants = [500.0]
  let out = @src.mc_calibrate(observed, calibrants, 0.01, "shift")
  assert_true((out[0] - 100.0).abs() < 0.001)
  assert_true((out[1] - 200.0).abs() < 0.001)
}

// ============================================================================
// Imputation: mc_is_missing and mc_impute
// ============================================================================

test "mc_is_missing_nan" {
  let nan = 0.0 / 0.0
  assert_true(@src.mc_is_missing(nan))
}

test "mc_is_missing_not_nan" {
  assert_true(!@src.mc_is_missing(5.0))
  assert_true(!@src.mc_is_missing(0.0))
  assert_true(!@src.mc_is_missing(-1.5))
}

test "mc_impute_zero" {
  let nan = 0.0 / 0.0
  let matrix = [
    [1.0, nan, 3.0],
    [nan, 5.0, 6.0],
  ]
  let out = @src.mc_impute(matrix, "zero", 0)
  assert_eq(out[0][1], 0.0)
  assert_eq(out[1][0], 0.0)
  // Non-missing values preserved
  assert_eq(out[0][0], 1.0)
  assert_eq(out[1][2], 6.0)
}

test "mc_impute_half_min" {
  let nan = 0.0 / 0.0
  // Row 0: non-missing values [2.0, 4.0], min = 2.0, half = 1.0
  let matrix = [[2.0, nan, 4.0]]
  let out = @src.mc_impute(matrix, "halfMin", 0)
  assert_true((out[0][1] - 1.0).abs() < 0.001)
}

test "mc_impute_mean" {
  let nan = 0.0 / 0.0
  // Row 0: non-missing [1.0, 4.0], mean = 2.5
  let matrix = [[1.0, nan, 4.0]]
  let out = @src.mc_impute(matrix, "mean", 0)
  assert_true((out[0][1] - 2.5).abs() < 0.001)
}

test "mc_impute_median" {
  let nan = 0.0 / 0.0
  // Row 0: non-missing [1.0, 4.0], median = 2.5
  let matrix = [[1.0, nan, 4.0]]
  let out = @src.mc_impute(matrix, "median", 0)
  assert_true((out[0][1] - 2.5).abs() < 0.001)
}

test "mc_impute_knn" {
  let nan = 0.0 / 0.0
  // Row 1 missing at col 1; nearest row (by cols 0 and 2) is row 0.
  // k=1 -> impute from row 0 col 1 = 2.0
  let matrix = [
    [1.0, 2.0, 3.0],
    [1.0, nan, 3.0],
    [5.0, 4.0, 3.0],
  ]
  let out = @src.mc_impute(matrix, "knn", 1)
  assert_true((out[1][1] - 2.0).abs() < 0.001)
}

test "mc_impute_knn_k2" {
  let nan = 0.0 / 0.0
  // k=2 -> mean of row 0 col 1 (2.0) and row 2 col 1 (4.0) = 3.0
  let matrix = [
    [1.0, 2.0, 3.0],
    [1.0, nan, 3.0],
    [5.0, 4.0, 3.0],
  ]
  let out = @src.mc_impute(matrix, "knn", 2)
  assert_true((out[1][1] - 3.0).abs() < 0.001)
}

// ============================================================================
// medianPolish: mc_median_polish
// ============================================================================

test "mc_median_polish_basic" {
  // Additive matrix: x[i][j] = overall + row_eff[i] + col_eff[j]
  // Expected: overall=4, row_eff=[-2, 2], col_eff=[-1, 1], residuals all 0
  let matrix = [
    [1.0, 3.0],
    [5.0, 7.0],
  ]
  let (r, row_eff, col_eff, overall) = @src.mc_median_polish(matrix, 10, 0.0001)
  assert_true((overall - 4.0).abs() < 0.001)
  assert_true((row_eff[0] - (-2.0)).abs() < 0.001)
  assert_true((row_eff[1] - 2.0).abs() < 0.001)
  assert_true((col_eff[0] - (-1.0)).abs() < 0.001)
  assert_true((col_eff[1] - 1.0).abs() < 0.001)
  // Residuals should be ~0
  for i in 0..<2 {
    for j in 0..<2 {
      assert_true(r[i][j].abs() < 0.001)
    }
  }
}

test "mc_median_polish_empty" {
  let empty : Array[Array[Double]] = []
  let (r, row_eff, col_eff, overall) = @src.mc_median_polish(empty, 10, 0.0001)
  assert_eq(r.length(), 0)
  assert_eq(row_eff.length(), 0)
  assert_eq(col_eff.length(), 0)
  assert_eq(overall, 0.0)
}

// ============================================================================
// robustSummary: mc_robust_summary
// ============================================================================

test "mc_robust_summary_basic" {
  // Per-column median:
  // col 0: median([1,3,5]) = 3.0
  // col 1: median([10,20,30]) = 20.0
  let matrix = [
    [1.0, 10.0],
    [3.0, 20.0],
    [5.0, 30.0],
  ]
  let out = @src.mc_robust_summary(matrix)
  assert_eq(out.length(), 2)
  assert_eq(out[0], 3.0)
  assert_eq(out[1], 20.0)
}

test "mc_robust_summary_empty" {
  let empty : Array[Array[Double]] = []
  let out = @src.mc_robust_summary(empty)
  assert_eq(out.length(), 0)
}

// ============================================================================
// MAD: mc_mad
// ============================================================================

test "mc_mad_basic" {
  // median = 3, deviations = [2,1,0,1,2], median of devs = 1.0
  // mad = 1.0 * 1.4826 = 1.4826
  let m = @src.mc_mad([1.0, 2.0, 3.0, 4.0, 5.0])
  assert_true((m - 1.4826).abs() < 0.001)
}

test "mc_mad_constant" {
  let m = @src.mc_mad([7.0, 7.0, 7.0])
  assert_eq(m, 0.0)
}

// ============================================================================
// Validity: mc_valid_peak_list, mc_has_missing
// ============================================================================

test "mc_valid_peak_list_valid" {
  let mz = [100.0, 200.0, 300.0]
  let intensity = [1.0, 2.0, 3.0]
  assert_true(@src.mc_valid_peak_list(mz, intensity))
}

test "mc_valid_peak_list_not_increasing" {
  // m/z 200 == 200 is not strictly increasing
  let mz = [100.0, 200.0, 200.0]
  let intensity = [1.0, 2.0, 3.0]
  assert_true(!@src.mc_valid_peak_list(mz, intensity))
}

test "mc_valid_peak_list_length_mismatch" {
  let mz = [100.0, 200.0]
  let intensity = [1.0, 2.0, 3.0]
  assert_true(!@src.mc_valid_peak_list(mz, intensity))
}

test "mc_valid_peak_list_negative_intensity" {
  let mz = [100.0, 200.0]
  let intensity = [-1.0, 2.0]
  assert_true(!@src.mc_valid_peak_list(mz, intensity))
}

test "mc_has_missing_clean" {
  let matrix = [
    [1.0, 2.0],
    [3.0, 4.0],
  ]
  assert_true(!@src.mc_has_missing(matrix))
}

test "mc_has_missing_with_nan" {
  let nan = 0.0 / 0.0
  let matrix = [
    [1.0, nan],
    [3.0, 4.0],
  ]
  assert_true(@src.mc_has_missing(matrix))
}

// ============================================================================
// Aggregation: mc_aggregate_rows
// ============================================================================

test "mc_aggregate_rows_sum" {
  let matrix = [
    [1.0, 10.0],
    [2.0, 20.0],
    [3.0, 30.0],
  ]
  let groups = [0, 0, 1]
  let out = @src.mc_aggregate_rows(matrix, groups, "sum")
  assert_eq(out.length(), 2)
  assert_true((out[0][0] - 3.0).abs() < 0.001)
  assert_true((out[0][1] - 30.0).abs() < 0.001)
  assert_true((out[1][0] - 3.0).abs() < 0.001)
  assert_true((out[1][1] - 30.0).abs() < 0.001)
}

test "mc_aggregate_rows_mean" {
  let matrix = [
    [1.0, 10.0],
    [2.0, 20.0],
    [3.0, 30.0],
  ]
  let groups = [0, 0, 1]
  let out = @src.mc_aggregate_rows(matrix, groups, "mean")
  assert_eq(out.length(), 2)
  assert_true((out[0][0] - 1.5).abs() < 0.001)
  assert_true((out[0][1] - 15.0).abs() < 0.001)
  assert_true((out[1][0] - 3.0).abs() < 0.001)
  assert_true((out[1][1] - 30.0).abs() < 0.001)
}

test "mc_aggregate_rows_median" {
  let matrix = [
    [1.0, 10.0],
    [2.0, 20.0],
    [3.0, 30.0],
  ]
  let groups = [0, 0, 1]
  let out = @src.mc_aggregate_rows(matrix, groups, "median")
  assert_eq(out.length(), 2)
  assert_true((out[0][0] - 1.5).abs() < 0.001)
  assert_true((out[0][1] - 15.0).abs() < 0.001)
  assert_true((out[1][0] - 3.0).abs() < 0.001)
  assert_true((out[1][1] - 30.0).abs() < 0.001)
}

// ============================================================================
// Normalization: mc_normalize_tic
// ============================================================================

test "mc_normalize_tic_basic" {
  // total = 10, out = v / total * 100
  let out = @src.mc_normalize_tic([1.0, 2.0, 3.0, 4.0])
  assert_eq(out.length(), 4)
  assert_true((out[0] - 10.0).abs() < 0.001)
  assert_true((out[1] - 20.0).abs() < 0.001)
  assert_true((out[2] - 30.0).abs() < 0.001)
  assert_true((out[3] - 40.0).abs() < 0.001)
}

test "mc_normalize_tic_zero_total" {
  // total = 0 -> all zeros
  let out = @src.mc_normalize_tic([0.0, 0.0, 0.0])
  assert_eq(out.length(), 3)
  assert_eq(out[0], 0.0)
  assert_eq(out[1], 0.0)
  assert_eq(out[2], 0.0)
}

test "mc_normalize_tic_sums_to_100" {
  let out = @src.mc_normalize_tic([10.0, 20.0, 30.0, 40.0])
  let total = @src.mc_sum(out)
  assert_true((total - 100.0).abs() < 0.001)
}
