///|
/// MsCoreUtils mass-spectrometry signal processing demo.
///
/// Demonstrates centroiding, peak detection, smoothing, baseline estimation,
/// calibration, imputation, and aggregation on synthetic MS data.

fn main() -> Unit {
  println("=== MsCoreUtils Signal Processing Demo ===\n")

  // 1. Peak detection and centroiding
  println("--- 1. Peak Detection & Centroiding ---\n")

  // Simulate profile-mode MS data with 3 peaks
  let mz : Array[Double] = Array::new()
  let intensity : Array[Double] = Array::new()
  for i in 0..<300 {
    let x = 100.0 + i.to_double() * 0.5
    mz.push(x)
    // Peak 1: centered at 120 (intensity 1000)
    let p1 = 1000.0 * @math.exp(-((x - 120.0) * (x - 120.0)) / 2.0)
    // Peak 2: centered at 180 (intensity 500)
    let p2 = 500.0 * @math.exp(-((x - 180.0) * (x - 180.0)) / 2.0)
    // Peak 3: centered at 220 (intensity 800)
    let p3 = 800.0 * @math.exp(-((x - 220.0) * (x - 220.0)) / 2.0)
    intensity.push(p1 + p2 + p3 + 10.0) // +10 = noise floor
  }

  let peaks = @bio.mc_local_maxima(intensity)
  println("Detected " + peaks.length().to_string() + " local maxima at m/z indices:")
  for p in peaks {
    println("  index " + p.to_string() + " -> m/z=" + mz[p].to_string() + ", intensity=" + intensity[p].to_string())
  }

  // Refine centroids
  let (centroid_mz, centroid_int) = @bio.mc_refine_centroids(mz, intensity, peaks, 3)
  println("\nRefined centroids:")
  for i in 0..<centroid_mz.length() {
    println("  m/z=" + centroid_mz[i].to_string() + ", intensity=" + centroid_int[i].to_string())
  }

  // 2. Smoothing
  println("\n--- 2. Signal Smoothing ---\n")

  // Create noisy signal
  let noisy : Array[Double] = [10.0, 15.0, 8.0, 12.0, 20.0, 18.0, 25.0, 22.0, 30.0, 28.0]
  let smoothed_ma = @bio.mc_smooth_moving_average(noisy, 2)
  let smoothed_sg = @bio.mc_smooth_savitzky_golay(noisy, 2)

  println("Original:     " + arr_to_string(noisy))
  println("Moving avg:   " + arr_to_string(smoothed_ma))
  println("Savitzky-Gol: " + arr_to_string(smoothed_sg))

  // 3. Baseline estimation (SNIP)
  println("\n--- 3. Baseline Estimation (SNIP) ---\n")

  let signal : Array[Double] = [10.0, 12.0, 100.0, 500.0, 100.0, 15.0, 10.0, 200.0, 800.0, 200.0, 12.0, 10.0]
  let baseline = @bio.mc_baseline_snip(signal, 4)
  println("Signal:    " + arr_to_string(signal))
  println("Baseline:  " + arr_to_string(baseline))
  let corrected = Array::new()
  for i in 0..<signal.length() {
    let c = signal[i] - baseline[i]
    corrected.push(if c > 0.0 { c } else { 0.0 })
  }
  println("Corrected: " + arr_to_string(corrected))

  // 4. Calibration
  println("\n--- 4. m/z Calibration ---\n")

  let observed : Array[Double] = [120.05, 180.10, 220.15, 300.20]
  let calibrants : Array[Double] = [120.00, 180.00, 220.00, 300.00]
  let calibrated_shift = @bio.mc_calibrate(observed, calibrants, 1.0, "shift")
  let calibrated_linear = @bio.mc_calibrate(observed, calibrants, 1.0, "linear")
  println("Observed:           " + arr_to_string(observed))
  println("Calibrants:         " + arr_to_string(calibrants))
  println("Calibrated (shift): " + arr_to_string(calibrated_shift))
  println("Calibrated (linear):" + arr_to_string(calibrated_linear))

  // 5. Peak matching
  println("\n--- 5. Peak Matching ---\n")

  let peaks_a : Array[Double] = [100.0, 200.0, 300.05, 400.1]
  let peaks_b : Array[Double] = [100.01, 200.02, 300.0, 400.5]
  let matched_tol = @bio.mc_join_peaks(peaks_a, peaks_b, 0.1, false)
  let matched_ppm = @bio.mc_join_peaks(peaks_a, peaks_b, 500.0, true)
  println("Peak list A: " + arr_to_string(peaks_a))
  println("Peak list B: " + arr_to_string(peaks_b))
  println("Matched (tol=0.1 Da): " + matched_tol.length().to_string() + " pairs")
  for i in 0..<matched_tol.length() {
    let (a, b) = matched_tol[i]
    println("  A[" + a.to_string() + "]=" + peaks_a[a].to_string() + " <-> B[" + b.to_string() + "]=" + peaks_b[b].to_string())
  }
  println("Matched (tol=500 ppm): " + matched_ppm.length().to_string() + " pairs")

  // 6. Imputation
  println("\n--- 6. Missing Value Imputation ---\n")

  let nan_val = 0.0 / 0.0
  let matrix : Array[Array[Double]] = [
    [1.0, 2.0, nan_val],
    [4.0, nan_val, 6.0],
    [7.0, 8.0, 9.0],
  ]
  println("Original matrix (NaN = missing):")
  print_matrix(matrix)

  let imputed_mean = @bio.mc_impute(matrix, "mean", 3)
  println("\nImputed (mean):")
  print_matrix(imputed_mean)

  let imputed_median = @bio.mc_impute(matrix, "median", 3)
  println("\nImputed (median):")
  print_matrix(imputed_median)

  let imputed_knn = @bio.mc_impute(matrix, "knn", 2)
  println("\nImputed (knn, k=2):")
  print_matrix(imputed_knn)

  // 7. Median polish aggregation
  println("\n--- 7. Median Polish Aggregation ---\n")

  let expr_matrix : Array[Array[Double]] = [
    [10.0, 12.0, 11.0],
    [20.0, 22.0, 21.0],
    [30.0, 32.0, 31.0],
  ]
  let (polished, row_eff, col_eff, overall) = @bio.mc_median_polish(expr_matrix, 10, 0.001)
  println("Input matrix:")
  print_matrix(expr_matrix)
  println("Overall effect: " + overall.to_string())
  println("Row effects: " + arr_to_string(row_eff))
  println("Col effects: " + arr_to_string(col_eff))
  println("Residuals:")
  print_matrix(polished)

  // 8. TIC normalization
  println("\n--- 8. TIC Normalization ---\n")

  let raw_intensities : Array[Double] = [100.0, 200.0, 300.0, 400.0]
  let normalized = @bio.mc_normalize_tic(raw_intensities)
  println("Raw:        " + arr_to_string(raw_intensities))
  println("Normalized: " + arr_to_string(normalized))

  println("\n=== Demo Complete ===")
}

fn arr_to_string(arr : Array[Double]) -> String {
  let mut s = "["
  for i in 0..<arr.length() {
    if i > 0 {
      s = s + ", "
    }
    s = s + arr[i].to_string()
  }
  s + "]"
}

fn print_matrix(m : Array[Array[Double]]) -> Unit {
  for row in m {
    println("  " + arr_to_string(row))
  }
}
