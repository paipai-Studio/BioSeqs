///|
/// Demo for xcms_peaks module: xcms-style chromatographic peak detection.

fn main {
  println("=== xcms_peaks Demo: Chromatographic Peak Detection ===")
  println("")

  // Example 1: Create a chromatogram with Gaussian peaks
  println("--- Step 1: Create Chromatogram ---")
  let rt : Array[Double] = Array::make(100, 0.0)
  let intensity : Array[Double] = Array::make(100, 0.0)
  for i in 0..<100 {
    rt[i] = i.to_double()
    // Two Gaussian peaks at 25 and 70
    let x1 = (i - 25).to_double()
    let x2 = (i - 70).to_double()
    let peak1 = 100.0 * @math.exp(-0.5 * x1 * x1 / 9.0)
    let peak2 = 80.0 * @math.exp(-0.5 * x2 * x2 / 4.0)
    intensity[i] = peak1 + peak2 + 1.0  // + baseline
  }
  let chrom = @src.xcp_chrom_new(rt, intensity, label="test")
  println("Chromatogram: " + rt.length().to_string() + " points")
  println("")

  // Example 2: Noise and SNR
  println("--- Step 2: Noise Estimation ---")
  let noise = @src.xcp_noise_mad(intensity)
  let snr = @src.xcp_snr(intensity)
  println("Noise (MAD): " + noise.to_string())
  println("SNR: " + snr.to_string())
  println("")

  // Example 3: centWave peak detection
  println("--- Step 3: centWave Peak Detection ---")
  let peaks = @src.xcp_centwave(chrom, snr_threshold=3.0)
  println("Peaks found: " + peaks.length().to_string())
  for p in peaks {
    println("  rt=" + p.rt.to_string() + " intensity=" + p.intensity.to_string() +
            " snr=" + p.snr.to_string() + " fwhm=" + p.fwhm.to_string())
  }
  println("")

  // Example 4: matchedFilter peak detection
  println("--- Step 4: matchedFilter Peak Detection ---")
  let peaks_mf = @src.xcp_matched_filter(chrom, window_size=5, snr_threshold=3.0)
  println("Peaks found: " + peaks_mf.length().to_string())
  for p in peaks_mf {
    println("  rt=" + p.rt.to_string() + " intensity=" + p.intensity.to_string() +
            " snr=" + p.snr.to_string())
  }
  println("")

  // Example 5: ASCII plot
  println("--- Step 5: ASCII Plot ---")
  let plot = @src.xcp_ascii_plot(chrom, peaks, width=50, height=10)
  println(plot)

  // Example 6: Export peaks
  println("--- Step 6: Export Peaks ---")
  let tsv = @src.xcp_export_peaks(peaks)
  println(tsv)
  println("")

  // Example 7: TIC statistics
  println("--- Step 7: TIC Statistics ---")
  let (total, mx, mean) = @src.xcp_tic_stats(intensity)
  println("Total intensity: " + total.to_string())
  println("Max intensity: " + mx.to_string())
  println("Mean intensity: " + mean.to_string())
}
