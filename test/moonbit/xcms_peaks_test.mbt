///|
/// Tests for xcms_peaks module (xcms-style chromatographic peak detection).

test "xcp_chrom_new" {
  let rt = [1.0, 2.0, 3.0, 4.0, 5.0]
  let intensity = [10.0, 20.0, 50.0, 20.0, 10.0]
  let chrom = @src.xcp_chrom_new(rt, intensity, label="test")
  assert_eq(chrom.rt.length(), 5)
  assert_eq(chrom.intensity[2], 50.0)
  assert_eq(chrom.label, "test")
}

test "xcp_noise_mad" {
  let intensity = [1.0, 2.0, 3.0, 4.0, 5.0, 100.0]  // outlier
  let noise = @src.xcp_noise_mad(intensity)
  assert_true(noise > 0.0)
  assert_true(noise < 50.0)  // MAD is robust to outlier
}

test "xcp_snr" {
  let intensity = [1.0, 2.0, 3.0, 100.0, 3.0, 2.0, 1.0]
  let snr = @src.xcp_snr(intensity)
  assert_true(snr > 1.0)
}

test "xcp_baseline_rolling_min" {
  let intensity = [10.0, 5.0, 20.0, 15.0, 30.0, 25.0]
  let base = @src.xcp_baseline_rolling_min(intensity, 3)
  assert_eq(base.length(), 6)
  assert_true(base[2] <= 10.0)
}

test "xcp_centwave_single_peak" {
  // Simulate a single Gaussian peak
  let rt : Array[Double] = Array::make(50, 0.0)
  let intensity : Array[Double] = Array::make(50, 0.0)
  for i in 0..<50 {
    rt[i] = i.to_double()
    // Gaussian centered at 25, sigma=3
    let x = (i - 25).to_double()
    intensity[i] = 100.0 * @math.exp(-0.5 * x * x / 9.0) + 1.0
  }
  let chrom = @src.xcp_chrom_new(rt, intensity)
  let peaks = @src.xcp_centwave(chrom, snr_threshold=2.0)
  assert_true(peaks.length() >= 1)
  assert_true(peaks[0].rt > 20.0 && peaks[0].rt < 30.0)
  assert_true(peaks[0].snr > 1.0)
}

test "xcp_centwave_no_peak" {
  let rt = [1.0, 2.0, 3.0, 4.0, 5.0]
  let intensity = [1.0, 1.0, 1.0, 1.0, 1.0]  // flat
  let chrom = @src.xcp_chrom_new(rt, intensity)
  let peaks = @src.xcp_centwave(chrom, snr_threshold=3.0)
  assert_eq(peaks.length(), 0)
}

test "xcp_matched_filter_single_peak" {
  let rt : Array[Double] = Array::make(50, 0.0)
  let intensity : Array[Double] = Array::make(50, 0.0)
  for i in 0..<50 {
    rt[i] = i.to_double()
    let x = (i - 25).to_double()
    intensity[i] = 50.0 * @math.exp(-0.5 * x * x / 4.0) + 1.0
  }
  let chrom = @src.xcp_chrom_new(rt, intensity)
  let peaks = @src.xcp_matched_filter(chrom, window_size=5, snr_threshold=2.0)
  assert_true(peaks.length() >= 1)
  assert_true(peaks[0].rt > 20.0 && peaks[0].rt < 30.0)
}

test "xcp_merge_peaks" {
  let peaks = [
    @src.XcpPeak::{mz: 0.0, rt: 25.0, rtmin: 20.0, rtmax: 30.0, intensity: 100.0, area: 500.0, snr: 10.0, fwhm: 5.0, peak_method: "centWave"},
    @src.XcpPeak::{mz: 0.0, rt: 26.0, rtmin: 22.0, rtmax: 32.0, intensity: 80.0, area: 400.0, snr: 8.0, fwhm: 5.0, peak_method: "matchedFilter"},
  ]
  let merged = @src.xcp_merge_peaks(peaks)
  assert_true(merged.length() >= 1)
}

test "xcp_export_peaks" {
  let peaks = [
    @src.XcpPeak::{mz: 100.5, rt: 25.0, rtmin: 20.0, rtmax: 30.0, intensity: 100.0, area: 500.0, snr: 10.0, fwhm: 5.0, peak_method: "centWave"},
  ]
  let tsv = @src.xcp_export_peaks(peaks)
  assert_true(tsv.contains("mz"))
  assert_true(tsv.contains("100.5"))
  assert_true(tsv.contains("centWave"))
}

test "xcp_tic_stats" {
  let intensity = [10.0, 20.0, 30.0, 40.0, 50.0]
  let (total, mx, mean) = @src.xcp_tic_stats(intensity)
  assert_eq(total, 150.0)
  assert_eq(mx, 50.0)
  assert_eq(mean, 30.0)
}

test "xcp_ascii_plot" {
  let rt = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let intensity = [1.0, 2.0, 5.0, 10.0, 5.0, 2.0, 1.0, 0.5]
  let chrom = @src.xcp_chrom_new(rt, intensity)
  let peaks = [@src.XcpPeak::{mz: 0.0, rt: 4.0, rtmin: 3.0, rtmax: 5.0, intensity: 10.0, area: 30.0, snr: 5.0, fwhm: 2.0, peak_method: "test"}]
  let plot = @src.xcp_ascii_plot(chrom, peaks, width=20, height=10)
  assert_true(plot.contains("+"))
  assert_true(plot.contains("|"))
}

test "xcp_gaussian" {
  let val = @src.xcp_gaussian(0.0, 0.0, 1.0)
  assert_true(val > 0.39 && val < 0.40)  // 1/sqrt(2pi) = 0.3989
}
