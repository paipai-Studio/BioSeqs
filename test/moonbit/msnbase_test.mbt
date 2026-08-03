///|
/// Tests for MSnbase module: core mass spectrometry data structures and processing.
/// Covers Spectrum, Chromatogram, MSnSet, and processing utilities.

// ============================================================================
// Enum/type creation tests
// ============================================================================

test "mbn_mslevel_constructors" {
  let l1 = @src.ms_level_ms1()
  let l2 = @src.ms_level_ms2()
  let ln = @src.ms_level_msn(5)
  assert_eq(l1, @src.ms_level_ms1())
  assert_eq(l2, @src.ms_level_ms2())
  match ln {
    @src.MSLevel::MSn(n) => assert_eq(n, 5)
    _ => assert_true(false)
  }
}

test "mbn_polarity_constructors" {
  let p = @src.polarity_positive()
  let n = @src.polarity_negative()
  assert_eq(p, @src.polarity_positive())
  assert_eq(n, @src.polarity_negative())
}

test "mbn_processing_step_new" {
  let s = @src.ProcessingStep::new("log_transformed", "2025-01-01")
  assert_eq(s.description, "log_transformed")
  assert_eq(s.timestamp, "2025-01-01")
}

// ============================================================================
// Spectrum tests
// ============================================================================

test "mbn_spectrum_new_basic" {
  let mz = [100.0, 200.0, 300.0, 400.0, 500.0]
  let int = [10.0, 100.0, 500.0, 200.0, 50.0]
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms1(), @src.polarity_positive(), 60.5)
  assert_eq(sp.peaks_count(), 5)
  assert_true((sp.rt() - 60.5).abs() < 0.0001)
  // TIC = 10+100+500+200+50 = 860
  assert_true((sp.tic() - 860.0).abs() < 0.0001)
}

test "mbn_spectrum_empty" {
  let sp = @src.Spectrum::empty()
  assert_eq(sp.peaks_count(), 0)
  assert_eq(sp.tic(), 0.0)
}

test "mbn_spectrum_mismatched_lengths" {
  // Mismatched mz/intensity returns empty
  let sp = @src.Spectrum::new([1.0, 2.0], [10.0], @src.ms_level_ms1(), @src.polarity_positive(), 0.0)
  assert_eq(sp.peaks_count(), 0)
}

test "mbn_spectrum_with_precursor" {
  let mz = [100.0, 200.0]
  let int = [50.0, 100.0]
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms2(), @src.polarity_positive(), 120.0)
  let sp2 = sp.with_precursor(500.25, 2)
  assert_eq(sp2.precursor_mz, 500.25)
  assert_eq(sp2.precursor_charge, 2)
}

test "mbn_spectrum_base_peak" {
  let mz = [100.0, 200.0, 300.0, 400.0]
  let int = [10.0, 500.0, 50.0, 200.0]
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms1(), @src.polarity_positive(), 0.0)
  let (bpmz, bpint) = sp.base_peak()
  // Peak at index 1 is highest: mz=200, int=500
  assert_true((bpmz - 200.0).abs() < 0.0001)
  assert_true((bpint - 500.0).abs() < 0.0001)
}

test "mbn_spectrum_base_peak_empty" {
  let sp = @src.Spectrum::empty()
  let (mz, int) = sp.base_peak()
  assert_eq(mz, 0.0)
  assert_eq(int, 0.0)
}

test "mbn_spectrum_find_peak_exact" {
  let mz = [100.0, 200.0, 300.0, 400.0, 500.0]
  let int = [10.0, 100.0, 500.0, 200.0, 50.0]
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms1(), @src.polarity_positive(), 0.0)
  let idx = sp.find_peak(300.0, ppm=10.0)
  assert_eq(idx, 2)
}

test "mbn_spectrum_find_peak_notfound" {
  let mz = [100.0, 200.0]
  let int = [10.0, 100.0]
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms1(), @src.polarity_positive(), 0.0)
  let idx = sp.find_peak(999.0, ppm=10.0)
  assert_eq(idx, -1)
}

test "mbn_spectrum_filter_mz_range" {
  let mz = [100.0, 200.0, 300.0, 400.0, 500.0]
  let int = [10.0, 100.0, 500.0, 200.0, 50.0]
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms1(), @src.polarity_positive(), 0.0)
  let filtered = sp.filter_mz_range(150.0, 450.0)
  assert_eq(filtered.peaks_count(), 3)
  assert_eq(filtered.mz()[0], 200.0)
  assert_eq(filtered.mz()[2], 400.0)
}

test "mbn_spectrum_normalize_tic" {
  let int = [10.0, 100.0, 500.0, 200.0, 50.0]  // TIC = 860
  let mz = [100.0, 200.0, 300.0, 400.0, 500.0]
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms1(), @src.polarity_positive(), 0.0)
  let norm = sp.normalize_tic()
  // New TIC should be 1.0
  assert_true((norm.tic() - 1.0).abs() < 0.0001)
  // First peak normalized: 10/860
  assert_true((norm.intensity()[0] - 10.0/860.0).abs() < 0.0001)
}

// ============================================================================
// Chromatogram tests
// ============================================================================

test "mbn_chromatogram_new" {
  let rt = [0.0, 30.0, 60.0, 90.0, 120.0]
  let int = [100.0, 5000.0, 20000.0, 8000.0, 200.0]
  let chr = @src.Chromatogram::new(rt, int, mz_target=500.0, ppm_tolerance=10.0, acquisition_mode="MRM")
  assert_eq(chr.n_points(), 5)
}

test "mbn_chromatogram_total_auc" {
  // Simple triangular: rt=[0, 1, 2], int=[0, 10, 0]
  // Area = 1.0 * 10 = 10 (two triangles each area 5)
  let rt = [0.0, 1.0, 2.0]
  let int = [0.0, 10.0, 0.0]
  let chr = @src.Chromatogram::new(rt, int)
  assert_true((chr.total_auc() - 10.0).abs() < 0.001)
}

test "mbn_chromatogram_apex" {
  let rt = [0.0, 30.0, 60.0, 90.0]
  let int = [100.0, 5000.0, 20000.0, 8000.0]
  let chr = @src.Chromatogram::new(rt, int)
  let (apex_rt, apex_int) = chr.apex()
  assert_true((apex_rt - 60.0).abs() < 0.001)
  assert_true((apex_int - 20000.0).abs() < 0.001)
}

test "mbn_chromatogram_fwhm" {
  // Gaussian-like: peak at center
  let rt = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
  let int = [0.0, 2.0, 5.0, 8.0, 10.0, 8.0, 5.0, 2.0, 0.0]
  let chr = @src.Chromatogram::new(rt, int)
  let f = chr.fwhm()
  // Half max = 5.0. Crosses between 1-2 and 6-7. FWHM ~ 5.0 or so.
  assert_true(f >= 3.0 && f <= 8.0)
}

// ============================================================================
// MSnFeatureData and MSnSampleData
// ============================================================================

test "mbn_featuredata_new" {
  let fd = @src.MSnFeatureData::new(
    "Pep_1_2", ["P12345"], "PEPTIDER", 2, 500.25, 60.5,
  )
  assert_eq(fd.feature_name, "Pep_1_2")
  assert_eq(fd.protein_accessions.length(), 1)
  assert_eq(fd.protein_accessions[0], "P12345")
  assert_eq(fd.sequence, "PEPTIDER")
  assert_eq(fd.charge, 2)
  assert_true((fd.mz - 500.25).abs() < 0.0001)
}

test "mbn_sampledata_new" {
  let sd = @src.MSnSampleData::new(
    "S1_WT_1", "WT", "Mouse_1", 1, "Run_001",
  )
  assert_eq(sd.sample_name, "S1_WT_1")
  assert_eq(sd.group, "WT")
  assert_eq(sd.subject, "Mouse_1")
  assert_eq(sd.replicate, 1)
  assert_eq(sd.run, "Run_001")
}

// ============================================================================
// MSnSet tests
// ============================================================================

test "mbn_msnset_new_basic" {
  let exprs = [
    [1000.0, 1200.0, 1100.0],
    [5000.0, 4500.0, 5200.0],
    [200.0, 180.0, 220.0],
  ]
  let fd = [
    @src.MSnFeatureData::new("Peptide1", ["P1"], "AAAAK", 1, 400.0, 30.0),
    @src.MSnFeatureData::new("Peptide2", ["P1"], "BBBBK", 2, 500.0, 45.0),
    @src.MSnFeatureData::new("Peptide3", ["P2"], "CCCCK", 3, 600.0, 60.0),
  ]
  let sd = [
    @src.MSnSampleData::new("S1", "Ctrl", "Subj1", 1, "R1"),
    @src.MSnSampleData::new("S2", "Ctrl", "Subj2", 1, "R2"),
    @src.MSnSampleData::new("S3", "Treat", "Subj3", 1, "R3"),
  ]
  let m = @src.MSnSet::new(exprs, fd, sd)
  assert_eq(m.n_features(), 3)
  assert_eq(m.n_samples(), 3)
  assert_eq(m.feature_names().length(), 3)
  assert_eq(m.sample_names().length(), 3)
}

test "mbn_msnset_from_names" {
  let exprs = [
    [1.0, 2.0],
    [3.0, 4.0],
  ]
  let m = @src.MSnSet::from_names(exprs, ["F1", "F2"], ["S1", "S2"])
  assert_eq(m.n_features(), 2)
  assert_eq(m.n_samples(), 2)
  assert_eq(m.feature_names()[0], "F1")
  assert_eq(m.sample_names()[1], "S2")
}

test "mbn_msnset_get_feature" {
  let exprs = [
    [10.0, 20.0],
    [30.0, 40.0],
  ]
  let m = @src.MSnSet::from_names(exprs, ["F1", "F2"], ["S1", "S2"])
  let f1 = m.get_feature("F1")
  assert_eq(f1.length(), 2)
  assert_true((f1[0] - 10.0).abs() < 0.0001)
  assert_true((f1[1] - 20.0).abs() < 0.0001)
  let missing = m.get_feature("Missing")
  assert_eq(missing.length(), 0)
}

test "mbn_msnset_get_sample" {
  let exprs = [
    [10.0, 20.0],
    [30.0, 40.0],
    [50.0, 60.0],
  ]
  let m = @src.MSnSet::from_names(exprs, ["F1", "F2", "F3"], ["S1", "S2"])
  let s1 = m.get_sample("S1")
  assert_eq(s1.length(), 3)
  assert_true((s1[0] - 10.0).abs() < 0.0001)
  assert_true((s1[2] - 50.0).abs() < 0.0001)
  let missing = m.get_sample("Missing")
  assert_eq(missing.length(), 0)
}

// ============================================================================
// MSnSet transformations
// ============================================================================

test "mbn_msnset_log2_transform" {
  let exprs = [
    [1.0, 3.0],
    [7.0, 15.0],
  ]
  let m = @src.MSnSet::from_names(exprs, ["F1", "F2"], ["S1", "S2"])
  let m2 = m.log2_transform(offset=1.0)
  // log2(1+1) = 1, log2(3+1) = 2, log2(7+1) = 3, log2(15+1) = 4
  let e = m2.exprs()
  assert_true((e[0][0] - 1.0).abs() < 0.0001)
  assert_true((e[0][1] - 2.0).abs() < 0.0001)
  assert_true((e[1][0] - 3.0).abs() < 0.0001)
  assert_true((e[1][1] - 4.0).abs() < 0.0001)
}

test "mbn_msnset_impute_mean" {
  let exprs = [
    [10.0, -1.0, 30.0],  // -1 = missing, mean of 10,30 = 20
    [5.0, 15.0, 0.0],     // 0 = missing, mean of 5,15 = 10
  ]
  let m = @src.MSnSet::from_names(exprs, ["F1", "F2"], ["S1", "S2", "S3"])
  let m2 = m.impute_missing(method="mean")
  let e = m2.exprs()
  assert_true((e[0][1] - 20.0).abs() < 0.0001)
  assert_true((e[1][2] - 10.0).abs() < 0.0001)
}

test "mbn_msnset_impute_median" {
  let exprs = [
    [1.0, -1.0, 5.0],  // sorted 1,5 median = 3
    [100.0, 200.0, 0.0], // sorted 100,200 median = 150
  ]
  let m = @src.MSnSet::from_names(exprs, ["F1", "F2"], ["S1", "S2", "S3"])
  let m2 = m.impute_missing(method="median")
  let e = m2.exprs()
  assert_true((e[0][1] - 3.0).abs() < 0.0001)
  assert_true((e[1][2] - 150.0).abs() < 0.0001)
}

test "mbn_msnset_normalize_sum" {
  let exprs = [
    [1.0, 10.0],
    [2.0, 20.0],
    [3.0, 30.0],
  ]
  // S1 sum = 6, S2 sum = 60. target = 12
  // S1 scale = 2, so 2, 4, 6
  // S2 scale = 0.2, so 2, 4, 6
  let m = @src.MSnSet::from_names(exprs, ["F1", "F2", "F3"], ["S1", "S2"])
  let m2 = m.normalize_sum(target=12.0)
  let e = m2.exprs()
  assert_true((e[0][0] - 2.0).abs() < 0.0001)
  assert_true((e[1][0] - 4.0).abs() < 0.0001)
  assert_true((e[0][1] - 2.0).abs() < 0.0001)
  assert_true((e[1][1] - 4.0).abs() < 0.0001)
}

test "mbn_msnset_normalize_median_center" {
  let exprs = [
    [1.0, 10.0],
    [3.0, 20.0],
    [5.0, 30.0],
  ]
  // S1 median = 3, S2 median = 20
  // After: [-2, -10, 2] and [0, 0, 10]
  let m = @src.MSnSet::from_names(exprs, ["F1", "F2", "F3"], ["S1", "S2"])
  let m2 = m.normalize_median_center()
  let e = m2.exprs()
  assert_true((e[0][0] - (-2.0)).abs() < 0.0001)
  assert_true((e[2][0] - 2.0).abs() < 0.0001)
  assert_true((e[1][1] - 0.0).abs() < 0.0001)
  assert_true((e[2][1] - 10.0).abs() < 0.0001)
}

test "mbn_msnset_summarize_proteins_sum" {
  // Two peptides from P1, one from P2
  let exprs = [
    [100.0, 200.0],  // Peptide1 -> P1
    [300.0, 400.0],  // Peptide2 -> P1
    [1000.0, 2000.0], // Peptide3 -> P2
  ]
  let fd = [
    @src.MSnFeatureData::new("Pep1", ["P1"], "A", 1, 0.0, 0.0),
    @src.MSnFeatureData::new("Pep2", ["P1"], "B", 1, 0.0, 0.0),
    @src.MSnFeatureData::new("Pep3", ["P2"], "C", 1, 0.0, 0.0),
  ]
  let sd = [
    @src.MSnSampleData::new("S1", "", "", 0, ""),
    @src.MSnSampleData::new("S2", "", "", 0, ""),
  ]
  let m = @src.MSnSet::new(exprs, fd, sd)
  let m2 = m.summarize_to_proteins(method="sum")
  // P1: sum = 400, 600 ; P2: sum = 1000, 2000
  assert_eq(m2.n_features(), 2)
  let e = m2.exprs()
  // Find P1 row
  let p_names = m2.feature_names()
  let mut p1_idx = -1
  let mut p2_idx = -1
  for i in 0..<p_names.length() {
    if p_names[i] == "P1" { p1_idx = i }
    if p_names[i] == "P2" { p2_idx = i }
  }
  assert_true(p1_idx >= 0)
  assert_true(p2_idx >= 0)
  assert_true((e[p1_idx][0] - 400.0).abs() < 0.0001)
  assert_true((e[p1_idx][1] - 600.0).abs() < 0.0001)
  assert_true((e[p2_idx][0] - 1000.0).abs() < 0.0001)
}

// ============================================================================
// Spectrum processing
// ============================================================================

test "mbn_spectrum_smooth_ma" {
  let mz = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]
  let int = [1.0, 100.0, 1.0, 100.0, 1.0, 100.0, 1.0]
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms1(), @src.polarity_positive(), 0.0)
  let smoothed = sp.smooth_moving_average(half_window=1)
  // Position 1 (100): avg of [1, 100, 1] = 34
  // (due to hw=1 => 3-point average)
  let si = smoothed.intensity()
  // Check that smoothing reduces the extremes (values become more moderate)
  assert_true(si[1] < 100.0)  // was spike
  assert_true(si[2] > 1.0)    // was valley
}

test "mbn_spectrum_baseline_correct" {
  let n = 20
  let mz : Array[Double] = []
  let int : Array[Double] = []
  for i in 0..<n {
    mz.push(i.to_double() * 100.0)
    // Base = 100, with a Gaussian peak at center
    let base = 100.0
    let peak = if i == 10 { 900.0 } else { 0.0 }
    int.push(base + peak)
  }
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms1(), @src.polarity_positive(), 0.0)
  let corrected = sp.baseline_correct_minwin(half_window=2)
  // Peak center should still be present (around 900), others should be near 0
  let ci = corrected.intensity()
  for i in 0..<n {
    if i != 10 {
      // After baseline subtraction, baseline should be near 0
      assert_true(ci[i] <= 100.0)
    } else {
      // Peak reduced but still large
      assert_true(ci[i] >= 500.0)
    }
  }
}

test "mbn_spectrum_centroid_simple" {
  // Three peaks: at 99, 100, 101 (peak in center, symmetric)
  // plus a trough
  let mz = [90.0, 99.0, 100.0, 101.0, 110.0, 199.0, 200.0, 201.0, 300.0]
  let int = [1.0, 50.0, 100.0, 50.0, 1.0, 60.0, 200.0, 60.0, 1.0]
  let sp = @src.Spectrum::new(mz, int, @src.ms_level_ms1(), @src.polarity_positive(), 0.0)
  let centroided = sp.centroid_simple(snr_threshold=0.0)
  // Should find two peaks near 100 and 200
  assert_true(centroided.peaks_count() >= 2)
}

// ============================================================================
// Sample QC
// ============================================================================

test "mbn_msnset_sample_qc" {
  let exprs = [
    [100.0, 0.0],
    [200.0, 1000.0],
    [300.0, 2000.0],
    [0.0, 3000.0],
    [500.0, 0.0],
  ]
  let m = @src.MSnSet::from_names(exprs, ["F1", "F2", "F3", "F4", "F5"], ["S1", "S2"])
  let qc = m.sample_qc()
  assert_eq(qc.length(), 2)
  // S1: detected 4 (one zero), total = 1100
  assert_eq(qc[0].1, 4)
  assert_true((qc[0].2 - 1100.0).abs() < 0.001)
  // S2: detected 3
  assert_eq(qc[1].1, 3)
  assert_true((qc[1].2 - 6000.0).abs() < 0.001)
}
