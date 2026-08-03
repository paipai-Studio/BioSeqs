///|
/// Tests for the PCD mass spectrometry format parser (pcd.mbt).

// ============================================================================
// PcdSpectrum Tests
// ============================================================================

test "pcd_spectrum_construction" {
  let spec = @src.PcdSpectrum::new(1, rt=12.5)
  assert_eq(spec.scan_number, 1)
  assert_true(spec.rt > 0.0)
  assert_eq(spec.num_peaks(), 0)
}

test "pcd_spectrum_add_peak" {
  let spec = @src.PcdSpectrum::new(1)
  spec.add_peak(500.5, 10000.0)
  spec.add_peak(600.5, 20000.0)
  assert_eq(spec.num_peaks(), 2)
  assert_true(spec.total_intensity() > 0.0)
}

test "pcd_spectrum_precursor" {
  let spec = @src.PcdSpectrum::new(1)
  spec.set_precursor_mz(Some(500.5))
  spec.set_precursor_charge(Some(2))
  match spec.precursor_mz {
    Some(mz) => assert_true(mz > 500.0 && mz < 501.0)
    None => assert_true(false)
  }
  match spec.precursor_charge {
    Some(z) => assert_eq(z, 2)
    None => assert_true(false)
  }
}

// ============================================================================
// PcdFile Tests
// ============================================================================

test "pcd_file_construction" {
  let file = @src.PcdFile::new()
  assert_eq(file.num_spectra, 0)
  assert_eq(file.num_peaks, 0)
}

test "pcd_file_add_spectrum" {
  let file = @src.PcdFile::new()
  let spec = @src.PcdSpectrum::new(1, rt=10.0)
  spec.add_peak(500.0, 1000.0)
  file.add_spectrum(spec)
  file.update_counts()
  assert_eq(file.num_spectra, 1)
  assert_eq(file.num_peaks, 1)
}

test "pcd_file_get_spectrum" {
  let file = @src.PcdFile::new()
  let spec1 = @src.PcdSpectrum::new(1)
  let spec2 = @src.PcdSpectrum::new(2)
  file.add_spectrum(spec1)
  file.add_spectrum(spec2)
  match file.get_spectrum(1) {
    Some(s) => assert_eq(s.scan_number, 1)
    None => assert_true(false)
  }
  match file.get_spectrum(3) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ============================================================================
// PCD Parsing Tests
// ============================================================================

test "pcd_parse_empty" {
  let file = @src.pcd_parse("")
  assert_eq(file.num_spectra, 0)
  assert_eq(file.num_peaks, 0)
}

test "pcd_parse_sample" {
  let content = @src.pcd_sample()
  let file = @src.pcd_parse(content)
  assert_true(file.num_spectra >= 1)
  assert_true(file.num_peaks >= 1)
}

test "pcd_parse_with_precursor" {
  let content = @src.pcd_sample()
  let file = @src.pcd_parse(content)
  assert_true(file.num_spectra >= 1)
}

// ============================================================================
// PCD Query Tests
// ============================================================================

test "pcd_total_peaks" {
  let file = @src.PcdFile::new()
  let spec1 = @src.PcdSpectrum::new(1)
  spec1.add_peak(500.0, 1000.0)
  spec1.add_peak(600.0, 2000.0)
  let spec2 = @src.PcdSpectrum::new(2)
  spec2.add_peak(700.0, 3000.0)
  file.add_spectrum(spec1)
  file.add_spectrum(spec2)
  file.update_counts()
  assert_eq(@src.pcd_total_peaks(file), 3)
}

test "pcd_bpc" {
  let file = @src.PcdFile::new()
  let spec1 = @src.PcdSpectrum::new(1, rt=10.0)
  spec1.add_peak(500.0, 1000.0)
  spec1.add_peak(600.0, 5000.0)
  let spec2 = @src.PcdSpectrum::new(2, rt=20.0)
  spec2.add_peak(700.0, 3000.0)
  spec2.add_peak(800.0, 4000.0)
  file.add_spectrum(spec1)
  file.add_spectrum(spec2)
  file.update_counts()
  let bpc = @src.pcd_bpc(file)
  assert_eq(bpc.length(), 2)
}

test "pcd_tic" {
  let file = @src.PcdFile::new()
  let spec1 = @src.PcdSpectrum::new(1, rt=10.0)
  spec1.add_peak(500.0, 1000.0)
  spec1.add_peak(600.0, 2000.0)
  let spec2 = @src.PcdSpectrum::new(2, rt=20.0)
  spec2.add_peak(700.0, 3000.0)
  file.add_spectrum(spec1)
  file.add_spectrum(spec2)
  file.update_counts()
  let tic = @src.pcd_tic(file)
  assert_eq(tic.length(), 2)
}

test "pcd_spectrum_by_mz" {
  let spec = @src.PcdSpectrum::new(1)
  spec.add_peak(500.0, 1000.0)
  spec.add_peak(600.0, 2000.0)
  spec.add_peak(700.0, 3000.0)
  let filtered = @src.pcd_spectrum_by_mz(spec, 550.0, 650.0)
  assert_eq(filtered.num_peaks(), 1)
}

// ============================================================================
// PCD Serialization Tests
// ============================================================================

test "pcd_write_roundtrip" {
  let content = @src.pcd_sample()
  let file = @src.pcd_parse(content)
  let written = @src.pcd_write(file)
  let file2 = @src.pcd_parse(written)
  assert_eq(file2.num_spectra, file.num_spectra)
}

// ============================================================================
// PCD Summary Tests
// ============================================================================

test "pcd_summary" {
  let file = @src.PcdFile::new()
  let spec = @src.PcdSpectrum::new(1)
  spec.add_peak(500.0, 1000.0)
  file.add_spectrum(spec)
  file.update_counts()
  let summary = @src.pcd_summary(file)
  assert_true(summary.contains("PCD File Summary"))
}
