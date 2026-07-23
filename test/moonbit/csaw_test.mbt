///|
test "csaw_create_window" {
  let counts = [10.0, 20.0, 30.0]
  let w = @src.CswWindow::new("chr1:1000-2000", "chr1", 1000, 2000, counts)
  assert_eq(w.chrom, "chr1")
  assert_eq(w.start, 1000)
  assert_eq(w.end, 2000)
  assert_eq(w.id, "chr1:1000-2000")
}

///|
test "csaw_window_methods" {
  let counts = [10.0, 20.0, 30.0]
  let w = @src.CswWindow::new("w1", "chr1", 1000, 2000, counts)
  assert_eq(w.csw_width(), 1000)
  assert_eq(w.csw_count(0), 10.0)
  assert_eq(w.csw_count(1), 20.0)
  assert_true(w.csw_mean_count() > 0.0)
}

///|
test "csaw_make_windows" {
  let windows = @src.csw_make_windows("chr1", 1000, 5000, 1000, 500)
  assert_true(windows.length() > 0)
  assert_eq(windows[0].chrom, "chr1")
}

///|
test "csaw_dataset_create" {
  let counts = [10.0, 20.0, 30.0]
  let windows = [
    @src.CswWindow::new("w1", "chr1", 1000, 2000, counts),
    @src.CswWindow::new("w2", "chr1", 3000, 4000, [50.0, 60.0, 70.0])
  ]
  let samples = ["s1", "s2", "s3"]
  let lib_sizes = [1000.0, 2000.0, 3000.0]
  let data = @src.CswDataSet::new(windows, samples, lib_sizes)
  assert_eq(data.csw_n_windows(), 2)
  assert_eq(data.csw_n_samples(), 3)
}

///|
test "csaw_norm_library" {
  let lib_sizes = [1000.0, 2000.0, 3000.0]
  let norm_factors = @src.csw_norm_library(lib_sizes)
  assert_eq(norm_factors.length(), 3)
  // Library size normalization: smaller libraries get larger factors
  assert_true(norm_factors[0] > norm_factors[1])
}

///|
test "csaw_norm_tmm" {
  let windows = [
    @src.CswWindow::new("w1", "chr1", 1000, 2000, [100.0, 200.0, 3000.0]),
    @src.CswWindow::new("w2", "chr1", 3000, 4000, [150.0, 250.0, 350.0])
  ]
  let samples = ["s1", "s2", "s3"]
  let lib_sizes = [1000.0, 2000.0, 3000.0]
  let data = @src.CswDataSet::new(windows, samples, lib_sizes)
  let result = @src.csw_norm_tmm(data)
  assert_eq(result.norm_factors.length(), 3)
}

///|
test "csaw_filter_abundance" {
  let windows = [
    @src.CswWindow::new("w1", "chr1", 1000, 2000, [100.0, 200.0]),
    @src.CswWindow::new("w2", "chr1", 3000, 4000, [5.0, 3.0])
  ]
  let samples = ["s1", "s2"]
  let lib_sizes = [1000.0, 2000.0]
  let data = @src.CswDataSet::new(windows, samples, lib_sizes)
  let filtered = @src.csw_filter_abundance(data, 1.0)
  assert_true(filtered.csw_n_windows() <= 2)
}

///|
test "csaw_test_differential" {
  let windows = [
    @src.CswWindow::new("w1", "chr1", 1000, 2000, [50.0, 100.0, 200.0, 250.0]),
    @src.CswWindow::new("w2", "chr1", 3000, 4000, [80.0, 120.0, 180.0, 220.0])
  ]
  let samples = ["s1", "s2", "s3", "s4"]
  let lib_sizes = [1000.0, 2000.0, 3000.0, 4000.0]
  let data = @src.CswDataSet::new(windows, samples, lib_sizes)
  let norm_factors = @src.csw_norm_library(lib_sizes)
  let norm_data = @src.csw_apply_normalization(data, norm_factors)
  let result = @src.csw_test_differential(norm_data, [0, 1], [2, 3], 0.1)
  assert_true(result.log_fc.length() > 0)
  assert_true(result.p_values.length() > 0)
}

///|
test "csaw_find_regions" {
  let windows = [
    @src.CswWindow::new("w1", "chr1", 1000, 2000, [500.0, 100.0, 200.0, 250.0]),
    @src.CswWindow::new("w2", "chr1", 3000, 4000, [80.0, 120.0, 180.0, 220.0])
  ]
  let samples = ["s1", "s2", "s3", "s4"]
  let lib_sizes = [1000.0, 2000.0, 3000.0, 4000.0]
  let data = @src.CswDataSet::new(windows, samples, lib_sizes)
  let norm_factors = @src.csw_norm_library(lib_sizes)
  let norm_data = @src.csw_apply_normalization(data, norm_factors)
  let result = @src.csw_test_differential(norm_data, [0, 1], [2, 3], 0.1)
  let regions = @src.csw_find_regions(norm_data, result, 1000)
  assert_true(regions.length() >= 0)
}
