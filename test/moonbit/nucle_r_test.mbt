///|
/// nucleR module tests - Nucleosome positioning analysis

///|
test "nuc_create_example_track" {
  let track = @src.nuc_create_example_track()
  assert_eq(track.nuc_get_chrom(), "chr1")
  assert_true(track.nuc_get_length() > 0)

  let pos100 = track.nuc_get_signal(100)
  assert_true(!pos100.is_nan())

  let pos_outside = track.nuc_get_signal(-1)
  assert_eq(pos_outside, 0.0)
}

///|
test "nuc_create_example_result" {
  let track = @src.nuc_create_example_track()
  let params = @src.NucCallParams::new()
  let result = @src.nuc_call_nucleosomes(track, params)
  assert_true(result.nuc_n_count() >= 0)
  assert_true(result.nuc_mean_spacing() >= 0.0)
  assert_true(result.nuc_mean_occupancy() >= 0.0)
  assert_true(result.nuc_frac_well_positioned() >= 0.0 && result.nuc_frac_well_positioned() <= 1.0)
}

///|
test "nuc_smooth_signal" {
  let track = @src.nuc_create_example_track()
  let smoothed = @src.nuc_smooth_signal(track, 50)

  assert_eq(smoothed.nuc_get_length(), track.nuc_get_length())
  assert_eq(smoothed.nuc_get_chrom(), track.nuc_get_chrom())

  let orig_val = track.nuc_get_signal(250)
  let smooth_val = smoothed.nuc_get_signal(250)
  assert_true(!smooth_val.is_nan())
  assert_true(!orig_val.is_nan())
}

///|
test "nuc_find_peaks" {
  let track = @src.nuc_create_example_track()
  let smoothed = @src.nuc_smooth_signal(track, 30)
  let peaks = @src.nuc_find_peaks(smoothed, 0.3, 30)
  assert_true(peaks.length() > 0)
}

///|
test "nuc_call_nucleosomes_with_params" {
  let track = @src.nuc_create_example_track()
  let params = @src.NucCallParams::new()
  let p1 = params.nuc_set_smooth_window(30)
  let p2 = p1.nuc_set_peak_threshold(0.3)
  let p3 = p2.nuc_set_min_distance(30)
  let result = @src.nuc_call_nucleosomes(track, p3)

  assert_true(result.nuc_n_count() > 0)
  assert_true(result.nuc_mean_spacing() >= 0.0)
  assert_true(result.nuc_mean_occupancy() >= 0.0)
}

///|
test "nuc_call_result_methods" {
  let track = @src.nuc_create_example_track()
  let params = @src.NucCallParams::new()
  let p1 = params.nuc_set_smooth_window(30)
  let p2 = p1.nuc_set_peak_threshold(0.3)
  let p3 = p2.nuc_set_min_distance(30)
  let result = @src.nuc_call_nucleosomes(track, p3)
  let n = result.nuc_n_count()
  assert_true(n > 0)

  for i = 0; i < n; i = i + 1 {
    let pos = result.nuc_get(i)
    assert_true(pos.nuc_chrom().length() > 0)
    assert_true(pos.nuc_start() < pos.nuc_end())
    assert_true(pos.nuc_center() > 0.0)
    assert_true(pos.nuc_width() >= 0.0)
    assert_true(pos.nuc_occupancy() >= 0.0 && pos.nuc_occupancy() <= 1.0)
    assert_true(pos.nuc_positioning_score() >= 0.0)
  }
}

///|
test "nuc_compare_positioning_identical" {
  let track = @src.nuc_create_example_track()
  let params = @src.NucCallParams::new()
  let p1 = params.nuc_set_smooth_window(30)
  let p2 = p1.nuc_set_peak_threshold(0.3)
  let p3 = p2.nuc_set_min_distance(30)
  let result = @src.nuc_call_nucleosomes(track, p3)
  let dynamic = @src.nuc_compare_positioning(result, result, 10)

  assert_true(dynamic.nuc_n_shared() >= 0)
  assert_eq(dynamic.nuc_n_gained(), 0)
  assert_eq(dynamic.nuc_n_lost(), 0)
  assert_true(dynamic.nuc_frac_changed() >= 0.0)
  assert_eq(dynamic.nuc_direction(), 0)
}

///|
test "nuc_position_methods" {
  let pos = @src.NucPosition::new(
    "chr1", 100, 300, 200.0, 200.0, 0.8, 0.01, true
  )
  assert_eq(pos.nuc_chrom(), "chr1")
  assert_eq(pos.nuc_start(), 100)
  assert_eq(pos.nuc_end(), 300)
  assert_eq(pos.nuc_center(), 200.0)
  assert_eq(pos.nuc_width(), 200.0)
  assert_eq(pos.nuc_occupancy(), 0.8)
  assert_eq(pos.nuc_positioning_score(), 0.01)
  assert_true(pos.nuc_is_well_positioned())

  let pos2 = @src.NucPosition::new(
    "chr2", 500, 800, 650.0, 300.0, 0.5, 0.003, false
  )
  assert_eq(pos2.nuc_chrom(), "chr2")
  assert_true(!pos2.nuc_is_well_positioned())
}

///|
test "nuc_filter_by_occupancy" {
  let track = @src.nuc_create_example_track()
  let params = @src.NucCallParams::new()
  let p1 = params.nuc_set_smooth_window(30)
  let p2 = p1.nuc_set_peak_threshold(0.3)
  let p3 = p2.nuc_set_min_distance(30)
  let result = @src.nuc_call_nucleosomes(track, p3)
  let filtered = result.nuc_filter_by_occupancy(0.1)
  assert_true(filtered.nuc_n_count() <= result.nuc_n_count())
}

///|
test "nuc_signal_track_methods" {
  let track = @src.nuc_create_example_track()
  assert_eq(track.nuc_get_chrom(), "chr1")
  assert_true(track.nuc_get_length() > 0)

  let pos = track.nuc_get_position(100)
  assert_true(pos >= 0)

  let signal = track.nuc_get_signal(pos)
  assert_true(!signal.is_nan())
}

///|
test "nuc_call_params_new" {
  let params = @src.NucCallParams::new()
  assert_eq(params.smooth_window, 100)
  assert_true(params.peak_threshold > 0.0)
  assert_eq(params.min_distance, 100)
  assert_eq(params.max_distance, 300)
  assert_true(params.occupancy_threshold > 0.0)
  assert_true(params.width_threshold > 0.0)
}

///|
test "nuc_call_params_setters" {
  let params = @src.NucCallParams::new()
  let p1 = params.nuc_set_smooth_window(50)
  assert_eq(p1.smooth_window, 50)

  let p2 = params.nuc_set_peak_threshold(0.3)
  assert_true((p2.peak_threshold - 0.3).abs() < 0.001)

  let p3 = params.nuc_set_min_distance(80)
  assert_eq(p3.min_distance, 80)

  let p4 = params.nuc_set_occupancy_threshold(0.5)
  assert_true((p4.occupancy_threshold - 0.5).abs() < 0.001)
}

///|
test "nuc_dynamic_result_methods" {
  let track = @src.nuc_create_example_track()
  let params = @src.NucCallParams::new()
  let p1 = params.nuc_set_smooth_window(30)
  let p2 = p1.nuc_set_peak_threshold(0.3)
  let p3 = p2.nuc_set_min_distance(30)
  let result = @src.nuc_call_nucleosomes(track, p3)
  let dynamic = @src.nuc_compare_positioning(result, result, 10)

  assert_true(dynamic.nuc_n_shared() >= 0)
  assert_eq(dynamic.nuc_n_gained(), 0)
  assert_eq(dynamic.nuc_n_lost(), 0)
  assert_true(dynamic.nuc_frac_changed() >= 0.0 && dynamic.nuc_frac_changed() <= 1.0)
  assert_eq(dynamic.nuc_direction(), 0)
}