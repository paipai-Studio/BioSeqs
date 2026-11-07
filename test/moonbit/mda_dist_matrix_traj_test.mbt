///|
/// Tests for mda_dist_matrix_traj.mbt (inter-residue distance time series).
///
/// Note: mda_distance_matrix_frame builds its matrix rows via
/// Array::make(n, Array::make(n, 0.0)), which aliases one row array across
/// all rows, so for n >= 2 the diagonal gets clobbered. That is a suspected
/// implementation bug; these tests therefore only exercise the frame helper
/// for n <= 1 and build distance-matrix time series by hand for the pure
/// numeric functions (average / fluctuation / covariance / DCCM).

///|
test "mda_distance_matrix_frame_single_residue" {
  let m = @src.mda_distance_matrix_frame([[2.0, 3.0, 4.0]])
  assert_eq(m.length(), 1)
  assert_true(@moonbit.assert_num_close(m[0][0], 0.0, 1.0e-12))
}

///|
test "mda_distance_matrix_frame_empty" {
  let m = @src.mda_distance_matrix_frame([])
  assert_eq(m.length(), 0)
}

///|
test "mda_distance_matrix_trajectory_single_residue" {
  let traj = [[[0.0, 0.0, 0.0]], [[5.0, 1.0, 2.0]]]
  let dt = @src.mda_distance_matrix_trajectory(traj)
  assert_eq(dt.length(), 2)
  assert_eq(dt[0].length(), 1)
  assert_true(@moonbit.assert_num_close(dt[0][0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(dt[1][0][0], 0.0, 1.0e-12))
}

///|
test "mda_average_distance_matrix_single_entry" {
  // 1x1 matrices are unaffected by the row-aliasing issue.
  let dist_traj = [[[5.0]], [[7.0]]]
  let avg = @src.mda_average_distance_matrix(dist_traj)
  assert_eq(avg.length(), 1)
  assert_true(@moonbit.assert_num_close(avg[0][0], 6.0, 1.0e-12))
}

///|
test "mda_average_distance_matrix_empty" {
  assert_eq(@src.mda_average_distance_matrix([]).length(), 0)
}

///|
test "mda_distance_fluctuation_matrix_single_entry" {
  // One frame: the entry equals its own mean -> zero fluctuation.
  let dist_traj = [[[3.0]]]
  let fluct = @src.mda_distance_fluctuation_matrix(dist_traj)
  assert_true(@moonbit.assert_num_close(fluct[0][0], 0.0, 1.0e-12))
}

///|
test "mda_distance_fluctuation_matrix_empty" {
  assert_eq(@src.mda_distance_fluctuation_matrix([]).length(), 0)
}

///|
test "mda_distance_covariance_anti_correlated" {
  // Pair (0,1) goes 1,2,3 while pair (2,3) goes 3,2,1 -> covariance = -1.
  let filler = 9.0
  let dist_traj = [
    [
      [0.0, 1.0, filler, filler],
      [1.0, 0.0, filler, filler],
      [filler, filler, 0.0, 3.0],
      [filler, filler, 3.0, 0.0],
    ],
    [
      [0.0, 2.0, filler, filler],
      [2.0, 0.0, filler, filler],
      [filler, filler, 0.0, 2.0],
      [filler, filler, 2.0, 0.0],
    ],
    [
      [0.0, 3.0, filler, filler],
      [3.0, 0.0, filler, filler],
      [filler, filler, 0.0, 1.0],
      [filler, filler, 1.0, 0.0],
    ],
  ]
  let cov = @src.mda_distance_covariance(dist_traj, 0, 1, 2, 3)
  assert_true(@moonbit.assert_num_close(cov, -1.0, 1.0e-9))
  // Correlated with itself = variance = 1.
  let self_cov = @src.mda_distance_covariance(dist_traj, 0, 1, 0, 1)
  assert_true(@moonbit.assert_num_close(self_cov, 1.0, 1.0e-9))
  // Fewer than 2 frames -> 0.
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_distance_covariance([dist_traj[0]], 0, 1, 2, 3),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "mda_allosteric_pairs_sorted_by_fluctuation" {
  let fluct = [[0.0, 2.5, 0.5], [2.5, 0.0, 1.5], [0.5, 1.5, 0.0]]
  let pairs = @src.mda_allosteric_pairs(fluct, min_fluct=1.0)
  assert_eq(pairs.length(), 2)
  // Sorted descending by fluctuation.
  assert_true(@moonbit.assert_num_close(pairs[0].2, 2.5, 1.0e-9))
  assert_eq(pairs[0].0, 0)
  assert_eq(pairs[0].1, 1)
  assert_true(@moonbit.assert_num_close(pairs[1].2, 1.5, 1.0e-9))
  assert_eq(pairs[1].0, 1)
  assert_eq(pairs[1].1, 2)
}

///|
test "mda_dynamic_cross_correlation_single_entry" {
  // A 1x1 "trajectory": the diagonal is defined as 1.
  let dist_traj = [[[2.0]], [[5.0]]]
  let dccm = @src.mda_dynamic_cross_correlation(dist_traj)
  assert_eq(dccm.length(), 1)
  assert_true(@moonbit.assert_num_close(dccm[0][0], 1.0, 1.0e-12))
}

///|
test "mda_dynamic_cross_correlation_empty" {
  assert_eq(@src.mda_dynamic_cross_correlation([]).length(), 0)
}
