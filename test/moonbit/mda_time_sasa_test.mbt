///|
/// Tests for mda_time_sasa.mbt (Shrake-Rupley SASA over trajectory).

///|
test "mda_time_dependent_sasa_isolated_atom_full_area" {
  // A lone atom exposes every sphere point -> SASA = 4*pi*r^2.
  let traj = [[[1.0, 2.0, 3.0]]]
  let sasa = @src.mda_time_dependent_sasa(traj, [2.0], n_points=10)
  assert_eq(sasa.length(), 1)
  assert_eq(sasa[0].length(), 1)
  assert_true(
    @moonbit.assert_num_close(sasa[0][0], 4.0 * 3.141592653589793 * 4.0, 1.0e-9),
  )
}

///|
test "mda_time_dependent_sasa_partially_buried" {
  // Two radius-2 atoms: B sits exactly on A's test point in the +y
  // direction (n_points=2 places both points at the poles). Each atom
  // has one buried and one exposed point -> half of 4*pi*r^2.
  let traj = [[[0.0, 0.0, 0.0], [0.0, 2.0, 0.0]]]
  let sasa = @src.mda_time_dependent_sasa(traj, [2.0, 2.0], n_points=2)
  let half = 0.5 * 4.0 * 3.141592653589793 * 4.0
  assert_true(@moonbit.assert_num_close(sasa[0][0], half, 1.0e-9))
  assert_true(@moonbit.assert_num_close(sasa[0][1], half, 1.0e-9))
}

///|
test "mda_time_dependent_sasa_empty_frame" {
  let sasa = @src.mda_time_dependent_sasa([[]], [1.0])
  assert_eq(sasa.length(), 1)
  assert_eq(sasa[0].length(), 0)
}

///|
test "mda_time_dependent_sasa_multi_frame_shapes" {
  // Two frames of one isolated atom: full area in both frames.
  let traj = [[[0.0, 0.0, 0.0]], [[5.0, 5.0, 5.0]]]
  let sasa = @src.mda_time_dependent_sasa(traj, [1.0], n_points=4)
  assert_eq(sasa.length(), 2)
  let full = 4.0 * 3.141592653589793
  assert_true(@moonbit.assert_num_close(sasa[0][0], full, 1.0e-9))
  assert_true(@moonbit.assert_num_close(sasa[1][0], full, 1.0e-9))
}

///|
test "mda_total_sasa_trajectory_sums" {
  let sasa_traj = [[1.0, 2.0, 3.0], [0.5, 0.5], []]
  let totals = @src.mda_total_sasa_trajectory(sasa_traj)
  assert_eq(totals.length(), 3)
  assert_true(@moonbit.assert_num_close(totals[0], 6.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(totals[1], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(totals[2], 0.0, 1.0e-12))
}

///|
test "mda_burial_exposure_events_transitions" {
  // SASA 20 -> 5 -> 15 -> 8 with threshold 10:
  // Burial at frame 1, Exposure at frame 2, Burial at frame 3.
  let events = @src.mda_burial_exposure_events(
    [20.0, 5.0, 15.0, 8.0],
    threshold=10.0,
  )
  assert_eq(events.length(), 3)
  assert_eq(events[0].0, 1)
  assert_eq(events[0].1, "Burial")
  assert_eq(events[1].0, 2)
  assert_eq(events[1].1, "Exposure")
  assert_eq(events[2].0, 3)
  assert_eq(events[2].1, "Burial")
}

///|
test "mda_burial_exposure_events_no_change" {
  // Constantly exposed series produces no events.
  assert_eq(
    @src.mda_burial_exposure_events([20.0, 25.0, 30.0], threshold=10.0).length(),
    0,
  )
  // Short series: no events.
  assert_eq(@src.mda_burial_exposure_events([20.0], threshold=10.0).length(), 0)
}

///|
test "mda_sasa_fluctuation_stddev" {
  // Atom 0 alternates 0 and 10 (population std = 5); atom 1 is constant.
  let sasa_traj = [[0.0, 10.0], [10.0, 10.0]]
  let fl = @src.mda_sasa_fluctuation(sasa_traj)
  assert_eq(fl.length(), 2)
  assert_true(@moonbit.assert_num_close(fl[0], 5.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(fl[1], 0.0, 1.0e-12))
}

///|
test "mda_sasa_fluctuation_empty" {
  assert_eq(@src.mda_sasa_fluctuation([]).length(), 0)
}
