///|
/// Tests for mda_water_analysis.mbt (water bridges, residence, occupancy).

///|
test "mda_find_water_bridges_between_two_residues" {
  // Protein atoms in residues 0 and 1 near a single water; a third atom
  // belongs to residue 1 again. The water bridges residues 0 and 1 via
  // two atom pairs (pairs within the same residue are skipped).
  let protein = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.5, 0.0, 0.0]]
  let water = [[2.0, 0.0, 0.0]]
  let res_of = [0, 1, 1]
  let bridges = @src.mda_find_water_bridges(protein, water, res_of, cutoff=3.0)
  assert_eq(bridges.length(), 2)
  assert_eq(bridges[0].residue1, 0)
  assert_eq(bridges[0].residue2, 1)
  assert_eq(bridges[0].water_idx, 0)
  assert_true(@moonbit.assert_num_close(bridges[0].distance1, 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(bridges[0].distance2, 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(bridges[1].distance2, 2.5, 1.0e-9))
  // One unique residue pair.
  assert_eq(@src.mda_water_bridge_pairs(bridges), 1)
}

///|
test "mda_find_water_bridges_same_residue_skipped" {
  // Both nearby atoms belong to the same residue: no bridge.
  let protein = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let water = [[0.5, 0.0, 0.0]]
  let bridges = @src.mda_find_water_bridges(protein, water, [0, 0], cutoff=3.0)
  assert_eq(bridges.length(), 0)
}

///|
test "mda_find_water_bridges_water_too_far" {
  // Water is 5 A away from everything: beyond the 3 A cutoff.
  let protein = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
  let water = [[5.0, 0.0, 0.0]]
  assert_eq(@src.mda_find_water_bridges(protein, water, [0, 1]).length(), 0)
}

///|
test "mda_find_water_bridges_three_residues_dedup_pairs" {
  // One water sits within cutoff of atoms of residues 0, 1, 2 -> three
  // bridges per water; duplicating the water doubles the bridge count
  // but the unique residue-pair count stays 3.
  let protein = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]]
  let water = [[1.5, 1.0, 0.0], [1.5, 1.0, 0.0]]
  let bridges = @src.mda_find_water_bridges(protein, water, [0, 1, 2])
  assert_eq(bridges.length(), 6)
  assert_eq(@src.mda_water_bridge_pairs(bridges), 3)
}

///|
test "mda_water_residence_time_two_visits" {
  // Water at x = 1, 1, 5, 1, 1 relative to a residue at the origin with
  // cutoff 2: two residence episodes of 2 frames each -> average 2.
  let traj_water = [
    [[1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [[5.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
  ]
  let residue = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
  ]
  let rt = @src.mda_water_residence_time(traj_water, residue, cutoff=2.0)
  assert_true(@moonbit.assert_num_close(rt, 2.0, 1.0e-9))
}

///|
test "mda_water_residence_time_degenerate" {
  // Water never enters the cutoff -> 0; empty trajectory -> 0.
  let traj_water = [[[9.0, 0.0, 0.0]], [[9.0, 0.0, 0.0]]]
  let residue = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_water_residence_time(traj_water, residue),
      0.0,
      1.0e-12,
    ),
  )
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_water_residence_time([], []),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "mda_water_residence_time_continuous_stay" {
  // A water that stays inside for the whole trajectory: one event whose
  // duration equals the trajectory length.
  let traj_water = [[[1.0, 0.0, 0.0]], [[1.5, 0.0, 0.0]], [[1.0, 0.0, 0.0]]]
  let residue = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let rt = @src.mda_water_residence_time(traj_water, residue, cutoff=2.0)
  assert_true(@moonbit.assert_num_close(rt, 3.0, 1.0e-9))
}

///|
test "mda_water_occupancy_fractions" {
  // Water 0 is near for 3 of 4 frames (occupancy 0.75); water 1 never.
  let traj_water = [
    [[1.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
    [[5.0, 0.0, 0.0], [9.0, 0.0, 0.0]],
  ]
  let residue = [
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
  ]
  let occ = @src.mda_water_occupancy(traj_water, residue, cutoff=2.0)
  assert_eq(occ.length(), 2)
  assert_true(@moonbit.assert_num_close(occ[0], 0.75, 1.0e-12))
  assert_true(@moonbit.assert_num_close(occ[1], 0.0, 1.0e-12))
}

///|
test "mda_water_occupancy_empty" {
  let occ = @src.mda_water_occupancy([], [])
  assert_eq(occ.length(), 0)
}
