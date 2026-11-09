///|
/// Tests for mda_interactions.mbt (salt bridges, contacts, stacking, coordination).

///|
test "mda_salt_bridges_lys_glu_pair" {
  // LYS-GLU at 3 A (within 4) forms a bridge; other pairs are too far
  // or not oppositely charged.
  let coords = [
    (0.0, 0.0, 0.0),
    (3.0, 0.0, 0.0),
    (10.0, 0.0, 0.0),
    (0.0, 0.0, 3.5),
    (20.0, 0.0, 0.0),
  ]
  let res_ids = [100, 101, 102, 103, 104]
  let res_names = ["LYS", "GLU", "ALA", "ARG", "ASP"]
  let bridges = @src.mda_salt_bridges(coords, res_ids, res_names)
  assert_eq(bridges.length(), 1)
  assert_eq(bridges[0].res1, 100)
  assert_eq(bridges[0].res2, 101)
  assert_true(@moonbit.assert_num_close(bridges[0].distance, 3.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(bridges[0].occupancy, 1.0, 1.0e-12))
}

///|
test "mda_salt_bridges_cutoff_excluded" {
  // GLU-ARG distance sqrt(21.25) ~ 4.61 > 4 -> no bridges.
  let coords = [(0.0, 0.0, 0.0), (0.0, 0.0, 4.61)]
  let res_ids = [1, 2]
  let res_names = ["GLU", "ARG"]
  assert_eq(@src.mda_salt_bridges(coords, res_ids, res_names).length(), 0)
}

///|
test "mda_salt_bridges_traj_occupancy" {
  // Frame 0: LYS-GLU at 1 A (contact); frame 1: separated by 10 A.
  // Occupancy over 2 frames = 0.5.
  let frames = [
    [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0)],
  ]
  let res_ids = [7, 8]
  let res_names = ["LYS", "GLU"]
  let bridges = @src.mda_salt_bridges_traj(frames, res_ids, res_names)
  assert_eq(bridges.length(), 1)
  assert_eq(bridges[0].res1, 7)
  assert_eq(bridges[0].res2, 8)
  assert_true(@moonbit.assert_num_close(bridges[0].occupancy, 0.5, 1.0e-12))
}

///|
test "mda_salt_bridges_traj_empty" {
  let bridges = @src.mda_salt_bridges_traj([], [], [])
  assert_eq(bridges.length(), 0)
}

///|
test "mda_hydrophobic_contacts_ala_val_leu" {
  // ALA-VAL at 2 A and VAL-LEU at 4 A (cutoff 5); GLY is not hydrophobic.
  let coords = [
    (0.0, 0.0, 0.0),
    (2.0, 0.0, 0.0),
    (6.0, 0.0, 0.0),
    (0.5, 0.0, 0.0),
  ]
  let res_ids = [10, 20, 30, 40]
  let res_names = ["ALA", "VAL", "LEU", "GLY"]
  let contacts = @src.mda_hydrophobic_contacts(coords, res_ids, res_names)
  assert_eq(contacts.length(), 2)
  assert_eq(contacts[0].res1, 10)
  assert_eq(contacts[0].res2, 20)
  assert_true(@moonbit.assert_num_close(contacts[0].distance, 2.0, 1.0e-9))
  assert_eq(contacts[1].res1, 20)
  assert_eq(contacts[1].res2, 30)
  assert_true(@moonbit.assert_num_close(contacts[1].distance, 4.0, 1.0e-9))
}

///|
test "mda_pi_stacking_face_to_face" {
  // PHE-TYR centroids 3 A apart -> face-to-face (angle hardcoded 0).
  // TRP at 10 A exceeds the 6 A cutoff.
  let coords = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (10.0, 0.0, 0.0)]
  let res_ids = [5, 6, 7]
  let res_names = ["PHE", "TYR", "TRP"]
  let stackings = @src.mda_pi_stacking(coords, res_ids, res_names)
  assert_eq(stackings.length(), 1)
  assert_eq(stackings[0].ring1_res, 5)
  assert_eq(stackings[0].ring2_res, 6)
  assert_true(@moonbit.assert_num_close(stackings[0].distance, 3.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(stackings[0].angle, 0.0, 1.0e-12))
  assert_eq(stackings[0].stacking_type, "face-to-face")
}

///|
test "mda_coordination_number_linear_chain" {
  // Atoms at x = 0, 1, 2, 5 with cutoff 1.5.
  let coords = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (2.0, 0.0, 0.0),
    (5.0, 0.0, 0.0),
  ]
  let result = @src.mda_coordination_number(coords, cutoff=1.5)
  assert_eq(result.length(), 4)
  assert_eq(result[0].atom_idx, 0)
  assert_eq(result[0].n_neighbors, 1)
  assert_eq(result[0].neighbor_indices[0], 1)
  assert_eq(result[1].n_neighbors, 2)
  assert_eq(result[1].neighbor_indices[0], 0)
  assert_eq(result[1].neighbor_indices[1], 2)
  assert_eq(result[2].n_neighbors, 1)
  assert_eq(result[3].n_neighbors, 0)
  assert_eq(result[3].neighbor_indices.length(), 0)
}

///|
test "mda_contact_map_trajectory_frequencies" {
  // 3 atoms, 2 frames. d01 = 1 then 7 (both <= 8) -> freq 1.0.
  // d12 = 9 then 3 -> freq 0.5. d02 = 10, 10 -> 0.
  let frames = [
    [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (10.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0), (7.0, 0.0, 0.0), (10.0, 0.0, 0.0)],
  ]
  let freq = @src.mda_contact_map_trajectory(frames, cutoff=8.0)
  assert_eq(freq.length(), 3)
  assert_true(@moonbit.assert_num_close(freq[0][1], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(freq[1][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(freq[1][2], 0.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(freq[0][2], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(freq[2][2], 0.0, 1.0e-12))
}

///|
test "mda_contact_map_trajectory_empty" {
  assert_eq(@src.mda_contact_map_trajectory([]).length(), 0)
}

///|
test "mda_contact_lifetimes_interrupted_contact" {
  // Pair (0,1): contact in frames 0,1 and 3 (distance 1, 1, 10, 1
  // with cutoff 8): total contact frames = 3 -> avg lifetime 3/4 = 0.75.
  // Atom 2 stays at 100 A: no lifetime entry.
  let frames = [
    [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (100.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (100.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (100.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (100.0, 0.0, 0.0)],
  ]
  let lifetimes = @src.mda_contact_lifetimes(frames, cutoff=8.0)
  assert_eq(lifetimes.length(), 1)
  assert_eq(lifetimes[0].0, 0)
  assert_eq(lifetimes[0].1, 1)
  assert_true(@moonbit.assert_num_close(lifetimes[0].2, 0.75, 1.0e-12))
}

///|
test "mda_contact_lifetimes_empty" {
  assert_eq(@src.mda_contact_lifetimes([]).length(), 0)
}
