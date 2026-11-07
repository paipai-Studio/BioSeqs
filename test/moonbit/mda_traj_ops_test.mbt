///|
/// Tests for mda_traj_ops.mbt (trajectory merge/slice/strip/recenter,
/// wrap/unwrap, mean structure, alignment, B-factors).

///|
test "mda_traj_merge_concatenates" {
  let f1 = [(0.0, 0.0, 0.0)]
  let f2 = [(1.0, 1.0, 1.0)]
  let f3 = [(2.0, 2.0, 2.0)]
  let merged = @src.mda_traj_merge([[f1], [f2, f3], []])
  assert_eq(merged.length(), 3)
  assert_true(@moonbit.assert_num_close(merged[0][0].0, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(merged[1][0].0, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(merged[2][0].0, 2.0, 1.0e-12))
  // Merging nothing yields nothing.
  assert_eq(@src.mda_traj_merge([]).length(), 0)
}

///|
test "mda_traj_slice_ranges" {
  let frames = [
    [(0.0, 0.0, 0.0)],
    [(1.0, 0.0, 0.0)],
    [(2.0, 0.0, 0.0)],
    [(3.0, 0.0, 0.0)],
    [(4.0, 0.0, 0.0)],
  ]
  // Contiguous slice [1, 4).
  let s1 = @src.mda_traj_slice(frames, 1, 4)
  assert_eq(s1.length(), 3)
  assert_true(@moonbit.assert_num_close(s1[0][0].0, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(s1[2][0].0, 3.0, 1.0e-12))
  // Stepped slice.
  let s2 = @src.mda_traj_slice(frames, 0, 5, step=2)
  assert_eq(s2.length(), 3)
  assert_true(@moonbit.assert_num_close(s2[1][0].0, 2.0, 1.0e-12))
  // Clamping: negative start -> 0; stop beyond the end -> n.
  let s3 = @src.mda_traj_slice(frames, -5, 99)
  assert_eq(s3.length(), 5)
  // Empty window.
  assert_eq(@src.mda_traj_slice(frames, 3, 1).length(), 0)
}

///|
test "mda_traj_strip_solvent_default_names" {
  let frames = [
    [(0.0, 0.0, 0.0), (9.0, 9.0, 9.0), (1.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0), (8.0, 8.0, 8.0), (1.0, 0.0, 0.0)],
  ]
  let stripped = @src.mda_traj_strip_solvent(frames, ["CA", "HOH", "N"])
  assert_eq(stripped.length(), 2)
  assert_eq(stripped[0].length(), 2)
  assert_true(@moonbit.assert_num_close(stripped[0][1].0, 1.0, 1.0e-12))
  // Custom solvent list.
  let stripped2 = @src.mda_traj_strip_solvent(frames, ["CA", "ION", "N"], solvent_names=[
    "ION",
  ])
  assert_eq(stripped2[0].length(), 2)
}

///|
test "mda_traj_recenter_to_origin" {
  let frames = [[(1.0, 1.0, 1.0), (3.0, 1.0, 1.0)]]
  let rec = @src.mda_traj_recenter(frames)
  assert_eq(rec.length(), 1)
  assert_true(@moonbit.assert_num_close(rec[0][0].0, -1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(rec[0][1].0, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(rec[0][0].1, 0.0, 1.0e-9))
  // Empty frames pass through untouched.
  let rec2 = @src.mda_traj_recenter([[]])
  assert_eq(rec2.length(), 1)
  assert_eq(rec2[0].length(), 0)
}

///|
test "mda_wrap_coords_into_box" {
  let wrapped = @src.mda_wrap_coords(
    [(12.0, -3.0, 5.0), (0.0, 10.0, -10.0)],
    (10.0, 10.0, 10.0),
  )
  assert_eq(wrapped.length(), 2)
  assert_true(@moonbit.assert_num_close(wrapped[0].0, 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(wrapped[0].1, 7.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(wrapped[0].2, 5.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(wrapped[1].0, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(wrapped[1].1, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(wrapped[1].2, 0.0, 1.0e-9))
}

///|
test "mda_unwrap_coords_minimum_image" {
  // Displacements of +6 and +11 in a 10-box map to -4 and +1 relative to
  // the reference at x = 5, giving unwrapped x of 1 and 6.
  let ref_struct = [(5.0, 5.0, 5.0), (5.0, 5.0, 5.0)]
  let coords = [(11.0, 5.0, 5.0), (16.0, 8.0, 5.0)]
  let unwrapped = @src.mda_unwrap_coords(coords, ref_struct, (10.0, 10.0, 10.0))
  assert_true(@moonbit.assert_num_close(unwrapped[0].0, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(unwrapped[0].1, 5.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(unwrapped[1].0, 6.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(unwrapped[1].1, 8.0, 1.0e-9))
}

///|
test "mda_mean_structure_averages_frames" {
  let frames = [
    [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
    [(2.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
  ]
  let mean = @src.mda_mean_structure(frames)
  assert_eq(mean.length(), 2)
  assert_true(@moonbit.assert_num_close(mean[0].0, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(mean[1].2, 0.0, 1.0e-9))
  assert_eq(@src.mda_mean_structure([]).length(), 0)
}

///|
test "mda_rmsd_to_mean_symmetric" {
  // Single atom at x = 0, 1, 2: mean 1, RMSDs [1, 0, 1].
  let frames = [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)], [(2.0, 0.0, 0.0)]]
  let rmsds = @src.mda_rmsd_to_mean(frames)
  assert_eq(rmsds.length(), 3)
  assert_true(@moonbit.assert_num_close(rmsds[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(rmsds[1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(rmsds[2], 1.0, 1.0e-9))
}

///|
test "mda_kabsch_align_returns_identity" {
  // The simplified implementation always returns the identity rotation.
  let rot = @src.mda_kabsch_align([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)], [
    (5.0, 5.0, 5.0),
    (6.0, 5.0, 5.0),
  ])
  assert_eq(rot.length(), 3)
  assert_true(@moonbit.assert_num_close(rot[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(rot[1][1], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(rot[2][2], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(rot[0][1], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(rot[1][2], 0.0, 1.0e-12))
}

///|
test "mda_bfactor_analysis_conversion" {
  // RMSF [1, 2, 3]: B_i = (8*pi^2/3) * rmsf_i^2; z-scores of [1,2,3].
  let results = @src.mda_bfactor_analysis([1.0, 2.0, 3.0])
  assert_eq(results.length(), 3)
  let scale = 8.0 * 3.141592653589793 * 3.141592653589793 / 3.0
  assert_eq(results[0].atom_idx, 0)
  assert_eq(results[2].atom_idx, 2)
  assert_true(@moonbit.assert_num_close(results[0].rmsf, 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(results[0].bfactor, scale, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(results[2].bfactor, scale * 9.0, 1.0e-9),
  )
  let std = (2.0 / 3.0).sqrt()
  assert_true(
    @moonbit.assert_num_close(results[0].normalized, -1.0 / std, 1.0e-9),
  )
  assert_true(@moonbit.assert_num_close(results[1].normalized, 0.0, 1.0e-12))
  assert_true(
    @moonbit.assert_num_close(results[2].normalized, 1.0 / std, 1.0e-9),
  )
}

///|
test "mda_bfactor_analysis_empty" {
  assert_eq(@src.mda_bfactor_analysis([]).length(), 0)
}
