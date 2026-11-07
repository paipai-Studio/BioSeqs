///|
/// Tests for mda_rmsd_rmsf.mbt (RMSD, RMSF, RMSD matrix, Kabsch centering).

///|
test "mda_rmsd_direct_distance" {
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_rmsd([(0.0, 0.0, 0.0)], [(3.0, 4.0, 0.0)]),
      5.0,
      1.0e-9,
    ),
  )
}

///|
test "mda_rmsd_not_translation_invariant" {
  // Plain RMSD does NOT superpose: a shift of 5 gives RMSD 5.
  let a = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
  let b = [(5.0, 0.0, 0.0), (6.0, 0.0, 0.0)]
  assert_true(@moonbit.assert_num_close(@src.mda_rmsd(a, b), 5.0, 1.0e-9))
  // Identical structures -> 0.
  assert_true(@moonbit.assert_num_close(@src.mda_rmsd(a, a), 0.0, 1.0e-12))
}

///|
test "mda_rmsd_degenerate_inputs" {
  // Mismatched lengths and empty inputs return 0.
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_rmsd([(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]),
      0.0,
      1.0e-12,
    ),
  )
  assert_true(@moonbit.assert_num_close(@src.mda_rmsd([], []), 0.0, 1.0e-12))
}

///|
test "mda_rmsd_aligned_translation_invariant" {
  // Aligned RMSD removes the centroid offset -> 0 for pure translation.
  let a = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
  let b = [(5.0, 5.0, 5.0), (6.0, 5.0, 5.0)]
  assert_true(
    @moonbit.assert_num_close(@src.mda_rmsd_aligned(a, b), 0.0, 1.0e-9),
  )
}

///|
test "mda_rmsd_aligned_perpendicular_vectors" {
  // a centered: (-1,0,0),(1,0,0); b centered: (0,-1,0),(0,1,0).
  // Per-atom squared distance 2 -> RMSD = sqrt(2).
  let a = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
  let b = [(0.0, 0.0, 0.0), (0.0, 2.0, 0.0)]
  assert_true(
    @moonbit.assert_num_close(@src.mda_rmsd_aligned(a, b), 2.0.sqrt(), 1.0e-9),
  )
}

///|
test "mda_rmsf_with_explicit_reference" {
  // Atom oscillates between 0 and 2; reference at 0 -> RMSF = sqrt(2).
  let frames = [[(0.0, 0.0, 0.0)], [(2.0, 0.0, 0.0)]]
  let rmsf = @src.mda_rmsf(frames, ref_coords=[(0.0, 0.0, 0.0)])
  assert_eq(rmsf.length(), 1)
  assert_true(@moonbit.assert_num_close(rmsf[0], 2.0.sqrt(), 1.0e-9))
}

///|
test "mda_rmsf_default_average_reference" {
  // Without a reference the average structure (x=1) is used -> RMSF = 1.
  let frames = [[(0.0, 0.0, 0.0)], [(2.0, 0.0, 0.0)]]
  let rmsf = @src.mda_rmsf(frames)
  assert_true(@moonbit.assert_num_close(rmsf[0], 1.0, 1.0e-9))
}

///|
test "mda_rmsf_empty" {
  assert_eq(@src.mda_rmsf([]).length(), 0)
  assert_eq(@src.mda_rmsf([[]]).length(), 0)
}

///|
test "mda_rmsd_matrix_unaligned" {
  // Single-atom frames at x = 0, 1, 2: pairwise distances 1, 2, 1.
  let frames = [[(0.0, 0.0, 0.0)], [(1.0, 0.0, 0.0)], [(2.0, 0.0, 0.0)]]
  let m = @src.mda_rmsd_matrix(frames)
  assert_eq(m.length(), 3)
  assert_true(@moonbit.assert_num_close(m[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(m[0][1], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(m[0][2], 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(m[1][2], 1.0, 1.0e-9))
  // Symmetric with zero diagonal.
  assert_true(@moonbit.assert_num_close(m[1][0], m[0][1], 1.0e-12))
  assert_true(@moonbit.assert_num_close(m[2][0], m[0][2], 1.0e-12))
}

///|
test "mda_rmsd_matrix_aligned_vs_unaligned" {
  // f1 is f0 translated by 5; f2 is f0 with the second atom pulled to 3.
  // Aligned RMSD: d(f0,f1) = 0 (same centered shape), d(f0,f2) = 1.
  // Unaligned: d(f0,f1) = 5, d(f0,f2) = sqrt(2).
  let f0 = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
  let f1 = [(5.0, 0.0, 0.0), (6.0, 0.0, 0.0)]
  let f2 = [(0.0, 0.0, 0.0), (3.0, 0.0, 0.0)]
  let ma = @src.mda_rmsd_matrix([f0, f1, f2], aligned=true)
  assert_true(@moonbit.assert_num_close(ma[0][1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(ma[0][2], 1.0, 1.0e-9))
  let mu = @src.mda_rmsd_matrix([f0, f1, f2])
  assert_true(@moonbit.assert_num_close(mu[0][1], 5.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(mu[0][2], 2.0.sqrt(), 1.0e-9))
}

///|
test "mda_kabsch_translation_only_fit" {
  // b is a translated by (10, -2, 3): the centroid translation is
  // (-10, 2, -3), rotation stays identity, RMSD 0.
  let a = [(0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 4.0, 0.0)]
  let b = [(10.0, -2.0, 3.0), (12.0, -2.0, 3.0), (10.0, 2.0, 3.0)]
  let k = @src.mda_kabsch(a, b)
  assert_true(@moonbit.assert_num_close(k.rmsd, 0.0, 1.0e-9))
  let (tx, ty, tz) = k.translation
  assert_true(@moonbit.assert_num_close(tx, -10.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(ty, 2.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(tz, -3.0, 1.0e-9))
  // Identity rotation.
  assert_true(@moonbit.assert_num_close(k.rotation[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(k.rotation[1][1], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(k.rotation[0][1], 0.0, 1.0e-12))
}

///|
test "mda_kabsch_degenerate_inputs" {
  // Empty / mismatched inputs return the identity result.
  let k = @src.mda_kabsch([], [])
  assert_true(@moonbit.assert_num_close(k.rmsd, 0.0, 1.0e-12))
  let (tx, ty, tz) = k.translation
  assert_true(@moonbit.assert_num_close(tx, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(ty, 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(tz, 0.0, 1.0e-12))
  let k2 = @src.mda_kabsch([(0.0, 0.0, 0.0)], [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0)])
  assert_true(@moonbit.assert_num_close(k2.rmsd, 0.0, 1.0e-12))
}
