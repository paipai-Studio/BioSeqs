///|
/// Tests for mda_analysis_pca.mbt (PCA, TICA, collective variables).

///|
test "mda_pca_empty_input" {
  let r = @src.mda_pca([])
  assert_eq(r.n_components, 0)
  assert_eq(r.eigenvalues.length(), 0)
  assert_eq(r.explained_variance_ratio.length(), 0)
  assert_eq(r.mean_coords.length(), 0)
  assert_eq(r.projections.length(), 0)
}

///|
test "mda_pca_empty_atoms" {
  let r = @src.mda_pca([[]])
  assert_eq(r.n_components, 0)
  assert_eq(r.mean_coords.length(), 0)
}

///|
test "mda_pca_single_axis_variance" {
  // 1 atom, 2 frames: x moves 0 -> 2, y/z fixed.
  // Var(x) = 1, Var(y) = Var(z) = 0, mean = (1, 0, 0).
  let frames = [[(0.0, 0.0, 0.0)], [(2.0, 0.0, 0.0)]]
  let r = @src.mda_pca(frames)
  assert_eq(r.n_components, 3)
  assert_eq(r.mean_coords.length(), 1)
  assert_true(@moonbit.assert_num_close(r.mean_coords[0].0, 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.mean_coords[0].1, 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.eigenvalues[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.eigenvalues[1], 0.0, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(r.explained_variance_ratio[0], 1.0, 1.0e-9),
  )
  // Projections are the centered coordinates of the top components.
  assert_eq(r.projections.length(), 2)
  assert_eq(r.projections[0].length(), 3)
  assert_true(@moonbit.assert_num_close(r.projections[0][0], -1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.projections[1][0], 1.0, 1.0e-9))
}

///|
test "mda_pca_n_components_limit" {
  // 2 atoms x 3 frames; only atom0.x varies (0, 2, 4) -> Var = 8/3.
  let frames = [
    [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
    [(2.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
    [(4.0, 0.0, 0.0), (0.0, 0.0, 0.0)],
  ]
  let r = @src.mda_pca(frames, n_components=2)
  assert_eq(r.n_components, 2)
  assert_eq(r.eigenvalues.length(), 2)
  assert_true(@moonbit.assert_num_close(r.eigenvalues[0], 8.0 / 3.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r.eigenvalues[1], 0.0, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(r.explained_variance_ratio[0], 1.0, 1.0e-9),
  )
  assert_eq(r.projections.length(), 3)
  assert_eq(r.projections[0].length(), 2)
  // Frame 0: x centered = -2 for atom0.
  assert_true(@moonbit.assert_num_close(r.projections[0][0], -2.0, 1.0e-9))
}

///|
test "mda_tica_fallback_to_pca" {
  // n_frames (3) <= lag_time (5) -> falls back to plain PCA.
  let frames = [[(0.0, 0.0, 0.0)], [(2.0, 0.0, 0.0)], [(4.0, 0.0, 0.0)]]
  let r = @src.mda_tica(frames, lag_time=5)
  assert_eq(r.n_components, 3)
  assert_true(@moonbit.assert_num_close(r.eigenvalues[0], 8.0 / 3.0, 1.0e-9))
}

///|
test "mda_tica_lagged_variance" {
  // 1 atom, x alternates 0,1,0,1,0; mean = 0.4; lag = 2 -> count = 3.
  // vx = ((-0.4)(-0.4) + (0.6)(0.6) + (-0.4)(-0.4)) / 3 = 0.68 / 3.
  let frames = [
    [(0.0, 0.0, 0.0)],
    [(1.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0)],
    [(1.0, 0.0, 0.0)],
    [(0.0, 0.0, 0.0)],
  ]
  let r = @src.mda_tica(frames, lag_time=2, n_components=1)
  assert_eq(r.n_components, 1)
  assert_true(@moonbit.assert_num_close(r.eigenvalues[0], 0.68 / 3.0, 1.0e-9))
  assert_true(
    @moonbit.assert_num_close(r.explained_variance_ratio[0], 1.0, 1.0e-9),
  )
  assert_eq(r.projections.length(), 5)
  assert_eq(r.projections[0].length(), 1)
}

///|
test "mda_cv_rgyr_square" {
  // 4 atoms forming a unit square in the xy-plane: Rg = sqrt(0.5).
  let frames = [
    [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (1.0, 1.0, 0.0)],
  ]
  let rgyr = @src.mda_cv_rgyr(frames)
  assert_eq(rgyr.length(), 1)
  assert_true(@moonbit.assert_num_close(rgyr[0], 0.5.sqrt(), 1.0e-9))
}

///|
test "mda_cv_rgyr_empty_frame" {
  let rgyr = @src.mda_cv_rgyr([[]])
  assert_eq(rgyr.length(), 1)
  assert_true(@moonbit.assert_num_close(rgyr[0], 0.0, 1.0e-9))
}

///|
test "mda_cv_end_to_end_345_triangle" {
  // First atom at origin, last atom at (3, 4, 0) -> distance 5.
  let frames = [[(0.0, 0.0, 0.0), (3.0, 4.0, 0.0)], [(1.0, 1.0, 1.0)], []]
  let r = @src.mda_cv_end_to_end(frames)
  assert_eq(r.length(), 3)
  assert_true(@moonbit.assert_num_close(r[0], 5.0, 1.0e-9))
  // Fewer than 2 atoms -> 0.
  assert_true(@moonbit.assert_num_close(r[1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(r[2], 0.0, 1.0e-9))
}
