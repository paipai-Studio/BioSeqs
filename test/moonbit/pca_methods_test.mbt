///|
/// Tests for the pcaMethods PPCA + BPCA port. Expected values are the output of
/// an independent numpy transcription of the pcaMethods equations on the same
/// synthetic matrix. Eigenvector signs are arbitrary, so loadings/scores are
/// compared by magnitude; scalars (tau, ss, imputed values) are compared
/// exactly.

///|
fn pm_nan() -> Double {
  0.0 / 0.0
}

///|
fn pm_data() -> Array[Array[Double]] {
  [
    [0.50792157, 0.24200472, -0.14188786, 0.31972680, 0.01227674],
    [-0.23431834, 0.46850265, pm_nan(), -0.50059657, -0.75143204],
    [0.75932131, 0.02362821, -0.69841887, 0.70615025, 0.51956590],
    [-1.05498977, -0.72517763, 0.14259525, -0.41810797, 0.41062713],
    [0.45841031, 0.74960856, 0.32468309, -0.04976495, -0.64565736],
    [pm_nan(), -0.93862724, -0.35686552, 0.04848338, 0.94033076],
    [1.11964842, 0.37134488, -0.66788357, 0.70307611, 0.17581264],
    [-0.86030953, 0.00690872, 0.80482140, -0.77911980, -0.51281908],
  ]
}

///|
fn pm_close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() < tol
}

///|
/// Sign-agnostic comparison (compares magnitudes).
fn pm_close_mag(a : Double, b : Double, tol : Double) -> Bool {
  (a.abs() - b.abs()).abs() < tol
}

///|
test "ppca_matches_reference" {
  let r = @src.ppca(matrix=pm_data(), n_pcs=2, max_iterations=2000, threshold=1.0e-5)
  // noise variance is invariant
  assert_true(pm_close(r.sigma2, 0.0009643751, 1.0e-6))
  // loadings (D×k) — magnitudes match the numpy reference
  assert_true(pm_close_mag(r.loadings[0][0], 0.708898, 1.0e-3))
  assert_true(pm_close_mag(r.loadings[0][1], 0.292039, 1.0e-3))
  assert_true(pm_close_mag(r.loadings[2][0], 0.435174, 1.0e-3))
  assert_true(pm_close_mag(r.loadings[4][1], 0.642863, 1.0e-3))
  // scores first row
  assert_true(pm_close_mag(r.scores[0][0], 0.629881, 1.0e-3))
  assert_true(pm_close_mag(r.scores[0][1], 0.18862, 1.0e-3))
}

///|
test "bpca_matches_reference" {
  let r = @src.bpca(matrix=pm_data(), n_pcs=2, max_steps=2000, threshold=1.0e-4)
  assert_true(pm_close(r.tau, 197.7274, 1.0e-2))
  // loadings (PA, D×k) — magnitudes
  assert_true(pm_close_mag(r.loadings[0][0], 0.698969, 1.0e-2))
  assert_true(pm_close_mag(r.loadings[0][1], 0.299191, 1.0e-2))
  assert_true(pm_close_mag(r.loadings[2][0], 0.48304, 1.0e-2))
  assert_true(pm_close_mag(r.loadings[4][1], 0.517957, 1.0e-2))
  // imputed values (exact — signs are meaningful here)
  assert_true(pm_close(r.imputed[1][2], 0.598401, 1.0e-3))
  assert_true(pm_close(r.imputed[5][0], -0.578966, 1.0e-3))
}

///|
test "bpca_imputation_is_reasonable" {
  // held-out truth values: [1,2]=0.664, [5,0]=-0.769. The imputation should be
  // in the right ballpark (within ~0.3).
  let r = @src.bpca(matrix=pm_data(), n_pcs=2, max_steps=2000, threshold=1.0e-4)
  assert_true(pm_close(r.imputed[1][2], 0.664209, 0.3))
  assert_true(pm_close(r.imputed[5][0], -0.768631, 0.3))
}

///|
test "ppca_bpca_shapes" {
  let r = @src.ppca(matrix=pm_data(), n_pcs=2, max_iterations=2000, threshold=1.0e-5)
  assert_eq(r.scores.length(), 8)
  assert_eq(r.scores[0].length(), 2)
  assert_eq(r.loadings.length(), 5)
  assert_eq(r.loadings[0].length(), 2)
  assert_eq(r.r2cum.length(), 2)
}