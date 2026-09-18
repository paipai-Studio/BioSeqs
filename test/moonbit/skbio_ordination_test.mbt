///|
/// Tests for the scikit-bio ordination + mantel port. Expected values come from
/// scikit-bio 0.7.3. CA score signs are arbitrary, so scores are compared by
/// magnitude; eigenvalues, proportions, centering, and the Mantel correlation
/// coefficient are compared exactly.

///|
fn skbio_ord_close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() < tol
}

///|
fn skbio_ord_close_mag(a : Double, b : Double, tol : Double) -> Bool {
  (a.abs() - b.abs()).abs() < tol
}

///|
test "skbio_center_distance_matrix" {
  let dm = [[0.0, 1.0, 2.0], [1.0, 0.0, 1.0], [2.0, 1.0, 0.0]]
  let c = @src.skbio_center_distance_matrix(dm)
  assert_true(skbio_ord_close(c[0][0], 1.0, 1.0e-9))
  assert_true(skbio_ord_close(c[0][1], 0.0, 1.0e-9))
  assert_true(skbio_ord_close(c[0][2], -1.0, 1.0e-9))
  assert_true(skbio_ord_close(c[1][0], 0.0, 1.0e-9))
  assert_true(skbio_ord_close(c[1][1], 0.0, 1.0e-9))
  assert_true(skbio_ord_close(c[2][0], -1.0, 1.0e-9))
  assert_true(skbio_ord_close(c[2][2], 1.0, 1.0e-9))
}

///|
test "skbio_ca_scaling1" {
  let x = [
    [1.0, 0.0, 2.0, 1.0],
    [0.0, 3.0, 1.0, 0.0],
    [2.0, 1.0, 0.0, 1.0],
  ]
  let r = @src.skbio_ca(x, 1)
  assert_eq(r.eigvals.length(), 2)
  assert_true(skbio_ord_close(r.eigvals[0], 0.470162, 1.0e-5))
  assert_true(skbio_ord_close(r.eigvals[1], 0.238172, 1.0e-5))
  assert_true(skbio_ord_close(r.proportion_explained[0], 0.663758, 1.0e-5))
  assert_true(skbio_ord_close(r.proportion_explained[1], 0.336242, 1.0e-5))
  // samples (magnitude)
  assert_true(skbio_ord_close_mag(r.sample_scores[0][0], 0.548885, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.sample_scores[0][1], 0.568969, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.sample_scores[1][0], 0.966748, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.sample_scores[2][1], 0.622809, 1.0e-4))
  // features (magnitude)
  assert_true(skbio_ord_close_mag(r.feature_scores[0][0], 0.981655, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.feature_scores[1][0], 1.319961, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.feature_scores[2][1], 1.667953, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.feature_scores[3][0], 1.028101, 1.0e-4))
}

///|
test "skbio_ca_scaling2" {
  let x = [
    [1.0, 0.0, 2.0, 1.0],
    [0.0, 3.0, 1.0, 0.0],
    [2.0, 1.0, 0.0, 1.0],
  ]
  let r = @src.skbio_ca(x, 2)
  assert_true(skbio_ord_close(r.eigvals[0], 0.470162, 1.0e-5))
  assert_true(skbio_ord_close(r.eigvals[1], 0.238172, 1.0e-5))
  assert_true(skbio_ord_close_mag(r.sample_scores[0][0], 0.800493, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.sample_scores[0][1], 1.165852, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.sample_scores[1][0], 1.409904, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.sample_scores[2][1], 1.276173, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.feature_scores[0][0], 0.673105, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.feature_scores[1][1], 0.236302, 1.0e-4))
  assert_true(skbio_ord_close_mag(r.feature_scores[2][1], 0.814008, 1.0e-4))
}

///|
test "skbio_mantel_correlation" {
  let x = [[0.0, 1.0, 2.0], [1.0, 0.0, 1.0], [2.0, 1.0, 0.0]]
  let y = [[0.0, 2.0, 7.0], [2.0, 0.0, 5.0], [7.0, 5.0, 0.0]]
  let r = @src.skbio_mantel(x, y, 999, 42)
  assert_true(skbio_ord_close(r.corr_coeff, 0.8029550685, 1.0e-7))
  assert_eq(r.n, 3)
  // p-value is permutation-based and must be a valid probability
  assert_true(r.p_value > 0.0 && r.p_value <= 1.0)
}