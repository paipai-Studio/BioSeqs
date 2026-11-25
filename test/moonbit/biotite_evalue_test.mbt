///|
/// Tests for biotite_evalue.mbt — Karlin-Altschul statistics.
test "biotite_evalue_defaults_and_monotonicity" {
  let est = @src.BiotiteEValueEstimator::new()
  assert_true(est.lambda > 0.0)
  assert_true(est.k > 0.0)
  let low = est.evalue(40.0, 100, 100)
  let high = est.evalue(60.0, 100, 100)
  assert_true(high < low)
  assert_true(low > 0.0)
}

///|
test "biotite_evalue_pvalue_range" {
  let est = @src.BiotiteEValueEstimator::new()
  let p = est.pvalue(50.0, 100, 100)
  assert_true(p >= 0.0)
  assert_true(p <= 1.0)
}

///|
test "biotite_evalue_bit_score_and_inverse" {
  let est = @src.BiotiteEValueEstimator::new()
  let e = est.evalue(50.0, 120, 80)
  let recovered = est.score_from_evalue(e, 120, 80)
  assert_true(@moonbit.assert_num_close(recovered, 50.0, 1.0e-6))
  assert_true(est.bit_score(50.0) > 0.0)
}

///|
test "biotite_evalue_estimate" {
  let matrix : Array[Array[Double]] = [[4.0, -1.0], [-1.0, 4.0]]
  let freqs : Array[Double] = [0.5, 0.5]
  let est = @src.BiotiteEValueEstimator::estimate(matrix, freqs)
  assert_true(est.lambda > 0.0)
  assert_true(est.k > 0.0)
}
