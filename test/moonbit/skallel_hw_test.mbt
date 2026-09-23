///|
/// Tests for scikit-allel Hardy–Weinberg statistics (array granularity).

///|
test "skallel_hw" {
  let g = [
    [[0, 0], [0, 0], [0, 0]],
    [[0, 0], [0, 1], [1, 1]],
    [[0, 0], [1, 1], [2, 2]],
    [[1, 1], [1, 2], [-1, -1]],
  ]
  // observed heterozygosity: 0, 1/3, 0, 0.5
  let ho = @src.skallel_heterozygosity_observed(g, -1.0)
  assert_true(assert_close(ho[0], 0.0, 1.0e-9))
  assert_true(assert_close(ho[1], 1.0 / 3.0, 1.0e-9))
  assert_true(assert_close(ho[2], 0.0, 1.0e-9))
  assert_true(assert_close(ho[3], 0.5, 1.0e-9))

  // expected heterozygosity: 0, 0.5, 2/3, 0.375
  let af = [
    [1.0],
    [0.5, 0.5],
    [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
    [0.75, 0.25],
  ]
  let he = @src.skallel_heterozygosity_expected(af, 2, -1.0)
  assert_true(assert_close(he[0], 0.0, 1.0e-9))
  assert_true(assert_close(he[1], 0.5, 1.0e-9))
  assert_true(assert_close(he[2], 2.0 / 3.0, 1.0e-9))
  assert_true(assert_close(he[3], 0.375, 1.0e-9))

  // inbreeding coefficient: nan(→ -999), 1/3, 1, -1/3
  let f = @src.skallel_inbreeding_coefficient(g, -999.0)
  assert_eq(f[0], -999.0)
  assert_true(assert_close(f[1], 1.0 / 3.0, 1.0e-9))
  assert_true(assert_close(f[2], 1.0, 1.0e-9))
  assert_true(assert_close(f[3], -1.0 / 3.0, 1.0e-9))
}