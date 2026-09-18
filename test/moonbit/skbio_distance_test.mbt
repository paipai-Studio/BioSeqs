///|
/// Tests for the scikit-bio ANOSIM + PERMANOVA port. Expected statistics come
/// from scikit-bio 0.7.3 (R=0.25 / F=2.0 for 2 groups, R=0.75 / F=8.0 for 3
/// groups). The permutation p-value is not bit-identical to numpy's RNG, so it
/// is only range-checked.

///|
test "skbio_anosim_statistic_2groups" {
  let dm = [[0.0, 1.0, 1.0, 4.0], [1.0, 0.0, 3.0, 2.0], [1.0, 3.0, 0.0, 3.0], [4.0, 2.0, 3.0, 0.0]]
  let g = [0, 0, 1, 1]
  let r = @src.skbio_anosim(dm, g, 99, 42)
  assert_true(assert_close(r.statistic, 0.25, 1.0e-9))
  assert_eq(r.n, 4)
  assert_eq(r.num_groups, 2)
  assert_true(r.p_value > 0.0 && r.p_value <= 1.0)
}

///|
test "skbio_permanova_statistic_2groups" {
  let dm = [[0.0, 1.0, 1.0, 4.0], [1.0, 0.0, 3.0, 2.0], [1.0, 3.0, 0.0, 3.0], [4.0, 2.0, 3.0, 0.0]]
  let g = [0, 0, 1, 1]
  let r = @src.skbio_permanova(dm, g, 99, 42)
  assert_true(assert_close(r.statistic, 2.0, 1.0e-9))
  assert_eq(r.n, 4)
  assert_eq(r.num_groups, 2)
  assert_true(r.p_value > 0.0 && r.p_value <= 1.0)
}

///|
test "skbio_anosim_permanova_3groups" {
  let dm = [
    [0.0, 1.0, 2.0, 3.0, 2.0, 1.0],
    [1.0, 0.0, 1.0, 2.0, 3.0, 2.0],
    [2.0, 1.0, 0.0, 1.0, 2.0, 3.0],
    [3.0, 2.0, 1.0, 0.0, 1.0, 2.0],
    [2.0, 3.0, 2.0, 1.0, 0.0, 1.0],
    [1.0, 2.0, 3.0, 2.0, 1.0, 0.0],
  ]
  let g = [0, 0, 1, 1, 2, 2]
  let a = @src.skbio_anosim(dm, g, 99, 42)
  assert_true(assert_close(a.statistic, 0.75, 1.0e-9))
  assert_eq(a.num_groups, 3)
  let p = @src.skbio_permanova(dm, g, 99, 42)
  assert_true(assert_close(p.statistic, 8.0, 1.0e-9))
  assert_eq(p.num_groups, 3)
  assert_true(p.p_value > 0.0 && p.p_value <= 1.0)
}