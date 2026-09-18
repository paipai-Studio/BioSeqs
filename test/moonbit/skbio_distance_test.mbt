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
///|
test "skbio_permdisp_statistic" {
  let dm = [
    [0.0, 0.223607, 0.111803, 2.828427, 3.162278, 3.590265],
    [0.223607, 0.0, 0.304138, 2.617250, 3.008322, 3.373426],
    [0.111803, 0.304138, 0.0, 2.865746, 3.148412, 3.677295],
    [2.828427, 2.617250, 2.865746, 0.0, 1.414214, 1.920937],
    [3.162278, 3.008322, 3.148412, 1.414214, 0.0, 3.330165],
    [3.590265, 3.373426, 3.677295, 1.920937, 3.330165, 0.0],
  ]
  let g = [0, 0, 0, 1, 1, 1]
  let r = @src.skbio_permdisp(dm, g, 99, 42)
  // scikit-bio permdisp(test='centroid') F = 4.410991776452
  assert_true(assert_close(r.statistic, 4.410991776452, 1.0e-5))
  assert_eq(r.n, 6)
  assert_eq(r.num_groups, 2)
  assert_true(r.p_value > 0.0 && r.p_value <= 1.0)
}
