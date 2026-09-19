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

///|
test "skbio_bioenv" {
  // tie-free dataset: Spearman depends only on ranks, robust to floating-point
  let dm = [
    [0.0, 1.9797207674, 2.2744954552, 3.3592500938, 2.6075188868, 2.4045200192],
    [1.9797207674, 0.0, 3.3362897722, 3.9784043017, 2.6069437758, 2.7906919854],
    [2.2744954552, 3.3362897722, 0.0, 2.9737874206, 3.2048016382, 3.4569544240],
    [3.3592500938, 3.9784043017, 2.9737874206, 0.0, 2.7403664856, 2.9079243854],
    [2.6075188868, 2.6069437758, 3.2048016382, 2.7403664856, 0.0, 2.7291111269],
    [2.4045200192, 2.7906919854, 3.4569544240, 2.9079243854, 2.7291111269, 0.0],
  ]
  // columns: v1, v2, v3
  let data = [
    [-0.0325217049, 0.8843898674, -0.5836004327],
    [-0.1117019496, 0.1104641432, 0.0637817743],
    [-1.2250558264, 0.0761402304, 1.3588234217],
    [-1.5471446781, 0.8593826880, 0.1193540257],
    [-0.6414703941, 2.0004165463, 0.7622597121],
    [-1.1992889021, 0.0745162288, 0.5766895837],
  ]
  let r = @src.skbio_bioenv(dm, data)
  assert_eq(r.correlations.length(), 3)
  assert_true(assert_close(r.correlations[0], 0.132142857, 1.0e-6))
  assert_true(assert_close(r.correlations[1], 0.075, 1.0e-6))
  assert_true(assert_close(r.correlations[2], -0.339285714, 1.0e-6))
  // best subset: size1 = {v1} (col 0); size2 = {v1, v3} (cols 0, 2)
  assert_eq(r.best_subsets[0][0], 0)
  assert_eq(r.best_subsets[1][0], 0)
  assert_eq(r.best_subsets[1][1], 2)
}
