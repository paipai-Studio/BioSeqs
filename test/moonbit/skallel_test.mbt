///|
/// Tests for the scikit-allel port. Expected values are taken verbatim from
/// the scikit-allel docstring examples for `hudson_fst` and `rogers_huff_r`.

///|
test "skallel_hudson_fst" {
  let ac1 = [[4, 0], [2, 2], [4, 0], [1, 2, 1], [2, 2]]
  let ac2 = [[0, 4], [2, 2], [4, 0], [0, 2, 2], [1, 1]]
  let (num, den) = @src.skallel_hudson_fst(ac1, ac2)
  let exp_num = [1.0, -0.16666667, 0.0, -0.125, -0.33333333]
  let exp_den = [1.0, 0.5, 0.0, 0.625, 0.5]
  assert_eq(num.length(), 5)
  for i = 0; i < 5; i = i + 1 {
    assert_true(assert_close(num[i], exp_num[i], 1.0e-5))
    assert_true(assert_close(den[i], exp_den[i], 1.0e-5))
  }
  assert_true(assert_close(@src.skallel_hudson_fst_block(ac1, ac2), 0.1428571428571429, 1.0e-5))
}

///|
test "skallel_rogers_huff_r" {
  let gn = [[0, 2, 0], [0, 2, 0], [2, 0, 2], [0, 1, -1]]
  let r = @src.skallel_rogers_huff_r(gn)
  let exp = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
  assert_eq(r.length(), 6)
  for i = 0; i < 6; i = i + 1 {
    assert_true(assert_close(r[i], exp[i], 1.0e-5))
  }
}

///|
test "skallel_pca" {
  let gn = [[0, 2], [0, 2], [2, 0], [0, 1]]
  let pca = @src.skallel_pca(gn, 2)
  assert_eq(pca.coords.length(), 2)
  assert_eq(pca.coords[0].length(), 2)
  assert_eq(pca.explained_variance_ratio.length(), 2)
  let s = pca.explained_variance_ratio[0] + pca.explained_variance_ratio[1]
  assert_true(assert_close(s, 1.0, 1.0e-6))
}

// ============================================================
// 单点 tajima_d / patterson_f4 / LD between
// ============================================================

test "skallel_tajima_d_single_point" {
  let ac = [[2, 2], [3, 1], [1, 3], [2, 2]]
  let r = @src.skallel_tajima_d(ac, 4)
  assert_eq(r.test_name, "Tajima's D")
  assert_true(r.statistic == r.statistic)  // 有限（非 NaN）
}

test "skallel_patterson_f4_identical" {
  // 四个群体完全相同时，f4 = Σ(p-p)(p-p) = 0
  let ac = [[2, 2], [3, 1], [1, 3]]
  let f = @src.skallel_patterson_f4(ac, ac, ac, ac)
  assert_true(assert_close(f, 0.0, 1.0e-12))
}

test "skallel_rogers_huff_r_between" {
  let gn = [[0, 2, 0], [0, 2, 0], [2, 0, 2]]
  let r = @src.skallel_rogers_huff_r_between(gn, gn)
  assert_eq(r.length(), 9)  // 3 × 3
  // 对角线（同一位点自身）r² = 1
  assert_true(assert_close(r[0], 1.0, 1.0e-6))
  assert_true(assert_close(r[4], 1.0, 1.0e-6))
  assert_true(assert_close(r[8], 1.0, 1.0e-6))
}
