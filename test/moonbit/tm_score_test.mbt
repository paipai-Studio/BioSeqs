///|
/// Tests for TM-score (tm_score.mbt).

test "tm_score_d0" {
  // d0 = 1.24 * cbrt(L-15) - 1.8 for L >= 22
  assert_true(assert_close(@src.tm_score_d0(50), 2.256122224633851, 1.0e-9))
  assert_true(assert_close(@src.tm_score_d0(22), 0.572034666637762, 1.0e-9))
  // L < 22 → 0.5
  assert_true(assert_close(@src.tm_score_d0(21), 0.5, 1.0e-9))
  assert_true(assert_close(@src.tm_score_d0(10), 0.5, 1.0e-9))
}

test "tm_score_aligned_identical" {
  let a = [@src.QsPoint::new(0.0, 0.0, 0.0)]
  let b = [@src.QsPoint::new(0.0, 0.0, 0.0)]
  // 完全重合 → TM-score = 1.0
  assert_true(assert_close(@src.tm_score_aligned(a, b, 1), 1.0, 1.0e-9))
}

test "tm_score_translation_invariant" {
  // 平移不影响叠加后的 TM-score（最优叠加后完全重合）
  let c1 = [
    @src.QsPoint::new(0.0, 0.0, 0.0),
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let c2 = [
    @src.QsPoint::new(5.0, 5.0, 5.0),
    @src.QsPoint::new(6.0, 5.0, 5.0),
    @src.QsPoint::new(5.0, 6.0, 5.0),
    @src.QsPoint::new(5.0, 5.0, 6.0),
  ]
  assert_true(assert_close(@src.tm_score(c1, c2), 1.0, 1.0e-6))
}

test "tm_score_rotation_invariant" {
  // 旋转 90° 也不影响 TM-score（非共面点，最优叠加后重合）
  let c1 = [
    @src.QsPoint::new(0.0, 0.0, 0.0),
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  // 绕 z 轴转 90°：(x,y,z) → (-y,x,z)
  let c2 = [
    @src.QsPoint::new(0.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(-1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  assert_true(assert_close(@src.tm_score(c1, c2), 1.0, 1.0e-6))
}

test "tm_score_in_range" {
  // 点扰动 → 分数在 (0, 1] 之间（接近 1 但 < 1）
  let c1 = [
    @src.QsPoint::new(0.0, 0.0, 0.0),
    @src.QsPoint::new(2.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 2.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 2.0),
  ]
  let c2 = [
    @src.QsPoint::new(0.0, 0.0, 0.0),
    @src.QsPoint::new(2.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 2.0, 0.0),
    @src.QsPoint::new(0.5, 0.5, 2.0),
  ]
  let s = @src.tm_score(c1, c2)
  assert_true(s > 0.0 && s <= 1.0 && s < 1.0)
}
