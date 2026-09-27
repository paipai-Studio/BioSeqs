///|
/// Tests for vegan rarefaction + decostand.

test "vegan_rarefy_two_species" {
  // 经典 vegan::rarefy 结果：两物种各 5 个个体，抽 1 个 → 期望 1 个物种
  assert_true(assert_close(@src.vegan_rarefy([5.0, 5.0], 1), 1.0, 1.0e-9))
  // 抽满 → 观察全部 2 个物种
  assert_true(assert_close(@src.vegan_rarefy([5.0, 5.0], 10), 2.0, 1.0e-9))
}

test "vegan_rarefy_single_species" {
  // 只有一个物种有丰度 → 任意深度都只有 1 个物种
  assert_true(assert_close(@src.vegan_rarefy([10.0, 0.0], 5), 1.0, 1.0e-9))
}

test "vegan_rarefaction_curve" {
  let curve = @src.vegan_rarefaction_curve([4.0, 4.0, 4.0], step=2)
  // 深度 2, 4, 6, 8, 10, 12（full）
  let last = curve[curve.length() - 1]
  assert_eq(last.0, 12)
  // 最深 = 观察物种数 3
  assert_true(assert_close(last.1, 3.0, 1.0e-9))
  // 曲线单调递增
  let mut ok = true
  let mut i = 1
  while i < curve.length() {
    if curve[i].1 < curve[i - 1].1 - 1.0e-9 {
      ok = false
    }
    i = i + 1
  }
  assert_true(ok)
}

test "vegan_decostand_total_max" {
  let m = [[2.0, 2.0], [3.0, 7.0]]
  let total = @src.vegan_decostand(m, "total")
  assert_true(assert_close(total[0][0], 0.5, 1.0e-9))
  assert_true(assert_close(total[1][1], 0.7, 1.0e-9))
  let mx = @src.vegan_decostand([[2.0, 8.0]], "max")
  assert_true(assert_close(mx[0][1], 1.0, 1.0e-9))
  assert_true(assert_close(mx[0][0], 0.25, 1.0e-9))
}

test "vegan_decostand_pa_hellinger_normalize" {
  let pa = @src.vegan_decostand([[0.0, 3.0], [5.0, 0.0]], "pa")
  assert_eq(pa[0][0], 0.0)
  assert_eq(pa[0][1], 1.0)
  assert_eq(pa[1][0], 1.0)
  assert_eq(pa[1][1], 0.0)

  let hell = @src.vegan_decostand([[1.0, 3.0]], "hellinger")
  assert_true(assert_close(hell[0][0], 0.5, 1.0e-9))           // sqrt(1/4)
  assert_true(assert_close(hell[0][1], 0.8660254038, 1.0e-9))  // sqrt(3/4)

  let norm = @src.vegan_decostand([[3.0, 4.0]], "normalize")
  assert_true(assert_close(norm[0][0], 0.6, 1.0e-9))
  assert_true(assert_close(norm[0][1], 0.8, 1.0e-9))
}
