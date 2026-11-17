// cobra3_robustness_test.mbt — 鲁棒性分析（单反应通量扫描）测试
// ===========================================================================

///|
fn cobra3_branch_two() -> @src.Cobra3Model {
  // 上游 EX_A 供给 10，下游两条支路各限 5，目标只走 BIO1：
  // 通量曲线在 [10,5] 摄入段有平台（支路容量饱和后生长不再增加）。
  let mut m = @src.cobra3_model_new("robust")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 5.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("A", -1.0), ("C", 1.0)], 0.0, 5.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO1", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO2", [("C", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO1", 1.0)
  m
}

///|
test "cobra3 robustness saturating branch has plateau" {
  let m = cobra3_branch_two()
  let curve = @src.cobra3_robustness(m, "EX_A", points=21)
  assert_eq(curve.length(), 21)
  // 摄入区间 [-10, 0]：两端点精确
  assert_true((curve[0].target_flux - -10.0).abs() < 1.0e-6)
  assert_true((curve[20].target_flux - 0.0).abs() < 1.0e-6)
  // 生长 = min(摄入, 5)：端点与平台段
  assert_true((curve[0].objective_max - 5.0).abs() < 1.0e-6)
  assert_true((curve[10].objective_max - 5.0).abs() < 1.0e-6)
  assert_true((curve[15].objective_max - 2.5).abs() < 1.0e-6)
  assert_true((curve[20].objective_max - 0.0).abs() < 1.0e-6)
  // 曲线单调不增（摄入越高、剩余生长空间越小）
  for i in 1..<curve.length() {
    assert_true(curve[i].objective_max <= curve[i - 1].objective_max + 1.0e-9)
  }
}

///|
test "cobra3 robustness forced chain is single valued" {
  let m = cobra3_toy()
  let curve = @src.cobra3_robustness(m, "EX_A", points=11)
  assert_eq(curve.length(), 11)
  // 线性链强制 v_R2 = -v_EX_A：max 与 min 重合，生长 = -通量
  for p in curve {
    assert_true((p.objective_max - p.objective_min).abs() < 1.0e-6)
    assert_true((p.objective_max + p.target_flux).abs() < 1.0e-6)
  }
  assert_true((curve[0].target_flux - -10.0).abs() < 1.0e-6)
  assert_true((curve[0].objective_max - 10.0).abs() < 1.0e-6)
  assert_true((curve[5].objective_max - 5.0).abs() < 1.0e-6)
  assert_true((curve[10].objective_max - 0.0).abs() < 1.0e-6)
}

///|
test "cobra3 robustness edge cases" {
  let m = cobra3_toy()
  // 未知反应 ID → 空曲线
  assert_eq(@src.cobra3_robustness(m, "NOPE").length(), 0)
  // points <= 0 → 空曲线
  assert_eq(@src.cobra3_robustness(m, "EX_A", points=0).length(), 0)
  // 单点扫描 → 只取范围下端
  let one = @src.cobra3_robustness(m, "EX_A", points=1)
  assert_eq(one.length(), 1)
  assert_true((one[0].target_flux - -10.0).abs() < 1.0e-6)
}
