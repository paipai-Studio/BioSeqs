// cobra3_cfba_test.mbt — cFBA 蛋白质配额约束测试

///|
/// cFBA 测试模型：EX_A（摄取 10）→ R1（A→B，不可逆）→ BIO（目标）。
fn cobra3_cfba_model() -> @src.Cobra3Model {
  let mut m = @src.cobra3_model_new("cfba")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  m
}

///|
test "cobra3 cfba protein quota limits growth" {
  let m = cobra3_cfba_model()
  let mw : Map[String, Double] = Map([("R1", 1.0), ("R2", 1.0)])
  // 配额充足（> 底物上限 10）→ 底物受限，生长 10
  let rich = @src.cobra3_cfba(m, mw, protein_budget=20.0, k_eff=1.0)
  assert_true(@moonbit.assert_num_close(rich.growth, 10.0, 1.0e-6))
  assert_true(rich.enzyme_load <= 10.0 + 1.0e-6)
  // 配额绑定时生长 = 配额，且酶占用打满
  let binding = @src.cobra3_cfba(m, mw, protein_budget=5.0, k_eff=1.0)
  assert_true(@moonbit.assert_num_close(binding.growth, 5.0, 1.0e-6))
  assert_true(@moonbit.assert_num_close(binding.enzyme_load, 5.0, 1.0e-6))
  assert_true(binding.fluxes.length() == 4)
  // 零配额 → 生长 0
  let zero = @src.cobra3_cfba(m, mw, protein_budget=0.0, k_eff=1.0)
  assert_true(zero.growth < 1.0e-9)
  assert_true(zero.enzyme_load < 1.0e-9)
}

///|
test "cobra3 cfba unassigned reactions bypass quota" {
  let m = cobra3_cfba_model()
  // 只给 R1 记 MW：LP 走免费的 R2 旁路 → 配额为 1 仍生长 10
  let mw : Map[String, Double] = Map([("R1", 1.0)])
  let bypass = @src.cobra3_cfba(m, mw, protein_budget=1.0, k_eff=1.0)
  assert_true(@moonbit.assert_num_close(bypass.growth, 10.0, 1.0e-6))
  assert_true(bypass.enzyme_load < 1.0e-9)
  // 空 MW 映射 → 等价普通 FBA
  let empty : Map[String, Double] = Map([])
  let plain = @src.cobra3_cfba(m, empty, protein_budget=0.001, k_eff=1.0)
  assert_true(@moonbit.assert_num_close(plain.growth, 10.0, 1.0e-6))
}

///|
test "cobra3 cfba sweep and edge cases" {
  let m = cobra3_cfba_model()
  let mw : Map[String, Double] = Map([("R1", 1.0), ("R2", 1.0)])
  // 生长 vs 配额曲线：线性上升后底物受限进入平台
  let sweep = @src.cobra3_cfba_sweep(
    m,
    mw,
    points=3,
    p_min=0.5,
    p_max=20.0,
    k_eff=1.0,
  )
  assert_eq(sweep.length(), 3)
  assert_true(@moonbit.assert_num_close(sweep[0].protein, 0.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(sweep[0].growth, 0.5, 1.0e-6))
  assert_true(@moonbit.assert_num_close(sweep[1].growth, 10.0, 1.0e-6))
  assert_true(@moonbit.assert_num_close(sweep[2].growth, 10.0, 1.0e-6))
  // points < 1（含负数）→ 空（负 capacity 构造已移到边界判断之后）
  assert_eq(@src.cobra3_cfba_sweep(m, mw, points=0, k_eff=1.0).length(), 0)
  assert_eq(@src.cobra3_cfba_sweep(m, mw, points=-1, k_eff=1.0).length(), 0)
  // 映射中未知反应 ID 与非正 MW 被忽略
  let weird : Map[String, Double] = Map([
    ("R1", 1.0),
    ("R_missing", 5.0),
    ("R2", -1.0),
  ])
  // R2 未计费 → 走 R2 旁路，配额 1 仍生长 10
  let ignored = @src.cobra3_cfba(m, weird, protein_budget=1.0, k_eff=1.0)
  assert_true(@moonbit.assert_num_close(ignored.growth, 10.0, 1.0e-6))
  // k_eff 缩放酶成本：MW/k_eff 减半 → 同配额下生长翻倍
  let mws : Map[String, Double] = Map([("R1", 100.0), ("R2", 100.0)])
  let slow = @src.cobra3_cfba(m, mws, protein_budget=5.0, k_eff=100.0)
  assert_true(@moonbit.assert_num_close(slow.growth, 5.0, 1.0e-6))
  let fast = @src.cobra3_cfba(m, mws, protein_budget=5.0, k_eff=200.0)
  assert_true(@moonbit.assert_num_close(fast.growth, 10.0, 1.0e-6))
}

///|
test "cobra3 cfba rejects nonpositive k_eff" {
  let m = cobra3_cfba_model()
  let mw : Map[String, Double] = Map([("R1", 1.0)])
  // k_eff = 0：非正比催化速率无意义 → 显式报错而非静默丢弃约束
  let raised0 = try {
    let _ = @src.cobra3_cfba(m, mw, protein_budget=1.0, k_eff=0.0)
    false
  } catch {
    _ => true
  }
  assert_true(raised0)
  let raised_neg = try {
    let _ = @src.cobra3_cfba(m, mw, protein_budget=1.0, k_eff=-1.0)
    false
  } catch {
    _ => true
  }
  assert_true(raised_neg)
  // sweep 同样透传报错
  let raised_sweep = try {
    let _ = @src.cobra3_cfba_sweep(m, mw, points=2, k_eff=0.0)
    false
  } catch {
    _ => true
  }
  assert_true(raised_sweep)
  // points<1 提前返回的分支同样校验 k_eff（入口统一）
  let raised_sweep0 = try {
    let _ = @src.cobra3_cfba_sweep(m, mw, points=0, k_eff=0.0)
    false
  } catch {
    _ => true
  }
  assert_true(raised_sweep0)
}
