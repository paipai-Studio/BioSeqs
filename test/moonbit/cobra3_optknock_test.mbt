// cobra3_optknock_test.mbt — OptKnock 菌株设计（bilevel 枚举式）测试
// ===========================================================================

///|
fn cobra3_optknock_model() -> @src.Cobra3Model {
  // 呼吸支路（A→B，无副产物）与发酵支路（A→B+P）并行：
  // 敲掉呼吸基因后 B 只能经发酵产生 → 产物 P 成为生长耦合。
  let mut m = @src.cobra3_model_new("optknock")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R_resp", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, "G_resp"),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn(
      "R_ferm",
      [("A", -1.0), ("B", 1.0), ("P", 1.0)],
      0.0,
      1000.0,
      "G_ferm",
    ),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("BIO", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_P", [("P", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "BIO", 1.0)
  m
}

///|
test "cobra3 optknock couples product to growth" {
  let reps = @src.cobra3_optknock_genes(
    cobra3_optknock_model(),
    "EX_P",
    ["G_resp", "G_ferm"],
    n_knockouts=1,
    min_growth_fraction=0.1,
  )
  assert_eq(reps.length(), 2)
  // 按产量降序：敲 G_resp（强制走发酵）排第一
  let best = reps[0]
  assert_eq(best.knockouts.length(), 1)
  assert_eq(best.knockouts[0], "G_resp")
  assert_true((best.growth - 10.0).abs() < 1.0e-6)
  assert_true((best.target_max - 10.0).abs() < 1.0e-6)
  // 生长下限 1.0 → 产物通量下界同值：生长耦合
  assert_true((best.target_min - 1.0).abs() < 1.0e-6)
  assert_true(best.coupled)
  // 敲 G_ferm：P 无来源，产量为 0、不耦合
  let alt = reps[1]
  assert_eq(alt.knockouts[0], "G_ferm")
  assert_true((alt.target_max - 0.0).abs() < 1.0e-9)
  assert_true((alt.target_min - 0.0).abs() < 1.0e-9)
  assert_true(!alt.coupled)
  assert_true((alt.growth - 10.0).abs() < 1.0e-6)
}

///|
test "cobra3 optknock double knockout skips lethal combo" {
  let reps = @src.cobra3_optknock_genes(
    cobra3_optknock_model(),
    "EX_P",
    ["G_resp", "G_ferm"],
    n_knockouts=2,
    min_growth_fraction=0.1,
  )
  // 双敲 {G_resp, G_ferm}：B 无来源 → 生长 0 < 下限 → 致死组合被剔除，
  // 只剩两个单敲结果
  assert_eq(reps.length(), 2)
  assert_eq(reps[0].knockouts[0], "G_resp")
  assert_eq(reps[1].knockouts[0], "G_ferm")
}

///|
test "cobra3 optknock strict growth floor" {
  // min_growth_fraction=1.0：要求敲除后生长不低于基线 10
  let reps = @src.cobra3_optknock_genes(
    cobra3_optknock_model(),
    "EX_P",
    ["G_resp", "G_ferm"],
    n_knockouts=1,
    min_growth_fraction=1.0,
  )
  assert_eq(reps.length(), 2)
  let best = reps[0]
  assert_eq(best.knockouts[0], "G_resp")
  // 生长必须 = 10 → 发酵通量 = 10 → 产物被完全锁定
  assert_true((best.target_min - 10.0).abs() < 1.0e-6)
  assert_true((best.target_max - 10.0).abs() < 1.0e-6)
  assert_true(best.coupled)
}

///|
test "cobra3 optknock edge cases" {
  let model = cobra3_optknock_model()
  // 未知目标反应
  assert_eq(@src.cobra3_optknock_genes(model, "NOPE", ["G_resp"]).length(), 0)
  // 空候选池
  assert_eq(@src.cobra3_optknock_genes(model, "EX_P", []).length(), 0)
  // n_knockouts = 0
  assert_eq(
    @src.cobra3_optknock_genes(model, "EX_P", ["G_resp"], n_knockouts=0).length(),
    0,
  )
  // 无目标（生物量）反应的模型 → 基线不生长
  let mut no_obj = @src.cobra3_model_new("noobj")
  no_obj = @src.cobra3_model_add_reaction(
    no_obj,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  no_obj = @src.cobra3_model_add_reaction(
    no_obj,
    cobra3_make_rxn("R1", [("A", -1.0)], 0.0, 1000.0, "G1"),
  )
  assert_eq(@src.cobra3_optknock_genes(no_obj, "EX_A", ["G1"]).length(), 0)
  // 与模型 GPR 无交集的候选基因：敲除无效（未失活任何反应）→ 被跳过
  let reps = @src.cobra3_optknock_genes(
    model,
    "EX_P",
    ["G_resp", "G_unrelated"],
    n_knockouts=1,
  )
  assert_eq(reps.length(), 1)
  assert_eq(reps[0].knockouts[0], "G_resp")
}
