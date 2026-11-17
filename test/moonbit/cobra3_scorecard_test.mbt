// cobra3_scorecard_test.mbt — Memote 式评分卡测试

///|
/// 评分卡用玩具模型：A 摄取 → B → 生物量，B 另有一条可逆支路 → C；
/// R2 缺 GPR、R3 缺 GPR 且 C 是死端（只产不耗）。
fn cobra3_scorecard_toy() -> @src.Cobra3Model {
  let mut m = @src.cobra3_model_new("draft")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, "g1"),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("B", -1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R3", [("B", -1.0), ("C", 1.0)], -1000.0, 1000.0, ""),
  )
  m = @src.cobra3_model_set_objective(m, "R2", 1.0)
  m
}

///|
test "cobra3 scorecard metrics on toy model" {
  let card = @src.cobra3_scorecard(cobra3_scorecard_toy())
  assert_true(card.model_id == "draft")
  // 6 个指标，权重和为 1
  assert_eq(card.metrics.length(), 6)
  let mut wsum = 0.0
  for m in card.metrics {
    wsum += m.weight
    assert_true(m.score >= 0.0 && m.score <= 1.0)
  }
  assert_true((wsum - 1.0).abs() < 1.0e-12)
  // 按指标 ID 查各项
  let get = fn(id : String) -> @src.Cobra3ScorecardMetric {
    let mut found = card.metrics[0]
    for m in card.metrics {
      if m.id == id {
        found = m
      }
    }
    found
  }
  // 无孤立代谢物（A/B/C 都参与反应）
  let cov = get("metabolite_coverage")
  assert_true(cov.score > 0.999)
  assert_eq(cov.passed, 3)
  assert_eq(cov.total, 3)
  // C 只产不耗 → 1 个死端（3 个有参与代谢物中的 2 个健康）
  let de = get("deadend_free")
  assert_eq(de.passed, 2)
  assert_eq(de.total, 3)
  assert_true((de.score - 2.0 / 3.0).abs() < 1.0e-9)
  // 内部反应 = R1/R2/R3；R2、R3 缺 GPR → 1/3
  let gpr = get("gpr_completeness")
  assert_eq(gpr.passed, 1)
  assert_eq(gpr.total, 3)
  // blocked：C 是死端（无出口）→ R3 不可载流，2/3 可载流
  let blk = get("blocked_free")
  assert_eq(blk.passed, 2)
  assert_eq(blk.total, 3)
  assert_true((blk.score - 2.0 / 3.0).abs() < 1.0e-9)
  // 总分与分档：缺 GPR + 死端 + blocked 拉低 → D 档
  assert_true(card.total_score > 0.4 && card.total_score < 0.6)
  assert_true(card.grade == "D")
}

///|
test "cobra3 scorecard penalizes missing gpr and deadend" {
  let good = @src.cobra3_scorecard(cobra3_scorecard_toy())
  // 参考（精修）模型：补全 GPR、加 D 消耗反应消除死端
  let mut refm = @src.cobra3_model_new("refined")
  refm = @src.cobra3_model_add_reaction(
    refm,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  refm = @src.cobra3_model_add_reaction(
    refm,
    cobra3_make_rxn("R1", [("A", -1.0), ("B", 1.0)], 0.0, 1000.0, "g1"),
  )
  refm = @src.cobra3_model_add_reaction(
    refm,
    cobra3_make_rxn("R2", [("B", -1.0)], 0.0, 1000.0, "g2"),
  )
  refm = @src.cobra3_model_add_reaction(
    refm,
    cobra3_make_rxn("R3", [("B", -1.0), ("C", 1.0)], -1000.0, 1000.0, "g3"),
  )
  refm = @src.cobra3_model_add_reaction(
    refm,
    cobra3_make_rxn("R4", [("C", -1.0), ("D", 1.0)], 0.0, 1000.0, "g4"),
  )
  refm = @src.cobra3_model_add_reaction(
    refm,
    cobra3_make_rxn("EX_D", [("D", -1.0)], -1000.0, 1000.0, ""),
  )
  refm = @src.cobra3_model_set_objective(refm, "R2", 1.0)
  let ref_card = @src.cobra3_scorecard(refm)
  // 参考模型无死端、GPR 完整
  let find = fn(card : @src.Cobra3Scorecard, id : String) -> Double {
    let mut s = 0.0
    for m in card.metrics {
      if m.id == id {
        s = m.score
      }
    }
    s
  }
  assert_true(find(ref_card, "deadend_free") > 0.999)
  assert_true(find(ref_card, "gpr_completeness") > 0.999)
  // 基准比较：草稿模型在 gpr/deadend 上更差，总分更低
  let cmp = @src.cobra3_scorecard_compare(good, ref_card)
  assert_true(cmp.total_score < cmp.reference_score)
  assert_true(cmp.delta < 0.0)
  assert_eq(cmp.diffs.length(), 6)
  // delta 升序（最差在前）→ 前两名是草稿模型缺失的两项
  assert_true(cmp.diffs[0].delta < 0.0 && cmp.diffs[1].delta < 0.0)
  assert_true(!cmp.diffs[0].improved)
  for d in cmp.diffs {
    assert_true((d.delta - (d.score - d.reference)).abs() < 1.0e-12)
  }
}

///|
test "cobra3 scorecard blocked and empty model" {
  // 口令反应：R1 的产物 D 无任何去路 → R1 被 FVA 判 blocked
  let mut m = @src.cobra3_model_new("blocked")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("A", -1.0), ("D", 1.0)], 0.0, 1000.0, "g1"),
  )
  m = @src.cobra3_model_set_objective(m, "R1", 1.0)
  let card = @src.cobra3_scorecard(m)
  let mut blk_score = 0.0
  let mut gpr_score = 0.0
  for x in card.metrics {
    if x.id == "blocked_free" {
      blk_score = x.score
      assert_eq(x.passed, 0)
      assert_eq(x.total, 1)
    }
    if x.id == "gpr_completeness" {
      gpr_score = x.score
    }
  }
  assert_true(blk_score < 1.0e-9)
  assert_true(gpr_score > 0.999)
  // 空模型：全部指标无统计对象 → 记 0 分、档位 E
  let empty = @src.cobra3_scorecard(@src.cobra3_model_new("empty"))
  assert_eq(empty.metrics.length(), 6)
  for x in empty.metrics {
    assert_true(x.score < 1.0e-12)
    assert_eq(x.total, 0)
  }
  assert_true(empty.total_score < 1.0e-12)
  assert_true(empty.grade == "E")
}
