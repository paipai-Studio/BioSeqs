// cobra3_reporter_test.mbt — Reporter metabolites（Patil & Nielsen 2005）测试
// ===========================================================================

///|
fn cobra3_reporter_model() -> @src.Cobra3Model {
  // 线性通路 A→B→C→D→E，中间穿插无 GPR 的交换/辅路：
  // 基因池 = {g1,g2,g3,g4,g6}（N=5），显著基因 {g1(3.5), g2(3.0)}（K=2）
  let mut m = @src.cobra3_model_new("reporter")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("EX_A", [("A_e", -1.0)], -10.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R_hot1", [("A_e", -1.0), ("B_e", 1.0)], 0.0, 1000.0, "g1"),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn(
      "R_hot2",
      [("B_e", -1.0), ("C_e", 1.0)],
      0.0,
      1000.0,
      "g1 and g2",
    ),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn(
      "R_cold",
      [("C_e", -1.0), ("D_e", 1.0)],
      0.0,
      1000.0,
      "g3 or g4",
    ),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn(
      "R_mid",
      [("D_e", -1.0), ("E_e", 1.0)],
      0.0,
      1000.0,
      "g1 and g6",
    ),
  )
  // 无 GPR 反应：X_e / Y_e 不该出现在输出中
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R_x", [("X_e", -1.0), ("Y_e", 1.0)], 0.0, 1000.0, ""),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R_x2", [("Y_e", -1.0)], 0.0, 1000.0, ""),
  )
  m
}

///|
fn cobra3_reporter_z() -> Map[String, Double] {
  Map([
    ("g1", 3.5),
    ("g2", 3.0),
    ("g3", 0.2),
    ("g4", 0.1),
    ("g5", 0.3),
    ("g6", 0.5),
  ])
}

///|
test "cobra3 reporter metabolites rank and statistics" {
  let reps = @src.cobra3_reporter_metabolites(
    cobra3_reporter_model(),
    cobra3_reporter_z(),
  )
  // 无 GPR 的 EX_A / R_x / R_x2 及其代谢物被排除，剩 A_e..E_e 共 5 个
  assert_eq(reps.length(), 5)
  // 输出按 p 升序：B_e(0.1) < A_e(0.4) < C_e(0.6) < E_e(0.7) < D_e(1.0)
  for i in 1..<reps.length() {
    assert_true(reps[i].p_value >= reps[i - 1].p_value)
  }
  let b = reps[0]
  assert_eq(b.metabolite_id, "B_e")
  assert_true((b.p_value - 0.1).abs() < 1.0e-9)
  // 反应分数：g1=3.5 支配 R_hot1/R_hot2
  assert_true((b.reporter_score - 3.5).abs() < 1.0e-12)
  assert_eq(b.n_reactions, 2)
  assert_eq(b.n_genes, 2)
  assert_eq(b.n_sig_genes, 2)
  // BH(m=5)：B_e 0.1*5 = 0.5（后续累计最小值不更低）
  assert_true((b.p_adj - 0.5).abs() < 1.0e-9)
  // A_e 只有 R_hot1 一个打分邻居
  let a = reps[1]
  assert_eq(a.metabolite_id, "A_e")
  assert_true((a.p_value - 0.4).abs() < 1.0e-9)
  assert_eq(a.n_reactions, 1)
  // C_e：邻域 4 基因、2 显著 → 3/5
  assert_eq(reps[2].metabolite_id, "C_e")
  assert_true((reps[2].p_value - 0.6).abs() < 1.0e-9)
  // D_e：邻域 {g3,g4,g1,g6} 1/4 显著，小池尾部饱和到 1.0
  let d = reps[4]
  assert_eq(d.metabolite_id, "D_e")
  assert_true(d.p_value > 0.99)
  // 反应分数取 |z| 最极端者：D_e 的 R_mid(3.5) 压过 R_cold(0.3)
  assert_true((d.reporter_score - 3.5).abs() < 1.0e-12)
  // p_adj ≥ p_value
  for r in reps {
    assert_true(r.p_adj >= r.p_value - 1.0e-12)
  }
}

///|
test "cobra3 reporter signed extreme and complex gpr" {
  // "g1 and g3"（复合体）与 "g8"（负变化 -4.0）的反应分数取 |z| 极端者的符号
  let mut m = @src.cobra3_model_new("reporter-signed")
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R1", [("M", -1.0), ("N", 1.0)], 0.0, 1000.0, "g1 and g3"),
  )
  m = @src.cobra3_model_add_reaction(
    m,
    cobra3_make_rxn("R2", [("N", -1.0), ("P", 1.0)], 0.0, 1000.0, "g8"),
  )
  let zs : Map[String, Double] = Map([("g1", 3.5), ("g3", 0.2), ("g8", -4.0)])
  let reps = @src.cobra3_reporter_metabolites(m, zs)
  assert_eq(reps.length(), 3)
  for r in reps {
    match r.metabolite_id {
      "M" => assert_true((r.reporter_score - 3.5).abs() < 1.0e-12)
      "N" => assert_true((r.reporter_score - -4.0).abs() < 1.0e-12)
      "P" => assert_true((r.reporter_score - -4.0).abs() < 1.0e-12)
      _ => assert_true(false)
    }
  }
}

///|
test "cobra3 reporter edge cases" {
  let model = cobra3_reporter_model()
  // 空 z 表 → 无打分反应 → 空输出
  assert_eq(@src.cobra3_reporter_metabolites(model, Map([])).length(), 0)
  // z 表基因与模型不交集 → 空输出
  let alien : Map[String, Double] = Map([("zz", 5.0)])
  assert_eq(@src.cobra3_reporter_metabolites(model, alien).length(), 0)
  // 提高 cutoff 到 4.0 → 全模型无显著基因，分数照常、p 全为 1
  let reps = @src.cobra3_reporter_metabolites(
    model,
    cobra3_reporter_z(),
    z_cutoff=4.0,
  )
  assert_eq(reps.length(), 5)
  for r in reps {
    assert_true(r.p_value > 0.99)
    assert_eq(r.n_sig_genes, 0)
  }
}
