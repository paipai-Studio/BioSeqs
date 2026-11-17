// kegg_fba_test.mbt — KEGG → FBA 草稿模型重构测试

///|
test "kegg reaction parse" {
  // 双向可逆式
  match @src.kegg_reaction_parse("C00031 + C00118 <=> C00197 + C00022") {
    Some((subs, prods)) => {
      assert_eq(subs.length(), 2)
      assert_eq(prods.length(), 2)
      assert_true(subs[0].0 == "C00031" && subs[1].0 == "C00118")
      assert_true(prods[0].0 == "C00197" && prods[1].0 == "C00022")
    }
    None => fail("parse failed")
  }
  // 数字系数 + 不可逆箭头
  match @src.kegg_reaction_parse("2 C00022 -> C00200") {
    Some((l, r)) => {
      assert_true(l.length() == 1 && l[0].1 == 2.0)
      assert_true(r.length() == 1 && r[0].0 == "C00200")
    }
    None => fail("parse failed")
  }
  // cpd: 前缀剥离
  match @src.kegg_reaction_parse("cpd:C00031 => cpd:C00022") {
    Some((l, _r)) => assert_true(l[0].0 == "C00031")
    None => fail("parse failed")
  }
  // 无分隔符 → None；右侧为空 → None
  assert_true(@src.kegg_reaction_parse("C00031") is None)
  assert_true(@src.kegg_reaction_parse("C00031 -> ") is None)
}

///|
test "kegg fba reconstruct from kgml" {
  let pathway = @src.parse_kgml(@src.kgml_sample_pathway()).unwrap()
  // 1. 闭环不开 exchange：无底物来源 → objective 0
  let closed = @src.kegg_fba_reconstruct(pathway, objective_reaction="R01786")
  assert_eq(closed.n_metabolites, 4)
  assert_eq(closed.n_reactions, 2)
  assert_eq(closed.n_exchanges, 0)
  assert_true(closed.issues.length() == 0)
  let sol0 = @src.fba_solve(closed.model)
  assert_true(sol0.status == "optimal")
  assert_true(sol0.objective_value < 1.0e-9)
  // 计量矩阵方向：S[C00118][R01786] = -1, S[C00197][R01786] = +1
  let mets = closed.model.metabolites
  let idx = fn(id : String) -> Int32 {
    let mut out = -1
    for i in 0..<mets.length() {
      if mets[i] == id {
        out = i
      }
    }
    out
  }
  let ridx = fn(id : String) -> Int32 {
    let mut out = -1
    for j in 0..<closed.model.reactions.length() {
      if closed.model.reactions[j].id == id {
        out = j
      }
    }
    out
  }
  assert_true(closed.model.stoichiometry[idx("C00118")][ridx("R01786")] == -1.0)
  assert_true(closed.model.stoichiometry[idx("C00197")][ridx("R01786")] == 1.0)
  // 可逆性 → 边界
  assert_true(closed.model.reactions[ridx("R01070")].lower_bound == -1000.0)
  assert_true(closed.model.reactions[ridx("R01786")].lower_bound == 0.0)
  // entry 基因关联（ko:Kxxxxx）
  assert_eq(closed.model.reactions[ridx("R01786")].genes, ["K00844"])
  assert_eq(closed.model.reactions[ridx("R01070")].genes, ["K01810"])
  // 2. 开全部 exchange：边界供料 → 生长
  let open = @src.kegg_fba_reconstruct(
    pathway,
    open_exchanges=true,
    objective_reaction="R01786",
  )
  assert_eq(open.n_exchanges, 4)
  let sol = @src.fba_solve(open.model)
  assert_true(sol.status == "optimal")
  assert_true(sol.objective_value > 0.0)
  assert_eq(sol.fluxes.length(), open.model.reactions.length())
}

///|
test "kegg fba mass balance check" {
  let m = @src.fba_model_new(["A", "B", "C", "D"])
  let m1 = m.add_reaction(
    {
      id: "R_empty",
      name: "",
      lower_bound: 0.0,
      upper_bound: 1000.0,
      objective: 0.0,
      genes: [],
    },
    [],
  )
  let m2 = m1.add_reaction(
    {
      id: "R_one",
      name: "",
      lower_bound: 0.0,
      upper_bound: 1000.0,
      objective: 0.0,
      genes: [],
    },
    [("A", 1.0)],
  )
  let m3 = m2.add_reaction(
    {
      id: "R_bal",
      name: "",
      lower_bound: 0.0,
      upper_bound: 1000.0,
      objective: 0.0,
      genes: [],
    },
    [("B", -1.0), ("C", 1.0)],
  )
  let m4 = m3.add_reaction(
    {
      id: "R_unbal",
      name: "",
      lower_bound: 0.0,
      upper_bound: 1000.0,
      objective: 0.0,
      genes: [],
    },
    [("C", -1.0), ("D", 1.0)],
  )
  let formulas : Map[String, String] = Map([
    ("A", "C1H4"),
    ("B", "C2H6O"),
    ("C", "C2H6O"),
    ("D", "C3H6O"),
  ])
  let issues = @src.kegg_fba_check_mass_balance(m4, formulas~)
  // empty / one_sided / imbalanced 各一条；R_bal 平衡无报告
  assert_eq(issues.length(), 3)
  let mut k_empty = 0
  let mut k_one = 0
  let mut k_imb = 0
  for x in issues {
    if x.kind == "empty" {
      k_empty = k_empty + 1
      assert_true(x.reaction_id == "R_empty")
    }
    if x.kind == "one_sided" {
      k_one = k_one + 1
      assert_true(x.reaction_id == "R_one")
    }
    if x.kind == "imbalanced" {
      k_imb = k_imb + 1
      assert_true(x.reaction_id == "R_unbal")
      // C2H6O -> C3H6O：净增 1 个 C
      assert_true(x.detail.contains("C:1"))
    }
  }
  assert_eq(k_empty, 1)
  assert_eq(k_one, 1)
  assert_eq(k_imb, 1)
  // 不提供 formulas：只做拓扑检查（R_bal/R_unbal 元素检查跳过）
  let topo = @src.kegg_fba_check_mass_balance(m4)
  assert_eq(topo.length(), 2)
}

///|
test "kegg fba gapfill closure" {
  let pathway = @src.parse_kgml(@src.kgml_sample_pathway()).unwrap()
  // 草稿：不开 exchange，仅手动开底物摄取 → 产物无外排出口，仍不生长
  let draft_model = @src.kegg_fba_reconstruct(
      pathway,
      objective_reaction="R01786",
    ).model
  let draft = @src.kegg_fba_add_exchange(draft_model, "C00118")
  let sol0 = @src.fba_solve(draft)
  assert_true(sol0.status == "optimal")
  assert_true(sol0.objective_value < 1.0e-9)
  // universal 池：一条产物外排 exchange（用 add_exchange 构造，因
  // "C00197 -> " 右侧为空不可由 kegg_reaction_parse 解析）
  let universal = @src.kegg_fba_add_exchange(
    @src.fba_model_new([]),
    "C00197",
    lower_bound=0.0,
  )
  assert_eq(universal.reactions.length(), 1)
  // gapfill 闭环：补上外排反应后恢复生长
  let results = @src.cobra_gapfill(draft, universal, lower_bound=0.05)
  assert_eq(results.length(), 1)
  assert_true(results[0].success)
  assert_eq(results[0].added_reactions, ["EX_C00197"])
  assert_true(results[0].objective_value > 0.05)
}
