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

///|
/// 构造最小 KGML：gene entry name 含空格分隔的两个 KO ID + 一条反应.
fn kegg_fba_multi_gene_xml(gene_name : String) -> String {
  let xml = "<?xml version=\"1.0\"?>\n"
  let s = xml +
    "<pathway name=\"path:tst00001\" org=\"tst\" number=\"00001\" title=\"Test\">\n" +
    "<entry id=\"1\" name=\"" +
    gene_name +
    "\" type=\"gene\" reaction=\"rn:R00001\" link=\"https://example.com/K\">"
  let s2 = s +
    "<graphics name=\"K\" fgcolor=\"#000000\" bgcolor=\"#BFFFBF\" type=\"rect\" x=\"100\" y=\"200\" width=\"46\" height=\"17\"/>"
  let s3 = s2 + "</entry>\n"
  let s4 = s3 +
    "<entry id=\"2\" name=\"cpd:C00001\" type=\"compound\" link=\"https://example.com/C1\">"
  let s5 = s4 +
    "<graphics name=\"C00001\" fgcolor=\"#000000\" bgcolor=\"#FFFFFF\" type=\"circle\" x=\"150\" y=\"250\" width=\"8\" height=\"8\"/>"
  let s6 = s5 + "</entry>\n"
  let s7 = s6 +
    "<entry id=\"3\" name=\"cpd:C00002\" type=\"compound\" link=\"https://example.com/C2\">"
  let s8 = s7 +
    "<graphics name=\"C00002\" fgcolor=\"#000000\" bgcolor=\"#FFFFFF\" type=\"circle\" x=\"160\" y=\"250\" width=\"8\" height=\"8\"/>"
  let s9 = s8 + "</entry>\n"
  let s10 = s9 +
    "<reaction name=\"rn:R00001\" type=\"irreversible\">\n" +
    "<substrate id=\"2\" name=\"cpd:C00001\"/>\n" +
    "<product id=\"3\" name=\"cpd:C00002\"/>\n" +
    "</reaction>\n" +
    "</pathway>"
  s10
}

///|
test "kegg fba multi gene entry name split" {
  // 空格分隔的多 KO ID 拆分为独立基因（KGML name 属性规范形式）
  let pw = @src.parse_kgml(kegg_fba_multi_gene_xml("ko:K00844 ko:K01623")).unwrap()
  let draft = @src.kegg_fba_reconstruct(pw, objective_reaction="R00001")
  assert_eq(draft.model.reactions[0].genes, ["K00844", "K01623"])
  // 逗号分隔形式同样拆分
  let pw2 = @src.parse_kgml(kegg_fba_multi_gene_xml("ko:K00844, ko:K01623")).unwrap()
  let draft2 = @src.kegg_fba_reconstruct(pw2, objective_reaction="R00001")
  assert_eq(draft2.model.reactions[0].genes, ["K00844", "K01623"])
  // 无空格逗号形式同样拆分
  let pw3 = @src.parse_kgml(kegg_fba_multi_gene_xml("ko:K00844,ko:K01623")).unwrap()
  let draft3 = @src.kegg_fba_reconstruct(pw3, objective_reaction="R00001")
  assert_eq(draft3.model.reactions[0].genes, ["K00844", "K01623"])
}

///|
test "kegg fba add equation duplicate compounds and ids" {
  let m = @src.fba_model_new([])
  // 同侧重复化合物系数合并：A + A -> 2B → A:-2, B:+2
  let m1 = @src.kegg_fba_add_equation(m, "R1", "A + A -> 2 B")
  assert_eq(m1.reactions.length(), 1)
  assert_eq(m1.metabolites.length(), 2)
  assert_true(m1.stoichiometry[0][0] == -2.0)
  assert_true(m1.stoichiometry[1][0] == 2.0)
  // 自转移方程净参与物为空 → 报错
  let zero = try {
    let _ = @src.kegg_fba_add_equation(m1, "R2", "A -> A")
    false
  } catch {
    _ => true
  }
  assert_true(zero)
  // 重复反应 ID → 报错
  let dup_eq = try {
    let _ = @src.kegg_fba_add_equation(m1, "R1", "A -> B")
    false
  } catch {
    _ => true
  }
  assert_true(dup_eq)
  // 重复 exchange ID → 报错
  let m2 = @src.kegg_fba_add_exchange(m1, "C")
  let dup_ex = try {
    let _ = @src.kegg_fba_add_exchange(m2, "C")
    false
  } catch {
    _ => true
  }
  assert_true(dup_ex)
}

///|
test "kegg fba reconstruct rejects duplicate pathway reactions" {
  // 复制 <reaction> 块制造重复反应 ID
  let base = kegg_fba_multi_gene_xml("ko:K00844")
  let xml2 = base.replace(
    old="</reaction>\n</pathway>",
    new="</reaction>\n<reaction name=\"rn:R00001\" type=\"irreversible\">\n<substrate id=\"2\" name=\"cpd:C00001\"/>\n<product id=\"3\" name=\"cpd:C00002\"/>\n</reaction>\n</pathway>",
  )
  let pw = @src.parse_kgml(xml2).unwrap()
  assert_eq(pw.reactions().length(), 2)
  let raised = try {
    let _ = @src.kegg_fba_reconstruct(pw)
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "kegg fba reconstruct merges duplicate participants" {
  // 同一化合物双侧同现：C00001 + C00002 -> C00001 + C00003
  // 正确净计量：C00002 -1、C00003 +1、C00001 净 0
  // （旧实现按列覆盖赋值会把 C00001 存成 +1，破坏质量守恒）
  let p = @src.KgmlPathway::new("path:ko00000", "ko", "00000", "merge test")
  let r = @src.KgmlReaction::new("rn:R99999", "irreversible")
  r.substrates().push("cpd:C00001")
  r.substrates().push("cpd:C00002")
  r.products().push("cpd:C00001")
  r.products().push("cpd:C00003")
  p.reactions().push(r)
  let draft = @src.kegg_fba_reconstruct(p)
  let mets = draft.model.metabolites
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
    for j in 0..<draft.model.reactions.length() {
      if draft.model.reactions[j].id == id {
        out = j
      }
    }
    out
  }
  assert_eq(draft.model.stoichiometry[idx("C00001")][ridx("R99999")], 0.0)
  assert_eq(draft.model.stoichiometry[idx("C00002")][ridx("R99999")], -1.0)
  assert_eq(draft.model.stoichiometry[idx("C00003")][ridx("R99999")], 1.0)
}

///|
test "kegg fba reconstruct rejects all-cancelling reaction" {
  // A -> A：净参与物为空，无法表达 → 报错
  let p = @src.KgmlPathway::new("path:ko00000", "ko", "00000", "zero test")
  let r = @src.KgmlReaction::new("rn:R99998", "irreversible")
  r.substrates().push("cpd:C00001")
  r.products().push("cpd:C00001")
  p.reactions().push(r)
  let raised = try {
    let _ = @src.kegg_fba_reconstruct(p)
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}
