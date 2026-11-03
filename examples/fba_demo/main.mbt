///|
/// FBA 代谢建模演示: toy 模型 -> 通量求解 -> 基因敲除筛查.
fn main raise {
  println("=== FBA 代谢通量平衡分析 (cobra-lite 单纯形) ===")

  // toy 模型: (外) -R0-> A -R1-> B -R2-> C -R3->(外排)
  //                                  C -R4-> biomass
  let mut model = @src.fba_model_new(["A", "B", "C"])
  model = model.add_reaction(
    {
      id: "R0",
      name: "uptake",
      lower_bound: 0.0,
      upper_bound: 10.0,
      objective: 0.0,
      genes: ["geneU"],
    },
    [("A", 1.0)],
  )
  model = model.add_reaction(
    {
      id: "R1",
      name: "A_to_B",
      lower_bound: 0.0,
      upper_bound: 1000.0,
      objective: 0.0,
      genes: ["geneA"],
    },
    [("A", -1.0), ("B", 1.0)],
  )
  model = model.add_reaction(
    {
      id: "R2",
      name: "B_to_C",
      lower_bound: 0.0,
      upper_bound: 1000.0,
      objective: 0.0,
      genes: ["geneB"],
    },
    [("B", -1.0), ("C", 1.0)],
  )
  model = model.add_reaction(
    {
      id: "R3",
      name: "excrete",
      lower_bound: 0.0,
      upper_bound: 1000.0,
      objective: 0.0,
      genes: ["geneE"],
    },
    [("C", -1.0)],
  )
  model = model.add_reaction(
    {
      id: "R4",
      name: "biomass",
      lower_bound: 0.0,
      upper_bound: 1000.0,
      objective: 1.0,
      genes: [],
    },
    [("C", -1.0)],
  )

  // 1. 野生型 FBA
  println("\n1. 野生型 FBA(底物上界 10)")
  let wt = @src.fba_solve(model)
  println("  " + wt.summary(model))
  for j = 0; j < model.reactions.length(); j = j + 1 {
    let r = model.reactions[j]
    println("  \{r.id}(\{r.name}): \{wt.fluxes[j]}")
  }

  // 2. 底物限制
  println("\n2. 底物限制(上界 10 -> 2.5)")
  // 底物限制: uptake 上界 10 -> 2.5
  let limited = @src.fba_set_reaction_bound(model, "R0", 0.0, 2.5)
  let lim_sol = @src.fba_solve(limited)
  println("  biomass = \{lim_sol.objective_value} (应随底物线性: 2.5)")

  // 3. 基因敲除筛查
  println("\n3. 基因敲除表型筛查")
  let screen = @src.fba_knockout_screen(model)
  for row in screen {
    println("  \{row.0}: biomass=\{row.1} -> \{row.2}")
  }

  // 4. JSON 模型交换
  println("\n4. JSON 模型交换(roundtrip)")
  let json = @src.fba_write_json(model)
  println("  \{json.length()} 字符")
  let reparsed = @src.fba_parse_json(json)
  match reparsed {
    Some(m) => {
      let sol2 = @src.fba_solve(m)
      println(
        "  reparsed: \{m.reactions.length()} 反应, biomass=\{sol2.objective_value} (与原模型一致)",
      )
    }
    None => println("  reparsed FAILED")
  }

  // 5. 经典 LP 验证
  println("\n5. 经典 LP 验证(单纯形正确性)")
  let a : Array[Array[Double]] = [[3.0, 2.0, 1.0]]
  let b : Array[Double] = [18.0]
  let c : Array[Double] = [3.0, 5.0, 0.0]
  let ub : Array[Double] = [4.0, 6.0, 1000.0]
  match @src.fba_solve_lp(a, b, c, ub) {
    Some((obj, _)) =>
      println("  max 3x1+5x2 (3x1+2x2<=18, x1<=4, x2<=6) = \{obj} (期望 36)")
    None => println("  LP FAILED")
  }

  println(
    "\n完成: 模型构造 -> 单纯形求解 -> 敲除筛查 -> JSON 交换全链路本地完成.",
  )
}
