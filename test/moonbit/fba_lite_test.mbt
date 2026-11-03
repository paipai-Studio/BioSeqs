// ===========================================================================
// fba_lite_test.mbt — FBA 代谢通量平衡测试(单纯形/toy 模型/敲除)
// ===========================================================================

///|
test "fba solve lp textbook max" {
  // max x1 + x2, s.t. x1 + 2*x2 = 4(拆成 <= 的等价形式), x1 <= 2
  // 标准型: x1 + 2*x2 + s1 = 4; x1 + s2 = 2
  // 这里直接用等式+上界形式: x1+2x2 <= 4 由上界与约束组合给出
  // 简化例子: max 3x1+5x2, x1<=4, 2x2<=12, 3x1+2x2<=18(经典)
  // 用等式约束行 + slack 已内置于 fba_solve_lp 的上界机制:
  // 改用显式等式: x1+x2=4 与 x1<=2 -> max x1+x2 的解 x1=2, x2=1, obj=4
  let a : Array[Array[Double]] = [[1.0, 2.0]]
  let b : Array[Double] = [4.0]
  let c : Array[Double] = [1.0, 1.0]
  let ub : Array[Double] = [2.0, 100.0]
  let (obj, x) = @src.fba_solve_lp(a, b, c, ub).unwrap()
  // x1+2x2=4, x1<=2 -> x1=2, x2=1, obj=3
  assert_true(obj > 2.99 && obj < 3.01)
  assert_true(x[0] > 1.99 && x[0] < 2.01)
  assert_true(x[1] > 0.99 && x[1] < 1.01)
}

///|
test "fba solve lp classic simplex" {
  // 经典: max 3x1+5x2, x1<=4, 2x2<=12, 3x1+2x2<=18 -> 最优 36
  // 等式化: 无等式约束, 仅上界+行约束 -> 用 a_eq 空 + 附加行
  // fba_solve_lp 的上界是逐变量的, 行约束需进 a_eq:
  // 3x1+2x2+s=18 -> 将 s 视为第三变量 c=0, ub=inf
  let a : Array[Array[Double]] = [[3.0, 2.0, 1.0]]
  let b : Array[Double] = [18.0]
  let c : Array[Double] = [3.0, 5.0, 0.0]
  let ub : Array[Double] = [4.0, 6.0, 1000.0]
  let (obj, _) = @src.fba_solve_lp(a, b, c, ub).unwrap()
  assert_true(obj > 35.99 && obj < 36.01)
}

///|
test "fba solve lp infeasible" {
  // x1 <= 2 但 x1 = 5
  let a : Array[Array[Double]] = [[1.0]]
  let b : Array[Double] = [5.0]
  let c : Array[Double] = [1.0]
  let ub : Array[Double] = [2.0]
  let r = @src.fba_solve_lp(a, b, c, ub)
  assert_true(r is None)
}

///|
test "fba solve lp unbounded" {
  // max x, 无约束行, 上界无穷
  let a : Array[Array[Double]] = []
  let b : Array[Double] = []
  let c : Array[Double] = [1.0]
  let ub : Array[Double] = [1.0e18]
  let r = @src.fba_solve_lp(a, b, c, ub)
  assert_true(r is None)
}

///|
test "fba solve lp negative rhs" {
  // rhs 负行: -x1 = -3 -> x1 = 3, max x1, ub=10
  let a : Array[Array[Double]] = [[-1.0]]
  let b : Array[Double] = [-3.0]
  let c : Array[Double] = [1.0]
  let ub : Array[Double] = [10.0]
  let (obj, x) = @src.fba_solve_lp(a, b, c, ub).unwrap()
  assert_true(obj > 2.99 && obj < 3.01)
  assert_true(x[0] > 2.99 && x[0] < 3.01)
}

///|
/// toy 代谢模型: A -(R1)-> B -(R2)-> C -(R3)-> 外排; 源 M -(R0)-> A; C -> biomass(R4)
fn fba_toy_model() -> @src.FbaModel raise {
  // metabolites: 内部代谢物 A, B, C(底物摄入为边界反应, 源不进稳态计量)
  let mut model = @src.fba_model_new(["A", "B", "C"])
  // R0: (外) -> A (底物摄入, 上界 10)
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
  // R1: A -> B
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
  // R2: B -> C
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
  // R3: C -> (外排)
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
  // R4: C -> biomass(目标)
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
  model
}

///|
test "fba toy model optimal growth" {
  let model = fba_toy_model() catch {
    _ => abort("toy model construction failed")
  }
  let sol = @src.fba_solve(model)
  assert_eq(sol.status, "optimal")
  // 底物上界 10 -> biomass = 10
  assert_true(sol.objective_value > 9.99 && sol.objective_value < 10.01)
  // 通量守恒: R0=R1=R2, R3+R4=R2
  assert_true(sol.fluxes[0] > 9.99 && sol.fluxes[0] < 10.01)
  assert_true(sol.fluxes[1] > 9.99 && sol.fluxes[1] < 10.01)
  let sum_out = sol.fluxes[3] + sol.fluxes[4]
  assert_true(sum_out > 9.99 && sum_out < 10.01)
}

///|
test "fba toy model steady state mass balance" {
  let model = fba_toy_model() catch {
    _ => abort("toy model construction failed")
  }
  let sol = @src.fba_solve(model)
  // S·v = 0 逐代谢物验证
  for i = 0; i < model.metabolites.length(); i = i + 1 {
    let mut bal = 0.0
    for j = 0; j < model.reactions.length(); j = j + 1 {
      bal = bal + model.stoichiometry[i][j] * sol.fluxes[j]
    }
    assert_true(bal > -1.0e-6 && bal < 1.0e-6)
  }
}

///|
test "fba knockout essential gene" {
  let model = fba_toy_model() catch {
    _ => abort("toy model construction failed")
  }
  // geneA(R1 A->B) 敲除 -> 通路断, biomass = 0
  let (ko_model, knocked) = @src.fba_knockout(model, ["geneA"])
  assert_eq(knocked, ["R1"])
  let sol = @src.fba_solve(ko_model)
  assert_true(sol.objective_value < 1.0e-6)
}

///|
test "fba knockout nonessential excretion" {
  let model = fba_toy_model() catch {
    _ => abort("toy model construction failed")
  }
  // geneE(R3 外排) 敲除 -> 全部 C 走 biomass, 生长不变
  let (ko_model, knocked) = @src.fba_knockout(model, ["geneE"])
  assert_eq(knocked, ["R3"])
  let sol = @src.fba_solve(ko_model)
  assert_true(sol.objective_value > 9.99 && sol.objective_value < 10.01)
  // R3 通量应为 0
  assert_true(sol.fluxes[3] < 1.0e-6 && sol.fluxes[3] > -1.0e-6)
}

///|
test "fba knockout screen phenotypes" {
  let model = fba_toy_model() catch {
    _ => abort("toy model construction failed")
  }
  let screen = @src.fba_knockout_screen(model)
  // geneU/geneA/geneB 敲断通路 -> non-growth
  // geneE(外排) -> growth
  assert_eq(screen.length(), 4)
  let results : Map[String, String] = Map([], capacity=0)
  for row in screen {
    results.set(row.0, row.2)
  }
  assert_eq(results.get("geneU").unwrap(), "non-growth")
  assert_eq(results.get("geneA").unwrap(), "non-growth")
  assert_eq(results.get("geneB").unwrap(), "non-growth")
  assert_eq(results.get("geneE").unwrap(), "growth")
}
