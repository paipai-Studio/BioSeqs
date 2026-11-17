// cobra3_cfba.mbt — cFBA 简化版：蛋白质配额下的资源分配约束
// （Sanchez et al. 2017, "Buying cells: metabolic enzyme budget"）
//
// 在标准 FBA 的稳态约束之上加一条线性"酶质量约束"：
//   Σ_j (MW_j / k_eff) · v_j ≤ P_total
// 其中 MW_j 为反应 j 的酶分子量（g/mol，来自 enzyme_mw 映射），k_eff 为
// 比催化速率（1/s），P_total 为蛋白质配额（g/gDW）。 LP 直接重解即可
// 得到配额下的最大生长。
//
// 简化口径（边界条件）：
// - 只对 enzyme_mw 中登记的反应计酶成本，未登记的反应免费（无酶旁路
//   可绕过配额，符合 cFBA 语义：只有有酶参与的通量消耗蛋白）
// - 酶成本按净正向通量线性计入：可逆反应负向通量会"退还"成本，因此
//   建议只给不可逆反应登记 MW；对可逆反应的严格双方向计费需要拆分
//   变量（简化版不做）
// - LP 不可行（如 P 不足以维持最小生长需求）时返回 growth = 0

///|
/// cFBA 参数错误.
pub suberror Cobra3Error {
  Cobra3Error(String)
}

///|
/// cFBA 求解结果.
pub(all) struct Cobra3CfbaResult {
  /// 配额下最大目标值（生物量通量）
  growth : Double
  /// 每反应通量（与模型反应同序；不可行时全 0）
  fluxes : Array[Double]
  /// 实际酶质量占用 Σ (MW/k_eff)·v（配额绑定时 = protein_budget）
  enzyme_load : Double
}

///|
/// 蛋白配额扫描点.
pub(all) struct Cobra3CfbaSweepPoint {
  /// 蛋白配额 P
  protein : Double
  /// 该配额下的最大生长
  growth : Double
}

///|
/// cFBA：在稳态约束 + 酶质量约束下最大化模型目标（生物量）。
/// `enzyme_mw`：反应 ID → 酶分子量（未登记反应不消耗蛋白配额）；
/// `protein_budget`：蛋白配额 P（g/gDW）；`k_eff`：统一比催化速率。
/// `protein_budget = 0` 且酶约束绑定时生长为 0；映射为空等价于普通 FBA。
/// `k_eff <= 0` 无物理意义（非正比催化速率会使约束失去意义），直接报错。
pub fn cobra3_cfba(
  model : Cobra3Model,
  enzyme_mw : Map[String, Double],
  protein_budget~ : Double = 1.0,
  k_eff~ : Double = 65.0,
) -> Cobra3CfbaResult raise Cobra3Error {
  if k_eff <= 0.0 {
    raise Cobra3Error("cFBA requires k_eff > 0, got \{k_eff}")
  }
  let n = model.reactions.length()
  let cons = cobra3_steady_state_constraints(model)
  // 酶质量约束：Σ (MW_j / k_eff) · v_j ≤ P
  let coeffs : Array[Double] = Array::make(n, 0.0)
  let mut n_assigned = 0
  for kv in enzyme_mw {
    let (rid, mw) = kv
    let idx = cobra3_model_reaction_index(model, rid)
    if idx >= 0 && mw > 0.0 && k_eff > 0.0 {
      coeffs[idx] = mw / k_eff
      n_assigned = n_assigned + 1
    }
  }
  if n_assigned > 0 {
    cons.push({ coeffs, rhs: protein_budget, sense: Cobra3Sense::LessEqual, })
  }
  let c = cobra3_model_objective_vector(model)
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  let sol = cobra3_solve_lp(cons, c, lb, ub)
  if sol.status is Optimal {
    let mut load = 0.0
    for j in 0..<n {
      load += coeffs[j] * sol.x[j]
    }
    { growth: sol.objective_value, fluxes: sol.x, enzyme_load: load, }
  } else {
    { growth: 0.0, fluxes: Array::make(n, 0.0), enzyme_load: 0.0, }
  }
}

///|
/// 蛋白配额扫描：在 [p_min, p_max] 上均匀取 `points` 个配额逐点求解，
/// 输出经典的"生长 vs 蛋白配额"曲线（配额小时线性上升，底物受限后
/// 进入平台期）。`k_eff` 透传给每次求解（入口统一校验，非正值报错）。
/// `points < 1` 返回空列表。
pub fn cobra3_cfba_sweep(
  model : Cobra3Model,
  enzyme_mw : Map[String, Double],
  points~ : Int32 = 10,
  p_min~ : Double = 0.2,
  p_max~ : Double = 2.0,
  k_eff~ : Double = 65.0,
) -> Array[Cobra3CfbaSweepPoint] raise Cobra3Error {
  // 参数校验与边界判断先于数组构造（负 capacity 会 panic）
  if k_eff <= 0.0 {
    raise Cobra3Error("cFBA requires k_eff > 0, got \{k_eff}")
  }
  if points < 1 {
    return []
  }
  let out : Array[Cobra3CfbaSweepPoint] = Array(capacity=points)
  for i in 0..<points {
    let t = if points == 1 {
      0.0
    } else {
      i.to_double() / (points - 1).to_double()
    }
    let p = p_min + t * (p_max - p_min)
    let r = cobra3_cfba(model, enzyme_mw, protein_budget=p, k_eff~)
    out.push({ protein: p, growth: r.growth, })
  }
  out
}
