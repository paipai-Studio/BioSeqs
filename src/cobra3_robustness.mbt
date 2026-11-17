// cobra3_robustness.mbt — 鲁棒性分析（单反应通量扫描）
//
// 经典 robustness analysis（Edwards & Palsson）：把某一反应的通量固定在
// 其可行范围内的扫描点上，逐点重新最大化生物量目标，得到"通量-生长"
// 曲线。曲线的平台段说明生长对该反应不敏感，下降段给出该反应的瓶颈
// 强度。与 cobra3_production_envelope（二维：目标 vs 产物通量）互补的
// 单维版本。

///|
pub(all) struct Cobra3RobustnessPoint {
  /// 扫描点上被固定的目标反应通量
  target_flux : Double
  /// 该点可行域内模型目标（生物量）的最大值
  objective_max : Double
  /// 该点可行域内模型目标的最小值（分析容错冗余用）
  objective_min : Double
}

///|
/// 鲁棒性曲线：把 `target_id` 的通量从其最小可行值扫到最大可行值（共
/// `points` 个等距点），每个点固定该通量后重新最大化/最小化模型目标。
/// 与 COBRApy 生态中 cameo 的 robustness analysis 对应。
pub fn cobra3_robustness(
  model : Cobra3Model,
  target_id : String,
  points~ : Int32 = 20,
) -> Array[Cobra3RobustnessPoint] {
  let idx = cobra3_model_reaction_index(model, target_id)
  let out : Array[Cobra3RobustnessPoint] = Array(capacity=points)
  if idx < 0 || points <= 0 {
    return out
  }
  let n = model.reactions.length()
  if n == 0 {
    return out
  }
  // 第一步：求目标反应自身的通量范围（两个独立 LP，与目标解耦）
  let base = cobra3_steady_state_constraints(model)
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  let e_max : Array[Double] = Array::make(n, 0.0)
  e_max[idx] = 1.0
  let sol_hi = cobra3_solve_lp(base, e_max, lb, ub)
  let e_min : Array[Double] = Array::make(n, 0.0)
  e_min[idx] = -1.0
  let sol_lo = cobra3_solve_lp(base, e_min, lb, ub)
  let mut lo = 0.0
  let mut hi = 0.0
  if sol_hi.status is Optimal {
    hi = sol_hi.objective_value
  }
  if sol_lo.status is Optimal {
    lo = -sol_lo.objective_value
  }
  if hi < lo {
    return out
  }
  // 第二步：逐点固定目标通量，重新优化模型目标
  let c_max : Array[Double] = Array::make(n, 0.0)
  let c_min : Array[Double] = Array::make(n, 0.0)
  for j in 0..<n {
    c_max[j] = model.reactions[j].objective
    c_min[j] = -model.reactions[j].objective
  }
  let fix : Array[Double] = Array::make(n, 0.0)
  fix[idx] = 1.0
  for p in 0..<points {
    let t = if points > 1 {
      p.to_double() / (points - 1).to_double()
    } else {
      0.0
    }
    let target_flux = lo + t * (hi - lo)
    let cons : Array[Cobra3LinConstraint] = Array(capacity=base.length() + 1)
    for c in base {
      cons.push(c)
    }
    cons.push({ coeffs: fix, rhs: target_flux, sense: Equal, })
    let sol_max = cobra3_solve_lp(cons, c_max, lb, ub)
    let sol_min = cobra3_solve_lp(cons, c_min, lb, ub)
    let mut obj_max = 0.0
    let mut obj_min = 0.0
    if sol_max.status is Optimal {
      obj_max = sol_max.objective_value
    }
    if sol_min.status is Optimal {
      obj_min = -sol_min.objective_value
    }
    out.push({ target_flux, objective_max: obj_max, objective_min: obj_min, })
  }
  out
}
