// cobra3_optknock.mbt — OptKnock 菌株设计：基因级敲除（bilevel 枚举式）
//
// 经典 OptKnock（Burgard et al. 2003, cameo `optknock`）：外层选择基因敲除
// 组合使目标产物产量最大化，内层细胞仍以最大生长为目标。完整 bilevel 需
// 带互补约束的 MILP；这里用工程上通用的枚举近似（同 cameo 的
// min_growth_fraction 语义）：
// 1. 完整模型 FBA 求基线生长 growth0（敲除只会移除通量能力，故任何敲除
//    组合的 growth <= growth0，即"max biomass ≤ 前值"自动成立）；
// 2. 枚举候选基因的所有 1..n_knockouts 组合（受 max_combinations 上限）；
// 3. 每个组合经 GPR 求值得到失活反应 → 固定其通量为 0；
// 4. 组合模型最大生长 < min_growth_fraction × growth0 的组合被剔除；
// 5. 在生物量下限约束下对目标反应做 FVA，target_max 即外层目标，
//    target_min > 0 表示产物成为生长耦合（growth-coupled）——敲不掉。
//
// 反应级敲除的旧版见 cobra3_eco.mbt 的 `cobra3_optknock`。

///|
pub(all) struct Cobra3OptKnockGeneResult {
  /// 敲除的基因 ID（按候选池出现顺序）
  knockouts : Array[String]
  /// 敲除后模型的最大生物量
  growth : Double
  /// 生物量下限约束下目标反应的最小通量
  target_min : Double
  /// 生物量下限约束下目标反应的最大通量（外层优化目标）
  target_max : Double
  /// target_min > 0：生长耦合产物，任何可行通量分布都必须产生产物
  coupled : Bool
}

///|
/// 基因级 OptKnock：枚举 `genes` 候选池的所有 1..`n_knockouts` 敲除组合
/// （经 GPR 求值映射到失活反应），返回按产物产量（target_max）降序排列的
/// 设计方案。`genes` 为空、目标反应不存在、模型无目标（生物量）或基线不
/// 生长时返回空列表。`max_combinations` 截断组合总数以控制组合爆炸。
pub fn cobra3_optknock_genes(
  model : Cobra3Model,
  target_id : String,
  genes : Array[String],
  n_knockouts~ : Int32 = 1,
  min_growth_fraction~ : Double = 0.1,
  max_combinations~ : Int32 = 200,
) -> Array[Cobra3OptKnockGeneResult] {
  let n = model.reactions.length()
  let t_idx = cobra3_model_reaction_index(model, target_id)
  if n == 0 || t_idx < 0 || genes.length() == 0 || n_knockouts <= 0 {
    return []
  }
  // 候选基因去重（保序）
  let seen : Map[String, Bool] = Map([])
  let pool : Array[String] = Array(capacity=genes.length())
  for g in genes {
    if not(seen.contains(g)) {
      seen.set(g, true)
      pool.push(g)
    }
  }
  let k_max = if n_knockouts > pool.length() {
    pool.length()
  } else {
    n_knockouts
  }
  // 目标（生物量）向量与基线生长
  let c_obj = cobra3_model_objective_vector(model)
  let base = cobra3_steady_state_constraints(model)
  let lb0 = cobra3_model_lower_bounds(model)
  let ub0 = cobra3_model_upper_bounds(model)
  let sol0 = cobra3_solve_lp(base, c_obj, lb0, ub0)
  if not(sol0.status is Optimal) {
    return []
  }
  let growth0 = sol0.objective_value
  if growth0 <= 0.0 {
    return []
  }
  let frac = if min_growth_fraction < 0.0 {
    0.0
  } else if min_growth_fraction > 1.0 {
    1.0
  } else {
    min_growth_fraction
  }
  let floor_growth = frac * growth0
  // 逐组合枚举
  let results : Array[Cobra3OptKnockGeneResult] = Array(capacity=0)
  let mut evaluated = 0
  for k in 1..<(k_max + 1) {
    // 组合下标生成器（字典序，next-combination 原地推进）
    let idx : Array[Int32] = Array(capacity=k)
    for i in 0..<k {
      idx.push(i)
    }
    let mut exhausted = false
    while not(exhausted) {
      if evaluated >= max_combinations {
        exhausted = true
        break
      }
      evaluated = evaluated + 1
      let ko_genes : Array[String] = Array(capacity=k)
      let ko_set : Map[String, Bool] = Map([])
      for i in 0..<k {
        let g = pool[idx[i]]
        ko_genes.push(g)
        ko_set.set(g, false)
      }
      // GPR 求值 → 失活反应。active map 必须覆盖模型全部 GPR 基因
      // （cobra3_gpr_eval 对缺失 key 返回 false，只放候选池基因会把
      // 池外基因控制的反应错误禁用）
      let active : Map[String, Bool] = Map([])
      for r in model.reactions {
        match r.gpr {
          Some(gpr) =>
            for g in cobra3_gpr_genes(gpr) {
              if not(active.contains(g)) {
                active.set(g, true)
              }
            }
          None => ()
        }
      }
      for g in ko_set.keys() {
        active.set(g, false)
      }
      let mut n_off = 0
      let lb = Array(capacity=n)
      let ub = Array(capacity=n)
      for j in 0..<n {
        lb.push(lb0[j])
        ub.push(ub0[j])
      }
      for r in model.reactions {
        match r.gpr {
          None => continue
          Some(gpr) =>
            if not(cobra3_gpr_eval(gpr, active)) {
              let j = cobra3_model_reaction_index(model, r.id)
              lb[j] = 0.0
              ub[j] = 0.0
              n_off = n_off + 1
            }
        }
      }
      // 未失活任何反应的组合与基线等价，跳过
      if n_off == 0 {
        exhausted = cobra3_optknock_genes_advance(idx, pool.length(), k)
        continue
      }
      // 组合模型最大生长
      let sol_g = cobra3_solve_lp(base, c_obj, lb, ub)
      let mut growth = 0.0
      if sol_g.status is Optimal {
        growth = sol_g.objective_value
      }
      if growth < floor_growth - 1.0e-9 {
        exhausted = cobra3_optknock_genes_advance(idx, pool.length(), k)
        continue
      }
      // 生物量下限约束下的目标产物 FVA
      let cons : Array[Cobra3LinConstraint] = Array(capacity=base.length() + 1)
      for c in base {
        cons.push(c)
      }
      cons.push({ coeffs: c_obj, rhs: floor_growth, sense: GreaterEqual, })
      let c_t : Array[Double] = Array::make(n, 0.0)
      c_t[t_idx] = 1.0
      let sol_max = cobra3_solve_lp(cons, c_t, lb, ub)
      c_t[t_idx] = -1.0
      let sol_min = cobra3_solve_lp(cons, c_t, lb, ub)
      let mut t_max = 0.0
      let mut t_min = 0.0
      if sol_max.status is Optimal {
        t_max = sol_max.objective_value
      }
      if sol_min.status is Optimal {
        t_min = -sol_min.objective_value
      }
      results.push({
        knockouts: ko_genes,
        growth,
        target_min: t_min,
        target_max: t_max,
        coupled: t_min > 1.0e-9,
      })
      exhausted = cobra3_optknock_genes_advance(idx, pool.length(), k)
    }
  }
  // 按 target_max 降序（并列按首个敲除基因 ID 升序，保证确定性）
  let order : Array[Int32] = Array(capacity=results.length())
  for i in 0..<results.length() {
    order.push(i)
  }
  for i in 1..<order.length() {
    let cur = order[i]
    let mut j = i - 1
    while j >= 0 {
      let a = results[order[j]]
      let b = results[cur]
      let worse = a.target_max < b.target_max - 1.0e-12 ||
        (
          a.target_max <= b.target_max + 1.0e-12 &&
          a.knockouts[0] > b.knockouts[0]
        )
      if worse {
        order[j + 1] = order[j]
        j = j - 1
      } else {
        break
      }
    }
    order[j + 1] = cur
  }
  let sorted : Array[Cobra3OptKnockGeneResult] = Array(capacity=order.length())
  for t in 0..<order.length() {
    sorted.push(results[order[t]])
  }
  sorted
}

///|
/// next-combination：把下标组合原地推进到字典序下一个；无下一个时返回
/// true。
fn cobra3_optknock_genes_advance(
  idx : Array[Int32],
  n : Int32,
  k : Int32,
) -> Bool {
  let mut i = k - 1
  while i >= 0 {
    if idx[i] < n - (k - i) {
      idx[i] = idx[i] + 1
      let s = i + 1
      for j in s..<k {
        idx[j] = idx[j - 1] + 1
      }
      return false
    }
    i = i - 1
  }
  true
}
