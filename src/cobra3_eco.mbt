///|
/// cobra3_eco.mbt — 生态扩展（非 COBRApy 核心，但常用）.
///
/// - EFM 枚举（零空间子集 + 支撑极小化，限制规模）
/// - iMAT / GIMME / CORDA-lite / GrowMatch-lite
/// - OptKnock-lite（单/双敲除枚举 + 内层 FBA）
/// - 群落 FBA（模型合并 + 共享 exchange）
/// - 热力学约束 FBA（给定 ΔG，禁止不利方向）
///
/// 这些都是简化实现，边界条件在函数注释中说明。

///|
/// EFM 枚举：枚举零空间基的非空子集，保留满足边界/方向一致且支撑极小的模式.
pub fn cobra3_elementary_flux_modes(
  model : Cobra3Model,
  max_modes~ : Int32 = 500,
  epsilon~ : Double = 1.0e-9,
) -> Array[Array[Double]] {
  let n = model.reactions.length()
  let basis = cobra3_nullspace_basis(model)
  let k = basis.length()
  let out : Array[Array[Double]] = Array(capacity=0)
  if n == 0 || k == 0 || k > 12 {
    return out
  }
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  let mut total = 1
  for _ in 0..<k {
    total = total * 2
  }
  let supports : Array[Array[Bool]] = Array(capacity=total)
  let vectors : Array[Array[Double]] = Array(capacity=total)
  for mask in 1..<total {
    let vec : Array[Double] = Array::make(n, 0.0)
    for b in 0..<k {
      if ((mask >> b) & 1) == 1 {
        for j in 0..<n {
          vec[j] = vec[j] + basis[b][j]
        }
      }
    }
    // 边界与方向一致性
    let mut ok = true
    let mut any_positive = false
    let mut any_negative = false
    for j in 0..<n {
      let v = vec[j]
      if v.abs() > epsilon {
        if v < lb[j] - 1.0e-9 || v > ub[j] + 1.0e-9 {
          ok = false
          break
        }
        // 方向一致性只约束内部反应；exchange 可正可负
        if !cobra3_model_is_exchange(model, model.reactions[j]) {
          if v > 0.0 {
            any_positive = true
          } else {
            any_negative = true
          }
        }
      }
    }
    // 内部反应不得正负混合（可逆反应允许，但简化只保留单向模式）
    if ok && any_positive && any_negative {
      ok = false
    }
    if ok {
      let support : Array[Bool] = Array::make(n, false)
      for j in 0..<n {
        support[j] = vec[j].abs() > epsilon
      }
      supports.push(support)
      vectors.push(vec)
    }
  }
  // 支撑极小化：若存在另一个模式的支撑是其真子集，则丢弃
  for i in 0..<supports.length() {
    let si = supports[i]
    let mut minimal = true
    for m in 0..<supports.length() {
      if m == i {
        continue
      }
      let sm = supports[m]
      let mut subset = true
      let mut strict = false
      for j in 0..<n {
        if sm[j] {
          if !si[j] {
            subset = false
            break
          }
        } else if si[j] {
          strict = true
        }
      }
      if subset && strict {
        minimal = false
        break
      }
    }
    if minimal {
      out.push(vectors[i])
      if out.length() >= max_modes {
        break
      }
    }
  }
  out
}

///|
fn cobra3_reaction_ids_for_genes(
  model : Cobra3Model,
  genes : Array[String],
) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  if genes.length() == 0 {
    return out
  }
  let active : Map[String, Bool] = Map([], capacity=genes.length())
  for g in genes {
    active[g] = true
  }
  for r in model.reactions {
    match r.gpr {
      Some(gpr) => if cobra3_gpr_eval(gpr, active) { out.push(r.id) }
      None => ()
    }
  }
  out
}

///|
/// 按模型反应边界估算 big-M，避免硬编码 1000 在更大通量范围下失效.
fn cobra3_model_big_m(model : Cobra3Model) -> Double {
  let mut big_m = 1000.0
  for r in model.reactions {
    let scale = if r.lower_bound.abs() > r.upper_bound.abs() {
      r.lower_bound.abs()
    } else {
      r.upper_bound.abs()
    }
    if 2.0 * scale > big_m {
      big_m = 2.0 * scale
    }
  }
  big_m
}

///|
/// iMAT：最大化高表达反应开启 + 低表达反应关闭的一致性，保持生长.
pub fn cobra3_imat(
  model : Cobra3Model,
  high_genes : Array[String],
  low_genes : Array[String],
  min_growth~ : Double = 0.1,
  max_nodes~ : Int32 = 3000,
) -> Cobra3FluxResult {
  let n = model.reactions.length()
  let high_ids = cobra3_reaction_ids_for_genes(model, high_genes)
  let low_ids = cobra3_reaction_ids_for_genes(model, low_genes)
  let nv = 2 * n
  let s = cobra3_model_stoichiometric_matrix(model)
  let cons : Array[Cobra3LinConstraint] = Array(
    capacity=s.length() + 1 + high_ids.length() * 2 + low_ids.length() * 2,
  )
  for i in 0..<s.length() {
    let coeffs : Array[Double] = Array::make(nv, 0.0)
    for j in 0..<n {
      coeffs[j] = s[i][j]
    }
    cons.push({ coeffs, rhs: 0.0, sense: Cobra3Sense::Equal, })
  }
  let obj_vec = cobra3_model_objective_vector(model)
  let grow : Array[Double] = Array::make(nv, 0.0)
  for j in 0..<n {
    grow[j] = obj_vec[j]
  }
  cons.push({ coeffs: grow, rhs: min_growth, sense: Cobra3Sense::GreaterEqual, })
  let lb : Array[Double] = Array::make(nv, 0.0)
  let ub : Array[Double] = Array::make(nv, 0.0)
  for j in 0..<n {
    lb[j] = model.reactions[j].lower_bound
    ub[j] = model.reactions[j].upper_bound
  }
  for j in 0..<n {
    lb[n + j] = 0.0
    ub[n + j] = 1.0
  }
  let big_m = cobra3_model_big_m(model)
  for id in high_ids {
    let idx = cobra3_model_reaction_index(model, id)
    if idx >= 0 {
      // y=1 活跃: v >= eps；y=0 不活跃: v >= eps - M
      let c1 : Array[Double] = Array::make(nv, 0.0)
      c1[idx] = 1.0
      c1[n + idx] = -big_m
      cons.push({
        coeffs: c1,
        rhs: 1.0e-6 - big_m,
        sense: Cobra3Sense::GreaterEqual,
      })
    }
  }
  for id in low_ids {
    let idx = cobra3_model_reaction_index(model, id)
    if idx >= 0 {
      // y=1 活跃: 不限制；y=0 不活跃: |v| <= eps
      let c1 : Array[Double] = Array::make(nv, 0.0)
      c1[idx] = 1.0
      c1[n + idx] = -big_m
      cons.push({ coeffs: c1, rhs: 1.0e-6, sense: Cobra3Sense::LessEqual, })
      let c2 : Array[Double] = Array::make(nv, 0.0)
      c2[idx] = -1.0
      c2[n + idx] = -big_m
      cons.push({ coeffs: c2, rhs: 1.0e-6, sense: Cobra3Sense::LessEqual, })
    }
  }
  let obj : Array[Double] = Array::make(nv, 0.0)
  for id in high_ids {
    let idx = cobra3_model_reaction_index(model, id)
    if idx >= 0 {
      obj[n + idx] = 1.0
    }
  }
  for id in low_ids {
    let idx = cobra3_model_reaction_index(model, id)
    if idx >= 0 {
      obj[n + idx] = -1.0
    }
  }
  let integer : Array[Int32] = Array(capacity=n)
  for j in 0..<n {
    integer.push(n + j)
  }
  let sol = cobra3_solve_milp(cons, obj, lb, ub, integer, max_nodes~)
  let fluxes : Array[Double] = Array(capacity=n)
  for j in 0..<n {
    fluxes.push(sol.x[j])
  }
  {
    status: sol.status,
    objective_value: sol.objective_value,
    fluxes,
    shadow_prices: Array::make(model.metabolites.length(), 0.0),
    reduced_costs: Array::make(n, 0.0),
  }
}

///|
/// GIMME：在保持最低生长下，最小化低表达反应的加权通量.
pub fn cobra3_gimme(
  model : Cobra3Model,
  low_genes : Array[String],
  min_growth~ : Double = 0.1,
) -> Cobra3FluxResult {
  let n = model.reactions.length()
  let low_ids = cobra3_reaction_ids_for_genes(model, low_genes)
  let base = cobra3_steady_state_constraints(model)
  let obj_vec = cobra3_model_objective_vector(model)
  if min_growth > 0.0 {
    base.push({
      coeffs: obj_vec,
      rhs: min_growth,
      sense: Cobra3Sense::GreaterEqual,
    })
  }
  let nv = 2 * n
  let lb : Array[Double] = Array::make(nv, 0.0)
  let ub : Array[Double] = Array::make(nv, 0.0)
  let obj : Array[Double] = Array::make(nv, 0.0)
  for j in 0..<n {
    lb[j] = model.reactions[j].lower_bound
    ub[j] = model.reactions[j].upper_bound
  }
  for j in 0..<n {
    lb[n + j] = 0.0
    ub[n + j] = 1.0e18
    let penalty = if low_ids.contains(model.reactions[j].id) {
      1.0
    } else {
      0.0
    }
    obj[n + j] = -penalty
  }
  let cons : Array[Cobra3LinConstraint] = Array(capacity=base.length() + 2 * n)
  for con in base {
    let coeffs : Array[Double] = Array::make(nv, 0.0)
    for j in 0..<con.coeffs.length() {
      coeffs[j] = con.coeffs[j]
    }
    cons.push({ coeffs, rhs: con.rhs, sense: con.sense, })
  }
  for j in 0..<n {
    let c1 : Array[Double] = Array::make(nv, 0.0)
    c1[j] = 1.0
    c1[n + j] = -1.0
    cons.push({ coeffs: c1, rhs: 0.0, sense: Cobra3Sense::LessEqual, })
    let c2 : Array[Double] = Array::make(nv, 0.0)
    c2[j] = -1.0
    c2[n + j] = -1.0
    cons.push({ coeffs: c2, rhs: 0.0, sense: Cobra3Sense::LessEqual, })
  }
  let sol = cobra3_solve_lp(cons, obj, lb, ub)
  let fluxes : Array[Double] = Array(capacity=n)
  for j in 0..<n {
    fluxes.push(sol.x[j])
  }
  {
    status: sol.status,
    objective_value: -sol.objective_value,
    fluxes,
    shadow_prices: Array::make(model.metabolites.length(), 0.0),
    reduced_costs: Array::make(n, 0.0),
  }
}

///|
pub(all) struct Cobra3CordaResult {
  context_model : Cobra3Model
  removed_reactions : Array[String]
}

///|
/// CORDA-lite：用 GIMME 结果删除低表达且非必需的反应，得到上下文模型.
pub fn cobra3_corda(
  model : Cobra3Model,
  low_genes : Array[String],
  min_growth~ : Double = 0.1,
) -> Cobra3CordaResult {
  let gimme = cobra3_gimme(model, low_genes, min_growth~)
  let removed : Array[String] = Array(capacity=0)
  let low_ids = cobra3_reaction_ids_for_genes(model, low_genes)
  let mut context = model
  for j in 0..<model.reactions.length() {
    if gimme.fluxes[j].abs() < 1.0e-8 && low_ids.contains(model.reactions[j].id) {
      context = cobra3_model_remove_reaction(
        context,
        model.reactions[j].id,
        remove_orphans=false,
      )
      removed.push(model.reactions[j].id)
    }
  }
  { context_model: context, removed_reactions: removed, }
}

///|
pub(all) struct Cobra3GrowMatchResult {
  added_reactions : Array[String]
  predicted_growth : Double
}

///|
/// GrowMatch-lite：用 gapfill 找到使模型恢复生长的最少补充反应.
pub fn cobra3_growmatch(
  model : Cobra3Model,
  universal : Cobra3Model,
  min_growth~ : Double = 0.05,
) -> Cobra3GrowMatchResult {
  let gap = cobra3_gapfill(model, universal, lower_bound=min_growth)
  {
    added_reactions: gap.added_reactions,
    predicted_growth: gap.objective_value,
  }
}

///|
pub(all) struct Cobra3OptKnockResult {
  knockout_ids : Array[String]
  growth : Double
  target_flux : Double
}

///|
/// OptKnock 内层问题：先取敲除模型最大生长，再在保持该生长下最大化目标产物.
fn cobra3_optknock_inner(
  ko : Cobra3Model,
  target_idx : Int32,
) -> (Double, Double)? {
  let growth_sol = cobra3_model_optimize(ko)
  if !(growth_sol.status is Optimal) {
    return None
  }
  let growth = growth_sol.objective_value
  let base = cobra3_steady_state_constraints(ko)
  if growth.abs() > 1.0e-9 {
    base.push({
      coeffs: cobra3_model_objective_vector(ko),
      rhs: cobra3_optimum_floor(growth, 1.0 - 1.0e-6),
      sense: Cobra3Sense::GreaterEqual,
    })
  }
  let n = ko.reactions.length()
  let c : Array[Double] = Array::make(n, 0.0)
  if target_idx >= 0 && target_idx < n {
    c[target_idx] = 1.0
  }
  let sol = cobra3_solve_lp(
    base,
    c,
    cobra3_model_lower_bounds(ko),
    cobra3_model_upper_bounds(ko),
  )
  let target = if sol.status is Optimal { sol.objective_value } else { 0.0 }
  Some((growth, target))
}

///|
/// 递归枚举尺寸为 `k` 的敲除组合；每测试一个组合 `budget[0]` 递减，
/// 归零即停止，避免组合爆炸.
fn cobra3_optknock_extend(
  model : Cobra3Model,
  target_idx : Int32,
  k : Int32,
  start : Int32,
  cur : Array[String],
  out : Array[Cobra3OptKnockResult],
  budget : Array[Int32],
) -> Unit {
  let n = model.reactions.length()
  let need = k - cur.length()
  let mut i = start
  while i < n && n - i >= need && budget[0] > 0 {
    cur.push(model.reactions[i].id)
    if cur.length() == k {
      budget[0] = budget[0] - 1
      let ko = cobra3_knockout_reactions(model, cur)
      match cobra3_optknock_inner(ko, target_idx) {
        Some(pair) => {
          let ids : Array[String] = Array(capacity=cur.length())
          for id in cur {
            ids.push(id)
          }
          out.push({ knockout_ids: ids, growth: pair.0, target_flux: pair.1, })
        }
        None => ()
      }
    } else {
      cobra3_optknock_extend(model, target_idx, k, i + 1, cur, out, budget)
    }
    ignore(cur.pop())
    i = i + 1
  }
}

///|
/// OptKnock：枚举单反应直到 `max_knockouts` 元敲除（默认 2，与旧行为一致）；
/// 内层先取最大生长，再在保持最大生长下最大化目标产物。`max_candidates`
/// 限制内层 LP 的求解次数，达到上限即停止枚举（结果可能不完整）。
pub fn cobra3_optknock(
  model : Cobra3Model,
  target_id : String,
  max_knockouts~ : Int32 = 2,
  max_candidates~ : Int32 = 5000,
) -> Array[Cobra3OptKnockResult] {
  let out : Array[Cobra3OptKnockResult] = Array(capacity=0)
  let n = model.reactions.length()
  let target_idx = cobra3_model_reaction_index(model, target_id)
  if target_idx < 0 || n > 80 {
    return out
  }
  let budget : Array[Int32] = [max_candidates]
  let mut k = 1
  while k <= max_knockouts && k <= n && budget[0] > 0 {
    let cur : Array[String] = Array(capacity=k)
    cobra3_optknock_extend(model, target_idx, k, 0, cur, out, budget)
    k = k + 1
  }
  out
}

///|
pub(all) struct Cobra3CommunityResult {
  model : Cobra3Model
  fluxes : Array[Double]
  members : Array[String]
  growth : Double
}

///|
fn cobra3_community_shared_metabolites(
  models : Array[Cobra3Model],
) -> Map[String, Bool] {
  let shared : Map[String, Bool] = Map([], capacity=8)
  for m in models {
    for r in m.reactions {
      if cobra3_model_is_exchange(m, r) {
        for met_id in r.metabolites.keys() {
          shared[met_id] = true
        }
      }
    }
  }
  shared
}

///|
fn cobra3_community_gpr(gpr : Cobra3Gpr, prefix : String) -> Cobra3Gpr {
  match gpr {
    Gene(id) => Cobra3Gpr::Gene(prefix + id)
    And(children) => {
      let out : Array[Cobra3Gpr] = Array(capacity=children.length())
      for c in children {
        out.push(cobra3_community_gpr(c, prefix))
      }
      Cobra3Gpr::And(out)
    }
    Or(children) => {
      let out : Array[Cobra3Gpr] = Array(capacity=children.length())
      for c in children {
        out.push(cobra3_community_gpr(c, prefix))
      }
      Cobra3Gpr::Or(out)
    }
  }
}

///|
/// 群落 FBA-lite：合并多个模型；仅共享同名 exchange 及其代谢物，
/// 内部反应/代谢物/基因加模型前缀，避免同名 biomass 互相覆盖.
pub fn cobra3_community_fba(
  models : Array[Cobra3Model],
) -> Cobra3CommunityResult {
  let empty = cobra3_model_new("community")
  let shared_mets = cobra3_community_shared_metabolites(models)
  let members : Array[String] = Array(capacity=models.length())
  let mets : Array[Cobra3Metabolite] = Array(capacity=8)
  let genes : Array[Cobra3Gene] = Array(capacity=8)
  let reactions : Array[Cobra3Reaction] = Array(capacity=8)
  let seen_mets : Map[String, Bool] = Map([], capacity=8)
  let seen_genes : Map[String, Bool] = Map([], capacity=8)
  let seen_rxns : Map[String, Bool] = Map([], capacity=8)
  for m in models {
    members.push(m.id)
    for met in m.metabolites {
      let nid = if shared_mets.contains(met.id) {
        met.id
      } else {
        m.id + ":" + met.id
      }
      if !seen_mets.contains(nid) {
        seen_mets[nid] = true
        mets.push({
          id: nid,
          name: met.name,
          formula: met.formula,
          charge: met.charge,
          compartment: met.compartment,
          boundary: met.boundary,
          annotation: met.annotation,
        })
      }
    }
    for g in m.genes {
      let nid = m.id + ":" + g.id
      if !seen_genes.contains(nid) {
        seen_genes[nid] = true
        genes.push({
          id: nid,
          name: g.name,
          annotation: g.annotation,
          functional: g.functional,
        })
      }
    }
    for r in m.reactions {
      let is_exchange = cobra3_model_is_exchange(m, r)
      let nid = if is_exchange { r.id } else { m.id + ":" + r.id }
      let participants : Map[String, Double] = Map(
        [],
        capacity=r.metabolites.length(),
      )
      for met_id in r.metabolites.keys() {
        let nk = if shared_mets.contains(met_id) {
          met_id
        } else {
          m.id + ":" + met_id
        }
        participants[nk] = r.metabolites[met_id]
        if !seen_mets.contains(nk) {
          seen_mets[nk] = true
          let met = cobra3_metabolite_new(nk)
          let comp = cobra3_metabolite_compartment_from_id(nk)
          mets.push({
            id: met.id,
            name: met.name,
            formula: met.formula,
            charge: met.charge,
            compartment: comp,
            boundary: comp != "c",
            annotation: met.annotation,
          })
        }
      }
      let gpr = match r.gpr {
        Some(g) => Some(cobra3_community_gpr(g, m.id + ":"))
        None => None
      }
      if seen_rxns.contains(nid) {
        // 共享 exchange：边界取并集，允许任一成员原本的摄取/分泌范围
        for idx in 0..<reactions.length() {
          if reactions[idx].id == nid {
            let lb = if r.lower_bound < reactions[idx].lower_bound {
              r.lower_bound
            } else {
              reactions[idx].lower_bound
            }
            let ub = if r.upper_bound > reactions[idx].upper_bound {
              r.upper_bound
            } else {
              reactions[idx].upper_bound
            }
            reactions[idx] = {
              id: reactions[idx].id,
              name: reactions[idx].name,
              metabolites: reactions[idx].metabolites,
              lower_bound: lb,
              upper_bound: ub,
              objective: reactions[idx].objective,
              subsystem: reactions[idx].subsystem,
              gpr: reactions[idx].gpr,
              annotation: reactions[idx].annotation,
            }
          }
        }
      } else {
        seen_rxns[nid] = true
        reactions.push({
          id: nid,
          name: r.name,
          metabolites: participants,
          lower_bound: r.lower_bound,
          upper_bound: r.upper_bound,
          objective: r.objective,
          subsystem: r.subsystem,
          gpr,
          annotation: r.annotation,
        })
      }
    }
  }
  let merged : Cobra3Model = {
    id: empty.id,
    name: empty.name,
    metabolites: mets,
    genes,
    reactions,
    compartments: empty.compartments,
    annotation: empty.annotation,
    undo_stack: Array(capacity=0),
  }
  let sol = cobra3_model_optimize(merged)
  { model: merged, fluxes: sol.fluxes, members, growth: sol.objective_value, }
}

///|
/// 热力学约束 FBA-lite：给定反应 ΔG°（kJ/mol），若 ΔG° < 0 禁止反向，> 0 禁止正向.
pub fn cobra3_thermo_fba(
  model : Cobra3Model,
  delta_g : Map[String, Double],
) -> Cobra3FluxResult {
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for r in model.reactions {
    match delta_g.get(r.id) {
      Some(dg) => {
        let mut lb = r.lower_bound
        let mut ub = r.upper_bound
        if dg > 1.0e-6 {
          // 吸能：禁止正向；若原本只能正向，则整个反应关闭而非留下 lb>ub
          ub = 0.0
          if lb > 0.0 {
            lb = 0.0
          }
        } else if dg < -1.0e-6 {
          // 放能：禁止反向；若原本只能反向，则整个反应关闭而非留下 lb>ub
          lb = 0.0
          if ub < 0.0 {
            ub = 0.0
          }
        }
        reactions.push({
          id: r.id,
          name: r.name,
          metabolites: r.metabolites,
          lower_bound: lb,
          upper_bound: ub,
          objective: r.objective,
          subsystem: r.subsystem,
          gpr: r.gpr,
          annotation: r.annotation,
        })
      }
      None => reactions.push(r)
    }
  }
  let thermo : Cobra3Model = {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: Array(capacity=0),
  }
  cobra3_model_optimize(thermo)
}
