///|
/// cobra3_analysis.mbt — 模型分析.
///
/// - 通量耦合分析（fully / partial / directional / none）
/// - 相平面（双反应通量网格）
/// - 模型/反应/代谢物摘要
/// - 敏感性（shadow price 排序）
/// - 模型差异比较

///|
pub(all) struct Cobra3CouplingPair {
  reaction_a : String
  reaction_b : String
  kind : String
}

///|
/// 通量耦合分析-lite：逐个强制 v_a >= eps，检查 v_b 的范围.
/// 反应数超过 max_reactions 时直接返回空数组（不做部分分析）.
pub fn cobra3_flux_coupling(
  model : Cobra3Model,
  epsilon~ : Double = 1.0e-3,
  max_reactions~ : Int32 = 40,
) -> Array[Cobra3CouplingPair] {
  let n = model.reactions.length()
  let out : Array[Cobra3CouplingPair] = Array(capacity=0)
  if n == 0 || n > max_reactions {
    return out
  }
  let base = cobra3_steady_state_constraints(model)
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  for i in 0..<n {
    let force : Array[Double] = Array::make(n, 0.0)
    force[i] = 1.0
    let cons : Array[Cobra3LinConstraint] = Array(capacity=base.length() + 1)
    for c in base {
      cons.push(c)
    }
    cons.push({ coeffs: force, rhs: epsilon, sense: Cobra3Sense::GreaterEqual, })
    for j in 0..<n {
      if i == j {
        continue
      }
      let id_a = model.reactions[i].id
      let id_b = model.reactions[j].id
      let c_max : Array[Double] = Array::make(n, 0.0)
      c_max[j] = 1.0
      let c_min : Array[Double] = Array::make(n, 0.0)
      c_min[j] = -1.0
      let smax = cobra3_solve_lp(cons, c_max, lb, ub)
      let smin = cobra3_solve_lp(cons, c_min, lb, ub)
      // 强制 a 正本身不可行（或未解出）时，任何 b 的耦合都无意义
      if !(smax.status is Optimal) || !(smin.status is Optimal) {
        continue
      }
      let mx = smax.objective_value
      let mn = -smin.objective_value
      let kind = if mn >= epsilon * 0.99 {
        // 强制 a 正时，b 必须正
        "fully_coupled"
      } else if mn >= -epsilon * 1.01 && mx >= epsilon * 0.99 {
        // b 可为正但不能显著为负
        "directionally_coupled"
      } else if mx <= epsilon * 1.01 && mn <= -epsilon * 0.99 {
        // b 可为负但不能显著为正
        "directionally_coupled"
      } else if mx - mn > epsilon {
        "partial_coupled"
      } else {
        "none"
      }
      if kind != "none" {
        out.push({ reaction_a: id_a, reaction_b: id_b, kind, })
      }
    }
  }
  out
}

///|
pub(all) struct Cobra3PhasePlanePoint {
  x_flux : Double
  y_flux : Double
  objective : Double
}

///|
/// 相平面：固定两个反应的通量，求目标函数值（不可行点的目标记为 0.0）.
pub fn cobra3_phase_plane(
  model : Cobra3Model,
  reaction_x : String,
  reaction_y : String,
  points~ : Int32 = 20,
) -> Array[Cobra3PhasePlanePoint] {
  let out : Array[Cobra3PhasePlanePoint] = Array(capacity=points * points)
  let ix = cobra3_model_reaction_index(model, reaction_x)
  let iy = cobra3_model_reaction_index(model, reaction_y)
  if ix < 0 || iy < 0 || points <= 0 {
    return out
  }
  let fva = cobra3_fva(model, fraction_of_optimum=0.0)
  let mut xmin = 0.0
  let mut xmax = 0.0
  let mut ymin = 0.0
  let mut ymax = 0.0
  for f in fva {
    if f.reaction_id == reaction_x {
      xmin = f.minimum
      xmax = f.maximum
    }
    if f.reaction_id == reaction_y {
      ymin = f.minimum
      ymax = f.maximum
    }
  }
  let n = model.reactions.length()
  let base = cobra3_steady_state_constraints(model)
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  let c = cobra3_model_objective_vector(model)
  for px in 0..<points {
    let tx = if points > 1 {
      px.to_double() / (points - 1).to_double()
    } else {
      0.0
    }
    let xv = xmin + tx * (xmax - xmin)
    for py in 0..<points {
      let ty = if points > 1 {
        py.to_double() / (points - 1).to_double()
      } else {
        0.0
      }
      let yv = ymin + ty * (ymax - ymin)
      let cons : Array[Cobra3LinConstraint] = Array(capacity=base.length() + 2)
      for con in base {
        cons.push(con)
      }
      let fx : Array[Double] = Array::make(n, 0.0)
      fx[ix] = 1.0
      cons.push({ coeffs: fx, rhs: xv, sense: Cobra3Sense::Equal, })
      let fy : Array[Double] = Array::make(n, 0.0)
      fy[iy] = 1.0
      cons.push({ coeffs: fy, rhs: yv, sense: Cobra3Sense::Equal, })
      let sol = cobra3_solve_lp(cons, c, lb, ub)
      let obj = if sol.status is Optimal { sol.objective_value } else { 0.0 }
      out.push({ x_flux: xv, y_flux: yv, objective: obj, })
    }
  }
  out
}

///|
pub(all) struct Cobra3ModelSummary {
  model_id : String
  n_metabolites : Int32
  n_reactions : Int32
  n_genes : Int32
  n_exchange : Int32
  objective : String
  growth : Double
  shadow_sum : Double
}

///|
pub fn cobra3_model_summary(model : Cobra3Model) -> Cobra3ModelSummary {
  let sol = cobra3_model_optimize(model)
  let mut shadow_sum = 0.0
  for v in sol.shadow_prices {
    shadow_sum = shadow_sum + v.abs()
  }
  let objective = match cobra3_model_objective_reaction(model) {
    Some(id) => id
    None => ""
  }
  {
    model_id: model.id,
    n_metabolites: model.metabolites.length(),
    n_reactions: model.reactions.length(),
    n_genes: model.genes.length(),
    n_exchange: cobra3_model_exchange_reactions(model).length(),
    objective,
    growth: sol.objective_value,
    shadow_sum,
  }
}

///|
pub fn cobra3_model_summary_string(model : Cobra3Model) -> String {
  let s = cobra3_model_summary(model)
  let sb = StringBuilder()
  let _ = sb.write("Model " + s.model_id + "\n")
  let _ = sb.write("  metabolites: " + s.n_metabolites.to_string() + "\n")
  let _ = sb.write("  reactions: " + s.n_reactions.to_string() + "\n")
  let _ = sb.write("  genes: " + s.n_genes.to_string() + "\n")
  let _ = sb.write("  exchanges: " + s.n_exchange.to_string() + "\n")
  let _ = sb.write("  objective: " + s.objective + "\n")
  let _ = sb.write("  growth: " + s.growth.to_string() + "\n")
  sb.to_string()
}

///|
pub(all) struct Cobra3ReactionSummary {
  reaction_id : String
  name : String
  flux : Double
  lower_bound : Double
  upper_bound : Double
  shadow_price_of_metabolites : Double
  gpr : String
  participants : Array[(String, Double)]
}

///|
pub fn cobra3_reaction_summary(
  model : Cobra3Model,
  reaction_id : String,
) -> Cobra3ReactionSummary? {
  let idx = cobra3_model_reaction_index(model, reaction_id)
  if idx < 0 {
    return None
  }
  let sol = cobra3_model_optimize(model)
  let r = model.reactions[idx]
  let gpr = match r.gpr {
    Some(g) => cobra3_gpr_to_string(g)
    None => ""
  }
  let mut shadow = 0.0
  for mp in r.metabolites.keys() {
    for i in 0..<model.metabolites.length() {
      if model.metabolites[i].id == mp && i < sol.shadow_prices.length() {
        shadow = shadow + sol.shadow_prices[i]
      }
    }
  }
  Some({
    reaction_id: r.id,
    name: r.name,
    flux: sol.fluxes[idx],
    lower_bound: r.lower_bound,
    upper_bound: r.upper_bound,
    shadow_price_of_metabolites: shadow,
    gpr,
    participants: cobra3_reaction_participants(r),
  })
}

///|
pub(all) struct Cobra3MetaboliteSummary {
  metabolite_id : String
  produced_by : Array[String]
  consumed_by : Array[String]
  shadow_price : Double
}

///|
pub fn cobra3_metabolite_summary(
  model : Cobra3Model,
  metabolite_id : String,
) -> Cobra3MetaboliteSummary? {
  let idx = cobra3_model_metabolite_index(model, metabolite_id)
  if idx < 0 {
    return None
  }
  let produced : Array[String] = Array(capacity=0)
  let consumed : Array[String] = Array(capacity=0)
  for r in model.reactions {
    match r.metabolites.get(metabolite_id) {
      Some(coeff) =>
        if coeff > 0.0 {
          produced.push(r.id)
        } else if coeff < 0.0 {
          consumed.push(r.id)
        }
      None => ()
    }
  }
  cobra3_sort_strings(produced)
  cobra3_sort_strings(consumed)
  let sol = cobra3_model_optimize(model)
  let shadow = if idx < sol.shadow_prices.length() {
    sol.shadow_prices[idx]
  } else {
    0.0
  }
  Some({
    metabolite_id,
    produced_by: produced,
    consumed_by: consumed,
    shadow_price: shadow,
  })
}

///|
pub(all) struct Cobra3SensitivityEntry {
  metabolite_id : String
  shadow_price : Double
}

///|
pub fn cobra3_sensitivity(model : Cobra3Model) -> Array[Cobra3SensitivityEntry] {
  let sol = cobra3_model_optimize(model)
  let out : Array[Cobra3SensitivityEntry] = Array(
    capacity=model.metabolites.length(),
  )
  for i in 0..<model.metabolites.length() {
    let sp = if i < sol.shadow_prices.length() {
      sol.shadow_prices[i]
    } else {
      0.0
    }
    out.push({ metabolite_id: model.metabolites[i].id, shadow_price: sp, })
  }
  for i in 1..<out.length() {
    let key = out[i]
    let mut j = i - 1
    while j >= 0 && out[j].shadow_price.abs() < key.shadow_price.abs() {
      out[j + 1] = out[j]
      j = j - 1
    }
    out[j + 1] = key
  }
  out
}

///|
pub(all) struct Cobra3ModelDiff {
  added_reactions : Array[String]
  removed_reactions : Array[String]
  changed_bounds : Array[String]
  changed_objective : Array[String]
  changed_stoichiometry : Array[String]
  added_metabolites : Array[String]
  removed_metabolites : Array[String]
}

///|

///|
fn cobra3_same_metabolites(
  a : Map[String, Double],
  b : Map[String, Double],
) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for k in a.keys() {
    match b.get(k) {
      Some(v) => if (a[k] - v).abs() > 1.0e-9 { return false }
      None => return false
    }
  }
  true
}

///|
pub fn cobra3_model_diff(a : Cobra3Model, b : Cobra3Model) -> Cobra3ModelDiff {
  let added_reactions : Array[String] = Array(capacity=0)
  let removed_reactions : Array[String] = Array(capacity=0)
  let changed_bounds : Array[String] = Array(capacity=0)
  let changed_objective : Array[String] = Array(capacity=0)
  let changed_stoichiometry : Array[String] = Array(capacity=0)
  for ra in a.reactions {
    match cobra3_model_get_reaction(b, ra.id) {
      None => removed_reactions.push(ra.id)
      Some(rb) => {
        if (ra.lower_bound - rb.lower_bound).abs() > 1.0e-9 ||
          (ra.upper_bound - rb.upper_bound).abs() > 1.0e-9 {
          changed_bounds.push(ra.id)
        }
        if (ra.objective - rb.objective).abs() > 1.0e-9 {
          changed_objective.push(ra.id)
        }
        if !cobra3_same_metabolites(ra.metabolites, rb.metabolites) {
          changed_stoichiometry.push(ra.id)
        }
      }
    }
  }
  for rb in b.reactions {
    if cobra3_model_reaction_index(a, rb.id) < 0 {
      added_reactions.push(rb.id)
    }
  }
  let added_metabolites : Array[String] = Array(capacity=0)
  let removed_metabolites : Array[String] = Array(capacity=0)
  for ma in a.metabolites {
    if cobra3_model_metabolite_index(b, ma.id) < 0 {
      removed_metabolites.push(ma.id)
    }
  }
  for mb in b.metabolites {
    if cobra3_model_metabolite_index(a, mb.id) < 0 {
      added_metabolites.push(mb.id)
    }
  }
  cobra3_sort_strings(added_reactions)
  cobra3_sort_strings(removed_reactions)
  cobra3_sort_strings(changed_bounds)
  cobra3_sort_strings(changed_objective)
  cobra3_sort_strings(changed_stoichiometry)
  cobra3_sort_strings(added_metabolites)
  cobra3_sort_strings(removed_metabolites)
  {
    added_reactions,
    removed_reactions,
    changed_bounds,
    changed_objective,
    changed_stoichiometry,
    added_metabolites,
    removed_metabolites,
  }
}
