///|
/// cobra3_flux.mbt — 在新对象模型上的通量分析.
///
/// - FBA / slim_optimize / shadow price / reduced cost
/// - pFBA（两阶段：先最优生长，再最小化总通量）
/// - FVA（fraction_of_optimum 真正生效）
/// - 单/双基因、单/双反应敲除（GPR 布尔求值）
/// - linear / quadratic MOMA、ROOM（真 MILP）、geometric FBA（QP）

///|
pub(all) struct Cobra3FluxResult {
  status : Cobra3Status
  objective_value : Double
  fluxes : Array[Double]
  shadow_prices : Array[Double]
  reduced_costs : Array[Double]
}

///|
fn cobra3_constraints_from_matrix(
  matrix : Array[Array[Double]],
  rhs : Array[Double],
) -> Array[Cobra3LinConstraint] {
  let out : Array[Cobra3LinConstraint] = Array(capacity=matrix.length())
  for i in 0..<matrix.length() {
    out.push({ coeffs: matrix[i], rhs: rhs[i], sense: Cobra3Sense::Equal, })
  }
  out
}

///|
fn cobra3_steady_state_constraints(
  model : Cobra3Model,
) -> Array[Cobra3LinConstraint] {
  let s = cobra3_model_stoichiometric_matrix(model)
  let rhs : Array[Double] = Array::make(s.length(), 0.0)
  cobra3_constraints_from_matrix(s, rhs)
}

///|

///|
fn cobra3_empty_flux(
  model : Cobra3Model,
  status : Cobra3Status,
) -> Cobra3FluxResult {
  {
    status,
    objective_value: 0.0,
    fluxes: Array::make(model.reactions.length(), 0.0),
    shadow_prices: Array::make(model.metabolites.length(), 0.0),
    reduced_costs: Array::make(model.reactions.length(), 0.0),
  }
}

///|
/// 标准 FBA：max c^T v，s.t. S v = 0，lb <= v <= ub.
pub fn cobra3_model_optimize(model : Cobra3Model) -> Cobra3FluxResult {
  let n = model.reactions.length()
  if n == 0 {
    return cobra3_empty_flux(model, Cobra3Status::Infeasible)
  }
  let constraints = cobra3_steady_state_constraints(model)
  let c = cobra3_model_objective_vector(model)
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  let sol = cobra3_solve_lp(constraints, c, lb, ub)
  {
    status: sol.status,
    objective_value: sol.objective_value,
    fluxes: sol.x,
    shadow_prices: sol.duals,
    reduced_costs: sol.reduced_costs,
  }
}

///|
pub fn cobra3_model_slim_optimize(model : Cobra3Model) -> Double {
  cobra3_model_optimize(model).objective_value
}

///|
/// 通量 -> id 映射.
pub fn cobra3_flux_map(
  model : Cobra3Model,
  result : Cobra3FluxResult,
) -> Map[String, Double] {
  let out : Map[String, Double] = Map([], capacity=model.reactions.length())
  for j in 0..<model.reactions.length() {
    out[model.reactions[j].id] = result.fluxes[j]
  }
  out
}

///|
/// pFBA：先求最大目标，再在目标不低于最优值 (1-eps) 下最小化 Σ|v|.
pub fn cobra3_pfba(model : Cobra3Model) -> Cobra3FluxResult {
  let n = model.reactions.length()
  if n == 0 {
    return cobra3_empty_flux(model, Cobra3Status::Infeasible)
  }
  let opt = cobra3_model_optimize(model)
  if !(opt.status is Optimal) {
    return opt
  }
  let constraints = cobra3_steady_state_constraints(model)
  // 目标约束：c·v >= opt*(1-1e-6)
  let c = cobra3_model_objective_vector(model)
  if opt.objective_value.abs() > 1.0e-9 {
    constraints.push({
      coeffs: c,
      rhs: opt.objective_value * (1.0 - 1.0e-6),
      sense: Cobra3Sense::GreaterEqual,
    })
  }
  // 变量：v (n) + t (n)，t_j >= |v_j|
  let nv = 2 * n
  let lb : Array[Double] = Array(capacity=nv)
  let ub : Array[Double] = Array(capacity=nv)
  let obj : Array[Double] = Array(capacity=nv)
  for j in 0..<n {
    lb.push(model.reactions[j].lower_bound)
    ub.push(model.reactions[j].upper_bound)
    obj.push(0.0)
  }
  for _ in 0..<n {
    lb.push(0.0)
    ub.push(1.0e18)
    obj.push(-1.0)
  }
  // 扩展约束：S v = 0 的 v 部分已有，t 部分补 0
  let full : Array[Cobra3LinConstraint] = Array(
    capacity=constraints.length() + 2 * n,
  )
  for con in constraints {
    let coeffs : Array[Double] = Array::make(nv, 0.0)
    for j in 0..<con.coeffs.length() {
      coerced_assign(coeffs, j, con.coeffs[j])
    }
    full.push({ coeffs, rhs: con.rhs, sense: con.sense, })
  }
  for j in 0..<n {
    // v_j - t_j <= 0
    let c1 : Array[Double] = Array::make(nv, 0.0)
    coerced_assign(c1, j, 1.0)
    coerced_assign(c1, n + j, -1.0)
    full.push({ coeffs: c1, rhs: 0.0, sense: Cobra3Sense::LessEqual, })
    // -v_j - t_j <= 0
    let c2 : Array[Double] = Array::make(nv, 0.0)
    coerced_assign(c2, j, -1.0)
    coerced_assign(c2, n + j, -1.0)
    full.push({ coeffs: c2, rhs: 0.0, sense: Cobra3Sense::LessEqual, })
  }
  let sol = cobra3_solve_lp(full, obj, lb, ub)
  if !(sol.status is Optimal) {
    return {
      status: sol.status,
      objective_value: opt.objective_value,
      fluxes: Array::make(n, 0.0),
      shadow_prices: Array::make(model.metabolites.length(), 0.0),
      reduced_costs: Array::make(n, 0.0),
    }
  }
  let fluxes : Array[Double] = Array(capacity=n)
  for j in 0..<n {
    fluxes.push(sol.x[j])
  }
  {
    status: sol.status,
    objective_value: opt.objective_value,
    fluxes,
    shadow_prices: Array::make(model.metabolites.length(), 0.0),
    reduced_costs: Array::make(n, 0.0),
  }
}

///|
fn coerced_assign(arr : Array[Double], idx : Int32, value : Double) -> Unit {
  arr[idx] = value
}

///|
pub(all) struct Cobra3FvaResult {
  reaction_id : String
  minimum : Double
  maximum : Double
}

///|
/// FVA：对每个反应求在稳态与目标 >= fraction*opt 下的通量范围.
/// processes 仅保留兼容参数（当前串行执行）；thresh 以下的结果归零.
/// 子问题非最优时该反应范围返回 0，调用方需自行区分不可行与真实 0.
pub fn cobra3_fva(
  model : Cobra3Model,
  fraction_of_optimum~ : Double = 1.0,
  thresh~ : Double = 1.0e-6,
  processes~ : Int32 = 1,
) -> Array[Cobra3FvaResult] {
  let n = model.reactions.length()
  let out : Array[Cobra3FvaResult] = Array(capacity=n)
  if n == 0 {
    return out
  }
  let _ = processes
  let opt = cobra3_model_optimize(model)
  if !(opt.status is Optimal) {
    for j in 0..<n {
      out.push({
        reaction_id: model.reactions[j].id,
        minimum: 0.0,
        maximum: 0.0,
      })
    }
    return out
  }
  let base = cobra3_steady_state_constraints(model)
  if fraction_of_optimum > 0.0 && opt.objective_value.abs() > 1.0e-9 {
    base.push({
      coeffs: cobra3_model_objective_vector(model),
      rhs: opt.objective_value * fraction_of_optimum,
      sense: Cobra3Sense::GreaterEqual,
    })
  }
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  for j in 0..<n {
    let c_min : Array[Double] = Array::make(n, 0.0)
    c_min[j] = -1.0
    let sol_min = cobra3_solve_lp(base, c_min, lb, ub)
    let c_max : Array[Double] = Array::make(n, 0.0)
    c_max[j] = 1.0
    let sol_max = cobra3_solve_lp(base, c_max, lb, ub)
    let mut mn = if sol_min.status is Optimal {
      -sol_min.objective_value
    } else {
      0.0
    }
    let mut mx = if sol_max.status is Optimal {
      sol_max.objective_value
    } else {
      0.0
    }
    if mn.abs() < thresh {
      mn = 0.0
    }
    if mx.abs() < thresh {
      mx = 0.0
    }
    out.push({ reaction_id: model.reactions[j].id, minimum: mn, maximum: mx, })
  }
  out
}

// ---------------------------------------------------------------------------
// 敲除
// ---------------------------------------------------------------------------

///|
pub(all) struct Cobra3DeletionResult {
  id : String
  growth : Double
  status : String
}

///|
fn cobra3_apply_gene_knockout(
  model : Cobra3Model,
  genes : Array[String],
) -> Cobra3Model {
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for r in model.reactions {
    let active = cobra3_reaction_active_in_model(model, r, genes)
    if active {
      reactions.push(r)
    } else {
      reactions.push({
        id: r.id,
        name: r.name,
        metabolites: r.metabolites,
        lower_bound: 0.0,
        upper_bound: 0.0,
        objective: r.objective,
        subsystem: r.subsystem,
        gpr: r.gpr,
        annotation: r.annotation,
      })
    }
  }
  {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: Array(capacity=0),
  }
}

///|
fn cobra3_apply_reaction_knockout(
  model : Cobra3Model,
  reaction_ids : Array[String],
) -> Cobra3Model {
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for r in model.reactions {
    if reaction_ids.contains(r.id) {
      reactions.push({
        id: r.id,
        name: r.name,
        metabolites: r.metabolites,
        lower_bound: 0.0,
        upper_bound: 0.0,
        objective: r.objective,
        subsystem: r.subsystem,
        gpr: r.gpr,
        annotation: r.annotation,
      })
    } else {
      reactions.push(r)
    }
  }
  {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: Array(capacity=0),
  }
}

///|
fn cobra3_gene_list(model : Cobra3Model) -> Array[String] {
  let ids : Array[String] = Array(capacity=model.genes.length())
  for g in model.genes {
    ids.push(g.id)
  }
  cobra3_sort_strings(ids)
  ids
}

///|
fn cobra3_status_text(status : Cobra3Status) -> String {
  match status {
    Optimal => "optimal"
    Infeasible => "infeasible"
    Unbounded => "unbounded"
    IterLimit => "iter_limit"
  }
}

///|
pub fn cobra3_single_gene_deletion(
  model : Cobra3Model,
) -> Array[Cobra3DeletionResult] {
  let out : Array[Cobra3DeletionResult] = Array(capacity=model.genes.length())
  for gene in cobra3_gene_list(model) {
    let ko = cobra3_apply_gene_knockout(model, [gene])
    let sol = cobra3_model_optimize(ko)
    let status = cobra3_status_text(sol.status)
    out.push({ id: gene, growth: sol.objective_value, status, })
  }
  out
}

///|
pub fn cobra3_single_reaction_deletion(
  model : Cobra3Model,
) -> Array[Cobra3DeletionResult] {
  let out : Array[Cobra3DeletionResult] = Array(
    capacity=model.reactions.length(),
  )
  for r in model.reactions {
    let ko = cobra3_apply_reaction_knockout(model, [r.id])
    let sol = cobra3_model_optimize(ko)
    let status = cobra3_status_text(sol.status)
    out.push({ id: r.id, growth: sol.objective_value, status, })
  }
  out
}

///|
pub fn cobra3_double_gene_deletion(
  model : Cobra3Model,
) -> Array[Cobra3DeletionResult] {
  let genes = cobra3_gene_list(model)
  let out : Array[Cobra3DeletionResult] = Array(
    capacity=genes.length() * genes.length(),
  )
  for i in 0..<genes.length() {
    for j in (i + 1)..<genes.length() {
      let ko = cobra3_apply_gene_knockout(model, [genes[i], genes[j]])
      let sol = cobra3_model_optimize(ko)
      let status = cobra3_status_text(sol.status)
      out.push({
        id: genes[i] + "|" + genes[j],
        growth: sol.objective_value,
        status,
      })
    }
  }
  out
}

///|
pub fn cobra3_double_reaction_deletion(
  model : Cobra3Model,
) -> Array[Cobra3DeletionResult] {
  let ids : Array[String] = Array(capacity=model.reactions.length())
  for r in model.reactions {
    ids.push(r.id)
  }
  let out : Array[Cobra3DeletionResult] = Array(
    capacity=ids.length() * ids.length(),
  )
  for i in 0..<ids.length() {
    for j in (i + 1)..<ids.length() {
      let ko = cobra3_apply_reaction_knockout(model, [ids[i], ids[j]])
      let sol = cobra3_model_optimize(ko)
      let status = cobra3_status_text(sol.status)
      out.push({
        id: ids[i] + "|" + ids[j],
        growth: sol.objective_value,
        status,
      })
    }
  }
  out
}

// ---------------------------------------------------------------------------
// MOMA / ROOM / geometric FBA
// ---------------------------------------------------------------------------

///|
pub(all) struct Cobra3MomaResult {
  growth : Double
  fluxes : Array[Double]
  distance : Double
  status : Cobra3Status
}

///|
/// MOMA：给定参考通量 v_wt 与敲除模型，求最接近参考的可行通量.
/// quadratic=true 用 QP（L2）；否则线性化（L1）.
pub fn cobra3_moma(
  model : Cobra3Model,
  knockout_genes : Array[String],
  quadratic~ : Bool = true,
) -> Cobra3MomaResult {
  let wt = cobra3_model_optimize(model)
  let ko = cobra3_apply_gene_knockout(model, knockout_genes)
  let n = model.reactions.length()
  if !(wt.status is Optimal) || n == 0 {
    return {
      growth: 0.0,
      fluxes: Array::make(n, 0.0),
      distance: 0.0,
      status: Cobra3Status::Infeasible,
    }
  }
  let lb = cobra3_model_lower_bounds(ko)
  let ub = cobra3_model_upper_bounds(ko)
  let s = cobra3_model_stoichiometric_matrix(ko)
  let zero : Array[Double] = Array::make(s.length(), 0.0)
  if quadratic {
    let q : Array[Array[Double]] = Array(capacity=n)
    for i in 0..<n {
      let row : Array[Double] = Array::make(n, 0.0)
      row[i] = 2.0
      q.push(row)
    }
    let c : Array[Double] = Array(capacity=n)
    for j in 0..<n {
      c.push(-2.0 * wt.fluxes[j])
    }
    let sol = cobra3_solve_qp(q, c, s, zero, lb, ub)
    let mut dist = 0.0
    for j in 0..<n {
      dist = dist +
        sol.x[j] * sol.x[j] -
        2.0 * wt.fluxes[j] * sol.x[j] +
        wt.fluxes[j] * wt.fluxes[j]
    }
    let growth = match cobra3_model_biomass_reaction(ko) {
      Some(id) =>
        match cobra3_model_reaction_index(ko, id) {
          i => if i >= 0 { sol.x[i] } else { 0.0 }
        }
      None => 0.0
    }
    { growth, fluxes: sol.x, distance: dist, status: sol.status, }
  } else {
    // L1: 变量 v(n) + d(n)，d_j >= |v_j - wt_j|
    let nv = 2 * n
    let lb2 : Array[Double] = Array(capacity=nv)
    let ub2 : Array[Double] = Array(capacity=nv)
    let obj : Array[Double] = Array(capacity=nv)
    for j in 0..<n {
      lb2.push(lb[j])
      ub2.push(ub[j])
      obj.push(0.0)
    }
    for _ in 0..<n {
      lb2.push(0.0)
      ub2.push(1.0e18)
      obj.push(-1.0)
    }
    let cons : Array[Cobra3LinConstraint] = Array(capacity=s.length() + 2 * n)
    for i in 0..<s.length() {
      let coeffs : Array[Double] = Array::make(nv, 0.0)
      for j in 0..<n {
        coeffs[j] = s[i][j]
      }
      cons.push({ coeffs, rhs: 0.0, sense: Cobra3Sense::Equal, })
    }
    for j in 0..<n {
      let c1 : Array[Double] = Array::make(nv, 0.0)
      c1[j] = 1.0
      c1[n + j] = -1.0
      cons.push({
        coeffs: c1,
        rhs: wt.fluxes[j],
        sense: Cobra3Sense::LessEqual,
      })
      let c2 : Array[Double] = Array::make(nv, 0.0)
      c2[j] = -1.0
      c2[n + j] = -1.0
      cons.push({
        coeffs: c2,
        rhs: -wt.fluxes[j],
        sense: Cobra3Sense::LessEqual,
      })
    }
    let sol = cobra3_solve_lp(cons, obj, lb2, ub2)
    let fluxes : Array[Double] = Array(capacity=n)
    for j in 0..<n {
      fluxes.push(sol.x[j])
    }
    let growth = match cobra3_model_biomass_reaction(ko) {
      Some(id) =>
        match cobra3_model_reaction_index(ko, id) {
          i => if i >= 0 { fluxes[i] } else { 0.0 }
        }
      None => 0.0
    }
    { growth, fluxes, distance: -sol.objective_value, status: sol.status, }
  }
}

///|
pub(all) struct Cobra3RoomResult {
  growth : Double
  fluxes : Array[Double]
  n_state_changes : Int32
  status : Cobra3Status
}

///|
/// ROOM：MILP，最小化相对参考解的"开关状态"改变数.
pub fn cobra3_room(
  model : Cobra3Model,
  knockout_genes : Array[String],
  delta~ : Double = 0.03,
  epsilon~ : Double = 0.001,
  max_nodes~ : Int32 = 2000,
) -> Cobra3RoomResult {
  let wt = cobra3_model_optimize(model)
  let ko = cobra3_apply_gene_knockout(model, knockout_genes)
  let n = model.reactions.length()
  if !(wt.status is Optimal) || n == 0 {
    return {
      growth: 0.0,
      fluxes: Array::make(n, 0.0),
      n_state_changes: 0,
      status: Cobra3Status::Infeasible,
    }
  }
  let s = cobra3_model_stoichiometric_matrix(ko)
  let lb0 = cobra3_model_lower_bounds(ko)
  let ub0 = cobra3_model_upper_bounds(ko)
  let nv = 2 * n
  let lb : Array[Double] = Array(capacity=nv)
  let ub : Array[Double] = Array(capacity=nv)
  for j in 0..<n {
    lb.push(lb0[j])
    ub.push(ub0[j])
  }
  for _ in 0..<n {
    lb.push(0.0)
    ub.push(1.0)
  }
  let obj : Array[Double] = Array::make(nv, 0.0)
  for j in 0..<n {
    obj[n + j] = -1.0
  }
  let cons : Array[Cobra3LinConstraint] = Array(capacity=s.length() + 2 * n)
  for i in 0..<s.length() {
    let coeffs : Array[Double] = Array::make(nv, 0.0)
    for j in 0..<n {
      coeffs[j] = s[i][j]
    }
    cons.push({ coeffs, rhs: 0.0, sense: Cobra3Sense::Equal, })
  }
  for j in 0..<n {
    let slack = delta * wt.fluxes[j].abs() + epsilon
    let scale = if lb0[j].abs() > ub0[j].abs() {
      lb0[j].abs()
    } else {
      ub0[j].abs()
    }
    let max_flux = if wt.fluxes[j].abs() > scale {
      wt.fluxes[j].abs()
    } else {
      scale
    }
    let m_j = 2.0 * (max_flux + slack) + 1.0
    let c1 : Array[Double] = Array::make(nv, 0.0)
    c1[j] = 1.0
    c1[n + j] = -m_j
    cons.push({
      coeffs: c1,
      rhs: wt.fluxes[j] + slack,
      sense: Cobra3Sense::LessEqual,
    })
    let c2 : Array[Double] = Array::make(nv, 0.0)
    c2[j] = -1.0
    c2[n + j] = -m_j
    cons.push({
      coeffs: c2,
      rhs: -wt.fluxes[j] + slack,
      sense: Cobra3Sense::LessEqual,
    })
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
  let growth = match cobra3_model_biomass_reaction(ko) {
    Some(id) =>
      match cobra3_model_reaction_index(ko, id) {
        i => if i >= 0 { fluxes[i] } else { 0.0 }
      }
    None => 0.0
  }
  {
    growth,
    fluxes,
    n_state_changes: (-sol.objective_value + 0.5).to_int(),
    status: sol.status,
  }
}

///|
/// geometric FBA：保持至少 99% 最优生长，最小化 Σ v^2.
pub fn cobra3_geometric_fba(model : Cobra3Model) -> Cobra3FluxResult {
  let n = model.reactions.length()
  if n == 0 {
    return cobra3_empty_flux(model, Cobra3Status::Infeasible)
  }
  let opt = cobra3_model_optimize(model)
  if !(opt.status is Optimal) {
    return opt
  }
  let s = cobra3_model_stoichiometric_matrix(model)
  let zero : Array[Double] = Array::make(s.length(), 0.0)
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  // 用生长反应下界表达"至少 99% 最优生长"
  match cobra3_model_biomass_reaction(model) {
    Some(id) => {
      let idx = cobra3_model_reaction_index(model, id)
      if idx >= 0 {
        let floor = opt.objective_value * 0.99
        if lb[idx] < floor {
          lb[idx] = floor
        }
      }
    }
    None => ()
  }
  let q : Array[Array[Double]] = Array(capacity=n)
  for i in 0..<n {
    let row : Array[Double] = Array::make(n, 0.0)
    row[i] = 2.0
    q.push(row)
  }
  let c : Array[Double] = Array::make(n, 0.0)
  let eq : Array[Array[Double]] = Array(capacity=s.length())
  for i in 0..<s.length() {
    eq.push(s[i])
  }
  let sol = cobra3_solve_qp(q, c, eq, zero, lb, ub)
  {
    status: sol.status,
    objective_value: opt.objective_value,
    fluxes: sol.x,
    shadow_prices: Array::make(model.metabolites.length(), 0.0),
    reduced_costs: Array::make(n, 0.0),
  }
}
