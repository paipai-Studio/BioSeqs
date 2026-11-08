///|
/// cobra3_design.mbt — 模型设计与采样.
///
/// - gapfill（MILP 最小添加反应数）
/// - minimal_medium（MILP 最小开启交换数）
/// - production_envelope（目标通量网格）
/// - MCS-lite（单/双反应最小割集）
/// - 通量采样（hit-and-run）+ sample_flux 统计 + 收敛诊断
/// - FastCC、loopless-lite（最小化内部循环通量）

///|
pub(all) struct Cobra3GapFillResult {
  added_reactions : Array[String]
  objective_value : Double
  status : Cobra3Status
}

///|
fn cobra3_merge_universal(
  model : Cobra3Model,
  universal : Cobra3Model,
) -> (Cobra3Model, Array[String]) {
  let mut merged = model
  let candidates : Array[String] = Array(capacity=0)
  for r in universal.reactions {
    if cobra3_model_reaction_index(merged, r.id) >= 0 {
      continue
    }
    merged = cobra3_model_add_reaction(merged, r)
    candidates.push(r.id)
  }
  (merged, candidates)
}

///|
/// gapfill：在 universal 模型中挑选最少反应，使合并模型目标 >= lower_bound.
pub fn cobra3_gapfill(
  model : Cobra3Model,
  universal : Cobra3Model,
  lower_bound~ : Double = 0.05,
  max_nodes~ : Int32 = 3000,
) -> Cobra3GapFillResult {
  let base = cobra3_model_optimize(model)
  if base.status is Optimal && base.objective_value >= lower_bound {
    return {
      added_reactions: Array(capacity=0),
      objective_value: base.objective_value,
      status: Cobra3Status::Optimal,
    }
  }
  let (merged, candidates) = cobra3_merge_universal(model, universal)
  let n = merged.reactions.length()
  let nc = candidates.length()
  if nc == 0 {
    return {
      added_reactions: Array(capacity=0),
      objective_value: base.objective_value,
      status: Cobra3Status::Infeasible,
    }
  }
  let nv = n + nc
  let s = cobra3_model_stoichiometric_matrix(merged)
  let cons : Array[Cobra3LinConstraint] = Array(
    capacity=s.length() + 1 + 2 * nc,
  )
  for i in 0..<s.length() {
    let coeffs : Array[Double] = Array::make(nv, 0.0)
    for j in 0..<n {
      coeffs[j] = s[i][j]
    }
    cons.push({ coeffs, rhs: 0.0, sense: Cobra3Sense::Equal, })
  }
  // 生长下界
  let obj_vec = cobra3_model_objective_vector(merged)
  let grow : Array[Double] = Array::make(nv, 0.0)
  for j in 0..<n {
    grow[j] = obj_vec[j]
  }
  cons.push({
    coeffs: grow,
    rhs: lower_bound,
    sense: Cobra3Sense::GreaterEqual,
  })
  let lb : Array[Double] = Array::make(nv, 0.0)
  let ub : Array[Double] = Array::make(nv, 0.0)
  for j in 0..<n {
    lb[j] = merged.reactions[j].lower_bound
    ub[j] = merged.reactions[j].upper_bound
  }
  for k in 0..<nc {
    let j = n + k
    lb[j] = 0.0
    ub[j] = 1.0
    let reaction_idx = cobra3_model_reaction_index(merged, candidates[k])
    let rlb = merged.reactions[reaction_idx].lower_bound
    let rub = merged.reactions[reaction_idx].upper_bound
    let scale = if rlb.abs() > rub.abs() { rlb.abs() } else { rub.abs() }
    let m = 2.0 * scale + 1.0
    // -M*y <= v <= M*y
    let c1 : Array[Double] = Array::make(nv, 0.0)
    c1[reaction_idx] = 1.0
    c1[j] = -m
    cons.push({ coeffs: c1, rhs: 0.0, sense: Cobra3Sense::LessEqual, })
    let c2 : Array[Double] = Array::make(nv, 0.0)
    c2[reaction_idx] = -1.0
    c2[j] = -m
    cons.push({ coeffs: c2, rhs: 0.0, sense: Cobra3Sense::LessEqual, })
  }
  let obj : Array[Double] = Array::make(nv, 0.0)
  for k in 0..<nc {
    obj[n + k] = -1.0
  }
  let integer : Array[Int32] = Array(capacity=nc)
  for k in 0..<nc {
    integer.push(n + k)
  }
  let sol = cobra3_solve_milp(cons, obj, lb, ub, integer, max_nodes~)
  if !(sol.status is Optimal) && !(sol.status is Feasible) {
    return {
      added_reactions: Array(capacity=0),
      objective_value: base.objective_value,
      status: sol.status,
    }
  }
  let added : Array[String] = Array(capacity=0)
  for k in 0..<nc {
    if sol.x[n + k] > 0.5 {
      added.push(candidates[k])
    }
  }
  // 用实际加入的反应重建模型并重解，避免返回固定的 lower_bound
  let mut completed = model
  for id in added {
    match cobra3_model_get_reaction(universal, id) {
      Some(r) => completed = cobra3_model_add_reaction(completed, r)
      None => ()
    }
  }
  let checked = cobra3_model_optimize(completed)
  let growth = if checked.status is Optimal {
    checked.objective_value
  } else {
    0.0
  }
  { added_reactions: added, objective_value: growth, status: sol.status, }
}

///|
pub(all) struct Cobra3MinimalMediumResult {
  medium : Map[String, Double]
  n_components : Int32
  objective_value : Double
  status : Cobra3Status
}

///|
/// minimal_medium：MILP 最小化开放的 exchange 数，同时保持生长 >= min_growth.
pub fn cobra3_minimal_medium(
  model : Cobra3Model,
  min_growth~ : Double = 0.1,
  max_import~ : Double = 10.0,
  max_nodes~ : Int32 = 3000,
) -> Cobra3MinimalMediumResult {
  let exchanges = cobra3_model_exchange_reactions(model)
  let n = model.reactions.length()
  let ne = exchanges.length()
  let nv = n + ne
  let s = cobra3_model_stoichiometric_matrix(model)
  let cons : Array[Cobra3LinConstraint] = Array(capacity=s.length() + 1 + ne)
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
  for k in 0..<ne {
    let ex_idx = cobra3_model_reaction_index(model, exchanges[k])
    lb[n + k] = 0.0
    ub[n + k] = 1.0
    // v_ex >= -max_import * y
    let c1 : Array[Double] = Array::make(nv, 0.0)
    c1[ex_idx] = -1.0
    c1[n + k] = -max_import
    cons.push({ coeffs: c1, rhs: 0.0, sense: Cobra3Sense::LessEqual, })
  }
  let obj : Array[Double] = Array::make(nv, 0.0)
  for k in 0..<ne {
    obj[n + k] = -1.0
  }
  let integer : Array[Int32] = Array(capacity=ne)
  for k in 0..<ne {
    integer.push(n + k)
  }
  let sol = cobra3_solve_milp(cons, obj, lb, ub, integer, max_nodes~)
  let medium : Map[String, Double] = Map([], capacity=ne)
  if !(sol.status is Optimal) && !(sol.status is Feasible) {
    return {
      medium,
      n_components: 0,
      objective_value: 0.0,
      status: sol.status,
    }
  }
  let mut n_components = 0
  for k in 0..<ne {
    if sol.x[n + k] > 0.5 {
      let ex_idx = cobra3_model_reaction_index(model, exchanges[k])
      let uptake = -sol.x[ex_idx]
      if uptake > 1.0e-9 {
        medium[exchanges[k]] = uptake
        n_components = n_components + 1
      }
    }
  }
  let mut growth = 0.0
  for j in 0..<n {
    growth = growth + obj_vec[j] * sol.x[j]
  }
  { medium, n_components, objective_value: growth, status: sol.status, }
}

///|
pub(all) struct Cobra3EnvelopePoint {
  target_flux : Double
  objective_max : Double
  objective_min : Double
  yield_max : Double
  yield_min : Double
}

///|
pub fn cobra3_production_envelope(
  model : Cobra3Model,
  target_id : String,
  points~ : Int32 = 20,
  carbon_source~ : String = "",
) -> Array[Cobra3EnvelopePoint] {
  let target_idx = cobra3_model_reaction_index(model, target_id)
  let out : Array[Cobra3EnvelopePoint] = Array(capacity=points)
  if target_idx < 0 || points <= 0 {
    return out
  }
  let fva = cobra3_fva(model, fraction_of_optimum=0.0)
  let mut lo = 0.0
  let mut hi = 0.0
  for f in fva {
    if f.reaction_id == target_id {
      lo = f.minimum
      hi = f.maximum
    }
  }
  if hi < lo {
    return out
  }
  let n = model.reactions.length()
  let carbon_idx = if carbon_source.length() == 0 {
    -1
  } else {
    cobra3_model_reaction_index(model, carbon_source)
  }
  let base = cobra3_steady_state_constraints(model)
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
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
    let fix : Array[Double] = Array::make(n, 0.0)
    fix[target_idx] = 1.0
    cons.push({ coeffs: fix, rhs: target_flux, sense: Cobra3Sense::Equal, })
    let c_max : Array[Double] = Array::make(n, 0.0)
    let c_min : Array[Double] = Array::make(n, 0.0)
    for j in 0..<n {
      c_max[j] = model.reactions[j].objective
      c_min[j] = -model.reactions[j].objective
    }
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
    let mut yield_max = 0.0
    let mut yield_min = 0.0
    if carbon_idx >= 0 {
      let uptake = if sol_max.status is Optimal {
        -sol_max.x[carbon_idx]
      } else {
        0.0
      }
      if uptake > 1.0e-9 {
        yield_max = obj_max / uptake
      }
      let uptake_min = if sol_min.status is Optimal {
        -sol_min.x[carbon_idx]
      } else {
        0.0
      }
      if uptake_min > 1.0e-9 {
        yield_min = obj_min / uptake_min
      }
    }
    out.push({
      target_flux,
      objective_max: obj_max,
      objective_min: obj_min,
      yield_max,
      yield_min,
    })
  }
  out
}

///|
pub(all) struct Cobra3McsResult {
  reaction_ids : Array[String]
  growth : Double
  size : Int32
}

///|
/// 判断 `ids` 是否包含 `sets` 中任一集合（即该组合不是极小割集）.
fn cobra3_mcs_contains_smaller(
  ids : Array[String],
  sets : Array[Array[String]],
) -> Bool {
  for set in sets {
    let mut contained = true
    for id in set {
      if !ids.contains(id) {
        contained = false
        break
      }
    }
    if contained {
      return true
    }
  }
  false
}

///|
/// 递归枚举尺寸为 `k` 的反应删除组合；组合被测试一次后 `budget[0]` 递减，
/// 归零即停止，避免大模型上组合爆炸.
fn cobra3_mcs_extend(
  model : Cobra3Model,
  threshold : Double,
  candidates : Array[Int32],
  k : Int32,
  start : Int32,
  cur : Array[String],
  found : Array[Array[String]],
  out : Array[Cobra3McsResult],
  budget : Array[Int32],
) -> Unit {
  let total = candidates.length()
  let need = k - cur.length()
  let mut i = start
  while i < total && total - i >= need && budget[0] > 0 {
    cur.push(model.reactions[candidates[i]].id)
    if cur.length() == k {
      budget[0] = budget[0] - 1
      // 包含已找到的更小割集 => 不是极小割集
      if !cobra3_mcs_contains_smaller(cur, found) {
        let ko = cobra3_knockout_reactions(model, cur)
        let sol = cobra3_model_optimize(ko)
        if sol.status is Optimal && sol.objective_value <= threshold {
          let ids : Array[String] = Array(capacity=cur.length())
          for id in cur {
            ids.push(id)
          }
          found.push(ids)
          out.push({ reaction_ids: ids, growth: sol.objective_value, size: k, })
        }
      }
    } else {
      cobra3_mcs_extend(
        model,
        threshold,
        candidates,
        k,
        i + 1,
        cur,
        found,
        out,
        budget,
      )
    }
    ignore(cur.pop())
    i = i + 1
  }
}

///|
/// 最小割集：按尺寸 1..max_size 递增枚举反应删除组合，跳过任何包含
/// 已知更小割集的组合，保证返回值都是极小割集（旧实现只支持 size<=2，
/// 且遇到单割集就漏掉不相交的双割集）。`max_candidates` 限制组合测试
/// 次数，达到上限即停止枚举。`exclude` 中的反应不会作为敲除候选
/// （例如保护 exchange / biomass 反应，对应 COBRApy 的 exclude 参数）。
pub fn cobra3_minimal_cut_sets(
  model : Cobra3Model,
  threshold~ : Double = 0.01,
  max_size~ : Int32 = 2,
  max_candidates~ : Int32 = 20000,
  exclude~ : Array[String] = [],
) -> Array[Cobra3McsResult] {
  let out : Array[Cobra3McsResult] = Array(capacity=0)
  let base = cobra3_model_optimize(model)
  if !(base.status is Optimal) || base.objective_value <= threshold {
    return out
  }
  let n = model.reactions.length()
  let candidates : Array[Int32] = Array(capacity=n)
  for j in 0..<n {
    if !exclude.contains(model.reactions[j].id) {
      candidates.push(j)
    }
  }
  let m = candidates.length()
  if m == 0 {
    return out
  }
  let found : Array[Array[String]] = Array(capacity=0)
  let budget : Array[Int32] = [max_candidates]
  let mut k = 1
  while k <= max_size && k <= m && budget[0] > 0 {
    let cur : Array[String] = Array(capacity=k)
    cobra3_mcs_extend(
      model, threshold, candidates, k, 0, cur, found, out, budget,
    )
    k = k + 1
  }
  out
}

///|
fn cobra3_knockout_reactions(
  model : Cobra3Model,
  ids : Array[String],
) -> Cobra3Model {
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for r in model.reactions {
    if ids.contains(r.id) {
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

// ---------------------------------------------------------------------------
// 采样
// ---------------------------------------------------------------------------

///|
/// 化学计量矩阵零空间基（返回 basis x n_reactions）.
pub fn cobra3_nullspace_basis(model : Cobra3Model) -> Array[Array[Double]] {
  let s = cobra3_model_stoichiometric_matrix(model)
  let m = s.length()
  let n = if m == 0 { 0 } else { s[0].length() }
  if n == 0 {
    return Array(capacity=0)
  }
  let a : Array[Array[Double]] = Array(capacity=m)
  for i in 0..<m {
    let row : Array[Double] = Array(capacity=n)
    for j in 0..<n {
      row.push(s[i][j])
    }
    a.push(row)
  }
  let pivot_cols : Array[Int32] = Array(capacity=0)
  let mut r = 0
  for col in 0..<n {
    let mut piv = r
    for i in r..<m {
      if a[i][col].abs() > a[piv][col].abs() {
        piv = i
      }
    }
    if r >= m || a[piv][col].abs() < 1.0e-10 {
      continue
    }
    if piv != r {
      let tmp = a[r]
      a[r] = a[piv]
      a[piv] = tmp
    }
    let d = a[r][col]
    for j in 0..<n {
      a[r][j] = a[r][j] / d
    }
    for i in 0..<m {
      if i != r {
        let factor = a[i][col]
        if factor != 0.0 {
          for j in 0..<n {
            a[i][j] = a[i][j] - factor * a[r][j]
          }
        }
      }
    }
    pivot_cols.push(col)
    r = r + 1
    if r == m {
      break
    }
  }
  let free : Array[Int32] = Array(capacity=n)
  for j in 0..<n {
    if !pivot_cols.contains(j) {
      free.push(j)
    }
  }
  let basis : Array[Array[Double]] = Array(capacity=free.length())
  for f in free {
    let vec : Array[Double] = Array::make(n, 0.0)
    vec[f] = 1.0
    for k in 0..<pivot_cols.length() {
      vec[pivot_cols[k]] = -a[k][f]
      if vec[pivot_cols[k]].abs() < 1.0e-12 {
        vec[pivot_cols[k]] = 0.0
      }
    }
    basis.push(vec)
  }
  basis
}

///|
fn cobra3_lcg_next(state : Int32) -> Int32 {
  let s = (state * 1103515245 + 12345) % 2147483647
  if s < 0 {
    -s
  } else {
    s
  }
}

///|
/// hit-and-run 通量采样；返回 samples x reactions.
pub fn cobra3_sample(
  model : Cobra3Model,
  n_samples~ : Int32 = 1000,
  seed~ : Int32 = 42,
  warmup~ : Int32 = 100,
) -> Array[Array[Double]] {
  let n = model.reactions.length()
  let out : Array[Array[Double]] = Array(capacity=n_samples)
  if n == 0 || n_samples <= 0 {
    return out
  }
  let zero_c : Array[Double] = Array::make(n, 0.0)
  let base = cobra3_solve_lp_eq(
    cobra3_model_stoichiometric_matrix(model),
    Array::make(model.metabolites.length(), 0.0),
    zero_c,
    cobra3_model_lower_bounds(model),
    cobra3_model_upper_bounds(model),
  )
  if !(base.status is Optimal) {
    return out
  }
  let x = cobra3_bounds_copy(base.x)
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  let null_basis = cobra3_nullspace_basis(model)
  if null_basis.length() == 0 {
    return out
  }
  let mut state = seed
  let mut step = 0
  let target_steps = warmup + n_samples
  while step < target_steps {
    // 随机方向：零空间基的随机线性组合
    let direction : Array[Double] = Array::make(n, 0.0)
    for b in null_basis {
      state = cobra3_lcg_next(state)
      let coef = state.to_double() / 2147483647.0 - 0.5
      for j in 0..<n {
        direction[j] = direction[j] + coef * b[j]
      }
    }
    // 归一化方向，避免不同基向量尺度主导步长
    let mut norm = 0.0
    for j in 0..<n {
      norm = norm + direction[j] * direction[j]
    }
    norm = @math.pow(norm, 0.5)
    if norm > 1.0e-12 {
      for j in 0..<n {
        direction[j] = direction[j] / norm
      }
    }
    // 计算可行步长区间
    let mut t_min = -1.0e18
    let mut t_max = 1.0e18
    for j in 0..<n {
      let d = direction[j]
      if d > 1.0e-12 {
        let hi = (ub[j] - x[j]) / d
        let lo = (lb[j] - x[j]) / d
        if hi < t_max {
          t_max = hi
        }
        if lo > t_min {
          t_min = lo
        }
      } else if d < -1.0e-12 {
        let hi = (lb[j] - x[j]) / d
        let lo = (ub[j] - x[j]) / d
        if hi < t_max {
          t_max = hi
        }
        if lo > t_min {
          t_min = lo
        }
      }
    }
    if t_max >= t_min {
      state = cobra3_lcg_next(state)
      let u = state.to_double() / 2147483647.0
      let t = t_min + u * (t_max - t_min)
      for j in 0..<n {
        let mut v = x[j] + t * direction[j]
        if v < lb[j] {
          v = lb[j]
        }
        if v > ub[j] {
          v = ub[j]
        }
        x[j] = v
      }
    }
    if step >= warmup && out.length() < n_samples {
      let sample : Array[Double] = Array(capacity=n)
      for j in 0..<n {
        sample.push(x[j])
      }
      out.push(sample)
    }
    step = step + 1
  }
  out
}

///|
/// 独立多链采样：每条链用不同 seed，返回 chain × sample × reaction.
pub fn cobra3_sample_chains(
  model : Cobra3Model,
  n_chains~ : Int32 = 4,
  n_samples~ : Int32 = 1000,
  seed~ : Int32 = 42,
  warmup~ : Int32 = 100,
) -> Array[Array[Array[Double]]] {
  let out : Array[Array[Array[Double]]] = Array(capacity=n_chains)
  for k in 0..<n_chains {
    out.push(cobra3_sample(model, n_samples~, seed=seed + k * 7919, warmup~))
  }
  out
}

///|
pub(all) struct Cobra3FluxDistribution {
  reaction_id : String
  mean : Double
  std : Double
  minimum : Double
  maximum : Double
  q05 : Double
  q95 : Double
}

///|
/// 从采样结果计算每个反应的通量分布.
pub fn cobra3_sample_flux(
  model : Cobra3Model,
  samples : Array[Array[Double]],
) -> Array[Cobra3FluxDistribution] {
  let out : Array[Cobra3FluxDistribution] = Array(
    capacity=model.reactions.length(),
  )
  if samples.length() == 0 {
    return out
  }
  for j in 0..<model.reactions.length() {
    let vals : Array[Double] = Array(capacity=samples.length())
    for s in samples {
      vals.push(s[j])
    }
    cobra3_sort_doubles(vals)
    let mut sum = 0.0
    for v in vals {
      sum = sum + v
    }
    let mean = sum / vals.length().to_double()
    let mut var_sum = 0.0
    for v in vals {
      var_sum = var_sum + (v - mean) * (v - mean)
    }
    let sd = @math.pow(var_sum / vals.length().to_double(), 0.5)
    let last = (vals.length() - 1).to_double()
    let q05 = vals[(last * 0.05).to_int()]
    let q95 = vals[(last * 0.95).to_int()]
    out.push({
      reaction_id: model.reactions[j].id,
      mean,
      std: sd,
      minimum: vals[0],
      maximum: vals[vals.length() - 1],
      q05,
      q95,
    })
  }
  out
}

///|
pub(all) struct Cobra3SamplingDiagnostics {
  r_hat : Double
  geweke_z : Double
  effective_sample_size : Double
}

///|
/// 多链采样的 R-hat 与单链 Geweke z / 有效样本数（简化诊断）.
pub fn cobra3_sampling_diagnostics(
  samples : Array[Double],
) -> Cobra3SamplingDiagnostics {
  let n = samples.length()
  if n < 10 {
    return { r_hat: 1.0, geweke_z: 0.0, effective_sample_size: 0.0, }
  }
  let half = n / 2
  let mean_all = cobra3_mean(samples, 0, n)
  let mean_a = cobra3_mean(samples, 0, half)
  let mean_b = cobra3_mean(samples, half, n)
  let var_a = cobra3_var(samples, 0, half, mean_a)
  let var_b = cobra3_var(samples, half, n, mean_b)
  let w = 0.5 * (var_a + var_b)
  let b = half.to_double() * (mean_a - mean_all) * (mean_a - mean_all) * 2.0
  let var_hat = (half.to_double() - 1.0) / half.to_double() * w +
    b / half.to_double()
  let r_hat = if w > 1.0e-15 { @math.pow(var_hat / w, 0.5) } else { 1.0 }
  let slit = n / 4
  let mean_front = cobra3_mean(samples, 0, slit)
  let mean_back = cobra3_mean(samples, n - slit, n)
  let var_front = cobra3_var(samples, 0, slit, mean_front)
  let var_back = cobra3_var(samples, n - slit, n, mean_back)
  let denom = @math.pow(
    var_front / slit.to_double() + var_back / slit.to_double(),
    0.5,
  )
  let geweke = if denom > 1.0e-15 {
    (mean_front - mean_back) / denom
  } else {
    0.0
  }
  // 用 lag-1 自相关估计 ESS
  let mut lag1 = 0.0
  let mut denom_ac = 0.0
  for i in 0..<n {
    denom_ac = denom_ac + (samples[i] - mean_all) * (samples[i] - mean_all)
  }
  let limit = n - 1
  for i in 0..<limit {
    lag1 = lag1 + (samples[i] - mean_all) * (samples[i + 1] - mean_all)
  }
  let rho = if denom_ac > 1.0e-15 { lag1 / denom_ac } else { 0.0 }
  let ess = if rho.abs() < 0.999 {
    n.to_double() * (1.0 - rho) / (1.0 + rho)
  } else {
    1.0
  }
  { r_hat, geweke_z: geweke, effective_sample_size: ess, }
}

///|
/// 多链诊断：R-hat 用链间/链内方差；Geweke z 取各链平均，ESS 取各链之和.
pub fn cobra3_sampling_diagnostics_chains(
  chains : Array[Array[Double]],
) -> Cobra3SamplingDiagnostics {
  let k = chains.length()
  if k == 0 {
    return { r_hat: 1.0, geweke_z: 0.0, effective_sample_size: 0.0, }
  }
  if k == 1 {
    return cobra3_sampling_diagnostics(chains[0])
  }
  let m = chains[0].length()
  let mut equal = true
  for c in chains {
    if c.length() != m {
      equal = false
    }
  }
  if !equal || m < 10 {
    return { r_hat: 1.0, geweke_z: 0.0, effective_sample_size: 0.0, }
  }
  let means : Array[Double] = Array(capacity=k)
  let vars : Array[Double] = Array(capacity=k)
  let mut grand = 0.0
  for c in chains {
    let mean = cobra3_mean(c, 0, m)
    means.push(mean)
    vars.push(cobra3_var(c, 0, m, mean))
    grand = grand + mean
  }
  grand = grand / k.to_double()
  let mut within = 0.0
  for i in 0..<k {
    within = within + vars[i]
  }
  within = within / k.to_double()
  let mut between = 0.0
  for i in 0..<k {
    between = between + (means[i] - grand) * (means[i] - grand)
  }
  between = between * m.to_double() / (k.to_double() - 1.0)
  let var_hat = (m.to_double() - 1.0) / m.to_double() * within +
    between / m.to_double()
  let mut r_hat = if within > 1.0e-15 {
    @math.pow(var_hat / within, 0.5)
  } else {
    1.0
  }
  if r_hat < 1.0 {
    r_hat = 1.0
  }
  let mut geweke = 0.0
  let mut ess = 0.0
  for c in chains {
    let d = cobra3_sampling_diagnostics(c)
    geweke = geweke + d.geweke_z
    ess = ess + d.effective_sample_size
  }
  { r_hat, geweke_z: geweke / k.to_double(), effective_sample_size: ess, }
}

///|
fn cobra3_mean(arr : Array[Double], start : Int32, stop : Int32) -> Double {
  let mut sum = 0.0
  let mut cnt = 0
  for i in start..<stop {
    sum = sum + arr[i]
    cnt = cnt + 1
  }
  if cnt == 0 {
    0.0
  } else {
    sum / cnt.to_double()
  }
}

///|
fn cobra3_var(
  arr : Array[Double],
  start : Int32,
  stop : Int32,
  mean : Double,
) -> Double {
  let mut sum = 0.0
  let mut cnt = 0
  for i in start..<stop {
    sum = sum + (arr[i] - mean) * (arr[i] - mean)
    cnt = cnt + 1
  }
  if cnt <= 1 {
    0.0
  } else {
    sum / (cnt - 1).to_double()
  }
}

// ---------------------------------------------------------------------------
// FastCC / loopless-lite
// ---------------------------------------------------------------------------

///|
pub(all) struct Cobra3ConsistencyResult {
  consistent_model : Cobra3Model
  removed_reactions : Array[String]
}

///|
fn cobra3_zero_objective_ranges(
  model : Cobra3Model,
) -> Array[(String, Double, Double)] {
  let n = model.reactions.length()
  let base = cobra3_steady_state_constraints(model)
  let lb = cobra3_model_lower_bounds(model)
  let ub = cobra3_model_upper_bounds(model)
  let out : Array[(String, Double, Double)] = Array(capacity=n)
  for j in 0..<n {
    let c_max : Array[Double] = Array::make(n, 0.0)
    c_max[j] = 1.0
    let c_min : Array[Double] = Array::make(n, 0.0)
    c_min[j] = -1.0
    let smax = cobra3_solve_lp(base, c_max, lb, ub)
    let smin = cobra3_solve_lp(base, c_min, lb, ub)
    let mx = if smax.status is Optimal { smax.objective_value } else { 0.0 }
    let mn = if smin.status is Optimal { -smin.objective_value } else { 0.0 }
    out.push((model.reactions[j].id, mn, mx))
  }
  out
}

///|
/// FastCC-lite：迭代移除无法承载通量的反应，得到一致子模型.
pub fn cobra3_fastcc(model : Cobra3Model) -> Cobra3ConsistencyResult {
  let mut current = model
  let removed : Array[String] = Array(capacity=0)
  let mut changed = true
  while changed {
    changed = false
    let ranges = cobra3_zero_objective_ranges(current)
    let to_remove : Array[String] = Array(capacity=0)
    for r in ranges {
      if r.1.abs() < 1.0e-9 && r.2.abs() < 1.0e-9 {
        to_remove.push(r.0)
      }
    }
    for id in to_remove {
      current = cobra3_model_remove_reaction(current, id, remove_orphans=false)
      removed.push(id)
      changed = true
    }
  }
  { consistent_model: current, removed_reactions: removed, }
}

///|
/// 检测可参与内部循环的反应（在 c=0 时通量范围跨过 0 且非零）.
pub fn cobra3_detect_loops(model : Cobra3Model) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for r in cobra3_zero_objective_ranges(model) {
    if r.2 > 1.0e-6 && r.1 < -1.0e-6 {
      out.push(r.0)
    }
  }
  out
}

///|
/// loopless-lite：先最大化目标，再最小化内部反应通量（抑制循环）.
pub fn cobra3_loopless_fba(model : Cobra3Model) -> Cobra3FluxResult {
  let opt = cobra3_model_optimize(model)
  if !(opt.status is Optimal) {
    return opt
  }
  let n = model.reactions.length()
  if n == 0 {
    return opt
  }
  let obj_vec = cobra3_model_objective_vector(model)
  let base = cobra3_steady_state_constraints(model)
  if opt.objective_value.abs() > 1.0e-9 {
    base.push({
      coeffs: obj_vec,
      rhs: cobra3_optimum_floor(opt.objective_value, 1.0 - 1.0e-6),
      sense: Cobra3Sense::GreaterEqual,
    })
  }
  let nv = 2 * n
  let lb : Array[Double] = Array::make(nv, 0.0)
  let ub : Array[Double] = Array::make(nv, 0.0)
  let obj : Array[Double] = Array::make(nv, 0.0)
  let integer : Array[Int32] = Array(capacity=n)
  for j in 0..<n {
    lb[j] = model.reactions[j].lower_bound
    ub[j] = model.reactions[j].upper_bound
    // exchange 不参与循环；内部反应用二元指示变量控制是否活跃
    if cobra3_model_is_exchange(model, model.reactions[j]) {
      lb[n + j] = 0.0
      ub[n + j] = 0.0
    } else {
      lb[n + j] = 0.0
      ub[n + j] = 1.0
      obj[n + j] = -1.0
      integer.push(n + j)
    }
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
    if cobra3_model_is_exchange(model, model.reactions[j]) {
      continue
    }
    let scale = if model.reactions[j].lower_bound.abs() >
      model.reactions[j].upper_bound.abs() {
      model.reactions[j].lower_bound.abs()
    } else {
      model.reactions[j].upper_bound.abs()
    }
    let m = 2.0 * scale + 1.0
    // z=0 -> v=0；z=1 -> |v| <= M
    let c1 : Array[Double] = Array::make(nv, 0.0)
    c1[j] = 1.0
    c1[n + j] = -m
    cons.push({ coeffs: c1, rhs: 0.0, sense: Cobra3Sense::LessEqual, })
    let c2 : Array[Double] = Array::make(nv, 0.0)
    c2[j] = -1.0
    c2[n + j] = -m
    cons.push({ coeffs: c2, rhs: 0.0, sense: Cobra3Sense::LessEqual, })
  }
  let sol = cobra3_solve_milp(cons, obj, lb, ub, integer)
  let fluxes : Array[Double] = Array(capacity=n)
  let mut growth = 0.0
  for j in 0..<n {
    fluxes.push(sol.x[j])
    growth = growth + obj_vec[j] * sol.x[j]
  }
  {
    status: sol.status,
    objective_value: growth,
    fluxes,
    shadow_prices: Array::make(model.metabolites.length(), 0.0),
    reduced_costs: Array::make(n, 0.0),
  }
}
