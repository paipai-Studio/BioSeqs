///|
/// cobra3_core.mbt — COBRApy 对象模型 v2（Metabolite / Gene / GPR / Reaction / Model）.
///
/// 与 fba_lite 的三数组平铺模型不同，这里按 COBRApy 的语义组织：
/// - Metabolite 带 formula / charge / compartment / annotation
/// - Gene 带 annotation / functional
/// - GPR 为结构化表达式（And/Or/Gene），敲除按布尔求值而非"含该基因即置零"
/// - Reaction 用 Map[metabolite_id, coefficient] 表示计量，带 subsystem / annotation
/// - Model 提供 by-id 索引、增删、快照/撤销、medium、质量平衡
///
/// 纯数据结构，不依赖求解器；求解入口见 cobra3_solver / cobra3_flux.

///|
pub(all) struct Cobra3Metabolite {
  id : String
  name : String
  formula : String
  charge : Double
  compartment : String
  boundary : Bool
  annotation : Map[String, String]
}

///|
pub(all) struct Cobra3Gene {
  id : String
  name : String
  annotation : Map[String, String]
  functional : Bool
}

///|
/// 结构化基因-蛋白-反应关联表达式.
pub enum Cobra3Gpr {
  Gene(String)
  And(Array[Cobra3Gpr])
  Or(Array[Cobra3Gpr])
}

///|
pub(all) struct Cobra3Reaction {
  id : String
  name : String
  metabolites : Map[String, Double]
  lower_bound : Double
  upper_bound : Double
  objective : Double
  subsystem : String
  gpr : Cobra3Gpr?
  annotation : Map[String, String]
}

///|
/// 模型快照，用于 undo / context 语义.
pub(all) struct Cobra3Snapshot {
  metabolites : Array[Cobra3Metabolite]
  genes : Array[Cobra3Gene]
  reactions : Array[Cobra3Reaction]
  compartments : Map[String, String]
}

///|
pub(all) struct Cobra3Model {
  id : String
  name : String
  metabolites : Array[Cobra3Metabolite]
  genes : Array[Cobra3Gene]
  reactions : Array[Cobra3Reaction]
  compartments : Map[String, String]
  annotation : Map[String, String]
  undo_stack : Array[Cobra3Snapshot]
}

///|
pub fn cobra3_metabolite_new(id : String) -> Cobra3Metabolite {
  {
    id,
    name: id,
    formula: "",
    charge: 0.0,
    compartment: "c",
    boundary: false,
    annotation: Map([], capacity=2),
  }
}

///|
pub fn cobra3_gene_new(id : String) -> Cobra3Gene {
  { id, name: id, annotation: Map([], capacity=2), functional: true, }
}

///|
pub fn cobra3_reaction_new(id : String) -> Cobra3Reaction {
  {
    id,
    name: id,
    metabolites: Map([], capacity=4),
    lower_bound: 0.0,
    upper_bound: 1000.0,
    objective: 0.0,
    subsystem: "",
    gpr: None,
    annotation: Map([], capacity=2),
  }
}

///|
pub fn cobra3_model_new(id : String) -> Cobra3Model {
  let compartments : Map[String, String] = Map([], capacity=2)
  compartments["c"] = "cytosol"
  {
    id,
    name: id,
    metabolites: Array(capacity=8),
    genes: Array(capacity=8),
    reactions: Array(capacity=8),
    compartments,
    annotation: Map([], capacity=2),
    undo_stack: Array(capacity=2),
  }
}

// ---------------------------------------------------------------------------
// 索引与查找
// ---------------------------------------------------------------------------

///|
pub fn cobra3_model_metabolite_index(model : Cobra3Model, id : String) -> Int32 {
  for i in 0..<model.metabolites.length() {
    if model.metabolites[i].id == id {
      return i
    }
  }
  -1
}

///|
pub fn cobra3_model_reaction_index(model : Cobra3Model, id : String) -> Int32 {
  for i in 0..<model.reactions.length() {
    if model.reactions[i].id == id {
      return i
    }
  }
  -1
}

///|
pub fn cobra3_model_gene_index(model : Cobra3Model, id : String) -> Int32 {
  for i in 0..<model.genes.length() {
    if model.genes[i].id == id {
      return i
    }
  }
  -1
}

///|
pub fn cobra3_model_get_metabolite(
  model : Cobra3Model,
  id : String,
) -> Cobra3Metabolite? {
  let idx = cobra3_model_metabolite_index(model, id)
  if idx < 0 {
    None
  } else {
    Some(model.metabolites[idx])
  }
}

///|
pub fn cobra3_model_get_reaction(
  model : Cobra3Model,
  id : String,
) -> Cobra3Reaction? {
  let idx = cobra3_model_reaction_index(model, id)
  if idx < 0 {
    None
  } else {
    Some(model.reactions[idx])
  }
}

///|
pub fn cobra3_model_get_gene(model : Cobra3Model, id : String) -> Cobra3Gene? {
  let idx = cobra3_model_gene_index(model, id)
  if idx < 0 {
    None
  } else {
    Some(model.genes[idx])
  }
}

///|
pub fn cobra3_model_has_reaction(model : Cobra3Model, id : String) -> Bool {
  cobra3_model_reaction_index(model, id) >= 0
}

// ---------------------------------------------------------------------------
// 模型编辑（函数式：所有修改返回新模型，并压入 undo 栈）
// ---------------------------------------------------------------------------

///|
pub fn cobra3_model_snapshot(model : Cobra3Model) -> Cobra3Snapshot {
  {
    metabolites: model.metabolites,
    genes: model.genes,
    reactions: model.reactions,
    compartments: model.compartments,
  }
}

///|
fn cobra3_model_with_history(
  model : Cobra3Model,
  snap : Cobra3Snapshot,
) -> Cobra3Model {
  let stack : Array[Cobra3Snapshot] = Array(capacity=0)
  // 最多保留最近 99 个快照 + 本次，避免反复编辑时历史无限增长
  let start = if model.undo_stack.length() >= 100 {
    model.undo_stack.length() - 99
  } else {
    0
  }
  for i in start..<model.undo_stack.length() {
    stack.push(model.undo_stack[i])
  }
  stack.push(snap)
  {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions: model.reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: stack,
  }
}

///|
/// 撤销最近一次编辑；无历史时原样返回.
pub fn cobra3_model_undo(model : Cobra3Model) -> Cobra3Model {
  let n = model.undo_stack.length()
  if n == 0 {
    return model
  }
  let prev = model.undo_stack[n - 1]
  let stack : Array[Cobra3Snapshot] = Array(capacity=0)
  for i in 0..<(n - 1) {
    stack.push(model.undo_stack[i])
  }
  {
    id: model.id,
    name: model.name,
    metabolites: prev.metabolites,
    genes: prev.genes,
    reactions: prev.reactions,
    compartments: prev.compartments,
    annotation: model.annotation,
    undo_stack: stack,
  }
}

///|
pub fn cobra3_model_has_history(model : Cobra3Model) -> Bool {
  model.undo_stack.length() > 0
}

///|
pub fn cobra3_model_add_metabolite(
  model : Cobra3Model,
  met : Cobra3Metabolite,
) -> Cobra3Model {
  if cobra3_model_metabolite_index(model, met.id) >= 0 {
    return model
  }
  let snapshot = cobra3_model_snapshot(model)
  let mets : Array[Cobra3Metabolite] = Array(capacity=0)
  for m in model.metabolites {
    mets.push(m)
  }
  mets.push(met)
  let next = {
    id: model.id,
    name: model.name,
    metabolites: mets,
    genes: model.genes,
    reactions: model.reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  cobra3_model_with_history(next, snapshot)
}

///|
pub fn cobra3_model_add_gene(
  model : Cobra3Model,
  gene : Cobra3Gene,
) -> Cobra3Model {
  if cobra3_model_gene_index(model, gene.id) >= 0 {
    return model
  }
  let snapshot = cobra3_model_snapshot(model)
  let genes : Array[Cobra3Gene] = Array(capacity=0)
  for g in model.genes {
    genes.push(g)
  }
  genes.push(gene)
  let next = {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes,
    reactions: model.reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  cobra3_model_with_history(next, snapshot)
}

///|
/// 添加反应；计量中的代谢物若不存在则自动创建，GPR 中基因若不存在则自动创建.
/// 自动创建代谢物/基因与添加反应属于同一次编辑，undo 一次即可整体回退.
pub fn cobra3_model_add_reaction(
  model : Cobra3Model,
  reaction : Cobra3Reaction,
) -> Cobra3Model {
  if cobra3_model_reaction_index(model, reaction.id) >= 0 {
    return model
  }
  let snapshot = cobra3_model_snapshot(model)
  let mets : Array[Cobra3Metabolite] = Array(
    capacity=model.metabolites.length() + reaction.metabolites.length(),
  )
  for met in model.metabolites {
    mets.push(met)
  }
  for met_id in reaction.metabolites.keys() {
    if cobra3_model_metabolite_index(model, met_id) < 0 {
      let met = cobra3_metabolite_new(met_id)
      let comp = cobra3_metabolite_compartment_from_id(met_id)
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
  let genes : Array[Cobra3Gene] = Array(capacity=model.genes.length() + 4)
  for gene in model.genes {
    genes.push(gene)
  }
  match reaction.gpr {
    Some(gpr) =>
      for gene_id in cobra3_gpr_genes(gpr) {
        if cobra3_model_gene_index(model, gene_id) < 0 {
          genes.push(cobra3_gene_new(gene_id))
        }
      }
    None => ()
  }
  let reactions : Array[Cobra3Reaction] = Array(
    capacity=model.reactions.length() + 1,
  )
  for r in model.reactions {
    reactions.push(r)
  }
  reactions.push(reaction)
  let next = {
    id: model.id,
    name: model.name,
    metabolites: mets,
    genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  cobra3_model_with_history(next, snapshot)
}

///|
pub fn cobra3_model_add_reaction_simple(
  model : Cobra3Model,
  id : String,
  participants : Array[(String, Double)],
  lower_bound : Double,
  upper_bound : Double,
) -> Cobra3Model {
  let met_map : Map[String, Double] = Map([], capacity=participants.length())
  for p in participants {
    let prev = match met_map.get(p.0) {
      Some(v) => v
      None => 0.0
    }
    met_map[p.0] = prev + p.1
  }
  let r : Cobra3Reaction = {
    id,
    name: id,
    metabolites: met_map,
    lower_bound,
    upper_bound,
    objective: 0.0,
    subsystem: "",
    gpr: None,
    annotation: Map([], capacity=2),
  }
  cobra3_model_add_reaction(model, r)
}

///|
/// 删除反应；remove_orphans 为 true 时清理不再参与任何反应的代谢物.
pub fn cobra3_model_remove_reaction(
  model : Cobra3Model,
  id : String,
  remove_orphans~ : Bool = true,
) -> Cobra3Model {
  let idx = cobra3_model_reaction_index(model, id)
  if idx < 0 {
    return model
  }
  let snapshot = cobra3_model_snapshot(model)
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for i in 0..<model.reactions.length() {
    if i != idx {
      reactions.push(model.reactions[i])
    }
  }
  let mut next = {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  if remove_orphans {
    next = cobra3_model_prune_orphan_metabolites(next)
  }
  cobra3_model_with_history(next, snapshot)
}

///|
pub fn cobra3_model_prune_orphan_metabolites(
  model : Cobra3Model,
) -> Cobra3Model {
  let used : Map[String, Bool] = Map([], capacity=model.metabolites.length())
  for r in model.reactions {
    for met_id in r.metabolites.keys() {
      used[met_id] = true
    }
  }
  let mets : Array[Cobra3Metabolite] = Array(capacity=0)
  for met in model.metabolites {
    if used.contains(met.id) {
      mets.push(met)
    }
  }
  {
    id: model.id,
    name: model.name,
    metabolites: mets,
    genes: model.genes,
    reactions: model.reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
}

///|
pub fn cobra3_model_prune_orphan_genes(model : Cobra3Model) -> Cobra3Model {
  let used : Map[String, Bool] = Map([], capacity=model.genes.length())
  for r in model.reactions {
    match r.gpr {
      Some(gpr) =>
        for g in cobra3_gpr_genes(gpr) {
          used[g] = true
        }
      None => ()
    }
  }
  let genes : Array[Cobra3Gene] = Array(capacity=0)
  for g in model.genes {
    if used.contains(g.id) {
      genes.push(g)
    }
  }
  {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes,
    reactions: model.reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
}

///|
pub fn cobra3_model_remove_metabolite(
  model : Cobra3Model,
  id : String,
) -> Cobra3Model {
  let idx = cobra3_model_metabolite_index(model, id)
  if idx < 0 {
    return model
  }
  let snapshot = cobra3_model_snapshot(model)
  let mets : Array[Cobra3Metabolite] = Array(capacity=0)
  for i in 0..<model.metabolites.length() {
    if i != idx {
      mets.push(model.metabolites[i])
    }
  }
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for r in model.reactions {
    let participants : Map[String, Double] = Map(
      [],
      capacity=r.metabolites.length(),
    )
    for k in r.metabolites.keys() {
      if k != id {
        participants[k] = r.metabolites[k]
      }
    }
    reactions.push({
      id: r.id,
      name: r.name,
      metabolites: participants,
      lower_bound: r.lower_bound,
      upper_bound: r.upper_bound,
      objective: r.objective,
      subsystem: r.subsystem,
      gpr: r.gpr,
      annotation: r.annotation,
    })
  }
  let next = {
    id: model.id,
    name: model.name,
    metabolites: mets,
    genes: model.genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  cobra3_model_with_history(next, snapshot)
}

///|
/// 修改反应边界（函数式）.
pub fn cobra3_model_set_reaction_bounds(
  model : Cobra3Model,
  id : String,
  lower_bound : Double,
  upper_bound : Double,
) -> Cobra3Model {
  let idx = cobra3_model_reaction_index(model, id)
  if idx < 0 {
    return model
  }
  let snapshot = cobra3_model_snapshot(model)
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    if i == idx {
      reactions.push({
        id: r.id,
        name: r.name,
        metabolites: r.metabolites,
        lower_bound,
        upper_bound,
        objective: r.objective,
        subsystem: r.subsystem,
        gpr: r.gpr,
        annotation: r.annotation,
      })
    } else {
      reactions.push(r)
    }
  }
  let next = {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  cobra3_model_with_history(next, snapshot)
}

///|
/// 设置单个目标系数；coeff 为 0 时等价于取消该反应的目标项.
pub fn cobra3_model_set_objective(
  model : Cobra3Model,
  id : String,
  coeff : Double,
) -> Cobra3Model {
  let idx = cobra3_model_reaction_index(model, id)
  if idx < 0 {
    return model
  }
  let snapshot = cobra3_model_snapshot(model)
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for i in 0..<model.reactions.length() {
    let r = model.reactions[i]
    if i == idx {
      reactions.push({
        id: r.id,
        name: r.name,
        metabolites: r.metabolites,
        lower_bound: r.lower_bound,
        upper_bound: r.upper_bound,
        objective: coeff,
        subsystem: r.subsystem,
        gpr: r.gpr,
        annotation: r.annotation,
      })
    } else {
      reactions.push(r)
    }
  }
  let next = {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  cobra3_model_with_history(next, snapshot)
}

///|
pub fn cobra3_model_clear_objective(model : Cobra3Model) -> Cobra3Model {
  let snapshot = cobra3_model_snapshot(model)
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for r in model.reactions {
    reactions.push({
      id: r.id,
      name: r.name,
      metabolites: r.metabolites,
      lower_bound: r.lower_bound,
      upper_bound: r.upper_bound,
      objective: 0.0,
      subsystem: r.subsystem,
      gpr: r.gpr,
      annotation: r.annotation,
    })
  }
  let next = {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  cobra3_model_with_history(next, snapshot)
}

///|
/// 主目标反应（objective 系数最大者，同分取先出现者）.
pub fn cobra3_model_objective_reaction(model : Cobra3Model) -> String? {
  let mut best = -1
  let mut best_val = 0.0
  for i in 0..<model.reactions.length() {
    if model.reactions[i].objective.abs() > best_val {
      best_val = model.reactions[i].objective.abs()
      best = i
    }
  }
  if best < 0 {
    None
  } else {
    Some(model.reactions[best].id)
  }
}

///|
/// 主 biomass 反应：优先 id 含 biomass 者，否则回退到目标系数最大的反应.
pub fn cobra3_model_biomass_reaction(model : Cobra3Model) -> String? {
  let mut found : String? = None
  for r in model.reactions {
    if found is None && r.id.to_lower().contains("biomass") {
      found = Some(r.id)
    }
  }
  match found {
    Some(id) => Some(id)
    None => cobra3_model_objective_reaction(model)
  }
}

// ---------------------------------------------------------------------------
// medium / exchange
// ---------------------------------------------------------------------------

///|
pub fn cobra3_model_is_exchange(
  model : Cobra3Model,
  reaction : Cobra3Reaction,
) -> Bool {
  if reaction.id.has_prefix("EX_") || reaction.id.has_prefix("Ex_") {
    return true
  }
  let mut is_boundary = false
  if reaction.metabolites.length() == 1 {
    for k in reaction.metabolites.keys() {
      if !is_boundary {
        match cobra3_model_get_metabolite(model, k) {
          Some(met) => if met.boundary { is_boundary = true }
          None => ()
        }
      }
    }
  }
  is_boundary
}

///|
pub fn cobra3_model_exchange_reactions(model : Cobra3Model) -> Array[String] {
  let ids : Array[String] = Array(capacity=0)
  for r in model.reactions {
    if cobra3_model_is_exchange(model, r) {
      ids.push(r.id)
    }
  }
  ids
}

///|
/// 当前培养基：exchange id -> 摄取上限（-lower_bound）.
pub fn cobra3_model_medium(model : Cobra3Model) -> Map[String, Double] {
  let medium : Map[String, Double] = Map([], capacity=8)
  for r in model.reactions {
    if cobra3_model_is_exchange(model, r) && r.lower_bound < 0.0 {
      medium[r.id] = -r.lower_bound
    }
  }
  medium
}

///|
/// 设置培养基：未列出的 exchange 关闭摄取；列出的按摄取上限打开.
pub fn cobra3_model_set_medium(
  model : Cobra3Model,
  medium : Map[String, Double],
) -> Cobra3Model {
  let snapshot = cobra3_model_snapshot(model)
  let reactions : Array[Cobra3Reaction] = Array(capacity=0)
  for r in model.reactions {
    if cobra3_model_is_exchange(model, r) {
      let lb = match medium.get(r.id) {
        Some(v) => -v
        None => 0.0
      }
      reactions.push({
        id: r.id,
        name: r.name,
        metabolites: r.metabolites,
        lower_bound: lb,
        upper_bound: r.upper_bound,
        objective: r.objective,
        subsystem: r.subsystem,
        gpr: r.gpr,
        annotation: r.annotation,
      })
    } else {
      reactions.push(r)
    }
  }
  let next = {
    id: model.id,
    name: model.name,
    metabolites: model.metabolites,
    genes: model.genes,
    reactions,
    compartments: model.compartments,
    annotation: model.annotation,
    undo_stack: model.undo_stack,
  }
  cobra3_model_with_history(next, snapshot)
}

///|
pub fn cobra3_model_set_medium_component(
  model : Cobra3Model,
  exchange_id : String,
  uptake : Double,
) -> Cobra3Model {
  let medium = cobra3_model_medium(model)
  medium[exchange_id] = uptake
  cobra3_model_set_medium(model, medium)
}

// ---------------------------------------------------------------------------
// 化学计量 / 转换
// ---------------------------------------------------------------------------

///|
/// 反应参与者按 id 排序输出，保证确定性.
pub fn cobra3_reaction_participants(
  reaction : Cobra3Reaction,
) -> Array[(String, Double)] {
  let ids : Array[String] = Array(capacity=0)
  for k in reaction.metabolites.keys() {
    ids.push(k)
  }
  cobra3_sort_strings(ids)
  let out : Array[(String, Double)] = Array(capacity=ids.length())
  for id in ids {
    match reaction.metabolites.get(id) {
      Some(v) => out.push((id, v))
      None => ()
    }
  }
  out
}

///|
/// 构造化学计量矩阵 S（metabolites x reactions）.
pub fn cobra3_model_stoichiometric_matrix(
  model : Cobra3Model,
) -> Array[Array[Double]] {
  let m = model.metabolites.length()
  let n = model.reactions.length()
  let s : Array[Array[Double]] = Array(capacity=m)
  for _ in 0..<m {
    s.push(Array::make(n, 0.0))
  }
  let index : Map[String, Int32] = Map([], capacity=m)
  for i in 0..<m {
    index[model.metabolites[i].id] = i
  }
  for j in 0..<n {
    let r = model.reactions[j]
    for met_id in r.metabolites.keys() {
      match index.get(met_id) {
        Some(i) => s[i][j] = r.metabolites[met_id]
        None => ()
      }
    }
  }
  s
}

///|
pub fn cobra3_model_objective_vector(model : Cobra3Model) -> Array[Double] {
  let c : Array[Double] = Array(capacity=model.reactions.length())
  for r in model.reactions {
    c.push(r.objective)
  }
  c
}

///|
pub fn cobra3_model_lower_bounds(model : Cobra3Model) -> Array[Double] {
  let lb : Array[Double] = Array(capacity=model.reactions.length())
  for r in model.reactions {
    lb.push(r.lower_bound)
  }
  lb
}

///|
pub fn cobra3_model_upper_bounds(model : Cobra3Model) -> Array[Double] {
  let ub : Array[Double] = Array(capacity=model.reactions.length())
  for r in model.reactions {
    ub.push(r.upper_bound)
  }
  ub
}

///|
fn cobra3_copy_string_map(m : Map[String, String]) -> Map[String, String] {
  let out : Map[String, String] = Map([], capacity=m.length())
  for k in m.keys() {
    out[k] = m[k]
  }
  out
}

///|
fn cobra3_copy_double_map(m : Map[String, Double]) -> Map[String, Double] {
  let out : Map[String, Double] = Map([], capacity=m.length())
  for k in m.keys() {
    out[k] = m[k]
  }
  out
}

///|
/// 深拷贝模型；Map 字段与 undo 历史不共享.
pub fn cobra3_model_copy(model : Cobra3Model) -> Cobra3Model {
  let mets : Array[Cobra3Metabolite] = Array(
    capacity=model.metabolites.length(),
  )
  for m in model.metabolites {
    mets.push({
      id: m.id,
      name: m.name,
      formula: m.formula,
      charge: m.charge,
      compartment: m.compartment,
      boundary: m.boundary,
      annotation: cobra3_copy_string_map(m.annotation),
    })
  }
  let genes : Array[Cobra3Gene] = Array(capacity=model.genes.length())
  for g in model.genes {
    genes.push({
      id: g.id,
      name: g.name,
      annotation: cobra3_copy_string_map(g.annotation),
      functional: g.functional,
    })
  }
  let reactions : Array[Cobra3Reaction] = Array(
    capacity=model.reactions.length(),
  )
  for r in model.reactions {
    reactions.push({
      id: r.id,
      name: r.name,
      metabolites: cobra3_copy_double_map(r.metabolites),
      lower_bound: r.lower_bound,
      upper_bound: r.upper_bound,
      objective: r.objective,
      subsystem: r.subsystem,
      gpr: r.gpr,
      annotation: cobra3_copy_string_map(r.annotation),
    })
  }
  {
    id: model.id,
    name: model.name,
    metabolites: mets,
    genes,
    reactions,
    compartments: cobra3_copy_string_map(model.compartments),
    annotation: cobra3_copy_string_map(model.annotation),
    undo_stack: Array(capacity=0),
  }
}

///|
/// 由 fba_lite 平铺模型转换（保留 S 矩阵符号约定）.
/// 一次性构造，避免逐条 add_reaction 带来的数组复制与 undo 栈开销.
pub fn cobra3_from_fba_model(model : FbaModel) -> Cobra3Model {
  let mets : Array[Cobra3Metabolite] = Array(
    capacity=model.metabolites.length(),
  )
  for met_id in model.metabolites {
    let met = cobra3_metabolite_new(met_id)
    let comp = cobra3_metabolite_compartment_from_id(met_id)
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
  let genes : Array[Cobra3Gene] = Array(capacity=8)
  let seen_genes : Map[String, Bool] = Map([], capacity=8)
  for r in model.reactions {
    for g in r.genes {
      if !seen_genes.contains(g) {
        seen_genes[g] = true
        genes.push(cobra3_gene_new(g))
      }
    }
  }
  let reactions : Array[Cobra3Reaction] = Array(
    capacity=model.reactions.length(),
  )
  for j in 0..<model.reactions.length() {
    let r = model.reactions[j]
    let participants : Map[String, Double] = Map([], capacity=4)
    for i in 0..<model.metabolites.length() {
      let v = model.stoichiometry[i][j]
      if v != 0.0 {
        participants[model.metabolites[i]] = v
      }
    }
    let nodes : Array[Cobra3Gpr] = Array(capacity=r.genes.length())
    for g in r.genes {
      nodes.push(Cobra3Gpr::Gene(g))
    }
    reactions.push({
      id: r.id,
      name: r.name,
      metabolites: participants,
      lower_bound: r.lower_bound,
      upper_bound: r.upper_bound,
      objective: r.objective,
      subsystem: "",
      gpr: if nodes.length() == 0 {
        None
      } else {
        Some(Cobra3Gpr::Or(nodes))
      },
      annotation: Map([], capacity=2),
    })
  }
  let compartments : Map[String, String] = Map([], capacity=2)
  compartments["c"] = "cytosol"
  {
    id: "model",
    name: "model",
    metabolites: mets,
    genes,
    reactions,
    compartments,
    annotation: Map([], capacity=2),
    undo_stack: Array(capacity=0),
  }
}

///|
/// 转回 fba_lite 平铺模型（供旧算法复用）.
/// FbaModel 只保存基因列表，GPR 的 AND/嵌套结构会退化为单层 OR；需要完整 GPR 时留在 Cobra3Model.
pub fn cobra3_to_fba_model(model : Cobra3Model) -> FbaModel {
  let met_ids : Array[String] = Array(capacity=0)
  for m in model.metabolites {
    met_ids.push(m.id)
  }
  let s = cobra3_model_stoichiometric_matrix(model)
  let reactions : Array[FbaReaction] = Array(capacity=0)
  for r in model.reactions {
    let genes : Array[String] = Array(capacity=0)
    match r.gpr {
      Some(gpr) =>
        for g in cobra3_gpr_genes(gpr) {
          genes.push(g)
        }
      None => ()
    }
    reactions.push({
      id: r.id,
      name: r.name,
      lower_bound: r.lower_bound,
      upper_bound: r.upper_bound,
      objective: r.objective,
      genes,
    })
  }
  { metabolites: met_ids, reactions, stoichiometry: s, }
}

///|
/// 与 cobra3_to_fba_model 相同，但额外返回反应的 GPR 字符串（key=reaction id），
/// 供需要保留 AND/嵌套结构的外部调用方使用.
pub fn cobra3_to_fba_model_with_gpr(
  model : Cobra3Model,
) -> (FbaModel, Map[String, String]) {
  let gprs : Map[String, String] = Map([], capacity=model.reactions.length())
  for r in model.reactions {
    match r.gpr {
      Some(g) => gprs[r.id] = cobra3_gpr_to_string(g)
      None => ()
    }
  }
  (cobra3_to_fba_model(model), gprs)
}

// ---------------------------------------------------------------------------
// GPR：解析、求值、基因集、化简、输出
// ---------------------------------------------------------------------------

///|
fn cobra3_is_gene_char(c : Int32) -> Bool {
  (c >= 48 && c <= 57) ||
  (c >= 65 && c <= 90) ||
  (c >= 97 && c <= 122) ||
  c == 95 ||
  c == 45 ||
  c == 46 ||
  c == 58
}

///|
fn cobra3_gpr_tokenize(text : String) -> Array[String]? {
  let tokens : Array[String] = Array(capacity=8)
  let mut i = 0
  while i < text.length() {
    let c = text[i].to_int()
    if c == 32 || c == 9 || c == 10 || c == 13 {
      i = i + 1
      continue
    }
    if c == 40 {
      tokens.push("(")
      i = i + 1
      continue
    }
    if c == 41 {
      tokens.push(")")
      i = i + 1
      continue
    }
    let start = i
    while i < text.length() && cobra3_is_gene_char(text[i].to_int()) {
      i = i + 1
    }
    if i == start {
      return None
    }
    tokens.push(text[start:i].to_owned())
  }
  Some(tokens)
}

///|
fn cobra3_gpr_parse_factor(
  tokens : Array[String],
  pos : Int32,
) -> (Cobra3Gpr?, Int32) {
  if pos >= tokens.length() {
    return (None, -1)
  }
  let t = tokens[pos]
  if t == "(" {
    let (node, next) = cobra3_gpr_parse_or(tokens, pos + 1)
    match node {
      None => (None, -1)
      Some(n) =>
        if next < 0 || next >= tokens.length() || tokens[next] != ")" {
          (None, -1)
        } else {
          (Some(n), next + 1)
        }
    }
  } else if t == ")" || t == "and" || t == "AND" || t == "or" || t == "OR" {
    (None, -1)
  } else {
    (Some(Cobra3Gpr::Gene(t)), pos + 1)
  }
}

///|
fn cobra3_gpr_parse_and(
  tokens : Array[String],
  pos : Int32,
) -> (Cobra3Gpr?, Int32) {
  let (first, next) = cobra3_gpr_parse_factor(tokens, pos)
  match first {
    None => (None, -1)
    Some(head) => {
      let nodes : Array[Cobra3Gpr] = Array(capacity=2)
      nodes.push(head)
      let mut p = next
      while p >= 0 &&
            p < tokens.length() &&
            (tokens[p] == "and" || tokens[p] == "AND") {
        let (child, next2) = cobra3_gpr_parse_factor(tokens, p + 1)
        match child {
          None => return (None, -1)
          Some(c) => {
            nodes.push(c)
            p = next2
          }
        }
      }
      if nodes.length() == 1 {
        (Some(nodes[0]), p)
      } else {
        (Some(Cobra3Gpr::And(nodes)), p)
      }
    }
  }
}

///|
fn cobra3_gpr_parse_or(
  tokens : Array[String],
  pos : Int32,
) -> (Cobra3Gpr?, Int32) {
  let (first, next) = cobra3_gpr_parse_and(tokens, pos)
  match first {
    None => (None, -1)
    Some(head) => {
      let nodes : Array[Cobra3Gpr] = Array(capacity=2)
      nodes.push(head)
      let mut p = next
      while p >= 0 &&
            p < tokens.length() &&
            (tokens[p] == "or" || tokens[p] == "OR") {
        let (child, next2) = cobra3_gpr_parse_and(tokens, p + 1)
        match child {
          None => return (None, -1)
          Some(c) => {
            nodes.push(c)
            p = next2
          }
        }
      }
      if nodes.length() == 1 {
        (Some(nodes[0]), p)
      } else {
        (Some(Cobra3Gpr::Or(nodes)), p)
      }
    }
  }
}

///|
/// 解析 COBRA 风格 GPR 字符串，如 "(A and B) or C".
pub fn cobra3_gpr_parse(text : String) -> Cobra3Gpr? {
  // 空表达式常量：与 cobra3_gpr_to_string 的输出保持一致
  if text == "true" {
    return Some(Cobra3Gpr::And([]))
  }
  if text == "false" {
    return Some(Cobra3Gpr::Or([]))
  }
  match cobra3_gpr_tokenize(text) {
    None => None
    Some(tokens) => {
      if tokens.length() == 0 {
        return None
      }
      let (node, next) = cobra3_gpr_parse_or(tokens, 0)
      match node {
        None => None
        Some(n) => if next == tokens.length() { Some(n) } else { None }
      }
    }
  }
}

///|
/// 按基因活性求值；未给出的基因视为不可用.
pub fn cobra3_gpr_eval(gpr : Cobra3Gpr, active : Map[String, Bool]) -> Bool {
  match gpr {
    Gene(id) =>
      match active.get(id) {
        Some(v) => v
        None => false
      }
    And(children) => {
      let mut ok = true
      for c in children {
        if !cobra3_gpr_eval(c, active) {
          ok = false
          break
        }
      }
      ok
    }
    Or(children) => {
      let mut ok = false
      for c in children {
        if cobra3_gpr_eval(c, active) {
          ok = true
          break
        }
      }
      ok
    }
  }
}

///|
pub fn cobra3_gpr_genes(gpr : Cobra3Gpr) -> Array[String] {
  let out : Array[String] = Array(capacity=4)
  cobra3_gpr_collect_genes(gpr, out)
  let dedup : Array[String] = Array(capacity=out.length())
  for g in out {
    if !dedup.contains(g) {
      dedup.push(g)
    }
  }
  cobra3_sort_strings(dedup)
  dedup
}

///|
fn cobra3_gpr_collect_genes(gpr : Cobra3Gpr, out : Array[String]) -> Unit {
  match gpr {
    Gene(id) => out.push(id)
    And(children) =>
      for c in children {
        cobra3_gpr_collect_genes(c, out)
      }
    Or(children) =>
      for c in children {
        cobra3_gpr_collect_genes(c, out)
      }
  }
}

///|
pub fn cobra3_gpr_to_string(gpr : Cobra3Gpr) -> String {
  match gpr {
    Gene(id) => id
    And(children) =>
      if children.length() == 0 {
        "true"
      } else {
        let sb = StringBuilder()
        let _ = sb.write("(")
        for i in 0..<children.length() {
          if i > 0 {
            let _ = sb.write(" and ")
          }
          let _ = sb.write(cobra3_gpr_to_string(children[i]))
        }
        let _ = sb.write(")")
        sb.to_string()
      }
    Or(children) =>
      if children.length() == 0 {
        "false"
      } else {
        let sb = StringBuilder()
        let _ = sb.write("(")
        for i in 0..<children.length() {
          if i > 0 {
            let _ = sb.write(" or ")
          }
          let _ = sb.write(cobra3_gpr_to_string(children[i]))
        }
        let _ = sb.write(")")
        sb.to_string()
      }
  }
}

///|
/// 基本化简：去重 + 单元素折叠 + 嵌套同型展平.
pub fn cobra3_gpr_simplify(gpr : Cobra3Gpr) -> Cobra3Gpr {
  match gpr {
    Gene(id) => Cobra3Gpr::Gene(id)
    And(children) => {
      let flat : Array[Cobra3Gpr] = Array(capacity=0)
      for c in children {
        match cobra3_gpr_simplify(c) {
          And(inner) =>
            for x in inner {
              flat.push(x)
            }
          other => flat.push(other)
        }
      }
      let uniq : Array[Cobra3Gpr] = Array(capacity=flat.length())
      let seen : Array[String] = Array(capacity=flat.length())
      for x in flat {
        let key = cobra3_gpr_to_string(x)
        if !seen.contains(key) {
          seen.push(key)
          uniq.push(x)
        }
      }
      if uniq.length() == 1 {
        uniq[0]
      } else {
        Cobra3Gpr::And(uniq)
      }
    }
    Or(children) => {
      let flat : Array[Cobra3Gpr] = Array(capacity=0)
      for c in children {
        match cobra3_gpr_simplify(c) {
          Or(inner) =>
            for x in inner {
              flat.push(x)
            }
          other => flat.push(other)
        }
      }
      let uniq : Array[Cobra3Gpr] = Array(capacity=flat.length())
      let seen : Array[String] = Array(capacity=flat.length())
      for x in flat {
        let key = cobra3_gpr_to_string(x)
        if !seen.contains(key) {
          seen.push(key)
          uniq.push(x)
        }
      }
      if uniq.length() == 1 {
        uniq[0]
      } else {
        Cobra3Gpr::Or(uniq)
      }
    }
  }
}

///|
/// 基因敲除后的 GPR：把指定基因替换为不可满足的 Or([]).
pub fn cobra3_gpr_knockout(gpr : Cobra3Gpr, gene_id : String) -> Cobra3Gpr {
  match gpr {
    Gene(id) =>
      if id == gene_id {
        Cobra3Gpr::Or([])
      } else {
        Cobra3Gpr::Gene(id)
      }
    And(children) => {
      let out : Array[Cobra3Gpr] = Array(capacity=0)
      for c in children {
        let k = cobra3_gpr_knockout(c, gene_id)
        if cobra3_gpr_is_false(k) {
          return Cobra3Gpr::Or([])
        }
        out.push(k)
      }
      Cobra3Gpr::And(out)
    }
    Or(children) => {
      let out : Array[Cobra3Gpr] = Array(capacity=0)
      for c in children {
        let k = cobra3_gpr_knockout(c, gene_id)
        if !cobra3_gpr_is_false(k) {
          out.push(k)
        }
      }
      if out.length() == 0 {
        Cobra3Gpr::Or([])
      } else if out.length() == 1 {
        out[0]
      } else {
        Cobra3Gpr::Or(out)
      }
    }
  }
}

///|
fn cobra3_gpr_is_false(gpr : Cobra3Gpr) -> Bool {
  match gpr {
    Or(children) => children.length() == 0
    _ => false
  }
}

///|
/// 给定敲除基因集合，反应是否仍可承载通量.
pub fn cobra3_reaction_active_after_knockout(
  reaction : Cobra3Reaction,
  knocked_out : Array[String],
) -> Bool {
  match reaction.gpr {
    None => true
    Some(gpr) => {
      let active : Map[String, Bool] = Map(
        [],
        capacity=cobra3_gpr_genes(gpr).length(),
      )
      for g in cobra3_gpr_genes(gpr) {
        active[g] = !knocked_out.contains(g)
      }
      cobra3_gpr_eval(gpr, active)
    }
  }
}

///|
/// 结合 Gene.functional 判断敲除后反应是否仍可承载通量.
pub fn cobra3_reaction_active_in_model(
  model : Cobra3Model,
  reaction : Cobra3Reaction,
  knocked_out : Array[String],
) -> Bool {
  match reaction.gpr {
    None => true
    Some(gpr) => {
      let active : Map[String, Bool] = Map(
        [],
        capacity=cobra3_gpr_genes(gpr).length(),
      )
      for g in cobra3_gpr_genes(gpr) {
        let functional = match cobra3_model_get_gene(model, g) {
          Some(gene) => gene.functional
          None => true
        }
        active[g] = functional && !knocked_out.contains(g)
      }
      cobra3_gpr_eval(gpr, active)
    }
  }
}

// ---------------------------------------------------------------------------
// 质量平衡 / formula
// ---------------------------------------------------------------------------

///|
fn cobra3_is_element_start(c : Int32) -> Bool {
  c >= 65 && c <= 90
}

///|
fn cobra3_is_element_char(c : Int32) -> Bool {
  // 元素符号首字母后只允许小写字母（Ca、Cl、Na），避免把 "OH" 当成元素
  c >= 97 && c <= 122
}

///|
/// 解析分子式元素计数，如 "C6H12O6"、"Ca(OH)2"；不支持电荷/同位素.

///|
/// 递归解析一组元素/括号，返回计数与下一位置.
fn cobra3_parse_formula_rec(
  formula : String,
  start : Int32,
  stop_at_paren : Bool,
) -> (Map[String, Double], Int32) {
  let out : Map[String, Double] = Map([], capacity=8)
  let mut i = start
  while i < formula.length() {
    let c = formula[i].to_int()
    if c == 41 {
      if stop_at_paren {
        return (out, i + 1)
      }
      i = i + 1
      continue
    }
    if c == 40 {
      let (inner, next) = cobra3_parse_formula_rec(formula, i + 1, true)
      i = next
      let mut mult = 1.0
      let num_start = i
      while i < formula.length() &&
            formula[i].to_int() >= 48 &&
            formula[i].to_int() <= 57 {
        i = i + 1
      }
      if i > num_start {
        mult = match parse_double(formula[num_start:i].to_owned()) {
          Some(v) => v
          None => 1.0
        }
      }
      for k in inner.keys() {
        let prev = match out.get(k) {
          Some(v) => v
          None => 0.0
        }
        out[k] = prev + inner[k] * mult
      }
      continue
    }
    if cobra3_is_element_start(c) {
      let el_start = i
      i = i + 1
      while i < formula.length() && cobra3_is_element_char(formula[i].to_int()) {
        i = i + 1
      }
      let element = formula[el_start:i].to_owned()
      let num_start = i
      while i < formula.length() &&
            formula[i].to_int() >= 48 &&
            formula[i].to_int() <= 57 {
        i = i + 1
      }
      let n = if i > num_start {
        match parse_double(formula[num_start:i].to_owned()) {
          Some(v) => v
          None => 1.0
        }
      } else {
        1.0
      }
      let prev = match out.get(element) {
        Some(v) => v
        None => 0.0
      }
      out[element] = prev + n
      continue
    }
    i = i + 1
  }
  (out, i)
}

///|
/// 解析分子式元素计数，如 "C6H12O6"、"Ca(OH)2"；不支持电荷/同位素.
pub fn cobra3_parse_formula(formula : String) -> Map[String, Double] {
  let (out, _) = cobra3_parse_formula_rec(formula, 0, false)
  out
}

///|
pub fn cobra3_metabolite_formula_weight(met : Cobra3Metabolite) -> Double {
  let counts = cobra3_parse_formula(met.formula)
  let mut weight = 0.0
  for element in counts.keys() {
    weight = weight + counts[element] * cobra3_element_weight(element)
  }
  weight
}

///|
fn cobra3_element_weight(element : String) -> Double {
  if element == "H" {
    1.008
  } else if element == "C" {
    12.011
  } else if element == "N" {
    14.007
  } else if element == "O" {
    15.999
  } else if element == "P" {
    30.974
  } else if element == "S" {
    32.06
  } else if element == "Na" {
    22.99
  } else if element == "K" {
    39.098
  } else if element == "Mg" {
    24.305
  } else if element == "Ca" {
    40.078
  } else if element == "Fe" {
    55.845
  } else if element == "Cl" {
    35.45
  } else if element == "Zn" {
    65.38
  } else if element == "Mn" {
    54.938
  } else if element == "Cu" {
    63.546
  } else if element == "Co" {
    58.933
  } else if element == "Se" {
    78.971
  } else {
    0.0
  }
}

///|
pub struct Cobra3MassBalance {
  reaction_id : String
  balanced : Bool
  detail : String
}

///|
/// 单反应质量/电荷平衡检查（只对提供 formula 的代谢物）.
pub fn cobra3_check_reaction_mass_balance(
  model : Cobra3Model,
  reaction_id : String,
) -> Cobra3MassBalance? {
  let idx = cobra3_model_reaction_index(model, reaction_id)
  if idx < 0 {
    return None
  }
  let r = model.reactions[idx]
  let elements : Map[String, Double] = Map([], capacity=8)
  let mut charge = 0.0
  let mut missing = 0
  for met_id in r.metabolites.keys() {
    let coeff = r.metabolites[met_id]
    match cobra3_model_get_metabolite(model, met_id) {
      None => missing = missing + 1
      Some(met) =>
        if met.formula.length() == 0 {
          missing = missing + 1
        } else {
          let counts = cobra3_parse_formula(met.formula)
          for el in counts.keys() {
            let prev = match elements.get(el) {
              Some(v) => v
              None => 0.0
            }
            elements[el] = prev + coeff * counts[el]
          }
          charge = charge + coeff * met.charge
        }
    }
  }
  let mut balanced = true
  for el in elements.keys() {
    let v = elements[el]
    if v > 1.0e-6 || v < -1.0e-6 {
      balanced = false
    }
  }
  if charge > 1.0e-6 || charge < -1.0e-6 {
    balanced = false
  }
  if missing > 0 {
    balanced = false
  }
  let detail = if missing > 0 {
    "missing formula for " + missing.to_string() + " metabolite(s)"
  } else if balanced {
    "balanced"
  } else {
    "element/charge imbalance"
  }
  Some({ reaction_id, balanced, detail, })
}

///|
pub fn cobra3_check_all_mass_balance(
  model : Cobra3Model,
) -> Array[Cobra3MassBalance] {
  let out : Array[Cobra3MassBalance] = Array(capacity=model.reactions.length())
  for r in model.reactions {
    match cobra3_check_reaction_mass_balance(model, r.id) {
      Some(b) => out.push(b)
      None => ()
    }
  }
  out
}

// ---------------------------------------------------------------------------
// 工具
// ---------------------------------------------------------------------------

///|
fn cobra3_metabolite_compartment_from_id(id : String) -> String {
  let mut bracket_start = -1
  let mut bracket_end = -1
  for i in 0..<id.length() {
    let c = id[i].to_int()
    if c == 91 {
      bracket_start = i
    }
    if c == 93 {
      bracket_end = i
    }
  }
  if bracket_start >= 0 && bracket_end > bracket_start + 1 {
    id[bracket_start + 1:bracket_end].to_owned()
  } else {
    "c"
  }
}

///|
/// 字典序插入排序（小数组足够）.
pub fn cobra3_sort_strings(arr : Array[String]) -> Unit {
  for i in 1..<arr.length() {
    let key = arr[i]
    let mut j = i - 1
    while j >= 0 && arr[j] > key {
      arr[j + 1] = arr[j]
      j = j - 1
    }
    arr[j + 1] = key
  }
}

///|
pub fn cobra3_sort_doubles(arr : Array[Double]) -> Unit {
  for i in 1..<arr.length() {
    let key = arr[i]
    let mut j = i - 1
    while j >= 0 && arr[j] > key {
      arr[j + 1] = arr[j]
      j = j - 1
    }
    arr[j + 1] = key
  }
}
