// kegg_fba.mbt — KEGG → FBA 草稿模型重构（ModelSEED / CarveME 思路）
//
// 打通 kegg/kgml 与 cobra 的最大断链：
// 1. 反应式解析：KEGG 方程字符串（"2 C00031 + C00118 <=> ..."）与 KGML
//    pathway 的 <reaction> 集合两种来源
// 2. 编译为 FbaModel：化合物 → 代谢物池（去重排序），反应 → 计量列
//    （可逆性 → 上下界），entry 基因（ko:Kxxxxx）→ GPR 简化的 genes
//    列表；可选开边界 exchange（ModelSEED 式边界）
// 3. 质量平衡检查：拓扑健全性（空列 / 单侧列）+ 可选 compound formula
//    元素配平
// 4. gapfill 闭环：草稿模型即 FbaModel，直接接入 cobra_gapfill 补缺失
//    反应（如产物外排 exchange）恢复生长
//
// 简化口径：KGML 不带化学计量系数（全部按 1 处理）；KEGG 方程系数支持
// 数字前缀；元素配平只统计参与者 formula 齐全的反应。

///|
/// 质量平衡问题条目.
pub(all) struct KeggFbaBalanceIssue {
  reaction_id : String
  /// "empty"（无参与物）| "one_sided"（单侧反应）| "imbalanced"（元素不守恒）
  kind : String
  /// 人读明细
  detail : String
}

///|
/// 草稿重建结果.
pub(all) struct KeggFbaDraft {
  model : FbaModel
  issues : Array[KeggFbaBalanceIssue]
  n_reactions : Int32
  n_metabolites : Int32
  n_exchanges : Int32
}

///|
/// 去掉 KEGG 常见前缀（cpd: / ko: / rn:）.
fn kegg_fba_strip_prefix(id : String) -> String {
  for p in ["cpd:", "ko:", "rn:"] {
    if id.has_prefix(p) {
      return id[p.length():].to_owned()
    }
  }
  id
}

///|
/// 拆分 KGML gene entry 的 name 属性为基因 ID 列表。name 可为空格
/// 或逗号（含无空格逗号形式）分隔的多个 KEGG ID，如
/// "ko:K00844 ko:K01623" / "ko:K00844,ko:K01623"；逐个剥离前缀后返回。
fn kegg_fba_gene_ids(name : String) -> Array[String] {
  let out : Array[String] = []
  for tok in kegg_fba_tokens(name) {
    // 同一 token 内的逗号也是分隔符（"A,B" 无空格形式）：二次切分
    let mut start = 0
    let n = tok.length()
    for i = 0; i < n; i = i + 1 {
      if tok[i:i + 1].to_owned() == "," {
        if i > start {
          out.push(kegg_fba_strip_prefix(tok[start:i].to_owned()))
        }
        start = i + 1
      }
    }
    if n > start {
      out.push(kegg_fba_strip_prefix(tok[start:n].to_owned()))
    }
  }
  out
}

///|
/// 合并底物/产物为净计量：同一化合物两侧同现或同侧重复时系数累加，
/// 按首次出现顺序返回剔除净零项后的 (化合物, 系数) 列表。
/// 底层 FbaModel::add_reaction 对重复参与者是覆盖赋值而非累加，
/// 因此所有构造参与者列表的路径都必须先经过本函数合并。
fn kegg_fba_merge_participants(
  subs : Array[(String, Double)],
  prods : Array[(String, Double)],
) -> Array[(String, Double)] {
  let order : Array[String] = Array(capacity=subs.length() + prods.length())
  let coeff : Map[String, Double] = Map([])
  for p in subs {
    let v = -p.1
    match coeff.get(p.0) {
      Some(prev) => coeff[p.0] = prev + v
      None => {
        coeff[p.0] = v
        order.push(p.0)
      }
    }
  }
  for p in prods {
    let v = p.1
    match coeff.get(p.0) {
      Some(prev) => coeff[p.0] = prev + v
      None => {
        coeff[p.0] = v
        order.push(p.0)
      }
    }
  }
  let out : Array[(String, Double)] = Array(capacity=order.length())
  for k in order {
    let v = coeff[k]
    if v.abs() > 1.0e-9 {
      out.push((k, v))
    }
  }
  out
}

///|
/// 检查反应 ID 是否已存在（防重复列）.
fn kegg_fba_has_reaction(model : FbaModel, id : String) -> Bool {
  for r in model.reactions {
    if r.id == id {
      return true
    }
  }
  false
}

///|
/// 按空白切分字符串（跳过连续空格）.
fn kegg_fba_tokens(s : String) -> Array[String] {
  let t = s.trim().to_owned()
  let out : Array[String] = Array(capacity=2)
  let mut pos = 0
  while pos < t.length() {
    while pos < t.length() && t[pos:pos + 1].to_owned() == " " {
      pos = pos + 1
    }
    if pos >= t.length() {
      break
    }
    match kgml_find_from(t, " ", pos) {
      Some(sp) => {
        out.push(t[pos:sp].to_owned())
        pos = sp + 1
      }
      None => {
        out.push(t[pos:].to_owned())
        break
      }
    }
  }
  out
}

///|
/// 解析单侧反应式："2 C00031 + C00118" → [("C00031", 2), ("C00118", 1)].
/// 空侧或出现非法片段返回 None。
fn kegg_fba_parse_side(s : String) -> Array[(String, Double)]? {
  let trimmed = s.trim().to_owned()
  if trimmed.length() == 0 {
    return None
  }
  let out : Array[(String, Double)] = Array(capacity=4)
  let mut pos = 0
  let mut done = false
  while !done {
    match kgml_find_from(trimmed, "+", pos) {
      Some(idx) => {
        match kegg_fba_parse_term(trimmed[pos:idx].to_owned()) {
          Some(t) => out.push(t)
          None => return None
        }
        pos = idx + 1
      }
      None => {
        match kegg_fba_parse_term(trimmed[pos:].to_owned()) {
          Some(t) => out.push(t)
          None => return None
        }
        done = true
      }
    }
  }
  if out.length() == 0 {
    None
  } else {
    Some(out)
  }
}

///|
/// 解析单个反应项："2 C00031" / "C00031" / "cpd:C00031" → (id, 系数).
fn kegg_fba_parse_term(term : String) -> (String, Double)? {
  let t = term.trim().to_owned()
  if t.length() == 0 {
    return None
  }
  match kgml_find_from(t, " ", 0) {
    Some(sp) => {
      let head = t[0:sp].trim().to_owned()
      match parse_double(head) {
        Some(c) => {
          let rest = t[sp:].trim().to_owned()
          if rest.length() == 0 {
            None
          } else {
            Some((kegg_fba_strip_prefix(rest), c))
          }
        }
        None => Some((kegg_fba_strip_prefix(t), 1.0))
      }
    }
    None => Some((kegg_fba_strip_prefix(t), 1.0))
  }
}

///|
/// 解析 KEGG 反应方程："2 C00031 + C00118 <=> C00197 + C00022"。
/// 分隔符按 "<=>" / "=>" / "->" / "=" 识别（最长优先），化合物前缀
/// "cpd:" 自动剥离；底物侧返回负系数由调用方处理。无法解析返回 None。
pub fn kegg_reaction_parse(
  equation : String,
) -> (Array[(String, Double)], Array[(String, Double)])? {
  let s = equation.trim().to_owned()
  for sep in ["<=>", "=>", "->", "="] {
    match kgml_find_from(s, sep, 0) {
      Some(i) => {
        let left = kegg_fba_parse_side(s[0:i].to_owned())
        let right = kegg_fba_parse_side(s[i + sep.length():].to_owned())
        return match (left, right) {
          (Some(l), Some(r)) => Some((l, r))
          _ => None
        }
      }
      None => continue
    }
  }
  None
}

///|
/// 扩展代谢物池（新代谢物追加零行）.
fn kegg_fba_extend_mets(model : FbaModel, extra : Array[String]) -> FbaModel {
  let mets : Array[String] = Array(
    capacity=model.metabolites.length() + extra.length(),
  )
  for m in model.metabolites {
    mets.push(m)
  }
  for m in extra {
    if !mets.contains(m) {
      mets.push(m)
    }
  }
  let n_r = model.reactions.length()
  let stoich : Array[Array[Double]] = Array(capacity=mets.length())
  for i in 0..<model.metabolites.length() {
    let row : Array[Double] = Array(capacity=n_r)
    for j in 0..<n_r {
      row.push(model.stoichiometry[i][j])
    }
    stoich.push(row)
  }
  for _i in 0..<(mets.length() - model.metabolites.length()) {
    stoich.push(Array::make(n_r, 0.0))
  }
  { metabolites: mets, reactions: model.reactions, stoichiometry: stoich, }
}

///|
/// 向模型追加一条 KEGG 方程反应（自动扩展新代谢物；底物取负系数）。
/// 同一化合物在方程中重复出现（同侧重复或两侧同现）时系数累加合并，
/// 合并后净参与物为空（如自转移 "A -> A"）则报错；反应 ID 与模型中
/// 已有反应重复时报错。
pub fn kegg_fba_add_equation(
  model : FbaModel,
  rxn_id : String,
  equation : String,
  reversible~ : Bool = false,
  objective~ : Double = 0.0,
  genes~ : Array[String] = [],
) -> FbaModel raise FbaError {
  if kegg_fba_has_reaction(model, rxn_id) {
    raise FbaError("reaction already exists: \{rxn_id}")
  }
  match kegg_reaction_parse(equation) {
    None => raise FbaError("invalid KEGG equation: \{equation}")
    Some((subs, prods)) => {
      // 合并重复化合物：系数累加，避免底层按列覆盖赋值丢失计量
      let participants = kegg_fba_merge_participants(subs, prods)
      if participants.length() == 0 {
        raise FbaError(
          "KEGG equation resolves to zero net participants: \{equation}",
        )
      }
      let extra : Array[String] = Array(capacity=participants.length())
      for p in participants {
        extra.push(p.0)
      }
      let m = kegg_fba_extend_mets(model, extra)
      let (lb, ub) = if reversible { (-1000.0, 1000.0) } else { (0.0, 1000.0) }
      m.add_reaction(
        {
          id: rxn_id,
          name: rxn_id,
          lower_bound: lb,
          upper_bound: ub,
          objective,
          genes,
        },
        participants,
      )
    }
  }
}

///|
/// 追加一条 exchange 反应（EX_<met>："met ->"，下界默认 -10 允许摄取），
/// 代谢物不存在时自动扩展。用于开边界与 gapfill universal 池构造。
/// 反应 ID（EX_<met>）与模型中已有反应重复时报错。
pub fn kegg_fba_add_exchange(
  model : FbaModel,
  met : String,
  lower_bound~ : Double = -10.0,
) -> FbaModel raise FbaError {
  let ex_id = "EX_" + met
  if kegg_fba_has_reaction(model, ex_id) {
    raise FbaError("reaction already exists: \{ex_id}")
  }
  let m = kegg_fba_extend_mets(model, [met])
  m.add_reaction(
    {
      id: ex_id,
      name: ex_id,
      lower_bound,
      upper_bound: 1000.0,
      objective: 0.0,
      genes: [],
    },
    [(met, -1.0)],
  )
}

///|
/// 质量平衡检查（ModelSEED 式）：对每条非 exchange 反应，
/// 1) 无参与物 → "empty"；2) 所有计量同号 → "one_sided"；
/// 3) 提供 formulas 时（Map 非空），参与者 formula 齐全的反应按元素
/// 守恒检查，不守恒 → "imbalanced"（detail 列出净增减元素）。
/// formula 缺失的反应跳过元素检查。
pub fn kegg_fba_check_mass_balance(
  model : FbaModel,
  formulas~ : Map[String, String] = Map([]),
) -> Array[KeggFbaBalanceIssue] {
  let issues : Array[KeggFbaBalanceIssue] = Array(capacity=0)
  for j in 0..<model.reactions.length() {
    let r = model.reactions[j]
    if r.id.has_prefix("EX_") {
      continue
    }
    let participants : Array[(String, Double)] = Array(capacity=4)
    for i in 0..<model.metabolites.length() {
      let v = model.stoichiometry[i][j]
      if v.abs() > 1.0e-9 {
        participants.push((model.metabolites[i], v))
      }
    }
    if participants.length() == 0 {
      issues.push({
        reaction_id: r.id,
        kind: "empty",
        detail: "reaction has no participants",
      })
      continue
    }
    let mut all_pos = true
    let mut all_neg = true
    for p in participants {
      if p.1 <= 0.0 {
        all_pos = false
      }
      if p.1 >= 0.0 {
        all_neg = false
      }
    }
    if all_pos || all_neg {
      issues.push({
        reaction_id: r.id,
        kind: "one_sided",
        detail: "all participants on one side",
      })
      continue
    }
    if formulas.length() == 0 {
      continue
    }
    let elements : Map[String, Double] = Map([])
    let mut covered = true
    for p in participants {
      match formulas.get(p.0) {
        Some(f) =>
          if f.length() == 0 {
            covered = false
          } else {
            let counts = cobra3_parse_formula(f)
            for el in counts.keys() {
              let prev = match elements.get(el) {
                Some(v) => v
                None => 0.0
              }
              elements[el] = prev + p.1 * counts[el]
            }
          }
        None => covered = false
      }
    }
    if !covered {
      continue
    }
    let mut imbalanced = false
    let mut detail = ""
    for el in elements.keys() {
      let v = elements[el]
      if v.abs() > 1.0e-6 {
        imbalanced = true
        detail = detail + el + ":" + v.to_string() + " "
      }
    }
    if imbalanced {
      issues.push({
        reaction_id: r.id,
        kind: "imbalanced",
        detail: detail.trim().to_owned(),
      })
    }
  }
  issues
}

///|
/// KEGG → FBA 草稿重建（pipeline 入口）：编译 KGML pathway 为 FbaModel
/// 并跑质量平衡检查。代谢物按 ID 排序保证确定性；基因从 etype=gene 的
/// entry 经其 reaction 字段关联（entry.name 的多个 KO ID 按空格/逗号
/// 拆分，ko:Kxxxxx）；可逆反应边界 ±1000，不可逆 [0, 1000]；
/// `objective_reaction` 置目标系数 1；`open_exchanges` 为全部化合物开
/// exchange（下界 `exchange_lb`）。pathway 中反应 ID 重复时报错。
pub fn kegg_fba_reconstruct(
  pathway : KgmlPathway,
  open_exchanges~ : Bool = false,
  exchange_lb~ : Double = -10.0,
  objective_reaction~ : String = "",
  formulas~ : Map[String, String] = Map([]),
) -> KeggFbaDraft raise FbaError {
  // 化合物池（去重 + 排序）
  let met_set : Map[String, Bool] = Map([])
  for rxn in pathway.reactions {
    for s in rxn.substrates {
      met_set[kegg_fba_strip_prefix(s)] = true
    }
    for p in rxn.products {
      met_set[kegg_fba_strip_prefix(p)] = true
    }
  }
  let mets : Array[String] = Array(capacity=met_set.length())
  for k in met_set.keys() {
    mets.push(k)
  }
  cobra3_sort_strings(mets)
  let mut model = fba_model_new(mets)
  // entry 基因 → reaction 关联（entry.reaction 可为空格分隔多个；
  // entry.name 同样可为空格/逗号分隔的多个 KO ID，逐个拆分）
  let gene_map : Map[String, Array[String]] = Map([])
  for e in pathway.get_entries_by_type("gene") {
    for gene in kegg_fba_gene_ids(e.name) {
      for tok in kegg_fba_tokens(e.reaction) {
        let rid = kegg_fba_strip_prefix(tok)
        if rid.length() == 0 {
          continue
        }
        match gene_map.get(rid) {
          Some(arr) => arr.push(gene)
          None => gene_map[rid] = [gene]
        }
      }
    }
  }
  // 编译反应
  for rxn in pathway.reactions {
    let rid = kegg_fba_strip_prefix(rxn.name)
    if kegg_fba_has_reaction(model, rid) {
      raise FbaError("duplicate reaction in pathway: \{rid}")
    }
    let reversible = rxn.rtype == "reversible"
    let (lb, ub) = if reversible { (-1000.0, 1000.0) } else { (0.0, 1000.0) }
    let raw_subs : Array[(String, Double)] = Array(
      capacity=rxn.substrates.length(),
    )
    for s in rxn.substrates {
      raw_subs.push((kegg_fba_strip_prefix(s), 1.0))
    }
    let raw_prods : Array[(String, Double)] = Array(
      capacity=rxn.products.length(),
    )
    for p in rxn.products {
      raw_prods.push((kegg_fba_strip_prefix(p), 1.0))
    }
    // 同一化合物双侧同现（如 A + B → A + C）时先合并，避免底层覆盖赋值
    let participants = kegg_fba_merge_participants(raw_subs, raw_prods)
    if participants.length() == 0 {
      raise FbaError("reaction resolves to zero net participants: \{rid}")
    }
    let genes = match gene_map.get(rid) {
      Some(g) => g
      None => []
    }
    model = model.add_reaction(
      {
        id: rid,
        name: rxn.name,
        lower_bound: lb,
        upper_bound: ub,
        objective: if rid == objective_reaction {
          1.0
        } else {
          0.0
        },
        genes,
      },
      participants,
    )
  }
  let mut n_ex = 0
  if open_exchanges {
    for met in mets {
      model = kegg_fba_add_exchange(model, met, lower_bound=exchange_lb)
      n_ex = n_ex + 1
    }
  }
  let issues = kegg_fba_check_mass_balance(model, formulas~)
  {
    model,
    issues,
    n_reactions: model.reactions.length(),
    n_metabolites: mets.length(),
    n_exchanges: n_ex,
  }
}
