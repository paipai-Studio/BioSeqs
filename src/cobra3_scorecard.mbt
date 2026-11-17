// cobra3_scorecard.mbt — Memote 式模型质量评分卡（memote scorecard 语义）
//
// 在 cobra_qc_check（逐条列问题）之上补"打分"层：对模型结构完备性按指标
// 打 0..1 分（score 越高越好），加权汇总为总分并分档，输出结构化
// scorecard；支持与参考模型评分卡（如手工精修重建 vs 草稿重建）逐指标
// 比较。
//
// 指标口径（简化实现，边界在注释中说明）：
// - metabolite_coverage  1 - 孤立代谢物比例（不参与任何反应）
// - connectivity         代谢物平均连接度（参与反应数），按阈值 8 归一
// - deadend_free         1 - 死端代谢物比例（拓扑上只产不耗或只耗不产，
//                        且不经 exchange 连接系统边界）
// - gpr_completeness     内部反应（非交换）带 GPR 注释的比例
// - mass_balance         内部反应质量/电荷平衡比例（依赖 formula 注释，
//                        缺 formula 的代谢物计入不平衡）
// - blocked_free         1 - 口令（blocked）反应比例（FVA 上下界全为 0、
//                        不可载流的内部反应；exchange 只统计有代谢物的，
//                        关闭的 uptake 属建模选择而非缺陷）

///|
/// 单指标评分条目.
pub(all) struct Cobra3ScorecardMetric {
  /// 指标 ID（如 "gpr_completeness"）
  id : String
  /// 所属板块（annotation / consistency / network）
  section : String
  /// 0..1 得分（越高越好；total=0 时记 0）
  score : Double
  /// 达标元素数（比例型指标）
  passed : Int32
  /// 参与统计的元素数（连续型指标为 0，均值在 detail）
  total : Int32
  /// 人读明细
  detail : String
  /// 加权总分时的权重
  weight : Double
}

///|
/// 结构化评分卡.
pub(all) struct Cobra3Scorecard {
  model_id : String
  /// 指标条目（按 weight 降序）
  metrics : Array[Cobra3ScorecardMetric]
  /// 加权总分（权重和归一后 0..1）
  total_score : Double
  /// 分档：A >= 0.9，B >= 0.75，C >= 0.6，D >= 0.4，否则 E
  grade : String
}

///|
/// 单指标比较结果.
pub(all) struct Cobra3ScorecardDiff {
  id : String
  /// 待评模型得分
  score : Double
  /// 参考模型得分
  reference : Double
  /// score - reference
  delta : Double
  /// delta > 0
  improved : Bool
}

///|
/// 评分卡比较（模型基准比较）.
pub(all) struct Cobra3ScorecardComparison {
  /// 两个评分卡共有指标的逐条比较（按 delta 升序，最差在前）
  diffs : Array[Cobra3ScorecardDiff]
  total_score : Double
  reference_score : Double
  /// total_score - reference_score
  delta : Double
}

///|
/// 简化版 memote 评分卡：对模型计算六个质量指标并加权汇总。空模型或全部
/// 指标无统计对象时各指标记 0 分。
pub fn cobra3_scorecard(model : Cobra3Model) -> Cobra3Scorecard {
  let n_rxn = model.reactions.length()
  // 代谢物拓扑统计：参与反应数 / 生产次数 / 消耗次数
  let part : Map[String, Int32] = Map([], capacity=model.metabolites.length())
  let produced : Map[String, Int32] = Map(
    [],
    capacity=model.metabolites.length(),
  )
  let consumed : Map[String, Int32] = Map(
    [],
    capacity=model.metabolites.length(),
  )
  for r in model.reactions {
    for met_id in r.metabolites.keys() {
      let coeff = r.metabolites[met_id]
      let prev = match part.get(met_id) {
        Some(v) => v
        None => 0
      }
      part[met_id] = prev + 1
      if coeff > 0.0 {
        let p = match produced.get(met_id) {
          Some(v) => v
          None => 0
        }
        produced[met_id] = p + 1
      } else if coeff < 0.0 {
        let c = match consumed.get(met_id) {
          Some(v) => v
          None => 0
        }
        consumed[met_id] = c + 1
      }
    }
  }
  let n_mets = model.metabolites.length()
  // 孤立代谢物
  let mut n_orphan = 0
  for m in model.metabolites {
    if not(part.contains(m.id)) {
      n_orphan = n_orphan + 1
    }
  }
  // 连接度（只统计有参与的代谢物）
  let mut sum_conn = 0.0
  for v in part.values() {
    sum_conn = sum_conn + v.to_double()
  }
  let n_connected = part.length()
  let avg_conn = if n_connected > 0 {
    sum_conn / n_connected.to_double()
  } else {
    0.0
  }
  // 死端代谢物：只产不耗或只耗不产，且不经 exchange 连接系统边界
  // （memote 语义：边界代谢物可摄取/分泌，不算死端）
  let ex_ids : Map[String, Bool] = Map([], capacity=n_rxn)
  for r in model.reactions {
    if cobra3_model_is_exchange(model, r) {
      ex_ids[r.id] = true
    }
  }
  let in_exchange : Map[String, Bool] = Map([], capacity=n_mets)
  for r in model.reactions {
    if ex_ids.contains(r.id) {
      for met_id in r.metabolites.keys() {
        in_exchange[met_id] = true
      }
    }
  }
  let mut n_deadend = 0
  for met_id in part.keys() {
    if in_exchange.contains(met_id) {
      continue
    }
    let np = match produced.get(met_id) {
      Some(v) => v
      None => 0
    }
    let nc = match consumed.get(met_id) {
      Some(v) => v
      None => 0
    }
    if (np > 0 && nc == 0) || (nc > 0 && np == 0) {
      n_deadend = n_deadend + 1
    }
  }
  // 内部反应（非 exchange、非空）：GPR / 质量平衡 / blocked
  let mut n_internal = 0
  let mut n_with_gpr = 0
  let mut n_balanced = 0
  let mut n_blocked = 0
  let fva = cobra3_fva(model, fraction_of_optimum=0.0)
  let fva_min : Map[String, Double] = Map([], capacity=n_rxn)
  let fva_max : Map[String, Double] = Map([], capacity=n_rxn)
  for f in fva {
    fva_min[f.reaction_id] = f.minimum
    fva_max[f.reaction_id] = f.maximum
  }
  for r in model.reactions {
    if cobra3_model_is_exchange(model, r) || r.metabolites.length() == 0 {
      continue
    }
    n_internal = n_internal + 1
    match r.gpr {
      Some(_) => n_with_gpr = n_with_gpr + 1
      None => ()
    }
    match cobra3_check_reaction_mass_balance(model, r.id) {
      Some(b) => if b.balanced { n_balanced = n_balanced + 1 }
      None => ()
    }
    let mn = match fva_min.get(r.id) {
      Some(v) => v
      None => 0.0
    }
    let mx = match fva_max.get(r.id) {
      Some(v) => v
      None => 0.0
    }
    if mn.abs() <= 1.0e-9 && mx.abs() <= 1.0e-9 {
      n_blocked = n_blocked + 1
    }
  }
  // 组装指标（score 为 0 的口径：total=0 表示无统计对象）
  let ratio = fn(passed : Int32, total : Int32) -> Double {
    if total == 0 {
      0.0
    } else {
      passed.to_double() / total.to_double()
    }
  }
  let coverage = ratio(n_mets - n_orphan, n_mets)
  let deadend_total = n_connected
  let deadend_free = ratio(deadend_total - n_deadend, deadend_total)
  let gpr_score = ratio(n_with_gpr, n_internal)
  let mass_score = ratio(n_balanced, n_internal)
  let blocked_free = ratio(n_internal - n_blocked, n_internal)
  let conn_score = if n_connected == 0 {
    0.0
  } else {
    let s = avg_conn / 8.0
    if s > 1.0 {
      1.0
    } else {
      s
    }
  }
  let metrics : Array[Cobra3ScorecardMetric] = Array(capacity=6)
  metrics.push({
    id: "gpr_completeness",
    section: "annotation",
    score: gpr_score,
    passed: n_with_gpr,
    total: n_internal,
    detail: "\{n_with_gpr}/\{n_internal} internal reactions have GPR",
    weight: 0.15,
  })
  metrics.push({
    id: "mass_balance",
    section: "consistency",
    score: mass_score,
    passed: n_balanced,
    total: n_internal,
    detail: "\{n_balanced}/\{n_internal} internal reactions mass/charge balanced (formula-dependent)",
    weight: 0.20,
  })
  metrics.push({
    id: "blocked_free",
    section: "consistency",
    score: blocked_free,
    passed: n_internal - n_blocked,
    total: n_internal,
    detail: "\{n_blocked} blocked internal reactions (FVA range = 0)",
    weight: 0.20,
  })
  metrics.push({
    id: "metabolite_coverage",
    section: "network",
    score: coverage,
    passed: n_mets - n_orphan,
    total: n_mets,
    detail: "\{n_orphan} orphan metabolites",
    weight: 0.15,
  })
  metrics.push({
    id: "deadend_free",
    section: "network",
    score: deadend_free,
    passed: deadend_total - n_deadend,
    total: deadend_total,
    detail: "\{n_deadend} dead-end metabolites (produced-only or consumed-only, exchange-connected excluded)",
    weight: 0.15,
  })
  metrics.push({
    id: "connectivity",
    section: "network",
    score: conn_score,
    passed: 0,
    total: 0,
    detail: "avg metabolite connectivity = \{avg_conn}",
    weight: 0.15,
  })
  // 加权总分（权重和恒为 1）
  let mut total = 0.0
  for m in metrics {
    total += m.score * m.weight
  }
  let grade = if total >= 0.9 {
    "A"
  } else if total >= 0.75 {
    "B"
  } else if total >= 0.6 {
    "C"
  } else if total >= 0.4 {
    "D"
  } else {
    "E"
  }
  { model_id: model.id, metrics, total_score: total, grade, }
}

///|
/// 模型基准比较：待评模型与参考模型的评分卡逐指标对比（只比较两者共有
/// 的指标），按 delta 升序（最差在前）；总分差一并列出。
pub fn cobra3_scorecard_compare(
  card : Cobra3Scorecard,
  reference : Cobra3Scorecard,
) -> Cobra3ScorecardComparison {
  let ref_score : Map[String, Double] = Map([], capacity=8)
  for m in reference.metrics {
    ref_score[m.id] = m.score
  }
  let diffs : Array[Cobra3ScorecardDiff] = Array(capacity=8)
  for m in card.metrics {
    match ref_score.get(m.id) {
      Some(rs) =>
        diffs.push({
          id: m.id,
          score: m.score,
          reference: rs,
          delta: m.score - rs,
          improved: m.score > rs + 1.0e-9,
        })
      None => ()
    }
  }
  // delta 升序（最差在前）
  for i in 1..<diffs.length() {
    let cur = diffs[i]
    let mut j = i - 1
    while j >= 0 && diffs[j].delta > cur.delta {
      diffs[j + 1] = diffs[j]
      j = j - 1
    }
    diffs[j + 1] = cur
  }
  {
    diffs,
    total_score: card.total_score,
    reference_score: reference.total_score,
    delta: card.total_score - reference.total_score,
  }
}
