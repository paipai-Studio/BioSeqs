// cobra3_reporter.mbt — Reporter metabolites（Patil & Nielsen 2005,
// cameo `reporter_metabolites`）
//
// 思路：基因表达 z-score 经 GPR 汇到反应（反应取其 GPR 成员基因中 |z| 最
// 大者的带符号 z）→ 每个代谢物聚合相邻反应的分数（取 |z| 最极端者的带符
// 号值）→ 显著性用超几何检验：代谢物邻域基因池中"显著基因"（|z| ≥
// z_cutoff）的富集程度相对全模型基因池的期望，BH 校正后输出热点排名。

///|
/// ln Γ(x)（Lanczos g=7 近似），供超几何检验的对数组合数使用。
fn cobra3_lgamma(x : Double) -> Double {
  let g = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
    1.5056327351493116e-7,
  ]
  if x < 0.5 {
    let s = @math.sin(@math.PI * x).abs()
    if s <= 0.0 {
      return 1.0e300
    }
    @math.ln(@math.PI / s) - cobra3_lgamma(1.0 - x)
  } else {
    let x1 = x - 1.0
    let mut a = g[0]
    for i in 1..<9 {
      a = a + g[i] / (x1 + i.to_double())
    }
    let t = x1 + 7.5
    0.5 * @math.ln(2.0 * @math.PI) + (x1 + 0.5) * @math.ln(t) - t + @math.ln(a)
  }
}

///|
/// ln C(n, k)；k 越界时返回 -1e300（对应概率 0）。
fn cobra3_ln_comb(n : Int32, k : Int32) -> Double {
  if k < 0 || k > n || n < 0 {
    return -1.0e300
  }
  cobra3_lgamma((n + 1).to_double()) -
  cobra3_lgamma((k + 1).to_double()) -
  cobra3_lgamma((n - k + 1).to_double())
}

///|
/// 超几何分布上侧尾概率 P(X >= k)，X ~ Hypergeom(N, K, n)。
fn cobra3_hyper_tail(
  n_total : Int32,
  k_sig : Int32,
  n_draw : Int32,
  k : Int32,
) -> Double {
  if k > n_draw || k > k_sig {
    return 0.0
  }
  if k <= 0 {
    return 1.0
  }
  let ln_denom = cobra3_ln_comb(n_total, n_draw)
  let upper = if n_draw < k_sig { n_draw } else { k_sig }
  let hi = upper + 1
  let mut p = 0.0
  for i in k..<hi {
    let ln_num = cobra3_ln_comb(k_sig, i) +
      cobra3_ln_comb(n_total - k_sig, n_draw - i)
    p = p + @math.exp(ln_num - ln_denom)
  }
  if p > 1.0 {
    return 1.0
  }
  p
}

///|
pub(all) struct Cobra3ReporterMetabolite {
  metabolite_id : String
  /// 相邻反应分数中 |z| 最极端者的带符号值
  reporter_score : Double
  /// 有基因证据的相邻反应数
  n_reactions : Int32
  /// 邻域基因池大小（相邻反应 GPR 成员基因并集中有 z 值者）
  n_genes : Int32
  /// 邻域基因池中 |z| ≥ z_cutoff 的基因数
  n_sig_genes : Int32
  p_value : Double
  /// Benjamini-Hochberg FDR 校正
  mut p_adj : Double
  /// 贡献分数的相邻反应 ID（按发现顺序）
  reactions : Array[String]
}

///|
/// Reporter metabolites：把 `gene_z`（基因 → 表达变化 z-score）经 GPR 汇到
/// 反应、再聚合到代谢物节点，用超几何检验评估每个代谢物邻域基因池中显著
/// 基因（|z| ≥ `z_cutoff`）的富集显著性。无 GPR 或 GPR 基因全部缺少 z 值
/// 的反应不参与打分；没有任何打分相邻反应的代谢物不出现在输出中。输出按
/// p_value 升序排列（并列按代谢物 ID）。
pub fn cobra3_reporter_metabolites(
  model : Cobra3Model,
  gene_z : Map[String, Double],
  z_cutoff~ : Double = 2.0,
) -> Array[Cobra3ReporterMetabolite] {
  // 第一步：反应分数 = GPR 成员基因（有 z 值者）中 |z| 最大者的带符号 z
  let rxn_score : Map[String, Double] = Map(
    [],
    capacity=model.reactions.length(),
  )
  let rxn_genes : Map[String, Array[String]] = Map(
    [],
    capacity=model.reactions.length(),
  )
  let gene_pool : Map[String, Bool] = Map([])
  let mut n_total = 0
  for r in model.reactions {
    match r.gpr {
      None => continue
      Some(gpr) => {
        let mut best = 0.0
        let mut found = false
        let scored : Array[String] = Array(capacity=0)
        for g in cobra3_gpr_genes(gpr) {
          match gene_z.get(g) {
            None => continue
            Some(z) => {
              scored.push(g)
              if not(gene_pool.contains(g)) {
                gene_pool.set(g, true)
                n_total = n_total + 1
              }
              if not(found) || z.abs() > best.abs() {
                best = z
                found = true
              }
            }
          }
        }
        if found {
          rxn_score.set(r.id, best)
          rxn_genes.set(r.id, scored)
        }
      }
    }
  }
  if rxn_score.is_empty() {
    return []
  }
  // 全模型打分基因池中的显著基因数
  let mut k_pool = 0
  for g in gene_pool.keys() {
    if gene_z.get(g).unwrap().abs() >= z_cutoff {
      k_pool = k_pool + 1
    }
  }
  // 第二步：代谢物聚合（相邻打分反应 + 邻域基因池）
  let met_ids : Map[String, Bool] = Map([])
  for r in model.reactions {
    if rxn_score.contains(r.id) {
      for kv in r.metabolites {
        met_ids.set(kv.0, true)
      }
    }
  }
  let results : Array[Cobra3ReporterMetabolite] = Array(capacity=0)
  for kv in met_ids {
    let met = kv.0
    let mut score = 0.0
    let mut has_score = false
    let mut n_rxn = 0
    let mut n_local = 0
    let genes : Map[String, Bool] = Map([])
    let rxns : Array[String] = Array(capacity=0)
    for r in model.reactions {
      match rxn_score.get(r.id) {
        None => continue
        Some(z) =>
          if r.metabolites.contains(met) {
            rxns.push(r.id)
            n_rxn = n_rxn + 1
            if not(has_score) || z.abs() > score.abs() {
              score = z
              has_score = true
            }
            for g in rxn_genes.get(r.id).unwrap() {
              if not(genes.contains(g)) {
                genes.set(g, true)
                n_local = n_local + 1
              }
            }
          }
      }
    }
    if n_rxn == 0 {
      continue
    }
    // 第三步：邻域超几何富集
    let mut k_local = 0
    for g in genes.keys() {
      if gene_z.get(g).unwrap().abs() >= z_cutoff {
        k_local = k_local + 1
      }
    }
    let p = cobra3_hyper_tail(n_total, k_pool, n_local, k_local)
    results.push({
      metabolite_id: met,
      reporter_score: score,
      n_reactions: n_rxn,
      n_genes: n_local,
      n_sig_genes: k_local,
      p_value: p,
      p_adj: 1.0,
      reactions: rxns,
    })
  }
  // 第四步：BH 校正 + 按 p_value 升序（并列按 ID）
  let order : Array[Int32] = Array(capacity=results.length())
  for i in 0..<results.length() {
    order.push(i)
  }
  // 简单插入排序：按 (p_value, metabolite_id)
  for i in 1..<order.length() {
    let cur = order[i]
    let mut j = i - 1
    while j >= 0 {
      let a = results[order[j]]
      let b = results[cur]
      let less = a.p_value > b.p_value ||
        (a.p_value == b.p_value && a.metabolite_id > b.metabolite_id)
      if less {
        order[j + 1] = order[j]
        j = j - 1
      } else {
        break
      }
    }
    order[j + 1] = cur
  }
  let m = order.length()
  let mut running_min = 1.0e300
  for t = 0; t < m; t = t + 1 {
    let idx = order[m - 1 - t]
    let rank = m - t
    let adj = results[idx].p_value * m.to_double() / rank.to_double()
    if adj < running_min {
      running_min = adj
    }
    results[idx].p_adj = running_min.min(1.0)
  }
  let sorted : Array[Cobra3ReporterMetabolite] = Array(capacity=m)
  for t in 0..<m {
    sorted.push(results[order[t]])
  }
  sorted
}
