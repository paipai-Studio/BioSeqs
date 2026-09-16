///|
/// TFBSTools：转录因子结合位点（TFBS）分析。
///
/// 对应 Bioconductor TFBSTools + PWMEnrich：JASPAR 矩阵访问、
/// PFM→PWM 转换、序列扫描（正/反向链）、motif 富集分析
/// （超几何检验 + BH-FDR 校正）、共识序列与信息含量。

/// JASPAR 矩阵记录。
pub(all) struct TfbsJasparMatrix {
  matrix_id : String
  name : String
  collection : String
  tax_group : String
  species : String
  pfm : Array[Array[Double]]
}

/// 扫描命中结果。
pub(all) struct TfbsHit {
  sequence_id : String
  start : Int
  end : Int
  strand : String
  score : Double
  sequence : String
}

/// 富集结果。
pub(all) struct TfbsEnrichResult {
  matrix_id : String
  tf_name : String
  n_targets_hit : Int
  n_background_hit : Int
  total_targets : Int
  total_background : Int
  fold_enrichment : Double
  p_value : Double
  adj_p : Double
}

/// 构造 JASPAR 矩阵。
pub fn TfbsJasparMatrix::new(
  matrix_id : String,
  name : String,
  pfm : Array[Array[Double]],
) -> TfbsJasparMatrix {
  TfbsJasparMatrix::{
    matrix_id: matrix_id, name: name, pfm: pfm,
    collection: "CORE", tax_group: "vertebrates", species: "",
  }
}

/// 矩阵宽度。
pub fn tfbs_matrix_width(pfm : Array[Array[Double]]) -> Int {
  pfm.length()
}

/// PFM 归一化（行归一化 + 伪计数）。
pub fn tfbs_pfm_normalize(
  counts : Array[Array[Double]],
  pseudocount? : Double = 0.8,
) -> Array[Array[Double]] {
  let pseudo = pseudocount
  let out : Array[Array[Double]] = Array::new()
  for i in 0..<counts.length() {
    let mut total = 0.0
    let col : Array[Double] = Array::new()
    for j in 0..<4 {
      col.push(counts[i][j] + pseudo)
      total = total + counts[i][j] + pseudo
    }
    for j in 0..<4 {
      col[j] = if total > 0.0 { col[j] / total } else { 0.25 }
    }
    out.push(col)
  }
  out
}

/// 模块内 ln2 常量。
let ln2 : Double = 0.6931471805599453

/// PFM → PWM（log-odds，背景默认均匀 0.25）。
pub fn tfbs_pfm_to_pwm(
  pfm : Array[Array[Double]],
  bg? : Array[Double],
) -> Array[Array[Double]] {
  let bg = if bg is Some(_) { bg.unwrap() } else { [0.25, 0.25, 0.25, 0.25] }
  let out : Array[Array[Double]] = Array::new()
  for i in 0..<pfm.length() {
    let col : Array[Double] = Array::new()
    for j in 0..<4 {
      let p = pfm[i][j]
      let ratio = if bg[j] > 0.0 {
        if p > 1.0e-9 { p / bg[j] } else { 1.0e-9 / bg[j] }
      } else {
        1.0
      }
      col.push(@math.ln(ratio) / ln2)
    }
    out.push(col)
  }
  out
}

/// 碱基索引 → 字符。
fn tfbs_base_char(idx : Int) -> Char {
  if idx == 0 { return 'A' }
  if idx == 1 { return 'C' }
  if idx == 2 { return 'G' }
  'T'
}

/// 共识序列。
pub fn tfbs_consensus(pfm : Array[Array[Double]]) -> String {
  let sb = StringBuilder()
  for i in 0..<pfm.length() {
    let mut best_idx = 0
    let mut best_val = pfm[i][0]
    for j in 1..<4 {
      if pfm[i][j] > best_val {
        best_val = pfm[i][j]
        best_idx = j
      }
    }
    sb.write_char(tfbs_base_char(best_idx))
  }
  sb.to_string()
}

/// 信息含量（每位置 bits + 总 bits）。
pub fn tfbs_info_content(pfm : Array[Array[Double]]) -> (Array[Double], Double) {
  let per_pos : Array[Double] = Array::new()
  let mut total = 0.0
  for i in 0..<pfm.length() {
    let mut ic = 0.0
    for j in 0..<4 {
      let p = pfm[i][j]
      if p > 1.0e-9 {
        ic = ic + p * @math.ln(p) / ln2
      }
    }
    let bits = 2.0 + ic
    per_pos.push(bits)
    total = total + bits
  }
  (per_pos, total)
}

/// 字符 → 碱基索引（A/C/G/T，否则 -1）。
fn tfbs_base_idx(b : Char) -> Int {
  if b == 'A' || b == 'a' { return 0 }
  if b == 'C' || b == 'c' { return 1 }
  if b == 'G' || b == 'g' { return 2 }
  if b == 'T' || b == 't' || b == 'U' || b == 'u' { return 3 }
  -1
}

/// 反向互补。
fn tfbs_reverse_complement(seq : String) -> String {
  let sb = StringBuilder()
  let mut i = seq.length() - 1
  while i >= 0 {
    let b = seq[i]
    if b == 'A' || b == 'a' { sb.write_char('T') }
    else if b == 'C' || b == 'c' { sb.write_char('G') }
    else if b == 'G' || b == 'g' { sb.write_char('C') }
    else if b == 'T' || b == 't' { sb.write_char('A') }
    else if b == 'U' || b == 'u' { sb.write_char('A') }
    else { sb.write_char('N') }
    i = i - 1
  }
  sb.to_string()
}

/// 在序列上滑动窗口打分（正向 + 可选反向链）。
pub fn tfbs_scan_sequence(
  sequence : String,
  pwm : Array[Array[Double]],
  threshold? : Double = 0.0,
  both_strands? : Bool = true,
) -> Array[TfbsHit] {
  let threshold = threshold
  let both = both_strands
  let width = pwm.length()
  let seq_len = sequence.length()
  let hits : Array[TfbsHit] = Array::new()
  if width == 0 || seq_len < width {
    return hits
  }
  // 正向链
  for start in 0..<(seq_len - width + 1) {
    let mut score = 0.0
    let end = start + width
    let sub = sequence[start:end].to_owned()
    for i2 in 0..<width {
      let b = sequence.unsafe_get(start + i2).unsafe_to_char()
      let idx = tfbs_base_idx(b)
      if idx >= 0 {
        score = score + pwm[i2][idx]
      }
    }
    if score >= threshold {
      hits.push(TfbsHit::{
        sequence_id: "seq", start: start, end: end,
        strand: "+", score: score, sequence: sub,
      })
    }
  }
  // 反向链
  if both {
    let rc = tfbs_reverse_complement(sequence)
    for start in 0..<(seq_len - width + 1) {
      let mut score = 0.0
      let end = start + width
      let sub = rc[start:end].to_owned()
      for i2 in 0..<width {
        let b = rc.unsafe_get(start + i2).unsafe_to_char()
        let idx = tfbs_base_idx(b)
        if idx >= 0 {
          score = score + pwm[i2][idx]
        }
      }
      if score >= threshold {
        let orig_start = seq_len - (start + width)
        let orig_end = seq_len - start
        hits.push(TfbsHit::{
          sequence_id: "seq", start: orig_start, end: orig_end,
          strand: "-", score: score, sequence: sub,
        })
      }
    }
  }
  hits
}

/// 批量扫描多序列。
pub fn tfbs_scan_batch(
  sequences : Array[(String, String)],
  pwm : Array[Array[Double]],
  threshold? : Double = 0.0,
) -> Array[TfbsHit] {
  let threshold = threshold
  let all_hits : Array[TfbsHit] = Array::new()
  for entry in sequences {
    let (seq_id, seq) = entry
    let hits = tfbs_scan_sequence(seq, pwm, threshold=threshold, both_strands=true)
    for h in hits {
      all_hits.push(TfbsHit::{
        sequence_id: seq_id, start: h.start, end: h.end,
        strand: h.strand, score: h.score, sequence: h.sequence,
      })
    }
  }
  all_hits
}

/// 标准正态右尾概率 P(Z >= z)。
fn tfbs_pnorm_tail(z : Double) -> Double {
  let abs_z = if z < 0.0 { -z } else { z }
  let t = 1.0 / (1.0 + 0.2316419 * abs_z)
  let d = 0.3989422804014327 * @math.exp(-0.5 * z * z)
  let c1 = 0.319381530
  let c2 = -0.356563782
  let c3 = 1.781477937
  let c4 = -1.821255978
  let c5 = 1.330274429
  let poly = c1 * t + c2 * t * t + c3 * t * t * t + c4 * t * t * t * t + c5 * t * t * t * t * t
  let p = d * poly
  if z > 0.0 { p } else { 1.0 - p }
}

/// motif 富集分析（超几何检验，正态近似）。
pub fn tfbs_enrich_hypergeometric(
  matrix_id : String,
  tf_name : String,
  n_targets_hit : Int,
  n_targets_total : Int,
  n_background_hit : Int,
  n_background_total : Int,
) -> TfbsEnrichResult {
  // 合并两组
  let total_pop = n_targets_total + n_background_total
  let total_hits = n_targets_hit + n_background_hit
  // 超几何均值/方差
  let p = if total_pop > 0 {
    total_hits.to_double() / total_pop.to_double()
  } else {
    0.0
  }
  let q = 1.0 - p
  let mean = n_targets_total.to_double() * p
  let variance = n_targets_total.to_double() * p * q *
    (total_pop.to_double() - n_targets_total.to_double()) /
    (if total_pop > 1 { total_pop.to_double() - 1.0 } else { 1.0 })
  let p_value = if variance > 0.0 {
    let z = (n_targets_hit.to_double() - 0.5 - mean) / variance.sqrt()
    tfbs_pnorm_tail(z)
  } else {
    1.0
  }
  // fold enrichment
  let target_frac = if n_targets_total > 0 {
    n_targets_hit.to_double() / n_targets_total.to_double()
  } else {
    0.0
  }
  let bg_frac = if n_background_total > 0 {
    n_background_hit.to_double() / n_background_total.to_double()
  } else {
    0.0
  }
  let fold_enrichment = if bg_frac > 0.0 { target_frac / bg_frac } else { 1.0 }
  TfbsEnrichResult::{
    matrix_id: matrix_id, tf_name: tf_name,
    n_targets_hit: n_targets_hit, n_background_hit: n_background_hit,
    total_targets: n_targets_total, total_background: n_background_total,
    fold_enrichment: fold_enrichment, p_value: p_value, adj_p: p_value,
  }
}

/// BH-FDR 多重检验校正。
pub fn tfbs_adjust_pvalues(results : Array[TfbsEnrichResult]) -> Array[TfbsEnrichResult] {
  let n = results.length()
  if n == 0 {
    return []
  }
  // 按 p 值升序
  let sorted = results.copy()
  for i in 0..<n {
    let mut min_idx = i
    let j_start = i + 1
    for j in j_start..<n {
      if sorted[j].p_value < sorted[min_idx].p_value {
        min_idx = j
      }
    }
    if min_idx != i {
      let tmp = sorted[i]
      sorted[i] = sorted[min_idx]
      sorted[min_idx] = tmp
    }
  }
  // BH 校正并构建新结果
   let adjusted : Array[TfbsEnrichResult] = Array::new()
  for i in 0..<n {
    let adj = sorted[i].p_value * n.to_double() / (i + 1).to_double()
    let adj_p = if adj > 1.0 { 1.0 } else { adj }
    adjusted.push(TfbsEnrichResult::{
      matrix_id: sorted[i].matrix_id, tf_name: sorted[i].tf_name,
      n_targets_hit: sorted[i].n_targets_hit, n_background_hit: sorted[i].n_background_hit,
      total_targets: sorted[i].total_targets, total_background: sorted[i].total_background,
      fold_enrichment: sorted[i].fold_enrichment, p_value: sorted[i].p_value, adj_p: adj_p,
    })
  }
  // 单调化（从后往前）
  let mut k = n - 2
  while k >= 0 {
    let prev_adj = adjusted[k + 1].adj_p
    let curr_adj = adjusted[k].adj_p
    if curr_adj > prev_adj {
      let item = adjusted[k]
      adjusted[k] = TfbsEnrichResult::{
        matrix_id: item.matrix_id, tf_name: item.tf_name,
        n_targets_hit: item.n_targets_hit, n_background_hit: item.n_background_hit,
        total_targets: item.total_targets, total_background: item.total_background,
        fold_enrichment: item.fold_enrichment, p_value: item.p_value, adj_p: prev_adj,
      }
    }
    k = k - 1
  }
  adjusted
}

/// 内置示例 JASPAR 矩阵（多条）。
pub fn tfbs_example_jaspar() -> Array[TfbsJasparMatrix] {
  [
    TfbsJasparMatrix::new("MA0001.1", "MA0001.1", [
      [0.8, 0.1, 0.05, 0.05],
      [0.05, 0.8, 0.1, 0.05],
      [0.05, 0.1, 0.8, 0.05],
      [0.05, 0.05, 0.1, 0.8],
    ]),
    TfbsJasparMatrix::new("MA0002.1", "MA0002.1", [
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
      [0.25, 0.25, 0.25, 0.25],
    ]),
  ]
}
