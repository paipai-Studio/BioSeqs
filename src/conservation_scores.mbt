///|
/// conservation_scores.mbt — GERP/phastCons 风格进化约束打分
/// 计算端(本地确定性实现).
///
/// 诚实声明(与 GERP++/phastCons 的差距):
/// - GERP++ 以"中性树期望替换数 vs 观测似然最优树替换数"之差
///   (rejected substitutions) 计算, 需替换模型似然; 本模块用
///   熵代理: RS-lite = tree_length x (1 - H_obs/H_max), 保守列
///   (低熵)得高分, 可变列(最大熵)得 0 分;
/// - phastCons 为 phylo-HMM 完整似然 + 参数标定; 本模块为两状态
///   Viterbi: 发射对数似然比(conserved 中心模型 vs 均匀中性模型)
///   + 固定转移概率, 输出逐位保守后验与保守区段.
///
/// 输入: 等长 MSA 行数组(碱基 ACGT/UN; gap '-' 或 '?' 与 N
/// 同义, 不参与列统计); tree_length 为中性树总分支长(替换/位点,
/// GERP-lite 用).
///
/// 数据不出本机.

///|
/// 一列的碱基计数(A/C/G/T, gap/N 排除).
priv struct CsColumn {
  /// 有效(非 gap)序列数.
  n : Int32
  /// 多数碱基计数.
  majority : Int32
  counts : Array[Int32]
}

///|
fn cs_base_index(b : Char) -> Int32 {
  // 小写归一化
  let c = if b >= 'a' && b <= 'z' {
    (b.to_int() - 32).to_char().unwrap_or(b)
  } else {
    b
  }
  match c {
    'A' => 0
    'C' => 1
    'G' => 2
    'T' => 3
    _ => -1
  }
}

///|
fn cs_is_gap(b : Char) -> Bool {
  b == '-' || b == '?' || b == '.' || b == '~'
}

///|
/// 统计 MSA 第 col 列(0-based).
fn cs_column(msa : Array[String], col : Int32) -> CsColumn {
  let counts : Array[Int32] = Array::make(4, 0)
  let mut n = 0
  for row in msa {
    if col >= row.length() {
      continue
    }
    let b = row.unsafe_get(col).unsafe_to_char()
    if cs_is_gap(b) || b == 'N' || b == 'n' {
      continue
    }
    let idx = cs_base_index(b)
    if idx >= 0 {
      counts[idx] = counts[idx] + 1
      n = n + 1
    }
  }
  let mut majority = 0
  for c in counts {
    if c > majority {
      majority = c
    }
  }
  CsColumn::{ n, majority, counts, }
}

///|
/// 列熵(以有效碱基频率计, ln 域).
fn cs_entropy(col : CsColumn) -> Double {
  if col.n <= 0 {
    return 0.0
  }
  let mut h = 0.0
  for c in col.counts {
    if c > 0 {
      let p = c.to_double() / col.n.to_double()
      h = h - p * @math.ln(p)
    }
  }
  h
}

// ---------------------------------------------------------------------------
// GERP-lite
// ---------------------------------------------------------------------------

///|
/// GERP-lite 拒绝打分(RS-lite): 逐列
/// tree_length x (1 - H/H_max), H_max = ln 4.
/// 保守列(全同) -> tree_length; 完全可变列 -> 0;
/// gap 占比过高的列(有效数 < min_effective)返回 0.
pub fn conservation_gerp(
  msa : Array[String],
  tree_length~ : Double = 1.0,
  min_effective~ : Int32 = 2,
) -> Array[Double] {
  if msa.length() == 0 {
    return []
  }
  let width = msa[0].length()
  let hmax = @math.ln(4.0)
  let out : Array[Double] = Array(capacity=width)
  for col = 0; col < width; col = col + 1 {
    let c = cs_column(msa, col)
    if c.n < min_effective {
      out.push(0.0)
      continue
    }
    let h = cs_entropy(c)
    let constrained = if hmax > 0.0 { 1.0 - h / hmax } else { 0.0 }
    out.push(tree_length * constrained)
  }
  out
}

// ---------------------------------------------------------------------------
// phastCons-lite 两状态 HMM
// ---------------------------------------------------------------------------

///|
/// 逐列发射对数似然比 ln P(col|conserved) - ln P(col|neutral).
/// conserved 模型: 多数碱基概率 p_keep(默认 0.95), 其余均分;
/// neutral 模型: 均匀 0.25.
fn cs_emission_llr(col : CsColumn, p_keep~ : Double = 0.95) -> Double {
  if col.n <= 0 {
    return 0.0
  }
  let p_rest = (1.0 - p_keep) / 3.0
  // conserved: 多数碱基以 p_keep, 其余 p_rest
  let ll_cons = col.majority.to_double() * @math.ln(p_keep) +
    (col.n - col.majority).to_double() * @math.ln(p_rest)
  // neutral: 均匀
  let ll_neu = col.n.to_double() * @math.ln(0.25)
  ll_cons - ll_neu
}

///|
/// 逐位保守后验(0-1): 发射 LLR 的 logistic 压缩
/// posterior = 1 / (1 + exp(-llr / scale)).
pub fn conservation_posterior(
  msa : Array[String],
  scale~ : Double = 2.0,
) -> Array[Double] {
  if msa.length() == 0 {
    return []
  }
  let width = msa[0].length()
  let out : Array[Double] = Array(capacity=width)
  for col = 0; col < width; col = col + 1 {
    let c = cs_column(msa, col)
    let llr = cs_emission_llr(c)
    let z = llr / scale
    // 数值稳定 logistic
    let p = if z >= 0.0 {
      1.0 / (1.0 + @math.exp(-z))
    } else {
      let e = @math.exp(z)
      e / (1.0 + e)
    }
    out.push(p)
  }
  out
}

///|
/// 保守区段(0-based 半开).
pub(all) struct CsSegment {
  start : Int32
  end_ : Int32
  /// 区段平均后验.
  mean_posterior : Double
}

///|
/// 两状态 Viterbi 解码保守区段.
///
/// 状态: 0 = neutral, 1 = conserved; 转移(固定):
/// self-neutral 0.99, ->conserved 0.01;
/// self-conserved 0.95, ->neutral 0.05.
/// 发射: neutral 对数均匀 0(常数), conserved 为 LLR.
pub fn conservation_viterbi_segments(
  msa : Array[String],
  posterior~ : Array[Double]? = None,
) -> Array[CsSegment] {
  if msa.length() == 0 {
    return []
  }
  let post = match posterior {
    Some(p) => p
    None => conservation_posterior(msa)
  }
  let width = post.length()
  if width == 0 {
    return []
  }
  // 转移对数
  let ln_self_n = @math.ln(0.99)
  let ln_to_c = @math.ln(0.01)
  let ln_self_c = @math.ln(0.95)
  let ln_to_n = @math.ln(0.05)
  // 发射对数比(conserved 相对 neutral): 后验 -> LLR 还原
  // llr = ln(p/(1-p)) * scale(默认 2.0); 这里直接用 ln(p/(1-p))
  // 的单调变换做 Viterbi 权重(等效的状态偏好).
  let llr = fn(p : Double) -> Double {
    if p <= 0.0 {
      return -20.0
    }
    if p >= 1.0 {
      return 20.0
    }
    @math.ln(p / (1.0 - p))
  }
  // Viterbi DP
  let v_n : Array[Double] = Array::make(width, 0.0)
  let v_c : Array[Double] = Array::make(width, 0.0)
  let ptr : Array[Int32] = Array::make(width, 0) // 每位的前驱状态
  v_n[0] = 0.0
  v_c[0] = @math.ln(0.05) + llr(post[0])
  for i = 1; i < width; i = i + 1 {
    let e = llr(post[i])
    // neutral
    let from_n = v_n[i - 1] + ln_self_n
    let from_c = v_c[i - 1] + ln_to_n
    if from_n >= from_c {
      v_n[i] = from_n
    } else {
      v_n[i] = from_c
    }
    // conserved
    let from_n2 = v_n[i - 1] + ln_to_c + e
    let from_c2 = v_c[i - 1] + ln_self_c + e
    if from_n2 >= from_c2 {
      v_c[i] = from_n2
    } else {
      v_c[i] = from_c2
    }
    // 记录最优前驱(用于回溯)
    let best_n : Double = if v_n[i] >= v_c[i] { v_n[i] } else { v_c[i] }
    // 前驱判定: 状态 i 取 neutral 时来自哪个
    let pred_n_from_n = v_n[i - 1] + ln_self_n
    let pred_n_from_c = v_c[i - 1] + ln_to_n
    if best_n == v_n[i] {
      ptr[i] = if pred_n_from_n >= pred_n_from_c { 0 } else { 1 }
    } else {
      let pred_c_from_n = v_n[i - 1] + ln_to_c
      let pred_c_from_c = v_c[i - 1] + ln_self_c
      ptr[i] = if pred_c_from_n >= pred_c_from_c { 0 } else { 1 }
    }
  }
  // 回溯终点状态
  let states : Array[Int32] = Array::make(width, 0)
  states[width - 1] = if v_c[width - 1] > v_n[width - 1] { 1 } else { 0 }
  let mut i = width - 1
  while i > 0 {
    states[i - 1] = if states[i] == 0 {
      ptr[i]
    } else {
      // conserved 态的前驱按其转移重算
      let from_n = v_n[i - 1] + ln_to_c
      let from_c = v_c[i - 1] + ln_self_c
      if from_n >= from_c {
        0
      } else {
        1
      }
    }
    i = i - 1
  }
  // 提取保守区段
  let out : Array[CsSegment] = Array(capacity=0)
  let mut seg_start : Int32? = None
  for k = 0; k < width; k = k + 1 {
    if states[k] == 1 {
      match seg_start {
        None => seg_start = Some(k)
        Some(_) => ()
      }
    } else {
      match seg_start {
        Some(s) => {
          let sum = cs_mean_range(post, s, k)
          out.push(CsSegment::{ start: s, end_: k, mean_posterior: sum, })
          seg_start = None
        }
        None => ()
      }
    }
  }
  match seg_start {
    Some(s) => {
      let sum = cs_mean_range(post, s, width)
      out.push(CsSegment::{ start: s, end_: width, mean_posterior: sum, })
    }
    None => ()
  }
  out
}

///|
fn cs_mean_range(post : Array[Double], start : Int32, end_ : Int32) -> Double {
  if end_ <= start {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = start
  while i < end_ {
    sum = sum + post[i]
    i = i + 1
  }
  sum / (end_ - start).to_double()
}

// ---------------------------------------------------------------------------
// GenomicScores 闭环与导出
// ---------------------------------------------------------------------------

///|
/// 逐列打分 -> GenomicScores 容器(GERP-lite 人口).
pub fn conservation_to_gscores(
  gs : GScores,
  chrom : String,
  scores : Array[Double],
) -> Unit raise SeqIOError {
  let positions : Array[Int32] = Array(capacity=scores.length())
  for i = 0; i < scores.length(); i = i + 1 {
    positions.push(i + 1)
  }
  gs.set_scores(chrom, "gerp_lite", positions, scores)
}

///|
/// 打分 -> TSV(chrom pos gerp posterior).
pub fn conservation_tsv(
  chrom : String,
  gerp : Array[Double],
  posterior : Array[Double],
) -> String {
  let sb = StringBuilder()
  sb.write_string("chrom\tpos\tgerp_lite\tposterior\n")
  for i = 0; i < gerp.length(); i = i + 1 {
    let p = if i < posterior.length() { posterior[i].to_string() } else { "." }
    sb.write_string(chrom)
    sb.write_char('\t')
    sb.write_string((i + 1).to_string())
    sb.write_char('\t')
    sb.write_string(gerp[i].to_string())
    sb.write_char('\t')
    sb.write_string(p)
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// 打分与区段 -> JSON 值.
pub fn conservation_value(
  gerp : Array[Double],
  posterior : Array[Double],
  segments : Array[CsSegment],
) -> JsonValue {
  let g : Array[JsonValue] = Array(capacity=gerp.length())
  for v in gerp {
    g.push(JsonNum(v))
  }
  let p : Array[JsonValue] = Array(capacity=posterior.length())
  for v in posterior {
    p.push(JsonNum(v))
  }
  let segs : Array[JsonValue] = Array(capacity=segments.length())
  for s in segments {
    segs.push(
      JsonObj(
        Map([
          ("start", JsonNum(s.start.to_double())),
          ("end", JsonNum(s.end_.to_double())),
          ("mean_posterior", JsonNum(s.mean_posterior)),
        ]),
      ),
    )
  }
  JsonObj(
    Map([
      ("gerp_lite", JsonArr(g)),
      ("posterior", JsonArr(p)),
      ("conserved_segments", JsonArr(segs)),
    ]),
  )
}

///|
/// JSON 文本.
pub fn conservation_json(
  gerp : Array[Double],
  posterior : Array[Double],
  segments : Array[CsSegment],
) -> String {
  json_to_string(conservation_value(gerp, posterior, segments))
}

///|
/// 摘要行.
pub fn conservation_summary_line(
  gerp : Array[Double],
  segments : Array[CsSegment],
) -> String {
  let mut max_g = 0.0
  let mut sum_g = 0.0
  for v in gerp {
    sum_g = sum_g + v
    if v > max_g {
      max_g = v
    }
  }
  let mean = if gerp.length() > 0 {
    sum_g / gerp.length().to_double()
  } else {
    0.0
  }
  let mut seg_bp = 0
  for s in segments {
    seg_bp = seg_bp + (s.end_ - s.start)
  }
  "conservation: cols=\{gerp.length()}, mean_gerp=\{mean}, max_gerp=\{max_g}, conserved_segments=\{segments.length()}, conserved_bp=\{seg_bp}"
}
