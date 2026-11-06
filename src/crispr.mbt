///|
/// crispr.mbt — CRISPR gRNA 设计器(PAM 扫描 + on-target 打分 + off-target 评估).
///
/// 对齐 CRISPOR / CHOPCHOP 的核心设计语义(本地筛查级):
/// - 候选生成: SpCas9(NGG) 20nt 间隔区扫描, 双链(反向互补侧 NGG -> CCN)
/// - on-target 打分: Doench 2016(Rule Set 2) 简化特征
///   (位置依赖碱基系数 + GC 区间 + 均聚物/TTTT 终止信号惩罚)
/// - off-target 评估: FM-index 近似搜索(复用 bwa.mbt), 统计全基因组
///   ≤3 错配位点数, MIT 风格位点得分(mismatch 位置加权), 聚合为
///   guide 特异性分
/// - 输出 JSON(衔接 Web/报告架构), 候选按总分降序
///
/// 内置 SpCas9(NGG)/SpG(SGG 缩短版兼容)/Cas12a(TTTV) 三 PAM 语义.
/// 数据不出本机.
pub(all) struct CrisprPam {
  /// 3' PAM 正则式语义: 用 IUPAC 简化集合表达, 如 ["N","G","G"]
  pattern : Array[String]
  /// PAM 位于间隔区 3'(Cas9) 或 5'(Cas12a) 端
  prime : String
}

///|
/// SpCas9 NGG PAM(20nt spacer + NGG).
pub fn crispr_pam_spcas9() -> CrisprPam {
  { pattern: ["N", "G", "G"], prime: "3", }
}

///|
/// SpG SGG PAM(缩短版, 兼容 NGG 语义的超集).
pub fn crispr_pam_spg() -> CrisprPam {
  { pattern: ["S", "G", "G"], prime: "3", }
}

///|
/// Cas12a TTTV PAM(PAM 在 5' 端, 23nt spacer).
pub fn crispr_pam_cas12a() -> CrisprPam {
  { pattern: ["T", "T", "T", "V"], prime: "5", }
}

///|
/// 一条 gRNA 候选.
pub(all) struct CrisprGuide {
  spacer : String
  pam : String
  strand : Int32
  /// 间隔区起点(正向链坐标, 0-based; 反向链候选已换算)
  start : Int32
  end : Int32
  gc_percent : Double
  mut on_target : Double
  mut specificity : Double
  mut total_score : Double
  mut off_targets : Array[CrisprOffTarget]
  mut warnings : Array[String]
}

///|
/// 一个 off-target 位点.
pub(all) struct CrisprOffTarget {
  pos : Int32
  strand : Int32
  mismatches : Int32
  site_score : Double
}

///|
/// PAM 单字符匹配(IUPAC 简化集合).
pub fn crispr_pam_char_matches(pattern_char : String, base : Char) -> Bool {
  match pattern_char {
    "N" => true
    "A" => base == 'A'
    "C" => base == 'C'
    "G" => base == 'G'
    "T" => base == 'T'
    "S" => base == 'G' || base == 'C'
    "W" => base == 'A' || base == 'T'
    "R" => base == 'A' || base == 'G'
    "Y" => base == 'C' || base == 'T'
    "V" => base != 'T'
    "B" => base != 'A'
    "H" => base != 'G'
    "D" => base != 'C'
    "K" => base == 'G' || base == 'T'
    "M" => base == 'A' || base == 'C'
    _ => false
  }
}

///|
/// 判断序列子串是否匹配 PAM.
fn crispr_pam_matches(pam : CrisprPam, seq : String, offset : Int32) -> Bool {
  let n = pam.pattern.length()
  if offset < 0 || offset + n > seq.length() {
    return false
  }
  for i = 0; i < n; i = i + 1 {
    let base = seq.unsafe_get(offset + i).unsafe_to_char()
    if !crispr_pam_char_matches(pam.pattern[i], base) {
      return false
    }
  }
  true
}

///|
/// Doench 2016 简化位置系数(20nt spacer, 0=P5 端).
/// 3' 端(PAM 近端)错配容忍度低 -> 系数高; 数据为筛查级近似.
fn crispr_position_weight(i : Int32) -> Double {
  // PAM 近端(位置 17-19)最关键
  if i >= 17 {
    1.0
  } else if i >= 13 {
    0.75
  } else if i >= 8 {
    0.5
  } else if i >= 4 {
    0.3
  } else {
    0.2
  }
}

///|
/// on-target 效力打分(Doench 2016 简化特征, [0,1]).
/// 特征: GC 区间(40-65% 最优) + TTTT 转录终止 + 均聚物 + 起始 G.
pub fn crispr_on_target_score(spacer : String) -> Double {
  let n = spacer.length()
  if n == 0 {
    return 0.0
  }
  let mut score = 0.4
  // GC 区间
  let mut gc = 0
  for i = 0; i < n; i = i + 1 {
    let c = spacer.unsafe_get(i).unsafe_to_char()
    if c == 'G' || c == 'C' {
      gc = gc + 1
    }
  }
  let gc_pct = gc.to_double() / n.to_double() * 100.0
  if gc_pct >= 40.0 && gc_pct <= 65.0 {
    score = score + 0.25
  } else {
    let dev = if gc_pct < 40.0 { 40.0 - gc_pct } else { gc_pct - 65.0 }
    let penalty = if dev > 20.0 { 0.2 } else { dev / 100.0 }
    score = score + 0.25 - penalty
  }
  // TTTT 转录终止信号(U6 启动子下致命)
  if spacer.contains("TTTT") {
    score = score - 0.25
  }
  // 均聚物(>=4 连续相同)
  let mut run = 1
  let mut max_run = 1
  for i = 1; i < n; i = i + 1 {
    if spacer.unsafe_get(i) == spacer.unsafe_get(i - 1) {
      run = run + 1
      if run > max_run {
        max_run = run
      }
    } else {
      run = 1
    }
  }
  if max_run >= 4 {
    score = score - 0.1
  }
  // 起始 G(U6 转录偏好)
  if spacer.unsafe_get(0).unsafe_to_char() == 'G' {
    score = score + 0.1
  }
  if score < 0.0 {
    score = 0.0
  }
  if score > 1.0 {
    score = 1.0
  }
  score
}

///|
/// MIT 风格 off-target 位点得分([0,1], 1=完全匹配).
/// mismatch 越多且越靠 PAM 近端(3' 端), 得分越低.
pub fn crispr_off_target_site_score(
  spacer : String,
  off_target_seq : String,
) -> Double {
  let n = spacer.length()
  if n == 0 || off_target_seq.length() < n {
    return 0.0
  }
  let mut mismatches = 0
  let mut penalty = 1.0
  for i = 0; i < n; i = i + 1 {
    if spacer.unsafe_get(i) != off_target_seq.unsafe_get(i) {
      mismatches = mismatches + 1
      let w = crispr_position_weight(i)
      penalty = penalty * (1.0 - 0.6 * w)
    }
  }
  if mismatches == 0 {
    return 1.0
  }
  // 平均错配距离因子
  let avg_dist = 1.0 / mismatches.to_double()
  penalty * avg_dist
}

///|
/// 扫描候选 gRNA(双链).
/// 返回候选(spacer/PAM/strand/坐标), 未打分.
pub fn crispr_scan_guides(
  sequence : String,
  pam~ : CrisprPam = crispr_pam_spcas9(),
  spacer_len~ : Int32 = 20,
) -> Array[CrisprGuide] {
  let n = sequence.length()
  let pam_len = pam.pattern.length()
  let out : Array[CrisprGuide] = Array(capacity=0)
  let is_cas9 = pam.prime == "3"
  for pos = 0; pos + pam_len <= n; pos = pos + 1 {
    if crispr_pam_matches(pam, sequence, pos) {
      let pam_seq = sequence[pos:pos + pam_len].to_owned()
      if is_cas9 {
        // 正向链: [spacer_start, pos) + PAM
        let s_start = pos - spacer_len
        if s_start >= 0 {
          let spacer = sequence[s_start:pos].to_owned()
          push_guide(out, spacer, pam_seq, 0, s_start, pos)
        }
      } else {
        // Cas12a: PAM + [pos+pam_len, pos+pam_len+spacer_len)
        let s_end = pos + pam_len + spacer_len
        if s_end <= n {
          let spacer = sequence[pos + pam_len:s_end].to_owned()
          push_guide(out, spacer, pam_seq, 0, pos + pam_len, s_end)
        }
      }
    }
    // 反向链: CCN(NGG 反向互补) / Cas12a 反向
    let rc_pam : Array[String] = []
    ignore(rc_pam)
    // 反向互补 PAM 序列逐位检测: 原5' 端出现互补模式
    if crispr_pam_matches_reverse(pam, sequence, pos) {
      let pam_seq = sequence[pos:pos + pam_len].to_owned()
      let rc = reverse_complement(pam_seq)
      if is_cas9 {
        // 反向链: PAM(原链 CCN) + spacer(其反向互补)
        let s_end = pos + pam_len + spacer_len
        if s_end <= n {
          let spacer_rc = sequence[pos + pam_len:s_end].to_owned()
          let spacer = reverse_complement(spacer_rc)
          push_guide(out, spacer, rc, 1, pos + pam_len, s_end)
        }
      } else {
        // Cas12a 反向: spacer 在 PAM 5' 侧
        let s_start = pos - spacer_len
        if s_start >= 0 {
          let spacer_rc = sequence[s_start:pos].to_owned()
          let spacer = reverse_complement(spacer_rc)
          push_guide(out, spacer, rc, 1, s_start, pos)
        }
      }
    }
  }
  out
}

///|
fn push_guide(
  out : Array[CrisprGuide],
  spacer : String,
  pam_seq : String,
  strand : Int32,
  start : Int32,
  end : Int32,
) -> Unit {
  let n = spacer.length()
  let mut gc = 0
  for i = 0; i < n; i = i + 1 {
    let c = spacer.unsafe_get(i).unsafe_to_char()
    if c == 'G' || c == 'C' {
      gc = gc + 1
    }
  }
  let gc_pct = if n > 0 { gc.to_double() / n.to_double() * 100.0 } else { 0.0 }
  out.push({
    spacer,
    pam: pam_seq,
    strand,
    start,
    end,
    gc_percent: gc_pct,
    on_target: 0.0,
    specificity: 1.0,
    total_score: 0.0,
    off_targets: [],
    warnings: [],
  })
}

///|
/// 反向互补语义的 PAM 匹配: 原链上出现 PAM 的反向互补模式.
/// (NGG 的反向互补为 CCN, 检测原链 pos..pos+len 是否为 CCN)
pub fn crispr_pam_matches_reverse(
  pam : CrisprPam,
  seq : String,
  pos : Int32,
) -> Bool {
  let n = pam.pattern.length()
  if pos < 0 || pos + n > seq.length() {
    return false
  }
  // PAM 反向互补后逐位匹配(从原链 pos 起正向读)
  for i = 0; i < n; i = i + 1 {
    // pam 反向互补的第 i 位 = complement(pam.pattern[n-1-i])
    let pc = pam.pattern[n - 1 - i]
    let base = seq.unsafe_get(pos + i).unsafe_to_char()
    let comp_ok = match pc {
      "N" => true
      "A" => base == 'T'
      "C" => base == 'G'
      "G" => base == 'C'
      "T" => base == 'A'
      "S" => base == 'C' || base == 'G'
      "W" => base == 'A' || base == 'T'
      "R" => base == 'C' || base == 'T'
      "Y" => base == 'A' || base == 'G'
      "V" => base == 'A'
      "B" => base == 'T'
      "H" => base == 'G'
      "D" => base == 'C'
      "K" => base == 'C' || base == 'A'
      "M" => base == 'T' || base == 'G'
      _ => false
    }
    if !comp_ok {
      return false
    }
  }
  true
}

///|
/// off-target 评估: FM-index 近似搜索(≤max_mismatches), 排除原位点.
/// 返回按位点得分降序的 off-target 列表与特异性分(0=多脱靶, 1=特异).
pub fn crispr_eval_off_targets(
  guide : CrisprGuide,
  reference : String,
  max_mismatches~ : Int32 = 3,
) -> (Array[CrisprOffTarget], Double) {
  let index = bwa_index_build("ref", reference)
  let hits = fm_inexact_search(index.fm, guide.spacer, max_mismatches)
  let offs : Array[CrisprOffTarget] = Array(capacity=0)
  let n = reference.length()
  let spacer_len = guide.spacer.length()
  let seen : Map[Int32, Bool] = Map([], capacity=0)
  for hit in hits {
    if hit.pos < 0 || hit.pos + spacer_len > n {
      continue
    }
    if seen.get(hit.pos) is Some(_) {
      continue
    }
    seen.set(hit.pos, true)
    // CRISPR gRNA-DNA 配对为 Hamming 型: 逐位比较重算错配数,
    // 过滤 FM 回溯的 indel 伪解释
    let site_seq = reference[hit.pos:hit.pos + spacer_len].to_owned()
    let mut mismatches = 0
    for i = 0; i < spacer_len; i = i + 1 {
      if guide.spacer.unsafe_get(i) != site_seq.unsafe_get(i) {
        mismatches = mismatches + 1
      }
    }
    if mismatches > max_mismatches {
      continue
    }
    // 原位点排除: 完全匹配且位置一致(正向)或为反向链原区段
    if mismatches == 0 {
      let is_original = if guide.strand == 0 {
        hit.pos == guide.start
      } else {
        // 反向链 spacer 对应原链 [start, end) 区段
        hit.pos == guide.start
      }
      if is_original {
        continue
      }
    }
    let score = crispr_off_target_site_score(guide.spacer, site_seq)
    offs.push({ pos: hit.pos, strand: 0, mismatches, site_score: score, })
  }
  offs.sort_by(fn(a, b) {
    if a.site_score > b.site_score {
      -1
    } else if a.site_score < b.site_score {
      1
    } else {
      compare(a.pos, b.pos)
    }
  })
  // 特异性: 1 - Σ site_score(截断到 0)
  let mut total_site = 0.0
  for o in offs {
    total_site = total_site + o.site_score
  }
  let specificity = if total_site >= 1.0 { 0.0 } else { 1.0 - total_site }
  (offs, specificity)
}

///|
/// 完整设计流程: 扫描 -> on-target 打分 -> off-target 评估 -> 排序.
pub fn crispr_design_guides(
  sequence : String,
  reference~ : String = "",
  pam~ : CrisprPam = crispr_pam_spcas9(),
  spacer_len~ : Int32 = 20,
  max_mismatches~ : Int32 = 3,
  top_n~ : Int32 = 20,
) -> Array[CrisprGuide] {
  let guides = crispr_scan_guides(sequence, pam~, spacer_len~)
  for guide in guides {
    guide.on_target = crispr_on_target_score(guide.spacer)
    if reference.length() >= guide.spacer.length() {
      let (offs, spec) = crispr_eval_off_targets(
        guide,
        reference,
        max_mismatches~,
      )
      guide.off_targets = offs
      guide.specificity = spec
    } else {
      guide.specificity = 1.0
    }
    // 警告
    if guide.spacer.contains("TTTT") {
      guide.warnings.push("TTTT: Pol III terminator")
    }
    if guide.gc_percent < 20.0 || guide.gc_percent > 80.0 {
      guide.warnings.push("extreme GC content")
    }
    if guide.off_targets.length() > 10 {
      guide.warnings.push("high off-target count")
    }
    // 总分: 0.5 * on-target + 0.5 * specificity
    guide.total_score = 0.5 * guide.on_target + 0.5 * guide.specificity
  }
  guides.sort_by(fn(a, b) {
    if a.total_score > b.total_score {
      -1
    } else if a.total_score < b.total_score {
      1
    } else {
      compare(a.start, b.start)
    }
  })
  if guides.length() > top_n {
    guides.truncate(top_n)
  }
  guides
}

///|
/// 设计摘要行.
pub fn crispr_guides_summary(guides : Array[CrisprGuide]) -> String {
  let mut n_off = 0
  for g in guides {
    n_off = n_off + g.off_targets.length()
  }
  "crispr: \{guides.length()} guides, \{n_off} total off-target sites"
}

///|
/// 候选 JSON 值(衔接报告/Web 架构).
pub fn CrisprGuide::to_json_value(self : CrisprGuide) -> JsonValue {
  let offs : Array[JsonValue] = Array(capacity=0)
  for o in self.off_targets {
    offs.push(
      JsonObj(
        Map([
          ("pos", JsonNum(o.pos.to_double())),
          ("strand", JsonNum(o.strand.to_double())),
          ("mismatches", JsonNum(o.mismatches.to_double())),
          ("site_score", JsonNum(o.site_score)),
        ]),
      ),
    )
  }
  let warns : Array[JsonValue] = Array(capacity=0)
  for w in self.warnings {
    warns.push(JsonStr(w))
  }
  JsonObj(
    Map([
      ("spacer", JsonStr(self.spacer)),
      ("pam", JsonStr(self.pam)),
      ("strand", JsonNum(self.strand.to_double())),
      ("start", JsonNum(self.start.to_double())),
      ("end", JsonNum(self.end.to_double())),
      ("gc_percent", JsonNum(self.gc_percent)),
      ("on_target", JsonNum(self.on_target)),
      ("specificity", JsonNum(self.specificity)),
      ("total_score", JsonNum(self.total_score)),
      ("off_targets", JsonArr(offs)),
      ("warnings", JsonArr(warns)),
    ]),
  )
}

///|
/// 设计结果 JSON 文本.
pub fn crispr_guides_json(guides : Array[CrisprGuide]) -> String {
  let arr : Array[JsonValue] = Array(capacity=0)
  for g in guides {
    arr.push(g.to_json_value())
  }
  json_to_string(JsonArr(arr))
}
