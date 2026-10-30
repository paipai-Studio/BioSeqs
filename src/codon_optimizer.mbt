///|
/// codon_optimizer.mbt — 蛋白 -> DNA 密码子优化(宿主偏好 + 序列约束).
///
/// 对齐 DNAchisel / codon-optimization 工具的核心语义(本地确定性实现):
/// - 按 NCBI 密码子表反向翻译蛋白序列
/// - 候选密码子按宿主 CAI 权重加权, 用确定性 RNG 打破平局
/// - 约束: GC 区间、禁止 motif(双链检查)、均聚物长度、可选规避 CpG
/// - 生成结果可做 GC/CAI/违规项报告, 数据不出本机

///|
/// 密码子优化错误.
pub suberror CodonOptimizerError {
  CodonOptimizerError(String)
}

///|
/// IUPAC 互补表.
let codon_opt_complement_map : Map[String, String] = Map(
  [
    ("A", "T"),
    ("T", "A"),
    ("C", "G"),
    ("G", "C"),
    ("N", "N"),
    ("R", "Y"),
    ("Y", "R"),
    ("S", "S"),
    ("W", "W"),
    ("K", "M"),
    ("M", "K"),
    ("B", "V"),
    ("V", "B"),
    ("D", "H"),
    ("H", "D"),
  ],
  capacity=15,
)

///|
/// IUPAC 简并碱基集合.
fn codon_opt_iupac_set(code : Int32) -> String {
  match code {
    65 => "A"
    67 => "C"
    71 => "G"
    84 => "T"
    85 => "T"
    78 => "ACGT"
    82 => "AG"
    89 => "CT"
    83 => "CG"
    87 => "AT"
    75 => "GT"
    77 => "AC"
    66 => "CGT"
    68 => "AGT"
    72 => "ACT"
    86 => "ACG"
    _ => ""
  }
}

///|
/// 单个碱基是否匹配 IUPAC 简并码.
pub fn codon_optimizer_base_matches(pattern : Char, base : Char) -> Bool {
  let set = codon_opt_iupac_set(pattern.to_int())
  if set.length() == 0 {
    return false
  }
  set.contains_char(base.to_ascii_uppercase())
}

///|
/// 在 `dna` 的 `start` 处是否匹配 IUPAC `motif`.
pub fn codon_optimizer_motif_matches(
  dna : String,
  start : Int32,
  motif : String,
) -> Bool {
  if start < 0 || start + motif.length() > dna.length() {
    return false
  }
  for i = 0; i < motif.length(); i = i + 1 {
    if !codon_optimizer_base_matches(
        motif.unsafe_get(i).unsafe_to_char(),
        dna.unsafe_get(start + i).unsafe_to_char(),
      ) {
      return false
    }
  }
  true
}

///|
/// 反向互补(IUPAC).
pub fn codon_optimizer_revcomp(seq : String) -> String {
  let sb = StringBuilder()
  for i = seq.length() - 1; i >= 0; i = i - 1 {
    let c = seq.unsafe_get(i).unsafe_to_char().to_ascii_uppercase()
    sb.write_string(codon_opt_complement_map.get(c.to_string()).unwrap_or("N"))
  }
  sb.to_string()
}

///|
/// 从限制酶名称构造禁止 motif 列表(位点 + 反向互补去重).
pub fn codon_optimizer_sites_for_enzymes(
  enzymes : Array[String],
) -> Array[String] {
  let out : Array[String] = Array(capacity=enzymes.length() * 2)
  for name in enzymes {
    match get_bio_restriction_enzyme(name) {
      Some(enzyme) => {
        let site = enzyme.site.to_upper()
        if site.length() > 0 && !out.contains(site) {
          out.push(site)
        }
        let rc = codon_optimizer_revcomp(site)
        if !out.contains(rc) {
          out.push(rc)
        }
      }
      None => ()
    }
  }
  out
}

///|
/// 一次优化使用的配置.
pub(all) struct CodonOptimizerConfig {
  table : Int32
  gc_min : Double
  gc_max : Double
  max_run : Int32
  avoid_cpg : Bool
  forbidden_motifs : Array[String]
  seed : UInt64
}

///|
/// 构造配置并校验参数.
pub fn codon_optimizer_config(
  table~ : Int32 = 1,
  gc_min~ : Double = 0.0,
  gc_max~ : Double = 1.0,
  max_run~ : Int32 = 0,
  avoid_cpg~ : Bool = false,
  forbidden_motifs~ : Array[String] = [],
  seed~ : UInt64 = 42UL,
) -> CodonOptimizerConfig raise CodonOptimizerError {
  if !codon_table_is_valid(table) {
    raise CodonOptimizerError("invalid NCBI codon table id \{table}")
  }
  if gc_min < 0.0 || gc_max > 1.0 || gc_min > gc_max {
    raise CodonOptimizerError("invalid GC range [\{gc_min}, \{gc_max}]")
  }
  {
    table,
    gc_min,
    gc_max,
    max_run,
    avoid_cpg,
    forbidden_motifs: forbidden_motifs.copy(),
    seed,
  }
}

///|
/// 序列中是否出现禁止 motif(正反链).
fn codon_opt_has_forbidden(dna : String, motifs : Array[String]) -> Bool {
  for motif in motifs {
    let m = motif.to_upper()
    if m.length() == 0 {
      continue
    }
    let rc = codon_optimizer_revcomp(m)
    for i = 0; i + m.length() <= dna.length(); i = i + 1 {
      if codon_optimizer_motif_matches(dna, i, m) ||
        codon_optimizer_motif_matches(dna, i, rc) {
        return true
      }
    }
  }
  false
}

///|
/// 统计禁止 motif 出现次数(正反链).
fn codon_opt_count_forbidden(dna : String, motifs : Array[String]) -> Int32 {
  let mut count = 0
  for motif in motifs {
    let m = motif.to_upper()
    if m.length() == 0 {
      continue
    }
    let rc = codon_optimizer_revcomp(m)
    for i = 0; i + m.length() <= dna.length(); i = i + 1 {
      if codon_optimizer_motif_matches(dna, i, m) ||
        codon_optimizer_motif_matches(dna, i, rc) {
        count = count + 1
      }
    }
  }
  count
}

///|
/// 最长均聚物长度.
fn codon_opt_max_homopolymer(dna : String) -> Int32 {
  if dna.length() == 0 {
    return 0
  }
  let mut best = 1
  let mut run = 1
  for i = 1; i < dna.length(); i = i + 1 {
    if dna.unsafe_get(i) == dna.unsafe_get(i - 1) {
      run = run + 1
      if run > best {
        best = run
      }
    } else {
      run = 1
    }
  }
  best
}

///|
/// GC 比例.
fn codon_opt_gc(dna : String) -> Double {
  if dna.length() == 0 {
    return 0.0
  }
  let mut gc = 0
  for i = 0; i < dna.length(); i = i + 1 {
    match dna.unsafe_get(i).to_int() {
      71 | 67 | 103 | 99 => gc = gc + 1
      _ => ()
    }
  }
  gc.to_double() / dna.length().to_double()
}

///|
/// 密码子是否包含 CpG.
fn codon_opt_has_cpg(codon : String) -> Bool {
  codon.to_upper().contains("CG")
}

///|
/// 在已生成序列的末端拼接新密码子后, 是否仍满足局部约束.
fn codon_opt_local_ok(
  prefix : String,
  codon : String,
  config : CodonOptimizerConfig,
) -> Bool {
  let trial = prefix + codon
  let join_start = prefix.length()
  for motif in config.forbidden_motifs {
    let m = motif.to_upper()
    if m.length() == 0 {
      continue
    }
    let rc = codon_optimizer_revcomp(m)
    let mut start = join_start - m.length() + 1
    if start < 0 {
      start = 0
    }
    while start + m.length() <= trial.length() {
      if codon_optimizer_motif_matches(trial, start, m) ||
        codon_optimizer_motif_matches(trial, start, rc) {
        return false
      }
      start = start + 1
    }
  }
  if config.max_run > 0 {
    let keep = config.max_run + codon.length() + 1
    let from = if prefix.length() > keep { prefix.length() - keep } else { 0 }
    if codon_opt_max_homopolymer(trial[from:].to_owned()) > config.max_run {
      return false
    }
  }
  if config.avoid_cpg {
    if codon_opt_has_cpg(codon) {
      return false
    }
    if prefix.length() > 0 {
      let last = prefix
        .unsafe_get(prefix.length() - 1)
        .unsafe_to_char()
        .to_ascii_uppercase()
      let first = codon.unsafe_get(0).unsafe_to_char().to_ascii_uppercase()
      if last == 'C' && first == 'G' {
        return false
      }
    }
  }
  true
}

///|
/// 按权重从候选密码子中选择一个.
fn codon_opt_select(
  candidates : Array[String],
  weights : Map[String, Double],
  rng : Rng,
) -> String {
  let mut total = 0.0
  let scores : Array[Double] = Array(capacity=candidates.length())
  for c in candidates {
    let mut w = weights.get(c).unwrap_or(0.0)
    if w <= 0.0 {
      w = 0.000001
    }
    scores.push(w)
    total = total + w
  }
  let mut r = rng.next_double() * total
  for i = 0; i < candidates.length(); i = i + 1 {
    r = r - scores[i]
    if r <= 0.0 {
      return candidates[i]
    }
  }
  candidates[candidates.length() - 1]
}

///|
/// DP 状态需要保留的末端长度(足够检测跨密码子 motif / 均聚物 / CpG).
fn codon_opt_suffix_len(config : CodonOptimizerConfig) -> Int32 {
  let mut len = 0
  if config.max_run > 0 {
    // 跨界均聚物 run 达 max_run+1 时, 左侧最多贡献 max_run 个碱基.
    len = config.max_run
  }
  for motif in config.forbidden_motifs {
    if motif.length() - 1 > len {
      len = motif.length() - 1
    }
  }
  if config.avoid_cpg && len < 1 {
    // 跨密码子 CpG(前缀末 C + 新密码子首 G)需要保留 1 个末端碱基.
    len = 1
  }
  len
}

///|
/// 取字符串末端 `len` 个碱基.
fn codon_opt_tail(s : String, len : Int32) -> String {
  if len <= 0 || s.length() == 0 {
    return ""
  }
  let from = if s.length() > len { s.length() - len } else { 0 }
  s[from:].to_owned()
}

///|
/// 反向翻译并做约束修复.
pub fn codon_optimize(
  protein : String,
  weights : Map[String, Double],
  config~ : CodonOptimizerConfig,
) -> String raise CodonOptimizerError {
  if protein.length() == 0 {
    raise CodonOptimizerError("protein sequence is empty")
  }
  let reverse = codon_table_reverse(table_id=config.table)
  let stops = codon_table_stop_codons(config.table)
  let suffix_len = codon_opt_suffix_len(config)
  let n = protein.length()

  // 逐位反向翻译 + 稀疏 DP: state = 已生成序列的末端 suffix,
  // 目标是最小化 -ln(宿主权重), 同时保证 motif/均聚物/CpG 局部合法.
  let mut states : Map[String, (Double, String, String)] = Map([], capacity=1)
  states[""] = (0.0, "", "")
  let history : Array[Map[String, (Double, String, String)]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let aa = protein.unsafe_get(i).unsafe_to_char().to_ascii_uppercase()
    let candidates = if aa == '*' {
      [if stops.length() > 0 { stops[0] } else { "TAA" }]
    } else if aa == 'X' {
      ["NNN"]
    } else {
      let list = reverse.get(aa.to_string()).unwrap_or([])
      if list.length() == 0 {
        raise CodonOptimizerError("unsupported amino acid `\{aa.to_string()}`")
      }
      list
    }
    let cur : Map[String, (Double, String, String)] = Map(
      [],
      capacity=states.length() * 2 + 1,
    )
    for st, info in states {
      for cand in candidates {
        let trial = st + cand
        if codon_opt_has_forbidden(trial, config.forbidden_motifs) ||
          (
            config.max_run > 0 &&
            codon_opt_max_homopolymer(trial) > config.max_run
          ) ||
          (config.avoid_cpg && codon_opt_has_cpg(trial)) {
          continue
        }
        let w = weights.get(cand).unwrap_or(0.0)
        let penalty = info.0 - (if w > 0.0 { @math.ln(w) } else { -13.8 })
        let ns = codon_opt_tail(trial, suffix_len)
        match cur.get(ns) {
          Some(old) => if penalty < old.0 { cur[ns] = (penalty, st, cand) }
          None => cur[ns] = (penalty, st, cand)
        }
      }
    }
    if cur.length() == 0 {
      // 约束不可满足时退化为最高权重候选, 避免整体失败.
      let mut best_state = ""
      let mut best_penalty = 1.0e300
      for st, info in states {
        if info.0 < best_penalty {
          best_penalty = info.0
          best_state = st
        }
      }
      let mut best_codon = candidates[0]
      let mut best_w = -1.0
      for cand in candidates {
        let w = weights.get(cand).unwrap_or(0.0)
        if w > best_w {
          best_w = w
          best_codon = cand
        }
      }
      cur[codon_opt_tail(best_state + best_codon, suffix_len)] = (
        best_penalty - (if best_w > 0.0 { @math.ln(best_w) } else { -13.8 }),
        best_state,
        best_codon,
      )
    }
    history.push(cur)
    states = cur
  }

  let mut st = ""
  let mut best_penalty = 1.0e300
  for s, info in states {
    if info.0 < best_penalty {
      best_penalty = info.0
      st = s
    }
  }
  let chosen : Array[String] = Array(capacity=n)
  for i = n - 1; i >= 0; i = i - 1 {
    match history[i].get(st) {
      Some(info) => {
        chosen.push(info.2)
        st = info.1
      }
      None => break
    }
  }
  let codons : Array[String] = Array(capacity=chosen.length())
  for i = chosen.length() - 1; i >= 0; i = i - 1 {
    codons.push(chosen[i])
  }
  let mut dna = codons.join("")
  // 全序列 GC 约束的局部修复: 逐次替换同义密码子朝目标区间移动.
  let mut repair_step = 0
  let max_steps = n * 4
  while repair_step < max_steps {
    let gc = codon_opt_gc(dna)
    let need_up = gc < config.gc_min
    let need_down = gc > config.gc_max
    let has_violation = codon_opt_has_forbidden(dna, config.forbidden_motifs) ||
      (config.max_run > 0 && codon_opt_max_homopolymer(dna) > config.max_run) ||
      (config.avoid_cpg && codon_opt_has_cpg(dna))
    if !need_up && !need_down && !has_violation {
      break
    }
    let mut changed = false
    for ci = 0; ci < codons.length(); ci = ci + 1 {
      let aa = protein
        .unsafe_get(ci)
        .unsafe_to_char()
        .to_ascii_uppercase()
        .to_string()
      let candidates = reverse.get(aa).unwrap_or([])
      for cand in candidates {
        if cand == codons[ci] {
          continue
        }
        let trial_codons = codons.copy()
        trial_codons[ci] = cand
        let trial = trial_codons.join("")
        let trial_gc = codon_opt_gc(trial)
        let gc_ok = if need_up {
          trial_gc > gc
        } else if need_down {
          trial_gc < gc
        } else {
          true
        }
        let viol_ok = !codon_opt_has_forbidden(trial, config.forbidden_motifs) &&
          !(config.max_run > 0 &&
          codon_opt_max_homopolymer(trial) > config.max_run) &&
          !(config.avoid_cpg && codon_opt_has_cpg(trial))
        if gc_ok && viol_ok {
          codons[ci] = cand
          dna = trial
          changed = true
          break
        }
      }
      if changed {
        break
      }
    }
    if !changed {
      break
    }
    repair_step = repair_step + 1
  }
  dna
}

///|
/// 优化结果报告.
pub(all) struct CodonOptimizationReport {
  dna : String
  gc : Double
  cai : Double
  max_run : Int32
  forbidden_count : Int32
}

///|
/// 统计 DNA 的 GC、CAI、最长均聚物与禁止 motif 数.
pub fn analyze_codon_optimization(
  dna : String,
  weights : Map[String, Double],
  table~ : Int32 = 1,
  forbidden_motifs~ : Array[String] = [],
) -> CodonOptimizationReport {
  let mut internal = dna.to_upper()
  for stop in codon_table_stop_codons(table) {
    if internal.has_suffix(stop) {
      internal = internal[0:internal.length() - stop.length()].to_owned()
      break
    }
  }
  {
    dna,
    gc: codon_opt_gc(dna),
    cai: codon_adaptation_index(internal, weights, table~),
    max_run: codon_opt_max_homopolymer(dna),
    forbidden_count: codon_opt_count_forbidden(dna, forbidden_motifs),
  }
}
