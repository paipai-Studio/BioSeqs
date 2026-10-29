///|
/// motif_discovery.mbt — de novo motif 发现(EM 与 Gibbs 采样, MEME/AlignACE 风格).
///
/// 对齐 MEME 的核心发现语义(本地确定性实现):
/// - EM(OOPS/ZOOPS): 期望-最大化交替, 位点后验 -> 加权 PWM 重估,
///   双链搜索, 伪计数平滑, 似然收敛判定
/// - Gibbs 采样(AlignACE 风格): 逐序列留一采样, 确定性 RNG(统一 Rng),
///   多重启取最优信息内容解
/// - 评分: 位点似然 / 信息含量 IC(相对背景), 输出 consensus 与
///   per-site 位置/打分
/// - 输出 Motif(与 pwm_from_counts/motif_scan/meme 序列化栈闭环)
///
/// 背景: 均匀 + 观测碱基组成混合. 数据不出本机.
pub suberror MotifDiscoveryError {
  MotifDiscoveryError(String)
}

///|
/// 发现结果(含 motif 与位点指认).
pub(all) struct DiscoveredMotif {
  name : String
  pwm : Array[Array[Double]]
  alphabet : Array[String]
  sites : Array[String]
  /// 每条序列的 (best_start, best_strand, score)
  site_positions : Array[(Int32, Int32, Double)]
  information_content : Double
  log_likelihood : Double
  width : Int32
  n_sequences : Int32
}

///|
fn motif_disco_alphabet() -> Array[String] {
  ["A", "C", "G", "T"]
}

///|
fn motif_disco_base_index(c : Char) -> Int32 {
  match c {
    'A' => 0
    'C' => 1
    'G' => 2
    'T' => 3
    _ => -1
  }
}

///|
/// 双链位点枚举: (start, strand) 全部 w-mer.
fn motif_disco_sites(seq : String, w : Int32) -> Array[(Int32, Int32, String)] {
  let out : Array[(Int32, Int32, String)] = Array(capacity=0)
  let n = seq.length()
  let mut pos = 0
  while pos + w <= n {
    out.push((pos, 0, seq[pos:pos + w].to_owned()))
    pos = pos + 1
  }
  // 反向链
  let rc = reverse_complement(seq)
  let mut p2 = 0
  while p2 + w <= n {
    // 反向链起点 start 对应原链坐标 [n - p2 - w, n - p2)
    out.push((n - p2 - w, 1, rc[p2:p2 + w].to_owned()))
    p2 = p2 + 1
  }
  out
}

///|
fn motif_disco_bg_freqs(
  sequences : Array[String],
  uniform_mix~ : Double = 0.5,
) -> Array[Double] {
  let counts : Array[Double] = [0.0, 0.0, 0.0, 0.0]
  let mut total = 0.0
  for seq in sequences {
    for i = 0; i < seq.length(); i = i + 1 {
      let bi = motif_disco_base_index(seq.unsafe_get(i).unsafe_to_char())
      if bi >= 0 {
        counts[bi] = counts[bi] + 1.0
        total = total + 1.0
      }
    }
  }
  if total == 0.0 {
    return [0.25, 0.25, 0.25, 0.25]
  }
  let mix = if uniform_mix < 0.0 {
    0.0
  } else if uniform_mix > 1.0 {
    1.0
  } else {
    uniform_mix
  }
  [0, 1, 2, 3].map(fn(i : Int32) -> Double {
    mix * 0.25 + (1.0 - mix) * counts[i] / total
  })
}

///|
/// w-mer 在 PWM 下的概率(相对背景归一化前的原始 P).
fn motif_disco_p_site(site : String, pwm : Array[Array[Double]]) -> Double {
  let mut p = 1.0
  for i = 0; i < site.length(); i = i + 1 {
    let bi = motif_disco_base_index(site.unsafe_get(i).unsafe_to_char())
    if bi < 0 {
      return 0.0
    }
    p = p * pwm[i][bi]
  }
  p
}

///|
/// 信息含量 IC(bits/position 总和, 相对背景).
fn motif_disco_ic(pwm : Array[Array[Double]], bg : Array[Double]) -> Double {
  let mut ic = 0.0
  for row in pwm {
    for a = 0; a < 4; a = a + 1 {
      let p = row[a]
      if p > 0.0 && bg[a] > 0.0 {
        ic = ic + p * @math.ln(p / bg[a]) / @math.ln(2.0)
      }
    }
  }
  if ic < 0.0 {
    0.0
  } else {
    ic
  }
}

///|
fn motif_disco_init_pwm(
  seed_site : String,
  pseudocount : Double,
) -> Array[Array[Double]] {
  let w = seed_site.length()
  let pwm : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < w; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for a = 0; a < 4; a = a + 1 {
      row.push(pseudocount)
    }
    let bi = motif_disco_base_index(seed_site.unsafe_get(i).unsafe_to_char())
    if bi >= 0 {
      row[bi] = row[bi] + 1.0
    }
    let total = 1.0 + pseudocount * 4.0
    let norm : Array[Double] = Array(capacity=0)
    for a = 0; a < 4; a = a + 1 {
      norm.push(row[a] / total)
    }
    pwm.push(norm)
  }
  pwm
}

///|
fn motif_disco_counts_to_pwm(
  counts : Array[Array[Double]],
  n_sites : Double,
  pseudocount : Double,
) -> Array[Array[Double]] {
  let w = counts.length()
  let pwm : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < w; i = i + 1 {
    let total = n_sites + pseudocount * 4.0
    let row : Array[Double] = Array(capacity=0)
    for a = 0; a < 4; a = a + 1 {
      row.push((counts[i][a] + pseudocount) / total)
    }
    pwm.push(row)
  }
  pwm
}

///|
/// 单轮 EM 主循环(从给定初始 PWM 迭代).
fn motif_em_run(
  all_sites : Array[Array[(Int32, Int32, String)]],
  pwm0 : Array[Array[Double]],
  bg : Array[Double],
  width : Int32,
  iterations : Int32,
  tolerance : Double,
  pseudocount : Double,
) -> (Array[Array[Double]], Double) {
  let n_seq = all_sites.length()
  let w = width
  let mut pwm = pwm0
  let mut prev_ll = -1.0e18
  let mut iter = 0
  let posteriors : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n_seq; i = i + 1 {
    posteriors.push(Array(capacity=0))
  }
  while iter < iterations {
    // E 步: 位点后验 ∝ P(site|pwm) / P(site|bg)
    let mut ll = 0.0
    for si = 0; si < n_seq; si = si + 1 {
      posteriors[si].truncate(0)
      let sites = all_sites[si]
      let mut zsum = 0.0
      for site in sites {
        let p_motif = motif_disco_p_site(site.2, pwm)
        let mut bgp = 1.0
        for c = 0; c < site.2.length(); c = c + 1 {
          let bi = motif_disco_base_index(site.2.unsafe_get(c).unsafe_to_char())
          if bi < 0 {
            bgp = 0.0
            break
          }
          bgp = bgp * bg[bi]
        }
        let z = if bgp > 0.0 { p_motif / bgp } else { 0.0 }
        posteriors[si].push(z)
        zsum = zsum + z
      }
      if zsum > 0.0 {
        for j = 0; j < posteriors[si].length(); j = j + 1 {
          posteriors[si][j] = posteriors[si][j] / zsum
        }
        ll = ll + @math.ln(zsum)
      }
    }
    // M 步: 加权计数 -> PWM
    let counts : Array[Array[Double]] = Array(capacity=0)
    for i = 0; i < w; i = i + 1 {
      counts.push([0.0, 0.0, 0.0, 0.0])
    }
    let mut n_weight = 0.0
    for si = 0; si < n_seq; si = si + 1 {
      let sites = all_sites[si]
      for j = 0; j < sites.length(); j = j + 1 {
        let weight = posteriors[si][j]
        if weight <= 0.0 {
          continue
        }
        n_weight = n_weight + weight
        let site = sites[j].2
        for c = 0; c < w; c = c + 1 {
          let bi = motif_disco_base_index(site.unsafe_get(c).unsafe_to_char())
          if bi >= 0 {
            counts[c][bi] = counts[c][bi] + weight
          }
        }
      }
    }
    if n_weight <= 0.0 {
      break
    }
    pwm = motif_disco_counts_to_pwm(counts, n_weight, pseudocount)
    if ll - prev_ll < tolerance && ll > -1.0e17 && iter > 0 {
      prev_ll = ll
      break
    }
    prev_ll = ll
    iter = iter + 1
  }
  (pwm, prev_ll)
}

///|
/// EM de novo motif 发现(MEME 风格多起点: 首序列 w-mer 种子枚举,
/// 短程 EM 择优后续跑到收敛). OOPS 每序列一位点.
pub fn motif_em(
  sequences : Array[String],
  width~ : Int32 = 10,
  iterations~ : Int32 = 100,
  tolerance~ : Double = 1.0e-6,
  pseudocount~ : Double = 0.25,
) -> DiscoveredMotif raise MotifDiscoveryError {
  if sequences.length() == 0 {
    raise MotifDiscoveryError("empty sequence set")
  }
  let w = if width >= 3 && width <= 30 { width } else { 10 }
  for seq in sequences {
    if seq.length() < w {
      raise MotifDiscoveryError(
        "sequence shorter than motif width (\{seq.length()} < \{w})",
      )
    }
  }
  let bg = motif_disco_bg_freqs(sequences)
  let all_sites : Array[Array[(Int32, Int32, String)]] = Array(capacity=0)
  for seq in sequences {
    all_sites.push(motif_disco_sites(seq, w))
  }
  // 多起点: 首序列正向链 w-mer 枚举(上限 15 个种子)
  let first = sequences[0]
  let n_seeds = if first.length() - w + 1 > 15 {
    15
  } else {
    first.length() - w + 1
  }
  let mut best_pwm : Array[Array[Double]]? = None
  let mut best_ll = -1.0e18
  for s = 0; s < n_seeds; s = s + 1 {
    let seed_site = first[s:s + w].to_owned()
    let pwm0 = motif_disco_init_pwm(seed_site, pseudocount)
    let (pwm, ll) = motif_em_run(
      all_sites, pwm0, bg, w, 15, tolerance, pseudocount,
    )
    if ll > best_ll {
      best_ll = ll
      best_pwm = Some(pwm)
    }
  }
  // 最优种子继续跑到收敛
  let pwm_final = match best_pwm {
    Some(p) => p
    None => motif_disco_init_pwm(first[0:w].to_owned(), pseudocount)
  }
  let (pwm, ll) = motif_em_run(
    all_sites, pwm_final, bg, w, iterations, tolerance, pseudocount,
  )
  motif_disco_finalize(sequences, all_sites, pwm, bg, ll, 15, "EM")
}

///|
/// Gibbs 采样 de novo motif 发现(AlignACE 风格, 确定性 RNG).
/// 多重启取 IC 最优解; 每轮逐序列留一采样新位点.
pub fn motif_gibbs(
  sequences : Array[String],
  width~ : Int32 = 10,
  restarts~ : Int32 = 3,
  sweeps~ : Int32 = 50,
  seed~ : UInt64 = 42,
  pseudocount~ : Double = 0.25,
) -> DiscoveredMotif raise MotifDiscoveryError {
  if sequences.length() == 0 {
    raise MotifDiscoveryError("empty sequence set")
  }
  let w = if width >= 3 && width <= 30 { width } else { 10 }
  for seq in sequences {
    if seq.length() < w {
      raise MotifDiscoveryError(
        "sequence shorter than motif width (\{seq.length()} < \{w})",
      )
    }
  }
  let bg = motif_disco_bg_freqs(sequences)
  let all_sites : Array[Array[(Int32, Int32, String)]] = Array(capacity=0)
  for seq in sequences {
    all_sites.push(motif_disco_sites(seq, w))
  }
  let n_seq = sequences.length()
  let rng = Rng::new(seed)
  let mut best : DiscoveredMotif? = None
  let mut best_ic = -1.0
  for restart = 0; restart < restarts; restart = restart + 1 {
    // 随机初始化: 每序列一个位点
    let starts : Array[Int32] = Array(capacity=0)
    for si = 0; si < n_seq; si = si + 1 {
      let sites = all_sites[si]
      let pick = rng.next_int(sites.length())
      starts.push(pick)
    }
    // Gibbs 扫描
    let mut sweep = 0
    while sweep < sweeps {
      for si = 0; si < n_seq; si = si + 1 {
        // 留一: 除 si 外的位点构建 PWM
        let counts : Array[Array[Double]] = Array(capacity=0)
        for i = 0; i < w; i = i + 1 {
          counts.push([0.0, 0.0, 0.0, 0.0])
        }
        let mut n_used = 0.0
        for sj = 0; sj < n_seq; sj = sj + 1 {
          if sj == si {
            continue
          }
          let site = all_sites[sj][starts[sj]].2
          for c = 0; c < w; c = c + 1 {
            let bi = motif_disco_base_index(site.unsafe_get(c).unsafe_to_char())
            if bi >= 0 {
              counts[c][bi] = counts[c][bi] + 1.0
            }
          }
          n_used = n_used + 1.0
        }
        if n_used <= 0.0 {
          continue
        }
        let pwm_loo = motif_disco_counts_to_pwm(counts, n_used, pseudocount)
        // 采样: 位点权重 ∝ P(site|pwm)/P(site|bg)
        let sites = all_sites[si]
        let weights : Array[Double] = Array(capacity=0)
        let mut wsum = 0.0
        for site in sites {
          let pm = motif_disco_p_site(site.2, pwm_loo)
          let mut bgp = 1.0
          for c = 0; c < site.2.length(); c = c + 1 {
            let bi = motif_disco_base_index(
              site.2.unsafe_get(c).unsafe_to_char(),
            )
            if bi < 0 {
              bgp = 0.0
              break
            }
            bgp = bgp * bg[bi]
          }
          let wt = if bgp > 0.0 { pm / bgp } else { 0.0 }
          weights.push(wt)
          wsum = wsum + wt
        }
        if wsum <= 0.0 {
          continue
        }
        // 轮盘采样(确定性 RNG)
        let target = rng.next_double() * wsum
        let mut acc = 0.0
        let mut chosen = sites.length() - 1
        for j = 0; j < sites.length(); j = j + 1 {
          acc = acc + weights[j]
          if acc >= target {
            chosen = j
            break
          }
        }
        starts[si] = chosen
      }
      sweep = sweep + 1
    }
    // 最终 PWM 与似然
    let counts : Array[Array[Double]] = Array(capacity=0)
    for i = 0; i < w; i = i + 1 {
      counts.push([0.0, 0.0, 0.0, 0.0])
    }
    for si = 0; si < n_seq; si = si + 1 {
      let site = all_sites[si][starts[si]].2
      for c = 0; c < w; c = c + 1 {
        let bi = motif_disco_base_index(site.unsafe_get(c).unsafe_to_char())
        if bi >= 0 {
          counts[c][bi] = counts[c][bi] + 1.0
        }
      }
    }
    let pwm_final = motif_disco_counts_to_pwm(
      counts,
      n_seq.to_double(),
      pseudocount,
    )
    let mut final_ll = 0.0
    for si = 0; si < n_seq; si = si + 1 {
      let site = all_sites[si][starts[si]].2
      final_ll = final_ll + @math.ln(motif_disco_p_site(site, pwm_final))
    }
    let result = motif_disco_finalize(
      sequences, all_sites, pwm_final, bg, final_ll, sweeps, "Gibbs",
    )
    let ic = result.information_content
    if ic > best_ic {
      best_ic = ic
      best = Some(result)
    }
  }
  match best {
    Some(b) => b
    None => raise MotifDiscoveryError("gibbs sampling failed")
  }
}

///|
/// 由最终 PWM 与 E 步后验取各位点最优指认, 组装结果.
fn motif_disco_finalize(
  sequences : Array[String],
  all_sites : Array[Array[(Int32, Int32, String)]],
  pwm : Array[Array[Double]],
  bg : Array[Double],
  ll : Double,
  iterations : Int32,
  algorithm : String,
) -> DiscoveredMotif {
  let n_seq = sequences.length()
  let positions : Array[(Int32, Int32, Double)] = Array(capacity=0)
  let site_strs : Array[String] = Array(capacity=0)
  for si = 0; si < n_seq; si = si + 1 {
    let mut best_score = -1.0
    let mut best : (Int32, Int32, String) = all_sites[si][0]
    for site in all_sites[si] {
      let s = motif_disco_p_site(site.2, pwm)
      if s > best_score {
        best_score = s
        best = site
      }
    }
    positions.push((best.0, best.1, best_score))
    site_strs.push(best.2)
  }
  let ic = motif_disco_ic(pwm, bg)
  {
    name: "discovered-\{algorithm}",
    pwm,
    alphabet: motif_disco_alphabet(),
    sites: site_strs,
    site_positions: positions,
    information_content: ic,
    log_likelihood: ll,
    width: pwm.length(),
    n_sequences: n_seq,
  }
}

///|
/// 发现结果 consensus(与 get_consensus 语义一致).
pub fn DiscoveredMotif::consensus(self : DiscoveredMotif) -> String {
  let sb = StringBuilder()
  for row in self.pwm {
    let mut best_i = 0
    let mut best_v = -1.0
    for a = 0; a < row.length(); a = a + 1 {
      if row[a] > best_v {
        best_v = row[a]
        best_i = a
      }
    }
    match self.alphabet[best_i] {
      "A" => sb.write_char('A')
      "C" => sb.write_char('C')
      "G" => sb.write_char('G')
      "T" => sb.write_char('T')
      _ => sb.write_char('N')
    }
  }
  sb.to_string()
}

///|
/// 结果摘要行.
pub fn DiscoveredMotif::summary(self : DiscoveredMotif) -> String {
  "motif-discovery(\{self.name}): width=\{self.width}, " +
  "sequences=\{self.n_sequences}, IC=\{self.information_content} bits, " +
  "logLik=\{self.log_likelihood}, consensus=\{self.consensus()}"
}

///|
/// 转换为库内 Motif(进入 motif_scan/meme 序列化栈).
pub fn DiscoveredMotif::to_motif(self : DiscoveredMotif) -> Motif {
  {
    name: self.name,
    pwm: self.pwm,
    alphabet: self.alphabet,
    sites: self.sites,
  }
}
