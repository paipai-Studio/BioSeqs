///|
/// kallisto_lite.mbt — 伪比对转录本定量(k-mer 等价类 + EM, kallisto 风格).
///
/// 对齐 kallisto (Bray et al. 2016) 的核心量化语义:
/// - 转录组 k-mer 索引: 每 k-mer -> 兼容转录本集合
/// - 伪比对: read 的全部 k-mer 的兼容集合取交集, 得到等价类
///   (kallisto 用 T-DBG + 严格 k-mer 兼容, 此处简化为集合交集, 语义等价)
/// - EM 估计: 等价类计数按 alpha_j/eff_len_j 权重分摊, 迭代至收敛
/// - 输出: est_counts / TPM / 有效长度, 生成 Salmon 风格 quant.sf
///
/// 与 tximport.parse_quant_sf 闭环: quantify -> to_quant_sf ->
/// parse_quant_sf -> tx_import 即可进入 DESeq2/edgeR 下游.
/// 数据不出本机.
pub suberror KallistoLiteError {
  KallistoLiteError(String)
}

///|
/// 转录组 k-mer 伪比对索引.
pub struct KallistoIndex {
  k : Int32
  tx_names : Array[String]
  tx_lens : Array[Int32]
  eff_lens : Array[Double]
  /// k-mer -> 兼容转录本 id 列表
  kmer_map : Map[String, Array[Int32]]
}

///|
/// 一个等价类(兼容转录本集合 + 观测 read 数).
pub(all) struct KallistoEc {
  tx_ids : Array[Int32]
  count : Int32
}

///|
/// 定量结果(每转录本一行).
pub(all) struct KallistoQuant {
  tx_names : Array[String]
  tx_lens : Array[Int32]
  eff_lens : Array[Double]
  est_counts : Array[Double]
  tpm : Array[Double]
  n_reads : Int32
  n_pseudoaligned : Int32
  n_ec : Int32
  iterations : Int32
}

///|
/// 构建转录组 k-mer 索引. k 取 [4, 31], 默认 11(教学级);
/// kallisto 生产默认 31, MoonBit 纯内存实现下小 k 更均衡.
pub fn kallisto_build_index(
  tx_names : Array[String],
  tx_seqs : Array[String],
  k~ : Int32 = 11,
) -> KallistoIndex raise KallistoLiteError {
  if tx_names.length() != tx_seqs.length() {
    raise KallistoLiteError(
      "names (\{tx_names.length()}) and sequences (\{tx_seqs.length()}) count mismatch",
    )
  }
  let kk = if k >= 4 && k <= 31 { k } else { 11 }
  let kmer_map : Map[String, Array[Int32]] = Map([], capacity=0)
  let eff_lens : Array[Double] = Array(capacity=0)
  let tx_lens : Array[Int32] = Array(capacity=0)
  for ti = 0; ti < tx_seqs.length(); ti = ti + 1 {
    let seq = tx_seqs[ti]
    let n = seq.length()
    tx_lens.push(n)
    let eff = if n - kk + 1 > 0 { (n - kk + 1).to_double() } else { 1.0 }
    eff_lens.push(eff)
    let mut p = 0
    while p + kk <= n {
      let key = seq[p:p + kk].to_owned()
      match kmer_map.get(key) {
        Some(ids) =>
          if ids.length() == 0 || ids[ids.length() - 1] != ti {
            ids.push(ti)
          }
        None => {
          let ids : Array[Int32] = Array(capacity=0)
          ids.push(ti)
          kmer_map.set(key, ids)
        }
      }
      p = p + 1
    }
  }
  { k: kk, tx_names, tx_lens, eff_lens, kmer_map, }
}

///|
/// 单条 read 伪比对: 全部 k-mer 兼容集合的交集.
/// 返回兼容转录本 id 列表(升序去重); 空表示未比对上.
pub fn kallisto_pseudoalign_read(
  index : KallistoIndex,
  read : String,
) -> Array[Int32] {
  let n = read.length()
  let k = index.k
  if n < k {
    return []
  }
  let mut result : Array[Int32]? = None
  let mut p = 0
  while p + k <= n {
    let key = read[p:p + k].to_owned()
    match index.kmer_map.get(key) {
      Some(ids) =>
        match result {
          Some(acc) => {
            // 交集(acc 与 ids 均升序)
            let merged : Array[Int32] = Array(capacity=0)
            let mut i = 0
            let mut j = 0
            while i < acc.length() && j < ids.length() {
              if acc[i] == ids[j] {
                merged.push(acc[i])
                i = i + 1
                j = j + 1
              } else if acc[i] < ids[j] {
                i = i + 1
              } else {
                j = j + 1
              }
            }
            result = Some(merged)
          }
          None => {
            let copy : Array[Int32] = Array(capacity=0)
            for x in ids {
              copy.push(x)
            }
            result = Some(copy)
          }
        }
      None => return []
    }
    match result {
      Some(acc) => if acc.length() == 0 { return [] }
      None => ()
    }
    p = p + 1
  }
  match result {
    Some(acc) => acc
    None => []
  }
}

///|
/// 批量伪比对: 聚合为等价类(兼容集合相同即同类).
pub fn kallisto_pseudoalign(
  index : KallistoIndex,
  reads : Array[String],
) -> Array[KallistoEc] {
  let ec_ids : Map[String, Array[Int32]] = Map([], capacity=0)
  let ec_counts : Map[String, Int32] = Map([], capacity=0)
  for read in reads {
    let ids = kallisto_pseudoalign_read(index, read)
    if ids.length() == 0 {
      continue
    }
    let sb = StringBuilder()
    for i = 0; i < ids.length(); i = i + 1 {
      if i > 0 {
        sb.write_char(',')
      }
      sb.write_string(ids[i].to_string())
    }
    let key = sb.to_string()
    match ec_counts.get(key) {
      Some(c) => ec_counts.set(key, c + 1)
      None => {
        ec_counts.set(key, 1)
        ec_ids.set(key, ids)
      }
    }
  }
  let out : Array[KallistoEc] = Array(capacity=0)
  for key in ec_counts.keys() {
    let count = ec_counts.get(key).unwrap()
    let ids = ec_ids.get(key).unwrap()
    out.push({ tx_ids: ids, count, })
  }
  out.sort_by(fn(a, b) {
    if a.tx_ids.length() != b.tx_ids.length() {
      a.tx_ids.length() - b.tx_ids.length()
    } else {
      a.tx_ids[0] - b.tx_ids[0]
    }
  })
  out
}

///|
/// EM 估计: 等价类计数按 alpha/eff_len 权重分摊至兼容转录本.
/// 迭代至 L1 变化 < tolerance 或达上限. 无兼容 read 的转录本计 0.
pub fn kallisto_em(
  index : KallistoIndex,
  ecs : Array[KallistoEc],
  iterations~ : Int32 = 200,
  tolerance~ : Double = 1.0e-6,
) -> (Array[Double], Int32) {
  let n_tx = index.tx_names.length()
  let alpha : Array[Double] = Array(capacity=0)
  for i = 0; i < n_tx; i = i + 1 {
    alpha.push(0.0)
  }
  // 初始化: 每类均匀分摊
  for ec in ecs {
    if ec.tx_ids.length() == 0 {
      continue
    }
    let share = ec.count.to_double() / ec.tx_ids.length().to_double()
    for id in ec.tx_ids {
      alpha[id] = alpha[id] + share
    }
  }
  // 预计算每类的权重基(1/eff_len)
  let inv_eff : Array[Double] = Array(capacity=0)
  for e in index.eff_lens {
    inv_eff.push(1.0 / e)
  }
  let mut iter = 0
  let max_iter = if iterations > 0 { iterations } else { 200 }
  while iter < max_iter {
    let next : Array[Double] = Array(capacity=0)
    for i = 0; i < n_tx; i = i + 1 {
      next.push(0.0)
    }
    for ec in ecs {
      let m = ec.tx_ids.length()
      if m == 0 {
        continue
      }
      let mut denom = 0.0
      for id in ec.tx_ids {
        denom = denom + alpha[id] * inv_eff[id]
      }
      if denom <= 0.0 {
        // 均匀回退
        let share = ec.count.to_double() / m.to_double()
        for id in ec.tx_ids {
          next[id] = next[id] + share
        }
        continue
      }
      for id in ec.tx_ids {
        next[id] = next[id] +
          ec.count.to_double() * (alpha[id] * inv_eff[id]) / denom
      }
    }
    // L1 收敛判定
    let mut delta = 0.0
    for i = 0; i < n_tx; i = i + 1 {
      let d = next[i] - alpha[i]
      let ad = if d >= 0.0 { d } else { -d }
      delta = delta + ad
    }
    for i = 0; i < n_tx; i = i + 1 {
      alpha[i] = next[i]
    }
    iter = iter + 1
    if delta < tolerance {
      break
    }
  }
  (alpha, iter)
}

///|
/// 全流程定量: 索引已建, 输入 read 序列数组, 输出每转录本计数与 TPM.
pub fn kallisto_quantify(
  index : KallistoIndex,
  reads : Array[String],
  iterations~ : Int32 = 200,
) -> KallistoQuant {
  let ecs = kallisto_pseudoalign(index, reads)
  let (alpha, iter) = kallisto_em(index, ecs, iterations~)
  // TPM: alpha/eff_len 归一化 * 1e6
  let tpm : Array[Double] = Array(capacity=0)
  let mut rate_sum = 0.0
  for i = 0; i < alpha.length(); i = i + 1 {
    rate_sum = rate_sum + alpha[i] / index.eff_lens[i]
  }
  for i = 0; i < alpha.length(); i = i + 1 {
    let v = if rate_sum > 0.0 {
      alpha[i] / index.eff_lens[i] / rate_sum * 1.0e6
    } else {
      0.0
    }
    tpm.push(v)
  }
  let mut n_aligned = 0
  for ec in ecs {
    n_aligned = n_aligned + ec.count
  }
  {
    tx_names: index.tx_names,
    tx_lens: index.tx_lens,
    eff_lens: index.eff_lens,
    est_counts: alpha,
    tpm,
    n_reads: reads.length(),
    n_pseudoaligned: n_aligned,
    n_ec: ecs.length(),
    iterations: iter,
  }
}

///|
/// 一步式便捷入口: 建索引并定量.
pub fn kallisto_run(
  tx_names : Array[String],
  tx_seqs : Array[String],
  reads : Array[String],
  k~ : Int32 = 11,
  iterations~ : Int32 = 200,
) -> KallistoQuant raise KallistoLiteError {
  let index = kallisto_build_index(tx_names, tx_seqs, k~)
  kallisto_quantify(index, reads, iterations~)
}

///|
/// 生成 Salmon 风格 quant.sf 文本(4 列: Name Length EffectiveLength TPM NumReads).
/// 与 tximport.parse_quant_sf 闭环.
pub fn KallistoQuant::to_quant_sf(self : KallistoQuant) -> String {
  let sb = StringBuilder()
  sb.write_string("Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
  for i = 0; i < self.tx_names.length(); i = i + 1 {
    sb.write_string(self.tx_names[i])
    sb.write_char('\t')
    sb.write_string(self.tx_lens[i].to_string())
    sb.write_char('\t')
    sb.write_string(kl_fmt(self.eff_lens[i]))
    sb.write_char('\t')
    sb.write_string(kl_fmt(self.tpm[i]))
    sb.write_char('\t')
    sb.write_string(kl_fmt(self.est_counts[i]))
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
fn kl_fmt(x : Double) -> String {
  let s = x.to_string()
  if s.contains(".") || s.contains("e") || s.contains("E") {
    s
  } else {
    s + ".0"
  }
}

///|
/// 定量摘要(用于报告).
pub fn KallistoQuant::summary(self : KallistoQuant) -> String {
  let n = self.tx_names.length()
  let mut tp = 0.0
  for t in self.tpm {
    tp = tp + t
  }
  let mut tc = 0.0
  for c in self.est_counts {
    tc = tc + c
  }
  let rate = if self.n_reads > 0 {
    self.n_pseudoaligned.to_double() / self.n_reads.to_double()
  } else {
    0.0
  }
  "kallisto-lite: \{n} transcripts, \{self.n_ec} equivalence classes, " +
  "\{self.n_pseudoaligned}/\{self.n_reads} reads pseudoaligned (" +
  (rate * 100.0).to_string() +
  "%), total TPM \{tp.to_string()}, total counts \{tc.to_string()}, EM converged in \{self.iterations} iterations"
}
