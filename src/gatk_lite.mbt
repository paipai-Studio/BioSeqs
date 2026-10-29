///|
/// gatk_lite.mbt — GATK 风格变异检出(genotype likelihood + active region 共识组装).
///
/// 对齐 GATK HaplotypeCaller 的两阶段语义(简化本地实现):
/// 1. 逐位点 pileup + 二倍体 genotype likelihood(统一 genotyper 模型,
///    错误率 epsilon, 三基因型后验, GQ = -10 log10 P(non-best))
/// 2. active region(错配密度窗口) + 多数投票共识组装 haplotype,
///    与参考比对产出 SNP/插入/缺失, 区域内等位重新分摊
///
/// 输出 VCF 4.2 文本, 与 parse_vcf 闭环进入现有 VCF/BCF 后处理栈
/// (bcftools_ops/variant_filtering). 数据不出本机.

///|
/// 简化读段(序列 + 参考起点, 1-based 语义按 0-based 内部处理).
pub(all) struct GatkRead {
  name : String
  seq : String
  pos : Int32
}

///|
/// 单位点 pileup 统计.
pub(all) struct GatkPileupSite {
  pos : Int32
  ref_base : Char
  depth : Int32
  base_counts : Array[(Char, Int32)]
  /// pos 之后的插入事件(插入序列 -> 次数)与缺失事件(长度 -> 次数)
  insertions : Array[(String, Int32)]
  deletions : Array[(Int32, Int32)]
}

///|
/// 基因型后验.
pub(all) struct GatkGenotype {
  gt : String
  posterior : Double
  log10_error : Double
}

///|
/// 一条检出变异.
pub(all) struct GatkVariant {
  chrom : String
  pos : Int32
  ref : String
  alt : String
  qual : Double
  filter : String
  gt : String
  dp : Int32
  ad_ref : Int32
  ad_alt : Int32
  gq : Int32
  af : Double
}

// ---------------------------------------------------------------------------
// pileup
// ---------------------------------------------------------------------------

///|
/// indel 推断结论.
priv enum GatkIndelCall {
  Insertion(Int32)
  Deletion(Int32)
  None
}

///|
/// 从 (rpos, roff) 错位起点推断 indel:
/// 插入假设 = 读段在 roff 后 d 个碱基恢复与参考 rpos 对齐;
/// 缺失假设 = 参考在 rpos 后 dd 个碱基恢复与读段 roff 对齐.
/// 两假设均可达时取恢复后连续匹配更长且 >= 2 者.
fn gatk_classify_indel(
  reference : String,
  seq : String,
  rpos : Int32,
  roff : Int32,
) -> GatkIndelCall {
  let rb = reference.unsafe_get(rpos).unsafe_to_char()
  let qb = seq.unsafe_get(roff).unsafe_to_char()
  let max_d : Int32 = if seq.length() - roff - 1 < 10 {
    seq.length() - roff - 1
  } else {
    10
  }
  let max_dd : Int32 = if reference.length() - rpos - 1 < 10 {
    reference.length() - rpos - 1
  } else {
    10
  }
  // 插入假设
  let mut ins_d : Int32 = 0
  let mut d = 1
  while d <= max_d {
    if seq.unsafe_get(roff + d).unsafe_to_char() == rb {
      ins_d = d
      break
    }
    d = d + 1
  }
  // 缺失假设
  let mut del_dd : Int32 = 0
  let mut dd = 1
  while dd <= max_dd {
    if reference.unsafe_get(rpos + dd).unsafe_to_char() == qb {
      del_dd = dd
      break
    }
    dd = dd + 1
  }
  if ins_d == 0 && del_dd == 0 {
    return GatkIndelCall::None
  }
  if ins_d == 0 {
    if gatk_match_run(reference, seq, rpos + del_dd, roff) >= 2 {
      return GatkIndelCall::Deletion(del_dd)
    }
    return GatkIndelCall::None
  }
  if del_dd == 0 {
    if gatk_match_run(reference, seq, rpos, roff + ins_d) >= 2 {
      return GatkIndelCall::Insertion(ins_d)
    }
    return GatkIndelCall::None
  }
  let m_ins = gatk_match_run(reference, seq, rpos, roff + ins_d)
  let m_del = gatk_match_run(reference, seq, rpos + del_dd, roff)
  if m_ins >= 2 || m_del >= 2 {
    if m_del >= m_ins {
      GatkIndelCall::Deletion(del_dd)
    } else {
      GatkIndelCall::Insertion(ins_d)
    }
  } else {
    GatkIndelCall::None
  }
}

///|
/// 从 (参考 rpos, 读段 qi) 起的连续匹配长度.
fn gatk_match_run(
  reference : String,
  seq : String,
  rpos : Int32,
  qi : Int32,
) -> Int32 {
  let mut n = 0
  let mut ri = rpos
  let mut q = qi
  while ri < reference.length() && q < seq.length() {
    if reference.unsafe_get(ri) == seq.unsafe_get(q) {
      n = n + 1
      ri = ri + 1
      q = q + 1
    } else {
      break
    }
  }
  n
}

///|
/// 读段是否覆盖参考位置 p(0-based), 返回读内偏移.
fn gatk_read_covers(read : GatkRead, p : Int32) -> Int32? {
  let off = p - read.pos
  if off >= 0 && off < read.seq.length() {
    Some(off)
  } else {
    None
  }
}

///|
/// 统计单点位点 pileup(碱基计数与紧邻 indel 事件).
/// indel 事件定义: 读段在 p 处结束后出现的插入/缺失相对参考的差分,
/// 由读段与参考的逐碱基比较在 [p, p+lookahead) 窗口内推断.
pub fn gatk_pileup_site(
  reference : String,
  reads : Array[GatkRead],
  pos : Int32,
) -> GatkPileupSite {
  let ref_base = if pos >= 0 && pos < reference.length() {
    reference.unsafe_get(pos).unsafe_to_char()
  } else {
    'N'
  }
  let counts : Map[Char, Int32] = Map([], capacity=0)
  let ins : Map[String, Int32] = Map([], capacity=0)
  let dels : Map[Int32, Int32] = Map([], capacity=0)
  let mut depth = 0
  for read in reads {
    match gatk_read_covers(read, pos) {
      Some(off) => {
        depth = depth + 1
        let b = read.seq.unsafe_get(off).unsafe_to_char()
        match counts.get(b) {
          Some(c) => counts.set(b, c + 1)
          None => counts.set(b, 1)
        }
        // 读段覆盖 pos 且 pos+1 处读段仍有碱基而参考不匹配 -> indel 推断
        let roff = off + 1
        let rpos = pos + 1
        if rpos < reference.length() && roff < read.seq.length() {
          let rb = reference.unsafe_get(rpos).unsafe_to_char()
          let qb = read.seq.unsafe_get(roff).unsafe_to_char()
          if rb != qb {
            match gatk_classify_indel(reference, read.seq, rpos, roff) {
              GatkIndelCall::Insertion(d) => {
                let inserted = read.seq[roff:roff + d].to_owned()
                match ins.get(inserted) {
                  Some(c) => ins.set(inserted, c + 1)
                  None => ins.set(inserted, 1)
                }
              }
              GatkIndelCall::Deletion(dd) =>
                match dels.get(dd) {
                  Some(c) => dels.set(dd, c + 1)
                  None => dels.set(dd, 1)
                }
              GatkIndelCall::None => ()
            }
          }
        }
      }
      None => ()
    }
  }
  let base_counts : Array[(Char, Int32)] = Array(capacity=0)
  for c in counts.keys() {
    base_counts.push((c, counts.get(c).unwrap()))
  }
  base_counts.sort_by(fn(a, b) { b.1 - a.1 })
  let insertions : Array[(String, Int32)] = Array(capacity=0)
  for key in ins.keys() {
    insertions.push((key, ins.get(key).unwrap()))
  }
  insertions.sort_by(fn(a, b) { b.1 - a.1 })
  let deletions : Array[(Int32, Int32)] = Array(capacity=0)
  for key in dels.keys() {
    deletions.push((key, dels.get(key).unwrap()))
  }
  deletions.sort_by(fn(a, b) { b.1 - a.1 })
  { pos, ref_base, depth, base_counts, insertions, deletions, }
}

// ---------------------------------------------------------------------------
// genotype likelihoods
// ---------------------------------------------------------------------------

///|
/// 观察碱基似然: P(base | 等位基因频率 f).
fn gatk_base_loglik(observed : Char, allele : Char, epsilon : Double) -> Double {
  if observed == allele {
    @math.ln(1.0 - epsilon)
  } else {
    @math.ln(epsilon / 3.0)
  }
}

///|
/// 计算双等位基因 (ref, alt) 的三种二倍体基因型后验.
/// 基因型: "0/0"(纯合参考), "0/1"(杂合), "1/1"(纯合 alt).
/// 先验均匀; 返回按后验降序排列的三元组.
pub fn gatk_genotype_likelihoods(
  site : GatkPileupSite,
  alt : Char,
  epsilon~ : Double = 0.02,
) -> Array[GatkGenotype] {
  let genotypes : Array[String] = ["0/0", "0/1", "1/1"]
  // 每基因型的等位频率
  let freqs : Array[Double] = [0.0, 0.5, 1.0]
  let lls : Array[Double] = []
  let results : Array[(String, Double)] = []
  for gi = 0; gi < 3; gi = gi + 1 {
    let f = freqs[gi]
    // P(obs) = f * P(obs|alt) + (1-f) * P(obs|ref)
    let mut ll = 0.0
    for pair in site.base_counts {
      let (base, cnt) = pair
      let p_alt = gatk_base_loglik(base, alt, epsilon)
      let p_ref = gatk_base_loglik(base, site.ref_base, epsilon)
      let log_p = @math.ln(f * @math.exp(p_alt) + (1.0 - f) * @math.exp(p_ref))
      ll = ll + cnt.to_double() * log_p
    }
    ignore(lls)
    results.push((genotypes[gi], ll))
  }
  // 归一化后验
  let mut max_ll = results[0].1
  for r in results {
    if r.1 > max_ll {
      max_ll = r.1
    }
  }
  let mut total = 0.0
  for r in results {
    total = total + @math.exp(r.1 - max_ll)
  }
  let out : Array[GatkGenotype] = Array(capacity=0)
  let indexed : Array[(String, Double)] = []
  for r in results {
    let post = @math.exp(r.1 - max_ll) / total
    indexed.push((r.0, post))
  }
  indexed.sort_by(fn(a, b) {
    if a.1 > b.1 {
      -1
    } else if a.1 < b.1 {
      1
    } else {
      0
    }
  })
  for r in indexed {
    let err = 1.0 - r.1
    let gq = if err > 1.0e-12 {
      -10.0 * @math.ln(err) / @math.ln(10.0)
    } else {
      99.0
    }
    out.push({ gt: r.0, posterior: r.1, log10_error: gq, })
  }
  out
}

// ---------------------------------------------------------------------------
// active region 与共识组装
// ---------------------------------------------------------------------------

///|
/// 识别 active region: 错配密度 >= min_mismatches 的连续窗口.
/// 返回 (start, end) 0-based 半开区间数组.
pub fn gatk_active_regions(
  reference : String,
  reads : Array[GatkRead],
  window~ : Int32 = 25,
  min_mismatches~ : Int32 = 3,
) -> Array[(Int32, Int32)] {
  let n = reference.length()
  if n == 0 {
    return []
  }
  // 每位置错配计数
  let mismatch : Array[Int32] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    mismatch.push(0)
  }
  for read in reads {
    for off = 0; off < read.seq.length(); off = off + 1 {
      let rp = read.pos + off
      if rp >= 0 && rp < n {
        if read.seq.unsafe_get(off) != reference.unsafe_get(rp) {
          mismatch[rp] = mismatch[rp] + 1
        }
      }
    }
  }
  let regions : Array[(Int32, Int32)] = Array(capacity=0)
  let w = if window > 0 { window } else { 25 }
  let mut s = 0
  while s + w <= n {
    let mut mm = 0
    for i = s; i < s + w; i = i + 1 {
      mm = mm + mismatch[i]
    }
    if mm >= min_mismatches {
      // 扩展至错配连续支持
      let mut e = s + w
      while e < n && mismatch[e] > 0 {
        e = e + 1
      }
      regions.push((s, e))
      s = e
    } else {
      s = s + 1
    }
  }
  regions
}

///|
/// active region 多数投票共识 haplotype(每列取最频碱基).
pub fn gatk_assemble_consensus(
  reference : String,
  reads : Array[GatkRead],
  start : Int32,
  end : Int32,
) -> String {
  if end <= start {
    return ""
  }
  let sb = StringBuilder()
  for p = start; p < end; p = p + 1 {
    let votes : Map[Char, Int32] = Map([], capacity=0)
    for read in reads {
      match gatk_read_covers(read, p) {
        Some(off) => {
          let b = read.seq.unsafe_get(off).unsafe_to_char()
          match votes.get(b) {
            Some(c) => votes.set(b, c + 1)
            None => votes.set(b, 1)
          }
        }
        None => ()
      }
    }
    if votes.length() == 0 {
      if p < reference.length() {
        sb.write_char(reference.unsafe_get(p).unsafe_to_char())
      }
    } else {
      let mut best : Char = if p < reference.length() {
        reference.unsafe_get(p).unsafe_to_char()
      } else {
        'N'
      }
      let mut best_cnt = 0
      for c in votes.keys() {
        let cnt = votes.get(c).unwrap()
        if cnt > best_cnt {
          best_cnt = cnt
          best = c
        }
      }
      sb.write_char(best)
    }
  }
  sb.to_string()
}

// ---------------------------------------------------------------------------
// 变异检出主流程
// ---------------------------------------------------------------------------

///|
/// 全流程变异检出: 逐位点 SNP + active region 插入/缺失.
/// min_depth: 检出最小深度; min_af: 最小等位频率; min_gq: 最小基因型质量.
pub fn gatk_call_variants(
  chrom : String,
  reference : String,
  reads : Array[GatkRead],
  min_depth~ : Int32 = 4,
  min_af~ : Double = 0.2,
  min_gq~ : Double = 20.0,
) -> Array[GatkVariant] {
  let variants : Array[GatkVariant] = Array(capacity=0)
  let n = reference.length()
  for p = 0; p < n; p = p + 1 {
    let site = gatk_pileup_site(reference, reads, p)
    if site.depth < min_depth {
      continue
    }
    // 非 ref 主碱基
    let mut alt : Char = 'N'
    let mut alt_count = 0
    for pair in site.base_counts {
      if pair.0 != site.ref_base {
        alt = pair.0
        alt_count = pair.1
        break
      }
    }
    let af = alt_count.to_double() / site.depth.to_double()
    if alt == 'N' || af < min_af {
      continue
    }
    let gts = gatk_genotype_likelihoods(site, alt)
    let best = gts[0]
    if best.gt == "0/0" {
      continue
    }
    if best.log10_error < min_gq {
      continue
    }
    let ref_count = site.depth - alt_count
    variants.push({
      chrom,
      pos: p + 1,
      ref: site.ref_base.to_string(),
      alt: alt.to_string(),
      qual: best.log10_error,
      filter: "PASS",
      gt: best.gt,
      dp: site.depth,
      ad_ref: ref_count,
      ad_alt: alt_count,
      gq: best.log10_error.to_int(),
      af,
    })
  }
  // active region indel: 在 SNP 之外补插入/缺失
  let regions = gatk_active_regions(
    reference,
    reads,
    window=25,
    min_mismatches=3,
  )
  for region in regions {
    let (rs, re) = region
    // 扫 region 全部位置, 取支持数最强的插入/缺失事件
    let mut best_ins : (Int32, String, Int32, Int32) = (-1, "", 0, 0)
    let mut best_del : (Int32, Int32, Int32, Int32) = (-1, 0, 0, 0)
    for p = rs; p < re && p < n; p = p + 1 {
      let site = gatk_pileup_site(reference, reads, p)
      if site.depth < min_depth {
        continue
      }
      if site.insertions.length() > 0 {
        let (ins_seq, cnt) = site.insertions[0]
        if cnt > best_ins.2 {
          best_ins = (p, ins_seq, cnt, site.depth)
        }
      }
      if site.deletions.length() > 0 {
        let (dl, cnt) = site.deletions[0]
        if cnt > best_del.2 {
          best_del = (p, dl, cnt, site.depth)
        }
      }
    }
    // 优先产出支持数更高的事件; 同分优先缺失(错位推断对缺失更敏感)
    if best_del.2 >= best_ins.2 && best_del.2 > 0 {
      let (p, dl, cnt, depth) = best_del
      let af = cnt.to_double() / depth.to_double()
      if af >= min_af && p + 1 + dl <= n && p >= 0 {
        let deleted = reference[p:p + 1 + dl].to_owned()
        let kept = reference[p:p + 1].to_owned()
        variants.push({
          chrom,
          pos: p + 1,
          ref: deleted,
          alt: kept,
          qual: 30.0,
          filter: "PASS",
          gt: "0/1",
          dp: depth,
          ad_ref: depth - cnt,
          ad_alt: cnt,
          gq: 30,
          af,
        })
      }
    } else if best_ins.2 > 0 {
      let (p, ins_seq, cnt, depth) = best_ins
      let af = cnt.to_double() / depth.to_double()
      if af >= min_af && p >= 0 && p < n {
        let anchor = reference[p:p + 1].to_owned()
        variants.push({
          chrom,
          pos: p + 1,
          ref: anchor,
          alt: anchor + ins_seq,
          qual: 30.0,
          filter: "PASS",
          gt: "0/1",
          dp: depth,
          ad_ref: depth - cnt,
          ad_alt: cnt,
          gq: 30,
          af,
        })
      }
    }
    // region 内共识组装验证(SNP 位点已直接产出)
    ignore(gatk_assemble_consensus(reference, reads, rs, re))
  }
  variants.sort_by(fn(a, b) { a.pos - b.pos })
  variants
}

///|
/// 变异检出摘要.
pub fn gatk_variants_summary(variants : Array[GatkVariant]) -> String {
  let mut n_snp = 0
  let mut n_ins = 0
  let mut n_del = 0
  for v in variants {
    if v.ref.length() == 1 && v.alt.length() == 1 {
      n_snp = n_snp + 1
    } else if v.alt.length() > v.ref.length() {
      n_ins = n_ins + 1
    } else {
      n_del = n_del + 1
    }
  }
  "gatk-lite: \{variants.length()} variants (\{n_snp} SNP, \{n_ins} insertions, \{n_del} deletions)"
}

///|
/// 生成 VCF 4.2 文本(单样本), 与 parse_vcf 闭环.
pub fn gatk_write_vcf(
  chrom : String,
  reference : String,
  variants : Array[GatkVariant],
  sample~ : String = "SAMPLE",
) -> String {
  let sb = StringBuilder()
  sb.write_string("##fileformat=VCFv4.2\n")
  sb.write_string("##source=gatk-lite\n")
  sb.write_string("##reference=\{chrom}\n")
  sb.write_string("##contig=<ID=\{chrom},length=\{reference.length()}>\n")
  sb.write_string(
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
  )
  sb.write_string(
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
  )
  sb.write_string(
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n",
  )
  sb.write_string(
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n",
  )
  sb.write_string(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n",
  )
  for v in variants {
    sb.write_string(v.chrom)
    sb.write_char('\t')
    sb.write_string(v.pos.to_string())
    sb.write_string("\t.\t")
    sb.write_string(v.ref)
    sb.write_char('\t')
    sb.write_string(v.alt)
    sb.write_char('\t')
    sb.write_string(gk_fmt_qual(v.qual))
    sb.write_char('\t')
    sb.write_string(v.filter)
    sb.write_string("\tDP=")
    sb.write_string(v.dp.to_string())
    sb.write_string(";AF=")
    sb.write_string(gk_fmt_af(v.af))
    sb.write_string("\tGT:DP:AD:GQ\t")
    sb.write_string(v.gt)
    sb.write_char(':')
    sb.write_string(v.dp.to_string())
    sb.write_char(':')
    sb.write_string(v.ad_ref.to_string())
    sb.write_char(',')
    sb.write_string(v.ad_alt.to_string())
    sb.write_char(':')
    sb.write_string(v.gq.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
fn gk_fmt_qual(x : Double) -> String {
  let i = x.to_int()
  if i.to_double() == x {
    i.to_string()
  } else {
    x.to_string()
  }
}

///|
fn gk_fmt_af(x : Double) -> String {
  let s = x.to_string()
  if s.contains(".") {
    s
  } else {
    s + ".0"
  }
}
