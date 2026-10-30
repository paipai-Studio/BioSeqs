///|
/// mutect_lite.mbt — Mutect2 风格肿瘤-正常体细胞变异检出(本地确定性实现).
///
/// 简化 Mutect2 / Strelka2 的 LOD 统计框架:
/// - tumor/normal 双 pileup(复用 gatk_pileup_site 的碱基与 indel 事件计数);
/// - tumor LOD = log10 P(k_t | n_t, f_hat) - log10 P(k_t | n_t, p_err),
///   备择假设为肿瘤携带 alt 等位(f_hat = k/n 二项式 MLE),
///   原假设为 alt 读数全由测序错误产生(p_err = epsilon/3);
/// - normal LOD: normal 中 alt 证据显著 -> germline 过滤;
/// - 人群 AF 过滤(pop_af 表, gnomAD 风格, panel of normals 语义);
/// - TMB(每兆碱基 PASS somatic 密度)与 SNV/indel 分类统计;
/// - somatic VCF 4.2 输出(tumor/normal 双样本列, TLOD/NLOD INFO 字段);
/// - PASS SNV -> SomaticMutation 衔接 mutational_patterns 的 96 通道
///   突变签名分析(注意: 该模块以 chromosome 字段作样本分组键).
///
/// 数据不出本机.

// ---------------------------------------------------------------------------
// 数据结构
// ---------------------------------------------------------------------------

///|
/// 一个样本(肿瘤或正常)的读段集合.
pub(all) struct MutectSample {
  name : String
  reads : Array[GatkRead]
}

///|
/// 检出参数(默认值对齐 Mutect2 常用量).
pub(all) struct MutectParams {
  /// tumor 中 alt 最少读数.
  min_alt_count : Int32
  /// tumor 最小 alt 频率.
  min_tumor_af : Double
  /// tumor LOD 阈值(Mutect2 经典 5.3).
  lod_t_threshold : Double
  /// normal LOD 阈值, 超过判 germline.
  lod_n_threshold : Double
  /// 人群 AF 过滤阈值(gnomAD 风格).
  pop_af_max : Double
  /// 测序错误率.
  epsilon : Double
}

///|
/// 默认参数.
pub fn mutect_params() -> MutectParams {
  {
    min_alt_count: 3,
    min_tumor_af: 0.02,
    lod_t_threshold: 5.3,
    lod_n_threshold: 2.2,
    pop_af_max: 0.001,
    epsilon: 0.02,
  }
}

///|
/// 一条体细胞候选(含被过滤的, 以 filter 标签区分).
pub(all) struct MutectSomaticVariant {
  chrom : String
  /// 1-based VCF 坐标.
  pos : Int32
  ref_allele : String
  alt : String
  tumor_ref : Int32
  tumor_alt : Int32
  normal_ref : Int32
  normal_alt : Int32
  tumor_af : Double
  normal_af : Double
  /// tumor log10 odds.
  lod_t : Double
  /// normal log10 odds(germline 证据).
  lod_n : Double
  /// PASS / germline / weak_evidence / pop_af.
  filter : String
}

///|
/// 汇总统计(含 TMB).
pub(all) struct MutectSummary {
  /// 全部候选位点数(含被过滤的).
  n_candidates : Int32
  n_pass : Int32
  n_germline : Int32
  n_weak_evidence : Int32
  n_pop_af : Int32
  /// PASS 中的 SNV 数.
  n_snv : Int32
  /// PASS 中的 indel 数.
  n_indel : Int32
  /// 每兆碱基 callable 区域的 PASS somatic 密度.
  tmb_per_mb : Double
}

// ---------------------------------------------------------------------------
// 统计核心
// ---------------------------------------------------------------------------

///|
/// 二项式 log PMF(对数域, lgamma 避免组合数溢出).
fn mutect_binom_log_pmf(k : Int32, n : Int32, p : Double) -> Double {
  if n <= 0 {
    return 0.0
  }
  if k < 0 || k > n {
    return -1.0e30
  }
  if p <= 0.0 {
    return if k == 0 { 0.0 } else { -1.0e30 }
  }
  if p >= 1.0 {
    return if k == n { 0.0 } else { -1.0e30 }
  }
  lgamma(n.to_double() + 1.0) -
  lgamma(k.to_double() + 1.0) -
  lgamma((n - k).to_double() + 1.0) +
  k.to_double() * @math.ln(p) +
  (n - k).to_double() * @math.ln(1.0 - p)
}

///|
/// LOD = log10 P(k|n, f_hat) - log10 P(k|n, p_err), f_hat = k/n 的 MLE.
/// k=0 或 n=0 时无 alt 证据, LOD 为 0.
fn mutect_lod(k : Int32, n : Int32, epsilon : Double) -> Double {
  if n <= 0 || k <= 0 {
    return 0.0
  }
  let f_hat = k.to_double() / n.to_double()
  let p_err = epsilon / 3.0
  (mutect_binom_log_pmf(k, n, f_hat) - mutect_binom_log_pmf(k, n, p_err)) /
  @math.ln(10.0)
}

///|
/// pileup 中指定碱基的计数.
fn mutect_base_count(site : GatkPileupSite, base : Char) -> Int32 {
  for pair in site.base_counts {
    if pair.0 == base {
      return pair.1
    }
  }
  0
}

///|
/// pileup 中与指定序列相同的插入事件计数.
fn mutect_insertion_count(site : GatkPileupSite, seq : String) -> Int32 {
  for pair in site.insertions {
    if pair.0 == seq {
      return pair.1
    }
  }
  0
}

///|
/// pileup 中指定长度的缺失事件计数.
fn mutect_deletion_count(site : GatkPileupSite, len : Int32) -> Int32 {
  for pair in site.deletions {
    if pair.0 == len {
      return pair.1
    }
  }
  0
}

///|
fn mutect_is_acgt(b : Char) -> Bool {
  b == 'A' || b == 'C' || b == 'G' || b == 'T'
}

///|
/// 构造单条候选并打过滤标签: pop_af -> germline -> weak_evidence -> PASS.
fn mutect_make_variant(
  chrom : String,
  pos0 : Int32,
  ref_allele : String,
  alt_allele : String,
  t_ref : Int32,
  t_alt : Int32,
  n_depth : Int32,
  n_alt : Int32,
  params : MutectParams,
  pop_af : Map[String, Double],
) -> MutectSomaticVariant {
  let n_ref : Int32 = if n_depth > n_alt { n_depth - n_alt } else { 0 }
  let t_depth = t_ref + t_alt
  let tumor_af : Double = if t_depth > 0 {
    t_alt.to_double() / t_depth.to_double()
  } else {
    0.0
  }
  let normal_af : Double = if n_depth > 0 {
    n_alt.to_double() / n_depth.to_double()
  } else {
    0.0
  }
  let lod_t = mutect_lod(t_alt, t_depth, params.epsilon)
  let lod_n = mutect_lod(n_alt, n_depth, params.epsilon)
  let key = "\{chrom}:\{pos0 + 1}:\{ref_allele}>\{alt_allele}"
  let pop = pop_af.get(key).unwrap_or(0.0)
  let filter = if pop > params.pop_af_max {
    "pop_af"
  } else if lod_n > params.lod_n_threshold {
    "germline"
  } else if lod_t < params.lod_t_threshold || tumor_af < params.min_tumor_af {
    "weak_evidence"
  } else {
    "PASS"
  }
  {
    chrom,
    pos: pos0 + 1,
    ref_allele,
    alt: alt_allele,
    tumor_ref: t_ref,
    tumor_alt: t_alt,
    normal_ref: n_ref,
    normal_alt: n_alt,
    tumor_af,
    normal_af,
    lod_t,
    lod_n,
    filter,
  }
}

// ---------------------------------------------------------------------------
// 主调用
// ---------------------------------------------------------------------------

///|
/// 肿瘤-正常体细胞变异检出主入口.
///
/// 逐参考位点对 tumor 做 pileup, 碱基错配/插入/缺失事件计数达到
/// min_alt_count 的进入候选; normal 同点 pileup 提供 germline 证据;
/// pop_af 表键格式为 "chrom:pos:ref>alt"(1-based).
/// normal 读段为空时即为 tumor-only 模式(lod_n 恒 0).
pub fn mutect_call(
  chrom : String,
  reference : String,
  tumor : MutectSample,
  normal : MutectSample,
  params~ : MutectParams = mutect_params(),
  pop_af~ : Map[String, Double] = Map([], capacity=0),
) -> Array[MutectSomaticVariant] {
  let out : Array[MutectSomaticVariant] = Array(capacity=0)
  let n_pos = reference.length()
  for p = 0; p < n_pos; p = p + 1 {
    let t_site = gatk_pileup_site(reference, tumor.reads, p)
    if t_site.depth <= 0 {
      continue
    }
    let n_site = gatk_pileup_site(reference, normal.reads, p)
    // SNV 候选: 非参考 ACGT 碱基
    for pair in t_site.base_counts {
      let (base, cnt) = pair
      if base == t_site.ref_base || mutect_is_acgt(base) == false {
        continue
      }
      if cnt < params.min_alt_count {
        continue
      }
      let ref_allele = t_site.ref_base.to_string()
      let alt_allele = base.to_string()
      out.push(
        mutect_make_variant(
          chrom,
          p,
          ref_allele,
          alt_allele,
          t_site.depth - cnt,
          cnt,
          n_site.depth,
          mutect_base_count(n_site, base),
          params,
          pop_af,
        ),
      )
    }
    // 插入候选: 锚定碱基 + 插入序列
    for ipair in t_site.insertions {
      let (ins_seq, cnt) = ipair
      if cnt < params.min_alt_count {
        continue
      }
      let ref_allele = t_site.ref_base.to_string()
      let alt_allele = ref_allele + ins_seq
      out.push(
        mutect_make_variant(
          chrom,
          p,
          ref_allele,
          alt_allele,
          t_site.depth - cnt,
          cnt,
          n_site.depth,
          mutect_insertion_count(n_site, ins_seq),
          params,
          pop_af,
        ),
      )
    }
    // 缺失候选: 锚定碱基 + 被删碱基
    for dpair in t_site.deletions {
      let (dl, cnt) = dpair
      if cnt < params.min_alt_count {
        continue
      }
      let end = p + 1 + dl
      if end >= n_pos {
        continue
      }
      let ref_allele = reference[p:end].to_owned()
      let alt_allele = t_site.ref_base.to_string()
      out.push(
        mutect_make_variant(
          chrom,
          p,
          ref_allele,
          alt_allele,
          t_site.depth - cnt,
          cnt,
          n_site.depth,
          mutect_deletion_count(n_site, dl),
          params,
          pop_af,
        ),
      )
    }
  }
  mutect_dedup_snv_artifacts(out, reference)
}

///|
/// SNV 优先清洗: gatk_classify_indel 的简化 indel 推断会把单碱基替换
/// 读段解释为相邻锚点处的伪插入(插入序列 = SNV alt 碱基 + 参考延续).
/// 判据: 插入候选 v(锚点 a, 插入序列 ins)与 SNV 候选 s(位点 a+1)满足
/// ins[0] == s.alt 且 ins[1:] == reference[a+2 : a+1+ins.length] 时,
/// v 是 s 的伪影, 移除之. 真插入不满足该参考子串判据, 不受影响.
fn mutect_dedup_snv_artifacts(
  variants : Array[MutectSomaticVariant],
  reference : String,
) -> Array[MutectSomaticVariant] {
  // 收集 SNV 候选的 1-based pos -> alt 碱基
  let snv_alts : Map[Int32, Char] = Map([], capacity=0)
  for v in variants {
    if v.ref_allele.length() == 1 && v.alt.length() == 1 {
      snv_alts.set(v.pos, v.alt.unsafe_get(0).unsafe_to_char())
    }
  }
  if snv_alts.length() == 0 {
    return variants
  }
  let out : Array[MutectSomaticVariant] = Array(capacity=variants.length())
  for v in variants {
    // 仅检查插入候选(alt = 锚定碱基 + 插入序列)
    let is_insertion = v.ref_allele.length() == 1 &&
      v.alt.length() > v.ref_allele.length()
    if is_insertion == false {
      out.push(v)
      continue
    }
    // 伪影判据: 锚点下一位置存在 SNV 候选, 且插入序列解释为
    // SNV alt + reference[a+2 ..]
    let anchor = v.pos // 1-based 锚点
    let ins = v.alt[1:].to_owned()
    let artifact = match snv_alts.get(anchor + 1) {
      Some(alt_char) =>
        if ins.length() == 0 {
          false
        } else {
          let head_match = ins.unsafe_get(0).unsafe_to_char() == alt_char
          // reference 0-based: 锚点 a0 = anchor-1, SNV 位点 p0 = anchor,
          // 插入序列后续应对应 reference[p0+1 : p0+ins.length]
          let start = anchor + 1
          let stop = anchor + ins.length()
          let tail_match = if stop <= reference.length() {
            ins[1:].to_owned() == reference[start:stop].to_owned()
          } else {
            false
          }
          head_match && tail_match
        }
      None => false
    }
    if artifact == false {
      out.push(v)
    }
  }
  out
}

///|
/// 仅保留 PASS 的体细胞变异.
pub fn mutect_pass(
  variants : Array[MutectSomaticVariant],
) -> Array[MutectSomaticVariant] {
  let out : Array[MutectSomaticVariant] = Array(capacity=variants.length())
  for v in variants {
    if v.filter == "PASS" {
      out.push(v)
    }
  }
  out
}

///|
/// 汇总统计与 TMB(callable_mb 为可调用区域兆碱基数).
pub fn mutect_summary(
  variants : Array[MutectSomaticVariant],
  callable_mb~ : Double = 1.0,
) -> MutectSummary {
  let mut n_pass = 0
  let mut n_germline = 0
  let mut n_weak = 0
  let mut n_pop = 0
  let mut n_snv = 0
  let mut n_indel = 0
  for v in variants {
    match v.filter {
      "PASS" => {
        n_pass = n_pass + 1
        if v.ref_allele.length() == 1 && v.alt.length() == 1 {
          n_snv = n_snv + 1
        } else {
          n_indel = n_indel + 1
        }
      }
      "germline" => n_germline = n_germline + 1
      "weak_evidence" => n_weak = n_weak + 1
      "pop_af" => n_pop = n_pop + 1
      _ => ()
    }
  }
  let tmb : Double = if callable_mb > 0.0 {
    n_pass.to_double() / callable_mb
  } else {
    0.0
  }
  {
    n_candidates: variants.length(),
    n_pass,
    n_germline,
    n_weak_evidence: n_weak,
    n_pop_af: n_pop,
    n_snv,
    n_indel,
    tmb_per_mb: tmb,
  }
}

// ---------------------------------------------------------------------------
// VCF 输出
// ---------------------------------------------------------------------------

///|
fn mutect_gt(af : Double, alt_count : Int32) -> String {
  if alt_count <= 0 {
    "0/0"
  } else if af >= 0.75 {
    "1/1"
  } else {
    "0/1"
  }
}

///|
/// 体细胞候选 -> somatic VCF 4.2 文本(tumor/normal 双样本列).
pub fn mutect_write_vcf(
  variants : Array[MutectSomaticVariant],
  tumor_name : String,
  normal_name : String,
) -> String {
  let sb = StringBuilder()
  sb.write_string("##fileformat=VCFv4.2\n")
  sb.write_string("##source=mutect_lite\n")
  sb.write_string(
    "##FILTER=<ID=germline,Description=\"Normal sample shows significant alt evidence (LOD_N above threshold)\">\n",
  )
  sb.write_string(
    "##FILTER=<ID=weak_evidence,Description=\"Tumor LOD below threshold or alt fraction too low\">\n",
  )
  sb.write_string(
    "##FILTER=<ID=pop_af,Description=\"Population allele frequency above threshold\">\n",
  )
  sb.write_string(
    "##INFO=<ID=TLOD,Number=1,Type=Float,Description=\"Log10 odds of variant in tumor\">\n",
  )
  sb.write_string(
    "##INFO=<ID=NLOD,Number=1,Type=Float,Description=\"Log10 odds of variant in normal (germline evidence)\">\n",
  )
  sb.write_string("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
  sb.write_string(tumor_name)
  sb.write_char('\t')
  sb.write_string(normal_name)
  sb.write_char('\n')
  for v in variants {
    sb.write_string(v.chrom)
    sb.write_char('\t')
    sb.write_string(v.pos.to_string())
    sb.write_string("\t.\t")
    sb.write_string(v.ref_allele)
    sb.write_char('\t')
    sb.write_string(v.alt)
    sb.write_char('\t')
    sb.write_string((v.lod_t * 10.0).to_string())
    sb.write_char('\t')
    sb.write_string(v.filter)
    sb.write_string("\tTLOD=")
    sb.write_string(v.lod_t.to_string())
    sb.write_string(";NLOD=")
    sb.write_string(v.lod_n.to_string())
    sb.write_string("\tGT:AD:AF\t")
    sb.write_string(mutect_gt(v.tumor_af, v.tumor_alt))
    sb.write_char(':')
    sb.write_string(v.tumor_ref.to_string())
    sb.write_char(',')
    sb.write_string(v.tumor_alt.to_string())
    sb.write_char(':')
    sb.write_string(v.tumor_af.to_string())
    sb.write_char('\t')
    sb.write_string(mutect_gt(v.normal_af, v.normal_alt))
    sb.write_char(':')
    sb.write_string(v.normal_ref.to_string())
    sb.write_char(',')
    sb.write_string(v.normal_alt.to_string())
    sb.write_char(':')
    sb.write_string(v.normal_af.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}

// ---------------------------------------------------------------------------
// 衔接 mutational_patterns
// ---------------------------------------------------------------------------

///|
/// PASS SNV -> SomaticMutation(供 mutational_patterns 的 96 通道
/// 突变签名分析). trinucleotide 上下文取参考序列 pos-1..pos+1,
/// 参考边界外以 N 填充.
/// 注意: mutational_patterns 以 chromosome 字段作样本分组键,
/// 故此处 chromosome 填 sample_name 而非真实染色体.
pub fn mutect_to_somatic_mutations(
  reference : String,
  variants : Array[MutectSomaticVariant],
  sample_name~ : String = "tumor",
) -> Array[SomaticMutation] {
  let out : Array[SomaticMutation] = Array(capacity=0)
  for v in variants {
    if v.filter != "PASS" {
      continue
    }
    if v.ref_allele.length() != 1 || v.alt.length() != 1 {
      continue
    }
    let p = v.pos - 1
    let left : Char = if p > 0 {
      reference.unsafe_get(p - 1).unsafe_to_char()
    } else {
      'N'
    }
    let right : Char = if p + 1 < reference.length() {
      reference.unsafe_get(p + 1).unsafe_to_char()
    } else {
      'N'
    }
    let ctx = left.to_string() + v.ref_allele + right.to_string()
    out.push(SomaticMutation::new(sample_name, v.pos, v.ref_allele, v.alt, ctx))
  }
  out
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// 体细胞候选 -> JSON 值.
pub fn mutect_value(variants : Array[MutectSomaticVariant]) -> JsonValue {
  let arr : Array[JsonValue] = Array(capacity=variants.length())
  for v in variants {
    arr.push(
      JsonObj(
        Map([
          ("chrom", JsonStr(v.chrom)),
          ("pos", JsonNum(v.pos.to_double())),
          ("ref", JsonStr(v.ref_allele)),
          ("alt", JsonStr(v.alt)),
          ("tumor_ref", JsonNum(v.tumor_ref.to_double())),
          ("tumor_alt", JsonNum(v.tumor_alt.to_double())),
          ("normal_ref", JsonNum(v.normal_ref.to_double())),
          ("normal_alt", JsonNum(v.normal_alt.to_double())),
          ("tumor_af", JsonNum(v.tumor_af)),
          ("normal_af", JsonNum(v.normal_af)),
          ("lod_t", JsonNum(v.lod_t)),
          ("lod_n", JsonNum(v.lod_n)),
          ("filter", JsonStr(v.filter)),
        ]),
      ),
    )
  }
  JsonArr(arr)
}

///|
/// 体细胞候选 -> JSON 文本.
pub fn mutect_json(variants : Array[MutectSomaticVariant]) -> String {
  json_to_string(mutect_value(variants))
}

///|
/// 汇总摘要行.
pub fn mutect_summary_line(summary : MutectSummary) -> String {
  "mutect-lite: candidates=\{summary.n_candidates}, pass=\{summary.n_pass}, germline=\{summary.n_germline}, weak=\{summary.n_weak_evidence}, pop_af=\{summary.n_pop_af}, snv=\{summary.n_snv}, indel=\{summary.n_indel}, tmb=\{summary.tmb_per_mb} per Mb"
}
