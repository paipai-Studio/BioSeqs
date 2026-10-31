///|
/// gvcf.mbt — gVCF 参考块表示 + 多样本合并 + 联合基因分型(简化本地实现).
///
/// 对齐 GATK 的 GVCF 语义(本地确定性实现, 数据不出本机):
/// - 变异位点记录携带 ALT(含 <NON_REF>)与逐样本基因型似然(PL);
/// - 无变异区间压缩为参考块记录(POS..END, ALT=<NON_REF>, 纯合参考);
/// - 多样本合并: 取各样本变异位点并集, 未被变异记录覆盖但落在参考块内的
///   样本判为纯合参考(0/0), 完全不覆盖的样本记为缺失(./.);
/// - 联合基因分型: 由各样本 PL 重新选择最优基因型(0/0, 0/1, 1/1)并重算 GQ.
///
/// 输入来自 gatk_lite 的检出变异, 输出多样本 VCF 4.2 文本, 与现有 VCF
/// 后处理栈(bcftools_ops/variant_filtering)衔接.

///|
/// gVCF 处理错误.
pub suberror GvcfError {
  GvcfError(String)
}

///|
/// 单个样本在某个记录上的基因型(含基因型似然).
pub(all) struct GvcfGenotype {
  /// 基因型, 如 "0/0" / "0/1" / "1/1" / "./.".
  gt : String
  /// phred 标度基因型似然, 顺序为 [0/0, 0/1, 1/1], 最小值最优; 空表示缺失.
  pl : Array[Int32]
  /// 覆盖深度.
  dp : Int32
  /// 基因型质量(GQ, 次小 PL 减最小 PL, 上限 99).
  gq : Int32
}

///|
/// gVCF 记录(变异记录或参考块记录).
pub(all) struct GvcfRecord {
  chrom : String
  /// 1-based 位置.
  pos : Int32
  /// 参考等位.
  ref_allele : String
  /// 替代等位(参考块为 ["<NON_REF>"]).
  alts : Array[String]
  /// 是否为参考块记录.
  is_block : Bool
  /// 参考块结束位置(1-based, 含); 变异记录等于 pos + ref.length() - 1.
  block_end : Int32
  /// 样本 -> 基因型.
  samples : Map[String, GvcfGenotype]
}

///|
/// 单染色体 gVCF 调用集.
pub(all) struct GvcfCallset {
  chrom : String
  samples : Array[String]
  records : Array[GvcfRecord]
}

// ============================================================================
// 基因型似然
// ============================================================================

///|
fn gvcf_pow(base : Double, exponent : Double) -> Double {
  if exponent == 0.0 {
    1.0
  } else if base <= 0.0 {
    0.0
  } else {
    @math.exp(exponent * @math.ln(base))
  }
}

///|
fn gvcf_log10(x : Double) -> Double {
  if x <= 0.0 {
    -300.0
  } else {
    @math.ln(x) / 2.302585092994046
  }
}

///|
/// 由等位深度 AD 计算基因型似然(二项/误差模型: e=0.01 测序错误率).
///
/// P(0/0) 正比 (1-e)^ref · e^alt; P(0/1) 正比 0.5^(ref+alt);
/// P(1/1) 正比 e^ref · (1-e)^alt. 转 phred 后减去最小值归一.
pub fn gvcf_pl_from_ad(ad_ref : Int32, ad_alt : Int32) -> Array[Int32] {
  let e = 0.01
  let r = ad_ref.to_double()
  let a = ad_alt.to_double()
  let p00 = gvcf_pow(1.0 - e, r) * gvcf_pow(e, a)
  let p01 = gvcf_pow(0.5, r + a)
  let p11 = gvcf_pow(e, r) * gvcf_pow(1.0 - e, a)
  let raw : Array[Double] = [
    -10.0 * gvcf_log10(p00),
    -10.0 * gvcf_log10(p01),
    -10.0 * gvcf_log10(p11),
  ]
  let mut m = raw[0]
  for x in raw {
    if x < m {
      m = x
    }
  }
  let out : Array[Int32] = Array(capacity=3)
  for x in raw {
    let v = (x - m).round().to_int()
    out.push(if v < 0 { 0 } else { v })
  }
  out
}

///|
/// 依据基因型字符串给出兜底 PL(深度不足时使用).
fn gvcf_default_pl_for_gt(gt : String) -> Array[Int32] {
  if gt == "0/1" {
    [30, 0, 30]
  } else if gt == "1/1" {
    [90, 30, 0]
  } else {
    [0, 60, 90]
  }
}

///|
/// 由 PL 计算 GQ(次小减最小, 上限 99).
pub fn gvcf_gq_from_pl(pl : Array[Int32]) -> Int32 {
  if pl.length() < 2 {
    return 0
  }
  let sorted = pl.copy()
  sorted.sort()
  let gq = sorted[1] - sorted[0]
  if gq > 99 {
    99
  } else {
    gq
  }
}

///|
/// 由 PL 选择最优基因型(最小值索引).
///
/// 支持多等位: PL 长度满足 n=(m+1)(m+2)/2 时按 VCF 基因型排序
/// (0/0, 0/1, 1/1, 0/2, 1/2, 2/2, ...) 解析最小项;
/// 长度不规范返回 "./.".
pub fn gvcf_gt_from_pl(pl : Array[Int32]) -> String {
  let n = pl.length()
  if n < 3 {
    return "./."
  }
  let mut m = 1
  while (m + 1) * (m + 2) / 2 < n {
    m = m + 1
  }
  if (m + 1) * (m + 2) / 2 != n {
    return "./."
  }
  let mut best = 0
  for idx = 1; idx < n; idx = idx + 1 {
    if pl[idx] < pl[best] {
      best = idx
    }
  }
  let mut j = 0
  while (j + 1) * (j + 2) / 2 <= best {
    j = j + 1
  }
  let i = best - j * (j + 1) / 2
  i.to_string() + "/" + j.to_string()
}

// ============================================================================
// 建库: 由 gatk_lite 变异构建单样本 gVCF(变异记录 + 参考块)
// ============================================================================

///|
fn gvcf_make_block(
  chrom : String,
  reference : String,
  start : Int32,
  end : Int32,
  sample : String,
) -> GvcfRecord {
  let ref_base = reference[start - 1:start].to_owned()
  let geno = GvcfGenotype::{ gt: "0/0", pl: [0, 60, 90], dp: 0, gq: 60, }
  GvcfRecord::{
    chrom,
    pos: start,
    ref_allele: ref_base,
    alts: ["<NON_REF>"],
    is_block: true,
    block_end: end,
    samples: Map([(sample, geno)]),
  }
}

///|
/// 由 gatk_lite 变异构建单样本 gVCF.
///
/// 变异之间的无变异区间压缩为参考块(END 标记), 变异位点输出为普通记录.
pub fn gvcf_build_from_variants(
  chrom : String,
  reference : String,
  variants : Array[GatkVariant],
  sample~ : String = "SAMPLE",
) -> GvcfCallset raise GvcfError {
  if reference.length() == 0 {
    raise GvcfError("reference sequence is empty")
  }
  let sorted = variants.copy()
  sorted.sort_by(fn(a, b) { a.pos - b.pos })
  let records : Array[GvcfRecord] = Array(capacity=sorted.length() * 2 + 1)
  let mut cursor = 1
  for v in sorted {
    let start = v.pos
    if start > cursor {
      records.push(gvcf_make_block(chrom, reference, cursor, start - 1, sample))
    }
    let pl = if v.ad_ref + v.ad_alt > 0 {
      gvcf_pl_from_ad(v.ad_ref, v.ad_alt)
    } else {
      gvcf_default_pl_for_gt(v.gt)
    }
    let geno = GvcfGenotype::{
      gt: v.gt,
      pl,
      dp: v.dp,
      gq: gvcf_gq_from_pl(pl),
    }
    let span = if v.ref.length() > 0 { v.ref.length() } else { 1 }
    records.push(GvcfRecord::{
      chrom,
      pos: start,
      ref_allele: v.ref,
      alts: [v.alt],
      is_block: false,
      block_end: start + span - 1,
      samples: Map([(sample, geno)]),
    })
    cursor = start + span
  }
  if cursor <= reference.length() {
    records.push(
      gvcf_make_block(chrom, reference, cursor, reference.length(), sample),
    )
  }
  GvcfCallset::{ chrom, samples: [sample], records, }
}

// ============================================================================
// 解析
// ============================================================================

///|
fn gvcf_parse_int(s : String) -> Int32 {
  let mut value = 0
  let mut negative = false
  let mut i = 0
  if s.length() > 0 && s.unsafe_get(0).unsafe_to_char() == '-' {
    negative = true
    i = 1
  }
  while i < s.length() {
    let c = s.unsafe_get(i).unsafe_to_char()
    if c >= '0' && c <= '9' {
      value = value * 10 + (c.to_int() - '0'.to_int())
    }
    i = i + 1
  }
  if negative {
    -value
  } else {
    value
  }
}

///|
fn gvcf_parse_pl(s : String) -> Array[Int32] {
  let out : Array[Int32] = []
  if s == "." || s == "" {
    return out
  }
  for part in s.split(",").to_array() {
    out.push(gvcf_parse_int(part.to_owned()))
  }
  out
}

///|
fn gvcf_split_alts(alt_field : String) -> Array[String] {
  let out : Array[String] = []
  if alt_field == "." || alt_field == "" {
    return out
  }
  for part in alt_field.split(",").to_array() {
    out.push(part.to_owned())
  }
  out
}

///|
/// 解析 INFO 中的 END=(参考块结束位置), 缺省返回 pos.
fn gvcf_parse_end(info : String, pos : Int32) -> Int32 {
  for part in info.split(";").to_array() {
    let token = part.to_owned()
    if token.has_prefix("END=") {
      return gvcf_parse_int(token[4:token.length()].to_owned())
    }
  }
  pos
}

///|
fn gvcf_parse_sample(field : String, keys : Array[String]) -> GvcfGenotype {
  let vals = field.split(":").to_array()
  let mut gt = "./."
  let mut pl : Array[Int32] = []
  let mut dp = 0
  let mut gq = 0
  for i, key in keys {
    if i < vals.length() {
      let value = vals[i].to_owned()
      if key == "GT" {
        gt = value
      } else if key == "PL" {
        pl = gvcf_parse_pl(value)
      } else if key == "DP" {
        dp = gvcf_parse_int(value)
      } else if key == "GQ" {
        gq = gvcf_parse_int(value)
      }
    }
  }
  GvcfGenotype::{ gt, pl, dp, gq, }
}

///|
/// 解析 gVCF 文本(变异记录 + 参考块记录).
pub fn gvcf_parse(text : String) -> GvcfCallset? {
  let lines = text.split("\n").to_array()
  let samples : Array[String] = []
  let records : Array[GvcfRecord] = []
  let mut chrom = ""
  for raw in lines {
    let line = raw.to_owned()
    if line.length() == 0 {
      continue
    }
    if line.unsafe_get(0).unsafe_to_char() == '#' {
      if line.has_prefix("#CHROM") {
        let cols = line.split("\t").to_array()
        let mut i = 9
        while i < cols.length() {
          samples.push(cols[i].to_owned())
          i = i + 1
        }
      }
      continue
    }
    let cols = line.split("\t").to_array()
    if cols.length() < 8 {
      continue
    }
    let rec_chrom = cols[0].to_owned()
    let pos = gvcf_parse_int(cols[1].to_owned())
    let ref_allele = cols[3].to_owned()
    let alts = gvcf_split_alts(cols[4].to_owned())
    let info = cols[7].to_owned()
    let block_end = gvcf_parse_end(info, pos)
    let mut only_non_ref = alts.length() > 0
    for alt in alts {
      if alt != "<NON_REF>" {
        only_non_ref = false
      }
    }
    let is_block = only_non_ref
    let keys : Array[String] = if cols.length() > 8 {
      let fmt = cols[8].to_owned()
      let out : Array[String] = []
      for part in fmt.split(":").to_array() {
        out.push(part.to_owned())
      }
      out
    } else {
      []
    }
    let sample_map : Map[String, GvcfGenotype] = Map(
      [],
      capacity=samples.length(),
    )
    for si = 0; si < samples.length(); si = si + 1 {
      let col = 9 + si
      if col < cols.length() && cols[col].to_owned() != "." {
        sample_map.set(
          samples[si],
          gvcf_parse_sample(cols[col].to_owned(), keys),
        )
      }
    }
    if chrom == "" {
      chrom = rec_chrom
    }
    records.push(GvcfRecord::{
      chrom: rec_chrom,
      pos,
      ref_allele,
      alts,
      is_block,
      block_end,
      samples: sample_map,
    })
  }
  if chrom == "" {
    None
  } else {
    Some(GvcfCallset::{ chrom, samples, records, })
  }
}

// ============================================================================
// 多样本合并 + 联合基因分型
// ============================================================================

///|
fn gvcf_pos_in_blocks(blocks : Array[(Int32, Int32)], pos : Int32) -> Bool {
  for interval in blocks {
    if pos >= interval.0 && pos <= interval.1 {
      return true
    }
  }
  false
}

///|
fn gvcf_merge_str_set(
  target : Array[String],
  seen : Map[String, Bool],
  value : String,
) -> Unit {
  if !seen.contains(value) {
    seen.set(value, true)
    target.push(value)
  }
}

///|
/// 在 ALT 列表中查找等位基因索引(0-based); 未找到返回 -1.
fn gvcf_alt_index(alts : Array[String], alt : String) -> Int32 {
  for i = 0; i < alts.length(); i = i + 1 {
    if alts[i] == alt {
      return i
    }
  }
  -1
}

///|
/// 解析 GT 中的单个等位索引数字; 非数字(含 '.')返回 None.
fn gvcf_parse_allele(token : String) -> Int32? {
  if token.length() == 0 {
    return None
  }
  let mut value = 0
  for i = 0; i < token.length(); i = i + 1 {
    let c = token.unsafe_get(i).unsafe_to_char()
    if c < '0' || c > '9' {
      return None
    }
    value = value * 10 + (c.to_int() - 48)
  }
  Some(value)
}

///|
/// 重映射单个等位索引: 0 保持 REF; a>0 映射到原 ALT 在合并 ALT 列表中的新索引.
fn gvcf_remap_allele(
  allele : Int32,
  orig_alts : Array[String],
  merged_alts : Array[String],
) -> Int32? {
  if allele == 0 {
    return Some(0)
  }
  if allele < 0 || allele > orig_alts.length() {
    return None
  }
  let idx = gvcf_alt_index(merged_alts, orig_alts[allele - 1])
  if idx < 0 {
    return None
  }
  Some(idx + 1)
}

///|
/// 重映射 GT 等位索引到合并后的 ALT 顺序(升序输出; 无法映射的等位记缺失).
fn gvcf_remap_gt(
  gt : String,
  orig_alts : Array[String],
  merged_alts : Array[String],
) -> String {
  let normalized = gt.replace_all(old="/", new="|")
  let parts = normalized.split("|").collect()
  if parts.length() != 2 {
    return gt
  }
  let alleles : Array[Int32] = Array(capacity=2)
  for part in parts {
    match gvcf_parse_allele(part.to_owned()) {
      Some(a) =>
        match gvcf_remap_allele(a, orig_alts, merged_alts) {
          Some(mapped) => alleles.push(mapped)
          None => ()
        }
      None => ()
    }
  }
  if alleles.length() != 2 {
    return "./."
  }
  let (a, b) = if alleles[0] <= alleles[1] {
    (alleles[0], alleles[1])
  } else {
    (alleles[1], alleles[0])
  }
  a.to_string() + "/" + b.to_string()
}

///|
/// 重映射 PL 到合并后的 ALT 顺序.
///
/// 仅当原 ALT 集合与合并 ALT 集合大小一致(双射重排)时输出完整 PL;
/// 否则返回空数组 —— 等位数变化时基因型似然空间无法一一对应,
/// 按 VCF 规范将整个 PL 字段置为缺失, 由下游保留重映射后的 GT.
fn gvcf_remap_pl(
  pl : Array[Int32],
  orig_alts : Array[String],
  merged_alts : Array[String],
) -> Array[Int32] {
  let m = orig_alts.length()
  let size = (m + 1) * (m + 2) / 2
  if m != merged_alts.length() || pl.length() != size {
    return []
  }
  let out : Array[Int32] = Array::make(size, 0)
  let filled : Array[Bool] = Array::make(size, false)
  for j = 0; j <= m; j = j + 1 {
    for i = 0; i <= j; i = i + 1 {
      let idx = j * (j + 1) / 2 + i
      match
        (
          gvcf_remap_allele(i, orig_alts, merged_alts),
          gvcf_remap_allele(j, orig_alts, merged_alts),
        ) {
        (Some(x), Some(y)) => {
          let (a, b) = if x <= y { (x, y) } else { (y, x) }
          let nidx = b * (b + 1) / 2 + a
          out[nidx] = pl[idx]
          filled[nidx] = true
        }
        _ => ()
      }
    }
  }
  for i = 0; i < size; i = i + 1 {
    if !filled[i] {
      return []
    }
  }
  out
}

///|
/// 合并多个(单样本)gVCF.
///
/// 取各样本变异位点的并集; 变异记录之间被参考块覆盖的样本记为纯合参考
/// (0/0, PL=[0,60,90]); 既非变异记录又无参考块覆盖的样本记缺失(./.).
/// 各样本 GT/PL 按其原 callset 的 ALT 顺序解释, 输出前重映射到
/// 合并后(排序)的 ALT 顺序, 避免同位点不同 ALT 时等位索引错位.
pub fn gvcf_merge(callsets : Array[GvcfCallset]) -> GvcfCallset raise GvcfError {
  if callsets.length() == 0 {
    raise GvcfError("no callsets to merge")
  }
  let chrom = callsets[0].chrom
  let samples : Array[String] = []
  let sample_seen : Map[String, Bool] = Map([], capacity=0)
  for cs in callsets {
    for s in cs.samples {
      gvcf_merge_str_set(samples, sample_seen, s)
    }
  }
  samples.sort()

  let sample_calls : Map[String, Map[Int32, GvcfGenotype]] = Map(
    [],
    capacity=samples.length(),
  )
  let sample_alts : Map[String, Map[Int32, Array[String]]] = Map(
    [],
    capacity=samples.length(),
  )
  let sample_blocks : Map[String, Array[(Int32, Int32)]] = Map(
    [],
    capacity=samples.length(),
  )
  for s in samples {
    sample_calls.set(s, Map([], capacity=0))
    sample_alts.set(s, Map([], capacity=0))
    sample_blocks.set(s, [])
  }
  let positions : Array[Int32] = []
  let pos_seen : Map[Int32, Bool] = Map([], capacity=0)
  let pos_ref : Map[Int32, String] = Map([], capacity=0)
  let pos_alts : Map[Int32, Array[String]] = Map([], capacity=0)
  let alt_seen : Map[Int32, Map[String, Bool]] = Map([], capacity=0)

  for cs in callsets {
    for record in cs.records {
      for s in cs.samples {
        if record.is_block {
          sample_blocks.get(s).unwrap().push((record.pos, record.block_end))
        } else {
          match record.samples.get(s) {
            Some(geno) => {
              sample_calls.get(s).unwrap().set(record.pos, geno)
              let record_alts : Array[String] = []
              for alt in record.alts {
                if alt != "<NON_REF>" {
                  record_alts.push(alt)
                }
              }
              sample_alts.get(s).unwrap().set(record.pos, record_alts)
            }
            None => ()
          }
        }
      }
      if !record.is_block {
        if !pos_seen.contains(record.pos) {
          pos_seen.set(record.pos, true)
          pos_ref.set(record.pos, record.ref_allele)
          pos_alts.set(record.pos, [])
          alt_seen.set(record.pos, Map([], capacity=0))
          positions.push(record.pos)
        }
        let alts = pos_alts.get(record.pos).unwrap()
        let seen = alt_seen.get(record.pos).unwrap()
        for alt in record.alts {
          if alt != "<NON_REF>" {
            gvcf_merge_str_set(alts, seen, alt)
          }
        }
      }
    }
  }
  positions.sort()

  let records : Array[GvcfRecord] = Array(capacity=positions.length())
  for pos in positions {
    let ref_allele = pos_ref.get(pos).unwrap()
    let alts = pos_alts.get(pos).unwrap()
    alts.sort()
    let sample_map : Map[String, GvcfGenotype] = Map(
      [],
      capacity=samples.length(),
    )
    for s in samples {
      let calls = sample_calls.get(s).unwrap()
      match calls.get(pos) {
        Some(geno) => {
          let orig_alts = match sample_alts.get(s).unwrap().get(pos) {
            Some(orig) => orig
            None => []
          }
          sample_map.set(s, GvcfGenotype::{
            gt: gvcf_remap_gt(geno.gt, orig_alts, alts),
            pl: gvcf_remap_pl(geno.pl, orig_alts, alts),
            dp: geno.dp,
            gq: geno.gq,
          })
        }
        None =>
          if gvcf_pos_in_blocks(sample_blocks.get(s).unwrap(), pos) {
            sample_map.set(s, GvcfGenotype::{
              gt: "0/0",
              pl: [0, 60, 90],
              dp: 0,
              gq: 60,
            })
          } else {
            sample_map.set(s, GvcfGenotype::{ gt: "./.", pl: [], dp: 0, gq: 0, })
          }
      }
    }
    let span = if ref_allele.length() > 0 { ref_allele.length() } else { 1 }
    records.push(GvcfRecord::{
      chrom,
      pos,
      ref_allele,
      alts,
      is_block: false,
      block_end: pos + span - 1,
      samples: sample_map,
    })
  }
  GvcfCallset::{ chrom, samples, records, }
}

///|
/// 联合基因分型: 对每个变异位点按各样本 PL 重新选择最优基因型并重算 GQ.
///
/// 缺失样本(PL 为空)保持 "./.". 返回新的调用集, 不修改输入.
pub fn gvcf_joint_genotype(callset : GvcfCallset) -> GvcfCallset {
  let records : Array[GvcfRecord] = Array(capacity=callset.records.length())
  for record in callset.records {
    let sample_map : Map[String, GvcfGenotype] = Map(
      [],
      capacity=record.samples.length(),
    )
    for s in callset.samples {
      match record.samples.get(s) {
        Some(geno) =>
          if geno.pl.length() >= 3 {
            sample_map.set(s, GvcfGenotype::{
              gt: gvcf_gt_from_pl(geno.pl),
              pl: geno.pl.copy(),
              dp: geno.dp,
              gq: gvcf_gq_from_pl(geno.pl),
            })
          } else {
            sample_map.set(s, geno)
          }
        None => ()
      }
    }
    records.push(GvcfRecord::{
      chrom: record.chrom,
      pos: record.pos,
      ref_allele: record.ref_allele,
      alts: record.alts.copy(),
      is_block: false,
      block_end: record.block_end,
      samples: sample_map,
    })
  }
  { chrom: callset.chrom, samples: callset.samples.copy(), records, }
}

// ============================================================================
// VCF 输出
// ============================================================================

///|
fn gvcf_join_ints(values : Array[Int32]) -> String {
  let sb = StringBuilder()
  for i, value in values {
    if i > 0 {
      sb.write_char(',')
    }
    sb.write_string(value.to_string())
  }
  sb.to_string()
}

///|
/// 调用集 -> 多样本 VCF 4.2 文本.
pub fn gvcf_write_vcf(callset : GvcfCallset) -> String {
  let sb = StringBuilder()
  sb.write_string("##fileformat=VCFv4.2\n")
  sb.write_string("##source=gvcf-lite\n")
  sb.write_string("##contig=<ID=\{callset.chrom}>\n")
  sb.write_string(
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
  )
  sb.write_string(
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n",
  )
  sb.write_string(
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n",
  )
  sb.write_string(
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype likelihoods\">\n",
  )
  sb.write_string("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
  for s in callset.samples {
    sb.write_char('\t')
    sb.write_string(s)
  }
  sb.write_char('\n')
  let sorted = callset.records.copy()
  sorted.sort_by(fn(a, b) { a.pos - b.pos })
  for record in sorted {
    sb.write_string(record.chrom)
    sb.write_char('\t')
    sb.write_string(record.pos.to_string())
    sb.write_string("\t.\t")
    sb.write_string(record.ref_allele)
    sb.write_char('\t')
    if record.alts.length() == 0 {
      sb.write_string(".")
    } else {
      for i, alt in record.alts {
        if i > 0 {
          sb.write_char(',')
        }
        sb.write_string(alt)
      }
    }
    sb.write_string("\t.\t.\t.\tGT:DP:GQ:PL")
    for s in callset.samples {
      sb.write_char('\t')
      match record.samples.get(s) {
        Some(geno) => {
          sb.write_string(geno.gt)
          sb.write_char(':')
          sb.write_string(geno.dp.to_string())
          sb.write_char(':')
          sb.write_string(geno.gq.to_string())
          sb.write_char(':')
          if geno.pl.length() == 0 {
            sb.write_char('.')
          } else {
            sb.write_string(gvcf_join_ints(geno.pl))
          }
        }
        None => sb.write_string("./.:0:0:.")
      }
    }
    sb.write_char('\n')
  }
  sb.to_string()
}

// ============================================================================
// JSON 序列化 / 反序列化
// ============================================================================

///|
fn gvcf_genotype_value(geno : GvcfGenotype) -> JsonValue {
  let pls : Array[JsonValue] = Array(capacity=geno.pl.length())
  for value in geno.pl {
    pls.push(JsonNum(value.to_double()))
  }
  JsonObj(
    Map([
      ("gt", JsonStr(geno.gt)),
      ("pl", JsonArr(pls)),
      ("dp", JsonNum(geno.dp.to_double())),
      ("gq", JsonNum(geno.gq.to_double())),
    ]),
  )
}

///|
fn gvcf_record_value(record : GvcfRecord) -> JsonValue {
  let alts : Array[JsonValue] = Array(capacity=record.alts.length())
  for alt in record.alts {
    alts.push(JsonStr(alt))
  }
  let names : Array[String] = []
  for s in record.samples.keys() {
    names.push(s)
  }
  names.sort()
  let pairs : Array[(String, JsonValue)] = Array(capacity=names.length())
  for s in names {
    pairs.push((s, gvcf_genotype_value(record.samples.get(s).unwrap())))
  }
  JsonObj(
    Map([
      ("chrom", JsonStr(record.chrom)),
      ("pos", JsonNum(record.pos.to_double())),
      ("ref", JsonStr(record.ref_allele)),
      ("alts", JsonArr(alts)),
      ("is_block", JsonBool(record.is_block)),
      ("block_end", JsonNum(record.block_end.to_double())),
      ("samples", JsonObj(Map(pairs))),
    ]),
  )
}

///|
/// 调用集 -> JSON 值.
pub fn gvcf_callset_value(callset : GvcfCallset) -> JsonValue {
  let records : Array[JsonValue] = Array(capacity=callset.records.length())
  for record in callset.records {
    records.push(gvcf_record_value(record))
  }
  let samples : Array[JsonValue] = Array(capacity=callset.samples.length())
  for s in callset.samples {
    samples.push(JsonStr(s))
  }
  JsonObj(
    Map([
      ("chrom", JsonStr(callset.chrom)),
      ("samples", JsonArr(samples)),
      ("records", JsonArr(records)),
    ]),
  )
}

///|
/// 调用集 -> JSON 文本.
pub fn gvcf_callset_to_json(callset : GvcfCallset) -> String {
  json_to_string(gvcf_callset_value(callset))
}

///|
fn gvcf_json_to_int(value : JsonValue) -> Int32? {
  match value {
    JsonNum(x) => Some(x.to_int())
    _ => None
  }
}

///|
fn gvcf_json_genotype(value : JsonValue) -> GvcfGenotype? {
  match value {
    JsonObj(fields) => {
      let mut gt = "./."
      let pl : Array[Int32] = []
      let mut dp = 0
      let mut gq = 0
      for field in fields {
        if field.0 == "gt" {
          match field.1 {
            JsonStr(s) => gt = s
            _ => ()
          }
        } else if field.0 == "dp" {
          match gvcf_json_to_int(field.1) {
            Some(x) => dp = x
            None => ()
          }
        } else if field.0 == "gq" {
          match gvcf_json_to_int(field.1) {
            Some(x) => gq = x
            None => ()
          }
        } else if field.0 == "pl" {
          match field.1 {
            JsonArr(arr) =>
              for entry in arr {
                match gvcf_json_to_int(entry) {
                  Some(x) => pl.push(x)
                  None => ()
                }
              }
            _ => ()
          }
        }
      }
      Some(GvcfGenotype::{ gt, pl, dp, gq, })
    }
    _ => None
  }
}

///|
fn gvcf_json_record(value : JsonValue) -> GvcfRecord? {
  match value {
    JsonObj(fields) => {
      let mut chrom = ""
      let mut pos = 0
      let mut ref_allele = ""
      let alts : Array[String] = []
      let mut is_block = false
      let mut block_end = 0
      let samples : Map[String, GvcfGenotype] = Map([], capacity=0)
      for field in fields {
        if field.0 == "chrom" {
          match field.1 {
            JsonStr(s) => chrom = s
            _ => ()
          }
        } else if field.0 == "pos" {
          match gvcf_json_to_int(field.1) {
            Some(x) => pos = x
            None => ()
          }
        } else if field.0 == "ref" {
          match field.1 {
            JsonStr(s) => ref_allele = s
            _ => ()
          }
        } else if field.0 == "is_block" {
          match field.1 {
            JsonBool(b) => is_block = b
            _ => ()
          }
        } else if field.0 == "block_end" {
          match gvcf_json_to_int(field.1) {
            Some(x) => block_end = x
            None => ()
          }
        } else if field.0 == "alts" {
          match field.1 {
            JsonArr(arr) =>
              for entry in arr {
                match entry {
                  JsonStr(s) => alts.push(s)
                  _ => ()
                }
              }
            _ => ()
          }
        } else if field.0 == "samples" {
          match field.1 {
            JsonObj(smap) =>
              for entry in smap {
                match gvcf_json_genotype(entry.1) {
                  Some(geno) => samples.set(entry.0, geno)
                  None => ()
                }
              }
            _ => ()
          }
        }
      }
      if chrom == "" {
        None
      } else {
        Some(GvcfRecord::{
          chrom,
          pos,
          ref_allele,
          alts,
          is_block,
          block_end,
          samples,
        })
      }
    }
    _ => None
  }
}

///|
/// JSON 文本 -> 调用集(roundtrip).
pub fn gvcf_callset_from_json(text : String) -> GvcfCallset? {
  match parse_json(text) {
    Some(JsonObj(fields)) => {
      let mut chrom = ""
      let samples : Array[String] = []
      let records : Array[GvcfRecord] = []
      for field in fields {
        if field.0 == "chrom" {
          match field.1 {
            JsonStr(s) => chrom = s
            _ => ()
          }
        } else if field.0 == "samples" {
          match field.1 {
            JsonArr(arr) =>
              for entry in arr {
                match entry {
                  JsonStr(s) => samples.push(s)
                  _ => ()
                }
              }
            _ => ()
          }
        } else if field.0 == "records" {
          match field.1 {
            JsonArr(arr) =>
              for entry in arr {
                match gvcf_json_record(entry) {
                  Some(record) => records.push(record)
                  None => ()
                }
              }
            _ => ()
          }
        }
      }
      if chrom == "" {
        return None
      }
      Some(GvcfCallset::{ chrom, samples, records, })
    }
    _ => None
  }
}

///|
/// 调用集摘要行.
pub fn gvcf_summary(callset : GvcfCallset) -> String {
  let mut n_variant = 0
  let mut n_block = 0
  for record in callset.records {
    if record.is_block {
      n_block = n_block + 1
    } else {
      n_variant = n_variant + 1
    }
  }
  "gvcf-lite: chrom=\{callset.chrom}, samples=\{callset.samples.length()}, variants=\{n_variant}, blocks=\{n_block}"
}
