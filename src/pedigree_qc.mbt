///|
/// pedigree_qc.mbt — 家系(家谱)与基因型一致性 QC(本地确定性实现).
///
/// 面向 VCF 的家系质控:
/// - 家谱解析(.fam 六列)与核心家系(trio)提取;
/// - 孟德尔错误检测: 逐变异位点检查 child 基因型是否能由双亲基因型组合得到;
/// - 样本间基因型一致性(concordance)统计;
/// - 纯合区段(ROH)检测: 连续纯合标记区段(支持最大间隔与最小标记数).
///
/// 与 snprelate.mbt 的亲缘系数/KING 互补: 后者给群体水平相关度, 本模块给
/// 家系结构与逐位点遗传一致性. 数据不出本机.

///|
/// 家系 QC 错误.
pub suberror PedigreeError {
  PedigreeError(String)
}

///|
/// 家系个体(.fam 行).
pub(all) struct PedigreeIndividual {
  family_id : String
  id : String
  father : String?
  mother : String?
  sex : Int32
  phenotype : String
}

///|
/// 家谱(个体集合).
pub(all) struct Pedigree {
  individuals : Array[PedigreeIndividual]
}

///|
/// 核心家系(子 + 双亲).
pub(all) struct Trio {
  child : String
  father : String
  mother : String
}

///|
/// 一个孟德尔错误位点.
pub(all) struct MendelianError {
  chrom : String
  pos : Int32
  child : String
  father_gt : String
  mother_gt : String
  child_gt : String
}

///|
/// 单个 trio 的孟德尔检查结果.
pub(all) struct TrioQc {
  trio : Trio
  sites_checked : Int32
  errors : Array[MendelianError]
}

///|
/// 两个样本的基因型一致性.
pub(all) struct ConcordanceResult {
  sample_a : String
  sample_b : String
  sites_compared : Int32
  concordant : Int32
  discordant : Int32
  concordance : Double
}

///|
/// 一个纯合区段.
pub(all) struct RohSegment {
  chrom : String
  start : Int32
  end : Int32
  n_markers : Int32
}

///|
/// 单个样本的 ROH 结果.
pub(all) struct RohResult {
  sample : String
  segments : Array[RohSegment]
  total_markers : Int32
}

// ============================================================================
// 基础解析
// ============================================================================

///|
fn pq_parse_int(s : String) -> Int32 {
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
/// 按空白(空格/制表符)切分并丢弃空字段.
fn pq_whitespace_tokens(line : String) -> Array[String] {
  let out : Array[String] = []
  for part in line.replace(old="\t", new=" ").split(" ").to_array() {
    let token = part.to_owned()
    if token.length() > 0 {
      out.push(token)
    }
  }
  out
}

///|
/// 父/母编号: "0" 或空视为缺失.
fn pq_parent_id(token : String) -> String? {
  if token == "0" || token == "" || token == "-9" {
    None
  } else {
    Some(token)
  }
}

///|
/// 解析 .fam 家谱文本(六列: FID IID PID MID Sex Phenotype).
pub fn pedigree_parse(text : String) -> Pedigree {
  let individuals : Array[PedigreeIndividual] = []
  for raw in text.split("\n").to_array() {
    let line = raw.trim().to_owned()
    if line.length() == 0 || line.has_prefix("#") {
      continue
    }
    let fields = pq_whitespace_tokens(line)
    if fields.length() < 6 {
      continue
    }
    individuals.push(PedigreeIndividual::{
      family_id: fields[0],
      id: fields[1],
      father: pq_parent_id(fields[2]),
      mother: pq_parent_id(fields[3]),
      sex: pq_parse_int(fields[4]),
      phenotype: fields[5],
    })
  }
  { individuals, }
}

///|
fn pq_has_individual(pedigree : Pedigree, id : String) -> Bool {
  for individual in pedigree.individuals {
    if individual.id == id {
      return true
    }
  }
  false
}

///|
/// 提取核心家系: 双亲均在谱中且非缺失的个体.
pub fn pedigree_trios(pedigree : Pedigree) -> Array[Trio] {
  let trios : Array[Trio] = []
  for individual in pedigree.individuals {
    match (individual.father, individual.mother) {
      (Some(father), Some(mother)) =>
        if pq_has_individual(pedigree, father) &&
          pq_has_individual(pedigree, mother) {
          trios.push(Trio::{ child: individual.id, father, mother, })
        }
      _ => ()
    }
  }
  trios
}

///|
/// 校验家谱: 个体编号唯一, 父/母引用均存在.
pub fn pedigree_validate(pedigree : Pedigree) -> Unit raise PedigreeError {
  let seen : Map[String, Bool] = Map([], capacity=pedigree.individuals.length())
  for individual in pedigree.individuals {
    if seen.contains(individual.id) {
      raise PedigreeError("duplicate individual: \{individual.id}")
    }
    seen.set(individual.id, true)
  }
  for individual in pedigree.individuals {
    match individual.father {
      Some(father) =>
        if !seen.contains(father) {
          raise PedigreeError("unknown father: \{father}")
        }
      None => ()
    }
    match individual.mother {
      Some(mother) =>
        if !seen.contains(mother) {
          raise PedigreeError("unknown mother: \{mother}")
        }
      None => ()
    }
  }
}

// ============================================================================
// 基因型解析
// ============================================================================

///|
/// 从 VCF 记录读取样本字段串(samples 值为整列字符串).
fn pq_sample_field(record : VcfRecord, sample : String) -> String? {
  match record.samples.get(sample) {
    Some(fields) => if fields.length() == 0 { None } else { Some(fields[0]) }
    None => None
  }
}

///|
/// 从 VCF 记录读取样本 GT 字符串(按 FORMAT 中 GT 位置取字段).
fn pq_gt(record : VcfRecord, sample : String) -> String? {
  match pq_sample_field(record, sample) {
    Some(field) => {
      if field == "." || field == "" {
        return None
      }
      let mut gt_idx = 0
      match record.format {
        Some(fmt) =>
          for i, key in fmt {
            if key == "GT" {
              gt_idx = i
            }
          }
        None => ()
      }
      let parts = field.split(":").to_array()
      if gt_idx < parts.length() {
        let gt = parts[gt_idx].to_owned()
        if gt == "" || gt == "." {
          None
        } else {
          Some(gt)
        }
      } else {
        None
      }
    }
    None => None
  }
}

///|
/// 解析 GT 串为等位索引数组; 缺失返回 None.
fn pq_parse_gt(gt : String) -> Array[Int32]? {
  if gt == "." || gt.length() == 0 {
    return None
  }
  let alleles : Array[Int32] = []
  for slash_part in gt.split("/").to_array() {
    for pipe_part in slash_part.split("|").to_array() {
      let token = pipe_part.to_owned()
      if token == "." || token == "" {
        return None
      }
      alleles.push(pq_parse_int(token))
    }
  }
  if alleles.length() == 0 {
    None
  } else {
    Some(alleles)
  }
}

///|
/// 双亲基因型集合是否含指定等位.
fn pq_parent_has(alleles : Array[Int32], allele : Int32) -> Bool {
  for a in alleles {
    if a == allele {
      return true
    }
  }
  false
}

///|
/// 孟德尔一致性: child 每个等位分别来自父/母.
fn pq_mendel_consistent(
  father : Array[Int32],
  mother : Array[Int32],
  child : Array[Int32],
) -> Bool {
  if father.length() < 2 || mother.length() < 2 || child.length() < 2 {
    return true
  }
  let c0 = child[0]
  let c1 = child[1]
  (pq_parent_has(father, c0) && pq_parent_has(mother, c1)) ||
  (pq_parent_has(father, c1) && pq_parent_has(mother, c0))
}

///|
/// 等位多重集是否相等(排序后逐位比较).
fn pq_alleles_equal(a : Array[Int32], b : Array[Int32]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  let sa = a.copy()
  let sb = b.copy()
  sa.sort()
  sb.sort()
  for i = 0; i < sa.length(); i = i + 1 {
    if sa[i] != sb[i] {
      return false
    }
  }
  true
}

// ============================================================================
// 孟德尔错误检测
// ============================================================================

///|
/// 对单个 trio 逐位点检查孟德尔一致性.
pub fn pedigree_mendelian_check(vcf : VcfFile, trio : Trio) -> TrioQc {
  let errors : Array[MendelianError] = []
  let mut checked = 0
  for record in vcf.records {
    match
      (
        pq_gt(record, trio.father),
        pq_gt(record, trio.mother),
        pq_gt(record, trio.child),
      ) {
      (Some(fg), Some(mg), Some(cg)) =>
        match (pq_parse_gt(fg), pq_parse_gt(mg), pq_parse_gt(cg)) {
          (Some(fa), Some(ma), Some(ca)) =>
            if fa.length() >= 2 && ma.length() >= 2 && ca.length() >= 2 {
              checked = checked + 1
              if !pq_mendel_consistent(fa, ma, ca) {
                errors.push(MendelianError::{
                  chrom: record.chrom,
                  pos: record.pos,
                  child: trio.child,
                  father_gt: fg,
                  mother_gt: mg,
                  child_gt: cg,
                })
              }
            }
          _ => ()
        }
      _ => ()
    }
  }
  { trio, sites_checked: checked, errors, }
}

///|
/// 全家系孟德尔检查.
pub fn pedigree_check_all(vcf : VcfFile, pedigree : Pedigree) -> Array[TrioQc] {
  let out : Array[TrioQc] = []
  for trio in pedigree_trios(pedigree) {
    out.push(pedigree_mendelian_check(vcf, trio))
  }
  out
}

// ============================================================================
// 样本一致性
// ============================================================================

///|
/// 两样本基因型一致性.
pub fn pedigree_concordance(
  vcf : VcfFile,
  sample_a : String,
  sample_b : String,
) -> ConcordanceResult {
  let mut compared = 0
  let mut concordant = 0
  for record in vcf.records {
    match (pq_gt(record, sample_a), pq_gt(record, sample_b)) {
      (Some(a), Some(b)) =>
        match (pq_parse_gt(a), pq_parse_gt(b)) {
          (Some(aa), Some(bb)) =>
            if aa.length() >= 2 && bb.length() >= 2 {
              compared = compared + 1
              if pq_alleles_equal(aa, bb) {
                concordant = concordant + 1
              }
            }
          _ => ()
        }
      _ => ()
    }
  }
  let discordant = compared - concordant
  {
    sample_a,
    sample_b,
    sites_compared: compared,
    concordant,
    discordant,
    concordance: if compared > 0 {
      concordant.to_double() / compared.to_double()
    } else {
      0.0
    },
  }
}

// ============================================================================
// 纯合区段(ROH)
// ============================================================================

///|
/// 单样本 ROH 检测: 连续纯合标记组成区段, 允许区段内标记间隔 <= max_gap,
/// 区段至少含 min_markers 个标记.
pub fn pedigree_roh(
  vcf : VcfFile,
  sample~ : String,
  min_markers~ : Int32 = 5,
  max_gap~ : Int32 = 1000000,
) -> RohResult {
  let markers : Array[(String, Int32, Bool)] = []
  for record in vcf.records {
    match pq_gt(record, sample) {
      Some(gt) =>
        match pq_parse_gt(gt) {
          Some(alleles) =>
            if alleles.length() >= 2 {
              markers.push((record.chrom, record.pos, alleles[0] == alleles[1]))
            }
          None => ()
        }
      None => ()
    }
  }
  let segments : Array[RohSegment] = []
  let mut run_chrom = ""
  let mut run_start = 0
  let mut run_end = 0
  let mut run_count = 0
  for marker in markers {
    let (chrom, pos, hom) = marker
    if hom {
      if run_count > 0 && chrom == run_chrom && pos - run_end <= max_gap {
        run_end = pos
        run_count = run_count + 1
      } else {
        if run_count >= min_markers {
          segments.push(RohSegment::{
            chrom: run_chrom,
            start: run_start,
            end: run_end,
            n_markers: run_count,
          })
        }
        run_chrom = chrom
        run_start = pos
        run_end = pos
        run_count = 1
      }
    } else if run_count > 0 {
      if run_count >= min_markers {
        segments.push(RohSegment::{
          chrom: run_chrom,
          start: run_start,
          end: run_end,
          n_markers: run_count,
        })
      }
      run_count = 0
    }
  }
  if run_count >= min_markers {
    segments.push(RohSegment::{
      chrom: run_chrom,
      start: run_start,
      end: run_end,
      n_markers: run_count,
    })
  }
  { sample, segments, total_markers: markers.length(), }
}

// ============================================================================
// JSON 汇总
// ============================================================================

///|
fn pq_obj(pairs : Array[(String, JsonValue)]) -> JsonValue {
  let m : Map[String, JsonValue] = Map([], capacity=pairs.length())
  for pair in pairs {
    m.set(pair.0, pair.1)
  }
  JsonObj(m)
}

///|
/// 家系 QC 汇总 JSON 值: trio 孟德尔错误 + 各样本 ROH.
pub fn pedigree_qc_report_value(
  pedigree : Pedigree,
  vcf : VcfFile,
  min_markers~ : Int32 = 5,
  max_gap~ : Int32 = 1000000,
) -> JsonValue {
  let trio_checks = pedigree_check_all(vcf, pedigree)
  let trio_json : Array[JsonValue] = Array(capacity=trio_checks.length())
  for check in trio_checks {
    trio_json.push(
      pq_obj([
        ("child", JsonStr(check.trio.child)),
        ("father", JsonStr(check.trio.father)),
        ("mother", JsonStr(check.trio.mother)),
        ("sites_checked", JsonNum(check.sites_checked.to_double())),
        ("mendelian_errors", JsonNum(check.errors.length().to_double())),
      ]),
    )
  }
  let roh_json : Array[JsonValue] = Array(
    capacity=pedigree.individuals.length(),
  )
  for individual in pedigree.individuals {
    let roh = pedigree_roh(vcf, sample=individual.id, min_markers~, max_gap~)
    roh_json.push(
      pq_obj([
        ("sample", JsonStr(individual.id)),
        ("roh_segments", JsonNum(roh.segments.length().to_double())),
        ("roh_markers", JsonNum(roh.total_markers.to_double())),
      ]),
    )
  }
  pq_obj([
    ("individuals", JsonNum(pedigree.individuals.length().to_double())),
    ("trios", JsonArr(trio_json)),
    ("roh", JsonArr(roh_json)),
  ])
}

///|
/// 家系 QC 汇总 JSON 文本.
pub fn pedigree_qc_report_json(
  pedigree : Pedigree,
  vcf : VcfFile,
  min_markers~ : Int32 = 5,
  max_gap~ : Int32 = 1000000,
) -> String {
  json_to_string(
    pedigree_qc_report_value(pedigree, vcf, min_markers~, max_gap~),
  )
}

///|
/// 家系 QC 摘要行.
pub fn pedigree_qc_summary(vcf : VcfFile, pedigree : Pedigree) -> String {
  let trios = pedigree_check_all(vcf, pedigree)
  let mut total_errors = 0
  let mut total_sites = 0
  for check in trios {
    total_errors = total_errors + check.errors.length()
    total_sites = total_sites + check.sites_checked
  }
  "pedigree-qc: individuals=\{pedigree.individuals.length()}, trios=\{trios.length()}, sites=\{total_sites}, mendelian_errors=\{total_errors}"
}
