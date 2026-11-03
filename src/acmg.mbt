///|
/// acmg.mbt — ACMG 2015 (Richards et al.) 临床变异致病性自动化分级
/// 与 ClinVar 记录摄取(本地确定性实现, InterVar 风格).
///
/// 证据规则(可计算子集):
/// - PVS1(致病极强): LoF 后果(无义 stop_gained / 移码 frameshift /
///   剪接供受体 splice_donor|acceptor ±1-2 / 起始丢失 start_lost);
/// - PS1(致病极强): 与 ClinVar 已知致病变异同位点同 alt;
/// - PM2(致病中等): 人群最大等位频率 < 0.001;
/// - PM4(中等): 框内 indel 改变蛋白长度(非重复区语义省略);
/// - PP2(致病支持): 错义变异位于"仅错义致病"基因;
/// - PP3(支持): 计算预测有害(预测分数 >= 阈值);
/// - PP5(支持): ClinVar 收录(>=2 星)无冲突;
/// - PM6(中等): de novo(父母未携带);
/// - BA1(良性极强): 人群频率 > 0.05;
/// - BS1(良性强): 人群频率高于致病预期(> 0.01);
/// - BP1(良性支持): 错义位于"LoF 主导"基因;
/// - BP4(良性支持): 计算预测良性(分数 < 阈值);
/// - BP6(良性支持): ClinVar 收录良性(>=2 星)无冲突;
/// - BP7(良性支持): 同义变异且远离剪接区;
/// - PS3/BS3(功能实验): functional_score 外部输入;
///
/// 组合规则(Richards et al. 2015 表 5, InterVar 实现):
/// - Pathogenic: (PVS1 & >=1PS) | (PVS1 & >=2PM) | (>=1PS & >=3PM) |
///   (>=1PS & 2PM & >=2PP) | (>=1PS & >=4PP) | (>=2PS);
/// - Likely Pathogenic: (PVS1 & 1PM) | (>=1PS & 1..2PM) |
///   (>=1PS & >=2PP) | (>=3PM) | (2PM & >=2PP) | (1PM & >=4PP);
/// - Benign: >=1BA | >=2BS;
/// - Likely Benign: (1BS & >=1BP) | >=2BP;
/// - VUS: 其余; 同时满足致病与良性组合时报冲突(降级 VUS 并说明).
///
/// 数据不出本机.

// ---------------------------------------------------------------------------
// ClinVar 摄取
// ---------------------------------------------------------------------------

///|
/// 一条 ClinVar 记录(简化 VCF 风格 TSV:
/// chrom pos ref alt CLNSIG review_star).
pub(all) struct ClinVarRecord {
  chrom : String
  pos : Int32
  ref_allele : String
  alt : String
  /// Pathogenic / Likely_pathogenic / VUS / Likely_benign / Benign /
  /// Conflicting_classifications.
  clnsig : String
  /// 审查星级 0-4.
  review_star : Int32
}

///|
/// 解析 ClinVar TSV 文本(每行: chrom pos ref alt clnsig star;
/// 以 '#' 开头的行与表头行跳过).
pub fn clinvar_parse(content : String) -> Array[ClinVarRecord] {
  let out : Array[ClinVarRecord] = Array(capacity=0)
  let lines = content.split("\n").to_array()
  for raw in lines {
    let line = raw.to_owned()
    if line.length() == 0 || line.unsafe_get(0).unsafe_to_char() == '#' {
      continue
    }
    let parts = line.split("\t").to_array()
    if parts.length() < 5 {
      continue
    }
    let chrom = parts[0].to_owned()
    // 表头行(非数字 pos)跳过
    let pos = acmg_parse_int_or(parts[1].to_owned(), -1)
    if pos <= 0 {
      continue
    }
    let ref_allele = parts[2].to_owned()
    let alt = parts[3].to_owned()
    let clnsig = parts[4].to_owned()
    let star = if parts.length() >= 6 {
      acmg_parse_int_or(parts[5].to_owned(), 0)
    } else {
      0
    }
    out.push(ClinVarRecord::{
      chrom,
      pos,
      ref_allele,
      alt,
      clnsig,
      review_star: star,
    })
  }
  out
}

///|
fn acmg_is_digit(c : Char) -> Bool {
  c == '0' ||
  c == '1' ||
  c == '2' ||
  c == '3' ||
  c == '4' ||
  c == '5' ||
  c == '6' ||
  c == '7' ||
  c == '8' ||
  c == '9'
}

///|
fn acmg_parse_int_or(s : String, fallback : Int32) -> Int32 {
  if s.length() == 0 {
    return fallback
  }
  let mut ok = true
  for i = 0; i < s.length(); i = i + 1 {
    let c = s.unsafe_get(i).unsafe_to_char()
    if acmg_is_digit(c) == false && !(i == 0 && c == '-') {
      ok = false
      break
    }
  }
  if ok {
    parse_int(s)
  } else {
    fallback
  }
}

///|
/// ClinVar 记录查找(同染色体同位置同 ref/alt; 无则 None).
pub fn clinvar_lookup(
  records : Array[ClinVarRecord],
  chrom : String,
  pos : Int32,
  ref_allele : String,
  alt : String,
) -> ClinVarRecord? {
  for r in records {
    if r.chrom == chrom &&
      r.pos == pos &&
      r.ref_allele == ref_allele &&
      r.alt == alt {
      return Some(r)
    }
  }
  None
}

// ---------------------------------------------------------------------------
// 证据输入与推导
// ---------------------------------------------------------------------------

///|
/// ACMG 证据推导输入(可自动计算子集; 外部证据由调用方补充).
pub(all) struct AcmgInput {
  chrom : String
  pos : Int32
  ref_allele : String
  alt : String
  /// VEP/SnpEff 后果标签(stop_gained / frameshift_variant / ...).
  consequence : String
  /// 人群最大等位频率(0 无数据时视为缺失).
  popmax_af : Double
  /// de novo(父母未携带)标志.
  de_novo : Bool
  /// 基因属性: LoF 为主导致病机制.
  gene_lof_major : Bool
  /// 基因属性: 仅错义致病.
  gene_missense_only : Bool
  /// 计算致病性预测分数(如 CADD 风格; 无数据 None).
  in_silico_score : Double?
  /// 阈值: 分数 >= 该值判有害. 默认 20.
  in_silico_deleterious_threshold : Double
  /// 功能实验分数(阳性支持致病, 阴性支持良性; None 无数据).
  functional_score : Double?
}

///|
/// 构造 ACMG 输入(外部证据以可选参数提供).
pub fn acmg_make_input(
  chrom : String,
  pos : Int32,
  ref_allele : String,
  alt : String,
  consequence : String,
  popmax_af~ : Double = 0.0,
  de_novo~ : Bool = false,
  gene_lof_major~ : Bool = false,
  gene_missense_only~ : Bool = false,
  in_silico_score~ : Double? = None,
  in_silico_deleterious_threshold~ : Double = 20.0,
  functional_score~ : Double? = None,
) -> AcmgInput {
  {
    chrom,
    pos,
    ref_allele,
    alt,
    consequence,
    popmax_af,
    de_novo,
    gene_lof_major,
    gene_missense_only,
    in_silico_score,
    in_silico_deleterious_threshold,
    functional_score,
  }
}

///|
/// 默认输入(全部外部证据缺失).
pub fn acmg_default_input(
  chrom : String,
  pos : Int32,
  ref_allele : String,
  alt : String,
  consequence : String,
) -> AcmgInput {
  acmg_make_input(chrom, pos, ref_allele, alt, consequence)
}

///|
fn acmg_is_lof(consequence : String) -> Bool {
  consequence == "stop_gained" ||
  consequence == "frameshift_variant" ||
  consequence == "splice_donor_variant" ||
  consequence == "splice_acceptor_variant" ||
  consequence == "start_lost"
}

///|
fn acmg_is_missense(consequence : String) -> Bool {
  consequence == "missense_variant"
}

///|
fn acmg_is_synonymous(consequence : String) -> Bool {
  consequence == "synonymous_variant"
}

///|
fn acmg_is_inframe_indel(consequence : String) -> Bool {
  consequence == "inframe_insertion" ||
  consequence == "inframe_deletion" ||
  consequence == "inframe_variant"
}

///|
/// 推导可自动计算的 ACMG 证据代码集合(不含外部证据).
pub fn acmg_derive_evidence(input : AcmgInput) -> Array[String] {
  let ev : Array[String] = Array(capacity=8)
  // 良性证据优先计算(频率硬约束)
  if input.popmax_af > 0.05 {
    ev.push("BA1")
  } else if input.popmax_af > 0.01 {
    ev.push("BS1")
  }
  if input.popmax_af > 0.0 && input.popmax_af < 0.001 {
    ev.push("PM2")
  }
  if input.popmax_af == 0.0 {
    // 频率缺失时不计 PM2(区分 absent 与无数据由调用方以极小值表达)
    ignore(0)
  }
  // 后果类
  if acmg_is_lof(input.consequence) {
    ev.push("PVS1")
  }
  if acmg_is_inframe_indel(input.consequence) {
    ev.push("PM4")
  }
  if acmg_is_missense(input.consequence) {
    if input.gene_missense_only {
      ev.push("PP2")
    }
    if input.gene_lof_major {
      ev.push("BP1")
    }
  }
  if acmg_is_synonymous(input.consequence) {
    // BP7: 同义且远离剪接区(splice_region 被单独标注时不算)
    if input.consequence != "splice_region_variant" {
      ev.push("BP7")
    }
  }
  // de novo
  if input.de_novo {
    ev.push("PM6")
  }
  // 计算预测
  match input.in_silico_score {
    Some(s) =>
      if s >= input.in_silico_deleterious_threshold {
        ev.push("PP3")
      } else {
        ev.push("BP4")
      }
    None => ()
  }
  // 功能实验
  match input.functional_score {
    Some(f) => if f >= 0.0 { ev.push("PS3") } else { ev.push("BS3") }
    None => ()
  }
  ev
}

///|
/// ClinVar 收录补充证据(PS1 同位点致病变异 / PP5 / BP6; 冲突不计).
pub fn acmg_clinvar_evidence(
  evidence : Array[String],
  clinvar : ClinVarRecord?,
) -> Array[String] {
  let out = evidence.copy()
  match clinvar {
    Some(r) => {
      if r.clnsig == "Conflicting_classifications" {
        return out
      }
      if r.review_star >= 2 {
        if r.clnsig == "Pathogenic" {
          out.push("PS1")
          out.push("PP5")
        } else if r.clnsig == "Likely_pathogenic" {
          out.push("PP5")
        } else if r.clnsig == "Benign" {
          out.push("BP6")
        } else if r.clnsig == "Likely_benign" {
          out.push("BP6")
        }
      }
    }
    None => ()
  }
  out
}

// ---------------------------------------------------------------------------
// 分级
// ---------------------------------------------------------------------------

///|
fn acmg_count(evidence : Array[String], prefix : String) -> Int32 {
  let mut n = 0
  for e in evidence {
    if e.length() >= 2 && e[0:2].to_owned() == prefix {
      n = n + 1
    }
  }
  n
}

///|
/// ACMG 五级分类结果.
pub(all) struct AcmgClassification {
  /// Pathogenic / Likely_pathogenic / VUS / Likely_benign / Benign.
  classification : String
  /// 命中的证据代码.
  evidence : Array[String]
  /// 组合路径说明(为何得到该分类).
  reasons : Array[String]
}

///|
/// ACMG 2015 组合规则分级(良性证据优先; 致病与良性同时满足时
/// 报冲突并降级 VUS).
pub fn acmg_classify(evidence : Array[String]) -> AcmgClassification {
  let has_pvs1 = evidence.contains("PVS1")
  let n_ps = acmg_count(evidence, "PS")
  let n_pm = acmg_count(evidence, "PM")
  let n_pp = acmg_count(evidence, "PP")
  let n_ba = acmg_count(evidence, "BA")
  let n_bs = acmg_count(evidence, "BS")
  let n_bp = acmg_count(evidence, "BP")
  // 良性组合
  let benign = n_ba >= 1 || n_bs >= 2
  let likely_benign = (n_bs == 1 && n_bp >= 1) || n_bp >= 2
  // 致病组合
  let pathogenic = (has_pvs1 && (n_ps >= 1 || n_pm >= 2)) ||
    (n_ps >= 1 && n_pm >= 3) ||
    (n_ps >= 1 && n_pm >= 2 && n_pp >= 2) ||
    (n_ps >= 1 && n_pp >= 4) ||
    n_ps >= 2
  let likely_pathogenic = (has_pvs1 && n_pm == 1) ||
    (n_ps >= 1 && n_pm >= 1 && n_pm <= 2) ||
    (n_ps >= 1 && n_pp >= 2) ||
    n_pm >= 3 ||
    (n_pm == 2 && n_pp >= 2) ||
    (n_pm == 1 && n_pp >= 4)
  let reasons : Array[String] = Array(capacity=4)
  if benign {
    reasons.push(
      "benign combination: \{n_ba}BA + \{n_bs}BS (ACMG Table 5 benign rules)",
    )
  }
  if pathogenic {
    reasons.push(
      "pathogenic combination: PVS1=\{has_pvs1}, \{n_ps}PS, \{n_pm}PM, \{n_pp}PP (ACMG Table 5 pathogenic rules)",
    )
  }
  if likely_pathogenic && pathogenic == false {
    reasons.push(
      "likely-pathogenic combination: PVS1=\{has_pvs1}, \{n_ps}PS, \{n_pm}PM, \{n_pp}PP",
    )
  }
  let classification = if benign && (pathogenic || likely_pathogenic) {
    reasons.push(
      "conflicting pathogenic and benign criteria - downgraded to VUS per ACMG guidelines",
    )
    "VUS"
  } else if benign {
    "Benign"
  } else if pathogenic {
    "Pathogenic"
  } else if likely_pathogenic && likely_benign {
    reasons.push(
      "conflicting likely-pathogenic and likely-benign criteria - downgraded to VUS",
    )
    "VUS"
  } else if likely_pathogenic {
    "Likely_pathogenic"
  } else if likely_benign {
    "Likely_benign"
  } else {
    reasons.push("insufficient criteria for any class - VUS by default")
    "VUS"
  }
  { classification, evidence, reasons, }
}

// ---------------------------------------------------------------------------
// 批量注释
// ---------------------------------------------------------------------------

///|
/// 一条变异的完整 ACMG 解读结果.
pub(all) struct AcmgResult {
  chrom : String
  pos : Int32
  ref_allele : String
  alt : String
  consequence : String
  evidence : Array[String]
  classification : String
  reasons : Array[String]
  /// ClinVar 一致性(空串 = 无记录).
  clinvar_sig : String
}

///|
/// 批量解读: 输入列表 + ClinVar 库 -> 逐条证据推导 + 分级.
pub fn acmg_annotate(
  inputs : Array[AcmgInput],
  clinvar~ : Array[ClinVarRecord] = [],
) -> Array[AcmgResult] {
  let out : Array[AcmgResult] = Array(capacity=inputs.length())
  for input in inputs {
    let cv = clinvar_lookup(
      clinvar,
      input.chrom,
      input.pos,
      input.ref_allele,
      input.alt,
    )
    let ev0 = acmg_derive_evidence(input)
    let ev = acmg_clinvar_evidence(ev0, cv)
    let cls = acmg_classify(ev)
    let clinvar_sig = match cv {
      Some(r) => "\{r.clnsig}(\{r.review_star}star)"
      None => ""
    }
    out.push(AcmgResult::{
      chrom: input.chrom,
      pos: input.pos,
      ref_allele: input.ref_allele,
      alt: input.alt,
      consequence: input.consequence,
      evidence: ev,
      classification: cls.classification,
      reasons: cls.reasons,
      clinvar_sig,
    })
  }
  out
}

// ---------------------------------------------------------------------------
// 导出
// ---------------------------------------------------------------------------

///|
/// 解读结果 -> TSV 文本(chrom pos ref alt consequence evidence
/// classification clinvar).
pub fn acmg_report_tsv(results : Array[AcmgResult]) -> String {
  let sb = StringBuilder()
  sb.write_string(
    "chrom\tpos\tref\talt\tconsequence\tevidence\tclassification\tclinvar\n",
  )
  for r in results {
    let mut ev = ""
    for i = 0; i < r.evidence.length(); i = i + 1 {
      if i > 0 {
        ev = ev + "&"
      }
      ev = ev + r.evidence[i]
    }
    sb.write_string(r.chrom)
    sb.write_char('\t')
    sb.write_string(r.pos.to_string())
    sb.write_char('\t')
    sb.write_string(r.ref_allele)
    sb.write_char('\t')
    sb.write_string(r.alt)
    sb.write_char('\t')
    sb.write_string(r.consequence)
    sb.write_char('\t')
    sb.write_string(ev)
    sb.write_char('\t')
    sb.write_string(r.classification)
    sb.write_char('\t')
    sb.write_string(r.clinvar_sig)
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// 解读结果 -> JSON 值.
pub fn acmg_value(results : Array[AcmgResult]) -> JsonValue {
  let arr : Array[JsonValue] = Array(capacity=results.length())
  for r in results {
    let ev : Array[JsonValue] = Array(capacity=r.evidence.length())
    for e in r.evidence {
      ev.push(JsonStr(e))
    }
    let reasons : Array[JsonValue] = Array(capacity=r.reasons.length())
    for x in r.reasons {
      reasons.push(JsonStr(x))
    }
    arr.push(
      JsonObj(
        Map([
          ("chrom", JsonStr(r.chrom)),
          ("pos", JsonNum(r.pos.to_double())),
          ("ref", JsonStr(r.ref_allele)),
          ("alt", JsonStr(r.alt)),
          ("consequence", JsonStr(r.consequence)),
          ("evidence", JsonArr(ev)),
          ("classification", JsonStr(r.classification)),
          ("reasons", JsonArr(reasons)),
          ("clinvar", JsonStr(r.clinvar_sig)),
        ]),
      ),
    )
  }
  JsonArr(arr)
}

///|
/// 解读结果 -> JSON 文本.
pub fn acmg_json(results : Array[AcmgResult]) -> String {
  json_to_string(acmg_value(results))
}

///|
/// 分级计数摘要行.
pub fn acmg_summary_line(results : Array[AcmgResult]) -> String {
  let mut n_p = 0
  let mut n_lp = 0
  let mut n_vus = 0
  let mut n_lb = 0
  let mut n_b = 0
  for r in results {
    match r.classification {
      "Pathogenic" => n_p = n_p + 1
      "Likely_pathogenic" => n_lp = n_lp + 1
      "VUS" => n_vus = n_vus + 1
      "Likely_benign" => n_lb = n_lb + 1
      "Benign" => n_b = n_b + 1
      _ => ()
    }
  }
  "acmg: total=\{results.length()}, pathogenic=\{n_p}, likely_pathogenic=\{n_lp}, vus=\{n_vus}, likely_benign=\{n_lb}, benign=\{n_b}"
}
