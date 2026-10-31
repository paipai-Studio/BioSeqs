///|
/// variant_benchmark.mbt — truth vs query 变异检出基准评估(vcfeval / hap.py 风格).
///
/// 提供:
/// - VCF 记录归一化(拆分多等位 + 去除 REF/ALT 公共前后缀)
/// - (chrom, pos, ref, alt) 键匹配
/// - 总体与按变异类型的 TP/FP/FN、Precision/Recall/F1
/// - 基因型(alt 剂量)一致性统计
/// - 按质量阈值生成 PR 曲线数据点
/// - JSON 输出, 数据不出本机

///|
/// 基准评估错误.
pub suberror VariantBenchmarkError {
  VariantBenchmarkError(String)
}

///|
/// 变异类型.
pub(all) enum VariantKind {
  Snp
  Indel
  Mnv
  Other
} derive(Eq, Debug)

///|
/// 归一化后的变异记录.
pub(all) struct BenchVariant {
  chrom : String
  pos : Int32
  ref_allele : String
  alt_allele : String
  /// 原始 GT 字符串(如 "0/1"), 缺失为 None.
  gt : String?
  /// alt 等位基因剂量(0/1/2), GT 缺失为 None.
  alt_dosage : Int32?
  qual : Double?
  is_pass : Bool
  kind : VariantKind
}

///|
/// 比较结果.
pub(all) struct BenchStats {
  truth_total : Int32
  query_total : Int32
  tp : Int32
  fp : Int32
  fn_ : Int32
  genotype_mismatches : Int32
  precision : Double
  recall : Double
  f1 : Double
  snp_truth : Int32
  snp_query : Int32
  snp_tp : Int32
  indel_truth : Int32
  indel_query : Int32
  indel_tp : Int32
}

///|
/// PR 曲线数据点.
pub(all) struct BenchPrPoint {
  threshold : Double
  tp : Int32
  fp : Int32
  fn_ : Int32
  precision : Double
  recall : Double
  f1 : Double
}

///|
/// 从 FORMAT / 样本字段中取 GT.
fn vcf_bench_gt(rec : VcfRecord, sample : String) -> String? {
  match rec.format {
    Some(fmt) => {
      let mut gt_idx = -1
      for i = 0; i < fmt.length(); i = i + 1 {
        if fmt[i] == "GT" {
          gt_idx = i
          break
        }
      }
      if gt_idx < 0 {
        return None
      }
      match rec.samples.get(sample) {
        Some(data) =>
          if gt_idx < data.length() {
            Some(data[gt_idx])
          } else {
            None
          }
        None => None
      }
    }
    None => None
  }
}

///|
/// GT -> alt 剂量; 缺失返回 None.
fn geno_alt_dosage(gt : String?) -> Int32? raise {
  match gt {
    None => None
    Some(text) => {
      let t = text.trim()
      if t == "." || t == "./." || t == ".|." || t.length() == 0 {
        return None
      }
      let slash_parts = t.split("/").to_array()
      let parts = if slash_parts.length() == 1 {
        t.split("|").to_array()
      } else {
        slash_parts
      }
      let mut dosage = 0
      for p in parts {
        let allele = p.to_owned()
        if allele == "." || allele.length() == 0 {
          return None
        }
        if allele == "1" {
          dosage = dosage + 1
        }
      }
      Some(dosage)
    }
  }
}

///|
/// 判断变异类型.
fn variant_kind(ref_allele : String, alt_allele : String) -> VariantKind {
  if ref_allele.length() == 1 && alt_allele.length() == 1 {
    Snp
  } else if ref_allele.length() == alt_allele.length() {
    if ref_allele.length() > 1 {
      Mnv
    } else {
      Other
    }
  } else {
    Indel
  }
}

///|
/// 归一化单条 VCF 记录(取第一个 ALT 与指定样本).
pub fn normalize_variant(
  rec : VcfRecord,
  sample~ : String = "",
) -> BenchVariant? raise {
  if rec.alt.length() == 0 {
    return None
  }
  let mut ref_allele = rec.ref_seq.to_upper()
  let mut alt_allele = rec.alt[0].to_upper()
  if alt_allele == "." || alt_allele.length() == 0 {
    return None
  }
  let mut pos = rec.pos
  // 去除公共后缀(至少保留 1 个碱基).
  while ref_allele.length() > 1 &&
        alt_allele.length() > 1 &&
        ref_allele.unsafe_get(ref_allele.length() - 1) ==
        alt_allele.unsafe_get(alt_allele.length() - 1) {
    ref_allele = ref_allele[0:ref_allele.length() - 1].to_owned()
    alt_allele = alt_allele[0:alt_allele.length() - 1].to_owned()
  }
  // 去除公共前缀并前移坐标.
  while ref_allele.length() > 1 &&
        alt_allele.length() > 1 &&
        ref_allele.unsafe_get(0) == alt_allele.unsafe_get(0) {
    ref_allele = ref_allele[1:].to_owned()
    alt_allele = alt_allele[1:].to_owned()
    pos = pos + 1
  }
  let target_sample = if sample.length() > 0 {
    sample
  } else {
    match rec.samples.keys().to_array() {
      [] => ""
      keys => keys[0]
    }
  }
  let gt = if target_sample.length() > 0 {
    vcf_bench_gt(rec, target_sample)
  } else {
    None
  }
  Some(BenchVariant::{
    chrom: rec.chrom,
    pos,
    ref_allele,
    alt_allele,
    gt,
    alt_dosage: geno_alt_dosage(gt),
    qual: rec.qual,
    is_pass: rec.filter.length() == 0 ||
    rec.filter.contains("PASS") ||
    rec.filter.contains("."),
    kind: variant_kind(ref_allele, alt_allele),
  })
}

///|
/// 归一化整个 VCF.
pub fn load_variants(
  vcf : VcfFile,
  sample~ : String = "",
  pass_only~ : Bool = false,
) -> Array[BenchVariant] raise {
  let out : Array[BenchVariant] = Array(capacity=vcf.records.length())
  for rec in vcf.records {
    match normalize_variant(rec, sample~) {
      Some(variant) => if !pass_only || variant.is_pass { out.push(variant) }
      None => ()
    }
  }
  out
}

///|
/// 变异匹配键.
fn bench_variant_key(v : BenchVariant) -> String {
  "\{v.chrom}:\{v.pos.to_string()}:\{v.ref_allele}:\{v.alt_allele}"
}

///|
/// 比较 truth 与 query.
pub fn compare_variants(
  truth : Array[BenchVariant],
  query : Array[BenchVariant],
  pass_only~ : Bool = false,
) -> BenchStats {
  let truth_map : Map[String, BenchVariant] = Map([], capacity=truth.length())
  for v in truth {
    if !pass_only || v.is_pass {
      let key = bench_variant_key(v)
      if !truth_map.contains(key) {
        truth_map[key] = v
      }
    }
  }
  let query_map : Map[String, BenchVariant] = Map([], capacity=query.length())
  for v in query {
    if !pass_only || v.is_pass {
      let key = bench_variant_key(v)
      if !query_map.contains(key) {
        query_map[key] = v
      }
    }
  }
  let mut tp = 0
  let mut mismatches = 0
  let mut snp_truth = 0
  let mut snp_tp = 0
  let mut indel_truth = 0
  let mut indel_tp = 0
  for key, t in truth_map {
    match t.kind {
      Snp => snp_truth = snp_truth + 1
      Indel => indel_truth = indel_truth + 1
      _ => ()
    }
    match query_map.get(key) {
      Some(q) => {
        tp = tp + 1
        match t.kind {
          Snp => snp_tp = snp_tp + 1
          Indel => indel_tp = indel_tp + 1
          _ => ()
        }
        match (t.alt_dosage, q.alt_dosage) {
          (Some(a), Some(b)) => if a != b { mismatches = mismatches + 1 }
          _ => ()
        }
      }
      None => ()
    }
  }
  let mut snp_query = 0
  let mut indel_query = 0
  for _, v in query_map {
    match v.kind {
      Snp => snp_query = snp_query + 1
      Indel => indel_query = indel_query + 1
      _ => ()
    }
  }
  let truth_total = truth_map.length()
  let query_total = query_map.length()
  let fp = query_total - tp
  let fn_ = truth_total - tp
  let precision = if tp + fp > 0 {
    tp.to_double() / (tp + fp).to_double()
  } else {
    0.0
  }
  let recall = if tp + fn_ > 0 {
    tp.to_double() / (tp + fn_).to_double()
  } else {
    0.0
  }
  let f1 = if precision + recall > 0.0 {
    2.0 * precision * recall / (precision + recall)
  } else {
    0.0
  }
  {
    truth_total,
    query_total,
    tp,
    fp,
    fn_,
    genotype_mismatches: mismatches,
    precision,
    recall,
    f1,
    snp_truth,
    snp_query,
    snp_tp,
    indel_truth,
    indel_query,
    indel_tp,
  }
}

///|
/// 按 query 质量阈值生成 PR 曲线数据点.
pub fn benchmark_pr_curve(
  truth : Array[BenchVariant],
  query : Array[BenchVariant],
  thresholds~ : Array[Double] = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
) -> Array[BenchPrPoint] {
  let out : Array[BenchPrPoint] = Array(capacity=thresholds.length())
  for threshold in thresholds {
    let filtered : Array[BenchVariant] = Array(capacity=query.length())
    for v in query {
      let q = match v.qual {
        Some(value) => value
        None => 0.0
      }
      if q >= threshold {
        filtered.push(v)
      }
    }
    let stats = compare_variants(truth, filtered)
    out.push({
      threshold,
      tp: stats.tp,
      fp: stats.fp,
      fn_: stats.fn_,
      precision: stats.precision,
      recall: stats.recall,
      f1: stats.f1,
    })
  }
  out
}

///|
/// 统计结果序列化为 JSON.
pub fn benchmark_to_json(stats : BenchStats) -> String {
  let root : Map[String, JsonValue] = Map([], capacity=10)
  root["truth_total"] = JsonNum(stats.truth_total.to_double())
  root["query_total"] = JsonNum(stats.query_total.to_double())
  root["tp"] = JsonNum(stats.tp.to_double())
  root["fp"] = JsonNum(stats.fp.to_double())
  root["fn"] = JsonNum(stats.fn_.to_double())
  root["genotype_mismatches"] = JsonNum(stats.genotype_mismatches.to_double())
  root["precision"] = JsonNum(stats.precision)
  root["recall"] = JsonNum(stats.recall)
  root["f1"] = JsonNum(stats.f1)
  json_to_string(JsonObj(root))
}
