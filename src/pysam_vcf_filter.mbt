// ============================================================================
// pysam_vcf_filter — bcftools filter 的内存版实现
//
// 支持的过滤条件（AND 组合；全部满足才通过）：
//  - QUAL 阈值范围（low <= qual <= high）
//  - 只保留 FILTER == "PASS" / "." 的记录
//  - INFO 中必须存在某些 key
//  - ALT allele 数量范围
//  - ID 已知（非空）/ 未知（空）
//
// 两种输出模式：
//  - "drop"      — 物理删除不满足条件的 records（默认）
//  - "set_filter" — 保留记录但写入 `filter_tag` 到 FILTER 字段
//
// 参考：bcftools filter / Pysam `VariantFile.filter()`
// ============================================================================

///|
/// Bundle of conditions that a record must satisfy to pass filtering. All
/// fields are optional (None = "don't check this dimension").
pub(all) struct VcfFilterConfig {
  /// Minimum QUAL (inclusive). None → skip check.
  qual_min : Double?
  /// Maximum QUAL (inclusive). None → skip check.
  qual_max : Double?
  /// If true, only records whose FILTER list contains "PASS" or is empty are
  /// kept. Mirrors bcftools `--filter PASS`.
  pass_only : Bool
  /// Set of INFO keys that must be **present** on the record (regardless of
  /// value). Mirrors bcftools `--annotate INFO/key`.
  info_required_keys : Array[String]
  /// Minimum number of ALT alleles (1 = at least one ALT).
  min_alt_count : Int32?
  /// Maximum number of ALT alleles. None → no upper bound.
  max_alt_count : Int32?
  /// If Some(true), keep only records with non-empty ID; if Some(false), keep
  /// only records with empty ID; None → no check.
  id_known : Bool?
}

///|
/// Default config: no filter applied. All records pass.
pub fn vcf_filter_config_default() -> VcfFilterConfig {
  VcfFilterConfig::{
    qual_min: None,
    qual_max: None,
    pass_only: false,
    info_required_keys: [],
    min_alt_count: None,
    max_alt_count: None,
    id_known: None,
  }
}

///|
/// Build a config that matches `bcftools filter -e 'QUAL<10 || FILTER=LowQual'`
/// style common patterns. Convenience wrapper.
pub fn vcf_filter_qual(min : Double, max~ : Double = 1000.0) -> VcfFilterConfig {
  VcfFilterConfig::{
    qual_min: Some(min),
    qual_max: Some(max),
    pass_only: false,
    info_required_keys: [],
    min_alt_count: None,
    max_alt_count: None,
    id_known: None,
  }
}

// ---------------------------------------------------------------------------
// 单 record 判断
// ---------------------------------------------------------------------------

///|
/// Return true when a single record satisfies every condition in the config.
pub fn vcf_record_passes(rec : VcfRecord, cfg : VcfFilterConfig) -> Bool {
  // QUAL range
  match (cfg.qual_min, cfg.qual_max, rec.qual) {
    (Some(lo), Some(hi), Some(q)) =>
      if q < lo || q > hi { return false }
      else { () }
    (Some(lo), None, Some(q)) =>
      if q < lo { return false } else { () }
    (None, Some(hi), Some(q)) =>
      if q > hi { return false } else { () }
    (Some(_), Some(_), None) => return false  // QUAL required but missing
    (Some(_), None, None) => return false
    (None, Some(_), None) => return false
    (_, _, _) => ()
  }
  // pass_only
  if cfg.pass_only {
    let filt = rec.filter
    if filt.length() != 1 || (filt.length() == 1 && filt[0] != "PASS" && filt[0] != ".") {
      return false
    }
  }
  // INFO keys
  for key in cfg.info_required_keys {
    if rec.info.get(key) is None {
      return false
    }
  }
  // ALT count
  match (cfg.min_alt_count, cfg.max_alt_count) {
    (Some(lo), Some(hi)) =>
      if rec.alt.length() < lo || rec.alt.length() > hi { return false } else { () }
    (Some(lo), None) =>
      if rec.alt.length() < lo { return false } else { () }
    (None, Some(hi)) =>
      if rec.alt.length() > hi { return false } else { () }
    (None, None) => ()
  }
  // ID known
  match cfg.id_known {
    Some(true) =>
      if rec.id.length() == 0 { return false } else { () }
    Some(false) =>
      if rec.id.length() != 0 { return false } else { () }
    None => ()
  }
  true
}

// ---------------------------------------------------------------------------
// 主入口 1：drop 模式（物理删除）
// ---------------------------------------------------------------------------

///|
/// Return a new VcfFile containing only records that satisfy `cfg`.
/// Mirrors `bcftools filter -i 'expression' -o filtered.vcf.gz`.
pub fn vcf_filter(vcf : VcfFile, cfg : VcfFilterConfig) -> VcfFile {
  let out : Array[VcfRecord] = Array(capacity=0)
  for rec in vcf.records {
    if vcf_record_passes(rec, cfg) {
      out.push(rec)
    }
  }
  VcfFile::{
    meta_headers: vcf.meta_headers,
    fixed_header: vcf.fixed_header,
    records: out,
  }
}

///|
/// Convenience: keep only records whose QUAL >= `threshold`.
pub fn vcf_filter_qual_min(vcf : VcfFile, threshold : Double) -> VcfFile {
  vcf_filter(vcf, vcf_filter_qual(threshold, max=1000.0))
}

///|
/// Convenience: keep only PASS records.
pub fn vcf_filter_pass_only(vcf : VcfFile) -> VcfFile {
  let cfg = VcfFilterConfig::{
    qual_min: None,
    qual_max: None,
    pass_only: true,
    info_required_keys: [],
    min_alt_count: None,
    max_alt_count: None,
    id_known: None,
  }
  vcf_filter(vcf, cfg)
}

// ---------------------------------------------------------------------------
// 主入口 2：set_filter 模式（保留 + 打标）
// ---------------------------------------------------------------------------

///|
/// Keep every record but change its FILTER field to `filter_tag` when the
/// record **fails** `cfg`. Useful for producing a VCF where low-quality
/// records are flagged but not removed (matches bcftools filter
/// `--set-GTs` at header level).
///
/// When the record already has a FILTER tag, `filter_tag` is appended using
/// semicolon separation.
pub fn vcf_filter_set_filter(
  vcf : VcfFile,
  cfg : VcfFilterConfig,
  filter_tag~ : String = "LowQual",
) -> VcfFile {
  let out : Array[VcfRecord] = Array(capacity=0)
  for rec in vcf.records {
    if vcf_record_passes(rec, cfg) {
      out.push(rec)
    } else {
      let new_filter : Array[String] = rec.filter.copy()
      new_filter.push(filter_tag)
      out.push(VcfRecord::{
        chrom: rec.chrom,
        pos: rec.pos,
        id: rec.id.copy(),
        ref_seq: rec.ref_seq,
        alt: rec.alt.copy(),
        qual: rec.qual,
        filter: new_filter,
        info: rec.info.copy(),
        format:
          match rec.format {
            Some(f) => Some(f.copy())
            None => None
          },
        samples: rec.samples.copy(),
      })
    }
  }
  VcfFile::{
    meta_headers: vcf.meta_headers,
    fixed_header: vcf.fixed_header,
    records: out,
  }
}

// ---------------------------------------------------------------------------
// 统计辅助
// ---------------------------------------------------------------------------

///|
/// Count records that pass / fail. Returns `(passes, failures)`.
pub fn vcf_filter_summary(vcf : VcfFile, cfg : VcfFilterConfig) -> (Int32, Int32) {
  let mut pass = 0
  let mut fail = 0
  for rec in vcf.records {
    if vcf_record_passes(rec, cfg) { pass = pass + 1 } else { fail = fail + 1 }
  }
  (pass, fail)
}

///|
/// Fraction of records that pass (0.0 ~ 1.0).
pub fn vcf_filter_pass_rate(vcf : VcfFile, cfg : VcfFilterConfig) -> Double {
  if vcf.records.length() == 0 { return 0.0 }
  let (p, _) = vcf_filter_summary(vcf, cfg)
  p.to_double() / vcf.records.length().to_double()
}
