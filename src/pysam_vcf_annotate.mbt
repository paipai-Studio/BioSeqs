// ============================================================================
// pysam_vcf_annotate — bcftools annotate 的内存版实现
//
// 核心思路：
//   对每条 VCF record，从一个注释源（BED 行 / 另一个 VCF / 内存 Map）里
//   按位置精确匹配或区间重叠匹配查出注释值，然后写入 INFO / ID 字段。
//
// 注意：为避免与 bedtools.mbt 的 BedInterval / bcftools_ops.mbt 的
// vcf_annotate_* 系列同名冲突，本文件类型名统一加 "Annot" 前缀、
// 函数名保留 vcf_annotate_* 但避开已占名。
//
// 参考：bcftools annotate / Pysam `VariantFile.annotate()`
// ============================================================================

///|
/// How a variant record is matched against annotation source rows.
pub enum AnnotMatchMode {
  /// `(chrom, pos)` 必须完全相等。适合 exact-call annotations.
  Exact
  /// BED 区间 `[start, end)` 包含 variant pos 即算命中。
  Overlap
}

pub fn annot_match_mode_exact() -> AnnotMatchMode { AnnotMatchMode::Exact }
pub fn annot_match_mode_overlap() -> AnnotMatchMode { AnnotMatchMode::Overlap }

// ---------------------------------------------------------------------------
// 注释源数据结构（全部加 Annot 前缀避免冲突）
// ---------------------------------------------------------------------------

///|
/// One annotation row from a BED-like source. Coordinates are 0-based
/// half-open (BED convention). `value` is what gets written to the VCF.
///
/// 注意：bedtools.mbt 也定义了 BedInterval，我们这里用 AnnotBedInterval
/// 避免重复声明错误。
pub struct AnnotBedInterval {
  pub chrom : String
  pub start : Int32
  pub end : Int32
  pub value : String
}

///|
/// Annotation source. One of three flavours; all answer the same question:
/// "给我 (chrom, pos)，返回一个可选的注释字符串".
pub enum VcfAnnotSource {
  /// 内存 BED 行数组
  Bed(Array[AnnotBedInterval])
  /// 另一个 VcfFile — 从中复制指定 INFO keys
  /// Second element = keys to copy (from source INFO → dest INFO).
  Vcf(VcfFile, Array[String])
  /// 精确坐标 Map，key = "chrom\tpos"（用 PosMap 命名避免遮蔽 MoonBit `Map` 类型）
  PosMap(Map[String, String])
}

pub fn vcf_annot_source_from_bed(intervals : Array[AnnotBedInterval]) -> VcfAnnotSource {
  VcfAnnotSource::Bed(intervals)
}
pub fn vcf_annot_source_from_vcf(vcf : VcfFile, keys : Array[String]) -> VcfAnnotSource {
  VcfAnnotSource::Vcf(vcf, keys)
}
pub fn vcf_annot_source_from_posmap(m : Map[String, String]) -> VcfAnnotSource {
  VcfAnnotSource::PosMap(m)
}

// ---------------------------------------------------------------------------
// BED 文本 → AnnotBedInterval 数组
// ---------------------------------------------------------------------------

///|
/// Parse BED text (tab-separated lines) into `Array[AnnotBedInterval]`.
///
/// BED 列约定：
///   0  chrom
///   1  start (0-based)
///   2  end   (half-open)
///   3..  extra columns
///
/// `value_col` — **1-based** index of the column whose value becomes the
/// annotation value (default = 4, i.e. BED name column). Set to 0 to use
/// the whole interval as value ("chrom:start-end").
///
/// Lines starting with `#` or "track"/"browser" are skipped.
pub fn parse_bed_annot(text : String, value_col~ : Int32 = 4) -> Array[AnnotBedInterval] {
  let out : Array[AnnotBedInterval] = Array(capacity=0)
  let lines = text.split("\n").to_array()
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].trim()
    if line.length() == 0 || line.unsafe_get(0).unsafe_to_char() == '#' {
      i = i + 1
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 3 {
      i = i + 1
      continue
    }
    let chrom = fields[0].to_owned()
    let start = parse_int(fields[1].to_owned())
    let end = parse_int(fields[2].to_owned())
    let value = if value_col <= 0 {
      chrom + ":" + start.to_string() + "-" + end.to_string()
    } else if value_col - 1 < fields.length() {
      fields[value_col - 1].to_owned()
    } else {
      ""
    }
    out.push(AnnotBedInterval::{ chrom, start, end, value })
    i = i + 1
  }
  out
}

///|
/// Build a `Map[String, String]` keyed by `"chrom\tpos"` from an array of
/// `AnnotBedInterval`. Useful for fast exact `(chrom, pos)` lookups.
pub fn bed_intervals_to_map(intervals : Array[AnnotBedInterval]) -> Map[String, String] {
  let m : Map[String, String] = Map([], capacity=intervals.length())
  for iv in intervals {
    let key = iv.chrom + "\t" + iv.start.to_string()
    m.set(key, iv.value)
  }
  m
}

// ---------------------------------------------------------------------------
// 位置匹配（内部）
// ---------------------------------------------------------------------------

///|
/// Find the first AnnotBedInterval matching `(chrom, pos)` under `mode`.
/// Returns None when nothing matches.
fn find_bed_hit(
  intervals : Array[AnnotBedInterval],
  chrom : String,
  pos : Int32,
  mode : AnnotMatchMode,
) -> AnnotBedInterval? {
  let mut i = 0
  while i < intervals.length() {
    let iv = intervals[i]
    if iv.chrom != chrom {
      i = i + 1
      continue
    }
    match mode {
      Exact => {
        if iv.start == pos { return Some(iv) }
      }
      Overlap => {
        if pos >= iv.start && pos < iv.end { return Some(iv) }
      }
    }
    i = i + 1
  }
  None
}

// ---------------------------------------------------------------------------
// 内部：从任意 VcfAnnotSource 解析单个 record 的注释值
// ---------------------------------------------------------------------------

///|
/// Single point of truth — given a source + record, return the annotation
/// string (or None on miss).
fn resolve_annot_value(
  source : VcfAnnotSource,
  rec : VcfRecord,
  mode : AnnotMatchMode,
) -> String? {
  match source {
    Bed(intervals) => {
      match find_bed_hit(intervals, rec.chrom, rec.pos, mode) {
        Some(iv) => Some(iv.value)
        None => None
      }
    }
    PosMap(m) => {
      // PosMap 只支持精确键查找，忽略 match_mode
      m.get(rec.chrom + "\t" + rec.pos.to_string())
    }
    Vcf(anno, keys) => {
      // Find first anno record matching (chrom, pos)
      let mut i = 0
      while i < anno.records.length() {
        let a = anno.records[i]
        if a.chrom == rec.chrom && a.pos == rec.pos {
          // Return first available key's value
          for k in keys {
            match a.info.get(k) {
              Some(v) => { return Some(v) }
              None => ()
            }
          }
          return None
        }
        i = i + 1
      }
      None
    }
  }
}

// ---------------------------------------------------------------------------
// 主入口 1：从 BED 文本注释 → INFO
// ---------------------------------------------------------------------------

///|
/// Annotate every record in `vcf` using raw BED text lines: find matching
/// BED rows for each variant (by mode) and write the value into
/// `INFO/dest_key`. Non-matching records keep their original INFO intact.
///
/// Mirrors `bcftools annotate -a annotations.bed -c CHROM,FROM,TO,INFO/dest`.
///
/// - `value_col` — 1-based column whose value gets copied; 0 means use the
///   whole interval string. Default 4 (BED name column).
/// - `mode` — `Exact` (chrom+start match) or `Overlap` (interval contains pos).
pub fn vcf_annotate_bed(
  vcf : VcfFile,
  bed_text : String,
  dest_key : String,
  value_col~ : Int32 = 4,
  mode~ : AnnotMatchMode = Overlap,
) -> VcfFile {
  let intervals = parse_bed_annot(bed_text, value_col=value_col)
  let out : Array[VcfRecord] = Array(capacity=vcf.records.length())
  for rec in vcf.records {
    match mode {
      Exact => {
        // 精确匹配：用 Map 加速
        let m = bed_intervals_to_map(intervals)
        let key = rec.chrom + "\t" + rec.pos.to_string()
        match m.get(key) {
          Some(val) => out.push(rec.info_set(dest_key, val))
          None => out.push(rec)
        }
      }
      Overlap => {
        match find_bed_hit(intervals, rec.chrom, rec.pos, Overlap) {
          Some(iv) => out.push(rec.info_set(dest_key, iv.value))
          None => out.push(rec)
        }
      }
    }
  }
  VcfFile::{
    meta_headers: vcf.meta_headers,
    fixed_header: vcf.fixed_header,
    records: out,
  }
}

// ---------------------------------------------------------------------------
// 主入口 2：从另一个 VCF 注释（INFO tag 复制）
// ---------------------------------------------------------------------------

///|
/// For each record in `vcf`, find records in `anno_vcf` at the same
/// `(chrom, pos)` and copy the listed INFO keys over.
/// When no anno record matches, the dest record is emitted unchanged.
///
/// Mirrors `bcftools annotate -a anno.vcf.gz -c INFO/CSQ,INFO/IMPACT`.
///
/// `keys` — INFO keys to copy from anno records. Missing keys on the anno
/// record are silently skipped.
pub fn vcf_annotate_from_vcf(
  vcf : VcfFile,
  anno_vcf : VcfFile,
  keys : Array[String],
  allow_overlap~ : Bool = false,
) -> VcfFile {
  let out : Array[VcfRecord] = Array(capacity=vcf.records.length())
  for rec in vcf.records {
    let matched : Array[VcfRecord] = Array(capacity=0)
    let mut i = 0
    while i < anno_vcf.records.length() {
      let a = anno_vcf.records[i]
      if a.chrom == rec.chrom {
        let a_end : Int32 = a.pos + a.ref_seq.length()
        if allow_overlap {
          if rec.pos >= a.pos && rec.pos < a_end {
            matched.push(a)
          }
        } else {
          if a.pos == rec.pos {
            matched.push(a)
          }
        }
      }
      i = i + 1
    }
    if matched.length() == 0 {
      out.push(rec)
    } else {
      let first = matched[0]
      let mut new_rec = rec
      for k in keys {
        match first.info.get(k) {
          Some(v) => new_rec = new_rec.info_set(k, v)
          None => ()
        }
      }
      out.push(new_rec)
    }
  }
  VcfFile::{
    meta_headers: vcf.meta_headers,
    fixed_header: vcf.fixed_header,
    records: out,
  }
}

// ---------------------------------------------------------------------------
// 主入口 3：写入/覆盖 ID 列（从注释源取值 — 与 bcftools_ops.mbt 的
//           vcf_annotate_set_id(file, id : String) 区分）
// ---------------------------------------------------------------------------

///|
/// Set (or append to) the `ID` column of each record using the value from
/// the annotation source. Mirrors `bcftools annotate --set-id +'value'`.
///
/// `append` — when true and the record already has IDs, the anno value is
/// appended as an additional ID component; when false the old IDs are
/// replaced.
pub fn vcf_annotate_id_from_source(
  vcf : VcfFile,
  source : VcfAnnotSource,
  match_mode~ : AnnotMatchMode = Overlap,
  append~ : Bool = false,
) -> VcfFile {
  let out : Array[VcfRecord] = Array(capacity=vcf.records.length())
  for rec in vcf.records {
    let val = resolve_annot_value(source, rec, match_mode)
    match val {
      Some(v) => {
        let new_id = if append && rec.id.length() > 0 {
          let x = rec.id.copy()
          x.push(v)
          x
        } else {
          [v]
        }
        out.push(VcfRecord::{
          chrom: rec.chrom,
          pos: rec.pos,
          id: new_id,
          ref_seq: rec.ref_seq,
          alt: rec.alt.copy(),
          qual: rec.qual,
          filter: rec.filter.copy(),
          info: rec.info.copy(),
          format:
            match rec.format {
              Some(f) => Some(f.copy())
              None => None
            },
          samples: rec.samples.copy(),
        })
      }
      None => out.push(rec)
    }
  }
  VcfFile::{
    meta_headers: vcf.meta_headers,
    fixed_header: vcf.fixed_header,
    records: out,
  }
}

// ---------------------------------------------------------------------------
// 主入口 4：通用 INFO 注释入口（接受 VcfAnnotSource 枚举）
// ---------------------------------------------------------------------------

///|
/// Write `INFO/dest_key` on every record whose `(chrom, pos)` matches the
/// annotation source. Unmatched records are emitted unchanged.
pub fn vcf_annotate_info_from_source(
  vcf : VcfFile,
  source : VcfAnnotSource,
  dest_key : String,
  match_mode~ : AnnotMatchMode = Overlap,
) -> VcfFile {
  let out : Array[VcfRecord] = Array(capacity=vcf.records.length())
  for rec in vcf.records {
    match resolve_annot_value(source, rec, match_mode) {
      Some(v) => out.push(rec.info_set(dest_key, v))
      None => out.push(rec)
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
/// Count records that received an annotation and those that didn't, by
/// comparing INFO keys before/after. Returns `(annotated, unchanged)`.
pub fn vcf_annotate_summary(
  before : VcfFile,
  after : VcfFile,
  dest_key : String,
) -> (Int32, Int32) {
  let mut annotated = 0
  let mut unchanged = 0
  let n = before.records.length()
  let mut i = 0
  while i < n {
    let old_has = before.records[i].info.get(dest_key) is Some(_)
    let new_has = after.records[i].info.get(dest_key) is Some(_)
    if (!old_has && new_has) || (old_has && new_has) {
      annotated = annotated + 1
    } else {
      unchanged = unchanged + 1
    }
    i = i + 1
  }
  (annotated, unchanged)
}

///|
/// Fraction of records that now carry `dest_key` after annotation.
pub fn vcf_annotate_rate(after : VcfFile, dest_key : String) -> Double {
  if after.records.length() == 0 { return 0.0 }
  let mut hit = 0
  for r in after.records {
    if r.info.get(dest_key) is Some(_) { hit = hit + 1 }
  }
  hit.to_double() / after.records.length().to_double()
}
